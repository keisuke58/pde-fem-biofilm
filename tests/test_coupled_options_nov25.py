"""Options added to coupled_living_consumption.py on 25 Nov 2026 (inputs and comparisons, the equations are unchanged):
- --inherit-radius R: a newly reached cell takes the phi_3D-weighted mean of all started cells within R (radius_mean),
- --living-one: the uptake of c and a uses sum_i phibar_i = 1 (the field model of Klempt et al. 2024 alone),
- --c-init-zero: the nutrient starts at c = 0 outside the held cells and fills the cube through Eq. 35 with dc/dt.
The functions and the field solver are taken from the script with ast, so its command-line set-up does not run."""
import ast
from pathlib import Path

import numpy as np
import pytest

_SRC = Path(__file__).resolve().parents[1] / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"


def _take(names, ns):
    tree = ast.parse(_SRC.read_text(encoding="utf-8"))
    body = [f for f in tree.body if isinstance(f, (ast.FunctionDef, ast.ClassDef)) and f.name in names]
    assert len(body) == len(names), names
    exec(compile(ast.Module(body=body, type_ignores=[]), str(_SRC), "exec"), ns)
    return [ns[n] for n in names]


def _means(radius, h, fast=False):
    pytest.importorskip("scipy")
    ns = {"np": np, "INHERIT_RADIUS": radius, "H": h, "FAST": fast}
    return _take(["radius_mean", "neighbours_mean"], ns)


def _field(shape, seed=0):
    rng = np.random.default_rng(seed)
    w = rng.random(shape) * (rng.random(shape) > 0.4)      # weight phi_3D, 0 where not started
    G = rng.random(shape + (12,))
    return G, w


def test_radius_below_the_cell_size_is_the_six_neighbours():
    rad, nb = _means(0.05, 0.1)
    G, w = _field((6, 5, 7))
    n1, d1 = rad(G, w)
    n2, d2 = nb(G, w)
    assert np.array_equal(n1, n2) and np.array_equal(d1, d2)


def test_radius_of_one_cell_equals_the_six_neighbours():
    rad, nb = _means(0.1, 0.1)                             # the ball of radius h holds the six face neighbours only
    G, w = _field((6, 5, 7), 1)
    n1, d1 = rad(G, w)
    n2, d2 = nb(G, w)
    assert np.allclose(n1, n2, atol=1e-10) and np.allclose(d1, d2, atol=1e-10)


def test_radius_mean_against_a_direct_sum():
    h, R = 0.1, 0.25
    rad, _ = _means(R, h)
    G, w = _field((7, 6, 8), 2)
    num, den = rad(G, w)
    idx = np.argwhere(np.ones(w.shape, bool))
    for p in ([0, 0, 0], [3, 2, 4], [6, 5, 7], [1, 5, 0]):
        d2 = ((idx - np.array(p)) ** 2).sum(1) * h * h
        m = (d2 <= R * R + 1e-12) & (d2 > 0)                # walls: nothing beyond the faces; the cell itself excluded
        q = idx[m]
        wq = w[q[:, 0], q[:, 1], q[:, 2]]
        assert np.isclose(den[tuple(p)], wq.sum(), atol=1e-10)
        assert np.allclose(num[tuple(p)], (wq[:, None] * G[q[:, 0], q[:, 1], q[:, 2]]).sum(0), atol=1e-10)


def test_radius_mean_on_the_quarter_domain_equals_the_full_domain():
    """--fast: x = 0 and z = 0 (low side of axes 0 and 2) are mirror planes. The quarter result must equal the result on
    the full domain built by mirroring."""
    h, R = 0.1, 0.25
    Gq, wq = _field((5, 6, 4), 3)
    full = lambda u: np.concatenate([np.concatenate([u[::-1, :, ::-1], u[::-1]], 2),
                                     np.concatenate([u[:, :, ::-1], u], 2)], 0)
    Gf, wf = full(Gq), full(wq)
    nq, dq = _means(R, h, fast=True)[0](Gq, wq)
    nf, df = _means(R, h, fast=False)[0](Gf, wf)
    assert np.allclose(dq, df[5:, :, 4:], atol=1e-10) and np.allclose(nq, nf[5:, :, 4:], atol=1e-10)


def _reaction(shape, h, living_one=False, transient=None):
    sp = pytest.importorskip("scipy.sparse")
    import scipy.sparse.linalg as spla
    ns = {"np": np, "sp": sp, "spla": spla, "SHAPE": shape, "H": h, "FAST": False, "AMG": False, "OBSTACLE": None,
          "OBST": None, "LIVING_ONE": living_one, "TRANSIENT_D": transient}
    return _take(["Reaction"], ns)[0]


def _blob(shape):
    g = np.meshgrid(*[np.linspace(-1, 1, n) for n in shape], indexing="ij")
    return np.exp(-4.0 * sum(x * x for x in g))


def test_living_one_ignores_the_living_fraction():
    shape, h = (8, 8, 8), 0.25
    phi = _blob(shape)
    L = 0.2 + 0.5 * _blob(shape)
    on = _reaction(shape, h, living_one=True)(4.0, -1)
    off = _reaction(shape, h, living_one=False)(4.0, -1)
    c_one = on.solve(phi, L)
    assert np.allclose(c_one, on.solve(phi, 0.1 * L), atol=1e-12)            # independent of sum_i phibar_i
    assert np.allclose(c_one, off.solve(phi, np.ones(shape)), atol=1e-12)    # = g phi_3D c of Klempt et al. 2024
    assert (off.solve(phi, L) - c_one).min() > 0.0                           # less living biomass leaves more nutrient


def test_nutrient_from_zero_fills_to_the_quasi_static_field():
    """--c-init-zero: backward Euler steps of (1/d) dc/dt = lap c - (g/d) phi_3D sum_i phibar_i c from c = 0 outside the
    held layer rise monotonically to the solution without dc/dt, within about 1e-7 T* for d = 1e8 mm^2/T*."""
    shape, h = (8, 8, 8), 0.25
    phi, L = _blob(shape), 0.2 + 0.5 * _blob(shape)
    nut = _reaction(shape, h, transient=1.0e8)(4.0, 0)
    nut.fixed = np.zeros(shape, bool)
    nut.fixed[:, 0, :] = True                                                # the layer where c is held
    static = nut.solve(phi, L)
    c, t, prev_mean = np.where(nut.fixed, 1.0, 0.0), 0.0, 0.0
    at_1e7 = None
    for tn in np.logspace(-11, -5, 49):
        c = nut.solve(phi, L, u_prev=c, dtau=tn - t); t = tn
        assert c.mean() >= prev_mean - 1e-12 and c.min() >= -1e-12 and c.max() <= 1 + 1e-12
        prev_mean = c.mean()
        if at_1e7 is None and tn >= 1e-7:
            at_1e7 = np.abs(c - static).max()
    assert at_1e7 < 0.05 and np.abs(c - static).max() < 1e-3
