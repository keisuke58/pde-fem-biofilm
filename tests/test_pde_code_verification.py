"""Code verification for the PDE discretisations, against exact solutions.

These pin `JAXFEM/pde_code_verification.py`. See `PDE_VERIFICATION_FINDINGS.md`
for what they found and what it means.

The point of the file is that every one of these targets is *exact* -- a cosh
profile, a parabola, a Neumann eigenfunction, a conservation law -- so a
failure here is a defect in the discretisation and not a disagreement with
another implementation.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "JAXFEM"))

import jax.numpy as jnp                                            # noqa: E402

from JAXFEM import pde_code_verification as v                      # noqa: E402
from JAXFEM import core_hamilton_1d_nutrient as p1                 # noqa: E402
from JAXFEM import core_hamilton_2d_nutrient as p2                 # noqa: E402

GRIDS = (21, 41, 81)


def _order(rows):
    return [np.log(a[2] / b[2]) / np.log(a[1] / b[1])
            for a, b in zip(rows[:-1], rows[1:])]


# ---------------------------------------------------------------------------
# the variant really is the committed routine plus one coefficient
# ---------------------------------------------------------------------------

def test_the_variant_at_coefficient_one_is_the_committed_step_bit_for_bit():
    """Without this the comparison would be two opinions, not one change."""
    rng = np.random.default_rng(0)
    for N in (21, 41):
        params, _, _ = v._params(N, 1.0, 4.0e6, 1.0e6)
        for _ in range(5):
            c = jnp.asarray(rng.uniform(0.1, 1.0, N))
            a = np.asarray(p1.nutrient_step(c, 1.0, params))
            b = np.asarray(v.nutrient_step_asis(c, 1.0, params))
            assert np.array_equal(a, b)


# ---------------------------------------------------------------------------
# 1D: the committed Neumann row costs one order
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def steady_1d():
    return {name: v.study(step, grids=GRIDS)
            for name, step in (("as-is", None),
                               ("mirror", v.nutrient_step_neumann2))}


@pytest.mark.parametrize("case", list(v.CASES))
def test_the_committed_solver_is_only_first_order(steady_1d, case):
    """Measured against cosh (first-order kinetics) and the parabola that is
    Klempt 2024 Eq. 35's own steady state. Both come out ~1.0, not ~2.0."""
    for p in _order(steady_1d["as-is"][case])[-1:]:
        assert 0.8 < p < 1.3, f"{case}: observed order {p}"


@pytest.mark.parametrize("case", list(v.CASES))
def test_reflecting_about_the_node_restores_second_order(steady_1d, case):
    rows = steady_1d["mirror"][case]
    if case.startswith("zeroth"):
        # a second-order central difference is EXACT on a quadratic, so there
        # is no order to measure -- the error sits at the asymptotic floor.
        assert max(r[2] for r in rows) < 1.0e-7
    else:
        for p in _order(rows)[-1:]:
            assert 1.8 < p < 2.2, f"{case}: observed order {p}"


def test_the_committed_error_is_worst_at_the_neumann_wall(steady_1d):
    """Localises the defect: the sup norm is attained at x=0 every time."""
    for case, rows in steady_1d["as-is"].items():
        for N, h, sup, at_wall, interior, _res in rows:
            assert at_wall == pytest.approx(sup), (case, N)
            assert interior < sup


def test_the_zeroth_order_floor_is_the_asymptotic_not_the_discretisation():
    """err / k_monod is flat, so the discretisation contributes nothing."""
    ratios = []
    for k in (1.0e-6, 1.0e-8, 1.0e-10):
        x, c, _ = v.steady_state(41, 1.0, 1.0, k, 1.0,
                                 step=v.nutrient_step_neumann2)
        ce = v.exact_zeroth_order(x, 1.0, 1.0, 1.0)
        ratios.append(float(np.max(np.abs(c - ce))) / k)
    assert max(ratios) / min(ratios) < 1.1, ratios


# ---------------------------------------------------------------------------
# 2D: the same treatment on four walls, and there it does not converge at all
# ---------------------------------------------------------------------------

def test_the_2d_neumann_laplacian_does_not_converge_at_the_walls():
    """The wall rows return half the correct second difference, so the error
    is half the exact eigenvalue 2(pi/L)^2 = 19.74 and is h-independent."""
    rows = v.study_2d_operator(p2.laplacian_2d_neumann, grids=GRIDS)
    errs = [r[2] for r in rows]
    assert all(abs(e - 0.5 * 2.0 * np.pi ** 2) < 0.05 for e in errs), errs
    assert max(errs) / min(errs) < 1.01           # flat under refinement
    assert all(r[3] < r[2] / 100.0 for r in rows)  # the interior is fine


def test_mirroring_makes_the_2d_laplacian_second_order():
    rows = v.study_2d_operator(v.laplacian_2d_neumann_mirror, grids=GRIDS)
    for p in _order(rows)[-1:]:
        assert 1.8 < p < 2.2, p


def test_the_committed_2d_diffusion_leaks_mass_through_a_zero_flux_wall():
    """The consequence in the units the 2D study reports. `diffusion_step` at
    core_hamilton_2d_nutrient.py:434 uses this operator, and
    `condition_spread_2d.py:173` uses that -- so this leak is on the path the
    committed condition-separation number was computed on."""
    N, D = 41, 1.0
    h = 1.0 / (N - 1)
    dt = 0.2 * h * h / D
    g = np.linspace(0.0, 1.0, N)
    X, Y = np.meshgrid(g, g, indexing="ij")
    u0 = np.exp(-((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 0.01)
    w = np.ones((N, N))                     # finite-volume (trapezoidal) measure
    w[0, :] *= 0.5
    w[-1, :] *= 0.5
    w[:, 0] *= 0.5
    w[:, -1] *= 0.5

    def drift(op):
        u = jnp.asarray(u0)
        m0 = float(np.sum(np.asarray(u) * w))
        for _ in range(2000):
            u = u + dt * D * op(u, h, h)
        return (float(np.sum(np.asarray(u) * w)) - m0) / m0

    assert drift(p2.laplacian_2d_neumann) < -0.01        # loses several percent
    assert abs(drift(v.laplacian_2d_neumann_mirror)) < 1.0e-12
