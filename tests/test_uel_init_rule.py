"""WP4 (KEIO_PLAN.ja.md): the start of a newly colonised Gauss point in the unified user element
p3/uel/unified_uel.f90 (init_rule 0-3), on a tiny column (4 elements, mechanics fixed) through the
Abaqus entry point UEL (p3/uel/emulate_abaqus.py):

  - rule 0 (the default) gives the numbers of the element before the rules were added (16 Oct 2026);
  - rule 1 starts a new Gauss point from the composition of the nearest colonised Gauss point at the
    end of the previous increment, rule 2 from the phi-weighted mean within 1.5 element sizes, rule 3
    from the point model's steady state; each is checked against the same start advanced by the
    driver entry point uel_advance, computed from the neighbours in Python.

Needs gfortran and jax (the point-model case file); skipped without them.
"""
import ctypes as ct
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
sys.path[:0] = [str(ROOT), str(ROOT / "JAXFEM"), str(UEL)]
pytest.importorskip("jax")

N, L, H, NINC = 4, 0.4, 0.1, 12               # element 3 starts its first Gauss points in increment 7 (c_ref = 0) or 8 (c_ref = 1)
GPC = np.array([-1, 1]) / np.sqrt(3)


@pytest.fixture(scope="module")
def em(tmp_path_factory):
    """emulate_abaqus imported from a scratch copy of p3/uel with a fresh build of the element"""
    gf = shutil.which("gfortran")
    if gf is None:
        pytest.skip("gfortran not available")
    d = tmp_path_factory.mktemp("uel")
    subprocess.run([gf, "-O2", "-fPIC", "-c", str(ROOT / "ansys_usermat/coupling/ecology_native.f"), "-o", "eco.o"], cwd=d, check=True)
    subprocess.run([gf, "-O2", "-fPIC", "-ffree-line-length-none", "-c", str(UEL / "unified_uel.f90"), "-o", "uel.o"], cwd=d, check=True)
    subprocess.run([gf, "-shared", "-o", "libunified_uel.so", "uel.o", "eco.o"], cwd=d, check=True)
    for f in ("driver.py", "emulate_abaqus.py"):
        shutil.copy(UEL / f, d / f)
    from make_column_inp import eco_file
    eco_file(d / "uel_eco.txt")
    os.environ["BIOFILM_ECO_CASE"] = str(d / "uel_eco.txt")
    sys.path.insert(0, str(d))
    for name in ("driver", "emulate_abaqus"):
        spec = importlib.util.spec_from_file_location(name, d / f"{name}.py")
        mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod; spec.loader.exec_module(mod)
    em = sys.modules["emulate_abaqus"]
    em.drv = sys.modules["driver"]
    em.drv.lib.uel_set_init_rule.argtypes = [ct.c_int, ct.c_double]
    em.drv.lib.uel_steady_state.argtypes = [ct.c_double, em.drv.D_P, ct.POINTER(ct.c_int)]
    return em


def run_column(em, rule, ninc, cref=1.0):
    """the tiny column: seed phi = 0.9 at x <= 0.1, c = 1, mechanics fixed; returns nodes, elems, U and the
    SV after every growth increment (index 0: after INIT).  c_ref = 1: the point model sees the local
    nutrient, so the composition differs from point to point (with c_ref = 0 every colonised point of the
    same age has the same composition and the rules cannot be told apart)"""
    drv = em.drv
    drv.lib.uel_set_init_rule(rule, 0.0)
    nodes, elems, xs = drv.column(N, L); props = em.F(drv.props_base(s=1.0, c_ref=cref)); xn = nodes[:, 0]
    SV = [em.F(np.zeros(em.NSVT)) for _ in elems]
    U0 = np.zeros(drv.NDN * len(nodes)); U = U0.copy()
    U[3::drv.NDN] = np.where(xn <= 0.1 + 1e-12, 0.9, 1e-4); U[4::drv.NDN] = 1.0
    U, SV, _ = em.increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, 1e-8, 1)
    fixed = np.zeros(len(U), bool); fixed[0::drv.NDN] = fixed[1::drv.NDN] = fixed[2::drv.NDN] = True
    fixed[4::drv.NDN] = np.isclose(xn, L)
    hist = [(U.copy(), [s.copy() for s in SV])]
    for k in range(ninc):
        U, SV, _ = em.increment(U.copy(), U, fixed, nodes, elems, SV, props, 0.2375 / 60, k + 1, ttot=1e-8 + k * 0.2375 / 60)
        hist.append((U.copy(), [s.copy() for s in SV]))
    return nodes, elems, props, hist


def comp(SV, e, em):
    return SV[e][:em.NSV * em.NGP].reshape(em.NGP, em.NSV)


def gauss_points(nodes, elems):
    """reference coordinates of the Gauss points, (nel, 8, 3), in the element's Gauss order, and N (8, 8)"""
    s = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
    Nq = []
    for k in range(2):
        for j in range(2):
            for i in range(2):
                g = np.array([GPC[i], GPC[j], GPC[k]])
                Nq.append(0.125 * np.prod(1 + s * g, axis=1))
    Nq = np.array(Nq)                                             # (gp, node)
    return np.einsum("qa,eai->eqi", Nq, nodes[elems]), Nq


def advance_from(em, g, e, nodes, elems, props, Uprev, Ucur, qs):
    """uel_advance (driver entry point, no registry) on element e with the start g at the Gauss points qs"""
    drv = em.drv
    ed = np.array([drv.NDN * n + d for n in elems[e] for d in range(drv.NDN)])
    sv = np.zeros((em.NSV, em.NGP))
    for q, gq in zip(qs, g):
        sv[0:12, q] = gq; sv[13, q] = 1.0
    sv = em.F(sv); ok = ct.c_int(0)
    drv.lib.uel_advance(em.F(Ucur[ed]), em.F(Uprev[ed]), em.F(nodes[elems[e]].T), sv, props, 0.2375 / 60, ct.byref(ok))
    assert ok.value == 1
    return np.asarray(sv)[0:12, qs].T


def new_starts(em, hist, e=2):
    """the first increment k in which Gauss points of element e start, and those points; returns
    k, qs, (U, SV) after k - 1, (U, SV) after k"""
    for k in range(1, len(hist)):
        before = comp(hist[k - 1][1], e, em)[:, 13]; after = comp(hist[k][1], e, em)[:, 13]
        qs = np.where((before < 0.5) & (after > 0.5))[0]
        if len(qs):
            assert len(qs) == 4, qs
            return k, qs, hist[k - 1], hist[k]
    raise AssertionError("element %d never started" % (e + 1))


def test_default_unchanged(em):
    """rule 0: the numbers of the element before WP4 (original build, 16 Oct 2026), to 1e-12"""
    nodes, elems, props, hist = run_column(em, 0, 10, cref=0.0)
    U, SV = hist[-1]
    phi = U[3::em.NDN].reshape(N + 1, 4)[:, 0]; c = U[4::em.NDN].reshape(N + 1, 4)[:, 0]
    assert np.allclose(phi, [0.7962488877106484, 0.6685861198372602, 0.09620625166514436, -0.015330426584631214,
                             0.004036566909464562], rtol=0, atol=1e-12)
    assert np.allclose(c, [0.9818473142528422, 0.9841592114612876, 0.9897722397062748, 0.9953224663875951, 1.0], rtol=0, atol=1e-12)
    g3 = comp(SV, 2, em)[0]
    assert np.allclose(g3[:12], [0.6537333574481498, 0.0515785781931356, 0.12891681189725754, 0.047198937307942494,
                                 0.03703634247847215, 0.08153597267504253, 0.79391600698931, 0.250752677389439,
                                 0.9180455644438681, 0.041852128830940444, 0.04532407687265921, 0.9320136330210733], rtol=0, atol=1e-12)
    assert abs(g3[12] - 0.0011570860508330805) < 1e-15
    codes = np.concatenate([comp(SV, e, em)[:, 13] for e in range(N)])
    assert set(np.rint(codes).astype(int)) <= {0, 1}             # started flag 1 only: the default start everywhere
    assert abs(np.array(SV).sum() - 99.62505261901687) < 1e-9 and abs(U.sum() - 26.003394525383545) < 1e-12


def test_seed_starts_from_default_under_every_rule(em):
    """the seed has no colonised neighbour at its start: rules 1-3 fall back to the default state there (code 1)"""
    ref = comp(run_column(em, 0, 1)[3][1][1], 0, em)
    for rule in (1, 2):
        sv = comp(run_column(em, rule, 1)[3][1][1], 0, em)
        assert np.array_equal(sv[:, 13], np.ones(8)) and np.array_equal(sv[:, :12], ref[:, :12])


def test_rule_nearest(em):
    nodes, elems, props, hist = run_column(em, 1, NINC)
    k, qs, (Uprev, SVprev), (Ucur, SVcur) = new_starts(em, hist)
    X, _ = gauss_points(nodes, elems)
    # colonised Gauss points at the end of the previous increment and their composition
    cols = [(X[e, q], comp(SVprev, e, em)[q, :12]) for e in range(N) for q in range(8) if comp(SVprev, e, em)[q, 13] > 0.5]
    assert len(cols) == 16
    g = []
    for q in qs:
        d = [np.sum((x - X[2, q]) ** 2) for x, _ in cols]
        g.append(cols[int(np.argmin(d))][1])
    expect = advance_from(em, g, 2, nodes, elems, props, Uprev, Ucur, qs)
    got = comp(SVcur, 2, em)
    assert np.array_equal(got[qs, 13], 2 * np.ones(4))
    assert np.allclose(got[qs, :12], expect, rtol=0, atol=1e-14)
    # and it is not the default start
    assert np.abs(got[qs, :12] - comp(run_column(em, 0, NINC)[3][k][1], 2, em)[qs, :12]).max() > 1e-3


def test_rule_weighted_mean(em):
    nodes, elems, props, hist = run_column(em, 2, NINC)
    k, qs, (Uprev, SVprev), (Ucur, SVcur) = new_starts(em, hist)
    X, Nq = gauss_points(nodes, elems)
    phi_reg = Nq @ hist[k - 2][0][3::em.NDN][elems].T                 # phi of the field at the start of the previous increment
    radius = 1.5 * H                                                 # 1.5 times the cube root of the element's volume
    cols = [(X[e, q], max(phi_reg[q, e], 0.0), comp(SVprev, e, em)[q, :12]) for e in range(N) for q in range(8)
            if comp(SVprev, e, em)[q, 13] > 0.5]
    g = []
    for q in qs:
        inr = [(w, gg) for x, w, gg in cols if np.sum((x - X[2, q]) ** 2) <= radius ** 2]
        assert len(inr) > 1
        w = np.array([v for v, _ in inr]); gg = np.array([v for _, v in inr])
        g.append((w[:, None] * gg).sum(0) / w.sum())
    expect = advance_from(em, g, 2, nodes, elems, props, Uprev, Ucur, qs)
    got = comp(SVcur, 2, em)
    assert np.array_equal(got[qs, 13], 3 * np.ones(4))
    assert np.allclose(got[qs, :12], expect, rtol=0, atol=1e-14)
    near = comp(run_column(em, 1, NINC)[3][k][1], 2, em)[qs, :12]
    assert np.abs(got[qs, :12] - near).max() > 1e-5                 # more than one neighbour in the radius, different compositions


def test_rule_steady_state(em):
    drv = em.drv
    g = np.zeros(12, order="F"); ok = ct.c_int(0)
    drv.lib.uel_steady_state(-1.0, g, ct.byref(ok))                 # under c*
    assert ok.value == 1 and abs(g[:6].sum() - 1) < 1e-12 and (g[:11] > 0).all() and (g[:11] < 1).all()
    nodes, elems, props, hist = run_column(em, 3, NINC, cref=0.0)    # c_ref = 0: every point sees c*, one steady state
    k, qs, (Uprev, SVprev), (Ucur, SVcur) = new_starts(em, hist)
    expect = advance_from(em, [g] * 4, 2, nodes, elems, props, Uprev, Ucur, qs)
    got = comp(SVcur, 2, em)
    assert np.array_equal(got[qs, 13], 4 * np.ones(4))
    assert np.allclose(got[qs, :12], expect, rtol=0, atol=1e-14)
    assert np.abs(expect - g).max() < 1e-7                          # a fixed point of the point model
