"""The implicit homeostatic gate of the unified user element (p3/uel/unified_uel.f90, PROPS(20) = form + 10;
p3/uel/homeostatic_test.py): the fully constrained cube is backward Euler on the ODE alpha_dot = k_alpha g(p(alpha))
to round-off, first order in dt, and reaches p_h from below without overshoot; the consistent tangent matches central
differences with and without the mean dilatation; the explicit gate (the default) still gives the stored results.

Needs gfortran (builds libunified_uel.so in p3/uel if it is missing or older than the source) and scipy; skipped without them. About 15 s.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
sys.path.insert(0, str(UEL))
pytest.importorskip("scipy")
pytestmark = pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")


@pytest.fixture(scope="module")
def ht(tmp_path_factory):
    so = UEL / "libunified_uel.so"
    srcs = (UEL / "unified_uel.f90", ROOT / "ansys_usermat" / "coupling" / "ecology_native.f")
    if not so.exists() or max(x.stat().st_mtime for x in srcs) > so.stat().st_mtime:   # rebuild when stale
        eco = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"
        subprocess.run(["gfortran", "-O2", "-fPIC", "-c", str(eco), "-o", "ecology_native.o"], cwd=UEL, check=True)
        subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-line-length-none", "-c", "unified_uel.f90", "-o", "unified_uel.o"],
                       cwd=UEL, check=True)
        subprocess.run(["gfortran", "-shared", "-o", "libunified_uel.so", "unified_uel.o", "ecology_native.o"], cwd=UEL, check=True)
    if not os.environ.get("BIOFILM_ECO_CASE"):
        from make_column_inp import eco_file
        p = tmp_path_factory.mktemp("eco") / "uel_eco.txt"
        eco_file(p)
        os.environ["BIOFILM_ECO_CASE"] = str(p)
    import homeostatic_test
    return homeostatic_test


K_ALPHA, T, P_REF = 1e-3, 1.1, 0.25          # as tests/test_uel_homeostatic.py


def test_constrained_cube_is_backward_euler_first_order(ht):
    """all u = 0, p = (phi^2 + f) K (1 - alpha^-3): the element with the implicit gate equals the scalar backward-Euler
    sequence on the ODE to round-off, never exceeds p_h, and its error against the ODE (scipy) halves with dt"""
    a_ode = ht.ode_constrained(K_ALPHA, T, P_REF, 0)[-1]
    errs = []
    for steps in (11, 22):
        xi, p, _, _ = ht.run(K_ALPHA, T, steps, "all", p_ref=P_REF, form=0, implicit=True)
        be = ht.backward_euler_constrained(K_ALPHA, T, steps, P_REF, 0)
        assert np.abs(xi + 1 - be).max() < 1e-9
        assert p[-1] == pytest.approx(ht.p_closed(1 + xi[-1]), rel=1e-9)
        assert p.max() <= P_REF and np.all(np.diff(xi) > 0)
        errs.append(abs(xi[-1] + 1 - a_ode))
    assert 1.9 < errs[0] / errs[1] < 2.1
    assert errs[1] < 1e-5 and a_ode - 1 < K_ALPHA * T


def test_reaches_the_homeostatic_pressure_without_overshoot(ht):
    """form 0, long run: p approaches p_h from below geometrically and ends within 1e-6 of it, alpha at alpha_h;
    the explicit gate on the same run overshoots p_h (by at most one increment of growth)"""
    xi, p, _, _ = ht.run(K_ALPHA, 16.0, 40, "all", p_ref=P_REF, form=0, implicit=True)
    assert p.max() <= P_REF and 0 <= 1 - p[-1] / P_REF < 1e-6
    assert xi[-1] == pytest.approx(ht.alpha_h(P_REF) - 1, rel=1e-6)
    assert np.all(np.diff(xi) >= 0) and xi[-1] - xi[-2] < 1e-6 * xi[0]
    xe, pe, _, _ = ht.run(K_ALPHA, 16.0, 40, "all", p_ref=P_REF, form=0, implicit=False)
    assert pe.max() > P_REF
    # the mean dilatation couples the eight Gauss points of the element through the volume term: same limit
    xm, pm, _, _ = ht.run(K_ALPHA, 16.0, 40, "all", p_ref=P_REF, form=0, implicit=True, md=True)
    assert pm.max() <= P_REF and 0 <= 1 - pm[-1] / P_REF < 1e-6
    # forms 1 and 2 (hyperbolic, no homeostatic value): backward Euler as well
    for form in (1, 2):
        xi, p, _, _ = ht.run(K_ALPHA, T, 22, "all", p_ref=P_REF, form=form, implicit=True)
        assert np.abs(xi + 1 - ht.backward_euler_constrained(K_ALPHA, T, 22, P_REF, form)).max() < 1e-9


def test_consistent_tangent(ht):
    """AMATRX against central differences on a distorted element with 0 < g < 1 at the end of the increment, with
    and without the mean dilatation and in small strain; the gate changes the tangent by order one"""
    for md, small, tol in ((False, False, 1e-6), (True, False, 1e-6), (False, True, 1e-7)):
        rel, g, changed = ht.tangent_check(p_ref=0.3, form=0, implicit=True, md=md, small=small, h=1e-5)
        assert 0.05 < g < 0.95 and rel < tol and changed > 1.0, (md, small, rel, g)


def test_explicit_default_unchanged(ht):
    """the explicit gate still equals the stored 2x6x8 result of the implant sector (p_h/E = 0.05, free outer face)
    and forward Euler on the ODE; the deck generator keeps PROPS(20) = form unless --gate-implicit"""
    import make_implant_uel_inp as mi
    xi, p, _, _ = ht.run(K_ALPHA, T, 11, "all", p_ref=P_REF, form=0)
    assert np.abs(xi + 1 - ht.euler_constrained(K_ALPHA, T, 11, P_REF, 0)).max() < 1e-9
    ref = json.load(open(UEL / "results" / "implant_2x6x8_ph0.05.json"))
    nodes, elems, layer, ijk = mi.sector(2, 6, 8)
    r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, p_ref=0.05)
    for k in ("alpha_minus_1", "vM_over_E", "p_over_E"):
        assert r["elements"][k] == ref["elements"][k], k
    assert not r["gate_implicit"] and mi.element_props(0.49, 0.05, 0)[19] == 0.0
    assert mi.element_props(0.49, 0.05, 0, gate_implicit=True)[19] == 10.0 and len(mi.element_props(0.49)) == 16
