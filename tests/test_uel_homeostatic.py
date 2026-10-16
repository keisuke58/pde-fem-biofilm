"""The homeostatic-pressure growth law in the unified user element (p3/uel/unified_uel.f90, PROPS(19) p_ref and
PROPS(20) form; p3/uel/homeostatic_test.py): the default is unchanged to the last bit, free growth is Eq. 36 of
Klempt et al. 2024 with the law on, the fully constrained cube follows the ODE alpha_dot = k_alpha g(p(alpha)) with
first order in dt, growth stops at p_h, and the tangent with the gate active matches central differences.

Needs gfortran (builds libunified_uel.so in p3/uel if it is missing) and scipy; skipped without them. About 10 s.
"""
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
    if not so.exists():
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


K_ALPHA, T, P_REF = 1e-3, 1.1, 0.25          # Klempt et al. 2024 Table 2 rate; p_h half the free-growth pressure


def test_default_unchanged(ht):
    """16 properties, and 20 with p_ref = 0, give the same bits; a gate of 1 changes nothing"""
    xi0, p0, _, _ = ht.run(K_ALPHA, T, 11, "all")
    orig = ht.props

    def props20(k, p_ref=0.0, form=0):
        return ht.F(np.r_[orig(k), 1e30, 0.0, 0.0, 0.0])
    ht.props = props20
    try:
        xi1, p1, _, _ = ht.run(K_ALPHA, T, 11, "all")
    finally:
        ht.props = orig
    assert np.array_equal(xi0, xi1) and np.array_equal(p0, p1)
    assert xi0[-1] == pytest.approx(K_ALPHA * T, rel=2e-6)              # Eq. 36 (1e-6: the frozen point model's drift)


def test_free_growth_is_eq36_with_the_law_on(ht):
    """unconstrained cube: p = 0, so every form gives Eq. 36, alpha - 1 = k_alpha t, and the law-off result"""
    ref = ht.run(K_ALPHA, T, 11, "free")[0]
    for form in (0, 1, 2):
        xi, p, _, _ = ht.run(K_ALPHA, T, 11, "free", p_ref=P_REF, form=form)
        assert np.abs(p).max() < 1e-10
        assert xi[-1] == pytest.approx(K_ALPHA * T, rel=2e-6)
        assert np.abs(xi - ref).max() < 1e-14


def test_constrained_cube_follows_the_ode_first_order(ht):
    """all u = 0: p = (phi^2 + f) K (1 - alpha^-3); the element equals forward Euler on the ODE to round-off and
    its error against the ODE (scipy) halves with dt"""
    a_ode = ht.ode_constrained(K_ALPHA, T, P_REF, 0)[-1]
    errs = []
    for steps in (11, 22):
        xi, p, _, _ = ht.run(K_ALPHA, T, steps, "all", p_ref=P_REF, form=0)
        eul = ht.euler_constrained(K_ALPHA, T, steps, P_REF, 0)
        assert np.abs(xi + 1 - eul).max() < 1e-9
        assert p[-1] == pytest.approx(ht.p_closed(1 + xi[-1]), rel=1e-9)
        assert 0 < p[-1] < P_REF                                        # approached from below (dt k dp/dalpha < p_h)
        errs.append(abs(xi[-1] + 1 - a_ode))
    assert 1.9 < errs[0] / errs[1] < 2.1
    assert errs[1] < 1e-5 and a_ode - 1 < K_ALPHA * T                   # the law slowed the growth


def test_growth_stops_at_the_homeostatic_pressure(ht):
    """form 0, T = 4 (free growth would give 4e-3): p -> p_h, alpha -> alpha_h, growth stops"""
    xi, p, _, _ = ht.run(K_ALPHA, 4.0, 40, "all", p_ref=P_REF, form=0)
    dpda = 3 * (1 + ht.FV) * ht.K
    assert P_REF * (1 - 1e-3) < p[-1] <= P_REF + 0.1 * K_ALPHA * dpda   # within the one-increment overshoot bound
    assert xi[-1] == pytest.approx(ht.alpha_h(P_REF) - 1, rel=1e-3)
    assert xi[-1] - xi[-2] < 1e-3 * xi[0]
    # form 2 has no homeostatic value: growth continues
    xi2, p2, _, _ = ht.run(K_ALPHA, 4.0, 40, "all", p_ref=P_REF, form=2)
    assert p2[-1] > P_REF and xi2[-1] - xi2[-2] > 0.1 * xi2[0]


def test_tangent_with_the_gate_active(ht):
    rel, g, changed = ht.tangent_check(p_ref=1.0, form=0)
    assert g == 0.5 and rel < 1e-7 and changed > 1e-3
