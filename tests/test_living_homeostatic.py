"""The homeostatic-pressure growth law in the living user element (abaqus_living/living_uel.f90, PROPS(22) p_ref and
PROPS(23) form; abaqus_living/homeostatic_driver.f90, test_tangent.f90), in the style of tests/test_uel_homeostatic.py:
the default is unchanged (law off and law on with p = 0 give the same growth), free growth is Eq. 36 of Klempt et al.
2024, the fully constrained cube follows the ODE alpha_dot = k_alpha g(p(alpha)) (forward Euler for the explicit gate,
backward Euler for the implicit one, first order against scipy), growth stops at p_h, and the tangent with the gate
active matches central differences of the residual.

Needs gfortran (builds into a temporary directory) and scipy; skipped without them. About 5 s.
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "abaqus_living"
ECO = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"
pytest.importorskip("scipy")
pytestmark = pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")

E, NU, FV, K_ALPHA = 10.0, 0.49, 1e-3, 1e-3            # Klempt et al. 2024 Table 2
DPDXI = E * (1 + FV) / (1 - 2 * NU)                    # fully constrained, small strain: p = E (phi^2 + f) xi / (1 - 2 nu)
P_REF = 0.25                                           # p_h; xi_h = p_h / DPDXI = 5.0e-4, below the free growth over T = 1.1


@pytest.fixture(scope="module")
def exe(tmp_path_factory):
    d = tmp_path_factory.mktemp("living_uel")
    ff = ["-O2", "-ffree-line-length-none", "-J", str(d), "-I", str(d)]
    subprocess.run(["gfortran", "-O2", "-J", str(d), "-c", str(ECO), "-o", str(d / "eco.o")], check=True)
    subprocess.run(["gfortran", "-cpp", *ff, "-c", str(SRC / "living_uel.f90"), "-o", str(d / "uel.o")], check=True)
    subprocess.run(["gfortran", *ff, "-c", str(SRC / "mms_src.f90"), "-o", str(d / "mms.o")], check=True)
    objs = [str(d / o) for o in ("uel.o", "mms.o", "eco.o")]
    subprocess.run(["gfortran", *ff, str(SRC / "homeostatic_driver.f90"), *objs, "-o", str(d / "hd")], check=True)
    subprocess.run(["gfortran", *ff, str(SRC / "test_tangent.f90"), *objs, "-o", str(d / "tt")], check=True)
    return d


def run(exe, support, p_ref, form, steps, T, raw=False):
    out = subprocess.run([str(exe / "hd"), support, repr(p_ref), str(form), str(steps), repr(T)],
                         capture_output=True, text=True, check=True).stdout
    if raw:
        return out
    a = np.array([[float(v) for v in line.split()] for line in out.strip().splitlines()])
    return a[:, 2], a[:, 3]                                       # alpha - 1, p (Gauss-point means) after every increment


def p_closed(xi):
    return DPDXI * xi


def gate(p, p_ref, form):
    if p_ref <= 0:
        return 1.0
    if form == 2:
        return 1.0 / (1.0 + max(p, 0.0) / p_ref)
    ph = p_ref * (1.0 + FV) if form == 1 else p_ref
    return min(1.0, max(0.0, 1.0 - p / ph))


def ode(T, p_ref, form):
    from scipy.integrate import solve_ivp
    s = solve_ivp(lambda t, x: [K_ALPHA * gate(p_closed(x[0]), p_ref, form)], (0.0, T), [0.0], method="DOP853",
                  rtol=1e-12, atol=1e-16)
    return s.y[0, -1]


def euler(T, steps, p_ref, form):
    x, dt, out = 0.0, T / steps, []
    for _ in range(steps):
        x = x + dt * K_ALPHA * gate(p_closed(x), p_ref, form); out.append(x)
    return np.array(out)


def backward_euler(T, steps, p_ref, form):
    x, dt, out = 0.0, T / steps, []
    for _ in range(steps):
        lo, hi = x, x + dt * K_ALPHA
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid - x - dt * K_ALPHA * gate(p_closed(mid), p_ref, form) > 0:
                hi = mid
            else:
                lo = mid
        x = 0.5 * (lo + hi); out.append(x)
    return np.array(out)


def test_default_unchanged_and_free_growth_is_eq36(exe):
    """law off (p_ref = 0) and the implicit gate with p_ref < 0 print the same bytes; free growth (p = 0) gives
    alpha - 1 = k_alpha t with the law off and on, every form"""
    assert run(exe, "all", 0.0, 0, 11, 1.1, raw=True) == run(exe, "all", -1.0, 10, 11, 1.1, raw=True)
    xi0, p0 = run(exe, "free", 0.0, 0, 11, 1.1)
    assert np.abs(p0).max() < 1e-14 and xi0[-1] == pytest.approx(K_ALPHA * 1.1, rel=1e-13)
    for form in (0, 1, 2, 10, 11, 12):
        xi, p = run(exe, "free", P_REF, form, 11, 1.1)
        assert np.abs(p).max() < 1e-12                   # the Newton residual (1e-13) leaves p of 1e-14 Pa
        assert np.abs(xi - xi0).max() < 1e-15            # g = 1 up to the round-off of p


def test_constrained_explicit_gate_is_forward_euler_first_order(exe):
    xi_ode = ode(1.1, P_REF, 0)
    errs = []
    for steps in (22, 44):
        xi, p = run(exe, "all", P_REF, 0, steps, 1.1)
        assert np.abs(xi - euler(1.1, steps, P_REF, 0)).max() < 1e-12 * xi[-1]
        assert p[-1] == pytest.approx(p_closed(xi[-1]), rel=1e-10)
        errs.append(abs(xi[-1] - xi_ode))
    assert 1.8 < errs[0] / errs[1] < 2.2 and errs[1] < 1e-5 and xi_ode < K_ALPHA * 1.1
    for form in (1, 2):
        xi, p = run(exe, "all", P_REF, form, 22, 1.1)
        assert np.abs(xi - euler(1.1, 22, P_REF, form)).max() < 1e-12 * xi[-1]


def test_constrained_implicit_gate_is_backward_euler_first_order(exe):
    xi_ode = ode(1.1, P_REF, 0)
    errs = []
    for steps in (22, 44):
        xi, p = run(exe, "all", P_REF, 10, steps, 1.1)
        assert np.abs(xi - backward_euler(1.1, steps, P_REF, 0)).max() < 1e-12 * xi[-1]
        errs.append(abs(xi[-1] - xi_ode))
    assert 1.8 < errs[0] / errs[1] < 2.2 and errs[1] < 1e-5
    for form in (1, 2):
        xi, p = run(exe, "all", P_REF, 10 + form, 22, 1.1)
        assert np.abs(xi - backward_euler(1.1, 22, P_REF, form)).max() < 1e-12 * xi[-1]


def test_growth_stops_at_the_homeostatic_pressure(exe):
    xi_h = P_REF / DPDXI
    xi, p = run(exe, "all", P_REF, 0, 40, 4.0)                   # explicit: overshoot at most one increment of growth
    assert P_REF * (1 - 1e-3) < p[-1] <= P_REF + 0.1 * K_ALPHA * DPDXI
    assert xi[-1] == pytest.approx(xi_h, rel=1e-3) and xi[-1] - xi[-2] < 1e-3 * xi[0]
    xi, p = run(exe, "all", P_REF, 10, 40, 4.0)                  # implicit: from below, no overshoot (geometric approach)
    assert p.max() <= P_REF * (1 + 1e-12) and p[-1] > P_REF * (1 - 1e-3)
    assert xi[-1] == pytest.approx(xi_h, rel=1e-3) and xi[-1] - xi[-2] < 1e-3 * xi[0]
    xi2, p2 = run(exe, "all", P_REF, 2, 40, 4.0)                 # form 2 has no homeostatic value
    assert p2[-1] > P_REF and xi2[-1] - xi2[-2] > 0.1 * xi2[0]


@pytest.mark.parametrize("mode", ["gate", "gate1", "gate2", "igate", "igate1", "igate2", "supg igate"])
def test_tangent_with_the_gate_active(exe, mode):
    out = subprocess.run([str(exe / "tt"), *mode.split()], capture_output=True, text=True, check=True).stdout
    lines = out.strip().splitlines()
    err = float(lines[0].split("=")[1].split(",")[0])
    g = np.array([float(v) for v in lines[-2].split(":")[1].split()])
    changed = float(lines[-1].split(":")[1])
    assert err < 1e-8
    assert ((g > 0) & (g < 1)).sum() >= 4                        # the gate is active at most Gauss points
    assert changed > (1e-3 if "igate" in mode else 1e-7)         # and the tangent differs from the law-off one
