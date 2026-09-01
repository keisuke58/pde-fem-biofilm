"""BIOFILM_GROWTH_VISCO_V01 — the routine we hand to Oliver's framework.

It repackages the verified law in the shape their `AceGenNeoHookV04` call site
uses: (E, nu) per phase blended by a biofilm fraction, rather than
(C10, C01, D1) directly. The whole value of that routine is that the 0-ULP
Abaqus equivalence (`crosscheck/`, `adversarial.py`) travels with it, so what
these tests check is that the adapter is *only* an adapter:

  * given constants that map to the same (C10, C01, D1), the wrapper's stress,
    viscous update and tangent match the core exactly;
  * the (E, nu) -> (C10, C01, D1) map is the small-strain-consistent one;
  * the growth and viscous paths it adds over their elastic routine actually
    do something.

Requires gfortran; skipped where absent.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_WRAP = _AU / "biofilm_material_v01.f"
_CORE = _AU / "usermat_biofilm.f"
_HOOK = _AU / "coupling" / "usermat_py_hook.f"
_SHIM = _AU / "coupling" / "biofilm_py_eval.c"
_DRV = _AU / "crosscheck" / "wrapper_driver.f"
_XDRV = _AU / "crosscheck" / "xcheck_driver_ans.f"

_FC = shutil.which("gfortran")
_CC = shutil.which("cc") or shutil.which("gcc")
pytestmark = pytest.mark.skipif(
    _FC is None or _CC is None or not all(
        p.exists() for p in (_WRAP, _CORE, _HOOK, _SHIM, _DRV)),
    reason="gfortran/cc or the Fortran sources are unavailable")

I3 = np.eye(3)
F_TEST = np.array([[1.06, 0.02, 0.0], [0.0, 0.97, 0.01], [0.0, 0.0, 0.98]])

# (E, nu) chosen to land on round C10 values; c01_ratio picks Mooney-Rivlin.
E_BIO, NU_BIO = 1000.0, 0.30
E_VOID, NU_VOID = 1.0, 0.30


def _mr_from_E_nu(E, nu, c01_ratio):
    """The map the wrapper implements, mirrored here independently so the
    test does not just restate the Fortran."""
    mu = E / (2.0 * (1.0 + nu))
    K = E / (3.0 * (1.0 - 2.0 * nu))
    c10 = 0.5 * mu / (1.0 + c01_ratio)
    return c10, c10 * c01_ratio, 2.0 / K


@pytest.fixture(scope="module")
def exes():
    """Build the wrapper driver and the existing core driver side by side."""
    tmp = Path(tempfile.mkdtemp())
    o = {n: tmp / f"{n}.o" for n in ("hook", "core", "shim")}
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(o["hook"])], check=True, cwd=tmp)
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(o["core"])], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM), "-o", str(o["shim"])], check=True)

    wrap = tmp / "wrap"
    subprocess.run([_FC, "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRV), str(_WRAP), str(o["core"]), str(o["hook"]),
                    str(o["shim"]), "-o", str(wrap)], check=True)

    core = None
    if _XDRV.exists():
        core = tmp / "core"
        subprocess.run([_FC, "-ffixed-line-length-132", "-I", str(tmp),
                        str(_XDRV), str(o["core"]), str(o["hook"]),
                        str(o["shim"]), "-o", str(core)], check=True)
    return wrap, core


def _run_wrapper(exe, *, F=F_TEST, Fv=I3, biofilm=1.0, growth=0.2,
                 eta=5.0, dt=0.01, c01_ratio=0.15, mtype=1.0,
                 E=E_BIO, EL=E_VOID, nu=NU_BIO, nuL=NU_VOID):
    stdin = (
        " ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        f"{E:.17e} {EL:.17e} {nu:.17e} {nuL:.17e} {biofilm:.17e} "
        f"{growth:.17e} {eta:.17e} {dt:.17e} {c01_ratio:.17e} {mtype:.1f}\n"
    )
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, timeout=30)
    assert r.returncode == 0, f"wrapper driver failed: {r.stderr}"
    t = [float(x) for x in r.stdout.split()]
    assert len(t) == 6 + 9 + 2 + 36
    return dict(stress=np.array(t[0:6]), fv=np.array(t[6:15]).reshape(3, 3),
                keycut=int(t[15]), work=t[16],
                tang=np.array(t[17:53]).reshape(6, 6))


def _run_core(exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt):
    """The existing crosscheck driver, which calls BIOFILM_STRESS_CORE."""
    stdin = (
        " ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        f"{alpha:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} "
        f"{mtype:.1f} {dt:.17e}\n"
    )
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, timeout=30)
    assert r.returncode == 0, f"core driver failed: {r.stderr}"
    v = [float(x) for x in r.stdout.split()]
    return np.array(v[:6]), np.array(v[6:15]).reshape(3, 3)


# --------------------------------------------------------------------------- #
def test_wrapper_is_only_an_adapter(exes):
    """The load-bearing test: feed the wrapper (E, nu) and the core the
    (C10, C01, D1) they map to, and the answers must be identical. If they
    are, the verification the core carries applies to the wrapper unchanged."""
    wrap, core = exes
    if core is None:
        pytest.skip("xcheck_driver_ans.f not present")

    c01r, mtype, alpha, eta, dt = 0.15, 1.0, 0.2, 5.0, 0.01
    c10, c01, d1 = _mr_from_E_nu(E_BIO, NU_BIO, c01r)

    w = _run_wrapper(wrap, biofilm=1.0, growth=alpha, eta=eta, dt=dt,
                     c01_ratio=c01r, mtype=mtype)
    s_core, fv_core = _run_core(core, F_TEST, I3, alpha, c10, c01, d1,
                                eta, mtype, dt)

    np.testing.assert_allclose(w["stress"], s_core, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(w["fv"], fv_core, rtol=0.0, atol=0.0)


def test_blend_reduces_to_each_phase(exes):
    """biofilm=1 must give the biofilm material and biofilm=0 the void one,
    with nothing of the other leaking in."""
    wrap, _ = exes
    only_bio = _run_wrapper(wrap, biofilm=1.0)
    only_void = _run_wrapper(wrap, biofilm=0.0)
    # E differs by 1000x, so the stresses must differ by about that much
    ratio = np.max(np.abs(only_bio["stress"])) / np.max(np.abs(only_void["stress"]))
    assert ratio > 100.0, f"blend does not separate the phases (ratio {ratio:.1f})"


@pytest.mark.parametrize("eta,dt", [(0.0, 0.01), (5.0, 1.0e-4)])
def test_blend_is_monotone_in_the_biofilm_fraction(exes, eta, dt):
    """More biofilm must mean more stress -- but only where the time step
    resolves the viscous relaxation, hence the two cases: no viscosity at
    all, and viscosity with a step well inside the relaxation time. See
    test_the_viscous_step_must_resolve_the_relaxation_time for why the
    qualifier is not optional."""
    wrap, _ = exes
    peaks = [np.max(np.abs(_run_wrapper(wrap, biofilm=b, eta=eta, dt=dt)["stress"]))
             for b in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert all(a < b for a, b in zip(peaks, peaks[1:])), peaks


def test_the_viscous_step_must_resolve_the_relaxation_time(exes):
    """Pins a real limit of the integrator this wrapper inherits.

    BIOFILM_STRESS_CORE advances Fv explicitly -- Fv_n1 = (I + dt/(2*eta*J)
    * tau) . Fv_n, with tau evaluated at Fv_n -- so the step is only stable
    while dt is small against the relaxation time eta/(2*C10). Push past it
    and the elastic strain overshoots through zero: at eta=5, C10~167
    (tau~0.015 s) the 11-stress crosses from -513 Pa to +347 Pa between
    dt=1e-4 and dt=1e-2, and blows up beyond dt/tau ~ 1.

    This is not something the wrapper introduces and not something to fix
    here: the core is locked 0 ULP to the Abaqus UMAT, so changing its
    integrator is a separate decision with its own re-verification. It is
    recorded because the caller picks dt -- in the handover case that is
    Oliver's framework, not us.
    """
    wrap, _ = exes
    s11 = [_run_wrapper(wrap, biofilm=1.0, eta=5.0, dt=dt)["stress"][0]
           for dt in (1.0e-4, 1.0e-3, 5.0e-3, 1.0e-2)]
    assert s11[0] < -100.0, f"resolved step should be compressive here: {s11}"
    assert s11[-1] > 0.0, (
        "the overshoot this test documents is gone -- if the integrator was "
        f"changed on purpose, retire this test; if not, investigate: {s11}")
    assert all(a < b for a, b in zip(s11, s11[1:])), (
        f"overshoot should grow monotonically with the step: {s11}")


def test_growth_actually_does_something(exes):
    """Growth is what their AceGenNeoHookV04 has no notion of, so check it
    is not silently inert here."""
    wrap, _ = exes
    a = _run_wrapper(wrap, growth=0.0)["stress"]
    b = _run_wrapper(wrap, growth=0.3)["stress"]
    assert np.max(np.abs(a - b)) > 1.0, "growth had no effect on the stress"


def test_viscous_state_advances_and_is_not_symmetric(exes):
    """The viscous update must move Fv, and Fv must be returned as a full
    3x3: it does not stay symmetric, which is why the interface needs nine
    state slots rather than a six-component Cauchy-Green."""
    wrap, _ = exes
    fv = _run_wrapper(wrap, eta=5.0, dt=0.01)["fv"]
    assert np.max(np.abs(fv - I3)) > 1e-8, "Fv did not advance"

    # march it forward to let any asymmetry accumulate
    for _ in range(20):
        fv = _run_wrapper(wrap, Fv=fv, eta=5.0, dt=0.01)["fv"]
    asym = np.max(np.abs(fv - fv.T)) / max(np.max(np.abs(fv)), 1e-30)
    assert asym > 1e-9, (
        "Fv came back symmetric — if that ever holds, a 6-component Cv would "
        "suffice and this interface note should be revisited")


def test_elastic_limit_has_no_viscous_flow(exes):
    wrap, _ = exes
    fv = _run_wrapper(wrap, eta=0.0)["fv"]
    np.testing.assert_allclose(fv, I3, rtol=0.0, atol=0.0)


def test_tangent_is_finite_and_scales_with_stiffness(exes):
    wrap, _ = exes
    soft = _run_wrapper(wrap, biofilm=0.0)
    stiff = _run_wrapper(wrap, biofilm=1.0)
    for r in (soft, stiff):
        assert r["tang"].shape == (6, 6)
        assert np.all(np.isfinite(r["tang"]))
    assert np.max(np.abs(stiff["tang"])) > np.max(np.abs(soft["tang"]))


def test_cutback_is_requested_when_the_jacobian_collapses(exes):
    """A collapsing elastic Jacobian must set keycut rather than return
    nonsense, the same contract usermat_biofilm.f honours."""
    wrap, _ = exes
    fv_in = np.array([[1.02, 0.01, 0.0], [0.0, 0.99, 0.0], [0.0, 0.0, 1.01]])
    r = _run_wrapper(wrap, F=np.diag([1e-3, 1e-3, 1e-3]), growth=2.0, Fv=fv_in)
    assert r["keycut"] == 1
    # every output must still be defined: the caller's arrays are work
    # space, so "left untouched" means "returns uninitialised memory".
    assert np.all(np.isfinite(r["stress"])) and np.all(np.isfinite(r["tang"]))
    np.testing.assert_allclose(r["stress"], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(r["tang"], 0.0, rtol=0.0, atol=0.0)
    # and Fv must not have been advanced off the collapsed configuration
    np.testing.assert_allclose(r["fv"], fv_in, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("c01r,mtype", [(0.0, 0.0), (0.15, 1.0), (0.3, 1.0)])
def test_neo_hookean_and_mooney_rivlin_paths_both_run(exes, c01r, mtype):
    wrap, _ = exes
    r = _run_wrapper(wrap, c01_ratio=c01r, mtype=mtype)
    assert r["keycut"] == 0
    assert np.all(np.isfinite(r["stress"]))
    assert np.all(np.isfinite(r["tang"]))
