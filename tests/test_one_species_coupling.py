"""The one-species coupling, end to end at one Gauss point, without ANSYS.

local phi -> BIOFILM_ALPHA_FROM_PHI (Klempt 2024 Eq. 36) ->
BIOFILM_GROWTH_VISCO_V01 (the routine handed to the partner) -> stress.

Everything is held against something that does not share the code: the
growth against the sum Eq. 36 integrates to, the stress against the closed
forms in ansys_usermat/apdl/closed_form_reference.py. Viscosity off and
neo-Hookean throughout -- the Klempt-2024 setting chosen on 2026-10-01.

Requires gfortran and a C compiler; skipped where absent.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
sys.path.insert(0, str(_AU))
sys.path.insert(0, str(_AU / "apdl"))

import one_species_reference as ref          # noqa: E402
import closed_form_reference as cf           # noqa: E402

_GFP = _AU / "growth_from_phi.f"
_WRAP = _AU / "biofilm_material_v01.f"
_CORE = _AU / "usermat_biofilm.f"
_HOOK = _AU / "coupling" / "usermat_py_hook.f"
_SHIM = _AU / "coupling" / "biofilm_py_eval.c"
_DRV = _AU / "crosscheck" / "one_species_driver.f"

_FC = shutil.which("gfortran")
_CC = shutil.which("cc") or shutil.which("gcc")
_WS2 = ["-lws2_32"] if sys.platform == "win32" else []
pytestmark = pytest.mark.skipif(
    _FC is None or _CC is None or not all(
        p.exists() for p in (_GFP, _WRAP, _CORE, _HOOK, _SHIM, _DRV)),
    reason="gfortran/cc or the Fortran sources are unavailable")

E_BIO, NU_BIO = 1000.0, 0.30
E_VOID, NU_VOID = 1.0, 0.30
I3 = np.eye(3)


@pytest.fixture(scope="module")
def exe():
    tmp = Path(tempfile.mkdtemp())
    o = {n: tmp / f"{n}.o" for n in ("hook", "core", "shim")}
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(o["hook"])], check=True, cwd=tmp)
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(o["core"])], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM), "-o", str(o["shim"])],
                   check=True)
    out = tmp / "one_species"
    subprocess.run([_FC, "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRV), str(_GFP), str(_WRAP), str(o["core"]),
                    str(o["hook"]), str(o["shim"]), "-o", str(out)] + _WS2,
                   check=True)
    return out


def _run(exe, phi, *, F=I3, k_alpha=50.0, dt=1.0e-3, biofilm=1.0):
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    stdin = (" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3))
             + "\n" + f"{E_BIO:.17e} {E_VOID:.17e} {NU_BIO:.17e} "
             f"{NU_VOID:.17e} {biofilm:.17e} {k_alpha:.17e} {dt:.17e} "
             f"{phi.size}\n" + "\n".join(f"{p:.17e}" for p in phi) + "\n")
    r = subprocess.run([str(exe)], input=stdin, capture_output=True,
                       text=True, timeout=30)
    assert r.returncode == 0, r.stderr
    t = r.stdout.split()
    n = phi.size
    return (np.array([float(x) for x in t[:n]]),
            np.array([float(x) for x in t[n:n + 6]]), int(t[n + 6]))


def _d1(E=E_BIO, nu=NU_BIO):
    return 2.0 / (E / (3.0 * (1.0 - 2.0 * nu)))


def _c10(E=E_BIO, nu=NU_BIO):
    return 0.5 * E / (2.0 * (1.0 + nu))


# --- growth: Eq. 36 and nothing else ---------------------------------------

def test_growth_is_eq36_for_constant_phi(exe):
    phi, k, dt, n = 0.4, 50.0, 1.0e-3, 25
    alpha, _, _ = _run(exe, [phi] * n, k_alpha=k, dt=dt)
    np.testing.assert_allclose(alpha[-1], k * phi * n * dt, rtol=1e-13)


def test_growth_matches_the_reference_bit_for_bit_on_a_varying_phi(exe):
    """A phi history like a front passing the point: rises, then saturates."""
    t = np.linspace(0.0, 1.0, 40)
    phi = 1.0 / (1.0 + np.exp(-12.0 * (t - 0.4)))
    alpha, _, _ = _run(exe, phi, k_alpha=7.0, dt=2.5e-3)
    np.testing.assert_array_equal(alpha, ref.alpha_history(phi, 7.0, 2.5e-3))


def test_no_biofilm_no_growth(exe):
    alpha, sig, kc = _run(exe, [0.0] * 10)
    assert np.all(alpha == 0.0) and np.all(sig == 0.0) and kc == 0


def test_nothing_is_added_to_eq36(exe):
    """No cap, no clamp, no gate: phi outside [0,1] and alpha past any cap
    the USERMAT carries go straight through, because Eq. 36 says so."""
    alpha, _, _ = _run(exe, [1.2] * 100, k_alpha=50.0, dt=1.0e-3)
    np.testing.assert_allclose(alpha[-1], 50.0 * 1.2 * 0.1, rtol=1e-13)
    alpha, _, _ = _run(exe, [-0.1] * 10, k_alpha=50.0, dt=1.0e-3)
    assert alpha[-1] < 0.0


def test_the_klempt_convention_is_a_constant_shift():
    a = ref.alpha_history([0.3] * 10, 5.0, 0.01)
    np.testing.assert_allclose([ref.alpha_klempt(x) for x in a], 1.0 + a)


# --- stress: through the handed-over routine, against closed form ----------

def test_free_growth_gives_zero_stress(exe):
    """F = (1+alpha) I is pure growth: Fe = I, no stress at all."""
    phi, k, dt, n = 0.5, 50.0, 1.0e-3, 8
    a = k * phi * n * dt
    _, sig, kc = _run(exe, [phi] * n, F=(1.0 + a) * I3, k_alpha=k, dt=dt)
    np.testing.assert_allclose(sig, cf.free_growth_stress(), atol=1e-10)
    assert kc == 0


def test_constrained_growth_matches_closed_form(exe):
    """F = I: all growth is resisted. Compared against the closed form plus
    the documented spherical term (DEVIATOR_SCALING_FINDING.md), exactly as
    tests/test_closed_form_reference.py does for the core itself."""
    phi, k, dt, n = 0.5, 50.0, 1.0e-3, 4
    alpha, sig, kc = _run(exe, [phi] * n, k_alpha=k, dt=dt)
    a = alpha[-1]
    want = cf.constrained_stress(a, _d1()) + cf.spurious_term(a, _c10())
    np.testing.assert_allclose(sig, want, rtol=1e-10)
    assert sig[0] < 0.0 and kc == 0          # compressive, and not refused
