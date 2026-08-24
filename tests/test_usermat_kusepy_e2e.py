"""End-to-end test of the KUSEPY branch through the REAL `usermat()` entry point.

Unlike test_coupling_shim.py (which exercises biofilm_py_eval.c directly) and
test_coupling_vs_fortran.py (which cross-checks the Python material core against
BIOFILM_STRESS_CORE, bypassing usermat() entirely), this test compiles and runs
the actual ANSYS-facing `usermat()` subroutine in usermat_biofilm.f, toggling
`kUsePy` via prop(6). It proves the whole chain end to end:

    usermat() -> biofilm_py_bridge (usermat_py_hook.f) -> biofilm_py_eval.c
              -> socket -> material_server.py -> back through the Abaqus->ANSYS
              Voigt reindex (MAP6) -> stress(6)/dsdePl(6x6)/ustatev

kUsePy=1 (the Python hook, live server) must match kUsePy=0 (the verified
inline Fortran core) to numerical precision for stress and the updated viscous
state, and to finite-difference-noise precision for the consistent tangent
dsdePl (both sides compute it via the same F-perturbation scheme, so the
residual mismatch is floating-point rounding between the Fortran and NumPy
evaluations, not a modelling difference).

Regression guard: earlier drafts of usermat_py_hook.f built `dsde` via a plain
`reshape(d36, [6, 6])`. d36 is a row-major (C-order) flatten of the Python
side's 6x6 Jacobian, but Fortran's RESHAPE fills column-major, so that line
silently returned the TRANSPOSE of the intended matrix -- invisible for
near-symmetric elastic cases but a large, sign-flipping discrepancy for
viscous/Mooney-Rivlin cases. The fix adds an explicit `transpose(...)`; the
dsdePl comparisons below are what would have caught it.

Requires gfortran, a C compiler, and the Fortran core/hook/driver sources;
skipped automatically if any are unavailable.
"""
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_COUP = _AU / "coupling"
_CORE = _AU / "usermat_biofilm.f"
_HOOK = _COUP / "usermat_py_hook.f"
_DRIVER = _COUP / "usermat_endtoend_driver.f"
_SHIM_C = _COUP / "biofilm_py_eval.c"

sys.path.insert(0, str(_COUP))
import material_server as ms                      # noqa: E402

_CC = shutil.which("cc") or shutil.which("gcc")
_FC = shutil.which("gfortran")
pytestmark = pytest.mark.skipif(
    _FC is None or _CC is None or not all(
        p.exists() for p in (_CORE, _HOOK, _DRIVER, _SHIM_C)),
    reason="gfortran/cc or the usermat sources are unavailable")


def _cases():
    """(name, F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt) battery covering
    elastic, viscous, and Mooney-Rivlin regimes, all near-identity so the
    inline core's numerical tangent stays well away from any clamp."""
    I = np.eye(3)
    return [
        ("elastic",
         np.diag([1.05, 0.98, 0.97]), I, 0.2, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0),
        ("elastic+shear",
         I + np.array([[0.05, 0.02, 0.0], [0.0, -0.02, 0.01], [0.0, 0.0, -0.03]]),
         I, 0.2, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0),
        ("viscous",
         np.array([[1.05, 0.02, 0.0], [0.0, 0.98, 0.01], [0.0, 0.0, 0.97]]),
         I, 0.3, 1.0, 0.0, 0.01, 5.0, 0.0, 0.1),
        ("mooney-rivlin",
         np.array([[1.08, 0.0, 0.03], [0.01, 0.95, 0.0], [0.0, 0.02, 1.02]]),
         I, 0.1, 0.8, 0.3, 0.02, 0.0, 1.0, 1.0),
        ("viscous+mooney-rivlin",
         np.array([[1.10, 0.04, 0.02], [0.0, 0.92, 0.03], [0.01, 0.0, 0.99]]),
         I, 0.4, 0.6, 0.2, 0.015, 3.0, 1.0, 0.05),
        ("prior Fv sheared",
         np.diag([1.1, 1.0, 0.9]),
         I + 0.05 * np.outer([1, 0, 0], [0, 1, 0]), 0.2, 1.0, 0.0, 0.01, 2.0, 0.0, 0.2),
    ]


@pytest.fixture(scope="module")
def e2e_exe():
    tmp = Path(tempfile.mkdtemp())
    hook_o, core_o, driver_o, shim_o = (tmp / f"{n}.o" for n in
                                        ("hook", "core", "driver", "shim"))
    exe = tmp / "e2e"
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(hook_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(core_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRIVER), "-o", str(driver_o)], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM_C), "-o", str(shim_o)], check=True)
    subprocess.run(["gfortran", "-o", str(exe), str(driver_o), str(hook_o),
                    str(core_o), str(shim_o)], check=True)
    return exe


@pytest.fixture(scope="module")
def server():
    try:
        srv = ms.socketserver.TCPServer(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv.server_address
    srv.shutdown()
    srv.server_close()


def _run(exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt, kusepy,
         host=None, port=None, timeout=20, state_mat=None):
    """state_mat: None leaves the per-IP material path disabled (prop(7)=0);
    otherwise a (C10, C01, D1, eta) tuple written to ustatev(11:14) with
    prop(7)=1, so the USERMAT takes its constants from state instead of prop."""
    if state_mat is None:
        kstmat, sm = 0.0, (0.0, 0.0, 0.0, 0.0)
    else:
        kstmat, sm = 1.0, state_mat
    stdin = (
        " ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        f"{alpha:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} "
        f"{mtype:.1f} {dt:.17e} {kusepy:.1f}\n" +
        f"{kstmat:.1f} " + " ".join(f"{v:.17e}" for v in sm) + "\n"
    )
    env = {"PATH": "/usr/bin:/bin"}
    if host is not None:
        env["BIOFILM_PY_HOST"] = host
        env["BIOFILM_PY_PORT"] = str(port)
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       env=env, timeout=timeout)
    assert r.returncode == 0, f"driver failed rc={r.returncode}: {r.stderr}"
    toks = [float(x) for x in r.stdout.split()]
    assert len(toks) == 6 + 9 + 2 + 36
    stress = np.array(toks[0:6])
    ustatev = np.array(toks[6:15])
    keycut = int(toks[15])
    dsde = np.array(toks[17:53]).reshape(6, 6)
    return stress, ustatev, keycut, dsde


@pytest.mark.parametrize("case", _cases(), ids=[c[0] for c in _cases()])
def test_kusepy_matches_inline_core(e2e_exe, server, case):
    _, F, Fv, alpha, c10, c01, d1, eta, mtype, dt = case
    host, port = server
    s0, u0, k0, d0 = _run(e2e_exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt, 0.0)
    s1, u1, k1, d1_ = _run(e2e_exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt, 1.0,
                           host=host, port=port)

    assert k0 == 0 and k1 == 0, "unexpected cut-back in a well-conditioned case"
    np.testing.assert_allclose(s1, s0, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(u1, u0, rtol=1e-9, atol=1e-12)
    # dsdePl: both paths use the same F-perturbation scheme (PERT=1e-7), so the
    # residual is Fortran-vs-NumPy floating-point noise, not a modelling gap.
    np.testing.assert_allclose(d1_, d0, rtol=1e-4, atol=1e-4)


def test_kusepy_falls_back_when_server_unreachable():
    """With no server listening, PYOK must come back false and usermat() must
    fall through to the verified inline core -- not crash or return garbage."""
    with socket.socket() as s:                     # grab a port, then free it
        s.bind(("127.0.0.1", 0))
        dead_port = s.getsockname()[1]
    _, F, Fv, alpha, c10, c01, d1, eta, mtype, dt = _cases()[0]

    tmp = Path(tempfile.mkdtemp())
    hook_o, core_o, driver_o, shim_o = (tmp / f"{n}.o" for n in
                                        ("hook", "core", "driver", "shim"))
    exe = tmp / "e2e"
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(hook_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(core_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRIVER), "-o", str(driver_o)], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM_C), "-o", str(shim_o)], check=True)
    subprocess.run(["gfortran", "-o", str(exe), str(driver_o), str(hook_o),
                    str(core_o), str(shim_o)], check=True)

    s0, u0, k0, d0 = _run(exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt, 0.0)
    s1, u1, k1, d1_ = _run(exe, F, Fv, alpha, c10, c01, d1, eta, mtype, dt, 1.0,
                           host="127.0.0.1", port=dead_port)
    assert k1 == 0
    np.testing.assert_allclose(s1, s0, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(u1, u0, rtol=1e-9, atol=1e-12)
