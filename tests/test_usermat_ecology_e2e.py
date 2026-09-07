"""End-to-end test of the ecology (0D Hamilton ODE) branch through the REAL
`usermat()` entry point -- the counterpart of test_usermat_kusepy_e2e.py for
coupling/README.md's next-steps #4: a live per-Gauss-point ecology step
driving the growth input alpha, instead of a precomputed field.

Compiles usermat_ecology_e2e_driver.f + usermat_py_hook.f +
../usermat_biofilm.f + biofilm_py_eval.c and drives the actual `usermat()`
subroutine with prop(8)=kUseEcology=1, checking:
  - ustatev(15:26) after one increment matches ecology_jax.ecology_step
    (the same Newton-per-step 0D Hamilton integrator) bit-for-bit through
    the wire,
  - ustatev(10) (alpha) advances by exactly dTime*k_alpha*phi_tot(g_new),
  - an all-zero incoming ecology state is seeded by INIT_ECO_IF_ZERO to
    match ecology_jax.default_initial_state() before stepping,
  - kUseEcology=0 leaves ustatev(10) and ustatev(15:26) untouched, and
  - a dead server falls back cleanly (no crash, alpha/state unchanged, no
    cut-back) rather than aborting the solve.

Requires gfortran, a C compiler, and the ecology driver source; skipped
automatically if any are unavailable.
"""
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

_WS2 = ["-lws2_32"] if sys.platform == "win32" else []

import numpy as np
import pytest

pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_COUP = _AU / "coupling"
_CORE = _AU / "usermat_biofilm.f"
_HOOK = _COUP / "usermat_py_hook.f"
_DRIVER = _COUP / "usermat_ecology_e2e_driver.f"
_SHIM_C = _COUP / "biofilm_py_eval.c"

sys.path.insert(0, str(_COUP))
import ecology_jax                                # noqa: E402
import material_server as ms                       # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402

_CC = shutil.which("cc") or shutil.which("gcc")
_FC = shutil.which("gfortran")
pytestmark = pytest.mark.skipif(
    _FC is None or _CC is None or not all(
        p.exists() for p in (_CORE, _HOOK, _DRIVER, _SHIM_C)),
    reason="gfortran/cc or the usermat/ecology sources are unavailable")

G0 = np.array([0.12, 0.10, 0.05, 0.02, 0.0, 0.71, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0])
K_ALPHA = 50.0
DT = 1.0e-5
NSTATEV = 26


@pytest.fixture(scope="module")
def e2e_exe():
    tmp = Path(tempfile.mkdtemp())
    hook_o, core_o, driver_o, shim_o = (tmp / f"{n}.o" for n in
                                        ("hook", "core", "driver", "shim"))
    exe = tmp / "eco_e2e"
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(hook_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_CORE), "-o", str(core_o)], check=True, cwd=tmp)
    subprocess.run(["gfortran", "-c", "-ffixed-line-length-132", "-I", str(tmp),
                    str(_DRIVER), "-o", str(driver_o)], check=True, cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM_C), "-o", str(shim_o)], check=True)
    subprocess.run(["gfortran", "-o", str(exe), str(driver_o), str(hook_o),
                    str(core_o), str(shim_o)] + _WS2, check=True)
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


def _run(exe, alpha, kuseeco, g_old, k_alpha=K_ALPHA, dt=DT,
         host=None, port=None, timeout=20):
    I = np.eye(3)
    stdin = (
        " ".join(f"{I[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        " ".join(f"{I[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
        f"{alpha:.17e} 1.0 0.0 0.01 0.0 0.0 {dt:.17e}\n" +
        f"{kuseeco:.1f} {k_alpha:.17e}\n" +
        " ".join(f"{v:.17e}" for v in THETA_DEMO.tolist()) + "\n" +
        " ".join(f"{v:.17e}" for v in g_old) + "\n"
    )
    env = dict(os.environ)                        # see test_usermat_kusepy_e2e.py
    if host is not None:
        env["BIOFILM_PY_HOST"] = host
        env["BIOFILM_PY_PORT"] = str(port)
    r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True,
                       env=env, timeout=timeout)
    assert r.returncode == 0, f"driver failed rc={r.returncode}: {r.stderr}"
    toks = [float(x) for x in r.stdout.split()]
    assert len(toks) == 6 + NSTATEV + 2
    ustatev = np.array(toks[6:6 + NSTATEV])
    keycut = int(toks[6 + NSTATEV])
    return ustatev, keycut


def test_ecology_step_matches_python_reference(e2e_exe, server):
    host, port = server
    ustatev, keycut = _run(e2e_exe, alpha=0.0, kuseeco=1.0, g_old=G0,
                           host=host, port=port)
    assert keycut == 0

    g_new_ref = np.asarray(ecology_jax.ecology_step(G0, THETA_DEMO, DT))
    np.testing.assert_allclose(ustatev[14:26], g_new_ref, rtol=1e-9, atol=1e-12)

    phi_tot = ecology_jax.living_fraction_total(g_new_ref)
    alpha_ref = 0.0 + DT * K_ALPHA * phi_tot
    assert ustatev[9] == pytest.approx(alpha_ref, rel=1e-9, abs=1e-12)


def test_ecology_seeds_default_state_when_zero(e2e_exe, server):
    host, port = server
    zero = np.zeros(12)
    ustatev, keycut = _run(e2e_exe, alpha=0.0, kuseeco=1.0, g_old=zero,
                           host=host, port=port)
    assert keycut == 0

    g_default = np.asarray(ecology_jax.default_initial_state())
    g_new_ref = np.asarray(ecology_jax.ecology_step(g_default, THETA_DEMO, DT))
    np.testing.assert_allclose(ustatev[14:26], g_new_ref, rtol=1e-6, atol=1e-9)


def test_ecology_disabled_leaves_alpha_and_state_untouched(e2e_exe):
    ustatev, keycut = _run(e2e_exe, alpha=0.3, kuseeco=0.0, g_old=G0)
    assert keycut == 0
    assert ustatev[9] == pytest.approx(0.3)
    np.testing.assert_allclose(ustatev[14:26], G0)


def test_ecology_falls_back_when_server_unreachable(e2e_exe):
    with socket.socket() as s:                    # grab a port, then free it
        s.bind(("127.0.0.1", 0))
        dead_port = s.getsockname()[1]
    ustatev, keycut = _run(e2e_exe, alpha=0.3, kuseeco=1.0, g_old=G0,
                           host="127.0.0.1", port=dead_port)
    assert keycut == 0
    assert ustatev[9] == pytest.approx(0.3)         # alpha stays at input
    np.testing.assert_allclose(ustatev[14:26], G0)  # state stays at input
