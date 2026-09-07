"""End-to-end test of the C shim's ecology path: Fortran-side ABI -> socket ->
Python 0D Hamilton ODE step. Mirrors test_coupling_shim.py for the material
bridge, but exercises biofilm_ecology_eval() (biofilm_py_eval.c) via
test_shim_ecology_main.c instead of biofilm_py_eval().

Requires a C compiler; skipped automatically if absent.
"""
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
_COUP = _ROOT / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))

import ecology_jax                                 # noqa: E402
import material_server as ms                       # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402

_CC = shutil.which("cc") or shutil.which("gcc")
_WS2 = ["-lws2_32"] if sys.platform == "win32" else []
pytestmark = pytest.mark.skipif(
    _CC is None or not (_COUP / "biofilm_py_eval.c").exists()
    or not (_COUP / "test_shim_ecology_main.c").exists(),
    reason="C compiler or shim source unavailable")

G0 = [0.12, 0.10, 0.05, 0.02, 0.0, 0.71, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0]
DT_H = 1.0e-5


@pytest.fixture(scope="module")
def shim_exe():
    tmp = tempfile.mkdtemp()
    exe = Path(tmp) / "test_shim_ecology"
    subprocess.run([_CC, str(_COUP / "test_shim_ecology_main.c"),
                    str(_COUP / "biofilm_py_eval.c"), "-o", str(exe)] + _WS2,
                   check=True)
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


def _run_shim(exe, host, port, timeout=20):
    stdin = " ".join(f"{v:.17g}" for v in (*G0, *THETA_DEMO.tolist(), DT_H)) + "\n"
    env = dict(os.environ)                        # see test_coupling_shim.py
    env["BIOFILM_PY_HOST"] = host
    env["BIOFILM_PY_PORT"] = str(port)
    return subprocess.run([str(exe)], input=stdin, capture_output=True,
                          text=True, env=env, timeout=timeout)


def test_ecology_shim_roundtrip_matches_python(shim_exe, server):
    host, port = server
    r = _run_shim(shim_exe, host, port)
    assert r.returncode == 0, f"shim failed rc={r.returncode}: {r.stderr}"

    g_new = [float(x) for x in r.stdout.split()]
    assert len(g_new) == 12

    ref = np.asarray(ecology_jax.ecology_step(G0, THETA_DEMO, DT_H))
    assert np.allclose(g_new, ref, rtol=1e-13, atol=1e-18)


def test_ecology_shim_reports_failure_when_server_absent(shim_exe):
    with socket.socket() as s:                    # grab a port, then free it
        s.bind(("127.0.0.1", 0))
        dead_port = s.getsockname()[1]
    r = _run_shim(shim_exe, "127.0.0.1", dead_port, timeout=20)
    assert r.returncode != 0
