"""End-to-end test of the C shim: Fortran-side ABI → socket → Python material model.

Compiles `biofilm_py_eval.c` (the symbol the Fortran USERMAT hook declares via
ISO_C_BINDING) together with a tiny C driver, starts material_server.py, and
checks the values that come back through the wire equal the in-process Python
evaluation. This proves the whole bridge — the piece an FE solver will actually
call — works, not just the Python half.

Requires a C compiler; skipped automatically if absent.
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
_COUP = _ROOT / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))

import material_server as ms                      # noqa: E402

_CC = shutil.which("cc") or shutil.which("gcc")
pytestmark = pytest.mark.skipif(
    _CC is None or not (_COUP / "biofilm_py_eval.c").exists(),
    reason="C compiler or shim source unavailable")

F = [1.15, 0.03, 0.0, -0.02, 0.97, 0.01, 0.0, 0.0, 1.02]
FV = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
PARAMS = (0.2, 2.0e-4, 5.0e-5, 5000.0, 8e-3, 1.0, 5.0)   # alpha,C10,C01,D1,eta,mtype,dt


@pytest.fixture(scope="module")
def shim_exe():
    tmp = tempfile.mkdtemp()
    exe = Path(tmp) / "test_shim"
    subprocess.run([_CC, str(_COUP / "test_shim_main.c"),
                    str(_COUP / "biofilm_py_eval.c"), "-o", str(exe)], check=True)
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
    stdin = " ".join(f"{v:.17g}" for v in (*F, *FV, *PARAMS)) + "\n"
    env = {"BIOFILM_PY_HOST": host, "BIOFILM_PY_PORT": str(port), "PATH": "/usr/bin:/bin"}
    return subprocess.run([str(exe)], input=stdin, capture_output=True,
                          text=True, env=env, timeout=timeout)


def test_shim_roundtrip_matches_python(shim_exe, server):
    host, port = server
    r = _run_shim(shim_exe, host, port)
    assert r.returncode == 0, f"shim failed rc={r.returncode}: {r.stderr}"

    vals = [float(x) for x in r.stdout.split()]
    assert len(vals) == 6 + 9 + 36
    stress, fv_new, dsde = vals[:6], vals[6:15], vals[15:]

    sv, fvn, _ = ms.stress_core(np.array(F).reshape(3, 3),
                                np.array(FV).reshape(3, 3), *PARAMS)
    D = ms.dsde_perturbation(np.array(F).reshape(3, 3),
                             np.array(FV).reshape(3, 3), PARAMS)

    # %.17g round-trips a double exactly; the wire adds no error of its own.
    assert np.allclose(stress, sv, rtol=1e-13, atol=1e-18)
    assert np.allclose(fv_new, fvn.reshape(9), rtol=1e-13, atol=1e-18)
    assert np.allclose(dsde, D.reshape(36), rtol=1e-9, atol=1e-14)


def test_shim_reports_failure_when_server_absent(shim_exe):
    """With no server listening the shim must fail cleanly (nonzero rc), so the
    Fortran side can fall back to the verified inline core instead of hanging."""
    with socket.socket() as s:                    # grab a port, then free it
        s.bind(("127.0.0.1", 0))
        dead_port = s.getsockname()[1]
    r = _run_shim(shim_exe, "127.0.0.1", dead_port, timeout=20)
    assert r.returncode != 0
