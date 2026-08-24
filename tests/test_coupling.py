"""End-to-end test of the Gauss-point material bridge skeleton.

Exercises the Python side exactly as the Fortran USERMAT will: marshal one
deformation state through the wire protocol to material_server and check the
response is a valid material evaluation. No Abaqus/ANSYS needed.
"""
import socket
import sys
import threading
from pathlib import Path

import numpy as np
import pytest

_COUP = Path(__file__).resolve().parents[1] / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))

import material_server as ms          # noqa: E402
import protocol                        # noqa: E402

PARAMS = dict(alpha=0.2, C10=2.0e-4, C01=5.0e-5, D1=5000.0, eta=8e-3, mtype=1.0, dt=5.0)
F = [1.15, 0.03, 0.0, -0.02, 0.97, 0.01, 0.0, 0.0, 1.02]
FV = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]


def _in_process():
    sv, fvn, det = ms.stress_core(
        np.array(F).reshape(3, 3), np.array(FV).reshape(3, 3),
        PARAMS["alpha"], PARAMS["C10"], PARAMS["C01"], PARAMS["D1"],
        PARAMS["eta"], PARAMS["mtype"], PARAMS["dt"])
    return sv, fvn, det


def test_core_is_finite_and_reasonable():
    sv, fvn, det = _in_process()
    assert np.all(np.isfinite(sv)) and np.all(np.isfinite(fvn))
    assert det > 0
    assert np.linalg.det(fvn) > 0            # viscous update stays invertible


def test_tangent_shape_and_approximate_symmetry():
    params = (PARAMS["alpha"], PARAMS["C10"], PARAMS["C01"], PARAMS["D1"],
              PARAMS["eta"], PARAMS["mtype"], PARAMS["dt"])
    D = ms.dsde_perturbation(np.array(F).reshape(3, 3), np.array(FV).reshape(3, 3), params)
    assert D.shape == (6, 6)
    assert np.all(np.isfinite(D))
    # dsde_perturbation deliberately does NOT symmetrise (see its docstring):
    # it mirrors usermat_biofilm.f's own unsymmetrised forward-difference
    # tangent, which was the one confirmed to converge under SOLID185/
    # NLGEOM,ON. So only approximate (truncation-level) symmetry holds here;
    # exact agreement with the Fortran core -- including this same asymmetry
    # -- is what tests/test_usermat_kusepy_e2e.py checks end to end.
    assert np.max(np.abs(D - D.T)) < 0.2 * max(np.max(np.abs(D)), 1e-12)


def test_socket_roundtrip():
    """Start the server, send one request as the Fortran client would, and
    check the wire response matches the in-process evaluation."""
    try:
        srv = ms.socketserver.TCPServer(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    host, port = srv.server_address
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        with socket.create_connection((host, port), timeout=5) as c:
            c.sendall(protocol.encode_request(F, FV, PARAMS["alpha"], PARAMS["C10"],
                      PARAMS["C01"], PARAMS["D1"], PARAMS["eta"], PARAMS["mtype"], PARAMS["dt"]))
            resp = protocol.decode_response(c.makefile("rb").readline())
    finally:
        srv.shutdown()
        srv.server_close()
    assert "error" not in resp, resp
    sv, fvn, det = _in_process()
    assert np.allclose(resp["stress"], sv, rtol=1e-12, atol=1e-14)
    assert np.allclose(resp["Fv_new"], fvn.reshape(9), rtol=1e-12, atol=1e-14)
    assert len(resp["dsdePl"]) == 36
