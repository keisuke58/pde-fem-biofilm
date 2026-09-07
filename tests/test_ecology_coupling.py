"""Python-side round-trip test for the 0D Hamilton ecology Gauss-point bridge
(coupling/ecology_jax.py, coupling/README.md next-steps #4).

Mirrors test_coupling.py's structure for the material bridge: exercises the
wire protocol and the in-process socket server exactly as the Fortran
USERMAT's ecology hook will, with no ANSYS needed. ecology_jax.ecology_step
is a thin wrapper around jax_hamilton_0d_5species_demo.newton_step_jit, so
these tests also guard against that wrapper silently drifting from the
verified 0D reference.
"""
import socket
import sys
import threading
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
_COUP = _ROOT / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))
sys.path.insert(0, str(_ROOT))

import ecology_jax                                     # noqa: E402
import material_server as ms                           # noqa: E402
import protocol                                         # noqa: E402
from jax_hamilton_0d_5species_demo import (             # noqa: E402
    THETA_DEMO, newton_step_jit, theta_to_matrices,
)

G0 = [0.12, 0.10, 0.05, 0.02, 0.0, 0.71, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0]
DT_H = 1.0e-5


def _reference_step(g, theta, dt_h):
    A, b_diag = theta_to_matrices(np.asarray(theta, dtype=np.float64))
    params = ecology_jax.default_hparams(dt_h)
    params["A"] = A
    params["b_diag"] = b_diag
    return np.asarray(newton_step_jit(np.asarray(g, dtype=np.float64), params))


def test_ecology_step_matches_demo_newton_step():
    """ecology_jax.ecology_step must not drift from the verified 0D reference
    it wraps -- rebuilds the same call independently rather than reusing the
    wrapper's own theta_to_matrices/default_hparams call sequence."""
    g_new = ecology_jax.ecology_step(G0, THETA_DEMO, DT_H)
    ref = _reference_step(G0, THETA_DEMO, DT_H)
    np.testing.assert_allclose(np.asarray(g_new), ref, rtol=1e-13, atol=1e-15)


def test_step_conserves_phi_partition():
    """g[5] (phi0) plus sum(phi) must stay at 1 -- the residual's own closure
    constraint (jax_hamilton_0d_5species_demo.residual, Q[11])."""
    g_new = np.asarray(ecology_jax.ecology_step(G0, THETA_DEMO, DT_H))
    assert np.isclose(np.sum(g_new[0:5]) + g_new[5], 1.0, atol=1e-8)


def test_living_fraction_total_in_range():
    g_new = ecology_jax.ecology_step(G0, THETA_DEMO, DT_H)
    phi_tot = ecology_jax.living_fraction_total(g_new)
    assert 0.0 <= phi_tot <= 1.0 + 1e-8


def test_default_initial_state_is_valid_composition():
    g0 = np.asarray(ecology_jax.default_initial_state())
    assert np.isclose(np.sum(g0[0:5]) + g0[5], 1.0, atol=1e-6)
    assert np.all(g0[0:5] >= 0.0) and np.all(g0[6:11] > 0.0)


def test_protocol_roundtrip():
    req = protocol.decode_ecology_request(
        protocol.encode_ecology_request(G0, THETA_DEMO, DT_H))
    assert req["kind"] == "ecology"
    resp = protocol.decode_response(ms.evaluate_ecology(req))
    assert "error" not in resp, resp
    np.testing.assert_allclose(
        resp["g_new"], np.asarray(ecology_jax.ecology_step(G0, THETA_DEMO, DT_H)),
        rtol=1e-12, atol=1e-14)


def test_server_dispatches_ecology_by_kind():
    """A request carrying kind="ecology" must be routed to evaluate_ecology,
    not the default material path, and a plain material request (no "kind")
    must still work unchanged -- the backward-compatibility guarantee
    protocol.py's docstring makes."""
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    host, port = srv.server_address
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        with socket.create_connection((host, port), timeout=5) as c, \
             c.makefile("rb") as rfile:
            # NB: c.makefile() holds its own reference to the socket fd, so
            # closing only `c` does not send the TCP FIN and the server's
            # `for line in self.rfile:` loop never sees EOF -- both must be
            # closed (the nested `with` above) or shutdown() below deadlocks
            # waiting for the still-running handler thread to return.
            c.sendall(protocol.encode_ecology_request(G0, THETA_DEMO, DT_H))
            eco_resp = protocol.decode_response(rfile.readline())

            c.sendall(protocol.encode_request(
                [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0], [1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0],
                0.1, 1.0, 0.0, 0.01, 0.0, 0.0, 1.0))
            mat_resp = protocol.decode_response(rfile.readline())
    finally:
        srv.shutdown()
        srv.server_close()

    assert "error" not in eco_resp and "g_new" in eco_resp
    assert "error" not in mat_resp and "stress" in mat_resp
    np.testing.assert_allclose(
        eco_resp["g_new"], np.asarray(ecology_jax.ecology_step(G0, THETA_DEMO, DT_H)),
        rtol=1e-12, atol=1e-14)
