#!/usr/bin/env python3
"""material_server.py — the Python side of the Gauss-point material bridge.

This is the SKELETON for the thesis' core deliverable: letting a Fortran
UMAT/USERMAT call a Python material model at each integration point. It is
deliberately runnable and testable end-to-end *without* Abaqus/ANSYS — the
Fortran side is replaced by a Python client in the tests.

`stress_core()` is a NumPy reference implementation that mirrors the verified
Fortran core `BIOFILM_STRESS_CORE` (Neo-Hookean/Mooney-Rivlin deviator + D1
pressure, backward-Euler viscous update). `dsde_perturbation()` builds the 6x6
material Jacobian by the same F-perturbation the UMAT uses. In production this
module would instead call the calibrated JAX material model (material_models.py
/ JAXFEM/); the interface is identical, so that swap is local.

  python material_server.py [--host 127.0.0.1 --port 8765]

Voigt order: Abaqus 11,22,33,12,13,23.
"""
from __future__ import annotations

import argparse
import json
import socketserver
import numpy as np

from protocol import (decode_request, encode_response, encode_error,
                      decode_ecology_request, encode_ecology_response)

I3 = np.eye(3)
VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]   # Abaqus order


def _sigma_and_fv(F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt):
    """Cauchy stress (3x3), updated Fv (3x3), detFe — mirrors the Fortran core."""
    Fg_inv = I3 / max(1.0 + alpha, 1e-15)
    Ftrial = F @ Fg_inv

    def elastic(Fv):
        Fe = Ftrial @ np.linalg.inv(Fv)
        Be = Fe @ Fe.T
        detFe = max(np.linalg.det(Fe), 1e-15)
        tmp1 = detFe ** (-2.0 / 3.0)
        i1b = tmp1 * np.trace(Be)
        return Be, detFe, tmp1, i1b

    def mr_extra(Be, tmp1, i1b):                     # C01 deviatoric contribution
        T3 = i1b * tmp1 * Be - tmp1 ** 2 * (Be @ Be)
        return T3 - (np.trace(T3) / 3.0) * I3

    # trial state → viscous flow driver (deviatoric Kirchhoff)
    Be, detFe, tmp1, i1b = elastic(Fv_old)
    tau = 2.0 * C10 * tmp1 * (Be - (i1b / 3.0) * I3)
    if mtype > 0.5:
        tau = tau + 2.0 * C01 * mr_extra(Be, tmp1, i1b)

    if eta > 1e-20:
        Fv_new = (I3 + dt / (2.0 * eta * detFe) * tau) @ Fv_old
    else:
        Fv_new = Fv_old.copy()

    # recompute with updated Fv → Cauchy stress
    Be, detFe, tmp1, i1b = elastic(Fv_new)
    press = (2.0 / D1) * (detFe - 1.0) * detFe
    sig = (2.0 * C10 * tmp1 * (Be - (i1b / 3.0) * I3) + press * I3) / detFe
    if mtype > 0.5:
        sig = sig + 2.0 * C01 * mr_extra(Be, tmp1, i1b) / detFe
    return sig, Fv_new, detFe


def stress_core(F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt):
    """Return (stress[6] Abaqus Voigt, Fv_new[3,3], detFe)."""
    sig, Fv_new, detFe = _sigma_and_fv(F, Fv_old, alpha, C10, C01, D1, eta, mtype, dt)
    sv = np.array([sig[i, j] for i, j in VOIGT])
    return sv, Fv_new, detFe


def dsde_perturbation(F, Fv_old, params, h=1.0e-7):
    """6x6 material Jacobian dσ/dε via symmetric spatial perturbations of F
    (the same F-perturbation scheme the verified UMAT uses, same step size
    PERT=1.0d-7 as usermat_biofilm.f, and deliberately NOT symmetrised --
    the inline Fortran core's unsymmetrised numerical tangent is the one
    that was confirmed to converge under SOLID185/NLGEOM,ON, so this path
    reproduces that behaviour rather than a different-but-plausible one)."""
    sv0, _, _ = stress_core(F, Fv_old, *params)
    D = np.zeros((6, 6))
    for k, (a, b) in enumerate(VOIGT):
        dE = np.zeros((3, 3))
        dE[a, b] += h / 2.0
        dE[b, a] += h / 2.0
        sv, _, _ = stress_core((I3 + dE) @ F, Fv_old, *params)
        D[:, k] = (sv - sv0) / h
    return D


# Tangent backend. "fd" is the finite difference above, matching
# usermat_biofilm.f's PERT=1e-7 exactly -- that equality is what makes
# kUsePy=1 vs kUsePy=0 a clean equivalence proof, so it stays the default.
# "jax" serves the exact AD tangent from material_jax instead (opt-in, since
# it deliberately breaks that bit-level equality at the ~3e-8 truncation
# level the FD tangent carries).
TANGENT_BACKEND = "fd"


def set_tangent_backend(name: str) -> None:
    global TANGENT_BACKEND
    if name not in ("fd", "jax"):
        raise ValueError(f"unknown tangent backend {name!r} (fd|jax)")
    if name == "jax":
        import material_jax  # noqa: F401  -- fail here, not per request
    TANGENT_BACKEND = name


def _tangent(F, Fv, params):
    if TANGENT_BACKEND == "jax":
        import material_jax
        return np.asarray(material_jax.dsde_exact(F, Fv, params), dtype=float)
    return dsde_perturbation(F, Fv, params)


def evaluate(req: dict) -> bytes:
    F = np.asarray(req["F"], float).reshape(3, 3)
    Fv = np.asarray(req["Fv"], float).reshape(3, 3)
    params = (req["alpha"], req["C10"], req["C01"], req["D1"],
              req["eta"], req["mtype"], req["dt"])
    sv, Fv_new, detFe = stress_core(F, Fv, *params)
    D = _tangent(F, Fv, params)
    return encode_response(sv, Fv_new.reshape(9), detFe, D.reshape(36))


def evaluate_ecology(req: dict) -> bytes:
    """0D Hamilton ecology ODE step -- see ecology_jax.py. Imported lazily so
    a plain material-bridge deployment (kUsePy=1, ecology unused) does not
    need jax installed, the same idiom set_tangent_backend uses for the
    optional 'jax' dsdePl backend."""
    import ecology_jax
    g_new = ecology_jax.ecology_step(req["g"], req["theta"], req["dt_h"])
    return encode_ecology_response(g_new)


class _Handler(socketserver.StreamRequestHandler):
    def handle(self):
        for line in self.rfile:
            line = line.strip()
            if not line:
                continue
            try:
                kind = json.loads(line).get("kind")
            except Exception as exc:                    # keep the server alive
                self.wfile.write(encode_error(exc))
                continue
            try:
                if kind == "ecology":
                    self.wfile.write(evaluate_ecology(decode_ecology_request(line)))
                else:
                    self.wfile.write(evaluate(decode_request(line)))
            except Exception as exc:                    # keep the server alive
                self.wfile.write(encode_error(exc))


class _Server(socketserver.ThreadingTCPServer):
    """Threaded so multiple persistent client connections can be serviced
    concurrently, not just the first one ever accepted.

    Found 2026-09-07 via a multi-element real-ANSYS run: ANSYS parallelises
    with MPI, not (only) OpenMP threads within one process -- a multi-
    element solve spawned 4 separate ANSYS.exe ranks, each opening its own
    persistent connection to this server (biofilm_py_eval.c's "one
    connection, reused for the whole run" design). Plain socketserver.
    TCPServer services one accepted connection's entire lifetime (its
    handler's `for line in self.rfile:` loop) before ever accepting the
    next -- so rank 0's connection monopolised the server for the rest of
    the run, and ranks 1-3 hung forever waiting for a response that could
    never come. This was mistaken at first for a client-side data race
    (fixed separately, and still a real, necessary fix for same-process
    OpenMP-style concurrency -- see biofilm_py_eval.c) because the
    observed symptom was so similar; it took reading the MPI "BAD
    TERMINATION ... RANK 0/1/2/3" messages from a killed hung run to find
    the real cause here.

    Python's GIL still serialises the actual evaluate()/evaluate_ecology()
    computation across threads, so this adds no real parallelism -- it
    only stops the server from permanently starving every client after
    the first."""
    daemon_threads = True
    allow_reuse_address = True


def serve(host="127.0.0.1", port=8765):
    with _Server((host, port), _Handler) as srv:
        print(f"material_server listening on {host}:{port}")
        srv.serve_forever()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--tangent", choices=("fd", "jax"), default="fd",
                    help="dsdePl backend: fd = finite difference matching the "
                         "USERMAT's PERT=1e-7 (default, keeps kUsePy=1 vs "
                         "kUsePy=0 an exact equivalence check); jax = exact "
                         "forward-mode AD (requires jax)")
    a = ap.parse_args()
    set_tangent_backend(a.tangent)
    serve(a.host, a.port)
