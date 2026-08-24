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
import socketserver
import numpy as np

from protocol import decode_request, encode_response, encode_error

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


def evaluate(req: dict) -> bytes:
    F = np.asarray(req["F"], float).reshape(3, 3)
    Fv = np.asarray(req["Fv"], float).reshape(3, 3)
    params = (req["alpha"], req["C10"], req["C01"], req["D1"],
              req["eta"], req["mtype"], req["dt"])
    sv, Fv_new, detFe = stress_core(F, Fv, *params)
    D = dsde_perturbation(F, Fv, params)
    return encode_response(sv, Fv_new.reshape(9), detFe, D.reshape(36))


class _Handler(socketserver.StreamRequestHandler):
    def handle(self):
        for line in self.rfile:
            line = line.strip()
            if not line:
                continue
            try:
                self.wfile.write(evaluate(decode_request(line)))
            except Exception as exc:                    # keep the server alive
                self.wfile.write(encode_error(exc))


def serve(host="127.0.0.1", port=8765):
    with socketserver.TCPServer((host, port), _Handler) as srv:
        print(f"material_server listening on {host}:{port}")
        srv.serve_forever()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args()
    serve(a.host, a.port)
