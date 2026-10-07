#!/usr/bin/env python3
"""Keio WP2, step 4: the homeostatic-pressure growth law for a biofilm on a
rigid substrate (tooth or implant surface), 2D plane-strain cross-section.

Domain 2 mm wide (x) and 1 mm high (y), substrate at y = 0 (u = 0), symmetry
on the sides (u_x = 0), free top. Biofilm at the start (phi = 1):
  film   y < 0.2 mm over the whole width (laterally confined by the substrate)
  patch  y < 0.2 mm and |x| < 0.5 mm (a colony that can also spread sideways)
Elsewhere phi = 0 with stiffness floor f = 1e-3. Same equations as
homeostatic.py:

    phi_dot   = beta lap(phi) + k_alpha alpha               (Klempt 2024 Eq. 34, no front/nutrient)
    alpha_dot = k_alpha phi max(0, 1 - p / p_h)              (Eq. 36 with a homeostatic pressure)
    p = -K tr(eps_e)  (compression > 0), plane strain, E (phi^2 + f), E = 10 Pa, nu = 0.49

p_h = P_h E k_alpha T*, P_h in {none, 1, 0.5, 0.25, 0.1}. The scale is E, not the bulk
modulus K: with nu = 0.49, K = 17 E, and a layer with a free top relieves the
hydrostatic part by thickening, so its pressure stays of order E (alpha - 1). For a thin film on a
rigid substrate the in-plane growth is blocked, so the stress is set by the
substrate and not by the surroundings, unlike the seed in a soft void
(homeostatic_confinement.py).

    python keio_wp2/homeostatic_layer.py -> keio_wp2/results_layer.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
LX, LY, H0, W = 2.0, 1.0, 0.2, 0.5
E0, NU, FLOOR, BETA, K_ALPHA = 10.0, 0.49, 1e-3, 0.02, 1e-3
K_BIO = E0 / (3 * (1 - 2 * NU))


class Mesh:
    def __init__(self, nx, ny, geom):
        self.nx, self.ny = nx, ny
        self.h = LX / nx
        assert abs(LY / ny - self.h) < 1e-12
        xc = (np.arange(nx) + 0.5) * self.h - LX / 2
        yc = (np.arange(ny) + 0.5) * self.h
        X, Y = np.meshgrid(xc, yc, indexing="ij")
        ins = Y < H0 if geom == "film" else (Y < H0) & (np.abs(X) < W)
        self.inside = ins.ravel()
        self.bottom = (Y < self.h).ravel()                # first element row on the substrate
        nnx, nny = nx + 1, ny + 1
        idx = np.arange(nnx * nny).reshape(nnx, nny)
        self.conn = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
        self.ndof = 2 * nnx * nny
        h = self.h
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)

        def bmat(xi, eta):
            dN = np.array([[cx * (1 + eta * cy), cy * (1 + xi * cx)] for cx, cy in corners]) / 4 * 2 / h
            b = np.zeros((3, 8))
            b[0, 0::2] = dN[:, 0]; b[1, 1::2] = dN[:, 1]
            b[2, 0::2] = dN[:, 1]; b[2, 1::2] = dN[:, 0]
            return b
        g = 1 / np.sqrt(3)
        Bg = [bmat(a, b) for a in (-g, g) for b in (-g, g)]
        self.B0 = bmat(0.0, 0.0)
        Ddev = np.array([[4 / 3, -2 / 3, 0], [-2 / 3, 4 / 3, 0], [0, 0, 1]])
        Dvol = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]], float)
        w = (h / 2) ** 2
        self.Kdev = sum(b.T @ Ddev @ b * w for b in Bg)
        self.Kvol = self.B0.T @ Dvol @ self.B0 * (4 * w)
        self.fvol = self.B0.T @ np.array([1.0, 1.0, 0.0]) * (4 * w)
        self.dofs = np.stack([2 * self.conn, 2 * self.conn + 1], -1).reshape(-1, 8)
        self.rows = np.repeat(self.dofs, 8, axis=1).ravel()
        self.cols = np.tile(self.dofs, (1, 8)).ravel()
        I, J = np.meshgrid(np.arange(nnx), np.arange(nny), indexing="ij")
        nid = idx.ravel(); I = I.ravel(); J = J.ravel()
        fixed = np.concatenate([2 * nid[J == 0], 2 * nid[J == 0] + 1,             # substrate
                                2 * nid[(I == 0) | (I == nx)]])                     # symmetry sides, u_x = 0
        self.free = np.setdiff1d(np.arange(self.ndof), np.unique(fixed))

    def lap(self, f):
        p = np.pad(f.reshape(self.nx, self.ny), 1, mode="edge")
        return ((p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4 * p[1:-1, 1:-1]) / self.h ** 2).ravel()


def mechanics(M, phi, delta):
    Ee = E0 * (phi ** 2 + FLOOR)
    mu = Ee / (2 * (1 + NU)); K = Ee / (3 * (1 - 2 * NU))
    vals = (mu[:, None, None] * M.Kdev + K[:, None, None] * M.Kvol).ravel()
    Kg = sp.csr_matrix((vals, (M.rows, M.cols)), shape=(M.ndof, M.ndof))
    f = np.zeros(M.ndof)
    np.add.at(f, M.dofs.ravel(), ((3 * K * delta)[:, None] * M.fvol).ravel())
    u = np.zeros(M.ndof)
    u[M.free] = spla.spsolve(Kg[M.free][:, M.free].tocsc(), f[M.free])
    e = (M.B0 @ u[M.dofs].T).T
    tr_e = e[:, 0] + e[:, 1] - 3 * delta                 # eps_zz = 0 (plane strain)
    exx, eyy, ezz = e[:, 0] - delta, e[:, 1] - delta, -delta
    d = np.stack([exx - tr_e / 3, eyy - tr_e / 3, ezz - tr_e / 3], 1)
    vm = np.sqrt(1.5) * 2 * mu * np.sqrt((d ** 2).sum(1) + 2 * (e[:, 2] / 2) ** 2)
    sxx = 2 * mu * d[:, 0] + K * tr_e
    return -K * tr_e, vm, sxx, u


def run(geom, P_h, nx=64, ny=32, dt=0.02, t_end=1.0):
    M = Mesh(nx, ny, geom)
    p_h = None if P_h is None else P_h * E0 * K_ALPHA
    phi = M.inside.astype(float)
    alpha = np.ones_like(phi)
    hist = []
    for k in range(int(round(t_end / dt))):
        p, vm, sxx, u = mechanics(M, phi, alpha - 1)
        g = np.ones_like(p) if p_h is None else np.clip(1 - p / p_h, 0.0, 1.0)
        phi = np.clip(phi + dt * (BETA * M.lap(phi) + K_ALPHA * alpha), 0, 1)
        alpha = alpha + dt * K_ALPHA * phi * g
        ins, bot = M.inside, M.inside & M.bottom
        uy = u[1::2].reshape(nx + 1, ny + 1)
        hist.append(((k + 1) * dt, float((alpha - 1)[ins].mean()), float(phi[ins].mean()),
                     float(p[ins].mean()), float(p[ins].max()), float(vm[ins].mean()),
                     float(sxx[bot].mean()), float(g[ins].mean()), float(uy.max())))
    return dict(geom=geom, P_h=P_h, hist=hist)


def main():
    t0 = time.time()
    out = {"K_bio": K_BIO, "scale_E_k_alpha": E0 * K_ALPHA, "runs": []}
    for geom in ("film", "patch"):
        for P_h in (None, 1.0, 0.5, 0.25, 0.1):
            r = run(geom, P_h)
            out["runs"].append(r)
            t, d, ph, pm, px, vm, sxx, g, uy = r["hist"][-1]
            print(f"{geom:5s} P_h={P_h}: alpha-1 {d:.3e}, phi {ph:.3f}, p mean {pm:.3e} max {px:.3e} Pa, "
                  f"vM {vm:.3e}, sxx at substrate {sxx:.3e}, g {g:.2f}, top lift {uy:.3e} mm "
                  f"[{time.time() - t0:.0f} s]", flush=True)
    (HERE / "results_layer.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
