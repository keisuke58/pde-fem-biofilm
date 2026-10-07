#!/usr/bin/env python3
"""Stress into growth in 3D: the 2D prototype (stress_to_growth_2d.py) on one
octant of the 2 mm cube, with symmetry planes x = y = z = 0, for a spherical
seed and a cubic seed (with corners).

Same model as in 2D (Klempt 2024 Eq. 30 and 36 divided by eta, no nutrient,
no front term):

    phi_dot = beta lap(phi) + k_alpha alpha - mu_star phi m,  alpha_dot = k_alpha phi,
    m = I : C_e_iso - 3 ~= 2 |dev eps_e|^2,

semi-implicit in the new term, mechanics from the step before (staggered).
Hex8, small strain, selective reduced integration, E (phi^2 + f), E = 10 Pa,
nu = 0.49. The question: in 2D plane strain the seed has deviatoric strain
inside (c = m/(alpha-1)^2 = 3.4 from the z constraint); a sphere in 3D should
have almost none inside, so the term would act only at the surface.

    python keio_wp2/stress_to_growth_3d.py [--n 16] -> keio_wp2/results_3d.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
L, E0, NU, FLOOR = 1.0, 10.0, 0.49, 1e-3      # octant edge [mm]
BETA, K_ALPHA = 0.02, 1e-3
SEEDS = {"sphere": 0.3, "cube": 0.25}          # radius / half edge [mm]
MU_TABLE2 = E0 / (2 * (1 + NU)) / 1e-10


class Mesh:
    def __init__(self, n, seed):
        self.n, self.h = n, L / n
        xc = (np.arange(n) + 0.5) * self.h
        X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
        if seed == "sphere":
            ins = X ** 2 + Y ** 2 + Z ** 2 <= SEEDS["sphere"] ** 2
        else:
            a = SEEDS["cube"]
            ins = (X <= a) & (Y <= a) & (Z <= a)
        self.inside = ins.ravel()
        # interior of the seed: elements whose 26 neighbours are all inside
        p = np.pad(ins, 1, mode="edge")
        core = np.ones_like(ins)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    core &= p[1 + dx:n + 1 + dx, 1 + dy:n + 1 + dy, 1 + dz:n + 1 + dz]
        self.core = core.ravel()
        self.rim = self.inside & ~self.core
        nn = n + 1
        idx = np.arange(nn ** 3).reshape(nn, nn, nn)
        c = [idx[:-1, :-1, :-1], idx[1:, :-1, :-1], idx[1:, 1:, :-1], idx[:-1, 1:, :-1],
             idx[:-1, :-1, 1:], idx[1:, :-1, 1:], idx[1:, 1:, 1:], idx[:-1, 1:, 1:]]
        self.conn = np.stack(c, -1).reshape(-1, 8)
        self.ndof = 3 * nn ** 3
        h = self.h
        corners = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)

        def bmat(xi):
            dN = np.array([[cx * (1 + xi[1] * cy) * (1 + xi[2] * cz),
                            cy * (1 + xi[0] * cx) * (1 + xi[2] * cz),
                            cz * (1 + xi[0] * cx) * (1 + xi[1] * cy)] for cx, cy, cz in corners]) / 8 * 2 / h
            b = np.zeros((6, 24))
            for a in range(8):
                dx, dy, dz = dN[a]; i = 3 * a
                b[0, i] = dx; b[1, i + 1] = dy; b[2, i + 2] = dz
                b[3, i] = dy; b[3, i + 1] = dx
                b[4, i + 1] = dz; b[4, i + 2] = dy
                b[5, i] = dz; b[5, i + 2] = dx
            return b
        g = 1 / np.sqrt(3)
        Bg = [bmat(np.array([a, b, cc])) for a in (-g, g) for b in (-g, g) for cc in (-g, g)]
        self.B0 = bmat(np.zeros(3))
        Ddev = np.zeros((6, 6))
        Ddev[:3, :3] = 2 * (np.eye(3) - 1 / 3)
        Ddev[3:, 3:] = np.eye(3)
        Dvol = np.zeros((6, 6)); Dvol[:3, :3] = 1.0
        w = (h / 2) ** 3
        self.Kdev = sum(b.T @ Ddev @ b * w for b in Bg)
        self.Kvol = self.B0.T @ Dvol @ self.B0 * (8 * w)
        self.fvol = self.B0.T @ np.array([1, 1, 1, 0, 0, 0.0]) * (8 * w)
        self.dofs = np.stack([3 * self.conn + k for k in range(3)], -1).reshape(-1, 24)
        self.rows = np.repeat(self.dofs, 24, axis=1).ravel()
        self.cols = np.tile(self.dofs, (1, 24)).ravel()
        I, J, K = np.meshgrid(np.arange(nn), np.arange(nn), np.arange(nn), indexing="ij")
        nid = idx.ravel()
        fixed = np.concatenate([3 * nid[I.ravel() == 0], 3 * nid[J.ravel() == 0] + 1, 3 * nid[K.ravel() == 0] + 2])
        self.free = np.setdiff1d(np.arange(self.ndof), fixed)

    def lap(self, f):
        n = self.n
        p = np.pad(f.reshape(n, n, n), 1, mode="edge")
        c = p[1:-1, 1:-1, 1:-1]
        lap = (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1]
               + p[1:-1, 1:-1, 2:] + p[1:-1, 1:-1, :-2] - 6 * c) / self.h ** 2
        return lap.ravel()


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
    exx, eyy, ezz = e[:, 0] - delta, e[:, 1] - delta, e[:, 2] - delta
    tr = exx + eyy + ezz
    dev2 = ((exx - tr / 3) ** 2 + (eyy - tr / 3) ** 2 + (ezz - tr / 3) ** 2
            + 2 * ((e[:, 3] / 2) ** 2 + (e[:, 4] / 2) ** 2 + (e[:, 5] / 2) ** 2))
    return 2 * dev2


def run(seed, n, mu_star, dt=0.01, t_end=1.0):
    M = Mesh(n, seed)
    phi = M.inside.astype(float)
    alpha = np.ones_like(phi)
    hist = []
    for k in range(int(round(t_end / dt))):
        m = mechanics(M, phi, alpha - 1)
        src = BETA * M.lap(phi) + K_ALPHA * alpha
        phi = np.clip((phi + dt * src) / (1 + dt * mu_star * m), 0.0, 1.0)
        alpha = alpha + dt * K_ALPHA * phi
        d = alpha - 1
        hist.append(((k + 1) * dt, float(d[M.inside].mean()), float(phi[M.inside].mean()),
                     float(phi[M.core].mean()) if M.core.any() else float("nan"),
                     float(phi[M.rim].mean()),
                     float(m[M.core].mean() / d[M.core].mean() ** 2) if M.core.any() else float("nan"),
                     float(m[M.rim].mean() / d[M.rim].mean() ** 2)))
    return dict(seed=seed, n=n, mu_star=mu_star, dt=dt, n_core=int(M.core.sum()), n_rim=int(M.rim.sum()),
                hist=hist)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    a = ap.parse_args(argv)
    t0 = time.time()
    out = {"mu_table2": MU_TABLE2, "n": a.n, "runs": []}
    for seed in ("sphere", "cube"):
        for ms in [0.0, 1e5, 1e6, 1e7, MU_TABLE2]:
            r = run(seed, a.n, ms, dt=0.02)
            out["runs"].append(r)
            t, d, p, pc, pr, cc, cr = r["hist"][-1]
            print(f"{seed:6s} mu*={ms:9.3g}: alpha-1 {d:.3e}, phi seed {p:.3f} core {pc:.3f} rim {pr:.3f}, "
                  f"c core {cc:.3f} rim {cr:.3f} ({r['n_core']}/{r['n_rim']} el.)  [{time.time() - t0:.0f} s]",
                  flush=True)
    (HERE / f"results_3d_n{a.n}.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
