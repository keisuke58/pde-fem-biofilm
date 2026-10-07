#!/usr/bin/env python3
"""Keio WP2, step 5 (towards real shapes, P3): the homeostatic-pressure
growth law on the biofilm layer around an implant collar, axisymmetric.

Geometry as abaqus_composition/make_implant_inp.py (values chosen for this
work, not from a paper): collar radius r_i = 2.05 mm (a 4.1 mm implant),
biofilm layer 0.25 mm thick (r_i ... r_i + 0.25) and 2.0 mm high (a sulcus).
The inner face is bonded to the titanium, taken as rigid (u = 0); the bottom
(z = 0, bone level) has u_z = 0; the outer face and the top are free.
phi = 1 in the inner quarter of the layer at the start (the first colonisers
on the titanium), 0 elsewhere (stiffness floor f = 1e-3); phi has no flux
through any face.

    phi_dot   = beta lap(phi) + k_alpha alpha                 (Klempt 2024 Eq. 34, no front term)
    alpha_dot = k_alpha phi max(0, 1 - p / p_h)                (Eq. 36 with a homeostatic pressure)
    p = -K tr(eps_e), axisymmetric small strain, E (phi^2 + f), E = 10 Pa, nu = 0.49

p_h = P_h E k_alpha T*, P_h in {none, 1, 0.5, 0.25, 0.1}; second set with p_h
proportional to the local stiffness, p_h = P_h E (phi^2 + f) k_alpha T*. The nutrient is left
out: without the front term it does not enter the growth. Reported: mean
alpha - 1 and phi in the layer, the largest pressure and von Mises stress,
and the hoop and shear stress on the titanium.

    python keio_wp2/homeostatic_implant.py -> keio_wp2/results_implant.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
RI, T_LAYER, HZ = 2.05, 0.25, 2.0
E0, NU, FLOOR, BETA, K_ALPHA = 10.0, 0.49, 1e-3, 0.02, 1e-3


class Mesh:
    def __init__(self, nr=8, nz=64):
        self.nr, self.nz = nr, nz
        self.hr, self.hz = T_LAYER / nr, HZ / nz
        rc = RI + (np.arange(nr) + 0.5) * self.hr
        zc = (np.arange(nz) + 0.5) * self.hz
        R, Z = np.meshgrid(rc, zc, indexing="ij")
        self.rc = R.ravel()
        self.seed = (R < RI + T_LAYER / 4).ravel()
        self.wall = (R < RI + self.hr).ravel()             # element ring on the titanium
        nnr, nnz = nr + 1, nz + 1
        idx = np.arange(nnr * nnz).reshape(nnr, nnz)
        self.conn = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
        rn = RI + np.arange(nnr) * self.hr
        zn = np.arange(nnz) * self.hz
        RN, ZN = np.meshgrid(rn, zn, indexing="ij")
        self.xn = np.c_[RN.ravel(), ZN.ravel()]
        self.ndof = 2 * nnr * nnz
        self.dofs = np.stack([2 * self.conn, 2 * self.conn + 1], -1).reshape(-1, 8)
        self.rows = np.repeat(self.dofs, 8, axis=1).ravel()
        self.cols = np.tile(self.dofs, (1, 8)).ravel()
        I, J = np.meshgrid(np.arange(nnr), np.arange(nnz), indexing="ij")
        nid = idx.ravel(); I = I.ravel(); J = J.ravel()
        fixed = np.concatenate([2 * nid[I == 0], 2 * nid[I == 0] + 1, 2 * nid[J == 0] + 1])
        self.free = np.setdiff1d(np.arange(self.ndof), np.unique(fixed))
        # element matrices per element (depend on r): unit mu and unit K parts, eigenstrain load
        g = 1 / np.sqrt(3)
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
        Ddev = np.zeros((4, 4)); Ddev[:3, :3] = 2 * (np.eye(3) - 1 / 3); Ddev[3, 3] = 1.0
        Dvol = np.zeros((4, 4)); Dvol[:3, :3] = 1.0

        def bmat(e, xi, eta):
            N = (1 + xi * corners[:, 0]) * (1 + eta * corners[:, 1]) / 4
            dN = np.c_[corners[:, 0] * (1 + eta * corners[:, 1]) / 4 * 2 / self.hr,
                       corners[:, 1] * (1 + xi * corners[:, 0]) / 4 * 2 / self.hz]
            r = N @ self.xn[self.conn[e], 0]
            b = np.zeros((4, 8))
            b[0, 0::2] = dN[:, 0]; b[1, 1::2] = dN[:, 1]; b[2, 0::2] = N / r
            b[3, 0::2] = dN[:, 1]; b[3, 1::2] = dN[:, 0]
            return b, r
        ne = len(self.conn)
        w = self.hr * self.hz / 4
        self.Kdev = np.zeros((ne, 8, 8)); self.Kvol = np.zeros((ne, 8, 8)); self.fv = np.zeros((ne, 8))
        self.B0 = np.zeros((ne, 4, 8))
        for e in range(ne):
            for a in (-g, g):
                for b_ in (-g, g):
                    B, r = bmat(e, a, b_)
                    self.Kdev[e] += B.T @ Ddev @ B * w * r
            B0, r0 = bmat(e, 0.0, 0.0)
            self.B0[e] = B0
            self.Kvol[e] = B0.T @ Dvol @ B0 * 4 * w * r0
            self.fv[e] = B0.T @ np.array([1, 1, 1, 0.0]) * 4 * w * r0

    def lap(self, f):
        """axisymmetric Laplacian on the element centres, zero flux on every face."""
        F = f.reshape(self.nr, self.nz)
        P = np.pad(F, 1, mode="edge")
        r = self.rc.reshape(self.nr, self.nz)
        rp, rm = r + self.hr / 2, r - self.hr / 2
        lr = (rp * (P[2:, 1:-1] - F) - rm * (F - P[:-2, 1:-1])) / (r * self.hr ** 2)
        lz = (P[1:-1, 2:] - 2 * F + P[1:-1, :-2]) / self.hz ** 2
        return (lr + lz).ravel()


def mechanics(M, phi, delta):
    Ee = E0 * (phi ** 2 + FLOOR)
    mu = Ee / (2 * (1 + NU)); K = Ee / (3 * (1 - 2 * NU))
    vals = (mu[:, None, None] * M.Kdev + K[:, None, None] * M.Kvol).ravel()
    Kg = sp.csr_matrix((vals, (M.rows, M.cols)), shape=(M.ndof, M.ndof))
    f = np.zeros(M.ndof)
    np.add.at(f, M.dofs.ravel(), ((3 * K * delta)[:, None] * M.fv).ravel())
    u = np.zeros(M.ndof)
    u[M.free] = spla.spsolve(Kg[M.free][:, M.free].tocsc(), f[M.free])
    e = np.einsum("eij,ej->ei", M.B0, u[M.dofs])          # rr, zz, tt, rz
    ee = e[:, :3] - delta[:, None]
    tr = ee.sum(1)
    d = ee - tr[:, None] / 3
    s = 2 * mu[:, None] * d + K[:, None] * tr[:, None]
    srz = mu * e[:, 3]
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2) + 3 * srz ** 2)
    return -K * tr, vm, s[:, 2], srz


def run(P_h, local=False, dt=0.02, t_end=1.0):
    M = Mesh()
    phi = M.seed.astype(float)
    alpha = np.ones_like(phi)
    # explicit diffusion: sub-steps below h^2 / (4 beta) (dt = 0.02 alone is unstable on this mesh)
    nsub = int(np.ceil(dt / (0.8 * min(M.hr, M.hz) ** 2 / (4 * BETA))))
    hist = []
    for k in range(int(round(t_end / dt))):
        p, vm, stt, srz = mechanics(M, phi, alpha - 1)
        if P_h is None:
            g = np.ones_like(p)
        else:
            p_h = P_h * E0 * K_ALPHA * ((phi ** 2 + FLOOR) if local else 1.0)
            g = np.clip(1 - p / p_h, 0.0, 1.0)
        for _ in range(nsub):
            phi = np.clip(phi + dt / nsub * (BETA * M.lap(phi) + K_ALPHA * alpha), 0, 1)
        alpha = alpha + dt * K_ALPHA * phi * g
        hist.append(((k + 1) * dt, float((alpha - 1).mean()), float(phi.mean()), float(p.max()), float(vm.max()),
                     float(stt[M.wall].min()), float(np.abs(srz[M.wall]).max()), float(g.mean())))
    return dict(P_h=P_h, local=local, hist=hist)


def main():
    t0 = time.time()
    out = {"scale_E_k_alpha": E0 * K_ALPHA, "runs": []}
    for local, P_h in [(False, None)] + [(lc, x) for lc in (False, True) for x in (1.0, 0.5, 0.25, 0.1)]:
        r = run(P_h, local)
        out["runs"].append(r)
        t, d, ph, px, vm, stt, srz, g = r["hist"][-1]
        print(f"{'local' if local else 'E    '} P_h={P_h}: alpha-1 {d:.3e}, phi {ph:.3f}, p max {px:.3e}, vM max {vm:.3e}, "
              f"hoop on Ti {stt:.3e}, shear on Ti {srz:.3e} Pa, g {g:.2f} [{time.time() - t0:.0f} s]", flush=True)
    (HERE / "results_implant.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
