#!/usr/bin/env python3
"""Keio WP2, step 6: the homeostatic-pressure growth law on a tooth crown,
axisymmetric, compared with the implant collar of step 5.

Geometry as abaqus_composition/make_implant_inp.py --bulge b (values chosen
for this work, not from a paper): the biofilm layer, 0.25 mm thick and 2 mm
high, lies on a surface of radius r_in(z) = r_i + b sin(pi z / 2h), r_i = 2.05 mm.
b = 0 is the implant collar (step 5); b = 1 mm is a crown that widens from
the cervical line to the height of contour (cf. tier2b_real/mesh_crown.py,
2.5 -> 3.8 mm over 2 mm). The inner face is bonded to the enamel or the
titanium, both taken as rigid (u = 0); the bottom has u_z = 0; the outer face
and the top are free. phi = 1 on the inner quarter of the layer at the start.
The nutrient is left out (no front term), as in step 5.

The tooth mesh is skewed, so unlike homeostatic_implant.py phi is nodal and
diffuses with linear finite elements (lumped mass, explicit sub-steps below
the Gershgorin bound); alpha lives on the elements and the element phi is the
mean of its nodes. Axisymmetric Q4 with selective reduced integration for the
mechanics. Equations as in step 5:

    phi_dot   = beta lap(phi) + k_alpha alpha,
    alpha_dot = k_alpha phi max(0, 1 - p / p_h),   p = -K tr(eps_e),
    p_h = P_h E k_alpha T*  ("E")   or   P_h E (phi^2 + f) k_alpha T*  ("local").

On the bonded face the traction is split into the normal part (tension
pulls the biofilm off the surface) and the shear along the surface.

    python keio_wp2/homeostatic_tooth.py -> keio_wp2/results_tooth.json
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
E0, NU, FLOOR, K_ALPHA = 10.0, 0.49, 1e-3, 1e-3
G = 1 / np.sqrt(3)
CORNERS = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)


def shape(xi, eta):
    N = (1 + xi * CORNERS[:, 0]) * (1 + eta * CORNERS[:, 1]) / 4
    dN = np.c_[CORNERS[:, 0] * (1 + eta * CORNERS[:, 1]), CORNERS[:, 1] * (1 + xi * CORNERS[:, 0])] / 4
    return N, dN


class Mesh:
    def __init__(self, bulge, nr=8, nz=64):
        self.nr, self.nz, self.bulge = nr, nz, bulge
        nnr, nnz = nr + 1, nz + 1
        I, K = np.meshgrid(np.arange(nnr), np.arange(nnz), indexing="ij")
        z = HZ * K / nz
        self.rin = lambda zz: RI + bulge * np.sin(0.5 * np.pi * zz / HZ)  # noqa: E731
        r = self.rin(z) + T_LAYER * I / nr
        self.xn = np.c_[r.ravel(), z.ravel()]
        idx = np.arange(nnr * nnz).reshape(nnr, nnz)
        self.conn = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
        ne, nn = len(self.conn), nnr * nnz
        self.nn, self.ndof = nn, 2 * nn
        self.dofs = np.stack([2 * self.conn, 2 * self.conn + 1], -1).reshape(-1, 8)
        self.rows = np.repeat(self.dofs, 8, axis=1).ravel()
        self.cols = np.tile(self.dofs, (1, 8)).ravel()
        Ir, Kr = I.ravel(), K.ravel()
        fixed = np.concatenate([2 * idx.ravel()[Ir == 0], 2 * idx.ravel()[Ir == 0] + 1, 2 * idx.ravel()[Kr == 0] + 1])
        self.free = np.setdiff1d(np.arange(self.ndof), np.unique(fixed))
        self.seed_n = (Ir <= 2).astype(float)          # nodes within the inner quarter (0.0625 mm)
        self.wall = np.repeat(np.arange(nr) == 0, nz)  # element ring on the bonded face
        # wall direction in (r, z) at the element centres of the first ring
        zc = (np.arange(nz) + 0.5) * HZ / nz
        slope = bulge * 0.5 * np.pi / HZ * np.cos(0.5 * np.pi * zc / HZ)
        nrm = np.sqrt(1 + slope ** 2)
        self.t_wall = np.c_[slope, np.ones(nz)] / nrm[:, None]      # along the surface, upwards
        self.n_wall = np.c_[np.ones(nz), -slope] / nrm[:, None]     # into the biofilm
        Ddev = np.zeros((4, 4)); Ddev[:3, :3] = 2 * (np.eye(3) - 1 / 3); Ddev[3, 3] = 1.0
        Dvol = np.zeros((4, 4)); Dvol[:3, :3] = 1.0
        self.Kdev = np.zeros((ne, 8, 8)); self.Kvol = np.zeros((ne, 8, 8)); self.fv = np.zeros((ne, 8))
        self.B0 = np.zeros((ne, 4, 8)); self.vol = np.zeros(ne)
        kd, md = np.zeros((ne, 4, 4)), np.zeros((ne, 4))
        for e in range(ne):
            X = self.xn[self.conn[e]]

            def bmat(xi, eta):
                N, dN = shape(xi, eta)
                J = dN.T @ X
                dx = dN @ np.linalg.inv(J).T
                rr = N @ X[:, 0]
                b = np.zeros((4, 8))
                b[0, 0::2] = dx[:, 0]; b[1, 1::2] = dx[:, 1]; b[2, 0::2] = N / rr
                b[3, 0::2] = dx[:, 1]; b[3, 1::2] = dx[:, 0]
                return b, rr * np.linalg.det(J), N, dx
            for a in (-G, G):
                for c in (-G, G):
                    B, w, N, dx = bmat(a, c)
                    self.Kdev[e] += B.T @ Ddev @ B * w
                    kd[e] += dx @ dx.T * w
                    md[e] += N * w
                    self.vol[e] += w
            B0, w0, _, _ = bmat(0.0, 0.0)
            self.B0[e] = B0
            self.Kvol[e] = B0.T @ Dvol @ B0 * 4 * w0
            self.fv[e] = B0.T @ np.array([1, 1, 1, 0.0]) * 4 * w0
        r4 = np.repeat(self.conn, 4, axis=1).ravel(); c4 = np.tile(self.conn, (1, 4)).ravel()
        self.Kphi = sp.csr_matrix((kd.ravel(), (r4, c4)), shape=(nn, nn))
        self.Mphi = np.zeros(nn); np.add.at(self.Mphi, self.conn.ravel(), md.ravel())
        self.md = md
        self.lam = np.max(np.abs(self.Kphi).sum(1).A1 / self.Mphi)   # Gershgorin bound of M^-1 K

    def phi_e(self, phin):
        return phin[self.conn].mean(1)

    def source(self, alpha):
        f = np.zeros(self.nn); np.add.at(f, self.conn.ravel(), (self.md * alpha[:, None]).ravel())
        return f / self.Mphi


def mechanics(M, phi, delta):
    Ee = E0 * (phi ** 2 + FLOOR)
    mu = Ee / (2 * (1 + NU)); K = Ee / (3 * (1 - 2 * NU))
    vals = (mu[:, None, None] * M.Kdev + K[:, None, None] * M.Kvol).ravel()
    Kg = sp.csr_matrix((vals, (M.rows, M.cols)), shape=(M.ndof, M.ndof))
    f = np.zeros(M.ndof)
    np.add.at(f, M.dofs.ravel(), ((3 * K * delta)[:, None] * M.fv).ravel())
    u = np.zeros(M.ndof)
    u[M.free] = spla.spsolve(Kg[M.free][:, M.free].tocsc(), f[M.free])
    e = np.einsum("eij,ej->ei", M.B0, u[M.dofs])          # rr, zz, tt, rz (engineering)
    ee = e[:, :3] - delta[:, None]
    tr = ee.sum(1)
    s = 2 * mu[:, None] * (ee - tr[:, None] / 3) + K[:, None] * tr[:, None]
    srz = mu * e[:, 3]
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2) + 3 * srz ** 2)
    # traction on the bonded face (first element ring, ordered in z)
    w = M.wall
    S = np.stack([np.stack([s[w, 0], srz[w]], -1), np.stack([srz[w], s[w, 1]], -1)], -2)
    tn = np.einsum("ei,eij,ej->e", M.n_wall, S, M.n_wall)
    ts = np.einsum("ei,eij,ej->e", M.t_wall, S, M.n_wall)
    return -K * tr, vm, s[w, 2], tn, ts


def run(bulge, P_h, local=False, beta=0.02, dt=0.02, t_end=1.0):
    M = Mesh(bulge)
    phin = M.seed_n.copy()
    alpha = np.ones(len(M.conn))
    nsub = int(np.ceil(dt / (1.8 / (beta * M.lam))))
    V = M.vol / M.vol.sum()
    hist = []
    for k in range(int(round(t_end / dt))):
        phi = M.phi_e(phin)
        p, vm, stt, tn, ts = mechanics(M, phi, alpha - 1)
        if P_h is None:
            g = np.ones_like(p)
        else:
            p_h = P_h * E0 * K_ALPHA * ((phi ** 2 + FLOOR) if local else 1.0)
            g = np.clip(1 - p / p_h, 0.0, 1.0)
        src = K_ALPHA * M.source(alpha)
        for _ in range(nsub):
            phin = np.clip(phin + dt / nsub * (-beta * (M.Kphi @ phin) / M.Mphi + src), 0, 1)
        alpha = alpha + dt * K_ALPHA * M.phi_e(phin) * g
        hist.append(((k + 1) * dt, float(V @ (alpha - 1)), float(V @ M.phi_e(phin)), float(p.max()),
                     float(vm.max()), float(stt.min()), float(tn.max()), float(np.abs(ts).max()), float(V @ g)))
    zc = (np.arange(M.nz) + 0.5) * HZ / M.nz
    prof = dict(z=zc.tolist(), alpha=(alpha - 1).reshape(M.nr, M.nz)[0].tolist(),
                tn=tn.tolist(), ts=ts.tolist(), stt=stt.tolist())
    return dict(bulge=bulge, P_h=P_h, local=local, beta=beta, hist=hist, wall=prof)


def main():
    t0 = time.time()
    out = {"scale_E_k_alpha": E0 * K_ALPHA, "runs": []}
    for beta in (0.02, 1e-3):
        for bulge in (0.0, 1.0):
            base = None
            for local, P_h in ((False, None), (False, 0.1), (True, 0.1)):
                r = run(bulge, P_h, local, beta)
                out["runs"].append(r)
                t, d, ph, px, vm, stt, tn, ts, g = r["hist"][-1]
                base = base or d
                print(f"beta={beta:g} {'tooth  ' if bulge else 'implant'} {'local' if local else 'E    '} "
                      f"P_h={P_h}: alpha-1 {d:.3e} ({d / base:.2f}), phi {ph:.3f}, p max {px:.2e}, "
                      f"vM {vm:.2e}, hoop {stt:.2e}, normal (tension) {tn:.2e}, shear {ts:.2e} Pa, g {g:.2f} "
                      f"[{time.time() - t0:.0f} s]", flush=True)
    (HERE / "results_tooth.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
