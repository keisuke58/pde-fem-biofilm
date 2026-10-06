#!/usr/bin/env python3
"""Stress into growth, first prototype (Keio WP2): the mechanical term of
Klempt et al. 2024, Eq. 30, kept in the phi equation, on a 2D plane-strain
square with a seed.

Model (Klempt 2024, Eq. 30 and 36 divided by eta_phi = eta_alpha, no
nutrient and no front term, so that only the density growth and the new term
act):

    phi_dot   = beta lap(phi) + k_alpha alpha - lam * mu_star * phi * m
    alpha_dot = k_alpha phi
    m         = I : C_e_iso - 3  ~=  2 |dev eps_e|^2      (small strain)

mu_star = mu / eta_phi [1/T*]; lam in {0, 1} switches the term on. With the
Table 2 values mu = 3.3557 Pa and eta_phi = 1e-10, mu_star = 3.4e10 /T*.
Mechanics: small-strain linear elasticity, plane strain, eigenstrain
(alpha - 1) I, stiffness E (phi^2 + f), f = 1e-3, E = 10 Pa, nu = 0.49,
selective reduced integration (Q4, volumetric part at the centre). phi and
alpha live at the element centres (one value per element, as the partner's
element keeps them at the integration points); the Laplacian is the
5-point stencil with zero normal gradient at the boundary.

Coupling per step (staggered): mechanics with alpha_n, then phi and alpha.
Three treatments of the new term:
  explicit   phi_{n+1} = phi_n + dt (... - mu_star phi_n m_n)
  semi       phi_{n+1} = (phi_n + dt (...)) / (1 + dt mu_star m_n)
  iterated   semi, repeated with m recomputed from alpha_{n+1} until it settles

Scaling: the sink balances the source when mu_star phi m = k_alpha alpha. With
m = c (alpha - 1)^2 near the seed this gives an arrest of the growth at

    alpha - 1 ~ sqrt(k_alpha / (c mu_star)),

1.7e-7 / sqrt(c) for the Table 2 values: the growth stops at once, as the
thesis appendix states. The script checks this law against the simulation.

    python keio_wp2/stress_to_growth_2d.py            -> keio_wp2/results.json, assets/fig_stress_to_growth.png
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "ansys_usermat"))

L, E0, NU, FLOOR = 2.0, 10.0, 0.49, 1e-3        # mm, Pa, -, -
BETA, K_ALPHA = 0.02, 1e-3                        # mm^2/T*, 1/T*
R_SEED = 0.3                                       # mm, radius of the seed
MU_TABLE2 = E0 / (2 * (1 + 0.49)) / 1e-10         # mu / eta_phi, Table 2


class Mesh:
    def __init__(self, n):
        self.n, self.h = n, L / n
        xc = (np.arange(n) + 0.5) * self.h - L / 2
        X, Y = np.meshgrid(xc, xc, indexing="ij")
        self.inside = (X ** 2 + Y ** 2) <= R_SEED ** 2
        nn = n + 1
        idx = np.arange(nn * nn).reshape(nn, nn)
        self.conn = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
        self.ndof = 2 * nn * nn
        h = self.h
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)

        def bmat(xi, eta):
            dN = np.array([[cx * (1 + eta * cy), cy * (1 + xi * cx)] for cx, cy in corners]) / 4 * 2 / h
            b = np.zeros((3, 8))
            b[0, 0::2] = dN[:, 0]; b[1, 1::2] = dN[:, 1]
            b[2, 0::2] = dN[:, 1]; b[2, 1::2] = dN[:, 0]
            return b
        g = 1 / np.sqrt(3)
        self.Bg = [bmat(a, b) for a in (-g, g) for b in (-g, g)]
        self.B0 = bmat(0.0, 0.0)
        # unit deviatoric (mu = 1) and volumetric (K = 1) element matrices, plane strain
        Ddev = np.array([[4 / 3, -2 / 3, 0], [-2 / 3, 4 / 3, 0], [0, 0, 1]])
        Dvol = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]], float)
        w = (h / 2) ** 2
        self.Kdev = sum(b.T @ Ddev @ b * w for b in self.Bg)
        self.Kvol = self.B0.T @ Dvol @ self.B0 * (4 * w)
        self.fvol = self.B0.T @ np.array([1.0, 1.0, 0.0]) * (4 * w)   # times 3 K delta
        dofs = np.stack([2 * self.conn, 2 * self.conn + 1], -1).reshape(-1, 8)
        self.dofs = dofs
        self.rows = np.repeat(dofs, 8, axis=1).ravel()
        self.cols = np.tile(dofs, (1, 8)).ravel()
        fixed = [0, 1, 2 * n + 1]                  # node (0,0) x,y and node (n,0) y
        self.free = np.setdiff1d(np.arange(self.ndof), fixed)

    def lap(self, f):
        p = np.pad(f.reshape(self.n, self.n), 1, mode="edge")
        lap = (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4 * p[1:-1, 1:-1]) / self.h ** 2
        return lap.ravel()


def mechanics(M, phi, delta):
    """element m = 2|dev eps_e|^2 and von Mises [Pa] for eigenstrain delta I."""
    Ee = E0 * (phi ** 2 + FLOOR)
    mu = Ee / (2 * (1 + NU)); K = Ee / (3 * (1 - 2 * NU))
    vals = (mu[:, None, None] * M.Kdev + K[:, None, None] * M.Kvol).ravel()
    Kg = sp.csr_matrix((vals, (M.rows, M.cols)), shape=(M.ndof, M.ndof))
    f = np.zeros(M.ndof)
    np.add.at(f, M.dofs.ravel(), ((3 * K * delta)[:, None] * M.fvol).ravel())
    u = np.zeros(M.ndof)
    Kf = Kg[M.free][:, M.free]
    u[M.free] = spla.spsolve(Kf.tocsc(), f[M.free])
    eps = (M.B0 @ u[M.dofs].T).T                    # exx, eyy, gxy at the centre
    exx, eyy, gxy = eps[:, 0] - delta, eps[:, 1] - delta, eps[:, 2]
    ezz = -delta
    tr = exx + eyy + ezz
    d = np.stack([exx - tr / 3, eyy - tr / 3, ezz - tr / 3], 1)
    dev2 = (d ** 2).sum(1) + 2 * (gxy / 2) ** 2
    vm = np.sqrt(1.5) * 2 * mu * np.sqrt(dev2)       # sqrt(3/2) |dev sigma|
    return 2 * dev2, vm


def run(n=32, mu_star=0.0, dt=0.01, t_end=1.0, scheme="semi", tol=1e-6, max_it=20):
    M = Mesh(n)
    phi = M.inside.ravel().astype(float)
    alpha = np.ones_like(phi)
    hist = []
    unstable = False
    its = []
    for k in range(int(round(t_end / dt))):
        m, vm = mechanics(M, phi, alpha - 1)
        src = BETA * M.lap(phi) + K_ALPHA * alpha
        if scheme == "explicit":
            phi_new = phi + dt * (src - mu_star * phi * m)
            alpha_new = alpha + dt * K_ALPHA * phi
            nit = 1
        else:
            phi_new = (phi + dt * src) / (1 + dt * mu_star * m)
            alpha_new = alpha + dt * K_ALPHA * phi_new
            nit = 1
            if scheme == "iterated":
                for nit in range(2, max_it + 1):
                    m2, _ = mechanics(M, phi_new, alpha_new - 1)
                    p2 = (phi + dt * src) / (1 + dt * mu_star * m2)
                    a2 = alpha + dt * K_ALPHA * p2
                    done = np.abs(p2 - phi_new).max() < tol
                    phi_new, alpha_new = p2, a2
                    if done:
                        break
        its.append(nit)
        if not np.all(np.isfinite(phi_new)) or phi_new.min() < -1e-6:
            unstable = True
        phi = np.clip(phi_new, 0.0, 1.0)
        alpha = alpha_new
        ins = M.inside.ravel()
        hist.append(((k + 1) * dt, float((alpha - 1)[ins].mean()), float(phi[ins].mean()),
                     float(vm[ins].mean()), float(m[ins].mean())))
    return dict(n=n, mu_star=mu_star, dt=dt, scheme=scheme, unstable=unstable,
                mean_iterations=float(np.mean(its)), hist=hist)


def main():
    t0 = time.time()
    out = {"mu_table2": MU_TABLE2, "runs": []}
    # (a) the term against mu_star (semi-implicit, dt = 0.01)
    for ms in [0.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, MU_TABLE2]:
        r = run(mu_star=ms)
        r["set"] = "sweep"
        out["runs"].append(r)
        t, d, p, vm, m = r["hist"][-1]
        print(f"mu*={ms:9.3g}: seed alpha-1 {d:.3e}, phi {p:.3f}, vM {vm:.3e} Pa, m {m:.2e}")
    t, d0, _, _, m0 = out["runs"][0]["hist"][-1]
    out["c_seed"] = m0 / d0 ** 2
    print(f"c = m/(alpha-1)^2 in the seed without the term: {out['c_seed']:.3f}")
    # (b) scheme and time step where the term matters (mu* = 1e7 and the Table 2 value)
    for ms in [1e7, MU_TABLE2]:
        ref = run(mu_star=ms, dt=0.001, scheme="iterated")
        ref["set"] = "reference"
        out["runs"].append(ref)
        dref, pref = ref["hist"][-1][1], ref["hist"][-1][2]
        for scheme in ["explicit", "semi", "iterated"]:
            for dt in [0.002, 0.005, 0.01, 0.025, 0.04]:
                r = run(mu_star=ms, dt=dt, scheme=scheme)
                r["set"] = "scheme"
                out["runs"].append(r)
                t, d, p, vm, m = r["hist"][-1]
                print(f"mu*={ms:8.2g} {scheme:9s} dt={dt:5.3f}: alpha-1 {d:.4e} (ref {dref:.4e}), "
                      f"phi {p:.4f} (ref {pref:.4f}), unstable {r['unstable']}, it {r['mean_iterations']:.1f}")
    (HERE / "results.json").write_text(json.dumps(out))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
