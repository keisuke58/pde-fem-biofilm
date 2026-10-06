#!/usr/bin/env python3
"""Keio P1 (numerics paper), first step: the growth field of Klempt et al. 2024
without front term and nutrient,

    phi_dot = beta lap(phi) + k_alpha alpha,   alpha_dot = k_alpha phi,
    zero normal gradient on the boundary, phi(0) = 1 in a disc of radius 0.3 mm,

on the 2 mm square, solved two ways:

  NEM       phi and alpha at the integration points (2x2 per square element, as
            in the partner's element), Laplacian from the neighbours by a
            second-order Taylor expansion and weighted least squares
            (Rudolf et al. 2025, Eq. 4-13; weight exp(-1/2 (4 d / beta*)^2),
            beta* = 2 h, 12 neighbours; boundary points get the row
            n . grad phi = 0), explicit Euler.
  Galerkin  phi nodal, bilinear Q4, consistent mass, backward Euler; alpha at
            the integration points (Eq. 36 as in the material routine).

Measured at T* = 1 against a fine Galerkin reference (n = 256, dt = 0.0025):
  - mean alpha - 1 in the seed (what the stress sees),
  - the step of alpha - 1 across the seed surface (largest jump between
    neighbouring points, normalised by the seed mean),
  - L2 difference of phi.
The time step of NEM is limited by the explicit update; Galerkin is not.

    python keio_p1/nem_vs_galerkin_2d.py -> keio_p1/results_2d.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
L, R, K_ALPHA, T_END = 2.0, 0.3, 1e-3, 1.0


def gauss_points(n):
    h = L / n
    g = h / 2 * (1 - 1 / np.sqrt(3))
    c = np.arange(n) * h
    xs = np.sort(np.concatenate([c + g, c + h - g]))
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    return np.c_[X.ravel(), Y.ravel()] - L / 2, h


def nem_operator(P, h, nb=12, bstar_fac=2.0):
    """sparse Laplacian with zero normal gradient on the boundary rows."""
    tree = cKDTree(P)
    _, idx = tree.query(P, nb + 1)
    rows, cols, vals = [], [], []
    lo, hi = P.min(0), P.max(0)
    for i, nbr in enumerate(idx):
        nbr = nbr[1:]
        d = P[nbr] - P[i]
        A = np.c_[d[:, 0] ** 2 / 2, d[:, 1] ** 2 / 2, d[:, 0] * d[:, 1], d[:, 0], d[:, 1]]
        w = np.exp(-0.5 * (4 * np.linalg.norm(d, axis=1) / (bstar_fac * h)) ** 2)
        nrm = np.zeros(2)
        nrm[0] = -1 if P[i, 0] <= lo[0] + 1e-12 else (1 if P[i, 0] >= hi[0] - 1e-12 else 0)
        nrm[1] = -1 if P[i, 1] <= lo[1] + 1e-12 else (1 if P[i, 1] >= hi[1] - 1e-12 else 0)
        Wt = w
        if nrm.any():                                   # Neumann row (Eq. 11-13)
            A = np.vstack([A, [0, 0, 0, nrm[0], nrm[1]]])
            Wt = np.r_[w, 10 * w.max()]
        D = np.linalg.solve(A.T @ (Wt[:, None] * A), A.T * Wt)
        Lrow = D[0] + D[1]
        coef = Lrow[:nb]                                # the Neumann row has b = 0
        rows += [i] * nb + [i]
        cols += list(nbr) + [i]
        vals += list(coef) + [-coef.sum()]
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(P), len(P)))


def run_nem(n, beta, dt):
    P, h = gauss_points(n)
    Lap = nem_operator(P, h)
    phi = (np.hypot(P[:, 0], P[:, 1]) <= R).astype(float)
    alpha = np.ones_like(phi)
    for _ in range(int(round(T_END / dt))):
        phi_new = phi + dt * (beta * (Lap @ phi) + K_ALPHA * alpha)
        alpha = alpha + dt * K_ALPHA * phi
        phi = phi_new
        if not np.isfinite(phi).all() or np.abs(phi).max() > 10:
            return None
    return P, phi, alpha


def q4(n):
    h = L / n
    nn = n + 1
    idx = np.arange(nn * nn).reshape(nn, nn)
    conn = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 4)
    xi = np.array([-1, 1]) / np.sqrt(3)
    corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    Me, Ke, Ng = np.zeros((4, 4)), np.zeros((4, 4)), []
    for a in xi:
        for b in xi:
            N = (1 + a * corners[:, 0]) * (1 + b * corners[:, 1]) / 4
            dN = np.c_[corners[:, 0] * (1 + b * corners[:, 1]), corners[:, 1] * (1 + a * corners[:, 0])] / 4 * 2 / h
            w = (h / 2) ** 2
            Me += np.outer(N, N) * w
            Ke += dN @ dN.T * w
            Ng.append(N)
    r = np.repeat(conn, 4, 1).ravel(); c = np.tile(conn, (1, 4)).ravel()
    M = sp.csr_matrix((np.tile(Me.ravel(), len(conn)), (r, c)), shape=(nn * nn, nn * nn))
    K = sp.csr_matrix((np.tile(Ke.ravel(), len(conn)), (r, c)), shape=(nn * nn, nn * nn))
    Fe = sum(N * (h / 2) ** 2 for N in Ng)               # integral of N_a
    F = np.zeros(nn * nn)
    np.add.at(F, conn.ravel(), np.tile(Fe, len(conn)))
    Ng = np.array(Ng)                                     # 4 gauss x 4 nodes
    xs = np.linspace(-L / 2, L / 2, nn)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    return dict(M=M, K=K, conn=conn, Ng=Ng, Xn=np.c_[X.ravel(), Y.ravel()], h=h, Fe=Fe)


def run_galerkin(n, beta, dt):
    S = q4(n)
    phi = (np.hypot(S["Xn"][:, 0], S["Xn"][:, 1]) <= R).astype(float)
    ne = len(S["conn"])
    alpha = np.ones((ne, 4))
    A = (S["M"] + dt * beta * S["K"]).tocsc()
    lu = spla.splu(A)
    w = (S["h"] / 2) ** 2
    for _ in range(int(round(T_END / dt))):
        # source k_alpha alpha integrated with the gauss values of alpha
        src = np.zeros(len(phi))
        np.add.at(src, S["conn"].ravel(), (K_ALPHA * alpha[:, :, None] * S["Ng"][None] * w).sum(1).ravel())
        phi = lu.solve(S["M"] @ phi + dt * src)
        phig = phi[S["conn"]] @ S["Ng"].T                 # ne x 4 gauss values
        alpha = alpha + dt * K_ALPHA * phig
    # integration point coordinates
    xi = np.array([-1, 1]) / np.sqrt(3)
    gp = []
    for a in xi:
        for b in xi:
            N = (1 + a * np.array([-1, 1, 1, -1])) * (1 + b * np.array([-1, -1, 1, 1])) / 4
            gp.append(S["Xn"][S["conn"]].transpose(0, 2, 1) @ N)
    P = np.stack(gp, 1).reshape(-1, 2)
    return P, phi, alpha.reshape(-1), S


def metrics(P, phi_p, alpha_p, ref):
    """phi_p, alpha_p at the points P; ref: (tree, phi, alpha) of the reference points."""
    seed = np.hypot(P[:, 0], P[:, 1]) <= R
    d = alpha_p - 1
    tree = cKDTree(P)
    pairs = tree.query_pairs(1.01 * np.min(cKDTree(P).query(P, 2)[0][:, 1]), output_type="ndarray")
    step = np.abs(d[pairs[:, 0]] - d[pairs[:, 1]]).max() / d[seed].mean()
    rtree, rphi, ralpha = ref
    _, j = rtree.query(P)
    l2 = np.sqrt(np.mean((phi_p - rphi[j]) ** 2))
    return dict(seed_alpha=float(d[seed].mean()), step=float(step), l2_phi=float(l2),
                err_alpha=float(abs(d[seed].mean() - (ralpha[j][seed] - 1).mean()) / (ralpha[j][seed] - 1).mean()))


def main():
    t0 = time.time()
    out = {"runs": []}
    for beta in (0.02, 1e-4):
        P, phi, alpha, S = run_galerkin(256, beta, 0.0025)
        phig = phi[S["conn"]] @ S["Ng"].T
        ref = (cKDTree(P), phig.reshape(-1), alpha)
        print(f"beta={beta}: reference done [{time.time() - t0:.0f} s]", flush=True)
        for n in (16, 32, 64, 128):
            h = L / n
            # NEM: the largest stable step in 2D is about h_gp^2/(4 beta); use 0.4 of it and 0.025 if larger
            hg = h / 2
            dt_lim = hg ** 2 / (4 * beta)
            for method, dt in (("nem", min(0.025, 0.4 * dt_lim)), ("nem_dt0.025", 0.025),
                               ("galerkin", 0.025)):
                t1 = time.time()
                if method.startswith("nem"):
                    if method == "nem_dt0.025" and 0.025 <= 0.4 * dt_lim:
                        continue
                    r = run_nem(n, beta, dt)
                    if r is None:
                        out["runs"].append(dict(beta=beta, n=n, method=method, dt=dt, unstable=True))
                        print(f"  n={n:3d} {method:12s} dt={dt:.2e}: unstable", flush=True)
                        continue
                    Pm, ph, al = r
                else:
                    Pm, phn, al, Sg = run_galerkin(n, beta, dt)
                    ph = (phn[Sg["conn"]] @ Sg["Ng"].T).reshape(-1)
                m = metrics(Pm, ph, al, ref)
                m.update(beta=beta, n=n, method=method, dt=dt, unstable=False, seconds=time.time() - t1)
                out["runs"].append(m)
                print(f"  n={n:3d} {method:12s} dt={dt:.2e}: seed alpha-1 {m['seed_alpha']:.4e} "
                      f"(err {m['err_alpha']:.2e}), step {m['step']:.3f}, L2 phi {m['l2_phi']:.2e}, "
                      f"{m['seconds']:.1f} s", flush=True)
    (HERE / "results_2d.json").write_text(json.dumps(out))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
