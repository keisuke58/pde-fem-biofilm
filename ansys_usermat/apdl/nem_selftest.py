#!/usr/bin/env python3
"""NEM operator self-test on a model of the 8^3 / 16^3 point set: 2x2x2 Gauss points
per hex + mesh nodes on the cube faces as surface points. Interior operator: weighted
least squares on the 9 Taylor terms (deck weights), stencil = nneigh nearest of all points.
Surface values: exact, IDW mean of GP neighbours (the 10/07 fix), or held at 0.
(CLOUD_TO_IKMHIWI03_2026-10-13 item 4; results in IKMHIWI03_TO_CLOUD.ja.md.)

    python ansys_usermat/apdl/nem_selftest.py
"""
import itertools, sys, numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nem_stability as N

def build(n, nneigh=30, bs=0.2, thr=1e-8):
    h = 2.0 / n
    ax = np.array([-1 + (i + o) * h for i in range(n) for o in N.OFF])
    G = np.array(list(itertools.product(ax, ax, ax)))
    nd = np.linspace(-1, 1, n + 1)
    S = np.array([p for p in itertools.product(nd, nd, nd) if np.any(np.abs(np.array(p)) == 1)])
    P = np.vstack([G, S]); nG = len(G)
    tree = cKDTree(P)
    _, idx = tree.query(P, nneigh + 1)
    ops = []
    for i in range(nG):
        nb = idx[i, 1:]; D = P[nb] - P[i]
        w = np.exp(-0.5 * (4 * np.linalg.norm(D, axis=1) / bs) ** 2); w[w <= thr] = 0
        dx, dy, dz = D.T
        A = np.column_stack([0.5*dx**2, 0.5*dy**2, 0.5*dz**2, dx*dy, dx*dz, dy*dz, dx, dy, dz])
        Dm = np.linalg.solve(A.T @ (w[:, None] * A), A.T * w)
        ops.append((nb, Dm[0] + Dm[1] + Dm[2], Dm[6:9]))
    # IDW surface update (GP neighbours only)
    idw = []
    for j in range(nG, len(P)):
        nb = idx[j, 1:]; nb = nb[nb < nG]
        d2 = np.sum((P[nb] - P[j]) ** 2, axis=1); idw.append((nb, 1 / d2 / np.sum(1 / d2)))
    # distance of each GP to the boundary
    dist = 1 - np.max(np.abs(G), axis=1)
    near = np.array([np.any(nb >= nG) for nb, _, _ in ops])
    return P, nG, ops, idw, near, h

def test(n, nneigh=30):
    P, nG, ops, idw, near, h = build(n, nneigh)
    x, y, z = P.T
    fields = {"linear u=x+2y-z": (x + 2*y - z, lambda q: np.tile([1, 2, -1], (len(q), 1)), 0.0),
              "quadratic u=x^2+yz": (x**2 + y*z, lambda q: np.column_stack([2*q[:, 0], q[:, 2], q[:, 1]]), 2.0)}
    print(f"== {n}^3, {nneigh} neighbours: {nG} GPs, {len(P)-nG} surface points, "
          f"{near.sum()} GPs with a surface point in the stencil")
    for name, (u, grad, lap) in fields.items():
        for mode in ("exact", "IDW", "zero"):
            v = u.copy()
            if mode == "IDW":
                for k, (nb, w) in enumerate(idw): v[nG + k] = w @ u[nb]
            if mode == "zero": v[nG:] = 0
            eL = np.array([c @ (v[nb] - v[i]) - lap for i, (nb, c, g) in enumerate(ops)])
            eG = np.array([np.max(np.abs(g @ (v[nb] - v[i]) - grad(P[i:i+1])[0])) for i, (nb, c, g) in enumerate(ops)])
            sb = np.max(np.abs(v[nG:] - u[nG:]))
            print(f"  {name:20s} surface {mode:5s}: max|u_b err| {sb:.2e} | lap err interior {np.abs(eL[~near]).max():.2e}"
                  f" near-boundary {np.abs(eL[near]).max():.2e} | grad err interior {eG[~near].max():.2e} near {eG[near].max():.2e}")

for n in (8, 16):
    test(n)
