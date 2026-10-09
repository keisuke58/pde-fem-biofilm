#!/usr/bin/env python3
"""NEM stability on the finite cube and the 8^3 seed runs (Appendix D).

nem_stability.py gives the limit on the infinite integration-point lattice.
This script builds the operator on the finite 2 mm cube (every point fits its
stencil over the 30 nearest points that exist, deck weights), takes its
spectral radius, and marches the one-species seed problem of the element's
8^3 limit runs (beta = 0.1 mm^2/T*, k_alpha = 1e-3, penalty 5, T* = 1.1) with
explicit Euler at their time steps. On 8^3 the limit is 0.063, but phi starts
constant in each element and its component along the unstable mode is about
1e-6, so within the 18 to 22 steps of a run the failure shows only above the
limit (ANSYS: smooth at 0.088, stopped at 0.096).

    python ansys_usermat/apdl/nem_finite_cube.py
"""
import sys, itertools, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nem_stability as N
from scipy.spatial import cKDTree
import scipy.sparse as sp, scipy.sparse.linalg as spla

def finite(n, nneigh=30, bs=0.2, thr=1e-8):
    h = 2.0/n
    ax = []
    for i in range(n):
        for o in N.OFF: ax.append(-1 + (i + o)*h)
    ax = np.array(ax)
    P = np.array(list(itertools.product(ax, ax, ax)))
    tree = cKDTree(P)
    rows, cols, vals = [], [], []
    d_, idx = tree.query(P, nneigh + 1)
    for i in range(len(P)):
        nb = idx[i, 1:]
        D = P[nb] - P[i]
        w = np.exp(-0.5*(4*np.linalg.norm(D, axis=1)/bs)**2); w[w <= thr] = 0
        dx, dy, dz = D.T
        A = np.column_stack([0.5*dx**2, 0.5*dy**2, 0.5*dz**2, dx*dy, dx*dz, dy*dz, dx, dy, dz])
        M = A.T @ (w[:, None]*A)
        try:
            Dm = np.linalg.solve(M, A.T*w)
        except np.linalg.LinAlgError:
            Dm = np.linalg.lstsq(M, A.T*w, rcond=None)[0]
        c = Dm[0] + Dm[1] + Dm[2]
        rows += [i]*len(nb) + [i]; cols += list(nb) + [i]; vals += list(c) + [-c.sum()]
    L = sp.csr_matrix((vals, (rows, cols)), shape=(len(P),)*2)
    ev = spla.eigs(L, k=6, which='LM', return_eigenvectors=False, maxiter=20000)
    rho = np.abs(ev).max()
    lead = ev[np.argmax(np.abs(ev))]
    return 2/(rho*h*h), lead*h*h, P, L

for n in (8, 16):
    lm, lead, P, L = finite(n)
    print(f"{n}^3 finite cube: lambda_max {lm:.4f}  (leading eigenvalue * h^2 = {lead:.3f})", flush=True)

# explicit Euler on the finite 8^3 operator, seed problem with beta = 0.1, as the ANSYS 8^3 limit runs
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp3_reference as W
cells8 = W.seed_cells_8(Path(__file__).resolve().parent / 'results/2026-10-wp3_fix/wp3c_n8_dt00125.json')
lm, lead, P, L = finite(8)
h = 0.25
cell = np.floor((P + 1)/h).astype(int)
seed = np.zeros((8, 8, 8), bool); seed[cells8[:, 0], cells8[:, 1], cells8[:, 2]] = True
inseed = seed[cell[:, 0], cell[:, 1], cell[:, 2]]
vals, vecs = spla.eigs(L, k=1, which='LM')
v = np.real(vecs[:, 0]); v /= np.linalg.norm(v)
phi0 = inseed.astype(float)
print("projection of the initial phi on the leading mode:", abs(v @ phi0)/np.linalg.norm(phi0))
for dt in (0.03, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11):
    beta = 0.1; lam_ = beta*dt/h**2
    phi = phi0.copy(); a = np.ones_like(phi); nst = int(round(1.1/dt)); ok = True
    for s in range(nst):
        pen = 5*(np.maximum(phi-1, 0) + np.minimum(phi, 0))
        phi, a = phi + dt*(beta*(L @ phi) + 1e-3*a - pen), a + dt*1e-3*phi
        if not np.isfinite(phi).all() or np.abs(phi).max() > 1e6: ok = False; break
    print(f"dt {dt:.3f} lambda {lam_:.3f}: steps {s+1}/{nst}, max|phi| {np.abs(phi).max():.3e}, min alpha-1 {(a-1).min():+.2e}")
