#!/usr/bin/env python3
"""Stability limit of the explicit field update of the IKM biofilm element.

The element stores phi at the 2 x 2 x 2 integration points of SOLID185 and
updates it explicitly, phi_{n+1} = phi_n + dt (beta L phi_n + ...), with L the
NEM Laplacian: a second-order Taylor polynomial fitted by weighted least
squares over the NEIGHBOR_CNT nearest integration points, Gaussian weight
w = exp(-(4 d / BETA_STAR)^2 / 2) cut at WMAT_THRESHOLD (OLIVER_MODEL_NOTES.md,
"What the NEM operator actually is"; the deck values are NEIGHBOR_CNT = 30,
BETA_STAR = 0.2 mm, WMAT_THRESHOLD = 1e-8).

On a uniform hexahedral mesh of size h the integration points form a lattice
with two points per cell and axis, at h (1/2 -/+ 1/(2 sqrt 3)), so the spacing
alternates between h/sqrt(3) inside an element and h (1 - 1/sqrt 3) across the
element boundary. The eight points of a cell are equivalent under the cubic
reflections, so the stencil of one of them gives the others by mirroring.

This script
  1. builds the WLS stencil of one integration point on the infinite lattice,
  2. forms the 8 x 8 Bloch symbol of L over the Brillouin zone of the cell
     lattice and finds its spectral radius rho,
  3. gives the explicit-Euler limit  lambda_max = 2 / (rho h^2),
     lambda = beta dt / h^2  (the 7-point stencil gives 1/6, the 1-D 3-point
     stencil 1/2, which is Eq. 3 of Rudolf et al. 2025),
  4. checks the limit by time stepping on a periodic block of cells.

    python ansys_usermat/apdl/nem_stability.py            -> table on stdout
    python ansys_usermat/apdl/nem_stability.py --json out.json --fig out.png
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

S3 = np.sqrt(3.0)
OFF = (0.5 - 1.0 / (2.0 * S3), 0.5 + 1.0 / (2.0 * S3))   # Gauss points of a unit cell, per axis


def lattice_points(h, ncell):
    """Integration points of a (2 ncell + 1)^3 block of cells around the origin cell."""
    pts = []
    cells = range(-ncell, ncell + 1)
    for cx, cy, cz in itertools.product(cells, repeat=3):
        for ox, oy, oz in itertools.product(OFF, repeat=3):
            pts.append(((cx + ox) * h, (cy + oy) * h, (cz + oz) * h))
    return np.array(pts)


def wls_stencil(x0, pts, nneigh, beta_star, threshold, tie="symmetric"):
    """Laplacian weights c_i of the WLS second-order fit at x0 over the nearest
    nneigh points. Returns (neighbour offsets, weights c_i); the self weight is
    -sum(c_i). tie: how to cut a distance shell that straddles nneigh:
    'symmetric' takes the whole shell, 'first' the first points in array order."""
    d = pts - x0
    r = np.linalg.norm(d, axis=1)
    order = np.argsort(r, kind="stable")
    order = order[r[order] > 1e-12 * max(1.0, r.max())]   # drop the point itself
    if tie == "symmetric":
        rc = r[order[nneigh - 1]]
        sel = order[r[order] <= rc * (1 + 1e-9)]
    else:
        sel = order[:nneigh]
    D = d[sel]
    w = np.exp(-0.5 * (4.0 * np.linalg.norm(D, axis=1) / beta_star) ** 2)
    w[w <= threshold] = 0.0
    dx, dy, dz = D.T
    A = np.column_stack([0.5 * dx**2, 0.5 * dy**2, 0.5 * dz**2,
                         dx * dy, dx * dz, dy * dz, dx, dy, dz])
    W = np.diag(w)
    M = A.T @ W @ A
    Dm = np.linalg.solve(M, A.T @ W)          # 9 x n: derivatives from (v_i - v_0)
    c = Dm[0] + Dm[1] + Dm[2]                 # Laplacian row
    return D, c, w


def bloch_radius(h, D, c, nk=24):
    """Spectral radius of the NEM Laplacian on the infinite lattice.
    The eight points of a cell are mirror images; the stencil of the reference
    point (low, low, low) is mirrored to the others. The symbol is 8 x 8."""
    subs = list(itertools.product(range(2), repeat=3))            # 0 = low, 1 = high per axis
    sub_index = {s: i for i, s in enumerate(subs)}
    off = np.array(OFF) * h
    x_ref = np.array([off[0]] * 3)
    ks = np.linspace(-np.pi / h, np.pi / h, nk, endpoint=False)
    worst = 0.0
    kworst = None
    for kx, ky, kz in itertools.product(ks, repeat=3):
        k = np.array([kx, ky, kz])
        Lk = np.zeros((8, 8), dtype=complex)
        for s in subs:
            i = sub_index[s]
            sgn = np.array([1 if si == 0 else -1 for si in s])   # mirror for 'high' points
            x_s = np.array([off[si] for si in s])
            for dvec, ci in zip(D, c):
                dm = dvec * sgn                                   # offset from point s
                x_n = x_s + dm
                # which sub-lattice point does x_n land on?
                cell = np.floor(x_n / h)
                frac = x_n / h - cell
                t = tuple(0 if abs(f - OFF[0]) < 1e-6 else 1 for f in frac)
                assert all(abs(f - OFF[ti]) < 1e-6 for f, ti in zip(frac, t))
                j = sub_index[t]
                Lk[i, j] += ci * np.exp(1j * k @ dm)
            Lk[i, i] += -np.sum(c)
        ev = np.linalg.eigvals(Lk)
        m = np.max(np.abs(ev))
        if m > worst:
            worst, kworst = m, k
    return worst, kworst


def march_check(h, D, c, lam, ncell=4, nstep=400, seed=0):
    """Explicit Euler on a periodic ncell^3 block; returns the growth factor per step."""
    pts = []
    subs = list(itertools.product(range(2), repeat=3))
    for cell in itertools.product(range(ncell), repeat=3):
        for s in subs:
            pts.append((cell, s))
    idx = {p: i for i, p in enumerate(pts)}
    n = len(pts)
    L = np.zeros((n, n))
    off = np.array(OFF)
    for (cell, s), i in idx.items():
        sgn = np.array([1 if si == 0 else -1 for si in s])
        x_s = (np.array(cell) + np.array([off[si] for si in s])) * h
        for dvec, ci in zip(D, c):
            x_n = x_s + dvec * sgn
            cn = np.floor(x_n / h + 1e-9)
            frac = x_n / h - cn
            t = tuple(0 if abs(f - OFF[0]) < 1e-6 else 1 for f in frac)
            j = idx[(tuple(int(v) % ncell for v in cn), t)]
            L[i, j] += ci
        L[i, i] -= np.sum(c)
    rng = np.random.default_rng(seed)
    u = rng.standard_normal(n)
    u -= u.mean()
    dt_beta = lam * h * h
    g = []
    for _ in range(nstep):
        v = u + dt_beta * (L @ u)
        g.append(np.linalg.norm(v) / np.linalg.norm(u))
        u = v / np.linalg.norm(v)
    return float(np.exp(np.mean(np.log(g[-50:]))))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path)
    ap.add_argument("--fig", type=Path)
    ap.add_argument("--nk", type=int, default=24)
    a = ap.parse_args(argv)

    cases = []
    meshes = [(8, 0.25), (16, 0.125), (24, 2.0 / 24)]
    for n, h in meshes:
        pts = lattice_points(h, 3)
        x0 = np.array([OFF[0] * h] * 3)
        for nneigh, tie in [(26, "symmetric"), (30, "first"), (32, "symmetric"), (56, "symmetric")]:
            for beta_star in [0.2, 1e9]:
                D, c, w = wls_stencil(x0, pts, nneigh, beta_star, 1e-8, tie)
                rho, kw = bloch_radius(h, D, c, a.nk)
                lam = 2.0 / (rho * h * h)
                cases.append(dict(mesh=n, h=h, nneigh=int(len(c)), tie=tie,
                                  beta_star=beta_star, lam_max=lam,
                                  k_worst=(np.array(kw) * h / np.pi).round(3).tolist()))
                print(f"{n:3d}^3 h={h:.4f} n={len(c):2d} {tie:9s} beta*={beta_star:<6g} "
                      f"lambda_max = {lam:.4f}  (k h/pi = {np.round(np.array(kw)*h/np.pi,2)})")
    # reference: 7-point stencil on a simple cubic lattice
    print("7-point finite differences: lambda_max = 1/6 =", 1 / 6)

    # time-stepping check for the deck values on each mesh
    print("\nmarch check (deck values: 30 neighbours, beta* = 0.2 mm), growth per step:")
    checks = []
    for n, h in meshes:
        pts = lattice_points(h, 3)
        x0 = np.array([OFF[0] * h] * 3)
        D, c, w = wls_stencil(x0, pts, 30, 0.2, 1e-8, "first")
        rho, _ = bloch_radius(h, D, c, a.nk)
        lam = 2.0 / (rho * h * h)
        row = dict(mesh=n, lam_max=lam, growth={})
        for f in (0.9, 0.98, 1.02, 1.1):
            g = march_check(h, D, c, f * lam)
            row["growth"][f] = g
            print(f"  {n:3d}^3  lambda = {f:.2f} lambda_max = {f*lam:.4f}: growth {g:.4f}")
        checks.append(row)

    out = dict(cases=cases, march=checks,
               note="lambda = beta dt / h^2; explicit Euler stable iff lambda <= lambda_max")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
        print("wrote", a.json)
    if a.fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import figstyle
        figstyle.apply(11)
        fig, ax = plt.subplots(figsize=(5.2, 3.0))
        for nn, mk in [(26, "o"), (30, "s"), (32, "^"), (56, "D")]:
            xs = [cs["mesh"] for cs in cases if cs["nneigh"] == nn and cs["beta_star"] == 0.2]
            ys = [cs["lam_max"] for cs in cases if cs["nneigh"] == nn and cs["beta_star"] == 0.2]
            if xs:
                ax.plot(xs, ys, marker=mk, label=f"{nn} neighbours")
        ax.axhline(1 / 6, color="k", ls="--", lw=0.9)
        ax.text(24.5, 1 / 6, "1/6", va="bottom", ha="right")
        ax.axhspan(0.144, 0.160, color="0.85", lw=0)
        ax.text(8.2, 0.152, "observed: smooth 0.144, oscillating 0.160", fontsize=8, va="center")
        ax.set(xlabel="elements per edge", ylabel=r"$\lambda_{\max}=\beta\Delta t/h^2$",
               xticks=[8, 16, 24])
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(a.fig, dpi=200)
        print("wrote", a.fig)


if __name__ == "__main__":
    main()
