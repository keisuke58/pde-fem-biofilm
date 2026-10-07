#!/usr/bin/env python3
"""The NEM Laplacian in a few lines: second-order Taylor expansion over the
nearest neighbours, solved by weighted least squares (Blaszczyk et al. 2022,
Rudolf et al. 2025, Eq. 4-10), on regular and jittered 3D point clouds.

    python ansys_usermat/nem_demo.py -> assets/fig_nem_demo.png

Panel (a): the Laplacian coefficients of one point of a jittered 2D cloud
(same construction, 5 unknown derivatives). Panel (b): the error of the 3D
Laplacian of phi = sin(pi x) sin(pi y) sin(pi z) on the unit cube at the
interior points, against the point spacing h, with 30 neighbours and the
Gaussian weight exp(-1/2 (4 d / beta*)^2), beta* = 4 h. Teaching figure; not
the partner element's code.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_nem_demo.png"
NB = 30


def rows(dx):
    """Taylor rows: second derivatives (unmixed, mixed), then first derivatives."""
    if dx.shape[1] == 3:
        x, y, z = dx.T
        return np.c_[x * x / 2, y * y / 2, z * z / 2, x * y, x * z, y * z, x, y, z], 3
    x, y = dx.T
    return np.c_[x * x / 2, y * y / 2, x * y, x, y], 2


def lap_coeffs(pts, i, idx, bstar):
    """coefficients L_j with lap(phi)(x_i) ~ sum_j L_j (phi_j - phi_i)."""
    dx = pts[idx] - pts[i]
    A, nd = rows(dx)
    w = np.exp(-0.5 * (4 * np.linalg.norm(dx, axis=1) / bstar) ** 2)
    D = np.linalg.solve(A.T @ (w[:, None] * A), A.T * w)    # m = D b
    return D[:nd].sum(axis=0)                                 # rows 1..nd: unmixed second derivatives


def cloud(n, jitter, dim, rng):
    g = (np.arange(n) + 0.5) / n
    P = np.stack(np.meshgrid(*([g] * dim), indexing="ij"), -1).reshape(-1, dim)
    return P + jitter * (rng.random(P.shape) - 0.5) / n


def error(n, jitter, rng):
    P = cloud(n, jitter, 3, rng)
    h = 1.0 / n
    f = lambda p: np.prod(np.sin(np.pi * p), axis=1)
    phi = f(P)
    inner = np.where(np.all((P > 0.2) & (P < 0.8), axis=1))[0]
    _, nbr = cKDTree(P).query(P[inner], NB + 1)
    err = []
    for i, idx in zip(inner, nbr):
        idx = idx[1:]
        L = lap_coeffs(P, i, idx, 4 * h)
        err.append(L @ (phi[idx] - phi[i]) - (-3 * np.pi ** 2 * phi[i]))
    return np.sqrt(np.mean(np.square(err))) / (3 * np.pi ** 2)


def main():
    figstyle.apply(size=10)
    rng = np.random.default_rng(1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))

    P = cloud(9, 0.6, 2, rng); i = len(P) // 2
    _, idx = cKDTree(P).query(P[i], 13); idx = idx[1:]
    L = lap_coeffs(P, i, idx, 4 / 9)
    ax[0].scatter(*P.T, s=8, c="0.75")
    sc = ax[0].scatter(*P[idx].T, s=60, c=L, cmap="Reds", vmin=0,
                       edgecolors="k", linewidths=0.4)
    ax[0].scatter(*P[i], s=90, marker="*", c="k")
    ax[0].set_aspect("equal"); ax[0].set_xticks([]); ax[0].set_yticks([])
    ax[0].set_title("(a) Laplacian weights of one point (2D, 12 neighbours)", fontsize=10)
    fig.colorbar(sc, ax=ax[0], fraction=0.046, label=r"weight $L_j$")

    ns = [8, 12, 16, 24, 32]
    for jit, mk, lab in ((0.0, "o-", "regular grid"), (0.6, "s-", "jittered points")):
        e = [error(n, jit, rng) for n in ns]
        ax[1].loglog([1 / n for n in ns], e, mk, label=lab, ms=5)
        print(lab, ["%.2e" % v for v in e])
    hh = np.array([1 / 32, 1 / 8])
    ax[1].loglog(hh, 0.5 * hh, "k:", lw=1, label=r"slope 1")
    ax[1].loglog(hh, 2.0 * hh ** 2, "k--", lw=1, label=r"slope 2")
    ax[1].set(xlabel=r"point spacing $h$", ylabel=r"relative RMS error of $\nabla^2\phi$")
    ax[1].set_title("(b) 3D, 30 neighbours, unit cube", fontsize=10)
    ax[1].legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(OUT, dpi=200); print("wrote", OUT)


if __name__ == "__main__":
    main()
