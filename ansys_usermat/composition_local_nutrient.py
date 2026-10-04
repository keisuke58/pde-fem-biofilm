#!/usr/bin/env python3
"""Composition that varies in space without spreading: the local nutrient in the
point model (step 1 of ROADMAP_TWO_WAY.md), tried in Python on the ANSYS model.

In the coupled ANSYS runs the point model uses a constant nutrient level c*, so
inside the packed seed every Gauss point has the same composition. Klempt et
al. 2026 treat c* as a given, possibly time-dependent, input; replacing it by
the local nutrient c of the field needs no new physics, only a fourth input to
the bridge. This script estimates what that would show, on the partner's
8^3 model (2 mm cube, seed of 32 elements, c = 1 held in the element layer at
y = -1 mm), without ANSYS:

  - steady nutrient d lap c = g phi c (first-order consumption in the seed,
    zero flux elsewhere), cell-centred on the 8^3 mesh; its strength is the
    Thiele number Lambda^2 = g L^2 / d with L = 1 mm. The partner's d and g are
    example inputs, not paper values, so Lambda is scanned;
  - in every seed element the point model (Klempt et al. 2026, cases 3 and 6)
    runs the coupled scheme of the ANSYS runs (amount phi_cap = 0.9, 10 coupling
    steps of 0.1, s = 0.15) with c* = c*_0 * c(x).

    python ansys_usermat/composition_local_nutrient.py --n 32 --dt 0.025
        -> assets/fig_composition_local_nutrient.png

Result (4 Oct 2026), share phi_1/(phi_1+phi_2) in the seed, converged
(32^3 grid, coupling step 0.025; the figure is made with these settings):
  Lambda = 3: case 3 0.652-0.665, case 6 0.024-0.219
  Lambda = 4: case 3 0.627-0.664, case 6 0.025-0.438
Convergence: from 16^3 to 32^3 the case-6 range at Lambda = 4 stays at
0.03-0.44; the coupling step matters more for case 3 (0.60 at 0.1, 0.66 at
0.025, the converged value of the ANSYS runs being 0.667). On the 8^3 ANSYS
mesh with step 0.1 the case-6 range is 0.12-0.44: the interior is resolved,
the edge facing the nutrient is not.
Case 3 (coexistence) depends little on the nutrient level. In case 6 the
takeover by species 2 is slower where the nutrient is low, so the seed's
interior, far from the nutrient face, keeps more of species 1: a composition
that varies in space without any spreading, once the consumption is strong
enough (Lambda >= 3). Lambda depends on the partner's d and g (example input,
not paper values); c* = c*_0 c is an assumption (c normalised to the held 1).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import composition_transport_check as T  # noqa: E402
import figstyle  # noqa: E402
import model_setup_fig as M  # noqa: E402

OUT = HERE.parent / "assets" / "fig_composition_local_nutrient.png"
N, H = 8, 0.25                                   # set by --n
LAMBDAS = [1.0, 2.0, 3.0, 4.0]


def grid_masks():
    """Seed and nutrient layer of the 8^3 deck, mapped by position onto an N^3 grid
    (a cell belongs to them when its centre lies in one of their 8^3 elements)."""
    _, cen, seed, nut, _ = M.read_deck(M.DECK)
    ax = -1 + H * (np.arange(N) + 0.5)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    P = np.stack([X, Y, Z], -1)
    def inside(ids):
        m = np.zeros((N,) * 3, bool)
        for e in ids:
            m |= np.all(np.abs(P - np.asarray(cen[e])) < 0.125 + 1e-9, axis=-1)
        return m
    return inside(seed), inside(nut)


def nutrient(S, D, lam):
    """d lap c = g S c, c = 1 on D, zero flux on the walls; d = 1, g = lam^2 (L = 1 mm)."""
    n = N ** 3
    ix = np.arange(n).reshape((N,) * 3)
    rows, cols, vals = [], [], []
    b = np.zeros(n)
    g = lam ** 2
    for i in range(N):
        for j in range(N):
            for k in range(N):
                r = ix[i, j, k]
                if D[i, j, k]:
                    rows.append(r); cols.append(r); vals.append(1.0); b[r] = 1.0
                    continue
                diag = g * H * H * S[i, j, k]
                for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                    a, bb, cc = i + di, j + dj, k + dk
                    if 0 <= a < N and 0 <= bb < N and 0 <= cc < N:
                        rows.append(r); cols.append(ix[a, bb, cc]); vals.append(-1.0); diag += 1.0
                rows.append(r); cols.append(r); vals.append(diag)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return spla.spsolve(A.tocsc(), b).reshape((N,) * 3)


def shares(case, cvals):
    """phi_1/(phi_1+phi_2) at T* = 1 for every nutrient level in cvals, in one
    vectorised call: c* = c*_0 * c per point, otherwise the coupled scheme of
    composition_transport_check (rescale, advance s*DT, rescale)."""
    import math
    import jax
    import jax.numpy as jnp
    eco = T.eco
    T.ms.set_case(case)
    th, hp = T.ms.ECOLOGY_CASE["theta"], dict(T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    dt_pm = T.S * T.DT
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    mask = eco.active_mask(2)
    tht = jnp.asarray(th, dtype=jnp.float64)
    eta = jnp.asarray(hp["eta"], dtype=jnp.float64)
    steps = jnp.arange(n_sub)
    run = jax.jit(jax.vmap(lambda g, c: eco._substep_scan(
        g, tht, dt_pm / n_sub, steps, mask, c, jnp.float64(hp["alpha"]), eta)[0]))
    cv = np.asarray(cvals, dtype=float)
    G = T.seed_state(cv.size)
    p = np.full(cv.size, T.PHI_CAP)
    cst = jnp.asarray(hp["c"] * cv)
    for _ in range(int(round(1.0 / T.DT))):
        G = T.rescale(G, p)
        G = T.rescale(np.asarray(run(jnp.asarray(G), cst)), p)
    return G[:, 0] / (G[:, 0] + G[:, 1])


def main(argv=None):
    import argparse
    global N, H
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="cells per edge (8 = the ANSYS mesh)")
    ap.add_argument("--dt", type=float, default=T.DT, help="coupling step (default 0.1)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--thesis", action="store_true", help="the thesis figure instead")
    args = ap.parse_args(argv)
    N, H = args.n, 2.0 / args.n
    T.DT = args.dt
    if args.thesis:
        return thesis_figure()
    figstyle.apply(size=11)
    S, D = grid_masks()
    cmap, norm = figstyle.klempt_cmap()
    fig, ax = plt.subplots(3, len(LAMBDAS), figsize=(13.5, 10))
    ext = (-1, 1, -1, 1)
    kx = N // 2 - 1                                  # section just below x = 0 (through the seed)
    for col, lam in enumerate(LAMBDAS):
        c = nutrient(S, D, lam)
        cs = c[S]
        print(f"Lambda = {lam}: c in the seed {cs.min():.2f} .. {cs.max():.2f}")
        for row, case in ((1, "2sp_case3"), (2, "2sp_case6")):
            sh = np.full((N,) * 3, np.nan)
            sh[S] = shares(case, cs)
            print(f"   {case}: share phi1/(phi1+phi2) {np.nanmin(sh):.3f} .. {np.nanmax(sh):.3f}")
            im = ax[row, col].imshow(sh[kx], origin="lower", extent=ext, cmap=cmap, norm=norm)
            ax[row, col].set_facecolor("#d9dde2")
        imc = ax[0, col].imshow(c[kx], origin="lower", extent=ext, cmap=cmap, norm=norm)
        ax[0, col].set_title(rf"$\Lambda={lam:g}$")
        for r in range(3):
            a = ax[r, col]; a.grid(False)
            figstyle.element_grid(a, np.linspace(-1, 1, 9), np.linspace(-1, 1, 9), alpha=0.25)
            a.set_xticks([-1, 0, 1]); a.set_yticks([-1, 0, 1])
            if col: a.set_yticklabels([])
            if r < 2: a.set_xticklabels([])
    ax[0, 0].set_ylabel("nutrient $c$\n$y$ [mm]")
    ax[1, 0].set_ylabel("case 3 (coexistence)\n$y$ [mm]")
    ax[2, 0].set_ylabel("case 6 (one species wins)\n$y$ [mm]")
    for a in ax[2]:
        a.set_xlabel(r"$z$ [mm]")
    fig.colorbar(imc, ax=ax[0, :], label="$c$ [-]", shrink=0.9, pad=0.015)
    fig.colorbar(im, ax=ax[1:, :], label=r"$\phi_1/(\phi_1+\phi_2)$", shrink=0.9, pad=0.015)
    fig.suptitle(f"Local nutrient in the point model (Python, ANSYS model 2 mm on a {N}$^3$ grid, coupling step {T.DT:g}, "
                 f"section $x={-H / 2:g}$ mm);\n"
                 r"nutrient held at $y=-1$ mm, consumed in the seed; $\Lambda^2=gL^2/d$; grey: no biofilm",
                 fontsize=11.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


THESIS_OUT = HERE.parent / "assets" / "fig_composition_local_nutrient_thesis.png"


def thesis_figure(lams=(2.0, 3.0, 4.0)):
    """Thesis version: nutrient c and the case-6 share on the section through the
    seed, for three Thiele numbers. No title (the caption carries it)."""
    figstyle.apply(size=11)
    S, D = grid_masks()
    cmap, norm = figstyle.klempt_cmap()
    kx = N // 2 - 1
    ext = (-1, 1, -1, 1)
    fig, ax = plt.subplots(2, len(lams), figsize=(10, 6.6))
    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.09, top=0.94, wspace=0.08, hspace=0.12)
    g = np.linspace(-1 + H / 2, 1 - H / 2, N)
    for col, lam in enumerate(lams):
        c = nutrient(S, D, lam)
        sh = np.full((N,) * 3, np.nan)
        sh[S] = shares("2sp_case6", c[S])
        imc = ax[0, col].imshow(c[kx], origin="lower", extent=ext, cmap=cmap, norm=norm)
        ims = ax[1, col].imshow(sh[kx], origin="lower", extent=ext, cmap=cmap, norm=norm)
        ax[1, col].set_facecolor("#d9dde2")
        ax[0, col].contour(g, g, S[kx].astype(float), [0.5], colors="white", linewidths=1.0)
        ax[0, col].set_title(rf"$\Lambda={lam:g}$")
        for r in range(2):
            a = ax[r, col]; a.grid(False)
            a.set_xticks([-1, 0, 1]); a.set_yticks([-1, 0, 1])
            if col: a.set_yticklabels([])
            if r == 0: a.set_xticklabels([])
        ax[1, col].set_xlabel(r"$z$ [mm]")
    ax[0, 0].set_ylabel("$y$ [mm]"); ax[1, 0].set_ylabel("$y$ [mm]")
    c1 = fig.add_axes([0.88, 0.54, 0.015, 0.40]); fig.colorbar(imc, cax=c1, label="nutrient $c$ [-]")
    c2 = fig.add_axes([0.88, 0.09, 0.015, 0.40]); fig.colorbar(ims, cax=c2, label=r"$\phi_1/(\phi_1+\phi_2)$, case 6")
    fig.savefig(THESIS_OUT, dpi=300)
    print("wrote", THESIS_OUT)


if __name__ == "__main__":
    main()
