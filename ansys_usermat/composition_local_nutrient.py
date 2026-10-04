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

    python ansys_usermat/composition_local_nutrient.py
        -> assets/fig_composition_local_nutrient.png

Result (4 Oct 2026), share phi_1/(phi_1+phi_2) in the 32 seed elements:
  Lambda = 1 (c in the seed 0.86-0.91): case 3 0.608, case 6 0.014-0.015
  Lambda = 2 (0.60-0.73):               case 3 0.604-0.606, case 6 0.015-0.017
  Lambda = 3 (0.38-0.57):               case 3 0.600-0.604, case 6 0.017-0.231
  Lambda = 4 (0.23-0.44):               case 3 0.596-0.601, case 6 0.118-0.437
Case 3 (coexistence) hardly depends on the nutrient level. In case 6 the
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
N, H = 8, 0.25
LAMBDAS = [1.0, 2.0, 3.0, 4.0]


def grid_masks():
    _, cen, seed, nut, _ = M.read_deck(M.DECK)
    idx = lambda p: tuple(int(round((v + 1) / H - 0.5)) for v in p)
    S = np.zeros((N,) * 3, bool); D = np.zeros((N,) * 3, bool)
    for e in seed:
        S[idx(cen[e])] = True
    for e in nut:
        D[idx(cen[e])] = True
    return S, D


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
    T.ms.set_case(case)
    th, hp = T.ms.ECOLOGY_CASE["theta"], dict(T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    out = []
    for c in cvals:
        h = dict(hp); h["c"] = hp["c"] * c
        adv = T.stepper(th, h)
        G = T.seed_state(1); p = np.array([T.PHI_CAP])
        for _ in range(10):
            G = T.rescale(G, p); G = T.rescale(adv(G), p)
        out.append(G[0, 0] / (G[0, 0] + G[0, 1]))
    return np.array(out)


def main():
    figstyle.apply(size=11)
    S, D = grid_masks()
    cmap, norm = figstyle.klempt_cmap()
    fig, ax = plt.subplots(3, len(LAMBDAS), figsize=(13.5, 10))
    ext = (-1, 1, -1, 1)
    kx = 3                                            # section x = -0.125 mm (through the seed)
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
            figstyle.element_grid(a, np.linspace(-1, 1, N + 1), np.linspace(-1, 1, N + 1), alpha=0.25)
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
    fig.suptitle("Local nutrient in the point model (Python, ANSYS model 2 mm, section $x=-0.125$ mm);\n"
                 r"nutrient held at $y=-1$ mm, consumed in the seed; $\Lambda^2=gL^2/d$; grey: no biofilm",
                 fontsize=11.5)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
