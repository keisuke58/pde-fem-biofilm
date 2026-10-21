#!/usr/bin/env python3
"""fig_maze.py -- maze runs of maze_coupled_2d.py, laid out as Table 7 of Klempt et al. 2024 (rows: times).

Columns: (a) phi_3D and (b) nutrient c without antibiotic; (c) living volume fraction sum_i phibar_i and
(d) antibiotic a with the antibiotic held at the entrance from T_END/3. Walls dark grey, cells with phi_3D < 0.1
light grey in (c).

    python ansys_usermat/apdl/fig_maze.py [--tend 8] -> assets/fig_maze.pdf/.png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

T_END = float(sys.argv[sys.argv.index("--tend") + 1]) if "--tend" in sys.argv else 8.0


def main():
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d0 = np.load(ROOT / "figdata" / f"maze_h2.5_tend{T_END:g}.npz")
    d1 = np.load(ROOT / "figdata" / f"maze_h2.5_tend{T_END:g}_ab.npz")
    ts = sorted({float(k.split("_")[1]) for k in d0.files if k.startswith("phi_")})
    pick = [ts[i] for i in (0, 2, 4, 6, 9)]
    op = d0["open"]
    X, Y = d0["X"], d0["Y"]
    h = X[1, 0] - X[0, 0]
    e = np.r_[X[:, 0] - h / 2, X[-1, 0] + h / 2]
    cols = [(d0, "phi", "viridis", r"(a) $\phi_{3D}$"), (d0, "c", "Blues", "(b) nutrient $c$"),
            (d1, "L", "Greens", r"(c) $\Sigma_i\bar\phi_i$, antibiotic"), (d1, "a", "Reds", "(d) antibiotic $a$")]
    fig, axs = plt.subplots(len(pick), 4, figsize=(6.3, 1.62 * len(pick)), constrained_layout=True)
    for r, T in enumerate(pick):
        for c_, (d, fld, cm, title) in enumerate(cols):
            a = axs[r, c_]
            a.set_facecolor("0.25")
            u = d[f"{fld}_{T}"].copy()
            u = np.where(op, u, np.nan)
            if fld == "L":
                u = np.where(d[f"phi_{T}"] >= 0.1, u, np.nan)
                a.pcolormesh(e, e, np.ma.masked_where(~op, np.ones_like(u)).T, cmap="Greys", vmin=0, vmax=5,
                             shading="flat", rasterized=True)
            im = a.pcolormesh(e, e, np.ma.masked_invalid(u).T, cmap=cm, vmin=0, vmax=1, shading="flat",
                              rasterized=True)
            a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
            for s in a.spines.values():
                s.set_visible(False)
            if r == 0:
                a.set_title(title, fontsize=7.5)
            if c_ == 0:
                a.set_ylabel(f"$T^*={T:g}$", fontsize=7.5)
            if r == len(pick) - 1:
                cb = fig.colorbar(im, ax=a, orientation="horizontal", fraction=0.06, pad=0.02)
                cb.ax.tick_params(labelsize=6)
    out = ROOT / "assets" / "fig_maze"
    fig.savefig(out.with_suffix(".pdf"), dpi=300); fig.savefig(out.with_suffix(".png"), dpi=300)
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    main()
