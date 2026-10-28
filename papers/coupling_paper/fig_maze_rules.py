#!/usr/bin/env python3
"""fig_maze_rules.py -- Figure 7 of the coupling paper: the maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11; read off the
figure, 4 um cells) with the agreed coupling and the antibiotic held at the entrance from T* = 20 (alpha*_0 = 100),
at T* = 60, for the two rules for newly reached cells (A-73, A-74): (a) new cells take the phi-weighted mean state of
their colonised neighbours, (b) new cells take the colony's start state. Colour: living volume fraction sum_i phibar_i
where phi_3D >= 0.05; walls grey, nutrient source (c = 1) outlined in blue, antibiotic source (entrance) in red.
(c) living volume fraction averaged over the colonised cells against T* for both rules; the antibiotic starts at the
dotted line. The Abaqus element with the start-state rule gives the same (A-73).

Data: ansys_usermat/apdl/maze_coupled_2d.py --h 4 --tend 60 --ab --maze-file abaqus_living/maze_k24.txt [--fresh-start]
-> figdata/maze_h4_tend60_ab_k24fig{,_fresh}.npz (fields, gitignored; rerun to regenerate) and
results/2026-10-maze/maze_h4_tend60_ab_k24fig{,_fresh}.json (key hist.living_mean). The fields are from the reruns of
27 Oct 2026 (after the inherit fix of A-73); their histories agree with the committed JSONs of 24 Oct to 5e-5.

    python papers/coupling_paper/fig_maze_rules.py -> papers/coupling_paper/fig_maze_rules.pdf
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

T_SHOW = 60.0
RULES = (("inherit", "maze_h4_tend60_ab_k24fig", "neighbours' state"),
         ("fresh", "maze_h4_tend60_ab_k24fig_fresh", "start state"))


def main():
    figstyle.apply(size=8.5)
    fig = plt.figure(figsize=(6.3, 2.45), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.045, 1.05], wspace=0.04)
    axs = [fig.add_subplot(gs[0, i]) for i in (0, 1)]
    cax = fig.add_subplot(gs[0, 2])
    axh = fig.add_subplot(gs[0, 3])
    m = None
    for ax, (rule, tag, name), letter in zip(axs, RULES, "ab"):
        z = np.load(ROOT / "figdata" / f"{tag}.npz")
        X, Y = z["X"], z["Y"]
        OPEN, NUT, ABS = (z[k].astype(bool) for k in ("open", "nut", "ab_src"))
        ext = [X.min(), X.max(), Y.min(), Y.max()]
        phi, L = z[f"phi_{T_SHOW}"], z[f"L_{T_SHOW}"]
        img = np.where(OPEN & (phi >= 0.05), L, np.nan)
        m = ax.imshow(img.T, origin="lower", extent=ext, cmap="Greens", vmin=0, vmax=1, interpolation="nearest")
        ax.imshow(np.where(OPEN, np.nan, 1.0).T, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=4,
                  interpolation="nearest")
        x1, y1 = X[:, 0], Y[0, :]
        ax.contour(x1, y1, NUT.T.astype(float), [0.5], colors="tab:blue", linewidths=0.6)
        ax.contour(x1, y1, ABS.T.astype(float), [0.5], colors="tab:red", linewidths=0.6)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_aspect("equal")
        ax.text(-0.02, 1.02, f"({letter})", transform=ax.transAxes, va="bottom", ha="left")
        col = ax.get_images()[0]  # the living fraction
        print(f"{rule:8s} T* = {T_SHOW:g}: colonised (phi >= 0.05) {int((OPEN & (phi >= 0.05)).sum())} cells, "
              f"living fraction mean {np.nanmean(img):.3f}, min {np.nanmin(img):.3f}, max {np.nanmax(img):.3f}")
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r"living volume fraction $\sum_i\bar\phi_i$")
    cb.ax.tick_params(labelsize=7)
    # (c) living fraction of the colonised cells over time
    for (rule, tag, name), col, ls in zip(RULES, ("C2", "C3"), ("-", "--")):
        h = json.load(open(ROOT / "ansys_usermat/apdl/results/2026-10-maze" / f"{tag}.json"))["hist"]
        t, liv = np.asarray(h["t"]), np.asarray(h["living_mean"])
        n = min(len(t), len(liv))
        t, liv = t[-n:] if len(t) > n else t, liv[-n:] if len(liv) > n else liv
        axh.plot(t, liv, color=col, ls=ls, lw=1.1, label=name)
        i_ab = np.searchsorted(t, 20.0)
        print(f"{rule:8s} living_mean: T* = 20 {liv[max(i_ab - 1, 0)]:.3f}, min after {liv[i_ab:].min():.3f} "
              f"(T* = {t[i_ab:][liv[i_ab:].argmin()]:.1f}), end {liv[-1]:.3f}")
    axh.axvline(20.0, color="0.4", ls=":", lw=0.8)
    axh.set_xlabel("$T^*$")
    axh.set_ylabel(r"mean $\sum_i\bar\phi_i$ of the colonised cells")
    axh.set_xlim(0, 60)
    axh.set_ylim(0, None)
    axh.legend(frameon=False, fontsize=7, loc="upper right")
    axh.text(-0.3, 1.02, "(c)", transform=axh.transAxes, va="bottom", ha="left")
    out = HERE / "fig_maze_rules.pdf"
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
