#!/usr/bin/env python3
"""The bonded two-layer cylinder under increasing prescribed growth alpha.

Before the interface was bonded (VGLUE; see README.md, cylinder correction)
this deck stopped converging between alpha = 0.01 and 0.015. Bonded, alpha =
0.01 ... 0.2 all reach the end time with 0 errors. Left: maximum radial
displacement of the outer surface against alpha, with the straight line
through the alpha = 0.01 run. Right: the radial displacement across the arc
divided by alpha -- identical curves would mean a purely linear response.

    python ansys_usermat/apdl/plot_cylinder_alpha_sweep.py \
        0.01=F:/biofilm_upf_kusepy/gcr_g0.txt 0.02=.../asw_0p02.txt ... \
        -o assets/growth_cylinder_alpha_sweep.png
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cylinder_mesh_convergence import THETA, reduce  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
BLUE = "#2a78d6"
# sequential blue, light -> dark (steps 250 .. 650), one per alpha
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="alpha=path to the result listing")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/growth_cylinder_alpha_sweep.png"))
    a = ap.parse_args(argv)

    runs = []
    for arg in a.runs:
        al, path = arg.split("=", 1)
        r = reduce(path)
        runs.append((float(al), r["u_r"]))
        print(f"alpha={al}: max u_r {r['u_r'].max():.5e}")
    runs.sort()
    alphas = np.array([al for al, _ in runs])
    umax = np.array([u.max() for _, u in runs])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import figstyle
    figstyle.apply(11)
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.0, 3.9),
                                   gridspec_kw={"width_ratios": [1, 1.25]})

    # ---- left: max u_r vs alpha, against the line through the first run ----
    xs = np.linspace(0, alphas.max() * 1.05, 50)
    axL.plot(xs, umax[0] / alphas[0] * xs, ls=(0, (4, 3)), lw=1.2, color=INK2,
             zorder=1)
    xl = alphas.max() * 0.72
    axL.text(xl + 0.004, umax[0] / alphas[0] * xl, "linear from $\\alpha-1$ = 0.01",
             ha="left", va="top", fontsize=8, color=INK2)
    axL.axvspan(0.01, 0.015, color=GRID, zorder=0)
    axL.text(0.0135, umax.max() * 0.55, "old limit\n(floating\nlayer)",
             ha="left", va="center", fontsize=7.5, color=INK2)
    axL.plot(alphas, umax, "o-", color=BLUE, lw=2, ms=6, zorder=3,
             markeredgecolor="white", markeredgewidth=1.5)
    for al, u in zip(alphas, umax):
        if al >= 0.05:
            axL.annotate(f"{u:.2e}", (al, u), xytext=(-6, 6),
                         textcoords="offset points", ha="right", fontsize=7.5,
                         color=INK2)
    axL.set_xlabel("prescribed growth $\\alpha-1$")
    axL.set_ylabel("max radial displacement $u_r$")
    axL.set_xlim(0, alphas.max() * 1.08)
    axL.set_ylim(0, umax.max() * 1.15)
    axL.set_title("all runs reach the end time, $\\alpha-1\\leq0.2$",
                  fontsize=10, color=INK)

    # ---- right: u_r / alpha across the arc ----
    for (al, u), col in zip(runs, SEQ):
        axR.plot(THETA, u / al, lw=2, color=col, label=f"$\\alpha-1$ = {al:g}")
    axR.set_xlabel("$\\theta$ across the arc  [deg]")
    axR.set_ylabel("$u_r / (\\alpha-1)$")
    axR.set_title("radial displacement across the arc",
                  fontsize=10, color=INK)
    axR.legend(frameon=False, fontsize=8, loc="lower center", ncol=5,
               handlelength=1.6, columnspacing=1.0)
    lo, hi = min((u / al).min() for al, u in runs), max((u / al).max() for al, u in runs)
    axR.set_ylim(lo - 0.35 * (hi - lo), hi + 0.1 * (hi - lo))

    for ax in (axL, axR):
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=200, facecolor="#fcfcfb")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
