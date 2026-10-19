#!/usr/bin/env python3
"""The bonded two-layer cylinder under increasing prescribed growth alpha.

Before the interface was bonded (VGLUE; see README.md, cylinder correction)
this deck stopped converging between alpha = 0.01 and 0.015. Bonded, alpha =
0.01 ... 0.2 all reach the end time with 0 errors. Left: maximum radial
displacement of the outer surface against alpha, with the straight line
through the alpha = 0.01 run. Right: the radial displacement across the arc
divided by alpha -- identical curves would mean a purely linear response.

    IKMHIWI03 (listings on F:):
    python ansys_usermat/apdl/plot_cylinder_alpha_sweep.py \
        0.01=F:/biofilm_upf_kusepy/gcr_g0.txt 0.02=.../asw_0p02.txt ...
        -> results/cylinder_alpha_sweep.json (reduced data, pushed) and the figure
    anywhere else, from the reduced data:
    python ansys_usermat/apdl/plot_cylinder_alpha_sweep.py -> assets/growth_cylinder_alpha_sweep.png
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


DATA = Path(__file__).resolve().parent / "results" / "cylinder_alpha_sweep.json"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", help="alpha=path to the result listing (IKMHIWI03); "
                    "without runs the reduced data of --data are plotted")
    ap.add_argument("--data", type=Path, default=DATA,
                    help="reduced data (alpha, u_r across the arc): written when runs are given, read otherwise")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/growth_cylinder_alpha_sweep.png"))
    a = ap.parse_args(argv)

    import json
    if a.runs:
        runs = []
        for arg in a.runs:
            al, path = arg.split("=", 1)
            r = reduce(path)
            runs.append((float(al), np.asarray(r["u_r"], float)))
            print(f"alpha={al}: max u_r {r['u_r'].max():.5e}")
        a.data.parent.mkdir(parents=True, exist_ok=True)
        a.data.write_text(json.dumps({"theta_deg": [float(t) for t in THETA],
                                      "runs": [{"alpha_minus_1": al, "u_r": [float(v) for v in u]}
                                               for al, u in sorted(runs)]}, indent=1))
        print(f"wrote {a.data}")
    else:
        d = json.loads(a.data.read_text())
        runs = [(r["alpha_minus_1"], np.asarray(r["u_r"], float)) for r in d["runs"]]
    runs.sort()
    alphas = np.array([al for al, _ in runs])
    umax = np.array([u.max() for _, u in runs])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import figstyle
    figstyle.apply(8.5)          # drawn at the thesis text width (6.3 in): text prints at 7.5-9 pt
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.3, 2.6),
                                   gridspec_kw={"width_ratios": [1, 1.25]})

    # ---- (a) max u_r vs alpha, against the line through the first run ----
    xs = np.linspace(0, alphas.max() * 1.05, 50)
    axL.plot(xs, umax[0] / alphas[0] * xs, ls=(0, (4, 3)), lw=0.9, color=INK2, zorder=1,
             label="linear from $\\alpha-1=0.01$")
    axL.plot(alphas, umax, "o-", color=BLUE, lw=1.3, ms=4, zorder=3, markeredgecolor="white",
             markeredgewidth=0.8, label="ANSYS")
    axL.set_xlabel("$\\alpha-1$")
    axL.set_ylabel("max $u_r$")
    axL.set_xlim(0, alphas.max() * 1.08)
    axL.set_ylim(0, umax.max() * 1.15)
    axL.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axL.set_title("(a) largest radial displacement", fontsize=9)
    axL.legend(fontsize=7.5, frameon=False, loc="upper left")

    # ---- (b) u_r / alpha across the arc ----
    for (al, u), col in zip(runs, SEQ):
        axR.plot(THETA, u / al, lw=1.2, color=col, label=f"{al:g}")
    axR.set_xlabel("$\\theta$ across the arc [deg]")
    axR.set_ylabel("$u_r / (\\alpha-1)$")
    axR.set_title("(b) displacement across the arc", fontsize=9)
    axR.legend(frameon=False, fontsize=7.5, loc="lower center", ncol=5, title="$\\alpha-1$",
               title_fontsize=7.5, handlelength=1.2, columnspacing=0.8)
    lo, hi = min((u / al).min() for al, u in runs), max((u / al).max() for al, u in runs)
    axR.set_ylim(lo - 0.45 * (hi - lo), hi + 0.1 * (hi - lo))

    for ax in (axL, axR):
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    fig.tight_layout(w_pad=1.5)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=250)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
