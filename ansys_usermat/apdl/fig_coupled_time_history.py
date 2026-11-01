#!/usr/bin/env python3
"""Seed means over time for the two coupled examples (separate and mixed start, equal amounts), with and without
antibiotic: (a) living volume fraction sum_i phibar_i, (b) share phi_1/(phi_1+phi_2), (c) nutrient c. The runs of
thesis Table tab:ch5-two-species (Python, 64^3, outer step 0.025, results/2026-10-3d). Printed width 16 cm, Times.

    python ansys_usermat/apdl/fig_coupled_time_history.py   -> assets/fig_coupled_time_history.png / .pdf
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

RES = HERE / "results" / "2026-10-3d"
P = "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_{}a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_{}.json"
CASES = [("separate", "", "separate", "#3b5bdb", "-"), ("separate, antibiotic", "", "separate_ab", "#3b5bdb", "--"),
         ("mixed", "mixeq_", "mixed", "#7048a8", "-"), ("mixed, antibiotic", "mixeq_", "mixed_ab", "#7048a8", "--")]
OUT = ROOT / "assets" / "fig_coupled_time_history.png"


def main():
    figstyle.apply()
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), constrained_layout=True)
    keys = [("living_seed", r"$\Sigma_i\,\bar\phi_i$", (0, 1.02)), ("share_a_seed", r"$\phi_1/(\phi_1+\phi_2)$", (0.4, 1.02)),
            ("c_seed", r"$c$", (0.4, 1.02))]
    for lab, mix, ex, col, ls in CASES:
        r = json.loads((RES / P.format(mix, ex)).read_text())["runs"][ex]
        for ax, (k, _, _) in zip(axs, keys):
            ax.plot(r["t"], r[k], color=col, ls=ls, lw=1.3, label=lab)
    for i, (ax, (_, yl, lim)) in enumerate(zip(axs, keys)):
        ax.axvline(1 / 3, color="0.55", lw=0.7, ls=":")
        ax.set_xlim(0, 1); ax.set_ylim(*lim); ax.set_xlabel(r"$T^*$"); ax.set_ylabel(yl)
        ax.text(0.03, 0.04 if i else 0.93, f"({'abc'[i]})", transform=ax.transAxes, va="bottom" if i else "top")
    axs[1].text(1 / 3 + 0.02, 0.42, "antibiotic\nfrom $T^*=1/3$", fontsize=7, color="0.35")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False)
    fig.savefig(OUT, dpi=400); fig.savefig(OUT.with_suffix(".pdf"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
