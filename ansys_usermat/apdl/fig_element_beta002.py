#!/usr/bin/env python3
"""The seeded element over time at beta = 0.02 mm^2/T* (16^3, dt = 0.025),
Eq. 36 in the material routine against the partner's growth variable; the
beta = 0.02 version of figs_1005.fig_elem, drawn from the run JSON (no F:).

    python ansys_usermat/apdl/fig_element_beta002.py -> assets/fig_element_beta002.png

Runs: results/2026-10-05_ansys/ds16_beta002_dt4.json (Eq. 36) and
results/2026-10-05_week/w16_1sppartner_b002.json (the partner's variable),
same deck otherwise. The post element is the seeded element at the 8^3
element 220 position (2151 on 16^3). Stresses MPa -> Pa; alpha is the
material routine's variable, alpha - 1 in the thesis notation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402

RES = HERE / "results"
RUNS = (("Eq. 36 in the material routine", RES / "2026-10-05_ansys" / "ds16_beta002_dt4.json", "-"),
        ("partner's growth variable", RES / "2026-10-05_week" / "w16_1sppartner_b002.json", "--"))
OUT = HERE.parents[1] / "assets" / "fig_element_beta002.png"


def main():
    figstyle.apply(size=8.5)     # drawn at the thesis text width (6.3 in); line style = run, see the caption
    fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.5))
    for k, (name, f, ls) in enumerate(RUNS):
        e = {k2: np.asarray(v, float) for k2, v in json.load(open(f))["elem_stress"].items()}
        t = e["time"][1:]
        ax[0].plot(t, e["seqv"][1:] * 1e6, "C0" + ls, lw=1.1, label="seed element" if k == 0 else None)
        ax[0].plot(t, e["seqv_nbr_max"][1:] * 1e6, "C3" + ls, lw=1.1, label="largest neighbour" if k == 0 else None)
        ax[1].plot(t, e["alpha"][1:], "C2" + ls, lw=1.1, label=("Eq. 36 in the material routine" if k == 0
                                                             else "element's own growth variable"))
        print(f"{name}: T* {t[-1]:.2f}, vM {e['seqv'][-1] * 1e6:.3e} Pa, neighbour max "
              f"{e['seqv_nbr_max'][-1] * 1e6:.3e} Pa, alpha-1 {e['alpha'][-1]:.4e}")
    ax[0].set(xlabel=r"$T^*$", ylabel="von Mises [Pa]", title="(a) von Mises stress")
    ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0)); ax[0].legend(fontsize=7.5, frameon=False)
    ax[1].set(xlabel=r"$T^*$", ylabel=r"$\alpha-1$", title=r"(b) growth of the seed element")
    ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0)); ax[1].legend(fontsize=7.5, frameon=False)
    fig.tight_layout(w_pad=1.5); fig.savefig(OUT, dpi=300); fig.savefig(OUT.with_suffix(".pdf")); print("wrote", OUT)


if __name__ == "__main__":
    main()
