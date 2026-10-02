#!/usr/bin/env python3
"""Slide figure: the point model used at every Gauss point against Klempt et
al. 2026, the two cases the coupling runs (Table 1, cases 3 and 6).

Lines: this repo's integrator (klempt2026_reproduction.simulate, variant
"repo"), living fractions phi_i psi_i and phi0. Markers: values read off the
paper's figures by eye (klempt2026_cases.py): the end values of phi_i psi_i
and the step at which phi0 reaches about 0.

    python JAXFEM/klempt2026_slide_fig.py   -> assets/fig_point_model_klempt2026.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from klempt2026_cases import CASES  # noqa: E402
from klempt2026_reproduction import simulate  # noqa: E402

OUT = HERE.parent / "assets" / "fig_point_model_klempt2026.png"
SPECIES = ["#2a78d6", "#eb6834"]
INK, MUTED = "#1f2933", "#7b8794"
TITLES = {"2sp_case3": "case 3 (Fig. 3): coexistence",
          "2sp_case6": "case 6 (Fig. 6): species 2 takes over"}


def main():
    sys.path.insert(0, str(HERE.parent / "ansys_usermat"))
    import figstyle
    figstyle.apply()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, name in zip(axs, TITLES):
        case = CASES[name]
        r = simulate(case)
        k = np.arange(r["phi"].shape[0])
        pb = r["phi"] * r["psi"]
        for i in range(2):
            ax.plot(k, pb[:, i], color=SPECIES[i], lw=2.2,
                    label=rf"species {i + 1}, $\phi_{i + 1}\psi_{i + 1}$")
            ax.plot(k[-1], case["paper"]["phibar_end"][i], "D", ms=9, mfc="none",
                    mec=SPECIES[i], mew=2)
        ax.plot(k, r["phi0"], color=INK, lw=1.6, ls="--", label=r"$\phi_0$ (auxiliary)")
        z = case["paper"]["phi0_zero_step"]
        ax.plot(z, 0.0, "D", ms=9, mfc="none", mec=INK, mew=2)
        ax.set(xlabel=r"time step ($\Delta t = 10^{-4}$)", ylim=(-0.04, 1.04), title=TITLES[name])
        ax.set_xlim(-0.02 * k[-1], 1.06 * k[-1])
    axs[0].set_ylabel("volume fraction")
    axs[0].plot([], [], "D", mfc="none", mec=MUTED, mew=2, label="read off the paper's figure")
    axs[0].legend(fontsize=9, loc="upper right")
    fig.suptitle("Point model of Klempt et al. 2026 as implemented here (lines) "
                 "against the published figures (markers)", fontsize=12)
    fig.tight_layout()
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
