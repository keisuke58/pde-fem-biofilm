"""fig_thesis_cal5.py -- thesis Ch. 4: the five-species point model calibrated in Ch. 3 (CH, MAP of the
measured-data calibration) in the element with the one-way coupling of the chapter (Eq. 36, point model scaled
to phi_3D, phi_cap = 0.9, s = 0.15), 8^3, at a seed Gauss point: (a) shares phi_i / sum_j phi_j,
(b) living volume fractions phi_i psi_i and their sum.

    python ansys_usermat/apdl/fig_thesis_cal5.py   -> assets/fig_thesis_cal5.pdf, .png
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat")); sys.path.insert(0, str(ROOT / "ansys_usermat/apdl"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from fig_cal5_appendix import NAMES, trace  # noqa: E402

if __name__ == "__main__":
    t, phi, psi = trace("w8_CH_g1_s015_raw_pm.json")
    sh = phi / phi.sum(1, keepdims=True); liv = phi * psi
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.1))
    for i in range(5):
        l, = ax[0].plot(t, sh[:, i], lw=1.5, label=NAMES[i])
        ax[1].plot(t, liv[:, i], lw=1.5, color=l.get_color())
    ax[1].plot(t, liv.sum(1), color="#222", lw=1.8, label="$\\sum_i\\phi_i\\psi_i$")
    ax[0].set(xlabel="$T^*$", ylabel="$\\phi_i/\\sum_j\\phi_j$", title="(a) composition", ylim=(0, 0.45), xlim=(0, 1))
    ax[1].set(xlabel="$T^*$", ylabel="$\\phi_i\\psi_i$", title="(b) living volume fractions", xlim=(0, 1))
    ax[0].legend(ncol=3, fontsize=8, frameon=False, loc="upper left")
    ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    out = ROOT / "assets" / "fig_thesis_cal5"
    fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
    print("final shares", sh[-1].round(3), "sum phi", phi.sum(1)[[0, -1]].round(3), "living", liv.sum(1)[[0, -1]].round(3))
