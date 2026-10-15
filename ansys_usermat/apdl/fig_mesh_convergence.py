"""fig_mesh_convergence.py -- mesh convergence of the seed means (CH, five species, T* = 1) for Eq. 36 and for
growth from the living biomass, 8^3, 16^3, 24^3: von Mises / E, alpha - 1 and the difference between the two.
Earlier calibration: all three meshes; calibration to the measured data: the meshes run so far.

    python ansys_usermat/apdl/fig_mesh_convergence.py   -> assets/fig_mesh_convergence.pdf, .png
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
E = 10.0   # Pa
SETS = (("earlier calibration", "w{n}_CH_g1_s015.json", "w{n}_CH_cal_k0882.json", "-"),
        ("measured data", "w{n}_CH_g1_s015_raw.json", "w{n}_CH_cal_k0814_raw.json", "--"))


def seed(fn):
    p = R / fn
    if not p.exists():
        return None
    r = json.load(open(p)); a = r["all_stress"]
    s = np.isin(np.array(a["elem"]), r["seed_BIOFILM1"])
    return np.mean(np.array(a["alpha"])[s]), np.mean(np.array(a["seqv"])[s]) * 1e6 / E


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for lab, fa, fc, ls in SETS:
        for k, (f, name, col) in enumerate(((fa, "Eq. 36", "#4c72b0"), (fc, "living biomass", "#dd8452"))):
            pts = [(n, seed(f.format(n=n))) for n in (8, 16, 24)]
            pts = [(n, v) for n, v in pts if v is not None]
            n = [p[0] for p in pts]
            ax[0].plot(n, [p[1][1] * 1e5 for p in pts], ls, marker="o", color=col, label=f"{name}, {lab}")
            ax[1].plot(n, [p[1][0] * 1e3 for p in pts], ls, marker="o", color=col)
        both = [(n, seed(fa.format(n=n)), seed(fc.format(n=n))) for n in (8, 16, 24)]
        both = [(n, a, c) for n, a, c in both if a is not None and c is not None]
        ax[2].plot([b[0] for b in both], [100 * (b[2][1] / b[1][1] - 1) for b in both], ls, marker="s", color="#333",
                   label=lab)
        for b in both:
            print(lab, b[0], f"vM/E {b[1][1]:.4e} {b[2][1]:.4e} alpha-1 {b[1][0]:.4e} {b[2][0]:.4e}")
    for a in ax:
        a.set_xticks([8, 16, 24], ["$8^3$", "$16^3$", "$24^3$"]); a.set_xlabel("mesh")
    ax[0].set(ylabel="seed von Mises$/E$ ($10^{-5}$)", title="(a) seed stress")
    ax[1].set(ylabel="seed $\\alpha-1$ ($10^{-3}$)", title="(b) seed growth", ylim=(0.6, 0.7))
    ax[2].set(ylabel="living biomass vs. Eq. 36 (%)", title="(c) difference in stress", ylim=(-10, 0))
    ax[0].legend(fontsize=7, frameon=False, loc="lower right")
    ax[2].legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    out = ROOT / "assets" / "fig_mesh_convergence"
    fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"), dpi=170)
