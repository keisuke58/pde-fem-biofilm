"""Figures for the one-page summary of the meeting of 19 Oct 2026
(meeting_1019_summary.tex): a sketch of the two coupling schemes and of growth
from the living biomass, and the calibrated five-species run in the element
(CH, 8^3) against Eq. 36.

    python ansys_usermat/apdl/meeting1019_figs.py
writes assets/fig_meeting1019_schemes.pdf and assets/fig_meeting1019_cal5.pdf
"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle
figstyle.apply()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = ROOT / "assets"


def box(ax, x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=fc, ec="0.3", lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.5)


def arrow(ax, a, b, text="", dy=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=10, lw=0.9, color="0.2"))
    if text:
        ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + dy, text, ha="center", va="bottom", fontsize=8.5)


def schemes():
    fig, axs = plt.subplots(1, 3, figsize=(10, 2.9))
    titles = ["(a) scaled: amount from the field",
              "(b) unscaled: $\\phi_0$ stays auxiliary",
              "(c) growth from the living biomass"]
    for ax, t in zip(axs, titles):
        ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1.08); ax.axis("off"); ax.set_title(t, fontsize=9)
    fld, pm, mech = "#dbe9f6", "#fde3c8", "#e2efd9"
    # (a)
    ax = axs[0]
    box(ax, 0.05, 0.68, 0.38, 0.22, "field\n$\\phi_{3D}$ (NEM)", fld)
    box(ax, 0.57, 0.68, 0.38, 0.22, "point model\n$\\phi_i,\\ \\psi_i$", pm)
    box(ax, 0.05, 0.12, 0.38, 0.22, "growth\n$\\dot\\alpha=k_\\alpha\\phi_{3D}$", mech)
    box(ax, 0.57, 0.12, 0.38, 0.22, "composition\n$\\phi_i/\\sum_j\\phi_j$", pm)
    arrow(ax, (0.43, 0.79), (0.57, 0.79), "scale $\\sum\\phi_i=\\phi_{3D}$", 0.13)
    arrow(ax, (0.24, 0.68), (0.24, 0.34))
    arrow(ax, (0.76, 0.68), (0.76, 0.34))
    ax.text(0.5, 0.0, "$\\phi_0=1-\\phi_{3D}$: physical void", ha="center", fontsize=8)
    # (b)
    ax = axs[1]
    box(ax, 0.05, 0.68, 0.38, 0.22, "field\n$\\phi_{3D}$ (NEM)", fld)
    box(ax, 0.57, 0.68, 0.38, 0.22, "point model\nown $\\phi_0$", pm)
    box(ax, 0.05, 0.12, 0.38, 0.22, "growth\n$\\dot\\alpha=k_\\alpha\\phi_{3D}$", mech)
    box(ax, 0.57, 0.12, 0.38, 0.22, "composition", pm)
    arrow(ax, (0.43, 0.79), (0.57, 0.79), "start at $\\phi_{\\min}$", 0.13)
    arrow(ax, (0.24, 0.68), (0.24, 0.34))
    arrow(ax, (0.76, 0.68), (0.76, 0.34))
    ax.text(0.5, 0.0, "composition depends on the local age", ha="center", fontsize=8)
    # (c)
    ax = axs[2]
    box(ax, 0.05, 0.68, 0.38, 0.22, "field\n$\\phi_{3D}$ (NEM)", fld)
    box(ax, 0.57, 0.68, 0.38, 0.22, "point model\n$\\psi_i$ (lysis)", pm)
    box(ax, 0.30, 0.12, 0.40, 0.22, "growth\n$\\dot\\alpha=\\kappa\\sum_i\\phi_i\\psi_i$", mech)
    arrow(ax, (0.43, 0.79), (0.57, 0.79), "scale", 0.13)
    arrow(ax, (0.76, 0.68), (0.60, 0.34), "")
    ax.text(0.5, 0.0, "lower $\\psi_i$ lowers growth from inside", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_meeting1019_schemes.pdf")


def cal5():
    R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
    pm = json.load(open(R / "w8_CH_cal_k0882_pm.json"))
    cols = pm["columns"]; C = {c: i for i, c in enumerate(cols)}
    rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]]]
    t = np.arange(1, len(rows) + 1) / len(rows)
    phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
    psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
    names = ["So", "An", "Vei", "Fn", "Pg"]

    def seed(fn):
        r = json.load(open(R / fn)); a = r["all_stress"]
        s = np.isin(np.array(a["elem"]), r["seed_BIOFILM1"])
        return np.array(a["alpha"])[s], np.array(a["seqv"])[s] * 1e6   # MPa -> Pa

    ab, vb = seed("w8_CH_g1_s015.json"); am, vm = seed("w8_CH_cal_k0882.json")
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
    for i in range(5):
        ax[0].plot(t, phi[:, i], label=names[i])
        ax[1].plot(t, phi[:, i] * psi[:, i])
    ax[0].set(xlabel="$T^*$", ylabel="$\\phi_i$", title="composition, seed Gauss point")
    ax[0].legend(fontsize=7.5, frameon=False, ncol=2)
    ax[1].set(xlabel="$T^*$", ylabel="$\\phi_i\\psi_i$", title="living fractions")
    ax[2].scatter(ab * 1e3, vb * 1e4, s=10, label="Eq. 36 ($\\phi$ of the field)")
    ax[2].scatter(am * 1e3, vm * 1e4, s=12, marker="x", label="living biomass")
    ax[2].set(xlabel="$\\alpha-1$ ($10^{-3}$), seed elements", ylabel="von Mises ($10^{-4}$ Pa)",
              title="seed growth and stress")
    ax[2].legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_meeting1019_cal5.pdf")


if __name__ == "__main__":
    schemes()
    cal5()
    print("wrote", OUT / "fig_meeting1019_schemes.pdf", OUT / "fig_meeting1019_cal5.pdf")
