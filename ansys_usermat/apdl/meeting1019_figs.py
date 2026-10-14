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


def last_per_step(rows, C):
    """One row per substep (the last, i.e. the converged one) and its time T* = substep/n."""
    d = {}
    for r in rows:
        d[int(r[C["isubst"]])] = r
    k = sorted(d)
    return [d[i] for i in k], np.array(k) / max(k)


def cal5():
    R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
    pm = json.load(open(R / "w8_CH_cal_k0882_pm.json"))
    cols = pm["columns"]; C = {c: i for i, c in enumerate(cols)}
    rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
    rows, t = last_per_step(rows, C)
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


def lysis():
    """Antibiotic in the point model, CH 8^3, T* = 1: seed means from the IKMHIWI03 report
    (IKMHIWI03_TO_CLOUD.ja.md, 14 Oct, table A-1/A-2; scheme_c_table.py)."""
    lab = ["(a)", "(a)\n+ antibiotic", "(c)", "(c)\n+ antibiotic"]
    a1 = np.array([6.7517e-4, 6.7517e-4, 6.7518e-4, 5.0066e-4])
    vm = np.array([4.6838e-4, 4.6838e-4, 4.3192e-4, 2.9397e-4])
    sp = np.array([2.616, 1.743, 2.132, 0.784])
    col = ["#7a8fa6", "#7a8fa6", "#c0603a", "#c0603a"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.2))
    for k, (y, yl, t) in enumerate([(a1 * 1e4, "$\\alpha-1$ ($10^{-4}$)", "seed growth"),
                                     (vm * 1e4, "von Mises ($10^{-4}$ Pa)", "seed stress"),
                                     (sp, "$\\sum_i\\psi_i$", "living fractions")]):
        b = ax[k].bar(range(4), y, color=col, width=0.65)
        for i in (1, 3):
            b[i].set_hatch("//"); b[i].set_edgecolor("white")
        ax[k].set_xticks(range(4), lab, fontsize=7.5)
        ax[k].set(ylabel=yl, title=t)
    for i, d in ((3, "$-26\\,\\%$"), ):
        ax[0].annotate(d, (i, a1[i] * 1e4), ha="center", va="bottom", fontsize=8)
    ax[1].annotate("$-32\\,\\%$", (3, vm[3] * 1e4), ha="center", va="bottom", fontsize=8)
    ax[1].annotate("same", (1, vm[1] * 1e4), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_meeting1019_lysis.pdf")


def lysis_trace():
    """Point model at the first seed Gauss point, scheme (c), without and with the antibiotic."""
    R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
    names = ["So", "An", "Vei", "Fn", "Pg"]

    def tr(fn):
        pm = json.load(open(R / fn)); C = {c: i for i, c in enumerate(pm["columns"])}
        rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
        rows, t = last_per_step(rows, C)
        phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
        psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
        return t, phi, psi

    t0, f0, p0 = tr("w8_CH_cal_k0882_pm.json"); t1, f1, p1 = tr("w8_CH_cal_k0882_ab_pm.json")
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.7))
    for i in range(5):
        l, = ax[0].plot(t0, f0[:, i] / f0.sum(1), label=names[i])
        ax[0].plot(t1, f1[:, i] / f1.sum(1), "--", color=l.get_color())
        ax[1].plot(t0, p0[:, i], color=l.get_color())
        ax[1].plot(t1, p1[:, i], "--", color=l.get_color())
    ax[2].plot(t0, (f0 * p0).sum(1), "k-", label="without antibiotic")
    ax[2].plot(t1, (f1 * p1).sum(1), "k--", label="with antibiotic")
    ax[0].set(xlabel="$T^*$", ylabel="$\\phi_i/\\sum_j\\phi_j$", title="composition")
    ax[0].legend(fontsize=7.5, frameon=False, ncol=2)
    ax[1].set(xlabel="$T^*$", ylabel="$\\psi_i$", title="living fraction of each species")
    ax[2].set(xlabel="$T^*$", ylabel="$\\sum_i\\phi_i\\psi_i$", title="living biomass (drives growth)")
    ax[2].legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_meeting1019_lysis_trace.pdf")


def switch_on():
    """Antibiotic switched on at T* = 0.5 (as Fritsch et al. 2025, Fig. 18), scheme (c), CH 8^3,
    first seed Gauss point; seed stress from the IKMHIWI03 report of 14 Oct (A-6 table)."""
    R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
    runs = [("w8_CH_cal_k0882_pm.json", "none", "k", "-"),
            ("w8_CH_cal_k0882_ab05_pm.json", "all $b_i=1$", "#4c72b0", "--"),
            ("w8_CH_cal_k0882_ab05_a065_pm.json", "all $b_i=1$, stronger", "#4c72b0", ":"),
            ("w8_CH_cal_k0882_ab05_So_pm.json", "So sensitive", "#dd8452", "--"),
            ("w8_CH_cal_k0882_ab05_ana_pm.json", "anaerobes sensitive", "#c44e52", "--")]
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.7), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for fn, lab, c, ls in runs:
        pm = json.load(open(R / fn)); C = {k: i for i, k in enumerate(pm["columns"])}
        rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
        rows, t = last_per_step(rows, C)
        phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
        psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
        ax[0].plot(t, (phi * psi).sum(1), color=c, ls=ls, label=lab)
        ax[1].plot(t, np.array([float(r[C["alpha_new"]]) for r in rows]) * 1e4, color=c, ls=ls)
    for a in ax[:2]:
        a.axvspan(0, 0.5, color="#4c72b0", alpha=0.06, lw=0); a.axvspan(0.5, 1, color="#c44e52", alpha=0.06, lw=0)
        a.axvline(0.5, color="0.4", lw=0.8)
    ax[0].set(xlabel="$T^*$", ylabel="$\\sum_i\\phi_i\\psi_i$", title="living biomass, seed Gauss point")
    ax[0].legend(fontsize=7, frameon=False, title="antibiotic from $T^*=0.5$", title_fontsize=7)
    ax[1].set(xlabel="$T^*$", ylabel="$\\alpha-1$ ($10^{-4}$)", title="growth, seed Gauss point")
    lab = ["(c), no\nantibiotic", "all\n$b_i=1$", "all,\nstronger", "So", "anaerobes"]
    d = np.array([-7.8, -16.1, -27.6, -18.7, -33.2])
    col = ["0.55", "#4c72b0", "#4c72b0", "#dd8452", "#c44e52"]
    ax[2].bar(range(5), d, color=col, width=0.65)
    for i, v in enumerate(d):
        ax[2].annotate(f"{v:.0f} %", (i, v), ha="center", va="top", fontsize=7, xytext=(0, -2), textcoords="offset points")
    ax[2].set_xticks(range(5), lab, fontsize=7)
    ax[2].set(ylabel="seed von Mises vs. (a) [%]", title="seed stress, $T^*=1$", ylim=(-40, 0))
    ax[2].axhline(0, color="0.3", lw=0.8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_meeting1019_switch.pdf")


if __name__ == "__main__":
    switch_on()
    lysis_trace()
    lysis()
    schemes()
    cal5()
    print("wrote", OUT / "fig_meeting1019_schemes.pdf", OUT / "fig_meeting1019_cal5.pdf")
