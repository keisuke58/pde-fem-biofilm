"""plot_klempt_depth.py
======================
Read klempt_extract_{geom}_{cond}.csv (from extract_klempt_odb.py)
and produce a 3x2 comparison figure for thesis §5.2.

Layout:
    col 0: tooth (P23 natural tooth)
    col 1: implant (Ti helical screw)
    row 0: alpha  (growth)
    row 1: E_gated [MPa]
    row 2: sigma_Mises [MPa]

x-axis: phi_gate (depth proxy, inverted: inner=left, outer=right)

Usage:
    PATH=/home/nishioka/texlive/2025/bin/x86_64-linux:$PATH python plot_klempt_depth.py
    python plot_klempt_depth.py --out thesis_fig_5_2.pdf  # (no LaTeX)
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── thesis_style ─────────────────────────────────────────────────────────────
THESIS_STYLE_DIR = os.path.expanduser("~/LUHsummer26/0820_Stochastic_FEM_oral_TBD")
sys.path.insert(0, THESIS_STYLE_DIR)
try:
    from thesis_style import use as ts_use, clean_ax, PALETTE
    FIGSIZE = ts_use(width_frac=1.0, aspect=0.88)
    HAS_THESIS = True
except ImportError:
    HAS_THESIS = False
    PALETTE = {
        "ch": "#2196A6", "dh": "#D95F02",
        "cs": "#5ab4ac", "ds": "#f1a340",
    }
    FIGSIZE = (7.0, 6.5)
    def clean_ax(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

HERE = os.path.dirname(os.path.abspath(__file__))

# 4 solid colors (no dashes) — static vs hobic distinguished by lightness
CONDS = {
    "commensal_hobic":  ("CH", PALETTE["ch"], "-",  2.0),
    "dysbiotic_hobic":  ("DH", PALETTE["dh"], "-",  2.0),
    "commensal_static": ("CS", PALETTE["cs"], "--", 1.2),
    "dysbiotic_static": ("DS", PALETTE["ds"], "--", 1.2),
}

GEOMS  = ["tooth", "implant"]
N_BINS = 20


def load_csv(geom, cond):
    path = os.path.join(HERE, "klempt_extract_%s_%s.csv" % (geom, cond))
    if not os.path.exists(path):
        return None
    return np.genfromtxt(path, delimiter=",", names=True)


def bin_by_phi_gate(data, n_bins=N_BINS):
    edges   = np.linspace(0.0, 1.0, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    alpha_m  = np.full(n_bins, np.nan)
    egated_m = np.full(n_bins, np.nan)
    mises_m  = np.full(n_bins, np.nan)
    for k in range(n_bins):
        mask = (data["phi_gate"] >= edges[k]) & (data["phi_gate"] < edges[k + 1])
        if mask.sum() == 0:
            continue
        alpha_m[k]  = np.nanmean(data["alpha"][mask])
        egated_m[k] = np.nanmean(data["E_gated_MPa"][mask])
        mises_m[k]  = np.nanmean(data["sigma_mises_MPa"][mask])
    return centers, alpha_m, egated_m, mises_m


def sigma_ratio_inner(geom, q=0.7):
    """CH/DH sigma_Mises ratio at innermost elements (phi_gate > q)."""
    d_ch = load_csv(geom, "commensal_hobic")
    d_dh = load_csv(geom, "dysbiotic_hobic")
    if d_ch is None or d_dh is None:
        return None
    m_ch = d_ch["phi_gate"] > q
    m_dh = d_dh["phi_gate"] > q
    if m_ch.sum() == 0 or m_dh.sum() == 0:
        return None
    denom = np.nanmean(d_dh["sigma_mises_MPa"][m_dh])
    if denom < 1e-20:
        return None
    return np.nanmean(d_ch["sigma_mises_MPa"][m_ch]) / denom


def make_figure(out_path):
    fig, axes = plt.subplots(3, 2, figsize=FIGSIZE, sharex=True)

    row_labels = [
        r"$\alpha$ (growth)",
        r"$E_\mathrm{gated}$ [MPa]",
        r"$\sigma_\mathrm{Mises}$ [MPa]",
    ]
    col_titles = ["Natural tooth (P23)", "Ti implant screw"]

    for col, geom in enumerate(GEOMS):
        for cond, (label, color, ls, lw) in CONDS.items():
            data = load_csv(geom, cond)
            if data is None:
                print("Missing: %s %s" % (geom, cond))
                continue
            x, alpha, egated, mises = bin_by_phi_gate(data)

            kw_line = dict(color=color, ls=ls, lw=lw, zorder=2)
            kw_mark = dict(color=color, ls="none", marker="o", ms=4,
                           markeredgewidth=0.4, markeredgecolor="white", zorder=3)
            for row, y in enumerate([alpha, egated, mises]):
                axes[row, col].plot(x, y, **kw_line)
                axes[row, col].plot(x, y, **kw_mark)

        axes[0, col].set_title(col_titles[col], pad=4)
        for row in range(3):
            clean_ax(axes[row, col])

    # y-axis labels (left column only)
    for row, label in enumerate(row_labels):
        axes[row, 0].set_ylabel(label)

    # x-axis label + ticks (bottom row)
    for col in range(2):
        axes[2, col].set_xlabel(r"$\phi_\mathrm{gate}$ (depth proxy)")
        axes[2, col].set_xticks([0.0, 0.5, 1.0])
        axes[2, col].set_xticklabels(
            ["outer\n(planktonic)", r"$0.5$", "inner\n(surface)"])

    # invert x-axis once (shared axis)
    axes[0, 0].invert_xaxis()

    # legend: custom handles so both line style AND marker are shown
    legend_handles = [
        Line2D([0,1], [0,0], color=PALETTE["ch"], lw=2.0, ls="-",
               marker="o", ms=4, markeredgewidth=0.4, markeredgecolor="white",
               label="CH"),
        Line2D([0,1], [0,0], color=PALETTE["dh"], lw=2.0, ls="-",
               marker="o", ms=4, markeredgewidth=0.4, markeredgecolor="white",
               label="DH"),
        Line2D([0,1], [0,0], color=PALETTE["cs"], lw=1.2, ls="--",
               marker="o", ms=4, markeredgewidth=0.4, markeredgecolor="white",
               label="CS"),
        Line2D([0,1], [0,0], color=PALETTE["ds"], lw=1.2, ls="--",
               marker="o", ms=4, markeredgewidth=0.4, markeredgecolor="white",
               label="DS"),
    ]
    axes[0, 1].legend(handles=legend_handles, loc="upper right", frameon=False,
                      handlelength=2.8, handletextpad=0.5)

    # σ^CH/σ^DH annotation — both geometries
    for col, geom in enumerate(GEOMS):
        r = sigma_ratio_inner(geom, q=0.7)
        if r is not None:
            axes[2, col].text(
                0.97, 0.95,
                r"$\sigma^\mathrm{CH}/\sigma^\mathrm{DH}=%.1f\times$" % r,
                ha="right", va="top",
                transform=axes[2, col].transAxes,
                fontsize=7 if HAS_THESIS else 7,
                color="0.4")

    fig.tight_layout(h_pad=0.6, w_pad=1.2)
    fig.savefig(out_path, bbox_inches="tight")
    print("Saved: %s" % out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "klempt_depth_comparison.pdf"))
    args = ap.parse_args()
    make_figure(args.out)


if __name__ == "__main__":
    main()
