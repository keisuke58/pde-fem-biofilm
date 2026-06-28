"""
plot_at2_comparison.py
======================
Thesis figure: AT2 phase-field fracture in CH vs DH biofilm conditions.

Left : d_max(t) curves for 3 conditions + t_crit vertical markers.
Right: Horizontal bar chart of t_crit with k_ratio annotation.

Output: fig_at2_comparison.pdf
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from thesis_style import use, PALETTE, clean_ax
from phase_field_at2_2d_sparse import (
    run, t_crit as find_tcrit,
    K_EFF_CH, K_EFF_DH, K_RATIO, N_STEPS,
)

# ── Run 3 conditions ──────────────────────────────────────────────────────────
NX, NZ   = 80, 50
N_STEPS2 = 200
PURPLE   = "#8B6BB1"

CASES = [
    (K_EFF_CH, "uniform",      "CH",      PALETTE["ch"],   "-"),
    (K_EFF_DH, "uniform",      "DH-unif", PALETTE["dh"],   "-"),
    (K_EFF_DH, "pg_substrate", "DH-Pg",   PURPLE,          "--"),
]

all_res = {}
for k, prof, lbl, col, ls in CASES:
    print(f"  Running {lbl} ...", flush=True)
    all_res[lbl] = run(k, prof, NX, NZ, n_steps=N_STEPS2, label="")

tc_ch  = find_tcrit(all_res["CH"]["results"])
tc_dhu = find_tcrit(all_res["DH-unif"]["results"])
tc_dhp = find_tcrit(all_res["DH-Pg"]["results"])

print(f"\n  t_crit: CH={tc_ch:.1f} s  DH-unif={tc_dhu:.1f} s  DH-Pg={tc_dhp:.1f} s")
print(f"  t_crit(DH-unif)/t_crit(CH) = {tc_dhu/tc_ch:.3f}  (theory k_ratio={K_RATIO:.2f})")

# ── Figure layout ─────────────────────────────────────────────────────────────
figsize = use(width_frac=1.0, aspect=0.42)
fig = plt.figure(figsize=figsize)
gs  = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[3, 2], wspace=0.38)
ax_curve = fig.add_subplot(gs[0])
ax_bar   = fig.add_subplot(gs[1])

# ── LEFT: d_max vs time ───────────────────────────────────────────────────────
THRESH = 0.5
for k, prof, lbl, col, ls in CASES:
    res = all_res[lbl]["results"]
    ts  = [r["t"]     for r in res]
    ds  = [r["d_max"] for r in res]
    ax_curve.plot(ts, ds, color=col, ls=ls, lw=1.3, label=lbl)

    tc = find_tcrit(res)
    if tc is not None:
        ax_curve.axvline(tc, color=col, ls=":", lw=0.9, alpha=0.65)

ax_curve.axhline(THRESH, color="0.4", ls="--", lw=0.7, label="$d=0.5$")

# Double-arrow for the 2.00 ratio
y_arr = 0.82
ax_curve.annotate("",
    xy=(tc_dhu, y_arr), xytext=(tc_ch, y_arr),
    arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.9, mutation_scale=8))
ax_curve.text((tc_ch + tc_dhu) / 2, y_arr + 0.025,
              r"$\times 2.00$", ha="center", va="bottom", fontsize=8, color="0.25")

# t_crit labels
for tc, lbl, col in [(tc_ch, "CH", PALETTE["ch"]),
                     (tc_dhu, "DH-unif", PALETTE["dh"]),
                     (tc_dhp, "DH-Pg", PURPLE)]:
    if tc is not None:
        ax_curve.text(tc - 2, 0.56,
                      f"{tc:.0f}" + r"\,s", ha="right", fontsize=7, color=col)

ax_curve.set_xlabel(r"Time $t$ [s]")
ax_curve.set_ylabel(r"$d_{\max}(t)$")
ax_curve.set_xlim(0, max(tc_dhu * 1.1, 180))
ax_curve.set_ylim(0, 1.05)
ax_curve.legend(loc="upper left", frameon=False, fontsize=7)
clean_ax(ax_curve)
ax_curve.text(-0.08, 1.04, r"\textbf{(a)}", transform=ax_curve.transAxes,
              fontsize=9, va="top")

# ── RIGHT: horizontal bar chart of t_crit ─────────────────────────────────────
bar_labels = ["DH-Pg",    "CH",        "DH-unif"]
bar_vals   = [tc_dhp,      tc_ch,       tc_dhu]
bar_colors = [PURPLE,      PALETTE["ch"], PALETTE["dh"]]
y_pos      = [2,           1,             0]

ax_bar.barh(y_pos, bar_vals, color=bar_colors, height=0.5,
            edgecolor="none", alpha=0.88)

for yp, tv, col in zip(y_pos, bar_vals, bar_colors):
    ax_bar.text(tv + 2, yp, f"{tv:.1f}" + r"\,s",
                va="center", fontsize=8, color=col)

ax_bar.text((tc_ch + tc_dhu) / 2, -0.18,
            r"$t_{\rm crit}^{\rm DH}/t_{\rm crit}^{\rm CH}=2.00$",
            ha="center", va="top", fontsize=7.5, color="0.25",
            transform=ax_bar.get_xaxis_transform())

ax_bar.set_yticks(y_pos)
ax_bar.set_yticklabels(bar_labels, fontsize=8)
ax_bar.set_xlabel(r"$t_{\rm crit}$ [s]")
ax_bar.set_xlim(0, tc_dhu * 1.25)
clean_ax(ax_bar)
ax_bar.text(-0.18, 1.04, r"\textbf{(b)}", transform=ax_bar.transAxes,
            fontsize=9, va="top")

outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "fig_at2_comparison.pdf")
fig.savefig(outpath)
print(f"\nSaved: {outpath}")
