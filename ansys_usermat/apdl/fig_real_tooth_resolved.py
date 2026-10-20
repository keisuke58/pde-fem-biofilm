"""fig_real_tooth_resolved.py -- summary figure of the real-molar runs with the converged stress: (a) growth alpha - 1
averaged over the film by height, (b) nutrient c by height (nutrient runs), (c) film mean von Mises of the end state
from the film submodel (h = 0.083 mm), relative to the run without nutrient and gate (REAL_TOOTH_RESOLVED.md).

    python ansys_usermat/apdl/fig_real_tooth_resolved.py   -> assets/fig_real_tooth_resolved.pdf, .png
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(9)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

R = HERE / "results" / "2026-10-3d"
P = json.load(open(R / "real_tooth_physics.json"))["runs"]
S = json.load(open(R / "real_tooth_physics_resolved.json"))
RUNS = [("base", "no nutrient, no gate", "k", "-"), ("nutrient", "nutrient", "#1f6fb4", "-"),
        ("gate", "gate", "#c0392b", "--"), ("nutrient_gate", "nutrient + gate", "#2a9d3f", "-."),
        ("nutrient_t2", "nutrient, Table 2 ratio", "#1f6fb4", ":")]

fig = plt.figure(figsize=(6.3, 2.9))
ax = [fig.add_axes([0.08, 0.32, 0.22, 0.6]), fig.add_axes([0.38, 0.32, 0.22, 0.6]),
      fig.add_axes([0.70, 0.32, 0.28, 0.6])]
for tag, lab, col, ls in RUNS:
    pr = P[tag]["profiles"]
    ax[0].plot(np.array(pr["am1"]) * 1e3, pr["z_field_mm"], color=col, ls=ls, lw=1.3, label=lab)
    if P[tag]["nutrient_on"]:
        ax[1].plot(pr["c"], pr["z_field_mm"], color=col, ls=ls, lw=1.3)
for a in ax[:2]:
    a.axhspan(-2.1, -1.1, color="0.88", lw=0, zorder=0)
    a.set_ylim(-3.2, -0.9)
ax[0].set_ylabel("height $z$ (mm)")
ax[0].text(0.30, -1.25, "cervical band", fontsize=7, va="center", color="0.35", transform=ax[0].get_yaxis_transform())
ax[0].set_xlabel(r"$\alpha-1$ ($10^{-3}$), film mean")
ax[1].set_xlabel("nutrient $c$, film mean")
ax[1].set_xlim(0, 1.0)
b = S["base"]["measures"]["film vM mean"]
vals = [S[t]["measures"]["film vM mean"] / b for t, *_ in RUNS]
y = np.arange(len(RUNS))[::-1]
ax[2].barh(y, vals, color=[c for *_, c, _ in RUNS], height=0.6)
for yi, v in zip(y, vals):
    ax[2].text(v + 0.02, yi, f"{v:.2f}", va="center", fontsize=8)
ax[2].set_yticks(y, ["base", "N", "G", "N + G", "N, Table 2"], fontsize=7.5)
ax[2].set_xlim(0, 1.25)
ax[2].set_xlabel("film mean von Mises / base")
for a in ax:
    a.grid(False)
for a, l in zip(ax, "abc"):
    a.text(0.02, 1.03, f"({l})", transform=a.transAxes, fontsize=9)
h, l = ax[0].get_legend_handles_labels()
l = ["no nutrient, no gate (base)", "nutrient (N)", "gate (G)", "nutrient + gate", "nutrient, Table 2 ratio"]
fig.legend(h, l, fontsize=7.5, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01))
out = ROOT / "assets" / "fig_real_tooth_resolved"
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
print("ratios", [round(v, 2) for v in vals])
