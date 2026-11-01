#!/usr/bin/env python3
"""fig_validation_volume.py -- model biofilm volume V(t)/V(day 1) of the flat layer with the calibrated five-species
point model (CH, CS, DH, DS), for the front speed of the paper, one tenth of it, and one tenth tied to the living
fraction. Model results only (ansys_usermat/apdl/results/2026-10-validation/); the measured volumes of Heine et al. 2025
are described in the text, not plotted here.

    python papers/coupling_paper/fig_validation_volume.py -> papers/coupling_paper/fig_validation_volume.pdf
"""
import glob
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

figstyle.apply()
RES = ROOT / "ansys_usermat/apdl/results/2026-10-validation"
PANELS = [("front1", False, r"(a) $r$ of the paper"), ("front0.1", False, r"(b) $r/10$"),
          ("front0.1", True, r"(c) $r/10$, speed $\times\sum_i\bar\phi_i$")]
COL = {"CH": "#1f4e8c", "CS": "#6b9bd1", "DH": "#c2571a", "DS": "#e3a26f"}
fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), sharey=True, constrained_layout=True)
for ax, (fr, fl, title) in zip(axs, PANELS):
    for c in ("CH", "CS", "DH", "DS"):
        pat = f"coupled_living_consumption_{fr}_ecocal5_{c}_dtc0.0025_phimin0.001_fortran_gc1*_tvd_wgmirror_fast_layer.json"
        fs = [f for f in glob.glob(str(RES / pat)) if ("frontliv" in f) == fl and "growliv" not in f and "_h0" not in f]
        r = list(json.load(open(fs[0]))["runs"].values())[0]
        t = [1 + 20 * x for x in r["t"]]
        V = r["phi_total_mm3"]
        ax.plot(t, [v / V[0] for v in V], color=COL[c], lw=1.2, label=c)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("day", fontsize=8)
    ax.set_xticks([1, 6, 11, 16, 21])
    ax.tick_params(labelsize=7.5)
axs[0].set_ylabel(r"$V/V(\mathrm{day}\ 1)$", fontsize=8)
axs[0].legend(fontsize=7, frameon=False)
fig.savefig(HERE / "fig_validation_volume.pdf")
print("wrote", HERE / "fig_validation_volume.pdf")
