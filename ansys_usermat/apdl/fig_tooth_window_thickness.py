"""fig_tooth_window_thickness.py -- film thickness of the cervical window (TOOTH_WINDOW.md, variant 2: saliva well mixed)
over T*, against the plaque thickness measured in situ on enamel (Al-Ahmad et al. 2007: about 15 um after one day and
50 um after seven days), to read off a time scale: how many days one T* corresponds to.

    python ansys_usermat/apdl/fig_tooth_window_thickness.py -> assets/fig_tooth_window_thickness.pdf/.png, prints the fit

Thickness = mean over the tooth surface of the film (phi_3D > 0.5) thickness, hist "film" -> "thickness_mean_mm" (code mm;
1 code mm = 10 um). The time scale is read from the slope of the late, nearly linear part (T* = 0.5 to 1) against the
mean plaque growth rate (50 - 15) / 6 um per day; both are linear fits, so this is an order of magnitude.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply()
RES = HERE / "results" / "2026-10-3d"
TAG = "front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_superbee_wgmirror_obststl_win{w}_stlx100_salmix0.05_box30_amg_{case}"
RUNS = [(128, "mixed", "window_ab", "mixed, antibiotic", "C0", "-"), (192, "mixed", "window_ab", "mixed, antibiotic (192³)", "C0", ":"),
        (128, "layered", "window_ab", "layered, antibiotic", "C1", "-"), (128, "mixed", "window", "mixed, no antibiotic", "C2", "-")]
PLAQUE_RATE = (50.0 - 15.0) / 6.0          # um per day, Al-Ahmad et al. 2007 (1 and 7 days)


def series(nf, w, case):
    d = json.load(open(RES / f"coupled_living_consumption_nf{nf}_{TAG.format(w=w, case=case)}.json"))
    r = list(d["runs"].values())[0]
    t = np.asarray(r["t"], float)
    h = 10.0 * np.asarray(r["film"]["thickness_mean_mm"], float)   # um
    n = min(len(t), len(h))
    return t[:n], h[:n]


fig, ax = plt.subplots(figsize=(3.3, 2.4))
out = {}
for nf, w, case, lab, col, ls in RUNS:
    try:
        t, h = series(nf, w, case)
    except FileNotFoundError:
        continue
    ax.plot(t, h, color=col, ls=ls, lw=1.3, label=lab)
    m = t >= 0.5
    k, b = np.polyfit(t[m], h[m], 1)
    out[f"{lab} {nf}^3"] = dict(slope_um_per_T=round(k, 1), days_per_T=round(k / PLAQUE_RATE, 1),
                                 T_at_50um=round(float(np.interp(50.0, h, t)), 3), h_T1=round(float(h[-1]), 1))
for hh, dd in ((15, "1 day"), (50, "7 days")):
    ax.axhline(hh, color="0.5", lw=0.6, ls="--")
    ax.text(0.99, hh + 1.0, f"plaque, {dd}", ha="right", va="bottom", fontsize=7, color="0.35", transform=ax.get_yaxis_transform())
ax.axvline(1 / 3, color="0.6", lw=0.6, ls=":")
ax.set_xlim(0, 1); ax.set_ylim(0, 75)
ax.set_xlabel(r"$T^*$", fontsize=9); ax.set_ylabel(r"film thickness [$\mu$m]", fontsize=9)
ax.tick_params(labelsize=8); ax.legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout()
A = HERE.parents[1] / "assets"
fig.savefig(A / "fig_tooth_window_thickness.pdf"); fig.savefig(A / "fig_tooth_window_thickness.png", dpi=400)
print(json.dumps(out, indent=1))
