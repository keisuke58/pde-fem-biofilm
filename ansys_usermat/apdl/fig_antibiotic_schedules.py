"""fig_antibiotic_schedules.py -- how the antibiotic is given in time: constant, pulses, sine (the same amount from T* = 1/3,
three periods of 2/9), separate start, 64^3, Fritsch et al. 2025 case I, alpha*_0 = 100 (Keio continuation, not the thesis).

(a) held antibiotic value a_b(T*): 1; 2 in the first half of each period and 0 in the second; 1 + sin(2 pi (T* - 1/3)/P).
(b) seed living volume fraction sum_i phibar_i; (c) seed antibiotic a. Growth law of Eq. 36 (solid) and growth from the
living biomass (--growth-living-energy, dashed; the living fraction is the same for both laws). Seed means at T* = 1 of
alpha - 1 and von Mises for both laws are printed (table of papers/coupling_paper/SINE_SCHEDULES.ja.md).

    python ansys_usermat/apdl/fig_antibiotic_schedules.py -> assets/fig_antibiotic_schedules.pdf/.png
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402

figstyle.apply()
R = HERE / "results" / "2026-10-3d"
P0 = "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv"
RUNS = {  # (schedule, law) -> file; the Eq. 36 pulse and sine runs are without dc/dt (within 0.1 %, CLAUDE.md)
    ("constant", "eq36"): P0 + "_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_wgmirror_fast_separate_ab.json",
    ("pulses", "eq36"): P0 + "_a100_gab1_gc1_tab0.333333_abpulse0.2222_tvd_wgmirror_fast_separate_ab.json",
    ("sine", "eq36"): P0 + "_a100_gab1_gc1_tab0.333333_absine0.2222_tvd_wgmirror_fast_separate_ab.json",
    ("constant", "living"): P0 + "_growlivE_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_wgmirror_fast_separate_ab.json",
    ("pulses", "living"): P0 + "_growlivE_a100_gab1_gc1_tab0.333333_transc1e+08_abpulse0.2222_tvd_wgmirror_fast_separate_ab.json",
    ("sine", "living"): P0 + "_growlivE_a100_gab1_gc1_tab0.333333_transc1e+08_absine0.2222_tvd_wgmirror_fast_separate_ab.json",
}
COL = {"constant": "0.2", "pulses": "#e8590c", "sine": "#1c7ed6"}
T_AB, P = 1 / 3, 2 / 9


def a_b(T, sched):
    if T <= T_AB:
        return 0.0
    s = (T - T_AB) / P
    return {"constant": 1.0, "pulses": 2.0 if s - math.floor(s) < 0.5 else 0.0, "sine": 1.0 + math.sin(2 * math.pi * s)}[sched]


fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), constrained_layout=True)
tt = np.linspace(0, 1, 2001)
for sched in ("constant", "pulses", "sine"):
    axs[0].plot(tt, [a_b(T, sched) for T in tt], color=COL[sched], lw=1.1, label=sched)
    for law, ls in (("eq36", "-"), ("living", "--")):
        r = json.load(open(R / RUNS[(sched, law)]))["runs"]["separate_ab"]
        t = np.asarray(r["t"], float)
        if law == "eq36":
            axs[1].plot(t, r["living_seed"], color=COL[sched], lw=1.1, marker="o", ms=2)
            axs[2].plot(t, r["a_seed"], color=COL[sched], lw=1.1, marker="o", ms=2)
        print(f"{sched:9s} {law:7s} living end {r['living_seed'][-1]:.3f} min {min(r['living_seed']):.3f} "
              f"alpha-1 {r['seed_am1']*1e4:.2f} vM {r['seed_vm_Pa']*1e4:.2f}")
for ax, (lab, yl) in zip(axs, (("(a) held antibiotic $a_b$", 2.1), (r"(b) seed $\Sigma_i\bar\phi_i$", 1.02), ("(c) seed antibiotic $a$", 2.1))):
    ax.axvline(T_AB, color="0.6", lw=0.7, ls=(0, (1, 2)))
    ax.set_xlim(0, 1); ax.set_ylim(0, yl); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title(lab, fontsize=8, loc="left")
axs[0].legend(frameon=False, fontsize=6.5, loc="upper left")
A = ROOT / "assets"
fig.savefig(A / "fig_antibiotic_schedules.pdf"); fig.savefig(A / "fig_antibiotic_schedules.png", dpi=400)
print("wrote", A / "fig_antibiotic_schedules.png")
