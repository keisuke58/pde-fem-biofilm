"""fig_nutrient_sine.py -- a nutrient supply that oscillates in time (the form of Klempt et al. 2026, Table 2, case 3),
with the two species of Fritsch et al. 2025, case I, and no antibiotic.

Top row, the point model alone (no field), 1000 steps of 1e-4, T* = 1 the whole run (s = 0.1):
  (a) c* = 100 (the coupled examples), (b) c* = 50 + 50 sin(500 t), t the point-model time, (c) c* = 50, the mean of (b).
  Solid: phi_i; dashed: phibar_i = phi_i psi_i; dotted: sum_i phibar_i.
Bottom row, the coupled model, separate start, 64^3, held nutrient 1, 0.5 + 0.5 sin(500 t) and 0.5 (c* = 100 c), seed
means over time: (d) living volume fraction, (e) nutrient, (f) share of species A. The held value is updated once per
outer step of 0.025 T* (five values per period of 0.126 T*), so the coupled "sine" is a sampled one.

    python ansys_usermat/apdl/fig_nutrient_sine.py -> assets/fig_nutrient_sine.pdf/.png
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
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(ROOT / "JAXFEM"))
import figstyle  # noqa: E402
import klempt2026_reproduction as KR  # noqa: E402
from fritsch2025_cases import CASES  # noqa: E402

figstyle.apply()
base = dict(CASES["fritsch_2sp"], alpha_star=0.0)
KR.DT = base["dt"]
STEPS = base["steps"]
SUP = [("constant, $c^*=100$", 100.0, "0.2"),
       (r"$c^*=50+50\sin(500\,t)$", lambda k: 50.0 + 50.0 * math.sin(500.0 * k * base["dt"]), "#c2255c"),
       ("constant, $c^*=50$ (its mean)", 50.0, "#1c7ed6")]
COL = ("#1c6fb0", "#c0392b")
R3 = HERE / "results" / "2026-10-3d-case3"
R0 = HERE / "results" / "2026-10-3d"
COUPLED = [R0 / "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_separate.json",
           R3 / "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_csine_tvd_wgmirror_fast_separate.json",
           R3 / "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_cconst0.5_tvd_wgmirror_fast_separate.json"]

fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.4), constrained_layout=True)
for j, (lab, cs, _) in enumerate(SUP):
    r = KR.simulate(dict(base, c_star=cs))
    t = np.arange(r["phi"].shape[0]) / STEPS
    pb = r["phi"] * r["psi"]
    ax = axs[0, j]
    for i in range(2):
        ax.plot(t, r["phi"][:, i], color=COL[i], lw=1.2, label=rf"$\phi_{i+1}$")
        ax.plot(t, pb[:, i], color=COL[i], lw=1.0, ls="--", label=rf"$\bar\phi_{i+1}$")
    ax.plot(t, pb.sum(1), color="0.25", lw=0.9, ls=":", label=r"$\Sigma_i\bar\phi_i$")
    ax.set_title(f"({'abc'[j]}) " + lab, fontsize=7.5, loc="left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
    print("0D", lab, "end phi", np.round(r["phi"][-1], 3), "sum phibar", round(float(pb[-1].sum()), 3))
axs[0, 0].set_ylabel("point model alone", fontsize=8)
axs[0, 0].legend(frameon=False, fontsize=6, loc="center right")
keys = (("living_seed", r"seed $\Sigma_i\bar\phi_i$"), ("c_seed", "seed nutrient $c$"), ("share_a_seed", r"seed $\phi_1/(\phi_1+\phi_2)$"))
for (lab, _, col), f in zip(SUP, COUPLED):
    r = json.load(open(f))["runs"]["separate"]
    t = np.asarray(r["t"], float)
    for j, (k, yl) in enumerate(keys):
        axs[1, j].plot(t, r[k], color=col, lw=1.1, marker="o" if "sin" in lab else None, ms=2.2, label=lab.replace("$c^*=", "$c^*=").replace(" (its mean)", ""))
    print("coupled", lab, "alpha-1", round(r["seed_am1"] * 1e4, 3), "vM", round(r["seed_vm_Pa"] * 1e4, 3), "living", round(r["living_seed"][-1], 3),
          "c mean", round(float(np.mean(r["c_seed"][1:])), 3))
for j, (k, yl) in enumerate(keys):
    ax = axs[1, j]
    ax.set_title(f"({'def'[j]}) coupled, " + yl, fontsize=7.5, loc="left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
axs[1, 2].legend(frameon=False, fontsize=6.3, loc="lower left")
A = ROOT / "assets"
fig.savefig(A / "fig_nutrient_sine.pdf"); fig.savefig(A / "fig_nutrient_sine.png", dpi=400)
print("wrote", A / "fig_nutrient_sine.png")
