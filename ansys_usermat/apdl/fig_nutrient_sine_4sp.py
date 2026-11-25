"""fig_nutrient_sine_4sp.py -- the oscillating nutrient supply with the four species of Klempt et al. 2026, Table 2,
case 3 (A without the 1/2, see JAXFEM/klempt2026_cases.py), no antibiotic. Same layout as fig_nutrient_sine.py
(Keio continuation, the four-species example is kept back until about June 2027).

Top row, the point model alone, 1500 steps of 1e-4, T* = 1 the whole run (s = 0.15):
  (a) c* = 100, (b) c* = 50 + 50 sin(500 t) (the paper's case 3), (c) c* = 50, the mean of (b).
  Solid: phi_i; dashed: phibar_i = phi_i psi_i; dotted: sum_i phibar_i.
Bottom row, the coupled model, separate start (species 1, 2 in the seed, 3, 4 in the corners), 64^3, held nutrient 1,
0.5 + 0.5 sin(500 t) and 0.5 (c* = 100 c), seed means over time: (d) living volume fraction, (e) nutrient,
(f) (phi_1 + phi_2) / sum_i phi_i. The held value is updated once per outer step of 0.005 T*.

    python ansys_usermat/apdl/fig_nutrient_sine_4sp.py -> assets/fig_nutrient_sine_4sp.pdf/.png
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
from klempt2026_cases import CASES, DT  # noqa: E402

figstyle.apply()
base = dict(CASES["4sp_case3"], alpha_star=0.0)
STEPS = base["steps"]
SUP = [("constant, $c^*=100$", 100.0, "0.2", "_transc1e+08_k26c3A1_tvd"),
       (r"$c^*=50+50\sin(500\,t)$", base["c_star"], "#c2255c", "_transc1e+08_k26c3A1_csine_tvd"),
       ("constant, $c^*=50$ (its mean)", 50.0, "#1c7ed6", "_transc1e+08_k26c3A1_cconst0.5_tvd")]
COL = ("#1c6fb0", "#4dabf7", "#c0392b", "#f08c00")   # species 1, 2 (seed, blue), 3, 4 (corners, red/orange)
R3 = HERE / "results" / "2026-10-3d-case3"
P0 = "coupled_living_consumption_nf64_front1_dtc0.005_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333"


def coupled(tag):
    f = R3 / f"{P0}{tag}_wgmirror_fast_separate.json"
    if not f.exists():
        hits = sorted(R3.glob(f"{P0.split('_dtc')[0]}*{tag.split('_transc1e+08')[1].replace('_tvd', '')}*_separate.json"))
        hits = [h for h in hits if "k26c3A1" in h.name and "startpos" not in h.name]
        f = hits[0]
    return f


fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.4), constrained_layout=True)
for j, (lab, cs, _, _) in enumerate(SUP):
    r = KR.simulate(dict(base, c_star=cs))
    t = np.arange(r["phi"].shape[0]) / STEPS
    pb = r["phi"] * r["psi"]
    ax = axs[0, j]
    for i in range(4):
        ax.plot(t, r["phi"][:, i], color=COL[i], lw=1.1, label=rf"$\phi_{i+1}$")
        ax.plot(t, pb[:, i], color=COL[i], lw=0.9, ls="--")
    ax.plot(t, pb.sum(1), color="0.25", lw=0.9, ls=":", label=r"$\Sigma_i\bar\phi_i$")
    ax.set_title(f"({'abc'[j]}) " + lab, fontsize=7.5, loc="left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
    print("0D", lab, "end phi", np.round(r["phi"][-1], 3), "phibar", np.round(pb[-1], 3), "sum phibar", round(float(pb[-1].sum()), 3))
axs[0, 0].set_ylabel("point model alone", fontsize=8)
axs[0, 0].legend(frameon=False, fontsize=6, loc="upper right", ncol=1)
keys = (("living_seed", r"seed $\Sigma_i\bar\phi_i$"), ("c_seed", "seed nutrient $c$"),
        ("share_a_seed", r"seed $(\phi_1+\phi_2)/\Sigma_i\phi_i$"))
for lab, _, col, tag in SUP:
    f = coupled(tag)
    r = json.load(open(f))["runs"]["separate"]
    t = np.asarray(r["t"], float)
    for j, (k, _) in enumerate(keys):
        axs[1, j].plot(t, r[k], color=col, lw=1.1, label=lab.replace(" (its mean)", ""))
    print("coupled", f.name[27:], "\n   alpha-1", round(r["seed_am1"] * 1e4, 3), "vM", round(r["seed_vm_Pa"] * 1e4, 3),
          "living", round(r["living_seed"][-1], 3), "share", round(r["share_a_seed"][-1], 3),
          "c mean", round(float(np.mean(r["c_seed"][1:])), 3))
for j, (k, yl) in enumerate(keys):
    ax = axs[1, j]
    ax.set_title(f"({'def'[j]}) coupled, " + yl, fontsize=7.5, loc="left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.set_xlabel(r"$T^*$", fontsize=8); ax.tick_params(labelsize=7)
axs[1, 2].legend(frameon=False, fontsize=6.3, loc="lower left")
A = ROOT / "assets"
fig.savefig(A / "fig_nutrient_sine_4sp.pdf"); fig.savefig(A / "fig_nutrient_sine_4sp.png", dpi=400)
print("wrote", A / "fig_nutrient_sine_4sp.png")
