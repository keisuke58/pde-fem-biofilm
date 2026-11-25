"""fig_c_start_zero.py -- the nutrient started at c = 0 everywhere except the held layer (Oliver's element, USolBeg)
against the start used so far (c solved from Eq. 35 for the start state). Oliver, 25 Nov 2026: "c = 0.93 in the seed at
t = 0 is not possible". 64^3, the main settings, with dc/dt (--transient-c 1e8), --c-init-zero.
(a) How c fills the cube from c = 0 (separate start, the start state, small backward-Euler steps): mean and minimum
over the seed, and the minimum over the cube, against T* on a log axis. (b) Seed nutrient over the run for both starts
and both examples: the two starts lie on top of each other after the first outer step. Prints the seed values at T* = 1.

    python ansys_usermat/apdl/fig_c_start_zero.py -> assets/fig_c_start_zero.pdf/.png
"""
import json
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
P = "coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_{mx}a100_gab1_gc1_tab0.333333_transc1e+08{z}_tvd_wgmirror_fast_{k}.json"
EX = [("separate", "", "separate_ab", "#3b5bdb"), ("mixed", "mixeq_", "mixed_ab", "#7048a8")]
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5), constrained_layout=True)
for ex, mx, k, col in EX:
    fb = R / P.format(mx=mx, z="", k=k)
    if not fb.exists():                       # the mixed run with dc/dt at outer step 0.025 exists as the video run (snapshots only)
        fb = R / P.format(mx=mx, z="", k=k).replace("tab0.333333_transc", "tab0.333333_snapsvm_transc")
    base = json.load(open(fb))["runs"][k]
    zero = json.load(open(R / P.format(mx=mx, z="_c0zero", k=k)))["runs"][k]
    if ex == "separate":
        cf = zero["c_fill"]
        t = np.array(cf["t"][1:])
        axs[0].semilogx(t, cf["seed_mean"][1:], color="#1c7ed6", lw=1.3, label="seed mean")
        axs[0].semilogx(t, cf["seed_min"][1:], color="#1c7ed6", lw=1.0, ls="--", label="seed minimum")
        axs[0].semilogx(t, cf["cube_min"][1:], color="0.45", lw=1.0, ls=":", label="cube minimum")
        axs[0].axhline(cf["quasi_static_seed_mean"], color="#c92a2a", lw=0.8, ls=(0, (4, 2)), label="start used so far (seed mean)")
        i = next(j for j, v in enumerate(cf["seed_mean"]) if v > 0.99 * cf["quasi_static_seed_mean"])
        print(f"fill: seed mean within 1 % of the quasi-static value at T* = {cf['t'][i]:.2e}; quasi-static {cf['quasi_static_seed_mean']:.3f}")
    for r, ls, lab in ((base, "-", "c solved at the start"), (zero, "--", "c = 0 at the start")):
        axs[1].plot(r["t"], r["c_seed"], color=col, ls=ls, lw=1.2, label=f"{ex}, {lab}")
    print(f"{ex:8s} T*=1 seed: living {base['living_seed'][-1]:.4f} / {zero['living_seed'][-1]:.4f}, share {base['share_a_seed'][-1]:.4f} / "
          f"{zero['share_a_seed'][-1]:.4f}, c {base['c_seed'][-1]:.4f} / {zero['c_seed'][-1]:.4f}, alpha-1 {base['seed_am1']*1e4:.3f} / "
          f"{zero['seed_am1']*1e4:.3f}e-4, vM {base['seed_vm_Pa']*1e4:.3f} / {zero['seed_vm_Pa']*1e4:.3f}e-4 Pa; c at T*=0.025: "
          f"{base['c_seed'][1]:.3f} / {zero['c_seed'][1]:.3f}")
axs[0].set_xlim(1e-11, 1e-5); axs[0].set_ylim(0, 1.02); axs[0].set_xlabel(r"$T^*$ (from $c = 0$)", fontsize=8)
axs[0].set_ylabel("$c$", fontsize=8); axs[0].set_title("(a) filling from $c = 0$, separate start", fontsize=8, loc="left")
axs[0].legend(frameon=False, fontsize=6, loc="upper left")
axs[1].axvline(1 / 3, color="0.5", ls=(0, (1, 2)), lw=0.8)
axs[1].set_xlim(0, 1); axs[1].set_ylim(0, 1.02); axs[1].set_xlabel(r"$T^*$", fontsize=8)
axs[1].set_title("(b) seed nutrient over the run, with antibiotic", fontsize=8, loc="left")
axs[1].legend(frameon=False, fontsize=6, loc="lower right")
for ax in axs:
    ax.tick_params(labelsize=7)
A = ROOT / "assets"
fig.savefig(A / "fig_c_start_zero.pdf"); fig.savefig(A / "fig_c_start_zero.png", dpi=400)
print("wrote", A / "fig_c_start_zero.png")
