#!/usr/bin/env python3
"""fig_ab_schedule.py -- Figure 6 of the coupling paper: the antibiotic given over time in three ways with the same
total dose from T* = 1/3 (A-71): constant a_b = 1, pulsed (2 for the first half of each period, 0 for the second) and
sine 1 + sin(2 pi (T* - 1/3)/P), period P = 2/9 (three periods up to T* = 1). Separate example, main settings
(Python TVD, mirrored wall gradients, inherit rule, 64^3). (a) held value a_b at the wall y = +1 mm (the value at the seed
follows it within the outer step, d = 1e8 mm^2/T*); (b) living volume fraction sum_i phibar_i at the seed; (c) seed mean
of alpha - 1 and of the von Mises stress at T* = 1 for the three schedules (each axis spans +-5 % of its mean).

Data: results/2026-10-3d/coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1
_tab0.333333[_abpulse0.2222|_absine0.2222]_tvd_wgmirror_fast_separate_ab.json (runs.separate_ab: t, a_seed,
living_seed, share_a_seed, seed_am1, seed_vm_Pa), made with
    coupled_living_consumption.py --nf 64 --fast --tvd --wall-grad mirror --gc 1 --gab 1 --alpha0 100 --front 1
        --tab 0.333333 --phi-min 0.001 --v2 --fortran --ansys-bc --inherit-living --only separate_ab
        [--ab-schedule pulse|sine --ab-period 0.2222]
(27 Oct 2026; the constant case is the main run of README_coupled_living.md).

    python papers/coupling_paper/fig_ab_schedule.py [--nf 32] -> papers/coupling_paper/fig_ab_schedule.pdf
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

RES = ROOT / "ansys_usermat/apdl/results/2026-10-3d"
NF = int(sys.argv[sys.argv.index("--nf") + 1]) if "--nf" in sys.argv else 64
T_AB, P = 1.0 / 3.0, 0.2222
SCHED = (("const", "constant", "C0", "-"), ("pulse", "pulsed", "C1", "--"), ("sine", "sine", "C2", ":"))


def a_wall(T, sched):
    """the held antibiotic value, as ab_value() of coupled_living_consumption.py"""
    if T <= T_AB + 1e-12:
        return 0.0
    s = (T - T_AB) / P
    if sched == "pulse":
        return 2.0 if (s - math.floor(s)) < 0.5 - 1e-12 else 0.0
    if sched == "sine":
        return 1.0 + math.sin(2.0 * math.pi * s)
    return 1.0


def load(sched):
    pre = "" if NF == 32 else f"nf{NF}_"
    mid = "" if sched == "const" else f"_ab{sched}{P:.4g}"
    f = RES / f"coupled_living_consumption_{pre}front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333{mid}_tvd_wgmirror_fast_separate_ab.json"
    d = json.load(open(f))
    assert d["setting"]["field"] == f"{NF}^3", (f, d["setting"])
    return d["runs"]["separate_ab"]


def main():
    figstyle.apply(size=8.5)
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.25, 1.25, 0.8]))
    Tf = np.linspace(0, 1, 2001)
    out = {}
    for k, (sched, name, col, ls) in enumerate(SCHED):
        r = load(sched)
        t = np.asarray(r["t"])
        axs[0].plot(Tf, [a_wall(T, sched) for T in Tf], color=col, ls=ls, lw=1, label=name)
        axs[1].plot(t, r["living_seed"], color=col, ls=ls, lw=1)
        axs[2].plot([k], [r["seed_am1"] * 1e4], color=col, marker="o", ms=4, ls="none")
        liv = np.asarray(r["living_seed"])
        i_ab = np.searchsorted(t, T_AB)
        out[sched] = (liv[-1], liv[i_ab:].min(), r["share_a_seed"][-1], r["seed_am1"], r["seed_vm_Pa"])
        print(f"{name:9s} {NF}^3: seed living end {liv[-1]:.3f}, min after T*=1/3 {liv[i_ab:].min():.3f}, "
              f"share A {r['share_a_seed'][-1]:.3f}, alpha-1 {r['seed_am1']:.3e}, vM {r['seed_vm_Pa']:.3e} Pa")
    ax2 = axs[2].twinx()
    for k, (sched, name, col, ls) in enumerate(SCHED):
        ax2.plot([k], [out[sched][4] * 1e4], color=col, marker="s", mfc="white", ms=4, ls="none")
    am = [out[s][3] * 1e4 for s, *_ in SCHED]
    vm = [out[s][4] * 1e4 for s, *_ in SCHED]
    print(f"spread over the schedules: alpha-1 {100 * (max(am) - min(am)) / np.mean(am):.1f} %, vM {100 * (max(vm) - min(vm)) / np.mean(vm):.1f} %")
    axs[0].set_ylabel(r"held value $a_b$ at $y = +1$ mm")
    axs[1].set_ylabel(r"seed $\sum_i\bar\phi_i$")
    axs[2].set_ylabel(r"seed $\alpha-1$ [$10^{-4}$], filled circles")
    ax2.set_ylabel(r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa], open squares")
    ax2.grid(False)
    ax2.spines["top"].set_visible(False)
    axs[2].spines["right"].set_visible(True)
    axs[2].set_xticks(range(3))
    axs[2].set_xticklabels([s[1] for s in SCHED], rotation=30, ha="right")
    axs[2].set_xlim(-0.5, 2.5)
    axs[2].set_ylim(0.95 * np.mean(am), 1.05 * np.mean(am)); ax2.set_ylim(0.95 * np.mean(vm), 1.05 * np.mean(vm))
    for ax in axs[:2]:
        ax.set_xlabel("$T^*$")
        ax.set_xlim(0, 1)
        ax.axvline(T_AB, color="0.4", ls=":", lw=0.7)
    axs[0].set_ylim(-0.05, 2.6)
    axs[1].set_ylim(0, 1)
    axs[0].legend(frameon=False, fontsize=7, loc="upper left", ncol=3, columnspacing=0.8, handlelength=1.6)
    for ax, letter in zip(axs, "abc"):
        ax.text(-0.28, 1.03, f"({letter})", transform=ax.transAxes, va="bottom", ha="left")
    out_pdf = HERE / "fig_ab_schedule.pdf"
    fig.savefig(out_pdf)
    print("wrote", out_pdf)


if __name__ == "__main__":
    main()
