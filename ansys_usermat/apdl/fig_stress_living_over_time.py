"""fig_stress_living_over_time.py -- seed means over time of the living volume fraction and of the von Mises stress, with
antibiotic, separate and mixed start (64^3, the --snap-stress runs of the videos). For the meeting of 30 Nov 2026: the
stress keeps rising after the antibiotic while the living biomass falls, because the growth alpha follows phi_3D.

    python ansys_usermat/apdl/fig_stress_living_over_time.py   -> assets/fig_stress_living_over_time.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402

C8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
RUNS = [("separate", "separate_ab", "", "#3b5bdb"), ("mixed", "mixed_ab", "_mixeq", "#7048a8")]


def main():
    figstyle.apply()
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.3), constrained_layout=True)
    for lab, case, mix, col in RUNS:
        z = np.load(ROOT / "figdata" / (f"coupled_living_{case}_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv{mix}_a100"
                                        f"_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_{case}.npz"))
        ts = sorted({float(k.split("_", 1)[1]) for k in z.files if k.startswith("vm16_")})
        s64, s16 = W.seed_mask(64, C8), W.seed_mask(z[f"vm16_{ts[0]}"].shape[0], C8)
        key = lambda n, t: f"{n}_{t}" if f"{n}_{t}" in z.files else f"{n}_{t:g}"
        L = [float(np.nan_to_num(z[key("L", t)])[s64].mean()) for t in ts]
        V = [1e4 * float(np.nan_to_num(z[key("vm16", t)])[s16].mean()) for t in ts]
        axs[0].plot(ts, L, color=col, label=lab); axs[1].plot(ts, V, color=col, label=lab)
    for ax in axs:
        ax.axvline(1 / 3, color="0.5", ls=":", lw=0.8); ax.set_xlabel("$T^*$"); ax.set_xlim(0, 1)
    axs[0].text(0.62, 0.9, "antibiotic from $T^*=1/3$", fontsize=7, color="0.35")
    axs[0].set_ylabel(r"seed $\Sigma_i\bar\phi_i$"); axs[0].set_ylim(0, 1.02)
    axs[1].set_ylabel(r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]"); axs[1].set_ylim(bottom=0)
    axs[1].legend(frameon=False, fontsize=7.5)
    for ax, t in zip(axs, "ab"):
        ax.text(0.03, 0.96, f"({t})", transform=ax.transAxes, va="top")
    o = ROOT / "assets" / "fig_stress_living_over_time.png"
    fig.savefig(o, dpi=600)
    print("wrote", o)


if __name__ == "__main__":
    main()
