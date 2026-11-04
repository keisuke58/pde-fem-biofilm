"""fig_four_species_time_history.py -- seed means over time for the four-species case of Klempt et al. 2026 (case 4) in the
coupled cube examples, with antibiotic from T* = 1/3 (64^3, the --snap-stress runs of the videos): (a) living volume fraction
sum_i phibar_i, (b) nutrient c, (c) growth alpha - 1, (d) von Mises stress; separate start solid, mixed start dashed.

    python ansys_usermat/apdl/fig_four_species_time_history.py   -> assets/fig_four_species_time_history.png
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
RUNS = [("separate start", "separate_ab", "", "-"), ("mixed start", "mixed_ab", "_mixeq", "--")]
COL = "#e8590c"


def main():
    figstyle.apply()
    fig, axs = plt.subplots(1, 4, figsize=(6.3, 1.9), constrained_layout=True)
    for lab, case, mix, ls in RUNS:
        z = np.load(ROOT / "figdata" / (f"coupled_living_{case}_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv{mix}_a100"
                                        f"_gab1_gc1_tab0.333333_snapsvm_transc1e+08_k26c4_tvd_wgmirror_fast_{case}.npz"))
        ts = sorted({float(k.split("_", 1)[1]) for k in z.files if k.startswith("vm16_")})
        s64, s16 = W.seed_mask(64, C8), W.seed_mask(z[f"vm16_{ts[0]}"].shape[0], C8)
        key = lambda n, t: f"{n}_{t}" if f"{n}_{t}" in z.files else f"{n}_{t:g}"
        m = lambda n, t, s=s64, f=1.0: f * float(np.nan_to_num(z[key(n, t)])[s].mean())
        axs[0].plot(ts, [m("L", t) for t in ts], ls, color=COL, label=lab)
        axs[1].plot(ts, [m("c", t) for t in ts], ls, color=COL)
        axs[2].plot(ts, [1e4 * (m("alpha", t) - 1.0) for t in ts], ls, color=COL)
        axs[3].plot(ts, [m("vm16", t, s16, 1e4) for t in ts], ls, color=COL)
    labels = [r"seed $\Sigma_i\bar\phi_i$", r"seed $c$", r"seed $\alpha-1$ [$10^{-4}$]", r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]"]
    for ax, t, yl in zip(axs, "abcd", labels):
        ax.axvline(1 / 3, color="0.5", ls=":", lw=0.8); ax.set_xlabel("$T^*$"); ax.set_xlim(0, 1)
        ax.set_ylabel(yl, fontsize=8); ax.tick_params(labelsize=7.5); ax.set_ylim(bottom=0)
        ax.text(0.04, 0.96, f"({t})", transform=ax.transAxes, va="top", fontsize=8.5)
    axs[0].set_ylim(0, 1.02); axs[1].set_ylim(0, 1.02)
    axs[0].legend(frameon=False, fontsize=6.5, loc="lower right")
    o = ROOT / "assets" / "fig_four_species_time_history.png"
    fig.savefig(o, dpi=600)
    print("wrote", o)


if __name__ == "__main__":
    main()
