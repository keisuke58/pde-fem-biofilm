"""fig_growth_living_over_time.py -- growth from the living biomass against Eq. 36 of Klempt 2024, seed means over time.

Keio continuation, not the thesis model. Eq. 36, d alpha/dt = k_alpha phi_3D, against the law in which only living biomass
grows, d alpha/dt = k_alpha phi_3D sum_i phibar_i (--growth-living-energy: the coupling energy -k_alpha alpha phi of Eq. 20
multiplied by sum_i phibar_i, so the source of Eq. 34 changes with it). Two species, antibiotic from T* = 1/3, 64^3, superbee
off (minmod), dc/dt, outer step 0.025, the --snap-stress runs. Rows: living volume fraction, alpha - 1, von Mises stress
(seed means); columns: separate and mixed (equal shares) start.

    python ansys_usermat/apdl/fig_growth_living_over_time.py -> assets/fig_growth_living_over_time.pdf/.png
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
CASES = [("separate start", "separate_ab", ""), ("mixed start", "mixed_ab", "_mixeq")]
LAWS = [(r"Eq. 36: $\dot\alpha=k_\alpha\phi_{3D}$", "", "0.45", "--"),
        (r"living: $\dot\alpha=k_\alpha\phi_{3D}\,\Sigma_i\bar\phi_i$", "_growlivE", "#c2255c", "-")]


def series(case, mix, law):
    z = np.load(ROOT / "figdata" / (f"coupled_living_{case}_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv{law}{mix}_a100"
                                    f"_gab1_gc1_tab0.333333_snapsvm_transc1e+08_tvd_wgmirror_fast_{case}.npz"))
    ts = sorted({float(k.split("_", 1)[1]) for k in z.files if k.startswith("vm16_")})
    key = lambda n, t: f"{n}_{t}" if f"{n}_{t}" in z.files else f"{n}_{t:g}"
    n = z[key("L", ts[0])].shape
    s_f = W.seed_mask(64, C8)
    if n[0] != 64:                                   # --fast stores the computed quarter: x, z >= 0
        s_f = s_f[32:, :, 32:]
    s16 = W.seed_mask(z[key("vm16", ts[0])].shape[0], C8)
    L = [float(np.nan_to_num(z[key("L", t)])[s_f].mean()) for t in ts]
    A = [1e4 * float((z[key("alpha", t)][s_f] - 1).mean()) for t in ts]
    V = [1e4 * float(np.nan_to_num(z[key("vm16", t)])[s16].mean()) for t in ts]
    return ts, L, A, V


def main():
    figstyle.apply()
    fig, axs = plt.subplots(3, 2, figsize=(6.3, 4.6), sharex=True, constrained_layout=True)
    for j, (title, case, mix) in enumerate(CASES):
        for lab, law, col, ls in LAWS:
            ts, L, A, V = series(case, mix, law)
            for i, y in enumerate((L, A, V)):
                axs[i, j].plot(ts, y, color=col, ls=ls, lw=1.3, label=lab)
            print(case, law or "eq36", f"T*=1: L={L[-1]:.3f} am1={A[-1]:.3f}e-4 vm={V[-1]:.3f}e-4 Pa")
        axs[0, j].set_title(title, fontsize=8.5, loc="center")
        axs[2, j].set_xlabel("$T^*$")
    for ax, t in zip(axs.T.ravel(), "adbecf"[:0] or ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]):
        ax.axvline(1 / 3, color="0.5", ls=":", lw=0.8); ax.set_xlim(0, 1); ax.set_ylim(bottom=0)
        ax.text(0.03, 0.95, t, transform=ax.transAxes, va="top", fontsize=8)
    axs[0, 0].set_ylabel(r"seed $\Sigma_i\bar\phi_i$"); axs[1, 0].set_ylabel(r"seed $\alpha-1$ [$10^{-4}$]")
    axs[2, 0].set_ylabel(r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa]")
    axs[0, 1].text(0.36, 0.86, "antibiotic from $T^*=1/3$", fontsize=7, color="0.35", transform=axs[0, 1].transData)
    axs[1, 0].legend(frameon=False, fontsize=7, loc="lower right")
    A_ = ROOT / "assets"
    fig.savefig(A_ / "fig_growth_living_over_time.pdf"); fig.savefig(A_ / "fig_growth_living_over_time.png", dpi=400)
    print("wrote", A_ / "fig_growth_living_over_time.png")


if __name__ == "__main__":
    main()
