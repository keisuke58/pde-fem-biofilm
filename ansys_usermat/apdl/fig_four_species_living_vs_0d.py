"""fig_four_species_living_vs_0d.py -- the living volume fraction phibar_i = phi_i psi_i of the four species of Klempt et al.
2026, case 4, over time: the point model alone (0D, c* = 100, alpha* = 100 from one third of the run, as in the paper)
against the coupled cube example (64^3, mean over the biofilm and over the seed; histories phibar_sp_* of the JSON).
The 0D run uses the same Fortran point model and the same time link s = 0.15 as the coupled runs.

    python ansys_usermat/apdl/fig_four_species_living_vs_0d.py   -> assets/fig_four_species_living_vs_0d.png
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(ROOT / "ansys_usermat/coupling"))
import figstyle  # noqa: E402
import ecology_native_ctypes as FN  # noqa: E402

COL = ["#3b5bdb", "#74c0fc", "#c92a2a", "#f59f00"]
RES = HERE / "results/2026-10-3d-case3"
RUNS = [("separate start", RES / ("coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1"
                                  "_tab0.333333_transc1e+08_k26c4_tvd_wgmirror_fast_separate_ab.json")),
        ("mixed start", RES / ("coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1"
                               "_tab0.333333_transc1e+08_k26c4_tvd_wgmirror_fast_mixed_ab.json"))]


def zero_d():
    """the point model alone: start 0.02 each, psi 0.999, c* = 100 throughout, alpha* = 100 from T* = 1/3"""
    FN.init(ROOT / "abaqus_composition" / "eco_klempt2026_t2c4_a100.txt")
    g = np.zeros((1, 12)); g[0, :4] = 0.02; g[0, 5] = 0.92; g[0, 6:11] = 0.999
    dt, s, nsub = 0.025, 0.15, 38
    t, out = [0.0], [g[0, :4] * g[0, 6:10]]
    for k in range(40):
        a0 = 1.0 if k * dt >= 1 / 3 - 1e-9 else 0.0
        a1 = 1.0 if (k + 1) * dt > 1 / 3 + 1e-9 else 0.0
        g, ok = FN.step_lin(g, s * dt, nsub, np.ones(1), np.ones(1), np.array([a0]), np.array([a1]))
        assert ok.all()
        t.append((k + 1) * dt); out.append(g[0, :4] * g[0, 6:10])
    return np.array(t), np.array(out)


def main():
    figstyle.apply()
    t0, l0 = zero_d()
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3), sharey=True, constrained_layout=True)
    for i in range(4):
        axs[0].plot(t0, l0[:, i], color=COL[i], lw=1.5, label=f"species {i + 1}")
    axs[0].set_title("point model alone (0D)", fontsize=8.5)
    for ax, (lab, f) in zip(axs[1:], RUNS):
        h = list(json.loads(f.read_text())["runs"].values())[0]
        t = np.array(h["t"]); P = np.array([[np.nan if v is None else v for v in r] for r in h["phibar_sp_bio"]])
        Ps = np.array([[np.nan if v is None else v for v in r] for r in h["phibar_sp_seed"]])
        for i in range(4):
            ax.plot(t[-len(P):], P[:, i], color=COL[i], lw=1.5)
            ax.plot(t[-len(Ps):], Ps[:, i], color=COL[i], lw=0.9, ls="--")
        ax.set_title(f"coupled cube, {lab}", fontsize=8.5)
    for ax in axs:
        ax.axvline(1 / 3, color="0.5", ls=":", lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("$T^*$")
    axs[0].set_ylabel(r"$\bar\phi_i=\phi_i\psi_i$")
    axs[0].legend(frameon=False, fontsize=7)
    axs[2].plot([], [], "k-", lw=1.5, label="mean over the biofilm"); axs[2].plot([], [], "k--", lw=0.9, label="mean over the seed")
    axs[2].legend(frameon=False, fontsize=7, loc="upper right")
    o = ROOT / "assets" / "fig_four_species_living_vs_0d.png"
    fig.savefig(o, dpi=600)
    print("wrote", o, "| 0D end", np.round(l0[-1], 3))


if __name__ == "__main__":
    main()
