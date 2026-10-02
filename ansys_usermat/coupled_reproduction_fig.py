#!/usr/bin/env python3
"""Slide figure: does the coupled scheme give back Klempt et al. 2026?

Per case (Table 1, cases 3 and 6), composition chi_1 = phi_1/(phi_1+phi_2)
over the point model's own time (steps of 1e-4, as in the paper):

  line     the point model on its own, from the paper's start (0.2, 0.2);
  circles  scheme A (rescaled to the field's amount), driven by the amount
           the point model produces itself: the coupling must change
           nothing (tests/test_composition_reduces_to_point_model.py);
  star     ANSYS, scheme B (independent point model at the Gauss point),
           s = 0.15, dt = 0.1, T* = 1, i.e. 1500 point-model steps. Values
           from the IKMHIWI03 runs of 2026-10-01/02, which replay the
           stand-alone reference exactly (judge_comp_trace.py --age);
  diamond  case 3 only: the end value read off the paper's Fig. 3.

    python ansys_usermat/coupled_reproduction_fig.py
        -> assets/fig_coupled_reproduction.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "coupling"))
import composition_reference as cref  # noqa: E402
import ecology_jax as eco  # noqa: E402
import material_server as ms  # noqa: E402

OUT = HERE.parent / "assets" / "fig_coupled_reproduction.png"
DT, NSTEPS, NSUB = 0.01, 15, 100            # coupling step 0.01 = 100 paper steps
ANSYS_B = {"2sp_case3": 0.6137, "2sp_case6": 0.0138}
PAPER = {"2sp_case3": (500, 0.60 / (0.60 + 0.40))}
TITLES = {"2sp_case3": "case 3: coexistence", "2sp_case6": "case 6: species 2 takes over"}
C_LINE, C_A, C_B = "#2a78d6", "#eb6834", "#1f2933"


def curves(case):
    ms.set_case(case)
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    g = cref.rescale(cref.seed(0.5), 0.4)
    t, chi, amount = [0], [0.5], []
    for k in range(NSTEPS * NSUB):              # paper resolution for the line
        if k % NSUB == 0:
            amount.append(g[0] + g[1])
        g, _ = eco.ecology_substeps(g, th, 1.0e-4, 1, 2, hp)
        g = np.asarray(g)
        t.append(k + 1)
        chi.append(g[0] / (g[0] + g[1]))
    coupled = cref.reference(amount, th, hp, DT, chi1=0.5)
    ms.set_case(None)
    ms.set_active_species(5)
    tc = np.arange(1, NSTEPS + 1) * NSUB
    cc = np.array([r[0] / (r[0] + r[1]) for r in coupled])
    return np.array(t), np.array(chi), tc, cc


def main():
    import figstyle
    figstyle.apply()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
    for ax, case in zip(axs, TITLES):
        t, chi, tc, cc = curves(case)
        dev = float(np.max(np.abs(cc - chi[tc])))
        ax.plot(t, chi, color=C_LINE, lw=2.2, label="point model alone (Klempt et al. 2026)")
        ax.plot(tc, cc, "o", ms=7, mfc="none", mec=C_A, mew=2,
                label="coupled A: rescaled to the field's amount (reference)")
        ax.plot(1500, ANSYS_B[case], "*", ms=16, color=C_B,
                label="coupled in ANSYS, B: independent")
        if case in PAPER:
            ax.plot(*PAPER[case], "D", ms=9, mfc="none", mec="#7b8794", mew=2,
                    label="read off the paper's figure")
        ax.text(0.97, 0.05 if case == "2sp_case3" else 0.55,
                f"max |A - alone| = {dev:.1e}\nANSYS B = {ANSYS_B[case]:.4f}, "
                f"alone = {chi[1500]:.4f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
                bbox=dict(fc="white", ec="#cbd2d9", boxstyle="round"))
        ax.set(xlabel=r"point-model time step ($\Delta t = 10^{-4}$)",
               title=TITLES[case], ylim=(-0.03, 1.0))
        print(f"{case}: max |A - alone| {dev:.2e}, alone at 1500 {chi[1500]:.4f}, "
              f"ANSYS B {ANSYS_B[case]}")
    axs[0].set_ylabel(r"composition $\chi_1 = \varphi_1/(\varphi_1+\varphi_2)$")
    axs[1].plot([], [], "D", ms=9, mfc="none", mec="#7b8794", mew=2,
                label="read off the paper's Fig. 3 (case 3)")
    axs[1].legend(fontsize=9, loc="upper right")
    fig.suptitle("Coupled schemes against the point model of Klempt et al. 2026: "
                 "same amount, same composition", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
