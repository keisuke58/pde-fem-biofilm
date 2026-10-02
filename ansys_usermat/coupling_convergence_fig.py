#!/usr/bin/env python3
"""Convergence of the coupled composition in the coupling step.

The field hands the point model its amount phi_3D once per deck increment,
so the point model sees phi_3D as piecewise constant in time. Refining the
increment dt must make the coupled composition converge, to first order, to
the one obtained when phi_3D is exchanged at every step of the point model
(dt_pm = 1e-4, the paper's step). Inside the point model the step is 1e-4 in
every run, so only the coupling step changes.

Scheme A (rescaled), prescribed smooth amount
phi_3D(t) = 0.05 + 0.425 (1 - cos(pi t)), T* in [0, 1], s = 0.15, cases 3, 6.

Result (2026-10-02): first order in the largest error over time; at the
deck step 0.1 of the ANSYS runs the transient composition is off by up to
0.07 (case 3) and 0.29 (case 6), the value at T* = 1 by 3e-3 / 9e-3. In the
seeded element of those runs phi_3D is constant (capped), so this error does
not arise there; where the amount changes in time (a front), a deck step of
0.01 or below is needed for the transient.

    python ansys_usermat/coupling_convergence_fig.py
        -> assets/fig_coupling_convergence.png
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
import figstyle  # noqa: E402
import material_server as ms  # noqa: E402

OUT = HERE.parent / "assets" / "fig_coupling_convergence.png"
S = 0.15
# dt such that s*dt is a whole number of 1e-4 steps (inner step fixed at 1e-4)
N_INNER = [150, 75, 30, 15, 6, 3]                 # inner steps per increment
DT = [n * 1.0e-4 / S for n in N_INNER]            # 0.1 ... 0.002
DT_REF = 1.0e-4 / S                               # exchange at every inner step
CASES = {"2sp_case3": ("case 3", "#2a78d6"), "2sp_case6": ("case 6", "#eb6834")}


def phi3(t):
    return 0.05 + 0.425 * (1.0 - np.cos(np.pi * t))


def chi_at_tenths(case, dt):
    """chi_1 at T* = 0.1, 0.2, ..., 1."""
    ms.set_case(case)
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    n = int(round(1.0 / dt))
    series = [phi3((k + 1) * dt) for k in range(n)]
    st = cref.reference(series, th, hp, dt, s=S)
    ms.set_case(None)
    c = np.array([g[0] / (g[0] + g[1]) for g in st])
    return c[n // 10 - 1::n // 10]


def main():
    figstyle.apply(12)
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.3))
    t = np.linspace(0, 1, 200)
    ax[0].plot(t, phi3(t), color="#1f2933", lw=2, label=r"prescribed $\phi_{3D}(t)$")
    for d, c in ((0.1, "#eb6834"), (0.02, "#2a78d6")):
        tk = np.arange(0, 1 + 1e-9, d)
        ax[0].step(tk, phi3(np.minimum(tk + d, 1)), where="post", color=c, lw=1.2,
                   label=rf"seen by the point model, $\Delta t$ = {d:g}")
    ax[0].set(xlabel=r"time $T^*$", ylabel=r"amount $\phi_{3D}$",
              title="prescribed amount and what the point model sees")
    ax[0].legend(loc="upper left", fontsize=10)

    for case, (lab, col) in CASES.items():
        ref = chi_at_tenths(case, DT_REF)
        runs = [chi_at_tenths(case, d) - ref for d in DT]
        emax = np.array([np.abs(r).max() for r in runs])
        eend = np.array([abs(r[-1]) for r in runs])
        rate = np.polyfit(np.log(DT), np.log(emax), 1)[0]
        ax[1].loglog(DT, emax, "o-", color=col, lw=2, ms=7,
                     label=f"{lab}, largest over $T^*$: order {rate:.2f}")
        ax[1].loglog(DT, eend, "s--", color=col, lw=1.2, ms=5, alpha=0.8,
                     label=f"{lab}, at $T^*$ = 1")
        print(f"{case}: dt {[round(d, 4) for d in DT]}; max {[f'{e:.2e}' for e in emax]}; "
              f"end {[f'{e:.2e}' for e in eend]}; order {rate:.2f}")
    x = np.array([DT[-1], DT[0]])
    ax[1].loglog(x, 2.0 * x, "k:", lw=1.2, label="first order")
    ax[1].axvline(0.1, color="#9aa5b1", ls=":", lw=1)
    ax[1].text(0.093, 2e-6, "deck step of\nthe ANSYS runs", color="#616e7c", fontsize=10,
               ha="right", va="bottom")
    ax[1].set(xlabel=r"coupling step $\Delta t$ (deck increment)",
              ylabel=r"error in $\phi_1/(\phi_1+\phi_2)$",
              title="composition converges in the coupling step")
    ax[1].legend(fontsize=9.5, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.suptitle(r"Coupled composition (rescaled scheme), point model at its own step $10^{-4}$; "
                 r"reference: amount exchanged at every point-model step", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
