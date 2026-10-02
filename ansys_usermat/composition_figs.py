#!/usr/bin/env python3
"""Composition figures without ANSYS: the two composition schemes over time,
and the sensitivity to the point-model clock s.

Both come from the stand-alone references the fragments are checked against
(composition_reference.reference / reference_age), for a prescribed amount
phi_3D(t) as the 3D field (Klempt 2024) produces it in a gradient-free
region. Two-species cases 3 and 6 of Klempt et al. 2026, Table 1.

  scheme A (rescaled):    the point model is rescaled every substep so that
                          phi_1 + phi_2 = min(phi_3D, phi_cap); it decides
                          only the composition chi_1 = phi_1 / (phi_1 + phi_2).
  scheme B (independent): the point model starts from (0.2, 0.2) when phi_3D
                          first reaches phi_min and then runs on its own.

Assumptions (not from a paper): s = 0.15, phi_cap = 0.9, phi_min = 0.01.

    python ansys_usermat/composition_figs.py   -> assets/fig_composition_*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "coupling"))
import composition_reference as cref  # noqa: E402
import material_server as ms  # noqa: E402

OUT = HERE.parent / "assets"
DT, T, S, PHI_CAP, PHI_MIN = 0.1, 1.0, 0.15, 0.9, 1.0e-2
CASES = {"2sp_case3": "case 3", "2sp_case6": "case 6"}
SCENARIOS = {
    r"seed, $\phi_{3D} = 1$": lambda t: 1.0,
    r"interior, $\phi_{3D} = 0.4$": lambda t: 0.4,
    r"front arriving at $T^* = 0.3$": lambda t: float(np.clip((t - 0.3) / 0.5, 0.0, 1.0)),
}
import figstyle  # noqa: E402
figstyle.apply()


def chi1(states):
    """chi_1 per substep; nan where the point model has not started."""
    out = []
    for g in states:
        tot = g[0] + g[1]
        out.append(g[0] / tot if g[6] > 0.0 and tot > 0.0 else np.nan)
    return np.array(out)


def run(case, f, s, scheme):
    ms.set_case(case)
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    n = int(round(T / DT))
    series = [f(k * DT) for k in range(1, n + 1)]
    if scheme == "A":
        st = cref.reference(series, th, hp, DT, phi_min=PHI_MIN, s=s, phi_cap=PHI_CAP)
        # before phi_3D reaches phi_min the seed is only carried: not plotted
        st = [g if p >= PHI_MIN else np.zeros(12) for g, p in zip(st, series)]
    else:
        st = cref.reference_age(series, th, hp, DT, phi_min=PHI_MIN, s=s)
    ms.set_case(None)
    return np.arange(1, n + 1) * DT, chi1(st)


def fig_schemes():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.6), sharex=True)
    rows = []
    for i, (case, cname) in enumerate(CASES.items()):
        for j, (sname, f) in enumerate(SCENARIOS.items()):
            a = ax[i, j]
            for scheme, lab, st in (("A", "A: rescaled to $\\phi_{3D}$", "o-"),
                                    ("B", "B: independent from arrival", "s--")):
                t, c = run(case, f, S, scheme)
                a.plot(t, c, st, ms=4, label=lab)
                rows.append((cname, sname, scheme, c[-1]))
            tt = np.linspace(0, T, 101)
            a2 = a.twinx()
            a2.plot(tt, [f(x) for x in tt], color="0.6", lw=1)
            a2.set_ylim(0, 1.05)
            a2.grid(False)
            a2.spines["right"].set_visible(True)
            if j == 2:
                a2.set_ylabel(r"$\phi_{3D}$ (grey)", color="0.4")
            else:
                a2.set_yticklabels([])
            a.set_ylim(-0.02, 1.0)
            a.set_title(f"{cname}, {sname}", fontsize=10)
            if j == 0:
                a.set_ylabel(r"share of species 1, $\phi_1/(\phi_1+\phi_2)$")
            if i == 1:
                a.set_xlabel(r"time $T^*$")
    ax[0, 0].legend(fontsize=9, loc="lower left")
    fig.suptitle(rf"Composition of the two schemes; Klempt et al. 2026 cases; "
                 rf"assumed: $s$ = {S}, $\phi_{{cap}}$ = {PHI_CAP}, $\phi_{{min}}$ = {PHI_MIN}, "
                 rf"$\Delta t$ = {DT}", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig_composition_schemes.png")
    plt.close(fig)
    return rows


def fig_s():
    s_vals = [0.05, 0.1, 0.15, 0.25, 0.5]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    rows = []
    for a, (case, cname) in zip(ax, CASES.items()):
        for k, (sname, f) in enumerate(SCENARIOS.items()):
            y = [run(case, f, s, "A")[1][-1] for s in s_vals]
            a.plot(s_vals, y, "o-", color=f"C{k}", label=f"A, {sname}")
            rows += [(cname, sname, "A", s, v) for s, v in zip(s_vals, y)]
        y = [run(case, SCENARIOS[r"seed, $\phi_{3D} = 1$"], s, "B")[1][-1] for s in s_vals]
        a.plot(s_vals, y, "s--", color="0.3", label="B, seed")
        rows += [(cname, "seed", "B", s, v) for s, v in zip(s_vals, y)]
        a.axvline(S, color="0.6", lw=1, ls=":")
        a.set(xscale="log", xlabel=r"point-model clock $s$ (assumed)",
              ylabel=r"$\phi_1/(\phi_1+\phi_2)$ at $T^* = 1$", title=cname)
        a.set_xticks(s_vals)
        a.set_xticklabels([f"{s:g}" for s in s_vals])
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].legend(fontsize=8)
    fig.suptitle(r"Sensitivity of the composition to $s$ (dotted: $s$ = 0.15 used in the runs); "
                 r"growth and stress do not depend on $s$", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig_composition_s.png")
    plt.close(fig)
    return rows


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for r in fig_schemes():
        print("schemes  %-7s %-34s %s  chi1(T*=1) = %.4f" % r)
    for r in fig_s():
        print("s-sweep  %-7s %-34s %s  s = %-5g chi1 = %.4f" % r)
