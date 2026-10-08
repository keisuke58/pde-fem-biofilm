#!/usr/bin/env python3
"""wp3_figs.py -- the two WP3 figures of the seed-growth runs in the partner's
element with the surface-point fix (RUN_WP3_IKMHIWI03.md, blocks C and D), from
the run JSON files only (no ANSYS needed).

    python ansys_usermat/apdl/wp3_figs.py
        -> assets/fig_wp3_stability.png, assets/fig_wp3_convergence.png

Figure 1, time-step limit of the explicit phi update: every run against
lambda = beta dt / h^2. A run is "stopped" when ANSYS ended before T* = 1.1
(its JSON has no all_stress), "oscillating" when it reached T* = 1.1 but
alpha - 1 is negative somewhere (alpha - 1 >= 0 whenever phi >= 0), otherwise
"smooth". Lines: 1/6 (seven-point Laplacian in 3-D) and 1/2 (h^2 / (2 dt) >= beta,
Rudolf et al. 2025 Eq. 3).
Figure 2, mesh convergence at beta = 0.02 mm^2/T*: seed mean von Mises stress
and alpha - 1 at fixed positions (wp3_alpha_profile.py) against h, for
dt = 0.025 and 0.0125 and extrapolated to dt -> 0 (first order).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
from compare_surface_fix import measures  # noqa: E402
from wp3_alpha_profile import profile  # noqa: E402

RES = HERE / "results"
FIX = RES / "2026-10-wp3_fix"
SURF = RES / "2026-10-07_surface_fix"
ASSETS = HERE.parents[1] / "assets"
L = 2.0   # cube edge, mm


def deck_values(rec, name):
    """(n, dt, beta) from the run name and the deck header line of make_wired_deck.py."""
    n = int(re.search(r"n(\d+)_|^ds(\d+)_", name).group(1) or re.search(r"^ds(\d+)_", name).group(1))
    hd = rec.get("deck_header", "")
    hd = hd if isinstance(hd, str) else " ".join(hd)
    dt = re.search(r"deltim ([\d.eE+-]+)", hd)
    beta = re.search(r"MY_BETA1=([\d.eE+-]+)", hd)
    return n, float(dt.group(1)) if dt else 0.025, float(beta.group(1)) if beta else 0.02


def classify():
    rows = []
    for j in sorted(FIX.glob("wp3[cd]_*.json")) + [SURF / "ds8_beta002_dt4.json", SURF / "ds16_beta002_dt4.json"]:
        rec = json.loads(j.read_text())
        n, dt, beta = deck_values(rec, j.stem)
        lam = beta * dt / (L / n) ** 2
        if "all_stress" not in rec:
            state = "stopped"
        elif min(rec["all_stress"]["alpha"]) < 0:
            state = "oscillating"
        else:
            state = "smooth"
        rows.append((lam, n, state, j.stem))
    return rows


def fig_stability():
    rows = classify()
    fig, ax = plt.subplots(figsize=(6.4, 2.8))
    style = {"smooth": ("o", "C0"), "oscillating": ("X", "C1"), "stopped": ("^", "C3")}
    ypos = {8: 0, 16: 1, 24: 2}
    for state, (mk, col) in style.items():
        pts = [(lam, ypos[n]) for lam, n, s, _ in rows if s == state]
        if pts:
            x, y = zip(*pts)
            ax.scatter(x, y, marker=mk, color=col, s=46, label=state, zorder=3)
    ax.axvline(1 / 6, color="k", lw=0.9, ls="--")
    ax.text(1 / 6 * 0.95, 1.5, r"$1/6$", ha="right", va="center", rotation=90)
    ax.axvline(0.5, color="0.5", lw=0.9, ls=":")
    ax.text(0.5 * 0.95, 1.5, r"$1/2$ (Rudolf et al. 2025)", ha="right", va="center", rotation=90,
            color="0.4", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(2e-3, 0.8)
    ax.set_ylim(-0.5, 2.9)
    ax.set_yticks([0, 1, 2], [r"$8^3$", r"$16^3$", r"$24^3$"])
    ax.set_xlabel(r"$\lambda=\beta\,\Delta t/h^2$")
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    fig.tight_layout()
    out = ASSETS / "fig_wp3_stability.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out, rows


def fig_convergence():
    runs = {8: (SURF / "ds8_beta002_dt4", FIX / "wp3c_n8_dt00125"),
            16: (SURF / "ds16_beta002_dt4", FIX / "wp3c_n16_dt00125"),
            24: (FIX / "wp3c_n24_dt0025", FIX / "wp3c_n24_dt00125")}
    h = np.array([L / n for n in runs])
    q = {"vM": [], "seed mean": [], "centre": [], "edge": []}
    for n, pair in runs.items():
        vals = []
        for p in pair:
            rec = json.loads(p.with_suffix(".json").read_text())
            m, _ = measures(rec)
            pr = profile(rec)
            vals.append((m["seed vM"], pr["seed mean"], pr["line"][0], pr["line"][4]))
        vals = np.array(vals)
        ext = 2 * vals[1] - vals[0]
        for k, key in enumerate(q):
            q[key].append((vals[0, k], vals[1, k], ext[k]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.8))
    v = np.array(q["vM"])
    for col, lab, mk in ((0, r"$\Delta t=0.025$", "o"), (1, r"$\Delta t=0.0125$", "s"), (2, r"$\Delta t\to 0$", "D")):
        a1.plot(h, v[:, col] * 1e4, marker=mk, label=lab)
    a1.set_xlabel(r"$h$ (mm)")
    a1.set_ylabel(r"seed mean $\sigma_{vM}$ ($10^{-4}$ Pa)")
    a1.legend(frameon=False, fontsize=8)
    for key, lab, mk in (("seed mean", "seed mean", "o"), ("centre", "centre", "s"), ("edge", "0.4 mm from centre", "^")):
        e = np.array(q[key])[:, 2]
        a2.plot(h, e / e[-1], marker=mk, label=lab)
    a2.axhline(1, color="0.6", lw=0.8)
    a2.set_xlabel(r"$h$ (mm)")
    a2.set_ylabel(r"$\alpha-1$ / value on $24^3$ ($\Delta t\to 0$)")
    a2.legend(frameon=False, fontsize=8)
    for a in (a1, a2):
        a.invert_xaxis()
    fig.tight_layout()
    out = ASSETS / "fig_wp3_convergence.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main():
    figstyle.apply()
    out1, rows = fig_stability()
    for lam, n, s, name in sorted(rows):
        print(f"{lam:.3f}  {n:>2}^3  {s:11s} {name}")
    print(out1)
    print(fig_convergence())


if __name__ == "__main__":
    main()
