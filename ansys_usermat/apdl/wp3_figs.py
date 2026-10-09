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
    if beta:
        b = float(beta.group(1))
    else:
        # a deck derived from a beta deck (wp3e_n8_b01_dt*: only deltim in its header):
        # the run name carries beta as _b01 = 0.1, _b005 = 0.05, _b001 = 0.01
        m = re.search(r"_b(\d+)_", name + "_")
        b = float("0." + m.group(1)[1:]) if m else 0.02
    return n, float(dt.group(1)) if dt else 0.025, b


# runs that ended for a reason other than stability: the 24^3 beta 0.05 run hung for 17 h after
# ANSYS switched from PCG to the sparse solver at T* = 0.025 and was stopped by hand (lambda 0.09)
NOT_STABILITY = {"wp3d_n24_b005_dt00125"}


def classify():
    rows = []
    for j in sorted(FIX.glob("wp3[cde]_*.json")) + [SURF / "ds8_beta002_dt4.json", SURF / "ds16_beta002_dt4.json"]:
        if j.stem in NOT_STABILITY:
            continue
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
    # runs to T* = 5 (name ends in _T5) sit slightly above their mesh row with a black edge:
    # at T* = 1.1 an unstable mode can stay invisible above the limit (8^3, Appendix D)
    for state, (mk, col) in style.items():
        for long_run in (False, True):
            pts = [(lam, ypos[n] + (0.22 if long_run else 0.0)) for lam, n, s, name in rows
                   if s == state and name.endswith("_T5") == long_run]
            if pts:
                x, y = zip(*pts)
                ax.scatter(x, y, marker=mk, color=col, s=46, zorder=3,
                           edgecolors="k" if long_run else "none", linewidths=1.0,
                           label=None if long_run else state)
    if any(name.endswith("_T5") for *_, name in rows):
        ax.scatter([], [], marker="o", facecolors="none", edgecolors="k", s=46, label=r"run to $T^*=5$")
    # predicted limit of the NEM stencil on the integration-point lattice
    # (nem_stability.py: 30 or 32 neighbours, beta* = 0.2 mm)
    pred = FIX / "nem_stability.json"
    if pred.exists():
        cases = json.loads(pred.read_text())["cases"]
        first = True
        for nmesh, y in ypos.items():
            lams = [c["lam_max"] for c in cases
                    if c["mesh"] == nmesh and c["beta_star"] == 0.2 and c["nneigh"] in (30, 32)]
            if lams:
                lo, hi = min(lams) * 0.985, max(lams) * 1.015
                ax.fill_betweenx([y - 0.32, y + 0.32], lo, hi, color="k", alpha=0.18, lw=0,
                                 zorder=1, label="predicted limit" if first else None)
                first = False
    ax.axvline(1 / 6, color="k", lw=0.9, ls="--")
    ax.text(1 / 6 * 0.95, 1.5, r"$1/6$", ha="right", va="center", rotation=90)
    ax.axvline(0.5, color="0.5", lw=0.9, ls=":")
    ax.text(0.5 * 0.95, 1.5, r"$1/2$ (Rudolf et al. 2025)", ha="right", va="center", rotation=90,
            color="0.4", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(2e-3, 0.8)
    ax.set_ylim(-0.5, 3.55)
    ax.set_yticks([0, 1, 2], [r"$8^3$", r"$16^3$", r"$24^3$"])
    ax.set_xlabel(r"$\lambda=\beta\,\Delta t/h^2$")
    ax.legend(loc="upper left", frameon=False, fontsize=8, ncol=5, columnspacing=1.0, handletextpad=0.3)
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
    # right: error of alpha - 1 against the reference solution (wp3_reference.py),
    # averaged over the same elements, after extrapolation to dt -> 0
    ref = json.loads((FIX / "wp3_reference.json").read_text())["nem_error"]
    hh = np.array([L / int(n) for n in sorted(ref, key=int)])
    for key, lab, mk in (("err seed mean", "seed mean", "o"), ("err L2", r"$L_2$ norm, cube", "D"),
                         ("err line", "centre", "s"), ("err line4", "0.4 mm from centre", "^")):
        e = []
        for n in sorted(ref, key=int):
            row = ref[n]["dt->0"]
            e.append(row["err line"][0] if key == "err line" else row["err line"][4] if key == "err line4" else row[key])
        a2.loglog(hh, np.abs(e) * 100, marker=mk, label=lab)
    e0 = abs(ref[min(ref, key=int)]["dt->0"]["err seed mean"]) * 100
    a2.loglog(hh, e0 * (hh / hh[0]) ** 2, color="0.6", lw=0.8, ls="--")
    a2.text(hh[-1] * 1.03, e0 * (hh[-1] / hh[0]) ** 2, r"$h^2$", fontsize=8, color="0.4", va="center")
    a2.set_xlabel(r"$h$ (mm)")
    a2.set_ylabel(r"error of $\alpha-1$ (%), $\Delta t\to 0$")
    a2.set_xticks(hh)
    a2.set_xticklabels([f"{v:.3g}" for v in hh])
    a2.set_ylim(0.1, 20)
    a2.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    a2.set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"])
    a2.minorticks_off()
    a2.legend(frameon=False, fontsize=8, loc="upper right")
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
