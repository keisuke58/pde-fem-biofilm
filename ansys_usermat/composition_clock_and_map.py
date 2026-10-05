#!/usr/bin/env python3
"""Two views of the two-species composition (case 6 of Klempt et al. 2026) in
the coupled scheme of the ANSYS runs (amount from the field, phi_cap = 0.9,
coupling step 0.025, psi = 0.999 at the start):

A  the nutrient as the point model's clock: the share phi_1/(phi_1+phi_2) at
   the end, plotted against the effective time s T* c_rel. Points: every seed
   element of the 8^3 ANSYS runs of the week set (case 6, Klempt 2024
   stiffness, beta = 0.02; c_rel = the element's nutrient at the end, s and
   T* from the run name). Lines: a packed point in Python at constant c_rel.
B  which species wins: the share at T* = 1 against the local nutrient c_rel
   and the start share chi_0, for s = 0.15 and 0.25 (Python, packed point).

    python ansys_usermat/composition_clock_and_map.py
        -> assets/fig_composition_clock.png, assets/fig_composition_map.png

Result (6 Oct 2026, 8^3 week runs up to 50 of 148):
A  All seed elements of all runs (consumption 0 to 6, s 0.05 to 1, T* 1 and 2)
   fall on one curve against s T* c_rel: the share stays near the start value
   0.5 up to s T* c_rel ~ 0.03 and drops to below 0.1 by ~ 0.08, as the packed
   point in Python (lines). s, T* and the local nutrient enter only through
   this product; the nutrient sets the local clock of the point model. After
   the takeover the level scatters between 0.02 and 0.12 (the amount, not the
   clock, then matters; T* = 2 sits higher).
B  The winner is set by the start share alone, with the boundary at 0.51-0.52
   for every nutrient level; the nutrient only sets how fast the outcome is
   reached (below c_rel ~ 0.3 it is not reached by T* = 1 at s = 0.15).
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import composition_transport_check as T  # noqa: E402
import figstyle  # noqa: E402

ROOT = HERE.parent
WEEK = HERE / "apdl" / "results" / "2026-10-05_week"
DT, CAP = 0.025, 0.9


def packed(chi0, crel, s, t_end, history=False):
    """share at t_end of a packed point (amount CAP), arrays chi0 / crel;
    history=True: the share after every coupling step, shape (steps, n)."""
    eco = T.eco
    T.ms.set_case("2sp_case6")
    th, hp = T.ms.ECOLOGY_CASE["theta"], dict(T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    dt_pm = s * DT
    n_sub = max(1, math.ceil(dt_pm / 1e-4 - 1e-9))
    tht = jnp.asarray(th, dtype=jnp.float64)
    eta = jnp.asarray(hp["eta"], dtype=jnp.float64)
    run = jax.jit(jax.vmap(lambda g, c: eco._substep_scan(
        g, tht, dt_pm / n_sub, jnp.arange(n_sub), eco.active_mask(2), c,
        jnp.float64(hp["alpha"]), eta)[0]))
    chi0 = np.asarray(chi0, float).ravel(); crel = np.asarray(crel, float).ravel()
    G = np.zeros((chi0.size, 12)); G[:, 0] = chi0; G[:, 1] = 1 - chi0; G[:, 6:8] = 0.999
    p = np.full(chi0.size, CAP)
    cst = jnp.asarray(hp["c"] * crel)
    hist = []
    for _ in range(int(round(t_end / DT))):
        G = T.rescale(np.asarray(run(jnp.asarray(T.rescale(G, p)), cst)), p)
        hist.append(G[:, 0] / (G[:, 0] + G[:, 1]))
    return np.array(hist) if history else hist[-1]


def ansys_points():
    out = []
    for f in sorted(WEEK.glob("w8_c6_g*_s*.json")):
        m = re.fullmatch(r"w8_c6_g(\d+)_s(\d+)(_T(\d+))?", f.stem)
        if not m:
            continue                                  # variants (gw, cap, ...) left out
        g = int(m.group(1))
        s = float("0." + m.group(2)[1:]) if m.group(2).startswith("0") else float(m.group(2))
        t_end = float(m.group(4)) if m.group(4) else 1.0
        d = json.load(open(f))
        nf = d["nut_field"]; seed = set(d["seed_BIOFILM1"])
        for e, c, p1, p2 in zip(nf["elem"], nf["nut1"], nf["phi1"], nf["phi2"]):
            if int(e) in seed:
                out.append((g, s, t_end, min(max(c, 0.0), 1.0), p1 / (p1 + p2)))
    return np.array(out)


def fig_clock():
    pts = ansys_points()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for i, g in enumerate(sorted(set(pts[:, 0]))):
        q = pts[pts[:, 0] == g]
        for t_end, mk in ((1.0, "o"), (2.0, "^")):
            r = q[q[:, 2] == t_end]
            if len(r):
                ax.plot(r[:, 1] * r[:, 2] * r[:, 3], r[:, 4], mk, ms=3.5, color=f"C{i}", alpha=0.7,
                        label=f"consumption {g:g}" + (", $T^*=2$" if t_end == 2 else ""))
    t = DT * np.arange(1, int(round(2.0 / DT)) + 1)
    for crel, ls in ((1.0, "-"), (0.4, "--")):
        h = packed([0.5], [crel], 0.25, 2.0, history=True)[:, 0]
        ax.plot(0.25 * t * crel, h, ls, color="k", lw=1.2,
                label=f"packed point, Python, $c_{{rel}}={crel:g}$")
    ax.set_xlabel(r"effective time $s\,T^*\,c_{\rm rel}$")
    ax.set_ylabel(r"$\phi_1/(\phi_1+\phi_2)$ at the end")
    ax.set_xscale("log"); ax.legend(fontsize=7.5, frameon=False, ncol=2)
    out = ROOT / "assets" / "fig_composition_clock.png"
    fig.tight_layout(); fig.savefig(out, dpi=200); print("wrote", out)
    return pts


def fig_map():
    c = np.linspace(0.05, 1.0, 20); chi = np.linspace(0.02, 0.98, 25)
    C, X = np.meshgrid(c, chi)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for a, s in zip(ax, (0.15, 0.25)):
        sh = packed(X, C, s, 1.0).reshape(X.shape)
        im = a.pcolormesh(c, chi, sh, vmin=0, vmax=1, cmap="RdBu_r", shading="nearest")
        a.contour(c, chi, sh, [0.5], colors="k", linewidths=1.2)
        a.set_title(f"$s={s:g}$, $T^*=1$", fontsize=11); a.set_xlabel(r"local nutrient $c_{\rm rel}$")
        print(f"s={s}: share range {sh.min():.3f}-{sh.max():.3f}")
    ax[0].set_ylabel(r"start share $\phi_1/(\phi_1+\phi_2)$")
    fig.colorbar(im, ax=ax, label=r"$\phi_1/(\phi_1+\phi_2)$ at $T^*=1$")
    out = ROOT / "assets" / "fig_composition_map.png"
    fig.savefig(out, dpi=200, bbox_inches="tight"); print("wrote", out)


if __name__ == "__main__":
    figstyle.apply(size=10)
    fig_clock()
    fig_map()
