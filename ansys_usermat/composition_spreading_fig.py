#!/usr/bin/env python3
"""Composition along a spreading front, in Python: the result ANSYS cannot
show yet (the partner element's front term does not spread the seed).

The growth field is Klempt et al. 2024 test case 4.1 as reproduced in
JAXFEM/klempt2024_timescale.py (Table 2, edge-line nutrient source, time
scale 1.5; growth on both faces and consumption g phi c, the two departures
from Eq. 34/35 as printed that the paper's curves need). On it, the point
model of Klempt et al. 2026 (cases 3 and 6) runs at every grid node with the
same scheme as the coupled ANSYS run: amount from the field, composition from
the point model, rescaled each coupling step (0.1), s = 0.15, phi_cap = 0.9,
phi_min = 0.01; a newly reached node starts from the seed composition, as in
ANSYS. Mid-plane z = 10 um; nutrient on the edge x = y = 20 um.

    python ansys_usermat/composition_spreading_fig.py
        -> assets/fig_composition_spreading.png
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
import composition_transport_check as T  # noqa: E402
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_composition_spreading.png"
TIMES = [1, 3, 6, 10]                       # coupling steps: T* = 0.1, 0.3, 0.6, 1.0


def run(snaps, case):
    """chi_1 = phi_1/(phi_1+phi_2) on the grid after each coupling step."""
    T.ms.set_case(case)
    adv = T.stepper(T.ms.ECOLOGY_CASE["theta"], T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    n = snaps[0].size
    G = np.zeros((n, 12))
    started = np.zeros(n, bool)
    out = {}
    for k, phi in enumerate(snaps[1:], start=1):
        p = np.minimum(phi.ravel(), T.PHI_CAP)
        new = (~started) & (p >= T.PHI_MIN)
        G[new] = T.seed_state(new.sum())
        started |= new
        act = started & (p >= T.PHI_MIN)
        G[act] = T.rescale(G[act], p[act])
        G[act] = T.rescale(adv(G[act]), p[act])
        chi = np.where(started, G[:, 0] / np.maximum(G[:, 0] + G[:, 1], 1e-300), np.nan)
        out[k] = chi.reshape(phi.shape)
    return out


def main():
    figstyle.apply(size=12)
    snaps = T.macro_snapshots()
    chis = {name: run(snaps, case) for case, name in T.CASES.items()}
    kz = T.K.N // 2
    ext = (0, T.K.L, 0, T.K.L)
    fig, ax = plt.subplots(3, len(TIMES), figsize=(13, 9.6))
    for c, k in enumerate(TIMES):
        phi = snaps[k][:, :, kz].T
        im0 = ax[0, c].imshow(phi, origin="lower", extent=ext, vmin=0, vmax=1, cmap="Blues")
        ax[0, c].set_title(rf"$T^* = {k * T.DT:.1f}$")
        for r, name in enumerate(chis, start=1):
            chi = chis[name][k][:, :, kz].T
            chi = np.where(phi >= 0.05, chi, np.nan)     # show composition where there is biofilm
            im = ax[r, c].imshow(chi, origin="lower", extent=ext, vmin=0, vmax=1, cmap="PuOr_r")
            ax[r, c].set_facecolor("#eef1f4")
        for r in range(3):
            a = ax[r, c]
            a.grid(False)
            a.plot([T.K.L], [T.K.L], marker="s", ms=9, color="#3a9d5d", clip_on=False)
            a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
            if c:
                a.set_yticklabels([])
            if r < 2:
                a.set_xticklabels([])
    ax[0, 0].set_ylabel("amount $\\phi$\n$y$ [$\\mu$m]")
    ax[1, 0].set_ylabel("case 3 (coexistence)\n$y$ [$\\mu$m]")
    ax[2, 0].set_ylabel("case 6 (one species wins)\n$y$ [$\\mu$m]")
    for a in ax[2]:
        a.set_xlabel(r"$x$ [$\mu$m]")
    fig.colorbar(im0, ax=ax[0, :], label=r"$\phi$", shrink=0.9, pad=0.015)
    fig.colorbar(im, ax=ax[1:, :], label=r"$\phi_1/(\phi_1+\phi_2)$", shrink=0.9, pad=0.015)
    fig.suptitle("Composition along a spreading front (Python); Klempt et al. 2024 test case 4.1 field, "
                 "Klempt et al. 2026 cases 3 and 6;\nmid-plane $z = 10\\ \\mu$m, nutrient on the edge "
                 "at the top right (green), grey: no biofilm ($\\phi < 0.05$)", fontsize=12)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)
    for name, d in chis.items():
        w = np.where(snaps[10] >= T.PHI_MIN, snaps[10], 0.0)
        print(f"{name}: biomass-weighted mean at T* = 1: {np.nansum(w * d[10]) / w.sum():.3f}")


if __name__ == "__main__":
    main()
