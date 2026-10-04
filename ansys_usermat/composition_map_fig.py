#!/usr/bin/env python3
"""Composition map on a section for a prescribed spreading amount (demo).

The amount phi_3D(r, t) is PRESCRIBED, not computed by ANSYS: a front
spreading radially from the seed,
    phi_3D = 0.9 * 0.5 * (1 - tanh((r - R(t)) / w)),  R = 0.15 + 0.75 t mm,
w = 0.06 mm. At each radius the coupled schemes run as in the material
routine (composition_reference.reference / reference_age; s = 0.15,
phi_cap = 0.9, phi_min = 0.01, dt = 0.01). The field is radially symmetric,
so each radius is run once and mapped onto the section.

What it shows: in the coupled model the composition follows the history of
each point. Material reached early by the front (the centre) is further along
the point model's dynamics than material the front has just reached (the
rim): for case 6, species 2 has taken over at the centre while the rim is
still mixed. Neither scheme is told anything about space; the pattern comes
from when the field's amount arrives.

    python ansys_usermat/composition_map_fig.py   -> assets/fig_composition_map.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "coupling"))
import composition_reference as cref  # noqa: E402
import figstyle  # noqa: E402
import material_server as ms  # noqa: E402

OUT = HERE.parent / "assets" / "fig_composition_map.png"
DT, S, PHI_CAP, PHI_MIN, W = 0.01, 0.15, 0.9, 0.01, 0.06
TIMES = (0.3, 0.6, 1.0)
R = np.linspace(0.0, 1.0, 41)
SP1, SP2 = "#e8833a", "#7a5195"
CHI = LinearSegmentedColormap.from_list("chi", [SP2, "#b9a3c9", "#f7f3fa"])


def phi3(r, t):
    return PHI_CAP * 0.5 * (1.0 - np.tanh((r - (0.15 + 0.75 * t)) / W))


def radial(case, scheme):
    """chi_1(r) at TIMES; nan where the point model has not started."""
    ms.set_case(case)
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    n = int(round(1.0 / DT))
    idx = [int(round(t / DT)) - 1 for t in TIMES]
    out = np.full((len(TIMES), len(R)), np.nan)
    for j, r in enumerate(R):
        series = [phi3(r, (k + 1) * DT) for k in range(n)]
        if scheme == "A":
            st = cref.reference(series, th, hp, DT, phi_min=PHI_MIN, s=S, phi_cap=PHI_CAP)
            first = next((k for k, p in enumerate(series) if p >= PHI_MIN), n)
        else:
            st = cref.reference_age(series, th, hp, DT, phi_min=PHI_MIN, s=S)
            first = next((k for k, g in enumerate(st) if g[6] > 0), n)
        for i, k in enumerate(idx):
            if k >= first:
                g = st[k]
                out[i, j] = g[0] / (g[0] + g[1])
    ms.set_case(None)
    return out


def to_map(prof, x):
    rr = np.hypot(*np.meshgrid(x, x))
    m = np.interp(rr, R, np.nan_to_num(prof, nan=-1.0), right=-1.0)
    return np.where(m < 0, np.nan, m)


def main():
    figstyle.apply(12)
    x = np.linspace(-1.0, 1.0, 241)
    ext = (-1, 1, -1, 1)
    rows = [("amount", None, None)] + [(f"{lab}", c, sch) for c, lab in
                                       (("2sp_case6", "case 6"),) for sch in ("A", "B")]
    fig, axs = plt.subplots(len(rows), len(TIMES), figsize=(10.6, 10.4),
                            constrained_layout=True)
    data = {}
    for _, c, sch in rows[1:]:
        data[(c, sch)] = radial(c, sch)
        print(c, sch, "chi1 at r=0:", np.round(data[(c, sch)][:, 0], 4),
              " at the front (last started r):",
              [round(float(v[~np.isnan(v)][-1]), 4) for v in data[(c, sch)]])
    for i, (lab, c, sch) in enumerate(rows):
        for jt, t in enumerate(TIMES):
            ax = axs[i, jt]
            ax.grid(False)
            if c is None:
                im0 = ax.imshow(phi3(np.hypot(*np.meshgrid(x, x)), t), origin="lower",
                                extent=ext, cmap="Blues", vmin=0, vmax=1)
                ax.set_title(rf"$T^*$ = {t:g}")
            else:
                im1 = ax.imshow(to_map(data[(c, sch)][jt], x), origin="lower", extent=ext,
                                cmap=CHI, vmin=0, vmax=0.5)
            ax.set_xticks([-1, 0, 1])
            ax.set_yticks([-1, 0, 1])
            ax.set_facecolor("#c3c9d1")
            if jt == 0:
                ax.set_ylabel({None: r"amount $\phi_{3D}$" "\n(prescribed)\n\ny [mm]",
                               "A": "case 6, scheme A\n(rescaled)\n\ny [mm]",
                               "B": "case 6, scheme B\n(independent)\n\ny [mm]"}[sch])
            if i == len(rows) - 1:
                ax.set_xlabel("x [mm]")
    fig.colorbar(im0, ax=axs[0, :], shrink=0.85, label=r"$\phi_{3D}$")
    cb = fig.colorbar(im1, ax=axs[1:, :], shrink=0.6,
                      label=r"share of species 1, $\phi_1/(\phi_1+\phi_2)$")
    cb.set_ticks([0, 0.25, 0.5])
    cb.set_ticklabels(["0: species 2 only", "0.25", "0.5: start (equal)"])
    fig.suptitle("Composition map from the coupled model on a prescribed spreading front "
                 "(demo, not an ANSYS result)\nKlempt et al. 2026 case 6; grey: point model "
                 r"not started ($\phi_{3D} < \phi_{min}$)", fontsize=12.5)
    fig.savefig(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
