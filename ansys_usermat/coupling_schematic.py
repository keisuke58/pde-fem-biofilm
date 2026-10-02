#!/usr/bin/env python3
"""Schematic of the one-way micro-macro coupling, for slides.

    python ansys_usermat/coupling_schematic.py   -> assets/fig_coupling_schematic.png

Left: the 3D growth field of Klempt et al. 2024 (the partner's element)
decides the amount phi_3D at every Gauss point. Middle: at each Gauss point
the point model of Klempt et al. 2026 decides only the composition. Right:
growth and stress use phi_3D. The composition is an output; it does not act
back on the field (two-way coupling is the outlook).
The small chi_1(t) inset is case 3 of the rescaled scheme
(composition_reference.reference, phi_3D = 1, s = 0.15, phi_cap = 0.9).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "assets"
INK, MUTED = "#1f2933", "#52606d"
FILL = {"field": "#dbeafe", "point": "#fde7d2", "mech": "#e2e8f0"}


def box(ax, x, y, w, h, kind, title, lines):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=FILL[kind], ec=INK, lw=1.0))
    ax.text(x + w / 2, y + h - 0.035, title, ha="center", va="top",
            fontsize=12, weight="bold", color=INK)
    ax.text(x + 0.02, y + h - 0.11, "\n".join(lines), ha="left", va="top",
            fontsize=10, color=INK, linespacing=1.55)


def arrow(ax, p, q, label, dy=0.025, style="-|>", ls="-", color=INK):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=16,
                                 lw=1.4, ls=ls, color=color))
    ax.text((p[0] + q[0]) / 2, (p[1] + q[1]) / 2 + dy, label, ha="center",
            va="bottom", fontsize=10, color=color)


def cube_inset(fig, rect):
    ax = fig.add_axes(rect, projection="3d")
    n = 8
    filled = np.ones((n, n, n), bool)
    i, j, _ = np.indices((n, n, n))
    filled[(i >= 3) & (j < 3)] = False
    fc = np.empty((n, n, n, 4))
    fc[:] = matplotlib.colors.to_rgba("#bfd4ee")
    fc[3, 3, 3] = matplotlib.colors.to_rgba("#1d4e89")
    ax.voxels(filled, facecolors=fc, edgecolor=(0, 0, 0, 0.12), lw=0.2)
    ax.set_axis_off()
    ax.set_facecolor("none")
    ax.view_init(elev=22, azim=-60)
    ax.set_box_aspect((1, 1, 1))


def chi_inset(fig, rect):
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(HERE / "coupling"))
    import composition_reference as cref
    import material_server as ms
    ms.set_case("2sp_case3")
    th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    st = cref.reference([1.0] * 10, th, hp, 0.1, phi_min=1e-2, s=0.15, phi_cap=0.9)
    ms.set_case(None)
    t = np.arange(0, 11) * 0.1
    c1 = np.array([0.5] + [g[0] / (g[0] + g[1]) for g in st])
    ax = fig.add_axes(rect)
    ax.stackplot(t, c1, 1 - c1, colors=("#e8833a", "#7a5195"), alpha=0.85)
    ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[0, 1], yticks=[0, 1])
    ax.set_xlabel("$T^*$", fontsize=8, labelpad=-6)
    ax.set_ylabel(r"$\chi_1$, $\chi_2$", fontsize=8, labelpad=-4)
    ax.tick_params(labelsize=7)
    ax.text(0.5, 0.3, "species 1", ha="center", fontsize=7, color="white")
    ax.text(0.5, 0.8, "species 2", ha="center", fontsize=7, color="white")


def main():
    fig = plt.figure(figsize=(13, 5.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")

    box(ax, 0.02, 0.30, 0.28, 0.58, "field", "Macro: 3D growth field",
        ["Klempt et al. 2024, partner's element",
         r"$\dot\varphi = \beta\Delta\varphi + k_a\alpha_K$ + front term",
         r"$\dot\alpha_K = k_a\varphi$",
         r"gives the amount $\varphi_{3D}(x,t)$"])
    box(ax, 0.36, 0.30, 0.28, 0.58, "point", "Micro: point model",
        ["Klempt et al. 2026, Table 1 cases",
         "one copy at every Gauss point",
         r"gives only the composition $\chi_i = \varphi_i/\sum\varphi$",
         ])
    box(ax, 0.70, 0.30, 0.28, 0.58, "mech", "Growth and stress",
        [r"$\dot\alpha = k_a\,\varphi_{3D}$,  $F_g = (1+\alpha)\,I$",
         r"$E(\varphi) = (\varphi^2 + f)\,E_{bio}$",
         r"$E_{bio}$ = 10 Pa, $\nu$ = 0.49 (Klempt 2024)",
         "von Mises and mean stress",
         "",
         "ANSYS, one species: the seed is",
         "compressed; its neighbours carry",
         "about 25 times its von Mises stress"])

    arrow(ax, (0.30, 0.66), (0.36, 0.66), r"$\varphi_{3D}$")
    ax.plot([0.16, 0.16, 0.84], [0.88, 0.95, 0.95], color=INK, lw=1.4)
    arrow(ax, (0.84, 0.95), (0.84, 0.885), "", dy=0)
    ax.text(0.50, 0.958, r"$\varphi_{3D}$: growth uses the amount only", ha="center",
            va="bottom", fontsize=10, color=INK)
    arrow(ax, (0.56, 0.30), (0.56, 0.16), "")
    ax.text(0.56, 0.14, r"output: composition map $\chi_i(x,t)$", ha="center", va="top",
            fontsize=10.5, color=INK)
    ax.plot([0.42, 0.42, 0.16], [0.30, 0.22, 0.22], color=MUTED, lw=1.4, ls="--")
    arrow(ax, (0.16, 0.22), (0.16, 0.295), "", dy=0, ls="--", color=MUTED)
    ax.text(0.29, 0.20, "two-way (outlook): composition acts on the field", ha="center",
            va="top", fontsize=9.5, color=MUTED, style="italic")
    ax.text(0.50, 0.035,
            "Two ways to run the point model: A) rescaled every step to the amount "
            r"$\varphi_{3D}$;  B) started when $\varphi_{3D}$ first reaches $\varphi_{min}$, then independent."
            "\nAssumed (not from a paper): clock $s$ = 0.15, $\\varphi_{cap}$ = 0.9, "
            r"$\varphi_{min}$ = 0.01, each with a sensitivity study.",
            ha="center", va="bottom", fontsize=9.5, color=MUTED)

    cube_inset(fig, [0.07, 0.30, 0.18, 0.30])
    chi_inset(fig, [0.43, 0.41, 0.14, 0.17])
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "fig_coupling_schematic.png", dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", OUT / "fig_coupling_schematic.png")


if __name__ == "__main__":
    main()
