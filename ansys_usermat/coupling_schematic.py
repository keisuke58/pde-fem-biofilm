#!/usr/bin/env python3
"""Schematic of the one-way micro-macro coupling, for slides and the thesis.

    python ansys_usermat/coupling_schematic.py   -> assets/fig_coupling_schematic.png

(1) Macro: the 3D growth field of Klempt et al. 2024 (the partner's element)
decides the amount phi_3D at every Gauss point. (2) Micro: one copy of the
point model of Klempt et al. 2026 per Gauss point decides only the
composition. (3) Growth and stress use phi_3D. The composition is an output;
it does not act back on the field (two-way coupling is the outlook).

The chi(t) inset is computed, not drawn: case 3, rescaled scheme
(composition_reference.reference, phi_3D = 1, s = 0.15, phi_cap = 0.9).
The stress sketch in (3) is a sketch, not a result.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_coupling_schematic.png"
INK, MUTED, LINE = "#1f2933", "#616e7c", "#9aa5b1"
PANEL = {"macro": ("#e8f1fb", "#1d4e89"), "micro": ("#fdf0e6", "#b45309"),
         "mech": ("#eef1f4", "#3e4c59")}
SP1, SP2 = "#e8833a", "#7a5195"
SEED, CELL = "#1d4e89", "#c9daee"

# panel geometry (axes fraction): x0, width; shared y0, height
Y0, H, HEAD = 0.20, 0.66, 0.075
AR = 14 / 6.4                 # figure width / height: one x unit is AR y units
PX = {"macro": (0.015, 0.300), "micro": (0.350, 0.300), "mech": (0.685, 0.300)}


def panel(ax, key, num, title, subtitle):
    x, w = PX[key]
    fill, accent = PANEL[key]
    ax.add_patch(FancyBboxPatch((x, Y0), w, H, boxstyle="round,pad=0,rounding_size=0.015",
                                fc=fill, ec=accent, lw=1.2, zorder=1))
    ax.add_patch(FancyBboxPatch((x, Y0 + H - HEAD), w, HEAD,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc=accent, ec=accent, lw=1.2, zorder=2))
    ax.add_patch(Rectangle((x, Y0 + H - HEAD), w, HEAD / 2, fc=accent, ec="none", zorder=2))
    ax.plot(x + 0.028, Y0 + H - HEAD / 2, "o", ms=24, mfc="white", mec="none", zorder=3)
    ax.text(x + 0.028, Y0 + H - HEAD / 2, num, ha="center", va="center", fontsize=13,
            weight="bold", color=accent, zorder=4)
    ax.text(x + 0.056, Y0 + H - HEAD / 2, title, ha="left", va="center", fontsize=14.5,
            weight="bold", color="white", zorder=4)
    ax.text(x + w / 2, Y0 + H - HEAD - 0.03, subtitle, ha="center", va="center",
            fontsize=11, style="italic", color=MUTED, zorder=4)


def eqs(ax, x, y, lines, size=12.5, gap=0.058):
    for k, s in enumerate(lines):
        ax.text(x, y - k * gap, s, ha="left", va="center", fontsize=size, color=INK, zorder=4)


def arrow(ax, p, q, color=INK, ls="-", lw=1.8, conn="arc3", z=5):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=18, lw=lw, ls=ls,
                                 color=color, connectionstyle=conn, zorder=z,
                                 shrinkA=0, shrinkB=0))


def cube_inset(fig, rect):
    ax = fig.add_axes(rect, projection="3d")
    n = 8
    filled = np.ones((n, n, n), bool)
    i, j, _ = np.indices((n, n, n))
    filled[(i >= 3) & (j < 3)] = False
    fc = np.empty((n, n, n, 4))
    fc[:] = matplotlib.colors.to_rgba(CELL)
    fc[3, 3, 3] = matplotlib.colors.to_rgba(SEED)
    ax.voxels(filled, facecolors=fc, edgecolor=(0.1, 0.15, 0.2, 0.18), lw=0.25)
    ax.set_axis_off()
    ax.set_facecolor("none")
    ax.view_init(elev=22, azim=-60)
    ax.set_box_aspect((1, 1, 1))
    return ax


def element_with_gauss_points(ax, cx, cy, s):
    """Oblique sketch of one hexahedral element with its 2x2x2 Gauss points."""
    d = np.array([0.45, 0.32 * AR]) * s                  # depth direction
    front = np.array([[0, 0], [s, 0], [s, s * AR], [0, s * AR]]) + [cx, cy]
    back = front + d
    for poly, a in ((back, 0.5), (front, 1.0)):
        ax.add_patch(Polygon(poly, closed=True, fc="white", ec=INK, lw=1.0, alpha=a, zorder=4))
    for k in range(4):
        ax.plot(*zip(front[k], back[k]), color=INK, lw=1.0, zorder=4)
    g = 0.5 - 0.5 / np.sqrt(3)
    pts = []
    for zz in (g, 1 - g):
        for yy in (g, 1 - g):
            for xx in (g, 1 - g):
                pts.append(np.array([cx + xx * s, cy + yy * s * AR]) + zz * d)
    for k, p in enumerate(pts):
        hot = k == 1
        ax.plot(*p, "o", ms=10 if hot else 7, mfc=PANEL["micro"][1] if hot else LINE,
                mec="white", mew=0.8, zorder=6 if hot else 5)
    return pts[1]


def chi_inset(fig, rect):
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
    ax.stackplot(t, c1, 1 - c1, colors=(SP1, SP2), alpha=0.9, lw=0)
    ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    ax.set_xlabel(r"time $T^*$", fontsize=10, labelpad=1)
    ax.set_ylabel("share", fontsize=10, labelpad=1)
    ax.tick_params(labelsize=8.5, length=2, pad=1.5)
    ax.text(0.5, 0.27, "species 1", ha="center", va="center", fontsize=10, color="white")
    ax.text(0.5, 0.80, "species 2", ha="center", va="center", fontsize=10, color="white")
    return ax


def stress_sketch(ax, cx, cy, h):
    """Sketch: the seed grows, is held by its neighbours, ends in compression."""
    n = 5
    hv = h * AR
    x0, y0 = cx - n * h / 2, cy - n * hv / 2
    for i in range(n):
        for j in range(n):
            r = max(abs(i - 2), abs(j - 2))
            col = SEED if r == 0 else ("#8fb3dc" if r == 1 else "#dde7f3")
            ax.add_patch(Rectangle((x0 + i * h, y0 + j * hv), h, hv, fc=col, ec="white",
                                   lw=1.2, zorder=4))
    c = np.array([cx, cy])
    for v in ((1, 0), (-1, 0), (0, AR), (0, -AR)):
        v = np.array(v, float)
        arrow(ax, c + v * 0.62 * h, c + v * 1.45 * h, color=PANEL["micro"][1], lw=1.6, z=6)
        arrow(ax, c + v * 2.45 * h, c + v * 1.6 * h, color=INK, lw=1.4, z=6)


def main():
    figstyle.apply(12)
    fig = plt.figure(figsize=(14, 6.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.set_aspect("auto")

    # (1) macro --------------------------------------------------------------
    panel(ax, "macro", "1", "Macro: 3D growth field", "Klempt et al. 2024, partner's element")
    x, _ = PX["macro"]
    eqs(ax, x + 0.025, 0.695, [
        r"$\dot\phi = \beta\,\Delta\phi + k_\alpha\,\alpha + $ front term",
        r"$\dot\alpha = k_\alpha\,\phi$,   $\alpha(0) = 1$",
    ])
    ax.text(x + 0.025, 0.575, r"gives the amount $\phi_{3D}(\mathbf{x},t)$",
            fontsize=12.5, color=PANEL["macro"][1], weight="bold", va="center", zorder=4)
    cube_inset(fig, [x + 0.005, 0.205, 0.20, 0.35])
    ax.text(x + 0.215, 0.315, "2 mm cube\n512 elements\nseed (dark)", fontsize=10.5,
            color=MUTED, va="center", linespacing=1.3, zorder=4)

    # (2) micro --------------------------------------------------------------
    panel(ax, "micro", "2", "Micro: point model", "Klempt et al. 2026, Table 1 cases")
    x, w = PX["micro"]
    hot = element_with_gauss_points(ax, x + 0.03, 0.50, 0.075)
    ax.text(x + 0.165, 0.665, "one copy at every\nGauss point", fontsize=11.5, color=INK,
            va="center", linespacing=1.25, zorder=4)
    ax.text(x + 0.16, 0.585, r"decides only the" "\n" r"composition $\phi_i\,/\,(\phi_1 + \phi_2)$",
            fontsize=11, color=PANEL["micro"][1], va="center", linespacing=1.35, zorder=4)
    arrow(ax, (hot[0] + 0.004, hot[1] - 0.012), (x + 0.11, 0.43), color=PANEL["micro"][1],
          lw=1.3, conn="arc3,rad=0.25")
    chi_inset(fig, [x + 0.075, 0.255, 0.20, 0.175])

    # (3) growth and stress --------------------------------------------------
    panel(ax, "mech", "3", "Growth and stress", "partner's element, Klempt 2024 stiffness")
    x, w = PX["mech"]
    eqs(ax, x + 0.025, 0.695, [
        r"$\dot\alpha = k_\alpha\,\phi_{3D}$,   $\mathbf{F}_g = \alpha\,\mathbf{I}$",
        r"$E(\phi) = (\phi^2 + f)\,E_{bio}$",
        r"$E_{bio}$ = 10 Pa,  $\nu$ = 0.49",
    ])
    stress_sketch(ax, x + 0.085, 0.375, 0.026)
    ax.text(x + 0.17, 0.375, "the seed grows,\nits neighbours\nhold it back:\nseed in compression,\n"
            "neighbours sheared\n(sketch)", fontsize=10.5, color=MUTED, va="center", linespacing=1.3,
            zorder=4)

    # connections ------------------------------------------------------------
    y_mid = 0.555
    arrow(ax, (PX["macro"][0] + PX["macro"][1], y_mid), (PX["micro"][0], y_mid))
    ax.text(0.3325, y_mid + 0.03, r"$\phi_{3D}$", ha="center", fontsize=13, color=INK)
    xa = PX["macro"][0] + PX["macro"][1] / 2
    xc = PX["mech"][0] + PX["mech"][1] / 2
    ax.plot([xa, xa, xc], [Y0 + H, 0.945, 0.945], color=INK, lw=1.8, zorder=5,
            solid_capstyle="round")
    arrow(ax, (xc, 0.945), (xc, Y0 + H + 0.004))
    ax.text(0.5, 0.958, r"$\phi_{3D}$: growth uses the amount only", ha="center",
            va="bottom", fontsize=12.5, color=INK)
    xm = PX["micro"][0] + PX["micro"][1] * 0.72
    arrow(ax, (xm, Y0), (xm, 0.105))
    ax.text(xm, 0.085, r"output: species $\phi_i(\mathbf{x},t)$, $\psi_i(\mathbf{x},t)$", ha="center",
            va="top", fontsize=12.5, color=INK)
    xb = PX["micro"][0] + PX["micro"][1] * 0.22
    ax.plot([xb, xb, xa], [Y0, 0.13, 0.13], color=MUTED, lw=1.5, ls=(0, (4, 3)), zorder=5)
    arrow(ax, (xa, 0.13), (xa, Y0 - 0.004), color=MUTED, ls=(0, (4, 3)), lw=1.5)
    ax.text((xa + xb) / 2, 0.115, "two-way coupling (outlook)", ha="center", va="top",
            fontsize=11, color=MUTED, style="italic")

    ax.text(0.5, 0.012,
            r"Point model run two ways: A) rescaled every step to $\phi_{3D}$;  "
            r"B) started when $\phi_{3D}$ first reaches $\phi_{min}$, then independent.   "
            r"Assumed, not from a paper: $s$ = 0.15, $\phi_{cap}$ = 0.9, $\phi_{min}$ = 0.01.",
            ha="center", va="bottom", fontsize=10.5, color=MUTED)

    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=220, bbox_inches="tight", facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
