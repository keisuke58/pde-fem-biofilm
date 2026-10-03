#!/usr/bin/env python3
"""The one-way coupling written out as equations: one increment at one
Gauss point. A companion to coupling_schematic.py (the overview); this one
shows each step the coupled ANSYS run takes, with the equation it uses and
where that equation comes from.

    python ansys_usermat/coupling_equations_fig.py   -> assets/fig_coupling_equations.png

(1) The partner's element updates its nodal fields (Klempt et al. 2024,
Eq. 34/35). (2) At each Gauss point the material routine integrates Eq. 36
and computes the stress from F_g = alpha I. (3) For two species, the point
model of Klempt et al. 2026 is scaled to the field's amount, advanced by
s dt and scaled again; it only decides the composition. Nothing in (3) acts
back on (1) or (2). Values not taken from a paper are marked in the footer.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_coupling_equations.png"
INK, MUTED = "#1f2933", "#616e7c"
COL = {"field": ("#e8f1fb", "#1d4e89"), "mech": ("#eef1f4", "#3e4c59"),
       "micro": ("#fdf0e6", "#b45309")}
Y0, H, HEAD = 0.17, 0.70, 0.065
PX = {"field": (0.012, 0.300), "mech": (0.350, 0.300), "micro": (0.688, 0.300)}


def panel(ax, key, num, title, sub):
    x, w = PX[key]
    fill, acc = COL[key]
    ax.add_patch(FancyBboxPatch((x, Y0), w, H, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=fill, ec=acc, lw=1.2, zorder=1))
    ax.add_patch(FancyBboxPatch((x, Y0 + H - HEAD), w, HEAD,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc=acc, ec=acc, lw=1.2, zorder=2))
    ax.add_patch(Rectangle((x, Y0 + H - HEAD), w, HEAD / 2, fc=acc, ec="none", zorder=2))
    ax.plot(x + 0.024, Y0 + H - HEAD / 2, "o", ms=21, mfc="white", mec="none", zorder=3)
    ax.text(x + 0.024, Y0 + H - HEAD / 2, num, ha="center", va="center", fontsize=12,
            weight="bold", color=acc, zorder=4)
    ax.text(x + 0.048, Y0 + H - HEAD / 2, title, ha="left", va="center", fontsize=13.5,
            weight="bold", color="white", zorder=4)
    ax.text(x + w / 2, Y0 + H - HEAD - 0.027, sub, ha="center", va="center",
            fontsize=10.5, style="italic", color=MUTED, zorder=4)


def lines(ax, key, y, rows, dy=0.058):
    """rows: (kind, text). kind 'step' = numbered label, 'eq' = equation,
    'note' = small grey remark."""
    x, w = PX[key]
    acc = COL[key][1]
    for kind, txt in rows:
        if kind == "step":
            ax.text(x + 0.016, y, txt, ha="left", va="center", fontsize=11,
                    weight="bold", color=acc, zorder=4)
            y -= dy * 0.72
        elif kind == "eq":
            ax.text(x + 0.030, y, txt, ha="left", va="center", fontsize=12.5,
                    color=INK, zorder=4)
            y -= dy
        elif kind == "note":
            ax.text(x + 0.030, y + 0.012, txt, ha="left", va="center", fontsize=9.5,
                    color=MUTED, zorder=4)
            y -= dy * 0.62
        elif kind == "gap":
            y -= dy * 0.35
    return y


def arrow(ax, p, q, color=INK, ls="-", lw=1.6, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=15, color=color,
                                 lw=lw, ls=ls, connectionstyle=f"arc3,rad={rad}", zorder=5))


def main():
    figstyle.apply()
    fig = plt.figure(figsize=(16, 8.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.968, r"One time increment $\Delta t$ at one Gauss point of the coupled ANSYS run",
            ha="center", va="center", fontsize=15, color=INK, weight="bold")

    panel(ax, "field", "1", "Element: growth field", "Klempt et al. 2024, partner's element (unchanged)")
    panel(ax, "mech", "2", "Material routine: growth, stress", "added in this work")
    panel(ax, "micro", "3", "Point model: composition", "Klempt et al. 2026, added in this work")

    top = Y0 + H - HEAD - 0.075
    lines(ax, "field", top, [
        ("step", "a  nutrient, quasi-static (Eq. 35)"),
        ("eq", r"$d\,\Delta c = g\,\phi,\qquad c = 1$ on the source"),
        ("step", "b  biofilm field (Eq. 34)"),
        ("eq", r"$\dot\phi = \beta\,\Delta\phi + k_\alpha\,\alpha"
               r" - r\,\dfrac{c}{k+c}\,\nabla\phi\cdot\dfrac{\nabla c}{|\nabla c|}$"),
        ("note", "diffusion, local source, front growth towards the nutrient"),
        ("eq", r"$\dot\alpha = k_\alpha\,\phi$  (Eq. 36, inside the element)"),
        ("step", "c  value at the Gauss point"),
        ("eq", r"$\phi_{3D} = \sum_I N_I\,\phi_I,\qquad 0 \leq \phi_{3D} \leq 1$"),
        ("note", r"$N_I$: shape functions, $\phi_I$: nodal values"),
    ])

    lines(ax, "mech", top, [
        ("step", "a  growth (Eq. 36, explicit)"),
        ("eq", r"$\alpha_{n+1} = \alpha_n + k_\alpha\,\phi_{3D}\,\Delta t,\qquad \alpha(0) = 1$"),
        ("step", "b  kinematics"),
        ("eq", r"$\mathbf{F}_g = \alpha_{n+1}\,\mathbf{I},\qquad"
               r"\mathbf{F}_e = \mathbf{F}\,\mathbf{F}_g^{-1},\qquad J_e = \det\mathbf{F}_e$"),
        ("note", r"growth is $\alpha - 1$ (stored as $\alpha - 1$ in the code)"),
        ("step", "c  stiffness and stress"),
        ("eq", r"$E(\phi_{3D}) = (\phi_{3D}^{\,2} + f)\,E_{bio},\qquad \nu = 0.49$"),
        ("eq", r"$\boldsymbol{\sigma} = \hat{\boldsymbol{\sigma}}(\mathbf{F}_e;\,E,\nu)$  (neo-Hookean)"),
        ("step", "d  equilibrium (ANSYS Newton iteration)"),
        ("eq", r"$\int_\Omega \boldsymbol{\sigma} : \delta\boldsymbol{\varepsilon}\;dv = 0$"),
        ("note", "a free part grows without stress; constraint gives stress"),
        ("note", r"next increment: $\alpha_{n+1}$, $\phi_i$, $\psi_i$ carried as state variables"),
    ])

    lines(ax, "micro", top, [
        ("step", "a  amount the point model sees"),
        ("eq", r"$\hat\phi = \mathrm{min}(\phi_{3D},\,\phi_{cap})$"),
        ("step", "b  scale in"),
        ("eq", r"$\phi_i \leftarrow \phi_i\,\dfrac{\hat\phi}{\phi_1+\phi_2},\qquad \phi_0 = 1 - \hat\phi$"),
        ("step", r"c  advance by $s\,\Delta t$ ($n_{sub}$ implicit steps)"),
        ("eq", r"$(\eta_{\phi,i}+\eta_i\psi_i^2)\,\dot\phi_i + \eta_i\phi_i\psi_i\,\dot\psi_i"
               r" = c^*\psi_i\sum_j A_{ij}\phi_j\psi_j - \gamma$"),
        ("eq", r"$\eta_i(\phi_i\psi_i\,\dot\phi_i + \phi_i^2\,\dot\psi_i)"
               r" = c^*\phi_i\sum_j A_{ij}\phi_j\psi_j - \alpha^* b_i\psi_i$"),
        ("note", r"$n_{sub} = \lceil s\,\Delta t / 10^{-4}\rceil$; barrier terms omitted"),
        ("step", "d  scale out and store"),
        ("eq", r"$\phi_1 + \phi_2 = \hat\phi$;  $\phi_i,\ \psi_i$ as state variables"),
        ("eq", r"output: $\phi_1/(\phi_1+\phi_2)$,  $\psi_i$"),
    ])

    # flows
    xf, wf = PX["field"]; xm, wm = PX["mech"]; xc, wc = PX["micro"]
    ymid = Y0 + 0.30
    arrow(ax, (xf + wf, ymid), (xm, ymid))
    ax.text((xf + wf + xm) / 2, ymid + 0.025, r"$\phi_{3D}$", ha="center", fontsize=13, color=INK)
    ax.plot([xf + wf * 0.55, xf + wf * 0.55, xc + wc * 0.5], [Y0 + H, 0.902, 0.902],
            color=INK, lw=1.6, zorder=5)
    arrow(ax, (xc + wc * 0.5, 0.902), (xc + wc * 0.5, Y0 + H))
    ax.text(xm + wm / 2, 0.913, r"$\phi_{3D}$: the amount only", ha="center", fontsize=12, color=INK)
    # no feedback: a crossed, dashed line without arrowhead
    yb = Y0 - 0.045
    ax.plot([xc + wc * 0.5, xc + wc * 0.5, xf + wf * 0.5, xf + wf * 0.5], [Y0, yb, yb, Y0],
            color=MUTED, lw=1.3, ls="--", zorder=5)
    ax.plot(xm + wm * 0.5, yb, "o", ms=17, mfc="white", mec=MUTED, mew=1.3, zorder=6)
    ax.text(xm + wm * 0.5, yb, r"$\times$", ha="center", va="center", fontsize=15,
            color=MUTED, zorder=7)
    ax.text(xm + wm * 0.5, yb - 0.032, "composition does not act back on (1) or (2): "
            "one-way coupling (two-way is the outlook)", ha="center", fontsize=11,
            style="italic", color=MUTED)
    ax.text(0.5, 0.035,
            r"Klempt et al. 2024 Table 2: $k_\alpha = 10^{-3}$ per $T^*$, $E_{bio} = 10$ Pa, $\nu = 0.49$.   "
            r"Klempt et al. 2026: $c^* = 100$, $\alpha^* = 10$, $A_{ij}$, $b_i$, $\eta_i$ per Table 1 case.   "
            "\nNot from a paper: "
            r"$s = 0.15$, $\phi_{cap} = 0.9$, $f = 10^{-3}$ (void stiffness), "
            r"and $\beta$, $r$, $k$, $d$, $g$ from the partner's example input.",
            ha="center", va="center", fontsize=10.5, color=MUTED, linespacing=1.6)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
