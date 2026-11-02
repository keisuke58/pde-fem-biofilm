#!/usr/bin/env python3
"""Set-up of the separate example with antibiotic (thesis Fig. ch5-3d), drawn like fig_model_setup.png: species A in
the seed of the ANSYS deck, species B in four small blocks along the edges parallel to z at the seed's height, the nutrient
held on the element layer y < -0.75 mm, the antibiotic held on the outer element shell from T* = 1/3, and the three
constrained corner nodes (rigid-body constraints of the stress solve). Shown on the 8^3 element grid of the deck; the
runs use the same regions on 64^3 cells (coupled_living_consumption.py --v2 --ansys-bc).

    python ansys_usermat/apdl/fig_setup_separate_ab.py [mixed]   -> assets/fig_setup_{separate,mixed}_ab.png / .pdf
"mixed": both species in the seed in equal amounts (phi_1 = phi_2 = 0.275), no edge blocks; the rest as above.
Orientation as the 3D results and the videos (fig_coupled_3d_panels.py, anim_coupled_tvd_3d.py): y points up, camera
from (3.6, 2.6, 4.4), the octant x, y, z > 0 removed in (b) and (c). matplotlib keeps its third axis vertical, so the
data are drawn with the plot axes (z, x, y), a cyclic permutation that keeps the handedness.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
from model_setup_fig import DECK, read_deck, to_grid  # noqa: E402

MIXED = len(sys.argv) > 1 and sys.argv[1] == "mixed"
OUT = HERE.parents[1] / "assets" / ("fig_setup_mixed_ab.png" if MIXED else "fig_setup_separate_ab.png")
SPAB = "#7048a8"                                       # both species in the seed (mixed example)
SPA, SPB, NUTR, AB, MESH, INK, FIX = "#3b5bdb", "#c92a2a", "#3a9d5d", "#e8890c", "#9aa5b1", "#1f2933", "#7a1f1f"


VIEW = dict(elev=float(np.degrees(np.arctan2(2.6, np.hypot(3.6, 4.4)))), azim=float(np.degrees(np.arctan2(3.6, 4.4))))
P = lambda a: np.transpose(a, (2, 0, 1))         # data (x, y, z) -> plot axes (z, x, y): y vertical


def frame(ax):
    for a in (-1, 1):
        for b in (-1, 1):
            ax.plot([-1, 1], [a, a], [b, b], color=INK, lw=0.8)
            ax.plot([a, a], [-1, 1], [b, b], color=INK, lw=0.8)
            ax.plot([a, a], [b, b], [-1, 1], color=INK, lw=0.8)
    ax.set_xlabel("$z$ [mm]"); ax.set_ylabel("$x$ [mm]"); ax.set_zlabel("$y$ [mm]")
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(**VIEW)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("none")


def main():
    figstyle.apply(size=13)
    xyz, cen, seed_ids, nut_ids, corners = read_deck(DECK)
    seed, nut = to_grid(cen, seed_ids), to_grid(cen, nut_ids)
    c = np.linspace(-1, 1, 9)
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    m = 0.5 * (c[:-1] + c[1:])
    XC, YC, ZC = np.meshgrid(m, m, m, indexing="ij")
    zmax = np.abs(ZC[seed]).max() + 1e-9
    seedc = SPAB if MIXED else SPA
    spb = np.zeros_like(seed) if MIXED else (np.abs(XC) > 0.75) & (np.abs(YC) > 0.75) & (np.abs(ZC) <= zmax)
    shell = (np.abs(XC) > 0.75) | (np.abs(YC) > 0.75) | (np.abs(ZC) > 0.75)
    front = (XC > 0) & (YC > 0) & (ZC > 0)           # the octant towards the camera, removed in (b) and (c)
    near = 3.6 * XC + 2.6 * YC > 0                   # blocks of B on the camera side of the seed (camera at (3.6, 2.6, 4.4))
    spb_back, spb_near = spb & ~near, spb & near
    seed, nut, spb, shell, front, spb_back, spb_near = (P(v) for v in (seed, nut, spb, shell, front, spb_back, spb_near))
    X, Y, Z = P(Z), P(X), P(Y)                       # plot coordinates

    fig = plt.figure(figsize=(16.5, 6.6))
    titles = ("(a) start: A and B mixed in the seed, supports" if MIXED else "(a) start: species A and B, supports",
              "(b) nutrient: $c = 1$ held from $T^* = 0$",
              "(c) antibiotic: $a = 1$ held from $T^* = 1/3$")
    for p in range(3):
        ax = fig.add_subplot(1, 3, p + 1, projection="3d")
        frame(ax)                                    # first, so that with computed_zorder off it lies behind the voxels
        if p == 0:
            ax.computed_zorder = False               # supports drawn on top (the node at (-1, -1, -1) lies behind the seed)
            # drawn in depth order (computed_zorder is off for the supports): B behind the seed, the seed, B in front
            ax.voxels(X, Y, Z, spb_back, facecolors=SPB + "dd", edgecolors="white", linewidth=0.3)
            ax.voxels(X, Y, Z, seed, facecolors=seedc + "dd", edgecolors="white", linewidth=0.3)
            if spb_near.any():
                ax.voxels(X, Y, Z, spb_near, facecolors=SPB + "dd", edgecolors="white", linewidth=0.3)
            for (x, y, z), dofs in corners.items():
                x, y, z = z, x, y                    # plot coordinates
                ax.scatter([x], [y], [z], s=70, marker="^", color=FIX, depthshade=False, zorder=10)
                ax.text(x * 1.18, y * 1.18, z - 0.22, ", ".join(d.strip().lower() for d in dofs), fontsize=11,
                        color=FIX, ha="center", zorder=11,
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
        elif p == 1:
            ax.voxels(X, Y, Z, ~nut & ~front, facecolors="#ffffff10", edgecolors=MESH + "55", linewidth=0.2)
            ax.voxels(X, Y, Z, nut & ~front, facecolors=NUTR + "99", edgecolors=NUTR, linewidth=0.3)
            ax.voxels(X, Y, Z, seed & ~front, facecolors=seedc + "dd", edgecolors="white", linewidth=0.3)
        else:
            ax.voxels(X, Y, Z, shell & ~front, facecolors=AB + "1a", edgecolors=AB + "66", linewidth=0.25)
            ax.voxels(X, Y, Z, seed & ~front, facecolors=seedc + "dd", edgecolors="white", linewidth=0.3)
        ax.set_title(titles[p], fontsize=14, y=1.06)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=seedc), plt.Rectangle((0, 0), 1, 1, fc=SPB),
               plt.Rectangle((0, 0), 1, 1, fc=NUTR, alpha=0.6), plt.Rectangle((0, 0), 1, 1, fc=AB, alpha=0.3),
               plt.Line2D([], [], marker="^", ls="", color=FIX)]
    labels = ["species A: seed, $\\phi_{3D} = 1$ at the start",
              "species B: four blocks along the edges parallel to $z$",
              "nutrient held: element layer on $y = -1$ mm",
              "antibiotic held: outer element shell",
              "constrained corner nodes (rigid-body supports)"]
    if MIXED:                                         # no edge blocks: one seed entry instead of A and B
        handles.pop(1)
        labels[:2] = ["seed: $\\phi_{3D} = 1$, A and B in equal amounts ($\\phi_1 = \\phi_2 = 0.275$)"]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=12, frameon=False)
    fig.text(0.5, 0.155, "$y$ points up, as in the 3D results and the videos.  (b), (c): octant $x, y, z > 0$ removed.  Zero flux of $c$ and $a$ elsewhere; "
             "no other displacement constraints (free swelling).", ha="center", fontsize=12, color=INK)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.9, bottom=0.27, wspace=-0.05)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05)
    print("wrote", OUT, "| seed", int(seed.sum()), "| B", int(spb.sum()), "| nutrient", int(nut.sum()),
          "| shell", int(shell.sum()), "| corners", corners)


if __name__ == "__main__":
    main()
