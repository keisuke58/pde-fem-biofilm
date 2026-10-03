#!/usr/bin/env python3
"""3D picture of the ANSYS model of the coupled runs: mesh, seed, nutrient
face and the three constrained corner nodes, read from the partner's deck
(element components BIOFILM1 and NUTRIENT1, the corner node components with
displacement constraints).

    python ansys_usermat/model_setup_fig.py   -> assets/fig_model_setup.png
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
import figstyle  # noqa: E402

DECK = HERE / "apdl" / "ds_oliver_wired_n5trial_inter34only_baseline_at0p0006.dat"
OUT = HERE.parent / "assets" / "fig_model_setup.png"
SEED, NUTR, MESH, INK = "#1d4e89", "#3a9d5d", "#9aa5b1", "#1f2933"


def read_deck(path):
    lines = path.read_text(encoding="latin-1").splitlines()
    i = next(k for k, l in enumerate(lines) if l.lower().startswith("nblock")) + 2
    xyz = {}
    while len(p := lines[i].split()) >= 4:
        xyz[int(p[0])] = np.array([float(v) for v in p[1:4]])
        i += 1
    k = next(j for j, l in enumerate(lines) if l.lower().startswith("eblock")) + 2
    cen = {}
    while len(p := lines[k].split()) >= 19:
        cen[int(p[10])] = np.mean([xyz[int(n)] for n in p[11:19]], 0)
        k += 1

    def comp(name):
        j = next(j for j, l in enumerate(lines) if l.startswith(f"CMBLOCK,{name}"))
        n = int(lines[j].split(",")[3])
        raw, j = [], j + 2
        while len(raw) < n:
            raw += [int(v) for v in lines[j].split()]
            j += 1
        ids, prev = [], None
        for v in raw:                     # a negative entry closes a range
            ids += list(range(prev + 1, -v + 1)) if v < 0 else [v]
            prev = abs(v)
        return ids

    fixed = []
    for j, l in enumerate(lines):
        if l.lower().startswith("d,_cm"):
            name, dof = l.split(",")[1], l.split(",")[2]
            fixed.append((name, dof))
    corners = {}
    for name, dof in fixed:
        node = comp(name)[0]
        corners.setdefault(tuple(xyz[node]), []).append(dof)
    return xyz, cen, comp("BIOFILM1"), comp("NUTRIENT1"), corners


def to_grid(cen, ids, h=0.25, n=8):
    g = np.zeros((n, n, n), bool)
    for e in ids:
        i, j, k = ((cen[e] + 1.0) / h - 0.5).round().astype(int)
        g[i, j, k] = True
    return g


def main():
    figstyle.apply(size=13)
    xyz, cen, seed_ids, nut_ids, corners = read_deck(DECK)
    seed, nut = to_grid(cen, seed_ids), to_grid(cen, nut_ids)
    edges = np.linspace(-1, 1, 9)
    X, Y, Z = np.meshgrid(edges, edges, edges, indexing="ij")

    fig = plt.figure(figsize=(11, 5.4))
    for p, (title, cut) in enumerate((("whole model", False),
                                      ("front half ($x < 0$) removed: the seed inside", True))):
        ax = fig.add_subplot(1, 2, p + 1, projection="3d")
        keep = np.ones_like(seed)
        if cut:
            keep[:4, :, :] = False
        if cut:                      # drawn first: the remaining mesh, faint
            void = keep & ~seed & ~nut
            ax.voxels(X, Y, Z, void, facecolors="#ffffff10", edgecolors=MESH + "55", linewidth=0.2)
        ax.voxels(X, Y, Z, nut & keep, facecolors=NUTR + "88", edgecolors=NUTR, linewidth=0.3)
        ax.voxels(X, Y, Z, seed & keep, facecolors=SEED + "dd", edgecolors="white", linewidth=0.3)
        # outline of the cube
        for a in (-1, 1):
            for b in (-1, 1):
                ax.plot([-1, 1], [a, a], [b, b], color=INK, lw=0.8)
                ax.plot([a, a], [-1, 1], [b, b], color=INK, lw=0.8)
                ax.plot([a, a], [b, b], [-1, 1], color=INK, lw=0.8)
        for (x, y, z), dofs in corners.items():
            ax.scatter([x], [y], [z], s=60, marker="^", color="#c0392b", depthshade=False, zorder=10)
            lab = ", ".join(d.strip().lower() for d in dofs)
            ax.text(x, y, z - 0.3, lab, fontsize=11, color="#c0392b", ha="center")
        ax.set_xlabel("$x$ [mm]"); ax.set_ylabel("$y$ [mm]"); ax.set_zlabel("$z$ [mm]")
        ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=13, azim=128)
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
            axis.pane.set_edgecolor("none")
        ax.set_title(title, fontsize=14, y=0.97)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=SEED), plt.Rectangle((0, 0), 1, 1, fc=NUTR, alpha=0.6),
               plt.Line2D([], [], marker="^", ls="", color="#c0392b")]
    fig.legend(handles, [f"seed: $\\phi = 1$ at the start ({len(seed_ids)} elements)",
                         f"nutrient: $c = 1$ held (layer on $y = -1$ mm, {len(nut_ids)} elements)",
                         "constrained corner nodes (fixed displacement components)"],
               loc="lower center", ncol=3, fontsize=12, frameon=False)
    fig.suptitle(r"ANSYS model: cube of 2 mm, $8\times8\times8 = 512$ hexahedral elements of 0.25 mm",
                 fontsize=15)
    fig.subplots_adjust(left=-0.02, right=1.0, top=0.9, bottom=0.12, wspace=-0.12)
    fig.savefig(OUT, dpi=200, bbox_inches="tight", pad_inches=0.05)
    print("wrote", OUT, "| seed", len(seed_ids), "| nutrient", len(nut_ids), "| corners", corners)


if __name__ == "__main__":
    main()
