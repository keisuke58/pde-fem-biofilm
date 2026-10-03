#!/usr/bin/env python3
"""Klempt et al. 2024 Fig. 3 redrawn from this repository's reproduction:
the isosurface phi = 0.8 of test case 4.1 at 5, 25, 45 and 100 % of the
simulation time, inside the 20 um cube, next to the paper's panels.

Surface orange, its cut by the cube's faces red, as in the paper.
Reproduction: Table 2 unchanged, growth on both faces, consumption g phi c,
time scale and nutrient source as given on the command line (defaults: the
best fit of klempt2024_corner_scan, see KLEMPT2024_REPRODUCTION.md).

    python JAXFEM/klempt2024_fig3.py [--source corner2] [--s 2]
        -> assets/fig_klempt2024_fig3.png
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from PIL import Image  # noqa: E402
from skimage.measure import marching_cubes  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

OUT = ROOT / "assets" / "fig_klempt2024_fig3.png"
PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
TIMES = [0.05, 0.25, 0.45, 1.00]
AZIM = 150
ORANGE, RED = "#e8a317", "#c0000a"


def source(name):
    """'edge' = the edge line x = y = L; 'cornerA' = the nodes of the A um
    cube at the corner (L, L, L)."""
    if name == "edge":
        return B.nutrient_mask("edge")
    a = float(name[len("corner"):])
    return (K.X >= K.L - a - 1e-9) & (K.Y >= K.L - a - 1e-9) & (K.Z >= K.L - a - 1e-9)


def fields(src, s, growth="abs"):
    K.R, K.BETA, K.K_A = 100 * s, 2 * s, 1e-3 * s
    phi, _, _ = K.setup("fig4_edge")
    mask = source(src)
    out, t = {}, 0.0
    for t_next in TIMES:
        _, phi, _ = B.run_setup(phi, mask, 1e8, "ic", "first_order", growth,
                                t_end=round(t_next - t, 6))
        out[t_next] = phi.copy()
        t = t_next
    return out


def draw(ax, phi, level=0.8):
    # pad with zeros so that the surface is closed where it meets the cube
    pad = np.pad(phi, 1, constant_values=0.0)
    v, f, _, _ = marching_cubes(pad, level)
    v = (v - 1) * K.H                                  # back to um
    v = np.clip(v, 0, K.L)
    tri = v[f]
    on_face = np.zeros(len(f), bool)                   # triangles lying on a cube face
    for d in range(3):
        for b in (0.0, K.L):
            on_face |= np.all(np.isclose(tri[:, :, d], b, atol=1e-6), axis=1)
    # simple Lambert shading from one light
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    light = np.array([0.4, -0.6, 0.7]); light /= np.linalg.norm(light)
    shade = 0.45 + 0.55 * np.abs(n @ light)
    base = np.where(on_face[:, None], np.array(matplotlib.colors.to_rgb(RED)),
                    np.array(matplotlib.colors.to_rgb(ORANGE)))
    ax.add_collection3d(Poly3DCollection(tri, facecolors=base * shade[:, None], edgecolor="none"))
    for a in (0, K.L):
        for b in (0, K.L):
            ax.plot([0, K.L], [a, a], [b, b], color="0.6", lw=0.6)
            ax.plot([a, a], [0, K.L], [b, b], color="0.6", lw=0.6)
            ax.plot([a, a], [b, b], [0, K.L], color="0.6", lw=0.6)
    ax.set_xlim(0, K.L); ax.set_ylim(0, K.L); ax.set_zlim(0, K.L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=22, azim=AZIM)       # nutrient corner (L, L, L) at the top left


def paper_panels():
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "10", "-l", "10", "-r", "200", "-png", str(PDF),
                        str(Path(tmp) / "p")], check=True)
        im = Image.open(next(Path(tmp).glob("p*.png"))).convert("RGB")
    w, h = im.size
    # the four panels of Fig. 3 (fractions of the page, read off the render)
    boxes = [(0.305, 0.07, 0.62, 0.315), (0.625, 0.07, 0.94, 0.315),
             (0.305, 0.345, 0.62, 0.59), (0.625, 0.345, 0.94, 0.59)]
    return [np.asarray(im.crop((int(a * w), int(b * h), int(c * w), int(d * h)))) for a, b, c, d in boxes]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="corner2")
    ap.add_argument("--s", type=float, default=2.0)
    ap.add_argument("--growth", default="abs", help='"abs" or e.g. "blend0.5" (diagnostic)')
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    figstyle.apply(size=11)
    F = fields(a.source, a.s, a.growth)
    P = paper_panels()
    fig = plt.figure(figsize=(14, 7.2))
    for i, t in enumerate(TIMES):
        ax = fig.add_subplot(2, 4, i + 1, projection="3d")
        draw(ax, F[t])
        ax.set_title(f"this work, {round(t * 100)} %", fontsize=11)
        ax2 = fig.add_subplot(2, 4, i + 5)
        ax2.imshow(P[i]); ax2.set_axis_off()
        ax2.set_title(f"paper, {round(t * 100)} %", fontsize=11)
    fig.suptitle(r"Klempt et al. 2024, Fig. 3: isosurface $\phi = 0.8$ of test case 4.1. "
                 f"This work: Python reproduction (Table 2, nutrient {a.source}, time scale {a.s:g}, growth {a.growth}); "
                 "paper: Fig. 3 (CC BY 4.0)", fontsize=11.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.02, wspace=0.02, hspace=0.12)
    fig.savefig(a.out, dpi=200)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
