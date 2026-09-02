#!/usr/bin/env python3
"""3D view of the two-layer curved-shell run, from the solver's own listing.

Parses the MAPDL `.txt` listing produced by the v222 wrapper run --- nodal
coordinates, nodal displacements, and the per-element von Mises table --- and
draws the geometry that was actually solved.

Deliberately does not place SEQV in space. The listing carries no element
connectivity (`ELIST` was never printed), so there is no honest way to say
where a given element sits; the stress is reported as a distribution instead.
Everything drawn is a quantity the file actually contains at a coordinate the
file actually gives.

    python ansys_usermat/apdl/plot_cylinder_3d.py \
        ansys_usermat/apdl/growth_cylinder_wrapper_result.txt -o assets/
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np

_NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"


def parse(path):
    """Returns (coords Nx3, usum N, seqv M). Streamed: the file is ~25 MB."""
    coords, disp, seqv = {}, {}, []
    mode = None
    # MAPDL packs columns without separators when a value is wide, so split on
    # the numbers themselves rather than on whitespace.
    num = re.compile(_NUM)
    with open(path, errors="replace") as fh:
        for line in fh:
            if "NODAL SOLUTION PER NODE" in line:
                mode = "u"; continue
            if "LIST ALL SELECTED NODES" in line:
                mode = "x"; continue
            if "ELEM" in line and "SEQV" in line:
                mode = "s"; continue
            if "PRINT S " in line or "SVAR" in line:
                mode = None; continue          # stop before the full stress dump
            head = line.lstrip()
            if not head or not head[0].isdigit():
                continue
            # A data row is numbers and whitespace, nothing else. MAPDL's page
            # furniture starts with a digit too -- the version banner
            # ("00000000  VERSION=WINDOWS x64 ... CP= 8.031") and the ESEL
            # summary ("11040 ELEMENTS (OF 12240 DEFINED)...") both parse as
            # plausible rows otherwise, and did: they put a radius of 65 and a
            # von Mises of 1.2e4 into the first version of this figure.
            if num.sub("", line).strip():
                continue
            v = [float(t.replace("D", "E")) for t in num.findall(line)]
            if mode == "x" and len(v) >= 4:
                coords[int(v[0])] = v[1:4]
            elif mode == "u" and len(v) >= 5:
                disp[int(v[0])] = v[4]             # USUM
            elif mode == "s" and len(v) == 2:
                seqv.append(v[1])
    ids = sorted(set(coords) & set(disp))
    if not ids:
        raise SystemExit("no nodes with both coordinates and displacement found")
    return (np.array([coords[i] for i in ids]),
            np.array([disp[i] for i in ids]),
            np.array(seqv))


def plot(xyz, usum, seqv, out, r_interface=4.3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = np.hypot(xyz[:, 0], xyz[:, 1])
    # The interface is at R_MID, not the midpoint of the radius range: the
    # substrate (4.0-4.3) is three times the growth layer's thickness.
    inner, outer = r < r_interface, r >= r_interface

    fig = plt.figure(figsize=(13.5, 4.3))

    # -- 1. the two layers, said plainly -------------------------------
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    ax.scatter(*xyz[inner].T, c="0.62", s=1.2, linewidths=0,
               label=f"substrate, $\\alpha=0$  ({inner.sum()} nodes)")
    ax.scatter(*xyz[outer].T, c="#c0392b", s=1.2, linewidths=0,
               label=f"growth layer  ({outer.sum()} nodes)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("two bonded layers", fontsize=10)
    ax.set_box_aspect((np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]), np.ptp(xyz[:, 2])))
    ax.view_init(elev=24, azim=-118)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.02), fontsize=7,
              frameon=False, markerscale=6)
    ax.tick_params(labelsize=7, pad=0)

    # -- 2. a radial section, where the interface is unmistakable ------
    ax1 = fig.add_subplot(1, 3, 2)
    th = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    sl = np.abs(xyz[:, 2] - xyz[:, 2].mean()) < 0.06      # a slice at mid-length
    p = ax1.scatter(th[sl], r[sl], c=usum[sl], s=5, cmap="viridis", linewidths=0)
    ax1.axhline(r_interface, color="#c0392b", lw=1.2, ls="--")
    ax1.text(th[sl].min(), r_interface, " interface", color="#c0392b",
             va="bottom", fontsize=8)
    ax1.set_xlabel("angle around the arc [deg]")
    ax1.set_ylabel("radius")
    ax1.set_title("section at mid-length, coloured by displacement", fontsize=10)
    ax1.text(0.5, 0.5, "two node rows through each layer\n= two elements through thickness",
             transform=ax1.transAxes, ha="center", fontsize=8, style="italic",
             color="0.4")
    fig.colorbar(p, ax=ax1, pad=0.02, label="USUM")

    # -- 3. the stress distribution ------------------------------------
    ax2 = fig.add_subplot(1, 3, 3)
    ax2.hist(seqv, bins=60, color="0.35")
    ax2.set_xlabel("element von Mises stress")
    ax2.set_ylabel("elements")
    ax2.set_title(f"SEQV over {len(seqv)} elements", fontsize=10)
    ax2.text(0.97, 0.95,
             f"min {seqv.min():.4g}\nmean {seqv.mean():.4g}\nmax {seqv.max():.4g}",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8,
             family="monospace")
    ax2.text(0.5, 0.55,
             "bimodal: the substrate does not grow,\nso its elements sit near zero",
             transform=ax2.transAxes, ha="center", fontsize=8, style="italic",
             color="0.35")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("listing", type=Path)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/cylinder_wrapper_3d.png"))
    ap.add_argument("--r-interface", type=float, default=4.3,
                    help="substrate/growth-layer radius (R_MID in the deck)")
    a = ap.parse_args(argv)
    xyz, usum, seqv = parse(a.listing)
    print(f"parsed {len(xyz)} nodes, {len(seqv)} element SEQV values")
    print(f"  radius {np.hypot(xyz[:,0],xyz[:,1]).min():.3f}"
          f" .. {np.hypot(xyz[:,0],xyz[:,1]).max():.3f}")
    print(f"  USUM   {usum.min():.4g} .. {usum.max():.4g}")
    print(f"  SEQV   {seqv.min():.4g} .. {seqv.max():.4g}, mean {seqv.mean():.6g}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(xyz, usum, seqv, a.out, a.r_interface)
    return 0


if __name__ == "__main__":
    sys.exit(main())
