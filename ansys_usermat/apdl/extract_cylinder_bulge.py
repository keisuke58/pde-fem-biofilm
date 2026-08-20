#!/usr/bin/env python3
"""extract_cylinder_bulge.py --- pull the outer-surface radial displacement
out of t_growth_cylinder_shell.dat's growth_cylinder_result.txt (PRNSOL,U +
NLIST output) and plot it against the arc angle.

The full result file is large (~24 MB for this mesh) and not committed;
regenerate it with:
    ansys_usermat/apdl/run_apdl.ps1 -Deck t_growth_cylinder_shell.dat
which writes it into the F:-based working directory (see that script's
own header). Point this script at wherever it landed.

Usage:
    python extract_cylinder_bulge.py <path-to-growth_cylinder_result.txt> [out.png]
"""
import csv
import math
import re
import sys
from pathlib import Path

LEN = 3.0
R_OUT = 4.4


def parse(text: str) -> tuple[dict[int, tuple[float, float, float]], dict[int, tuple[float, float, float]]]:
    disp: dict[int, tuple[float, float, float]] = {}
    coord: dict[int, tuple[float, float, float]] = {}
    in_u_block = in_n_block = False
    for line in text.splitlines():
        if "POST1 NODAL DEGREE OF FREEDOM LISTING" in line:
            in_u_block, in_n_block = True, False
            continue
        if "LIST ALL SELECTED NODES" in line:
            in_u_block, in_n_block = False, True
            continue
        if "PRINT SVAR" in line or "STORE SEQV" in line:
            in_n_block = False
            continue

        if in_u_block:
            m = re.match(r"^\s*(\d+)\s+(.*)$", line)
            if m and re.search(r"[.\dE]", m.group(2)):
                nums = re.findall(r"-?\d+\.\d+(?:E[+-]\d+)?", m.group(2))
                if len(nums) == 4:
                    ux, uy, uz, _usum = map(float, nums)
                    disp[int(m.group(1))] = (ux, uy, uz)
        elif in_n_block:
            parts = line.split()
            if len(parts) == 4 and parts[0].isdigit():
                x, y, z = map(float, parts[1:])
                coord[int(parts[0])] = (x, y, z)
    return disp, coord


def outer_surface_slice(disp, coord, r_out=R_OUT, z_mid=LEN / 2,
                         r_tol=0.005, z_tol=0.03):
    """(theta_deg, u_radial) at the outer growth-layer surface, mid-length."""
    rows = []
    for nid in set(disp) & set(coord):
        x, y, z = coord[nid]
        r = math.hypot(x, y)
        if abs(r - r_out) < r_tol and abs(z - z_mid) < z_tol:
            ux, uy, _uz = disp[nid]
            u_r = (ux * x + uy * y) / r
            rows.append((math.degrees(math.atan2(y, x)), u_r))
    rows.sort()
    return rows


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    result_path = Path(sys.argv[1])
    out_png = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("cylinder_bulge.png")

    disp, coord = parse(result_path.read_text(encoding="utf-8", errors="replace"))
    rows = outer_surface_slice(disp, coord)
    if not rows:
        print("No nodes matched the outer-surface/mid-length slice -- geometry constants "
              "(R_OUT/LEN) at the top of this script may not match the deck.")
        sys.exit(1)

    csv_path = out_png.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta_deg", "u_radial"])
        w.writerows(rows)
    print(f"{len(rows)} points -> {csv_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available -- CSV written, skipping the plot.")
        return

    theta = [r[0] for r in rows]
    u_r = [r[1] * 1000 for r in rows]

    fig, ax = plt.subplots(figsize=(6.0, 4.2), dpi=200)
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.plot(theta, u_r, color="#2a78d6", linewidth=2.0, zorder=3)
    ax.scatter(theta, u_r, s=14, color="#2a78d6", zorder=4, edgecolor="#fcfcfb", linewidth=0.5)
    ax.set_xlabel(r"$\theta$ (deg, around the arc)", fontsize=10)
    ax.set_ylabel(r"outward radial displacement $u_r$ ($\times10^{-3}$)", fontsize=10)
    ax.grid(True, color="#e1e0d9", linewidth=0.7, zorder=0)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, facecolor="#fcfcfb", bbox_inches="tight")
    print(f"plot -> {out_png}")


if __name__ == "__main__":
    main()
