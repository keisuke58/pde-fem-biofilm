"""plot_section.py -- section maps of growth and stress from all_stress.csv.

    python ansys_usermat/apdl/plot_section.py all_stress_<job>.csv [--axis x]
           [--out assets/section_<job>.png] [--title "..."]

Reads the CSV written by callsite/post_all_stress.mac (with the centroid
columns cx, cy, cz added 2026-10-02), cuts the mesh through the element with
the largest growth alpha (the seed), normal to --axis, and draws three panels
on that section: alpha, von Mises stress and mean stress p = (SX+SY+SZ)/3, the
stresses in Pa (the deck is in MPa). alpha and von Mises use a sequential
colour scale; p uses a diverging one centred on zero, so compression and
tension read at a glance. The seed is outlined.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read(path, grid=None, size=None):
    """CSV with centroid columns; without them, --grid/--size assume a
    regular block (see plot_3d.py, which prints a check of that)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from plot_3d import read as read3d
    return read3d(path, grid, size)


def section(col, axis):
    ax_i = "xyz".index(axis)
    c = np.stack([col["cx"], col["cy"], col["cz"]], axis=1)
    seed = int(np.argmax(col["alpha"]))
    tol = 1e-6 * (np.ptp(c[:, ax_i]) or 1.0)
    on = np.abs(c[:, ax_i] - c[seed, ax_i]) < tol
    u_i, v_i = [i for i in range(3) if i != ax_i]
    us = np.unique(np.round(c[on, u_i], 9))
    vs = np.unique(np.round(c[on, v_i], 9))
    iu = np.searchsorted(us, np.round(c[on, u_i], 9))
    iv = np.searchsorted(vs, np.round(c[on, v_i], 9))
    grids = {}
    for k in ("alpha", "seqv", "p"):
        g = np.full((len(vs), len(us)), np.nan)
        g[iv, iu] = col[k][on]
        grids[k] = g
    seedmask = np.zeros((len(vs), len(us)), bool)
    a_on = col["alpha"][on]
    seedmask[iv, iu] = a_on > 0.5 * a_on.max()
    return grids, us, vs, "xyz"[u_i], "xyz"[v_i], c[seed, ax_i], seedmask


def main(argv=None):
    import matplotlib
    matplotlib.use("Agg")
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import figstyle
    figstyle.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--axis", default="x", choices=("x", "y", "z"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default="")
    ap.add_argument("--grid", type=int, help="old CSV without centroids: N per side")
    ap.add_argument("--size", type=float, help="old CSV without centroids: side [mm]")
    a = ap.parse_args(argv)
    col = read(a.csv, a.grid, a.size)
    g, us, vs, un, vn, pos, seed = section(col, a.axis)
    pa = 1e6                                   # MPa -> Pa
    du = (us[1] - us[0]) if len(us) > 1 else 1.0
    dv = (vs[1] - vs[0]) if len(vs) > 1 else 1.0
    ext = [us[0] - du / 2, us[-1] + du / 2, vs[0] - dv / 2, vs[-1] + dv / 2]
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0), constrained_layout=True)
    panels = [("alpha", "growth $\\alpha - 1$ [-]", 1.0, "Blues", None),
              ("seqv", "von Mises [Pa]", pa, "Blues", None),
              ("p", "mean stress $p$ [Pa] (− compression)", pa, "RdBu_r", "div")]
    for ax, (k, lab, sc, cmap, kind) in zip(axs, panels):
        d = g[k] * sc
        if kind == "div":
            lim = np.nanmax(np.abs(d)) or 1.0
            im = ax.imshow(d, origin="lower", extent=ext, cmap=cmap,
                           norm=TwoSlopeNorm(0.0, -lim, lim))
        else:
            im = ax.imshow(d, origin="lower", extent=ext, cmap=cmap,
                           vmin=0.0)
        ax.contour(np.where(seed, 1.0, 0.0), levels=[0.5], colors="k",
                   linewidths=1.0, origin="lower", extent=ext)
        ax.grid(False)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(f"{un} [mm]")
        ax.set_ylabel(f"{vn} [mm]")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle((a.title + "  " if a.title else "")
                 + f"section {a.axis} = {pos:.3f} mm through the seed "
                 "(outlined)", fontsize=10)
    out = Path(a.out or Path(a.csv).with_suffix(".section.png"))
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
