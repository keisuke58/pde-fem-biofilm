"""plot_3d.py -- 3D cut-away view of the stress around the seed, from all_stress.csv.

    python ansys_usermat/apdl/plot_3d.py all_stress_<job>.csv
           [--out assets/stress3d_<job>.png] [--title "..."]
    python ansys_usermat/apdl/plot_3d.py old_all_stress.csv --grid 8 --size 2.0

The element centroids come from the cx, cy, cz columns of
callsite/post_all_stress.mac. For CSVs written before those columns existed,
pass the deck the run used (--deck run.dat): the centroids are then computed
from its nblock/eblock. --grid N --size L instead ASSUMES a regular N^3 block
numbered x fastest, then y, then z; the partner's decks are not numbered that
way (5 Oct: 8^3 runs x reversed, then z, then y), so the picture comes out
mirrored and with axes swapped. Use --grid only without the deck.

The block x >= x_seed, y < y_seed is cut away so the seed and the elements
around it are visible. Left: von Mises in Pa on a log scale (the stress drops
by orders of magnitude away from the seed). Right: mean stress
p = (SX+SY+SZ)/3 in Pa, diverging scale centred on zero (blue compression,
red tension). The deck is in MPa.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def deck_centroids(deck):
    """element number -> centroid, from the deck's nblock and eblock (solid format)."""
    lines = Path(deck).read_text(encoding="latin-1").splitlines()
    nodes, cen, i = {}, {}, 0
    while i < len(lines):
        s = lines[i].strip().lower()
        if s.startswith("nblock"):
            i += 2
            while not lines[i].strip().startswith("-1") and not lines[i].strip().upper().startswith("N,"):
                f = lines[i].split()
                nodes[int(f[0])] = np.array([float(x) for x in f[-3:]])
                i += 1
        elif s.startswith("eblock"):
            i += 2
            while not lines[i].strip().startswith("-1"):
                f = [int(x) for x in lines[i].split()]
                cen[f[10]] = np.mean([nodes[n] for n in f[11:19]], axis=0)
                i += 1
        i += 1
    return cen


def read(path, grid=None, size=None, deck=None):
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r]
    head = [h.strip() for h in rows[0]]
    a = np.array([[float(x) for x in r] for r in rows[1:]])
    col = {h: a[:, i] for i, h in enumerate(head)}
    col["_centroids"] = "csv"
    if not {"cx", "cy", "cz"} <= set(head):
        if deck is not None:
            cen = deck_centroids(deck)
            c = np.array([cen[int(e)] for e in col["elem"]])
            col["cx"], col["cy"], col["cz"] = c[:, 0], c[:, 1], c[:, 2]
            col["_centroids"] = "deck"
        elif grid is not None:
            e = col["elem"].astype(int) - 1
            h = size / grid
            for k, name in enumerate(("cx", "cy", "cz")):
                col[name] = -size / 2 + h * ((e // grid ** k) % grid + 0.5)
            col["_centroids"] = "grid"
        else:
            raise SystemExit("no centroid columns: pass the run's deck (--deck), "
                             "or --grid/--size (assumed numbering)")
    col["p"] = (col["sx"] + col["sy"] + col["sz"]) / 3.0
    return col


def to_grid(col):
    c = np.stack([col["cx"], col["cy"], col["cz"]], axis=1)
    axes = [np.unique(np.round(c[:, k], 9)) for k in range(3)]
    idx = [np.searchsorted(axes[k], np.round(c[:, k], 9)) for k in range(3)]
    shape = tuple(len(x) for x in axes)
    if len(col["elem"]) != np.prod(shape):
        raise SystemExit(f"{len(col['elem'])} elements do not fill a "
                         f"{shape} block: this view needs a regular mesh")
    g = {}
    for k in ("seqv", "p", "alpha"):
        v = np.full(shape, np.nan)
        v[idx[0], idx[1], idx[2]] = col[k]
        g[k] = v
    return g, axes


def check(g):
    """Seed = largest alpha; its face neighbours should carry the most stress."""
    s = np.unravel_index(np.nanargmax(g["alpha"]), g["alpha"].shape)
    q = g["seqv"]
    nb = []
    for d in range(3):
        for sgn in (-1, 1):
            j = list(s)
            j[d] += sgn
            if 0 <= j[d] < q.shape[d]:
                nb.append(q[tuple(j)])
    top = np.sort(q.ravel())[::-1][:len(nb) + 1]
    ok = min(nb) >= top[-1]
    print(f"seed (largest alpha) at grid index {tuple(int(x) for x in s)}; "
          f"face neighbours von Mises {min(nb):.3g} .. {max(nb):.3g} "
          f"[MPa]; the {len(nb) + 1} largest are >= {top[-1]:.3g}: "
          f"{'neighbours are the most stressed, OK' if ok else 'CHECK FAILED'}")
    return s, ok


def main(argv=None):
    import matplotlib
    matplotlib.use("Agg")
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import figstyle
    figstyle.apply()
    import matplotlib.pyplot as plt
    from matplotlib import colors, cm

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv")
    ap.add_argument("--out")
    ap.add_argument("--title", default="")
    ap.add_argument("--grid", type=int)
    ap.add_argument("--size", type=float)
    ap.add_argument("--deck", help="the run's deck, for CSVs without centroid columns")
    a = ap.parse_args(argv)

    col = read(a.csv, a.grid, a.size, a.deck)
    g, axes = to_grid(col)
    if col["_centroids"] == "grid":
        s, _ = check(g)
    else:
        s = np.unravel_index(np.nanargmax(g["alpha"]), g["alpha"].shape)
    q, p = g["seqv"] * 1e6, g["p"] * 1e6          # MPa -> Pa
    n = q.shape
    edges = []
    for ax in axes:
        h = np.diff(ax).mean() if len(ax) > 1 else 1.0
        edges.append(np.concatenate([ax - h / 2, [ax[-1] + h / 2]]))
    X, Y, Z = np.meshgrid(*edges, indexing="ij")
    i, j, _k = np.indices(n)
    show = ~((i >= s[0]) & (j < s[1]))

    qmax = np.nanmax(q)
    qn = colors.LogNorm(vmin=qmax * 1e-3, vmax=qmax)
    pa = np.nanmax(np.abs(p))
    pn = colors.TwoSlopeNorm(vcenter=0.0, vmin=-pa, vmax=pa)
    panels = (("von Mises [Pa], log scale", np.clip(q, qmax * 1e-3, None), qn, "Blues"),
              (r"mean stress $p$ [Pa] ($-$ compression)", p, pn, "RdBu_r"))

    fig = plt.figure(figsize=(13, 6))
    for k, (title, v, norm, cmap) in enumerate(panels):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        fc = getattr(cm, cmap)(norm(v.ravel())).reshape(v.shape + (4,))
        fc[..., 3] = 1.0
        ax.voxels(X, Y, Z, show, facecolors=fc, edgecolor=(0, 0, 0, 0.15), linewidth=0.3)
        ax.set(xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]", title=title)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=22, azim=-60)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, pad=0.1)
    fig.suptitle((a.title + "  " if a.title else "")
                 + r"block $x \geq x_{seed}$, $y < y_{seed}$ cut away (seed face exposed)", fontsize=11)
    out = Path(a.out) if a.out else Path(a.csv).with_suffix(".3d.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
