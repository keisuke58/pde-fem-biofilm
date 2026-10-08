"""wp3_alpha_profile.py -- alpha - 1 of the WP3 seed-growth runs measured in a
way that does not change with the mesh (RUN_WP3_IKMHIWI03.md, block C: the
"surface" and "interior" means of compare_surface_fix.py are element layers,
so their region shrinks with h).

    python ansys_usermat/apdl/wp3_alpha_profile.py RUN.json [RUN.json ...]

Per run:
- seed mean: mean alpha - 1 over all seed elements (the seed is the same
  physical region on every mesh, refine_deck.py splits the elements);
- domain integral: sum of (alpha - 1) * element volume over the cube (mm^3);
- line: alpha - 1 at fixed points on the x axis through the seed centre,
  trilinear interpolation of the element-centroid values (a regular grid).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

OFFSETS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)   # mm from the seed centre along x


def profile(rec):
    a = rec["all_stress"]
    elem = np.array(a["elem"], int)
    c = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    al = np.array(a["alpha"])
    seed = np.isin(elem, rec.get("seed_BIOFILM1", []))
    axes = [np.unique(np.round(c[:, i], 9)) for i in range(3)]
    n = [len(x) for x in axes]
    grid = np.full(n, np.nan)
    idx = [np.searchsorted(axes[i], np.round(c[:, i], 9)) for i in range(3)]
    grid[idx[0], idx[1], idx[2]] = al
    h = axes[0][1] - axes[0][0]
    f = RegularGridInterpolator(axes, grid, bounds_error=False, fill_value=None)
    ctr = c[seed].mean(0)
    pts = np.array([[ctr[0] + d, ctr[1], ctr[2]] for d in OFFSETS])
    return {"n": n[0], "seed mean": al[seed].mean(), "integral": al.sum() * h ** 3,
            "line": f(pts)}


def main():
    print("run | n | seed mean | integral mm^3 | alpha - 1 at x - x_c = "
          + " / ".join(f"{d:g}" for d in OFFSETS) + " mm")
    for p in sys.argv[1:]:
        r = profile(json.loads(Path(p).read_text()))
        line = " / ".join(f"{v:.3e}" for v in r["line"])
        print(f"{Path(p).stem} | {r['n']} | {r['seed mean']:.4e} | {r['integral']:.4e} | {line}")


if __name__ == "__main__":
    main()
