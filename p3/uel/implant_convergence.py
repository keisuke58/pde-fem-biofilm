"""implant_convergence.py -- mesh-convergence table for results/implant_<nr>x<nt>x<nz>[_b1].json
(make_implant_uel_inp.py --run).

    python implant_convergence.py [--bulge 0|1]

Per mesh: seed (inner layer) means of alpha - 1 and vM/E; the same over the interior of the seed layer, i.e. without
the top row (z = 2.0), the bottom row (z = 0) and the two columns at the symmetry planes; the largest vM/E over
all elements, its (r, theta, z) and the element's row (k from the top: 0 = the top row). The element order is
that of make_implant_uel_inp.sector(): e = i + nr (j + nt k).
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MESHES = [(2, 6, 8), (4, 12, 16), (6, 18, 24)]


def row(nr, nt, nz, bulge):
    f = HERE / "results" / f"implant_{nr}x{nt}x{nz}{'_b1' if bulge else ''}.json"
    if not f.exists():
        return None
    d = json.load(open(f))
    e = d["elements"]
    a, v = np.array(e["alpha_minus_1"]), np.array(e["vM_over_E"])
    n = np.arange(len(a))
    i, j, k = n % nr, (n // nr) % nt, n // (nr * nt)
    seed = i == 0
    inner = seed & (k > 0) & (k < nz - 1) & (j > 0) & (j < nt - 1)
    im = int(v.argmax())
    x, y, z = (np.array(e[c]) for c in ("cx", "cy", "cz"))
    return dict(mesh=f"{nr}x{nt}x{nz}", n=len(a), seed_a=a[seed].mean(), seed_v=v[seed].mean(), int_a=a[inner].mean(),
                int_v=v[inner].mean(), max_v=v[im], r=float(np.hypot(x[im], y[im])), th=float(np.degrees(np.arctan2(y[im], x[im]))),
                z=float(z[im]), layer=int(i[im]), k_from_top=int(nz - 1 - k[im]), j=int(j[im]),
                max_inner_wo_edges=float(v[(k > 0) & (k < nz - 1) & (j > 0) & (j < nt - 1)].max()))


if __name__ == "__main__":
    for bulge in (0, 1) if "--bulge" not in sys.argv else (int(sys.argv[sys.argv.index("--bulge") + 1]),):
        print(f"bulge {bulge}")
        print("mesh      elems  seed a-1   seed vM/E  int a-1    int vM/E   max vM/E  (r, theta, z)            layer row_from_top  max w/o edge rows")
        for m in MESHES:
            r = row(*m, bulge)
            if r:
                print(f"{r['mesh']:9s} {r['n']:5d}  {r['seed_a']:.5f}   {r['seed_v']:.5f}   {r['int_a']:.5f}   {r['int_v']:.5f}   {r['max_v']:.4f}   "
                      f"({r['r']:.3f}, {r['th']:.1f}, {r['z']:.3f})  {r['layer']:3d}   {r['k_from_top']:3d}   {r['max_inner_wo_edges']:.4f}")
