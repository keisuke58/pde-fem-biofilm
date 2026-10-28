"""implant_convergence.py -- mesh convergence table of the implant sector from results/implant_*.json
(written by make_implant_uel_inp.py --run).

    python implant_convergence.py [--bulge 0|1] [--fixed-seed] [--md]

Levels are nr x nt x nz = f (2 x 6 x 8), f = 1, 2, 3, ...; the element size is h = h_1 / f. By default the seed
is the inner element layer (make_implant_inp.py; it thins with the mesh, so the problem changes with the level);
--fixed-seed takes the runs with --seed-layers f (seed 0.125 mm thick at every level). For every
quantity q the table gives the value per level, the change from the level before, and the observed order
p from three consecutive levels (q = q_inf + C h^p, p solved from the ratio of the two differences; needs
monotone differences, else "-"), plus the extrapolated q_inf of the last three levels.

The "interior" of the .json leaves out one row of elements at each edge, a band that thins with the mesh; the
table adds the seed core, a fixed region: seed elements whose centroid has 0.25 <= z <= 1.75 and
15 <= theta <= 75 degrees (computed here from the per-element data of the .json).
"""
import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
QUANT = [("seed_alpha_minus_1", "seed alpha - 1"), ("seed_vM_over_E", "seed vM/E"),
         ("seed_interior_alpha_minus_1", "interior alpha - 1"), ("seed_interior_vM_over_E", "interior vM/E"),
         ("core_alpha_minus_1", "core alpha - 1"), ("core_vM_over_E", "core vM/E"),
         ("max_vM_over_E", "max vM/E (all)"), ("max_interior_vM_over_E", "max vM/E (interior)")]
CORE_Z, CORE_TH = (0.25, 1.75), (15.0, 75.0)


def load(bulge, fixed_seed=False):
    rows = []
    for f in glob.glob(str(HERE / "results" / "implant_*.json")):
        m = re.match(r"implant_(\d+)x(\d+)x(\d+)(_b1)?(_sl\d+)?\.json$", Path(f).name)
        if not m or bool(m.group(4)) != bool(bulge):
            continue
        r = json.load(open(f)); r["f"] = int(m.group(1)) // 2
        r.update(core(r))
        if r.get("seed_layers", 1) != (r["f"] if fixed_seed else 1):
            continue
        rows.append(r)
    return sorted(rows, key=lambda r: r["f"])


def core(r):
    """means of alpha - 1 and vM/E over the seed core (fixed region, see the header)"""
    el = r["elements"]
    th = np.degrees(np.arctan2(el["cy"], el["cx"])); z = np.array(el["cz"])
    m = (np.array(el["layer"]) < r.get("seed_layers", 1)) & (z >= CORE_Z[0]) & (z <= CORE_Z[1]) & (th >= CORE_TH[0]) & (th <= CORE_TH[1])
    if not m.any():
        return dict(core_elements=0, core_alpha_minus_1=None, core_vM_over_E=None)
    return dict(core_elements=int(m.sum()), core_alpha_minus_1=float(np.array(el["alpha_minus_1"])[m].mean()),
                core_vM_over_E=float(np.array(el["vM_over_E"])[m].mean()))


def order(f1, f2, f3, q1, q2, q3):
    """p of q = q_inf + C f^-p from three levels; None unless the differences shrink with the same sign"""
    d12, d23 = q1 - q2, q2 - q3
    if d12 == 0 or d23 == 0 or d12 * d23 < 0 or abs(d23) >= abs(d12):
        return None, None
    lo, hi = 1e-3, 20.0
    g = lambda p: (f1 ** -p - f2 ** -p) / (f2 ** -p - f3 ** -p) - d12 / d23  # noqa: E731
    if g(lo) * g(hi) > 0:
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if g(lo) * g(mid) > 0 else (lo, mid)
    p = 0.5 * (lo + hi)
    qinf = q3 - d23 * f3 ** -p / (f2 ** -p - f3 ** -p)
    return p, qinf


def table(rows, md=False):
    out = []
    head = ["level", "elements", "time"] + [n for _, n in QUANT] + ["max at (r, theta, z)"]
    out.append(("| " + " | ".join(head) + " |") if md else "  ".join(head))
    if md:
        out.append("|" + "---|" * len(head))
    for r in rows:
        v = [f"{r['nr']}x{r['nt']}x{r['nz']}", str(r["elements_total"]), f"{r.get('wall_time_s', float('nan')):.0f} s"]
        v += [f"{r[k]:.5g}" if r.get(k) is not None else "-" for k, _ in QUANT]
        v += ["(%.3f, %.1f, %.3f)" % tuple(r["max_centroid_r_theta_z"])]
        out.append(("| " + " | ".join(v) + " |") if md else "  ".join(v))
    if len(rows) >= 2:
        out.append("")
        out.append("relative change from the level before, and observed order p (q_inf) from three levels:")
        for k, n in QUANT:
            q = [r.get(k) for r in rows]
            if any(x is None for x in q):
                continue
            ch = ["-"] + [f"{(q[i] - q[i - 1]) / q[i - 1]:+.2%}" for i in range(1, len(q))]
            ps = []
            for i in range(2, len(q)):
                p, qi = order(rows[i - 2]["f"], rows[i - 1]["f"], rows[i]["f"], q[i - 2], q[i - 1], q[i])
                ps.append("p = -" if p is None else f"p = {p:.2f} (q_inf {qi:.5g})")
            out.append(f"  {n:22s} " + ", ".join(ch) + ("; " + "; ".join(ps) if ps else ""))
    return "\n".join(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulge", type=int, default=0)
    ap.add_argument("--fixed-seed", action="store_true", help="runs with --seed-layers f (seed 0.125 mm at every level)")
    ap.add_argument("--md", action="store_true", help="Markdown table")
    a = ap.parse_args()
    rows = load(a.bulge, a.fixed_seed)
    print(f"implant sector, bulge {a.bulge}, seed = {'0.125 mm' if a.fixed_seed else 'one element layer'}, {len(rows)} levels")
    print(table(rows, a.md))
