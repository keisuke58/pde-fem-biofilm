"""judge_comp_trace.py -- judge a comp_trace.csv written by the composition
mode (prop(28) = 7) of the partner-element usermat.

    python ansys_usermat/apdl/judge_comp_trace.py F:\\biofilm_upf_wired\\comp_trace.csv \\
        [--case 2sp_case3] [--k-alpha 1e-3] [--phi-min 0.01] [--s S]

1. composition_reference.check_comp_trace (once per increment, carried,
   sum phi_i = phi_3D, one server call per point and substep, Eq. 36, replay);
2. every traced point's whole history against the stand-alone scheme
   composition_reference.reference driven by that point's own phi_3D(t).
Exit 1 if anything fails, so a script can stop on it.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "coupling"))
import composition_reference as cr  # noqa: E402
import material_server as ms  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trace")
    ap.add_argument("--case", default="2sp_case3")
    ap.add_argument("--k-alpha", type=float, default=1.0e-3)
    ap.add_argument("--phi-min", type=float, default=0.01)
    ap.add_argument("--age", action="store_true",
                    help="trace from prop(28) = 8 (age_trace.csv): independent point model, no scaling")
    ap.add_argument("--phi-init", type=float, nargs=2, default=(0.2, 0.2),
                    help="phi_1, phi_2 at the start, prop(30)/prop(32) (--age only)")
    ap.add_argument("--s", type=float, default=None,
                    help="point-model clock, prop(31); default: read from the trace (dt_pm / dtime)")
    a = ap.parse_args()

    ms.set_case(a.case)
    theta, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
    rows = cr.read_comp_trace(a.trace)
    if not rows:
        print("FAIL: empty trace")
        return 1
    if a.age:
        checks = cr.check_age_trace(rows, a.k_alpha, theta, hp)
    else:
        checks = cr.check_comp_trace(rows, a.k_alpha, theta, hp)
    print("checks:", checks)

    by_pt = defaultdict(dict)
    hits = defaultdict(list)
    for r in rows:
        key = (r["elem"], r["ip"])
        by_pt[key].setdefault((r["ldstep"], r["isubst"]), r)
        hits[key].append(r["hit"])
    worst = 0.0
    print(" elem ip steps  phi3 first..last   chi1 end  psi1    psi2    "
          "calls/cached/held  max diff")
    for key, subs in sorted(by_pt.items()):
        seq = [subs[k] for k in sorted(subs)]
        s_pt = a.s if a.s is not None else seq[0]["dt_pm"] / seq[0]["dtime"]
        if a.age:
            ref = cr.reference_age([r["phi3"] for r in seq], theta, hp,
                                   seq[0]["dtime"], phi_min=a.phi_min, s=s_pt,
                                   phi_init=tuple(a.phi_init))
        else:
            ref = cr.reference([r["phi3"] for r in seq], theta, hp,
                               seq[0]["dtime"], phi_min=a.phi_min, s=s_pt,
                               dt_pm_series=None if a.s is not None
                               else [r["dt_pm"] for r in seq])
        d = max(float(np.max(np.abs(np.asarray(r["g_new"]) - g)))
                for r, g in zip(seq, ref))
        worst = max(worst, d)
        g = seq[-1]["g_new"]
        s = g[0] + g[1]
        h = hits[key]
        print(f"{key[0]:5d} {key[1]:2d} {len(seq):5d}  {seq[0]['phi3']:.4f}..{seq[-1]['phi3']:.4f}"
              f"   {g[0] / s if s else float('nan'):.4f}    {g[6]:.4f}  {g[7]:.4f}  "
              f"{h.count(0)}/{h.count(1)}/{h.count(2)}   {d:.1e}")
    print("max |ANSYS - stand-alone scheme|:", worst)
    ok = all(v for k, v in checks.items() if k != "replay_worst") and worst == 0.0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
