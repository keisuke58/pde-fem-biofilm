"""A-83: the mixed example (equal amounts, antibiotic) with the minmod and the MC limiter of the TVD front term,
32^3 / 64^3 / 128^3: seed alpha - 1 and von Mises, and the amount of phi_3D the front term adds (--budget).

    python limiter_a83.py
"""
import glob
import json
import math
from pathlib import Path

D = Path(__file__).resolve().parents[1] / "results" / "2026-10-3d"


def run(n, lim, budget):
    pre = "" if n == 32 else f"nf{n}_"
    name = (f"coupled_living_consumption_{pre}front1_phimin0.001_v2_fortran_ansysbc_startpos{'_budget' if budget else ''}"
            f"_mixeq_a100_gab1_gc1_tab0.333333_tvd{lim}_wgmirror_fast_mixed_ab.json")
    f = glob.glob(str(D / name))
    return list(json.load(open(f[0]))["runs"].values())[0] if f else None


def rich(v):
    f3, f2, f1 = v
    e32, e21 = f2 - f3, f1 - f2
    p = math.log(abs(e32 / e21)) / math.log(2) if e21 and e32 / e21 > 0 else float("nan")
    pp = p if p == p and p > 0 else 1.0
    return p, f1 + (f1 - f2) / (2 ** pp - 1), 125 * abs((f1 - f2) / f1) / (2 ** pp - 1)


for lim, label in (("", "minmod"), ("_mc", "MC")):
    rs = [run(n, lim, lim == "_mc") or run(n, lim, False) for n in (32, 64, 128)]
    for q in ("seed_am1", "seed_vm_Pa"):
        v = [r[q] for r in rs]
        p, ext, g = rich(v)
        print(f"{label:6s} {q:10s} " + " ".join(f"{x:.3e}" for x in v) + f"  order {p:.2f}, extrapolated {ext:.3e}, GCI {g:.0f} %")
    amt = [(r.get("budget") or {}).get("phi_total_end") for r in rs]
    print(f"{label:6s} phi_3D total [mm^3]: {amt}")
