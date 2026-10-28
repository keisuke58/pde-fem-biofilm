"""A-84: grid convergence (32^3 / 64^3 / 128^3, Python TVD, minmod) of the seed alpha - 1 and von Mises for every
example of the thesis table, with the rule for newly reached points given as the argument (startpos: by position,
the A-75 runs; inhliv: the neighbours' living composition, A-82).

    python gci_examples_a84.py [startpos|inhliv]
"""
import glob
import json
import math
import sys
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results"
RULE = sys.argv[1] if len(sys.argv) > 1 else "startpos"
CASES = [("separate, antibiotic", "2026-10-3d", "", "separate_ab"),
         ("separate, no antibiotic", "2026-10-3d", "", "separate"),
         ("mixed (equal), antibiotic", "2026-10-3d", "_mixeq", "mixed_ab"),
         ("mixed (equal), no antibiotic", "2026-10-3d", "_mixeq", "mixed"),
         ("four species case 4, separate", "2026-10-3d-case3", "", "separate_ab"),
         ("four species case 4, mixed", "2026-10-3d-case3", "", "mixed_ab")]


def get(d, mix, ex, n, k26):
    pre = "" if n == 32 else f"nf{n}_"
    name = (f"coupled_living_consumption_{pre}front1_phimin0.001_v2*ansysbc_{RULE}{mix}_a100_gab1_gc1_tab0.333333"
            f"{'_k26c4' if k26 else ''}_tvd_wgmirror_fast_{ex}.json")
    f = sorted(glob.glob(str(R / d / name)))
    return list(json.load(open(f[0]))["runs"].values())[0] if f else None


def rich(v):
    f3, f2, f1 = v
    e32, e21 = f2 - f3, f1 - f2
    p = math.log(abs(e32 / e21)) / math.log(2) if e21 and e32 / e21 > 0 else float("nan")
    pp = p if p == p and p > 0 else 1.0
    return p, f1 + (f1 - f2) / (2 ** pp - 1), 125 * abs((f1 - f2) / f1) / (2 ** pp - 1)


print(f"rule: {RULE}")
print("| example | quantity | 32^3 | 64^3 | 128^3 | order | extrapolated | GCI |")
print("|---|---|---|---|---|---|---|---|")
for label, d, mix, ex in CASES:
    rs = [get(d, mix, ex, n, "four" in label) for n in (32, 64, 128)]
    if any(r is None for r in rs):
        print(f"| {label} | (missing: {[n for n, r in zip((32, 64, 128), rs) if r is None]}) |")
        continue
    for q, lab in (("seed_am1", "alpha-1"), ("seed_vm_Pa", "vM [Pa]")):
        v = [r[q] for r in rs]
        p, ext, g = rich(v)
        print(f"| {label} | {lab} | {v[0]:.3e} | {v[1]:.3e} | {v[2]:.3e} | {p:.2f} | {ext:.3e} | {g:.1f} % |")
