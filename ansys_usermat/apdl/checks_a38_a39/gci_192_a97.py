"""A-97: four grids (32^3 / 64^3 / 128^3 / 192^3, Python TVD, outer step 0.0025, inherit rule) of the seed averages at
T* = 1, with and without antibiotic; GCI from the three finest grids (64/128/192, ratios 2 and 1.5, Celik et al. 2008
iteration, Fs = 1.25), as in A-92/A-94.

    python gci_192_a97.py [--transc 1e8]     (--transc: the runs with dc/dt, file tag _transc1e8)
"""
import glob
import json
import math
import sys
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results" / "2026-10-3d"
TC = f"_transc{sys.argv[sys.argv.index('--transc') + 1]}" if "--transc" in sys.argv else ""
CASES = [("separate, antibiotic", "", "separate_ab"), ("separate, no antibiotic", "", "separate"),
         ("mixed (equal), antibiotic", "_mixeq", "mixed_ab"), ("mixed (equal), no antibiotic", "_mixeq", "mixed")]


def get(mix, ex, n):
    pre = "" if n == 32 else f"nf{n}_"
    tc = TC if n == 192 else ""
    name = (f"coupled_living_consumption_{pre}front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv{mix}"
            f"_a100_gab1_gc1_tab0.333333{tc}_tvd_wgmirror_fast_{ex}.json")
    f = sorted(glob.glob(str(R / name)))
    return list(json.load(open(f[0]))["runs"].values())[0] if f else None


def gci(f1, f2, f3, r21=1.5, r32=2.0):
    """f1 finest. Celik et al. 2008: fixed-point iteration for the order p with unequal ratios."""
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0 or e32 / e21 <= 0:
        return float("nan"), float("nan"), float("nan")
    s = 1.0
    p = abs(math.log(abs(e32 / e21))) / math.log(r21)
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
    ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    return p, ext, 125 * abs((f1 - f2) / f1) / (r21 ** p - 1)


print(f"files at 192^3: {'with dc/dt ' + TC if TC else 'quasi-static limit'}")
print("| example | quantity | 32^3 | 64^3 | 128^3 | 192^3 | p | extrapolated | GCI |")
print("|---|---|---|---|---|---|---|---|---|")
for label, mix, ex in CASES:
    rs = [get(mix, ex, n) for n in (32, 64, 128, 192)]
    if any(r is None for r in rs):
        print(f"| {label} | missing {[n for n, r in zip((32, 64, 128, 192), rs) if r is None]} |")
        continue
    for q, lab, sc in (("seed_am1", "alpha-1 [e-4]", 1e4), ("seed_vm_Pa", "vM [e-4 Pa]", 1e4)):
        v = [r[q] * sc for r in rs]
        p, ext, g = gci(v[3], v[2], v[1])
        print(f"| {label} | {lab} | {v[0]:.3f} | {v[1]:.3f} | {v[2]:.3f} | {v[3]:.3f} | {p:.2f} | {ext:.3f} | {g:.1f} % |")
    for q, lab in (("living_seed", "living"), ("share_a_seed", "share A"), ("c_seed", "c")):
        print(f"| {label} | {lab} (T*=1) | " + " | ".join(f"{r[q][-1]:.4f}" for r in rs) + " | | | |")
    print(f"| {label} | phi_3D max (domain, end) | " + " | ".join(f"{r['phi_max'][-1]:.4f}" for r in rs) + " | | | |")
