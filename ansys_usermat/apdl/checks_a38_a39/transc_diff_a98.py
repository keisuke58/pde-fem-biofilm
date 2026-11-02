"""A-98: 192^3 without antibiotic, with dc/dt (--transient-c 1e8, file tag _transc1e+08) against the quasi-static limit
(the thesis values), seed averages at T* = 1. Prints a Markdown table.

    python transc_diff_a98.py
"""
import glob
import json
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results" / "2026-10-3d"
print("| example | quantity | limit | with dc/dt | difference |")
print("|---|---|---|---|---|")
for label, ex, mix in (("separate", "separate", ""), ("mixed (equal)", "mixed", "_mixeq")):
    stem = f"coupled_living_consumption_nf192_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv{mix}_a100_gab1_gc1_tab0.333333"
    a = glob.glob(str(R / f"{stem}_tvd_wgmirror_fast_{ex}.json"))
    b = glob.glob(str(R / f"{stem}_transc1e+08_tvd_wgmirror_fast_{ex}.json"))
    if not (a and b):
        print(f"| {label} | missing (limit {len(a)}, dc/dt {len(b)}) | | | |")
        continue
    A = list(json.load(open(a[0]))["runs"].values())[0]
    B = list(json.load(open(b[0]))["runs"].values())[0]
    rows = [("alpha-1", A["seed_am1"], B["seed_am1"]), ("vM [Pa]", A["seed_vm_Pa"], B["seed_vm_Pa"])]
    rows += [(k, A[k][-1], B[k][-1]) for k in ("living_seed", "share_a_seed", "c_seed")]
    rows += [("phi_3D max (domain, all t)", max(A["phi_max"]), max(B["phi_max"]))]
    for q, x, y in rows:
        print(f"| {label} | {q} | {x:.6e} | {y:.6e} | {100 * (y / x - 1):+.3f} % |")
