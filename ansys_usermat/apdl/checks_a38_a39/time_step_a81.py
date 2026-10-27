"""A-81: time-step convergence of the main Python runs (64^3, TVD, separate example with antibiotic): outer step
(point model, c and a) Delta T* = 0.025 (main), 0.0125, 0.00625; observed order and the error of the main step against
the Richardson value. Also the grid convergence of the mixed example (32^3, 64^3, 128^3).

    python time_step_a81.py
"""
import glob
import json
import math
from pathlib import Path

D = Path(__file__).resolve().parents[1] / "results" / "2026-10-3d"
TAIL = "phimin0.001_v2_fortran_ansysbc_startpos{mix}_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast_{ex}.json"


def get(name):
    f = glob.glob(str(D / name))
    r = list(json.load(open(f[0]))["runs"].values())[0]
    return {"alpha-1": r["seed_am1"], "vM": r["seed_vm_Pa"], "sum phibar": r["living_seed"][-1]}


def richardson(v, ratio=2.0):
    f3, f2, f1 = v
    e32, e21 = f2 - f3, f1 - f2
    p = math.log(abs(e32 / e21)) / math.log(ratio) if e21 and e32 / e21 > 0 else float("nan")
    pp = p if p == p and p > 0 else 1.0
    ext = f1 + (f1 - f2) / (ratio ** pp - 1)
    return p, ext, 100 * 1.25 * abs((f1 - f2) / f1) / (ratio ** pp - 1)


sep = TAIL.format(mix="", ex="separate_ab")
runs = [get("coupled_living_consumption_nf64_front1_" + sep.replace("_fortran", "")),   # the main run (JAX point model, = Fortran to 1e-16)
        get("coupled_living_consumption_nf64_front1_dtc0.0125_" + sep),
        get("coupled_living_consumption_nf64_front1_dtc0.00625_" + sep)]
print("time step (separate, antibiotic, 64^3): Delta T* = 0.025 / 0.0125 / 0.00625")
for q in runs[0]:
    v = [r[q] for r in runs]
    p, ext, _ = richardson(v)
    print(f"  {q:10s} {v[0]:.4e} {v[1]:.4e} {v[2]:.4e}  order {p:.2f}, error of 0.025: {100 * abs(v[0] - ext) / abs(ext):.2f} %")
mix = TAIL.format(mix="_mixeq", ex="mixed_ab")
runs = [get("coupled_living_consumption_front1_" + mix), get("coupled_living_consumption_nf64_front1_" + mix),
        get("coupled_living_consumption_nf128_front1_" + mix)]
print("grid (mixed, equal amounts, antibiotic): 32^3 / 64^3 / 128^3")
for q in runs[0]:
    v = [r[q] for r in runs]
    p, ext, g = richardson(v)
    print(f"  {q:10s} {v[0]:.4e} {v[1]:.4e} {v[2]:.4e}  order {p:.2f}, extrapolated {ext:.4e}, GCI {g:.1f} %")
