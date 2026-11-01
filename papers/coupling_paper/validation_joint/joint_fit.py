#!/usr/bin/env python3
"""Joint fit of the volume comparison (Keio continuation, not the thesis model): one front speed and one lysis rate for
all four conditions of Heine et al. 2025, one detachment (V0, tau) for the two flow conditions (CH, DH), none for the
static ones (CS, DS). Grid search over the 32^3 layer runs of coupled_living_consumption.py (--lysis, --erosion,
--erosion-tau; papers/coupling_paper/VALIDATION_RESULTS.ja.md).

    python papers/coupling_paper/validation_joint/joint_fit.py run [--jobs 4]      # the missing runs
    python papers/coupling_paper/validation_joint/joint_fit.py fit --measured heine_volume_rel.json
The measured values stay outside this repository; fit prints the ranking and writes ranking.json (model errors only).
"""
import itertools
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOLVER = ROOT / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"
RES = ROOT / "ansys_usermat" / "apdl" / "results" / "2026-10-validation"
FRONTS = [0.1, 0.2]
LYSIS = [0, 1, 2, 4, 8]
EROS = [None, (4, 0.03), (8, 0.03), (8, 0.06), (12, 0.04)]
FLOW, STATIC = ["CH", "DH"], ["CS", "DS"]
DAYS = [1, 3, 6, 10, 15, 21]


def tag(front, cond, k, ero):
    t = f"coupled_living_consumption_front{front:g}_ecocal5_{cond}_dtc0.0025_phimin0.001_fortran"
    t += f"_lys{k:g}" if k else ""
    t += f"_ero{ero[0]:g}tau{ero[1]:g}" if ero else ""
    return RES / (t + "_gc1_transc1e+08_tvd_wgmirror_fast_layer.json")


def cases():
    for f, k in itertools.product(FRONTS, LYSIS):
        for c in STATIC:
            yield f, c, k, None
        for c, e in itertools.product(FLOW, EROS):
            yield f, c, k, e


def run_one(case):
    f, c, k, e = case
    if tag(*case).exists():
        return
    args = ["--nf", "32", "--fast", "--fortran", "--only", "layer", "--front", f"{f:g}", "--gc", "1", "--phi-min", "0.001",
            "--tvd", "--wall-grad", "mirror", "--transient-c", "1e8", "--dtc", "0.0025",
            "--eco-file", f"abaqus_composition/eco_cal5_{c}.txt"]
    if k:
        args += ["--lysis", f"{k:g}"]
    if e:
        args += ["--erosion", f"{e[0]:g}", "--erosion-tau", f"{e[1]:g}"]
    subprocess.run([sys.executable, str(SOLVER), *args], cwd=ROOT, capture_output=True)
    print("done", case, tag(*case).exists(), flush=True)


def rel(case):
    h = json.loads(tag(*case).read_text())["runs"]["layer"]
    v = np.interp([(d - 1) / 20 for d in DAYS], h["t"], h["phi_total_mm3"])
    return v / v[0]


def main():
    mode = sys.argv[1]
    if mode == "run":
        jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 4
        with ThreadPoolExecutor(jobs) as ex:
            list(ex.map(run_one, list(cases())))
        return
    meas = json.loads(Path(sys.argv[sys.argv.index("--measured") + 1]).read_text())
    m = {c: np.array([meas[c][str(d)]["rel"] for d in DAYS]) for c in FLOW + STATIC}
    err = lambda c, r: float(np.sqrt(np.mean(np.log(np.maximum(r[1:], 1e-6) / m[c][1:]) ** 2)))
    rows = []
    for f, k, e in itertools.product(FRONTS, LYSIS, EROS):
        per = {c: err(c, rel((f, c, k, None))) for c in STATIC}
        per.update({c: err(c, rel((f, c, k, e))) for c in FLOW})
        rows.append(dict(front=f, lysis=k, erosion=e, per=per, mean=float(np.mean(list(per.values())))))
    rows.sort(key=lambda r: r["mean"])
    for r in rows[:8]:
        print(f"front {r['front']:g} lysis {r['lysis']:g} erosion {r['erosion']}: mean {r['mean']:.2f} ",
              " ".join(f"{c} {v:.2f}" for c, v in r["per"].items()))
    (HERE / "ranking.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
