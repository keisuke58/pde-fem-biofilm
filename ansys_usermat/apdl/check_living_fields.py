#!/usr/bin/env python3
"""check_living_fields.py -- checks of the fixed ANSYS coupling runs (v2f_*, fsens_*, A-22) that need only the pushed
JSON, complementing check_living_trace.py (which checks the point model alone).

    python ansys_usermat/apdl/check_living_fields.py v2f_8_mix_ab [v2f_8_sep_ab ...] [--results DIR]

Per run (fields_<run>.json, <run>_pm.json):
  1. bounds: phi in [0, 1], c and a in [0, 1], L in [0, 1], share of A in [0, 1] at every element and time;
  2. symmetry: the examples are symmetric under x -> -x and z -> -z; largest difference of each field between
     mirrored elements (should be at solver tolerance; a larger value means an indexing error);
  3. antibiotic: a = 0 at every element before T* = 0.5;
  4. hand-over of the fields to the point model (pm trace): the inputs of a step continue those of the step before,
     c0(k+1) = c1(k) and a0(k+1) = a1(k), and the start state is the end state of the step before (consecutive
     steps only; the trace skips calls while phi is below its threshold; such a gap in v2f_8_sep_ab, element 25,
     (2, 4), was bridged by replaying the missing step with the Fortran point model: 8e-13);
  5. L passed back equals sum_i phi_i psi_i of the end state;
  6. pairs with and without antibiotic (run and run without "_ab"): identical up to T* = 0.5.
Exit 1 if a check fails.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL_SYM, TOL_HAND, TOL_L = 1e-6, 1e-12, 1e-12


def load(res, run):
    if not (res / f"fields_{run}.json").exists():         # superseded generations are in results/.../archive (22 Oct)
        res = res / "archive"
    f = json.loads((res / f"fields_{run}.json").read_text())
    p = res / f"{run}_pm.json"
    return f, (json.loads(p.read_text()) if p.exists() else None)


def mirror_index(f):
    x, y, z = (np.round(np.array(f[k]), 9) for k in ("cx", "cy", "cz"))
    key = {(a, b, c): i for i, (a, b, c) in enumerate(zip(x, y, z))}
    out = {}
    for name, m in (("x", (-1, 1, 1)), ("z", (1, 1, -1))):
        idx = [key.get((round(m[0] * a, 9), b, round(m[2] * c, 9))) for a, b, c in zip(x, y, z)]
        out[name] = None if any(i is None for i in idx) else np.array(idx)
    return out


def check_run(res, run):
    f, pm = load(res, run)
    ok, lines = True, []
    t = np.array(f["t"])
    fields = ("phi", "nut", "ab", "L", "shareA", "alpha", "seqv")
    # 1. bounds
    for k in ("phi", "nut", "ab", "L", "shareA"):
        u = np.array(f[k], float)
        lo, hi = np.nanmin(u), np.nanmax(u)
        bad = lo < -1e-9 or hi > 1 + 1e-9
        ok &= not bad
        lines.append(f"  bounds {k:6s} [{lo:.4f}, {hi:.4f}] {'FAIL' if bad else 'ok'}")
    # 2. symmetry
    mi = mirror_index(f)
    for ax, idx in mi.items():
        if idx is None:
            lines.append(f"  symmetry {ax}: mesh not mirrored, skipped")
            continue
        worst = []
        for k in fields:
            u = np.array(f[k], float)
            scale = max(np.nanmax(np.abs(u)), 1e-30)
            worst.append((np.nanmax(np.abs(u - u[:, idx])) / scale, k))
        d, k = max(worst)
        ok &= d < TOL_SYM
        lines.append(f"  symmetry {ax}: largest relative difference {d:.1e} ({k}) {'ok' if d < TOL_SYM else 'FAIL'}")
    # 3. antibiotic before 0.5
    if run.endswith("_ab") or "_ab_" in run or run.startswith("fsens"):
        a = np.array(f["ab"], float)
        pre = a[t < 0.5 - 1e-9]
        m = np.abs(pre).max() if pre.size else 0.0
        ok &= m == 0.0
        lines.append(f"  antibiotic before T*=0.5: max |a| {m:.1e} {'ok' if m == 0 else 'FAIL'}")
    # 4./5. pm trace
    if pm is not None:
        C = {k: i for i, k in enumerate(pm["columns"])}
        last = {}
        for r in pm["seed_trace"]:
            if str(r[C["hit"]]).strip() == "0":
                last[(str(r[C["elem"]]).strip(), str(r[C["ip"]]).strip(), int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
        gp = sorted({k[:2] for k in last})
        dc = da = ds = dl = 0.0
        nsteps = nmiss = 0
        for g in gp:
            ks = sorted(k for k in last if k[:2] == g)
            nsteps += len(ks)
            for k0, k1 in zip(ks, ks[1:]):
                r0, r1 = last[k0], last[k1]
                # the trace skips calls while phi is below its threshold: compare consecutive steps only
                if not (k1[2:] == (k0[2], k0[3] + 1) or (k1[2] == k0[2] + 1 and k1[3] == 1)):
                    nmiss += 1
                    continue
                dc = max(dc, abs(float(r1[C["c0"]]) - float(r0[C["c1"]])))
                da = max(da, abs(float(r1[C["a0"]]) - float(r0[C["a1"]])))
                ds = max(ds, max(abs(float(r1[C[f"gold_{i}"]]) - float(r0[C[f"gnew_{i}"]])) for i in range(1, 13)))
            for k in ks:
                r = last[k]
                phis = [float(r[C[f"gnew_{i}"]]) for i in (1, 2)]
                psis = [float(r[C[f"gnew_{i}"]]) for i in (7, 8)]
                dl = max(dl, abs(float(r[C["L"]]) - sum(p * s for p, s in zip(phis, psis))))
        for name, v, tol in (("c0(k+1) - c1(k)", dc, TOL_HAND), ("a0(k+1) - a1(k)", da, TOL_HAND),
                             ("start state - end state before", ds, TOL_HAND), ("L - sum phi_i psi_i", dl, TOL_L)):
            ok &= v < tol
            lines.append(f"  hand-over {name:30s} {v:.1e} {'ok' if v < tol else 'FAIL'}")
        lines.append(f"  ({len(gp)} traced Gauss points, {nsteps} steps, {nmiss} gaps in the trace skipped;"
                     " species columns gnew_1,2 / gnew_7,8 assume 2 species)")
    else:
        lines.append("  no pm trace")
    return ok, lines, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--results", default=str(HERE / "results/2026-10-living"))
    a = ap.parse_args()
    res = Path(a.results)
    allok = True
    loaded = {}
    for run in a.runs:
        ok, lines, f = check_run(res, run)
        loaded[run] = f
        allok &= ok
        print(f"{run}: {'PASS' if ok else 'FAIL'}")
        print("\n".join(lines))
    # 6. pairs
    for run in a.runs:
        if run.endswith("_ab"):
            base = run[:-3]
            if not (res / f"fields_{base}.json").exists() and (res / "archive" / f"fields_{base}.json").exists():
            res = res / "archive"
        if not (res / f"fields_{base}.json").exists():
                continue
            f1, f0 = loaded[run], json.loads((res / f"fields_{base}.json").read_text())
            t = np.array(f1["t"])
            ii = [i for i, tt in enumerate(t) if tt <= 0.5 + 1e-9 and tt in f0["t"]]
            d = max((np.nanmax(np.abs(np.array(f1[k][i], float) - np.array(f0[k][f0["t"].index(t[i])], float)))
                     for k in ("phi", "nut", "L", "alpha", "seqv") for i in ii), default=np.nan)
            ok = d == 0.0
            allok &= ok
            print(f"pair {run} / {base}: largest difference up to T*=0.5 {d:.1e} {'ok' if ok else 'FAIL'}")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
