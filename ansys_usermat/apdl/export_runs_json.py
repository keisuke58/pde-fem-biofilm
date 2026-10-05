#!/usr/bin/env python3
"""export_runs_json.py -- the CSV outputs of wired ANSYS runs as one JSON per run,
so that a session without F: (the cloud) can redraw figures from them.

    python ansys_usermat/apdl/export_runs_json.py OUTDIR run1 run2 ...

For each run name <run> it reads from F:\\biofilm_upf_wired (CSV files are
git-ignored):
  all_stress_<run>.csv   every element at the last result set: elem, seqv, sx,
                         sy, sz, alpha (SVAR 84), centroid cx, cy, cz [MPa, mm]
  elem_stress_<run>.csv  one element over time (the deck's --post-elem)
  nut_field_<run>.csv    if present: elem, Nut1 (SVAR 51), point-model phi_1,
                         phi_2 (SVAR 72, 73), alpha, centroid
  comp_trace_<run>.csv   if present (point-model runs): reduced to "share_history",
                         one row per substep: time T*, and for the seed (BIOFILM1)
                         and layers 1 and 2 around it (summarize_runs_json.layers)
                         the mean, min and max of phi_1/(phi_1+phi_2), the mean
                         amount phi_1+phi_2 and the number of points where the point
                         model ran (phi_3D >= 0.01). Traced points are integration
                         point 1 of each element the fragment traces. Nut1 is not in
                         the trace, so it is only in nut_field (last step).
and the deck's first line (how the deck was generated) plus the element
numbers of the seed component BIOFILM1. Values are rounded to 7 significant
digits. Missing files are skipped and listed under "missing".
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

W = Path(r"F:\biofilm_upf_wired")


def table(path: Path):
    rows = [r for r in csv.reader(open(path)) if r]
    head = [h.strip() for h in rows[0]]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        for h, x in zip(head, r):
            cols[h].append(float(f"{float(x):.7g}"))
    return cols


def seed_ids(deck: Path):
    lines = deck.read_text(encoding="latin-1").splitlines()
    i = next(k for k, l in enumerate(lines) if l.upper().startswith("CMBLOCK,BIOFILM1,"))
    n = int(lines[i].split(",")[3])
    raw, j = [], i + 2
    while len(raw) < n:
        raw += [int(x) for x in lines[j].split()]
        j += 1
    out: list[int] = []
    for v in raw:                                   # a negative entry means "up to"
        out += list(range(out[-1] + 1, -v + 1)) if v < 0 else [v]
    return out


def share_history(trace: Path, all_stress: dict, seed: list[int], phi_min=0.01):
    import numpy as np
    from summarize_runs_json import layers
    elem = np.array(all_stress["elem"], int)
    cen = np.stack([all_stress["cx"], all_stress["cy"], all_stress["cz"]], 1)
    lay = dict(zip(elem.tolist(), layers(cen, np.isin(elem, seed))[0].tolist()))
    last = {}                                       # (step, substep, elem) -> last call
    for r in csv.reader(open(trace)):
        if not r or not r[0].strip().isdigit() or r[1].strip() != "1":
            continue
        last[(int(r[2]), int(r[3]), int(r[0]))] = (float(r[6]), float(r[8]), float(r[24]), float(r[25]))
    dt = {k[:2]: v[0] for k, v in last.items()}
    steps = sorted(dt)
    t, dts = 0.0, {}
    for k in steps:                                 # T* at the end of each substep
        t += dt[k]
        dts[k] = t
    by_step = {k: [] for k in steps}
    for (ls, ss, e), (_, phi, p1, p2) in last.items():
        by_step[(ls, ss)].append((lay.get(e, -1), phi, p1, p2))
    # executables before 5 Oct evening traced only phi >= 0.5 (TRACE_PHI_MIN):
    # their seed/layer means are over the dense points only
    low = any(0.01 <= v[1] < 0.5 and e % 37 != 1 for (_, _, e), v in last.items())
    cols = {"points": "phi >= 0.01" if low else "phi >= 0.5 only (old trace rule): means biased to dense points",
            "t": []}
    for g in ("seed", "layer1", "layer2"):
        for q in ("mean", "min", "max", "amount", "n"):
            cols[f"{g}_{q}"] = []
    for k in steps:
        cols["t"].append(float(f"{dts[k]:.7g}"))
        rows = np.array(by_step[k]) if by_step[k] else np.zeros((0, 4))
        for g, L in (("seed", 0), ("layer1", 1), ("layer2", 2)):
            m = (rows[:, 0] == L) & (rows[:, 1] >= phi_min) if len(rows) else np.zeros(0, bool)
            s = rows[m, 2] + rows[m, 3]
            sh = rows[m, 2] / np.where(s > 0, s, 1)
            for q, v in (("mean", sh.mean() if m.any() else None), ("min", sh.min() if m.any() else None),
                         ("max", sh.max() if m.any() else None), ("amount", s.mean() if m.any() else None)):
                cols[f"{g}_{q}"].append(None if v is None else float(f"{v:.7g}"))
            cols[f"{g}_n"].append(int(m.sum()))
    return cols


def main(argv):
    outdir, runs = Path(argv[0]), argv[1:]
    outdir.mkdir(parents=True, exist_ok=True)
    for run in runs:
        rec = {"run": run, "missing": []}
        deck = W / f"{run}.dat"
        if deck.exists():
            rec["deck_header"] = deck.read_text(encoding="latin-1").splitlines()[0]
            rec["seed_BIOFILM1"] = seed_ids(deck)
        for key, name in (("all_stress", f"all_stress_{run}.csv"),
                          ("elem_stress", f"elem_stress_{run}.csv"),
                          ("nut_field", f"nut_field_{run}.csv")):
            p = W / name
            if p.exists():
                rec[key] = table(p)
            else:
                rec["missing"].append(name)
        tr = W / f"comp_trace_{run}.csv"
        if tr.exists() and "cx" in rec.get("all_stress", {}) and rec.get("seed_BIOFILM1"):
            rec["share_history"] = share_history(tr, rec["all_stress"], rec["seed_BIOFILM1"])
        (outdir / f"{run}.json").write_text(json.dumps(rec, separators=(",", ":")))
        print(f"{run}: {', '.join(k for k in ('all_stress', 'elem_stress', 'nut_field', 'share_history') if k in rec)}"
              + (f"  (missing {rec['missing']})" if rec["missing"] else ""))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1:])
