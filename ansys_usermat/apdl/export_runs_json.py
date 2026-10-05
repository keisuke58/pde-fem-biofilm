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
        (outdir / f"{run}.json").write_text(json.dumps(rec, separators=(",", ":")))
        print(f"{run}: {', '.join(k for k in ('all_stress', 'elem_stress', 'nut_field') if k in rec)}"
              + (f"  (missing {rec['missing']})" if rec["missing"] else ""))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1:])
