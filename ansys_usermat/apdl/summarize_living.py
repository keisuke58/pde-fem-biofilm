#!/usr/bin/env python3
"""summarize_living.py -- tables of the fixed coupling runs (A-22/A-23): for every fields_<job>.json of the given
prefixes in results/2026-10-living/, the seed and corner-column means at each stored time (phi, sum phibar = L,
share of A, c, a, alpha - 1, von Mises), the replay check of the point model (check_living_trace.py), and for the
sep_ab runs the seed's c and L over all increments (from the pm trace). Writes results/2026-10-living/A23_summary.md.

    python ansys_usermat/apdl/summarize_living.py [PREFIX ...]      (default: v2f_ fsens_)
The point-model file of each run: abaqus_composition/eco_fritsch_2sp_a50.txt, with alpha*_0 replaced for the
fsens_a<N> runs (alpha*_0 = N).
"""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = HERE / "results" / "2026-10-living" / "archive"   # the v2f_/fsens_ generations it summarises were archived (22 Oct)
ECO = ROOT / "abaqus_composition" / "eco_fritsch_2sp_a50.txt"
COLS = ("phi", "L", "shareA", "nut", "ab", "alpha", "seqv")


def eco_for(job, tmp):
    m = re.search(r"_a(\d+)$", job)
    if not m:
        return ECO
    lines = ECO.read_text().splitlines()
    lines[1] = f"100.0 {float(m.group(1)):.1f}"
    p = Path(tmp) / f"eco_{job}.txt"
    p.write_text("\n".join(lines) + "\n")
    return p


def check(job, tmp):
    if not (RES / f"{job}_pm.json").exists():
        return "no trace"
    r = subprocess.run([sys.executable, str(HERE / "check_living_trace.py"), job, "--eco", str(eco_for(job, tmp))],
                       capture_output=True, text=True)
    out = r.stdout.strip().splitlines()
    d = [l for l in out if "max |gnew" in l]
    return (out[-1] if out else r.stderr[-200:]) + (f" ({d[0].split('=')[-1].strip()})" if d else "")


def series(job):
    p = RES / f"{job}_pm.json"
    if not p.exists():
        return []
    d = json.load(open(p))
    C = {k: i for i, k in enumerate(d["columns"])}
    seed = str(d["picked"]["seed"])
    last = {}
    for r in d["seed_trace"]:
        if r[0] == seed and r[C["hit"]].strip() == "0":
            last[(int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
    return [(k, float(r[C["c1"]]), float(r[C["L"]]), float(r[C["phi_used"]])) for k, r in sorted(last.items())]


def main():
    prefixes = sys.argv[1:] or ["v2f_", "fsens_"]
    jobs = sorted(p.stem[len("fields_"):] for p in RES.glob("fields_*.json")
                  if not p.stem.endswith("_dense") and any(p.stem[len("fields_"):].startswith(x) for x in prefixes))
    out = ["# A-23 summary (summarize_living.py)", "",
           "Seed = the Chapter 4 seed elements, corners = elements with |x|, |y| > 0.75 mm; element means. "
           "L = sum phibar_i = sum phi_i psi_i, share A = phi_1/(phi_1 + phi_2), alpha - 1, vM in Pa.", ""]
    with tempfile.TemporaryDirectory() as tmp:
        for job in jobs:
            d = json.load(open(RES / f"fields_{job}.json"))
            el = np.array(d["elem"])
            seed = np.isin(el, d["seed"])
            corners = np.isin(el, d["edges"])
            out += [f"## {job}", "", f"trace check: {check(job, tmp)}", "",
                    "| T* | region | phi | L | share A | c | a | alpha - 1 | vM [Pa] |", "|---|---|---|---|---|---|---|---|---|"]
            for k, t in enumerate(d["t"]):
                for name, m in (("seed", seed), ("corners", corners)):
                    if not m.any():
                        continue
                    v = {c: float(np.array(d[c][k])[m].mean()) for c in COLS}
                    out.append(f"| {t:g} | {name} | {v['phi']:.3f} | {v['L']:.3f} | {v['shareA']:.3f} | {v['nut']:.3f} | "
                               f"{v['ab']:.3f} | {v['alpha']:.3e} | {v['seqv'] * 1e6:.3e} |")
            phi = np.array(d["phi"])
            out += ["", f"phi over all elements and stored times: [{phi.min():.3f}, {phi.max():.3f}]", ""]
            if "sep_ab" in job:
                s = series(job)
                if s:
                    out += [f"seed element {json.load(open(RES / f'{job}_pm.json'))['picked']['seed']}, "
                            "integration point 1, per increment (ldstep, isubst): c seen, L after the step, phi", "",
                            "| step | c | L | phi |", "|---|---|---|---|"]
                    out += [f"| {k} | {c:.4f} | {L:.4f} | {p:.4f} |" for k, c, L, p in s]
                    out.append("")
    p = RES / "A23_summary.md"
    p.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"wrote {p}: {len(jobs)} runs")


if __name__ == "__main__":
    main()
