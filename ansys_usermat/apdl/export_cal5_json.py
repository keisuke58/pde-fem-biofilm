#!/usr/bin/env python3
"""export_cal5_json.py -- the point-model trace of a calibrated five-species run
(research idea 7) as JSON, without rounding, so the cloud can compare it digit
by digit with the Python point model (ecology_jax).

    python ansys_usermat/apdl/export_cal5_json.py OUTDIR RUN --theta-json T.json
        [--source "Tmcmc202601@<sha>:docs/revision/generated/final_theta_MAP/CH.json"]

Reads from BIOFILM_WORKDIR (default F:\\biofilm_upf_ch5):
  comp_trace_<RUN>.csv  one row per point-model call: elem, integration point,
                        ldstep, isubst, n_sub, hit, dTime, coupling dt, phi_3D used,
                        amount the point model saw, alpha_n, alpha_new, g_old(12),
                        g_new(12) (g = phi_1..5, phi_0, psi_1..5, gamma)
  <RUN>.dat             first line (how the deck was made) and TBDATA (theta)
Writes OUTDIR/<RUN>_pm.json:
  seed_trace   every row of the element of the deck's --post-elem (integration
               point 1), the numbers as printed (strings, all digits)
  last_rows    the last row of every traced element and integration point 1
  theta        the 20 values of the MAP file and the deck's prop(8:27), as text
  constants    c*, alpha*, eta, n of the run (ecology_constants; the Fortran
               defaults the native build uses without BIOFILM_ECO_CASE)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "abaqus_composition")]

COLS = (["elem", "ip", "ldstep", "isubst", "n_sub", "hit", "dtime", "dt_couple",
         "phi_used", "phi_seen", "alpha_n", "alpha_new"]
        + [f"gold_{k}" for k in range(1, 13)] + [f"gnew_{k}" for k in range(1, 13)])


def deck_props(dat: Path):
    """prop(1..n) of TB,USER material 1 as the deck's text."""
    vals = {}
    for line in open(dat, errors="replace"):
        m = re.match(r"TBDATA,(\d+),(.*)", line.strip())
        if m:
            start = int(m.group(1))
            for k, v in enumerate(m.group(2).split(",")):
                if v.strip():
                    vals[start + k] = v.strip()
        elif vals and not line.startswith("TBDATA"):
            break
    return vals


def post_elem(dat: Path) -> str:
    for line in open(dat, errors="replace"):
        m = re.match(r"\*USE,post_elem_stress\.mac,(\d+)", line.strip())
        if m:
            return m.group(1)
    raise SystemExit(f"{dat}: no *USE,post_elem_stress.mac,<elem>")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("run")
    ap.add_argument("--theta-json", required=True)
    ap.add_argument("--source", default="")
    a = ap.parse_args()
    w = Path(os.environ.get("BIOFILM_WORKDIR") or r"F:\biofilm_upf_ch5")
    dat = w / f"{a.run}.dat"
    import ecology_constants as ec
    from write_eco_cfg import theta_json_cfg
    n, cs, als, eta, theta = theta_json_cfg(a.theta_json)
    props = deck_props(dat)
    deck_theta = [props.get(i, "0.0") for i in range(8, 28)]
    if [float(x) for x in deck_theta] != theta:
        raise SystemExit("deck prop(8:27) != MAP theta")
    seed = post_elem(dat)
    rows, last = [], {}
    with open(w / f"comp_trace_{a.run}.csv") as f:
        for r in csv.reader(f):
            r = [x.strip() for x in r]
            if len(r) != len(COLS) or r[1] != "1":
                continue
            if r[0] == seed:
                rows.append(r)
            last[r[0]] = r
    out = {
        "run": a.run,
        "deck_first_line": open(dat, errors="replace").readline().strip(),
        "theta_source": a.source,
        "theta": [repr(x) for x in theta],
        "deck_theta": deck_theta,
        "constants": {"n": n, "c_star": cs, "alpha_star": als, "eta": eta,
                      "K_hill": ec.K_HILL if hasattr(ec, "K_HILL") else 0.0},
        "columns": COLS,
        "seed_elem": int(seed),
        "seed_trace": rows,
        "last_rows": [last[k] for k in sorted(last, key=int)],
    }
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    p = Path(a.outdir) / f"{a.run}_pm.json"
    p.write_text(json.dumps(out, indent=0))
    print(f"wrote {p}: {len(rows)} seed rows, {len(last)} traced elements")


if __name__ == "__main__":
    main()
