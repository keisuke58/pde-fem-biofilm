"""compare_surface_fix.py -- what the surface-point fix (FRONT_TERM_FIX.md,
"phi = 0 at the surface points") changes in the thesis runs.

    python ansys_usermat/apdl/compare_surface_fix.py FIXDIR REF1DIR [REF2DIR ...] [--out file.md]

FIXDIR holds the JSON of the reruns with the fixed executable
(F:\\biofilm_upf_nativefix, export_runs_json.py); each run's reference JSON is
looked up in the REF dirs (the week chain / 5 Oct runs, executable without the
fix). Per run, reference / fixed / relative change of the measures of
summarize_runs_json.py (Pa): seed mean von Mises and p, layer 1 and 2 mean p,
alpha - 1 (seed interior, seed surface), and where the point model ran the seed
and layer-1 means of phi_1/(phi_1+phi_2). "max |d alpha|" is the largest change
of alpha - 1 of any element relative to the largest alpha - 1: the faces.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_runs_json import layers, surface  # noqa: E402


def measures(rec):
    a = rec["all_stress"]
    elem = np.array(a["elem"], int)
    cen = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    seed = np.isin(elem, rec.get("seed_BIOFILM1", []))
    lay, h = layers(cen, seed)
    vm = np.array(a["seqv"]) * 1e6
    p = (np.array(a["sx"]) + np.array(a["sy"]) + np.array(a["sz"])) / 3 * 1e6
    al = np.array(a["alpha"])
    surf = surface(cen, seed, h)
    inner = seed & ~surf
    m = {"seed vM": vm[seed].mean(), "seed p": p[seed].mean(),
         "layer 1 p": p[lay == 1].mean(), "layer 2 p": p[lay == 2].mean(),
         "alpha seed interior": al[inner].mean() if inner.any() else np.nan,
         "alpha seed surface": al[surf].mean()}
    nf = rec.get("nut_field")
    if nf:
        idx = {e: i for i, e in enumerate(elem)}
        order = np.array([idx[e] for e in nf["elem"]])
        p1, p2 = np.array(nf["phi1"]), np.array(nf["phi2"])
        s = p1 + p2
        ran = (s > 1e-6) & (s < 0.95)
        share = np.where(ran, p1 / np.where(ran, s, 1), np.nan)
        S, L = seed[order], lay[order]
        m["share seed"] = np.nanmean(share[S])
        if ((L == 1) & ran).any():
            m["share layer 1"] = np.nanmean(share[(L == 1) & ran])
    return m, dict(zip(elem.tolist(), al.tolist()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixdir")
    ap.add_argument("refdirs", nargs="+")
    ap.add_argument("--out")
    a = ap.parse_args()
    lines = ["# Surface-point fix: thesis runs before and after", "",
             "Reference: the runs as in the thesis (phi = 0 at the surface points). Fixed: the same decks with",
             "phi at the surface points from its neighbours (F:\\biofilm_upf_nativefix). Pa; change = (fixed - ref) / |ref|.", ""]
    worst = 0.0
    for fj in sorted(Path(a.fixdir).glob("*.json")):
        ref = next((Path(d) / fj.name for d in a.refdirs if (Path(d) / fj.name).exists()), None)
        if ref is None:
            continue
        rf, fx = json.loads(ref.read_text()), json.loads(fj.read_text())
        if "all_stress" not in rf or "all_stress" not in fx:
            lines.append(f"## {fj.stem}: no all_stress"); continue
        (mr, ar), (mf, af) = measures(rf), measures(fx)
        amax = max(abs(v) for v in ar.values())
        dmax = max(abs(af[e] - ar[e]) for e in ar if e in af) / amax
        lines += [f"## {fj.stem}", "", "| measure | reference | fixed | change |", "|---|---|---|---|"]
        for k in mr:
            if k in mf:
                ch = (mf[k] - mr[k]) / abs(mr[k]) if mr[k] else np.nan
                if k.startswith(("seed", "alpha", "share")):
                    worst = max(worst, abs(ch))
                lines.append(f"| {k} | {mr[k]:.4e} | {mf[k]:.4e} | {100 * ch:+.3f} % |")
        lines += [f"| max \\|d alpha\\| / max alpha (any element) | | | {100 * dmax:.3f} % |", ""]
    lines += [f"Largest change of a seed / alpha / share measure: {100 * worst:.3f} %", ""]
    txt = "\n".join(lines)
    print(txt)
    if a.out:
        Path(a.out).write_text(txt)


if __name__ == "__main__":
    main()
