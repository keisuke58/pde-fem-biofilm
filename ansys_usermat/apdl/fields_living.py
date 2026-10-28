"""fields_living.py -- collect fields_1..4.csv of post_fields_living.mac (A-19, the coupling with the living share)
into results/2026-10-living/fields_<run>.json and print the seed / edge-column summary.

    python ansys_usermat/apdl/fields_living.py RUN DIR_WITH_CSV DECK

JSON: {"run", "t": [4 set times], "elem", "cx", "cy", "cz", "seed", "edges", and per time [time][element]:
"phi", "nut", "ab", "L", "shareA", "alpha", "seqv"}, 7 significant digits. phi = SVAR 41, nut = SVAR 51 (c), ab =
SVAR 52 (a, the antibiotic), L = SVAR 87 (living share of the point model, 0 before it starts), shareA = SVAR 89
(phi_1/(phi_1 + phi_2), 0 before the start), alpha = SVAR 84 (alpha - 1), seqv in MPa (deck /units,MPA).
seed = the Chapter 4 seed elements (base component BIOFILM2 of the original deck is not used; the seed is the
deck's BIOFILM1 minus the edge columns), edges = elements with |x|, |y| > 0.75 mm.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from deck_components import read_components  # noqa: E402

COLS = ("phi", "nut", "ab", "L", "shareA", "alpha", "seqv")


def g7(x):
    return float(f"{x:.7g}")


def main(run, d, deck, n="4"):
    """n = 4: post_fields_living.mac; n = 40: post_fields_living_dense.mac (fields_<run>_dense.json, no trace);
    n = 11: post_fields_living_early.mac (fields_<run>_early.json, no trace)"""
    n = int(n)
    out = {"run": run, "t": []}
    for k in range(1, n + 1):
        rows = [{c.strip(): v for c, v in r.items()} for r in csv.DictReader(open(Path(d) / f"fields_{k}.csv"))]
        rows = [r for r in rows if any(float(r[c]) != 0 for c in ("cx", "cy", "cz"))]
        out["t"].append(g7(float(rows[0]["t"])))
        if k == 1:
            out["elem"] = [int(float(r["elem"])) for r in rows]
            for c in ("cx", "cy", "cz"):
                out[c] = [g7(float(r[c])) for r in rows]
        else:
            assert [int(float(r["elem"])) for r in rows] == out["elem"]
        for c in COLS:
            out.setdefault(c, []).append([g7(float(r[c])) for r in rows])
    cx, cy = np.array(out["cx"]), np.array(out["cy"])
    edges = (np.abs(cx) > 0.75) & (np.abs(cy) > 0.75)
    bio1 = set(read_components(deck)["BIOFILM1"])
    seed = np.array([e in bio1 for e in out["elem"]]) & ~edges
    out["seed"] = [int(e) for e, s in zip(out["elem"], seed) if s]
    out["edges"] = [int(e) for e, s in zip(out["elem"], edges) if s]
    res = HERE / "results" / "2026-10-living"
    res.mkdir(parents=True, exist_ok=True)
    p = res / {4: f"fields_{run}.json", 11: f"fields_{run}_early.json"}.get(n, f"fields_{run}_dense.json")
    p.write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote {p}: t = {out['t']}, {len(out['elem'])} elements, seed {seed.sum()}, edges {edges.sum()}")
    print("| T* | region | phi | L | share A | c | a | alpha-1 | vM [Pa] |")
    for k, t in enumerate(out["t"]):
        if n == 40 and (k + 1) % 10:
            continue
        for name, m in (("seed", seed), ("edges", edges)):
            v = {c: np.array(out[c][k])[m].mean() for c in COLS}
            print(f"| {t:g} | {name} | {v['phi']:.4f} | {v['L']:.4f} | {v['shareA']:.4f} | {v['nut']:.4f} | "
                  f"{v['ab']:.4f} | {v['alpha']:.4e} | {v['seqv'] * 1e6:.4e} |")
    phi = np.array(out["phi"])
    print(f"phi range over all elements and times: [{phi.min():.4f}, {phi.max():.4f}]")
    ab = np.array(out["ab"])
    print(f"a range over all elements and times: [{ab.min():.4g}, {ab.max():.4g}]")
    if n == 4:
        trace(run, d, out, seed, edges)


TRACE_COLS = (["elem", "ip", "ldstep", "isubst", "n_sub", "hit", "dtime", "dt_couple", "phi_used", "phi_seen",
               "alpha_n", "alpha_new"] + [f"gold_{i}" for i in range(1, 13)] + [f"gnew_{i}" for i in range(1, 13)]
              + ["c0", "c1", "a0", "a1", "L"])


def trace(run, d, out, seed, edges):
    """point-model trace (comp_trace.csv, integration point 1) of the seed element nearest the centre and of the
    edge element nearest the column's centre: every call, all digits, for the cloud's reproduction (check 1)"""
    p = Path(d) / "comp_trace.csv"
    if not p.exists():
        return
    cx, cy, cz = (np.array(out[c]) for c in ("cx", "cy", "cz"))
    el = np.array(out["elem"])
    pick = {}
    r2 = cx ** 2 + cy ** 2 + cz ** 2
    pick["seed"] = int(el[seed][np.argmin(r2[seed])])
    if edges.any():
        q = (np.abs(cx) - 0.875) ** 2 + (np.abs(cy) - 0.875) ** 2 + cz ** 2
        pick["edge"] = int(el[edges][np.argmin(q[edges])])
    want = {str(v) for v in pick.values()}
    rows = [r for r in csv.reader(open(p)) if r and r[0] in want and r[1] == "1"]
    # the layout check_living_trace.py reads (<job>_pm.json: columns, seed_trace)
    res = HERE / "results" / "2026-10-living" / f"{run}_pm.json"
    res.write_text(json.dumps({"run": run, "columns": TRACE_COLS, "picked": pick,
                               "note": "every call at integration point 1 (equilibrium iterations included; hit = 1 "
                                       "is a cached repeat, 2 = point model not started); c0, c1, a0, a1 the "
                                       "nutrient and antibiotic pairs interpolated over the inner steps, L = sum "
                                       "phi_i psi_i after the step",
                               "seed_trace": [[x.strip() for x in r] for r in rows]}, separators=(",", ":")))
    print(f"wrote {res}: {pick}, {len(rows)} rows")


if __name__ == "__main__":
    main(*sys.argv[1:5])
