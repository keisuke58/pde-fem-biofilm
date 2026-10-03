#!/usr/bin/env python3
"""One setting for all of Klempt 2024's curves: a joint search.

Fitting each figure on its own proves little: enough knobs fit anything. The
test here is stricter. A setting (k, r, g) must reproduce, at once, test case
4.1 (Fig. 4) and both sub-cases of 4.2 (Fig. 7, "high" and "low"), phi and c,
against the curves read off the PDF (klempt2024_digitize.py). Only the two
things the paper's figures leave open are chosen per case: the width of the
nutrient strip in 4.1 (Fig. 2 is a sketch) and the thickness of the initial
disk in 4.2 (the paper's mean phi at T* = 0.01 is ten times one node layer).

Model form, from the single-case studies (klempt2024_case1_bc.py,
klempt2024_case2.py): growth on both faces |grad phi . n_c| (the paper's text,
not Eq. 34's sign) and first-order consumption g phi c (Fig. 7's plateau),
seed as an initial value. Values searched (Table 2: k = 1, r = 100, g = 1e8 /
1e10 for 4.2 "high" / "low", 1e8 for 4.1):
    k in {1, 0.3, 0.1, 0.03},  r x {1, 2, 4},  g x {1, 2, 4}

Score per curve: root mean square difference over the compared times, phi and
c averaged. A setting's score is the worst of its three curves, so it has to
hold everywhere.

    python JAXFEM/klempt2024_variant_search.py      (about an hour)
        -> JAXFEM/klempt2024_results/variant_search.json
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

OUT = Path(__file__).resolve().parent / "klempt2024_results" / "variant_search.json"
KS = [1.0, 0.3, 0.1, 0.03]
RS = [1.0, 2.0, 4.0]
GS = [1.0, 2.0, 4.0]
STRIPS = ["band1", "band2"]
LAYERS = [1, 9]
T41 = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]


def digitized():
    d = json.loads((C.RES / "paper_curves_digitized.json").read_text())
    t = np.array(d["t"])

    def at(vals, ts):
        y = np.array([np.nan if v is None else v for v in vals])
        ok = ~np.isnan(y)
        return np.interp(ts, t[ok], y[ok])
    return {"4.1": {k: at(d["fig4"][k], T41) for k in ("phi", "c")},
            "4.2 high": {k: at(d["fig7"][f"{k}_high"], C.T_CMP) for k in ("phi", "c")},
            "4.2 low": {k: at(d["fig7"][f"{k}_low"], C.T_CMP) for k in ("phi", "c")}}


def score(rec, ts, target):
    tt = np.array(rec["t"])
    rms = {k: float(np.sqrt(np.mean((np.interp(ts, tt, rec[k]) - target[k]) ** 2)))
           for k in ("phi", "c")}
    return 0.5 * (rms["phi"] + rms["c"]), rms


def run41(strip, g):
    phi0, _, _ = K.setup("fig4_edge")
    rec, _, _ = B.run_setup(phi0, B.nutrient_mask(strip), g, "ic", "first_order", "abs")
    return rec


def run42(sub, layers, g):
    phi0, mask = C.setup(layers)
    rec, _, _ = B.run_setup(phi0, mask, g, "ic", "first_order", "abs")
    return rec


def main():
    tg = digitized()
    try:
        res = json.loads(OUT.read_text())
    except (OSError, ValueError):
        res = {}
    for k, rs, gs in itertools.product(KS, RS, GS):
        key = f"k={k:g}/r_x{rs:g}/g_x{gs:g}"
        if key in res:
            continue
        t0 = time.time()
        K.K_M, K.R = k, 100.0 * rs
        entry = {}
        best = None
        for strip in STRIPS:
            s, rms = score(run41(strip, 1e8 * gs), T41, tg["4.1"])
            if best is None or s < best[0]:
                best = (s, rms, strip)
        entry["4.1"] = {"score": best[0], "rms": best[1], "strip": best[2]}
        for sub, g in (("high", 1e8), ("low", 1e10)):
            best = None
            for layers in LAYERS:
                s, rms = score(run42(sub, layers, g * gs), C.T_CMP, tg[f"4.2 {sub}"])
                if best is None or s < best[0]:
                    best = (s, rms, layers)
            entry[f"4.2 {sub}"] = {"score": best[0], "rms": best[1], "layers": best[2]}
        entry["worst"] = max(entry[c]["score"] for c in ("4.1", "4.2 high", "4.2 low"))
        res[key] = entry
        OUT.write_text(json.dumps(res, indent=1))
        print(f"{key:24s} worst {entry['worst']:.3f}  "
              + "  ".join(f"{c} {entry[c]['score']:.3f}" for c in ("4.1", "4.2 high", "4.2 low"))
              + f"  ({time.time() - t0:.0f}s)", flush=True)
    K.K_M, K.R = 1.0, 100.0
    print("\nbest settings (worst-curve score):")
    for key, e in sorted(res.items(), key=lambda kv: kv[1]["worst"])[:10]:
        print(f"  {key:24s} {e['worst']:.3f}   4.1 {e['4.1']['score']:.3f} ({e['4.1']['strip']})"
              f"   high {e['4.2 high']['score']:.3f} (layers {e['4.2 high']['layers']})"
              f"   low {e['4.2 low']['score']:.3f} (layers {e['4.2 low']['layers']})")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
