#!/usr/bin/env python3
"""Fine scan of the isotropic share w of the front term (KLEMPT2024_REPRODUCTION.md
sec. 15): growth = r c/(k+c) |grad phi| [w + (1 - w) |n_phi . n_c|], DIAGNOSTIC.
Table 2 otherwise, consumption g phi c, 4.1 with the 2 um corner source; a
time scale per run. Scored against the digitised curves (RMS, phi and c averaged).

    python JAXFEM/klempt2024_blend_scan.py  -> JAXFEM/klempt2024_results/blend_scan.json
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_case2 as C  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402
import klempt2024_variant_search as VS  # noqa: E402

OUT = Path(__file__).resolve().parent / "klempt2024_results" / "blend_scan.json"
WS = [0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]
SCALES = {"4.1": [0.75, 1.0, 1.25], "4.2 high": [5.0, 6.0, 7.0], "4.2 low": [2.5, 3.0, 3.5]}


def corner(a=2.0):
    return (K.X >= K.L - a - 1e-9) & (K.Y >= K.L - a - 1e-9) & (K.Z >= K.L - a - 1e-9)


def one(job):
    w, case, s = job
    K.R, K.BETA, K.K_A = 100 * s, 2 * s, 1e-3 * s
    tg = VS.digitized()
    g = f"blend{w}"
    if case == "4.1":
        p, _, _ = K.setup("fig4_edge")
        rec, _, _ = B.run_setup(p, corner(), 1e8, "ic", "first_order", g)
        sc, rms = VS.score(rec, VS.T41, tg["4.1"])
    else:
        p, m = C.setup(1)
        rec, _, _ = B.run_setup(p, m, C.G[case.split()[1]], "ic", "first_order", g)
        sc, rms = VS.score(rec, C.T_CMP, tg[case])
    return w, case, s, sc, rms


def main():
    jobs = [(w, c, s) for w in WS for c, ss in SCALES.items() for s in ss]
    res = {}
    with Pool(4) as pool:
        for w, c, s, sc, rms in pool.imap_unordered(one, jobs):
            res.setdefault(f"{w:g}", {}).setdefault(c, {})[f"{s:g}"] = {"score": sc, **rms}
            print(f"w={w:g} {c:9s} s={s:<5g} score {sc:.3f}", flush=True)
            OUT.write_text(json.dumps(res, indent=1))
    print("\nbest per w (time scale per run):")
    for w in sorted(res, key=float):
        best = {c: min(v.items(), key=lambda kv: kv[1]["score"]) for c, v in res[w].items()}
        worst = max(b[1]["score"] for b in best.values())
        print(f"  w={w:5s} worst {worst:.3f}  " + "  ".join(
            f"{c} {b[1]['score']:.3f} (s={b[0]})" for c, b in sorted(best.items())))


if __name__ == "__main__":
    main()
