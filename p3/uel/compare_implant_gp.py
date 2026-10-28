"""compare_implant_gp.py -- compares an Abaqus run of the implant sector (the Gauss-point record uel_gp.txt
that UEXTERNALDB writes, last block = end of the analysis) with the Python emulation's results/implant_*.json
(make_implant_uel_inp.py --run), element by element.

    python compare_implant_gp.py <uel_gp.txt> <implant_<nr>x<nt>x<nz>[_b1].json> [--tol 1e-5]

Abaqus numbers the elements as the deck does (deck order = results order), so element e of the record is
element e - 1 of the .json. Prints the seed means of both, the largest vM/E of both and the largest
element-wise difference of alpha - 1 and vM/E; exit code 1 when a difference exceeds --tol (the .json keeps
7 significant digits, so differences below about 1e-7 are rounding).
"""
import argparse
import json
import sys

import numpy as np


def last_block(path):
    """element means (alpha - 1, von Mises) of the last block of uel_gp.txt, (E, 2), 0-based element order"""
    blocks, cur = [], None
    for line in open(path):
        if line.startswith("# step"):
            cur = []; blocks.append(cur)
        elif line.startswith("#") or not line.strip():
            continue
        elif cur is not None:
            cur.append(line.split())
    if not blocks:
        raise SystemExit(f"no record block in {path}")
    rec = np.array(blocks[-1], float)                                   # element gp x y z phi c alpha-1 living vM p
    ne = int(rec[:, 0].max())
    out = np.zeros((ne, 2)); cnt = np.zeros(ne)
    for row in rec:
        e = int(row[0]) - 1
        out[e] += row[[7, 9]]; cnt[e] += 1
    return out / cnt[:, None]


def compare(gp_path, json_path, tol=1e-5):
    r = json.load(open(json_path))
    el = r["elements"]
    ab = last_block(gp_path)
    py = np.c_[el["alpha_minus_1"], el["vM_over_E"]]
    if len(ab) != len(py):
        raise SystemExit(f"element counts differ: Abaqus {len(ab)}, Python {len(py)}")
    seed = np.array(el["layer"]) == 0
    d = np.abs(ab - py)
    lines = [f"{len(ab)} elements, seed {seed.sum()}",
             f"seed alpha - 1: Abaqus {ab[seed, 0].mean():.7g}, Python {py[seed, 0].mean():.7g}",
             f"seed vM/E:      Abaqus {ab[seed, 1].mean():.7g}, Python {py[seed, 1].mean():.7g}",
             f"max vM/E:       Abaqus {ab[:, 1].max():.7g} (element {ab[:, 1].argmax() + 1}), "
             f"Python {py[:, 1].max():.7g} (element {py[:, 1].argmax() + 1})",
             f"largest element-wise difference: alpha - 1 {d[:, 0].max():.2e}, vM/E {d[:, 1].max():.2e}"]
    ok = d.max() <= tol
    lines.append("PASS" if ok else f"FAIL (tol {tol:g})")
    return ok, "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("gp"); ap.add_argument("json"); ap.add_argument("--tol", type=float, default=1e-5)
    a = ap.parse_args()
    ok, text = compare(a.gp, a.json, a.tol)
    print(text)
    sys.exit(0 if ok else 1)
