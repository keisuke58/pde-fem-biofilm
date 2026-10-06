"""partner_elem_sets.py -- the nutrient and seed elements of a partner-style deck
as (i, j, k) indices on its n^3 grid, for make_klempt_inp.py --elem-sets, so an
Abaqus run uses the same regions as a run of the partner-style deck in ANSYS.

    python abaqus_composition/partner_elem_sets.py DECK.dat OUT.json [--case 41|42]

--case 42: the deck's NUTRIENT1 and BIOFILM1 components (the Fig. 7 set-up);
--case 41: the nutrient on the corner element at (max x, max y, max z) and a
seed of the elements whose centroid lies within L/4 of the centre (Klempt 2024
test case 4.1: a 5 um sphere in the 20 um cube), as built for the ANSYS runs.
Only the mesh and the component lists are read.
"""
import argparse
import json
import re
from pathlib import Path


def read(deck):
    lines = Path(deck).read_text(errors="replace").splitlines()
    nodes, elems, comps, i = {}, {}, {}, 0
    while i < len(lines):
        l = lines[i]
        if l.lower().startswith("nblock"):
            i += 2
            while not lines[i].startswith("-1"):
                p = lines[i].split()
                nodes[int(p[0])] = [float(x) for x in p[1:4]]
                i += 1
        elif l.lower().startswith("eblock"):
            i += 2
            while not lines[i].strip().startswith("-1"):
                p = [int(x) for x in lines[i].split()]
                if len(p) >= 19:
                    elems[p[10]] = p[11:19]
                i += 1
        elif l.upper().startswith("CMBLOCK") and ",ELEM" in l.upper():
            name, cnt = l.split(",")[1].strip().upper(), int(l.split(",")[3])
            ids, i = [], i + 2
            while len(ids) < cnt:
                ids += [int(x) for x in re.findall(r"-?\d+", lines[i])]
                i += 1
            full = []
            for v in ids:
                full += list(range(full[-1] + 1, -v + 1)) if v < 0 else [v]
            comps[name] = full
            continue
        i += 1
    return nodes, elems, comps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("deck")
    ap.add_argument("out")
    ap.add_argument("--case", choices=("41", "42"), default="42")
    a = ap.parse_args()
    nodes, elems, comps = read(a.deck)
    cen = {e: [sum(nodes[n][k] for n in c) / 8 for k in range(3)] for e, c in elems.items()}
    lo = [min(v[k] for v in nodes.values()) for k in range(3)]
    hi = [max(v[k] for v in nodes.values()) for k in range(3)]
    n = round(len(elems) ** (1 / 3))
    h = [(hi[k] - lo[k]) / n for k in range(3)]
    ijk = {e: [int((c[k] - lo[k]) // h[k]) for k in range(3)] for e, c in cen.items()}
    if a.case == "42":
        nut, seed = comps["NUTRIENT1"], comps["BIOFILM1"]
    else:
        L = hi[0] - lo[0]
        mid = [(p + q) / 2 for p, q in zip(lo, hi)]
        nut = [max(cen, key=lambda e: sum(cen[e]))]
        seed = [e for e, c in cen.items() if sum((c[k] - mid[k]) ** 2 for k in range(3)) <= (0.25 * L) ** 2 + 1e-12]
    out = {"n": n, "deck": str(a.deck), "case": a.case,
           "nut": sorted(ijk[e] for e in nut), "seed": sorted(ijk[e] for e in seed)}
    Path(a.out).write_text(json.dumps(out))
    print(f"{a.out}: n = {n}, nutrient {len(nut)} elements, seed {len(seed)} elements")


if __name__ == "__main__":
    main()
