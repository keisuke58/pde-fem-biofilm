#!/usr/bin/env python3
"""check_deck_symmetry.py -- is a run deck mirror-symmetric in x and z? (addendum 10, check S1)

The examples of the agreed coupling are symmetric under x -> -x and z -> -z, but the ANSYS results are not
(up to 10 % in phi on 8^3, ~1 % in the Chapter 4 runs on 16^3). This checks the deck itself: node coordinates,
every element component (seed, nutrient layer, corner blocks, antibiotic layer, ...), and the node order inside
each element (a mirrored element whose local node numbering is not mirrored is fine for standard elements, but an
element that builds its own neighbour stencil from the node order may then be direction-biased; reported only).

    python ansys_usermat/apdl/check_deck_symmetry.py DECK.dat
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deck_components import centroids, read_components, read_mesh  # noqa: E402


def main(deck):
    nodes, elems = read_mesh(deck)
    comps = read_components(deck)
    xyz = {n: np.round(np.array(c, float), 9) for n, c in nodes.items()}
    nkey = {tuple(c): n for n, c in xyz.items()}
    cen = {e: np.round(c, 9) for e, c in centroids(nodes, elems).items()}
    ekey = {tuple(c): e for e, c in cen.items()}
    ok = True
    print(f"{deck}: {len(nodes)} nodes, {len(elems)} elements, components {sorted(comps)}")
    for ax, m in (("x", np.array([-1, 1, 1])), ("z", np.array([1, 1, -1]))):
        miss_n = sum(tuple(np.round(m * c, 9)) not in nkey for c in xyz.values())
        emap = {e: ekey.get(tuple(np.round(m * c, 9))) for e, c in cen.items()}
        miss_e = sum(v is None for v in emap.values())
        ok &= miss_n == 0 and miss_e == 0
        print(f"  mirror {ax}: nodes without a mirror {miss_n}, elements without a mirror {miss_e}")
        for name, ids in sorted(comps.items()):
            s = set(ids)
            bad = [e for e in s if e in emap and emap[e] not in s]
            ok &= not bad
            print(f"    component {name:12s} {len(s):5d} elements, not mirrored: {len(bad)}"
                  + (f" (e.g. {bad[:5]})" if bad else ""))
        # local node order: mirror each node of e and compare with the mirror element's node list position by position
        same = diff = 0
        for e, ns in elems.items():
            f = emap.get(e)
            if f is None:
                continue
            mirrored = [nkey.get(tuple(np.round(m * xyz[n], 9))) for n in ns]
            if mirrored == list(elems[f]):
                same += 1
            else:
                diff += 1
        print(f"    local node order of mirrored element = mirror of the node order: {same}, different: {diff}"
              " (information only)")
    print("PASS" if ok else "FAIL: the deck is not mirror-symmetric")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
