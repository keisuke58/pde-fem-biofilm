"""deck_components.py -- element centroids of a wired deck and of its components (NUTRIENT1, BIOFILM1, ...), and
helpers to write new element components for the coupling runs of 20 Oct 2026 (A-19).

    python ansys_usermat/apdl/deck_components.py DECK            # where each component lies
"""
import re
import sys
from pathlib import Path

import numpy as np


def read_mesh(path):
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    nodes, elems = {}, {}
    i = 0
    while i < len(lines):
        s = lines[i].strip().lower()
        if s.startswith("nblock"):
            i += 2
            while not lines[i].strip().startswith("-1") and lines[i].strip() and not lines[i].strip().lower().startswith("n,"):
                t = lines[i].split()
                if len(t) < 4:
                    break
                nodes[int(t[0])] = [float(x) for x in t[1:4]]
                i += 1
        elif s.startswith("eblock"):
            i += 2
            while i < len(lines) and not lines[i].strip().startswith("-1"):
                t = [int(x) for x in lines[i].split()]
                if len(t) >= 19:
                    elems[t[10]] = t[11:19]
                i += 1
        i += 1
    return nodes, elems


def centroids(nodes, elems):
    return {e: np.mean([nodes[n] for n in ns], axis=0) for e, ns in elems.items()}


def read_components(path):
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    comps = {}
    for i, l in enumerate(lines):
        m = re.match(r"\s*CMBLOCK,\s*(\w+)\s*,\s*ELEM\s*,\s*(\d+)", l, re.I)
        if m:
            name, n = m.group(1).upper(), int(m.group(2))
            vals, j = [], i + 2
            while len(vals) < n:
                vals += [int(x) for x in lines[j].split()]
                j += 1
            comps[name] = vals[:n]
    return comps


def cmblock(name, ids):
    ids = sorted(ids)
    out = [f"CMBLOCK,{name},ELEM,{len(ids):9d}", "(8i10)"]
    for k in range(0, len(ids), 8):
        out.append("".join(f"{x:10d}" for x in ids[k:k + 8]))
    return "\n".join(out)


def main(deck):
    nodes, elems = read_mesh(deck)
    c = centroids(nodes, elems)
    print(f"{len(nodes)} nodes, {len(elems)} elements")
    for name, ids in read_components(deck).items():
        x = np.array([c[e] for e in ids if e in c])
        if len(x):
            print(f"{name:12s} {len(ids):4d}  x [{x[:, 0].min():.3f}, {x[:, 0].max():.3f}]  "
                  f"y [{x[:, 1].min():.3f}, {x[:, 1].max():.3f}]  z [{x[:, 2].min():.3f}, {x[:, 2].max():.3f}]")


if __name__ == "__main__":
    main(sys.argv[1])
