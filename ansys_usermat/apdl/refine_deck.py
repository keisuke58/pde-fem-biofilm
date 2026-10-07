#!/usr/bin/env python3
"""refine_deck.py -- the partner's Workbench deck on a finer hexahedral mesh.

    python ansys_usermat/apdl/refine_deck.py BASE.dat OUT.dat --n 16

The partner's model is a 2 mm cube of 8 x 8 x 8 SOLID185 elements. For a mesh
study this writes the same deck with n x n x n elements: the nblock and the
eblock are replaced, and every component (CMBLOCK) is mapped onto the new mesh
by its place in space, so the boundary conditions stay where they were:
  - an element component (seed BIOFILM1/2, nutrient layers, ...): a new element
    belongs to it when its centroid lies inside one of the component's old
    elements;
  - a node component is classified first and then rebuilt: all nodes, all
    surface nodes, the nodes of one face, a single node (the three corner
    nodes with displacement constraints: same coordinates), or, failing those,
    the new nodes inside the old nodes' bounding box (reported).
Element attributes (mat, type, real, secnum, esys) and the node order inside an
element are copied from the old element that contains the new one. Nothing
else in the deck changes; the partner's routines count nodes and elements with
*GET at run time. Every component's old and new size is printed, and old
element numbers named with --track are mapped to the new elements that cover
them (for --elem / --post-elem of make_wired_deck.py).

Check before trusting it: `--n 8` must reproduce the base deck's components
(same sizes, same places); this is what `--selftest` does.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

CM = re.compile(r"^CMBLOCK,([^,]+),(NODE|ELEM),\s*(\d+)", re.I)


def read(path):
    return Path(path).read_text(encoding="latin-1").splitlines()


def parse(lines):
    i_nb = next(k for k, l in enumerate(lines) if l.lower().startswith("nblock"))
    xyz, k = {}, i_nb + 2
    while lines[k].strip() != "-1" and len(lines[k].split()) >= 4:
        p = lines[k].split()
        xyz[int(p[0])] = np.array([float(v) for v in p[1:4]])
        k += 1
    nb_end = k                                            # the "-1" line
    i_eb = next(k for k, l in enumerate(lines) if l.lower().startswith("eblock"))
    elems, k = {}, i_eb + 2
    while lines[k].strip() != "-1":
        p = [int(v) for v in lines[k].split()]
        elems[p[10]] = (p[:10], p[11:19])
        k += 1
    eb_end = k
    comps = []
    for j, l in enumerate(lines):
        m = CM.match(l)
        if not m:
            continue
        n = int(m.group(3))
        raw, r = [], j + 2
        while len(raw) < n and r < len(lines) and re.match(r"^\s*-?\d", lines[r]):
            raw += [int(v) for v in lines[r].split()]
            r += 1
        ids, prev = [], None
        for v in raw:
            ids += list(range(prev + 1, -v + 1)) if v < 0 else [v]
            prev = abs(v)
        comps.append({"name": m.group(1).strip(), "kind": m.group(2).upper(), "line": j,
                      "end": r, "ids": ids})
    return i_nb, nb_end, i_eb, eb_end, xyz, elems, comps


def build(n, lo, hi):
    ax = [np.linspace(lo[d], hi[d], n + 1) for d in range(3)]
    idx = np.arange(1, (n + 1) ** 3 + 1).reshape(n + 1, n + 1, n + 1)
    nodes = {int(idx[i, j, k]): np.array([ax[0][i], ax[1][j], ax[2][k]])
             for i in range(n + 1) for j in range(n + 1) for k in range(n + 1)}
    return ax, idx, nodes


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--track", default="220", help="old element numbers to map, comma separated")
    ap.add_argument("--selftest", action="store_true", help="--n 8 against the base, nothing written")
    a = ap.parse_args(argv)
    lines = read(a.base)
    i_nb, nb_end, i_eb, eb_end, xyz, elems, comps = parse(lines)
    P = np.array(list(xyz.values()))
    lo, hi = P.min(0), P.max(0)
    n = 8 if a.selftest else a.n
    ax, idx, nodes = build(n, lo, hi)
    tol = 1e-6 * float((hi - lo).max())

    # old elements: bounds, attributes, and the corner pattern of their node order
    old = {}
    for e, (attr, nn) in elems.items():
        c = np.array([xyz[v] for v in nn])
        old[e] = {"lo": c.min(0), "hi": c.max(0), "attr": attr,
                  "pat": [tuple((c[q] > c.mean(0)).astype(int)) for q in range(8)]}
    olo = np.array([o["lo"] for o in old.values()])
    ohi = np.array([o["hi"] for o in old.values()])
    okeys = list(old)

    def owner(p):
        hit = np.nonzero(np.all((olo - tol <= p) & (p <= ohi + tol), axis=1))[0]
        return okeys[hit[0]]

    new = {}
    e = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                e += 1
                cen = np.array([(ax[0][i] + ax[0][i + 1]) / 2, (ax[1][j] + ax[1][j + 1]) / 2,
                                (ax[2][k] + ax[2][k + 1]) / 2])
                o = owner(cen)
                nn = [int(idx[i + s[0], j + s[1], k + s[2]]) for s in old[o]["pat"]]
                new[e] = {"cen": cen, "attr": old[o]["attr"], "nodes": nn, "old": o}

    # components
    surf = lambda p: any(abs(p[d] - lo[d]) < tol or abs(p[d] - hi[d]) < tol for d in range(3))
    allsurf_old = {v for v, p in xyz.items() if surf(p)}
    out_comps = {}
    for c in comps:
        if c["kind"] == "ELEM":
            members = set(c["ids"])
            out_comps[c["name"]] = [e for e, r in new.items() if r["old"] in members]
            rule = "elements inside"
        else:
            s = set(c["ids"])
            pts = np.array([xyz[v] for v in c["ids"]])
            if s == set(xyz):
                rule, keep = "all nodes", lambda p: True
            elif s == allsurf_old:
                rule, keep = "surface nodes", surf
            elif len(s) == 1:
                q = pts[0]
                rule, keep = "single node", lambda p, q=q: np.allclose(p, q, atol=tol)
            else:
                face = [(d, b) for d in range(3) for b in (lo[d], hi[d]) if np.all(np.abs(pts[:, d] - b) < tol)]
                if face:
                    d, b = face[0]
                    rule, keep = f"face {'xyz'[d]} = {b:g}", lambda p, d=d, b=b: abs(p[d] - b) < tol
                else:
                    blo, bhi = pts.min(0), pts.max(0)
                    rule = "bounding box (check!)"
                    keep = lambda p, blo=blo, bhi=bhi: bool(np.all((blo - tol <= p) & (p <= bhi + tol)))
            out_comps[c["name"]] = [v for v, p in nodes.items() if keep(p)]
        c["rule"] = rule
        print(f"{c['name']:18s} {c['kind']:4s} {len(c['ids']):6d} -> {len(out_comps[c['name']]):6d}   {rule}")

    for t in [int(x) for x in a.track.split(",") if x]:
        cover = [e for e, r in new.items() if r["old"] == t]
        print(f"old element {t} -> new elements {cover}")

    if a.selftest:
        bad = [c["name"] for c in comps if len(out_comps[c["name"]]) != len(c["ids"])]
        print("SELFTEST", "FAIL: " + ", ".join(bad) if bad else "PASS (all component sizes kept at n = 8)")
        return 1 if bad else 0
    if not a.out:
        sys.exit("give OUT.dat (or --selftest)")

    out = lines[:i_nb]
    out.append(f"nblock,3,,{len(nodes)}")
    out.append(lines[i_nb + 1])
    out += [f"{v:9d}{p[0]:20.9E}{p[1]:20.9E}{p[2]:20.9E}" for v, p in nodes.items()]
    out += lines[nb_end:i_eb]
    out.append(re.sub(r",\d+\s*$", f",{len(new)}", lines[i_eb]))
    out.append(lines[i_eb + 1])
    for e, r in new.items():
        rec = r["attr"] + [e] + r["nodes"]
        out.append("".join(f"{v:9d}" for v in rec))
    pos = eb_end
    for c in comps:
        out += lines[pos:c["line"]]
        ids = out_comps[c["name"]]
        out.append(f"CMBLOCK,{c['name']:8s},{c['kind']},{len(ids):9d}")
        out.append("(8i10)")
        out += ["".join(f"{v:10d}" for v in ids[q:q + 8]) for q in range(0, len(ids), 8)]
        pos = c["end"]
    out += lines[pos:]
    Path(a.out).write_text("\n".join(out) + "\n", encoding="latin-1")
    print(f"wrote {a.out}: {len(nodes)} nodes, {len(new)} elements")
    # read back
    _, _, _, _, xyz2, el2, comps2 = parse(read(a.out))
    assert len(xyz2) == len(nodes) and len(el2) == len(new)
    for c in comps2:
        assert len(c["ids"]) == len(out_comps[c["name"]]), c["name"]
    print("read back: nodes, elements and every component as written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
