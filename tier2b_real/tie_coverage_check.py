#!/usr/bin/env python3
"""How much of each *TIE interface is actually tied? Read it off the .inp.

The cylinder decks were bonded with NUMMRG, which merges only coincident
nodes; in the 4-condition deck that came to four corner nodes out of the whole
interface, and every stress those decks produced was wrong (see
ansys_usermat/apdl/README.md). *TIE is a different mechanism and does not have
that bug -- it constrains a slave surface to a master surface by equations, so
the two meshes need not match.

It has its own version of the same failure, though. With ADJUST=NO, a slave
node further from the master surface than POSITION TOLERANCE is **dropped from
the tie**, and Abaqus says so in a warning rather than stopping. A tolerance
set by hand (build_assembly.py uses 0.5, 0.6, 1.0 and 2.8 mm) can therefore
leave part of an interface unbonded while the job still runs and reports
stresses -- which is exactly how the cylinder bug stayed invisible.

The other direction matters too: a node well inside the tolerance but at a real
distance is tied rigidly *across that gap*, which stiffens the joint. So the
useful output is not a pass/fail but the gap distribution against the
tolerance, per tie.

This reads only the generated .inp -- nodes, elements, surfaces and the *TIE
cards with their tolerances -- so it needs no Abaqus, no solve and no .odb, and
it can be run before submitting a job.

    python tier2b_real/tie_coverage_check.py tier2b_real.inp

numpy only, deliberately: scipy is not installed everywhere this has to run.
"""
import argparse
import re
import sys
from collections import defaultdict

import numpy as np

# Abaqus face -> local node indices (1-based), corners first.
FACES = {
    "C3D4":  {1: (1, 2, 3), 2: (1, 4, 2), 3: (2, 4, 3), 4: (3, 4, 1)},
    "C3D10": {1: (1, 2, 3, 5, 6, 7), 2: (1, 4, 2, 8, 9, 5),
              3: (2, 4, 3, 9, 10, 6), 4: (3, 4, 1, 10, 8, 7)},
}


def parse_inp(path):
    nodes, elems, etype = {}, {}, {}
    surfaces = defaultdict(list)          # name -> [(eid, face_id), ...]
    ties = []                             # (name, slave, master, tol)
    mode, cur, cur_type = None, None, None
    with open(path, errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                kw = line.split(",")[0].strip().upper()
                up = line.upper()
                if kw == "*NODE":
                    mode = "node"
                elif kw == "*ELEMENT":
                    m = re.search(r"TYPE\s*=\s*([A-Z0-9]+)", up)
                    cur_type = m.group(1) if m else None
                    mode = "elem"
                elif kw == "*SURFACE":
                    m = re.search(r"NAME\s*=\s*([^,\s]+)", up)
                    cur = m.group(1) if m else None
                    mode = "surf" if "TYPE=ELEMENT" in up.replace(" ", "") else None
                elif kw == "*TIE":
                    m = re.search(r"NAME\s*=\s*([^,\s]+)", up)
                    t = re.search(r"POSITION\s+TOLERANCE\s*=\s*([0-9.eE+-]+)", up)
                    ties.append([m.group(1) if m else "?", None, None,
                                 float(t.group(1)) if t else None])
                    mode = "tie"
                else:
                    mode = None
                continue
            if mode == "node":
                p = [x for x in line.split(",") if x.strip()]
                if len(p) >= 4:
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif mode == "elem":
                p = [x for x in line.split(",") if x.strip()]
                if len(p) >= 2:
                    eid = int(p[0])
                    elems[eid] = [int(x) for x in p[1:]]
                    etype[eid] = cur_type
            elif mode == "surf" and cur:
                p = [x.strip() for x in line.split(",") if x.strip()]
                if len(p) >= 2 and p[1].upper().startswith("S"):
                    try:
                        surfaces[cur].append((int(p[0]), int(p[1][1:])))
                    except ValueError:
                        pass
            elif mode == "tie" and ties and ties[-1][1] is None:
                p = [x.strip().upper() for x in line.split(",") if x.strip()]
                if len(p) >= 2:
                    ties[-1][1], ties[-1][2] = p[0], p[1]
                    mode = None
    return nodes, elems, etype, surfaces, ties


def face_nodes(surf, elems, etype, corners_only):
    """Node ids on a surface; corners only when the geometry is what's wanted."""
    out = []
    for eid, fid in surf:
        conn, t = elems.get(eid), etype.get(eid)
        loc = FACES.get(t, {}).get(fid)
        if conn is None or loc is None:
            continue
        loc = loc[:3] if corners_only else loc
        out.extend(conn[i - 1] for i in loc if i - 1 < len(conn))
    return sorted(set(out))


def face_triangles(surf, elems, etype, nodes):
    tri = []
    for eid, fid in surf:
        conn, t = elems.get(eid), etype.get(eid)
        loc = FACES.get(t, {}).get(fid)
        if conn is None or loc is None:
            continue
        ids = [conn[i - 1] for i in loc[:3]]
        if all(i in nodes for i in ids):
            tri.append([nodes[i] for i in ids])
    return np.asarray(tri, dtype=float)


def _tri_dist(p, V):
    """Exact distance from one point to each triangle (Ericson, RTCD 5.1.5).

    Vertex distance is not good enough here. A node sitting 0.1 mm above the
    middle of a face is 0.1 mm from the surface but can be several mm from any
    of its vertices, and this routine's whole job is to say how far a node is
    from the surface it is meant to be tied to. Using the vertex distance would
    report interfaces as unbonded that are perfectly fine.
    """
    a, b, c = V[:, 0], V[:, 1], V[:, 2]
    ab, ac, ap_ = b - a, c - a, p - a
    d1 = (ab * ap_).sum(1); d2 = (ac * ap_).sum(1)
    bp = p - b
    d3 = (ab * bp).sum(1); d4 = (ac * bp).sum(1)
    cp = p - c
    d5 = (ab * cp).sum(1); d6 = (ac * cp).sum(1)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    den = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v_ab = np.clip(np.where(d1 - d3 != 0, d1 / (d1 - d3), 0.0), 0, 1)
        v_ac = np.clip(np.where(d2 - d6 != 0, d2 / (d2 - d6), 0.0), 0, 1)
        den_bc = (d4 - d3) + (d5 - d6)
        v_bc = np.clip(np.where(den_bc != 0, (d4 - d3) / den_bc, 0.0), 0, 1)
        inv = np.where(den != 0, 1.0 / np.where(den == 0, 1.0, den), 0.0)
        vv, ww = vb * inv, vc * inv

    Q = a + vv[:, None] * ab + ww[:, None] * ac           # interior case
    Q = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, Q)
    Q = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, Q)
    Q = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, Q)
    Q = np.where(((vc <= 0) & (d1 >= 0) & (d3 <= 0))[:, None],
                 a + v_ab[:, None] * ab, Q)
    Q = np.where(((vb <= 0) & (d2 >= 0) & (d6 <= 0))[:, None],
                 a + v_ac[:, None] * ac, Q)
    Q = np.where(((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0))[:, None],
                 b + v_bc[:, None] * (c - b), Q)
    return np.sqrt(((p - Q) ** 2).sum(1))


def point_tri_dist(P, T):
    """Distance from each point to the nearest triangle, bucketed by centroid.

    Checking every point against every triangle is O(len(P)*len(T)); the
    buckets cut it to the triangles that could plausibly be nearest, then
    _tri_dist gives the exact distance to each of those.
    """
    if len(T) == 0:
        return np.full(len(P), np.inf)
    cen = T.mean(axis=1)
    lo, hi = cen.min(axis=0), cen.max(axis=0)
    span = np.maximum(hi - lo, 1e-9)
    n = max(1, int(round(len(T) ** (1 / 3))))
    key = lambda X: np.clip(((X - lo) / span * n).astype(int), 0, n - 1)
    grid = defaultdict(list)
    for i, k in enumerate(key(cen)):
        grid[tuple(k)].append(i)
    out = np.empty(len(P))
    for j, (p, k) in enumerate(zip(P, key(P))):
        cand, r = [], 0
        while not cand and r <= n:
            for dx in range(-r, r + 1):
                for dy in range(-r, r + 1):
                    for dz in range(-r, r + 1):
                        if max(abs(dx), abs(dy), abs(dz)) != r:
                            continue
                        cand.extend(grid.get((k[0] + dx, k[1] + dy, k[2] + dz), []))
            r += 1
        V = T[cand] if cand else T
        out[j] = _tri_dist(p, V).min()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("--tol", type=float, default=None,
                    help="override every tie's POSITION TOLERANCE")
    a = ap.parse_args(argv)

    nodes, elems, etype, surfaces, ties = parse_inp(a.inp)
    print(f"{a.inp}: {len(nodes)} nodes, {len(elems)} elements, "
          f"{len(surfaces)} surfaces, {len(ties)} ties\n")
    if not ties:
        print("no *TIE cards -- nothing to check")
        return 0

    worst = 0
    for name, slave, master, tol in ties:
        tol = a.tol if a.tol is not None else tol
        if slave not in surfaces or master not in surfaces:
            print(f"{name}: surfaces {slave}/{master} not found as element "
                  f"surfaces -- skipped (node-based surface?)")
            continue
        sn = face_nodes(surfaces[slave], elems, etype, corners_only=False)
        P = np.array([nodes[i] for i in sn if i in nodes], dtype=float)
        T = face_triangles(surfaces[master], elems, etype, nodes)
        d = point_tri_dist(P, T)

        print(f"{name}:  {slave} -> {master}   tolerance {tol}")
        print(f"  slave nodes {len(P)},  master faces {len(T)}")
        print(f"  gap  min {d.min():.4f}  median {np.median(d):.4f}  "
              f"p95 {np.percentile(d, 95):.4f}  max {d.max():.4f}")
        if tol:
            out = int((d > tol).sum())
            pct = 100.0 * out / max(len(d), 1)
            print(f"  beyond tolerance: {out} of {len(d)} ({pct:.2f}%) "
                  f"-> DROPPED from the tie")
            far = int((d > 0.5 * tol).sum())
            print(f"  tied across more than half the tolerance: {far} "
                  f"({100.0*far/max(len(d),1):.2f}%) -> stiffened, not dropped")
            worst = max(worst, pct)
            if pct > 1.0:
                print("  ** more than 1% of this interface is not tied **")
        print()

    print("Distances are exact point-to-triangle against the master faces --")
    print("the same quantity Abaqus compares POSITION TOLERANCE against.")
    return 1 if worst > 1.0 else 0


if __name__ == "__main__":
    sys.exit(main())
