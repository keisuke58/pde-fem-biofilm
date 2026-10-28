#!/usr/bin/env python3
"""nem_getnhgb_counts.py -- what the element's neighbour search actually did, read from an ANSYS output written with
DEBUG_KEY_NEM = 100 (addendum 15 (1) and 16).

With that deck parameter the element prints, for every point (Gauss points, then surface nodes), one record before
and one after the enlargement of the search radius: the number of neighbours found inside the initial radius, the
length of the extended list, the neighbour offsets (dx, dy, dz) and IDs, the number of enlargements and the final
radius. Rule of the element, as read in its source (summary only, the code is the partner's): the points inside the
initial radius are taken in element order up to NEIGHBOR_CNT; if fewer, the radius grows by 10 % at a time and points
of the extended list are added, again in list (element) order, until there are exactly NEIGHBOR_CNT.

    python ansys_usermat/apdl/nem_getnhgb_counts.py OUT.txt [--deck DECK] [--json results/2026-10-living/X.json]

Output: per point the count inside the initial radius, the enlargements, the final radius and, with the point
positions rebuilt from the deck (nem_neighbour_symmetry.build_points, checked against the printed offsets), the
number of points strictly inside the final radius, i.e. how many candidates of the last 10 % shell were left out
because of the element order.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def records(path):
    """Yield (kind, text) for every 'DEBUG getnhgb before/after' record."""
    buf, kind = [], None
    with open(path, errors="replace") as fh:
        for line in fh:
            if "DEBUG getnhgb" in line:
                if kind:
                    yield kind, " ".join(buf)
                kind = "before" if "before" in line else "after"
                buf = [line]
            elif kind:
                buf.append(line)
                if len(buf) > 3000:
                    yield kind, " ".join(buf)
                    kind, buf = None, []
    if kind:
        yield kind, " ".join(buf)


def take(text, label, n):
    i = text.find(label)
    if i < 0:
        return None
    toks = re.findall(NUM, text[i + len(label):])
    return [float(t.replace("D", "E").replace("d", "E")) for t in toks[:n]]


def parse(text, cap):
    sid = int(take(text, "sID:", 1)[0])
    nlive = int(take(text, "sNeighCntLive:", 1)[0])
    next_ = int(take(text, "sExtNeighCntLive:", 1)[0])
    r = dict(sid=sid, n=nlive, ext=next_)
    r["dx"] = take(text, "vdx:", nlive)
    r["dy"] = take(text, "vdy:", nlive)
    r["dz"] = take(text, "vdz:", nlive)
    col = take(text, "vLi_ColID:", cap + 1)
    r["col"] = [int(c) for c in col] if col else None
    rr = take(text, "Refinment of the search radius:", 1)
    if rr:
        r["rr"] = int(rr[0])
        r["radius"] = take(text, "sLdp_Radius:", 1)[0]
    return r


def emulate(deck, cap):
    """The element's rule rebuilt in Python on the deck's mesh (checked against the printed neighbour IDs).

    Initial radius 0.8 (3 V N / (4 pi n_points))^(1/3) with V the box of the surface nodes; candidates are the Gauss
    points and surface nodes of every element whose first Gauss point lies within 5 radii, in element order (each
    surface node once); inside the radius they are taken in that order up to N, the rest form the extended list; then
    the radius grows by 10 % until N are found, again in list order. Returns per point: neighbour IDs (1-based),
    final radius, number of points strictly inside the final radius."""
    sys.path.insert(0, str(HERE))
    from deck_components import read_mesh
    from nem_neighbour_symmetry import build_points
    nodes, elems = read_mesh(deck)
    eids, _, gp, _, surf = build_points(deck)
    P = np.vstack([gp, surf])
    ngp = len(gp)
    nid = sorted(nodes)
    surf_ids = [n for n in nid if np.isclose(np.abs(nodes[n]), 1.0).any()]
    sid_of_node = {n: ngp + k + 1 for k, n in enumerate(surf_ids)}
    box = surf.max(0) - surf.min(0)
    r0 = 0.8 * (3.0 * np.prod(box) * cap / (4.0 * np.pi * len(P))) ** (1.0 / 3.0)
    # candidate order per element: its 8 Gauss points, then its surface nodes
    first_gp = np.array([P[k * 8] for k in range(len(eids))])
    out = []
    for s in range(1, len(P) + 1):
        x0 = P[s - 1]
        near = np.flatnonzero(np.linalg.norm(first_gp - x0, axis=1) < 5.0 * r0)
        seen, order = set(), []
        for k in near:
            order += list(range(k * 8 + 1, k * 8 + 9))
            for n in elems[eids[k]][:8]:
                if n in sid_of_node and n not in seen:
                    seen.add(n)
                    order.append(sid_of_node[n])
        order = np.array([o for o in order if o != s])
        d = np.linalg.norm(P[order - 1] - x0, axis=1)
        inside = d < r0
        sel = list(order[inside][:cap])
        # extended list: every other candidate, in pass order (also those inside the radius but over the cap)
        chosen = set(sel)
        ext = [(o, dd) for o, dd in zip(order, d) if o not in chosen]
        R = r0
        taken = np.zeros(len(ext), bool)
        while len(sel) < cap:
            R *= 1.1
            for j, (o, dd) in enumerate(ext):
                if not taken[j] and dd < R and len(sel) < cap:
                    sel.append(o)
                    taken[j] = True
        d_all = np.linalg.norm(P - x0, axis=1)
        d_all[s - 1] = np.inf
        out.append(dict(sid=s, col=[int(c) for c in sel], radius=float(R), n_inside_final=int((d_all < R).sum())))
    return r0, ngp, out


def main_emulate(a):
    rows = []
    for cap in a.emulate:
        r0, ngp, out = emulate(a.deck, cap)
        tr = [o for o in out if o["n_inside_final"] > cap]
        gp_tr = sum(1 for o in tr if o["sid"] <= ngp)
        sn_tr = len(tr) - gp_tr
        nin = [o["n_inside_final"] for o in out if o["sid"] <= ngp]
        rows.append(dict(cap=cap, r0=r0, gauss_points_truncated=gp_tr, surface_nodes_truncated=sn_tr,
                         n_points=len(out), inside_final_gp=(min(nin), max(nin))))
        print(f"  cap {cap}: initial radius {r0:.4f} mm; last shell cut in element order at {gp_tr} of {ngp} Gauss "
              f"points and {sn_tr} of {len(out) - ngp} surface nodes; points inside the final radius {min(nin)}-{max(nin)}")
        if a.check_against and cap == a.cap:
            bef_aft = {}
            for kind, text in records(a.check_against):
                if kind == "after":
                    try:
                        rec = parse(text, cap)
                    except (TypeError, IndexError, ValueError):
                        continue
                    bef_aft.setdefault(rec["sid"], rec)
            same = sum(1 for o in out if o["sid"] in bef_aft and bef_aft[o["sid"]]["col"][1:1 + cap] == o["col"])
            print(f"    neighbour lists identical to the ANSYS output at {same} of {len(bef_aft)} points")
            rows[-1]["identical_to_ansys"] = [same, len(bef_aft)]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emulate", type=int, nargs="*", default=None,
                    help="caps for the Python rebuild of the rule (needs --deck); no ANSYS output needed")
    ap.add_argument("--check-against", default=None, help="ANSYS output (DEBUG_KEY_NEM = 100) to compare with")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--cap", type=int, default=30)
    ap.add_argument("--r0", type=float, default=None, help="initial radius (default: read from the output)")
    ap.add_argument("--deck", default=None, help="deck for the point positions (default: none, counts only)")
    ap.add_argument("--ngp", type=int, default=4096, help="number of Gauss points (8^3: 512 x 8)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    if a.emulate:
        rows = main_emulate(a)
        if a.json:
            Path(a.json).write_text(json.dumps(dict(deck=str(a.deck), emulate=rows), indent=1))
            print("wrote", a.json)
        return

    r0 = a.r0
    if r0 is None:
        with open(a.out, errors="replace") as fh:
            for line in fh:
                if "Estimated Search Radius" in line:
                    r0 = float(re.findall(NUM, line.split(":")[-1])[0])
                    break
    bef, aft = {}, {}
    for kind, text in records(a.out):
        try:
            rec = parse(text, a.cap)
        except (TypeError, IndexError, ValueError):
            continue
        (bef if kind == "before" else aft).setdefault(rec["sid"], rec)
    sids = sorted(set(bef) & set(aft))
    print(f"{a.out}: initial radius {r0}, cap {a.cap}; {len(bef)} before / {len(aft)} after records, {len(sids)} paired")

    P = None
    if a.deck:
        sys.path.insert(0, str(HERE))
        from nem_neighbour_symmetry import build_points
        _, _, gp, _, surf = build_points(a.deck)
        P = np.vstack([gp, surf])
    rows = []
    worst_check = 0.0
    for s in sids:
        b, f = bef[s], aft[s]
        dist = np.sqrt(np.array(f["dx"]) ** 2 + np.array(f["dy"]) ** 2 + np.array(f["dz"]) ** 2)
        row = dict(sid=s, n_initial=b["n"], ext=b["ext"], n_final=f["n"], enlargements=f["rr"],
                   radius_final=f["radius"], dmax_selected=float(dist.max()) if dist.size else 0.0)
        if P is not None:
            x0 = P[s - 1]
            ids = np.array(f["col"][1:1 + f["n"]]) - 1
            off = P[ids] - x0
            worst_check = max(worst_check, float(np.abs(off - np.column_stack([f["dx"], f["dy"], f["dz"]])).max()))
            d_all = np.linalg.norm(P - x0, axis=1)
            d_all[s - 1] = np.inf
            row["n_inside_initial"] = int((d_all < r0).sum())
            row["n_inside_final"] = int((d_all < f["radius"]).sum())
            # the shell added last: (radius/1.1, radius)
            prev = f["radius"] / 1.1 if f["rr"] > 0 else 0.0
            row["n_last_shell"] = int(((d_all >= prev) & (d_all < f["radius"])).sum())
            row["taken_last_shell"] = int(((dist >= prev) & (dist < f["radius"])).sum())
        rows.append(row)

    def stat(key, sel):
        v = [r[key] for r in rows if sel(r)]
        return (min(v), max(v)) if v else None

    gp_sel = lambda r: r["sid"] <= a.ngp
    sn_sel = lambda r: r["sid"] > a.ngp
    summ = {}
    for name, sel in (("gauss_points", gp_sel), ("surface_nodes", sn_sel)):
        summ[name] = {k: stat(k, sel) for k in rows[0] if k not in ("sid",)}
        summ[name]["count"] = sum(1 for r in rows if sel(r))
        if P is not None:
            summ[name]["points_truncated_in_last_shell"] = sum(
                1 for r in rows if sel(r) and r["n_inside_final"] > r["n_final"])
    for name, s in summ.items():
        print(f"  {name} ({s['count']}):")
        for k, v in s.items():
            if k != "count":
                print(f"    {k}: {v}")
    if P is not None:
        print(f"  offsets rebuilt from the deck vs printed: largest difference {worst_check:.1e}")
    if a.json:
        Path(a.json).write_text(json.dumps(dict(out=str(a.out), r0=r0, cap=a.cap, summary=summ,
                                                offset_check=worst_check if P is not None else None,
                                                points=rows), indent=1))
        print("wrote", a.json)


if __name__ == "__main__":
    main()
