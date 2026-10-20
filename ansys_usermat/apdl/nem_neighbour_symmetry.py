#!/usr/bin/env python3
"""nem_neighbour_symmetry.py -- does the neighbour search of the NEM operator break the mirror symmetry? (A-23)

IKMHIWI03 found (A-23) that results of mirror-symmetric set-ups are not symmetric, already with only the diffusion
of phi switched on, and traced it to the neighbour search of the element (GetNhgb): it goes through the points in
element order and keeps the first NEIGHBOR_CNT points found inside the search radius. Near the surface, where the
surface nodes are points too, the radius holds more points than the cap, so which ones are dropped depends on the
numbering. This script rebuilds that operator in Python on the mesh of the 8^3 decks (element and node numbering
from the deck in this folder, the same mesh as the coupling runs) and marches the diffusion of phi from the seed:

    phi_dot = beta Lap(phi),  beta = 0.02 mm^2/T*,  T* = 0.25,  phi = 1 in the 32 seed elements, 0 elsewhere,

for the rule of the element ("first": first N found in element order) with N = 30, 60, 120, and for a
distance-sorted search ("sorted": the N nearest, a distance shell that straddles N taken whole). Output: the largest
difference of the element-mean phi between mirrored elements (x -> -x, z -> -z), as check_living_fields.py reports
for the ANSYS runs (sym_a0: N = 30 x 8.5e-3, z 6.4e-3; 60: 2.5e-3, 4.2e-3; 120: 4.6e-4, 2.0e-3).

Assumptions (the element's code is not in this repository): Gauss points of an element at centroid + (node -
centroid)/sqrt(3) in the element's node order; the points are the Gauss points in element order followed by the
surface nodes in node order; surface-node values are the inverse-distance-squared mean of the Gauss points within
the radius (zero normal gradient, approximately); weighted least squares of second order with the Gaussian weight
exp(-(4 r / beta*)^2 / 2), beta* = 0.2 mm, threshold 1e-8; search radius 0.187 mm (the element's estimate on 8^3,
A-23) and, since that radius holds at most 16 points here (so the cap is never reached), also 0.25 and 0.3 mm.
The result is meant to show the mechanism and its size, not to reproduce ANSYS digit by digit.

Result (21 Oct): radius 0.187 mm: symmetric to 4e-15 for every rule (at most 16 points, the cap is not reached).
Radius 0.25 mm (up to 32 points): first 30 gives x 2.3e-3, z 5.9e-3 (ANSYS sym_a0 8.5e-3, 6.4e-3); first 120 and
sorted 30 give 4e-15. Radius 0.3 mm (up to 62 points): first 30 gives 0.12 and 0.22. With surface points placed at
the projections of the boundary Gauss points instead of the nodes the numbers are the same. So the truncation in
element order alone produces an asymmetry of the size seen in ANSYS, and either a cap above the point count or a
distance-sorted search removes it to round-off. ANSYS with 120 still shows 4.6e-4 and 2.0e-3: either its radius
holds more than 120 points, or there is a second source (to be settled on IKMHIWI03, addendum 15).

    python ansys_usermat/apdl/nem_neighbour_symmetry.py [--deck DECK] [--radius 0.187]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from deck_components import centroids, read_mesh  # noqa: E402

BETA, T_END, BS, THR = 0.02, 0.25, 0.2, 1e-8


def build_points(deck):
    nodes, elems = read_mesh(deck)
    cen = centroids(nodes, elems)
    eids = sorted(elems)
    gp, gp_elem = [], []
    for e in eids:
        c = np.asarray(cen[e])
        for n in elems[e][:8]:
            gp.append(c + (np.asarray(nodes[n]) - c) / np.sqrt(3.0))
            gp_elem.append(e)
    gp = np.array(gp)
    nid = sorted(nodes)
    surf = np.array([nodes[n] for n in nid if np.isclose(np.abs(nodes[n]), 1.0).any()], float)
    return eids, cen, gp, np.array(gp_elem), surf


def wls_row(D, nneigh_used):
    r = np.linalg.norm(D, axis=1)
    w = np.exp(-0.5 * (4.0 * r / BS) ** 2)
    w[w <= THR] = 0.0
    dx, dy, dz = D.T
    A = np.column_stack([0.5 * dx**2, 0.5 * dy**2, 0.5 * dz**2, dx * dy, dx * dz, dy * dz, dx, dy, dz])
    M = A.T @ (w[:, None] * A)
    Dm = np.linalg.lstsq(M, A.T * w, rcond=None)[0]
    return Dm[0] + Dm[1] + Dm[2]


def operator(gp, surf, radius, cap, rule):
    """Laplacian at the Gauss points: L = Lg (on Gauss points) + Ls (on surface nodes) with surface = S gp."""
    P = np.vstack([gp, surf])
    ng = len(gp)
    rows, cols, vals = [], [], []
    count_max = 0
    for i in range(ng):
        d = np.linalg.norm(P - gp[i], axis=1)
        inside = np.flatnonzero((d > 1e-12) & (d <= radius * (1 + 1e-12)))
        count_max = max(count_max, len(inside))
        if rule == "first":
            sel = inside[:cap]                         # element order, then node order
        else:
            o = inside[np.argsort(d[inside], kind="stable")]
            if len(o) > cap:
                rc = d[o[cap - 1]]
                o = o[d[o] <= rc * (1 + 1e-9)]          # a shell straddling the cap is taken whole
            sel = o
        c = wls_row(P[sel] - gp[i], len(sel))
        rows += [i] * len(sel) + [i]
        cols += list(sel) + [i]
        vals += list(c) + [-c.sum()]
    L = sp.csr_matrix((vals, (rows, cols)), shape=(ng, len(P)))
    # surface nodes: inverse-distance-squared mean of the Gauss points within the radius
    srows, scols, svals = [], [], []
    for j, s in enumerate(surf):
        d = np.linalg.norm(gp - s, axis=1)
        nb = np.flatnonzero(d <= radius)
        w = 1.0 / d[nb] ** 2
        srows += [j] * len(nb); scols += list(nb); svals += list(w / w.sum())
    S = sp.csr_matrix((svals, (srows, scols)), shape=(len(surf), ng))
    return (L[:, :ng] + L[:, ng:] @ S).tocsr(), count_max


def mirror(cen, eids, m):
    key = {tuple(np.round(cen[e], 6)): k for k, e in enumerate(eids)}
    return np.array([key[tuple(np.round(np.asarray(cen[e]) * m, 6))] for e in eids])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deck", default=str(HERE / "ds_oliver_wired_n5trial_inter34only_baseline_at0p0006.dat"))
    ap.add_argument("--radius", type=float, nargs="+", default=[0.187, 0.25, 0.3])
    ap.add_argument("--seed-from", default=str(HERE / "results/2026-10-living/fields_v2f_8_mix_ab.json"))
    a = ap.parse_args()
    eids, cen, gp, gp_elem, surf = build_points(a.deck)
    seed = set(json.loads(Path(a.seed_from).read_text())["seed"])
    phi0 = np.array([1.0 if e in seed else 0.0 for e in gp_elem])
    ix, iz = mirror(cen, eids, np.array([-1, 1, 1])), mirror(cen, eids, np.array([1, 1, -1]))
    out = {}
    print(f"{len(gp)} Gauss points, {len(surf)} surface nodes")
    for radius, rule, cap in [(r, ru, c) for r in a.radius
                              for ru, c in [("first", 30), ("first", 60), ("first", 120), ("sorted", 30)]]:
        L, nmax = operator(gp, surf, radius, cap, rule)
        lam = np.abs(sp.linalg.eigs(L, k=1, which="LM", return_eigenvectors=False)[0])
        nstep = int(np.ceil(T_END * BETA * lam / 0.5)) + 1
        dt = T_END / nstep
        phi = phi0.copy()
        for _ in range(nstep):
            phi = phi + dt * BETA * (L @ phi)
        pe = phi.reshape(-1, 8).mean(1)
        dx, dz = float(np.abs(pe - pe[ix]).max()), float(np.abs(pe - pe[iz]).max())
        key = f"r{radius}_{rule}_{cap}"
        out[key] = dict(radius=radius, rule=rule, cap=cap, x=dx, z=dz, points_in_radius_max=nmax, steps=nstep,
                        seed_mean=float(pe[[k for k, e in enumerate(eids) if e in seed]].mean()))
        print(f"  radius {radius:.3f} {rule:6s} N = {cap:3d}: mirror difference x {dx:.1e}, z {dz:.1e}; "
              f"points in the radius up to {nmax}; seed mean phi {out[key]['seed_mean']:.4f}", flush=True)
    p = HERE / "results/2026-10-living/nem_neighbour_symmetry.json"
    p.write_text(json.dumps(dict(beta=BETA, t_end=T_END, runs=out,
                                 ansys_sym_a0=dict(first_30=[8.5e-3, 6.4e-3], first_60=[2.5e-3, 4.2e-3],
                                                   first_120=[4.6e-4, 2.0e-3])), indent=1))
    print("wrote", p)


if __name__ == "__main__":
    main()
