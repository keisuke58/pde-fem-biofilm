#!/usr/bin/env python3
"""Sign of the stress around the growing seed, predicted in Python.

Klempt et al. 2024 (and PAMM 2023) report pressure inside the biofilm and a
ring of tension at its edge. The ANSYS runs have not yet output the mean
stress of the seed's neighbours (open item, read on IKMHIWI03). This script
uses the seed problem of mesh_study_seed.py (B-bar, E = 10 Pa, nu = 0.49,
growth 1.1e-3 in the 32 seed elements, void floor 1e-3) and reports, for the
elements that share a face with the seed (first layer) and the next layer:
the mean stress p, the largest principal stress and the hoop-like stress
tangential to the seed surface, on 8^3 and 16^3 (layers of the 8^3 size).

    python ansys_usermat/apdl/neighbour_sign.py
    python ansys_usermat/apdl/neighbour_sign.py --csv all_stress_<job>.csv [--grid 8]

With --csv the same layer table is made from an ANSYS run (the CSV of
post_all_stress.mac, MPa -> Pa): the seed is every element whose alpha column
exceeds half its maximum. A CSV without centroid columns needs --grid 8
(regular 2 mm block, element numbering as plot_3d.py assumes; check its OK). That CSV has no shear stresses, so only the mean
stress p is compared, not the hoop and radial parts.

Result (3 Oct, Python): the seed is compressed (p = -3.8e-5 Pa on 8^3), the
first layer around it is in tension on average (p = +5.8e-6 Pa, 56 % of its
elements; 16^3: +4.1e-6, 70 %), with tension along the seed surface (hoop
+9.2e-6) and slight compression across it (radial -1e-6). This is the ring
of tension of Klempt et al. 2024 and PAMM 2023, in mean, not in every
element: the elements at the seed's corners and edges carry both signs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_study_seed as S  # noqa: E402


def stresses(n, boxes):
    """centroid stress tensors (Voigt) and the seed mask, re-solving as in S.solve."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    ax = -1 + h * (np.arange(n) + 0.5)
    cen = np.stack([c.ravel() for c in np.meshgrid(ax, ax, ax, indexing="ij")], 1)
    ins = np.zeros(len(cen), bool)
    for b in boxes:
        ins |= np.all(np.abs(cen - b) < 0.125 + 1e-9, axis=1)
    Ee = S.E * np.where(ins, 1.0 + S.FLOOR, S.FLOOR)
    eps0 = np.where(ins, S.GROWTH, 0.0)
    rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(3 * (n + 1) ** 3,) * 2)
    F = np.zeros(3 * (n + 1) ** 3)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    node = lambda p: idx[tuple(int(round((c + 1) / h)) for c in p)]
    fix = [3 * node(S.CORNERS[0]) + d for d in (0, 1, 2)] + [3 * node(S.CORNERS[1]) + d for d in (1, 2)] \
        + [3 * node(S.CORNERS[2]) + 2]
    free = np.setdiff1d(np.arange(K.shape[0]), fix)
    u = np.zeros(K.shape[0])
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    return cen, ins, Ee[:, None] * (eps @ D.T)


def tensor(v):
    return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]])


def from_csv(path, grid=None, size=2.0, deck=None):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from plot_3d import read
    col = read(path, grid, size, deck)
    cen = np.stack([col["cx"], col["cy"], col["cz"]], 1)
    ins = col["alpha"] > 0.5 * col["alpha"].max()
    h = np.min(np.diff(np.unique(np.round(cen[:, 0], 9))))
    d = np.full(len(cen), np.inf)
    for b in cen[ins]:
        d = np.minimum(d, np.max(np.maximum(np.abs(cen - b) - h / 2, 0), axis=1))
    layer = np.ceil(d / 0.25 - 1e-9).astype(int)
    p = col["p"] * 1e6
    print(f"{path}: {len(cen)} elements, seed {ins.sum()} (alpha > half its max); p in Pa")
    print(f"{'region':>10} {'elements':>8} {'p mean':>10} {'p>0 share':>9}")
    for name, sel in (("seed", ins), ("layer 1", (~ins) & (layer == 1)),
                      ("layer 2", (~ins) & (layer == 2)), ("farther", (~ins) & (layer >= 3))):
        print(f"{name:>10} {sel.sum():8d} {p[sel].mean():10.2e} {np.mean(p[sel] > 0):9.2f}")


def main():
    if "--csv" in sys.argv:
        g = int(sys.argv[sys.argv.index("--grid") + 1]) if "--grid" in sys.argv else None
        dk = sys.argv[sys.argv.index("--deck") + 1] if "--deck" in sys.argv else None
        return from_csv(sys.argv[sys.argv.index("--csv") + 1], g, deck=dk)
    boxes = S.seed_boxes()
    print("stresses in Pa; layer = distance from the seed in 8^3 element sizes (0.25 mm)")
    for n in (8, 16):
        cen, ins, sig = stresses(n, boxes)
        # distance (max-norm, in units of 0.25 mm) from each element centroid to the seed region
        d = np.full(len(cen), np.inf)
        for b in boxes:
            d = np.minimum(d, np.max(np.maximum(np.abs(cen - b) - 0.125, 0), axis=1))
        layer = np.ceil(d / 0.25 - 1e-9).astype(int)
        seedc = cen[ins].mean(0)
        print(f"\n{n}^3   {'region':>16} {'elements':>8} {'p mean':>10} {'p>0 share':>9}"
              f" {'s1 mean':>10} {'hoop mean':>10} {'radial mean':>11}")
        for name, sel in (("seed", ins), ("layer 1", (~ins) & (layer == 1)),
                          ("layer 2", (~ins) & (layer == 2)), ("farther", (~ins) & (layer >= 3))):
            p = sig[sel, :3].mean(1)
            s1, hoop, rad = [], [], []
            for c, v in zip(cen[sel], sig[sel]):
                T = tensor(v)
                s1.append(np.linalg.eigvalsh(T)[-1])
                r = c - seedc
                r = r / max(np.linalg.norm(r), 1e-12)
                rad.append(r @ T @ r)
                hoop.append((np.trace(T) - r @ T @ r) / 2)
            print(f"{'':6}{name:>16} {sel.sum():8d} {p.mean():10.2e} {np.mean(p > 0):9.2f}"
                  f" {np.mean(s1):10.2e} {np.mean(hoop):10.2e} {np.mean(rad):11.2e}")


if __name__ == "__main__":
    main()
