#!/usr/bin/env python3
"""Mesh study of the seed stress (Check 3 of the 5 Oct slides), in Python.

The ANSYS model is the partner's 2 mm cube with 8 x 8 x 8 = 512 hexahedra.
Its stresses have not been checked against a finer mesh. This script solves
the same mechanical problem with a small-strain linear-elastic hex8 solver on
8^3, 16^3 and 32^3 meshes:
  - the seed = the 32 elements of BIOFILM1 in the partner's deck, as a region
    in space (a fine element belongs to it when its centroid lies in a seed
    element of the 8^3 mesh);
  - growth: eigenstrain (alpha - 1) I = 1.1e-3 I in the seed, 0 elsewhere
    (Check 3: k_alpha phi t at T* = 1.1);
  - stiffness E (phi^2 + f), phi = 1 in the seed, 0 elsewhere, f = 1e-3
    (Klempt 2024 Eq. 20 weighting, the floor of the runs), E = 10 Pa, nu = 0.49,
    selective reduced integration (volumetric part at the centre point);
  - only rigid-body motion constrained, at the three corner nodes of the deck.
Reported per mesh: von Mises and mean stress in the seed (mean and max) and
the largest von Mises outside it. The quantity of Check 3, "neighbours carry
about 25x the seed's von Mises", is the ratio of the last to the first.

    python ansys_usermat/apdl/mesh_study_seed.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import model_setup_fig as M  # noqa: E402  (reads the partner's deck)

E, NU, FLOOR, GROWTH = 10.0, 0.49, 1e-3, 1.1e-3
CORNERS = [(-1, -1, -1), (1, -1, -1), (-1, 1, -1)]


def element_matrices(h, nu=NU, full=False):
    """full=True: volumetric part with 2x2x2 points too (no B-bar), for the locking check."""
    g = np.array([-1, 1]) / np.sqrt(3)
    corners = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                        [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)

    def bmat(xi):
        dN = np.array([[cx * (1 + xi[1] * cy) * (1 + xi[2] * cz),
                        cy * (1 + xi[0] * cx) * (1 + xi[2] * cz),
                        cz * (1 + xi[0] * cx) * (1 + xi[1] * cy)] for cx, cy, cz in corners]) / 8 * 2 / h
        b = np.zeros((6, 24))
        for a in range(8):
            dx, dy, dz = dN[a]
            i = 3 * a
            b[0, i] = dx; b[1, i + 1] = dy; b[2, i + 2] = dz
            b[3, i] = dy; b[3, i + 1] = dx
            b[4, i + 1] = dz; b[4, i + 2] = dy
            b[5, i] = dz; b[5, i + 2] = dx
        return b
    lam = nu / ((1 + nu) * (1 - 2 * nu))
    mu = 1 / (2 * (1 + nu))
    Dd = np.diag([2 * mu] * 3 + [mu] * 3)
    Dv = np.zeros((6, 6)); Dv[:3, :3] = lam
    m = np.array([1, 1, 1, 0, 0, 0.])
    Kd = np.zeros((24, 24)); fd = np.zeros(24)
    for a in g:
        for b_ in g:
            for c in g:
                B = bmat((a, b_, c))
                Kd += B.T @ Dd @ B * (h / 2) ** 3
                fd += B.T @ Dd @ m * (h / 2) ** 3
                if full:
                    Kd += B.T @ Dv @ B * (h / 2) ** 3
                    fd += B.T @ Dv @ m * (h / 2) ** 3
    B0 = bmat((0, 0, 0))
    Kv = 0 if full else B0.T @ Dv @ B0 * h ** 3
    fv = 0 if full else B0.T @ Dv @ m * h ** 3
    return Kd + Kv, fd + fv, B0, Dd + Dv, m


def seed_boxes():
    xyz, cen, seed_ids, _, _ = M.read_deck(M.DECK)
    return np.array([cen[e] for e in seed_ids])            # centroids of the 8^3 seed elements


def solve(n, boxes, nu=NU, full=False):
    h = 2.0 / n
    Ke0, fe0, B0, D, m = element_matrices(h, nu, full)
    ne = n
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + ne, j:j + ne, k:k + ne].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    ax = -1 + h * (np.arange(ne) + 0.5)
    cx, cy, cz = np.meshgrid(ax, ax, ax, indexing="ij")
    cen = np.stack([cx.ravel(), cy.ravel(), cz.ravel()], 1)
    inseed = np.zeros(len(cen), bool)
    for b in boxes:                                        # 0.25 mm boxes of the 8^3 seed
        inseed |= np.all(np.abs(cen - b) < 0.125 + 1e-9, axis=1)
    Ee = E * np.where(inseed, 1.0 + FLOOR, FLOOR)
    eps0 = np.where(inseed, GROWTH, 0.0)
    rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(3 * (n + 1) ** 3,) * 2)
    F = np.zeros(3 * (n + 1) ** 3)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    node = lambda p: idx[tuple(int(round((c + 1) / h)) for c in p)]
    fix = [3 * node(CORNERS[0]) + d for d in (0, 1, 2)] + [3 * node(CORNERS[1]) + d for d in (1, 2)] \
        + [3 * node(CORNERS[2]) + 2]
    free = np.setdiff1d(np.arange(K.shape[0]), fix)
    u = np.zeros(K.shape[0])
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = sig[:, :3].mean(1)
    s = sig[:, :3] - p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return inseed, vm, p


def main():
    boxes = seed_boxes()
    print(f"seed = {len(boxes)} elements of the 8^3 deck; growth {GROWTH:g}, E = {E:g} Pa, nu = {NU}")
    print(f"{'mesh':>6} {'elements':>9} {'seed vM mean':>13} {'seed vM max':>12} {'seed p mean':>12}"
          f" {'outside vM max':>15} {'ratio':>6}   [Pa]")
    for n in (8, 16, 32):
        t0 = time.time()
        ins, vm, p = solve(n, boxes)
        print(f"{n:>4}^3 {n ** 3:>9} {vm[ins].mean():13.3e} {vm[ins].max():12.3e} {p[ins].mean():12.3e}"
              f" {vm[~ins].max():15.3e} {vm[~ins].max() / vm[ins].mean():6.1f}   ({time.time() - t0:.0f} s)",
              flush=True)


if __name__ == "__main__":
    main()
