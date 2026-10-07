#!/usr/bin/env python3
"""Two-way coupling, step 3 of ROADMAP_TWO_WAY.md, estimate in Python:
the composition acts on the stiffness.

E = (phi^2 + f) * (chi_1 E_1 + chi_2 E_2), chi_i = phi_i/(phi_1+phi_2), on the
partner's 8^3 model with the seed problem of apdl/mesh_study_seed.py (growth
1.1e-3 in the 32 seed elements, nu = 0.49, B-bar). The composition in the
seed is the one of composition_local_nutrient.py (local nutrient in the point
model, case 6, Thiele number 4), the only result so far in which the
composition varies across the seed. Species-specific moduli are not known
(Boel et al. 2013); E_2/E_1 is a sensitivity parameter, with the mean
modulus kept at 10 Pa.

    python ansys_usermat/two_way_step3.py

Result (4 Oct 2026): seed composition chi_1 = 0.12-0.44; seed averages
  E2/E1 = 1: von Mises 5.703e-5 Pa, mean stress -3.768e-5 Pa
  E2/E1 = 2: +0.1 % / +0.1 %
  E2/E1 = 5: +0.2 % / +0.1 %
The seed's stress hardly depends on its own stiffness: the seed is about 1000
times stiffer than the void around it, so the soft surroundings set the
stress (a stiff inclusion in a soft matrix). Composition-dependent stiffness
matters only where the biofilm itself carries the load, e.g. a contiguous
biofilm on a substrate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "apdl"))
import composition_local_nutrient as L  # noqa: E402
import mesh_study_seed as S  # noqa: E402

RATIOS = [1.0, 2.0, 5.0]


def solve(Ee, ins):
    n = 8; h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
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
    sig = Ee[:, None] * (eps @ D.T)
    p = sig[:, :3].mean(1); s = sig[:, :3] - p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm, p


def main():
    Sm, Dm = L.grid_masks()
    c = L.nutrient(Sm, Dm, 4.0)
    chi = np.full((8,) * 3, 0.5)
    chi[Sm] = L.shares("2sp_case6", c[Sm])
    # element order of mesh_study_seed: i (x) slowest, k (z) fastest -> ravel C order
    ins = Sm.ravel(); chi_e = chi.ravel()
    print(f"seed composition chi_1: {chi_e[ins].min():.3f} .. {chi_e[ins].max():.3f}")
    base = None
    for r in RATIOS:
        e1 = 2 * S.E / (1 + r); e2 = r * e1                  # mean (e1 + e2)/2 = E
        Eb = chi_e * e1 + (1 - chi_e) * e2
        Ee = np.where(ins, (1 + S.FLOOR) * Eb, S.FLOOR * S.E)
        vm, p = solve(Ee, ins)
        out = (vm[ins].mean(), p[ins].mean(), vm[ins].std() / vm[ins].mean())
        if base is None:
            base = out
        print(f"E2/E1 = {r:g}: seed vM mean {out[0]:.3e} Pa ({out[0] / base[0] - 1:+.1%}), "
              f"p mean {out[1]:.3e} Pa ({out[1] / base[1] - 1:+.1%}), vM spread {out[2]:.2f}")


if __name__ == "__main__":
    main()
