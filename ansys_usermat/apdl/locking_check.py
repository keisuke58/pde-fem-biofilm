#!/usr/bin/env python3
"""Volumetric locking check for the seed problem (incompressibility, nu = 0.49).

Klempt et al. 2024 use nu = 0.49, close to incompressible. Linear hexahedra
with full 2x2x2 integration lock in that limit: the volumetric part is too
stiff and the stresses are wrong. ANSYS SOLID185 with the default KEYOPT(2) = 0
(the partner's deck sets no KEYOPT) uses B-bar, i.e. the volumetric part at
the centre point, as mesh_study_seed.py does. This script solves the seed
problem of mesh_study_seed.py with both integrations for several nu and two
meshes, so that one can see

  - whether full integration would change the result at nu = 0.49 (locking),
  - whether B-bar stays stable as nu -> 0.5.

    python ansys_usermat/apdl/locking_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_study_seed as S  # noqa: E402

NUS = [0.3, 0.45, 0.49, 0.499]
MESHES = [8, 16]


def main():
    boxes = S.seed_boxes()
    print(f"seed {len(boxes)} elements, growth {S.GROWTH:g}, E = {S.E:g} Pa; stresses in Pa")
    print(f"{'mesh':>5} {'nu':>6} {'integration':>12} {'seed vM mean':>13} {'seed p mean':>12}"
          f" {'outside vM max':>15}")
    rows = []
    for n in MESHES:
        for nu in NUS:
            for full in (False, True):
                ins, vm, p = S.solve(n, boxes, nu=nu, full=full)
                r = (n, nu, "full 2x2x2" if full else "B-bar", vm[ins].mean(), p[ins].mean(),
                     vm[~ins].max())
                rows.append(r)
                print(f"{n:>3}^3 {nu:6.3f} {r[2]:>12} {r[3]:13.3e} {r[4]:12.3e} {r[5]:15.3e}",
                      flush=True)
    return rows


if __name__ == "__main__":
    main()
