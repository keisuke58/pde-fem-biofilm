"""A-87: mesh error of the mechanics alone. The coupled runs solve the stress on a 16^3 hex8 mesh whatever the field
grid (block means of phi_3D and alpha - 1). Here the end state of a 128^3 run (phi, alpha saved in the snapshot
file) is put on 16^3, 32^3 and 64^3 stress meshes, and the seed mean of the von Mises stress is compared.

    python stress_mesh_a87.py <figdata npz of a 128^3 run with alpha_1.0>
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402

Z = np.load(sys.argv[1])
phi, alpha = Z["phi_1.0"], Z["alpha_1.0"]
cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
for nm in (16, 32, 64):
    vm, p = S.solve_fields(nm, W.block_mean(phi, nm), W.block_mean(alpha - 1, nm))
    seed = W.seed_mask(nm, cells8)
    print(f"stress mesh {nm}^3: seed von Mises {vm[seed].mean():.4e} Pa, seed mean stress {p[seed].mean():.4e} Pa, "
          f"max vM {np.nanmax(vm):.3e}", flush=True)
