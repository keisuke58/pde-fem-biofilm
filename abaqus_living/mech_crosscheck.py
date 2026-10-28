"""mech_crosscheck.py -- the mechanics of the Abaqus element against the Python stress solver on the same field (A-51):
the last Gauss-point record of an Abaqus run (phi, alpha - 1 at the Gauss points) is averaged per element and solved
by wp3_stress_reference.solve_fields (hex8, selective reduced integration, stress at the centroid, corner constraints);
the seed means of von Mises are compared with the Abaqus record (mean over the Gauss points).

    python mech_crosscheck.py F:\\abaqus_work\\lv16_sep_ab_dt000125\\living_gp.txt 16
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat/apdl"))
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402

gp, n = Path(sys.argv[1]), int(sys.argv[2])
blocks, cur = [], []
for line in gp.read_text().splitlines():
    if line.startswith("# step"):
        if cur:
            blocks.append(np.array(cur))
        cur = []
    elif line and not line.startswith("#"):
        cur.append([float(x) for x in line.split()])
A = np.array(cur) if cur else blocks[-1]
h = 2.0 / n
i, j, k = (np.clip(np.floor((A[:, c] + 1.0) / h).astype(int), 0, n - 1) for c in (2, 3, 4))
cnt = np.zeros((n,) * 3)
np.add.at(cnt, (i, j, k), 1.0)


def cellmean(col):
    acc = np.zeros((n,) * 3); np.add.at(acc, (i, j, k), A[:, col]); return acc / np.maximum(cnt, 1)


phi, am1, vm_ab = cellmean(5), cellmean(8), cellmean(10)
phi2 = np.zeros((n,) * 3); np.add.at(phi2, (i, j, k), A[:, 5] ** 2); phi2 /= np.maximum(cnt, 1)
seed = W.seed_mask(n, W.seed_cells_8(HERE.parent / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json"))
BOX = ([((-1.0, -0.125, -0.75), (1.0, 0.125, -0.5)), ((-1.0, -0.125, 0.5), (1.0, 0.125, 0.75)),
        ((-0.125, -0.75, -1.0), (0.125, -0.5, 1.0)), ((-0.125, 0.5, -1.0), (0.125, 0.75, 1.0))] if "--obstacle" in sys.argv else None)
vm_py, p_py = S.solve_fields(n, phi, am1, obstacles=BOX)
vm_py2, _ = S.solve_fields(n, np.sqrt(phi2), am1, obstacles=BOX)        # stiffness from the Gauss-point mean of phi^2
seed = seed & np.isfinite(vm_py)
print(f"{gp}: seed von Mises [Pa]: Abaqus (Gauss-point mean) {vm_ab[seed].mean():.4e}; "
      f"Python solver on the Abaqus field: E((mean phi)^2) {vm_py[seed].mean():.4e}, E(mean phi^2) {vm_py2[seed].mean():.4e}")
