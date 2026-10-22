"""compare_obstacle.py -- Python (cells, figdata npz snapshots) against Abaqus (Gauss-point records, cell means) for the
obstacle columns at 16^3 (A-49): phi and c per region (seed, cells next to an obstacle face, the rest) at T* = 0.25 ... 1.

    python compare_obstacle.py <python npz> <abaqus living_gp.txt> 16
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat/apdl"))
import wp3_reference as W  # noqa: E402

npz, gp, n = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
h = 2.0 / n
ax = -1 + (np.arange(n) + 0.5) * h
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
obst = np.zeros((n,) * 3, bool)
for lo, hi in (((-1.0, -0.125, -0.75), (1.0, 0.125, -0.5)), ((-1.0, -0.125, 0.5), (1.0, 0.125, 0.75)),
               ((-0.125, -0.75, -1.0), (0.125, -0.5, 1.0)), ((-0.125, 0.5, -1.0), (0.125, 0.75, 1.0))):
    obst |= (X > lo[0]) & (X < hi[0]) & (Y > lo[1]) & (Y < hi[1]) & (Z > lo[2]) & (Z < hi[2])
seed = W.seed_mask(n, W.seed_cells_8(HERE.parent / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json"))
near = np.zeros_like(obst)
for a in range(3):
    for s in (1, -1):
        near |= np.roll(obst, s, a)
near &= ~obst
regs = (("seed", seed & ~obst), ("next to obstacle", near & ~seed), ("rest", ~obst & ~near & ~seed))
z = np.load(npz)
blocks, cur, t = [], [], None
for line in gp.read_text().splitlines():
    if line.startswith("# step"):
        if cur:
            blocks.append((t, np.array(cur)))
        t = float(line.split()[-1]); cur = []
    elif line and not line.startswith("#"):
        cur.append([float(x) for x in line.split()])
if cur:
    blocks.append((t, np.array(cur)))
for t, A in blocks:
    T = round(t, 2)
    if f"phi_{T}" not in z.files:
        continue
    i, j, k = (np.clip(np.floor((A[:, c] + 1.0) / h).astype(int), 0, n - 1) for c in (2, 3, 4))
    acc = {f: np.zeros((n,) * 3) for f in ("phi", "c")}; cnt = np.zeros((n,) * 3)
    np.add.at(cnt, (i, j, k), 1.0)
    np.add.at(acc["phi"], (i, j, k), A[:, 5]); np.add.at(acc["c"], (i, j, k), A[:, 6])
    ab = {f: acc[f] / np.maximum(cnt, 1) for f in acc}
    py = {"phi": z[f"phi_{T}"], "c": z[f"c_{T}"]}
    print(f"T* {T}: " + "; ".join(f"{nm} phi Py {py['phi'][m].mean():.3f} Ab {ab['phi'][m].mean():.3f}, "
                                  f"c Py {py['c'][m].mean():.3f} Ab {ab['c'][m].mean():.3f}" for nm, m in regs))
