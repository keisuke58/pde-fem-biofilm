"""initial_seed_check.py -- initial phi of the deck as the element means of the nodal values (what the Galerkin element
starts from) against the cell values 1 / 0 of the Python reference and the ANSYS Gauss points (A-46).

    python initial_seed_check.py deck.inp n
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat/apdl"))
import wp3_reference as W  # noqa: E402

deck, n = Path(sys.argv[1]), int(sys.argv[2])
init = deck.read_text().split("*STEP, NAME=INIT")[1].split("*END STEP")[0]
vals = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"^(\d+), 11, 11, ([0-9.eE+-]+)$", init, re.M)}
nid = lambda i, j, k: 1 + i + (n + 1) * (j + (n + 1) * k)  # noqa: E731
off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
seed = W.seed_mask(n, W.seed_cells_8(HERE.parent / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json"))
em = np.array([[[np.mean([vals[nid(i + a, j + b, k + c)] for a, b, c in off]) for k in range(n)]
                for j in range(n)] for i in range(n)])
ext = [int(np.ptp(np.argwhere(seed)[:, d]) + 1) for d in range(3)]
print(f"seed: {int(seed.sum())} cells, extent {ext} cells")
print(f"initial phi, mean over the seed elements: Galerkin (nodal) {em[seed].mean():.3f}; cells / Gauss points 1.000")
print(f"initial phi, cube mean: Galerkin {em.mean():.4f}")
