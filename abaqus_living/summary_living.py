"""summary_living.py -- the Gauss-point records of living_uel (living_gp.txt) per time: mean phi over the cube, and
over the seed elements (wp3 seed cells) phi, living, c, alpha - 1, von Mises (Pa), as the ANSYS / Python tables.

    python summary_living.py living_gp.txt n [--json out.json]
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ansys_usermat/apdl"))
import wp3_reference as W  # noqa: E402

path, n = Path(sys.argv[1]), int(sys.argv[2])
blocks, cur, t = [], [], None
for line in path.read_text().splitlines():
    if line.startswith("# step"):
        if cur:
            blocks.append((t, np.array(cur)))
        t = float(line.split()[-1]); cur = []
    elif line and not line.startswith("#"):
        cur.append([float(x) for x in line.split()])
if cur:
    blocks.append((t, np.array(cur)))
cells8 = W.seed_cells_8(HERE.parent / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
seed = W.seed_mask(n, cells8)
h = 2.0 / n
rows = []
for t, A in blocks:
    # the cell of each Gauss point from its position (element numbers skip the obstacle cells of --obstacle decks)
    i, j, k = (np.clip(np.floor((A[:, c] + 1.0) / h).astype(int), 0, n - 1) for c in (2, 3, 4))
    s = seed[i, j, k]
    m = {"t": t, "phi_mean": float(A[:, 5].mean()), "phi_max": float(A[:, 5].max()),
         "seed_phi": float(A[s, 5].mean()), "seed_living": float(A[s, 9].mean()), "seed_c": float(A[s, 6].mean()),
         "seed_a": float(A[s, 7].mean()), "seed_am1": float(A[s, 8].mean()), "seed_vm_Pa": float(A[s, 10].mean())}
    rows.append(m)
    print(f"T* {t:.4f}: phi mean {m['phi_mean']:.4f} max {m['phi_max']:.3f}; seed phi {m['seed_phi']:.3f} "
          f"L {m['seed_living']:.3f} c {m['seed_c']:.3f} a {m['seed_a']:.3f} alpha-1 {m['seed_am1']:.3e} "
          f"vM {m['seed_vm_Pa']:.3e} Pa")
if "--json" in sys.argv:
    Path(sys.argv[sys.argv.index("--json") + 1]).write_text(json.dumps(rows, indent=1))
