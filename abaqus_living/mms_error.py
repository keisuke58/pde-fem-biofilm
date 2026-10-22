"""mms_error.py -- L2 and max error of the Abaqus element against the manufactured solution (mms.py) at the last
Gauss-point record of living_gp.txt (A-50). Equal weights: the meshes are uniform.

    python mms_error.py F:\\abaqus_work\\mms16\\living_gp.txt [...]
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mms as M  # noqa: E402

for p in sys.argv[1:]:
    blocks, cur, t = [], [], None
    for line in Path(p).read_text().splitlines():
        if line.startswith("# step"):
            if cur:
                blocks.append((t, np.array(cur)))
            t = float(line.split()[-1]); cur = []
        elif line and not line.startswith("#"):
            cur.append([float(x) for x in line.split()])
    if cur:
        blocks.append((t, np.array(cur)))
    t, A = blocks[-1]
    T = t - 1e-8                                   # total time minus the 1e-8 start step = step time
    x, y, z = A[:, 2], A[:, 3], A[:, 4]
    ep = A[:, 5] - M.phi_m(x, y, z, T); ec = A[:, 6] - M.c_m(x, y, z, T)
    n = round((len(A) / 8) ** (1 / 3))
    print(f"{p}: {n}^3, T* {T:.4g}: L2 error phi {np.sqrt((ep ** 2).mean()):.3e} (max {np.abs(ep).max():.3e}), "
          f"c {np.sqrt((ec ** 2).mean()):.3e} (max {np.abs(ec).max():.3e})")
