"""ansys_seed_table.py -- mean phi per stored time and the seed means at T* = 1 of ANSYS runs (fields_<run>.json),
for the side-by-side tables with the Abaqus element (A-46).

    python ansys_seed_table.py run [run ...]
"""
import json
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent.parent / "ansys_usermat/apdl/results/2026-10-living"
for run in sys.argv[1:]:
    p = R / f"fields_{run}.json"
    if not p.exists():
        print(f"{run}: not yet")
        continue
    f = json.loads(p.read_text())
    s = np.isin(np.array(f["elem"]), f["seed"])
    print(f"{run}: mean phi {[round(float(np.mean(x)), 3) for x in f['phi']]}; seed L {np.array(f['L'][-1])[s].mean():.3f} "
          f"c {np.array(f['nut'][-1])[s].mean():.3f} alpha-1 {np.array(f['alpha'][-1])[s].mean():.3e} "
          f"vM {np.array(f['seqv'][-1])[s].mean() * 1e6:.3e} Pa")
