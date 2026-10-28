"""A-71: the antibiotic in time (--ab-schedule const|pulse|sine, the same dose after T* = 1/3) on 64^3, TVD, mirrored
walls: seed sum phibar, share of A, alpha - 1 and von Mises at T* = 1, and the lowest sum phibar after the onset.

  python ab_schedule_a71.py
"""
import glob
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parents[1] / "results" / "2026-10-3d"
print(f"{'schedule':10s} {'start':5s} {'sum phibar':>10s} {'min after':>10s} {'share A':>8s} {'alpha-1':>10s} {'vM [Pa]':>10s}")
for ex in ("separate_ab", "mixed_ab"):
    for s in ("const", "pulse", "sine"):
        f = glob.glob(str(RES / f"*nf64*_ab{s}0.2222*_{ex}.json"))
        if not f:
            continue
        r = json.load(open(f[0]))["runs"][ex]
        t, L = np.array(r["t"]), np.array(r["living_seed"])
        print(f"{s:10s} {ex[:3]:5s} {L[-1]:10.3f} {L[t > 1 / 3].min():10.3f} {r['share_a_seed'][-1]:8.3f} "
              f"{r['seed_am1']:10.3e} {r['seed_vm_Pa']:10.3e}")
