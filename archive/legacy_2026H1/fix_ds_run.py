#!/usr/bin/env python3
"""Run the DS (dysbiotic_static) tooth UMAT with the CORRECTED composition
(raw CLSM, V.dispar-dominant) vs the buggy So-dominant one. Exact sigma."""
import sys, subprocess
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "JAXFEM"))
import propagate_posterior_umat as P  # reuses make_inp/launch/extract, N_STEPS=2000

OUT = P.OUT
# buggy (current) vs corrected (raw CLSM dysbiotic-static, representative D10)
CASES = {
    "DS_buggy_Sodom":   np.array([0.944, 0.011, 0.011, 0.011, 0.011]),
    "DS_corrected_D6":  np.array([0.064, 0.054, 0.585, 0.107, 0.190]),
    "DS_corrected_D10": np.array([0.036, 0.057, 0.568, 0.129, 0.209]),
    "DS_corrected_D21": np.array([0.035, 0.130, 0.368, 0.203, 0.265]),
}
items = [("dysbiotic_static", phi / phi.sum(), tag) for tag, phi in CASES.items()]
sigma, meta = P.run_batch(items, n_par=4)
print("\n=== DS composition fix: max von Mises ===")
for tag in CASES:
    s = sigma.get(tag)
    print("  %-18s sigma=%s kPa  (E_voigt=%.3e, alpha_max=%.3f)"
          % (tag, ("%.2f" % (s * 1e3)) if s else "FAIL",
             float(CASES[tag] / CASES[tag].sum() @ np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])),
             meta[tag]["alpha_max"]))
print("\nheadline MAP_SIGMA[DS] (buggy) = 13.63 kPa")
