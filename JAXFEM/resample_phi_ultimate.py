#!/usr/bin/env python3
"""resample_phi_ultimate.py
============================
Re-run 0D ODE on ultimate_10000p TMCMC posterior (K=50 uniform draws).

Fix for C1: posterior_ci_0d.py used a wrong path (not ultimate_10000p)
and a basin filter that dropped 96% of CS samples.
This script uses the correct ultimate_10000p posterior with NO basin filter.

Saves: _ci_0d_results/{cond}/samples_0d_ultimate.json
       (same schema as samples_0d.json — plug-in replacement for CI script)

Usage
-----
  python JAXFEM/resample_phi_ultimate.py
  python JAXFEM/resample_phi_ultimate.py --n-samples 100 --seed 0
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE   = Path(__file__).resolve().parent
_FEM    = _HERE.parent
_ULTIMATE = Path("/home/nishioka/IKM_Hiwi/Tmcmc202601/data_5species/_runs/ultimate_10000p")
_CI0D   = _FEM / "_ci_0d_results"

sys.path.insert(0, str(_FEM))   # for material_models, JAXFEM

CONDITIONS = ["commensal_hobic", "dysbiotic_hobic", "commensal_static", "dysbiotic_static"]
COND_TO_ULTIMATE = {
    "commensal_hobic":  "commensal_hobic",
    "dysbiotic_hobic":  "dh_baseline",
    "commensal_static": "commensal_static",
    "dysbiotic_static": "dysbiotic_static",
}


def uniform_subsample(samples: np.ndarray, k: int, seed: int = 42) -> np.ndarray:
    """Draw K samples uniformly spaced across the full chain (not random)."""
    n = len(samples)
    if n <= k:
        return samples
    idx = np.linspace(0, n - 1, k, dtype=int)
    return samples[idx]


def run_0d_batch(cond: str, thetas: np.ndarray) -> list[dict]:
    """Run solve_0d_single for each row in thetas. No basin filter."""
    from posterior_ci_0d import solve_0d_single   # in _FEM/
    results = []
    n = len(thetas)
    for i, theta in enumerate(thetas):
        t0 = time.time()
        try:
            r = solve_0d_single(theta)
            phi = np.array(r["phi_final"])
            if phi.sum() < 0.01:
                continue   # degenerate (all-zero); discard
            results.append({
                "sample_idx": int(i),
                "phi_final":  phi.tolist(),
                "phi_So": float(phi[0]),
                "phi_An": float(phi[1]),
                "phi_Vd": float(phi[2]),
                "phi_Fn": float(phi[3]),
                "phi_Pg": float(phi[4]),
                "di_0d":  r.get("di_0d", None),
            })
        except Exception as e:
            print(f"    [{cond}] sample {i} failed: {e}")
        if (i + 1) % 10 == 0:
            print(f"    [{cond}] {i+1}/{n}  t={time.time()-t0:.1f}s/sample")
    return results


def main(n_samples: int = 50, seed: int = 42, conditions: list[str] | None = None):
    conds = conditions or CONDITIONS
    all_ok = True
    for cond in conds:
        up_name = COND_TO_ULTIMATE[cond]
        samples_path = _ULTIMATE / up_name / "samples.npy"
        if not samples_path.exists():
            print(f"[{cond}] ultimate_10000p not found: {samples_path}")
            all_ok = False
            continue

        all_samples = np.load(samples_path)
        thetas = uniform_subsample(all_samples, n_samples, seed)
        print(f"\n[{cond}] {len(thetas)}/{len(all_samples)} samples (uniform stride)")

        results = run_0d_batch(cond, thetas)
        print(f"  → {len(results)}/{len(thetas)} succeeded")

        # Save
        out_dir = _CI0D / cond
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "samples_0d_ultimate.json"
        json.dump(results, open(out, "w"), indent=2)
        print(f"  Saved: {out}")

        if len(results) < 5:
            print(f"  ⚠️  only {len(results)} valid samples — CI will be unreliable")
            all_ok = False

    return all_ok


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=50)
    parser.add_argument("--seed",      type=int, default=42)
    parser.add_argument("--conditions", nargs="+", default=None,
                        choices=CONDITIONS)
    args = parser.parse_args()
    ok = main(args.n_samples, args.seed, args.conditions)
    sys.exit(0 if ok else 1)
