#!/usr/bin/env python3
"""How much does t_growth_cylinder_ecology.dat's result depend on THETA_DEMO
specifically, versus being a structural property of the two seed
compositions (region A vs region B in ecology_cylinder_reference.py)?

THETA_DEMO is jax_hamilton_0d_5species_demo.py's placeholder, not a
TMCMC-calibrated value (coupling/README.md next-steps #4). There is no
"theta-free" mode -- the 0D Hamilton reaction term structurally needs some
20-vector (15 independent A_ij + 5 b_i, theta_to_matrices). The choices are:
demo placeholder, TMCMC-calibrated (deferred), or a hand-picked value. This
script checks whether the *qualitative* result reference.py reports (region
A's all-zero/default seed grows faster than region B's higher-initial-
biomass seed, ratio B/A ~ 0.63) is a robust feature of the two seeds, or an
artifact of THETA_DEMO's specific numbers -- by rerunning the same 10-step
chained ecology_step trajectory under several theta variants:

  demo      THETA_DEMO itself (baseline, matches ecology_cylinder_reference.py)
  weak      A scaled x0.3 (much weaker species interactions)
  strong    A scaled x3.0 (much stronger species interactions)
  sign_flip off-diagonal A entries negated (cooperative <-> competitive),
            diagonal and b left alone
  random_N  5 random theta draws, same magnitude/sign range as THETA_DEMO
            component-wise (fixed seed, reproducible), b left at THETA_DEMO

If the B/A ordering (which region grows faster) flips across variants, that
result is theta-dependent and should not be described as a physical
prediction before calibration -- only as evidence the plumbing responds to
theta at all. If it's stable across every variant tried, that is at least
weak evidence the ordering comes from the seed compositions themselves
(region B saturates psi faster, etc.) rather than from THETA_DEMO's
particular numbers -- still not a calibrated result, but a different kind
of claim than "this one demo theta happens to give this order."
"""
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coupling"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import ecology_jax
from jax_hamilton_0d_5species_demo import THETA_DEMO, theta_to_matrices

DT = 1.0e-5
N_SUBSTEPS = 10
K_ALPHA = 50.0

REGION_B_SEED = np.array([
    0.05, 0.05, 0.15, 0.25, 0.20,
    0.30,
    0.999, 0.999, 0.999, 0.999, 0.999,
    0.0,
])

# theta layout (theta_to_matrices): indices 0,1,2 -> A[0,0],A[0,1]=A[1,0],A[1,1];
# 3,4 -> b[0],b[1]; 5,6,7 -> A[2,2],A[2,3]=A[3,2],A[3,3]; 8,9 -> b[2],b[3];
# 10,11,12,13 -> A[0,2],A[0,3],A[1,2],A[1,3] (off-diagonal cross terms);
# 14 -> A[4,4]; 15 -> b[4]; 16,17,18,19 -> A[0,4],A[1,4],A[2,4],A[3,4].
DIAG_IDX = [0, 2, 5, 7, 14]                  # A[i,i] entries
OFFDIAG_IDX = [1, 6, 10, 11, 12, 13, 16, 17, 18, 19]  # A[i,j], i!=j entries
B_IDX = [3, 4, 8, 9, 15]                     # b entries


def make_variant(name, base=THETA_DEMO):
    t = np.asarray(base, dtype=np.float64).copy()
    if name == "demo":
        return t
    if name == "weak":
        t[DIAG_IDX] *= 0.3
        t[OFFDIAG_IDX] *= 0.3
        return t
    if name == "strong":
        t[DIAG_IDX] *= 3.0
        t[OFFDIAG_IDX] *= 3.0
        return t
    if name == "sign_flip":
        t[OFFDIAG_IDX] *= -1.0
        return t
    raise ValueError(name)


def random_variant(rng, base=THETA_DEMO):
    t = np.asarray(base, dtype=np.float64).copy()
    lo, hi = float(np.min(t[DIAG_IDX + OFFDIAG_IDX])), float(np.max(t[DIAG_IDX + OFFDIAG_IDX]))
    idx = DIAG_IDX + OFFDIAG_IDX
    t[idx] = rng.uniform(lo, hi, size=len(idx))
    return t


def chained_alpha(g0, theta):
    g = np.asarray(g0, dtype=np.float64)
    alpha = 0.0
    for _ in range(N_SUBSTEPS):
        g = np.asarray(ecology_jax.ecology_step(g, theta, DT))
        alpha += DT * K_ALPHA * ecology_jax.living_fraction_total(g)
    return alpha


def report(label, theta):
    a_a = chained_alpha(ecology_jax.default_initial_state(), theta)
    a_b = chained_alpha(REGION_B_SEED, theta)
    ratio = a_b / a_a
    faster = "A" if ratio < 1.0 else "B"
    print(f"{label:12s}  alpha_A={a_a:.6e}  alpha_B={a_b:.6e}  "
          f"B/A={ratio:.4f}  faster-growing region: {faster}")
    return ratio


if __name__ == "__main__":
    print(f"dt={DT}, N_SUBSTEPS={N_SUBSTEPS}, k_alpha={K_ALPHA}\n")
    ratios = {}
    for name in ("demo", "weak", "strong", "sign_flip"):
        ratios[name] = report(name, make_variant(name))

    rng = np.random.default_rng(0)
    for i in range(5):
        ratios[f"random_{i}"] = report(f"random_{i}", random_variant(rng))

    faster_regions = {("A" if r < 1.0 else "B") for r in ratios.values()}
    print()
    if len(faster_regions) == 1:
        print(f"STABLE across all {len(ratios)} variants tried: region "
              f"{faster_regions.pop()} always grows faster.")
    else:
        print(f"NOT STABLE: the faster-growing region flips across variants "
              f"({sorted(faster_regions)}) -- the B/A ordering is theta-"
              f"dependent, not a fixed property of the two seeds.")
