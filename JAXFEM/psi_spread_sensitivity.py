#!/usr/bin/env python3
"""How different would per-species viability have to be for composition to
reach the stress at all?

The four clinical conditions come out within 0.4% of each other in stress, and
the reason is structural: the condition reaches growth only through
alpha = k_alpha * integral(phi_tot), and phi_tot = sum_i phi_i psi_i is a
single total. The CLSM compositions are fractions, so every condition starts
at sum_i phi_i = 1 and the differences cancel in that sum.

That cancellation is exact only because psi is species-independent. The model
seeds it that way -- ecology_jax.default_initial_state and every 4-region
reference use psi_i = 0.999 for all five -- and psi is precisely the quantity
this repository has declined to substitute, because the workbook's "living
cells" ratio exceeds 1 and is not the model's psi in [0, 1].

So "composition does not reach the stress" is a statement about the model as
seeded, not a property of the physics. This script measures what it would take
to overturn it: give the species different viabilities, sweep how different,
and watch the spread in alpha across the four conditions come back.

The pattern of a real per-species psi is unknown, so nothing here assumes one.
At each spread s the script samples many random patterns and reports the
distribution, with s defined so that it reads directly:

    psi_i = PSI_MAX * (1 - s * v_i),   v_i in [0, 1], min 0, max 1

so s is exactly (psi_max - psi_min) / psi_max, and psi never exceeds PSI_MAX.
One pattern is drawn per sample and applied to all four conditions, since
viability is a property of a species, not of the condition it sits in.

The output is the answer to "we have no data": not a number for psi, but the
precision any future measurement of it would have to reach before it could
change the conclusion.

    python JAXFEM/psi_spread_sensitivity.py
    python JAXFEM/psi_spread_sensitivity.py --samples 100 --horizon deck
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "ansys_usermat" / "coupling"))
sys.path.insert(0, str(_REPO / "ansys_usermat" / "apdl"))

import ecology_jax                                           # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO         # noqa: E402
from plot_heine_phi_psi import load                          # noqa: E402

XLSX = _REPO / "data" / "heine_species_distribution_biofilm.xlsx"
PSI_MAX = 0.999
K_ALPHA = 50.0

CONDITIONS = [
    ("CS", "Commensal", "Static all cells"),
    ("CH", "Commensal", "HOBIC all cells"),
    ("DS", "Dysbiotic", "Static all cells"),
    ("DH", "Dysbiotic", "HOBIC all cells"),
]

# The deck's own increment, and the time the alpha spread was found to peak at.
HORIZONS = {"deck": (1.0e-4, 10), "peak": (1.0e-2, 1000)}

SPREADS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50]


def day1_compositions():
    """Day-1 (Tag 1) CLSM composition per condition, as the references seed it."""
    import openpyxl
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    out = {}
    for name, sheet, key in CONDITIONS:
        rows = load(wb[sheet])[key][1]
        phi = np.array([np.nanmean(np.array(r, dtype=float)) for r in rows]) / 100.0
        out[name] = phi * (0.999999 / phi.sum())
    return out


def alpha(phi, psi, dt_h, n_sub):
    g0 = np.concatenate([phi, [1.0 - phi.sum()], psi, [0.0]])
    _, phi_int = ecology_jax.ecology_substeps(g0, THETA_DEMO, dt_h, n_sub)
    return K_ALPHA * phi_int


def pattern(rng):
    """A random per-species viability shape, rescaled to span [0, 1] exactly."""
    v = rng.random(5)
    lo, hi = v.min(), v.max()
    return (v - lo) / (hi - lo) if hi > lo else np.zeros(5)


def _crossing(ss, med, target):
    """The smallest spread at which the median reaches `target`.

    Not np.interp: that assumes the x it searches is increasing, and the
    median is not. Over the longer horizon a small psi spread slightly
    *reduces* the condition spread before it grows, so interpolating against
    the whole curve is reading a non-monotonic function backwards. Scan for
    the first crossing instead and interpolate only inside the segment that
    brackets it, and say plainly when the baseline is already above the
    target rather than reporting s = 0 for it.
    """
    if med[0] >= target:
        return f"already exceeded with psi uniform ({med[0]:.2%} at s = 0)"
    for i in range(1, len(med)):
        if med[i] >= target:
            lo, hi = med[i - 1], med[i]
            frac = 0.0 if hi == lo else (target - lo) / (hi - lo)
            return f"s ~ {ss[i - 1] + frac * (ss[i] - ss[i - 1]):.1%}"
    return f"not reached by s = {ss[-1]:.0%}"


def spread_of(alphas):
    a = np.asarray(alphas, dtype=float)
    return float((a.max() - a.min()) / a.mean())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=40)
    ap.add_argument("--horizon", choices=list(HORIZONS) + ["both"], default="both")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("JAXFEM/klempt2024_results/psi_spread.json"))
    a = ap.parse_args(argv)

    comps = day1_compositions()
    rng = np.random.default_rng(a.seed)
    horizons = list(HORIZONS) if a.horizon == "both" else [a.horizon]
    results = {}

    for hname in horizons:
        dt_h, n_sub = HORIZONS[hname]
        print(f"\nhorizon {hname}: dt_h={dt_h}, n_sub={n_sub}, "
              f"k_alpha={K_ALPHA}, psi_max={PSI_MAX}")
        print(f"{'psi spread s':>13}  {'median':>9}  {'p90':>9}  {'max':>9}"
              "   alpha spread across CS/CH/DS/DH")
        rows = []
        for s in SPREADS:
            n = 1 if s == 0.0 else a.samples      # s=0 has only one pattern
            sp = []
            for _ in range(n):
                psi = PSI_MAX * (1.0 - s * pattern(rng))
                sp.append(spread_of([alpha(comps[c], psi, dt_h, n_sub)
                                     for c, _, _ in CONDITIONS]))
            sp = np.asarray(sp)
            rows.append({"s": s, "n": n, "median": float(np.median(sp)),
                         "p90": float(np.percentile(sp, 90)),
                         "max": float(sp.max())})
            print(f"{s:13.3f}  {np.median(sp):8.3%}  {np.percentile(sp,90):8.3%}"
                  f"  {sp.max():8.3%}", flush=True)
        results[hname] = rows

        med = np.array([r["median"] for r in rows])
        ss = np.array([r["s"] for r in rows])
        print("  psi spread needed for an alpha spread of:")
        for target in (0.01, 0.05, 0.091):
            print(f"    {target:.1%}: {_crossing(ss, med, target)}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(
        {"psi_max": PSI_MAX, "k_alpha": K_ALPHA, "samples": a.samples,
         "horizons": {k: dict(dt_h=HORIZONS[k][0], n_sub=HORIZONS[k][1])
                      for k in horizons},
         "results": results}, indent=1))
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
