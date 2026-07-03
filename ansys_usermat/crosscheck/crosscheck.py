#!/usr/bin/env python3
"""crosscheck.py — Abaqus UMAT core  vs  ANSYS USERMAT core equivalence.

Compiles the two constitutive cores (the *actual Fortran*, not a re-implementation)
and drives both over a battery of deformation states, verifying that they return
the **same Cauchy stress tensor, viscous state, and Je to machine precision** —
the key claim for handing the ANSYS USERMAT to Felix.

  Abaqus  umat_biofilm_visco.f   BIOFILM_STRESS_CORE   Voigt order 11,22,33,12,13,23
  ANSYS   usermat_biofilm.f      BIOFILM_STRESS_CORE   Voigt order 11,22,33,12,23,13

The two files use different component orderings on purpose (Abaqus vs ANSYS
convention); this harness reconstructs the full 3x3 tensor from each and compares.
Because the tangent in both top-level routines is a finite-difference of this same
core, core equivalence implies tangent equivalence.

Run:
    python ansys_usermat/crosscheck/crosscheck.py          # report
    python -m pytest ansys_usermat/crosscheck/crosscheck.py # as a test
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]                       # repo root
_ABQ_SRC = _ROOT / "umat_biofilm_visco.f"
_ANS_SRC = _ROOT / "ansys_usermat" / "usermat_biofilm.f"
_ABA_INC = _ROOT / "umat_tangent_test"         # holds ABA_PARAM.INC

# Voigt (0-based (i,j)) reconstruction maps
ABQ_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]   # 11,22,33,12,13,23
ANS_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]   # 11,22,33,12,23,13


def _build(tmp: Path) -> tuple[Path, Path]:
    xabq, xans = tmp / "xabq", tmp / "xans"
    subprocess.run(
        ["gfortran", "-ffixed-line-length-132", f"-I{_ABA_INC}",
         str(_HERE / "xcheck_driver_abq.f"), str(_ABQ_SRC), "-o", str(xabq)],
        check=True)
    subprocess.run(
        ["gfortran", "-ffixed-line-length-132",
         str(_HERE / "xcheck_driver_ans.f"), str(_ANS_SRC), "-o", str(xans)],
        check=True)
    return xabq, xans


def _fmt(F, Fv, params) -> str:
    a, c10, c01, d1, eta, mtype, dt = params
    rows = []
    rows.append(" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)))
    rows.append(" ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)))
    rows.append(f"{a:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} {mtype:.1f} {dt:.17e}")
    return "\n".join(rows) + "\n"


def _run(exe: Path, F, Fv, params):
    out = subprocess.run([str(exe)], input=_fmt(F, Fv, params),
                         capture_output=True, text=True, check=True).stdout.split()
    vals = np.array([float(x.replace("E", "e").replace("D", "e")) for x in out])
    sv, fvn, det = vals[:6], vals[6:15].reshape(3, 3), vals[15]
    return sv, fvn, det


def _tensor(sv, vmap):
    s = np.zeros((3, 3))
    for k, (i, j) in enumerate(vmap):
        s[i, j] = s[j, i] = sv[k]
    return s


def _cases():
    """Deterministic battery covering identity, uniaxial, shear, growth, elastic/
    viscous limits, and Mooney-Rivlin."""
    rng = np.random.default_rng(0)
    I = np.eye(3)
    base = dict(c10=0.2e-3, c01=0.0, d1=1.0 / (0.2e-3), eta=8e-3, mtype=0.0, dt=5.0)

    def P(F, Fv=I, **kw):
        p = {**base, **kw}
        return (F, Fv, (p["alpha"], p["c10"], p["c01"], p["d1"], p["eta"],
                        p["mtype"], p["dt"]))

    cases = []
    cases.append(("identity+growth", *P(I, alpha=0.12)))
    cases.append(("uniaxial", *P(np.diag([1.2, 0.95, 0.95]), alpha=0.3)))
    cases.append(("simple shear", *P(I + 0.15 * np.outer([1, 0, 0], [0, 1, 0]), alpha=0.05)))
    cases.append(("elastic (eta=0)", *P(np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=0.0)))
    cases.append(("frozen (eta huge)", *P(np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=1e6)))
    cases.append(("large growth", *P(I, alpha=2.0)))
    cases.append(("no growth", *P(np.diag([1.05, 0.98, 1.02]), alpha=0.0)))
    cases.append(("Mooney-Rivlin", *P(np.diag([1.15, 0.95, 0.92]), alpha=0.15,
                                       c01=0.05e-3, mtype=1.0)))
    # random finite-strain + random prior viscous state
    for n in range(12):
        F = I + 0.2 * (rng.random((3, 3)) - 0.5)
        Fv = I + 0.05 * (rng.random((3, 3)) - 0.5)
        cases.append((f"random#{n}", *P(F, Fv=Fv, alpha=float(rng.uniform(0, 1.5)),
                                        eta=float(10 ** rng.uniform(-4, -1)))))
    return cases


def run_all(verbose=True):
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        xabq, xans = _build(Path(td))
        worst = 0.0
        rows = []
        for name, F, Fv, params in _cases():
            sa, fva, da = _run(xabq, F, Fv, params)
            sb, fvb, db = _run(xans, F, Fv, params)
            d_sig = np.abs(_tensor(sa, ABQ_MAP) - _tensor(sb, ANS_MAP)).max()
            d_fv = np.abs(fva - fvb).max()
            d_det = abs(da - db)
            d = max(d_sig, d_fv, d_det)
            worst = max(worst, d)
            rows.append((name, d_sig, d_fv, d_det))
        if verbose:
            print("Abaqus UMAT core  vs  ANSYS USERMAT core  (Cauchy stress tensor)")
            print(f"  {'case':<18} {'|Δσ|':>10} {'|ΔFv|':>10} {'|ΔJe|':>10}")
            print(f"  {'-'*50}")
            for name, ds, df, dd in rows:
                print(f"  {name:<18} {ds:>10.2e} {df:>10.2e} {dd:>10.2e}")
            print(f"  {'-'*50}")
            print(f"  worst discrepancy over {len(rows)} cases: {worst:.2e}")
            print("  => the ANSYS USERMAT reproduces the verified Abaqus law "
                  + ("to machine precision." if worst < 1e-11 else "WITH DISCREPANCY."))
        return worst


def test_ansys_matches_abaqus():
    assert run_all(verbose=False) < 1e-11


if __name__ == "__main__":
    sys.exit(0 if run_all() < 1e-11 else 1)
