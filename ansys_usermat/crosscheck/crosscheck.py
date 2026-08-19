#!/usr/bin/env python3
"""crosscheck.py — constitutive-core equivalence between two UMAT/USERMAT sources.

Compiles two constitutive cores (the *actual Fortran*, not a re-implementation)
and drives both over a battery of deformation states, verifying that they return
the **same Cauchy stress tensor, viscous state, and Je to machine precision**.

Defaults to the pair this repo has verified so far:

  Abaqus  umat_biofilm_visco.f   BIOFILM_STRESS_CORE   Voigt order 11,22,33,12,13,23
  ANSYS   usermat_biofilm.f      BIOFILM_STRESS_CORE   Voigt order 11,22,33,12,23,13

Pass --right-src/--right-driver/--right-voigt (etc.) to compare against a third
implementation instead — e.g. Felix's UserElement/UserMat core, once it arrives
(see THESIS_ASSIGNMENT.md §4.2). Writing the driver for a new source is the
part that can't be automated (it has to call whatever that source's entry
point actually is) — copy xcheck_driver_template.f and fill it in.

The left/right files may use different Voigt component orderings on purpose;
this harness reconstructs the full 3x3 tensor from each and compares. Because
the tangent in both top-level routines is a finite-difference of this same
core, core equivalence implies tangent equivalence.

Run:
    python ansys_usermat/crosscheck/crosscheck.py          # Abaqus vs ANSYS (default)
    python ansys_usermat/crosscheck/crosscheck.py \\
        --right-src path/to/felix_core.f \\
        --right-driver path/to/xcheck_driver_felix.f \\
        --right-voigt 11,22,33,12,13,23 --right-name felix   # ANSYS vs Felix
    python -m pytest ansys_usermat/crosscheck/crosscheck.py # as a test (default pair)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]                       # repo root
_ABQ_SRC = _ROOT / "umat_biofilm_visco.f"
_ANS_SRC = _ROOT / "ansys_usermat" / "usermat_biofilm.f"
_ABA_INC = _ROOT / "umat_tangent_test"         # holds ABA_PARAM.INC

# Voigt (0-based (i,j)) reconstruction maps, keyed by the same order strings
# used on the CLI ("11,22,33,12,13,23" style).
ABQ_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]   # 11,22,33,12,13,23
ANS_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]   # 11,22,33,12,23,13
_VOIGT_PAIR_TO_MAP = {"12,13,23": ABQ_MAP, "12,23,13": ANS_MAP}


def _parse_voigt(spec: str):
    """'11,22,33,12,13,23' -> the (i,j) map, reusing ABQ_MAP/ANS_MAP for the two
    orderings this repo has needed so far; extend _VOIGT_PAIR_TO_MAP if a third
    implementation uses yet another shear ordering."""
    tail = ",".join(spec.split(",")[3:])  # just the shear part distinguishes them
    try:
        return _VOIGT_PAIR_TO_MAP[tail]
    except KeyError:
        raise SystemExit(
            f"Unrecognised Voigt order {spec!r} (shear part {tail!r}). "
            f"Known shear orderings: {sorted(_VOIGT_PAIR_TO_MAP)}. Add a new "
            f"entry to _VOIGT_PAIR_TO_MAP if this is a genuinely new convention."
        )


def _build_pair(tmp: Path, left_driver: Path, left_src: Path,
                 right_driver: Path, right_src: Path,
                 left_inc: Path | None) -> tuple[Path, Path]:
    xleft, xright = tmp / "xleft", tmp / "xright"
    left_cmd = ["gfortran", "-ffixed-line-length-132"]
    if left_inc is not None:
        left_cmd += [f"-I{left_inc}"]
    left_cmd += [str(left_driver), str(left_src), "-o", str(xleft)]
    subprocess.run(left_cmd, check=True)
    subprocess.run(
        ["gfortran", "-ffixed-line-length-132",
         str(right_driver), str(right_src), "-o", str(xright)],
        check=True)
    return xleft, xright


def _build(tmp: Path) -> tuple[Path, Path]:
    """Back-compat wrapper: the original Abaqus-vs-ANSYS pair."""
    return _build_pair(
        tmp,
        _HERE / "xcheck_driver_abq.f", _ABQ_SRC,
        _HERE / "xcheck_driver_ans.f", _ANS_SRC,
        _ABA_INC,
    )


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


def run_all(verbose=True, *, left_driver=None, left_src=None, left_inc=None,
            left_map=ABQ_MAP, left_name="Abaqus UMAT",
            right_driver=None, right_src=None, right_map=ANS_MAP,
            right_name="ANSYS USERMAT"):
    left_driver = left_driver or (_HERE / "xcheck_driver_abq.f")
    left_src = left_src or _ABQ_SRC
    left_inc = _ABA_INC if left_inc is None and left_driver == _HERE / "xcheck_driver_abq.f" else left_inc
    right_driver = right_driver or (_HERE / "xcheck_driver_ans.f")
    right_src = right_src or _ANS_SRC

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        xleft, xright = _build_pair(Path(td), left_driver, left_src,
                                     right_driver, right_src, left_inc)
        worst = 0.0
        rows = []
        for name, F, Fv, params in _cases():
            sa, fva, da = _run(xleft, F, Fv, params)
            sb, fvb, db = _run(xright, F, Fv, params)
            d_sig = np.abs(_tensor(sa, left_map) - _tensor(sb, right_map)).max()
            d_fv = np.abs(fva - fvb).max()
            d_det = abs(da - db)
            d = max(d_sig, d_fv, d_det)
            worst = max(worst, d)
            rows.append((name, d_sig, d_fv, d_det))
        if verbose:
            print(f"{left_name}  vs  {right_name}  (Cauchy stress tensor)")
            print(f"  {'case':<18} {'|Δσ|':>10} {'|ΔFv|':>10} {'|ΔJe|':>10}")
            print(f"  {'-'*50}")
            for name, ds, df, dd in rows:
                print(f"  {name:<18} {ds:>10.2e} {df:>10.2e} {dd:>10.2e}")
            print(f"  {'-'*50}")
            print(f"  worst discrepancy over {len(rows)} cases: {worst:.2e}")
            print(f"  => {right_name} reproduces {left_name} "
                  + ("to machine precision." if worst < 1e-11 else "WITH DISCREPANCY."))
        return worst


def test_ansys_matches_abaqus():
    assert run_all(verbose=False) < 1e-11


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--left-src", type=Path, default=_ABQ_SRC)
    p.add_argument("--left-driver", type=Path, default=_HERE / "xcheck_driver_abq.f")
    p.add_argument("--left-inc", type=Path, default=_ABA_INC)
    p.add_argument("--left-voigt", default="11,22,33,12,13,23",
                    help="Voigt shear order for --left-src, e.g. 11,22,33,12,13,23")
    p.add_argument("--left-name", default="Abaqus UMAT")
    p.add_argument("--right-src", type=Path, default=_ANS_SRC)
    p.add_argument("--right-driver", type=Path, default=_HERE / "xcheck_driver_ans.f")
    p.add_argument("--right-voigt", default="11,22,33,12,23,13",
                    help="Voigt shear order for --right-src, e.g. 11,22,33,12,23,13")
    p.add_argument("--right-name", default="ANSYS USERMAT")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    worst = run_all(
        left_driver=args.left_driver, left_src=args.left_src, left_inc=args.left_inc,
        left_map=_parse_voigt(args.left_voigt), left_name=args.left_name,
        right_driver=args.right_driver, right_src=args.right_src,
        right_map=_parse_voigt(args.right_voigt), right_name=args.right_name,
    )
    sys.exit(0 if worst < 1e-11 else 1)
