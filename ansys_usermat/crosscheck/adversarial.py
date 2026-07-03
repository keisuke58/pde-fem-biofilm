#!/usr/bin/env python3
"""adversarial.py — try to BREAK the Abaqus-UMAT / ANSYS-USERMAT equivalence.

`crosscheck.py` proves the two constitutive cores agree on 20 curated states.
This harness is adversarial: instead of confirming a handful of nice cases it
*hunts* for a deformation / parameter state where the ANSYS USERMAT core and the
verified Abaqus UMAT core diverge — across a wide random sweep and a battery of
pathological states that deliberately drive every branch in the core:

  * near-incompressible volumetric branch      (D1 -> 0, PRESS ~ (2/D1)(J-1)J)
  * near-singular elastic Fe / detFe clamp     (|F| small, growth-driven collapse)
  * Mooney-Rivlin C01 branch                    (mtype=1, the term a real bug hid in)
  * elastic vs frozen viscous limits           (eta -> 0, eta -> huge)
  * shrinkage / extreme growth                  (alpha near -1, alpha large)
  * tiny / huge time increments                 (dt across decades)
  * heavy rotation & simple shear in each plane (Voigt 5<->6 swap integrity)

It also runs a *physical* correctness probe independent of the port: Cauchy-stress
frame indifference, sigma(Q F) == Q sigma(F) Q^T, on each core separately.

Because both cores execute the same float64 algebra, agreement is expected to be
*exact* (0 ULP) — as `crosscheck.py` already sees. The point of the sweep is that
"exact on 20 states" is a much weaker claim than "exact on thousands of states
including the nasty ones"; the latter is what the thesis verification chapter wants.

    python ansys_usermat/crosscheck/adversarial.py            # report
    python ansys_usermat/crosscheck/adversarial.py -n 20000   # bigger sweep
    python -m pytest ansys_usermat/crosscheck/adversarial.py   # as a test
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]
_ABQ_SRC = _ROOT / "umat_biofilm_visco.f"
_ANS_SRC = _ROOT / "ansys_usermat" / "usermat_biofilm.f"
_ABA_INC = _ROOT / "umat_tangent_test"          # holds ABA_PARAM.INC

# Voigt (0-based (i,j)) reconstruction maps
ABQ_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]   # 11,22,33,12,13,23
ANS_MAP = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]   # 11,22,33,12,23,13

EXACT_TOL = 1e-11        # equivalence must be at/below this (expected: exactly 0)
OBJ_TOL = 1e-7           # frame-indifference residual (finite-strain, float64)


# --------------------------------------------------------------------------- #
# build + batched execution
# --------------------------------------------------------------------------- #
def build(tmp: Path) -> tuple[Path, Path]:
    xabq, xans = tmp / "fabq", tmp / "fans"
    subprocess.run(
        ["gfortran", "-ffixed-line-length-132", f"-I{_ABA_INC}",
         str(_HERE / "fuzz_driver_abq.f"), str(_ABQ_SRC), "-o", str(xabq)],
        check=True)
    subprocess.run(
        ["gfortran", "-ffixed-line-length-132",
         str(_HERE / "fuzz_driver_ans.f"), str(_ANS_SRC), "-o", str(xans)],
        check=True)
    return xabq, xans


def _case_block(F, Fv, p) -> str:
    a, c10, c01, d1, eta, mtype, dt = p
    r = [" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)),
         " ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)),
         f"{a:.17e} {c10:.17e} {c01:.17e} {d1:.17e} {eta:.17e} {mtype:.1f} {dt:.17e}"]
    return "\n".join(r)


def run_batch(exe: Path, cases) -> np.ndarray:
    """Feed every case through one process; return (N,16) output array."""
    stdin = "\n".join(_case_block(F, Fv, p) for F, Fv, p in cases) + "\n"
    out = subprocess.run([str(exe)], input=stdin, capture_output=True,
                         text=True, check=True).stdout.split()
    vals = np.array([float(x.replace("E", "e").replace("D", "e")) for x in out])
    return vals.reshape(-1, 16)


def _tensor(sv, vmap):
    s = np.zeros((3, 3))
    for k, (i, j) in enumerate(vmap):
        s[i, j] = s[j, i] = sv[k]
    return s


def compare(rows_abq, rows_ans):
    """Per-case max discrepancy in Cauchy stress tensor, Fv_new, detFe."""
    d = np.empty(len(rows_abq))
    for k, (ra, rb) in enumerate(zip(rows_abq, rows_ans)):
        d_sig = np.abs(_tensor(ra[:6], ABQ_MAP) - _tensor(rb[:6], ANS_MAP)).max()
        d_fv = np.abs(ra[6:15] - rb[6:15]).max()
        d_det = abs(ra[15] - rb[15])
        d[k] = max(d_sig, d_fv, d_det)
    return d


# --------------------------------------------------------------------------- #
# case generators
# --------------------------------------------------------------------------- #
def rand_F(rng, scale):
    """Random finite-strain gradient with a guaranteed-positive determinant."""
    F = np.eye(3) + scale * (rng.random((3, 3)) - 0.5)
    # keep det well away from 0 so we test physics, not just the clamp
    if np.linalg.det(F) < 0.05:
        F = np.eye(3) + 0.5 * scale * (rng.random((3, 3)) - 0.5)
    return F


def fuzz_cases(rng, n):
    """Wide random sweep across deformation and every material parameter."""
    cases = []
    for _ in range(n):
        F = rand_F(rng, scale=float(rng.uniform(0.05, 1.6)))
        Fv = np.eye(3) + 0.15 * (rng.random((3, 3)) - 0.5)
        c10 = float(10 ** rng.uniform(-5, -1))
        mtype = float(rng.integers(0, 2))
        c01 = c10 * float(rng.uniform(0, 0.8)) if mtype else 0.0
        d1 = float(10 ** rng.uniform(-2, 6))            # incompressible..soft
        eta = 0.0 if rng.random() < 0.15 else float(10 ** rng.uniform(-6, 2))
        alpha = float(rng.uniform(-0.8, 3.0))           # shrink..strong growth
        dt = float(10 ** rng.uniform(-6, 3))
        cases.append((F, Fv, (alpha, c10, c01, d1, eta, mtype, dt)))
    return cases


def pathological_cases():
    """Curated corner cases that each pin one branch to an extreme."""
    I = np.eye(3)
    b = dict(c10=0.2e-3, c01=0.0, d1=5000.0, eta=8e-3, mtype=0.0, dt=5.0)

    def C(name, F, Fv=None, **kw):
        p = {**b, **kw}
        return (name, F, I if Fv is None else Fv,
                (p["alpha"], p["c10"], p["c01"], p["d1"], p["eta"],
                 p["mtype"], p["dt"]))

    S12 = I + 0.9 * np.outer([1, 0, 0], [0, 1, 0])
    S13 = I + 0.9 * np.outer([1, 0, 0], [0, 0, 1])
    S23 = I + 0.9 * np.outer([0, 1, 0], [0, 0, 1])
    th = 0.7
    Q = np.array([[np.cos(th), -np.sin(th), 0],
                  [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    return [
        C("near-incompressible", np.diag([1.2, 0.95, 0.9]), alpha=0.2, d1=1e-3),
        C("very-soft-volum.", np.diag([1.3, 0.9, 0.85]), alpha=0.4, d1=1e6),
        C("near-singular-F", np.diag([0.06, 0.06, 0.06]), alpha=0.0),
        C("shrinkage alpha->-1", np.diag([1.1, 1.0, 0.9]), alpha=-0.75),
        C("extreme growth", I, alpha=8.0),
        C("MR near-incompr.", np.diag([1.25, 0.9, 0.9]), alpha=0.2,
          c01=0.15e-3, mtype=1.0, d1=1e-3),
        C("MR + shear12", S12, alpha=0.15, c01=0.1e-3, mtype=1.0),
        C("shear12 large", S12, alpha=0.1),
        C("shear13 large", S13, alpha=0.1),
        C("shear23 large", S23, alpha=0.1),
        C("pure rotation", Q, alpha=0.0),
        C("rotation+growth", Q @ np.diag([1.1, 1.0, 0.9]), alpha=0.3),
        C("tiny dt", np.diag([1.1, 1.0, 0.95]), alpha=0.2, dt=1e-6, eta=1e-3),
        C("huge dt", np.diag([1.1, 1.0, 0.95]), alpha=0.2, dt=1e3, eta=1e-3),
        C("eta->0 elastic", np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=0.0),
        C("eta huge frozen", np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=1e8),
        C("prior Fv sheared", np.diag([1.1, 1.0, 0.9]),
          Fv=I + 0.2 * np.outer([1, 0, 0], [0, 1, 0]), alpha=0.2),
    ]


# --------------------------------------------------------------------------- #
# physical probe: Cauchy-stress frame indifference (per core)
# --------------------------------------------------------------------------- #
def objectivity(exe, rng, n=200):
    """Frame indifference of the Cauchy stress: sigma(Q F) == Q sigma(F) Q^T.

    Tested in the elastic limit (eta=0, Fv=I) where the model reduces to an
    isotropic hyperelastic response in Be = Fe Fe^T, Fe = F Fg^{-1} (Fg^{-1}
    isotropic). There Be(QF) = Q Be(F) Q^T *exactly*, so an isotropic Cauchy
    stress must co-rotate to float64 precision. (With eta>0 the two-step
    trial/update carries the intermediate-config internal variable Fv, whose
    superposed-rotation transformation is subtler and not the port's concern.)
    Still adversarial: any sign flip or index transposition in the stress
    assembly — including the Voigt 5<->6 map — breaks co-rotation.
    """
    vmap = ABQ_MAP if exe.name == "fabq" else ANS_MAP
    I = np.eye(3)
    base, rot, Qs = [], [], []
    for _ in range(n):
        # moderate, well-conditioned F so detFe stays away from the clamp
        F = I + 0.6 * (rng.random((3, 3)) - 0.5)
        if np.linalg.det(F) < 0.3:
            F = I + 0.3 * (rng.random((3, 3)) - 0.5)
        p = (float(rng.uniform(0.0, 1.0)), 0.2e-3,
             (0.1e-3 if rng.random() < 0.5 else 0.0), 5000.0,
             0.0,                                  # eta = 0 (elastic)
             (1.0 if rng.random() < 0.5 else 0.0), 5.0)
        Q, _ = np.linalg.qr(rng.random((3, 3)))    # random rotation via QR
        if np.linalg.det(Q) < 0:
            Q[:, 0] = -Q[:, 0]
        Qs.append(Q)
        base.append((F, I, p))
        rot.append((Q @ F, I, p))                  # Fv=I fixed; rotate F only
    rb = run_batch(exe, base)
    rr = run_batch(exe, rot)
    worst = 0.0
    for k in range(n):
        s0 = _tensor(rb[k][:6], vmap)
        s1 = _tensor(rr[k][:6], vmap)
        worst = max(worst, np.abs(Qs[k] @ s0 @ Qs[k].T - s1).max())
    return worst


# --------------------------------------------------------------------------- #
def run_all(n=8000, seed=12345, verbose=True):
    rng = np.random.default_rng(seed)
    with tempfile.TemporaryDirectory() as td:
        xabq, xans = build(Path(td))

        # 1) pathological battery (named) --------------------------------
        path = pathological_cases()
        pa = run_batch(xabq, [(F, Fv, p) for _, F, Fv, p in path])
        pn = run_batch(xans, [(F, Fv, p) for _, F, Fv, p in path])
        pd = compare(pa, pn)

        # 2) wide random fuzz --------------------------------------------
        fz = fuzz_cases(rng, n)
        fa = run_batch(xabq, fz)
        fn = run_batch(xans, fz)
        fd = compare(fa, fn)

        # any non-finite output is itself a failure
        nonfinite = int((~np.isfinite(fa)).sum() + (~np.isfinite(fn)).sum())

        # 3) frame indifference (physical, per core) ---------------------
        obj_abq = objectivity(xabq, np.random.default_rng(seed + 1))
        obj_ans = objectivity(xans, np.random.default_rng(seed + 1))

        worst_equiv = max(pd.max(), fd.max())
        ipath = int(np.argmax(pd))
        ifz = int(np.argmax(fd))

        if verbose:
            print("Adversarial Abaqus-UMAT vs ANSYS-USERMAT equivalence hunt")
            print("=" * 62)
            print(f"[1] pathological battery ({len(path)} cases)")
            for (name, *_), d in zip(path, pd):
                mark = "ok" if d < EXACT_TOL else "!!"
                print(f"      {mark} {name:<22} |Δ| = {d:.2e}")
            print(f"    worst: {path[ipath][0]}  ({pd[ipath]:.2e})")
            print(f"[2] random fuzz ({n} cases, wide param+strain sweep)")
            print(f"      worst |Δ| = {fd[ifz]:.2e}   median = "
                  f"{np.median(fd):.2e}   non-finite outputs = {nonfinite}")
            if fd[ifz] >= EXACT_TOL:
                print("      offending case params "
                      "(alpha,c10,c01,d1,eta,mtype,dt):")
                print(f"        {tuple(round(v, 6) for v in fz[ifz][2])}")
            print(f"[3] Cauchy-stress frame indifference (per core, physical)")
            print(f"      Abaqus core  max |σ(QF) − Qσ(F)Qᵀ| = {obj_abq:.2e}")
            print(f"      ANSYS  core  max |σ(QF) − Qσ(F)Qᵀ| = {obj_ans:.2e}")
            print("=" * 62)
            ok = (worst_equiv < EXACT_TOL and nonfinite == 0
                  and max(obj_abq, obj_ans) < OBJ_TOL)
            print(f"  equivalence worst over {len(path)+n} cases: "
                  f"{worst_equiv:.2e}")
            print("  => " + ("PASS — ANSYS USERMAT core is bit-identical to the "
                             "verified Abaqus law across the full adversarial "
                             "sweep, and both are frame-indifferent."
                             if ok else "FAIL — see cases marked above."))
        return dict(worst_equiv=worst_equiv, nonfinite=nonfinite,
                    obj_abq=obj_abq, obj_ans=obj_ans)


def test_adversarial_equivalence():
    r = run_all(n=4000, verbose=False)
    assert r["worst_equiv"] < EXACT_TOL
    assert r["nonfinite"] == 0
    assert max(r["obj_abq"], r["obj_ans"]) < OBJ_TOL


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=8000, help="random fuzz case count")
    ap.add_argument("--seed", type=int, default=12345)
    a = ap.parse_args()
    r = run_all(n=a.n, seed=a.seed)
    ok = (r["worst_equiv"] < EXACT_TOL and r["nonfinite"] == 0
          and max(r["obj_abq"], r["obj_ans"]) < OBJ_TOL)
    sys.exit(0 if ok else 1)
