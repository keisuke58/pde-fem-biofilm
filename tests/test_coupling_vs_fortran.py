"""Cross-validate the coupling's Python material core against the REAL Fortran core.

`ansys_usermat/coupling/material_server.py` claims to be a line-by-line mirror of
the verified Abaqus core `BIOFILM_STRESS_CORE` (umat_biofilm_visco.f). This test
proves it: it compiles the actual Fortran core (via the existing batch driver used
by the adversarial harness), drives both implementations over the same battery of
deformation states, and compares Cauchy stress, the updated viscous state Fv, and
detFe.

Closing this loop means the verification chain is complete:

    Abaqus UMAT  ==  ANSYS USERMAT  ==  Python material core

so swapping the Fortran law for the Python/JAX model at the Gauss point provably
does not change the physics.

Requires gfortran; skipped automatically if it is absent.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_XCHK = _ROOT / "ansys_usermat" / "crosscheck"
_ABQ_SRC = _ROOT / "umat_biofilm_visco.f"
_ABA_INC = _ROOT / "umat_tangent_test"          # holds ABA_PARAM.INC
_DRIVER = _XCHK / "fuzz_driver_abq.f"

sys.path.insert(0, str(_ROOT / "ansys_usermat" / "coupling"))
import material_server as ms                     # noqa: E402

pytestmark = pytest.mark.skipif(
    shutil.which("gfortran") is None or not _ABQ_SRC.exists() or not _DRIVER.exists(),
    reason="gfortran or the Fortran core/driver is unavailable")

# Python core returns Abaqus Voigt order 11,22,33,12,13,23 — same as the driver.
#
# Agreement is asserted RELATIVE to the magnitude of the Fortran result. That is
# the meaningful comparison here: the random battery deliberately includes
# degenerate states (a large viscous step drives Fv toward singular, the detFe
# clamp fires, and the returned stress reaches ~1e30 — non-physical but a valid
# code path). An absolute tolerance is meaningless against 1e30; both
# implementations agree to ~1e-15 relative across the whole battery, including
# those clamped states.
RTOL = 1e-11
TINY = 1e-30


def _cases():
    """(F, Fv, params) battery: named corner cases + random finite strains."""
    rng = np.random.default_rng(20260804)
    I = np.eye(3)
    base = dict(c10=2.0e-4, c01=0.0, d1=5000.0, eta=8e-3, mtype=0.0, dt=5.0)

    def P(F, Fv=I, **kw):
        p = {**base, **kw}
        return (F, Fv, (p["alpha"], p["c10"], p["c01"], p["d1"], p["eta"],
                        p["mtype"], p["dt"]))

    out = [
        ("identity+growth", *P(I, alpha=0.12)),
        ("uniaxial", *P(np.diag([1.2, 0.95, 0.95]), alpha=0.3)),
        ("shear12", *P(I + 0.15 * np.outer([1, 0, 0], [0, 1, 0]), alpha=0.05)),
        ("shear13", *P(I + 0.15 * np.outer([1, 0, 0], [0, 0, 1]), alpha=0.05)),
        ("shear23", *P(I + 0.15 * np.outer([0, 1, 0], [0, 0, 1]), alpha=0.05)),
        ("elastic eta=0", *P(np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=0.0)),
        ("frozen eta huge", *P(np.diag([1.1, 1.0, 0.9]), alpha=0.2, eta=1e6)),
        ("large growth", *P(I, alpha=2.0)),
        ("no growth", *P(np.diag([1.05, 0.98, 1.02]), alpha=0.0)),
        ("Mooney-Rivlin", *P(np.diag([1.15, 0.95, 0.92]), alpha=0.15,
                             c01=0.5e-4, mtype=1.0)),
        ("MR + shear", *P(I + 0.2 * np.outer([1, 0, 0], [0, 1, 0]), alpha=0.1,
                          c01=0.4e-4, mtype=1.0)),
        ("near-incompressible", *P(np.diag([1.2, 0.95, 0.9]), alpha=0.2, d1=1e-3)),
        ("prior Fv sheared", *P(np.diag([1.1, 1.0, 0.9]),
                                Fv=I + 0.2 * np.outer([1, 0, 0], [0, 1, 0]), alpha=0.2)),
    ]
    for n in range(15):
        F = I + 0.25 * (rng.random((3, 3)) - 0.5)
        Fv = I + 0.06 * (rng.random((3, 3)) - 0.5)
        mtype = float(rng.integers(0, 2))
        out.append((f"random#{n}", *P(F, Fv=Fv,
                                      alpha=float(rng.uniform(0, 1.5)),
                                      c01=(1.0e-4 * rng.random() if mtype else 0.0),
                                      mtype=mtype,
                                      eta=float(10 ** rng.uniform(-4, -1)))))
    return out


def _fortran_results(cases, exe):
    """Run every case through the compiled Fortran core; return (N,16) array."""
    blocks = []
    for _, F, Fv, p in cases:
        a, c10, c01, d1, eta, mtype, dt = p
        blocks.append(" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)))
        blocks.append(" ".join(f"{Fv[i, j]:.17e}" for i in range(3) for j in range(3)))
        blocks.append(f"{a:.17e} {c10:.17e} {c01:.17e} {d1:.17e} "
                      f"{eta:.17e} {mtype:.1f} {dt:.17e}")
    out = subprocess.run([str(exe)], input="\n".join(blocks) + "\n",
                         capture_output=True, text=True, check=True).stdout.split()
    vals = np.array([float(x.replace("D", "e").replace("E", "e")) for x in out])
    return vals.reshape(len(cases), 16)


@pytest.fixture(scope="module")
def comparison():
    cases = _cases()
    with tempfile.TemporaryDirectory() as td:
        exe = Path(td) / "fabq"
        subprocess.run(
            ["gfortran", "-ffixed-line-length-132", f"-I{_ABA_INC}",
             str(_DRIVER), str(_ABQ_SRC), "-o", str(exe)], check=True)
        fort = _fortran_results(cases, exe)
    def rel(py, fo):
        """max |py - fo| scaled by the magnitude of the Fortran result."""
        py, fo = np.atleast_1d(py), np.atleast_1d(fo)
        return float(np.abs(py - fo).max() / max(np.abs(fo).max(), TINY))

    rows = []
    for k, (name, F, Fv, p) in enumerate(cases):
        sv, fvn, det = ms.stress_core(F, Fv, *p)
        degenerate = fort[k][15] <= 1.0000001e-15      # the detFe clamp fired
        rows.append((name, rel(sv, fort[k][:6]), rel(fvn.reshape(9), fort[k][6:15]),
                     rel(det, fort[k][15]), degenerate))
    return rows


@pytest.mark.parametrize("idx", range(len(_cases())))
def test_python_core_matches_fortran(comparison, idx):
    name, r_sig, r_fv, r_det, _ = comparison[idx]
    assert r_sig < RTOL, f"{name}: relative |Δσ| = {r_sig:.3e}"
    assert r_fv < RTOL, f"{name}: relative |ΔFv| = {r_fv:.3e}"
    assert r_det < RTOL, f"{name}: relative |ΔdetFe| = {r_det:.3e}"


def test_worst_case_summary(comparison):
    worst = max(max(r[1], r[2], r[3]) for r in comparison)
    assert worst < RTOL, (
        f"worst relative discrepancy over {len(comparison)} cases: {worst:.3e}")


def test_battery_covers_both_regimes(comparison):
    """The battery must exercise well-conditioned states *and* the degenerate
    detFe-clamp path — otherwise the equivalence claim is only half-tested."""
    degenerate = sum(1 for r in comparison if r[4])
    assert degenerate >= 1, "no degenerate (clamped) state exercised"
    assert len(comparison) - degenerate >= 10, "too few well-conditioned states"
