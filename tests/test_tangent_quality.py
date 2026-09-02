"""The delivered tangent, held against exact AD across the parameter range.

BIOFILM_GROWTH_VISCO_V01 returns a consistent tangent obtained by
finite-difference perturbation of F. That choice is the most likely cause of
the failure mode that would land on us during integration: a tangent that is
subtly wrong does not produce a wrong answer, it produces a Newton iteration
that will not converge -- and the partner group, seeing their solve stall after
swapping in our routine, has no way to tell whether the fault is the tangent,
the field solve, or their own wiring.

`material_jax.dsde_exact` differentiates the same perturbation exactly
(`jax.jacfwd`), so it is an independent reference rather than a second
finite-difference. Previously the two had been compared at a single state;
this sweeps the range the routine is actually expected to see: the biofilm
fraction across three decades of stiffness, Poisson ratios up to
near-incompressible, growth up to 0.35, both material paths, elastic and
viscous, and steps up to the guard threshold.

Requires gfortran and jax; skipped where either is absent.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("JAX_ENABLE_X64", "1")

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_FC = shutil.which("gfortran")
_CC = shutil.which("cc") or shutil.which("gcc")
jax = pytest.importorskip("jax")
pytestmark = pytest.mark.skipif(_FC is None or _CC is None,
                                reason="gfortran/cc unavailable")

sys.path.insert(0, str(_AU / "coupling"))
import material_jax as mj  # noqa: E402

I3 = np.eye(3)

# material_jax orders shear 12,13,23 (Abaqus); the Fortran orders it 12,23,13
# (ANSYS). Comparing the 6x6 tangents without this permutation compares
# mismatched rows and columns, and a stress check on a diagonal deformation
# will not catch it because the shear components are zero there.
_TO_ANSYS = [0, 1, 2, 3, 5, 4]


def _ansys(D):
    return np.asarray(D, dtype=float)[np.ix_(_TO_ANSYS, _TO_ANSYS)]


def _mr(E, nu, c01r):
    mu = E / (2.0 * (1.0 + nu))
    K = E / (3.0 * (1.0 - 2.0 * nu))
    c10 = 0.5 * mu / (1.0 + c01r)
    return c10, c10 * c01r, 2.0 / K


@pytest.fixture(scope="module")
def wrap():
    tmp = Path(tempfile.mkdtemp())

    def sh(*a):
        subprocess.run([str(x) for x in a], check=True, cwd=tmp)

    sh(_FC, "-c", "-ffixed-line-length-132", "-J", tmp,
       _AU / "coupling/usermat_py_hook.f", "-o", tmp / "h.o")
    sh(_FC, "-c", "-ffixed-line-length-132", "-I", tmp,
       _AU / "usermat_biofilm.f", "-o", tmp / "c.o")
    sh(_CC, "-c", "-fPIC", _AU / "coupling/biofilm_py_eval.c", "-o", tmp / "s.o")
    exe = tmp / "wrap"
    sh(_FC, "-ffixed-line-length-132", "-I", tmp,
       _AU / "crosscheck/wrapper_driver.f", _AU / "biofilm_material_v01.f",
       tmp / "c.o", tmp / "h.o", tmp / "s.o", "-o", exe)
    return exe


def _run(exe, F, bio, nu, alpha, eta, dt, c01r, mt):
    txt = (" ".join(f"{F[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
           " ".join(f"{I3[i, j]:.17e}" for i in range(3) for j in range(3)) + "\n" +
           f"{1000.0:.17e} {1.0:.17e} {nu:.17e} {nu:.17e} {bio:.17e} "
           f"{alpha:.17e} {eta:.17e} {dt:.17e} {c01r:.17e} {mt:.1f}\n")
    r = subprocess.run([str(exe)], input=txt, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"}, timeout=60)
    assert r.returncode == 0, r.stderr
    v = [float(x) for x in r.stdout.split()]
    return int(v[15]), np.array(v[17:53]).reshape(6, 6)


F_STATES = {
    "near-identity": np.array([[1.01, 0.0, 0.0], [0.0, 0.995, 0.0], [0.0, 0.0, 1.0]]),
    "shear+stretch": np.array([[1.06, 0.02, 0.0], [0.0, 0.97, 0.01], [0.0, 0.0, 0.98]]),
    "strong shear": np.array([[1.10, 0.12, 0.03], [0.05, 0.92, 0.0], [0.02, 0.04, 1.05]]),
    "large stretch": np.array([[1.35, 0.0, 0.0], [0.0, 0.88, 0.0], [0.0, 0.0, 0.90]]),
}

# The measured bound is ~1.1e-6 over this grid. 1e-5 leaves room for a
# different gfortran without letting a real regression through: a tangent that
# had genuinely gone wrong would miss by orders of magnitude, not by 10x.
TOL = 1.0e-5


def test_the_voigt_orders_really_do_differ():
    """Guards the permutation above. If material_jax ever switches to ANSYS
    order, the sweep would silently start comparing permuted matrices and
    still pass, because the permutation would then be the error."""
    assert mj.VOIGT[4] == (0, 2) and mj.VOIGT[5] == (1, 2), (
        "material_jax is no longer in Abaqus shear order; _TO_ANSYS is wrong")


@pytest.mark.parametrize("state", list(F_STATES))
def test_finite_difference_tangent_matches_exact_ad(wrap, state):
    F = F_STATES[state]
    worst, where, n = 0.0, None, 0
    for bio in (0.0, 0.25, 1.0):
        for nu in (0.30, 0.45, 0.49):
            for alpha in (0.0, 0.20, 0.35):
                for c01r, mt in ((0.0, 0.0), (0.15, 1.0)):
                    E = bio * 1000.0 + (1.0 - bio) * 1.0
                    c10, c01, d1 = _mr(E, nu, c01r)
                    for eta, frac in ((0.0, None), (5.0, 0.3), (5.0, 0.45)):
                        dt = 1.0e-4 if eta == 0.0 else frac * eta / (2.0 * c10)
                        kc, dfd = _run(wrap, F, bio, nu, alpha, eta, dt, c01r, mt)
                        if kc:
                            continue          # refused by the step guard
                        dex = _ansys(mj.dsde_exact(
                            F, I3, (alpha, c10, c01, d1, eta, mt, dt)))
                        rel = np.max(np.abs(dfd - dex)) / np.max(np.abs(dex))
                        n += 1
                        if rel > worst:
                            worst, where = rel, (bio, nu, alpha, mt, eta)
    assert n > 0
    assert worst < TOL, (
        f"tangent error {worst:.2e} at (biofilm, nu, alpha, mtype, eta)={where}; "
        "a tangent this far from exact is a Newton-convergence problem for "
        "whoever calls this routine")
