"""The shape-function gradient of the unified user element (p3/uel/unified_uel.f90, shape_grad) on elements that are
not axis-aligned bricks: until 17 Oct 2026 it was Jm^-T dN/dxi instead of Jm^-1 dN/dxi, right only for a symmetric
Jm. Checks: free growth u = xi X leaves no mechanics residual on a rotated and on a skewed hexahedron (as on the
aligned one); a uniform strain gives the same nodal-force norm on the rotated element as on the aligned one (the
element is objective); the discrete gradient of a linear field is exact on a skewed element; the closed-form tangent
still matches central differences there. With and without the mean dilatation.

Needs gfortran (builds libunified_uel.so in p3/uel if it is missing or older than the source) and scipy; skipped
without them. About 5 s.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
sys.path.insert(0, str(UEL))
pytest.importorskip("scipy")
pytestmark = pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")


@pytest.fixture(scope="module")
def dv(tmp_path_factory):
    so = UEL / "libunified_uel.so"
    srcs = (UEL / "unified_uel.f90", ROOT / "ansys_usermat" / "coupling" / "ecology_native.f")
    if not so.exists() or max(x.stat().st_mtime for x in srcs) > so.stat().st_mtime:   # rebuild when stale
        eco = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"
        subprocess.run(["gfortran", "-O2", "-fPIC", "-c", str(eco), "-o", "ecology_native.o"], cwd=UEL, check=True)
        subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-line-length-none", "-c", "unified_uel.f90", "-o", "unified_uel.o"],
                       cwd=UEL, check=True)
        subprocess.run(["gfortran", "-shared", "-o", "libunified_uel.so", "unified_uel.o", "ecology_native.o"], cwd=UEL, check=True)
    if not os.environ.get("BIOFILM_ECO_CASE"):
        from make_column_inp import eco_file
        p = tmp_path_factory.mktemp("eco") / "uel_eco.txt"
        eco_file(p)
        os.environ["BIOFILM_ECO_CASE"] = str(p)
    import driver
    driver.setup(np.zeros(20))
    return driver


X0 = np.array([[a, b, c] for a in (0.0, 1.0) for b in (0.0, 1.0) for c in (0.0, 1.0)])
EL = np.array([0, 4, 6, 2, 1, 5, 7, 3])                                   # node order of driver.column's hexahedron
cz, sz = np.cos(np.pi / 4), np.sin(np.pi / 4); cx, sx = np.cos(np.pi / 6), np.sin(np.pi / 6)
ROT = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]]) @ np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
SKEW = np.array([[1.0, 0.3, 0.0], [0.0, 1.0, 0.2], [0.1, 0.0, 1.0]])
SHAPES = {"aligned": np.eye(3), "rotated": ROT, "skewed": SKEW}


def element(dv, X, u_of_x, xi, md):
    """residual and tangent of one element with phi = 0.8, c = 1, growth xi stored, u = u_of_x(X) at the nodes"""
    from driver import F, NDN, NDOF, NGP, NSV, lib, props_base
    lib.uel_set_finite(1, 1 if md else 0)
    props = F(props_base(nu=0.3))
    sv = F(np.zeros((NSV, NGP))); sv[12, :] = xi
    U = np.zeros(NDOF); U[3::NDN] = 0.8; U[4::NDN] = 1.0
    for a in range(8):
        U[NDN * a:NDN * a + 3] = u_of_x(X[a])
    R = F(np.zeros(NDOF)); K = F(np.zeros((NDOF, NDOF)))
    lib.uel_element(F(U), F(U), F(X.T), sv, props, 1e-30, R, K)
    return np.array(R), np.array(K), U, sv, props


@pytest.mark.parametrize("md", [False, True])
def test_free_growth_has_no_residual_on_rotated_and_skewed_elements(dv, md):
    for name, A in SHAPES.items():
        X = X0[EL] @ A.T
        R, _, _, _, _ = element(dv, X, lambda x: 0.1 * x, 0.1, md)
        assert np.abs(R.reshape(8, 5)[:, :3]).max() < 1e-14, name


def test_uniform_strain_force_is_objective(dv):
    """the nodal forces of the rotated element under the rotated strain have the norm of the aligned element's"""
    B = np.array([[0.02, 0.01, 0.0], [0.0, -0.01, 0.005], [0.003, 0.0, 0.015]])
    norms = {}
    for name in ("aligned", "rotated"):
        A = SHAPES[name]; X = X0[EL] @ A.T
        R, _, _, _, _ = element(dv, X, lambda x, A=A: A @ B @ A.T @ x, 0.0, True)
        norms[name] = np.linalg.norm(R.reshape(8, 5)[:, :3])
    assert norms["aligned"] > 1e-3
    assert norms["rotated"] == pytest.approx(norms["aligned"], rel=1e-12)


def test_tangent_matches_differences_on_skewed_element(dv):
    from driver import F, NDN, NDOF, lib
    X = X0[EL] @ SKEW.T
    B = np.array([[0.05, 0.02, 0.0], [0.0, -0.03, 0.01], [0.01, 0.0, 0.04]])
    R, K, U, sv, props = element(dv, X, lambda x: B @ x, 0.05, True)
    Kfd = np.zeros((NDOF, NDOF))
    for m in range(NDOF):
        h = 1e-6 * max(1.0, abs(U[m]))
        Rp = F(np.zeros(NDOF)); Rm = F(np.zeros(NDOF)); Kd = F(np.zeros((NDOF, NDOF)))
        Up = U.copy(); Up[m] += h; Um = U.copy(); Um[m] -= h
        lib.uel_element(F(Up), F(U), F(X.T), sv, props, 1e-30, Rp, Kd)
        lib.uel_element(F(Um), F(U), F(X.T), sv, props, 1e-30, Rm, Kd)
        Kfd[:, m] = (np.array(Rp) - np.array(Rm)) / (2 * h)
    mech = np.arange(NDOF) % NDN < 3
    assert np.abs(K[np.ix_(mech, mech)] - Kfd[np.ix_(mech, mech)]).max() < 1e-7 * np.abs(K[np.ix_(mech, mech)]).max()
