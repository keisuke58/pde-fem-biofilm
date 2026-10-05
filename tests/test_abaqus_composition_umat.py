"""The Abaqus composition UMAT (abaqus_composition/make_umat.py) without Abaqus.

The generated file is compiled with gfortran next to a small driver that calls
UMAT the way Abaqus does (two equilibrium iterations per increment from the
converged state), linked against the real C shim, and run against a live
material server (case 6 of Klempt et al. 2026). Checked, one integration
point, phi = 1 (capped at 0.9 in the point model), F = I, 40 increments of 0.025:
  - the share phi_1/(phi_1+phi_2) is the ANSYS value 0.0209 (elem 220, 1 Oct);
  - alpha - 1 = k_alpha phi T* (Eq. 36);
  - the stress is the closed form of constrained growth plus the spherical
    term of DEVIATOR_SCALING_FINDING.md (the stress core is shared with ANSYS);
  - c = 0.3 with c_ref = 1 changes the composition (local nutrient reaches it).
On IKMHIWI03 the same file ran in Abaqus 2024 (abaqus_composition/run_comp.ps1).
"""
import os
import shutil
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat"),
                str(ROOT / "abaqus_composition")]
pytest.importorskip("jax")
import material_server as ms  # noqa: E402
import make_umat  # noqa: E402

FC, CC = shutil.which("gfortran"), shutil.which("cc") or shutil.which("gcc")
pytestmark = pytest.mark.skipif(FC is None or CC is None, reason="gfortran/cc unavailable")

DRIVER = """\
      PROGRAM DRIVE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION STRESS(6), STATEV(100), SV(100), DDSDDE(6,6),
     1 DDSDDT(6), DRPLDE(6), STRAN(6), DSTRAN(6), TIME(2),
     2 PREDEF(2), DPRED(2), PROPS(46), COORDS(3), DROT(3,3),
     3 DFGRD0(3,3), DFGRD1(3,3)
      INTEGER JSTEP(4)
      CHARACTER*80 CMNAME
      READ(*,*) DT, NINC, PHI, C
      READ(*,*) (PROPS(I), I = 1, 46)
      SV = 0.0D0
      DFGRD1 = 0.0D0
      DO I = 1, 3
        DFGRD1(I,I) = 1.0D0
      END DO
      DFGRD0 = DFGRD1
      PREDEF(1) = PHI
      PREDEF(2) = C
      DPRED = 0.0D0
      JSTEP = 1
      DO K = 1, NINC
        DO IT = 1, 2
          STATEV = SV
          PNEWDT = 1.0D0
          CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1      DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DT, TEMP,
     2      DTEMP, PREDEF, DPRED, CMNAME, 3, 3, 6, 100, PROPS, 46,
     3      COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, 1, 1, 0,
     4      0, JSTEP, K)
        END DO
        SV = STATEV
      END DO
      PRINT '(6(1X,ES24.16E3))', STRESS
      PRINT '(3(1X,ES24.16E3))', SV(72), SV(73), SV(84)
      END

      SUBROUTINE XIT
      STOP 1
      END

      SUBROUTINE GETOUTDIR(D, L)
      CHARACTER*(*) D
      INTEGER L
      D = '.'
      L = 1
      END
"""


@pytest.fixture(scope="module")
def build():
    tmp = Path(tempfile.mkdtemp())
    (tmp / "umat_comp.for").write_text(make_umat.build(), encoding="latin-1")
    (tmp / "ABA_PARAM.INC").write_text("      IMPLICIT REAL*8(A-H,O-Z)\n")
    (tmp / "drive.f").write_text(DRIVER)
    subprocess.run([FC, "-c", "-ffixed-form", "-ffixed-line-length-132", "-x", "f77-cpp-input",
                    "umat_comp.for", "-o", "umat_comp.o"], cwd=tmp, check=True)
    subprocess.run([CC, "-c", str(ROOT / "ansys_usermat" / "coupling" / "biofilm_py_eval.c"),
                    "-o", "shim.o"], cwd=tmp, check=True)
    ws2 = ["-lws2_32"] if sys.platform == "win32" else []
    r = subprocess.run([FC, "-ffixed-line-length-132", "drive.f", "umat_comp.o", "shim.o",
                        "-o", "drive"] + ws2, cwd=tmp, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return tmp


@pytest.fixture(scope="module")
def server():
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    ms.set_case("2sp_case6")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        yield srv
    finally:
        srv.shutdown()
        srv.server_close()
        ms.set_case(None)
        ms.set_active_species(5)


def run(build, server, c=1.0, cref=0.0):
    p = [0.0] * 46
    p[0], p[6] = 1.0, 1e-3
    p[7:27] = list(ms.ECOLOGY_CASE["theta"])
    p[27], p[28], p[29], p[30], p[31], p[32] = 7.0, 0.01, 0.5, 0.15, 0.9, cref
    p[42], p[43], p[44], p[45] = 1e-5, -1e-3, 0.49, 0.3
    inp = f"0.025 40 1.0 {c}\n" + " ".join(f"{x:.17e}" for x in p) + "\n"
    exe = build / ("drive.exe" if sys.platform == "win32" else "drive")
    r = subprocess.run([str(exe)], input=inp, capture_output=True, text=True, cwd=build, timeout=600,
                       env=dict(os.environ, BIOFILM_PY_PORT=str(server.server_address[1])))
    assert r.returncode == 0, r.stderr
    v = [float(x) for x in r.stdout.split()]
    return v[:6], v[6:]


def test_share_alpha_and_stress(build, server):
    s, (p1, p2, alpha) = run(build, server)
    assert abs(p1 / (p1 + p2) - 0.0209) < 5e-4
    assert abs(alpha - 1e-3) < 1e-15
    E, nu = (1.0 + 1e-3) * 1e-5, 0.49           # Klempt branch: (phi^2 + f) E
    K, C10 = E / (3 * (1 - 2 * nu)), E / (4 * (1 + nu))
    J = (1 + alpha) ** -3
    closed = K * (J - 1)
    spurious = 2 * C10 * (1 - (1 + alpha) ** 2) / J      # DEVIATOR_SCALING_FINDING.md sec. 8
    for i in range(3):
        assert s[i] == pytest.approx(closed + spurious, rel=1e-10)
    assert max(abs(x) for x in s[3:]) < 1e-20


def test_local_nutrient_reaches_the_point_model(build, server):
    _, (a1, a2, _) = run(build, server)
    _, (b1, b2, _) = run(build, server, c=0.3, cref=1.0)
    assert b1 / (b1 + b2) - a1 / (a1 + a2) > 0.1
