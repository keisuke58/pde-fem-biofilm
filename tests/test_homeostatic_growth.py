"""Homeostatic-pressure growth law (prop(49:51) of phi_mode_exec.inc) in the
Abaqus composition UMAT, without Abaqus.

alpha_dot = k_alpha phi max(0, 1 - p/p_h), p = -tr(sigma)/3 at the start of the
increment (keio_wp2/README.ja.md; an assumption of this work). One integration
point, one species (prop(28) = 1, no material server), phi = 1, F = I (a fully
constrained cube). The stress of constrained growth has a closed form (the one
test_abaqus_composition_umat.py checks), so p(alpha) is known and:
  - prop(49) = 0 leaves Eq. 36 unchanged: alpha = k_alpha phi T;
  - with p_h > 0, alpha follows the recurrence
        alpha_{n+1} = alpha_n + k_alpha phi dt max(0, 1 - p(alpha_n)/p_h)
    exactly, and stops at alpha* with p(alpha*) = p_h;
  - the local-stiffness form p_h = p_ref (phi^2 + f) does the same with its p_h;
  - the driver hands every equilibrium iteration the stress at the start of
    the increment, as Abaqus and ANSYS do.
The same fragment is pasted into the partner's ANSYS usermat (paste_fragments.py).
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat"),
                str(ROOT / "abaqus_composition")]
import make_umat  # noqa: E402

FC, CC = shutil.which("gfortran"), shutil.which("cc") or shutil.which("gcc")
pytestmark = pytest.mark.skipif(FC is None or CC is None, reason="gfortran/cc unavailable")

DRIVER = """\
      PROGRAM DRIVE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION STRESS(6), SS(6), STATEV(100), SV(100), DDSDDE(6,6),
     1 DDSDDT(6), DRPLDE(6), STRAN(6), DSTRAN(6), TIME(2),
     2 PREDEF(2), DPRED(2), PROPS(51), COORDS(3), DROT(3,3),
     3 DFGRD0(3,3), DFGRD1(3,3)
      INTEGER JSTEP(4)
      CHARACTER*80 CMNAME
      READ(*,*) DT, NINC, PHI
      READ(*,*) (PROPS(I), I = 1, 51)
      SV = 0.0D0
      SS = 0.0D0
      DFGRD1 = 0.0D0
      DO I = 1, 3
        DFGRD1(I,I) = 1.0D0
      END DO
      DFGRD0 = DFGRD1
      PREDEF(1) = PHI
      PREDEF(2) = -1.0D30
      DPRED = 0.0D0
      JSTEP = 1
      DO K = 1, NINC
        DO IT = 1, 2
          STATEV = SV
          STRESS = SS
          PNEWDT = 1.0D0
          CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1      DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DT, TEMP,
     2      DTEMP, PREDEF, DPRED, CMNAME, 3, 3, 6, 100, PROPS, 51,
     3      COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, 1, 1, 0,
     4      0, JSTEP, K)
        END DO
        SV = STATEV
        SS = STRESS
        PRINT '(4(1X,ES24.16E3))', SV(84), SS(1), SV(95), SV(96)
      END DO
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

KALPHA, DT, NINC = 1e-3, 0.025, 400
E0, F_VOID, NU = 1e-5, 1e-3, 0.49


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


def run(build, pref=0.0, form=0, f=0.0, phi=1.0):
    p = [0.0] * 51
    p[0], p[6], p[27] = 1.0, KALPHA, 1.0              # prop(7) k_alpha, prop(28) = 1
    p[42], p[43], p[44], p[45] = E0, -F_VOID, NU, 0.3  # Klempt 2024 stiffness branch
    p[48], p[49], p[50] = pref, form, f
    inp = f"{DT} {NINC} {phi}\n" + " ".join(f"{x:.17e}" for x in p) + "\n"
    exe = build / ("drive.exe" if sys.platform == "win32" else "drive")
    r = subprocess.run([str(exe)], input=inp, capture_output=True, text=True, cwd=build, timeout=600)
    assert r.returncode == 0, r.stderr
    return [[float(x) for x in line.split()] for line in r.stdout.splitlines() if line.strip()]


def p_closed(alpha, phi=1.0):
    """-sigma_ii of constrained growth (test_abaqus_composition_umat.py)."""
    E = (phi * phi + F_VOID) * E0
    K, C10 = E / (3 * (1 - 2 * NU)), E / (4 * (1 + NU))
    J = (1 + alpha) ** -3
    return -(K * (J - 1) + 2 * C10 * (1 - (1 + alpha) ** 2) / J)


def alpha_star(ph):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if p_closed(mid) < ph else (lo, mid)
    return (lo + hi) / 2


def recurrence(ph):
    a, out = 0.0, []
    for _ in range(NINC):
        a = a + KALPHA * DT * max(0.0, 1 - p_closed(a) / ph)
        out.append(a)
    return out


def test_off_is_eq36(build):
    rows = run(build)
    assert rows[-1][0] == pytest.approx(KALPHA * DT * NINC, rel=1e-12)
    assert rows[-1][2] == 0.0 and rows[-1][3] == 0.0      # state 95/96 untouched


def test_stress_is_the_closed_form(build):
    for a, s11, *_ in run(build)[::50]:
        assert -s11 == pytest.approx(p_closed(a), rel=1e-10)


@pytest.mark.parametrize("form", [0, 1])
def test_growth_stops_at_p_h(build, form):
    pref = 2.5e-7
    ph = pref * (1.0 + F_VOID) if form == 1 else pref
    rows = run(build, pref=pref, form=form, f=F_VOID)
    ref = recurrence(ph)
    for (a, _, p, g), r in zip(rows, ref):
        assert a == pytest.approx(r, rel=1e-9, abs=1e-15)
        assert 0.0 <= g <= 1.0
    a_end, s_end, p_end, g_end = rows[-1]
    assert a_end == pytest.approx(alpha_star(ph), rel=1e-6)
    assert -s_end == pytest.approx(ph, rel=1e-5)          # p reaches p_h ...
    assert max(-r[1] for r in rows) <= ph * (1 + 1e-9)    # ... from below
    assert g_end < 1e-5
    assert a_end < 0.1 * KALPHA * DT * NINC               # far below Eq. 36's growth
