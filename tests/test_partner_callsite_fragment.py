"""Pre-flight for the call-site fragments before they meet the partner's file.

The fragments in ansys_usermat/apdl/callsite/ are pasted into the partner's
usermat on IKMHIWI03, where they cannot be tested here. So they are compiled
into a mock routine that has the same variables in scope, with STRICT
72-column fixed form and IMPLICIT NONE -- what ifort assumes for a .F file --
and driven the way ANSYS drives a usermat: several equilibrium iterations per
increment, each starting from the last CONVERGED state, which is committed
only after the last iteration. The trace they write is then read by the very
checker the real run will be judged with.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
_CS = _AU / "apdl" / "callsite"
sys.path.insert(0, str(_AU))
import one_species_reference as ref   # noqa: E402

_FC = shutil.which("gfortran")
pytestmark = pytest.mark.skipif(_FC is None, reason="gfortran unavailable")

MOCK = """\
      SUBROUTINE MOCKMAT(prop, ustatev, dTime, elemId, kDomIntPt,
     &   ldstep, isubst, Sdp_bio1_n, Sdp_locbio1_n, Sdp_bio2_n,
     &   Sdp_locbio2_n, sGrowth)
      IMPLICIT NONE
      DOUBLE PRECISION prop(30), ustatev(100), dTime, sGrowth
      DOUBLE PRECISION Sdp_bio1_n, Sdp_locbio1_n
      DOUBLE PRECISION Sdp_bio2_n, Sdp_locbio2_n
      INTEGER elemId, kDomIntPt, ldstep, isubst
      INCLUDE 'phi_mode_decl.inc'
      sGrowth = prop(5)
      INCLUDE 'phi_mode_exec.inc'
      END

      PROGRAM DRIVE
      IMPLICIT NONE
      DOUBLE PRECISION prop(30), ust(100,3), work(100), sg, dt
      DOUBLE PRECISION b1, b2
      INTEGER e, ie, s, it, ielem(3), mode
      DATA ielem /1, 2, 998/
      prop = 0.0D0
      READ(*,*) mode, prop(7), dt
      prop(1) = 1.0D0
      prop(6) = DBLE(mode)
      ust = 0.0D0
      DO s = 1, 5
        DO ie = 1, 3
          e = ielem(ie)
          b1 = 0.1D0*ie + 0.05D0*s
          b2 = 0.02D0*s
          DO it = 1, 3
            work = ust(:,ie)
            CALL MOCKMAT(prop, work, dt, e, 1, 1, s, b1, 0.5D0*b1,
     &                   b2, 0.5D0*b2, sg)
          END DO
          ust(:,ie) = work
        END DO
      END DO
      PRINT '(3ES24.16E3)', ust(84,1), ust(84,2), ust(84,3)
      END
"""


def _build_and_run(mode, k_alpha=0.5, dt=0.1):
    tmp = Path(tempfile.mkdtemp())
    for f in ("phi_mode_decl.inc", "phi_mode_exec.inc"):
        shutil.copy(_CS / f, tmp / f)
    (tmp / "mock.f").write_text(MOCK)
    exe = tmp / "mock"
    r = subprocess.run([_FC, "-ffixed-form", "-ffixed-line-length-72",
                        "-Werror=line-truncation", "-I", str(tmp),
                        str(tmp / "mock.f"), str(_AU / "growth_from_phi.f"),
                        "-o", str(exe)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    line = f"{mode} {k_alpha} {dt}\n"
    r = subprocess.run([str(exe)], input=line, capture_output=True,
                       text=True, cwd=tmp, timeout=30)
    assert r.returncode == 0, r.stderr
    alphas = [float(x) for x in r.stdout.split()]
    return tmp, alphas


def test_the_fragments_compile_as_strict_72_column_fixed_form():
    _build_and_run(2)


@pytest.mark.parametrize("mode", [2, 3])
def test_one_and_two_species_pass_the_bring_up_checks(mode):
    tmp, _ = _build_and_run(mode)
    rows = ref.read_trace(tmp / "phi_trace.csv")
    c = ref.check_trace(rows, 0.5)
    assert c["once_per_increment"], "growth counted per iteration"
    assert c["eq36"] and c["carried"]
    for r in rows:
        phi = r["bio1"] if mode == 2 else r["bio1"] + r["bio2"]
        assert abs(r["phi_used"] - phi) < 1e-15


def test_only_the_sampled_elements_are_traced():
    tmp, _ = _build_and_run(2)
    elems = {r["elem"] for r in ref.read_trace(tmp / "phi_trace.csv")}
    assert elems == {1, 998}                     # element 2 is not sampled


def test_alpha_matches_eq36_summed_over_increments():
    _, alphas = _build_and_run(3, k_alpha=0.5, dt=0.1)
    # element 1: b1 = 0.1 + 0.05 s, b2 = 0.02 s, s = 1..5
    want = sum(0.5 * (0.1 + 0.05 * s + 0.02 * s) * 0.1 for s in range(1, 6))
    assert abs(alphas[0] - want) < 1e-14


def test_mode_one_and_below_leave_growth_alone():
    _, alphas = _build_and_run(1)
    assert alphas == [0.0, 0.0, 0.0]
