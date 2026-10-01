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
_CC = shutil.which("cc") or shutil.which("gcc")
_HOOK = _AU / "coupling" / "usermat_py_hook.f"
_SHIM = _AU / "coupling" / "biofilm_py_eval.c"
_WS2 = ["-lws2_32"] if sys.platform == "win32" else []
pytestmark = pytest.mark.skipif(_FC is None or _CC is None,
                                reason="gfortran/cc unavailable")

MOCK = """\
      SUBROUTINE MOCKMAT(nProp, prop, ustatev, dTime, elemId,
     &   kDomIntPt, ldstep, isubst, Sdp_bio1_n, Sdp_locbio1_n,
     &   Sdp_bio2_n, Sdp_locbio2_n, Sdp_sumLocal, Sbio_GrowthConst,
     &   keycut)
      USE biofilm_py_bridge
      USE biofilm_split
      IMPLICIT NONE
      INTEGER keycut
      INTEGER nProp
      DOUBLE PRECISION prop(nProp), ustatev(100), dTime
      DOUBLE PRECISION Sbio_GrowthConst, Sdp_sumLocal
      DOUBLE PRECISION Sdp_bio1_n, Sdp_locbio1_n
      DOUBLE PRECISION Sdp_bio2_n, Sdp_locbio2_n
      INTEGER elemId, kDomIntPt, ldstep, isubst
      INCLUDE 'phi_mode_decl.inc'
      Sbio_GrowthConst = prop(5)
      INCLUDE 'phi_mode_exec.inc'
      END

      PROGRAM DRIVE
      IMPLICIT NONE
      DOUBLE PRECISION prop(28), ust(100,3), work(100), sg, dt
      DOUBLE PRECISION b1, b2
      INTEGER e, ie, s, it, ielem(3), mode, np, kc
      DATA ielem /1, 2, 38/
      prop = 0.0D0
      READ(*,*) mode, prop(7), dt, np
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      kc = 0
      prop(1) = 1.0D0
      prop(28) = DBLE(mode)
      ust = 0.0D0
      DO s = 1, 5
        DO ie = 1, 3
          e = ielem(ie)
          b1 = 0.1D0*ie + 0.05D0*s
          b2 = 0.02D0*s
          DO it = 1, 3
            work = ust(:,ie)
            CALL MOCKMAT(np, prop, work, dt, e, 1, 1, s, b1,
     &                   1.0D0+0.01D0*b1, b2, 1.0D0+0.01D0*b2,
     &                   1.0D0+0.005D0*(b1+b2), sg, kc)
          END DO
          ust(:,ie) = work
        END DO
      END DO
      PRINT '(3ES24.16E3)', ust(84,1), ust(84,2), ust(84,3)
      END
"""


_SPLIT = _CS / "split_rates.f"


def build_mock(source=None):
    """Compile the fragments into the mock, linked against the real bridge
    module and C shim -- the same objects the ANSYS build links -- and the
    split hand-over module. `source` replaces the default mock program."""
    tmp = Path(tempfile.mkdtemp())
    for f in ("phi_mode_decl.inc", "phi_mode_exec.inc"):
        shutil.copy(_CS / f, tmp / f)
    (tmp / "mock.f").write_text(source if source is not None else MOCK)
    subprocess.run([_FC, "-c", "-ffixed-line-length-132", "-J", str(tmp),
                    str(_HOOK), "-o", str(tmp / "hook.o")], check=True,
                   cwd=tmp)
    subprocess.run([_FC, "-c", "-ffixed-form", "-ffixed-line-length-72",
                    "-Werror=line-truncation", "-J", str(tmp),
                    str(_SPLIT), "-o", str(tmp / "split.o")], check=True,
                   cwd=tmp)
    subprocess.run([_CC, "-c", "-fPIC", str(_SHIM), "-o",
                    str(tmp / "shim.o")], check=True)
    exe = tmp / "mock"
    r = subprocess.run([_FC, "-ffixed-form", "-ffixed-line-length-72",
                        "-fcheck=bounds",
                        "-Werror=line-truncation", "-I", str(tmp),
                        str(tmp / "mock.f"), str(_AU / "growth_from_phi.f"),
                        str(tmp / "hook.o"), str(tmp / "split.o"),
                        str(tmp / "shim.o"),
                        "-o", str(exe)] + _WS2, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return tmp, exe


def run_mock(tmp, exe, mode, k_alpha, dt, nprop=28, theta=None, env=None,
             timeout=30):
    th = theta if theta is not None else [0.0] * 20
    line = f"{mode} {k_alpha} {dt} {nprop}\n" + "".join(
        f"{x:.17e}\n" for x in th)
    r = subprocess.run([str(exe)], input=line, capture_output=True,
                       text=True, cwd=tmp, timeout=timeout, env=env)
    return r


def _build_and_run(mode, k_alpha=0.5, dt=0.1, nprop=28):
    tmp, exe = build_mock()
    r = run_mock(tmp, exe, mode, k_alpha, dt, nprop)
    assert r.returncode == 0, r.stderr
    alphas = [float(x) for x in r.stdout.split()]
    return tmp, alphas


def test_the_fragments_compile_as_strict_72_column_fixed_form():
    _build_and_run(1)


@pytest.mark.parametrize("mode", [1, 2])
def test_one_and_two_species_pass_the_bring_up_checks(mode):
    tmp, _ = _build_and_run(mode)
    rows = ref.read_trace(tmp / "phi_trace.csv")
    c = ref.check_trace(rows, 0.5)
    assert c["once_per_increment"], "growth counted per iteration"
    assert c["eq36"] and c["carried"]
    for r in rows:
        phi = r["bio1"] if mode == 1 else r["bio1"] + r["bio2"]
        assert abs(r["phi_used"] - phi) < 1e-15


def test_only_the_sampled_elements_are_traced():
    tmp, _ = _build_and_run(1)
    elems = {r["elem"] for r in ref.read_trace(tmp / "phi_trace.csv")}
    assert elems == {1, 38}                      # stride 37: 2 is not sampled


def test_alpha_matches_eq36_summed_over_increments():
    _, alphas = _build_and_run(2, k_alpha=0.5, dt=0.1)
    # element 1: b1 = 0.1 + 0.05 s, b2 = 0.02 s, s = 1..5
    want = sum(0.5 * (0.1 + 0.05 * s + 0.02 * s) * 0.1 for s in range(1, 6))
    assert abs(alphas[0] - want) < 1e-14


def test_phi_mode_off_leaves_growth_alone():
    _, alphas = _build_and_run(0)
    assert alphas == [0.0, 0.0, 0.0]


def test_a_27_constant_deck_never_reads_prop_28():
    """The ecology decks pass 27 constants; the guard must keep the
    fragment from reading past the end of prop there. Bounds checking turns
    an out-of-range read into a runtime failure, so a pass means no read."""
    _, alphas = _build_and_run(1, nprop=27)
    assert alphas == [0.0, 0.0, 0.0]


def test_mode_3_takes_the_partners_alpha():
    """prop(28) = 3: alpha is the partner's own Sdp_sumLocal - 1, whatever
    k_alpha says -- their USSFin already integrates Eq. 36 as locbio."""
    tmp, alphas = _build_and_run(3, k_alpha=123.0)
    rows = ref.read_trace(tmp / "phi_trace.csv")
    assert rows and ref.partner_alpha_gap(rows) < 1e-15
    # element 1, last increment: b1 = 0.15 + 0.05*4 ... sumLocal - 1
    s = 5
    b1, b2 = 0.1 + 0.05 * s, 0.02 * s
    assert abs(alphas[0] - 0.005 * (b1 + b2)) < 1e-15


def test_points_past_the_phi_threshold_are_traced_even_off_stride():
    """Element 2 is not on the stride; in two-species mode its phi reaches
    0.55 by the last increment and must then appear."""
    tmp, _ = _build_and_run(2)
    rows = ref.read_trace(tmp / "phi_trace.csv")
    e2 = [r for r in rows if r["elem"] == 2]
    assert e2 and all(r["phi_used"] >= 0.5 for r in e2)
