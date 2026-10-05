"""prop(35): start share for points the biofilm reaches later, in the composition
mode (prop(28) = 7) of the partner-element fragment.

Two Gauss points through the real fragment (mock material call, real bridge and
C shim, live material server, case 6 of Klempt et al. 2026): point A holds
biofilm from the start (phi = 1), point B has none for the first steps
(phi = 0 < phi_min) and is reached by the biofilm afterwards.
Checked:
  - without prop(35) both start from prop(30), as before;
  - with prop(35) > 0, B starts from prop(35) and A from prop(30);
  - A's state is bit for bit the same with and without prop(35).
"""
import os
import subprocess
import sys
import threading
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
sys.path.insert(0, str(_AU / "coupling"))
sys.path.insert(0, str(_AU))
sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip("jax")
import material_server as ms                            # noqa: E402
import test_partner_callsite_fragment as tp             # noqa: E402

pytestmark = pytest.mark.skipif(tp._FC is None or tp._CC is None,
                                reason="gfortran/cc unavailable")

_SUB = tp.MOCK[:tp.MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM LATEDRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(35), wa(100), wb(100), sg, dt, phib
      INTEGER it, s, kc, nsteps, nlate
      prop = 0.0D0
      READ(*,*) dt, nsteps, nlate, prop(30), prop(35)
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(7) = 1.0D-3
      prop(28) = 7.0D0
      prop(29) = 1.0D-2
      prop(31) = 0.5D0
      prop(32) = 0.9D0
      kc = 0
      wa = 0.0D0
      wb = 0.0D0
      DO s = 1, nsteps
        CALL MOCKMAT(35, prop, wa, dt, 1, 1, 1, s, 1.0D0,
     &               1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
        phib = 0.0D0
        IF (s .GT. nlate) phib = 1.0D0
        CALL MOCKMAT(35, prop, wb, dt, 2, 1, 1, s, phib,
     &               1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
        IF (s .EQ. 1) PRINT '(2ES26.17E3)', wb(72), wb(73)
      END DO
      PRINT '(12ES26.17E3)', (wa(71 + it), it = 1, 12)
      PRINT '(12ES26.17E3)', (wb(71 + it), it = 1, 12)
      END
"""

DT, NSTEPS, NLATE = 0.01, 3, 1


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


@pytest.fixture(scope="module")
def build():
    return tp.build_mock(DRIVER)


def _run(server, build, chi0, chi_late):
    tmp, exe = build
    th = ms.ECOLOGY_CASE["theta"]
    line = f"{DT} {NSTEPS} {NLATE} {chi0} {chi_late}\n" + "".join(f"{x:.17e}\n" for x in th)
    r = subprocess.run([str(exe)], input=line, capture_output=True, text=True, cwd=tmp,
                       timeout=600, env=dict(os.environ, BIOFILM_PY_PORT=str(server.server_address[1])))
    assert r.returncode == 0, r.stdout + r.stderr
    rows = [np.array([float(x) for x in l.split()]) for l in r.stdout.strip().splitlines()]
    return rows[0], rows[1], rows[2]          # B's stored phi_1, phi_2 after step 1; A; B at the end


def test_without_prop35_both_start_from_prop30(server, build):
    b1, a, b = _run(server, build, 0.5, 0.0)
    assert b1[0] / (b1[0] + b1[1]) == pytest.approx(0.5, abs=1e-15)


def test_prop35_sets_only_the_late_points(server, build):
    b1, a, b = _run(server, build, 0.5, 0.2)
    assert b1[0] / (b1[0] + b1[1]) == pytest.approx(0.2, abs=1e-15)
    _, a0, b0 = _run(server, build, 0.5, 0.0)
    assert np.array_equal(a, a0)              # the seed point is untouched
    assert not np.array_equal(b, b0)          # the late point starts elsewhere
