"""Pre-flight for prop(28) = 8: composition from an independent point model.

Klempt et al. 2026 formulate a material point model and state that phi0 is a
numerical auxiliary, not a physical void. In this mode the point model is
therefore NOT scaled to the field: at each Gauss point it runs on its own from
the substep in which the field first reaches phi_min, and supplies only the
composition; the amount stays the field's. Mock usermat + live server
(--case 2sp_case3), three points: one growing from 0.35, one in the void that
is filled from substep 3, one in the seed (phi_3D = 1).
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
import composition_reference as cref                    # noqa: E402
import material_server as ms                            # noqa: E402
from test_partner_callsite_fragment import (            # noqa: E402
    _CC, _FC, MOCK, build_mock)

pytestmark = pytest.mark.skipif(_FC is None or _CC is None,
                                reason="gfortran/cc unavailable")

DT, NSTEPS, K_ALPHA, PHI_MIN, S = 0.01, 6, 1.0e-3, 0.01, 0.5
PHI1, PHI2 = 0.2, 0.2
ELEMS = (1, 38, 75)


def phi3(e, s):
    if e == 1:
        return 0.3 + 0.05 * s
    if e == 38:
        return 0.0 if s < 3 else 0.2
    return 1.0


_SUB = MOCK[:MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM AGEDRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(32), ust(100,3), work(100), sg, dt, b1
      INTEGER it, s, ie, e, kc, nsteps, ielem(3)
      DATA ielem /1, 38, 75/
      prop = 0.0D0
      READ(*,*) prop(7), dt, nsteps, prop(29), prop(30), prop(31),
     &          prop(32)
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(28) = 8.0D0
      kc = 0
      ust = 0.0D0
      DO s = 1, nsteps
        DO ie = 1, 3
          e = ielem(ie)
          IF (e .EQ. 1) b1 = 0.3D0 + 0.05D0 * s
          IF (e .EQ. 38) THEN
            b1 = 0.0D0
            IF (s .GE. 3) b1 = 0.2D0
          END IF
          IF (e .EQ. 75) b1 = 1.0D0
          DO it = 1, 3
            work = ust(:,ie)
            CALL MOCKMAT(32, prop, work, dt, e, 1, 1, s, b1,
     &                   1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
          END DO
          ust(:,ie) = work
        END DO
      END DO
      PRINT '(3ES24.16E3)', ust(84,1), ust(84,2), ust(84,3)
      END
"""


@pytest.fixture(scope="module")
def run():
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    ms.set_case("2sp_case3")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        tmp, exe = build_mock(DRIVER)
        th = ms.ECOLOGY_CASE["theta"]
        line = (f"{K_ALPHA} {DT} {NSTEPS} {PHI_MIN} {PHI1} {S} {PHI2}\n"
                + "".join(f"{x:.17e}\n" for x in th))
        r = subprocess.run([str(exe)], input=line, capture_output=True,
                           text=True, cwd=tmp, timeout=600,
                           env=dict(os.environ,
                                    BIOFILM_PY_PORT=str(srv.server_address[1])))
        assert r.returncode == 0, r.stdout + r.stderr
        yield (cref.read_comp_trace(tmp / "age_trace.csv"), th,
               dict(ms.ECOLOGY_CASE["hp"]))
    finally:
        srv.shutdown()
        srv.server_close()
        ms.set_case(None)
        ms.set_active_species(5)


def _last(rows, e, s):
    return np.asarray([r for r in rows if r["elem"] == e
                       and r["isubst"] == s][-1]["g_new"])


def test_trace_passes_every_check(run):
    rows, th, hp = run
    assert {r["elem"] for r in rows} == set(ELEMS)
    c = cref.check_age_trace(rows, K_ALPHA, th, hp)
    for k in ("once_per_increment", "carried", "held", "one_call", "eq36",
              "replay"):
        assert c[k], (k, c)


@pytest.mark.parametrize("e", ELEMS)
def test_each_point_equals_the_stand_alone_scheme(run, e):
    rows, th, hp = run
    want = cref.reference_age([phi3(e, s) for s in range(1, NSTEPS + 1)],
                              th, hp, DT, phi_min=PHI_MIN, s=S,
                              phi_init=(PHI1, PHI2))
    for s in range(1, NSTEPS + 1):
        assert np.array_equal(want[s - 1], _last(rows, e, s)), (e, s)


def test_the_amount_does_not_change_the_composition(run):
    """Points 1 and 75 start together; their amounts differ (0.35.. vs 1)
    but the point model is not scaled, so their states are identical."""
    rows, *_ = run
    for s in range(1, NSTEPS + 1):
        assert np.array_equal(_last(rows, 1, s), _last(rows, 75, s))


def test_a_later_arrival_follows_the_same_path_later(run):
    """Point 38 is filled at substep 3: its state then is point 1's at 1."""
    rows, *_ = run
    assert not _last(rows, 38, 2).any()                 # not started
    for s in range(3, NSTEPS + 1):
        assert np.array_equal(_last(rows, 38, s), _last(rows, 1, s - 2))
