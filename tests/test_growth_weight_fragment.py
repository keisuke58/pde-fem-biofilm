"""Species-weighted growth law in the composition mode (prop(28) = 7,
prop(36) = d; two_way_growth_rate.py). One packed Gauss point, case 6 of
Klempt et al. 2026, mock usermat with a live material server. Checked:
  - without prop(36), or with d = 0, alpha is Eq. 36 bit for bit;
  - with d != 0 the alpha increment of a step is Eq. 36's times
    1 + d (2 chi_1 - 1), chi_1 from the point model's state after the step;
  - the composition itself does not change.
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
sys.path.insert(0, str(_AU / "apdl"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip("jax")
import material_server as ms                            # noqa: E402
import test_partner_callsite_fragment as tp             # noqa: E402

pytestmark = pytest.mark.skipif(tp._FC is None or tp._CC is None,
                                reason="gfortran/cc unavailable")

DT, S, CAP = 0.01, 0.5, 0.9

_SUB = tp.MOCK[:tp.MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM GWDRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(36), work(100), sg, dt
      INTEGER it, s, kc, nsteps, np
      prop = 0.0D0
      READ(*,*) dt, nsteps, np, prop(31), prop(32), prop(36)
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(7) = 1.0D-3
      prop(28) = 7.0D0
      prop(29) = 1.0D-2
      prop(30) = 0.5D0
      kc = 0
      work = 0.0D0
      DO s = 1, nsteps
        CALL MOCKMAT(np, prop, work, dt, 1, 1, 1, s, 1.0D0,
     &               1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
      END DO
      PRINT '(13ES26.17E3)', (work(71 + it), it = 1, 13)
      END
"""


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


def _run(server, build, nsteps, np_, d):
    tmp, exe = build
    line = (f"{DT} {nsteps} {np_} {S} {CAP} {d}\n"
            + "".join(f"{x:.17e}\n" for x in ms.ECOLOGY_CASE["theta"]))
    r = subprocess.run([str(exe)], input=line, capture_output=True, text=True,
                       cwd=tmp, timeout=600,
                       env=dict(os.environ,
                                BIOFILM_PY_PORT=str(server.server_address[1])))
    assert r.returncode == 0, r.stdout + r.stderr
    return np.array([float(x) for x in r.stdout.split()])   # g(1:12), alpha


def test_d_zero_is_eq36(server, build):
    base = _run(server, build, 3, 35, 0.0)
    assert np.array_equal(_run(server, build, 3, 36, 0.0), base)
    assert np.array_equal(_run(server, build, 3, 35, 0.5), base)   # prop(36) not passed


def test_one_step_factor(server, build):
    base = _run(server, build, 1, 36, 0.0)
    for d in (1 / 3, 0.5, -0.4):
        w = _run(server, build, 1, 36, d)
        assert np.array_equal(w[:12], base[:12])                    # composition unchanged
        chi = w[0] / (w[0] + w[1])
        assert w[12] == pytest.approx(base[12] * (1 + d * (2 * chi - 1)), rel=1e-12)
        assert w[12] != base[12]
