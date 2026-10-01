"""Pre-flight for the composition mode (prop(28) = 7).

The 3D field decides the amount phi_3D = bio1 + bio2 and is not changed;
the point model (Klempt et al. 2026, case 2sp_case3 via the server's --case)
decides only the composition. The fragments run in the mock usermat against
a live server; three points get prescribed phi_3D(t): a growing one, one in
the void that is filled later (below phi_min first), and one in the seed
region (phi_3D = 1, clamped). Judged by composition_reference: the trace
checks, a bit-exact replay of every server call, and the stand-alone
rescaled scheme (the exact target for a gradient-free region).
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

DT, NSTEPS, K_ALPHA, PHI_MIN, CHI1 = 0.01, 6, 1.0e-3, 0.01, 0.5
ELEMS = (1, 38, 75)                      # all on the trace stride


def phi3(e, s):
    if e == 1:
        return 0.3 + 0.05 * s            # growing
    if e == 38:
        return 0.0 if s < 3 else 0.2     # void, then filled
    return 1.0                           # seed region, clamped


_SUB = MOCK[:MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM COMPDRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(30), ust(100,3), work(100), sg, dt, b1
      INTEGER it, s, ie, e, kc, nsteps, ielem(3)
      DATA ielem /1, 38, 75/
      prop = 0.0D0
      READ(*,*) prop(7), dt, nsteps, prop(29), prop(30)
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(28) = 7.0D0
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
            CALL MOCKMAT(30, prop, work, dt, e, 1, 1, s, b1,
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
        line = (f"{K_ALPHA} {DT} {NSTEPS} {PHI_MIN} {CHI1}\n"
                + "".join(f"{x:.17e}\n" for x in th))
        r = subprocess.run([str(exe)], input=line, capture_output=True,
                           text=True, cwd=tmp, timeout=600,
                           env=dict(os.environ,
                                    BIOFILM_PY_PORT=str(srv.server_address[1])))
        assert r.returncode == 0, r.stdout + r.stderr
        yield (cref.read_comp_trace(tmp / "comp_trace.csv"), th,
               dict(ms.ECOLOGY_CASE["hp"]),
               [float(x) for x in r.stdout.split()])
    finally:
        srv.shutdown()
        srv.server_close()
        ms.set_case(None)
        ms.set_active_species(5)


def test_trace_passes_every_check(run):
    rows, th, hp, _ = run
    assert {r["elem"] for r in rows} == set(ELEMS)
    c = cref.check_comp_trace(rows, K_ALPHA, th, hp)
    for k in ("once_per_increment", "carried", "amount", "one_call", "eq36",
              "replay"):
        assert c[k], (k, c)


@pytest.mark.parametrize("e", ELEMS)
def test_each_point_equals_the_stand_alone_scheme(run, e):
    rows, th, hp, _ = run
    want = cref.reference([phi3(e, s) for s in range(1, NSTEPS + 1)], th,
                          hp, DT, chi1=CHI1, phi_min=PHI_MIN)
    got = [[r for r in rows if r["elem"] == e and r["isubst"] == s][-1]
           for s in range(1, NSTEPS + 1)]
    for w, g in zip(want, got):
        assert np.array_equal(w, np.asarray(g["g_new"])), (e, g["isubst"])


def test_the_void_point_holds_until_filled(run):
    rows, *_ = run
    e38 = [r for r in rows if r["elem"] == 38]
    assert all(r["hit"] == 2 for r in e38 if r["isubst"] < 3)
    assert all(r["hit"] != 2 for r in e38 if r["isubst"] >= 3)
    first = [r for r in e38 if r["isubst"] == 1][0]
    assert first["g_new"][0] == CHI1 and first["g_new"][1] == 1 - CHI1


def test_the_amount_is_the_3d_field_and_growth_is_eq36(run):
    rows, _, _, alphas = run
    for r in rows:
        if r["hit"] != 2:
            assert abs(r["g_new"][0] + r["g_new"][1] - r["phi3"]) <= 4e-16
            assert r["g_new"][5] == 1.0 - r["phi3"]
    want = sum(K_ALPHA * phi3(1, s) * DT for s in range(1, NSTEPS + 1))
    assert abs(alphas[0] - want) < 1e-16


def test_server_refuses_a_deck_for_another_case():
    ms.set_case("2sp_case3")
    try:
        bad = list(ms.ECOLOGY_CASE["theta"])
        bad[1] = -1.0                                    # case 6's a12
        with pytest.raises(ValueError, match="does not match case"):
            ms.evaluate_ecology({"g": list(cref.seed()), "theta": bad,
                                 "dt_h": 0.01, "n_sub": 1})
    finally:
        ms.set_case(None)
        ms.set_active_species(5)
