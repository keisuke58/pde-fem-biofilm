"""Pre-flight for the transport/reaction split (prop(28) = 5, 6).

The usermat fragments run the point model at every Gauss point and hand the
reaction rate R_i to USSFin through biofilm_split (split_rates.f). Here the
fragments are compiled into the mock usermat (strict 72-column fixed form,
bounds checking), linked against the real bridge hook and C shim and a live
material_server, and driven together with a stand-in for USSFin: an explicit
zero-flux diffusion of bio_i on a line of Gauss points plus dt * R_i -- the
update the partner's USSFin gets, with its local source replaced by R_i.

  (a) a spatially uniform field reproduces the stand-alone point model;
  (b) mass: transport conserves, the reaction adds exactly dt * sum R_i;
  (c) with the reaction off it is transport only and growth is Eq. 36;
  (d) below phi_min there is no reaction (and no server call);
  (e) one server call per point and substep, not per equilibrium iteration.
"""
import os
import sys
import threading
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
sys.path.insert(0, str(_AU / "coupling"))
sys.path.insert(0, str(_AU))
sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip("jax")
import ecology_jax as eco                               # noqa: E402
import material_server as ms                            # noqa: E402
import one_species_reference as ref                     # noqa: E402
from test_partner_callsite_fragment import (            # noqa: E402
    _CC, _FC, MOCK, build_mock, run_mock)

pytestmark = pytest.mark.skipif(_FC is None or _CC is None,
                                reason="gfortran/cc unavailable")

NE, GPPE = 2, 8
NP = NE * GPPE
DT, NSTEPS, K_ALPHA = 0.1, 4, 0.5

_SUB = MOCK[:MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + f"""\
      PROGRAM SPLITDRV
      USE biofilm_split
      IMPLICIT NONE
      INTEGER NP
      PARAMETER (NP = {NP})
      DOUBLE PRECISION prop(29), ust(100,NP), work(100), sg, dt
      DOUBLE PRECISION B1(NP), B2(NP), N1(NP), N2(NP), RR(2,NP)
      DOUBLE PRECISION beta, phimin, b10, b20
      INTEGER mode, nsteps, it, s, i, e, k, kc, ld, is, init, np2
      prop = 0.0D0
      READ(*,*) mode, prop(7), dt, np2, beta, phimin, init, b10, b20,
     &          nsteps
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(28) = DBLE(mode)
      prop(29) = phimin
      kc = 0
      ust = 0.0D0
      DO i = 1, NP
        B1(i) = b10
        B2(i) = 0.0D0
        IF (mode .EQ. 6) B2(i) = b20
        IF (init .EQ. 1 .AND. i .GT. NP / 2) THEN
          B1(i) = b20
          IF (mode .EQ. 6) B2(i) = b10
        END IF
      END DO
      OPEN(20, FILE='field.csv', STATUS='REPLACE')
      DO i = 1, NP
        WRITE(20,'(I0,",",I0,4(",",ES24.16E3))') 0, i, B1(i), B2(i),
     &        0.0D0, 0.0D0
      END DO
      DO s = 1, nsteps
C       usermat: every Gauss point, three iterations, commit the last
        DO i = 1, NP
          e = (i - 1) / 8 + 1
          k = MOD(i - 1, 8) + 1
          DO it = 1, 3
            work = ust(:,i)
            CALL MOCKMAT(np2, prop, work, dt, e, k, 1, s, B1(i),
     &                   1.0D0, B2(i), 1.0D0, 1.0D0, sg, kc)
          END DO
          ust(:,i) = work
        END DO
C       USSFin stand-in: explicit zero-flux diffusion + dt * R_i
        DO i = 1, NP
          CALL SPLIT_GET(i, RR(1,i), RR(2,i), ld, is)
          IF (ld .NE. 1 .OR. is .NE. s) STOP 2
        END DO
        DO i = 1, NP
          N1(i) = B1(i) + dt * (RR(1,i) + beta * LAPL(B1, i))
          N2(i) = B2(i) + dt * (RR(2,i) + beta * LAPL(B2, i))
        END DO
        B1 = N1
        B2 = N2
        DO i = 1, NP
          WRITE(20,'(I0,",",I0,4(",",ES24.16E3))') s, i, B1(i), B2(i),
     &          RR(1,i), RR(2,i)
        END DO
      END DO
      CLOSE(20)
      PRINT '(ES24.16E3)', ust(84,1)
      CONTAINS
      DOUBLE PRECISION FUNCTION LAPL(B, i)
      DOUBLE PRECISION B(NP)
      INTEGER i
      DOUBLE PRECISION L, R
      L = B(i)
      R = B(i)
      IF (i .GT. 1) L = B(i - 1)
      IF (i .LT. NP) R = B(i + 1)
      LAPL = L - 2.0D0 * B(i) + R
      END FUNCTION
      END
"""


def _theta(n):
    th = [0.0] * 20
    th[0], th[3] = 1.34, 0.32                 # a11, b1
    if n == 2:
        th[1], th[2], th[4] = -0.18, 1.79, 1.49   # a12, a22, b2
    return th


@pytest.fixture(scope="module")
def port():
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv.server_address[1]
    srv.shutdown()
    srv.server_close()
    ms.set_active_species(5)


@pytest.fixture(scope="module")
def built():
    return build_mock(DRIVER)


def _run(built, port, n, init, b10, b20, beta, phimin, nsteps=NSTEPS):
    tmp, exe = built
    for f in ("field.csv", "split_trace.csv", "phi_trace.csv"):
        (tmp / f).unlink(missing_ok=True)
    ms.set_active_species(n)
    th = _theta(n)
    line = (f"{4 + n} {K_ALPHA} {DT} 29 {beta} {phimin} {init} "
            f"{b10} {b20} {nsteps}\n" + "".join(f"{x:.17e}\n" for x in th))
    import subprocess
    r = subprocess.run([str(exe)], input=line, capture_output=True,
                       text=True, cwd=tmp, timeout=600,
                       env=dict(os.environ, BIOFILM_PY_PORT=str(port)))
    assert r.returncode == 0, r.stdout + r.stderr
    field = defaultdict(dict)
    for row in (tmp / "field.csv").read_text().splitlines():
        v = [x.strip() for x in row.split(",")]
        field[int(v[0])][int(v[1])] = [float(x) for x in v[2:]]
    trace = []
    p = tmp / "split_trace.csv"
    if p.exists():
        for row in p.read_text().splitlines():
            v = [x.strip() for x in row.split(",")]
            trace.append({"elem": int(v[0]), "ip": int(v[1]),
                          "isubst": int(v[3]), "nsub": int(v[4]),
                          "hit": int(v[5]), "r1": float(v[9])})
    phi = ref.read_trace(tmp / "phi_trace.csv") \
        if (tmp / "phi_trace.csv").exists() else []
    return field, trace, phi, th


@pytest.mark.parametrize("n", [1, 2])
def test_a_uniform_field_reproduces_the_point_model(built, port, n):
    b10, b20 = 0.3, 0.2
    field, _, _, th = _run(built, port, n, 0, b10, b20, 0.5, 0.0)
    g = np.zeros(12)
    g[0] = b10
    if n == 2:
        g[1] = b20
    g[5] = 1.0 - g[:5].sum()
    g[6:6 + n] = 0.999
    worst = 0.0
    for s in range(1, NSTEPS + 1):
        vals = {tuple(field[s][i][:2]) for i in range(1, NP + 1)}
        assert len(vals) == 1, "uniform field did not stay uniform"
        g, _ = eco.ecology_substeps(g, th, DT, int(DT / 1e-4), n)
        g = np.asarray(g)
        worst = max(worst, abs(field[s][1][0] - g[0]),
                    abs(field[s][1][1] - g[1]))
    assert worst < 1e-12, worst
    assert field[NSTEPS][1][0] > b10          # it did react


@pytest.mark.parametrize("n", [1, 2])
def test_b_transport_conserves_and_reaction_adds_dt_R(built, port, n):
    field, _, _, _ = _run(built, port, n, 1, 0.4, 0.05, 0.5, 0.0)
    for s in range(1, NSTEPS + 1):
        for c in range(n):
            m_old = sum(field[s - 1][i][c] for i in range(1, NP + 1))
            m_new = sum(field[s][i][c] for i in range(1, NP + 1))
            added = DT * sum(field[s][i][2 + c] for i in range(1, NP + 1))
            assert abs(m_new - m_old - added) < 1e-13, (s, c)
            assert added != 0.0


def test_c_reaction_off_is_transport_and_eq36(built, port):
    field, trace, phi, _ = _run(built, port, 1, 1, 0.4, 0.05, 0.5, 2.0)
    for s in range(1, NSTEPS + 1):
        assert all(field[s][i][2] == 0.0 for i in range(1, NP + 1))
        m_old = sum(field[s - 1][i][0] for i in range(1, NP + 1))
        m_new = sum(field[s][i][0] for i in range(1, NP + 1))
        assert abs(m_new - m_old) < 1e-14
    assert field[NSTEPS][NP // 2][0] < 0.4    # it diffused
    assert trace and all(r["hit"] == 2 for r in trace)   # no server call
    c = ref.check_trace(phi, K_ALPHA)
    assert c["eq36"] and c["once_per_increment"] and c["carried"]


def test_d_no_reaction_below_phi_min(built, port):
    field, _, _, _ = _run(built, port, 1, 1, 0.4, 0.0, 0.0, 0.01, 1)
    left = [field[1][i][2] for i in range(1, NP // 2 + 1)]
    right = [field[1][i][2] for i in range(NP // 2 + 1, NP + 1)]
    assert all(r != 0.0 for r in left)
    assert all(r == 0.0 for r in right)


def test_e_one_server_call_per_point_and_substep(built, port):
    _, trace, _, _ = _run(built, port, 1, 0, 0.3, 0.0, 0.0, 0.0)
    calls = defaultdict(list)
    for r in trace:
        calls[(r["elem"], r["ip"], r["isubst"])].append(r["hit"])
    assert calls and all(h == [0, 1, 1] for h in calls.values()), calls
    assert all(r["nsub"] == 1000 for r in trace)
