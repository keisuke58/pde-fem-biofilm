"""Four species in the composition mode (prop(28) = 7, prop(37) = 4; 5 Oct).

The server runs case 4sp_case1 of Klempt et al. 2026 (Eq. 19-20); the
fragment carries phi_1..phi_4, psi_1..psi_4. Checked:
  - the server's four-species case is the paper reproduction's model
    (klempt2026_reproduction.py, which matches Fig. 7);
  - three points of the mock usermat (growing, void then filled, packed seed)
    equal composition_reference.reference(n = 4) bit for bit;
  - the species-weighted growth law with f_i = prop(37 + i) multiplies the
    alpha increment by sum_i f_i chi_i, and leaves the composition alone.
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
sys.path.insert(0, str(_ROOT / "JAXFEM"))
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip("jax")
import composition_reference as cref                    # noqa: E402
import material_server as ms                            # noqa: E402
from test_partner_callsite_fragment import (            # noqa: E402
    _CC, _FC, MOCK, build_mock)

DT, NSTEPS, K_ALPHA, PHI_MIN, S, CAP = 0.01, 6, 1.0e-3, 0.01, 0.5, 0.9
ELEMS = (1, 38, 75)
F = (1.25, 1.0, 0.75, 0.5)


def phi3(e, s):
    if e == 1:
        return 0.3 + 0.05 * s
    if e == 38:
        return 0.0 if s < 3 else 0.2
    return 1.0


def test_server_case_is_the_paper_model():
    import ecology_jax as eco
    import klempt2026_reproduction as R
    from klempt2026_cases import CASES
    case = CASES["4sp_case1"]
    want = R.simulate(case)
    ms.set_case("4sp_case1")
    try:
        g0, _ = R.build(case)
        g, _ = eco.ecology_substeps(np.asarray(g0), ms.ECOLOGY_CASE["theta"], 1e-4 * case["steps"],
                                    case["steps"], 4, ms.ECOLOGY_CASE["hp"])
    finally:
        ms.set_case(None)
        ms.set_active_species(5)
    assert np.allclose(np.asarray(g)[:4], want["phi"][-1], rtol=1e-10, atol=1e-14)
    assert np.allclose(np.asarray(g)[6:10], want["psi"][-1], rtol=1e-10, atol=1e-14)


_SUB = MOCK[:MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM C4DRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(41), ust(100,3), work(100), sg, dt, b1
      INTEGER it, s, ie, e, kc, nsteps, ielem(3)
      DATA ielem /1, 38, 75/
      prop = 0.0D0
      READ(*,*) prop(7), dt, nsteps, prop(29), prop(31), prop(32),
     &          prop(36), prop(38), prop(39), prop(40), prop(41)
      DO it = 1, 20
        READ(*,*) prop(7 + it)
      END DO
      prop(1) = 1.0D0
      prop(28) = 7.0D0
      prop(37) = 4.0D0
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
            CALL MOCKMAT(41, prop, work, dt, e, 1, 1, s, b1,
     &                   1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
          END DO
          ust(:,ie) = work
        END DO
      END DO
      PRINT '(3ES26.17E3)', ust(84,1), ust(84,2), ust(84,3)
      END
"""


@pytest.fixture(scope="module")
def runs():
    if _FC is None or _CC is None:
        pytest.skip("gfortran/cc unavailable")
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    ms.set_case("4sp_case1")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        tmp, exe = build_mock(DRIVER)
        th = ms.ECOLOGY_CASE["theta"]
        out = {}
        for d in (0.0, 1.0):
            f = F if d else (0.0,) * 4
            line = (f"{K_ALPHA} {DT} {NSTEPS} {PHI_MIN} {S} {CAP} {d} "
                    + " ".join(str(x) for x in f) + "\n"
                    + "".join(f"{x:.17e}\n" for x in th))
            trace = tmp / "comp_trace.csv"
            if trace.exists():
                trace.unlink()
            r = subprocess.run([str(exe)], input=line, capture_output=True,
                               text=True, cwd=tmp, timeout=600,
                               env=dict(os.environ,
                                        BIOFILM_PY_PORT=str(srv.server_address[1])))
            assert r.returncode == 0, r.stdout + r.stderr
            out[d] = (cref.read_comp_trace(trace), [float(x) for x in r.stdout.split()])
        yield out, th, dict(ms.ECOLOGY_CASE["hp"])
    finally:
        srv.shutdown()
        srv.server_close()
        ms.set_case(None)
        ms.set_active_species(5)


@pytest.mark.parametrize("e", ELEMS)
def test_each_point_equals_the_reference(runs, e):
    out, th, hp = runs
    rows = out[0.0][0]
    want = cref.reference([phi3(e, s) for s in range(1, NSTEPS + 1)], th, hp,
                          DT, phi_min=PHI_MIN, s=S, phi_cap=CAP, n=4)
    got = [[r for r in rows if r["elem"] == e and r["isubst"] == s][-1]
           for s in range(1, NSTEPS + 1)]
    for w, g in zip(want, got):
        assert np.array_equal(w, np.asarray(g["g_new"])), (e, g["isubst"])
    assert np.all(np.asarray(got[-1]["g_new"])[:4] > 0)           # all four carried


def test_weighted_growth(runs):
    out, _, _ = runs
    base, wtd = out[0.0], out[1.0]
    for e in ELEMS:
        rb = [r for r in base[0] if r["elem"] == e]
        rw = [r for r in wtd[0] if r["elem"] == e]
        for a, b in zip(rb, rw):
            assert np.array_equal(np.asarray(a["g_new"]), np.asarray(b["g_new"]))
            g = np.asarray(b["g_new"]); sh = g[:4] / g[:4].sum()
            fac = sum(f * x for f, x in zip(F, sh))
            inc_b = a["alpha_new"] - a["alpha_n"]
            inc_w = b["alpha_new"] - b["alpha_n"]
            assert inc_w == pytest.approx(fac * inc_b, rel=1e-12, abs=1e-300)
