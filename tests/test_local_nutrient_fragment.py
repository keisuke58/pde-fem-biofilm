"""Pre-flight for the local nutrient in the composition mode (prop(28) = 7,
prop(33) = c_ref > 0; two-way step 1 of ROADMAP_TWO_WAY.md).

The fragment's nutrient source line (CM_NUT = -1.0D30) is replaced, as
paste_fragments.py --nut-var does, by CM_NUT = prop(34) so the test can
prescribe the nutrient. One packed Gauss point, case 6 of Klempt et al. 2026.
Checked:
  - as shipped (CM_NUT = -1D30) prop(33) has no effect;
  - c_rel = 1 gives the plain composition mode bit for bit;
  - c_rel < 1 reaches the server and changes the composition, by exactly
    what the server gives for that c_rel;
  - paste_fragments.set_nut_source refuses what would not compile.
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
import paste_fragments as pf                            # noqa: E402
import test_partner_callsite_fragment as tp             # noqa: E402

pytestmark = pytest.mark.skipif(tp._FC is None or tp._CC is None,
                                reason="gfortran/cc unavailable")

DT, NSTEPS, S, CAP = 0.01, 3, 0.5, 0.9

_SUB = tp.MOCK[:tp.MOCK.index("      PROGRAM DRIVE")]
DRIVER = _SUB + """\
      PROGRAM NUTDRV
      IMPLICIT NONE
      DOUBLE PRECISION prop(34), work(100), sg, dt
      INTEGER it, s, kc, nsteps
      prop = 0.0D0
      READ(*,*) dt, nsteps, prop(31), prop(32), prop(33), prop(34)
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
        CALL MOCKMAT(34, prop, work, dt, 1, 1, 1, s, 1.0D0,
     &               1.0D0, 0.0D0, 1.0D0, 1.0D0, sg, kc)
      END DO
      PRINT '(12ES26.17E3)', (work(71 + it), it = 1, 12)
      END
"""


def _build(patched):
    exec_lines = pf.read_lines(tp._CS / "phi_mode_exec.inc")
    if patched:
        exec_lines = pf.set_nut_source(exec_lines, "prop(34)")
    src = tp._CS
    tmpcs = Path(__import__("tempfile").mkdtemp())
    for f in ("phi_mode_decl.inc", "split_rates.f"):
        if (src / f).exists():
            (tmpcs / f).write_bytes((src / f).read_bytes())
    (tmpcs / "phi_mode_exec.inc").write_text("\n".join(exec_lines) + "\n")
    old = tp._CS
    tp._CS = tmpcs
    try:
        return tp.build_mock(DRIVER)
    finally:
        tp._CS = old


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
def builds():
    return {"shipped": _build(False), "patched": _build(True)}


def _run(server, build, c_ref, nut):
    tmp, exe = build
    th = ms.ECOLOGY_CASE["theta"]
    line = (f"{DT} {NSTEPS} {S} {CAP} {c_ref} {nut}\n"
            + "".join(f"{x:.17e}\n" for x in th))
    r = subprocess.run([str(exe)], input=line, capture_output=True, text=True,
                       cwd=tmp, timeout=600,
                       env=dict(os.environ,
                                BIOFILM_PY_PORT=str(server.server_address[1])))
    assert r.returncode == 0, r.stdout + r.stderr
    return np.array([float(x) for x in r.stdout.split()])


def test_as_shipped_prop33_has_no_effect(server, builds):
    base = _run(server, builds["shipped"], 0.0, 0.0)
    assert np.array_equal(_run(server, builds["shipped"], 1.0, 0.3), base)


def test_full_nutrient_equals_the_plain_mode(server, builds):
    base = _run(server, builds["patched"], 0.0, 0.0)
    assert np.array_equal(_run(server, builds["patched"], 1.0, 1.0), base)
    assert np.array_equal(_run(server, builds["patched"], 2.0, 3.0), base)  # clamped


def test_low_nutrient_reaches_the_server(server, builds):
    base = _run(server, builds["patched"], 0.0, 0.0)
    low = _run(server, builds["patched"], 1.0, 0.3)
    assert not np.array_equal(low, base)
    share = lambda g: g[0] / (g[0] + g[1])
    assert abs(share(low) - share(base)) > 1e-6
    # same as c_ref = 2, nut = 0.6 (only the ratio enters)
    assert np.array_equal(_run(server, builds["patched"], 2.0, 0.6), low)


def test_negative_nutrient_is_starved_not_full(server, builds):
    # zero-order consumption lets the field's c go below 0 (5 Oct, consumption > 6);
    # such a point must see c_rel = 0, not fall back to the plain mode (c_rel = 1)
    base = _run(server, builds["patched"], 0.0, 0.0)
    zero = _run(server, builds["patched"], 1.0, 0.0)
    assert not np.array_equal(zero, base)
    assert np.array_equal(_run(server, builds["patched"], 1.0, -0.3), zero)


def test_set_nut_source_refuses_bad_names():
    lines = pf.read_lines(tp._CS / "phi_mode_exec.inc")
    assert sum(l == pf.NUT_LINE for l in lines) == 1
    out = pf.set_nut_source(lines, "vGdp_Nut1_n(ID)")
    assert "          CM_NUT = vGdp_Nut1_n(ID)" in out
    for bad in ("1abc", "x; y", "a" * 70):
        with pytest.raises(SystemExit):
            pf.set_nut_source(lines, bad)
