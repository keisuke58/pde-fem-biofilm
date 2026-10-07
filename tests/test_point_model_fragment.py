"""Pre-flight for prop(28) = 4, the point-model bridge, before IKMHIWI03.

The call-site fragments are compiled into the same mock usermat as in
test_partner_callsite_fragment.py, linked against the real bridge hook and C
shim, and pointed at a real material_server run with one active species.
With dt = 0.1 the hook asks for 1000 inner steps per increment. The trace
it writes is judged by the same checker the ANSYS run will be judged by:
growth once per increment, inner state carried between increments, the
inner step count, the inactive species exactly zero, and a replay of every
increment with ecology_substeps.
"""
import os
import sys
import threading
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_AU = _ROOT / "ansys_usermat"
sys.path.insert(0, str(_AU / "coupling"))
sys.path.insert(0, str(_AU))
sys.path.insert(0, str(Path(__file__).resolve().parent))
pytest.importorskip("jax")
import material_server as ms                            # noqa: E402
import one_species_reference as ref                     # noqa: E402
from test_partner_callsite_fragment import (            # noqa: E402
    _CC, _FC, build_mock, run_mock)
from hamilton_pde_jaxfem import THETA_DEMO              # noqa: E402

pytestmark = pytest.mark.skipif(_FC is None or _CC is None,
                                reason="gfortran/cc unavailable")


def _theta_one_species():
    """Only a11 and b1 non-zero: the n = 1 prop(8:27) layout
    (theta_to_matrices: a11 = theta[0], b1 = theta[3])."""
    th = [0.0] * 20
    th[0], th[3] = float(THETA_DEMO[0]), float(THETA_DEMO[3])
    return th


@pytest.fixture(scope="module")
def server():
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    ms.set_active_species(1)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv.server_address[1]
    srv.shutdown()
    srv.server_close()
    ms.set_active_species(5)


@pytest.fixture(scope="module")
def run(server):
    tmp, exe = build_mock()
    th = _theta_one_species()
    env = dict(os.environ, BIOFILM_PY_PORT=str(server))
    r = run_mock(tmp, exe, 4, 0.5, 0.1, theta=th, env=env, timeout=300)
    assert r.returncode == 0, r.stderr
    rows = ref.read_pm_trace(tmp / "pm_trace.csv")
    return rows, [float(x) for x in r.stdout.split()], th


def test_every_sampled_point_gets_one_row_per_iteration(run):
    rows, _, _ = run
    assert {r["elem"] for r in rows} == {1, 38}  # stride 37: 2 not sampled
    assert len(rows) == 2 * 5 * 3                # elements x incr x iter
    assert all(r["nsub"] == 1000 for r in rows)


def test_point_model_trace_passes_every_check(run):
    rows, _, th = run
    c = ref.check_pm_trace(rows, 0.5, 1, theta=th)
    for key in ("once_per_increment", "carried", "nsub", "inactive_zero"):
        assert c[key], key
    assert c["replay"], c
    assert c["replay_worst_g"] == 0.0, c


def test_alpha_grows_and_matches_the_stored_value(run):
    rows, alphas, _ = run
    last = [r for r in rows if r["elem"] == 1][-1]
    assert last["alpha_new"] > 0.0
    assert alphas[0] == last["alpha_new"]
    # every element starts from the same seed and sees the same theta, so
    # the point model ignores the mock's per-element bio values
    assert alphas[0] == alphas[1] == alphas[2]

