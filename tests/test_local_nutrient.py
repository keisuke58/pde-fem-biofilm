"""Local nutrient in the point model (ROADMAP_TWO_WAY.md step 1, 2026-10-04).

An optional c_rel scales the nutrient level of one point-model call,
c* = c*_0 * c_rel. Pinned here:
  - without c_rel, and with c_rel = 1, the result is the old one bit for bit
    (with and without a Klempt et al. 2026 case loaded);
  - c_rel = 0.5 equals ecology_substeps with c* halved;
  - negative c_rel is refused by the server;
  - the C shim's biofilm_ecology_eval_c sends c_rel, and a negative value
    sends nothing (equal to biofilm_ecology_eval).
"""
import os
import shutil
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import json
import numpy as np
import pytest

pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
_COUP = _ROOT / "ansys_usermat" / "coupling"
sys.path.insert(0, str(_COUP))

import ecology_jax                                 # noqa: E402
import material_server as ms                       # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO  # noqa: E402

G0 = [0.12, 0.10, 0.05, 0.02, 0.0, 0.71, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0]
DT = 1.0e-3
NSUB = 10


def _call(c_rel=None, theta=None):
    req = {"kind": "ecology", "g": G0, "theta": list(theta if theta is not None else THETA_DEMO.tolist()),
           "dt_h": DT, "n_sub": NSUB}
    if c_rel is not None:
        req["c_rel"] = c_rel
    return json.loads(ms.evaluate_ecology(req))


@pytest.fixture
def no_case():
    ms.set_case(None)
    n0 = ms.ECOLOGY_ACTIVE
    ms.set_active_species(5)
    yield
    ms.set_active_species(n0)


def test_c_rel_one_equals_absent(no_case):
    assert _call() == _call(1.0)


def test_c_rel_one_equals_absent_with_case():
    ms.set_case("2sp_case6")
    try:
        th = ms.ECOLOGY_CASE["theta"]
        g0 = [0.3, 0.3, 0, 0, 0, 0.4, 0.999, 0.999, 0, 0, 0, 0]
        base = {"kind": "ecology", "g": g0, "theta": th, "dt_h": DT, "n_sub": NSUB}
        a = json.loads(ms.evaluate_ecology(dict(base)))
        b = json.loads(ms.evaluate_ecology(dict(base, c_rel=1.0)))
        assert a == b
        c = json.loads(ms.evaluate_ecology(dict(base, c_rel=0.5)))
        hp = dict(ms.ECOLOGY_CASE["hp"]); hp["c"] *= 0.5
        g, pint = ecology_jax.ecology_substeps(g0, th, DT, NSUB, 2, hp)
        assert np.array_equal(np.asarray(c["g_new"]), np.asarray(g))
        assert c["g_new"] != a["g_new"]
    finally:
        ms.set_case(None)
        ms.set_active_species(5)


def test_c_rel_half_scales_c_star(no_case):
    r = _call(0.5)
    hp = {"c": ecology_jax.C_STAR * 0.5, "alpha": ecology_jax.ALPHA_STAR, "eta": [1.0] * 5}
    g, pint = ecology_jax.ecology_substeps(G0, THETA_DEMO.tolist(), DT, NSUB, 5, hp)
    assert np.array_equal(np.asarray(r["g_new"]), np.asarray(g))
    assert r["g_new"] != _call()["g_new"]


def test_negative_c_rel_refused(no_case):
    with pytest.raises(ValueError):
        _call(-0.1)


_CC = shutil.which("cc") or shutil.which("gcc")
_DRIVER = r"""
#include <stdio.h>
int biofilm_ecology_eval_c(const double*, const double*, double, int, double, double*, double*);
void biofilm_py_close(void);
int main(int argc, char **argv) {
    double g[12], th[20], dt, c, gn[12], pi; int i, rc;
    for (i = 0; i < 12; i++) if (scanf("%lf", &g[i]) != 1) return 100;
    for (i = 0; i < 20; i++) if (scanf("%lf", &th[i]) != 1) return 100;
    if (scanf("%lf %lf", &dt, &c) != 2) return 100;
    rc = biofilm_ecology_eval_c(g, th, dt, 10, c, gn, &pi);
    if (rc) { biofilm_py_close(); return rc; }
    for (i = 0; i < 12; i++) printf("%.17g\n", gn[i]);
    biofilm_py_close(); return 0;
}
"""


@pytest.mark.skipif(_CC is None, reason="no C compiler")
def test_c_shim_sends_c_rel(no_case):
    tmp = Path(tempfile.mkdtemp())
    (tmp / "drv.c").write_text(_DRIVER)
    exe = tmp / "drv"
    ws2 = ["-lws2_32"] if sys.platform == "win32" else []
    subprocess.run([_CC, str(tmp / "drv.c"), str(_COUP / "biofilm_py_eval.c"), "-o", str(exe)] + ws2,
                   check=True)
    try:
        srv = ms._Server(("127.0.0.1", 0), ms._Handler)
    except OSError:
        pytest.skip("cannot bind a local socket in this environment")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    host, port = srv.server_address
    env = dict(os.environ, BIOFILM_PY_HOST=host, BIOFILM_PY_PORT=str(port))
    try:
        def run(c):
            stdin = " ".join(f"{v:.17g}" for v in (*G0, *THETA_DEMO.tolist(), DT, c)) + "\n"
            r = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, env=env, timeout=60)
            assert r.returncode == 0, r.stderr
            return [float(x) for x in r.stdout.split()]
        assert run(0.5) == _call(0.5)["g_new"]
        assert run(-1.0) == _call()["g_new"]
    finally:
        srv.shutdown(); srv.server_close()
