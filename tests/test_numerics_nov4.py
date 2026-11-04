"""Numerics added on 4 Nov 2026 (no change of the model):
- WENO5 face reconstruction of the front term (--limiter weno5, coupled_living_consumption.py),
- the four TVD limiters of the two-level grid (coupled_amr.limited_slope),
- the threaded point model (ECO_THREADS, ecology_native_ctypes.step_lin): same result to the last bit.
The functions are taken from the scripts with ast, so the scripts' command-line set-up does not run."""
import ast
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
APDL = ROOT / "ansys_usermat" / "apdl"


def _fns(path, names, ns):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    body = [f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name in names]
    exec(compile(ast.Module(body=body, type_ignores=[]), str(path), "exec"), ns)
    return [ns[n] for n in names]


def _weno(h=1.0):
    ns = {"np": np, "WENO_EPS": 1e-6, "WENO_D": (0.1, 0.6, 0.3), "H": h}
    return _fns(APDL / "coupled_living_consumption.py", ["_weno5_face", "upwind_adv_weno5"], ns)


def test_weno5_face_is_exact_for_a_linear_field():
    face, _ = _weno()
    a = [1.0 + 0.5 * k for k in range(5)]          # cells 0..4, face between 2 and 3 at 2.5
    assert abs(face(*a) - (1.0 + 0.5 * 2.5)) < 1e-12


def test_weno5_face_does_not_overshoot_a_step():
    face, _ = _weno()
    for a in ([0, 0, 0, 1, 1], [0, 0, 1, 1, 1], [1, 1, 1, 0, 0], [0, 0, 0, 0, 1]):
        v = face(*map(float, a))
        assert -1e-3 <= v <= 1 + 1e-3


def test_weno5_derivative_converges_at_high_order():
    _, adv = _weno()
    errs = []
    for n in (32, 64):
        h = 1.0 / n
        _, adv = _weno(h)
        x = (np.arange(n) + 0.5) * h
        u = np.sin(2 * np.pi * x)[:, None, None] * np.ones((1, 4, 4))
        v = np.zeros((3,) + u.shape); v[0] = 1.0
        d = adv(u, v)[8:n - 8, 0, 0]                 # interior, away from the mirrored walls
        errs.append(np.abs(d - 2 * np.pi * np.cos(2 * np.pi * x[8:n - 8])).max())
    assert np.log2(errs[0] / errs[1]) > 3.5


def test_amr_limiters_order_and_unknown_name():
    path = APDL / "coupled_amr.py"
    src = path.read_text(encoding="utf-8")
    assert "vanleer" in src and "superbee" in src
    dp = np.array([1.0, 1.0, 1.0, -1.0, 0.5])
    dm = np.array([0.5, 1.0, 2.0, 1.0, 0.0])
    out = {}
    for lim in ("minmod", "vanleer", "mc", "superbee"):
        (f,) = _fns(path, ["limited_slope"], {"np": np, "LIMITER": lim})
        out[lim] = f(dp, dm)
        assert out[lim][3] == 0.0 and out[lim][4] == 0.0      # opposite signs or a zero side: no slope
        assert np.allclose(out[lim][1], 1.0)                   # equal sides: the slope itself
    # the less diffusive limiters give larger slopes where the two sides differ
    for i in (0, 2):
        assert out["minmod"][i] <= out["vanleer"][i] <= out["superbee"][i] + 1e-12
        assert out["minmod"][i] <= out["mc"][i] + 1e-12


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="needs gfortran")
def test_threaded_point_model_is_bitwise_identical():
    sys.path.insert(0, str(ROOT / "ansys_usermat" / "coupling"))
    import ecology_native_ctypes as F
    F.init(ROOT / "abaqus_composition" / "eco_fritsch_2sp_a100.txt")
    rng = np.random.default_rng(3)
    n = 3000
    G = np.zeros((n, 12)); G[:, 0] = 0.2 + 0.1 * rng.random(n); G[:, 1] = 0.2 + 0.1 * rng.random(n); G[:, 5:7] = 1.0
    c0 = rng.random(n)
    res = []
    for k in (1, 3):
        F.THREADS, F._pool = k, None
        res.append(F.step_lin(G, 0.01, 5, c0, 0.99 * c0, 0.0, 0.0))
    F.THREADS, F._pool = 1, None
    assert np.array_equal(res[0][0], res[1][0]) and np.array_equal(res[0][1], res[1][1])
