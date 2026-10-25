"""Keio-continuation options of coupled_living_consumption.py (off by default, not in the thesis):
--front-living (front speed times the living volume fraction of the cell and its neighbours) and
--stiff-species (stiffness factor from the species shares)."""
import ast
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"


def _fn(name):
    tree = ast.parse(_SRC.read_text(encoding="utf-8"))
    fn = next(f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == name)
    ns = {"np": np}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(_SRC), "exec"), ns)
    return ns[name]


def test_living_speed_takes_the_neighbours_maximum():
    f = _fn("living_speed")
    L = np.zeros((5, 5, 5))
    L[2, 2, 2] = 0.8
    out = f(L)
    assert out[2, 2, 2] == 0.8
    for d in [(1, 2, 2), (3, 2, 2), (2, 1, 2), (2, 3, 2), (2, 2, 1), (2, 2, 3)]:
        assert out[d] == 0.8
    assert out[1, 1, 2] == 0.0 and out[0, 0, 0] == 0.0


def test_living_speed_is_one_for_fully_living():
    f = _fn("living_speed")
    assert np.allclose(f(np.ones((4, 4, 4))), 1.0)


def test_stiffness_factor_equal_moduli_is_one():
    f = _fn("stiffness_factor")
    G = np.zeros((3, 12))
    G[0, :2] = (0.3, 0.2)
    G[1, :2] = (0.0, 0.5)
    assert np.allclose(f(G, [1.0, 1.0]), 1.0)


def test_stiffness_factor_weights_by_share():
    f = _fn("stiffness_factor")
    G = np.zeros((2, 12))
    G[0, :2] = (0.3, 0.1)          # shares 0.75, 0.25
    out = f(G, [2.0, 1.0])
    assert np.isclose(out[0], 1.75)
    assert out[1] == 1.0           # point model not started: factor 1
