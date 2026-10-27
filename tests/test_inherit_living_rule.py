"""The rule for newly reached points that inherit the living composition (A-82,
coupled_living_consumption.living_start): the amount of the neighbours is kept, it is shared among the species as their
living biomass phi_i psi_i, psi starts fresh (0.999) and gamma at 0; without living biomass the shares of phi_i."""
import ast
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"


def _living_start(n_sp):
    tree = ast.parse(_SRC.read_text(encoding="utf-8"))
    fn = next(f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == "living_start")
    ns = {"np": np, "N_SP": n_sp}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(_SRC), "exec"), ns)
    return ns["living_start"]


def _state(phi, psi, n_sp):
    g = np.zeros(12)
    g[:n_sp] = phi
    g[5] = 1.0 - sum(phi)
    g[6:6 + n_sp] = psi
    g[11] = 0.3
    return g


def test_amount_kept_and_shared_by_living_biomass():
    f = _living_start(2)
    g = _state([0.4, 0.2], [0.1, 0.9], 2)
    out = f(g[None])[0]
    assert np.isclose(out[:2].sum(), 0.6)
    lv = np.array([0.4 * 0.1, 0.2 * 0.9])
    assert np.allclose(out[:2], 0.6 * lv / lv.sum())
    assert np.allclose(out[6:11], 0.999) and out[11] == 0.0


def test_four_species_and_dead_neighbours():
    f = _living_start(4)
    g = _state([0.1, 0.2, 0.3, 0.1], [0.0, 0.0, 0.0, 0.0], 4)
    out = f(g[None])[0]
    assert np.allclose(out[:4], [0.1, 0.2, 0.3, 0.1])     # no living biomass: the shares of phi_i
    g = _state([0.1, 0.2, 0.3, 0.1], [1.0, 0.0, 0.5, 0.0], 4)
    out = f(g[None])[0]
    assert np.isclose(out[:4].sum(), 0.7) and out[1] == 0.0 and out[3] == 0.0
