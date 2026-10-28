"""The homeostatic-pressure growth law as an option of the Python reference of the agreed coupling
(ansys_usermat/apdl/coupled_living_consumption.py, --homeostatic P_H_OVER_E [--gate-form 0|1|2]; Keio continuation, off
by default, not in the thesis), in the style of tests/test_keio_options.py: the gate g(p) of the elements, the
cell <-> stress-mesh mapping, and the fully constrained cube: the explicit gate of the Python side (the hex8 element
matrices of the stress solve, forward Euler in the outer step) against the ODE alpha_dot = k_alpha g(p(alpha)) solved by
scipy (first order) and, when gfortran is there, against the Abaqus element driven without Abaqus
(abaqus_living/homeostatic_driver.f90, the same cube, the same p_h, step and T): the two Euler sequences must agree
to round-off. About 1 s without gfortran, 6 s with it (one compilation)."""
import ast
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_SRC = ROOT / "ansys_usermat" / "apdl" / "coupled_living_consumption.py"
E, NU, FV, K_ALPHA, P_REF = 10.0, 0.49, 1e-3, 1e-3, 0.25          # as tests/test_living_homeostatic.py


def _fn(name):
    tree = ast.parse(_SRC.read_text(encoding="utf-8"))
    fn = next(f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == name)
    ns = {"np": np}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(_SRC), "exec"), ns)
    return ns[name]


gate_fn, prolong = _fn("gate_fn"), _fn("prolong")


def test_gate_forms():
    assert gate_fn(0.0, 1.0, 0.5, 0) == 1.0 and gate_fn(-0.3, 1.0, 0.5, 0) == 1.0     # tension: no gate
    assert gate_fn(0.25, 1.0, 0.5, 0) == 0.5 and gate_fn(0.7, 1.0, 0.5, 0) == 0.0     # clipped at p_h
    assert np.isclose(gate_fn(0.25, 0.5, 1.0, 1), 1.0 - 0.25 / (1.0 * (0.25 + 1e-3)))  # p_h = p_ref (phi^2 + f)
    assert np.isclose(gate_fn(0.5, 1.0, 0.5, 2), 0.5) and gate_fn(-1.0, 1.0, 0.5, 2) == 1.0
    g = gate_fn(np.array([0.0, 0.1, 1.0]), np.ones(3), 0.5, 0)
    assert np.allclose(g, [1.0, 0.8, 0.0])


def test_prolong_is_the_right_inverse_of_block_mean():
    sys.path.insert(0, str(_SRC.parent))
    from wp3_reference import block_mean
    u = np.random.default_rng(0).standard_normal((4, 4, 4))
    v = prolong(u, 12)
    assert v.shape == (12, 12, 12) and np.array_equal(v[3:6, 6:9, 9:12], np.full((3, 3, 3), u[1, 2, 3]))
    assert np.abs(block_mean(v, 4) - u).max() < 1e-15                   # the mean of 27 equal values, to round-off


def _dpdxi():
    """p = -tr(sigma)/3 of the fully constrained cube per unit alpha - 1 from the hex8 matrices of the stress solve
    (mesh_study_seed.element_matrices, the B-bar element of wp3_stress_reference.solve_fields): u = 0, small strain"""
    sys.path.insert(0, str(_SRC.parent))
    import mesh_study_seed as M
    _, _, _, D, m = M.element_matrices(1.0, NU, False)
    return E * (1.0 + FV) * float((D @ m)[:3].mean())


def euler(T, steps, form, dpdxi):
    """the explicit gate of the script on the constrained cube: p of the start of the step, alpha += dt k_alpha phi g"""
    x, dt, out = 0.0, T / steps, []
    for _ in range(steps):
        x = x + dt * K_ALPHA * 1.0 * float(gate_fn(dpdxi * x, 1.0, P_REF, form, FV)); out.append(x)
    return np.array(out)


def test_constrained_pressure_is_the_closed_form():
    assert np.isclose(_dpdxi(), E * (1.0 + FV) / (1.0 - 2.0 * NU), rtol=1e-12)      # 500.5 Pa


def test_constrained_cube_first_order_against_the_ode():
    from scipy.integrate import solve_ivp
    dpdxi = _dpdxi()
    s = solve_ivp(lambda t, x: [K_ALPHA * float(gate_fn(dpdxi * x[0], 1.0, P_REF, 0, FV))], (0.0, 1.1), [0.0],
                  method="DOP853", rtol=1e-12, atol=1e-16)
    xi_ode = s.y[0, -1]
    errs = [abs(euler(1.1, n, 0, dpdxi)[-1] - xi_ode) for n in (22, 44)]
    assert 1.8 < errs[0] / errs[1] < 2.2 and errs[1] < 1e-5 and xi_ode < K_ALPHA * 1.1
    xi = euler(4.0, 40, 0, dpdxi)                                 # growth stops at p_h (free growth would give 4e-3)
    assert xi[-1] == pytest.approx(P_REF / dpdxi, rel=1e-3) and xi[-1] - xi[-2] < 1e-3 * xi[0]


@pytest.fixture(scope="module")
def driver(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    d = tmp_path_factory.mktemp("living_uel")
    src = ROOT / "abaqus_living"
    ff = ["-O2", "-ffree-line-length-none", "-J", str(d), "-I", str(d)]
    subprocess.run(["gfortran", "-O2", "-J", str(d), "-c", str(ROOT / "ansys_usermat/coupling/ecology_native.f"),
                    "-o", str(d / "eco.o")], check=True)
    subprocess.run(["gfortran", "-cpp", *ff, "-c", str(src / "living_uel.f90"), "-o", str(d / "uel.o")], check=True)
    subprocess.run(["gfortran", *ff, "-c", str(src / "mms_src.f90"), "-o", str(d / "mms.o")], check=True)
    subprocess.run(["gfortran", *ff, str(src / "homeostatic_driver.f90"), *(str(d / o) for o in ("uel.o", "mms.o", "eco.o")),
                    "-o", str(d / "hd")], check=True)
    return d / "hd"


@pytest.mark.parametrize("form,steps,T", [(0, 22, 1.1), (0, 44, 1.1), (1, 22, 1.1), (2, 22, 1.1), (0, 40, 4.0)])
def test_python_gate_equals_the_abaqus_element_on_the_constrained_cube(driver, form, steps, T):
    """abaqus_living/homeostatic_driver.f90 (all u = 0, phi = c = 1, no point model, E = 10 Pa, nu = 0.49, f = 1e-3,
    k_alpha = 1e-3, p_h = 0.25 Pa, explicit gate): alpha - 1 and p after every increment against the Python sequence"""
    out = subprocess.run([str(driver), "all", repr(P_REF), str(form), str(steps), repr(T)],
                         capture_output=True, text=True, check=True).stdout
    a = np.array([[float(v) for v in line.split()] for line in out.strip().splitlines()])
    dpdxi = _dpdxi()
    xi = euler(T, steps, form, dpdxi)
    assert np.abs(a[:, 2] - xi).max() < 1e-13 * xi[-1]
    assert np.abs(a[:, 3] - dpdxi * xi).max() < 1e-9 * P_REF           # the element's Newton residual (1e-13)
