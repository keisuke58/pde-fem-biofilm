"""Code verification for the PDE discretisations, against exact solutions.

These pin `JAXFEM/pde_code_verification.py`. See `PDE_VERIFICATION_FINDINGS.md`
for what it found and `JAXFEM/boundary_fix_impact.py` for what the fix did to
the reported 2D condition spread.

Every target here is *exact* -- a cosh profile, a parabola, a Neumann
eigenfunction, a conservation law -- so a failure is a defect in the
discretisation and not a disagreement with another implementation.

Both sides are pinned deliberately: that the committed wall is second order and
conservative, and that the pre-2026-09-30 wall was not. Keeping the second
means the measurement that motivated the change cannot quietly stop being
reproducible.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("jax")

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "JAXFEM"))

import jax.numpy as jnp                                            # noqa: E402

from JAXFEM import pde_code_verification as v                      # noqa: E402
from JAXFEM import core_hamilton_1d_nutrient as p1                 # noqa: E402
from JAXFEM import core_hamilton_2d_nutrient as p2                 # noqa: E402

GRIDS = (21, 41, 81)


def _order(rows):
    return [np.log(a[2] / b[2]) / np.log(a[1] / b[1])
            for a, b in zip(rows[:-1], rows[1:])]


# ---------------------------------------------------------------------------
# the variants really are the committed routines with one coefficient changed
# ---------------------------------------------------------------------------

def test_the_1d_mirror_variant_is_the_committed_step_bit_for_bit():
    """Without this the comparisons below would be two opinions, not one
    change."""
    rng = np.random.default_rng(0)
    for N in (21, 41):
        params, _, _ = v._params(N, 1.0, 4.0e6, 1.0e6)
        for _ in range(5):
            c = jnp.asarray(rng.uniform(0.1, 1.0, N))
            a = np.asarray(p1.nutrient_step(c, 1.0, params))
            b = np.asarray(v.nutrient_step_mirror(c, 1.0, params))
            assert np.array_equal(a, b)


def test_the_2d_mirror_variant_is_the_committed_operator_bit_for_bit():
    rng = np.random.default_rng(1)
    for N in (11, 21):
        u = jnp.asarray(rng.uniform(0.0, 1.0, (N, N)))
        h = 1.0 / (N - 1)
        a = np.asarray(p2.laplacian_2d_neumann(u, h, h))
        b = np.asarray(v.laplacian_2d_neumann_mirror(u, h, h))
        assert np.array_equal(a, b)


# ---------------------------------------------------------------------------
# 1D: second order now, first order before
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def steady_1d():
    return {name: v.study(step, grids=GRIDS)
            for name, step in (("committed", None),
                               ("pre-fix", v.nutrient_step_half))}


@pytest.mark.parametrize("case", list(v.CASES))
def test_the_committed_solver_is_second_order(steady_1d, case):
    """Measured against cosh (first-order kinetics) and the parabola that is
    Klempt 2024 Eq. 35's own steady state."""
    rows = steady_1d["committed"][case]
    if case.startswith("zeroth"):
        # a second-order central difference is EXACT on a quadratic, so there
        # is no order to measure -- the error sits at the asymptotic floor.
        assert max(r[2] for r in rows) < 1.0e-7
    else:
        for p in _order(rows)[-1:]:
            assert 1.8 < p < 2.2, f"{case}: observed order {p}"


@pytest.mark.parametrize("case", list(v.CASES))
def test_the_pre_fix_wall_was_only_first_order(steady_1d, case):
    for p in _order(steady_1d["pre-fix"][case])[-1:]:
        assert 0.8 < p < 1.3, f"{case}: observed order {p}"


def test_the_pre_fix_error_was_worst_at_the_neumann_wall(steady_1d):
    """Localises the defect that was fixed: the sup norm was attained at x=0
    on every grid in every case."""
    for case, rows in steady_1d["pre-fix"].items():
        for N, h, sup, at_wall, interior, _res in rows:
            assert at_wall == pytest.approx(sup), (case, N)
            assert interior < sup


def test_the_committed_solver_beats_the_pre_fix_one_on_every_grid(steady_1d):
    for case in v.CASES:
        for new, old in zip(steady_1d["committed"][case],
                            steady_1d["pre-fix"][case]):
            assert new[2] < old[2], (case, new[0], new[2], old[2])


def test_the_zeroth_order_floor_is_the_asymptotic_not_the_discretisation():
    """err / k_monod is flat, so the discretisation contributes nothing: the
    committed solver reproduces Klempt Eq. 35's steady state exactly, to the
    limit of the c >> k_monod approximation."""
    ratios = []
    for k in (1.0e-6, 1.0e-8, 1.0e-10):
        x, c, _ = v.steady_state(41, 1.0, 1.0, k, 1.0)
        ce = v.exact_zeroth_order(x, 1.0, 1.0, 1.0)
        ratios.append(float(np.max(np.abs(c - ce))) / k)
    assert max(ratios) / min(ratios) < 1.1, ratios


# ---------------------------------------------------------------------------
# 2D: the pre-fix operator did not converge at the walls at all
# ---------------------------------------------------------------------------

def test_the_committed_2d_laplacian_is_second_order():
    rows = v.study_2d_operator(p2.laplacian_2d_neumann, grids=GRIDS)
    for p in _order(rows)[-1:]:
        assert 1.8 < p < 2.2, p


def test_the_pre_fix_2d_laplacian_did_not_converge_at_the_walls():
    """The wall rows returned half the correct second difference, so the error
    was half the exact eigenvalue 2(pi/L)^2 = 19.74 and h-independent."""
    rows = v.study_2d_operator(v.laplacian_2d_neumann_half, grids=GRIDS)
    errs = [r[2] for r in rows]
    assert all(abs(e - 0.5 * 2.0 * np.pi ** 2) < 0.05 for e in errs), errs
    assert max(errs) / min(errs) < 1.01           # flat under refinement
    assert all(r[3] < r[2] / 100.0 for r in rows)  # the interior was fine


def test_the_committed_2d_diffusion_conserves_mass_through_a_zero_flux_wall():
    """The consequence in the units the 2D study reports.
    `diffusion_step_species_2d` (core_hamilton_2d_nutrient.py:426) calls this
    operator at line 434, and `condition_spread_2d.py:173` calls that -- so
    this is on the path the condition-separation number is computed on."""
    N, D = 41, 1.0
    h = 1.0 / (N - 1)
    dt = 0.2 * h * h / D
    g = np.linspace(0.0, 1.0, N)
    X, Y = np.meshgrid(g, g, indexing="ij")
    u0 = np.exp(-((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 0.01)
    w = np.ones((N, N))                     # finite-volume (trapezoidal) measure
    w[0, :] *= 0.5
    w[-1, :] *= 0.5
    w[:, 0] *= 0.5
    w[:, -1] *= 0.5

    def drift(op):
        u = jnp.asarray(u0)
        m0 = float(np.sum(np.asarray(u) * w))
        for _ in range(2000):
            u = u + dt * D * op(u, h, h)
        return (float(np.sum(np.asarray(u) * w)) - m0) / m0

    assert abs(drift(p2.laplacian_2d_neumann)) < 1.0e-12
    assert drift(v.laplacian_2d_neumann_half) < -0.01   # lost several percent


# ---------------------------------------------------------------------------
# the AT2 damage solver -- the counter-example: it was already right
# ---------------------------------------------------------------------------

def test_the_at2_damage_solver_is_second_order():
    """By manufactured solution, `d = 0.5 + 0.2 cos(n pi z/L)`, which satisfies
    both Neumann conditions exactly. A uniform driving field would make d''
    vanish and prove nothing -- the same blindness that let the nutrient
    defect survive test_pde_uniform_consistency."""
    rows = v.study_at2(grids=(21, 41, 81), n_mode=2)
    assert all(r[3] > 0.0 for r in rows), "manufactured H must stay positive"
    for p in _order(rows)[-1:]:
        assert 1.9 < p < 2.1, p


def test_the_at2_solver_reflects_about_the_node():
    """The reason it passes: its Neumann rows use ghost d[-1] = d[1] AND double
    the coefficient of d[1]. The same repository held both the right and the
    wrong version of the same boundary, in different files."""
    rows = v.study_at2(grids=(41, 81), n_mode=4)
    assert rows[-1][2] < rows[0][2] / 3.0, rows     # ~4x per halving


# ---------------------------------------------------------------------------
# the Allen-Cahn interface against its exact tanh profile
# ---------------------------------------------------------------------------

def test_the_allen_cahn_interface_matches_its_exact_tanh_profile():
    """phi = (1 + tanh(x/(2 xi)))/2 with xi = sqrt(beta/Gamma) is the exact
    stationary solution of Klempt Eq. 34's double well against its Laplacian,
    and that xi is the interface width felix_complete_reproduction.py
    documents. The residual of the module's own right-hand side on it
    converges at second order."""
    rows = v.study_tanh_interface(grids=(41, 81, 161))
    for p in _order(rows)[-1:]:
        assert 1.9 < p < 2.1, p


# ---------------------------------------------------------------------------
# census of every zero-flux wall in the repository
# ---------------------------------------------------------------------------

def test_the_reported_klempt_numbers_do_not_inherit_the_boundary_defect():
    """The one that matters most: `klempt2024_quantitative.lap` is where the
    reported Klempt 2024 figures come from. It uses numpy's `reflect` padding,
    which puts the ghost at node 1 -- the correct mirror -- so it agrees with
    a hand-built mirrored 7-point Laplacian exactly."""
    assert v.klempt2024_lap_matches_mirror() == 0.0


def test_the_census_of_zero_flux_walls_is_what_it_is_recorded_to_be():
    """Pins which operators are right and which are not, so a new instance of
    the `mode="edge"` idiom (or a silent fix of one of these) shows up here.
    The whole defect is one word: numpy `pad(mode="edge")` puts the ghost at
    the boundary node, `mode="reflect"` puts it at node 1."""
    ops = v.census_operators()
    for name in v.CENSUS_CORRECT:
        if name not in ops:
            continue
        rows = v.study_2d_operator(ops[name], grids=(21, 41, 81))
        assert min(_order(rows)) > 1.8, (name, rows)
    for name in v.CENSUS_DEFECTIVE:
        rows = v.study_2d_operator(ops[name], grids=(21, 41, 81))
        errs = [r[2] for r in rows]
        assert max(errs) / min(errs) < 1.01, (name, errs)   # flat: no convergence
        assert all(r[3] < r[2] / 100.0 for r in rows), (name, rows)  # interior fine
