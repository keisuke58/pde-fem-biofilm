"""The composition coupling reduces to the point model of Klempt et al. 2026.

In the composition mode the 3D field decides the amount and the point model
only the composition: its state is scaled to the field's amount before and
after every step. If the field's amount is the one the point model would
produce on its own, the scaling must change nothing and the coupled scheme
must give back the paper's composition. This is the property that answers
"how is the point model connected to the field, and is that justified": the
coupled result can differ from the paper only where the field's amount
differs from the point model's.

Measured (2026-10-02): case 3 agrees to 1.2e-5, case 6 to 6e-17, at coupling
steps 0.01 down to 0.0025.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

_AU = Path(__file__).resolve().parents[1] / "ansys_usermat"
sys.path.insert(0, str(_AU / "coupling"))
sys.path.insert(0, str(_AU))
pytest.importorskip("jax")
import composition_reference as cref                    # noqa: E402
import ecology_jax as eco                               # noqa: E402
import material_server as ms                            # noqa: E402

DT, NSTEPS, NSUB = 0.01, 15, 100          # 1500 steps of 1e-4, T = 0.15


@pytest.mark.parametrize("case, tol", [("2sp_case3", 5e-5),
                                       ("2sp_case6", 1e-12)])
def test_same_amount_gives_the_papers_composition(case, tol):
    ms.set_case(case)
    try:
        th, hp = ms.ECOLOGY_CASE["theta"], ms.ECOLOGY_CASE["hp"]
        # the point model on its own, from the paper's start 0.2 / 0.2
        g = cref.rescale(cref.seed(0.5), 0.4)
        amount, chi_alone = [], []
        for _ in range(NSTEPS):
            amount.append(g[0] + g[1])                  # at the step's start
            g, _ = eco.ecology_substeps(g, th, DT, NSUB, 2, hp)
            g = np.asarray(g)
            chi_alone.append(g[0] / (g[0] + g[1]))
        # the coupled scheme, driven by that same amount
        coupled = cref.reference(amount, th, hp, DT, chi1=0.5)
        chi_coupled = [r[0] / (r[0] + r[1]) for r in coupled]
        worst = max(abs(a - b) for a, b in zip(chi_alone, chi_coupled))
        assert worst < tol, worst
        # and the amount is the field's, not the point model's
        assert all(abs(r[0] + r[1] - a) < 1e-15
                   for r, a in zip(coupled, amount))
    finally:
        ms.set_case(None)
        ms.set_active_species(5)
