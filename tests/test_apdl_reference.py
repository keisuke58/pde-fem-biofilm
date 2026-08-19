"""Guard the closed-form reference values used by the ANSYS growth deck.

`ansys_usermat/apdl/reference_values.json` is what an ANSYS run is compared
against. If the constitutive core ever changes, the committed reference must
change with it in the same commit — otherwise the ANSYS check silently starts
validating against a stale target.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_APDL = _ROOT / "ansys_usermat" / "apdl"
_REF = _APDL / "reference_values.json"

sys.path.insert(0, str(_ROOT / "ansys_usermat" / "coupling"))
from material_server import stress_core  # noqa: E402


@pytest.fixture(scope="module")
def ref():
    return json.loads(_REF.read_text())


def test_reference_matches_the_constitutive_core(ref):
    """Regenerating must be a no-op — the committed file is not hand-edited."""
    rc = subprocess.run(
        [sys.executable, str(_APDL / "make_reference.py"), "--check"],
        capture_output=True, text=True,
    )
    assert rc.returncode == 0, rc.stderr


def test_constrained_growth_is_pure_hydrostatic_compression(ref):
    """F = I under isotropic growth can only produce uniform pressure.

    A nonzero shear component here would mean the ANSYS Voigt map (VI/VJ) is
    mis-wired — the exact porting trap this deck exists to catch.
    """
    for label, case in ref["cases"].items():
        s = case["stress"]
        assert s[0] == pytest.approx(s[1]) == pytest.approx(s[2]), label
        assert s[3] == s[4] == s[5] == 0.0, f"{label}: growth produced shear"
        assert s[0] < 0, f"{label}: constrained growth must compress, not stretch"


def test_growth_compresses_more_as_alpha_rises(ref):
    elastic = {c["alpha"]: c["stress"][0] for c in ref["cases"].values() if c["eta"] == 0}
    assert elastic[0.20] < elastic[0.05] < 0


def test_viscosity_relaxes_the_growth_stress(ref):
    """A finite eta lets Fv absorb part of the growth, so |sigma| drops."""
    by = {(c["alpha"], c["eta"]): c["stress"][0] for c in ref["cases"].values()}
    for alpha in (0.05, 0.20):
        assert abs(by[(alpha, 8e-3)]) < abs(by[(alpha, 0.0)]), alpha


def test_reference_reproduces_from_F_identity(ref):
    """Recompute independently of make_reference.py's own bookkeeping."""
    p = ref["properties"]
    I3 = np.eye(3)
    for label, case in ref["cases"].items():
        sv, _, je = stress_core(I3, I3, case["alpha"], p["C10"], p["C01"],
                                p["D1"], case["eta"], p["mtype"], p["dt"])
        assert je == pytest.approx(case["Je"], rel=1e-12), label
        np.testing.assert_allclose(sv, case["stress"], rtol=1e-12, err_msg=label)


def test_deck_time_matches_the_viscous_reference_dt():
    """The viscous rows are only reproducible if the deck's TIME equals dt."""
    deck = (_APDL / "t_growth_constrained.dat").read_text()
    dt = json.loads(_REF.read_text())["properties"]["dt"]
    assert f"TIME,{dt}" in deck, f"deck TIME must be {dt} to match the reference"
