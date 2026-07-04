"""Golden-value regression guard for the headline stress results (#3).

Freezes the thesis's headline numbers so a pipeline change that silently moves a
result trips CI instead of slipping into the paper. If a change is intended,
update the frozen constants in the same commit — that makes the result change
explicit and reviewable.

Data-provenance note (see DS_composition_fix.md): a dysbiotic-static (DS)
composition copy-bug once made sigma_DS ~8x too high (13.63 kPa). It was
corrected to ~1.55 kPa (tooth) / ~1.06 kPa (implant). The corrected values live
in the stress-CI artifacts (JAXFEM/_posterior_ci/klempt_stress_ci_*.json), which
this test treats as the source of truth. `tooth_klempt_comparison.json`'s
`_flat_golden` block was NOT regenerated after that fix, so its DS entries are
still the stale ~8x values — tracked below by `test_comparison_json_DS_is_stale`
so the issue self-clears when the artifact is regenerated. The headline
sigma_CH/sigma_DH = 6.44x is DS-independent and unaffected.
"""
import json
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_CMP = _ROOT / "tooth_klempt_comparison.json"
_CI_TOOTH = _ROOT / "JAXFEM" / "_posterior_ci" / "klempt_stress_ci_tooth.json"
_CI_IMPLANT = _ROOT / "JAXFEM" / "_posterior_ci" / "klempt_stress_ci_implant.json"

CONDS = ("commensal_hobic", "commensal_static", "dysbiotic_hobic", "dysbiotic_static")

# --- corrected MAP von Mises stress (stress-CI artifacts, post DS-fix) ---------
GOLDEN_SIGMA_TOOTH = {
    "commensal_hobic": 0.01372, "commensal_static": 0.00877,
    "dysbiotic_hobic": 0.00213, "dysbiotic_static": 0.00155,
}
GOLDEN_SIGMA_IMPLANT = {
    "commensal_hobic": 0.00613, "commensal_static": 0.00428,
    "dysbiotic_hobic": 0.00130, "dysbiotic_static": 0.00106,
}
HEADLINE_RATIO = 6.44   # sigma_CH / sigma_DH (tooth MAP)

# --- trustworthy _flat_golden entries (all conditions EXCEPT the stale DS) -----
GOLDEN_FLAT = {
    "tooth_commensal_hobic_mises_max": 0.013715936802327633,
    "tooth_commensal_hobic_mises_mean": 0.0023347710769373383,
    "tooth_commensal_hobic_alpha_max": 1.8026622533798218,
    "tooth_commensal_static_mises_max": 0.008766991086304188,
    "tooth_commensal_static_mises_mean": 0.0015552624695879757,
    "tooth_commensal_static_alpha_max": 1.5454940795898438,
    "tooth_dysbiotic_hobic_mises_max": 0.002130083041265607,
    "tooth_dysbiotic_hobic_mises_mean": 0.0004195389666694968,
    "tooth_dysbiotic_hobic_alpha_max": 0.8981017470359802,
    "implant_commensal_hobic_mises_max": 0.0061302389949560165,
    "implant_commensal_hobic_mises_mean": 0.001409647256472921,
    "implant_commensal_hobic_alpha_max": 1.682484745979309,
    "implant_commensal_static_mises_max": 0.0042835925705730915,
    "implant_commensal_static_mises_mean": 0.0009916232741983006,
    "implant_commensal_static_alpha_max": 1.4424611330032349,
    "implant_dysbiotic_hobic_mises_max": 0.0012987378286197782,
    "implant_dysbiotic_hobic_mises_mean": 0.0003112205559087755,
    "implant_dysbiotic_hobic_alpha_max": 0.8382282853126526,
}
# Stale (pre-fix) DS entries still sitting in _flat_golden, vs the corrected value.
STALE_DS = {
    "tooth_dysbiotic_static_mises_max": (0.013625433668494225, 0.00155),
    "implant_dysbiotic_static_mises_max": (0.00610277010127902, 0.00106),
}

_have_ci = _CI_TOOTH.exists() and _CI_IMPLANT.exists()
pytestmark = pytest.mark.skipif(not _have_ci, reason="stress-CI artifacts absent")


@pytest.fixture(scope="module")
def ci():
    return (json.loads(_CI_TOOTH.read_text()), json.loads(_CI_IMPLANT.read_text()))


@pytest.mark.parametrize("cond", CONDS)
def test_tooth_sigma_frozen(ci, cond):
    tooth, _ = ci
    assert tooth[cond]["sigma_map"] == pytest.approx(GOLDEN_SIGMA_TOOTH[cond], rel=1e-4)


@pytest.mark.parametrize("cond", CONDS)
def test_implant_sigma_frozen(ci, cond):
    _, implant = ci
    assert implant[cond]["sigma_map"] == pytest.approx(GOLDEN_SIGMA_IMPLANT[cond], rel=1e-4)


def test_headline_ratio(ci):
    tooth, _ = ci
    ratio = tooth["commensal_hobic"]["sigma_map"] / tooth["dysbiotic_hobic"]["sigma_map"]
    assert ratio == pytest.approx(HEADLINE_RATIO, abs=0.02), ratio


def test_corrected_DS_not_regressed_to_bug(ci):
    """The DS fix must hold: corrected DS stress ~1.55/1.06 kPa, well below the
    ~8x buggy 13.6/6.1 that the copy-bug produced."""
    tooth, implant = ci
    assert tooth["dysbiotic_static"]["sigma_map"] < 0.004      # not the ~0.0136 bug
    assert implant["dysbiotic_static"]["sigma_map"] < 0.003    # not the ~0.0061 bug


@pytest.mark.skipif(not _CMP.exists(), reason="comparison json absent")
@pytest.mark.parametrize("key", sorted(GOLDEN_FLAT))
def test_flat_golden_frozen(key):
    flat = json.loads(_CMP.read_text())["_flat_golden"]
    assert flat[key] == pytest.approx(GOLDEN_FLAT[key], rel=1e-9), (
        f"{key} drifted — update GOLDEN_FLAT in the same commit if intended.")


@pytest.mark.skipif(not _CMP.exists(), reason="comparison json absent")
@pytest.mark.parametrize("key", sorted(STALE_DS))
@pytest.mark.xfail(strict=True,
                   reason="tooth_klempt_comparison.json _flat_golden DS entries "
                          "are stale (pre DS-fix, ~8x). Regenerate the artifact "
                          "post-fix; then this xpasses -> move DS into GOLDEN_FLAT "
                          "and drop STALE_DS.")
def test_comparison_json_DS_is_stale(key):
    """Self-clearing tracker: asserts the comparison-json DS entry equals the
    corrected value. It currently does NOT (stale) -> xfail. When the artifact
    is regenerated it will -> xpass -> strict failure prompts cleanup."""
    flat = json.loads(_CMP.read_text())["_flat_golden"]
    _stale, corrected = STALE_DS[key]
    assert flat[key] == pytest.approx(corrected, rel=0.05)
