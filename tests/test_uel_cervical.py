"""The cervical-line measure of the tooth-like sector (p3/uel/cervical_study.py): the band at the bottom of the
bonded face and the seed core are non-empty, fixed regions at every level, the concentration ratio is finite, and
with the homeostatic-pressure law on the band's stress is lower than without it. Coarsest level, 2x6x8, bulge 1,
bottom u_z = 0, a third of the GROW step (nstop = 20) so the two runs take about 6 s together.

Needs gfortran (builds libunified_uel.so in p3/uel if it is missing or older than the source) and scipy; skipped without them.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
UEL = ROOT / "p3" / "uel"
sys.path.insert(0, str(UEL))
pytest.importorskip("scipy")
pytestmark = pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")


@pytest.fixture(scope="module")
def lib(tmp_path_factory):
    so = UEL / "libunified_uel.so"
    srcs = (UEL / "unified_uel.f90", ROOT / "ansys_usermat" / "coupling" / "ecology_native.f")
    if not so.exists() or max(x.stat().st_mtime for x in srcs) > so.stat().st_mtime:   # rebuild when stale
        eco = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"
        subprocess.run(["gfortran", "-O2", "-fPIC", "-c", str(eco), "-o", "ecology_native.o"], cwd=UEL, check=True)
        subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-line-length-none", "-c", "unified_uel.f90", "-o", "unified_uel.o"],
                       cwd=UEL, check=True)
        subprocess.run(["gfortran", "-shared", "-o", "libunified_uel.so", "unified_uel.o", "ecology_native.o"], cwd=UEL, check=True)
    if not os.environ.get("BIOFILM_ECO_CASE"):
        from make_column_inp import eco_file
        p = tmp_path_factory.mktemp("eco") / "uel_eco.txt"
        eco_file(p)
        os.environ["BIOFILM_ECO_CASE"] = str(p)
    return so


def test_masks_on_committed_results():
    """band = the seed elements (one layer at 2x6x8, two at 4x12x16) of the bottom rows within 0.25 mm and 15-75
    degrees (4 and 32 elements); the core as implant_convergence (24 and 192); the two do not overlap and the regions
    are the same in z at both levels"""
    import cervical_study as cs
    for name, nband, ncore in (("implant_2x6x8_b1.json", 4, 24), ("implant_4x12x16_b1_sl2.json", 32, 192)):
        r = json.load(open(UEL / "results" / name))
        band, core, seed = cs.masks(r)
        assert band.sum() == nband and core.sum() == ncore and not (band & core).any()
        assert band.sum() <= seed.sum() and core.sum() <= seed.sum()
        z = np.array(r["elements"]["cz"])
        assert z[band].max() <= cs.BAND_Z[1] and z[core].min() >= cs.CORE_Z[0]
        s = cs.measures(r)
        assert np.isfinite(s["ratio"]) and s["ratio"] > 0 and s["band_max_vM_over_E"] >= s["band_vM_over_E"]
    assert cs.tag_of(2, "uz", "fixed", 0.05) == "4x12x16_b1_sl2_ph0.05_of_bz"
    assert cs.tag_of(1, "fixed", "free", 0.0) == "2x6x8_b1_bf"


def test_law_lowers_the_cervical_band(lib):
    import cervical_study as cs
    import make_implant_uel_inp as mi
    nodes, elems, layer, ijk = mi.sector(2, 6, 8, 1.0)
    kw = dict(nu=0.49, md=True, finite=True, seed_layers=1, bottom="uz", outer="free", nstop=20)
    off = cs.measures(dict(mi.run(nodes, elems, layer, ijk, **kw), nr=2, nt=6, nz=8, bulge=1.0))
    on = cs.measures(dict(mi.run(nodes, elems, layer, ijk, p_ref=0.05, **kw), nr=2, nt=6, nz=8, bulge=1.0))
    for s in (off, on):
        assert s["band_elements"] == 4 and s["core_elements"] == 24 and s["bottom"] == "uz"
        assert np.isfinite(s["ratio"]) and s["ratio"] > 0 and np.isfinite(s["band_max_vM_over_E"])
        assert s["band_max_at_r_theta_z"][2] == pytest.approx(0.125) and 15 <= s["band_max_at_r_theta_z"][1] <= 75
    # the gate only reduces growth, so the band's stress is lower with the law (its own alpha - 1 is 1e-6 at this
    # time and set by the nutrient, not by the gate, so it is not compared)
    assert on["band_vM_over_E"] < off["band_vM_over_E"] and on["band_max_vM_over_E"] < off["band_max_vM_over_E"]
    assert on["core_vM_over_E"] < off["core_vM_over_E"] and on["core_alpha_minus_1"] < off["core_alpha_minus_1"]
    assert on["seed_max_p_over_E"] < off["seed_max_p_over_E"]
    # the table and the spread lines render
    text = cs.table([off, on]); assert "| 2x6x8 | uz | free | off |" in text and "| 0.05 |" in text
    assert "ratio" in cs.spread([off, on] + [dict(on, outer="fixed")])
