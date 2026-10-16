"""The unified user element on the implant sector (p3/uel/make_implant_uel_inp.py) at the coarsest level,
2 x 6 x 8, run through the Python emulation of Abaqus (emulate_abaqus.py): the seed means and the largest
von Mises stress of results/implant_2x6x8.json (the values Abaqus reproduced to 4e-7 on fifa, 16 Oct 2026),
plus the convergence table (implant_convergence.py) and the Abaqus comparison (compare_implant_gp.py) on
synthetic data.

Needs gfortran (builds libunified_uel.so in p3/uel if it is missing) and scipy; skipped without them.
The run takes about 10 s.
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
def lib():
    so = UEL / "libunified_uel.so"
    if not so.exists():
        eco = ROOT / "ansys_usermat" / "coupling" / "ecology_native.f"
        subprocess.run(["gfortran", "-O2", "-fPIC", "-c", str(eco), "-o", "ecology_native.o"], cwd=UEL, check=True)
        subprocess.run(["gfortran", "-O2", "-fPIC", "-ffree-line-length-none", "-c", "unified_uel.f90", "-o", "unified_uel.o"],
                       cwd=UEL, check=True)
        subprocess.run(["gfortran", "-shared", "-o", "libunified_uel.so", "unified_uel.o", "ecology_native.o"], cwd=UEL, check=True)
    return so


@pytest.fixture(scope="module")
def eco_case(tmp_path_factory):
    from make_column_inp import eco_file
    p = tmp_path_factory.mktemp("eco") / "uel_eco.txt"
    eco_file(p)
    os.environ["BIOFILM_ECO_CASE"] = str(p)        # emulate_abaqus asserts it at import
    return p


def test_implant_2x6x8_matches_reference(lib, eco_case):
    import make_implant_uel_inp as mi
    ref = json.load(open(UEL / "results" / "implant_2x6x8.json"))
    nodes, elems, layer, ijk = mi.sector(2, 6, 8)
    assert len(elems) == 96 and len(nodes) == 189
    r = mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True)
    for k in ("seed_alpha_minus_1", "seed_vM_over_E", "seed_interior_alpha_minus_1", "seed_interior_vM_over_E",
              "max_vM_over_E", "max_interior_vM_over_E"):
        assert r[k] == pytest.approx(ref[k], rel=1e-6), k
    assert r["seed_elements"] == 48 and r["seed_interior_elements"] == 24
    assert r["max_centroid_r_theta_z"] == pytest.approx(ref["max_centroid_r_theta_z"], abs=1e-6)
    # the interior excludes exactly the rows at z = 0, z = 2 and theta = 0, 90 degrees
    jk = np.array([mi.element_jk(e, 2, 6) for e in range(96)])
    inner = np.array(r["elements"]["interior"], bool)
    assert np.array_equal(inner, (jk[:, 0] > 0) & (jk[:, 0] < 5) & (jk[:, 1] > 0) & (jk[:, 1] < 7))
    # growth reported as alpha - 1 with alpha(0) = 1: positive, largest near the nutrient (top row)
    xi = np.array(r["elements"]["alpha_minus_1"])
    assert xi.min() >= 0.0 and xi[(layer == 0) & (jk[:, 1] == 7)].mean() > xi[(layer == 0) & (jk[:, 1] == 0)].mean()


def test_compare_gp_record_roundtrip(tmp_path):
    """a uel_gp.txt whose Gauss-point values equal the element means of the .json compares to 0"""
    from compare_implant_gp import compare
    ref = json.load(open(UEL / "results" / "implant_2x6x8.json"))
    el = ref["elements"]
    lines = ["# step 1 increment 5 total time 1.0E-8", "# element gp x y z phi c alpha-1 living von_Mises pressure"]
    for e in range(len(el["layer"])):
        for q in range(8):
            lines.append(f"{e + 1:8d}{q + 1:3d}" + "".join(f"{v:16.8E}" for v in (0, 0, 0, 0.5, 0.5, 0.0, 0.3, 0.0, 0.0)))
    lines += ["# step 2 increment 60 total time 2.3750001E-01", "# element gp x y z phi c alpha-1 living von_Mises pressure"]
    for e in range(len(el["layer"])):
        for q in range(8):
            vals = (el["cx"][e], el["cy"][e], el["cz"][e], 0.5, 0.5, el["alpha_minus_1"][e], 0.3, el["vM_over_E"][e], 0.0)
            lines.append(f"{e + 1:8d}{q + 1:3d}" + "".join(f"{v:16.8E}" for v in vals))
    gp = tmp_path / "uel_gp.txt"; gp.write_text("\n".join(lines) + "\n")
    ok, text = compare(gp, UEL / "results" / "implant_2x6x8.json", tol=1e-7)
    assert ok, text
    assert f"{ref['seed_vM_over_E']:.7g}"[:6] in text


def test_convergence_order_recovers_p():
    from implant_convergence import order
    q = lambda f, p: 1.0 + 0.5 * f ** -p  # noqa: E731
    for p0 in (1.0, 2.0):
        p, qinf = order(1, 2, 3, q(1, p0), q(2, p0), q(3, p0))
        assert p == pytest.approx(p0, abs=1e-6) and qinf == pytest.approx(1.0, abs=1e-6)
    assert order(1, 2, 3, 1.0, 1.1, 1.3)[0] is None     # differences growing: no order
