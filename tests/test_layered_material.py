"""The alpha(depth) -> layered-material generator.

The ANSYS side reads the growth driver from a per-material state slot, so it
takes alpha piecewise-uniformly; the Abaqus side takes a continuous field
through the temperature slot. Representing a depth-resolved growth field on the
ANSYS route therefore means binning it into layers, and the point of these
tests is that the binning is a controlled approximation rather than an
unexamined one -- including that what it emits actually passes the deck
pre-flight, which is the closed loop that matters.
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "ansys_usermat" / "apdl"))
import make_layered_material as mlm  # noqa: E402
_GEN = _ROOT / "ansys_usermat" / "apdl" / "make_layered_material.py"
_CHECK = _ROOT / "ansys_usermat" / "apdl" / "check_deck.py"


@pytest.fixture
def field(tmp_path):
    """A nutrient-limited profile: growth concentrated at the fluid side."""
    x = np.linspace(0.0, 0.2, 30)
    a = 0.05 * np.exp(-((0.2 - x) / 0.07) ** 2)
    p = tmp_path / "alpha_field_1d.csv"
    np.savetxt(p, np.c_[x / 0.2, x, np.ones_like(x), np.ones_like(x), a, a / 3],
               delimiter=",",
               header="x_norm,x_mm,c_final,phi_total_final,alpha_x,eps_growth_x",
               comments="")
    return p


def _gen(field, *args):
    r = subprocess.run([sys.executable, str(_GEN), str(field), *args],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_binning_error_falls_as_layers_are_added(field):
    """The trade the generator exists to make, stated as a number rather than
    a guess. One layer is the spatially averaged alpha the 0-D route uses."""
    out = _gen(field, "--convergence")
    errs = {}
    for line in out.splitlines():
        parts = line.replace("!", "").split()
        if len(parts) == 2 and parts[0].isdigit():
            errs[int(parts[0])] = float(parts[1])
    assert set((1, 2, 4, 8, 16)) <= set(errs)
    ns = sorted(errs)
    assert all(errs[a] > errs[b] for a, b in zip(ns, ns[1:])), errs
    assert errs[1] > 0.5, "a single layer should be a poor fit to a graded field"
    assert errs[16] < 0.1


def test_the_emitted_block_matches_the_routine_s_property_layout(field):
    """prop is (E, nu, eta, C01/C10, mtype) and the state block is ten slots
    with alpha in the tenth -- the layout usermat_wrapper_v01_smoketest.f
    unpacks. Getting either wrong is silent."""
    out = _gen(field, "-n", "2")
    assert out.count("TB,USER,") == 2
    assert all(f"TB,USER,{n},1,5" in out for n in (1, 2))
    assert out.count("TB,STATE,") == 2
    assert "TBDATA,10," in out


def test_no_TBDATA_call_exceeds_six_values(field):
    """APDL drops the seventh value onward without complaining."""
    out = _gen(field, "-n", "6")
    for line in out.splitlines():
        if line.strip().upper().startswith("TBDATA,"):
            vals = [v for v in line.split("!")[0].split(",")[2:] if v.strip()]
            assert len(vals) <= 6, line


def test_a_deck_built_from_the_output_passes_the_pre_flight(field, tmp_path):
    """The closed loop: what the generator emits must survive the checker that
    every deck is supposed to go through before its numbers are reported."""
    block = _gen(field, "-n", "3", "--substrate")
    deck = tmp_path / "generated.dat"
    deck.write_text("/PREP7\nET,1,SOLID185\n" + block +
                    "\nTIME,1.0\nNSUBST,10,50,10\nSOLVE\n")
    r = subprocess.run([sys.executable, str(_CHECK), str(deck)],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stdout + r.stderr
    for bad in ("PURELY ELASTIC", "past the array", "no growth", "drops the rest"):
        assert bad not in r.stdout, r.stdout


def test_the_substrate_layer_does_not_grow(field):
    out = _gen(field, "-n", "2", "--substrate")
    first = out.split("! ---- material 2")[0]
    assert "TBDATA,10,0  " in first or "TBDATA,10,0 " in first, first


# ---------------------------------------------------------------------------
# The producer/consumer contract on the CSV columns.
#
# This is where a real bug lived. The generator's candidate list was
# ("alpha_x", "alpha"), its docstring said the extractor writes `alpha_x`, and
# the fixture above builds a CSV with `alpha_x` -- so everything agreed with
# everything, and none of it agreed with `extract_alpha_field_1d.py`, which
# writes `alpha_x_monod`. The generator could not read a single real file it
# was written for, and the tests could not see it because they asserted against
# an invented schema rather than the producer's.
#
# So these read the header out of the extractor's source. A rename on either
# side now fails here instead of at the machine.
# ---------------------------------------------------------------------------
_EXTRACTOR = _ROOT / "JAXFEM" / "extract_alpha_field_1d.py"
_MACRO = _ROOT / "generate_abaqus_eigenstrain.py"


def _extractor_header():
    m = re.search(r'^\s*header\s*=\s*"([^"]+)"', _EXTRACTOR.read_text(),
                  re.MULTILINE)
    assert m, f"no header literal found in {_EXTRACTOR.name}"
    return m.group(1).split(",")


def test_generator_reads_the_extractors_actual_columns(tmp_path):
    cols = _extractor_header()
    n = len(cols)
    data = np.tile(np.linspace(0.1, 0.9, 12)[:, None], (1, n))
    p = tmp_path / "from_extractor.csv"
    np.savetxt(p, data, delimiter=",", header=",".join(cols), comments="")

    # Exactly one alpha column, or the generator would rightly refuse to guess.
    alpha_cols = [c for c in cols if c in mlm.ALPHA_COLS]
    assert alpha_cols == [c for c in cols if c.startswith("alpha")], cols
    assert len(alpha_cols) == 1, alpha_cols
    # Depth is resolved by the generator's own preference order, not by the
    # order the columns happen to appear in the file.
    depth_col = next(c for c in mlm.DEPTH_COLS if c in cols)

    out = _gen(p, "--convergence")
    assert f"({depth_col} vs {alpha_cols[0]})" in out, out


def test_two_alpha_columns_are_refused_rather_than_guessed(tmp_path):
    """`alpha_x` and `alpha_x_monod` are different physical fields -- the
    unweighted and Monod-weighted growth drivers. Picking one by list order
    would be a silent choice between them, so it is an error instead."""
    cols = ["x_mm", "alpha_x", "alpha_x_monod"]
    p = tmp_path / "ambiguous.csv"
    np.savetxt(p, np.tile(np.linspace(0.1, 0.9, 8)[:, None], (1, 3)),
               delimiter=",", header=",".join(cols), comments="")

    r = subprocess.run([sys.executable, str(_GEN), str(p), "--convergence"],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode != 0
    assert "--alpha-col" in r.stdout + r.stderr

    # ...and naming one resolves it.
    assert "vs alpha_x_monod" in _gen(p, "--convergence",
                                      "--alpha-col", "alpha_x_monod")


def test_generator_reads_the_per_condition_macro_columns(tmp_path):
    """The other producer, and the one the results chapter actually needs.

    `macro_eigenstrain_<condition>_hybrid.csv` carries the per-condition growth
    field. It spells its columns `depth_mm` / `alpha_monod`, where the 1D
    extractor spells them `x_mm` / `alpha_x_monod`, and the generator had
    neither -- so it could not read either producer. The names are taken from
    the consumer that defines the schema rather than restated here.
    """
    src = _MACRO.read_text()
    block = src[src.index("result = {"):src.index('"path": path')]
    cols = re.findall(r'"(\w+)":\s*d\[', block)
    assert {"depth_mm", "alpha_monod"} <= set(cols), cols

    p = tmp_path / "macro_eigenstrain_commensal_static_hybrid.csv"
    body = np.tile(np.linspace(0.1, 0.9, 12)[:, None], (1, len(cols)))
    with p.open("w") as f:
        f.write("# generated by a test, mimicking the real file's comment block\n")
        f.write(",".join(cols) + "\n")
        for row in body:
            f.write(",".join(f"{v:.17g}" for v in row) + "\n")

    depth_col = next(c for c in mlm.DEPTH_COLS if c in cols)
    alpha_cols = [c for c in mlm.ALPHA_COLS if c in cols]
    assert len(alpha_cols) == 1, alpha_cols
    assert f"({depth_col} vs {alpha_cols[0]})" in _gen(p, "--convergence")
