"""The alpha(depth) -> layered-material generator.

The ANSYS side reads the growth driver from a per-material state slot, so it
takes alpha piecewise-uniformly; the Abaqus side takes a continuous field
through the temperature slot. Representing a depth-resolved growth field on the
ANSYS route therefore means binning it into layers, and the point of these
tests is that the binning is a controlled approximation rather than an
unexamined one -- including that what it emits actually passes the deck
pre-flight, which is the closed loop that matters.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
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
