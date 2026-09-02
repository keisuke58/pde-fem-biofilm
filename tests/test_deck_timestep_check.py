"""The deck time-step checker itself, not the decks.

Kept separate from a check of the repository's own decks on purpose. One of
those decks currently permits a step that costs ~30% of the von Mises stress,
and whether to change it is a decision about that deck's stored comparison,
not something a test should force at import time. What must not rot is the
checker: if it stopped detecting the case it was written for, nobody would
notice, because its output is only read when someone remembers to run it.
"""
import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_CHECK = _ROOT / "ansys_usermat" / "apdl" / "check_deck_timestep.py"


def _run(deck_text, tmp_path, name="d.dat"):
    p = tmp_path / name
    p.write_text(deck_text)
    r = subprocess.run([sys.executable, str(_CHECK), str(p)],
                       capture_output=True, text=True, timeout=60)
    return r.returncode, r.stdout


ELASTIC = "TBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\nTIME,5.0\nNSUBST,4,20,1\n"
VISCOUS_COARSE = "TBDATA,1,0.2E-3,0.0,5.0E3,0.008,0.0,0.0\nTIME,5.0\nNSUBST,4,20,1\n"
VISCOUS_FINE = "TBDATA,1,0.2E-3,0.0,5.0E3,0.008,0.0,0.0\nTIME,5.0\nNSUBST,200,500,200\n"
VISCOUS_REFUSED = "TBDATA,1,0.2E-3,0.0,5.0E3,0.008,0.0,0.0\nTIME,50.0\nNSUBST,1,1,1\n"


def test_an_elastic_deck_is_never_constrained(tmp_path):
    code, out = _run(ELASTIC, tmp_path)
    assert code == 0
    assert "no relaxation time" in out


def test_a_coarse_viscous_deck_is_flagged(tmp_path):
    """The case the checker exists for: within the routine's guard, so nothing
    refuses it, but far too coarse for a number that gets reported."""
    code, out = _run(VISCOUS_COARSE, tmp_path)
    assert code == 1, out
    assert "too coarse to report" in out
    assert "tau = 20" in out


def test_it_reads_the_third_NSUBST_argument_as_the_largest_step(tmp_path):
    """NSBMN is the fewest substeps, so it sets the biggest step AUTOTS may
    take -- the opposite of the intuitive reading, and the whole reason the
    coarse case above is easy to miss."""
    _, out = _run(VISCOUS_COARSE, tmp_path)
    line = [l for l in out.splitlines() if "coarsest" in l]
    assert line and "dt =         5" in line[0], out


def test_a_resolved_viscous_deck_passes(tmp_path):
    code, out = _run(VISCOUS_FINE, tmp_path)
    assert code == 0, out
    assert "too coarse" not in out


def test_a_step_past_the_guard_is_reported_as_refused(tmp_path):
    code, out = _run(VISCOUS_REFUSED, tmp_path)
    assert code == 1
    assert "REFUSED" in out
