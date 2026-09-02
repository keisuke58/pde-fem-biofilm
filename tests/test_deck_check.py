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
_CHECK = _ROOT / "ansys_usermat" / "apdl" / "check_deck.py"


def _run(deck_text, tmp_path, name="d.dat"):
    p = tmp_path / name
    p.write_text(deck_text)
    r = subprocess.run([sys.executable, str(_CHECK), str(p)],
                       capture_output=True, text=True, timeout=60)
    return r.returncode, r.stdout


# A well-formed deck: 10 state slots, slot 10 written, TBDATA split at six.
_STATE_OK = "TB,STATE,1,,10\nTBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0\nTBDATA,7,0.0,0.0,1.0\nTBDATA,10,0.05\n"

ELASTIC = ("TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\n" + _STATE_OK +
           "TIME,5.0\nNSUBST,4,20,1\n")
_VISC = "TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.008,0.0,0.0\n" + _STATE_OK
VISCOUS_COARSE = _VISC + "TIME,5.0\nNSUBST,4,20,1\n"
VISCOUS_FINE = _VISC + "TIME,5.0\nNSUBST,200,500,200\n"
VISCOUS_REFUSED = _VISC + "TIME,50.0\nNSUBST,1,1,1\n"

# The silent-failure decks. Each of these solves without complaint in ANSYS
# and returns a stress field that looks entirely normal.
TOO_FEW_STATEV = ("TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\n"
                  "TB,STATE,1,,9\nTBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0\n"
                  "TBDATA,7,0.0,0.0,1.0\nTIME,5.0\nNSUBST,4,20,1\n")
FV_OUT_OF_BOUNDS = ("TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\n"
                    "TB,STATE,1,,6\nTBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0\n"
                    "TIME,5.0\nNSUBST,4,20,1\n")
ALPHA_NEVER_WRITTEN = ("TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\n"
                       "TB,STATE,1,,10\nTBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0\n"
                       "TBDATA,7,0.0,0.0,1.0\nTIME,5.0\nNSUBST,4,20,1\n")
TBDATA_OVERLONG = ("TB,USER,1,1,6\nTBDATA,1,0.2E-3,0.0,5.0E3,0.0,0.0,0.0\n"
                   "TB,STATE,1,,10\n"
                   "TBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0\n"
                   "TBDATA,10,0.05\nTIME,5.0\nNSUBST,4,20,1\n")


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


# --------------------------------------------------------------------------- #
# The silent failures. Each solves cleanly in ANSYS and looks right.

def test_too_few_state_slots_silently_disables_growth(tmp_path):
    """The worst of them. usermat reads alpha only if nStatev >= 10, so nine
    slots means Fg = I -- a purely elastic solve, reported as a growth one, in
    a thesis about growth-induced stress."""
    code, out = _run(TOO_FEW_STATEV, tmp_path)
    assert code == 1, out
    assert "PURELY ELASTIC" in out


def test_far_too_few_state_slots_is_reported_as_out_of_bounds(tmp_path):
    code, out = _run(FV_OUT_OF_BOUNDS, tmp_path)
    assert code == 1, out
    assert "past the array" in out


def test_declaring_the_slot_without_writing_it_is_still_no_growth(tmp_path):
    """Ten slots declared, none of them filled with alpha: same outcome, and
    even easier to miss because the declaration looks correct."""
    code, out = _run(ALPHA_NEVER_WRITTEN, tmp_path)
    assert code == 1, out
    assert "no growth" in out


def test_an_overlong_TBDATA_is_caught(tmp_path):
    """APDL takes six values and drops the rest without complaint -- a real
    deck bug in this repository's history."""
    code, out = _run(TBDATA_OVERLONG, tmp_path)
    assert code == 1, out
    assert "drops the rest" in out


def test_a_well_formed_deck_passes_all_of_the_above(tmp_path):
    code, out = _run(ELASTIC, tmp_path)
    assert code == 0, out
    for phrase in ("PURELY ELASTIC", "past the array", "no growth", "drops the rest"):
        assert phrase not in out


# --------------------------------------------------------------------------- #
# Abaqus decks. Same silent failure, different route: growth arrives as
# temperature (ALPHA_GROWTH = TEMP + DTEMP), not as a state variable.

_ABQ_HEAD = "*MATERIAL, NAME=BIO\n*USER MATERIAL, CONSTANTS=5\n2e-4,0.,5e3,8e-3,0.\n*DEPVAR\n9\n"
ABQ_OK = _ABQ_HEAD + "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nALL,0.0\n*STEP\n*TEMPERATURE\nALL,0.05\n*END STEP\n"
ABQ_NO_TEMPERATURE = _ABQ_HEAD + "*STEP\n*STATIC\n*END STEP\n"
ABQ_INITIAL_ONLY = _ABQ_HEAD + "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nALL,0.05\n*STEP\n*STATIC\n*END STEP\n"
ABQ_NOT_OURS = "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000.,0.3\n"


def test_an_abaqus_deck_without_a_temperature_field_grows_by_zero(tmp_path):
    """The Abaqus counterpart of too-few state slots: the UMAT reads growth
    from TEMP, so a job with no temperature field is purely elastic and says
    nothing about it."""
    code, out = _run(ABQ_NO_TEMPERATURE, tmp_path, name="d.inp")
    assert code == 1, out
    assert "PURELY ELASTIC" in out


def test_a_well_formed_abaqus_deck_passes(tmp_path):
    code, out = _run(ABQ_OK, tmp_path, name="d.inp")
    assert code == 0, out


def test_initial_temperature_alone_is_called_out(tmp_path):
    """Not an error -- a constant alpha is a legitimate job -- but worth
    saying, because it is indistinguishable from a missing step field until
    the results come out flat."""
    code, out = _run(ABQ_INITIAL_ONLY, tmp_path, name="d.inp")
    assert code == 0, out
    assert "never changes" in out


def test_a_deck_that_is_not_ours_is_not_flagged(tmp_path):
    """False positives are how a checker stops being read. A deck with no
    *USER MATERIAL is somebody else's problem."""
    code, out = _run(ABQ_NOT_OURS, tmp_path, name="d.inp")
    assert code == 0, out
    assert "not one of ours" in out
