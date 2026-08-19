"""Tests for the Abaqus cost extractor.

The real .sta/.msg/.dat files live on the author's workstation and are not in
the repository, so the parsers are exercised against fixtures written in the
exact column layout Abaqus/Standard emits. That is the point of the test: the
script has to be right *before* it is pointed at data that cannot be replayed
here.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from extract_abaqus_cost import collect_job, parse_dat, parse_msg, parse_sta, to_markdown

_ROOT = Path(__file__).resolve().parents[1]

# Increment 3 is attempted twice (attempt 1 cut back, attempt 2 converged),
# which is the case a naive line-count gets wrong.
_STA = """\

                               STEP     INC ATT SEVERE EQUIL TOTAL  TOTAL    STEP       INC OF
                                                DISCON ITERS ITERS  TIME/    TIME/LPF  TIME/LPF
                                                ITERS               FREQ
     1     1   1     0      2     2   0.100      0.100     0.100
     1     2   1     0      3     3   0.300      0.300     0.200
     1     3   1     1      5     6   0.300      0.300     0.200
     1     3   2     0      4     4   0.500      0.500     0.100
     1     4   1     0      2     2   1.000      1.000     0.500
"""

_MSG = """\
 ***WARNING: THE SYSTEM MATRIX HAS 3 NEGATIVE EIGENVALUES.

                              JOB TIME SUMMARY
  USER TIME (SEC)      =        1200.0
  SYSTEM TIME (SEC)    =          46.8
  TOTAL CPU TIME (SEC) =        1246.8
  WALLCLOCK TIME (SEC) =           456
"""

_DAT = """\
 NUMBER OF ELEMENTS IS 43080
 NUMBER OF NODES IS 8985
 TOTAL NUMBER OF VARIABLES IN THE MODEL: 26955
"""


@pytest.fixture
def job(tmp_path):
    stem = tmp_path / "p23_klempt_A_commensal_hobic"
    stem.with_suffix(".sta").write_text(_STA)
    stem.with_suffix(".msg").write_text(_MSG)
    stem.with_suffix(".dat").write_text(_DAT)
    return stem


def test_sta_counts_converged_increments_not_attempts(job):
    sta = parse_sta(job.with_suffix(".sta"))
    assert sta["increments"] == 4, "increment 3 was retried; it is still one increment"
    assert sta["attempts"] == 5
    assert sta["cutbacks"] == 1
    assert sta["total_iters"] == 2 + 3 + 6 + 4 + 2
    assert sta["severe_discon_iters"] == 1
    assert sta["final_step_time"] == pytest.approx(1.0)


def test_sta_ignores_header_and_banner_lines(job):
    """The header rows contain words, not six leading integers."""
    rows = parse_sta(job.with_suffix(".sta"))["rows"]
    assert len(rows) == 5
    assert all(isinstance(r["step"], int) for r in rows)


def test_msg_timing_and_diagnostics(job):
    msg = parse_msg(job.with_suffix(".msg"))
    assert msg["wallclock_s"] == pytest.approx(456.0)
    assert msg["total_cpu_s"] == pytest.approx(1246.8)
    assert msg["user_cpu_s"] == pytest.approx(1200.0)
    assert msg["warnings"] == 1
    assert msg["errors"] == 0


def test_dat_model_size(job):
    dat = parse_dat(job.with_suffix(".dat"))
    assert dat["n_elements"] == 43080
    assert dat["n_nodes"] == 8985
    assert dat["n_dof"] == 26955


def test_derived_cost_metrics(job):
    j = collect_job(job)
    assert j["iters_per_increment"] == pytest.approx(17 / 4, abs=0.01)
    assert j["s_per_iteration"] == pytest.approx(456 / 17, abs=0.01)
    # 456 s over 43080 elements x 17 iterations
    assert j["us_per_element_iteration"] == pytest.approx(
        1e6 * 456 / (43080 * 17), rel=1e-3
    )


def test_missing_files_degrade_without_error(tmp_path):
    """A job with only a .sta must still report iteration counts."""
    stem = tmp_path / "sta_only"
    stem.with_suffix(".sta").write_text(_STA)
    j = collect_job(stem)
    assert j["increments"] == 4
    assert j["wallclock_s"] is None
    assert j["s_per_iteration"] is None, "must not invent a timing"
    assert j["n_elements"] is None


def test_markdown_renders_missing_values(tmp_path):
    stem = tmp_path / "sta_only"
    stem.with_suffix(".sta").write_text(_STA)
    md = to_markdown([collect_job(stem)])
    assert "| Job |" in md and "sta_only" in md
    assert "—" in md, "absent values render as a dash, never as 0"


def test_cli_scans_a_directory(job, tmp_path):
    out = tmp_path / "cost.json"
    rc = subprocess.run(
        [sys.executable, str(_ROOT / "extract_abaqus_cost.py"),
         str(tmp_path), "-o", str(out), "--markdown"],
        capture_output=True, text=True,
    )
    assert rc.returncode == 0, rc.stderr
    payload = json.loads(out.read_text())
    assert len(payload["jobs"]) == 1
    assert payload["jobs"][0]["job"] == "p23_klempt_A_commensal_hobic"
    assert payload["jobs"][0]["wallclock_s"] == pytest.approx(456.0)


def test_cli_reports_when_nothing_found(tmp_path):
    rc = subprocess.run(
        [sys.executable, str(_ROOT / "extract_abaqus_cost.py"), str(tmp_path)],
        capture_output=True, text=True,
    )
    assert rc.returncode == 1
    assert "no Abaqus job files" in rc.stderr
