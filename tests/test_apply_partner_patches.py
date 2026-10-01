"""The partner-file patcher, on synthetic fixed-form snippets only -- the
partner's real source is never in this repository."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys_usermat" / "apdl"))
import apply_partner_patches as ap   # noqa: E402

SNIPPET = (
    "      Sdp_sumBio   = Sdp_bio1_n + Sdp_bio1_n\n"
    "      Sdp_sumLocal = (Sdp_locbio1_n + Sdp_locbio2_n)/2\n")


def test_fixes_the_typo_and_leaves_everything_else():
    new, status = ap.patch_text(SNIPPET)
    assert status == "patched"
    assert "Sdp_sumBio = Sdp_bio1_n + Sdp_bio2_n" in new
    assert "Sdp_bio1_n + Sdp_bio1_n" not in new
    assert "Sdp_sumLocal = (Sdp_locbio1_n + Sdp_locbio2_n)/2" in new


def test_is_idempotent():
    once, _ = ap.patch_text(SNIPPET)
    twice, status = ap.patch_text(once)
    assert status == "already" and twice == once


def test_respects_fixed_form():
    new, _ = ap.patch_text(SNIPPET)
    for ln in new.splitlines():
        assert ln[:1] in ("C", " "), ln              # comment col 1 or code
        assert len(ln.split("!")[0].rstrip()) <= 72 or ln.startswith("C")


def test_keeps_crlf_line_endings():
    new, _ = ap.patch_text(SNIPPET.replace("\n", "\r\n"))
    assert "\r\n" in new and "\n\n" not in new.replace("\r\n", "")


def test_refuses_when_the_line_is_absent():
    with pytest.raises(ap.PatchError):
        ap.patch_text("      Sdp_sumBio = something_else\n")


def test_refuses_when_the_line_is_ambiguous():
    with pytest.raises(ap.PatchError):
        ap.patch_text(SNIPPET + SNIPPET)


def test_end_to_end_writes_a_backup(tmp_path):
    f = tmp_path / "Usermat_test.F"
    f.write_text(SNIPPET)
    assert ap.main([str(f)]) == 0
    assert "Sdp_bio2_n" in f.read_text()
    assert len(list(tmp_path.glob("Usermat_test.F.orig-*"))) == 1
    assert ap.main([str(f)]) == 0                    # second run: no-op
    assert len(list(tmp_path.glob("Usermat_test.F.orig-*"))) == 1
