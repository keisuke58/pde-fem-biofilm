"""Every top-level document must be findable from the guided tour.

REPO_MAP.md is the entry point a new reader (or a co-supervisor, or the
author in three months) is pointed at. Documents get added far more often
than the map gets updated, and nothing fails when they diverge -- on
2026-09-01 ten new files had accumulated, none of them listed anywhere. By
then the map is not just incomplete, it is misleading: it reads as a complete
inventory.

This is the cheapest possible fix. It does not check that the description is
good, only that the file is mentioned somewhere a reader would look.
"""
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_INDEXES = ("REPO_MAP.md", "DOCS.md", "README.md")

# Documents that are deliberately not in the tour. Keep this short, and give
# each one a reason -- an exemption without a reason is how the check rots.
_EXEMPT = {
    "CLAUDE.md": "instructions for the assistant, not repository documentation",
    "REPO_MAP.md": "is the index",
    "DOCS.md": "is an index",
    "README.md": "is the front page",
}


def _indexed_text():
    out = []
    for name in _INDEXES:
        p = _ROOT / name
        if p.exists():
            out.append(p.read_text(encoding="utf-8", errors="replace"))
    assert out, "no index document found at all"
    return "\n".join(out)


def test_every_top_level_document_is_reachable_from_the_index():
    text = _indexed_text()
    missing = sorted(
        p.name for p in _ROOT.glob("*.md")
        if p.name not in _EXEMPT and p.name not in text)
    assert not missing, (
        "these top-level documents are not mentioned in "
        f"{', '.join(_INDEXES)}: {missing}\n"
        "Add a line to REPO_MAP.md, or add the file to _EXEMPT here with a "
        "reason. An index that silently goes stale is worse than no index, "
        "because it still reads as complete.")


def test_every_top_level_directory_of_substance_is_reachable():
    """Same argument for directories that hold work rather than artefacts."""
    text = _indexed_text()
    skip = {".github", ".claude", "__pycache__", "assets", "data",
            "docs", "runs", "scratch", "node_modules", ".pytest_cache"}

    def of_substance(d):
        if not d.is_dir() or d.name.startswith(".") or d.name in skip:
            return False
        # holds work, not just outputs
        return any(d.rglob("*.py")) or any(d.rglob("*.f")) or any(d.glob("*.md"))

    missing = sorted(d.name for d in _ROOT.iterdir()
                     if of_substance(d) and d.name not in text)
    assert not missing, (
        f"these directories are not mentioned in {', '.join(_INDEXES)}: "
        f"{missing}")


@pytest.mark.parametrize("name", sorted(_EXEMPT))
def test_exemptions_still_exist(name):
    """An exemption for a deleted file quietly widens the check's blind spot."""
    assert (_ROOT / name).exists(), (
        f"{name} is exempt from the index check but no longer exists; "
        "remove it from _EXEMPT")
