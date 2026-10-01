#!/usr/bin/env python3
"""Apply our agreed edits to the partner's call-site file, on IKMHIWI03.

The partner's source never enters this public repository, so the edits travel
as this script instead of as a diff. Run it against the working copy:

    python ansys_usermat/apdl/apply_partner_patches.py ^
        F:\\biofilm_upf_wired\\Usermat_P21-V21_Conection_Test.F

What it does, and nothing else:

  * **typo** -- ``Sdp_sumBio = Sdp_bio1_n + Sdp_bio1_n`` becomes
    ``... + Sdp_bio2_n``. Confirmed as a typo with Oliver on 2026-10-01, and
    he agreed we fix it in our working copy.

It refuses rather than guesses: the target line must occur exactly once, or
the file is left untouched and the script exits non-zero. Re-running on an
already-patched file is a no-op. A timestamped backup is written before any
change. Fixed-form Fortran is respected (comment in column 1, code from
column 7, nothing past column 72).

The one-species branch (``prop(6) = 2``, see ONE_SPECIES_COUPLING.md) is NOT
inserted here yet: where it goes depends on code around the existing ecology
call that this repository has not seen, and on which pool variable holds the
Gauss-point phi -- both still open.
"""
from __future__ import annotations

import re
import shutil
import sys
import time
from pathlib import Path

TYPO = re.compile(r"^(?P<ind>[ \t]*)Sdp_sumBio\s*=\s*Sdp_bio1_n\s*\+\s*Sdp_bio1_n"
                  r"(?P<tail>\s*(!.*)?)$", re.IGNORECASE)
FIXED = re.compile(r"^\s*Sdp_sumBio\s*=\s*Sdp_bio1_n\s*\+\s*Sdp_bio2_n\b",
                   re.IGNORECASE)
NOTE = ("C     typo fixed (was bio1_n + bio1_n): confirmed with Oliver and "
        "agreed\nC     2026-10-01, see pde-fem-biofilm ONE_SPECIES_COUPLING.md\n")


class PatchError(RuntimeError):
    pass


def patch_text(text: str) -> tuple[str, str]:
    """Return (new_text, status) where status is 'patched' or 'already'."""
    lines = text.splitlines(keepends=True)
    hits = [i for i, ln in enumerate(lines) if TYPO.match(ln.rstrip("\r\n"))]
    done = [i for i, ln in enumerate(lines) if FIXED.match(ln)]
    if not hits and len(done) == 1:
        return text, "already"
    if len(hits) != 1:
        raise PatchError(f"expected exactly one 'Sdp_sumBio = Sdp_bio1_n + "
                         f"Sdp_bio1_n' line, found {len(hits)} -- not touching "
                         "the file")
    i = hits[0]
    eol = "\r\n" if lines[i].endswith("\r\n") else "\n"
    m = TYPO.match(lines[i].rstrip("\r\n"))
    new = f"{m.group('ind')}Sdp_sumBio = Sdp_bio1_n + Sdp_bio2_n{m.group('tail')}"
    code = new.split("!")[0].rstrip()
    if len(code.expandtabs(6)) > 72:
        raise PatchError("patched line would pass column 72")
    lines[i] = NOTE.replace("\n", eol) + new + eol
    return "".join(lines), "patched"


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    path = Path(argv[0])
    text = path.read_text(encoding="latin-1")
    try:
        new, status = patch_text(text)
    except PatchError as e:
        print(f"REFUSED: {e}")
        return 1
    if status == "already":
        print("typo: already fixed, nothing to do")
        return 0
    backup = path.with_name(path.name + time.strftime(".orig-%Y%m%d-%H%M%S"))
    shutil.copy2(path, backup)
    path.write_text(new, encoding="latin-1", newline="")
    print(f"typo: fixed. backup -> {backup}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
