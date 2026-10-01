"""paste_fragments.py -- put the repo's phi-mode fragments into the partner's
usermat (a local working copy that is never committed), the same way every time.

    python ansys_usermat/apdl/paste_fragments.py F:\\biofilm_upf_wired\\Usermat_P21-V21_v222.F

What it does:
  1. backs the target up as <target>.prepaste-<timestamp>;
  2. makes sure "use biofilm_split" follows "use biofilm_py_bridge";
  3. replaces the declaration block and the executable block with
     callsite/phi_mode_decl.inc and callsite/phi_mode_exec.inc, wrapped in
     BEGIN/END marker lines so the next run finds them exactly;
  4. copies callsite/split_rates.f next to the target (module biofilm_split);
  5. reads the result back and checks each block equals its fragment.

Without markers (a target pasted by hand) the old blocks are found from their
first comment line; the declaration block ends before the partner's
"!-----" separator, the executable block before CALL BIOFILM_GROWTH_VISCO_V01.
Anything ambiguous (zero or several matches) stops the script and leaves the
target untouched.
"""
from __future__ import annotations

import re
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CALLSITE = HERE / "callsite"

MARK = {
    "decl": ("C     >>> pde-fem-biofilm fragment BEGIN decl",
             "C     <<< pde-fem-biofilm fragment END decl"),
    "exec": ("C     >>> pde-fem-biofilm fragment BEGIN exec",
             "C     <<< pde-fem-biofilm fragment END exec"),
}
USE_BRIDGE = re.compile(r"^\s*use\s+biofilm_py_bridge\b", re.I)
USE_SPLIT = re.compile(r"^\s*use\s+biofilm_split\b", re.I)


def fail(msg: str) -> None:
    sys.exit(f"STOP (target unchanged): {msg}")


def read_lines(p: Path) -> list[str]:
    return p.read_text(encoding="latin-1").replace("\r\n", "\n").split("\n")


def one(lines: list[str], pred, what: str, start: int = 0) -> int:
    hits = [i for i in range(start, len(lines)) if pred(lines[i])]
    if len(hits) != 1:
        fail(f"{what}: {len(hits)} matches, expected 1")
    return hits[0]


def find_block(lines, kind):
    b, e = MARK[kind]
    if any(l == b for l in lines):
        i = one(lines, lambda l: l == b, f"{kind} BEGIN marker")
        j = one(lines, lambda l: l == e, f"{kind} END marker")
        if j < i:
            fail(f"{kind} END marker before BEGIN")
        return i, j + 1
    if kind == "decl":
        i = one(lines, lambda l: l.startswith("C     --- phi mode, declarations"), "decl start")
        j = next((k for k in range(i, len(lines)) if re.match(r"^\s*!-{20,}", lines[k])), None)
        if j is None:
            fail("decl end (partner's !----- separator) not found")
        if "PM_DTMAX" not in "".join(lines[i:j]) and "CM_" not in "".join(lines[i:j]):
            fail("decl block does not look like a previous paste")
        return i, j
    i = one(lines, lambda l: l.startswith("C     --- phi mode, executable part"), "exec start")
    j = one(lines, lambda l: "CALL BIOFILM_GROWTH_VISCO_V01(" in l.upper(), "growth call", i)
    while j > i and lines[j - 1].strip() == "":
        j -= 1
    return i, j


def main(target: Path) -> None:
    lines = read_lines(target)
    frag = {k: read_lines(CALLSITE / f"phi_mode_{k}.inc") for k in ("decl", "exec")}
    for k in frag:
        while frag[k] and frag[k][-1] == "":
            frag[k].pop()

    # locate both blocks on the original, replace the later one first
    d0, d1 = find_block(lines, "decl")
    x0, x1 = find_block(lines, "exec")
    if not d1 <= x0:
        fail("decl block is not before exec block")
    out = (lines[:d0] + [MARK["decl"][0]] + frag["decl"] + [MARK["decl"][1]]
           + lines[d1:x0] + [MARK["exec"][0]] + frag["exec"] + [MARK["exec"][1]]
           + lines[x1:])

    if not any(USE_SPLIT.match(l) for l in out):
        u = one(out, USE_BRIDGE.match, "use biofilm_py_bridge")
        out.insert(u + 1, "      use biofilm_split")

    stamp = time.strftime("%Y%m%d-%H%M%S")
    backup = target.with_name(target.name + f".prepaste-{stamp}")
    shutil.copy2(target, backup)
    target.write_bytes("\n".join(out).encode("latin-1"))
    shutil.copy2(CALLSITE / "split_rates.f", target.parent / "split_rates.f")

    # read back and verify
    back = read_lines(target)
    for k in ("decl", "exec"):
        i, j = find_block(back, k)
        if back[i + 1:j - 1] != frag[k]:
            fail(f"read-back of {k} differs from the fragment (backup: {backup})")
    if sum(bool(USE_SPLIT.match(l)) for l in back) != 1:
        fail("use biofilm_split not present exactly once")
    print(f"OK: {target}\n  backup  {backup.name}\n  decl {len(frag['decl'])} lines, "
          f"exec {len(frag['exec'])} lines, split_rates.f copied")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]))
