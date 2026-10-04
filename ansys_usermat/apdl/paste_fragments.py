"""paste_fragments.py -- put the repo's phi-mode fragments into the partner's
usermat (a local working copy that is never committed), the same way every time.

    python ansys_usermat/apdl/paste_fragments.py F:\\biofilm_upf_wired\\Usermat_P21-V21_v222.F
    python ansys_usermat/apdl/paste_fragments.py <target> --nut-var NAME

--nut-var NAME (optional): the partner's Gauss-point nutrient variable, e.g.
the pool value of Nut1. It replaces the fragment's line "CM_NUT = -1.0D0" by
"CM_NUT = NAME", which switches on the local nutrient for prop(33) > 0
(two-way step 1). Without it the pasted code is the same as before.

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


NUT_LINE = "          CM_NUT = -1.0D0"
NUT_VAR = re.compile(r"^[A-Za-z][A-Za-z0-9_]*(\([A-Za-z0-9_, ]+\))?$")


def set_nut_source(exec_lines: list[str], name: str) -> list[str]:
    """Replace the fragment's nutrient source line by CM_NUT = name."""
    if not NUT_VAR.match(name):
        fail(f"--nut-var {name!r} is not a Fortran variable or array element")
    new = f"          CM_NUT = {name}"
    if len(new) > 72:
        fail(f"--nut-var {name!r} makes the line longer than 72 columns")
    hits = [i for i, l in enumerate(exec_lines) if l == NUT_LINE]
    if len(hits) != 1:
        fail(f"nutrient source line: {len(hits)} matches, expected 1")
    out = list(exec_lines)
    out[hits[0]] = new
    return out


def main(target: Path, nut_var: str | None = None) -> None:
    lines = read_lines(target)
    frag = {k: read_lines(CALLSITE / f"phi_mode_{k}.inc") for k in ("decl", "exec")}
    for k in frag:
        while frag[k] and frag[k][-1] == "":
            frag[k].pop()
    if nut_var:
        frag["exec"] = set_nut_source(frag["exec"], nut_var)

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
          f"exec {len(frag['exec'])} lines, split_rates.f copied"
          + (f"\n  nutrient source: CM_NUT = {nut_var}" if nut_var else ""))


if __name__ == "__main__":
    args = sys.argv[1:]
    nut = None
    if "--nut-var" in args:
        i = args.index("--nut-var")
        if i + 1 >= len(args):
            sys.exit(__doc__)
        nut = args[i + 1]
        del args[i:i + 2]
    if len(args) != 1:
        sys.exit(__doc__)
    main(Path(args[0]), nut)
