#!/usr/bin/env python3
"""patch_usermat_to_v222.py — retarget Oliver's USERMAT from ANSYS 2024 R2 to v222.

Oliver's UPF pool is written against the **2024 R2** `usermat` interface, which
takes 41 arguments. IKMHIWI03 has **v222**, whose interface takes 42: the two
reserved slots `var1`/`var2` are named `pVolDer(3)` and `hrmflg` in 2024 R2,
and `var8` was dropped there. See `V222_PORT_INSTRUCTIONS.md` §1.1.

This script performs that retarget mechanically, so the edit is reproducible
when Oliver sends a newer version rather than something to redo by hand.

    python patch_usermat_to_v222.py Usermat_P21-V21_Conection_Test.F -o Usermat_v222.F

Four changes, all verified necessary against the delivered source:

  1. the subroutine argument list  -> the v222 form;
  2. drop the `pVolDer (3),` entry from the DOUBLE PRECISION block;
  3. drop the standalone `DOUBLE PRECISION hrmflg`;
  4. **drop `data var1/0.0d0/` and `data var2/0.0d0/`.**

(4) is the one that is easy to miss and does not fail obviously. Under 2024 R2
`var1`/`var2` are not arguments, so they are locals and initialising them with
DATA is legal. Under v222 they *become dummy arguments*, and a DATA statement
for a dummy argument is a hard compile error. Leaving them in turns a working
port into an error message that points at the wrong thing.

`pVolDer` and `hrmflg` are only declared, never used in the routine body — the
script checks that and refuses to proceed if a future version starts using
them, since dropping them would then change behaviour.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# The v222 argument list, taken from ansys_usermat/usermat_biofilm.f — the
# version verified in-solver on this machine's ANSYS.
V222_ARGS = """      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, var1, var2, var3, var4,
     &                   var5, var6, var7, var8)
"""


def _body_after_declarations(src: str) -> str:
    """Everything after the declaration block, for the usage check."""
    m = re.search(r"^\s*DOUBLE PRECISION\s+var0.*?$", src, re.M | re.I)
    return src[m.end():] if m else src


def patch(src: str) -> tuple[str, list[str]]:
    """Return (patched source, list of changes applied)."""
    changes: list[str] = []

    body = _body_after_declarations(src)
    used = [n for n in ("pVolDer", "hrmflg") if re.search(rf"\b{n}\b", body)]
    if used:
        raise SystemExit(
            f"refusing to patch: {', '.join(used)} is USED in the routine body, "
            "not merely declared. Dropping it would change behaviour — port it "
            "deliberately instead."
        )

    # 1. the argument list
    new, n = re.subn(r"      subroutine\s+usermat\s*\(.*?\)\s*\n",
                     V222_ARGS, src, count=1, flags=re.S | re.I)
    if not n:
        raise SystemExit("could not find the usermat subroutine statement")
    src = new
    changes.append("argument list -> v222 (42 args: var1, var2, ..., var8)")

    # 2. pVolDer out of the DOUBLE PRECISION block
    src, n = re.subn(r"^\s*&\s*pVolDer\s*\(\s*3\s*\)\s*,\s*\n", "", src,
                     count=1, flags=re.M | re.I)
    if n:
        changes.append("removed pVolDer(3) from the DOUBLE PRECISION block")

    # 3. the standalone hrmflg declaration
    src, n = re.subn(r"^\s*DOUBLE PRECISION\s+hrmflg\s*\n", "", src,
                     count=1, flags=re.M | re.I)
    if n:
        changes.append("removed the DOUBLE PRECISION hrmflg declaration")

    # 4. var8 must exist; var1/var2 must lose their DATA initialisers
    src, n = re.subn(r"(DOUBLE PRECISION\s+var0.*?var6\s*,\s*var7)(?!\s*,\s*var8)",
                     r"\1, var8", src, count=1, flags=re.S | re.I)
    if n:
        changes.append("added var8 to the var declaration")

    for v in ("var1", "var2"):
        src, n = re.subn(rf"^\s*data\s+{v}\s*/[^/]*/\s*\n", "", src,
                         count=1, flags=re.M | re.I)
        if n:
            changes.append(f"removed 'data {v}/../' "
                           f"({v} is a dummy argument under v222; DATA on a "
                           f"dummy argument is a compile error)")

    return src, changes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path, help="Oliver's Usermat_*.F")
    ap.add_argument("-o", "--out", type=Path, required=True, help="patched output")
    a = ap.parse_args(argv)

    patched, changes = patch(a.source.read_text(errors="replace"))
    a.out.write_text(patched)

    print(f"{a.source.name} -> {a.out.name}")
    for c in changes:
        print(f"  - {c}")
    if len(changes) < 5:
        print("\n  NOTE: fewer changes than expected. The source may already be "
              "patched, or its layout may have moved — check the output before "
              "building.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
