"""Material-override wrapper around build_assembly.py for the OAT sensitivity sweep.

Same CLI as build_assembly.py (imp_cache, job_name, [crown_E]) but optionally honours a
``MATS_OVERRIDE`` env var containing a JSON dict ``{mat_name: [E_MPa, nu]}``. The override is applied
to ``build_assembly.MATS`` AFTER all module-level definitions but BEFORE its main path executes, so
the headline behaviour with the env var UNSET is byte-identical to ``python build_assembly.py``.

The headline build_assembly.py is preserved untouched (per the project convention of new variants
as new files); the sensitivity sweep run_crown_sensitivity.sh calls THIS wrapper so the headline
.inp generator is never modified.

Usage (run from /home/nishioka/IKM_Hiwi/FEM/tier2b_real):

    MATS_OVERRIDE='{"GINGIVA": [1.5, 0.45]}' python build_assembly_override.py \
        cache_implant.npz tier2b_crown_sens_gingiva_low

When ``MATS_OVERRIDE`` is unset the script is equivalent to plain build_assembly.py.
"""
from __future__ import annotations

import json
import os
import runpy
import sys

import build_assembly


def main() -> int:
    override = os.environ.get("MATS_OVERRIDE")
    if override:
        try:
            data = json.loads(override)
        except json.JSONDecodeError as e:
            sys.exit(f"[build_assembly_override] invalid MATS_OVERRIDE JSON: {e}")
        for mat, val in data.items():
            if mat not in build_assembly.MATS:
                print(f"[build_assembly_override] warning: unknown material '{mat}', "
                      "adding (may be ignored by the assembler).")
            build_assembly.MATS[mat] = (float(val[0]), float(val[1]))
        print(f"[build_assembly_override] applied {len(data)} material override(s): {data}")

    # Re-enter build_assembly's CLI path as if it had been invoked directly. runpy preserves
    # sys.argv, __name__ == "__main__" semantics, and the module-level MATS dict we just mutated.
    runpy.run_module("build_assembly", run_name="__main__", alter_sys=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
