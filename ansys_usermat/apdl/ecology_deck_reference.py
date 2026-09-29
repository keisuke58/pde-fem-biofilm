#!/usr/bin/env python3
"""Independent EXPECTED values for the single-point ecology decks
(t_growth_ecology*.dat), read straight from each deck.

Every node of these decks is fixed, so F = I and the result follows from the
ecology alone: parse TB,USER / TB,STATE of material 1 and TIME, advance the
ecology state over TIME in ceil(TIME / 1e-4) sub-steps (usermat_biofilm.f's
DT_ECO_MAX), accumulate alpha += k_alpha * phi_int, and evaluate the stress
with material_server.stress_core. An all-zero ecology state is seeded with
ecology_jax.default_initial_state(), as INIT_ECO_IF_ZERO does.

Uses the constants in ecology_constants.py (c* = 25, no Hill gate,
alpha* = 0 since 2026-09-29).

    python ansys_usermat/apdl/ecology_deck_reference.py [deck.dat ...]
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "ansys_usermat" / "coupling"))

import ecology_jax  # noqa: E402
from material_server import stress_core  # noqa: E402

DT_ECO_MAX = 1.0e-4
DECKS = ["t_growth_ecology.dat", "t_growth_ecology_clsm_phi.dat",
         "t_growth_ecology_multi.dat", "t_growth_ecology_substep.dat"]


def table(lines, kind, mat=1):
    """Values of the first TB,<kind>,<mat> block, by TBDATA start index."""
    vals, inside = {}, False
    for ln in lines:
        code = ln.split("!")[0].strip()
        if not code:
            continue
        head = [t.strip().upper() for t in code.split(",")]
        if head[0] == "TB":
            if inside:
                break
            inside = head[1] == kind and head[2] == str(mat)
            continue
        if inside and head[0] == "TBDATA":
            start = int(head[1])
            for k, v in enumerate(head[2:]):
                if v:
                    vals[start + k] = float(v)
        elif inside:
            break
    n = max(vals) if vals else 0
    return np.array([vals.get(i, 0.0) for i in range(1, n + 1)])


def reference(deck):
    lines = Path(deck).read_text(errors="replace").splitlines()
    prop = table(lines, "USER")
    state = np.zeros(26)
    s = table(lines, "STATE")
    state[:len(s)] = s
    t = [float(m.group(1)) for ln in lines
         if (m := re.match(r"\s*TIME\s*,\s*([-+0-9.Ee]+)", ln))][0]
    C10, C01, D1, eta, mtype = prop[0:5]
    k_alpha, theta = prop[8], prop[9:29]

    g0 = state[14:26]
    if not np.any(g0):
        g0 = np.asarray(ecology_jax.default_initial_state())
    n_sub = math.ceil(t / DT_ECO_MAX) if t > DT_ECO_MAX else 1
    g, phi_int = ecology_jax.ecology_substeps(g0, theta, t, n_sub)
    alpha = state[9] + k_alpha * float(phi_int)
    Fv = state[0:9].reshape(3, 3)
    sv, _, _ = stress_core(np.eye(3), Fv, alpha, C10, C01, D1, eta, mtype, t)
    return {"TIME": t, "n_sub": n_sub, "alpha": alpha, "g": np.asarray(g),
            "S": sv}


def main(argv):
    decks = argv or [str(_HERE / d) for d in DECKS]
    for d in decks:
        r = reference(d)
        g = r["g"]
        print(f"--- {Path(d).name}  (TIME={r['TIME']:g}, n_sub={r['n_sub']})")
        print("  g_new(phi)   = [" + ", ".join(f"{v:.9g}" for v in g[0:5]) + "]")
        print(f"  g_new(phi0)  = {g[5]:.9g}")
        print("  g_new(psi)   = [" + ", ".join(f"{v:.9g}" for v in g[6:11]) + "]")
        print(f"  g_new(gamma) = {g[11]:.9g}")
        print(f"  alpha_new    = {r['alpha']:.9e}")
        print(f"  SX=SY=SZ     = {r['S'][0]:.8e}   (shear max "
              f"{np.abs(r['S'][3:]).max():.1e})")


if __name__ == "__main__":
    main(sys.argv[1:])
