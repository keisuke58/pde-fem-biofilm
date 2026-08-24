#!/usr/bin/env python3
"""composition_to_material.py — CLSM composition phi -> per-integration-point
material constants (C10, C01, D1, eta) for the USERMAT's kStateMat=1 path.

This is the "stiffness E(phi)" leg of the model (RESEARCH_MODEL.md sec.3), the
one that runs *alongside* the growth field alpha rather than through it. Until
now the USERMAT took C10/C01/D1/eta as constants in prop(1:4), so every Gauss
point had the same stiffness and the only thing distinguishing the four
clinical conditions was alpha. That leaves out the largest mechanical
difference in the whole study: E spans roughly 995 Pa (commensal) to 32 Pa
(dysbiotic), a ~30x contrast.

Composition is CLSM-*measured* input, not something the solve evolves, so
these constants are known before the FE run starts. They are therefore
computed once here and delivered as initial state (ustatev(11:14)) -- there is
deliberately no per-increment Python call on this path. The socket bridge
(material_server.py / kUsePy=1) solves a different problem: swapping the
constitutive *law*, not its coefficients.

    python composition_to_material.py --phi 0.2,0.2,0.2,0.2,0.2
    python composition_to_material.py --E 995 --di 0.05 --matid 1

Chain:  phi -> E(phi)          material_models.compute_E_composite
             -> DI(phi)        material_models.compute_di
             -> C10,C01,D1     material_models.compute_mooney_rivlin_params
             -> eta(DI)        material_models.compute_viscosity_di
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import material_models as mm  # noqa: E402

# The USERMAT reads these four from ustatev(11:14), in this order.
STATE_SLOTS = ("C10", "C01", "D1", "eta")
STATE_BASE = 11                      # ustatev index of the first of them
NSTATEV_WITH_MATERIAL = 14           # what TB,STATE must be sized to


def material_constants_from_E(E, di, nu=0.30,
                              c01_ratio=mm.C01_RATIO_DEFAULT):
    """(E [Pa], DI) -> {C10, C01, D1, eta} in USERMAT units (Pa, 1/Pa, Pa*s).

    Kept separate from the phi entry point so a condition-level E that was
    calibrated directly (the E_map values in run_ve_twin_experiment.py, say)
    can be used without round-tripping through a composition that may not
    reproduce it.
    """
    E = np.asarray(E, dtype=np.float64)
    di = np.asarray(di, dtype=np.float64)
    mr = mm.compute_mooney_rivlin_params(E, nu=nu, c01_ratio=c01_ratio)
    return {
        "C10": mr["C10"],
        "C01": mr["C01"],
        "D1": mr["D1"],
        "eta": mm.compute_viscosity_di(di),
        "E": E,
        "DI": di,
    }


def material_constants_from_composition(phi_species, nu=0.30,
                                        c01_ratio=mm.C01_RATIO_DEFAULT):
    """(N,5) species volume fractions -> {C10, C01, D1, eta, E, DI}.

    phi rows need not sum to 1; every downstream model normalises internally.
    """
    phi = np.asarray(phi_species, dtype=np.float64)
    if phi.shape[-1] != 5:
        raise ValueError(
            f"expected 5 species in the last axis, got shape {phi.shape}. "
            "Order is [S. oralis, A. naeslundii, Veillonella, F. nucleatum, "
            "P. gingivalis] (RESEARCH_MODEL.md sec.1).")
    return material_constants_from_E(
        mm.compute_E_composite(phi), mm.compute_di(phi),
        nu=nu, c01_ratio=c01_ratio)


def _tbdata_lines(start_index: int, values) -> list[str]:
    """TBDATA in chunks of at most 6 values.

    ANSYS silently drops anything past the 6th value in a single TBDATA call.
    That has already cost this project one debugging session (see
    apdl/t_growth_free.dat and THESIS_PLAYBOOK.md), so the chunking lives here
    rather than in every caller.
    """
    vals = list(values)
    lines = []
    for off in range(0, len(vals), 6):
        chunk = vals[off:off + 6]
        body = ",".join(f"{v:.9G}" for v in chunk)
        lines.append(f"TBDATA,{start_index + off},{body}")
    return lines


def apdl_state_block(const: dict, matid: int = 1, alpha: float = 0.0,
                     index: int | None = None) -> str:
    """Emit the TB,STATE block initialising one material's integration points.

    `const` is a result dict from either entry point above; `index` picks a row
    when the arrays are vectorised. ANSYS applies TB,STATE per *material*, so
    spatially varying composition means one material per composition bin --
    call this once per bin and give each its own matid.
    """
    def pick(key):
        arr = np.atleast_1d(np.asarray(const[key], dtype=np.float64))
        return float(arr[0] if index is None else arr[index])

    c10, c01, d1, eta = (pick(k) for k in STATE_SLOTS)
    mtype = 1.0 if c01 > 0.0 else 0.0        # prop(5): C01=0 degenerates to Neo-Hookean
    lines = [
        f"! --- material {matid}: E={pick('E'):.6G} Pa, DI={pick('DI'):.4f} ---",
        f"TB,USER,{matid},1,7",
        f"TBDATA,1,{c10:.9G},{c01:.9G},{d1:.9G},{eta:.9G},{mtype:.1f},0.0",
        "TBDATA,7,1.0                             ! prop(7)=kStateMat -> read ustatev(11:14)",
        f"TB,STATE,{matid},,{NSTATEV_WITH_MATERIAL}",
        "TBDATA,1,1.0,0.0,0.0,0.0,1.0,0.0         ! Fv(1,1..3), Fv(2,1..3)",
        "TBDATA,7,0.0,0.0,1.0                     ! Fv(3,1..3)",
        f"TBDATA,10,{alpha:.9G}   ! alpha (growth driver)",
    ]
    lines += [
        _tbdata_lines(STATE_BASE, (c10, c01, d1, eta))[0]
        + "   ! C10,C01,D1,eta per IP"
    ]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--phi", help="5 comma-separated species fractions")
    g.add_argument("--E", type=float, help="Young's modulus [Pa] directly")
    ap.add_argument("--di", type=float, default=None,
                    help="dysbiosis index (required with --E)")
    ap.add_argument("--nu", type=float, default=0.30)
    ap.add_argument("--c01-ratio", type=float, default=mm.C01_RATIO_DEFAULT)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--matid", type=int, default=1)
    ap.add_argument("--apdl", action="store_true", help="emit the TB,STATE block")
    a = ap.parse_args(argv)

    if a.phi is not None:
        phi = np.array([float(x) for x in a.phi.split(",")])
        const = material_constants_from_composition(
            phi, nu=a.nu, c01_ratio=a.c01_ratio)
    else:
        if a.di is None:
            ap.error("--di is required with --E (viscosity is a function of DI)")
        const = material_constants_from_E(
            a.E, a.di, nu=a.nu, c01_ratio=a.c01_ratio)

    if a.apdl:
        print(apdl_state_block(const, matid=a.matid, alpha=a.alpha))
    else:
        for k in ("E", "DI", *STATE_SLOTS):
            print(f"{k:>5s} = {float(np.atleast_1d(const[k])[0]):.6G}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
