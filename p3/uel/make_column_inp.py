"""make_column_inp.py -- Abaqus input for check 2 of driver.py (a column of the unified user element),
so the element can be run inside Abaqus (Keio server or IKMHIWI03) and compared with the driver.

    python make_column_inp.py [n] [outdir]        (default n = 20, outdir = .)

Writes:
  uel_column.inp   n hexahedra along x on [0, 1], cross-section 1/n; user element U1 (dofs 1, 2, 3, 11, 12);
  uel_eco.txt      the point model for ecology_native.f (BIOFILM_ECO_CASE): the CH calibration case of
                   unified_1d.ch_case with b_i = 1, the Day-1 shares scaled to 0.9 as in driver.py.
Run (Linux):
  export BIOFILM_ECO_CASE=$PWD/uel_eco.txt
  python make_abaqus_source.py            (one file: the fixed-form point model + the free-form element)
  abaqus job=uel_column user=unified_uel_abaqus.for cpus=1 interactive
  (cpus=1: the point-model constants are module variables set once)

Steps:
  1. 1e-8 long: phi (dof 11) and c (dof 12) are prescribed at every node (seed 0.9 at x <= 0.25, 1e-4
     elsewhere; c = 1), the only way to set initial values of these dofs in Abaqus;
  2. T = 0.2375 in 60 fixed increments: phi and c free except c = 1 at x = 1; mechanics held fixed, as in
     driver.py check 2. Expected (driver.py, 20 elements): phi(0) 0.776, c(0) 0.936.
Nodal output: NT11 = phi, NT12 = c (*NODE PRINT).

Procedure (15 Oct 2026, A-11): Abaqus does not allow boundary conditions on DOF 11 and 12 in *STATIC, so both
steps are a transient *COUPLED TEMPERATURE-DISPLACEMENT without DELTMX, i.e. with fixed increments (the steady-state
form tried on IKMHIWI03 chose its own increments and gave wrong values). UNSYMM: the element tangent is not symmetric.
The transient procedure needs a specific heat: dummy C3D8T elements on the same nodes with 1e-9 material (overlay).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "JAXFEM"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from driver import column, props_base  # noqa: E402


def eco_file(path):
    import unified_1d as u1
    g0, th = u1.ch_case()
    th = np.array(th, float); th[[3, 4, 8, 9, 15]] = 1.0
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum()
    lines = ["5", "25.0 0.0", " ".join(["1.0"] * 5), " ".join(f"{v:.17g}" for v in th),
             "phi_init " + " ".join(f"{v:.17g}" for v in g[0:5]), "newton 60 1e-26"]
    path.write_text("\n".join(lines) + "\n")



def overlay(elems, offset=100000):
    """dummy coupled elements on the same nodes: a transient *COUPLED TEMPERATURE-DISPLACEMENT needs a specific
    heat, and the model has no material otherwise (IKMHIWI03, third try of A-11). Stiffness, conductivity and
    specific heat 1e-9 of the element's, so they change nothing (relative 1e-9)."""
    out = ["*ELEMENT, TYPE=C3D8T, ELSET=DUMMY"]
    out += [f"{offset + e + 1}, " + ", ".join(str(k + 1) for k in el) for e, el in enumerate(elems)]
    out += ["*SOLID SECTION, ELSET=DUMMY, MATERIAL=DUMMY", "", "*MATERIAL, NAME=DUMMY", "*ELASTIC", "1e-9, 0.3",
            "*CONDUCTIVITY", "1e-9", "*SPECIFIC HEAT", "1e-9", "*DENSITY", "1.0"]
    return out

def deck(n, path):
    nodes, elems, xs = column(n)
    P = props_base(s=1.0)
    seed = nodes[:, 0] <= 0.25 + 1e-12
    top = np.isclose(nodes[:, 0], 1.0)
    out = ["*HEADING", f"unified two-level model, user element, column of {n} elements (p3/uel/driver.py check 2)",
           "*NODE"]
    out += [f"{i + 1}, {x:.10g}, {y:.10g}, {z:.10g}" for i, (x, y, z) in enumerate(nodes)]
    out += ["*USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES=16, I PROPERTIES=3, COORDINATES=3, VARIABLES=128", "1, 2, 3, 11, 12",
            "*ELEMENT, TYPE=U1, ELSET=COL"]
    out += [f"{e + 1}, " + ", ".join(str(k + 1) for k in el) for e, el in enumerate(elems)]
    vals = [f"{v:.10g}" for v in P]
    out += ["*UEL PROPERTY, ELSET=COL", ", ".join(vals[0:8]), ", ".join(vals[8:16]), "1, 0, 0"]   # JPROPS: finite strain, no mean dilatation, nutrient in the reference configuration
    out += overlay(elems)

    def nset(name, mask):
        ids = [str(i + 1) for i in np.where(mask)[0]]
        out.append(f"*NSET, NSET={name}")
        out.extend(", ".join(ids[k:k + 16]) for k in range(0, len(ids), 16))
    nset("ALL", np.ones(len(nodes), bool)); nset("SEED", seed); nset("REST", ~seed); nset("TOP", top)
    nset("X0", np.isclose(nodes[:, 0], 0.0))
    out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", "1e-8, 1e-8",
            "*BOUNDARY", "ALL, 1, 3, 0.0", "SEED, 11, 11, 0.9", "REST, 11, 11, 1e-4", "ALL, 12, 12, 1.0",
            "*END STEP",
            "*STEP, NAME=GROW, INC=1000", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{0.2375 / 60:.10g}, 0.2375",
            "*BOUNDARY, OP=NEW", "ALL, 1, 3, 0.0", "TOP, 12, 12, 1.0",
            "*NODE PRINT, NSET=X0, FREQUENCY=10", "NT", "*NODE PRINT, NSET=ALL, FREQUENCY=60", "NT",
            "*END STEP"]
    path.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    od = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
    od.mkdir(parents=True, exist_ok=True)
    deck(n, od / "uel_column.inp"); eco_file(od / "uel_eco.txt")
    print(f"wrote {od / 'uel_column.inp'} and {od / 'uel_eco.txt'}")
