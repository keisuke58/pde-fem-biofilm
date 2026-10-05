"""make_inp.py -- one-element Abaqus check of the composition UMAT (make_umat.py).

    python abaqus_composition/make_inp.py OUT.inp [--case 2sp_case6] [--phi 1.0]
        [--c 1.0] [--cref 0] [--s 0.15] [--cap 0.9] [--kalpha 1e-3] [--dt 0.025]
        [--T 1.0] [--bc clamped|free] [--nsp 2]

A unit cube of one C3D8 element (8 integration points), phi and c held
uniform through field variables 1 and 2. --bc clamped fixes every node, so
F = I and the stress has the closed form of fully constrained growth (with the
pressure term of DEVIATOR_SCALING_FINDING.md); --bc free holds only rigid-body
motion, so uniform growth gives no stress. Material constants: Klempt 2024
(E = 10 Pa as 1e-5 MPa, nu = 0.49; E_void = -1e-3 selects that branch), the
point model's theta from the material server's case, the rest as the ANSYS
composition runs (prop(29) = 0.01, prop(30) = 0.5).
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat")]


def props(a):
    import material_server as ms
    ms.set_case(a.case)
    p = [0.0] * 46
    p[0] = 1.0
    p[6] = a.kalpha
    p[7:27] = list(ms.ECOLOGY_CASE["theta"])
    p[27] = 7.0
    p[28], p[29], p[30], p[31], p[32] = 0.01, 0.5, a.s, a.cap, a.cref
    p[36] = float(a.nsp)
    p[42], p[43], p[44], p[45] = 1e-5, -1e-3, 0.49, 0.3
    return p


def inp(a):
    p = props(a)
    rows = [", ".join(f"{x:.17g}" for x in p[i:i + 8]) for i in range(0, len(p), 8)]
    bc = ("NALL, 1, 3\n" if a.bc == "clamped" else "1, 1, 3\n2, 2, 3\n4, 1, 1\n4, 3, 3\n5, 1, 2\n")
    return f"""*HEADING
 One-element check of the composition UMAT ({a.case}, phi {a.phi}, c {a.c}, {a.bc})
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=NALL, GENERATE
1, 8, 1
*SOLID SECTION, ELSET=EALL, MATERIAL=BIOFILM
*MATERIAL, NAME=BIOFILM
*USER MATERIAL, CONSTANTS={len(p)}
{chr(10).join(rows)}
*DEPVAR
100
*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1
NALL, {a.phi}
*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=2
NALL, {a.c}
*STEP, NLGEOM=YES, INC=100000
*STATIC, DIRECT
{a.dt}, {a.T}
*BOUNDARY
{bc}*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=INTEGRATION POINTS, SUMMARY=NO, TOTALS=NO
S
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=INTEGRATION POINTS, SUMMARY=NO, TOTALS=NO
SDV72, SDV73, SDV74, SDV75, SDV84
*OUTPUT, FIELD
*ELEMENT OUTPUT
S, SDV
*END STEP
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--case", default="2sp_case6")
    ap.add_argument("--phi", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=1.0)
    ap.add_argument("--cref", type=float, default=0.0)
    ap.add_argument("--s", type=float, default=0.15)
    ap.add_argument("--cap", type=float, default=0.9)
    ap.add_argument("--kalpha", type=float, default=1e-3)
    ap.add_argument("--dt", type=float, default=0.025)
    ap.add_argument("--T", type=float, default=1.0)
    ap.add_argument("--bc", choices=("clamped", "free"), default="clamped")
    ap.add_argument("--nsp", type=int, default=2)
    a = ap.parse_args()
    Path(a.out).write_text(inp(a))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
