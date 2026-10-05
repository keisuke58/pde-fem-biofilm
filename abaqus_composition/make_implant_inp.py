"""make_implant_inp.py -- the composition model on an implant: a biofilm layer
around the transmucosal collar of a dental implant, phi, nutrient and
composition solved by Abaqus as on the cube (make_cube_inp.py --case --cons).

    python abaqus_composition/make_implant_inp.py OUT.inp [--nr 6] [--nt 24] [--nz 20]
        [--case 2sp_case6] [--cons 6] [--T 1.1] [--dt 0.025]

Geometry (mm), all values chosen for this work, not from a paper: implant
collar radius 2.05 (a 4.1 mm implant), biofilm layer 0.25 thick, 2.0 high (a
sulcus depth); a 90 degree sector with symmetry planes (x = 0, y = 0).
  - the inner face (r = 2.05) is bonded to the titanium, taken as rigid
    (E_Ti / E_biofilm ~ 1e10): u = 0 there;
  - the nutrient enters through the top face (z = 2.0, the gingival margin):
    c = 1 there; Eq. 35 quasi-static with d = 1 and consumption g = --cons;
  - phi = 1 in the inner element ring (the first colonisers on the titanium),
    as the nodal share of seed elements (make_cube_inp.py --ic fraction);
  - phi has no flux through any face (Eq. 34's natural condition).
Material and growth as on the cube: beta = 0.02 mm^2/T*, k_alpha = 1e-3,
Klempt 2024 stiffness (E = 10 Pa, nu = 0.49), units mm, MPa; the composition
of --case (Klempt et al. 2026) through the point model (prop(28) = 7).
Writes OUT.json with each element's centroid (r, theta, z) and the seed.
"""
import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--nr", type=int, default=6)
    ap.add_argument("--nt", type=int, default=24)
    ap.add_argument("--nz", type=int, default=20)
    ap.add_argument("--ri", type=float, default=2.05)
    ap.add_argument("--thick", type=float, default=0.25)
    ap.add_argument("--height", type=float, default=2.0)
    ap.add_argument("--case", default="2sp_case6")
    ap.add_argument("--nsp", type=int, default=2)
    ap.add_argument("--cons", type=float, default=6.0)
    ap.add_argument("--diff", type=float, default=1.0)
    ap.add_argument("--cref", type=float, default=1.0)
    ap.add_argument("--beta", type=float, default=0.02)
    ap.add_argument("--kalpha", type=float, default=1e-3)
    ap.add_argument("--pen", type=float, default=5.0)
    ap.add_argument("--s", type=float, default=0.15)
    ap.add_argument("--cap", type=float, default=0.9)
    ap.add_argument("--chi0", type=float, default=0.5)
    ap.add_argument("--dt", type=float, default=0.025)
    ap.add_argument("--T", type=float, default=1.1)
    a = ap.parse_args()
    nr, nt, nz = a.nr, a.nt, a.nz
    ro = a.ri + a.thick
    nid = lambda i, j, k: 1 + i + (nr + 1) * j + (nr + 1) * (nt + 1) * k  # noqa: E731
    nodes = []
    for k in range(nz + 1):
        for j in range(nt + 1):
            th = 0.5 * math.pi * j / nt
            for i in range(nr + 1):
                r = a.ri + a.thick * i / nr
                nodes.append((nid(i, j, k), r * math.cos(th), r * math.sin(th), a.height * k / nz))
    elems, seed, cen = [], [], {}
    around = {}
    e = 0
    for k in range(nz):
        for j in range(nt):
            for i in range(nr):
                e += 1
                con = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                       nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                elems.append((e, con))
                cen[e] = (a.ri + a.thick * (i + 0.5) / nr, 0.5 * math.pi * (j + 0.5) / nt,
                          a.height * (k + 0.5) / nz, i, j, k)
                for q in con:
                    around.setdefault(q, [0, 0])[0] += 1
                if i == 0:
                    seed.append(e)
                    for q in con:
                        around[q][1] += 1
    tol = 1e-9
    inner = [q for q, x, y, z in nodes if abs(math.hypot(x, y) - a.ri) < tol]
    top = [q for q, x, y, z in nodes if abs(z - a.height) < tol]
    symx = [q for q, x, y, z in nodes if abs(x) < tol]       # theta = 90 deg: u_x = 0
    symy = [q for q, x, y, z in nodes if abs(y) < tol]       # theta = 0: u_y = 0

    sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat")]
    import material_server as ms
    ms.set_case(a.case)
    p = [0.0] * 47
    p[0], p[6] = 1.0, a.kalpha
    p[7:27] = list(ms.ECOLOGY_CASE["theta"])
    p[27], p[28], p[29], p[30], p[31] = 7.0, 0.01, a.chi0, a.s, a.cap
    p[32] = a.cref
    p[36] = float(a.nsp)
    p[42], p[43], p[44], p[45] = 1e-5, -1e-3, 0.49, 0.3
    p[46] = 1.0                                      # phi = temperature
    p += [1.0]                                       # c from the UEL
    OFF = 1000000
    rows = lambda v: "\n".join(", ".join(f"{x:.17g}" for x in v[i:i + 8]) for i in range(0, len(v), 8))  # noqa: E731
    lst = lambda v: "\n".join(", ".join(str(x) for x in v[i:i + 16]) for i in range(0, len(v), 16))  # noqa: E731
    ic = "\n".join(f"{q}, {around[q][1] / around[q][0]:.10g}" for q in sorted(around) if around[q][1])
    txt = f"""*HEADING
 Biofilm on an implant collar: r {a.ri}-{ro} mm, h {a.height} mm, 90 deg sector, {nr}x{nt}x{nz} C3D8T, {a.case}, nutrient g {a.cons} from the top
*NODE
{chr(10).join(f"{i}, {x:.12g}, {y:.12g}, {z:.12g}" for i, x, y, z in nodes)}
*ELEMENT, TYPE=C3D8T, ELSET=EALL
{chr(10).join(f"{i}, " + ", ".join(map(str, c)) for i, c in elems)}
*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=3, COORDINATES=3, VARIABLES=1, UNSYMM
11, 12
*ELEMENT, TYPE=U1, ELSET=NUTEL
{chr(10).join(f"{i + OFF}, " + ", ".join(map(str, c)) for i, c in elems)}
*UEL PROPERTY, ELSET=NUTEL
{a.diff}, {a.cons}, {OFF}
*NSET, NSET=NALL, GENERATE
1, {len(nodes)}, 1
*NSET, NSET=NINNER
{lst(inner)}
*NSET, NSET=NTOP
{lst(top)}
*NSET, NSET=NSYMX
{lst(symx)}
*NSET, NSET=NSYMY
{lst(symy)}
*ELSET, ELSET=SEED
{lst(seed)}
*SOLID SECTION, ELSET=EALL, MATERIAL=BIOFILM
*MATERIAL, NAME=BIOFILM
*DENSITY
1.0
*USER MATERIAL, CONSTANTS={len(p)}, TYPE=MECHANICAL
{rows(p)}
*USER MATERIAL, CONSTANTS=3, TYPE=THERMAL
{a.beta}, {a.kalpha}, {a.pen}
*DEPVAR
100
*INITIAL CONDITIONS, TYPE=TEMPERATURE
NALL, 0.0
{ic}
*STEP, NLGEOM=YES, INC=100000, UNSYMM=YES
*COUPLED TEMPERATURE-DISPLACEMENT
{a.dt}, {a.T}
*BOUNDARY
NINNER, 1, 3
NSYMX, 1, 1
NSYMY, 2, 2
NTOP, 12, 12, 1.0
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
S
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
SDV84, TEMP, SDV72, SDV73, SDV51
*OUTPUT, FIELD
*ELEMENT OUTPUT
S, SDV
*NODE OUTPUT
NT, U
*END STEP
"""
    Path(a.out).write_text(txt)
    Path(a.out).with_suffix(".json").write_text(json.dumps(
        {"ri": a.ri, "ro": ro, "height": a.height, "nr": nr, "nt": nt, "nz": nz, "seed": seed,
         "cen": {str(e): v for e, v in cen.items()}}))
    print(f"wrote {a.out}: {len(elems)} elements, seed {len(seed)}, inner {len(inner)} nodes, top {len(top)} nodes")


if __name__ == "__main__":
    main()
