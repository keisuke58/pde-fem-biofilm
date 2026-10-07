"""make_implant_inp.py -- the composition model on an implant: a biofilm layer
around the transmucosal collar of a dental implant, phi, nutrient and
composition solved by Abaqus as on the cube (make_cube_inp.py --case --cons).

    python abaqus_composition/make_implant_inp.py OUT.inp [--nr 6] [--nt 24] [--nz 20]
        [--case 2sp_case6] [--cons 6] [--T 1.1] [--dt 0.025] [--gw 0]

With zero-order consumption c goes negative where g phi H^2 / (2 d) > 1
(H the depth below the margin): --cons 6 on 2 mm gives c down to -0.9 (the
point model then counts the point as starved); --cons 1 (Klempt 2024 Table 2,
as the ANSYS runs) keeps c > 0.6.

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

Options added 6 Oct 2026:
  --front r [--blend w] [--kmono k] [--hstab h] [--eps e]: the front term of
    Eq. 34 with the nutrient of the UEL, as make_cube_inp.py (without
    --blend: Eq. 34 as printed; with --blend: growth on every face, my
    modification);
  --bulge b --nut outer: a tooth instead of the implant. The inner radius
    grows from --ri at the base to --ri + b at the top, as the crown of a
    tooth widens from the cervical line to the height of contour, and the
    nutrient (saliva) enters through the outer face of the layer instead of
    the top. The tooth's enamel is taken as rigid (E ~ 80 GPa). Values are
    mine, not from a paper.
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
    ap.add_argument("--gw", type=float, default=0.0,
                    help="prop(36) = d: species-weighted growth, alpha_dot = k_alpha phi (1 + d (2 chi_1 - 1)) "
                         "(not from a paper; 0 = Eq. 36)")
    ap.add_argument("--front", type=float, default=None, help="r (mm/T*) of the front term of Eq. 34")
    ap.add_argument("--kmono", type=float, default=1.0)
    ap.add_argument("--hstab", type=float, default=0.0)
    ap.add_argument("--eps", type=float, default=0.0)
    ap.add_argument("--blend", type=float, default=None)
    ap.add_argument("--bulge", type=float, default=0.0, help="inner radius ri + bulge sin(pi z / 2h): a tooth crown")
    ap.add_argument("--nut", choices=("top", "outer"), default="top")
    ap.add_argument("--first-order", action="store_true",
                    help="consumption g phi c (UEL property 4 = 1) instead of the paper's g phi")
    ap.add_argument("--ph", type=float, default=None,
                    help="P_h of the homeostatic growth law (constants 49-51): p_h = P_h E k_alpha T*")
    ap.add_argument("--ph-local", action="store_true",
                    help="p_h = P_h E (phi^2 + f) k_alpha T* (local stiffness, f = 1e-3)")
    ap.add_argument("--base", choices=("free", "uz", "clamped"), default="free",
                    help="bottom face z = 0: free (default, as before), u_z = 0, or u = 0")
    a = ap.parse_args()
    nr, nt, nz = a.nr, a.nt, a.nz
    ro = a.ri + a.thick
    rin = lambda z: a.ri + a.bulge * math.sin(0.5 * math.pi * z / a.height)  # noqa: E731
    nid = lambda i, j, k: 1 + i + (nr + 1) * j + (nr + 1) * (nt + 1) * k  # noqa: E731
    nodes, ii = [], {}
    for k in range(nz + 1):
        z = a.height * k / nz
        for j in range(nt + 1):
            th = 0.5 * math.pi * j / nt
            for i in range(nr + 1):
                r = rin(z) + a.thick * i / nr
                nodes.append((nid(i, j, k), r * math.cos(th), r * math.sin(th), z))
                ii[nid(i, j, k)] = i
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
                zc = a.height * (k + 0.5) / nz
                cen[e] = (rin(zc) + a.thick * (i + 0.5) / nr, 0.5 * math.pi * (j + 0.5) / nt, zc, i, j, k)
                for q in con:
                    around.setdefault(q, [0, 0])[0] += 1
                if i == 0:
                    seed.append(e)
                    for q in con:
                        around[q][1] += 1
    tol = 1e-9
    inner = [q for q, x, y, z in nodes if ii[q] == 0]
    if a.nut == "top":
        top = [q for q, x, y, z in nodes if abs(z - a.height) < tol]
    else:
        top = [q for q, x, y, z in nodes if ii[q] == nr]
    symx = [q for q, x, y, z in nodes if abs(x) < tol]       # theta = 90 deg: u_x = 0
    symy = [q for q, x, y, z in nodes if abs(y) < tol]       # theta = 0: u_y = 0
    base = [q for q, x, y, z in nodes if abs(z) < tol]       # bottom face (--base)

    sys.path[:0] = [str(ROOT / "ansys_usermat" / "coupling"), str(ROOT / "ansys_usermat")]
    import material_server as ms
    ms.set_case(a.case)
    p = [0.0] * 47
    p[0], p[6] = 1.0, a.kalpha
    p[7:27] = list(ms.ECOLOGY_CASE["theta"])
    p[27], p[28], p[29], p[30], p[31] = 7.0, 0.01, a.chi0, a.s, a.cap
    p[32] = a.cref
    p[35] = a.gw
    p[36] = float(a.nsp)
    p[42], p[43], p[44], p[45] = 1e-5, -1e-3, 0.49, 0.3
    p[46] = 1.0                                      # phi = temperature
    p += [1.0]                                       # c from the UEL
    if a.ph is not None:
        # constants 49-51: homeostatic-pressure growth law (phi_mode_exec.inc);
        # p_ref = P_h E k_alpha T* with T* = 1 and E = constant 43 (keio_wp2)
        p += [a.ph * p[42] * a.kalpha, 1.0 if a.ph_local else 0.0, -p[43]]
    OFF = 1000000
    uprop = [a.diff, a.cons, OFF] + ([1.0] if a.first_order else [])
    bc = {"free": "", "uz": "NBASE, 3, 3\n", "clamped": "NBASE, 1, 3\n"}[a.base]
    rows = lambda v: "\n".join(", ".join(f"{x:.17g}" for x in v[i:i + 8]) for i in range(0, len(v), 8))  # noqa: E731
    lst = lambda v: "\n".join(", ".join(str(x) for x in v[i:i + 16]) for i in range(0, len(v), 16))  # noqa: E731
    ic = "\n".join(f"{q}, {around[q][1] / around[q][0]:.10g}" for q in sorted(around) if around[q][1])
    therm = [a.beta, a.kalpha, a.pen]
    if a.front is not None:
        therm += [a.front, a.kmono, a.hstab, a.eps] + ([a.blend] if a.blend is not None else [])
    what = "a tooth's crown" if a.bulge or a.nut == "outer" else "an implant collar"
    txt = f"""*HEADING
 Biofilm on {what}: r {a.ri}(+{a.bulge})-{ro} mm, h {a.height} mm, 90 deg sector, {nr}x{nt}x{nz} C3D8T, {a.case}, nutrient g {a.cons}{" phi c" if a.first_order else " phi"} from the {a.nut} face, front {a.front} blend {a.blend}, base {a.base}, P_h {a.ph}{" local" if a.ph_local else ""}
*NODE
{chr(10).join(f"{i}, {x:.12g}, {y:.12g}, {z:.12g}" for i, x, y, z in nodes)}
*ELEMENT, TYPE=C3D8T, ELSET=EALL
{chr(10).join(f"{i}, " + ", ".join(map(str, c)) for i, c in elems)}
*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES={len(uprop)}, COORDINATES=3, VARIABLES=1, UNSYMM
11, 12
*ELEMENT, TYPE=U1, ELSET=NUTEL
{chr(10).join(f"{i + OFF}, " + ", ".join(map(str, c)) for i, c in elems)}
*UEL PROPERTY, ELSET=NUTEL
{", ".join(str(x) for x in uprop)}
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
*NSET, NSET=NBASE
{lst(base)}
*ELSET, ELSET=SEED
{lst(seed)}
*SOLID SECTION, ELSET=EALL, MATERIAL=BIOFILM
*MATERIAL, NAME=BIOFILM
*DENSITY
1.0
*USER MATERIAL, CONSTANTS={len(p)}, TYPE=MECHANICAL
{rows(p)}
*USER MATERIAL, CONSTANTS={len(therm)}, TYPE=THERMAL
{", ".join(f"{x:.10g}" for x in therm)}
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
{bc}NTOP, 12, 12, 1.0
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
S
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
SDV84, TEMP, SDV72, SDV73, SDV51{", SDV95, SDV96" if a.ph is not None else ""}
*OUTPUT, FIELD
*ELEMENT OUTPUT
S, SDV
*NODE OUTPUT
NT, U
*END STEP
"""
    Path(a.out).write_text(txt)
    Path(a.out).with_suffix(".json").write_text(json.dumps(
        {"ri": a.ri, "ro": ro, "bulge": a.bulge, "nut": a.nut, "height": a.height, "nr": nr, "nt": nt, "nz": nz, "seed": seed,
         "cen": {str(e): v for e, v in cen.items()}}))
    print(f"wrote {a.out}: {len(elems)} elements, seed {len(seed)}, inner {len(inner)} nodes, top {len(top)} nodes")


if __name__ == "__main__":
    main()
