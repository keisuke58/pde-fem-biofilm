"""make_klempt_inp.py -- Klempt et al. 2024's own test cases in Abaqus, with the
front term of Eq. 34 as printed (the partner's element cannot run it), to be
compared with the independent finite-difference reproduction
JAXFEM/klempt2024_quantitative.py (growth="printed", zero-order consumption).

    python abaqus_composition/make_klempt_inp.py OUT.inp [--case fig7_high|fig4_edge]
        [--n 20] [--dt 0.002] [--T 0.2] [--pen 100] [--every 10] [--eps 0]

The paper's 20 um cube in um (n^3 C3D8T, nodes on the reproduction's grid
for n = 20), Table 2 values: beta = 2 um^2/T*, k_alpha = 1e-3, r = 100 um/T*,
k = 1, g/d = 1e8/1e10 (only the ratio enters the quasi-static Eq. 35).
fig7_high: phi = 1 at the nodes 1 um above the bottom within 2.5 um of the
axis, c = 1 on the bottom face. fig4_edge: phi = 1 at the nodes within 5 um of
the centre, c = 1 on the edge x = y = 20 um. The reproduction clips phi to
[0, 1]; here a stiff penalty (P = 100) does that. Mechanics as in the
one-species runs (it does not act back on phi). Prints the centroid phi
(TEMP) and c (SDV51) every --every increments. --eps (1/um): regularised
n_c = grad c / sqrt(|grad c|^2 + eps^2), as eps in the reproduction's run().
"""
import argparse
import numpy as np
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--case", choices=("fig7_high", "fig7_low", "fig4_edge", "fig4_corner", "advect"), default="fig7_high")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--dt", type=float, default=0.002)
    ap.add_argument("--T", type=float, default=0.2)
    ap.add_argument("--pen", type=float, default=100.0)
    ap.add_argument("--every", type=int, default=10)
    ap.add_argument("--no-front", action="store_true")
    ap.add_argument("--eps", type=float, default=0.0)
    ap.add_argument("--seed-field", action="store_true",
                    help="initial phi = the paper-grid seed as its trilinear interpolant (seed_field in the "
                         "reproduction): the same continuous seed on any n (grid studies)")
    ap.add_argument("--blend", type=float, default=None,
                    help="w: growth on every face, r c/(k+c) (w |grad phi| + (1-w) |n_c . grad phi|), in place of "
                         "the printed front term (KLEMPT2024_REPRODUCTION.md sec. 15; a hypothesis)")
    ap.add_argument("--first-order", action="store_true", help="consumption g phi c (not the paper's g phi)")
    ap.add_argument("--felix", action="store_true",
                    help="front term as in the first author's current implementation: added only where positive, "
                         "times (1 - phi) below phi = 1 (UMATHT constant 9; ansys_usermat/FELIX_VS_PARTNER_DIFF.md)")
    ap.add_argument("--scale", type=float, default=1.0, help="time scale s: beta, k_alpha and r times s (sec. 12)")
    ap.add_argument("--E", type=float, default=1e-5, help="YOUNG_BIO (default 1e-5: E = 10 Pa in MPa)")
    ap.add_argument("--s-every", type=int, default=0, help="also print S and SDV84 every this many increments")
    ap.add_argument("--hstab", type=float, default=None,
                    help="length in the artificial diffusion v hstab/2 of the front term (default: the element "
                         "size; 0 = none: with the growth form --blend it slowed the front, 4.2 high)")
    a = ap.parse_args()
    L, n = 20.0, a.n
    h = L / n
    nid = lambda i, j, k: 1 + i + (n + 1) * j + (n + 1) ** 2 * k  # noqa: E731
    nodes = [(nid(i, j, k), h * i, h * j, h * k)
             for k in range(n + 1) for j in range(n + 1) for i in range(n + 1)]
    elems = []
    e = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                e += 1
                elems.append((e, [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                                  nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                                  nid(i, j + 1, k + 1)]))
    tol = 1e-9
    top = []
    if a.case == "advect":                           # check of the front term: c = 1 / 0 on bottom / top,
        seed = [q for q, x, y, z in nodes if (x - 10) ** 2 + (y - 10) ** 2 + (z - 10) ** 2 <= 9 + tol]
        src = [q for q, x, y, z in nodes if abs(z) < tol]          # no consumption: c linear in z,
        top = [q for q, x, y, z in nodes if abs(z - L) < tol]      # grad c the same everywhere
    elif a.case == "fig4_corner":                    # sec. 14: the nutrient in a 2 um block at the corner
        seed = [q for q, x, y, z in nodes if (x - 10) ** 2 + (y - 10) ** 2 + (z - 10) ** 2 <= 25 + tol]
        src = [q for q, x, y, z in nodes if min(x, y, z) >= L - 2 - tol]
    elif a.case == "fig4_edge":
        seed = [q for q, x, y, z in nodes if (x - 10) ** 2 + (y - 10) ** 2 + (z - 10) ** 2 <= 25 + tol]
        src = [q for q, x, y, z in nodes if abs(x - L) < tol and abs(y - L) < tol]
    else:
        seed = [q for q, x, y, z in nodes if abs(z - h) < tol and (x - 10) ** 2 + (y - 10) ** 2 <= 6.25 + tol]
        src = [q for q, x, y, z in nodes if abs(z) < tol]
    ic = "NSEED, 1.0"
    if a.seed_field and a.case != "advect":
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "JAXFEM"))
        import klempt2024_quantitative as kq
        val = kq.seed_field(a.case.replace("fig4_corner", "fig4_edge"),*(np.array([q[i] for q in nodes]) for i in (1, 2, 3)))
        ic = chr(10).join(f"{q[0]}, {v:.10g}" for q, v in zip(nodes, val) if v > 1e-12)
    g = 1e10 if a.case == "fig7_low" else (0.0 if a.case == "advect" else 1e8)
    lst = lambda v: "\n".join(", ".join(str(x) for x in v[i:i + 16]) for i in range(0, len(v), 16))  # noqa: E731
    p = [0.0] * 47
    s = a.scale
    p[0], p[6], p[27] = 1.0, 1e-3 * s, 1.0
    p[42], p[43], p[44], p[45] = a.E, -1e-3, 0.49, 0.3
    p[46] = 1.0
    p += [1.0]                                       # constant 48: c from the UEL
    rows = "\n".join(", ".join(f"{x:.17g}" for x in p[i:i + 8]) for i in range(0, len(p), 8))
    OFF = 1000000
    r = 0.0 if a.no_front else 100.0 * s
    therm = [2.0 * s, 1e-3 * s, a.pen, r, 1.0, h if a.hstab is None else a.hstab, a.eps] + ([a.blend] if a.blend is not None else [])
    if a.felix:
        assert a.blend is None, "--felix replaces --blend"
        therm += [-1.0, 1.0]
    uprop = [1.0, g / 1e10, OFF] + ([1.0] if a.first_order else [])
    sprint = (f"*EL PRINT, ELSET=EALL, FREQUENCY={a.s_every}, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO\n"
              "S, SDV84\n") if a.s_every else ""
    txt = f"""*HEADING
 Klempt 2024 {a.case} in Abaqus: 20 um cube, {n}^3 C3D8T, Eq. 34 as printed (front term r = {r}), Eq. 35 zero order
*NODE
{chr(10).join(f"{i}, {x:.10g}, {y:.10g}, {z:.10g}" for i, x, y, z in nodes)}
*ELEMENT, TYPE=C3D8T, ELSET=EALL
{chr(10).join(f"{i}, " + ", ".join(map(str, c)) for i, c in elems)}
*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES={len(uprop)}, COORDINATES=3, VARIABLES=1, UNSYMM
11, 12
*ELEMENT, TYPE=U1, ELSET=NUTEL
{chr(10).join(f"{i + OFF}, " + ", ".join(map(str, c)) for i, c in elems)}
*UEL PROPERTY, ELSET=NUTEL
{", ".join(f"{x:.10g}" for x in uprop)}
*NSET, NSET=NALL, GENERATE
1, {(n + 1) ** 3}, 1
*NSET, NSET=NSEED
{lst(seed)}
*NSET, NSET=NSRC
{lst(src)}
{("*NSET, NSET=NTOP" + chr(10) + lst(top)) if top else "** no top face"}
*SOLID SECTION, ELSET=EALL, MATERIAL=BIOFILM
*MATERIAL, NAME=BIOFILM
*DENSITY
1.0
*USER MATERIAL, CONSTANTS={len(p)}, TYPE=MECHANICAL
{rows}
*USER MATERIAL, CONSTANTS={len(therm)}, TYPE=THERMAL
{chr(10).join(", ".join(f"{x:.10g}" for x in therm[i:i + 8]) for i in range(0, len(therm), 8))}
*DEPVAR
100
*INITIAL CONDITIONS, TYPE=TEMPERATURE
NALL, 0.0
{ic}
*STEP, NLGEOM=YES, INC=100000, UNSYMM=YES
*COUPLED TEMPERATURE-DISPLACEMENT
{a.dt}, {a.T}
*BOUNDARY
{nid(0, 0, 0)}, 1, 3
{nid(n, 0, 0)}, 2, 3
{nid(0, n, 0)}, 3, 3
NSRC, 12, 12, 1.0
{"NTOP, 12, 12, 0.0" + chr(10) if top else ""}*EL PRINT, ELSET=EALL, FREQUENCY={a.every}, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
TEMP, SDV51
{sprint}*OUTPUT, FIELD, FREQUENCY={a.every}
*NODE OUTPUT
NT
*END STEP
"""
    Path(a.out).write_text(txt)
    print(f"wrote {a.out}: {len(elems)} elements, seed {len(seed)} nodes, nutrient source {len(src)} nodes")


if __name__ == "__main__":
    main()
