"""make_cube_inp.py -- the partner's ANSYS model in Abaqus, with phi solved by
Abaqus itself (coupled temperature-displacement, phi = temperature, UMATHT).

    python abaqus_composition/make_cube_inp.py OUT.inp --json RUN.json [--n 8]
        [--beta 0.02] [--kalpha 1e-3] [--pen 5] [--dt 0.025] [--T 1.1]

The 2 mm cube (-1..1 mm) of n^3 C3D8T elements; the seed = the elements whose
centroids are those of BIOFILM1 in an ANSYS run's JSON (export_runs_json.py),
phi = 1 at their nodes and 0 elsewhere (MY_BIOSTART1 = 1); only rigid-body
motion held (3-2-1 at three corners); no flux of phi through the faces
(Eq. 34's boundary condition, natural here). One species, Eq. 36 in the UMAT
(prop(28) = 1), Klempt 2024 stiffness (E = 10 Pa, nu = 0.49), units mm, MPa.
Writes OUT.seed.txt with the seed's Abaqus element numbers.
"""
import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--json", required=True)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--beta", type=float, default=0.02)
    ap.add_argument("--kalpha", type=float, default=1e-3)
    ap.add_argument("--pen", type=float, default=5.0)
    ap.add_argument("--dt", type=float, default=0.025)
    ap.add_argument("--T", type=float, default=1.1)
    ap.add_argument("--ic", choices=("nodes", "fraction"), default="fraction",
                    help="phi at t = 0: 1 at every node of a seed element (nodes), or the share of "
                         "seed elements among the elements around the node (fraction; closer to the "
                         "partner element, whose phi lives at the integration points: 1 at the seed's, 0 elsewhere)")
    a = ap.parse_args()
    n, h = a.n, 2.0 / a.n
    rec = json.loads(Path(a.json).read_text())
    st = rec["all_stress"]
    cen = {e: (x, y, z) for e, x, y, z in zip(st["elem"], st["cx"], st["cy"], st["cz"])}
    seed_c = [cen[e] for e in rec["seed_BIOFILM1"]]

    nid = lambda i, j, k: 1 + i + (n + 1) * j + (n + 1) ** 2 * k  # noqa: E731
    nodes = [(nid(i, j, k), -1 + h * i, -1 + h * j, -1 + h * k)
             for k in range(n + 1) for j in range(n + 1) for i in range(n + 1)]
    elems, seed, seed_nodes = [], [], set()
    around = {}                                      # node -> [elements around it, seed elements]
    e = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                e += 1
                con = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                       nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                elems.append((e, con))
                for q in con:
                    around.setdefault(q, [0, 0])[0] += 1
                c = (-1 + h * (i + 0.5), -1 + h * (j + 0.5), -1 + h * (k + 0.5))
                if any(max(abs(c[q] - s[q]) for q in range(3)) < 1e-6 for s in seed_c):
                    seed.append(e)
                    seed_nodes.update(con)
                    for q in con:
                        around[q][1] += 1
    assert len(seed) == len(seed_c), f"seed: {len(seed)} of {len(seed_c)} centroids found"

    p = [0.0] * 47
    p[0], p[6], p[27] = 1.0, a.kalpha, 1.0          # prop(28) = 1: one species, Eq. 36
    p[42], p[43], p[44], p[45] = 1e-5, -1e-3, 0.49, 0.3
    p[46] = 1.0                                      # phi = temperature
    rows = lambda v: "\n".join(", ".join(f"{x:.17g}" for x in v[i:i + 8]) for i in range(0, len(v), 8))  # noqa: E731
    lst = lambda v: "\n".join(", ".join(str(x) for x in v[i:i + 16]) for i in range(0, len(v), 16))  # noqa: E731
    c1, c2, c3 = nid(0, 0, 0), nid(n, 0, 0), nid(0, n, 0)
    if a.ic == "nodes":
        ic = "NSEED, 1.0"
    else:
        ic = "\n".join(f"{q}, {around[q][1] / around[q][0]:.10g}" for q in sorted(seed_nodes))
    txt = f"""*HEADING
 Partner's 2 mm cube, {n}^3 C3D8T, phi solved by Abaqus (UMATHT), beta {a.beta}, initial phi {a.ic}, from {Path(a.json).name}
*NODE
{chr(10).join(f"{i}, {x:.10g}, {y:.10g}, {z:.10g}" for i, x, y, z in nodes)}
*ELEMENT, TYPE=C3D8T, ELSET=EALL
{chr(10).join(f"{i}, " + ", ".join(map(str, c)) for i, c in elems)}
*NSET, NSET=NALL, GENERATE
1, {(n + 1) ** 3}, 1
*NSET, NSET=NSEED
{lst(sorted(seed_nodes))}
*ELSET, ELSET=SEED
{lst(seed)}
*SOLID SECTION, ELSET=EALL, MATERIAL=BIOFILM
*MATERIAL, NAME=BIOFILM
*DENSITY
1.0
*USER MATERIAL, CONSTANTS=47, TYPE=MECHANICAL
{rows(p)}
*USER MATERIAL, CONSTANTS=3, TYPE=THERMAL
{a.beta}, {a.kalpha}, {a.pen}
*DEPVAR
100
*INITIAL CONDITIONS, TYPE=TEMPERATURE
NALL, 0.0
{ic}
*STEP, NLGEOM=YES, INC=100000
*COUPLED TEMPERATURE-DISPLACEMENT
{a.dt}, {a.T}
*BOUNDARY
{c1}, 1, 3
{c2}, 2, 3
{c3}, 3, 3
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
S
*EL PRINT, ELSET=EALL, FREQUENCY=100000, POSITION=CENTROIDAL, SUMMARY=NO, TOTALS=NO
SDV84, TEMP
*OUTPUT, FIELD
*ELEMENT OUTPUT
S, SDV
*NODE OUTPUT
NT, U
*END STEP
"""
    Path(a.out).write_text(txt)
    Path(a.out).with_suffix(".seed.txt").write_text(" ".join(map(str, seed)))
    print(f"wrote {a.out}: {len(elems)} elements, seed {len(seed)} ({len(seed_nodes)} nodes)")


if __name__ == "__main__":
    main()
