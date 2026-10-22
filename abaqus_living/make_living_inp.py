"""make_living_inp.py -- Abaqus input for the agreed coupling with the user element living_uel.f90 (A-46), on the
n^3 cube [-1, 1]^3 mm of the ANSYS decks and the Python reference, with their examples:

    python make_living_inp.py n layout out.inp [--dt 0.005] [--tend 1] [--t-ab 0.333333] [--no-ab]
        layout: sep (seed species A, the four corner blocks species B) or mixeq (seed only, equal start shares)

Set-up as coupled_living_consumption.py --v2 --ansys-bc --start-position --front 1 --phi-min 0.001 --alpha0 100
--gab 1 --gc 1 --tab 0.333333 --artdiff: front r = 10 mm/T*, K = 1, beta = 0.02 + r h/4, k_alpha = 1e-3 /T*,
penalty 5, g l^2/d = 1 for the nutrient and the antibiotic, phi_min = 0.001, s = 0.1, n_sub as the reference
(ceil(s dt / 1e-4)), E_bio = 10 Pa, f = 1e-3, nu = 0.49. Nutrient c = 1 on the nodes of the element layer y < -0.75
(the ANSYS layer), antibiotic a = 1 on the nodes of the outer element shell from T* = t_ab (0 before); phi: no flux
(nothing to do in a Galerkin element); mechanics: the three corner nodes of wp3_stress_reference.py.

Initial phi at a node: the mean of the cell values (1 in the seed and blocks, 0 elsewhere) of its elements.

Procedure (p3/uel, A-11): a 1e-8 step that prescribes phi, c, a at every node (the only way to set initial values of
these dofs), then the growth step, both transient *COUPLED TEMPERATURE-DISPLACEMENT with fixed increments; UNSYMM;
dummy C3D8T elements with 1e-9 material on the same nodes for the specific heat. Gauss-point records go to the file
BIOFILM_UEL_GP (living_uel.f90, UEXTERNALDB).
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "ansys_usermat/apdl"))
import wp3_reference as W  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("n", type=int)
ap.add_argument("layout", choices=["sep", "mixeq"])
ap.add_argument("out")
ap.add_argument("--dt", type=float, default=0.005)
ap.add_argument("--tend", type=float, default=1.0)
ap.add_argument("--t-ab", type=float, default=0.333333)
ap.add_argument("--no-ab", action="store_true")
ap.add_argument("--obstacle", choices=["columns"], help="the four rigid bars of coupled_living_consumption.py --obstacle columns (Klempt 2024 Sect. 4.3): their elements are left out (no biofilm, no nutrient, zero flux through their faces); the mechanics has no contact yet (A-47)")
ap.add_argument("--front", type=float, default=1.0, help="front term factor (0: off, beta = 0.02; as the Python --front)")
ap.add_argument("--supg", action="store_true", help="A-54: SUPG for the front term (PROPS(20) = 1) instead of the added diffusion r h/4 (beta = beta_0)")
ap.add_argument("--quarter", action="store_true", help="A-53: only x > 0, z > 0 (the examples are symmetric in x and z, as the Python --fast); "
                "symmetry u_x = 0 on x = 0, u_z = 0 on z = 0 and u_y = 0 at (0, -1, 0) instead of the corner constraints")
ap.add_argument("--mms", action="store_true", help="manufactured solution (A-50, mms.py): point model off (L = 1), k_alpha = 0, "
                "phi and c start from the solution, c = 1 on y = -1 only, no antibiotic; PROPS(19) = 1")
a = ap.parse_args()
if a.quarter and (a.obstacle or a.mms):
    sys.exit("--quarter only without --obstacle and --mms")

n = a.n
h = 2.0 / n
r_front, K_monod, s_link = 10.0 * a.front, 1.0, 0.1
beta = W.BETA + (0.0 if a.supg else r_front * h / 4.0)
n_sub = max(1, math.ceil(s_link * a.dt / 1e-4 - 1e-9))
tot = 0.55                                   # Fritsch et al. 2025 case I start, 0.25 + 0.30
props = [10.0, 0.49, 1e-3, beta, W.K_ALPHA, W.PENALTY, r_front, K_monod, 1.0, 1.0, 0.001, s_link, n_sub,
         1 if a.layout == "sep" else 2, 0.5625, tot, tot / 2, tot / 2]
if a.mms:
    props[4] = 0.0
    props.append(1.0)
if a.supg:
    props += [0.0] * (19 - len(props)) + [1.0]
    a.no_ab = True

# cells: seed (wp3 reference) and the four corner blocks at the seed's height (v2)
cells8 = W.seed_cells_8(ROOT / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
seed = W.seed_mask(n, cells8)
ax = -1 + (np.arange(n) + 0.5) * h
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
blocks = (np.abs(X) > 0.75) & (np.abs(Y) > 0.75) & (np.abs(Z) <= np.abs(Z[seed]).max() + 1e-9)
phi_cell = (seed | (blocks if a.layout == "sep" else False)).astype(float)
obst = np.zeros_like(seed)
BOXES = ((((-1.0, -0.125, -0.75), (1.0, 0.125, -0.5)), ((-1.0, -0.125, 0.5), (1.0, 0.125, 0.75)),
          ((-0.125, -0.75, -1.0), (0.125, -0.5, 1.0)), ((-0.125, 0.5, -1.0), (0.125, 0.75, 1.0)))
         if a.obstacle else ())
if a.obstacle:
    for lo, hi in BOXES:
        obst |= (X > lo[0]) & (X < hi[0]) & (Y > lo[1]) & (Y < hi[1]) & (Z > lo[2]) & (Z < hi[2])
    phi_cell[obst] = 0.0

xn = -1 + np.arange(n + 1) * h
nid = lambda i, j, k: 1 + i + (n + 1) * (j + (n + 1) * k)  # noqa: E731
off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
elems = []
for k in range(n):
    for j in range(n):
        for i in range(n):
            if not obst[i, j, k] and (not a.quarter or (X[i, j, k] > 0 and Z[i, j, k] > 0)):
                elems.append(((i, j, k), [nid(i + di, j + dj, k + dk) for di, dj, dk in off]))
phi_node = np.zeros((n + 1,) * 3); cnt = np.zeros((n + 1,) * 3)
for (i, j, k), _ in elems:
    for di, dj, dk in off:
        phi_node[i + di, j + dj, k + dk] += phi_cell[i, j, k]; cnt[i + di, j + dj, k + dk] += 1
used = cnt > 0
phi_node = np.where(used, phi_node / np.maximum(cnt, 1), 0.0)
c_node = np.ones_like(phi_node)
if a.mms:
    sys.path.insert(0, str(HERE))
    import mms as M
    XN, YN, ZN = np.meshgrid(xn, xn, xn, indexing="ij")
    phi_node = M.phi_m(XN, YN, ZN, 0.0); c_node = M.c_m(XN, YN, ZN, 0.0)

out = ["*HEADING", f"agreed coupling (living share), user element living_uel, {n}^3, {a.layout}", "*NODE"]
for k in range(n + 1):
    for j in range(n + 1):
        for i in range(n + 1):
            if used[i, j, k]:
                out.append(f"{nid(i, j, k)}, {xn[i]:.10g}, {xn[j]:.10g}, {xn[k]:.10g}")
out += [f"*USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES={len(props)}, COORDINATES=3, VARIABLES=128", "1, 2, 3, 11, 12, 13",
        "*ELEMENT, TYPE=U1, ELSET=ALLE"]
out += [f"{e + 1}, " + ", ".join(map(str, nodes)) for e, (_, nodes) in enumerate(elems)]
vals = [f"{v:.12g}" for v in props]
out += ["*UEL PROPERTY, ELSET=ALLE"] + [", ".join(vals[k:k + 8]) for k in range(0, len(vals), 8)]
# dummy coupled elements for the specific heat of the transient procedure
out += ["*ELEMENT, TYPE=C3D8T, ELSET=DUMMY"]
out += [f"{100000 + e + 1}, " + ", ".join(map(str, nodes)) for e, (_, nodes) in enumerate(elems)]
out += ["*SOLID SECTION, ELSET=DUMMY, MATERIAL=DUMMY", "", "*MATERIAL, NAME=DUMMY", "*ELASTIC", "1e-9, 0.3",
        "*CONDUCTIVITY", "1e-9", "*SPECIFIC HEAT", "1e-9", "*DENSITY", "1.0"]


def nset(name, ids):
    ids = sorted(ids)
    return [f"*NSET, NSET={name}"] + [", ".join(map(str, ids[k:k + 16])) for k in range(0, len(ids), 16)]


eps = 1e-9
allN = [(i, j, k) for k in range(n + 1) for j in range(n + 1) for i in range(n + 1) if used[i, j, k]]
supply = [nid(*p) for p in allN if xn[p[1]] <= (-1.0 if a.mms else -0.75) + eps]
shell = [nid(*p) for p in allN if max(abs(xn[p[0]]), abs(xn[p[1]]), abs(xn[p[2]])) >= 0.75 - eps]
out += nset("ALLN", [nid(*p) for p in allN]) + nset("SUPPLY", supply) + nset("SHELL", shell)
c0, c1, c2 = nid(0, 0, 0), nid(n, 0, 0), nid(0, n, 0)          # (-1,-1,-1), (1,-1,-1), (-1,1,-1)
if not a.quarter:
    out += [f"*NSET, NSET=CORNER0", f"{c0}", f"*NSET, NSET=CORNER1", f"{c1}", f"*NSET, NSET=CORNER2", f"{c2}"]
else:
    out += nset("SYMX", [nid(*p) for p in allN if abs(xn[p[0]]) < eps]) + nset("SYMZ", [nid(*p) for p in allN if abs(xn[p[2]]) < eps])
    out += ["*NSET, NSET=YFIX", f"{nid(n // 2, 0, n // 2)}"]
# frictionless contact with the rigid obstacles (as wp3_stress_reference.solve_fields with obstacles): every box face
# that is not on the cube boundary is a rigid quad (R3D4, normal out of the obstacle) and the mesh nodes on that face a
# node-based surface; hard, frictionless contact. With obstacles the contact holds the body (no corner constraints).
contact, nr = [], 900000
for b, (lo, hi) in enumerate(BOXES):
    for ax in range(3):
        bb, cc = (ax + 1) % 3, (ax + 2) % 3
        for face, sgn in ((lo[ax], -1.0), (hi[ax], 1.0)):
            if abs(abs(face) - 1.0) < eps:
                continue
            on = [nid(*p) for p in allN if abs(xn[p[ax]] - face) < eps
                  and lo[bb] - eps <= xn[p[bb]] <= hi[bb] + eps and lo[cc] - eps <= xn[p[cc]] <= hi[cc] + eps]
            if not on:
                continue
            tag = f"B{b}A{ax}{'P' if sgn > 0 else 'M'}"
            pad = h / 2
            corners = [(lo[bb] - pad, lo[cc] - pad), (hi[bb] + pad, lo[cc] - pad), (hi[bb] + pad, hi[cc] + pad),
                       (lo[bb] - pad, hi[cc] + pad)]
            if sgn < 0:
                corners = corners[::-1]
            ids = []
            for u, v in corners:
                x3 = [0.0, 0.0, 0.0]; x3[ax] = face; x3[bb] = u; x3[cc] = v
                nr += 1; ids.append(nr)
                contact.append(f"{nr}, {x3[0]:.10g}, {x3[1]:.10g}, {x3[2]:.10g}")
            nr += 1; ref = nr
            x3 = [0.0, 0.0, 0.0]; x3[ax] = face
            contact.append(f"{ref}, {x3[0]:.10g}, {x3[1]:.10g}, {x3[2]:.10g}")
            contact.append(f"__ELEM__{tag}__{ids[0]}, {ids[1]}, {ids[2]}, {ids[3]}__REF__{ref}__SLAVE__" + " ".join(map(str, on)))
if contact:
    nodes_c = [c for c in contact if not c.startswith("__ELEM__")]
    out += ["*NODE"] + nodes_c
    ne = 950000
    for c in contact:
        if not c.startswith("__ELEM__"):
            continue
        tag, rest = c[len("__ELEM__"):].split("__", 1)
        quad, rest = rest.split("__REF__")
        ref, slave = rest.split("__SLAVE__")
        ne += 1
        out += [f"*ELEMENT, TYPE=R3D4, ELSET=R{tag}", f"{ne}, {quad}", f"*NSET, NSET=RP{tag}", ref,
                f"*RIGID BODY, REF NODE=RP{tag}, ELSET=R{tag}"]
        out += nset(f"S{tag}", [int(s) for s in slave.split()])
        out += [f"*SURFACE, NAME=SS{tag}, TYPE=NODE", f"S{tag}, 1.0", f"*SURFACE, NAME=MS{tag}", f"R{tag}, SPOS",
                f"*CONTACT PAIR, INTERACTION=FRICLESS, TYPE=NODE TO SURFACE", f"SS{tag}, MS{tag}"]
    # linear penalty instead of hard contact (A-49: the hard contact chattered at T* = 0.17 on 16^3 and the fixed increment
    # could not be cut back); stiffness about 100 E_bio / h, so the penetration is about 1 % of the growth displacement
    kpen = 100.0 * props[0] / h
    out += ["*SURFACE INTERACTION, NAME=FRICLESS", "*FRICTION", "0.0", "*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR",
            f"{kpen:.6g}"]
    out += nset("RPALL", [int(c[len("__ELEM__"):].split("__REF__")[1].split("__")[0]) for c in contact if c.startswith("__ELEM__")])
t_ab = 10.0 if a.no_ab else a.t_ab
out += ["*AMPLITUDE, NAME=ABON, TIME=STEP TIME", f"0.0, 0.0, {t_ab - 1e-9:.10g}, 0.0, {t_ab:.10g}, 1.0, 100.0, 1.0"]
CTRL = ["*CONTROLS, PARAMETERS=FIELD, FIELD=GLOBAL", "1e-8, 1e-8", "*CONTROLS, PARAMETERS=TIME INCREMENTATION", "20, 30"]
if BOXES:
    CTRL = ["*CONTROLS, ANALYSIS=DISCONTINUOUS"] + CTRL
t0 = 1e-8
out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{t0:g}, {t0:g}", *CTRL,
        "*BOUNDARY", "ALLN, 1, 3, 0.0", "ALLN, 13, 13, 0.0"] + ([] if a.mms else ["ALLN, 12, 12, 1.0"]) + (["RPALL, ENCASTRE"] if contact else [])
for p in allN:
    out.append(f"{nid(*p)}, 11, 11, {phi_node[p]:.12g}")
    if a.mms:
        out.append(f"{nid(*p)}, 12, 12, {c_node[p]:.12g}")
out += ["*END STEP"]
ninc = int(round(a.tend / a.dt))
out += [f"*STEP, NAME=GROW, INC={ninc + 10}, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT",
        f"{a.dt:.10g}, {a.tend:.10g}", *CTRL,
        "*BOUNDARY, OP=NEW", *(["RPALL, ENCASTRE"] if contact else ["SYMX, 1, 1, 0.0", "SYMZ, 3, 3, 0.0", "YFIX, 2, 2, 0.0"] if a.quarter else ["CORNER0, 1, 3, 0.0", "CORNER1, 2, 3, 0.0", "CORNER2, 3, 3, 0.0"]),
        "SUPPLY, 12, 12, 1.0",
        "*BOUNDARY, OP=NEW, AMPLITUDE=ABON", "SHELL, 13, 13, 1.0",
        # results come from the Gauss-point file; the .odb only at the end (by default every increment: 270 MB at 16^3)
        f"*OUTPUT, FIELD, FREQUENCY={ninc}", "*NODE OUTPUT", "U, NT", "*END STEP"]
Path(a.out).write_text("\n".join(out) + "\n")
print(f"wrote {a.out}: {n}^3, {a.layout}, dt {a.dt}, {ninc} increments, n_sub {n_sub}, beta {beta:g}, "
      f"seed {int(seed.sum())} cells, blocks {int(blocks.sum()) if a.layout == 'sep' else 0}")
