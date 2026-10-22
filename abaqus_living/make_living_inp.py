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
a = ap.parse_args()

n = a.n
h = 2.0 / n
r_front, K_monod, s_link = 10.0, 1.0, 0.1
beta = W.BETA + r_front * h / 4.0
n_sub = max(1, math.ceil(s_link * a.dt / 1e-4 - 1e-9))
tot = 0.55                                   # Fritsch et al. 2025 case I start, 0.25 + 0.30
props = [10.0, 0.49, 1e-3, beta, W.K_ALPHA, W.PENALTY, r_front, K_monod, 1.0, 1.0, 0.001, s_link, n_sub,
         1 if a.layout == "sep" else 2, 0.5625, tot, tot / 2, tot / 2]

# cells: seed (wp3 reference) and the four corner blocks at the seed's height (v2)
cells8 = W.seed_cells_8(ROOT / "ansys_usermat/apdl/results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
seed = W.seed_mask(n, cells8)
ax = -1 + (np.arange(n) + 0.5) * h
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
blocks = (np.abs(X) > 0.75) & (np.abs(Y) > 0.75) & (np.abs(Z) <= np.abs(Z[seed]).max() + 1e-9)
phi_cell = (seed | (blocks if a.layout == "sep" else False)).astype(float)

xn = -1 + np.arange(n + 1) * h
nid = lambda i, j, k: 1 + i + (n + 1) * (j + (n + 1) * k)  # noqa: E731
off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
elems = []
for k in range(n):
    for j in range(n):
        for i in range(n):
            elems.append(((i, j, k), [nid(i + di, j + dj, k + dk) for di, dj, dk in off]))
phi_node = np.zeros((n + 1,) * 3); cnt = np.zeros((n + 1,) * 3)
for (i, j, k), _ in elems:
    for di, dj, dk in off:
        phi_node[i + di, j + dj, k + dk] += phi_cell[i, j, k]; cnt[i + di, j + dj, k + dk] += 1
phi_node /= cnt

out = ["*HEADING", f"agreed coupling (living share), user element living_uel, {n}^3, {a.layout}", "*NODE"]
for k in range(n + 1):
    for j in range(n + 1):
        for i in range(n + 1):
            out.append(f"{nid(i, j, k)}, {xn[i]:.10g}, {xn[j]:.10g}, {xn[k]:.10g}")
out += ["*USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES=18, COORDINATES=3, VARIABLES=128", "1, 2, 3, 11, 12, 13",
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
allN = [(i, j, k) for k in range(n + 1) for j in range(n + 1) for i in range(n + 1)]
supply = [nid(*p) for p in allN if xn[p[1]] <= -0.75 + eps]
shell = [nid(*p) for p in allN if max(abs(xn[p[0]]), abs(xn[p[1]]), abs(xn[p[2]])) >= 0.75 - eps]
out += nset("ALLN", [nid(*p) for p in allN]) + nset("SUPPLY", supply) + nset("SHELL", shell)
c0, c1, c2 = nid(0, 0, 0), nid(n, 0, 0), nid(0, n, 0)          # (-1,-1,-1), (1,-1,-1), (-1,1,-1)
out += [f"*NSET, NSET=CORNER0", f"{c0}", f"*NSET, NSET=CORNER1", f"{c1}", f"*NSET, NSET=CORNER2", f"{c2}"]
t_ab = 10.0 if a.no_ab else a.t_ab
out += ["*AMPLITUDE, NAME=ABON, TIME=STEP TIME", f"0.0, 0.0, {t_ab - 1e-9:.10g}, 0.0, {t_ab:.10g}, 1.0, 100.0, 1.0"]
CTRL = ["*CONTROLS, PARAMETERS=FIELD, FIELD=GLOBAL", "1e-8, 1e-8", "*CONTROLS, PARAMETERS=TIME INCREMENTATION", "20, 30"]
t0 = 1e-8
out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{t0:g}, {t0:g}", *CTRL,
        "*BOUNDARY", "ALLN, 1, 3, 0.0", "ALLN, 12, 12, 1.0", "ALLN, 13, 13, 0.0"]
for p in allN:
    out.append(f"{nid(*p)}, 11, 11, {phi_node[p]:.12g}")
out += ["*END STEP"]
ninc = int(round(a.tend / a.dt))
out += [f"*STEP, NAME=GROW, INC={ninc + 10}, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT",
        f"{a.dt:.10g}, {a.tend:.10g}", *CTRL,
        "*BOUNDARY, OP=NEW", "CORNER0, 1, 3, 0.0", "CORNER1, 2, 3, 0.0", "CORNER2, 3, 3, 0.0", "SUPPLY, 12, 12, 1.0",
        "*BOUNDARY, OP=NEW, AMPLITUDE=ABON", "SHELL, 13, 13, 1.0", "*END STEP"]
Path(a.out).write_text("\n".join(out) + "\n")
print(f"wrote {a.out}: {n}^3, {a.layout}, dt {a.dt}, {ninc} increments, n_sub {n_sub}, beta {beta:g}, "
      f"seed {int(seed.sum())} cells, blocks {int(blocks.sum()) if a.layout == 'sep' else 0}")
