"""make_maze_inp.py -- Abaqus input for the maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11; the mask read off the figure,
abaqus_living/maze_k24.txt) with the agreed coupling in the user element living_uel.f90 (A-67): one layer of elements
(semi-2D as the paper), the walls not discretised, u_n = 0 on every wall face (the paper lets the biofilm slide along
the walls), u_z = 0 on both faces of the layer, nutrient c = 1 on the nodes of the source cells, the colony phi = 1 at
the entrance, the antibiotic on the entrance nodes from T* = t_ab (--ab).

    python make_maze_inp.py out.inp [--cell 4] [--split 1] [--dt 0.025] [--tend 20] [--ab] [--t-ab 20/3] [--supg]
        --split k: each 4 um cell into k x k x 1 elements

Units mm, T*. Parameters of Klempt 2024 Table 6 converted to mm: beta = 10 um^2/T* = 1e-5 mm^2/T*, r = 200 um/T* =
0.2 mm/T*, k_alpha = 1e-6 /T*, K = 1, q = g/d = 1e8/1e13 /um^2 = 10 /mm^2 (the element solves lap c - q phi L c = 0);
the rest as the cube examples: E = 10 Pa, nu = 0.49, f = 1e-3, penalty 5, phi_min 0.001, s = 0.1, start (0.25, 0.30)
of Fritsch case I in the colony, mixed start for newly reached points (layout 2). Stabilisation of the front term:
SUPG (--supg) or beta + r h/4.
"""
import argparse
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--mask", default=str(HERE / "maze_k24.txt"))
ap.add_argument("--cell", type=float, default=4.0, help="cell size of the mask [um]")
ap.add_argument("--split", type=int, default=1)
ap.add_argument("--dt", type=float, default=0.025)
ap.add_argument("--tend", type=float, default=20.0)
ap.add_argument("--ab", action="store_true")
ap.add_argument("--t-ab", type=float, default=None)
ap.add_argument("--supg", action="store_true")
a = ap.parse_args()

m = np.loadtxt(a.mask, dtype=str, comments="#")
m = np.array([[int(ch) for ch in row] for row in m])          # rows = y of the figure (top first), columns = x
m = np.kron(m, np.ones((a.split, a.split), int))[::-1]              # rows upward (y = +y at the entrance)
h = a.cell / a.split * 1e-3                                       # element size [mm]
ny, nx = m.shape
beta0, r_front, K, q, ka, s_link = 1e-5, 0.2, 1.0, 10.0, 1e-6, 0.1
n_sub = max(1, math.ceil(s_link * a.dt / 1e-4 - 1e-9))
beta = beta0 + (0.0 if a.supg else r_front * h / 4.0)
tot = 0.55
props = [10.0, 0.49, 1e-3, beta, ka, 5.0, r_front, K, q, q, 0.001, s_link, n_sub, 2, 0.0, tot, 0.25, 0.30]
if a.supg:
    props += [0.0, 1.0]

# nodes: (i, j, k) with i along x (columns), j along y (rows, figure top = j = 0 ... keep as is), k = 0, 1
nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)  # noqa: E731
off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
elems, used = [], set()
for j in range(ny):
    for i in range(nx):
        if m[j, i] > 0:
            nodes = [nid(i + di, j + dj, dk) for di, dj, dk in off]
            elems.append(((i, j), nodes)); used.update(nodes)
# nodal phi: mean of the cell values (colony cells 1) over the elements of the node
phi_node = {}; cnt = {}
for (i, j), nodes in elems:
    v = 1.0 if m[j, i] == 2 else 0.0
    for nd in nodes:
        phi_node[nd] = phi_node.get(nd, 0.0) + v; cnt[nd] = cnt.get(nd, 0) + 1
phi_node = {nd: phi_node[nd] / cnt[nd] for nd in phi_node}
# node sets
def cell_nodes(code):
    out = set()
    for (i, j), nodes in elems:
        if m[j, i] == code:
            out.update(nodes)
    return sorted(out)
supply, colony = cell_nodes(3), cell_nodes(2)
# wall faces: a path cell whose neighbour (in x or y) is a wall or outside -> the two (x4 nodes) of that face get u_n = 0
fix_x, fix_y = set(), set()
for (i, j), nodes in elems:
    for di, dj, axis in ((1, 0, "x"), (-1, 0, "x"), (0, 1, "y"), (0, -1, "y")):
        ii, jj = i + di, j + dj
        if not (0 <= ii < nx and 0 <= jj < ny) or m[jj, ii] == 0:
            # explicit: the face nodes
            if axis == "x":
                i0 = i + (1 if di == 1 else 0)
                fs = [nid(i0, j, 0), nid(i0, j + 1, 0), nid(i0, j, 1), nid(i0, j + 1, 1)]
                fix_x.update(fs)
            else:
                j0 = j + (1 if dj == 1 else 0)
                fs = [nid(i, j0, 0), nid(i + 1, j0, 0), nid(i, j0, 1), nid(i + 1, j0, 1)]
                fix_y.update(fs)

out = ["*HEADING", f"maze of Klempt 2024 Fig. 11, agreed coupling, user element living_uel, {h * 1e3:g} um cells, one layer",
       "*NODE"]
for k in range(2):
    for j in range(ny + 1):
        for i in range(nx + 1):
            nd = nid(i, j, k)
            if nd in used:
                out.append(f"{nd}, {(i - nx / 2) * h:.8g}, {(j - ny / 2) * h:.8g}, {k * h:.8g}")
out += [f"*USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES={len(props)}, COORDINATES=3, VARIABLES=128", "1, 2, 3, 11, 12, 13",
        "*ELEMENT, TYPE=U1, ELSET=ALLE"]
out += [f"{e + 1}, " + ", ".join(map(str, nodes)) for e, (_, nodes) in enumerate(elems)]
vals = [f"{v:.12g}" for v in props]
out += ["*UEL PROPERTY, ELSET=ALLE"] + [", ".join(vals[k:k + 8]) for k in range(0, len(vals), 8)]
out += ["*ELEMENT, TYPE=C3D8T, ELSET=DUMMY"]
out += [f"{100000 + e + 1}, " + ", ".join(map(str, nodes)) for e, (_, nodes) in enumerate(elems)]
out += ["*SOLID SECTION, ELSET=DUMMY, MATERIAL=DUMMY", "", "*MATERIAL, NAME=DUMMY", "*ELASTIC", "1e-9, 0.3",
        "*CONDUCTIVITY", "1e-9", "*SPECIFIC HEAT", "1e-9", "*DENSITY", "1.0"]


def nset(name, ids):
    ids = sorted(ids)
    return [f"*NSET, NSET={name}"] + [", ".join(map(str, ids[k:k + 16])) for k in range(0, len(ids), 16)]


out += nset("ALLN", used) + nset("SUPPLY", supply) + nset("COLONY", colony) + nset("FIXX", fix_x) + nset("FIXY", fix_y)
t_ab = a.t_ab if a.t_ab is not None else a.tend / 3.0
if not a.ab:
    t_ab = 10.0 * a.tend
out += ["*AMPLITUDE, NAME=ABON, TIME=STEP TIME", f"0.0, 0.0, {t_ab - 1e-9:.10g}, 0.0, {t_ab:.10g}, 1.0, {100 * a.tend:.10g}, 1.0"]
CTRL = ["*CONTROLS, PARAMETERS=FIELD, FIELD=GLOBAL", "1e-8, 1e-8", "*CONTROLS, PARAMETERS=TIME INCREMENTATION", "20, 30"]
t0 = 1e-8
out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{t0:g}, {t0:g}", *CTRL,
        "*BOUNDARY", "ALLN, 1, 3, 0.0", "ALLN, 13, 13, 0.0", "ALLN, 12, 12, 1.0"]
for nd in sorted(used):
    out.append(f"{nd}, 11, 11, {phi_node.get(nd, 0.0):.12g}")
out += ["*END STEP"]
ninc = int(round(a.tend / a.dt))
out += [f"*STEP, NAME=GROW, INC={ninc + 10}, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT",
        f"{a.dt:.10g}, {a.tend:.10g}", *CTRL,
        "*BOUNDARY, OP=NEW", "ALLN, 3, 3, 0.0", "FIXX, 1, 1, 0.0", "FIXY, 2, 2, 0.0", "SUPPLY, 12, 12, 1.0",
        "*BOUNDARY, OP=NEW, AMPLITUDE=ABON", "COLONY, 13, 13, 1.0",
        f"*OUTPUT, FIELD, FREQUENCY={ninc}", "*NODE OUTPUT", "U, NT", "*END STEP"]
Path(a.out).write_text("\n".join(out) + "\n")
print(f"wrote {a.out}: {len(elems)} elements ({h * 1e3:g} um), {len(used)} nodes, supply {len(supply)}, colony {len(colony)}, "
      f"wall nodes x {len(fix_x)} y {len(fix_y)}, dt {a.dt}, {ninc} increments, n_sub {n_sub}, beta {beta:g}")
