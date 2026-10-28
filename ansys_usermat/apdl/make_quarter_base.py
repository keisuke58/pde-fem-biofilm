"""make_quarter_base.py -- A-55: the quarter x > 0, z > 0 of a base deck (base_w<n>_c6_g1_s015.dat) for the separate
and mixed examples, which are symmetric in x and z. The nodes and elements are renumbered 1..N (the partner's APDL
indexes arrays of size NODE_CNT by node number), the element and node components are carried over, and the boundary
node components follow the quarter box: XMIN is the plane x = 0, ZMIN the plane z = 0 (their outward normals are -x and
-z, as for a wall), and ALLSURFACENODES, WALL_BC, IC_PHI_NODE are all nodes on its six faces. The symmetry planes are
thus zero-flux walls for phi, c and a, as the mirror symmetry requires. make_living_deck.py --quarter then replaces the
corner constraints by the symmetry constraints.

    python make_quarter_base.py base_w16_c6_g1_s015.dat base_w16q_c6_g1_s015.dat
"""
import re
import sys
from pathlib import Path

import numpy as np

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
text = src.read_text(encoding="latin-1")
nl = "\r\n" if "\r\n" in text else "\n"
lines = text.split(nl)

# --- mesh
i_nb = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("nblock"))
nodes, k = {}, i_nb + 2
while not lines[k].strip().startswith("-1"):
    t = lines[k].split()
    nodes[int(t[0])] = np.array([float(x) for x in t[1:4]])
    k += 1
nb_end = k
i_eb = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("eblock"))
elems, k = {}, i_eb + 2
erows = {}
while not lines[k].strip().startswith("-1"):
    t = [int(x) for x in lines[k].split()]
    elems[t[10]] = t[11:19]; erows[t[10]] = t
    k += 1
eb_end = k
cen = {e: np.mean([nodes[n] for n in ns], axis=0) for e, ns in elems.items()}
keep_e = [e for e in sorted(elems) if cen[e][0] > 0 and cen[e][2] > 0]
keep_n = sorted({n for e in keep_e for n in elems[e]})
emap = {e: i + 1 for i, e in enumerate(keep_e)}
nmap = {n: i + 1 for i, n in enumerate(keep_n)}
xyz = {nmap[n]: nodes[n] for n in keep_n}
h = min(abs(a - b) for a in {round(v[0], 9) for v in nodes.values()} for b in {round(v[0], 9) for v in nodes.values()} if a != b)
eps = 1e-6

new_nb = [lines[i_nb].split(",")[0] + ",3,," + str(len(keep_n)), lines[i_nb + 1]]
new_nb += [f"{nmap[n]:9d}" + "".join(f"{v:20.9E}" for v in nodes[n]) for n in keep_n]
new_eb = [f"eblock,19,solid,,{len(keep_e)}", lines[i_eb + 1]]
for e in keep_e:
    t = list(erows[e]); t[10] = emap[e]; t[11:19] = [nmap[n] for n in t[11:19]]
    new_eb.append("".join(f"{v:9d}" for v in t))

# --- components
def block(name, kind, ids):
    ids = sorted(ids)
    out = [f"CMBLOCK,{name:8s},{kind},{len(ids):9d}", "(8i10)"]
    out += ["".join(f"{x:10d}" for x in ids[j:j + 8]) for j in range(0, len(ids), 8)]
    return out

on = lambda p, ax, v: abs(p[ax] - v) < eps                                  # noqa: E731
face = {"XMIN": (0, 0.0), "XMAX": (0, 1.0), "YMIN": (1, -1.0), "YMAX": (1, 1.0), "ZMIN": (2, 0.0), "ZMAX": (2, 1.0)}
surf = [n for n, p in xyz.items() if any(on(p, ax, v) for ax, v in face.values())]
cent_new = {emap[e]: cen[e] for e in keep_e}


def nearest_mirror(old_ids):
    """kept elements nearest to the mirror images (|x|, y, |z|) of the old elements"""
    out = set()
    for e in old_ids:
        c = cen[e] * np.array([-1 if cen[e][0] < 0 else 1, 1, -1 if cen[e][2] < 0 else 1])
        out.add(min(cent_new, key=lambda q: np.sum((cent_new[q] - c) ** 2)))
    return out


rx = re.compile(r"^\s*CMBLOCK,\s*(\w+)\s*,\s*(NODE|ELEM)\s*,\s*(\d+)", re.I)
out, i, report = [], 0, []
while i < len(lines):
    if i == i_nb:
        out += new_nb; i = nb_end; continue
    if i == i_eb:
        out += new_eb; i = eb_end; continue
    m = rx.match(lines[i])
    if not m:
        out.append(lines[i]); i += 1; continue
    name, kind, n = m.group(1).upper(), m.group(2).upper(), int(m.group(3))
    vals, j = [], i + 2
    while len(vals) < n:
        vals += [int(x) for x in lines[j].split()]
        j += 1
    if kind == "ELEM":
        ids = {emap[e] for e in vals if e in emap}
        if name == "IC_PHI_ELEM":
            ids = set(emap.values())
        if not ids and name in ("HEATUP_ELEM", "EXPO_TEMP_ELEM", "NUTRIENT2", "BIOFILM2"):   # NUTRIENT2: replaced by make_living_deck
            ids = nearest_mirror(vals)
    else:
        ids = {nmap[x] for x in vals if x in nmap}
        if name in face:
            ax, v = face[name]
            ids = {q for q, p in xyz.items() if on(p, ax, v)}
        elif name in ("ALLSURFACENODES", "WALL_BC", "IC_PHI_NODE"):
            ids = set(surf)
        elif name == "ALL_NODES":
            ids = set(xyz)
    report.append(f"{name} {kind}: {n} -> {len(ids)}")
    if ids:
        out += block(name, kind, ids)
    i = j
dst.write_text(nl.join(out), encoding="latin-1", newline="")
print(f"wrote {dst}: {len(keep_n)} of {len(nodes)} nodes, {len(keep_e)} of {len(elems)} elements, h = {h:g}")
print("; ".join(report))
