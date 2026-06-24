"""Global C3D4 -> C3D10 quadratic-tet converter for the tier2b coupled implant assembly.

The headline crowned model (tier2b_crown, ~275 k elements) is meshed with LINEAR tets (C3D4). Linear
tets under-predict the curved thread-root / crestal stress concentration by ~9-15% (the very reason the
ISO-14801 coupons were already run at C3D10). This script promotes the WHOLE assembly to quadratic tets
so the crown moment-arm result can be checked at second order, WITHOUT touching build_assembly.py.

Correctness points handled here:
  * ONE GLOBAL edge -> mid-node map over every element, so the conforming PDL-inner / dentin interface
    (which SHARES corner node ids) also shares its mid-edge nodes -> the interface stays watertight at
    quadratic order (a per-body conversion would tear it).
  * NSET=ALLN is extended with EVERY new mid-node, so the growth step's `*TEMPERATURE ALLN, 1.0` drives
    a UNIFORM eigenstrain inside each biofilm element (a mid-node left out of ALLN would default to T=0
    and create a spurious temperature gradient -> wrong growth strain).
  * NSET=FIXED is extended with mid-nodes whose BOTH corner endpoints are fixed (i.e. the edge lies on a
    planar crop face) -> the Dirichlet boundary is not left under-constrained at the new mid-nodes.
  * *CLOAD node lists (TOOTOP / IMPTOP) are kept corner-only: the 100 N resultant is preserved and, by
    St-Venant, the crestal quantity of interest is insensitive to the occlusal load's nodal lumping.
  * *SURFACE element-face (Sn) defs, *TIE, *EMBEDDED ELEMENT, materials and steps pass through verbatim
    -- all valid for C3D10 (Abaqus uses the mid-face nodes automatically).

C3D10 node order: corners n1..n4, then mid-edge n5(1-2) n6(2-3) n7(3-1) n8(1-4) n9(2-4) n10(3-4).

Run:  python convert_c3d4_to_c3d10.py <in.inp> <out.inp>
"""
import sys
from pathlib import Path

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "tier2b_crown.inp")
DST = Path(sys.argv[2] if len(sys.argv) > 2 else "tier2b_crown_q.inp")

# C3D10 mid-edge ordering: (local corner a, local corner b) for mid-nodes 5..10 (0-based corner idx)
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]

lines = SRC.read_text().splitlines()

# ---- split into blocks (header line starting with '*' + its data lines) ----
blocks = []  # (header_line, [data_lines])
cur = None
for ln in lines:
    if ln.startswith("*"):
        if cur is not None:
            blocks.append(cur)
        cur = (ln, [])
    else:
        if cur is None:
            cur = ("", [])
        cur[1].append(ln)
if cur is not None:
    blocks.append(cur)

# ---- pass 1: parse nodes + C3D4 elements ----
node_xyz = {}            # id -> (x,y,z)
max_nid = 0
elem_blocks = []         # index into blocks list for each C3D4 element block -> parsed conn
for bi, (hdr, data) in enumerate(blocks):
    h = hdr.upper().replace(" ", "")
    if h.startswith("*NODE") and not h.startswith("*NODEOUTPUT"):
        for d in data:
            if not d.strip():
                continue
            p = d.split(",")
            nid = int(p[0])
            node_xyz[nid] = (float(p[1]), float(p[2]), float(p[3]))
            if nid > max_nid:
                max_nid = nid
    elif h.startswith("*ELEMENT") and "TYPE=C3D4" in h:
        conn = []
        for d in data:
            if not d.strip():
                continue
            p = [int(x) for x in d.split(",") if x.strip()]
            conn.append((p[0], p[1:5]))     # (eid, [n1,n2,n3,n4])
        elem_blocks.append((bi, conn))

# ---- build ONE global edge -> mid-node map ----
edge_mid = {}
next_nid = max_nid + 1
for _bi, conn in elem_blocks:
    for _eid, c in conn:
        for a, b in EDGES:
            ka, kb = c[a], c[b]
            key = (ka, kb) if ka < kb else (kb, ka)
            if key not in edge_mid:
                edge_mid[key] = next_nid
                xa, xb = node_xyz[key[0]], node_xyz[key[1]]
                node_xyz[next_nid] = (0.5 * (xa[0] + xb[0]), 0.5 * (xa[1] + xb[1]), 0.5 * (xa[2] + xb[2]))
                next_nid += 1
new_mids = list(range(max_nid + 1, next_nid))
print("corners=%d  new mid-nodes=%d  total nodes=%d  elements=%d"
      % (max_nid, len(new_mids), next_nid - 1, sum(len(c) for _, c in elem_blocks)))


def mid(a, b):
    return edge_mid[(a, b) if a < b else (b, a)]


# ---- parse FIXED / ALLN sets so we can augment them ----
def parse_set(data):
    s = set()
    for d in data:
        for x in d.split(","):
            x = x.strip()
            if x:
                s.add(int(x))
    return s


fixed_set = None
for hdr, data in blocks:
    if hdr.upper().replace(" ", "").startswith("*NSET,NSET=FIXED"):
        fixed_set = parse_set(data)
        break

# mid-nodes on a fully-fixed edge (both endpoints fixed) -> also fixed
fixed_mids = []
if fixed_set is not None:
    for (a, b), m in edge_mid.items():
        if a in fixed_set and b in fixed_set:
            fixed_mids.append(m)
    print("FIXED corners=%d  +fixed mid-nodes=%d" % (len(fixed_set), len(fixed_mids)))


def fmt_ids(ids, per=16):
    out = []
    ids = list(ids)
    for k in range(0, len(ids), per):
        out.append(" " + ",".join(str(x) for x in ids[k:k + per]))
    return out


# ---- pass 2: re-emit ----
out = []
for bi, (hdr, data) in enumerate(blocks):
    h = hdr.upper().replace(" ", "")
    if h.startswith("*NODE") and not h.startswith("*NODEOUTPUT"):
        out.append(hdr)
        # original corner nodes (verbatim formatting) + appended mid-nodes
        for d in data:
            out.append(d)
        for m in new_mids:
            x, y, z = node_xyz[m]
            out.append(" %d, %.5f, %.5f, %.5f" % (m, x, y, z))
    elif h.startswith("*ELEMENT") and "TYPE=C3D4" in h:
        out.append(hdr.replace("C3D4", "C3D10").replace("c3d4", "C3D10"))
        # find this block's parsed conn
        conn = next(c for (j, c) in elem_blocks if j == bi)
        for eid, c in conn:
            m5, m6, m7 = mid(c[0], c[1]), mid(c[1], c[2]), mid(c[2], c[0])
            m8, m9, m10 = mid(c[0], c[3]), mid(c[1], c[3]), mid(c[2], c[3])
            out.append(" %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d"
                       % (eid, c[0], c[1], c[2], c[3], m5, m6, m7, m8, m9, m10))
    elif h.startswith("*NSET,NSET=FIXED"):
        out.append(hdr)
        out.extend(data)
        out.extend(fmt_ids(fixed_mids))
    elif h.startswith("*NSET,NSET=ALLN"):
        out.append(hdr)
        out.extend(data)
        out.extend(fmt_ids(new_mids))      # every mid-node gets the growth temperature
    else:
        out.append(hdr)
        out.extend(data)

DST.write_text("\n".join(out) + "\n")
print("wrote", DST)
