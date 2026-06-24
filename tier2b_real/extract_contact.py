"""Micromotion of a debonded implant. Run: abaqus python extract_contact.py <tag> <z_bone_top>
Relative interface displacement (implant node vs nearest socket node) = micromotion; compare 150 um."""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb
TAG = sys.argv[1]; ZBT = float(sys.argv[2]); D = 4.1; JOB = "cc_%s" % TAG

# parse nodes + element-set membership from inp
ncoord = {}; tiN = set(); boN = set(); mode = None; cur = None
for line in open(JOB + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper()
        if u.startswith("*NODE") and "OUTPUT" not in u:
            mode = "n"
        elif u.startswith("*ELEMENT"):
            mode = "e"
            for tok in s.split(","):
                if tok.strip().upper().startswith("ELSET"):
                    cur = tok.split("=")[1].strip()
        else:
            mode = None
        continue
    if not s or s.startswith("**"):
        continue
    p = s.split(",")
    if mode == "n" and len(p) >= 4:
        ncoord[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    elif mode == "e" and len(p) >= 5:
        tgt = tiN if cur == "TI" else boN
        for v in p[1:5]:
            tgt.add(int(v))

o = openOdb(JOB + ".odb"); fr = o.steps.values()[-1].frames[-1]
U = {}
for v in fr.fieldOutputs["U"].values:
    U[v.nodeLabel] = (float(v.data[0]), float(v.data[1]), float(v.data[2]))
o.close()

# implant interface nodes: embedded (z<ZBT) and near flank radius
imp_if = [n for n in tiN if n in ncoord and ncoord[n][2] < ZBT and
          math.hypot(ncoord[n][0], ncoord[n][1]) > D / 2 - 0.6]
bon_if = [n for n in boN if n in ncoord and ncoord[n][2] < ZBT and
          math.hypot(ncoord[n][0], ncoord[n][1]) > D / 2 - 0.6]
# nearest bone node to each implant interface node -> relative displacement
rel = []
for ni in imp_if:
    xi = ncoord[ni]
    best = None; bd = 1e9
    for nb in bon_if:
        xb = ncoord[nb]
        dd = (xi[0]-xb[0])**2 + (xi[1]-xb[1])**2 + (xi[2]-xb[2])**2
        if dd < bd:
            bd = dd; best = nb
    if best is None or ni not in U or best not in U:
        continue
    du = [U[ni][k] - U[best][k] for k in range(3)]
    rel.append(math.sqrt(du[0]**2 + du[1]**2 + du[2]**2) * 1000.0)   # um
rel.sort()
mm_max = rel[-1] if rel else 0.0
mm_p95 = rel[int(0.95*(len(rel)-1))] if rel else 0.0
rec = {"tag": TAG, "z_bone_top": ZBT, "micromotion_max_um": mm_max, "micromotion_p95_um": mm_p95,
       "brunski_150um_exceeded": mm_max > 150.0, "n_pairs": len(rel)}
open("contact_results.jsonl", "a").write(json.dumps(rec) + "\n")
print(json.dumps(rec))
