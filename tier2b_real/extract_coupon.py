"""Extract implant-coupon design metrics. Run: abaqus python extract_coupon.py <tag> <L>
Appends one JSON line to coupon_results.jsonl:
  thread_vM  : p99 von Mises in the embedded TITANIUM threaded body (thread-root stress concentration)
  ti_max     : max von Mises anywhere in Ti
  bone_p95   : p95 von Mises in bone
  crest_p95  : p95 von Mises in bone near the crestal (bone-top) ring
  disp_um    : max nodal displacement magnitude (um); stiffness = 100 N / disp
"""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb

TAG = sys.argv[1]
L = float(sys.argv[2])
JOB = "coupon_%s" % TAG
Z_BONE_TOP = L - 3.0

# centroids + material from INP
nodes = {}; elems = {}; emat = {}; mode = None; curset = None
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
                    curset = tok.split("=")[1].strip()
        else:
            mode = None
        continue
    if not s or s.startswith("**"):
        continue
    p = s.split(",")
    try:
        if mode == "n" and len(p) >= 4:
            nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
        elif mode == "e" and len(p) >= 5:
            e = int(p[0]); emat[e] = curset
    except ValueError:
        pass

o = openOdb(JOB + ".odb")
fr = o.steps.values()[-1].frames[-1]
Z = {v.elementLabel: v.data[2] for v in fr.fieldOutputs["COORD"].values}
vm = {}
for v in fr.fieldOutputs["S"].values:
    s = v.data
    vm[v.elementLabel] = math.sqrt(0.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*(s[3]**2+s[4]**2+s[5]**2))
umax = 0.0
for v in fr.fieldOutputs["U"].values:
    d = v.data
    umax = max(umax, math.sqrt(d[0]**2 + d[1]**2 + d[2]**2))
o.close()


def pct(vals, q):
    if not vals:
        return 0.0
    vals = sorted(vals); k = int(q * (len(vals) - 1))
    return vals[k]

thread = [vm[e] for e in emat if emat[e] == "TI" and e in Z and 0.5 <= Z[e] <= L - 0.3]
ti_all = [vm[e] for e in emat if emat[e] == "TI"]
bone = [vm[e] for e in emat if emat[e] == "BONE"]
crest = [vm[e] for e in emat if emat[e] == "BONE" and e in Z and abs(Z[e] - Z_BONE_TOP) <= 1.2]
rec = {"tag": TAG, "L": L, "thread_vM": pct(thread, 0.99), "ti_max": max(ti_all) if ti_all else 0,
       "bone_p95": pct(bone, 0.95), "crest_p95": pct(crest, 0.95),
       "disp_um": umax * 1000.0, "stiffness_N_per_um": (100.0 / umax / 1000.0) if umax else 0}
open("coupon_results.jsonl", "a").write(json.dumps(rec) + "\n")
print(json.dumps(rec))
