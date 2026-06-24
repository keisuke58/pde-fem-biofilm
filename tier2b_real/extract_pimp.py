"""Peri-implantitis coupon metrics. Run: abaqus python extract_pimp.py <tag> <z_bone_top>
Appends JSON to pimp_results.jsonl:
  crest_p95 : p95 vM in bone within 1.2 mm of the current bone level (resorption-relevant crestal stress)
  bone_max  : max vM in bone
  thread_exp: max vM in Ti over the EXPOSED thread (bone_top < z < platform L=10)
  ti_max    : max vM in Ti
  disp_um   : max nodal displacement; stiffness = 100/disp
"""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb

TAG = sys.argv[1]; ZBT = float(sys.argv[2]); L = 10.0
JOB = "pimp_%s" % TAG
emat = {}; mode = None; cur = None
for line in open(JOB + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper()
        mode = "e" if u.startswith("*ELEMENT") else None
        if mode:
            for tok in s.split(","):
                if tok.strip().upper().startswith("ELSET"):
                    cur = tok.split("=")[1].strip()
        continue
    if not s or s.startswith("**"):
        continue
    if mode == "e":
        try:
            emat[int(s.split(",")[0])] = cur
        except ValueError:
            pass

o = openOdb(JOB + ".odb"); fr = o.steps.values()[-1].frames[-1]
Z = {v.elementLabel: v.data[2] for v in fr.fieldOutputs["COORD"].values}
vm = {}
for v in fr.fieldOutputs["S"].values:
    s = v.data
    vm[v.elementLabel] = math.sqrt(0.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*(s[3]**2+s[4]**2+s[5]**2))
umax = 0.0
for v in fr.fieldOutputs["U"].values:
    d = v.data; umax = max(umax, math.sqrt(d[0]**2+d[1]**2+d[2]**2))
o.close()


def pct(a, q):
    if not a:
        return 0.0
    a = sorted(a); return a[int(q*(len(a)-1))]

bone = [vm[e] for e in emat if emat[e] == "BONE"]
crest = [vm[e] for e in emat if emat[e] == "BONE" and e in Z and abs(Z[e]-ZBT) <= 1.2]
texp = [vm[e] for e in emat if emat[e] == "TI" and e in Z and ZBT <= Z[e] <= L]
ti = [vm[e] for e in emat if emat[e] == "TI"]
rec = {"tag": TAG, "z_bone_top": ZBT, "crest_p95": pct(crest, 0.95),
       "bone_max": max(bone) if bone else 0, "thread_exp": max(texp) if texp else 0,
       "ti_max": max(ti) if ti else 0, "disp_um": umax*1000.0,
       "stiffness_N_per_um": (100.0/umax/1000.0) if umax else 0}
open("pimp_results.jsonl", "a").write(json.dumps(rec) + "\n")
print(json.dumps(rec))
