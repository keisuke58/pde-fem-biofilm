"""Transmucosal axisymmetric field. Run: abaqus python extract_tm.py <tag>. Writes tm_<tag>_field.json."""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb
TAG = sys.argv[1]; JOB = "tm_%s" % TAG
emat = {}; mode = None; cur = None
for line in open(JOB + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper(); mode = "e" if u.startswith("*ELEMENT") else None
        if mode:
            for t in s.split(","):
                if t.strip().upper().startswith("ELSET"):
                    cur = t.split("=")[1].strip()
        continue
    if mode == "e" and s and not s.startswith("**"):
        try:
            emat[int(s.split(",")[0])] = cur
        except ValueError:
            pass
o = openOdb(JOB + ".odb"); fr = o.steps.values()[-1].frames[-1]
C = {v.elementLabel: v.data for v in fr.fieldOutputs["COORD"].values}
out = []
for v in fr.fieldOutputs["S"].values:
    e = v.elementLabel; s = v.data
    if e not in C:
        continue
    s11, s22, s33, s12 = float(s[0]), float(s[1]), float(s[2]), float(s[3])
    vm = math.sqrt(0.5*((s11-s22)**2+(s22-s33)**2+(s33-s11)**2)+3*s12**2)
    out.append({"r": float(C[e][0]), "z": float(C[e][1]), "mat": emat.get(e), "vm": vm})
o.close()
json.dump(out, open("tm_%s_field.json" % TAG, "w"))
print("wrote tm_%s_field.json (%d elems)" % (TAG, len(out)))
