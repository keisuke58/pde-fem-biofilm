"""Dump bone-element (r, z, vM) field for the mechanostat map. Run: abaqus python extract_pimp_field.py <tag>"""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb
TAG = sys.argv[1]; JOB = "pimp_%s" % TAG
emat = {}; mode = None; cur = None
for line in open(JOB + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper(); mode = "e" if u.startswith("*ELEMENT") else None
        if mode:
            for tok in s.split(","):
                if tok.strip().upper().startswith("ELSET"):
                    cur = tok.split("=")[1].strip()
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
    e = v.elementLabel
    if emat.get(e) != "BONE" or e not in C:
        continue
    s = v.data
    vm = math.sqrt(0.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*(s[3]**2+s[4]**2+s[5]**2))
    x, y, z = float(C[e][0]), float(C[e][1]), float(C[e][2])
    out.append({"r": math.hypot(x, y), "x": x, "y": y, "z": z, "vm": float(vm)})
o.close()
json.dump(out, open("pimp_%s_bonefield.json" % TAG, "w"))
print("wrote pimp_%s_bonefield.json (%d bone elems)" % (TAG, len(out)))
