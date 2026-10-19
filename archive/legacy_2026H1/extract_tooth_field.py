"""Element centroid (from INP) + von Mises (from odb) for a conformal biofilm tooth.
Run: abaqus python extract_tooth_field.py <job_basename>  (reads <job>.inp and <job>.odb)"""
from __future__ import annotations
import sys, json, math
from odbAccess import openOdb
job = sys.argv[1]
# parse INP for nodes + C3D4 elements
nodes = {}; elems = {}
mode = None
for line in open(job + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper()
        if u.startswith("*NODE"): mode = "n"
        elif u.startswith("*ELEMENT"): mode = "e" if "C3D4" in u else None
        else: mode = None
        continue
    if not s or s.startswith("**"): continue
    p = s.split(",")
    try:
        if mode == "n" and len(p) >= 4:
            nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
        elif mode == "e" and len(p) >= 5:
            elems[int(p[0])] = [int(x) for x in p[1:5]]
    except ValueError:
        pass
cent = {e: tuple(sum(nodes[n][d] for n in nl)/len(nl) for d in range(3)) for e, nl in elems.items() if all(n in nodes for n in nl)}
o = openOdb(job + ".odb"); fr = list(o.steps.values())[-1].frames[-1]
S = {v.elementLabel: v.data for v in fr.fieldOutputs["S"].values}
out = []
for el, s in S.items():
    if el not in cent: continue
    vm = math.sqrt(0.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*(s[3]**2+s[4]**2+s[5]**2))
    x, y, z = cent[el]
    out.append({"x": x, "y": y, "z": z, "vm": float(vm)})
json.dump({"job": job, "els": out}, open(job + "_3d.json", "w"))
print("%s: %d elems, vM max=%.4g, xyz extent x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f]" % (
    job, len(out), max(e["vm"] for e in out),
    min(e["x"] for e in out), max(e["x"] for e in out),
    min(e["y"] for e in out), max(e["y"] for e in out),
    min(e["z"] for e in out), max(e["z"] for e in out)))
o.close()
