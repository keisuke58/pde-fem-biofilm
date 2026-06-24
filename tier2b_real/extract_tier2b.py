"""Extract per-element centroid, material and von Mises (both steps) for the Tier-2(b) assembly.
Run: abaqus python extract_tier2b.py   (reads tier2b_real.inp + tier2b_real.odb)"""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb

JOB = sys.argv[1] if len(sys.argv) > 1 else "tier2b_real"
nodes = {}; elems = {}; emat = {}; mode = None; curset = None
for line in open(JOB + ".inp"):
    s = line.strip()
    if s.startswith("*"):
        u = s.upper()
        if u.startswith("*NODE") and "OUTPUT" not in u:
            mode = "n"
        elif u.startswith("*ELEMENT"):
            mode = "e" if "C3D4" in u else None
            curset = None
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
            e = int(p[0]); elems[e] = [int(x) for x in p[1:5]]; emat[e] = curset
    except ValueError:
        pass

cent = {}
for e, nl in elems.items():
    if all(n in nodes for n in nl):
        cent[e] = tuple(sum(nodes[n][d] for n in nl) / 4.0 for d in range(3))

o = openOdb(JOB + ".odb")
steps = list(o.steps.values())

def vm_map(step):
    fr = step.frames[-1]
    out = {}
    for v in fr.fieldOutputs["S"].values:
        s = v.data
        out[v.elementLabel] = math.sqrt(0.5 * ((s[0]-s[1])**2 + (s[1]-s[2])**2 + (s[2]-s[0])**2)
                                        + 3*(s[3]**2 + s[4]**2 + s[5]**2))
    return out

vmg = vm_map(steps[0])
vmo = vm_map(steps[-1])
out = []
for e in elems:
    if e not in cent:
        continue
    x, y, z = cent[e]
    out.append({"x": x, "y": y, "z": z, "mat": emat[e],
                "vmg": float(vmg.get(e, 0.0)), "vmo": float(vmo.get(e, 0.0))})
json.dump({"job": JOB, "els": out}, open(JOB + "_field.json", "w"))
o.close()

# quick per-material report
import collections
agg = collections.defaultdict(lambda: [0, 0.0, 0.0])
for r in out:
    a = agg[r["mat"]]; a[0] += 1; a[1] = max(a[1], r["vmg"]); a[2] = max(a[2], r["vmo"])
print("material   n      vM_growth_max  vM_occlusal_max")
for m, a in sorted(agg.items()):
    print("%-9s %7d  %12.3f  %12.3f" % (m, a[0], a[1], a[2]))
print("wrote %s_field.json (%d elems)" % (JOB, len(out)))
