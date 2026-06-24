"""Angular distribution of crestal peri-implant bone stress (study (3): eccentric load -> saucerisation
direction). Run: abaqus python extract_pimp_angular.py <tag> <z_bone_top> [nsec=12]
Bins BONE elements within 1.2 mm of the bone level by azimuth theta=atan2(y,x) and reports the mean &
p90 vM per sector + the peak sector -> where the saucer-shaped defect would start. Appends to
ci_ecc_results.jsonl."""
from __future__ import print_function
import sys, json, math
from odbAccess import openOdb

TAG = sys.argv[1]; ZBT = float(sys.argv[2]); NSEC = int(sys.argv[3]) if len(sys.argv) > 3 else 12
JOB = "pimp_%s" % TAG
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
vm = {}
for v in fr.fieldOutputs["S"].values:
    s = v.data
    vm[v.elementLabel] = math.sqrt(0.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*(s[3]**2+s[4]**2+s[5]**2))
o.close()

sec = [[] for _ in range(NSEC)]
for e in emat:
    if emat[e] != "BONE" or e not in C:
        continue
    x, y, z = C[e][0], C[e][1], C[e][2]
    if abs(z - ZBT) > 1.2:
        continue
    th = math.atan2(y, x)                       # +x = buccal (the eccentric-load side)
    k = int((th + math.pi) / (2*math.pi) * NSEC) % NSEC
    sec[k].append(vm[e])


def p(a, q):
    if not a:
        return 0.0
    a = sorted(a); return a[int(q*(len(a)-1))]


centers = [-math.pi + (k+0.5)*2*math.pi/NSEC for k in range(NSEC)]
mean = [(sum(a)/len(a) if a else 0.0) for a in sec]
p90 = [p(a, 0.90) for a in sec]
kmax = max(range(NSEC), key=lambda k: p90[k])
rec = {"tag": TAG, "z_bone_top": ZBT, "nsec": NSEC, "theta": centers, "mean": mean, "p90": p90,
       "peak_theta_deg": math.degrees(centers[kmax]), "peak_p90": p90[kmax],
       "asym": (max(p90) / max(min(p90), 1e-6))}
open("ci_ecc_results.jsonl", "a").write(json.dumps(rec) + "\n")
print("%s  peak@%.0fdeg  p90=%.1f  asym(max/min)=%.2f" %
      (TAG, rec["peak_theta_deg"], rec["peak_p90"], rec["asym"]))
