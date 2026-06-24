"""Crown moment-arm crestal-bone p95 under the EXACT thesis metric (verbatim from
masterarbeit_ansys_fem/extensions/fig_implant_crown_fem.py: peri_implant_bone_peak), applied to the
linear (C3D4) and quadratic (C3D10) jobs so the headline 18->27 / x1.5 can be re-stated at second order.

Metric (identical to the thesis figure): p95 of occlusal vM over CORTICAL/CANCELLOUS/BONE elements in the
crestal peri-implant DISK r<=3.0 mm of the implant axis T23=(-69.4,-41.0), z in [CREST-3, CREST+1.5] with
CREST=29.0. p95 (not max) to skip the crown-insensitive neck/TIE singularity.

Run:  python crestal_p95_thesis_metric.py
"""
import json
import numpy as np

HERE = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
T23 = np.array([-69.4, -41.0]); CREST = 29.0


def p95(job):
    d = json.load(open("%s/%s_field.json" % (HERE, job)))["els"]
    x = np.array([e["x"] for e in d]); y = np.array([e["y"] for e in d]); z = np.array([e["z"] for e in d])
    vm = np.array([e["vmo"] for e in d]); mat = np.array([e["mat"] for e in d])
    r = np.hypot(x - T23[0], y - T23[1])
    pib = np.isin(mat, ["CORTICAL", "CANCELLOUS", "BONE"]) & (r <= 3.0) & (z >= CREST - 3) & (z <= CREST + 1.5)
    return float(np.percentile(vm[pib], 95)), int(pib.sum())


jobs = [("crown C3D4", "tier2b_crown"), ("bare  C3D4", "tier2b_generic"),
        ("crown C3D10", "tier2b_crown_q"), ("bare  C3D10", "tier2b_generic_q")]
v = {}
for lab, j in jobs:
    p, n = p95(j)
    v[j] = p
    print("%-12s  crestal bone p95 = %5.1f MPa  (n=%d)" % (lab, p, n))

r4 = v["tier2b_crown"] / v["tier2b_generic"]
r10 = v["tier2b_crown_q"] / v["tier2b_generic_q"]
print("\nmoment-arm ratio (crown / bare):")
print("  C3D4  : %5.1f / %5.1f = x%.2f   (thesis headline ~18->27, x1.5)" %
      (v["tier2b_crown"], v["tier2b_generic"], r4))
print("  C3D10 : %5.1f / %5.1f = x%.2f" % (v["tier2b_crown_q"], v["tier2b_generic_q"], r10))
print("  crown delta %+.1f%%   bare delta %+.1f%%   ratio delta %+.1f%%" %
      (100 * (v["tier2b_crown_q"] - v["tier2b_crown"]) / v["tier2b_crown"],
       100 * (v["tier2b_generic_q"] - v["tier2b_generic"]) / v["tier2b_generic"],
       100 * (r10 - r4) / r4))
