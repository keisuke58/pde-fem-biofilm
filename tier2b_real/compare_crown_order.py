"""Compare crestal peri-implant BONE stress between the linear (C3D4) and quadratic (C3D10) crown
assemblies with ONE identical metric, to test whether the crown moment-arm result survives second order.

NOTE: this uses a SECONDARY annulus shell (r in [2.0,3.5], z in [crest-3,crest]). For the headline
moment-arm number that matches the thesis (27/17 MPa, x1.5), use `crestal_p95_thesis_metric.py`, which
replicates `fig_implant_crown_fem.py::peri_implant_bone_peak` (disk r<=3, z in [26,30.5]) verbatim. The
two shells agree on the CONCLUSION (ratio order-robust) but the annulus over-weights the near-crest
cortical ring, so its ratio Δ is a shell artefact, not physics.

Metric (same for both jobs): occlusal-step von Mises 95th-percentile over CORTICAL+CANCELLOUS bone
elements in the crestal peri-implant shell r in [2.0,3.5] mm of the implant axis, just below the platform
(z in [crest-3, crest]); far-crop-edge elements excluded. The implant axis (cx,cy) and the crest height
are derived from the TI element centroids in each job, so the two meshes are measured the same way.

Run:  python compare_crown_order.py tier2b_crown tier2b_crown_q
"""
import sys
import json
import numpy as np

HERE = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"


def load(job):
    d = json.load(open("%s/%s_field.json" % (HERE, job)))["els"]
    x = np.array([r["x"] for r in d]); y = np.array([r["y"] for r in d]); z = np.array([r["z"] for r in d])
    mat = np.array([r["mat"] for r in d]); vmo = np.array([r["vmo"] for r in d]); vmg = np.array([r["vmg"] for r in d])
    return x, y, z, mat, vmo, vmg


def crestal_p95(job, verbose=True):
    x, y, z, mat, vmo, vmg = load(job)
    ti = mat == "TI"
    # crest = top of the osseointegrated bone ~ platform; derive from TI just below abutment.
    # implant axis from TI centroids in the bony zone (z below the platform):
    zcrest = np.percentile(z[ti], 88)              # near the platform/crest height
    near = ti & (z < zcrest) & (z > zcrest - 8)
    cx, cy = float(x[near].mean()), float(y[near].mean())
    r = np.hypot(x - cx, y - cy)
    # far-edge crop exclusion: drop the outer 1 mm shell of the overall bounding box
    bx0, bx1, by0, by1, bz0 = x.min() + 1, x.max() - 1, y.min() + 1, y.max() - 1, z.min() + 1
    faredge = (x < bx0) | (x > bx1) | (y < by0) | (y > by1) | (z < bz0)
    bone = np.isin(mat, ["BONE", "CORTICAL", "CANCELLOUS"]) & ~faredge
    shell = bone & (r >= 2.0) & (r < 3.5) & (z < zcrest) & (z > zcrest - 3.0)
    n = int(shell.sum())
    p95 = float(np.percentile(vmo[shell], 95)) if n else 0.0
    mean = float(vmo[shell].mean()) if n else 0.0
    pmax = float(vmo[shell].max()) if n else 0.0
    if verbose:
        print("%-18s axis=(%.2f,%.2f) crest z=%.2f  shell n=%d  crestal bone p95=%.2f  mean=%.2f  max=%.2f MPa"
              % (job, cx, cy, zcrest, n, p95, mean, pmax))
    return p95, mean, n


if __name__ == "__main__":
    jobs = sys.argv[1:] or ["tier2b_crown", "tier2b_crown_q"]
    res = {j: crestal_p95(j) for j in jobs}
    # only auto-print the order verdict for a true linear/quadratic pair (X, X_q)
    if len(jobs) == 2 and jobs[1] == jobs[0] + "_q":
        a, b = jobs
        pa, pb = res[a][0], res[b][0]
        print("\nC3D4 -> C3D10 crestal bone p95:  %.2f -> %.2f MPa  (delta %+.1f%%)"
              % (pa, pb, 100 * (pb - pa) / pa if pa else 0))
        print("=> the crestal-bone result is %s to quadratic order"
              % ("ROBUST" if pa and abs(pb - pa) / pa < 0.15 else "SENSITIVE (order matters)"))
