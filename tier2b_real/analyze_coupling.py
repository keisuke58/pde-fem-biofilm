import sys, json, numpy as np
JOB = sys.argv[1] if len(sys.argv) > 1 else "tier2b_real"
d = json.load(open("/home/nishioka/IKM_Hiwi/FEM/tier2b_real/%s_field.json" % JOB))["els"]
x = np.array([r["x"] for r in d]); y = np.array([r["y"] for r in d]); z = np.array([r["z"] for r in d])
mat = np.array([r["mat"] for r in d]); vmo = np.array([r["vmo"] for r in d]); vmg = np.array([r["vmg"] for r in d])
T24 = np.array([-63.9, -41.2]); T23 = np.array([-69.4, -41.0])
r24 = np.hypot(x - T24[0], y - T24[1]); r23 = np.hypot(x - T23[0], y - T23[1])
X0, X1, Y0, Y1, Z0 = -73.0, -59.0, -47.0, -37.5, 15.0
faredge = (x < X0 + 1.0) | (x > X1 - 1.0) | (y < Y0 + 1.0) | (y > Y1 - 1.0) | (z < Z0 + 1.0)
bone = np.isin(mat, ["BONE", "CORTICAL", "CANCELLOUS"]) & ~faredge


def shell(rr, lo, hi, zmax=28.5):
    m = bone & (rr >= lo) & (rr < hi) & (z < zmax)
    if not m.sum():
        return 0, 0, 0
    return int(m.sum()), float(np.percentile(vmo[m], 95)), float(vmo[m].mean())


print("=== OCCLUSAL step: peri-implant vs peri-tooth BONE (interface shell r[2.0,3.5], below crest) ===")
ni, p95i, mni = shell(r23, 2.0, 3.5)
nt, p95t, mnt = shell(r24, 2.0, 3.5)
print("  peri-IMPLANT (tooth23, Ti, NO PDL): n=%d  p95 vM=%.1f  mean=%.1f MPa" % (ni, p95i, mni))
print("  peri-TOOTH   (tooth24, PDL):        n=%d  p95 vM=%.1f  mean=%.1f MPa" % (nt, p95t, mnt))
print("  ratio implant/tooth  p95=%.2f  mean=%.2f" % (p95i / max(p95t, 1e-9), mni / max(mnt, 1e-9)))

mi = bone & (r23 >= 2) & (r23 < 3.5) & (z < 28.5)
mt = bone & (r24 >= 2) & (r24 < 3.5) & (z < 28.5)
print("=== GROWTH step (dysbiotic biofilm) peri-bone p95 ===")
print("  peri-implant=%.2f  peri-tooth=%.2f MPa" % (np.percentile(vmg[mi], 95), np.percentile(vmg[mt], 95)))

# depth profile of occlusal bone stress along each column (load transmission depth)
print("=== occlusal bone p95 vM by depth band (load transmitted through bone) ===")
for zlo, zhi in [(26, 28.5), (23, 26), (20, 23), (16, 20)]:
    mi = bone & (r23 < 3.5) & (z >= zlo) & (z < zhi)
    mt = bone & (r24 < 3.5) & (z >= zlo) & (z < zhi)
    pi = np.percentile(vmo[mi], 95) if mi.sum() else 0
    pt = np.percentile(vmo[mt], 95) if mt.sum() else 0
    print("  z[%2d,%2d): implant=%6.1f  tooth=%6.1f MPa" % (zlo, zhi, pi, pt))
