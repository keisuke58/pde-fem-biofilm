"""phi totals by region: Python snapshot vs ANSYS early fields (8^3, separate_ab)."""
import glob, json, sys
import numpy as np

N, H = 8, 0.25
z = np.load(glob.glob(sys.argv[1])[0])
d = json.load(open(sys.argv[2]))
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]


def grid(v):
    g = np.full((N, N, N), np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


xc = -0.875 + H * np.arange(N)
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
seed = grid(np.isin(np.array(d["elem"]), d["seed"]).astype(float)) > 0.5
p0 = grid(d["phi"][0])
blk = (p0 > 0.99) & ~seed
wall = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)
print("blocks", blk.sum(), "seed", seed.sum(), "seed at wall", (seed & wall).sum(), "blocks at wall", (blk & wall).sum())
for T, k in ((0.025, 1), (0.05, 2)):
    pp = z[f"phi_{T}"]; ap = grid(d["phi"][k])
    for nm, m in (("seed", seed), ("blocks", blk), ("wall other", wall & ~blk & ~seed), ("interior other", ~wall & ~seed & ~blk)):
        print(f"T*={T} {nm:15s} n={m.sum():3d}  sum phi Py {pp[m].sum():7.2f} ANSYS {ap[m].sum():7.2f}"
              f"  max Py {pp[m].max():.3f} ANSYS {ap[m].max():.3f}")
    i = np.unravel_index(np.argmax(pp), pp.shape)
    print("  Py argmax at", [float(xc[q]) for q in i])
    print("  phi along x at y=0.875, z=0.125  Py", np.round(pp[:, 7, 4], 3), " ANSYS", np.round(ap[:, 7, 4], 3))
    print("  phi along y at x=0.875, z=0.125  Py", np.round(pp[7, :, 4], 3), " ANSYS", np.round(ap[7, :, 4], 3))
