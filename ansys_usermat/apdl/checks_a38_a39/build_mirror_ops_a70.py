"""A-70: MLS operators at the integration points with mirrored neighbours instead of surface points. For a point i
within the stencil radius of a face, the neighbours (and the point itself) are mirrored across that face (and across
two or three faces at edges and corners); the mirrored copies carry the value of their original, so zero flux through
the faces is built in and the fit stays well conditioned (the surface points of the element sit at 0.21 h from the
first Gauss points and their values are themselves fits of those points: the raw stencils there are inconsistent,
error 65 in lap(c_m) on 16^3, A-69). The fit is the element's own (9 terms, Gaussian weights 4 d / BETA_STAR,
BETA_STAR = 0.2), over the GP neighbours of the dumped stencil plus their mirrors; surface columns are dropped.

  python build_mirror_ops_a70.py <dir with nem_dump.txt> [--check] [--lap-cons]    -> <dir>\\mirror_ops.txt (and lap_cons.txt)

Format: nGP; then for each of lap, dx, dy, dz: nnz, row pointers (1-based, 10 per line), "column value" pairs (columns
<= nGP).  --check: residual of lap(c_m) and of grad(c_m) of the manufactured solution (mms.py) by distance to the faces.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

BETA_STAR = float(sys.argv[sys.argv.index("--beta-star") + 1]) if "--beta-star" in sys.argv else 0.2   # weight scale (the element: 0.2 mm)
WTHRES = 1e-8
d = Path(sys.argv[1])
with open(d / "nem_dump.txt") as f:
    nGP, nTot, nnz = (int(x) for x in f.readline().split())
    pos = np.array([f.readline().split()[1:] for _ in range(nTot)], float)
    rows = np.loadtxt(f, max_rows=nnz, usecols=(0, 1, 5))
r, c = rows[:, 0].astype(int) - 1, rows[:, 1].astype(int) - 1
# the stencil of the element: the columns with a non-zero Laplacian coefficient (the sorted GetNhgb keeps fewer than
# the 30 listed neighbours on 16^3: a tie group is dropped, A-44)
nb = [[] for _ in range(nTot)]
for i, j, v in zip(r, c, rows[:, 2]):
    if i != j and j < nGP and v != 0.0:
        nb[i].append(j)
# --cutoff R: the stencil = all integration points within R h (h = element size), symmetric by construction, instead
# of the element's list (whose dropped tie group breaks the reflection symmetry near the walls on 16^3)
if "--cutoff" in sys.argv:
    from scipy.spatial import cKDTree
    hel = 2.0 / round((nGP / 8) ** (1 / 3))
    Rc = float(sys.argv[sys.argv.index("--cutoff") + 1]) * hel
    tree = cKDTree(pos[:nGP])
    nb = [[j for j in tree.query_ball_point(pos[i], Rc) if j != i] for i in range(nGP)] + [[] for _ in range(nTot - nGP)]
    print(f"cutoff {Rc:.4f}: neighbours per point {np.mean([len(x) for x in nb[:nGP]]):.1f} (min {min(len(x) for x in nb[:nGP])}, max {max(len(x) for x in nb[:nGP])})")
lo, hi = pos[:nGP].min(0), pos[:nGP].max(0)
gp_h = (hi - lo).min() / 100                                   # not used; faces are the box faces at +-1
FACES = [(a, s) for a in range(3) for s in (-1.0, 1.0)]        # box [-1, 1]^3
radius = max(np.linalg.norm(pos[j] - pos[0]) for j in nb[0])   # stencil radius of point 0 (all stencils alike)
ops = {k: ([1], [], []) for k in ("lap", "dx", "dy", "dz")}
for i in range(nGP):
    js = nb[i]
    pts = [(pos[j], j) for j in js] + [(pos[i], i)]
    # mirror across every face within the radius (and the combinations: edges and corners)
    near = [(a, s) for a, s in FACES if abs(pos[i, a] - s) < radius]
    allpts = list(pts)
    for a, s in near:
        allpts += [(p * np.where(np.arange(3) == a, -1, 1) + np.where(np.arange(3) == a, 2 * s, 0), j) for p, j in allpts]
    P = np.array([p for p, _ in allpts]); J = np.array([j for _, j in allpts])
    dx = P - pos[i]
    keep = np.linalg.norm(dx, axis=1) > 1e-12
    dx, J = dx[keep], J[keep]
    A = np.c_[0.5 * dx[:, 0] ** 2, 0.5 * dx[:, 1] ** 2, 0.5 * dx[:, 2] ** 2, dx[:, 0] * dx[:, 1], dx[:, 0] * dx[:, 2], dx[:, 1] * dx[:, 2], dx]
    w = np.exp(-0.5 * (4 * np.linalg.norm(dx, axis=1) / BETA_STAR) ** 2)
    w[w <= WTHRES] = 0.0
    D = np.linalg.solve(A.T @ (w[:, None] * A), A.T * w)          # 9 x m
    for k, row in (("lap", D[0] + D[1] + D[2]), ("dx", D[6]), ("dy", D[7]), ("dz", D[8])):
        acc = {}
        for j, v in zip(J, row):
            if j != i:                         # a mirrored copy of the centre carries f_i - f_i = 0
                acc[j] = acc.get(j, 0.0) + v
        acc[i] = -sum(acc.values())           # the centre: minus the sum (the fit is in differences)
        ptr, col, val = ops[k]
        for j in sorted(acc):
            if abs(acc[j]) > 1e-15:
                col.append(j + 1); val.append(acc[j])
        ptr.append(len(col) + 1)
with open(d / "mirror_ops.txt", "w") as f:
    f.write(f"{nGP}\n")
    for k in ("lap", "dx", "dy", "dz"):
        ptr, col, val = ops[k]
        f.write(f"{len(col)}\n")
        for q in range(0, len(ptr), 10):
            f.write(" ".join(str(x) for x in ptr[q:q + 10]) + "\n")
        for j, v in zip(col, val):
            f.write(f"{j} {v:.17e}\n")
print("wrote", d / "mirror_ops.txt")
# --lap-cons: the conservative (symmetrised) mirror Laplacian as lap_cons.txt for LAP_FILE = 1
if "--lap-cons" in sys.argv:
    ptr, col, val = ops["lap"]
    Lm = sp.csr_matrix((np.array(val), np.array(col) - 1, np.array(ptr) - 1), shape=(nGP, nGP))
    Ls = (0.5 * (Lm + Lm.T)).tolil(); Ls.setdiag(0.0); Ls = Ls.tocsr()
    Lc = (Ls - sp.diags(np.asarray(Ls.sum(1)).ravel())).tocsr(); Lc.eliminate_zeros(); Lc.sort_indices()
    q = (pos[:nGP] ** 2).sum(1); v = Lc @ q; inner = np.abs(pos[:nGP]).max(1) < 0.5
    print(f"lap_cons from the mirror operator: inner lap(r^2) {v[inner].min():.4f} .. {v[inner].max():.4f}, |column sum| max {np.abs(Lc.sum(0)).max():.1e}")
    with open(d / "lap_cons.txt", "w") as f:
        f.write(f"{nGP} {Lc.nnz}\n")
        p1 = Lc.indptr + 1
        for k in range(0, len(p1), 10):
            f.write(" ".join(str(x) for x in p1[k:k + 10]) + "\n")
        for j, vv in zip(Lc.indices + 1, Lc.data):
            f.write(f"{j} {vv:.17e}\n")
    print("wrote", d / "lap_cons.txt")
# --lap-cons-plain: the plain mirror Laplacian as lap_cons.txt (consistent to the walls; its mass drift on 16^3 is 3e-5
# over T* = 1 at beta = 0.33, A-70, so the conservative fold is not needed)
if "--lap-cons-plain" in sys.argv:
    ptr, col, val = ops["lap"]
    with open(d / "lap_cons.txt", "w") as f:
        f.write(f"{nGP} {len(col)}\n")
        for k in range(0, len(ptr), 10):
            f.write(" ".join(str(x) for x in ptr[k:k + 10]) + "\n")
        for j, vv in zip(col, val):
            f.write(f"{j} {vv:.17e}\n")
    print("wrote", d / "lap_cons.txt", "(plain mirror Laplacian)")
if "--check" in sys.argv:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "abaqus_living"))
    import mms as M
    PI = np.pi
    x, y, z = pos[:nGP].T; t = 0.0
    cm = M.c_m(x, y, z, t)
    A_ = 0.3 * (1 + 0.5 * t); s_ = np.sin(PI * (y + 1) / 4) ** 2; g_ = 1 + 0.5 * np.cos(PI * (x + 1)) * np.cos(PI * (z + 1))
    s_yy = (PI / 4) ** 2 * 2 * np.cos(PI * (y + 1) / 2); g_xx = -0.5 * PI ** 2 * np.cos(PI * (x + 1)) * np.cos(PI * (z + 1))
    lap_ex = -A_ * (s_yy * g_ + s_ * 2 * g_xx)
    dy_ex = -A_ * (PI / 4) * np.sin(PI * (y + 1) / 2) * g_
    def mat(k):
        ptr, col, val = ops[k]
        return sp.csr_matrix((np.array(val), np.array(col) - 1, np.array(ptr) - 1), shape=(nGP, nGP))
    dist = np.minimum(pos[:nGP] + 1, 1 - pos[:nGP]).min(1)
    for k, ex in (("lap", lap_ex), ("dy", dy_ex)):
        res = mat(k) @ cm - ex
        print(f"mirror {k}(c_m) error:", " ".join(f"[{a:.2f},{b:.2f}): {np.abs(res[(dist >= a) & (dist < b)]).max():.2e}"
                                                 for a, b in ((0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 9))))
    L = mat("lap")
    print(f"row sums |max| {np.abs(L.sum(1)).max():.1e}, column sums |max| {np.abs(L.sum(0)).max():.1e}")
