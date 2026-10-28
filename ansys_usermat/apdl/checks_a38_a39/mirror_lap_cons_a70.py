"""A-70: the conservative Laplacian of phi with mirrored neighbours (zero flux built in), consistent to the walls:
the MLS rows of the integration points are built on the lattice extended by the mirror images across the faces
(as build_mirror_ops_a70.py); the ghost rows are the reflections of the real rows; the operator on the extended set
is symmetrised, its diagonal set to minus the row sum, and the ghost columns are folded onto their originals. By the
reflection symmetry the folded operator keeps zero column sums (mass conservation) and the consistency of the mirror
rows (A-69: the symmetric part of the unfolded mirror operator loses the accuracy next to the walls).

  python mirror_lap_cons_a70.py <dir with nem_dump.txt>   -> <dir>\\lap_cons.txt   (LAP_FILE = 1)
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

BETA_STAR, WTHRES = 0.2, 1e-8
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
radius = max(np.linalg.norm(pos[j] - pos[0]) for j in nb[0])
FACES = [(a, s) for a in range(3) for s in (-1.0, 1.0)]


def reflect(p, m):
    q = p.copy()
    for a in range(3):
        if m >> a & 1:
            q[a] = 2 * np.sign(q[a]) - q[a] if abs(q[a]) > 0 else q[a]     # across the face on q's side
    return q


# real rows on the extended lattice: entries keyed by (j, mask); mask bit a = reflected across the face of axis a
# (the face on the side of the point: all mirrored points are outside the box, so the side is the sign of pos)
real = {}
for i in range(nGP):
    near = [a for a, s in FACES if abs(pos[i, a] - s) < radius]
    base = [(pos[j], j, 0) for j in nb[i]] + [(pos[i], i, 0)]
    allp = list(base)
    for a in near:
        s = 1.0 if pos[i, a] > 0 else -1.0
        allp += [(np.where(np.arange(3) == a, 2 * s - p, p), j, m | (1 << a)) for p, j, m in allp]
    P = np.array([p for p, _, _ in allp]); keys = [(j, m) for _, j, m in allp]
    dx = P - pos[i]
    keep = np.linalg.norm(dx, axis=1) > 1e-12
    dx = dx[keep]; keys = [k for k, kp in zip(keys, keep) if kp]
    A = np.c_[0.5 * dx[:, 0] ** 2, 0.5 * dx[:, 1] ** 2, 0.5 * dx[:, 2] ** 2, dx[:, 0] * dx[:, 1], dx[:, 0] * dx[:, 2], dx[:, 1] * dx[:, 2], dx]
    w = np.exp(-0.5 * (4 * np.linalg.norm(dx, axis=1) / BETA_STAR) ** 2); w[w <= WTHRES] = 0.0
    D = np.linalg.solve(A.T @ (w[:, None] * A), A.T * w)
    row = D[0] + D[1] + D[2]
    acc = {}
    for k, v in zip(keys, row):
        if k != (i, 0):
            acc[k] = acc.get(k, 0.0) + v
    acc[(i, 0)] = -sum(acc.values())
    real[i] = acc
# the extended operator: real rows + ghost rows (reflections); symmetrise over all pairs that occur
ext = {}
for i, acc in real.items():
    for k, v in acc.items():
        ext[((i, 0), k)] = v
masks = set(k[1] for acc in real.values() for k in acc) - {0}
for i, acc in real.items():
    for m in masks:
        # the ghost (i, m) exists if some real row references it
        pass
refd = set(k for acc in real.values() for k in acc if k[1] != 0)
for (i, m) in refd:
    for (j, mj), v in real[i].items():
        ext[((i, m), (j, mj ^ m))] = v
sym = {}
for (a, b), v in ext.items():
    sym[(a, b)] = 0.5 * (v + ext.get((b, a), 0.0))
    if (b, a) not in ext:
        sym[(b, a)] = sym[(a, b)]
# fold the ghost columns of the real rows, diagonal = minus the row sum
L = sp.lil_matrix((nGP, nGP))
for (a, b), v in sym.items():
    if a[1] != 0 or a == b:
        continue
    L[a[0], b[0]] += v
L = L.tocsr(); L.setdiag(0.0); L.eliminate_zeros()
L = (L - sp.diags(np.asarray(L.sum(1)).ravel())).tocsr(); L.sort_indices()
q = (pos[:nGP] ** 2).sum(1); v = L @ q
dist = np.minimum(pos[:nGP] + 1, 1 - pos[:nGP]).min(1)
print(f"nGP {nGP}: nnz {L.nnz}; |row sum| max {np.abs(L.sum(1)).max():.1e}, |column sum| max {np.abs(L.sum(0)).max():.1e}")
print("lap(r^2) (exact 6; r^2 has a flux through the walls, so only the inner bands count):",
      " ".join(f"[{a:.2f},{b:.2f}): {v[(dist >= a) & (dist < b)].min():.3f}..{v[(dist >= a) & (dist < b)].max():.3f}" for a, b in ((0.2, 0.4), (0.4, 9))))
# the manufactured c (zero flux): consistency to the walls
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "abaqus_living"))
import mms as M  # noqa: E402
PI = np.pi
x, y, z = pos[:nGP].T
cm = M.c_m(x, y, z, 0.0)
A_ = 0.3; s_ = np.sin(PI * (y + 1) / 4) ** 2; g_ = 1 + 0.5 * np.cos(PI * (x + 1)) * np.cos(PI * (z + 1))
lap_ex = -A_ * ((PI / 4) ** 2 * 2 * np.cos(PI * (y + 1) / 2) * g_ + s_ * 2 * (-0.5 * PI ** 2 * np.cos(PI * (x + 1)) * np.cos(PI * (z + 1))))
res = L @ cm - lap_ex
print("lap(c_m) error:", " ".join(f"[{a:.2f},{b:.2f}): {np.abs(res[(dist >= a) & (dist < b)]).max():.2e}" for a, b in ((0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 9))))
with open(d / "lap_cons.txt", "w") as f:
    f.write(f"{nGP} {L.nnz}\n")
    p1 = L.indptr + 1
    for k in range(0, len(p1), 10):
        f.write(" ".join(str(x) for x in p1[k:k + 10]) + "\n")
    for j, vv in zip(L.indices + 1, L.data):
        f.write(f"{j} {vv:.17e}\n")
print("wrote", d / "lap_cons.txt")
