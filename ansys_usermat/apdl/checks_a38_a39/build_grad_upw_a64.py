"""A-64: one-sided (upwind) gradients of phi at the integration points from the NEM dump (nem_dump.txt), for the front
term in Oliver's element (FRONT_UPWIND = 1, build living11): for each axis a and each side s, a weighted least-squares
plane fit over the neighbours j of the NEM stencil with x_j,a - x_i,a on side s (and the point itself); the a-component
of that fit is the one-sided derivative. The element then takes the side opposite to the sign of the front velocity.
Surface points (nGP+1..nTot) carry the value of their nearest integration point (SURF_NEAREST), as in lap_cons.

  python build_grad_upw_a64.py <dir with nem_dump.txt>        -> <dir>\\grad_upw.txt

Format: line 1 nGP; line 2 the six nnz (x-, x+, y-, y+, z-, z+); then per operator the row pointers (1-based, 10 per
line) and one "column value" pair per entry (columns may be surface points, > nGP).
Check: the derivative of x_a (exact 1) and of x_a^2 (exact 2 x_a) at the inner points."""
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
with open(d / "nem_dump.txt") as f:
    nGP, nTot, nnz = (int(x) for x in f.readline().split())
    pos = np.array([f.readline().split()[1:] for _ in range(nTot)], float)
    rows = np.loadtxt(f, max_rows=nnz, usecols=(0, 1))
r, c = rows[:, 0].astype(int) - 1, rows[:, 1].astype(int) - 1
nb = [[] for _ in range(nTot)]
for i, j in zip(r, c):
    if i != j:
        nb[i].append(j)
h = np.min([np.linalg.norm(pos[j] - pos[0]) for j in nb[0]])      # smallest neighbour distance
ops = {}
names = ["xm", "xp", "ym", "yp", "zm", "zp"]
for a in range(3):
    for s, nm in ((-1, names[2 * a]), (1, names[2 * a + 1])):
        ptr, col, val = [1], [], []
        for i in range(nGP):
            js = [j for j in nb[i] if s * (pos[j, a] - pos[i, a]) > 0.1 * h]
            dx = pos[js] - pos[i]
            A = np.c_[dx, np.ones(len(js))]
            A = np.r_[A, [[0, 0, 0, 1.0]]]                       # the point itself
            dist = np.r_[np.linalg.norm(dx, axis=1), 0.0]
            w = np.exp(-(dist / (1.5 * h)) ** 2)
            D = np.linalg.solve(A.T @ (w[:, None] * A), (A.T * w))   # 4 x (len(js)+1)
            g = D[a]                                             # row of d/dx_a
            for j, v in zip(js + [i], g):
                if abs(v) > 1e-14:
                    col.append(j + 1); val.append(v)
            ptr.append(len(col) + 1)
        ops[nm] = (np.array(ptr), np.array(col), np.array(val))
# check at the inner points
inner = np.abs(pos[:nGP]).max(1) < 0.5
for nm, (ptr, col, val) in ops.items():
    a = "xyz".index(nm[0])
    for fn, ex, lab in ((pos[:, a], np.ones(nGP), "x_a"), (pos[:, a] ** 2, 2 * pos[:nGP, a], "x_a^2")):
        out = np.zeros(nGP)
        for i in range(nGP):
            sl = slice(ptr[i] - 1, ptr[i + 1] - 1)
            out[i] = np.dot(val[sl], fn[col[sl] - 1])
        err = np.abs(out - ex)[inner]
        print(f"{nm}: d({lab})/dx_{'xyz'[a]} at inner points, max error {err.max():.2e}, mean {err.mean():.2e}")
with open(d / "grad_upw.txt", "w") as f:
    f.write(f"{nGP}\n" + " ".join(str(len(ops[nm][1])) for nm in names) + "\n")
    for nm in names:
        ptr, col, val = ops[nm]
        for k in range(0, len(ptr), 10):
            f.write(" ".join(str(x) for x in ptr[k:k + 10]) + "\n")
        for j, v in zip(col, val):
            f.write(f"{j} {v:.17e}\n")
print("wrote", d / "grad_upw.txt")
