"""A-42: lap_cons.txt for LAP_FILE = 1 (build F:\\biofilm_upf_living7) from nem_dump.txt of a run on the same mesh:
the NEM Laplacian of the integration points with the surface points eliminated by the SURF_NEAREST rule, made
conservative (symmetric part, diagonal = -(sum of the off-diagonal row)), so rows and columns sum to zero.

  python build_lap_cons_a42.py <dir with nem_dump.txt>        -> <dir>\\lap_cons.txt

Format: nGP nnz / row pointers (1-based, nGP + 1 values, 10 per line) / one "column value" pair per entry."""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

d = Path(sys.argv[1])
with open(d / "nem_dump.txt") as f:
    nGP, nTot, nnz = (int(x) for x in f.readline().split())
    pos = np.array([f.readline().split()[1:] for _ in range(nTot)], float)
    rows = np.loadtxt(f, max_rows=nnz, usecols=(0, 1, 5))
r, c = rows[:, 0].astype(int) - 1, rows[:, 1].astype(int) - 1
LAP = sp.csr_matrix((rows[:, 2], (r, c)), shape=(nTot, nTot))
S = sp.lil_matrix((nTot - nGP, nGP))
for i in range(nGP, nTot):
    js = np.array([j for j in LAP.indices[LAP.indptr[i]:LAP.indptr[i + 1]] if j < nGP])
    d2 = ((pos[js] - pos[i]) ** 2).sum(1)
    w = (d2 <= d2.min() * 1.0001).astype(float)
    S[i - nGP, js] = w / w.sum()
Le = (LAP[:nGP, :nGP] + LAP[:nGP, nGP:] @ S.tocsr()).tocsr()
Ls = (0.5 * (Le + Le.T)).tolil()
Ls.setdiag(0.0)
Ls = Ls.tocsr()
Lc = (Ls - sp.diags(np.asarray(Ls.sum(1)).ravel())).tocsr()
Lc.eliminate_zeros(); Lc.sort_indices()
print(f"nGP {nGP}: nnz {Lc.nnz} (NEM {Le.nnz}); |row sum| max {np.abs(Lc.sum(1)).max():.1e}, "
      f"|column sum| max {np.abs(Lc.sum(0)).max():.1e}")
q = (pos[:nGP] ** 2).sum(1); inner = np.abs(pos[:nGP]).max(1) < 0.5
v = Lc @ q
print(f"lap(r^2) at inner points {v[inner].min():.4f} .. {v[inner].max():.4f} (exact 6)")
with open(d / "lap_cons.txt", "w") as f:
    f.write(f"{nGP} {Lc.nnz}\n")
    ptr = Lc.indptr + 1
    for k in range(0, len(ptr), 10):
        f.write(" ".join(str(x) for x in ptr[k:k + 10]) + "\n")
    for j, val in zip(Lc.indices + 1, Lc.data):
        f.write(f"{j} {val:.17e}\n")
print("wrote", d / "lap_cons.txt")
