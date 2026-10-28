"""A-44: why the symmetric conservative Laplacian loses accuracy at 24^3, and an alternative that keeps it.

L_e   = NEM Laplacian of the integration points, surface points eliminated (SURF_NEAREST rule); rows sum to 0.
sym   = symmetric part with the diagonal from the row sums (A-42): rows and columns sum to 0.
wrow  = diag(w) L_e with w the left null vector of L_e (w^T L_e = 0), scaled to 1 inside: rows and columns sum to 0,
        each row only scaled by w_i, so it stays as accurate as L_e where w = 1 (inside); rows near the faces are
        scaled by w_i (diffusion there weaker or stronger by |w_i - 1|).

  python lap_cons_check_a44.py <dir with nem_dump.txt> [--write wrow|sym]"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

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
x = pos[:nGP]
L1 = np.abs(x).max(1)
inner = L1 < 0.5
q = (x ** 2).sum(1)


def report(name, L):
    v = L @ q
    cs = np.asarray(L.sum(0)).ravel(); rs = np.asarray(L.sum(1)).ravel()
    asym = abs(L - L.T)
    print(f"{name:6s}: lap(r^2) inner {v[inner].min():.4f} .. {v[inner].max():.4f}, all {np.percentile(v, 1):.3f} .. "
          f"{np.percentile(v, 99):.3f} (1-99 %); |row sum| {np.abs(rs).max():.1e}, |col sum| {np.abs(cs).max():.1e}; "
          f"asymmetry inner rows {asym[inner].sum() / abs(L)[inner].sum():.2e}")


report("L_e", Le)
Ls = (0.5 * (Le + Le.T)).tolil(); Ls.setdiag(0.0); Ls = Ls.tocsr()
Lsym = (Ls - sp.diags(np.asarray(Ls.sum(1)).ravel())).tocsr()
report("sym", Lsym)
# left null vector: L_e^T w = 0
vals, vecs = sla.eigs(Le.T, k=1, sigma=1e-8, which="LM")
w = np.real(vecs[:, 0]); w = w / np.median(w[inner])   # 1 inside: the inner rows stay as they are
print(f"w: min {w.min():.4f}, max {w.max():.4f}; inner {w[inner].min():.4f} .. {w[inner].max():.4f}; "
      f"near faces (|x|max > 0.75) {w[L1 > 0.75].min():.4f} .. {w[L1 > 0.75].max():.4f}; |L_e^T w| {np.abs(Le.T @ w).max():.1e}")
Lw = (sp.diags(w) @ Le).tocsr()
report("wrow", Lw)
if "--write" in sys.argv:
    kind = sys.argv[sys.argv.index("--write") + 1]
    Lc = {"wrow": Lw, "sym": Lsym}[kind]
    Lc.eliminate_zeros(); Lc.sort_indices()
    with open(d / "lap_cons.txt", "w") as f:
        f.write(f"{nGP} {Lc.nnz}\n")
        ptr = Lc.indptr + 1
        for k in range(0, len(ptr), 10):
            f.write(" ".join(str(v) for v in ptr[k:k + 10]) + "\n")
        for j, val in zip(Lc.indices + 1, Lc.data):
            f.write(f"{j} {val:.17e}\n")
    print("wrote", d / "lap_cons.txt", "(", kind, ")")
