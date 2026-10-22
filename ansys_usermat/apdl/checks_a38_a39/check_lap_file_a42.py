"""A-42: does the run use lap_cons.txt? The Laplacian in front_dump.txt against lap_cons.txt applied to the dumped phi;
and the sum over the points (0 for the conservative operator).

  python check_lap_file_a42.py <dir with lap_cons.txt and front_dump.txt>"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

d = Path(sys.argv[1])
with open(d / "lap_cons.txt") as f:
    n, nnz = (int(x) for x in f.readline().split())
    ptr = []
    while len(ptr) < n + 1:
        ptr += [int(x) for x in f.readline().split()]
    cv = np.loadtxt(f, max_rows=nnz)
Lc = sp.csr_matrix((cv[:, 1], cv[:, 0].astype(int) - 1, np.array(ptr) - 1), shape=(n, n))
A = np.loadtxt(d / "front_dump.txt")
for k in range(len(A) // n):
    B = A[k * n:(k + 1) * n]; B = B[np.argsort(B[:, 0])]
    lap = Lc @ B[:, 2]
    print(f"call {k}: max |lap_cons phi - dumped lap| {np.abs(lap - B[:, 10]).max():.2e} (max |lap| {np.abs(B[:, 10]).max():.1f});"
          f" sum of dumped lap {B[:, 10].sum():+.3f}")
