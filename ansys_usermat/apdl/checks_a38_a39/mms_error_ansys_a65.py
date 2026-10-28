"""A-65: L2 error of an ANSYS run with MMS_SRC = 1 (build living12) against the manufactured solution of abaqus_living/mms.py,
from the element means of fields_<k>.csv (post_fields_living_t0.05.mac): the element mean of phi_m over its 8 Gauss
points (centroid +- h/(2 sqrt 3)) against the element mean of SVAR 41 (phi) and SVAR 51 (c).

    python mms_error_ansys_a65.py <csv dir> n [k=4]
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE / "abaqus_living"))
import mms as M  # noqa: E402

d, n = Path(sys.argv[1]), int(sys.argv[2])
k = int(sys.argv[3]) if len(sys.argv) > 3 else 4
f = d / f"fields_{k}.csv"
first = open(f).readline()
A = np.loadtxt(f, delimiter=",", skiprows=0 if first.strip()[0].isdigit() else 1)
t, phi, c = A[0, 1], A[:, 2], A[:, 3]
x, y, z = A[:, 9], A[:, 10], A[:, 11]
g = 2.0 / n / (2 * np.sqrt(3.0))
pm = np.zeros_like(phi); cm = np.zeros_like(phi)
for sx in (-1, 1):
    for sy in (-1, 1):
        for sz in (-1, 1):
            pm += M.phi_m(x + sx * g, y + sy * g, z + sz * g, t) / 8
            cm += M.c_m(x + sx * g, y + sy * g, z + sz * g, t) / 8
ep, ec = phi - pm, c - cm
print(f"{d}: {n}^3, T* {t:.4g}: L2 error phi {np.sqrt((ep ** 2).mean()):.3e} (max {np.abs(ep).max():.3e}), "
      f"c {np.sqrt((ec ** 2).mean()):.3e} (max {np.abs(ec).max():.3e})")
