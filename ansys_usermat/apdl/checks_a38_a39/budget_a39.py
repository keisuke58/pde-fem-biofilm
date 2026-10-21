"""Budget of the mean phi from front_dump.txt: d<phi>/dt split into beta lap, front, penalty, local growth; the
Laplacian sum over all points (0 for a zero-flux, conservative operator)."""
import sys
import numpy as np

BETA = float(sys.argv[2]) if len(sys.argv) > 2 else 0.645
cols = ["id", "t", "phi", "gpx", "gpy", "gpz", "gcx", "gcy", "gcz", "c", "lap", "ori", "growth", "pen", "dt", "phinew"]
A = np.loadtxt(sys.argv[1])
D = {k: A[:, i] for i, k in enumerate(cols)}
n = int(D["id"].max())
blocks = [slice(i, i + n) for i in range(0, len(A) - n + 1, n)]
print("rows", len(A), "calls", len(blocks))
print("call  mean phi   d<phi>/dt   beta*<lap>   <front>   <pen>   rest   sum lap / point (rows sorted by call)")
for k, s in enumerate(blocks):
    phi, new, dt = D["phi"][s], D["phinew"][s], D["dt"][s][0]
    rate = (new.mean() - phi.mean()) / dt
    lap = BETA * D["lap"][s].mean(); fr = (-D["ori"][s] * D["growth"][s]).mean(); pen = D["pen"][s].mean()
    print(f"{k:3d}  {phi.mean():.4f}   {rate:+.4f}     {lap:+.4f}     {fr:+.4f}   {pen:+.4f}  {rate - lap - fr - pen:+.4f}"
          f"   {D['lap'][s].mean():+.3f}")
