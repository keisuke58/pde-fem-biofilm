"""NEM operators (front_dump.txt) vs finite differences on the same element means: Laplacian and gradient,
interior and wall cells, per call (8^3, 8 Gauss points per element averaged)."""
import json, sys
import numpy as np

N, H = 8, 0.25
cols = ["id", "t", "phi", "gpx", "gpy", "gpz", "gcx", "gcy", "gcz", "c", "lap", "ori", "growth", "pen", "dt", "phinew"]
A = np.loadtxt(sys.argv[1])
D = {k: A[:, i] for i, k in enumerate(cols)}
d = json.load(open(sys.argv[2]))
el = np.array(d["elem"])
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]
n = int(D["id"].max()); gp = n // len(el)
xc = -0.875 + H * np.arange(N)
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
wall = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)


def grid(v):
    g = np.full((N,) * 3, np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


def lap7(u):
    out = np.zeros_like(u)
    for a in range(3):
        up = np.roll(u, -1, a) - u; dn = np.roll(u, 1, a) - u
        e = [slice(None)] * 3; e[a] = -1; up[tuple(e)] = 0
        e[a] = 0; dn[tuple(e)] = 0
        out += up + dn
    return out / H ** 2


def gmirror(u, a):
    g = np.gradient(u, H, axis=a)
    for j in (0, -1):
        e = [slice(None)] * 3; e[a] = j; g[tuple(e)] *= 0.5
    return g


print("call  <phi>   lap: slope NEM/FD interior, wall;  sum NEM, sum FD | grad phi: slope NEM/FD(mirror) interior, wall")
for k in range(len(A) // n):
    s = slice(k * n, (k + 1) * n)
    e = (D["id"][s].astype(int) - 1) // gp
    per = {}
    for c in ("phi", "lap", "gpx", "gpy", "gpz"):
        v = np.zeros(len(el)); np.add.at(v, e, D[c][s] / gp); per[c] = grid(v)
    L7 = lap7(per["phi"])
    fit = lambda x, y: float((x * y).sum() / (x * x).sum())
    li, lw = fit(L7[~wall], per["lap"][~wall]), fit(L7[wall], per["lap"][wall])
    gF = np.concatenate([gmirror(per["phi"], a)[None] for a in range(3)])
    gN = np.concatenate([per[c][None] for c in ("gpx", "gpy", "gpz")])
    gi = fit(gF[:, ~wall].ravel(), gN[:, ~wall].ravel()); gw = fit(gF[:, wall].ravel(), gN[:, wall].ravel())
    print(f"{k:3d}  {np.nanmean(per['phi']):.4f}   {li:.2f}, {lw:.2f};   {per['lap'].sum():+8.1f}, {L7.sum():+6.1f} |   {gi:.2f}, {gw:.2f}")
