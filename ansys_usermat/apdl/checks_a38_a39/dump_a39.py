"""front_dump.txt (USSFin, A-39) vs the Python operators on the same phi and c, element by element (8^3)."""
import json, sys
import numpy as np

N, H = 8, 0.25
R, K = 10.0, 1.0
cols = ["id", "t", "phi", "gpx", "gpy", "gpz", "gcx", "gcy", "gcz", "c", "lap", "ori", "growth", "pen", "dt", "phinew"]
A = np.loadtxt(sys.argv[1])
D = {k: A[:, i] for i, k in enumerate(cols)}
d = json.load(open(sys.argv[2]))                      # fields json: element centroids
el = np.array(d["elem"])
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]
times = np.unique(D["t"])
print("times", times[:12], "ids per time", [int((D["t"] == t).sum()) for t in times[:3]], "max id", int(D["id"].max()))
gp = int(D["id"].max()) // len(el)
print("points per element", gp)
xc = -0.875 + H * np.arange(N)
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
wall = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)


def grid(v):
    g = np.full((N,) * 3, np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


def pygrad(u, zero_wall):
    out = []
    for a in range(3):
        g = np.gradient(u, H, axis=a)
        if zero_wall:
            for j in (0, -1):
                e = [slice(None)] * 3; e[a] = j; g[tuple(e)] = 0.0
        out.append(g)
    return out


def pylap(u):
    out = np.zeros_like(u)
    for a in range(3):
        up = np.roll(u, -1, a) - u; dn = np.roll(u, 1, a) - u
        e = [slice(None)] * 3; e[a] = -1; up[tuple(e)] = 0
        e[a] = 0; dn[tuple(e)] = 0
        out += up + dn
    return out / H ** 2


for t in times[: int(sys.argv[3]) if len(sys.argv) > 3 else 3]:
    m = D["t"] == t
    ide = ((D["id"][m].astype(int) - 1) // gp)       # 0-based element index in USSFin order
    per = {k: np.zeros(len(el)) for k in cols}
    for k in cols:
        np.add.at(per[k], ide, D[k][m] / gp)
    phi, c = grid(per["phi"]), grid(per["c"])
    an = {k: grid(per[k]) for k in ("gpx", "gpy", "gpz", "gcx", "gcy", "gcz", "lap", "ori", "growth")}
    an_front = -an["ori"] * an["growth"]
    gpz_, gcz_ = pygrad(phi, True), pygrad(c, True)
    py_front = -sum(R * c / (K + c) * gcz_[a] / np.sqrt(sum(x ** 2 for x in gcz_) + 1e-6) * gpz_[a] for a in range(3))
    pl = pylap(phi)
    print(f"\nt = {t:.5f}: phi mean {np.nanmean(phi):.4f}, c min {np.nanmin(c):.4f} max {np.nanmax(c):.4f}")
    for nm, msk in (("wall", wall), ("interior", ~wall)):
        print(f"  {nm:8s}: sum front ANSYS {an_front[msk].sum():+9.2f}  Python(zero wall) {py_front[msk].sum():+9.2f};"
              f" sum lap ANSYS {an['lap'][msk].sum():+9.1f}  Python {pl[msk].sum():+9.1f}")
    # a corner block cell and a seed-edge cell
    for j in ((7, 7, 4), (0, 7, 3), (6, 7, 4), (4, 4, 5)):
        print(f"  cell {[float(xc[q]) for q in j]}: phi {phi[j]:.3f} c {c[j]:.4f} | grad phi AN "
              f"{np.round([an['gpx'][j], an['gpy'][j], an['gpz'][j]], 2)} Py {np.round([g[j] for g in pygrad(phi, False)], 2)}"
              f" | grad c AN {np.round([an['gcx'][j], an['gcy'][j], an['gcz'][j]], 4)} Py {np.round([g[j] for g in pygrad(c, False)], 4)}"
              f" | lap AN {an['lap'][j]:+.2f} Py {pl[j]:+.2f} | front AN {an_front[j]:+.3f} Py {py_front[j]:+.3f}")
