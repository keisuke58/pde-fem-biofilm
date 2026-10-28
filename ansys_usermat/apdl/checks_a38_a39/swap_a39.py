"""T* = 0.25 -> 0.5 with the Python operators (zero wall gradient), c held fixed: phi and c from Python or ANSYS."""
import glob, json, sys
import numpy as np

N, H = 8, 0.25
R, K, BETA, PEN = 10.0, 1.0, 0.02 + 10.0 * H / 4, 5.0
z = np.load(glob.glob(sys.argv[1])[0])
d = json.load(open(sys.argv[2]))
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]


def grid(v):
    g = np.full((N,) * 3, np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


def grad(u):
    out = []
    for a in range(3):
        g = np.gradient(u, H, axis=a)
        for j in (0, -1):
            e = [slice(None)] * 3; e[a] = j; g[tuple(e)] = 0.0
        out.append(g)
    return out


def lap(u):
    out = np.zeros_like(u)
    for a in range(3):
        up = np.roll(u, -1, a) - u; dn = np.roll(u, 1, a) - u
        e = [slice(None)] * 3; e[a] = -1; up[tuple(e)] = 0
        e[a] = 0; dn[tuple(e)] = 0
        out += up + dn
    return out / H ** 2


def front(phi, c):
    g = grad(c); gp = grad(phi)
    nrm = np.sqrt(sum(x ** 2 for x in g) + 1e-6)
    return -sum(R * c / (K + c) * g[a] / nrm * gp[a] for a in range(3))


def step(ph, c, T=0.25, n=20000):
    for _ in range(n):
        ph = ph + T / n * (BETA * lap(ph) + front(ph, c) - PEN * (np.maximum(ph - 1, 0) + np.minimum(ph, 0)))
    return ph


k = d["t"].index(0.25)
P = {"Py": z["phi_0.25"], "AN": grid(d["phi"][k])}
C = {"Py": z["c_0.25"], "AN": grid(d["nut"][k])}
for nm in ("Py", "AN"):
    print(f"{nm} at 0.25: phi mean {P[nm].mean():.4f} min {P[nm].min():.3f} max {P[nm].max():.3f}; "
          f"c mean {C[nm].mean():.4f} min {C[nm].min():.3f} max {C[nm].max():.3f}")
print("phi along y at x=z=0.125 Py", np.round(P["Py"][4, :, 4], 3), " AN", np.round(P["AN"][4, :, 4], 3))
print("c   along y at x=z=0.125 Py", np.round(C["Py"][4, :, 4], 3), " AN", np.round(C["AN"][4, :, 4], 3))
for p in ("Py", "AN"):
    for c in ("Py", "AN"):
        ph = step(P[p], C[c])
        print(f"phi {p}, c {c}: mean {P[p].mean():.4f} -> {ph.mean():.4f}; mean front term {front(P[p], C[c]).mean():+.4f} /T*")
