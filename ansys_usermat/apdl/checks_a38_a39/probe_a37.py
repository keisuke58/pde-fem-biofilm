"""Apply the Python operators (central front term, beta + r h/4 Laplacian, penalty) to the ANSYS fields at
T* = 0.0025 and integrate to 0.0275 with c held fixed; compare with the ANSYS field at 0.0275."""
import json, sys
import numpy as np

f = sys.argv[1]
d = json.load(open(f))
N, H = 8, 0.25
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]


def grid(v):
    g = np.full((N, N, N), np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


R, K, BETA, PEN = 10.0, 1.0, 0.02 + 10.0 * H / 4, 5.0
if "--front0" in sys.argv:
    R, BETA = 0.0, 0.02


def lap(u):
    out = np.zeros_like(u)
    for a in range(3):
        up = np.roll(u, -1, a) - u; dn = np.roll(u, 1, a) - u
        e = [slice(None)] * 3; e[a] = -1; up[tuple(e)] = 0
        e[a] = 0; dn[tuple(e)] = 0
        out += up + dn
    return out / H ** 2


def front(phi, c):
    g = np.gradient(c, H); gp = np.gradient(phi, H)
    nrm = np.sqrt(sum(x ** 2 for x in g) + 1e-6)
    return -sum(R * c / (K + c) * g[a] / nrm * gp[a] for a in range(3))


k0 = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 0; k1 = k0 + 1
phi0, phi1 = grid(d["phi"][k0]), grid(d["phi"][k1])
c0 = grid(d["nut"][k0])
dt_out = d["t"][k1] - d["t"][k0]
n = int(2000 * dt_out / 0.025) + 1
ph = phi0.copy()
for _ in range(n):
    rhs = BETA * lap(ph) + front(ph, c0) - PEN * (np.maximum(ph - 1, 0) + np.minimum(ph, 0))
    ph = ph + dt_out / n * rhs
seed = np.zeros(N ** 3, bool); el = np.array(d["elem"])
seedg = grid(np.isin(el, d["seed"]).astype(float)) > 0.5
print(f"T* {d['t'][k0]} -> {d['t'][k1]}, R = {R}, beta = {BETA}")
print(f"mean phi: start {phi0.mean():.4f}  ANSYS {phi1.mean():.4f}  Python ops on ANSYS field {ph.mean():.4f}")
print(f"seed mean: start {phi0[seedg].mean():.4f}  ANSYS {phi1[seedg].mean():.4f}  Python ops {ph[seedg].mean():.4f}")
print(f"max: ANSYS {phi1.max():.4f}  Python ops {ph.max():.4f}")
# one seed-edge element: the seed cell with the largest ANSYS loss, and its outer neighbour
g = np.gradient(c0, H); gp = np.gradient(phi0, H); nrm = np.sqrt(sum(x ** 2 for x in g) + 1e-6)
gpn = np.sqrt(sum(x ** 2 for x in gp))
for name, m in (("seed", seedg), ("outside", ~seedg & (phi0 < 0.01))):
    sel = np.argwhere(m)
    dd = (phi1 - phi0)[m]
    j = tuple(sel[np.argmax(np.abs(dd))])
    nd = sum(g[a][j] * gp[a][j] for a in range(3)) / max(nrm[j] * gpn[j], 1e-12)
    print(f"{name} cell {j} (x,y,z = {[-0.875 + H * q for q in j]}): phi {phi0[j]:.3f}->{phi1[j]:.3f} ANSYS, "
          f"{ph[j]:.3f} Python ops; |grad phi| {gpn[j]:.3f}, n_phi.n_c {nd:+.3f}, c {c0[j]:.4f}, c/(k+c) {c0[j]/(K+c0[j]):.4f}, "
          f"|grad c| {nrm[j]:.2e}, grad c {[f'{g[a][j]:+.2e}' for a in range(3)]}, lap {lap(phi0)[j]:+.2f}, "
          f"front {front(phi0, c0)[j]:+.2f}, beta lap {BETA * lap(phi0)[j]:+.2f}")
print("c range", np.nanmin(c0), np.nanmax(c0))
