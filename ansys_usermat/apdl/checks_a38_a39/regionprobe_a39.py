"""Python operators on the ANSYS early fields, one 0.025 interval, by region, for both wall gradients;
the term split (diffusion / front) in the corner blocks."""
import json, sys
import numpy as np

N, H = 8, 0.25
R, K, BETA, PEN = 10.0, 1.0, 0.02 + 10.0 * H / 4, 5.0
d = json.load(open(sys.argv[1]))
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]


def grid(v):
    g = np.full((N,) * 3, np.nan)
    g[ix[0], ix[1], ix[2]] = v
    return g


def grad(u, wall):
    out = []
    for a in range(3):
        g = np.gradient(u, H, axis=a)
        if wall == "zero":
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


def front(phi, c, wall):
    g = grad(c, wall); gp = grad(phi, wall)
    nrm = np.sqrt(sum(x ** 2 for x in g) + 1e-6)
    return -sum(R * c / (K + c) * g[a] / nrm * gp[a] for a in range(3))


xc = -0.875 + H * np.arange(N)
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
seed = grid(np.isin(np.array(d["elem"]), d["seed"]).astype(float)) > 0.5
blk = (grid(d["phi"][0]) > 0.99) & ~seed
wall = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)
regs = (("seed", seed), ("blocks", blk), ("wall other", wall & ~blk & ~seed), ("interior", ~wall & ~seed & ~blk))
ks = [int(a) for a in sys.argv[2:]] or [1, 2]
for k0 in ks:
    p0, p1 = grid(d["phi"][k0]), grid(d["phi"][k0 + 1])
    c0 = grid(d["nut"][k0 + (1 if k0 == 0 else 0)])      # k0 = 0: c of the next set (the start set is the IC)
    print(f"T* {d['t'][k0]} -> {d['t'][k0 + 1]}" + ("  (c from the next set)" if k0 == 0 else ""))
    print("   " + "  ".join(f"{nm}: AN {p1[m].sum():5.2f}" for nm, m in regs) + f"   (start {'  '.join(f'{p0[m].sum():5.2f}' for _, m in regs)})")
    for w in ("onesided", "zero"):
        ph = p0.copy(); n = 400
        for _ in range(n):
            ph = ph + (d["t"][k0 + 1] - d["t"][k0]) / n * (BETA * lap(ph) + front(ph, c0, w) - PEN * (np.maximum(ph - 1, 0) + np.minimum(ph, 0)))
        print(f"   {w:8s} " + "  ".join(f"{nm}: Py {ph[m].sum():5.2f}" for nm, m in regs)
              + f"   blocks: diffusion {(BETA * lap(p0))[blk].sum():+.1f}, front {front(p0, c0, w)[blk].sum():+.1f} /T*")
