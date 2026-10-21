"""Python 8^3 separate_ab with a different gradient at the walls, to test whether the Python/ANSYS gap is the
wall cells. variant: onesided (as the script, np.gradient), mirror (zero-flux ghost: (u1 - u0) / 2h), zero."""
import importlib.util, json, sys
import numpy as np

variant, repo, ansys_json = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ["x", "--nf", "8", "--front", "1", "--phi-min", "0.001", "--v2", "--ansys-bc", "--start-position",
            "--tend", "0.25", "--alpha0", "100", "--gab", "1", "--gc", "1", "--tab", "0.333333",
            "--snap-every", "0.025", "--artdiff"] + (["--dtc", "0.005"] if "DTC5" in __import__("os").environ else []) + [ "--only", "separate_ab"]
spec = importlib.util.spec_from_file_location("clc", f"{repo}/ansys_usermat/apdl/coupled_living_consumption.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
H = M.H


def grad_wall(u):
    out = []
    for a in range(3):
        g = np.gradient(u, H, axis=a)
        e0 = [slice(None)] * 3; e0[a] = 0
        e1 = [slice(None)] * 3; e1[a] = -1
        if variant == "mirror":
            g[tuple(e0)] *= 0.5; g[tuple(e1)] *= 0.5
        elif variant == "zero":
            g[tuple(e0)] = 0.0; g[tuple(e1)] = 0.0
        out.append(g)
    return out


if variant != "onesided":
    M.grad = grad_wall
M.CKPT_SEC = 1e9
hist, snaps = M.run("separate", True)
d = json.load(open(ansys_json))
N = 8
ix = [np.rint((np.array(d[k]) + 0.875) / H).astype(int) for k in ("cx", "cy", "cz")]
print(f"variant {variant}")
print("T*      Py mean   ANSYS mean   Py max  ANSYS max")
for k, T in enumerate(d["t"]):
    Tp = round(T - 0.0025, 3)
    if Tp not in snaps:
        continue
    g = np.full((N,) * 3, np.nan); g[ix[0], ix[1], ix[2]] = d["phi"][k]
    print(f"{Tp:.3f}  {snaps[Tp]['phi'].mean():.4f}    {np.nanmean(g):.4f}      {snaps[Tp]['phi'].max():.3f}   {np.nanmax(g):.3f}")
xc = -0.875 + H * np.arange(N)
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
gr = lambda v: (lambda g: (g.__setitem__((ix[0], ix[1], ix[2]), v), g)[1])(np.full((N,) * 3, np.nan))
seedm = gr(np.isin(np.array(d["elem"]), d["seed"]).astype(float)) > 0.5
p0 = gr(d["phi"][0]); blk = (p0 > 0.99) & ~seedm
wall = (np.abs(X) > 0.75) | (np.abs(Y) > 0.75) | (np.abs(Z) > 0.75)
regs = (("seed", seedm), ("blocks", blk), ("wall other", wall & ~blk & ~seedm), ("interior", ~wall & ~seedm & ~blk))
for k, T in enumerate(d["t"]):
    Tp = round(T - 0.0025, 3)
    if Tp not in (0.025, 0.1, 0.25):
        continue
    ap = gr(d["phi"][k]); pp = snaps[Tp]["phi"]
    print(f"T* {Tp}: " + "; ".join(f"{nm} Py {pp[m].sum():.1f} AN {ap[m].sum():.1f}" for nm, m in regs))
print("seed living", round(hist["living_seed"][-1], 4), "c seed", round(hist["c_seed"][-1], 4), "alpha-1 seed", hist["seed_am1"])
