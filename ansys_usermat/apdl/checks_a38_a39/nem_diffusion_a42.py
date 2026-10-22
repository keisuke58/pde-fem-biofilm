"""A-42: diffusion only, d phi/dt = beta L phi on the ANSYS integration points (8^3, 4096 points) from the separate
example's start (seed + corner blocks phi = 1), with the NEM operator and two surface rules, and with the NEM operator
made conservative (symmetric part, diagonal from the row sums). The mean over the points (equal Gauss weights) should
stay constant with zero-flux walls.

  python nem_diffusion_a42.py nem_dump.txt front_dump.txt"""
import sys
import numpy as np
import scipy.sparse as sp

sys.argv += [] if len(sys.argv) > 2 else []
exec(open(__file__.replace("nem_diffusion_a42.py", "nem_replay_a42.py")).read().split("# --- 1.")[0])
A = np.loadtxt(sys.argv[2])
B = A[:nGP]; B = B[np.argsort(B[:, 0])]
phi0 = B[:, 2]                                  # the start state (first call)
beta, T = 0.645, 0.25


def conservative(Le):
    Ls = (0.5 * (Le + Le.T)).tolil()
    Ls.setdiag(0.0)
    Ls = Ls.tocsr()
    return (Ls - sp.diags(np.asarray(Ls.sum(1)).ravel())).tocsr()


ops = {"NEM, inverse-square surface (old)": eff(LAP, surface_rule("inv2")),
       "NEM, nearest surface (SURF_NEAREST)": eff(LAP, surface_rule("nearest"))}
ops["NEM nearest, made conservative"] = conservative(ops["NEM, nearest surface (SURF_NEAREST)"])
lam = {k: float(np.abs(sp.linalg.eigs(v, k=1, which="LM", return_eigenvectors=False))[0]) for k, v in ops.items()}
for k, Lo in ops.items():
    dt = 0.5 / (beta * lam[k]); n = int(np.ceil(T / dt)); dt = T / n
    ph = phi0.copy(); means = [ph.mean()]
    for i in range(n):
        ph = ph + dt * beta * (Lo @ ph)
        if (i + 1) % (n // 5) == 0:
            means.append(ph.mean())
    cs = np.asarray(Lo.sum(0)).ravel()
    print(f"{k:38s}: |colsum| max {np.abs(cs).max():6.1f}; mean phi " + " -> ".join(f"{m:.4f}" for m in means)
          + f"; max {ph.max():.3f}")
# consistency of the conservative operator: Laplacian of x^2 + y^2 + z^2 (= 6) at points away from the walls
q = (pos[:nGP] ** 2).sum(1)
inner = np.abs(pos[:nGP]).max(1) < 0.5
for k, Lo in ops.items():
    v = Lo @ q
    print(f"{k:38s}: lap(r^2) inside mean {v[inner].mean():.3f}, range {v[inner].min():.3f} .. {v[inner].max():.3f} (exact 6)")
