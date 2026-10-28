"""A-42: replay of the ANSYS phi operator from its own NEM matrices (nem_dump.txt, NEM_DUMP = 1) on the values of
front_dump.txt (FRONT_DUMP), and the leak of the effective Laplacian through the walls for several surface-point rules.

  python nem_replay_a42.py nem_dump.txt front_dump.txt [beta]

Surface points (rows nGP+1 .. nTot) take phi from their integration-point neighbours (S); the integration points use
L_eff = L_GG + L_GS S. A zero-flux, conservative operator has zero column sums (equal point volumes)."""
import sys
import numpy as np
import scipy.sparse as sp

nem, fdump = sys.argv[1], sys.argv[2]
with open(nem) as f:
    nGP, nTot, nnz = (int(x) for x in f.readline().split())
    pos = np.array([f.readline().split()[1:] for _ in range(nTot)], float)
    rows = np.loadtxt(f, max_rows=nnz)
r, c = rows[:, 0].astype(int) - 1, rows[:, 1].astype(int) - 1
DX, DY, DZ, LAP = (sp.csr_matrix((rows[:, k], (r, c)), shape=(nTot, nTot)) for k in (2, 3, 4, 5))
print(f"nGP {nGP}, surface points {nTot - nGP}, nnz {nnz}; neighbours per row {nnz / nTot:.1f}")
G, S_ = slice(0, nGP), slice(nGP, nTot)


def surface_rule(kind):
    """S: (nTot - nGP) x nGP, phi_surface = S phi_GP"""
    out = sp.lil_matrix((nTot - nGP, nGP))
    for i in range(nGP, nTot):
        lo, hi = LAP.indptr[i], LAP.indptr[i + 1]
        js = np.array([j for j in LAP.indices[lo:hi] if j < nGP])
        d2 = ((pos[js] - pos[i]) ** 2).sum(1)
        if kind == "inv2":
            w = 1.0 / np.maximum(d2, 1e-30)
        elif kind == "nearest":
            w = (d2 <= d2.min() * 1.0001).astype(float)
        out[i - nGP, js] = w / w.sum()
    return out.tocsr()


def eff(Op, S):
    return (Op[G, G] + Op[G, S_] @ S).tocsr()


# --- 1. replay against the dump (SURF_NEAREST run)
cols = ["id", "t", "phi", "gpx", "gpy", "gpz", "gcx", "gcy", "gcz", "c", "lap", "ori", "growth", "pen", "dt", "phinew"]
A = np.loadtxt(fdump)
S = surface_rule("nearest")
Le, Xe, Ye, Ze = eff(LAP, S), eff(DX, S), eff(DY, S), eff(DZ, S)
beta = float(sys.argv[3]) if len(sys.argv) > 3 else 0.645
for k in range(len(A) // nGP):
    B = A[k * nGP:(k + 1) * nGP]
    order = np.argsort(B[:, 0]); B = B[order]
    phi = B[:, 2]
    lap = Le @ phi; gx = Xe @ phi
    rest = (B[:, 15] - phi) / B[0, 14] - (-B[:, 11] * B[:, 12] + beta * B[:, 10] + B[:, 13])
    print(f"call {k}: lap max |replay - ANSYS| {np.abs(lap - B[:, 10]).max():.2e} (max |lap| {np.abs(B[:, 10]).max():.1f}); "
          f"grad x {np.abs(gx - B[:, 3]).max():.2e}; phi update rest (k_local Bioloc) {rest.min():+.2e} .. {rest.max():+.2e}")

# --- 2. leak: column sums of L_eff for the surface rules
xyz = pos[:nGP]
L1 = np.abs(xyz).max(1)
for kind in ("inv2", "nearest"):
    Le = eff(LAP, surface_rule(kind))
    cs = np.asarray(Le.sum(0)).ravel()
    rs = np.asarray(Le.sum(1)).ravel()
    print(f"{kind:8s}: row sums |max| {np.abs(rs).max():.2e} (constants exact); column sums: total {cs.sum():+.1f}, "
          f"min {cs.min():+.1f}, max {cs.max():+.1f}; points with |colsum| > 1: {(np.abs(cs) > 1).sum()} of {nGP}")
    for nm, m in (("near a face", (L1 > 0.75)), ("inside", L1 <= 0.75)):
        print(f"           {nm:12s}: sum of column sums {cs[m].sum():+9.1f} over {m.sum()} points")
