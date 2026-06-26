#!/usr/bin/env python3
"""
tooth_pde3d.py -- ULTRA-RIGOROUS: solve the Klempt growth reaction-diffusion
PDE directly in 3D on the real tooth biofilm tet mesh (no 2D->3D mapping, no
imposed depth ramp). alpha(node) comes straight from the PDE -> UMAT -> sigma.

Physics (P1 FEM on the C3D4 shell, physical mm units):
  dphi/dt = D_PHI lap(phi) + mu_eff phi(1-phi) c/(K_M+c)
  dc/dt   = D_C   lap(c)   - GAMMA_C phi c/(K_M+c)   ; c=1 Dirichlet on OUTER (oral) surface
  dalpha/dt = k_alpha_eff phi c/(K_M+c)
Nutrient enters from the oral (outer, layer 4) surface; attachment (inner) no-flux.
Penetration length sqrt(D_C/GAMMA_C) ~ 0.045 mm < 0.2 mm shell -> nutrient-limited.
"""
import sys, time
import numpy as np
import scipy.sparse as sp

D_PHI, D_C = 5e-4, 1e-2          # mm^2 / t*   (same kinetics as klempt_pde_multispecies)
GAMMA_C, K_M = 5.0, 0.3
K_ALPHA = np.array([1.0, 0.8, 0.4, 0.6, 0.3])
MU_SP = np.array([2.0, 1.6, 1.2, 1.4, 0.8])
NV = 1797           # nodes per shell layer
N_LAYERS_NODE = 5   # node layers (0=inner/attachment ... 4=outer/oral)
N_STEPS = 2000


def parse_mesh(inp):
    nodes, tets = [], []
    rd = None
    for ln in open(inp):
        s = ln.strip()
        low = s.lower()
        if low.startswith("*node") and "output" not in low:
            rd = "n"; continue
        if low.startswith("*element") and "output" not in low:
            rd = "e"; continue
        if s.startswith("*"):
            rd = None; continue
        if not s:
            continue
        if rd == "n":
            p = s.split(",")
            try:
                nodes.append([float(p[1]), float(p[2]), float(p[3])])
            except (ValueError, IndexError):
                rd = None
        elif rd == "e":
            p = s.split(",")
            try:
                tets.append([int(x) for x in p[1:5]])
            except (ValueError, IndexError):
                rd = None
    return np.array(nodes), np.array(tets, dtype=int) - 1  # 0-indexed


def assemble(nodes, tets):
    """P1 FEM: lumped mass Ml (N,), stiffness K (NxN sparse), for D=1."""
    N = len(nodes)
    Ml = np.zeros(N)
    rows, cols, vals = [], [], []
    for tet in tets:
        p = nodes[tet]                       # (4,3)
        J = (p[1:] - p[0]).T                 # (3,3)
        detJ = np.linalg.det(J)
        V = abs(detJ) / 6.0
        if V < 1e-18:
            continue
        # grads of barycentric coords: g0 = -(g1+g2+g3); [g1,g2,g3] = inv(J)^T
        Ginv = np.linalg.inv(J).T            # rows = grad lambda_{1,2,3}
        g = np.zeros((4, 3))
        g[1:] = Ginv
        g[0] = -Ginv.sum(0)
        Ke = V * (g @ g.T)                   # (4,4)
        for a in range(4):
            Ml[tet[a]] += V / 4.0
            for b in range(4):
                rows.append(tet[a]); cols.append(tet[b]); vals.append(Ke[a, b])
    K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    return Ml, K


def solve(Ml, K, phi_vec, outer_nodes, dt=2e-3, n_steps=1000):
    """ESTABLISHED biofilm: the 0.2mm conformal shell IS a dense, fully-colonized
    biofilm, so the phase-field density phi = 1 uniformly (the growth term
    mu*phi*(1-phi) vanishes at phi=1 -> phi stays 1). Nutrient c diffuses from
    the oral (outer) surface (Dirichlet c=1) and is consumed; alpha accumulates
    dalpha/dt = k_alpha*phi*monod(c). Backward-Euler (implicit) diffusion ->
    unconditionally stable. t_end = dt*n_steps = 2.0 matches the 2D headline.
    Returns (alpha[node], phi[node]=1)."""
    from scipy.sparse.linalg import factorized
    N = len(Ml)
    k_alpha = float(phi_vec @ K_ALPHA)
    M = sp.diags(Ml)
    A_c = (M + dt * D_C * K).tolil()
    for nidx in outer_nodes:                           # Dirichlet c=1 on oral surface
        A_c.rows[nidx] = [nidx]; A_c.data[nidx] = [1.0]
    solve_c = factorized(A_c.tocsc())
    phi = np.ones(N)                                   # established dense biofilm
    c = np.ones(N)
    alpha = np.zeros(N)
    for _ in range(n_steps):
        monod = c / (K_M + c)
        rhs_c = Ml * (c - dt * GAMMA_C * phi * monod)
        rhs_c[outer_nodes] = 1.0
        c = solve_c(rhs_c)
        np.clip(c, 0.0, 1.0, out=c)
        alpha = alpha + dt * (k_alpha * phi * monod)
    return alpha, phi


if __name__ == "__main__":
    inp = sys.argv[1] if len(sys.argv) > 1 else "p23_biofilm_commensal_static.inp"
    nodes, tets = parse_mesh(inp)
    print("mesh: %d nodes, %d tets" % (len(nodes), len(tets)))
    t = time.time()
    Ml, K = assemble(nodes, tets)
    print("FE assembled in %.1fs" % (time.time() - t))
    outer = np.arange(4 * NV, 5 * NV)        # layer 4 = oral/outer
    np.save("_pde3d_Ml.npy", Ml)
    sp.save_npz("_pde3d_K.npz", K)
    np.save("_pde3d_outer.npy", outer)
    print("saved FE operators (_pde3d_*).")
