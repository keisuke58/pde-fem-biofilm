"""
phase_field_at2_2d_sparse.py
============================
2-D AT2 phase-field fracture using scipy.sparse (structured FD grid).

Variational energy functional (AT2, Miehe et al. 2010):
  E(u,d) = ∫ [g(d)ψ_e + G_c/(2ℓ)(d² + ℓ²|∇d|²)] dΩ
  g(d)   = (1-d)²  (quadratic degradation; +k_res regularisation)
  ψ_e    = ½ C:ε^m:ε^m  (plane-strain undegraded energy from actual strain field)

Euler–Lagrange PDE for d (with irreversibility via history field H):
  -G_c·ℓ·Δd + (G_c/ℓ + 2H)·d = 2H     Neumann BCs: ∇d·n=0
  H(x,z,t) = max_{s≤t} ψ_e(ε^m(u(s)))  (undegraded energy from mechanical solve)

Critical eigenstrain (analytic, prescribed-eigenstrain limit):
  α_c = sqrt(G_c(1-ν) / (2Eℓ))    from ψ_e(α_c) = G_c/(2ℓ)

Staggered scheme (per time step):
  1. Solve 2-D plane-strain elasticity with degraded stiffness g(d^{k-1})
       div(g · C : (ε - α I)) = 0   with u=0 at substrate, traction-free elsewhere
  2. Compute ψ_e from actual strain field (undegraded C, actual ε^m)
  3. Update H = max(H, ψ_e)
  4. Solve AT2 PDE for d^k; enforce irreversibility d^k >= d^{k-1}

Discretisation:
  - Mechanics: bilinear Q4 FEM on (Nx-1)x(Nz-1) cells (1-pt Gauss quadrature,
    avoids volumetric locking for quasi-incompressible nu=0.45)
  - Phase-field: 5-pt FD on Nx×Nz grid (Neumann ghost-node)
  - h <= ell/2 required for accurate phase-field diffuse crack width
"""
from __future__ import annotations

import argparse
import math
import time

import numpy as np
from scipy.sparse import coo_matrix, diags, eye, kron, lil_matrix
from scipy.sparse.linalg import spsolve

# ────────────────────────────────────────────────────────────────────────────
# Physical parameters (SI)
# ────────────────────────────────────────────────────────────────────────────
E_BIO  = 1.0e3      # Young's modulus [Pa]
NU     = 0.45        # Poisson ratio (quasi-incompressible biofilm)
G_C    = 1.0e-5     # Fracture energy [J/m2]
ELL    = 5.0e-6     # Phase-field length scale [m]
W      = 100.0e-6   # Domain width [m]
L      = 60.0e-6    # Film thickness [m]
K_RES  = 1.0e-8     # Residual stiffness (prevents singular K at d=1)

K_RATIO  = 6.44 ** (1.0 / 2.68)   # k_eff^b(CH)/k_eff^b(DH) ~ 2.00
# alpha_c from psi_e(alpha_c)=G_c/(2*ell), psi_e=E/(1-nu)*alpha^2
ALPHA_C  = math.sqrt(G_C * (1.0 - NU) / (2.0 * E_BIO * ELL))  # ~0.0235
N_STEPS  = 120
K_EFF_CH = ALPHA_C / N_STEPS * 1.6
K_EFF_DH = K_EFF_CH / K_RATIO


# ────────────────────────────────────────────────────────────────────────────
# Lame constants
# ────────────────────────────────────────────────────────────────────────────
def lame(E: float, nu: float) -> tuple[float, float]:
    lam = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu  = E / (2.0 * (1.0 + nu))
    return lam, mu


# ────────────────────────────────────────────────────────────────────────────
# Phase-field: 5-pt FD Laplacian (Neumann BCs via ghost nodes)
# ────────────────────────────────────────────────────────────────────────────
def _lap1d_neumann(n: int, h: float):
    mat = lil_matrix((n, n))
    for i in range(n):
        mat[i, i] = -2.0
        if i > 0:      mat[i, i - 1] = 1.0
        if i < n - 1:  mat[i, i + 1] = 1.0
    mat[0, 1]     += 1.0
    mat[n-1, n-2] += 1.0
    return mat.tocsr() / h**2


def build_laplacian(Nx: int, Nz: int, hx: float, hz: float):
    """2-D Laplacian on Nx x Nz grid (x outer, z inner)."""
    Lx = _lap1d_neumann(Nx, hx)
    Lz = _lap1d_neumann(Nz, hz)
    return (kron(eye(Nx, format="csr"), Lz, format="csr")
            + kron(Lx, eye(Nz, format="csr"), format="csr"))


# ────────────────────────────────────────────────────────────────────────────
# Mechanics: bilinear Q4 FEM assembly (vectorised)
# ────────────────────────────────────────────────────────────────────────────
def build_elasticity_system(
    Nx: int, Nz: int, hx: float, hz: float,
    g_2d: np.ndarray, alpha_2d: np.ndarray,
    lam: float, mu: float,
) -> tuple:
    """
    Assemble 2-D plane-strain elasticity K·u = F via bilinear Q4 elements,
    1-pt Gauss quadrature (reduced integration avoids volumetric locking).

    Node DOF: nid = i*Nz + j;  DOF_x = 2*nid,  DOF_z = 2*nid+1.
    Element nodes CCW: n0=(i,j) n1=(i+1,j) n2=(i+1,j+1) n3=(i,j+1).

    Plane-strain with isotropic eigenstrain alpha (eps_yy_total=0):
      eps_yy^m = -alpha  (constraint; contributes to pre-stress)
    Pre-stress: sigma^pre_Voigt = -(3*lam+2*mu)*alpha * [1,1,0]^T
    C_2D (Voigt [eps_xx, eps_zz, gamma_xz]):
      [[lam+2mu, lam, 0], [lam, lam+2mu, 0], [0, 0, mu]]
    """
    dNdx = np.array([-1.0, 1.0, 1.0, -1.0]) / (2.0 * hx)
    dNdz = np.array([-1.0, -1.0, 1.0, 1.0]) / (2.0 * hz)

    # B0 (3x8): Voigt [eps_xx, eps_zz, gamma_xz]
    # cols: [ux_n0, uz_n0, ux_n1, uz_n1, ux_n2, uz_n2, ux_n3, uz_n3]
    B0 = np.zeros((3, 8))
    for a in range(4):
        B0[0, 2*a]   = dNdx[a]          # eps_xx
        B0[1, 2*a+1] = dNdz[a]          # eps_zz
        B0[2, 2*a]   = dNdz[a]          # gamma_xz (u_x part)
        B0[2, 2*a+1] = dNdx[a]          # gamma_xz (u_z part)

    C2d = np.array([[lam + 2*mu, lam, 0.0],
                    [lam, lam + 2*mu, 0.0],
                    [0.0, 0.0,        mu ]])
    k0 = B0.T @ C2d @ B0               # (8x8) reference stiffness
    f0 = B0.T @ np.array([1.0, 1.0, 0.0])  # (8,) reference pre-stress shape

    # Element corner averages (vectorised over all (Nx-1)x(Nz-1) elements)
    Ni, Nj = Nx - 1, Nz - 1
    I_e = np.repeat(np.arange(Ni), Nj)
    J_e = np.tile(np.arange(Nj), Ni)

    g_e = 0.25 * (g_2d[I_e, J_e] + g_2d[I_e+1, J_e]
                  + g_2d[I_e+1, J_e+1] + g_2d[I_e, J_e+1])
    a_e = 0.25 * (alpha_2d[I_e, J_e] + alpha_2d[I_e+1, J_e]
                  + alpha_2d[I_e+1, J_e+1] + alpha_2d[I_e, J_e+1])

    nid = np.array([I_e*Nz+J_e, (I_e+1)*Nz+J_e,
                    (I_e+1)*Nz+(J_e+1), I_e*Nz+(J_e+1)])  # (4, Ne)
    dofs = np.empty((8, len(I_e)), dtype=np.int64)
    for a in range(4):
        dofs[2*a]   = 2 * nid[a]
        dofs[2*a+1] = 2 * nid[a] + 1

    area     = hx * hz
    ke_coeff = area * g_e                                    # (Ne,)
    # ke_vals[a, b, e] = ke_coeff[e] * k0[a, b]
    ke_vals  = k0[:, :, np.newaxis] * ke_coeff[np.newaxis, np.newaxis, :]  # (8,8,Ne)

    row_idx = np.broadcast_to(dofs[:, np.newaxis, :], (8, 8, len(I_e)))
    col_idx = np.broadcast_to(dofs[np.newaxis, :, :], (8, 8, len(I_e)))

    N_dof = 2 * Nx * Nz
    K = coo_matrix(
        (ke_vals.ravel(), (row_idx.ravel(), col_idx.ravel())),
        shape=(N_dof, N_dof),
    ).tocsr()

    # RHS from eigenstrain pre-stress: f_e = -(3lam+2mu)*area*g_e*a_e * f0
    fe_coeff = (3.0*lam + 2.0*mu) * area * g_e * a_e       # (Ne,) RHS = B^T * C * eps^0 > 0
    F = np.zeros(N_dof)
    for a in range(8):
        np.add.at(F, dofs[a], fe_coeff * f0[a])

    return K, F


def apply_dirichlet(K, F, fixed_dofs: np.ndarray):
    """
    Enforce u=0 at fixed_dofs (CSR in-place, no format conversion).
    Zero each fixed row, set diagonal=1, F[dof]=0.
    Column zeroing omitted because prescribed value = 0 (no RHS correction needed).
    """
    K2 = K.copy()
    F2 = F.copy()
    for dof in fixed_dofs:
        d = int(dof)
        rs, re = K2.indptr[d], K2.indptr[d + 1]
        K2.data[rs:re] = 0.0
        cols = K2.indices[rs:re]
        idx  = int(np.searchsorted(cols, d))
        if idx < (re - rs) and cols[idx] == d:
            K2.data[rs + idx] = 1.0
        F2[d] = 0.0
    return K2, F2


def compute_psi_e(
    ux_2d: np.ndarray, uz_2d: np.ndarray, alpha_2d: np.ndarray,
    hx: float, hz: float, lam: float, mu: float,
) -> np.ndarray:
    """
    Plane-strain undegraded elastic energy density: psi_e = 0.5 * C:eps^m:eps^m.

    Mechanical strains (total - eigenstrain):
      eps_xx^m = du_x/dx - alpha
      eps_zz^m = du_z/dz - alpha
      eps_yy^m = -alpha  (plane-strain: eps_yy_total = 0)
      gamma_xz = du_x/dz + du_z/dx

    Undegraded stresses (full C, no g(d)):
      sigma_xx^0 = lam*e_kk + 2mu*eps_xx^m
      sigma_yy^0 = lam*e_kk
      sigma_zz^0 = lam*e_kk + 2mu*eps_zz^m
      sigma_xz^0 = mu*gamma_xz

    psi_e = 0.5*(sigma_xx^0*eps_xx^m + sigma_yy^0*eps_yy^m
                 + sigma_zz^0*eps_zz^m + sigma_xz^0*gamma_xz)

    H uses this undegraded energy — standard AT2 staggered (Miehe 2010).
    """
    exx_m = np.gradient(ux_2d, hx, axis=0) - alpha_2d
    ezz_m = np.gradient(uz_2d, hz, axis=1) - alpha_2d
    eyy_m = -alpha_2d                                       # plane-strain
    gxz   = (np.gradient(ux_2d, hz, axis=1)
              + np.gradient(uz_2d, hx, axis=0))

    ekk  = exx_m + eyy_m + ezz_m
    sxx0 = lam * ekk + 2.0 * mu * exx_m
    syy0 = lam * ekk
    szz0 = lam * ekk + 2.0 * mu * ezz_m
    sxz0 = mu * gxz

    return 0.5 * (sxx0*exx_m + syy0*eyy_m + szz0*ezz_m + sxz0*gxz)


# ────────────────────────────────────────────────────────────────────────────
# Growth profile
# ────────────────────────────────────────────────────────────────────────────
def make_k_field(k_base: float, profile: str,
                 X: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """Spatially varying k_eff^b at each node. X, Z: (Nx, Nz). Returns (Nx*Nz,)."""
    z_norm    = Z / L
    z_profile = z_norm ** 0.8 + 0.01   # surface-biased; +0.01 avoids zero

    if profile == "z_gradient":
        phi = z_profile
    elif profile == "z_patchy":
        patch = 1.0 + 0.6 * np.sin(2.0 * np.pi * X / W * 2.5)
        phi = z_profile * np.clip(patch, 0.1, None)
    elif profile == "pg_substrate":
        phi = np.exp(-Z / (L / 4.0))
    elif profile == "pg_combo":
        phi = z_profile / K_RATIO + 0.8 * np.exp(-Z / (L / 5.0))
    elif profile == "uniform":
        phi = np.ones_like(Z)
    else:
        raise ValueError(profile)

    phi /= phi.mean()
    return (k_base * phi).ravel()


# ────────────────────────────────────────────────────────────────────────────
# Fully-coupled staggered AT2 solver
# ────────────────────────────────────────────────────────────────────────────
def run_coupled(
    k_base: float, profile: str,
    Nx: int, Nz: int, n_steps: int = N_STEPS,
    label: str = "",
) -> dict:
    """
    Full AT2 staggered with 2-D plane-strain mechanical solve at each step.

    Step k:
      g = (1-d^{k-1})^2 + K_RES
      K(g)*u = F(g,alpha)  ->  u_x, u_z  [Q4 FEM, 2*Nx*Nz DOFs]
      psi_e = 0.5*C:eps^m:eps^m  (undegraded, actual strain field)
      H^k = max(H^{k-1}, psi_e)
      AT2 PDE -> d^k;  irreversibility d^k = max(d^k, d^{k-1})
    """
    lam, mu = lame(E_BIO, NU)
    hx = W / (Nx - 1);  hz = L / (Nz - 1)
    X, Z = np.meshgrid(np.linspace(0, W, Nx), np.linspace(0, L, Nz), indexing="ij")
    k_field = make_k_field(k_base, profile, X, Z).reshape(Nx, Nz)

    # Substrate BC: u_x=u_z=0 at j=0 (z=0)
    j0_nids    = np.arange(Nx) * Nz
    fixed_dofs = np.concatenate([2*j0_nids, 2*j0_nids+1])

    LAP   = build_laplacian(Nx, Nz, hx, hz)
    A_pf0 = -G_C * ELL * LAP

    d_flat = np.zeros(Nx * Nz)
    H_flat = np.zeros(Nx * Nz)
    results = []
    sub_mask = (Z < hz * 2).ravel()

    for step in range(1, n_steps + 1):
        t        = step * (ALPHA_C * 1.6 / n_steps) / k_base
        alpha_2d = k_field * t

        d_2d = d_flat.reshape(Nx, Nz)
        g_2d = (1.0 - d_2d)**2 + K_RES

        # ── Mechanical solve ─────────────────────────────────────────────
        K_el, F_el = build_elasticity_system(
            Nx, Nz, hx, hz, g_2d, alpha_2d, lam, mu)
        K_el, F_el = apply_dirichlet(K_el, F_el, fixed_dofs)
        u_vec = spsolve(K_el, F_el)

        ux_2d = u_vec[0::2].reshape(Nx, Nz)
        uz_2d = u_vec[1::2].reshape(Nx, Nz)

        # ── Undegraded elastic energy from actual strain ──────────────────
        psi_e = np.maximum(
            compute_psi_e(ux_2d, uz_2d, alpha_2d, hx, hz, lam, mu),
            0.0)  # numerical noise guard

        # ── Phase-field update ────────────────────────────────────────────
        H_flat = np.maximum(H_flat, psi_e.ravel())
        diag_vals = G_C / ELL + 2.0 * H_flat
        A_pf = A_pf0 + diags(diag_vals, format="csr")
        d_new = spsolve(A_pf, 2.0 * H_flat)
        d_flat = np.clip(np.maximum(d_new, d_flat), 0.0, 1.0)

        d_max  = d_flat.max()
        d_base = d_flat[sub_mask].max() if sub_mask.any() else 0.0

        results.append({
            "step": step, "t": t,
            "alpha_mean": float(alpha_2d.mean()),
            "d_max":  d_max, "d_base": d_base,
            "d_2d":   d_flat.reshape(Nx, Nz).copy(),
            "psi_e":  psi_e.copy(),
        })

        if label and step % 20 == 0:
            print(f"  [{label}] step {step:3d}  "
                  f"a_mean={alpha_2d.mean():.4f}  "
                  f"d_max={d_max:.3f}  d_base={d_base:.3f}")

        if d_base > 0.95:
            if label:
                print(f"  [{label}] substrate saturated at step {step}")
            break

    return {"label": label, "results": results,
            "X": X, "Z": Z, "d_final": d_flat.reshape(Nx, Nz)}


# ────────────────────────────────────────────────────────────────────────────
# Prescribed-eigenstrain solver (fast reference, no mechanical solve)
# ────────────────────────────────────────────────────────────────────────────
def run(
    k_base: float, profile: str,
    Nx: int, Nz: int, n_steps: int = N_STEPS,
    label: str = "",
) -> dict:
    """
    AT2 with prescribed-eigenstrain approx: no mechanical solve.
    psi_e = E/(1-nu)*alpha^2 (biaxial plane-stress closed form).
    Fast reference / comparison for the coupled solver.
    """
    hx = W / (Nx - 1);  hz = L / (Nz - 1)
    X, Z = np.meshgrid(np.linspace(0, W, Nx), np.linspace(0, L, Nz), indexing="ij")
    k_field = make_k_field(k_base, profile, X, Z)

    LAP   = build_laplacian(Nx, Nz, hx, hz)
    A_pf0 = -G_C * ELL * LAP

    d_flat = np.zeros(Nx * Nz)
    H_flat = np.zeros(Nx * Nz)
    results = []
    sub_mask = (Z < hz * 2).ravel()

    for step in range(1, n_steps + 1):
        t     = step * (ALPHA_C * 1.6 / n_steps) / k_base
        alpha = k_field * t
        psi_e = E_BIO / (1.0 - NU) * alpha**2   # biaxial plane-stress limit

        H_flat = np.maximum(H_flat, psi_e)
        diag_vals = G_C / ELL + 2.0 * H_flat
        A_pf = A_pf0 + diags(diag_vals, format="csr")
        d_new = spsolve(A_pf, 2.0 * H_flat)
        d_flat = np.clip(np.maximum(d_new, d_flat), 0.0, 1.0)

        d_max  = d_flat.max()
        d_base = d_flat[sub_mask].max() if sub_mask.any() else 0.0

        results.append({
            "step": step, "t": t,
            "alpha_mean": float(alpha.mean()),
            "d_max": d_max, "d_base": d_base,
            "d_2d":  d_flat.reshape(Nx, Nz).copy(),
        })

        if label and step % 30 == 0:
            print(f"  [{label}] step {step:3d}  d_max={d_max:.3f}")

        if d_base > 0.95:
            break

    return {"label": label, "results": results,
            "X": X, "Z": Z, "d_final": d_flat.reshape(Nx, Nz)}


# ────────────────────────────────────────────────────────────────────────────
# Helper
# ────────────────────────────────────────────────────────────────────────────
def t_crit(results: list, threshold: float = 0.5) -> float | None:
    for r in results:
        if r["d_max"] > threshold:
            return r["t"]
    return None


# ────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot",    action="store_true")
    parser.add_argument("--coupled", action="store_true",
                        help="full 2-D mechanical coupling via Q4 FEM each step")
    parser.add_argument("--nx",      type=int,   default=120)
    parser.add_argument("--nz",      type=int,   default=80)
    parser.add_argument("--n-steps", type=int,   default=N_STEPS)
    parser.add_argument("--save",    default="phase_field_at2_2d_coupled.pdf")
    parser.add_argument("--width",   type=float, default=None,
                        help="domain width [m] (default 100e-6)")
    parser.add_argument("--length",  type=float, default=None,
                        help="film thickness [m] (default 60e-6)")
    args = parser.parse_args()

    global W, L
    if args.width:  W = args.width
    if args.length: L = args.length

    Nx, Nz = args.nx, args.nz
    hx = W / (Nx - 1);  hz = L / (Nz - 1)
    solver = run_coupled if args.coupled else run
    mode   = "fully-coupled Q4 FEM" if args.coupled else "prescribed-eigenstrain"
    lam_v, mu_v = lame(E_BIO, NU)

    print("=" * 70)
    print(f"AT2 Phase-Field 2-D  [{mode}]")
    print("=" * 70)
    print(f"  E={E_BIO:.0f} Pa  nu={NU}  G_c={G_C:.1e} J/m2  ell={ELL*1e6:.0f} um")
    print(f"  lam={lam_v:.1f} Pa  mu={mu_v:.1f} Pa  (3lam+2mu)={3*lam_v+2*mu_v:.1f} Pa")
    print(f"  Domain {W*1e6:.0f}x{L*1e6:.0f} um  grid {Nx}x{Nz}={Nx*Nz} nodes")
    print(f"  h_x={hx*1e6:.2f} um  h_z={hz*1e6:.2f} um  "
          f"h/ell_x={hx/ELL:.2f}  h/ell_z={hz/ELL:.2f}", end="")
    if max(hx, hz) > ELL / 2:
        print("  <- WARNING h > ell/2", end="")
    print()
    print(f"  alpha_c={ALPHA_C:.4f}  K_RATIO={K_RATIO:.3f}")
    if args.coupled:
        print(f"  Mechanical DOFs: {2*Nx*Nz}   Phase-field DOFs: {Nx*Nz}")
    print()

    cases = [
        (K_EFF_CH, "z_gradient",   "CH"),
        (K_EFF_DH, "z_gradient",   "DH-base"),
        (K_EFF_DH, "pg_substrate", "DH-Pg"),
    ]
    all_res = {}
    for k, prof, lbl in cases:
        print(f"Running {lbl:8s} ...")
        t0 = time.time()
        all_res[lbl] = solver(k, prof, Nx, Nz, args.n_steps, label=lbl)
        print(f"  done in {time.time()-t0:.1f} s\n")

    print("-" * 70)
    tc = {lbl: t_crit(out["results"]) for lbl, out in all_res.items()}
    for lbl, out in all_res.items():
        dm   = out["results"][-1]["d_max"]
        tc_s = f"{tc[lbl]:.1f} s" if tc[lbl] else "not reached"
        print(f"  {lbl:10s}  t_crit={tc_s:>12s}  d_max_final={dm:.3f}")
    print()
    if tc["CH"] and tc["DH-base"]:
        print(f"  Effect 1 (baseline): t_crit(DH-base)/t_crit(CH) = "
              f"{tc['DH-base']/tc['CH']:.3f}   theory ~ {K_RATIO:.3f}")
    if tc["CH"] and tc["DH-Pg"]:
        print(f"  Effect 2 (DH+Pg):    t_crit(DH-Pg)/t_crit(CH)  = "
              f"{tc['DH-Pg']/tc['CH']:.3f}")

    if not args.plot:
        return

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize

    ratio_str = (f"{tc['DH-base']/tc['CH']:.2f}"
                 if (tc["CH"] and tc["DH-base"]) else "N/A")
    pg_str = (f"{tc['DH-Pg']/tc['CH']:.2f}"
              if (tc["CH"] and tc["DH-Pg"]) else "N/A")

    x_um = np.linspace(0, W * 1e6, Nx)
    z_um = np.linspace(0, L * 1e6, Nz)

    fig = plt.figure(figsize=(15, 11))
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.30,
                            height_ratios=[1, 1, 0.85])
    norm_d    = Normalize(vmin=0, vmax=1)
    fracs     = [0.30, 0.65, 1.0]
    stages    = ["Early", "Mid", "Final"]
    map_cases = [("CH", 0), ("DH-Pg", 1)]

    for lbl, row in map_cases:
        res = all_res[lbl]["results"]
        for col, (frac, stage) in enumerate(zip(fracs, stages)):
            idx   = max(0, min(len(res)-1, int(len(res)*frac)))
            d_2d  = res[idx]["d_2d"]
            t_val = res[idx]["t"]
            ax    = fig.add_subplot(gs[row, col])
            im = ax.pcolormesh(x_um, z_um, d_2d.T, norm=norm_d,
                               cmap="hot_r", shading="gouraud", rasterized=True)
            try:
                ax.contour(x_um, z_um, d_2d.T, levels=[0.5],
                           colors="cyan", linewidths=0.7, linestyles="--")
            except Exception:
                pass
            ax.set_xlim(0, W*1e6);  ax.set_ylim(0, L*1e6)
            ax.set_xlabel("x [um]", fontsize=7)
            ax.set_ylabel("z [um]", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_title(f"{lbl} - {stage}  t={t_val:.0f}s  "
                         f"d_max={d_2d.max():.2f}", fontsize=8)
            cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
            cb.set_label("d", fontsize=7);  cb.ax.tick_params(labelsize=6)

    ax_l = fig.add_subplot(gs[2, :2])
    ax_r = fig.add_subplot(gs[2,  2])
    colors  = {"CH": "#1f77b4", "DH-base": "#2ca02c", "DH-Pg": "#d62728"}
    markers = {"CH": "o",       "DH-base": "s",        "DH-Pg": "^"}

    for lbl, out in all_res.items():
        res   = out["results"]
        ts    = [r["t"]      for r in res]
        dmaxs = [r["d_max"]  for r in res]
        dbas  = [r["d_base"] for r in res]
        c, m  = colors[lbl], markers[lbl]
        ax_l.plot(ts, dmaxs, color=c, marker=m, ms=2, lw=1.5, label=f"{lbl} d_max")
        ax_l.plot(ts, dbas,  color=c, marker=m, ms=2, lw=0.8,
                  ls="--", label=f"{lbl} d_substrate")

    ax_l.axhline(0.5, color="gray", ls=":", lw=0.8, label="d=0.5 (t_crit)")
    ax_l.set_xlabel("Time [s]", fontsize=8);  ax_l.set_ylabel("Damage d", fontsize=8)
    ax_l.set_title("d_max(t) and d_substrate(t)", fontsize=8)
    ax_l.legend(fontsize=6, ncol=2, loc="upper left");  ax_l.set_ylim(0, 1.08)

    lbls_b = list(all_res.keys())
    vals_b = [all_res[lb]["results"][-1]["d_max"] for lb in lbls_b]
    bars = ax_r.bar(lbls_b, vals_b,
                    color=[colors[lb] for lb in lbls_b],
                    alpha=0.8, edgecolor="k", lw=0.5)
    ax_r.set_ylim(0, 1.15);  ax_r.set_ylabel("d_max (final)", fontsize=8)
    ax_r.set_title("Final damage", fontsize=8);  ax_r.tick_params(labelsize=7)
    for bar, v in zip(bars, vals_b):
        ax_r.text(bar.get_x() + bar.get_width()/2, v + 0.03,
                  f"{v:.2f}", ha="center", va="bottom", fontsize=8)

    fig.suptitle(
        f"AT2 Phase-Field 2-D  [{mode}]  -  Biofilm growth fracture\n"
        f"CH (z-gradient, fast)  |  DH-base (2x slower)  |  DH-Pg (Pg@substrate)\n"
        f"E={E_BIO:.0f} Pa  nu={NU}  G_c={G_C:.0e} J/m2  ell={ELL*1e6:.0f} um  "
        f"grid {Nx}x{Nz}  h/ell={max(hx,hz)/ELL:.2f}\n"
        f"K_RATIO={K_RATIO:.2f}  t_crit(DH-base/CH)={ratio_str}  "
        f"t_crit(DH-Pg/CH)={pg_str}",
        fontsize=8,
    )
    plt.savefig(args.save, bbox_inches="tight", dpi=180)
    print(f"\nSaved: {args.save}")


if __name__ == "__main__":
    main()
