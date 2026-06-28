"""
phase_field_at2_2d.py
=====================
2-D AT2 phase-field fracture for biofilm growth on a rigid substrate.
Uses FEniCSx (dolfinx) for the variational formulation.

Geometry
--------
  Rectangular domain  Ω = [0, W] × [0, L]
    z=0 : substrate   (bottom, rigid attachment)
    z=L : free surface (top)
    x=0,W : periodic  (or symmetry)

  W = 100 µm (lateral width)
  L =  60 µm (thickness)

Physics
-------
  Biofilm grows laterally but substrate constrains expansion.
  In-plane eigenstress (biaxial):  σ₀ = −E·α(x,z,t)

  AT2 energy functional:
    𝒻(d) = ∫_Ω [ (1−d)² Ψ_e  +  G_c/(4ℓ) (d²/ℓ + ℓ |∇d|²) ] dΩ

  Staggered (alternating-minimisation) solve:
    Step 1 — Mechanics  : find u minimising 𝒻 for fixed d
    Step 2 — History    : H ← max(H, Ψ_e(u))
    Step 3 — Damage     : find d minimising 𝒻 for fixed H (irreversibility)

  In the biaxial-growth 1-field reduction (no displacement solve needed):
    Ψ_e(x,z,t) = ½ E α(x,z,t)²    (stored compression energy)
    Damage PDE:  −G_c ℓ Δd + G_c/ℓ d = 2(1−d)H

  This is exact for uniform lateral growth under full constraint.
  When the depth profile φ(x,z) varies spatially, H(x,z) drives localisation.

Growth model
------------
  α(x,z,t) = k_eff^b(x,z) · t
  k_eff^b(x,z) = Σ_i φ_i(x,z) K_i^b

  Three conditions (from thesis MAP):
    CH  : uniform Pg  →  uniform k_eff^b = K_CH
    DH-u: uniform Pg  →  uniform k_eff^b = K_DH < K_CH
    DH-c: Pg concentrated near z=0  →  exponential depth profile

  K_CH / K_DH = 6.44^(1/2.68) ≈ 2.00

Parameters (SI)
---------------
  E   = 1 kPa,  G_c = 10 µJ/m²,  ℓ = 5 µm,  W = 100 µm,  L = 60 µm

Usage
-----
  conda activate fenics
  python phase_field_at2_2d.py              # print results
  python phase_field_at2_2d.py --plot       # save PDF figures
  python phase_field_at2_2d.py --n 60      # coarser mesh (faster)
"""
from __future__ import annotations

import argparse
import math
import sys

import numpy as np

# --------------------------------------------------------------------------
# Physical parameters (SI)
# --------------------------------------------------------------------------
E_BIO = 1.0e3      # Young's modulus  [Pa]
G_C   = 1.0e-5     # Fracture energy  [J/m²]
ELL   = 5.0e-6     # Phase-field length scale  [m]
W     = 100.0e-6   # Domain width     [m]
L     = 60.0e-6    # Film thickness   [m]

K_RATIO   = 6.44 ** (1.0 / 2.68)   # k_eff^b(CH) / k_eff^b(DH) ≈ 2.00
ALPHA_C   = math.sqrt(G_C / (E_BIO * ELL))   # critical eigenstrain ≈ 0.045
N_STEPS   = 100                      # loading increments
K_EFF_CH  = ALPHA_C / N_STEPS * 1.6  # reach 1.6×α_c in N_STEPS steps
K_EFF_DH  = K_EFF_CH / K_RATIO


# --------------------------------------------------------------------------
# Depth profile factory
# --------------------------------------------------------------------------
def make_k_field(k_base: float, profile: str, coords: np.ndarray) -> np.ndarray:
    """
    Return k_eff^b at each DOF location.

    coords : (N, 2) array of (x, z) coordinates
    profile : 'uniform' | 'pg_substrate'
    """
    z = coords[:, 1]
    if profile == "uniform":
        phi_rel = np.ones_like(z)
    elif profile == "pg_substrate":
        # Pg concentrated near z=0 (substrate), characteristic depth L/4
        phi_rel = np.exp(-z / (L / 4.0))
    else:
        raise ValueError(f"Unknown profile: {profile}")
    # Normalise to mean=1 so k_base = mean(k_eff^b)
    phi_rel = phi_rel / phi_rel.mean()
    return k_base * phi_rel


# --------------------------------------------------------------------------
# FEniCSx implementation
# --------------------------------------------------------------------------
def run_fenics(k_base: float, profile: str, n_mesh: int,
               n_steps: int = N_STEPS, label: str = "") -> dict:
    """
    Staggered AT2 solve via FEniCSx.
    Returns dict with time series of d_max, d_base, t.
    """
    from mpi4py import MPI
    import dolfinx
    from dolfinx import mesh as dmesh, fem, io
    from dolfinx.fem import functionspace, Function, Constant, form
    from dolfinx.fem.petsc import LinearProblem
    import ufl
    from petsc4py import PETSc

    comm = MPI.COMM_WORLD

    # -- Mesh: rectangle [0,W]×[0,L]
    msh = dmesh.create_rectangle(
        comm,
        [np.array([0.0, 0.0]), np.array([W, L])],
        [n_mesh, int(n_mesh * L / W)],
        cell_type=dmesh.CellType.triangle,
    )

    # -- Function space (CG1 for scalar damage)
    V = functionspace(msh, ("Lagrange", 1))

    # -- Functions
    d     = Function(V, name="damage")
    d_old = Function(V, name="damage_old")   # for irreversibility
    H_fn  = Function(V, name="history")       # history field

    # -- Coordinate array for k_eff^b field
    coords = V.tabulate_dof_coordinates()[:, :2]
    k_field_np = make_k_field(k_base, profile, coords)
    k_fn = Function(V, name="k_eff_b")
    k_fn.x.array[:] = k_field_np

    # -- Trial/test functions for damage PDE
    v = ufl.TestFunction(V)
    d_trial = ufl.TrialFunction(V)

    # AT2 damage PDE (linearised): for fixed H
    # -G_c·ℓ·Δd + (G_c/ℓ + 2H)·d = 2H
    # Weak form: G_c·ℓ·∫∇d·∇v + (G_c/ℓ + 2H)·∫d·v = 2·∫H·v
    gc_ell = Constant(msh, PETSc.ScalarType(G_C * ELL))
    gc_over_ell = Constant(msh, PETSc.ScalarType(G_C / ELL))

    dx = ufl.Measure("dx", domain=msh)

    def assemble_damage_system(H: Function):
        a = (gc_ell * ufl.dot(ufl.grad(d_trial), ufl.grad(v))
             + (gc_over_ell + 2.0 * H) * d_trial * v) * dx
        L_ = 2.0 * H * v * dx
        return a, L_

    # -- Time loop
    results = []
    H_fn.x.array[:] = 0.0
    d.x.array[:] = 0.0

    for step in range(1, n_steps + 1):
        t = step * (ALPHA_C * 1.6 / n_steps) / k_base

        # Eigenstrain and elastic energy density at DOFs
        alpha_np = k_field_np * t
        psi_e_np = 0.5 * E_BIO * alpha_np**2

        # History update (irreversibility)
        H_fn.x.array[:] = np.maximum(H_fn.x.array, psi_e_np)

        # Damage solve (Neumann BCs: no essential BCs needed)
        d_old.x.array[:] = d.x.array.copy()
        a_d, L_d = assemble_damage_system(H_fn)
        prob = LinearProblem(a_d, L_d,
                             petsc_options_prefix="damage_",
                             petsc_options={"ksp_type": "cg",
                                            "pc_type": "hypre"})
        d_new = prob.solve()

        # Enforce irreversibility: d >= d_old
        d.x.array[:] = np.maximum(d_new.x.array, d_old.x.array)
        d.x.array[:] = np.clip(d.x.array, 0.0, 1.0)

        d_max  = float(d.x.array.max())
        # d at substrate: nodes near z=0
        z_coords = coords[:, 1]
        substrate_mask = z_coords < (L / n_mesh * 2)
        d_base = float(d.x.array[substrate_mask].max()) if substrate_mask.any() else 0.0

        results.append({
            "step": step, "t": t,
            "alpha_mean": float(alpha_np.mean()),
            "d_max": d_max,
            "d_base": d_base,
            "d_array": d.x.array.copy(),
            "coords": coords.copy(),
        })

        if label and step % 20 == 0:
            print(f"  [{label}] step {step:3d}/{n_steps}  "
                  f"ᾱ={alpha_np.mean():.3f}  d_max={d_max:.3f}  d_base={d_base:.3f}")

        if d_base > 0.95:
            if label:
                print(f"  [{label}] substrate fully damaged at step {step}")
            break

    return {"label": label, "results": results,
            "coords": coords, "d_final": d.x.array.copy()}


# --------------------------------------------------------------------------
# t_crit helper
# --------------------------------------------------------------------------
def t_crit(results: list, threshold: float = 0.5) -> float | None:
    for r in results:
        if r["d_max"] > threshold:
            return r["t"]
    return None


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("-n", "--n-mesh", type=int, default=40,
                        help="Mesh resolution in x (default 40)")
    parser.add_argument("--n-steps", type=int, default=N_STEPS)
    parser.add_argument("--save", type=str, default="phase_field_at2_2d.pdf")
    args = parser.parse_args()

    print("=" * 62)
    print("AT2 Phase-Field 2-D: Biofilm Growth Fracture  (FEniCSx)")
    print("=" * 62)
    print(f"  E={E_BIO:.0f} Pa  G_c={G_C:.1e} J/m²  ℓ={ELL*1e6:.0f} µm")
    print(f"  Domain: {W*1e6:.0f}×{L*1e6:.0f} µm  mesh: {args.n_mesh}×{int(args.n_mesh*L/W)}")
    print(f"  α_c={ALPHA_C:.3f}  k_ratio(CH/DH)={K_RATIO:.2f}")
    print()

    cases = [
        (K_EFF_CH, "uniform",       "CH"),
        (K_EFF_DH, "uniform",       "DH-unif"),
        (K_EFF_DH, "pg_substrate",  "DH-Pg"),
    ]
    all_results = {}
    for k, profile, label in cases:
        print(f"Running {label:8s} (k={k:.3e}, profile={profile}) ...")
        out = run_fenics(k, profile, args.n_mesh, args.n_steps, label=label)
        all_results[label] = out
        print()

    print("-" * 62)
    print("Summary:")
    for label, out in all_results.items():
        tc = t_crit(out["results"])
        dm = out["results"][-1]["d_max"]
        print(f"  {label:10s}: t_crit={tc:.1f}s" if tc else f"  {label:10s}: t_crit=not reached", end="")
        print(f"   d_max_final={dm:.3f}")

    res_ch  = all_results["CH"]["results"]
    res_dhu = all_results["DH-unif"]["results"]
    res_dhp = all_results["DH-Pg"]["results"]
    tc_ch  = t_crit(res_ch)
    tc_dhu = t_crit(res_dhu)
    tc_dhp = t_crit(res_dhp)

    print()
    print("  Effect 1 (uniform):  t_crit(DH-unif)/t_crit(CH) =",
          f"{tc_dhu/tc_ch:.2f}" if (tc_ch and tc_dhu) else "N/A",
          f"  (theory {K_RATIO:.2f})")
    print("  Effect 2 (Pg@sub):   t_crit(DH-Pg)/t_crit(CH)   =",
          f"{tc_dhp/tc_ch:.2f}" if (tc_ch and tc_dhp) else "N/A",
          "  (<1 → local substrate damage)")

    if not args.plot:
        return

    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    fig.subplots_adjust(hspace=0.45, wspace=0.35)

    fracs  = [0.35, 0.70, 1.0]
    titles = ["Early", "Mid", "Final"]

    for row, (label, out) in enumerate(all_results.items()):
        coords = out["coords"]
        xs = coords[:, 0] * 1e6
        zs = coords[:, 1] * 1e6
        tri = mtri.Triangulation(xs, zs)

        res = out["results"]
        for col, frac in enumerate(fracs):
            idx = max(0, min(len(res) - 1, int(len(res) * frac)))
            d_vals = res[idx]["d_array"]
            alpha_m = res[idx]["alpha_mean"]

            ax = axes[row, col]
            tc_ = axes[row, col].tripcolor(tri, d_vals, vmin=0, vmax=1,
                                            cmap="hot_r", shading="gouraud")
            ax.set_xlim(0, W * 1e6)
            ax.set_ylim(0, L * 1e6)
            ax.set_xlabel("x [µm]", fontsize=7)
            ax.set_ylabel("z [µm]", fontsize=7)
            ax.set_title(f"{label} — {titles[col]}\nᾱ={alpha_m:.3f}  d_max={d_vals.max():.2f}",
                         fontsize=8)
            plt.colorbar(tc_, ax=ax, label="d", fraction=0.04)

    fig.suptitle(
        "AT2 Phase-Field 2-D: Biofilm Growth Fracture\n"
        f"E={E_BIO:.0f} Pa, G_c={G_C:.1e} J/m², ℓ={ELL*1e6:.0f} µm  |  "
        f"k_ratio(CH/DH)={K_RATIO:.2f}",
        fontsize=11
    )
    plt.savefig(args.save, bbox_inches="tight", dpi=150)
    print(f"\nSaved {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
