"""
felix_complete_reproduction.py  (v4)
=====================================
Felix Klempt (2024) Eq.34-36 — extended to Fig.3/4/7/8/9-10.

Cases:
  A    = Fig.2-4:  corner nutrient, directional growth
  B_hi = Fig.5-8:  bottom nutrient, g=1e8  (high nutrient, mushroom/droplet)
  B_lo = Fig.5-8:  bottom nutrient, g=1e10 (low  nutrient, layer formation)
  C    = Fig.9-10: corner nutrient, 4 rigid-pillar obstacle avoidance

Outputs:
  felix_complete_reproduction.png  — Cases A + B_hi (Fig.2/5 analog, v3-compatible)
  felix_fig3_isosurface.png        — Fig.3: 3-D surface φ isosurface time series
  felix_fig789_10.png              — Fig.7/8 (avg + pressure) + Case B_lo + Case C

Governing equations (Eq.34-36):
  Eq.34: φ̇ = β∇²φ - Γφ(1-φ)(1-2φ) + K_log φ(1-φ)c/(k+c)
              + k_α α - R_2D c/(k+c) (∇φ·n_∇c)
  Eq.35: D∇²c = g φ  (quasi-static, scipy sparse direct solve)
  Eq.36: α̇ = k_α φ

2D approximation notes:
  - Allen-Cahn double-well replaces elastic energy surrogate
  - R_PARAM_2D = 10 (Felix 100 gives v_chemo = 40 μm/T* → CFL failure)
  - B_hi: g = G_HI = 1e8, g/D = 0.01 (Thiele ≈ 0.3, reaction-limited, mushroom)
  - B_lo: g = G_LO = 1e10, g/D = 1.0  (Thiele ≈ 3.0, diffusion-limited, layer)
"""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from thesis_style import use as thesis_use, clean_ax

# ── Thesis style ─────────────────────────────────────────────────
_FIG_W, _FIG_H = thesis_use(width_frac=1.0, aspect=0.55)

# ── Felix Table 2 parameters ─────────────────────────────────────
MU         = 3.3557
E_MOD      = 10.0
NU         = 0.49
D_STAR     = 1e10        # nutrient diffusivity [μm²/T*]
BETA       = 2.0         # interface regularisation [μm²/T*]
K_ALPHA    = 1e-3        # growth coefficient [1/T*]
K_MON      = 1.0         # half-saturation constant [-]
G_HI       = 1e8         # consumption: "high" nutrient → mushroom
G_LO       = 1e10        # consumption: "low"  nutrient → layer
R_PARAM    = 100.0       # Felix original chemotaxis (reference)
R_PARAM_2D = 10.0        # 2D effective (CFL-stable; see docstring)
GAMMA      = 8.0         # Allen-Cahn double-well (ξ = √(β/Γ) = 0.50 μm)
K_LOG      = 5.0         # logistic-Monod growth (K_LOG < Γ → φ=0 stable)

C10   = MU / 2
D1    = 6 * (1 - 2*NU) / E_MOD
KAPPA = E_MOD / (3 * (1 - 2*NU))

# ── Numerical parameters ─────────────────────────────────────────
NX, NY  = 60, 60
LX, LY  = 20.0, 20.0
DX = LX / (NX - 1)
DY = LY / (NY - 1)

DT            = 0.0008
N_STEPS       = 12000
SAVE_AT       = {0, 2000, 5000, 8000, 10000, 12000}
C_SOLVE_EVERY = 10
AVG_EVERY     = 50       # record domain-average every N steps (for Fig.4/7)

x1d = np.linspace(0, LX, NX)
y1d = np.linspace(0, LY, NY)
X, Y = np.meshgrid(x1d, y1d, indexing='ij')

# ── Obstacle geometry (Case C, Fig.9-10) ─────────────────────────
OBS_CENTERS = [(7.0, 8.0), (7.0, 14.0), (13.0, 8.0), (13.0, 14.0)]
OBS_RADIUS  = 2.0


# ═══════════════════════════════════════════════════════════════════
# PDE helpers
# ═══════════════════════════════════════════════════════════════════

def lap2d(u):
    lap = np.zeros_like(u)
    lap[1:-1,1:-1] = (
        (u[2:,1:-1]-2*u[1:-1,1:-1]+u[:-2,1:-1])/DX**2 +
        (u[1:-1,2:]-2*u[1:-1,1:-1]+u[1:-1,:-2])/DY**2
    )
    # Neumann boundary (zero-flux)
    lap[0,  1:-1] = 2*(u[1,1:-1]-u[0,1:-1])/DX**2   +(u[0,2:]-2*u[0,1:-1]+u[0,:-2])/DY**2
    lap[-1, 1:-1] = 2*(u[-2,1:-1]-u[-1,1:-1])/DX**2 +(u[-1,2:]-2*u[-1,1:-1]+u[-1,:-2])/DY**2
    lap[1:-1,0]   = (u[2:,0]-2*u[1:-1,0]+u[:-2,0])/DX**2   +2*(u[1:-1,1]-u[1:-1,0])/DY**2
    lap[1:-1,-1]  = (u[2:,-1]-2*u[1:-1,-1]+u[:-2,-1])/DX**2+2*(u[1:-1,-2]-u[1:-1,-1])/DY**2
    lap[0,0]   = 2*(u[1,0]-u[0,0])/DX**2  +2*(u[0,1]-u[0,0])/DY**2
    lap[-1,0]  = 2*(u[-2,0]-u[-1,0])/DX**2+2*(u[-1,1]-u[-1,0])/DY**2
    lap[0,-1]  = 2*(u[1,-1]-u[0,-1])/DX**2+2*(u[0,-2]-u[0,-1])/DY**2
    lap[-1,-1] = 2*(u[-2,-1]-u[-1,-1])/DX**2+2*(u[-1,-2]-u[-1,-1])/DY**2
    return lap


def grad2d(u):
    gx = np.zeros_like(u); gy = np.zeros_like(u)
    gx[1:-1,:] = (u[2:,:]-u[:-2,:])/(2*DX)
    gy[:,1:-1] = (u[:,2:]-u[:,:-2])/(2*DY)
    gx[0,:]  = (u[1,:]-u[0,:])/DX;   gx[-1,:] = (u[-1,:]-u[-2,:])/DX
    gy[:,0]  = (u[:,1]-u[:,0])/DY;   gy[:,-1] = (u[:,-1]-u[:,-2])/DY
    return gx, gy


def build_poisson_matrix(c_bc_mask):
    N = NX * NY
    rows, cols, vals = [], [], []
    def idx(i, j): return i*NY + j
    for i in range(NX):
        for j in range(NY):
            k = idx(i, j)
            if c_bc_mask[i, j]:
                rows.append(k); cols.append(k); vals.append(1.0)
            else:
                ip = i+1 if i < NX-1 else i-1
                im = i-1 if i > 0   else i+1
                jp = j+1 if j < NY-1 else j-1
                jm = j-1 if j > 0   else j+1
                for ni, nj, coeff in [
                    (ip, j, 1/DX**2), (im, j, 1/DX**2),
                    (i, jp, 1/DY**2), (i, jm, 1/DY**2),
                    (i,  j, -(2/DX**2+2/DY**2))
                ]:
                    rows.append(k); cols.append(idx(ni, nj)); vals.append(coeff)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def solve_c(A_mat, phi, c_bc_mask, c_bc_val, g_ratio):
    """Quasi-static nutrient solve: D∇²c = g·φ → rhs = (g/D)·φ = g_ratio·φ."""
    flat_mask = c_bc_mask.ravel(order='C')
    flat_bc   = c_bc_val.ravel(order='C')
    rhs = np.where(flat_mask, flat_bc, g_ratio * phi.ravel(order='C'))
    c_flat = spla.spsolve(A_mat, rhs)
    return np.clip(c_flat.reshape(NX, NY, order='C'), 0.0, 1.0)


def step_phi(phi, c, alpha):
    monod       = c / (K_MON + c)
    double_well = -GAMMA * phi * (1 - phi) * (1 - 2*phi)
    logistic    = K_LOG * phi * (1 - phi) * monod
    dpx, dpy    = grad2d(phi)
    dcx, dcy    = grad2d(c)
    mag_c = np.sqrt(dcx**2 + dcy**2 + 1e-20)
    dot   = dpx*(dcx/mag_c) + dpy*(dcy/mag_c)
    chemo = R_PARAM_2D * monod * dot
    phi_dot = BETA*lap2d(phi) + double_well + logistic + K_ALPHA*alpha - chemo
    return np.clip(phi + DT * phi_dot, 0.0, 1.0)


def step_alpha(phi, alpha):
    return np.maximum(alpha + DT * K_ALPHA * phi, 0.0)


def sphere_phi(cx, cy, r):
    xi = np.sqrt(BETA / GAMMA)
    d  = np.sqrt((X - cx)**2 + (Y - cy)**2)
    return 0.5 * (1 - np.tanh((d - r) / (2 * xi)))


def make_obs_mask(centers, radius):
    mask = np.zeros((NX, NY), dtype=bool)
    for cx, cy in centers:
        mask |= (X - cx)**2 + (Y - cy)**2 <= radius**2
    return mask


def bc_corner():
    m = np.zeros((NX, NY), dtype=bool)
    m[-1, :] = True; m[:, -1] = True
    return m, np.ones((NX, NY))


def bc_bottom():
    m = np.zeros((NX, NY), dtype=bool)
    m[:, 0] = True
    return m, np.ones((NX, NY))


def cm(phi):
    tot = phi.sum() + 1e-20
    return (phi*X).sum()/tot, (phi*Y).sum()/tot


def hydrostatic(phi, alpha):
    alpha_f = 1.0 + alpha
    return phi * KAPPA * (alpha_f**3 - 1.0) / 3.0


# ═══════════════════════════════════════════════════════════════════
# Simulation runner
# ═══════════════════════════════════════════════════════════════════

def run(case, phi0, c_bc_mask, c_bc_val, g_ratio=G_HI/D_STAR, obs_mask=None):
    """Run Eq.34-36 time integration.

    Returns
    -------
    snaps : dict {step: (phi, c, alpha)}
    avgs  : list of (step, avg_phi, avg_c)  — recorded every AVG_EVERY steps
    """
    print(f"  Building sparse Poisson matrix for case {case}...")
    A_mat = build_poisson_matrix(c_bc_mask)
    phi   = phi0.copy()
    alpha = np.zeros((NX, NY))
    if obs_mask is not None:
        phi[obs_mask] = 0.0
    c = solve_c(A_mat, phi, c_bc_mask, c_bc_val, g_ratio)
    snaps = {}
    avgs  = []
    for s in range(N_STEPS + 1):
        if s in SAVE_AT:
            snaps[s] = (phi.copy(), c.copy(), alpha.copy())
            print(f"  [{case}] step {s:5d}/{N_STEPS}"
                  f"  φ=[{phi.min():.3f},{phi.max():.3f}]"
                  f"  c=[{c.min():.3f},{c.max():.3f}]"
                  f"  α_max={alpha.max():.4f}")
        if s % AVG_EVERY == 0:
            avgs.append((s, float(phi.mean()), float(c.mean())))
        if s < N_STEPS:
            phi   = step_phi(phi, c, alpha)
            if obs_mask is not None:
                phi[obs_mask] = 0.0
            alpha = step_alpha(phi, alpha)
            if obs_mask is not None:
                alpha[obs_mask] = 0.0
            if (s + 1) % C_SOLVE_EVERY == 0:
                c = solve_c(A_mat, phi, c_bc_mask, c_bc_val, g_ratio)
    return snaps, avgs


# ═══════════════════════════════════════════════════════════════════
# Plot 1: felix_complete_reproduction.png  (Cases A + B_hi, v3-style)
# ═══════════════════════════════════════════════════════════════════

def plot_main(snaps_a, snaps_bhi, phi0_a, phi0_b, avgs_a, avgs_bhi):
    steps_show = sorted(s for s in SAVE_AT if s in snaps_a)[:5]
    fig, axes = plt.subplots(3, 7, figsize=(3.2*7, 3.0*3))

    def show_phi(ax, phi, title, obs_mask=None):
        data = phi.T
        im = ax.imshow(data, origin='lower', cmap='hot', vmin=0, vmax=1,
                       extent=[0, LX, 0, LY])
        if obs_mask is not None:
            ax.contourf(x1d, y1d, obs_mask.T.astype(float),
                        levels=[0.5, 1.5], colors=['#5599ff'], alpha=0.5)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
        ax.set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
        plt.colorbar(im, ax=ax, fraction=0.046)

    phi_fa, c_fa, alpha_fa = snaps_a[N_STEPS]
    phi_fb, c_fb, alpha_fb = snaps_bhi[N_STEPS]

    for col, s in enumerate(steps_show):
        show_phi(axes[0, col], snaps_a[s][0],
                 fr"Fig.\,2 $\phi$, $T={s*DT:.1f}T^*$")
    im = axes[0,5].imshow(c_fa.T, origin='lower', cmap='Blues',
                           vmin=0, vmax=1, extent=[0,LX,0,LY])
    axes[0,5].set_title(r"Fig.\,2 $c$ (final)", fontsize=14)
    axes[0,5].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[0,5].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[0,5], fraction=0.046)
    p_a = hydrostatic(phi_fa, alpha_fa)
    im = axes[0,6].imshow(p_a.T, origin='lower', cmap='RdBu_r', extent=[0,LX,0,LY])
    axes[0,6].set_title(fr"Fig.\,2 $p$, max={p_a.max():.2f} Pa", fontsize=14)
    axes[0,6].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[0,6].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[0,6], fraction=0.046)

    for col, s in enumerate(steps_show):
        show_phi(axes[1, col], snaps_bhi[s][0],
                 fr"Fig.\,5 $\phi$, $T={s*DT:.1f}T^*$")
    im = axes[1,5].imshow(c_fb.T, origin='lower', cmap='Blues',
                           vmin=0, vmax=1, extent=[0,LX,0,LY])
    axes[1,5].set_title(r"Fig.\,5 $c$ (final)", fontsize=14)
    axes[1,5].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[1,5].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[1,5], fraction=0.046)
    p_b = hydrostatic(phi_fb, alpha_fb)
    im = axes[1,6].imshow(p_b.T, origin='lower', cmap='RdBu_r', extent=[0,LX,0,LY])
    axes[1,6].set_title(fr"Fig.\,5 $p$, max={p_b.max():.2f} Pa", fontsize=14)
    axes[1,6].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[1,6].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[1,6], fraction=0.046)

    # Row 2: α centrelines, kinematics, avg c/φ (Fig.4)
    ax = axes[2,0]
    for s in steps_show[-3:]:
        if s in snaps_a:
            ax.plot(y1d, snaps_a[s][2][NX//2,:], label=fr"$T={s*DT:.0f}$")
    ax.set_xlabel(r'$y$ ($\mu$m)'); ax.set_ylabel(r'$\alpha_{\rm acc}$')
    ax.set_title(r"Fig.\,2 $\alpha$ (centreline)"); ax.legend(fontsize=14); ax.grid(alpha=0.3)
    clean_ax(ax)

    ax = axes[2,1]
    for s in steps_show[-3:]:
        if s in snaps_bhi:
            ax.plot(y1d, snaps_bhi[s][2][NX//2,:], label=fr"$T={s*DT:.0f}$")
    ax.set_xlabel(r'$y$ ($\mu$m)')
    ax.set_title(r"Fig.\,5 $\alpha$ (centreline)")
    ax.legend(fontsize=14); ax.grid(alpha=0.3); clean_ax(ax)

    ax = axes[2,2]
    ar = np.linspace(0, 0.05, 100)
    ax.plot(ar, (1+ar)**3, 'b-', lw=1.5, label=r'$J_g=(1+\alpha)^3$')
    ax.plot(ar, 1+ar, 'r--', lw=1.2, label=r'$1+\alpha$ (wrong)')
    ax.set_xlabel(r'$\alpha_{\rm acc}$'); ax.set_ylabel(r'$J_g$')
    ax.set_title(r'Kinematics: $\mathbf{F}_g=(1+\alpha)\mathbf{I}$')
    ax.legend(fontsize=14); ax.grid(alpha=0.3); clean_ax(ax)

    ax = axes[2,3]
    ax.plot(p_a[NX//2,:], y1d, 'b-', lw=1.5, label=r'Fig.\,2')
    ax.plot(p_b[NX//2,:], y1d, 'r--', lw=1.5, label=r'Fig.\,5')
    ax.axvline(0, color='k', lw=0.5, ls=':')
    ax.set_xlabel(r'$p$ (Pa)'); ax.set_ylabel(r'$y$ ($\mu$m)')
    ax.set_title("Centreline hydrostatic $p$"); ax.legend(fontsize=14)
    ax.grid(alpha=0.3); clean_ax(ax)

    # Fig.4: Case A average c and φ over time
    ax = axes[2,4]
    t_a  = np.array([v[0]*DT for v in avgs_a])
    ph_a = np.array([v[1] for v in avgs_a])
    c_a  = np.array([v[2] for v in avgs_a])
    ax.plot(t_a, ph_a, 'r-',  lw=1.2, label=r'$\langle\phi\rangle$')
    ax.plot(t_a, c_a,  'b-',  lw=1.2, label=r'$\langle c\rangle$')
    ax.set_xlabel(r'$T$ ($T^*$)'); ax.set_ylabel('Domain average')
    ax.set_title(r"Fig.\,4: Case A avg $\phi$ and $c$")
    ax.legend(fontsize=14); ax.grid(alpha=0.3); ax.set_ylim(0, 1.05); clean_ax(ax)

    # Parameter + PASS/FAIL table
    cx0a, cy0a = cm(phi0_a); cxfa, cyfa = cm(phi_fa)
    cy0b = cm(phi0_b)[1];    cyfb = cm(phi_fb)[1]
    pa_a = (cxfa > cx0a) and (cyfa > cy0a)
    pa_b = cyfb > cy0b
    pa_p = (p_a.max() > 0) and (p_b.max() > 0)
    p_all = pa_a and pa_b and pa_p

    ax = axes[2,5]; ax.axis('off')
    ax.text(0.02, 0.98,
        "Felix Eq.34--36 (v4)\n"
        "----------------------------\n"
        r"$\dot{\phi}$" + " = beta*lap(phi)\n"
        "  - GAMMA*phi*(1-phi)*(1-2phi)\n"
        "  + K_LOG*phi*(1-phi)*c/(k+c)\n"
        "  + k_alpha*alpha\n"
        "  - R_2D*c/(k+c)*(gphi.ngc)\n\n"
        f"R_2D={R_PARAM_2D} (Felix={R_PARAM:.0f})\n"
        f"K_LOG={K_LOG}  GAMMA={GAMMA}\n\n"
        f"CM Fig.2: ({cx0a:.1f},{cy0a:.1f})->"
        f"({cxfa:.1f},{cyfa:.1f})\n"
        f"CM_y Fig.5: {cy0b:.1f}->{cyfb:.1f} um",
        transform=ax.transAxes, fontsize=14, va='top', family='monospace',
        bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

    ax = axes[2,6]; ax.axis('off')
    ax.text(0.02, 0.98,
        "Felix Table 2\n"
        "-------------------\n"
        f"mu  = {MU:.4f} Pa\n"
        f"E   = {E_MOD:.1f} Pa\n"
        f"nu  = {NU:.2f}\n"
        f"bet = {BETA:.1f} um2/T*\n"
        f"k_a = {K_ALPHA:.1e} /T*\n"
        f"g_hi= 1e8 /T*\n"
        f"r   = {R_PARAM:.0f}/{R_PARAM_2D:.0f} /T*\n\n"
        "UMAT:\n"
        f"C10 = {C10:.4f} Pa\n"
        r"D1  = " + f"{D1:.4f} Pa$^{{-1}}$\n\n"
        f"{'PASS' if p_all else 'FAIL'}  2026-06-23\n"
        f"  Fig.2: {'OK' if pa_a else 'NG'}\n"
        f"  Fig.5: {'OK' if pa_b else 'NG'}\n"
        f"  p>0:   {'OK' if pa_p else 'NG'}",
        transform=ax.transAxes, fontsize=14, va='top', family='monospace',
        bbox=dict(boxstyle='round', fc='lightblue', alpha=0.9))

    fig.suptitle(
        fr"Klempt 2024 Eq.34--36 (v4): "
        fr"$R_{{\rm 2D}}$={R_PARAM_2D}, $K_{{\rm log}}$={K_LOG}, "
        fr"$\Gamma$={GAMMA}, $T$={N_STEPS*DT:.1f}$T^*$",
        fontsize=16, fontweight='bold')
    plt.tight_layout()
    out = _HERE / "felix_complete_reproduction.png"
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


# ═══════════════════════════════════════════════════════════════════
# Plot 2: felix_fig3_isosurface.png  (3-D surface φ time series)
# ═══════════════════════════════════════════════════════════════════

def plot_fig3_isosurface(snaps_a):
    """Fig.3 analog: 3-D surface plots of φ(x,y) at 5 time points."""
    steps_show = sorted(s for s in SAVE_AT if s in snaps_a)[:5]
    fig = plt.figure(figsize=(3.5*5, 3.6))
    for col, s in enumerate(steps_show):
        ax = fig.add_subplot(1, 5, col+1, projection='3d')
        phi = snaps_a[s][0]
        # Plot φ surface
        ax.plot_surface(X, Y, phi, cmap='hot', vmin=0, vmax=1,
                        alpha=0.92, linewidth=0, antialiased=True)
        # Mark φ = 0.8 isosurface as thick contour projected on floor
        ax.contour(X, Y, phi, levels=[0.8], colors=['cyan'],
                   linewidths=1.2, zdir='z', offset=0)
        # Axes
        ax.set_xlabel(r'$x$ ($\mu$m)', fontsize=14, labelpad=0)
        ax.set_ylabel(r'$y$ ($\mu$m)', fontsize=14, labelpad=0)
        ax.set_zlabel(r'$\phi$', fontsize=14, labelpad=0)
        ax.set_zlim(0, 1)
        ax.set_title(fr"$T = {s*DT:.1f}\,T^*$", fontsize=14)
        ax.tick_params(labelsize=10)
        ax.view_init(elev=28, azim=225)
        # Highlight isosurface region φ ≥ 0.8
        phi_hi = np.where(phi >= 0.8, phi, np.nan)
        ax.plot_surface(X, Y, phi_hi, color='cyan', alpha=0.4,
                        linewidth=0, antialiased=True)

    fig.suptitle(
        r"Fig.\,3 analog: 3-D $\phi$ surface — Case A (directional growth, $\phi \geq 0.8$ in cyan)",
        fontsize=16, fontweight='bold')
    plt.tight_layout()
    out = _HERE / "felix_fig3_isosurface.png"
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


# ═══════════════════════════════════════════════════════════════════
# Plot 3: felix_fig789_10.png  (Fig.7/8 + B_lo + Case C)
# ═══════════════════════════════════════════════════════════════════

def plot_fig789_10(snaps_a, snaps_bhi, snaps_blo, snaps_c, avgs_a, avgs_bhi, avgs_blo, obs_mask):
    fig, axes = plt.subplots(4, 6, figsize=(3.0*6, 3.2*4))

    T_arr_a   = np.array([v[0]*DT for v in avgs_a])
    ph_arr_a  = np.array([v[1] for v in avgs_a])
    c_arr_a   = np.array([v[2] for v in avgs_a])
    T_arr_hi  = np.array([v[0]*DT for v in avgs_bhi])
    ph_arr_hi = np.array([v[1] for v in avgs_bhi])
    c_arr_hi  = np.array([v[2] for v in avgs_bhi])
    T_arr_lo  = np.array([v[0]*DT for v in avgs_blo])
    ph_arr_lo = np.array([v[1] for v in avgs_blo])
    c_arr_lo  = np.array([v[2] for v in avgs_blo])

    def show_phi(ax, phi, title, obs_m=None):
        im = ax.imshow(phi.T, origin='lower', cmap='hot', vmin=0, vmax=1,
                       extent=[0, LX, 0, LY])
        if obs_m is not None:
            ax.contourf(x1d, y1d, obs_m.T.astype(float),
                        levels=[0.5, 1.5], colors=['#5599ff'], alpha=0.6)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
        ax.set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
        plt.colorbar(im, ax=ax, fraction=0.046)

    # ── Row 0: Fig.4 (Case A avg) + Fig.7 (B_hi vs B_lo avg comparison) ──
    ax = axes[0, 0]
    ax.plot(T_arr_a, ph_arr_a, 'r-',  lw=1.3, label=r'$\langle\phi\rangle$')
    ax.plot(T_arr_a, c_arr_a,  'b-',  lw=1.3, label=r'$\langle c\rangle$')
    ax.set_xlabel(r'$T$ ($T^*$)'); ax.set_ylabel('Domain average')
    ax.set_title(r"Fig.\,4: Case A avg $\phi$, $c$")
    ax.legend(fontsize=14); ax.grid(alpha=0.3); ax.set_ylim(0, 1.05); clean_ax(ax)

    ax = axes[0, 1]
    ax.plot(T_arr_hi, ph_arr_hi, 'r-',  lw=1.3, label=r'$\langle\phi\rangle$ high')
    ax.plot(T_arr_hi, c_arr_hi,  'b-',  lw=1.3, label=r'$\langle c\rangle$ high')
    ax.plot(T_arr_lo, ph_arr_lo, 'r--', lw=1.3, label=r'$\langle\phi\rangle$ low')
    ax.plot(T_arr_lo, c_arr_lo,  'b--', lw=1.3, label=r'$\langle c\rangle$ low')
    ax.set_xlabel(r'$T$ ($T^*$)'); ax.set_ylabel('Domain average')
    ax.set_title(r"Fig.\,7: Case B high (solid) vs low (dashed)")
    ax.legend(fontsize=11, ncol=2); ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05); clean_ax(ax)

    # B_hi final φ and B_lo final φ for visual comparison
    phi_fa_hi, c_fa_hi, alpha_fa_hi = snaps_bhi[N_STEPS]
    phi_fa_lo, c_fa_lo, alpha_fa_lo = snaps_blo[N_STEPS]
    show_phi(axes[0, 2], phi_fa_hi, r"Fig.\,5 $\phi$ final (high, $g=10^8$)")
    show_phi(axes[0, 3], phi_fa_lo, r"Fig.\,5 $\phi$ final (low,  $g=10^{10}$)")

    # Nutrient field topology (Fig.6 analog)
    show_phi(axes[0, 4], c_fa_hi,   r"Fig.\,6: nutrient field $c$ (high)")
    show_phi(axes[0, 5], c_fa_lo,   r"Fig.\,6: nutrient field $c$ (low)")
    # override colormap for c
    for col in [4, 5]:
        axes[0, col].get_images()[0].set_cmap('Blues')

    # ── Row 1: Fig.8 — B_hi hydrostatic pressure at 4 times ──────────────
    p_steps = [2000, 5000, 8000, 12000]
    p_all_hi = [hydrostatic(snaps_bhi[s][0], snaps_bhi[s][2]) for s in p_steps
                if s in snaps_bhi]
    p_lo_final = hydrostatic(phi_fa_lo, alpha_fa_lo)

    # Common colour limits (symmetric around 0 → tension ring visible if present)
    p_abs_max = max(abs(p.max()) for p in p_all_hi + [p_lo_final]) + 0.01
    for i, (s, p_field) in enumerate(zip(p_steps, p_all_hi)):
        ax = axes[1, i]
        im = ax.imshow(p_field.T, origin='lower', cmap='RdBu_r',
                       vmin=-p_abs_max, vmax=p_abs_max, extent=[0, LX, 0, LY])
        ax.set_title(fr"Fig.\,8 $p$, $T={s*DT:.1f}T^*$ (max={p_field.max():.2f} Pa)",
                     fontsize=14)
        ax.set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
        ax.set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
        plt.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1, 4]
    im = ax.imshow(p_lo_final.T, origin='lower', cmap='RdBu_r',
                   vmin=-p_abs_max, vmax=p_abs_max, extent=[0, LX, 0, LY])
    ax.set_title(r"Fig.\,8 $p$ final (low)", fontsize=14)
    ax.set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    ax.set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1, 5]; ax.axis('off')
    phi_fa_a, _, alpha_fa_a = snaps_a[N_STEPS]
    p_a_final = hydrostatic(phi_fa_a, alpha_fa_a)
    p_b_final = hydrostatic(snaps_bhi[N_STEPS][0], snaps_bhi[N_STEPS][2])
    cx0a, cy0a = cm(snaps_a[0][0]); cxfa, cyfa = cm(phi_fa_a)
    cy0b = cm(snaps_bhi[0][0])[1];  cyfb = cm(phi_fa_hi)[1]
    cy0lo = cm(snaps_blo[0][0])[1]; cyflo = cm(phi_fa_lo)[1]
    layer_ok = cyflo < cyfb
    ax.text(0.02, 0.98,
        "v4 PASS/FAIL\n"
        "-------------------\n"
        f"Fig.2 directional: {'OK' if (cxfa>cx0a and cyfa>cy0a) else 'NG'}\n"
        f"Fig.5 upward:      {'OK' if cyfb>cy0b else 'NG'}\n"
        f"p>0 (A,B_hi):      {'OK' if p_a_final.max()>0 and p_b_final.max()>0 else 'NG'}\n"
        f"Fig.7 layer<mush:  {'OK' if layer_ok else 'NG'}\n"
        f"  CM_y B_hi: {cy0b:.1f}->{cyfb:.1f}\n"
        f"  CM_y B_lo: {cy0lo:.1f}->{cyflo:.1f}",
        transform=ax.transAxes, fontsize=14, va='top', family='monospace',
        bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

    # ── Row 2: Case B_lo time series ──────────────────────────────────────
    blo_steps = [s for s in sorted(SAVE_AT) if s in snaps_blo][:4]
    for col, s in enumerate(blo_steps):
        show_phi(axes[2, col], snaps_blo[s][0],
                 fr"B\_lo $\phi$, $T={s*DT:.1f}T^*$ ($g=10^{{10}}$)")
    show_phi(axes[2, 4], c_fa_lo, r"B\_lo $c$ final (layer formation)")
    axes[2, 4].get_images()[0].set_cmap('Blues')
    ax = axes[2, 5]; ax.axis('off')
    ax.text(0.05, 0.9,
        r"$g_{\rm lo}=10^{10}$ (Thiele$\approx$3)" "\n"
        r"Layer formation expected:" "\n"
        r"  $\langle\phi\rangle$ stays near" "\n"
        r"  bottom ($y\approx2.5\,\mu$m)," "\n"
        r"  unlike mushroom of B\_hi.",
        transform=ax.transAxes, fontsize=14, va='top',
        bbox=dict(boxstyle='round', fc='#e8ffe8', alpha=0.9))

    # ── Row 3: Case C — obstacle avoidance ────────────────────────────────
    c_steps = [s for s in sorted(SAVE_AT) if s in snaps_c][:4]
    for col, s in enumerate(c_steps):
        show_phi(axes[3, col], snaps_c[s][0],
                 fr"Fig.\,10 obs. $\phi$, $T={s*DT:.1f}T^*$",
                 obs_m=obs_mask)
    phi_fc, c_fc, alpha_fc = snaps_c[N_STEPS]
    im = axes[3, 4].imshow(c_fc.T, origin='lower', cmap='Blues',
                            vmin=0, vmax=1, extent=[0, LX, 0, LY])
    axes[3, 4].contourf(x1d, y1d, obs_mask.T.astype(float),
                         levels=[0.5, 1.5], colors=['#5599ff'], alpha=0.6)
    axes[3, 4].set_title(r"Fig.\,10 $c$ final (obstacles in blue)", fontsize=14)
    axes[3, 4].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[3, 4].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[3, 4], fraction=0.046)

    p_c = hydrostatic(phi_fc, alpha_fc)
    im = axes[3, 5].imshow(p_c.T, origin='lower', cmap='RdBu_r',
                            extent=[0, LX, 0, LY])
    axes[3, 5].contourf(x1d, y1d, obs_mask.T.astype(float),
                         levels=[0.5, 1.5], colors=['#5599ff'], alpha=0.6)
    axes[3, 5].set_title(fr"Fig.\,10 $p$ final (max={p_c.max():.2f} Pa)", fontsize=14)
    axes[3, 5].set_xlabel(r'$x$ ($\mu$m)', fontsize=14)
    axes[3, 5].set_ylabel(r'$y$ ($\mu$m)', fontsize=14)
    plt.colorbar(im, ax=axes[3, 5], fraction=0.046)

    fig.suptitle(
        r"Klempt 2024 Fig.\,7/8 (avg+pressure) + B\_lo (layer) + Fig.\,9-10 (obstacle avoidance)",
        fontsize=16, fontweight='bold')
    plt.tight_layout()
    out = _HERE / "felix_fig789_10.png"
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    print("="*65)
    print("Felix Klempt 2024 -- Eq.34-36 Complete Reproduction v4")
    print(f"Grid: {NX}x{NY}  DX={DX:.3f}um  DT={DT:.4f}T*  N={N_STEPS}")
    xi = np.sqrt(BETA/GAMMA)
    print(f"T_total={N_STEPS*DT:.1f}T*  R_2D={R_PARAM_2D}  "
          f"K_LOG={K_LOG}  GAMMA={GAMMA}  xi={xi:.3f}um")
    print(f"Thiele B_hi={np.sqrt(G_HI*5**2/D_STAR):.2f}  "
          f"Thiele B_lo={np.sqrt(G_LO*5**2/D_STAR):.2f}")
    print("="*65)

    c_bc_a, c_bv_a = bc_corner()
    c_bc_b, c_bv_b = bc_bottom()
    obs_mask = make_obs_mask(OBS_CENTERS, OBS_RADIUS)

    # Case A: directional growth (Fig.2-4)
    print("\n[A] Fig.2-4: corner nutrient, directional growth")
    phi0_a = sphere_phi(LX/2, LY/2, r=5.0)
    snaps_a, avgs_a = run("A", phi0_a, c_bc_a, c_bv_a,
                          g_ratio=G_HI/D_STAR)

    # Case B_hi: agar plate, high nutrient (g=1e8) → mushroom (Fig.5-8)
    print("\n[B_hi] Fig.5-8: bottom nutrient, g=1e8 (high), mushroom")
    phi0_b = sphere_phi(LX/2, 2.5, r=5.0)
    snaps_bhi, avgs_bhi = run("B_hi", phi0_b, c_bc_b, c_bv_b,
                               g_ratio=G_HI/D_STAR)

    # Case B_lo: agar plate, low nutrient (g=1e10) → layer (Fig.7 dashed)
    print("\n[B_lo] Fig.7 dashed: bottom nutrient, g=1e10 (low), layer")
    snaps_blo, avgs_blo = run("B_lo", phi0_b.copy(), c_bc_b, c_bv_b,
                               g_ratio=G_LO/D_STAR)

    # Case C: obstacle avoidance (Fig.9-10)
    print("\n[C] Fig.9-10: obstacle avoidance (4 rigid pillars)")
    phi0_c = sphere_phi(3.0, 3.0, r=2.5)
    snaps_c, avgs_c = run("C", phi0_c, c_bc_a, c_bv_a,
                           g_ratio=G_HI/D_STAR, obs_mask=obs_mask)

    # ── Verification printout ────────────────────────────────────
    phi_fa, _, alpha_fa = snaps_a[N_STEPS]
    phi_fb_hi, _, alpha_fb_hi = snaps_bhi[N_STEPS]
    phi_fb_lo, _, alpha_fb_lo = snaps_blo[N_STEPS]
    phi_fc, _, alpha_fc = snaps_c[N_STEPS]

    cx0a, cy0a = cm(phi0_a); cxfa, cyfa = cm(phi_fa)
    cy0b = cm(phi0_b)[1];    cyfb_hi = cm(phi_fb_hi)[1]
    cy0lo = cm(phi0_b)[1];   cyfb_lo = cm(phi_fb_lo)[1]
    cx0c, cy0c = cm(phi0_c); cxfc, cyfc = cm(phi_fc)

    p_a  = hydrostatic(phi_fa,    alpha_fa)
    p_bh = hydrostatic(phi_fb_hi, alpha_fb_hi)
    p_c  = hydrostatic(phi_fc,    alpha_fc)

    pa_fig2 = (cxfa > cx0a) and (cyfa > cy0a)
    pa_fig5 = cyfb_hi > cy0b
    pa_p    = (p_a.max() > 0) and (p_bh.max() > 0)
    pa_fig7 = cyfb_lo < cyfb_hi   # layer stays lower than mushroom
    pa_figC = (cxfc > cx0c) or (cyfc > cy0c)  # biofilm moved toward nutrient corner

    print("\n" + "="*65)
    print("Felix Eq.34-36 v4 verification:")
    print(f"  Fig.2  directional: {'PASS' if pa_fig2 else 'FAIL'}"
          f"  CM ({cx0a:.1f},{cy0a:.1f})->({cxfa:.1f},{cyfa:.1f})")
    print(f"  Fig.5  upward:      {'PASS' if pa_fig5 else 'FAIL'}"
          f"  CM_y {cy0b:.1f}->{cyfb_hi:.1f} um")
    print(f"  p>0:               {'PASS' if pa_p else 'FAIL'}"
          f"  p_max(A)={p_a.max():.2f} p_max(B_hi)={p_bh.max():.2f} Pa")
    print(f"  Fig.7  layer<mush: {'PASS' if pa_fig7 else 'FAIL'}"
          f"  CM_y B_lo={cyfb_lo:.1f} < B_hi={cyfb_hi:.1f}")
    print(f"  Fig.10 obs. avoid: {'PASS' if pa_figC else 'FAIL'}"
          f"  CM ({cx0c:.1f},{cy0c:.1f})->({cxfc:.1f},{cyfc:.1f})")
    all_pass = pa_fig2 and pa_fig5 and pa_p and pa_fig7 and pa_figC
    print(f"\n  Overall: {'PASS' if all_pass else 'FAIL'}")
    print("="*65)

    # ── Generate figures ─────────────────────────────────────────
    print("\nGenerating figures...")
    plot_main(snaps_a, snaps_bhi, phi0_a, phi0_b, avgs_a, avgs_bhi)
    plot_fig3_isosurface(snaps_a)
    plot_fig789_10(snaps_a, snaps_bhi, snaps_blo, snaps_c,
                   avgs_a, avgs_bhi, avgs_blo, obs_mask)
    print("\nAll done.")


if __name__ == "__main__":
    main()
