"""anim_nutrient_3d.py -- the Chapter 4 seed demonstration with a nutrient field (demonstration, not a thesis result).

Field as anim_scenarios_3d.py (48^3 cells of the 2 mm cube, RK4), plus the nutrient c of Klempt et al. 2024:

    d phi/dt   = beta lap(phi) + k_alpha alpha - penalty (max(phi - 1, 0) + min(phi, 0))     Eq. 34 without the front term
    0          = d lap(c) - g phi c,   c = 1 on the top face z = +1, zero flux elsewhere        Eq. 35, quasi-static,
                                                                                                consumption g phi c
    d alpha/dt = k_alpha phi c                                                                  Eq. 36 with the factor c

beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, T* = 1.1. The nutrient enters through the ratio g/d only:
Table 2 (d = 1e10 um^2/T*, g = 1e8 /T*) converted to the 2 mm cube like beta (KLEMPT2024_REPRODUCTION.md,
klempt2024_cases.deck_values) gives g/d = 1 /mm^2, i.e. a penetration depth (d/g)^1/2 = 1 mm. The nutrient is
solved quasi-statically (the paper's time homogenisation, eta_c = 1e-10): with the converted d = 1e8 mm^2/T* an
explicit step would need dt < 4e-12 T*, so the steady equation is solved by conjugate gradients once per step.
Mechanics as anim_scenarios_3d.py: E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, small strain,
hex8 B-bar on 16^3, bottom face clamped, solved once per frame. Both runs, with the nutrient and with c = 1 held
(no limitation), are made in one go; the video shows the run with the nutrient.

  corner      seed x < 0, y < 0, z < -0.5 (the corner seed of the law case), g/d = 1 /mm^2
  layer       flat layer z < -0.75 over the whole bottom face, g/d = 1 /mm^2
  corner_g10  as corner with g x 10 (g = 1e9 /T*, between the paper's two Fig. 7 values 1e8 and 1e10)
  layer_g10   as layer with g x 10

    python ansys_usermat/apdl/anim_nutrient_3d.py corner|layer|corner_g10|layer_g10 [--no-video]
        -> assets/anim_nutrient_<case>_3d.mp4, .png; cache figdata/anim_nutrient_<case>_3d.npz;
           ansys_usermat/apdl/results/2026-10-3d/nutrient_<case>.json
"""
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")        # threaded BLAS makes the 110k-vector dot products of the CG 100x slower
import numpy as np  # noqa: E402
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

CASE = sys.argv[1]
NO_VIDEO = "--no-video" in sys.argv
NF, NFRAMES, NM = 48, 88, 16
BETA, KA, PEN, TEND = 0.02, 1e-3, 5.0, 1.1
E_BIO, FLOOR, NU = 10.0, 1e-3, 0.49
G_TABLE2, D_TABLE2, L_PAPER_UM, L_CUBE_MM = 1e8, 1e10, 20.0, 2.0          # Klempt 2024 Table 2, 20 um cube
G_OVER_D_TABLE2 = G_TABLE2 / (D_TABLE2 * (L_CUBE_MM * 1e3 / L_PAPER_UM) ** 2 * 1e-6)   # 1/mm^2 on the 2 mm cube
G_FACTOR = 10.0 if CASE.endswith("_g10") else 1.0
Q = G_OVER_D_TABLE2 * G_FACTOR                                           # g/d in 1/mm^2
CACHE = ROOT / "figdata" / f"anim_nutrient_{CASE}_3d.npz"
RESDIR = HERE / "results" / "2026-10-3d"


def up(a, n):
    f = n // a.shape[0]
    return np.repeat(np.repeat(np.repeat(a, f, 0), f, 1), f, 2)


def down(a, n):
    f = a.shape[0] // n
    return a.reshape(n, f, n, f, n, f).mean((1, 3, 5))


# ---- geometry and seeds on the field grid
ax = -1 + (np.arange(NF) + 0.5) * 2 / NF
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
base = CASE.replace("_g10", "")
if base == "corner":
    seed = (X < 0) & (Y < 0) & (Z < -0.5)
elif base == "layer":
    seed = Z < -0.75
else:
    raise SystemExit("case: corner, layer, corner_g10 or layer_g10")
seed_m = down(seed.astype(float), NM) > 0.5
H = 2.0 / NF


def lap(u):
    """7-point Laplacian with zero flux at all six walls (mirror cells)."""
    out = -6.0 * u
    out[1:] += u[:-1]; out[:-1] += u[1:]; out[0] += u[0]; out[-1] += u[-1]
    out[:, 1:] += u[:, :-1]; out[:, :-1] += u[:, 1:]; out[:, 0] += u[:, 0]; out[:, -1] += u[:, -1]
    out[:, :, 1:] += u[:, :, :-1]; out[:, :, :-1] += u[:, :, 1:]; out[:, :, 0] += u[:, :, 0]; out[:, :, -1] += u[:, :, -1]
    return out / (H * H)


# ---- quasi-static nutrient: -lap(c) + q phi c = 0, c = 1 on z = +1 (ghost cell 2 - c), zero flux elsewhere
def solve_c(phi, c0):
    q = Q * np.clip(phi, 0.0, 1.0)
    n3 = NF ** 3

    def matvec(v):
        u = v.reshape((NF,) * 3)
        out = -lap(u) + q * u
        out[:, :, -1] += 2.0 * u[:, :, -1] / (H * H)        # Dirichlet ghost: (2 cD - u) replaces the mirror value u
        return out.ravel()

    diag = (6.0 / (H * H) + q).ravel()
    diag.reshape((NF,) * 3)[:, :, -1] += 2.0 / (H * H)
    b = np.zeros((NF,) * 3); b[:, :, -1] = 2.0 / (H * H); b = b.ravel()
    A = spla.LinearOperator((n3, n3), matvec=matvec, dtype=float)
    M = spla.LinearOperator((n3, n3), matvec=lambda v: v / diag, dtype=float)
    it = [0]

    def cb(_):
        it[0] += 1
    c, info = spla.cg(A, b, x0=c0.ravel(), M=M, rtol=1e-10, maxiter=5000, callback=cb)
    res = np.linalg.norm(matvec(c) - b) / np.linalg.norm(b)
    return c.reshape((NF,) * 3), it[0], float(res)


def mech(phi_m, am1_m):
    n = NM
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, NU, False)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    Ee = E_BIO * (np.clip(phi_m.ravel(), 0, 1) ** 2 + FLOOR)
    eps0 = am1_m.ravel()
    rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    ndof = 3 * (n + 1) ** 3
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(ndof), fix)
    u = np.zeros(ndof)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = -sig[:, :3].mean(1)
    s = sig[:, :3] + p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape((n,) * 3), p.reshape((n,) * 3), u.reshape((n + 1,) * 3 + (3,))


def rhs(phi, alpha, c):
    dphi = BETA * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
    return dphi, KA * phi * c


def run(limited):
    """RK4 on the field grid, c quasi-static once per step (limited) or c = 1 (not limited).
    Returns the per-frame lists (stress once per frame)."""
    dt = 0.15 * H * H / BETA
    per = int(np.ceil(TEND / NFRAMES / dt))
    dt = TEND / (NFRAMES * per)
    phi = np.where(seed, 1.0, 0.0)
    alpha = np.ones_like(phi)
    c = np.ones_like(phi)
    out = dict(t=[], phi=[], am1=[], c=[], vm=[], p=[], u=[], cg_it=[], cg_res=[])
    t0 = time.time()
    for fr in range(NFRAMES + 1):
        if limited:
            c, it, res = solve_c(phi, c)
            out["cg_it"].append(it); out["cg_res"].append(res)
        if limited or fr == NFRAMES:
            vm, p, u = mech(down(phi, NM), down(alpha - 1, NM))
        else:
            vm = p = np.full((NM,) * 3, np.nan); u = np.full((NM + 1,) * 3 + (3,), np.nan)
        out["t"].append(fr * per * dt); out["phi"].append(phi.astype(np.float32))
        out["am1"].append((alpha - 1).astype(np.float32)); out["c"].append(c.astype(np.float32))
        out["vm"].append(vm.astype(np.float32)); out["p"].append(p.astype(np.float32)); out["u"].append(u.astype(np.float32))
        print(f"{CASE} {'nutrient' if limited else 'c = 1   '} frame {fr}/{NFRAMES}  "
              f"c in biofilm [{c[phi > 0.5].min():.3f}, {c[phi > 0.5].max():.3f}]  max vM {np.nanmax(vm):.3e}  "
              f"{time.time() - t0:.0f} s", flush=True)
        if fr == NFRAMES:
            break
        for _ in range(per):
            k1 = rhs(phi, alpha, c)
            k2 = rhs(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1], c)
            k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1], c)
            k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1], c)
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    out["dt"] = dt; out["per"] = per
    return out


if not CACHE.exists():
    L = run(True)
    N0 = run(False)
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, t=np.array(L["t"]), phi=np.array(L["phi"]), am1=np.array(L["am1"]), c=np.array(L["c"]),
                        vm=np.array(L["vm"]), p=np.array(L["p"]), u=np.array(L["u"]),
                        cg_it=np.array(L["cg_it"]), cg_res=np.array(L["cg_res"]), dt=L["dt"], per=L["per"],
                        phi0=N0["phi"][-1], am1_0=N0["am1"][-1], vm0=N0["vm"][-1], p0=N0["p"][-1], u0=N0["u"][-1],
                        q=Q)

D = np.load(CACHE)
t, PHI, AM1, C, VM, P, U = D["t"], D["phi"], D["am1"], D["c"], D["vm"], D["p"], D["u"]


# ---------------------------------------------------------------- numbers
def stats(phi_f, am1_f, c_f, vm, p):
    phi_m = down(np.asarray(phi_f, float), NM)
    bio_f = np.asarray(phi_f) > 0.5
    bio_m = phi_m > 0.5
    am1_m = down(np.asarray(am1_f, float), NM)
    kz = np.arange(NM)
    sel = seed_m
    imax = np.unravel_index(np.argmax(np.where(bio_m, vm, -1)), vm.shape)
    zc = -1 + (np.arange(NM) + 0.5) * 2 / NM
    out = {
        "seed_mean_alpha_minus_1": float(am1_m[sel].mean()),
        "seed_max_alpha_minus_1": float(am1_m[sel].max()),
        "seed_min_alpha_minus_1": float(am1_m[sel].min()),
        "seed_mean_vm_Pa": float(vm[sel].mean()),
        "seed_max_vm_Pa": float(vm[sel].max()),
        "seed_p95_vm_Pa": float(np.percentile(vm[sel], 95)),
        "seed_mean_p_Pa": float(p[sel].mean()),
        "seed_max_p_Pa": float(p[sel].max()),
        "biofilm_c_min": float(np.asarray(c_f)[bio_f].min()),
        "biofilm_c_max": float(np.asarray(c_f)[bio_f].max()),
        "biofilm_c_mean": float(np.asarray(c_f)[bio_f].mean()),
        "domain_c_mean": float(np.asarray(c_f).mean()),
        "biofilm_cells_phi_gt_0.5_field_grid": int(bio_f.sum()),
        "max_vm_in_biofilm_Pa": float(vm[imax]),
        "max_vm_location_mm": [float(zc[imax[0]]), float(zc[imax[1]]), float(zc[imax[2]])],
        # profiles over the stress-mesh layers k (z = -1 + (k + 0.5)/8 mm) within the seed
        "seed_layer_z_mm": [float(zc[k]) for k in kz if sel[:, :, k].any()],
        "seed_layer_mean_alpha_minus_1": [float(am1_m[:, :, k][sel[:, :, k]].mean()) for k in kz if sel[:, :, k].any()],
        "seed_layer_mean_vm_Pa": [float(vm[:, :, k][sel[:, :, k]].mean()) for k in kz if sel[:, :, k].any()],
        "seed_layer_mean_c": [float(down(np.asarray(c_f, float), NM)[:, :, k][sel[:, :, k]].mean())
                              for k in kz if sel[:, :, k].any()],
    }
    return out


with_n = stats(PHI[-1], AM1[-1], C[-1], VM[-1], P[-1])
without = stats(D["phi0"], D["am1_0"], np.ones((NF,) * 3), D["vm0"], D["p0"])
dt = float(D["dt"])
res = {
    "case": CASE,
    "model": {
        "phi": "d phi/dt = beta lap(phi) + k_alpha alpha - penalty (max(phi-1,0) + min(phi,0))  [Eq. 34, no front term]",
        "c": "0 = d lap(c) - g phi c, c = 1 on z = +1, zero flux on the other faces  [Eq. 35 quasi-static, g phi c]",
        "alpha": "d alpha/dt = k_alpha phi c  [Eq. 36 with the factor c]",
        "mechanics": "E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, hex8 B-bar, z = -1 clamped",
    },
    "parameters": {
        "beta_mm2_per_T": BETA, "k_alpha_per_T": KA, "penalty": PEN, "T_end": TEND,
        "g_over_d_per_mm2": float(Q), "g_over_d_table2_converted_per_mm2": float(G_OVER_D_TABLE2),
        "g_factor_vs_table2": G_FACTOR, "penetration_depth_mm": float(Q ** -0.5),
        "E_Pa": E_BIO, "stiffness_floor": FLOOR, "nu": NU,
        "field_grid": NF, "stress_grid": NM, "frames": NFRAMES, "dt_T": dt, "rk4_steps_per_frame": int(D["per"]),
    },
    "stability": {
        "phi_diffusion_number_beta_dt_over_h2": float(BETA * dt / H ** 2),
        "rk4_limit_7point": 2.785 / 12,
        "penalty_dt": float(PEN * dt),
        "c_explicit_diffusion_number_d_dt_over_h2": float(D_TABLE2 * (L_CUBE_MM * 1e3 / L_PAPER_UM) ** 2 * 1e-6 * dt / H ** 2),
        "c_explicit_dt_limit_T": float(2.785 / 12 * H ** 2 / (D_TABLE2 * (L_CUBE_MM * 1e3 / L_PAPER_UM) ** 2 * 1e-6)),
        "c_treatment": "quasi-static, conjugate gradients with Jacobi preconditioner, rtol 1e-10, once per RK4 step",
        "c_reaction_number_q_h2": float(Q * H ** 2),
        "cg_iterations_max": int(D["cg_it"].max()), "cg_relative_residual_max": float(D["cg_res"].max()),
    },
    "with_nutrient": with_n,
    "without_nutrient_c_fixed_1": without,
    "ratio_with_over_without": {k: with_n[k] / without[k] for k in
                                ("seed_mean_alpha_minus_1", "seed_mean_vm_Pa", "seed_max_vm_Pa", "seed_mean_p_Pa")
                                if without[k] != 0},
}
RESDIR.mkdir(parents=True, exist_ok=True)
(RESDIR / f"nutrient_{CASE}.json").write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps({k: res[k] for k in ("with_nutrient", "without_nutrient_c_fixed_1", "ratio_with_over_without")}, indent=1))
if NO_VIDEO:
    sys.exit(0)

# ---------------------------------------------------------------- render
import pyvista as pv  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402

pv.OFF_SCREEN = True
pv.global_theme.font.family = "times"
MAG = 50.0
_pts = {}


def disp(Uf, n):
    if n not in _pts:
        xn = np.linspace(-1, 1, n + 1)
        Xn, Yn, Zn = np.meshgrid(xn, xn, xn, indexing="ij")
        _pts[n] = np.stack([Xn.ravel(order="F"), Yn.ravel(order="F"), Zn.ravel(order="F")], 1)
    xm = np.linspace(-1, 1, Uf.shape[0])
    return np.stack([RegularGridInterpolator((xm,) * 3, Uf[..., c])(_pts[n]) for c in range(3)], 1)


def grid(field, Uf, keep):
    n = field.shape[0]
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1, -1, -1)).cast_to_structured_grid()
    g.points = g.points + MAG * disp(Uf, n)
    g.cell_data["v"] = np.nan_to_num(field).ravel(order="F")
    g.cell_data["keep"] = keep.ravel(order="F").astype(float)
    return g.threshold(0.5, scalars="keep") if not keep.all() else g


# no cut: both seeds touch the faces x = -1 and y = -1 that face the camera, so the walls show c inside the biofilm
all_f = np.ones((NF,) * 3, bool)
all_m = np.ones((NM,) * 3, bool)
pm_end = down(PHI[-1].astype(float), NM)
vtop = float(np.percentile(VM[-1][pm_end > 0.3], 95)) * 1e3
atop = float(np.nanmax(AM1[-1])) * 1e3                           # own colour scale; the note gives the ratio to c = 1
cmin = float(np.floor(C[-1].min() * 10) / 10)
PANELS = [(lambda i: (C[i], U[i], all_f), "nutrient c", "Blues", (cmin, 1.0), "%.2f"),
          (lambda i: (AM1[i] * 1e3, U[i], all_f), "alpha - 1 (1e-3)", "Oranges", (0, atop), "%.2f"),
          (lambda i: (VM[i] * 1e3, U[i], all_m), "von Mises (mPa)", "viridis", (0, vtop), "%.1f")]
CAM = [(-4.0, -3.3, 2.6), (0, 0, -0.1), (0, 0, 1)]
R = res["ratio_with_over_without"]
NOTE = (f"c = 1 on the top face, g/d = {Q:g} /mm^2 (Table 2 x {G_FACTOR:g}); seed mean alpha - 1, von Mises: "
        f"{R['seed_mean_alpha_minus_1']:.2f}, {R['seed_mean_vm_Pa']:.2f} x the run with c = 1; "
        f"bottom clamped; displacement x{MAG:g}")

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
acts = []
for j, (get, lab, cm, clim, fmt) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    a = pl.add_mesh(grid(*get(0)), scalars="v", cmap=cm, clim=clim, ambient=0.55, diffuse=0.5, specular=0.0,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                         width=0.7, height=0.07, label_font_size=24, fmt=fmt, n_labels=5,
                                         font_family="times", color="black"))
    pl.add_text(lab, position="upper_edge", font_size=14, font="times", color="black")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).extract_feature_edges(), color="#9a9a9a", line_width=1.0)
    pl.add_mesh(pv.Plane(center=(0, 0, -1.005), direction=(0, 0, 1), i_size=2.6, j_size=2.6), color="#b0b0b0")
    pl.camera_position = CAM
    pl.camera.zoom(0.82)
    acts.append((a, get))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.subplot(0, 2)
pl.add_text(NOTE, position="lower_right", font_size=7, font="times", color="#555555")
pl.link_views()


def frame(i, azim):
    for a, get in acts:
        a.mapper.SetInputData(grid(*get(i)))
    title.SetText(2, f"T* = {t[i]:.2f}")
    pl.camera.azimuth = azim


out = ROOT / "assets" / f"anim_nutrient_{CASE}_3d"
cam0 = pl.camera.azimuth
pl.open_movie(str(out.with_suffix(".mp4")), framerate=16, quality=9)
nf = len(t)
for i in range(nf):
    frame(i, cam0 + 30.0 * i / (nf - 1) - 15.0)
    pl.write_frame()
for k in range(24):
    frame(nf - 1, cam0 + 15.0 + 15.0 * (k + 1) / 24)
    pl.write_frame()
pl.screenshot(str(out.with_suffix(".png")))
pl.close()
print("wrote", out.with_suffix(".mp4"))
