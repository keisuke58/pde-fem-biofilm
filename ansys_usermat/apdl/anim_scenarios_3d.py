"""anim_scenarios_3d.py -- three 3D videos with the Chapter 4 seed model (demonstrations, not thesis results).

Model as wp3_reference (d phi/dt = beta lap phi + k_alpha alpha - penalty, beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*,
penalty 5, T* = 1.1), stiffness E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, small strain, hex8 B-bar,
bottom face z = -1 clamped. Field on 48^3 (RK4), stress on 16^3 (24^3 for the tooth), staggered: the stress is solved
once per frame (88 frames) and held over the frame's RK4 steps.

  law    seed in the corner (as anim_seed_offcentre_3d.py corner) with the homeostatic-pressure law of P3,
         alpha_dot = k_alpha phi g, g = min(1, max(0, 1 - p/p_h)), p = -tr(sigma)/3, explicit gate per element;
         p_h = 0.4 x the 95th percentile of p in the biofilm at the end of the run without the law.
         Video: phi with the law, von Mises without and with the law on one colour scale.
  two    two seeds of 0.5 x 0.5 x 0.5 mm on the bottom face, 1 mm apart; they spread and meet.
  tooth  a stiff half cylinder (radius 0.7 mm, axis along y on the bottom face, E = 1e4 Pa, no growth, no flux of phi
         into it) with a biofilm layer 0.25 mm thick on one flank of it.

    python ansys_usermat/apdl/anim_scenarios_3d.py law|two|tooth
        -> assets/anim_<case>_3d.mp4, .png; cache figdata/anim_<case>_3d.npz
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

CASE = sys.argv[1]
NF, NFRAMES = 48, 88
NM = 24 if CASE == "tooth" else 16
BETA, KA, PEN, TEND = 0.02, 1e-3, 5.0, 1.1
E_BIO, E_TOOTH, FLOOR, NU = 10.0, 1e4, 1e-3, 0.49
CACHE = ROOT / "figdata" / f"anim_{CASE}_3d.npz"


def up(a, n):
    f = n // a.shape[0]
    return np.repeat(np.repeat(np.repeat(a, f, 0), f, 1), f, 2)


def down(a, n):
    f = a.shape[0] // n
    return a.reshape(n, f, n, f, n, f).mean((1, 3, 5))


# ---- geometry and seeds on the field grid
ax = -1 + (np.arange(NF) + 0.5) * 2 / NF
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
tooth = np.zeros((NF,) * 3, bool)
if CASE == "law":
    seed = (X < 0) & (Y < 0) & (Z < -0.5)
elif CASE == "two":
    seed = (np.abs(Y) < 0.25) & (Z < -0.5) & ((np.abs(X + 0.5) < 0.25) | (np.abs(X - 0.5) < 0.25))
elif CASE == "tooth":
    r = np.sqrt(X ** 2 + (Z + 1) ** 2)
    tooth = r < 0.7
    seed = (r >= 0.7) & (r < 0.95) & (X < 0.15) & (np.abs(Y) < 0.6)
else:
    raise SystemExit("case: law, two or tooth")
tooth_m = down(tooth.astype(float), NM) > 0.5          # tooth elements of the stress mesh
OPEN = ~tooth


def lap(u):
    """7-point Laplacian with zero flux at the walls and at the tooth."""
    out = np.zeros_like(u)
    for axn in range(3):
        for sgn in (1, -1):
            un = np.roll(u, -sgn, axn)
            on = np.roll(OPEN, -sgn, axn)
            edge = [slice(None)] * 3
            edge[axn] = -1 if sgn == 1 else 0
            ok = OPEN & on
            ok[tuple(edge)] = False
            out += np.where(ok, un - u, 0.0)
    return out / (2 / NF) ** 2


def mech(phi_m, am1_m):
    n = NM
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, NU, False)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    Ee = np.where(tooth_m.ravel(), E_TOOTH, E_BIO * (np.clip(phi_m.ravel(), 0, 1) ** 2 + FLOOR))
    eps0 = np.where(tooth_m.ravel(), 0.0, am1_m.ravel())
    rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    ndof = 3 * (n + 1) ** 3
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(ndof), fix)
    u = np.zeros(ndof)
    Kff = K[free][:, free].tocsc()
    if n <= 16:
        u[free] = spla.spsolve(Kff, F[free])
    else:
        dinv = 1.0 / Kff.diagonal()
        u[free], info = spla.cg(Kff, F[free], M=spla.LinearOperator(Kff.shape, matvec=lambda x: dinv * x),
                                rtol=1e-10, maxiter=50000)
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = -sig[:, :3].mean(1)
    s = sig[:, :3] + p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape((n,) * 3), p.reshape((n,) * 3), u.reshape((n + 1,) * 3 + (3,))


P_H = None
if CASE == "law":
    c = np.load(ROOT / "figdata/anim_seed_corner_3d.npz")         # the corner run without the law
    _, p0, _ = mech(down(c["phi"][-1].astype(float), NM), down(c["am1"][-1].astype(float), NM))
    pm = down(c["phi"][-1].astype(float), NM)
    P_H = 0.4 * float(np.percentile(p0[pm > 0.3], 95))
    print("p_h =", P_H, "Pa")

if not CACHE.exists():
    h = 2 / NF
    dt = 0.15 * h * h / BETA
    per = int(np.ceil(TEND / NFRAMES / dt))
    dt = TEND / (NFRAMES * per)
    phi = np.where(seed, 1.0, 0.0)
    alpha = np.ones_like(phi)
    T, P, A, V, U, G = [], [], [], [], [], []

    def rhs(phi, alpha, g):
        dphi = BETA * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
        return np.where(OPEN, dphi, 0.0), np.where(OPEN, KA * phi * g, 0.0)

    for fr in range(NFRAMES + 1):
        vm, p, u = mech(down(phi, NM), down(np.where(OPEN, alpha - 1, 0), NM))
        g = np.ones((NM,) * 3) if P_H is None else np.clip(1 - p / P_H, 0, 1)
        T.append(fr * per * dt); P.append(np.where(OPEN, phi, np.nan).astype(np.float32))
        A.append(np.where(OPEN, alpha - 1, np.nan).astype(np.float32)); V.append(vm.astype(np.float32))
        U.append(u.astype(np.float32)); G.append(g.astype(np.float32))
        print(f"{CASE} frame {fr}/{NFRAMES}  max vM {vm.max():.3e}", flush=True)
        if fr == NFRAMES:
            break
        gf = up(g, NF)
        for _ in range(per):
            k1 = rhs(phi, alpha, gf)
            k2 = rhs(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1], gf)
            k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1], gf)
            k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1], gf)
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, t=np.array(T), phi=np.array(P), am1=np.array(A), vm=np.array(V), u=np.array(U),
                        g=np.array(G), p_h=np.array(P_H if P_H else np.nan))

D = np.load(CACHE)
t, PHI, AM1, VM, U = D["t"], D["phi"], D["am1"], D["vm"], D["u"]

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
    if CASE == "two":                                   # cut the front half away (y < 0) to show both seeds
        yc = -1 + (np.arange(n) + 0.5) * 2 / n
        keep = keep & (yc[None, :, None] > 0)
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1, -1, -1)).cast_to_structured_grid()
    g.points = g.points + MAG * disp(Uf, n)
    g.cell_data["v"] = np.nan_to_num(field).ravel(order="F")
    g.cell_data["keep"] = keep.ravel(order="F").astype(float)
    return g.threshold(0.5, scalars="keep") if not keep.all() else g


bio_f = OPEN
bio_m = ~tooth_m


def vmax_of(V, phi):
    pm = down(np.nan_to_num(phi), NM)
    return float(np.percentile(V[pm > 0.3], 95))


if CASE == "law":
    c = np.load(ROOT / "figdata/anim_seed_corner_3d.npz")
    VM0, U0 = c["vm"], c["u"]
    vtop = vmax_of(VM0[-1], c["phi"][-1]) * 1e3
    PANELS = [(lambda i: (PHI[i], U[i], bio_f), "phase phi (with the law)", "Greens", (0, 1), "%.1f"),
              (lambda i: (VM0[i] * 1e3, U0[i], bio_m), "von Mises, mPa (without the law)", "viridis", (0, vtop), "%.1f"),
              (lambda i: (VM[i] * 1e3, U[i], bio_m), "von Mises, mPa (with the law)", "viridis", (0, vtop), "%.1f")]
    CAM = [(-4.0, -3.3, 2.6), (0, 0, -0.1), (0, 0, 1)]
    NOTE = f"homeostatic pressure p_h = {float(D['p_h']) * 1e3:.1f} mPa; bottom face clamped; displacement x{MAG:g}"
else:
    vtop = vmax_of(VM[-1], PHI[-1]) * 1e3
    PANELS = [(lambda i: (PHI[i], U[i], bio_f), "phase phi", "Greens", (0, 1), "%.1f"),
              (lambda i: (AM1[i] * 1e3, U[i], bio_f), "alpha - 1 (1e-3)", "Oranges",
               (0, float(np.nanmax(AM1)) * 1e3), "%.1f"),
              (lambda i: (VM[i] * 1e3, U[i], bio_m), "von Mises (mPa)", "viridis", (0, vtop), "%.1f")]
    CAM = {"two": [(-1.8, -5.0, 2.4), (0, 0, -0.2), (0, 0, 1)],
           "tooth": [(-4.6, -3.2, 2.2), (0, 0, -0.3), (0, 0, 1)]}[CASE]
    NOTE = ("tooth (grey) E = 1e4 Pa; " if CASE == "tooth" else "") + f"bottom face clamped; displacement x{MAG:g}"

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
acts = []
for j, (get, lab, cm, clim, fmt) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    a = pl.add_mesh(grid(*get(0)), scalars="v", cmap=cm, clim=clim, ambient=0.55, diffuse=0.5, specular=0.0,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                         width=0.7, height=0.07, label_font_size=24, fmt=fmt, n_labels=5,
                                         font_family="times", color="black"))
    tm = None
    if CASE == "tooth":
        tm = pl.add_mesh(grid(np.zeros((NM,) * 3), U[0], tooth_m), color="#c8c3b8", ambient=0.5, diffuse=0.6)
    pl.add_text(lab, position="upper_edge", font_size=14, font="times", color="black")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).extract_feature_edges(), color="#9a9a9a", line_width=1.0)
    pl.add_mesh(pv.Plane(center=(0, 0, -1.005), direction=(0, 0, 1), i_size=2.6, j_size=2.6), color="#b0b0b0")
    pl.camera_position = CAM
    pl.camera.zoom(0.82)
    acts.append((a, tm, get))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.subplot(0, 2)
pl.add_text(NOTE, position="lower_right", font_size=8, font="times", color="#555555")
pl.link_views()


def frame(i, azim):
    for a, tm, get in acts:
        a.mapper.SetInputData(grid(*get(i)))
        if tm is not None:
            tm.mapper.SetInputData(grid(np.zeros((NM,) * 3), U[i], tooth_m))
    title.SetText(2, f"T* = {t[i]:.2f}")
    pl.camera.azimuth = azim


out = ROOT / "assets" / f"anim_{CASE}_3d"
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
