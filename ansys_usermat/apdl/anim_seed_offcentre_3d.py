"""anim_seed_offcentre_3d.py -- 3D video of the Chapter 4 seed problem with the seed moved off the centre, so that the
field and the stress are no longer symmetric. Same model and values as the verification problem (wp3_reference:
beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, E = (phi^2 + 1e-3) 10 Pa, nu = 0.49, T* = 1.1), same seed volume
(32 cells of the 8^3 mesh = 0.5 mm^3), but:
  corner: seed 1 x 1 x 0.5 mm in the corner x = -1, y = -1, z = -1 (touches three faces);
  edge:   seed 1 x 1 x 0.5 mm on the edge x = -1, z = -1, centred in y (touches two faces);
and the bottom face z = -1 is clamped (the tooth surface) instead of the three corner constraints.
Not a thesis result: a demonstration with the thesis model.

Field 48^3 (RK4), stress 16^3 (hex8, B-bar), displacement drawn magnified.

    python ansys_usermat/apdl/anim_seed_offcentre_3d.py corner|edge
        -> assets/anim_seed_<case>_3d.mp4, .png;  cache figdata/anim_seed_<case>_3d.npz
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
CASE = sys.argv[1] if len(sys.argv) > 1 else "corner"
CACHE = ROOT / "figdata" / f"anim_seed_{CASE}_3d.npz"
NF, NM, NFRAMES = 48, 16, 88
MAG = 50.0

m8 = np.zeros((8, 8, 8), bool)
if CASE == "corner":
    m8[0:4, 0:4, 0:2] = True
elif CASE == "edge":
    m8[0:4, 2:6, 0:2] = True
else:
    raise SystemExit("case: corner or edge")
assert m8.sum() == 32


def mechanics(n, phi, am1):
    import mesh_study_seed as S
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, 0.49, False)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    Ee = 10.0 * (np.clip(phi.ravel(), 0, 1) ** 2 + 1e-3)
    eps0 = am1.ravel()
    rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
    ndof = 3 * (n + 1) ** 3
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()       # bottom face z = -1 clamped
    free = np.setdiff1d(np.arange(ndof), fix)
    u = np.zeros(ndof)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    s = sig[:, :3] - sig[:, :3].mean(1, keepdims=True)
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape(n, n, n), u.reshape(n + 1, n + 1, n + 1, 3)


if not CACHE.exists():
    import wp3_reference as W
    f = NF // 8
    seed = np.repeat(np.repeat(np.repeat(m8, f, 0), f, 1), f, 2)
    h = W.L / NF
    dt = 0.15 * h * h / W.BETA
    per = int(np.ceil(W.T_END / NFRAMES / dt))
    dt = W.T_END / (NFRAMES * per)
    phi = seed.astype(float)
    alpha = np.ones_like(phi)
    lap = np.empty_like(phi)
    T, P, A, V, U = [], [], [], [], []
    for fr in range(NFRAMES + 1):
        if fr:
            for _ in range(per):
                k1p, k1a = W.rhs(phi, alpha, h, lap)
                k2p, k2a = W.rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a, h, lap)
                k3p, k3a = W.rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a, h, lap)
                k4p, k4a = W.rhs(phi + dt * k3p, alpha + dt * k3a, h, lap)
                phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
                alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        vm, u = mechanics(NM, W.block_mean(phi, NM), W.block_mean(alpha - 1, NM))
        T.append(fr * per * dt); P.append(phi.astype(np.float32)); A.append((alpha - 1).astype(np.float32))
        V.append(vm.astype(np.float32)); U.append(u.astype(np.float32))
        print(f"{CASE} frame {fr}/{NFRAMES}", flush=True)
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, t=np.array(T), phi=np.array(P), am1=np.array(A), vm=np.array(V), u=np.array(U))

D = np.load(CACHE)
t, PHI, AM1, VM, U = D["t"], D["phi"], D["am1"], D["vm"], D["u"]

import pyvista as pv  # noqa: E402

pv.OFF_SCREEN = True
pv.global_theme.font.family = "times"
pv.global_theme.font.color = "black"


def coarse_u(i, n):
    """Displacement at the nodes of an n^3 grid (n a multiple of NM): trilinear from the 16^3 nodes."""
    from scipy.interpolate import RegularGridInterpolator
    x16 = np.linspace(-1, 1, NM + 1)
    xn = np.linspace(-1, 1, n + 1)
    X, Y, Z = np.meshgrid(xn, xn, xn, indexing="ij")
    pts = np.stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")], 1)
    return np.stack([RegularGridInterpolator((x16,) * 3, U[i][..., c])(pts) for c in range(3)], 1)


def grid(n, field, i):
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1, -1, -1)).cast_to_structured_grid()
    g.points = g.points + MAG * coarse_u(i, n)
    g.cell_data["v"] = field.ravel(order="F")
    return g


# colour limit of the stress: 95th percentile over the biofilm elements at the end (the clamped edge under the seed
# is singular and would wash out the rest); higher values are drawn in the top colour
_pm = PHI[-1].reshape(NM, NF // NM, NM, NF // NM, NM, NF // NM).mean((1, 3, 5))
VMAX = float(np.percentile(VM[-1][_pm > 0.3], 95))
PANELS = [(PHI, NF, "phase phi", "Greens", (0, 1), 1.0, "%.1f"),
          (AM1, NF, "alpha - 1 (1e-3)", "Oranges", (0, float(AM1.max()) * 1e3), 1e3, "%.1f"),
          (VM, NM, "von Mises (mPa)", "viridis", (0, VMAX * 1e3), 1e3, "%.1f")]
CAM = {"corner": [(-4.0, -3.3, 2.6), (0, 0, -0.1), (0, 0, 1)],
       "edge": [(-5.0, -1.6, 2.4), (0, 0, -0.1), (0, 0, 1)]}[CASE]   # looking at the seed side

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 900), border=False)
acts = []
for j, (F, n, lab, cm, clim, sc, fmt) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    g = grid(n, F[0] * sc, 0)
    a = pl.add_mesh(g, scalars="v", cmap=cm, clim=clim, ambient=0.55, diffuse=0.5, specular=0.0,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                         width=0.7, height=0.07, label_font_size=24, fmt=fmt, n_labels=5,
                                         font_family="times", color="black"))
    pl.add_text(lab, position="upper_edge", font_size=15, font="times", color="black")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).extract_feature_edges(), color="#9a9a9a", line_width=1.0)
    pl.add_mesh(pv.Plane((0, 0, -1.005), (0, 0, 1), 2.6, 2.6), color="#b0b0b0")      # the clamped base
    pl.camera_position = CAM
    pl.camera.zoom(0.82)
    acts.append((a, F, n, sc))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.subplot(0, 2)
pl.add_text(f"bottom face clamped; displacement x{MAG:g}", position="lower_right", font_size=9, font="times",
            color="#555555")
pl.link_views()


def frame(i, azim):
    for a, F, n, sc in acts:
        a.mapper.SetInputData(grid(n, F[i] * sc, i))
    title.SetText(2, f"T* = {t[i]:.2f}")
    pl.camera.azimuth = azim


out = ROOT / "assets" / f"anim_seed_{CASE}_3d"
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
print("wrote", out.with_suffix(".mp4"), "max vM Pa", float(VM[-1].max()))
