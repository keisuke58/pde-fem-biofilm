"""anim_wp3_seed_3d.py -- 3D video of the verification problem of Chapter 4 (seed growth in the 2 mm cube): the cube
with one octant cut away, coloured by phi, alpha - 1 and von Mises, while the camera turns slowly.

Same data as anim_wp3_seed.py (reference field on 48^3 from wp3_reference, stress on 16^3 from
wp3_stress_reference.solve_fields); the fields are cached in figdata/anim_wp3_seed_3d.npz.

    python ansys_usermat/apdl/anim_wp3_seed_3d.py   -> assets/anim_wp3_seed_3d.mp4, .gif, .png (last frame)
Needs pyvista and imageio-ffmpeg (pip install pyvista imageio-ffmpeg); renders off screen.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
CACHE = ROOT / "figdata" / "anim_wp3_seed_3d.npz"
NF, NM, NFRAMES = 48, 16, 88

if not CACHE.exists():
    import wp3_reference as W
    import wp3_stress_reference as S

    cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    h = W.L / NF
    dt = 0.15 * h * h / W.BETA
    per = int(np.ceil(W.T_END / NFRAMES / dt))
    dt = W.T_END / (NFRAMES * per)
    phi = W.seed_mask(NF, cells8).astype(float)
    alpha = np.ones_like(phi)
    lap = np.empty_like(phi)
    T, P, A, V = [], [], [], []
    for f in range(NFRAMES + 1):
        if f:
            for _ in range(per):
                k1p, k1a = W.rhs(phi, alpha, h, lap)
                k2p, k2a = W.rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a, h, lap)
                k3p, k3a = W.rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a, h, lap)
                k4p, k4a = W.rhs(phi + dt * k3p, alpha + dt * k3a, h, lap)
                phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
                alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        vm, _ = S.solve_fields(NM, W.block_mean(phi, NM), W.block_mean(alpha - 1, NM))
        T.append(f * per * dt); P.append(phi.astype(np.float32)); A.append((alpha - 1).astype(np.float32))
        V.append(vm.astype(np.float32))
        print(f"frame {f}/{NFRAMES}", flush=True)
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, t=np.array(T), phi=np.array(P), am1=np.array(A), vm=np.array(V))

D = np.load(CACHE)
t, PHI, AM1, VM = D["t"], D["phi"], D["am1"], D["vm"]

import pyvista as pv  # noqa: E402

pv.OFF_SCREEN = True
pv.global_theme.font.family = "times"
pv.global_theme.font.color = "black"


def grid(n, field):
    g = pv.ImageData(dimensions=(n + 1,) * 3, spacing=(2.0 / n,) * 3, origin=(-1, -1, -1))
    g.cell_data["v"] = field.ravel(order="F")
    return g


CUT = (0, 1.01, 0, 1.01, 0, 1.01)                   # octant removed, facing the camera
PANELS = [(PHI, NF, "phase phi", "Greens", (0, 1), 1.0, "%.1f"),
          (AM1, NF, "alpha - 1 (1e-3)", "Oranges", (0, float(AM1.max()) * 1e3), 1e3, "%.1f"),
          (VM, NM, "von Mises (mPa)", "viridis", (0, float(VM.max()) * 1e3), 1e3, "%.2f")]

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 900), border=False)
actors = []
for j, (F, n, lab, cm, clim, sc, fmt) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    g = grid(n, F[0] * sc)
    cut = g.clip_box(CUT, invert=True)
    a = pl.add_mesh(cut, scalars="v", cmap=cm, clim=clim, show_edges=False, ambient=0.55, diffuse=0.5, specular=0.0,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03, width=0.7,
                                         height=0.07, title_font_size=30, label_font_size=24, fmt=fmt, n_labels=5,
                                         font_family="times", color="black"))
    pl.add_text(lab, position="upper_edge", font_size=15, font="times", color="black")
    pl.add_mesh(pv.Box((-1, 1, -1, 1, -1, 1)).extract_feature_edges(), color="#595959", line_width=1.5)
    pl.camera_position = [(4.2, 3.4, 3.0), (0, 0, -0.05), (0, 0, 1)]
    pl.camera.zoom(0.82)
    actors.append((g, a, F, sc))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.link_views()


def frame(i, azim):
    for g, a, F, sc in actors:
        g.cell_data["v"] = (F[i] * sc).ravel(order="F")
        a.mapper.SetInputData(g.clip_box(CUT, invert=True))
    title.SetText(2, f"T* = {t[i]:.2f}")
    pl.camera.azimuth = azim


out = ROOT / "assets" / "anim_wp3_seed_3d"
cam0 = pl.camera.azimuth
pl.open_movie(str(out.with_suffix(".mp4")), framerate=16, quality=9)
nf = len(t)
for i in range(nf):
    frame(i, cam0 + 40.0 * i / (nf - 1) - 20.0)
    pl.write_frame()
for k in range(24):                                  # hold the last state and keep turning
    frame(nf - 1, cam0 + 20.0 + 15.0 * (k + 1) / 24)
    pl.write_frame()
pl.screenshot(str(out.with_suffix(".png")))
pl.close()
print("wrote", out.with_suffix(".mp4"))
