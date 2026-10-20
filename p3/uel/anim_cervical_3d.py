"""anim_cervical_3d.py -- 3D video (synthetic: no experimental data; the four-species case of Fritsch et al. 2025, see synthetic_eco) of the tooth-like sector of the cervical study (P3 (e)): the biofilm layer around a
crown-shaped collar (bulge 1), 4 x 12 x 16, inner face bonded, bottom u_z = 0, outer face free, nutrient from the
outer face (saliva), seed of 0.125 mm; without the homeostatic-pressure law and with p_h/E = 0.05. The same runs as
cervical_study.py (mi.run with these options), now with the state stored at every increment. The 90 degree sector is
mirrored to the full ring for the picture (the problem is symmetric about both cut planes).

    python p3/uel/anim_cervical_3d.py            # runs both cases (cached in figdata/anim_cervical_<p>.npz), renders
    -> assets/anim_cervical_3d.mp4, .png
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCR = Path(os.environ.get("UQ_SCRATCH", "/tmp/anim_cervical"))
NR = 4


def cache(p):
    return ROOT / "figdata" / f"anim_cervical_f4_{p:g}.npz"


def synthetic_eco(path):
    """point model without experimental data: the four-species synthetic case of Fritsch et al. 2025 (Case II, M3:
    theta* = the "true mean" of their Fig. 16, c* = 25, alpha* = 0, phi_init 0.02 each; JAXFEM/fritsch2025_cases.py),
    put into the five-species slots of the element with species 5 absent (phi_init 0, a55 = 1, no interactions);
    same file format as make_column_inp.eco_file"""
    sys.path.insert(0, str(ROOT / "JAXFEM"))
    from fritsch2025_cases import CASES
    c = CASES["fritsch_4sp_M3"]
    A, b = c["A"], c["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    th[5], th[6], th[7], th[8], th[9] = A[2][2], A[2][3], A[3][3], b[2], b[3]
    th[10], th[11], th[12], th[13] = A[0][2], A[0][3], A[1][2], A[1][3]
    th[14], th[15] = 1.0, 1.0
    lines = ["5", f"{c['c_star']} 0.0", " ".join(["1.0"] * 5), " ".join(f"{v:.17g}" for v in th),
             "phi_init " + " ".join(f"{v:.17g}" for v in c["phi_init"] + [0.0]), "newton 60 1e-26"]
    path.write_text("\n".join(lines) + "\n")


if len(sys.argv) > 1 and sys.argv[1] == "--run":
    p_ref = float(sys.argv[2])
    sys.path.insert(0, str(HERE))
    import emulate_abaqus as ea
    import make_implant_uel_inp as mi
    from emulate_abaqus import NDN, NGP, NSV
    snaps = []
    orig = ea.increment

    def rec(U0, U, fixed, nodes, elems, SV, props, dt, kinc, ttot=0.0, **kw):
        Un, SVn, it = orig(U0, U, fixed, nodes, elems, SV, props, dt, kinc, ttot=ttot, **kw)
        xi = np.array([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv in SVn])
        st = np.array([sv[NSV * NGP:].reshape(NGP, 2).mean(0) for sv in SVn])
        snaps.append((ttot + dt, Un.reshape(-1, NDN)[:, :4].copy(), xi, st[:, 0], st[:, 1]))
        return Un, SVn, it

    ea.increment = rec
    nodes, elems, layer, ijk = mi.sector(NR, 3 * NR, 4 * NR, 1.0)
    mi.run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=NR // 2, p_ref=p_ref, outer="free",
           bottom="uz", nut="outer")
    snaps = snaps[1:]                                     # the first call is INIT
    cache(p_ref).parent.mkdir(exist_ok=True)
    np.savez_compressed(cache(p_ref), t=np.array([s[0] for s in snaps]), u=np.array([s[1] for s in snaps]),
                        xi=np.array([s[2] for s in snaps]), vm=np.array([s[3] for s in snaps]),
                        p=np.array([s[4] for s in snaps]), nodes=nodes, elems=elems, layer=layer)
    print("stored", len(snaps), "increments")
    sys.exit()

# ---- run both cases in subprocesses (the point model is read from BIOFILM_ECO_CASE once per process)
SCR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
eco = SCR / "uel_eco_synthetic.txt"
synthetic_eco(eco)
procs = [subprocess.Popen([sys.executable, __file__, "--run", str(p)], env=dict(os.environ, BIOFILM_ECO_CASE=str(eco)))
         for p in (0.0, 0.05) if not cache(p).exists()]
for pr in procs:
    assert pr.wait() == 0

A, B = np.load(cache(0.0)), np.load(cache(0.05))
nodes, elems = A["nodes"], A["elems"]
# common frames: the increments of the run without the law (fixed 60), the other run by nearest time
tA, tB = A["t"], B["t"]
iB = [int(np.abs(tB - t).argmin()) for t in tA]

import pyvista as pv  # noqa: E402

pv.OFF_SCREEN = True
pv.global_theme.font.family = "times"
MAG = max(1.0, float(round(0.1 / max(np.abs(A['u'][..., :3]).max(), np.abs(B['u'][..., :3]).max()), 1)))   # magnified only when small
MIRRORS = [(1, 1), (-1, 1), (-1, -1)]                    # three quarters of the ring; the front one cut away


def ring(Un, cell, kind="cell"):
    blocks = []
    for sx, sy in MIRRORS:
        pts = nodes + MAG * Un[:, :3]
        pts = pts * np.array([sx, sy, 1.0])
        el = elems if sx * sy > 0 else elems[:, [1, 0, 3, 2, 5, 4, 7, 6]]
        g = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: el}, pts)
        if kind == "cell":
            g.cell_data["v"] = cell
        else:
            g.point_data["v"] = cell
        blocks.append(g)
    return blocks[0].merge(blocks[1:])


def core(Un):
    """the collar under the biofilm (rigid, grey): the inner surface of the layer, closed to the axis"""
    nr = NR
    inner = nodes[::nr + 1]
    pts = []
    for sx, sy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
        pts.append(inner * np.array([sx, sy, 1.0]))
    return pv.PolyData(np.vstack(pts)).delaunay_3d().extract_surface()


vtop = float(np.percentile(np.concatenate([A["vm"][-1], B["vm"][-1]]), 98))
PANELS = [(lambda i: (A["u"][i], A["u"][i][:, 3], "point"), "phase phi (without the law)", "Greens", (0, 1)),
          (lambda i: (A["u"][i], A["vm"][i], "cell"), "von Mises / E (without the law)", "viridis", (0, vtop)),
          (lambda i: (B["u"][iB[i]], B["vm"][iB[i]], "cell"), "von Mises / E (with the law, p_h/E = 0.05)", "viridis",
           (0, vtop))]
CAM = [(12.0, -10.0, 4.5), (0, 0, 0.9), (0, 0, 1)]
COLLAR = core(None)

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
acts = []
for j, (get, lab, cm, clim) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    Un, v, kind = get(0)
    a = pl.add_mesh(ring(Un, v, kind), scalars="v", cmap=cm, clim=clim, ambient=0.5, diffuse=0.55, specular=0.0,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                         width=0.7, height=0.07, label_font_size=24, fmt="%.2f" if j == 0 else "%.3f",
                                         n_labels=5, font_family="times", color="black"))
    pl.add_mesh(COLLAR, color="#c8c3b8", ambient=0.45, diffuse=0.6)
    pl.add_text(lab, position="upper_edge", font_size=13, font="times", color="black")
    pl.camera_position = CAM
    pl.camera.zoom(1.0)
    acts.append((a, get))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.subplot(0, 2)
pl.add_text(f"four-species synthetic case of Fritsch et al. 2025 (no experimental data); collar (grey) rigid; nutrient from the outside; displacement x{MAG:g}", position="lower_right",
            font_size=8, font="times", color="#555555")
pl.link_views()


def frame(i, azim):
    for a, get in acts:
        a.mapper.SetInputData(ring(*get(i)))
    title.SetText(2, f"t = {tA[i]:.3f}")
    pl.camera.azimuth = azim


out = ROOT / "assets" / "anim_cervical_f4_3d"
cam0 = pl.camera.azimuth
pl.open_movie(str(out.with_suffix(".mp4")), framerate=12, quality=9)
nf = len(tA)
for i in range(nf):
    frame(i, cam0 + 40.0 * i / (nf - 1) - 20.0)
    pl.write_frame()
for k in range(24):
    frame(nf - 1, cam0 + 20.0 + 20.0 * (k + 1) / 24)
    pl.write_frame()
pl.screenshot(str(out.with_suffix(".png")))
pl.close()
print("wrote", out.with_suffix(".mp4"), "frames", nf)
