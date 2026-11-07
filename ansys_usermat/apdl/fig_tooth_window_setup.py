"""fig_tooth_window_setup.py -- set-up figure of the cervical window on the molar crown (papers/coupling_paper/
REAL_TOOTH_WINDOW.ja.md; the runs of coupled_living_consumption.py --obstacle stl --window mixed|layered --stl-scale 100
--box 30 --only window_ab, the film stress of real_tooth_coupled.py --window). The window geometry is taken from the
solver itself (the module is imported with the run's options; crown_mesh, SDF, OBST), not re-derived.

(a) the crown (tier2b_real/_crown_real_watertight.stl, Open-Full-Jaw, Gholamalizadeh et al. 2022) with the window
    at the cervical -x flank; (b) the window in 3D: the crown (voxel body of the field grid) in grey, the start film
    of the mixed start, the octant x < 0, y < 0, z > 0 removed, the faces with their boundary conditions; (c) the
    layered start (A on the inner half, B on the outer half); (d) the section y = 0 at the start with the boundary
    conditions on the edges. Lengths in um (1 code mm = 10 um). The solver's z is vertical and is drawn up: the
    camera of the cube examples ((3.6, 2.6, 4.4) with y up) with y and z exchanged and x, y mirrored, so that the
    open (saliva) side faces the viewer. Printed width 16 cm, Times (figstyle), 400 dpi; no title.

    python ansys_usermat/apdl/fig_tooth_window_setup.py [--nf 64]   -> assets/fig_tooth_window_setup.pdf / .png
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pyvista as pv

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(size=8)
NF = int(sys.argv[sys.argv.index("--nf") + 1]) if "--nf" in sys.argv else 64
# the solver with the window's options (no point model, no run): its crown placement, SDF and tooth mask
sys.argv = [sys.argv[0], "--nf", str(NF), "--front", "1", "--v2", "--tvd", "--obstacle", "stl", "--window", "mixed",
            "--stl-scale", "100", "--box", "30", "--only", "window_ab"]
sys.path.insert(0, str(HERE))
import coupled_living_consumption as CL  # noqa: E402

BOX, H, SDF, OBST, H0 = CL.BOX, CL.H, CL.SDF, CL.OBST, CL.FILM_H0
UM = 10.0                                               # um per code mm
CROWN = CL.crown_mesh()
FILM0 = (SDF >= 0) & (SDF < H0)
INNER = FILM0 & (SDF < 0.5 * H0)
COL = dict(tooth="#c9c4b8", mixed="#7b5ea7", A="#3b6fb6", B="#c0392b", saliva="#5b9bd5", gingiva="#c78f5a", cut="#6a994e")
LIT = dict(ambient=0.45, diffuse=0.6, specular=0.0)
S = BOX / 2
CAM = [(-3.6 * S, -4.4 * S, 2.6 * S), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]   # the cube camera, z up, open side in front


def voxels(mask):
    g = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-S,) * 3)
    g.cell_data["keep"] = mask.ravel(order="F").astype(float)
    return g.threshold(0.5, scalars="keep")


class Render:
    """an off-screen render and the projection of 3D points into the cropped image (vtk display coordinates)"""

    def __init__(self, add, cam, zoom=1.0, size=(1500, 1400)):
        import vtk
        pl = pv.Plotter(off_screen=True, window_size=size)
        pl.enable_anti_aliasing("ssaa")
        pl.set_background("white")
        add(pl)
        if isinstance(cam, dict):                       # view direction, the whole scene in view
            pl.view_vector(cam["dir"], viewup=cam["up"])
            pl.reset_camera()
        else:
            pl.camera_position = cam
        pl.camera.zoom(zoom)
        img = pl.screenshot(return_img=True)
        self.proj = {}
        c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld()

        def project(p):
            c.SetValue(*p)
            px, py = c.GetComputedDoubleDisplayValue(pl.renderer)
            return px, py
        self._project = project
        rows = np.where((img < 250).any(axis=(1, 2)))[0]
        cols = np.where((img < 250).any(axis=(0, 2)))[0]
        pad = 12
        self.r0, self.c0 = max(rows[0] - pad, 0), max(cols[0] - pad, 0)
        self.img = img[self.r0:rows[-1] + pad, self.c0:cols[-1] + pad]
        self.hw = img.shape[0]
        self.pts = {}
        self.pl = pl

    def xy(self, p):
        """image coordinates (column, row) of the world point p"""
        px, py = self._project(p)
        return px - self.c0, (self.hw - 1 - py) - self.r0

    def close(self):
        self.pl.close()


# ---------------------------------------------------------------- (a) the crown with the window
def add_overview(pl):
    pl.add_mesh(CROWN, color=COL["tooth"], smooth_shading=True, ambient=0.35, diffuse=0.7, specular=0.1)
    box = pv.Box((-S, S, -S, S, -S, S))
    pl.add_mesh(box, color="#c0392b", opacity=0.9)
    pl.add_mesh(box.outline(), color="#c0392b", line_width=3)


lo, hi = CROWN.bounds[::2], CROWN.bounds[1::2]
RA = Render(add_overview, dict(dir=(-1.0, -0.55, 0.35), up=(0, 0, 1)), zoom=1.25, size=(1300, 1300))
pa_win = RA.xy((-S, -S, 0))
pa_base = RA.xy((lo[0] + 0.35 * (hi[0] - lo[0]), lo[1] + 0.2 * (hi[1] - lo[1]), lo[2]))
RA.close()


# ---------------------------------------------------------------- (b), (c) the window in 3D
def cut(mask):
    ax_ = -S + (np.arange(NF) + 0.5) * H
    X, Y, Zc = np.meshgrid(ax_, ax_, ax_, indexing="ij")
    return mask & ~((X < 0) & (Y < 0) & (Zc > 0))


def add_window(layered):
    def add(pl):
        pl.add_mesh(pv.Box((-S, S, -S, S, -S, S)).outline(), color="grey", line_width=2)
        pl.add_mesh(voxels(cut(OBST)), color=COL["tooth"], **LIT)
        if layered:
            pl.add_mesh(voxels(cut(INNER)), color=COL["A"], **LIT)
            pl.add_mesh(voxels(cut(FILM0 & ~INNER)), color=COL["B"], **LIT)
        else:
            pl.add_mesh(voxels(cut(FILM0)), color=COL["mixed"], **LIT)
        # the faces as translucent planes on the part left after the cut (saliva: x = -W/2; gingiva: z = -W/2; cuts)
        d = 0.02 * S                                   # planes just outside the box (no z-fighting with the crown)
        pl.add_mesh(pv.Plane(center=(-S - d, S / 2, 0), direction=(1, 0, 0), i_size=BOX, j_size=S).triangulate(),
                    color=COL["saliva"], opacity=0.3)
        pl.add_mesh(pv.Plane(center=(-S - d, -S / 2, -S / 2), direction=(1, 0, 0), i_size=S, j_size=S).triangulate(),
                    color=COL["saliva"], opacity=0.3)
        pl.add_mesh(pv.Plane(center=(0, 0, -S - d), direction=(0, 0, 1), i_size=BOX, j_size=BOX).triangulate(),
                    color=COL["gingiva"], opacity=0.3)
        pl.add_mesh(pv.Plane(center=(0, -S - d, -S / 2), direction=(0, 1, 0), i_size=BOX, j_size=S).triangulate(),
                    color=COL["cut"], opacity=0.22)
        pl.add_mesh(pv.Plane(center=(S / 2, -S - d, S / 2), direction=(0, 1, 0), i_size=S, j_size=S).triangulate(),
                    color=COL["cut"], opacity=0.22)
    return add


RB = Render(add_window(False), CAM, zoom=1.0)
pb = {k: RB.xy(p) for k, p in dict(saliva=(-S, 0.6 * S, -0.5 * S), gingiva=(-0.3 * S, -0.5 * S, -S), cut_y=(0.1 * S, -S, -0.5 * S),
                                   top=(-0.4 * S, 0.7 * S, S), tooth=(0.6 * S, -0.5 * S, 0.5 * S), film=(-0.12 * S, 0.6 * S, 0.3 * S),
                                   ox=(-S, -S, -S), ex=(-0.5 * S, -S, -S), ey=(-S, -0.5 * S, -S), ez=(-S, -S, -0.5 * S)).items()}
RB.close()
RC_ = Render(add_window(True), CAM, zoom=1.0)
RC_.close()


# ---------------------------------------------------------------- the figure
fig = plt.figure(figsize=(6.3, 6.0))
# (a)
ax = fig.add_axes([0.01, 0.60, 0.33, 0.38])
ax.imshow(RA.img, interpolation="none"); ax.set_axis_off()
ax.text(0.0, 1.0, "(a)", transform=ax.transAxes, va="bottom", fontsize=8)
ax.annotate("window, 300 $\\mu$m,\nat the cervical $-x$ flank", xy=pa_win, xytext=(0.28, -0.05), textcoords="axes fraction",
            fontsize=7.5, ha="center", va="center", arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.8), color="#c0392b",
            annotation_clip=False)
ax.annotate("cervical base", xy=pa_base, xytext=(0.82, -0.05), textcoords="axes fraction", fontsize=7.5, ha="center", va="center",
            arrowprops=dict(arrowstyle="-|>", color="0.3", lw=0.6), color="0.3", annotation_clip=False)
ax.text(0.5, 0.99, "molar crown (6.7 mm wide)", transform=ax.transAxes, fontsize=7.5, ha="center", va="top", color="0.3")
# (b)
ax = fig.add_axes([0.46, 0.57, 0.44, 0.40])
ax.imshow(RB.img, interpolation="none"); ax.set_axis_off()
ax.text(-0.25, 1.0, "(b)", transform=ax.transAxes, va="bottom", fontsize=8)
lab = dict(fontsize=7, ha="center", va="center", textcoords="axes fraction", annotation_clip=False)
arr = dict(arrowstyle="-|>", lw=0.7, color="0.25", shrinkB=2)
ax.annotate("saliva face $x=-W/2$:\n$c=1$, $a=1$ from $T^*=1/3$;\nfree", xy=pb["saliva"], xytext=(-0.17, 0.38),
            arrowprops=arr, color=COL["saliva"], **lab)
ax.annotate("gingiva $z=-W/2$:\nzero flux, $u_z=0$", xy=pb["gingiva"], xytext=(0.0, -0.08), arrowprops=arr, color="#9c6a35", **lab)
ax.annotate("cut faces $y=\\pm W/2$:\nzero flux, roller ($u_y=0$)", xy=pb["cut_y"], xytext=(0.62, -0.08),
            arrowprops=arr, color=COL["cut"], **lab)
ax.annotate("cut face $z=+W/2$:\nzero flux, roller ($u_z=0$)", xy=pb["top"], xytext=(0.66, 1.06),
            arrowprops=arr, color=COL["cut"], **lab)
ax.annotate("crown: rigid ($u=0$),\nno biofilm, no flux", xy=pb["tooth"], xytext=(1.08, 0.55), arrowprops=arr, color="0.3", **lab)
ax.annotate("start film, 20 $\\mu$m\n(mixed start)", xy=pb["film"], xytext=(0.15, 1.06), arrowprops=arr, color=COL["mixed"], **lab)
for k, name, col in (("ex", "x", "#d9480f"), ("ey", "y", "#2b8a3e"), ("ez", "z", "#1c3faa")):
    ax.annotate("", xy=pb[k], xytext=pb["ox"], arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0, mutation_scale=7), annotation_clip=False)
    tx = pb["ox"][0] + 1.3 * (pb[k][0] - pb["ox"][0]); ty = pb["ox"][1] + 1.3 * (pb[k][1] - pb["ox"][1])
    ax.text(tx, ty, f"${name}$", color=col, fontsize=8, ha="center", va="center")
# (c)
ax = fig.add_axes([0.01, 0.12, 0.33, 0.38])
ax.imshow(RC_.img, interpolation="none"); ax.set_axis_off()
ax.text(0.0, 1.0, "(c)", transform=ax.transAxes, va="bottom", fontsize=8)
ax.text(0.5, -0.02, "layered start: A (blue) on the inner 10 $\\mu$m,\nB (red) on the outer 10 $\\mu$m",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
# (d) section y = 0
ax = fig.add_axes([0.53, 0.12, 0.32, 0.30])
j = NF // 2
e = (-S + np.arange(NF + 1) * H) * UM
ax.pcolormesh(e, e, OBST[:, j, :].T, cmap=ListedColormap(["white", COL["tooth"]]), vmin=0, vmax=1, shading="flat", rasterized=True)
film = np.ma.masked_where(~FILM0[:, j, :], np.where(INNER[:, j, :], 1.0, 2.0))
ax.pcolormesh(e, e, film.T, cmap=ListedColormap([COL["A"], COL["B"]]), vmin=1, vmax=2, shading="flat", rasterized=True)
ax.set_aspect("equal"); ax.grid(False)
ax.set_xlim(e[0], e[-1]); ax.set_ylim(e[0], e[-1])
ax.set_xlabel(r"$x$ [$\mu$m]", labelpad=1); ax.set_ylabel(r"$z$ [$\mu$m]", labelpad=1)
ax.tick_params(labelsize=7, pad=1)
ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
for sp, col in (("left", COL["saliva"]), ("bottom", "#9c6a35"), ("top", COL["cut"])):
    ax.spines[sp].set_color(col); ax.spines[sp].set_linewidth(2.0)
ax.text(0.0, 1.10, "(d)  section $y=0$, $T^*=0$", transform=ax.transAxes, va="bottom", fontsize=8)
ax.text(-0.30, 0.5, "$c=1$, $a=1$ ($T^*\\geq 1/3$); free", transform=ax.transAxes, rotation=90, ha="center", va="center",
        fontsize=7, color=COL["saliva"])
ax.text(0.5, -0.18, "gingiva: zero flux, $u_z=0$", transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#9c6a35")
ax.text(0.5, 1.02, "cut: zero flux, $u_z=0$", transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=COL["cut"])
ax.text(85, 60, "crown\nrigid, no flux", fontsize=7, ha="center", va="center", color="0.3")
ax.text(-95, 60, "saliva", fontsize=7, ha="center", va="center", color="0.3")
ax.text(-60, -95, "film\n20 $\\mu$m", fontsize=7, ha="right", va="center", color="0.2")
ax.annotate("", xy=(-10, -95), xytext=(-55, -95), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.2"))
ax.plot([-140, -40], [-138, -138], color="k", lw=1.5); ax.text(-90, -132, "100 $\\mu$m", fontsize=7, ha="center", va="bottom")
ax.text(0.5, -0.27, "section plane $y=0$; the faces $y=\\pm W/2$: zero flux, $u_y=0$", transform=ax.transAxes,
        ha="center", va="top", fontsize=7, color=COL["cut"])
out = ROOT / "assets" / "fig_tooth_window_setup.pdf"
fig.savefig(out, dpi=400)
fig.savefig(out.with_suffix(".png"), dpi=400)
print("wrote", out, "window placement", CL.WIN_INFO)
