"""anim_real_tooth_3d.py -- the Chapter 4 seed model on a real molar surface (demonstration, not a thesis result).

Geometry: lower right first molar (FDI 46) from the CC BY 4.0 dataset of Diaz, Ferro Sanchez and Campo (Mendeley Data,
doi:10.17632/xjsx7nfhj8.1), file Teeth/STL/D46.stl, bundled unchanged as assets/geometry/mendeley_xjsx7nfhj8_D46.stl
(see assets/geometry/README.md). The tooth is used at its real size
(crown about 10 x 12 mm): it is centred on the crown and placed as oriented in the jaw (occlusal plane horizontal) in a
cube of L = 16 mm so that the crown top sits at z = +5.5 mm and the bottom face z = -8 mm cuts the roots about 4 mm
below the furcation; the bottom face is clamped, which holds the tooth.

Diffusion coefficient. The thesis value beta = 0.02 mm^2/T* belongs to the 2 mm cube, where the seed is 0.25 mm thick
and spreads by sqrt(beta T) = 0.15 mm (60 % of its thickness) in T = 1.1 T*. At real size the film here is 0.5 mm thick,
so beta is scaled by (0.5/0.25)^2 = 4 to keep that ratio: beta = 0.08 mm^2/T*, spreading 0.3 mm. The two obvious
alternatives were tried and rejected: the thesis value unchanged spreads 0.15 mm, below one cell of any affordable grid;
scaling beta with the cube edge (x 8^2 = 1.28 mm^2/T*, the same dimensionless problem as the 2 mm cube with the tooth
shrunk 1:8) spreads 1.2 mm and dilutes the 0.5 mm film to about a quarter within the run (max alpha - 1 = 3.3e-4 instead
of about 1e-3), so the film turns into a cloud. k_alpha is unchanged, so alpha - 1 and the stresses are comparable with
the 2 mm demonstrations (the small-strain mechanics has no length scale).

Model as anim_scenarios_3d.py: d phi/dt = beta lap phi + k_alpha alpha - penalty, k_alpha = 1e-3 /T*, penalty 5,
T* = 1.1, RK4 on a 96^3 cell grid (h = 0.167 mm) with zero flux at the walls and at the tooth (masked Laplacian; the
tooth is the set of cells whose centre lies inside the STL). Seed: phi = 1 in the cells within three cells (0.5 mm) of
the tooth surface, on one flank (within 60 degrees of the -x direction) and in a 2 mm band at the cervical third of the crown
(z = -3 .. -1 mm, the cemento-enamel junction is near z = -1.6). Mechanics: small-strain hex8 B-bar on 32^3
(mesh_study_seed.element_matrices), E = 10 Pa (phi^2 + 1e-3) in the film and the space, E = 1e4 Pa in the tooth without
growth, nu = 0.49, eigenstrain (alpha - 1) I, bottom face clamped, AMG-preconditioned CG (pyamg). The stress does not act
back on the growth here, so it is solved after the field run, every third frame, and interpolated in between.

    python ansys_usermat/apdl/anim_real_tooth_3d.py [--png-only]
        -> assets/anim_real_tooth_3d.mp4, .png; cache figdata/anim_real_tooth_3d.npz;
           ansys_usermat/apdl/results/2026-10-3d/real_tooth.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

PNG_ONLY = "--png-only" in sys.argv
L = 16.0                                   # cube edge, mm
NF, NM = 96, 32                            # field cells, stress elements per edge
NFRAMES, STRESS_EVERY = 66, 3
BETA0, KA, PEN, TEND = 0.02, 1e-3, 5.0, 1.1
SEED_THICK = 0.5                           # film thickness, mm
BETA = BETA0 * (SEED_THICK / 0.25) ** 2    # 0.08 mm^2/T*, see the docstring
E_BIO, E_TOOTH, FLOOR, NU = 10.0, 1e4, 1e-3, 0.49
Z_TOP = 5.5                                # crown top in the cube, mm
SEED_Z, SEED_HALF_ANGLE = (-3.0, -1.0), 60.0
STL = ROOT / "assets" / "geometry" / "mendeley_xjsx7nfhj8_D46.stl"
CACHE = ROOT / "figdata" / "anim_real_tooth_3d.npz"
OUT = ROOT / "assets" / "anim_real_tooth_3d"
RES = HERE / "results" / "2026-10-3d" / "real_tooth.json"
H = L / NF
WALL = {}


def up(a, n):
    f = n // a.shape[0]
    return np.repeat(np.repeat(np.repeat(a, f, 0), f, 1), f, 2)


def down(a, n):
    f = a.shape[0] // n
    return a.reshape(n, f, n, f, n, f).mean((1, 3, 5))


# ---------------------------------------------------------------- geometry
import pyvista as pv  # noqa: E402


def load_tooth():
    """The molar, centred on its crown in x and y, crown top at Z_TOP."""
    m = pv.read(STL).triangulate().clean()
    P = np.asarray(m.points)
    top = P[:, 2].max()
    crown = P[P[:, 2] > top - 7.0]                         # the crown is about 7 mm high
    shift = np.array([crown[:, 0].mean(), crown[:, 1].mean(), top - Z_TOP])
    m.points = P - shift
    return m


tooth_mesh = load_tooth()
ax = -L / 2 + (np.arange(NF) + 0.5) * H
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
t0 = time.time()
cen = pv.ImageData(dimensions=(NF,) * 3, spacing=(H,) * 3, origin=(ax[0],) * 3)
SDF = np.asarray(cen.compute_implicit_distance(tooth_mesh)["implicit_distance"]).reshape((NF,) * 3, order="F")
WALL["sdf"] = time.time() - t0
tooth = SDF < 0
OPEN = ~tooth
ang = np.degrees(np.arctan2(Y, -X))                          # 0 at the -x flank
seed = OPEN & (SDF < SEED_THICK) & (Z > SEED_Z[0]) & (Z < SEED_Z[1]) & (np.abs(ang) < SEED_HALF_ANGLE)
tooth_m = down(tooth.astype(float), NM) > 0.5
print(f"tooth: {tooth.sum()} cells = {tooth.sum() * H ** 3:.0f} mm^3 inside the cube; seed {seed.sum()} cells "
      f"= {seed.sum() * H ** 3:.2f} mm^3; sdf {WALL['sdf']:.1f} s", flush=True)


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
    return out / H ** 2


# ---------------------------------------------------------------- mechanics
class Mech:
    def __init__(self, n):
        import pyamg
        self.pyamg = pyamg
        self.n = n
        h = L / n
        self.Ke0, self.fe0, self.B0, self.D, self.m = S.element_matrices(h, NU, False)
        idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
        off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
        self.dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
        self.rows = np.repeat(self.dofs, 24, 1).ravel()
        self.cols = np.tile(self.dofs, (1, 24)).ravel()
        self.ndof = 3 * (n + 1) ** 3
        fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
        self.free = np.setdiff1d(np.arange(self.ndof), fix)
        xn = np.linspace(-L / 2, L / 2, n + 1)
        Xn, Yn, Zn = np.meshgrid(xn, xn, xn, indexing="ij")
        P = np.stack([Xn.ravel(), Yn.ravel(), Zn.ravel()], 1)
        B = np.zeros((self.ndof, 6))                         # rigid body modes for the AMG near-nullspace
        B[0::3, 0] = 1; B[1::3, 1] = 1; B[2::3, 2] = 1
        B[0::3, 3] = -P[:, 1]; B[1::3, 3] = P[:, 0]
        B[1::3, 4] = -P[:, 2]; B[2::3, 4] = P[:, 1]
        B[0::3, 5] = P[:, 2]; B[2::3, 5] = -P[:, 0]
        self.Bf = B[self.free]
        self.u_prev = None

    def solve(self, phi_m, am1_m):
        n = self.n
        Ee = np.where(tooth_m.ravel(), E_TOOTH, E_BIO * (np.clip(phi_m.ravel(), 0, 1) ** 2 + FLOOR))
        eps0 = np.where(tooth_m.ravel(), 0.0, am1_m.ravel())
        K = sp.csr_matrix(((self.Ke0[None] * Ee[:, None, None]).ravel(), (self.rows, self.cols)),
                          shape=(self.ndof, self.ndof))
        F = np.zeros(self.ndof)
        np.add.at(F, self.dofs.ravel(), (self.fe0[None] * (Ee * eps0)[:, None]).ravel())
        free = self.free
        Kff = K[free][:, free].tocsr()
        u = np.zeros(self.ndof)
        if np.abs(F[free]).max() > 0:
            ml = self.pyamg.smoothed_aggregation_solver(Kff, B=self.Bf, max_coarse=1000)
            res = []
            u[free] = ml.solve(F[free], x0=self.u_prev, tol=1e-8, accel="cg", maxiter=600, residuals=res)
            self.u_prev = u[free].copy()
            self.last_iters = len(res)
        else:
            self.last_iters = 0
        eps = (self.B0 @ u[self.dofs].T).T - eps0[:, None] * self.m
        sig = Ee[:, None] * (eps @ self.D.T)
        p = -sig[:, :3].mean(1)
        s = sig[:, :3] + p[:, None]
        vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
        return vm.reshape((n,) * 3), p.reshape((n,) * 3), u.reshape((n + 1,) * 3 + (3,))


# ---------------------------------------------------------------- run
if not CACHE.exists():
    dt = 0.15 * H * H / BETA
    per = int(np.ceil(TEND / NFRAMES / dt))
    dt = TEND / (NFRAMES * per)
    print(f"field: dt = {dt:.4f} T*, {per} steps per frame, {NFRAMES * per} steps", flush=True)
    phi = np.where(seed, 1.0, 0.0)
    alpha = np.ones_like(phi)
    T, P, A = [], [], []

    def rhs(phi, alpha):
        dphi = BETA * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
        return np.where(OPEN, dphi, 0.0), np.where(OPEN, KA * phi, 0.0)

    t0 = time.time()
    for fr in range(NFRAMES + 1):
        T.append(fr * per * dt)
        P.append(phi.astype(np.float16)); A.append(np.where(OPEN, alpha - 1, 0).astype(np.float16))
        if fr == NFRAMES:
            break
        for _ in range(per):
            k1 = rhs(phi, alpha)
            k2 = rhs(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1])
            k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1])
            k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1])
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if fr % 11 == 0:
            print(f"field frame {fr}/{NFRAMES}  phi sum {phi.sum() * H ** 3:.3f} mm^3  max alpha-1 "
                  f"{(alpha - 1)[OPEN].max():.2e}", flush=True)
    WALL["field"] = time.time() - t0
    P, A = np.array(P), np.array(A)

    mech = Mech(NM)
    solved = list(range(0, NFRAMES + 1, STRESS_EVERY))
    if solved[-1] != NFRAMES:
        solved.append(NFRAMES)
    VMs, PRs, ITs = {}, {}, {}
    t0 = time.time()
    for fr in solved:
        t1 = time.time()
        vm, p, _ = mech.solve(down(P[fr].astype(float), NM), down(A[fr].astype(float), NM))
        VMs[fr], PRs[fr], ITs[fr] = vm, p, mech.last_iters
        print(f"stress frame {fr}/{NFRAMES}  max vM in film {vm[~tooth_m].max():.3e} Pa  {mech.last_iters} it "
              f"{time.time() - t1:.0f} s", flush=True)
    WALL["stress"] = time.time() - t0
    V = np.zeros((NFRAMES + 1, NM, NM, NM), np.float32)
    PR = np.zeros_like(V)
    for a, b in zip(solved[:-1], solved[1:]):               # linear in time between the solved frames
        for fr in range(a, b + 1):
            w = (fr - a) / (b - a)
            V[fr] = (1 - w) * VMs[a] + w * VMs[b]
            PR[fr] = (1 - w) * PRs[a] + w * PRs[b]
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, t=np.array(T), phi=P, am1=A, vm=V, p=PR, sdf=SDF.astype(np.float32),
                        solved=np.array(solved), iters=np.array([ITs[f] for f in solved]),
                        wall_sdf=WALL["sdf"], wall_field=WALL["field"], wall_stress=WALL["stress"])

D = np.load(CACHE)
t, PHI, AM1, VM = D["t"], D["phi"], D["am1"], D["vm"]
nf = len(t)

# ---------------------------------------------------------------- numbers
film_m = lambda i: down(PHI[i].astype(float), NM) > 0.3                      # noqa: E731
vm_end = VM[-1]
sel = film_m(-1) & ~tooth_m
vmax = float(vm_end[sel].max())
loc = np.unravel_index(np.argmax(np.where(sel, vm_end, -1)), vm_end.shape)
axm = -L / 2 + (np.arange(NM) + 0.5) * L / NM
loc_mm = [float(axm[k]) for k in loc]
vtop = float(np.percentile(vm_end[sel], 95))
phi_end = PHI[-1].astype(float)
results = dict(
    geometry=dict(source="Diaz, Ferro Sanchez, Campo (2024), Mendeley Data v1, doi:10.17632/xjsx7nfhj8.1, "
                         "Teeth/STL/D46.stl (lower right first molar), CC BY 4.0",
                  file=str(STL.relative_to(ROOT)), crown_top_z_mm=Z_TOP,
                  tooth_cells_in_cube=int(tooth.sum()), tooth_volume_in_cube_mm3=float(tooth.sum() * H ** 3)),
    settings=dict(L_mm=L, field_cells=NF, stress_elements=NM, h_field_mm=H, h_stress_mm=L / NM,
                  beta_mm2_per_Tstar=BETA, beta_thesis_2mm_cube=BETA0, beta_scale=(SEED_THICK / 0.25) ** 2, k_alpha=KA,
                  penalty=PEN, T_end=TEND, frames=NFRAMES, stress_every=STRESS_EVERY,
                  E_bio_Pa=E_BIO, E_tooth_Pa=E_TOOTH, floor=FLOOR, nu=NU,
                  seed=dict(thickness_mm=SEED_THICK, cells_thick=SEED_THICK / H, z_mm=list(SEED_Z),
                            half_angle_deg=SEED_HALF_ANGLE, flank="-x")),
    seed=dict(cells=int(seed.sum()), volume_mm3=float(seed.sum() * H ** 3),
              phi_mean_end=float(phi_end[seed].mean()), am1_mean_end=float(AM1[-1].astype(float)[seed].mean()),
              am1_max_end=float(AM1[-1].astype(float)[OPEN].max())),
    film=dict(phi_volume_mm3_start=float(PHI[0].astype(float).sum() * H ** 3),
              phi_volume_mm3_end=float(phi_end.sum() * H ** 3),
              cells_phi_gt_0p5_end=int((phi_end > 0.5).sum()), cells_phi_gt_0p05_end=int((phi_end > 0.05).sum())),
    von_mises=dict(max_Pa=vmax, max_location_mm=loc_mm, p95_Pa=vtop,
                   film_elements=int(sel.sum()), max_pressure_Pa=float(D["p"][-1][sel].max()),
                   amg_iterations=[int(k) for k in D["iters"]]),
    wall_s=dict(sdf=float(D["wall_sdf"]), field=float(D["wall_field"]), stress=float(D["wall_stress"])),
)
print(json.dumps(results, indent=1))

# ---------------------------------------------------------------- render
pv.OFF_SCREEN = True
pv.global_theme.font.family = "times"

# von Mises from the stress mesh onto the field cells: tooth elements blanked, nearest film value, then linear
def vm_fine(i):
    v = np.where(tooth_m, np.nan, VM[i])
    bad = np.isnan(v)
    if bad.any():
        ind = distance_transform_edt(bad, return_distances=False, return_indices=True)
        v = v[tuple(ind)]
    f = RegularGridInterpolator((axm,) * 3, v, bounds_error=False, fill_value=None)
    return f(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)


def film_surface(field, keep):
    """Envelope of the cells with phi > 0.05 (smoothed), coloured by the maximum of the field through the film:
    at each surface point the field is sampled along the inward normal over DEPTH mm and the maximum is taken, so the
    colour shows the film under the surface and not only the thin diffuse front at the envelope."""
    g = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-L / 2,) * 3)
    g.cell_data["keep"] = keep.ravel(order="F").astype(float)
    surf = g.threshold(0.5, scalars="keep").extract_surface(algorithm="dataset_surface").extract_largest()
    if surf.n_points == 0:
        surf.point_data["v"] = np.zeros(0)
        return surf
    surf = surf.smooth_taubin(n_iter=60, pass_band=0.05, feature_angle=150.0)
    surf = surf.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
    pts = np.asarray(surf.points)
    nrm = np.asarray(surf["Normals"])
    f = RegularGridInterpolator((ax,) * 3, np.nan_to_num(field), bounds_error=False, fill_value=0.0)
    depths = np.linspace(0.0, DEPTH, 10)
    vals = np.stack([f(pts - d * nrm) for d in depths], 0).max(0)
    surf.point_data["v"] = vals
    return surf


DEPTH = 1.5                                              # mm, sampling depth of the film colouring
# the tooth surface drawn half a cell inside its voxel set, so that it does not poke through the film
tooth_draw = tooth_mesh.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
tooth_draw.points = np.asarray(tooth_draw.points) - 0.5 * H * np.asarray(tooth_draw["Normals"])
tooth_draw = tooth_draw.clip_box([-L / 2, L / 2, -L / 2, L / 2, -L / 2, L / 2], invert=False)

KEEP = lambda i: (PHI[i].astype(float) > 0.05) & OPEN       # noqa: E731
PANELS = [(lambda i: (PHI[i].astype(float), KEEP(i)), "phase phi", "Greens", (0, 1), "%.1f"),
          (lambda i: (AM1[i].astype(float) * 1e3, KEEP(i)), "alpha - 1 (1e-3)", "Oranges",
           (0, float(AM1[-1].astype(float).max()) * 1e3), "%.1f"),
          (lambda i: (vm_fine(i) * 1e3, KEEP(i)), "von Mises (mPa), film", "viridis", (0, vtop * 1e3), "%.1f")]
CAM = [(-3.6 * L / 2, -1.9 * L / 2, 1.3 * L / 2), (0, 0, -0.7), (0, 0, 1)]
NOTE = (f"molar FDI 46 (grey, E = 1e4 Pa) at real size, cube {L:g} mm; film 0.5 mm, beta = {BETA:.2f} mm^2/T* "
        f"(4 x the 2 mm cube value, see the script); bottom face clamped; undeformed; "
        f"film surface = phi > 0.05, colour = maximum of the field through the film")

pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
acts = []
for j, (get, lab, cm, clim, fmt) in enumerate(PANELS):
    pl.subplot(0, j)
    pl.set_background("white")
    a = pl.add_mesh(film_surface(*get(0)), scalars="v", cmap=cm, clim=clim, ambient=0.45, diffuse=0.6,
                    specular=0.1, smooth_shading=True,
                    scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                         width=0.7, height=0.07, label_font_size=24, fmt=fmt, n_labels=5,
                                         font_family="times", color="black"))
    pl.add_mesh(tooth_draw, color="#d8d3c8", ambient=0.35, diffuse=0.7, specular=0.15, smooth_shading=True)
    pl.add_text(lab, position="upper_edge", font_size=14, font="times", color="black")
    pl.add_mesh(pv.Box((-L / 2, L / 2) * 3).extract_feature_edges(), color="#9a9a9a", line_width=1.0)
    pl.add_mesh(pv.Plane(center=(0, 0, -L / 2 - 0.01), direction=(0, 0, 1), i_size=1.3 * L, j_size=1.3 * L),
                color="#b8b8b8")
    pl.camera_position = CAM
    pl.camera.zoom(1.0)
    acts.append((a, get))
pl.subplot(0, 0)
title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
pl.subplot(0, 2)
pl.add_text(NOTE, position="lower_right", font_size=8, font="times", color="#555555")
pl.link_views()


def frame(i, azim):
    for a, get in acts:
        a.mapper.SetInputData(film_surface(*get(i)))
    title.SetText(2, f"T* = {t[i]:.2f}")
    pl.camera.azimuth = azim


OUT.parent.mkdir(exist_ok=True)
cam0 = pl.camera.azimuth
t0 = time.time()
if PNG_ONLY:
    frame(nf - 1, cam0)
    pl.screenshot(str(OUT.with_suffix(".png")))
else:
    pl.open_movie(str(OUT.with_suffix(".mp4")), framerate=16, quality=9)
    for i in range(nf):
        frame(i, cam0 + 30.0 * i / (nf - 1) - 15.0)
        pl.write_frame()
    for k in range(24):
        frame(nf - 1, cam0 + 15.0 + 15.0 * (k + 1) / 24)
        pl.write_frame()
    frame(nf - 1, cam0)
    pl.screenshot(str(OUT.with_suffix(".png")))
pl.close()
results["wall_s"]["render"] = time.time() - t0
RES.parent.mkdir(parents=True, exist_ok=True)
RES.write_text(json.dumps(results, indent=1) + "\n")
print("wrote", OUT.with_suffix(".mp4" if not PNG_ONLY else ".png"), "and", RES)
