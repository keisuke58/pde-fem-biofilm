"""real_tooth_model.py -- the seed model on the real molar with a nutrient field and the homeostatic-pressure gate,
each on or off (demonstration, not a thesis result).

Geometry, seed, field grid, stress mesh and rendering are those of anim_real_tooth_3d.py (copied, that script is not
changed): molar FDI 46 at real size in a 16 mm cube, 96^3 field cells, beta = 0.08 mm^2/T* (scaled, see there),
k_alpha = 1e-3 /T*, penalty 5, T = 1.1 T*, 0.5 mm film on the -x flank at z = -3 .. -1 mm, 32^3 hex8 B-bar stress mesh
with pyamg, E = 10 Pa (phi^2 + 1e-3), tooth 1e4 Pa without growth, nu = 0.49, bottom face clamped.

Added here, on top of d phi/dt = beta lap phi + k_alpha alpha - penalty (Eq. 34 of Klempt et al. 2024 without the
front term):

  nutrient   0 = d lap c - g phi c  (Eq. 35, quasi-static, consumption g phi c), c = 1 on the five outer faces of
             the box except the clamped bottom (saliva around the crown), zero flux at the bottom and at the tooth.
             Only g/d enters: 10 /mm^2 (ten times the Table 2 ratio converted as in anim_nutrient_3d.py, penetration
             depth 0.32 mm, so a gradient shows across the 0.5 mm film); the Table 2 ratio 1 /mm^2 is run for the
             numbers only. Solved once per RK4 step (sparse 7-point matrix, pyamg + CG, warm start).
  gate       d alpha/dt = k_alpha phi c g,  g = min(1, max(0, 1 - p/p_h)),  p = -tr(sigma)/3 from the stress at the
             start of the frame, held over the frame's RK4 steps (as the law case of anim_scenarios_3d.py). With the
             gate the stress is solved every frame. p_h = 0.5 p_0, p_0 = the mean pressure in the film (phi > 0.5 on
             the stress mesh) at the end of the run without nutrient and gate (the rule of p3/uel/README.md,
             "Threshold p_h = p_0").
  The factor c in d alpha/dt is an assumption (Eq. 36 as printed has no c), as in anim_nutrient_3d.py.

    python ansys_usermat/apdl/real_tooth_model.py base|nutrient|gate|nutrient_gate|nutrient_t2 [--no-video] [--png-only]
        -> figdata/real_tooth_<tag>.npz (base: the cache of anim_real_tooth_3d.py), assets/anim_real_tooth_<tag>_3d.mp4/.png
    python ansys_usermat/apdl/real_tooth_model.py --summary
        -> ansys_usermat/apdl/results/2026-10-3d/real_tooth_physics.json, assets/fig_real_tooth_physics.pdf/.png
"""
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")        # threaded BLAS slows the long-vector operations of the solvers
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402
from scipy.ndimage import distance_transform_edt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

L = 16.0                                   # cube edge, mm
NF, NM = 96, 32                            # field cells, stress elements per edge
NFRAMES, STRESS_EVERY = 66, 3
BETA0, KA, PEN, TEND = 0.02, 1e-3, 5.0, 1.1
SEED_THICK = 0.5                           # film thickness, mm
BETA = BETA0 * (SEED_THICK / 0.25) ** 2    # 0.08 mm^2/T*, see anim_real_tooth_3d.py
E_BIO, E_TOOTH, FLOOR, NU = 10.0, 1e4, 1e-3, 0.49
Z_TOP = 5.5                                # crown top in the cube, mm
SEED_Z, SEED_HALF_ANGLE = (-3.0, -1.0), 60.0
Z_CEJ, BAND_HALF = -1.6, 0.5               # cemento-enamel junction level and the half width of the cervical band, mm
G_TABLE2, D_TABLE2, L_PAPER_UM, L_CUBE_MM = 1e8, 1e10, 20.0, 2.0          # Klempt 2024 Table 2, 20 um cube
Q_TABLE2 = G_TABLE2 / (D_TABLE2 * (L_CUBE_MM * 1e3 / L_PAPER_UM) ** 2 * 1e-6)   # g/d = 1 /mm^2 (converted like beta)
P_H_FRACTION = 0.5                         # p_h = 0.5 p_0
STL = ROOT / "assets" / "geometry" / "mendeley_xjsx7nfhj8_D46.stl"
BASE_CACHE = ROOT / "figdata" / "anim_real_tooth_3d.npz"
RESDIR = HERE / "results" / "2026-10-3d"
H = L / NF
RUNS = {                                   # tag: (nutrient, gate, g/d in 1/mm^2)
    "base": (False, False, None),
    "nutrient": (True, False, 10 * Q_TABLE2),
    "gate": (False, True, None),
    "nutrient_gate": (True, True, 10 * Q_TABLE2),
    "nutrient_t2": (True, False, Q_TABLE2),
}


def cache_path(tag):
    return BASE_CACHE if tag == "base" else ROOT / "figdata" / f"real_tooth_{tag}.npz"


def up(a, n):
    f = n // a.shape[0]
    return np.repeat(np.repeat(np.repeat(a, f, 0), f, 1), f, 2)


def down(a, n):
    f = a.shape[0] // n
    return a.reshape(n, f, n, f, n, f).mean((1, 3, 5))


# ---------------------------------------------------------------- geometry (as anim_real_tooth_3d.py)
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
axm = -L / 2 + (np.arange(NM) + 0.5) * L / NM
Zm = np.meshgrid(axm, axm, axm, indexing="ij")[2]
if BASE_CACHE.exists():
    SDF = np.load(BASE_CACHE)["sdf"].astype(float)
else:
    cen = pv.ImageData(dimensions=(NF,) * 3, spacing=(H,) * 3, origin=(ax[0],) * 3)
    SDF = np.asarray(cen.compute_implicit_distance(tooth_mesh)["implicit_distance"]).reshape((NF,) * 3, order="F")
tooth = SDF < 0
OPEN = ~tooth
ang = np.degrees(np.arctan2(Y, -X))                          # 0 at the -x flank
seed = OPEN & (SDF < SEED_THICK) & (Z > SEED_Z[0]) & (Z < SEED_Z[1]) & (np.abs(ang) < SEED_HALF_ANGLE)
tooth_m = down(tooth.astype(float), NM) > 0.5
BAND_F = (Z > Z_CEJ - BAND_HALF) & (Z < Z_CEJ + BAND_HALF)   # cervical band on the field grid
BAND_M = (Zm > Z_CEJ - BAND_HALF) & (Zm < Z_CEJ + BAND_HALF)  # and on the stress mesh (centres z = -1.75, -1.25)


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


# ---------------------------------------------------------------- nutrient
class Nutrient:
    """Quasi-static c: -lap c + q phi c = 0 in the open cells, c = 1 (ghost cell) on the five outer faces except the
    bottom, zero flux at the bottom and at the tooth; the tooth cells are decoupled and carry c = 0."""

    def __init__(self, q):
        import pyamg
        self.pyamg = pyamg
        self.q = q
        n3 = NF ** 3
        idx = np.arange(n3).reshape((NF,) * 3)
        rows, cols, vals = [], [], []
        diag = np.zeros((NF,) * 3)
        for axn in range(3):
            for sgn in (1, -1):
                on = np.roll(OPEN, -sgn, axn)
                jn = np.roll(idx, -sgn, axn)
                edge = [slice(None)] * 3
                edge[axn] = -1 if sgn == 1 else 0
                ok = OPEN & on
                ok[tuple(edge)] = False
                diag += ok
                rows.append(idx[ok]); cols.append(jn[ok]); vals.append(-np.ones(ok.sum()))
                if not (axn == 2 and sgn == -1):               # Dirichlet ghost on the outer face (not the bottom)
                    face = np.zeros((NF,) * 3, bool)
                    face[tuple(edge)] = True
                    diag += 2.0 * (face & OPEN)
        self.b = np.zeros((NF,) * 3)
        for axn in range(3):
            for sgn in (1, -1):
                if axn == 2 and sgn == -1:
                    continue
                edge = [slice(None)] * 3
                edge[axn] = -1 if sgn == 1 else 0
                self.b[tuple(edge)] += 2.0 / H ** 2 * OPEN[tuple(edge)]
        rows.append(idx[tooth]); cols.append(idx[tooth]); vals.append(np.ones(tooth.sum()) * H ** 2)   # identity rows
        self.A0 = sp.csr_matrix((np.concatenate(vals) / H ** 2, (np.concatenate(rows), np.concatenate(cols))),
                                shape=(n3, n3))
        self.diag0 = (diag / H ** 2).ravel()
        self.c = np.where(OPEN, 1.0, 0.0)
        self.it = 0

    def solve(self, phi):
        qphi = np.where(OPEN, self.q * np.clip(phi, 0.0, 1.0), 0.0).ravel()
        A = self.A0 + sp.diags(qphi + self.diag0)
        b = self.b.ravel()
        ml = self.pyamg.smoothed_aggregation_solver(A, max_coarse=500)
        res = []
        c = ml.solve(b, x0=self.c.ravel(), tol=1e-9, accel="cg", maxiter=400, residuals=res)
        self.it = len(res)
        self.res = float(np.linalg.norm(A @ c - b) / np.linalg.norm(b))
        self.c = np.clip(c.reshape((NF,) * 3), 0.0, 1.0)
        return self.c


# ---------------------------------------------------------------- mechanics (as anim_real_tooth_3d.py)
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
        self.tooth_m = down(tooth.astype(float), n) > 0.5

    def solve(self, phi_m, am1_m):
        n = self.n
        tm = self.tooth_m.ravel()
        Ee = np.where(tm, E_TOOTH, E_BIO * (np.clip(phi_m.ravel(), 0, 1) ** 2 + FLOOR))
        eps0 = np.where(tm, 0.0, am1_m.ravel())
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


# ---------------------------------------------------------------- the threshold from the base run
def p_threshold():
    """p_0 = mean pressure in the film (phi > 0.5 on the stress mesh) at the end of the base run; p_h = 0.5 p_0."""
    B = np.load(BASE_CACHE)
    pm = down(B["phi"][-1].astype(float), NM)
    sel = (pm > 0.5) & ~tooth_m
    p0 = float(B["p"][-1][sel].mean())
    return p0, P_H_FRACTION * p0


# ---------------------------------------------------------------- run
def run(nutrient, gate, q=None, tag=None, stress_every=STRESS_EVERY, nm=NM, verbose=True):
    """Field run with the nutrient and/or the gate; returns the cache dict and writes figdata/real_tooth_<tag>.npz.
    Without the gate the stress is solved every stress_every frames and interpolated linearly in between; with the
    gate it is solved every frame (the gate needs the pressure at the start of each frame)."""
    tag = tag or ("nutrient_" if nutrient else "") + ("gate" if gate else "")
    cache = cache_path(tag)
    if cache.exists():
        return dict(np.load(cache))
    walls = {}
    dt = 0.15 * H * H / BETA
    per = int(np.ceil(TEND / NFRAMES / dt))
    dt = TEND / (NFRAMES * per)
    if verbose:
        print(f"[{tag}] field: dt = {dt:.4f} T*, {per} steps per frame; nutrient {nutrient} q = {q}; gate {gate}",
              flush=True)
    phi = np.where(seed, 1.0, 0.0)
    alpha = np.ones_like(phi)
    c = np.where(OPEN, 1.0, 0.0) if nutrient else None
    nut = Nutrient(q) if nutrient else None
    mech = Mech(nm)
    p0 = p_h = np.nan
    if gate:
        p0, p_h = p_threshold()
        if verbose:
            print(f"[{tag}] p_0 = {p0 * 1e3:.3f} mPa, p_h = {p_h * 1e3:.3f} mPa", flush=True)
    solved = set(range(0, NFRAMES + 1, 1 if gate else stress_every)) | {NFRAMES}
    T, P, A, C, VMs, PRs, ITs, GS, CGI = [], [], [], [], {}, {}, {}, [], []
    gf = np.ones((NF,) * 3)
    t_field = t_c = t_stress = 0.0
    t_run = time.time()

    def rhs(phi, alpha, cfac, gfac):
        dphi = BETA * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
        return np.where(OPEN, dphi, 0.0), np.where(OPEN, KA * phi * cfac * gfac, 0.0)

    for fr in range(NFRAMES + 1):
        t1 = time.time()
        if nutrient:
            c = nut.solve(phi)
            CGI.append(nut.it)
        t_c += time.time() - t1
        T.append(fr * per * dt)
        P.append(phi.astype(np.float16)); A.append(np.where(OPEN, alpha - 1, 0).astype(np.float16))
        if nutrient:
            C.append(c.astype(np.float16))
        if fr in solved:
            t1 = time.time()
            vm, p, _ = mech.solve(down(phi, nm), down(np.where(OPEN, alpha - 1, 0), nm))
            VMs[fr], PRs[fr], ITs[fr] = vm, p, mech.last_iters
            t_stress += time.time() - t1
            if gate:
                g = np.clip(1 - p / p_h, 0.0, 1.0)
                gf = up(g, NF)
                film_m = (down(phi, nm) > 0.5) & ~mech.tooth_m
                GS.append(float((p[film_m] >= p_h).mean()) if film_m.any() else 0.0)
        if verbose and (fr % 11 == 0 or fr in solved and gate and fr % 3 == 0):
            msg = f"[{tag}] frame {fr}/{NFRAMES}  phi sum {phi.sum() * H ** 3:.3f} mm^3  max alpha-1 " \
                  f"{(alpha - 1)[OPEN].max():.2e}"
            if nutrient:
                msg += f"  c in film [{c[phi > 0.5].min():.3f}, {c[phi > 0.5].max():.3f}] ({nut.it} it, {nut.res:.1e})"
            if fr in solved:
                msg += f"  max vM {VMs[fr][~mech.tooth_m].max():.3e} Pa ({ITs[fr]} it)"
            if gate:
                msg += f"  gate closed {GS[-1]:.2f}"
            print(msg + f"  {time.time() - t_run:.0f} s", flush=True)
        if fr == NFRAMES:
            break
        t1 = time.time()
        cf = c if nutrient else 1.0
        for _ in range(per):
            k1 = rhs(phi, alpha, cf, gf)
            k2 = rhs(phi + 0.5 * dt * k1[0], alpha + 0.5 * dt * k1[1], cf, gf)
            k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1], cf, gf)
            k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1], cf, gf)
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t_field += time.time() - t1
    solved = sorted(solved)
    V = np.zeros((NFRAMES + 1, nm, nm, nm), np.float32)
    PR = np.zeros_like(V)
    for a, b in zip(solved[:-1], solved[1:]):               # linear in time between the solved frames
        for fr in range(a, b + 1):
            w = (fr - a) / (b - a)
            V[fr] = (1 - w) * VMs[a] + w * VMs[b]
            PR[fr] = (1 - w) * PRs[a] + w * PRs[b]
    out = dict(t=np.array(T), phi=np.array(P), am1=np.array(A), vm=V, p=PR, solved=np.array(solved),
               iters=np.array([ITs[f] for f in solved]), nutrient=nutrient, gate=gate,
               q=np.nan if q is None else q, p_0=p0, p_h=p_h, nm=nm, dt=dt, per=per,
               wall_field=t_field, wall_c=t_c, wall_stress=t_stress)
    if nutrient:
        out["c"] = np.array(C)
        out["cg_iters"] = np.array(CGI)
    if gate:
        out["gate_closed_share"] = np.array(GS)
        out["g_end"] = np.clip(1 - PRs[NFRAMES] / p_h, 0.0, 1.0).astype(np.float32)
    cache.parent.mkdir(exist_ok=True)
    np.savez_compressed(cache, **out)
    if verbose:
        print(f"[{tag}] wrote {cache}; field {t_field:.0f} s, nutrient {t_c:.0f} s, stress {t_stress:.0f} s",
              flush=True)
    return out


# ---------------------------------------------------------------- numbers
def stats(D):
    """End-of-run numbers of one cache: seed and film means, nutrient in the film, cervical band versus the rest."""
    nm = int(D["nm"]) if "nm" in D else NM
    tm = down(tooth.astype(float), nm) > 0.5
    band_m = (np.meshgrid(*[-L / 2 + (np.arange(nm) + 0.5) * L / nm] * 3, indexing="ij")[2] > Z_CEJ - BAND_HALF) & \
             (np.meshgrid(*[-L / 2 + (np.arange(nm) + 0.5) * L / nm] * 3, indexing="ij")[2] < Z_CEJ + BAND_HALF)
    phi = D["phi"][-1].astype(float)
    am1 = D["am1"][-1].astype(float)
    vm, p = D["vm"][-1].astype(float), D["p"][-1].astype(float)
    c = D["c"][-1].astype(float) if "c" in D else np.ones_like(phi)
    film_f = (phi > 0.5) & OPEN
    pm = down(phi, nm)
    film_m = (pm > 0.5) & ~tm
    film_m3 = (pm > 0.3) & ~tm                               # the wider film of anim_real_tooth_3d.py (one element thick)
    sb, sr = film_m & band_m, film_m & ~band_m
    out = dict(
        seed=dict(phi_mean=float(phi[seed].mean()), am1_mean=float(am1[seed].mean()),
                  c_mean=float(c[seed].mean())),
        film=dict(cells=int(film_f.sum()), volume_phi_mm3=float(phi.sum() * H ** 3),
                  am1_mean=float(am1[film_f].mean()), am1_max=float(am1[OPEN].max()),
                  c_min=float(c[film_f].min()), c_mean=float(c[film_f].mean()), c_max=float(c[film_f].max()),
                  elements=int(film_m.sum()), vm_mean_Pa=float(vm[film_m].mean()), vm_max_Pa=float(vm[film_m].max()),
                  vm_p95_Pa=float(np.percentile(vm[film_m], 95)), p_mean_Pa=float(p[film_m].mean()),
                  p_max_Pa=float(p[film_m].max()),
                  elements_phi_gt_0p3=int(film_m3.sum()), vm_mean_phi_gt_0p3_Pa=float(vm[film_m3].mean()),
                  vm_max_phi_gt_0p3_Pa=float(vm[film_m3].max())),
        band=dict(z_mm=[Z_CEJ - BAND_HALF, Z_CEJ + BAND_HALF], elements=int(sb.sum()), rest_elements=int(sr.sum()),
                  vm_mean_Pa=float(vm[sb].mean()), vm_p95_Pa=float(np.percentile(vm[sb], 95)),
                  vm_max_Pa=float(vm[sb].max()), p_mean_Pa=float(p[sb].mean()),
                  rest_vm_mean_Pa=float(vm[sr].mean()), rest_vm_p95_Pa=float(np.percentile(vm[sr], 95)),
                  rest_vm_max_Pa=float(vm[sr].max()), rest_p_mean_Pa=float(p[sr].mean()),
                  ratio_vm_mean_band_over_rest=float(vm[sb].mean() / vm[sr].mean()),
                  ratio_vm_p95_band_over_rest=float(np.percentile(vm[sb], 95) / np.percentile(vm[sr], 95)),
                  ratio_vm_max_band_over_rest=float(vm[sb].max() / vm[sr].max()),
                  am1_mean=float(am1[film_f & BAND_F].mean()), rest_am1_mean=float(am1[film_f & ~BAND_F].mean())),
    )
    p0, p_h = p_threshold()                                   # the same threshold for every run (from the base run)
    out["film"]["share_p_ge_p_h_end"] = float((p[film_m] >= p_h).mean())
    out["film"]["share_p_ge_p_h_end_phi_gt_0p3"] = float((p[film_m3] >= p_h).mean())
    out["band"]["share_p_ge_p_h_end"] = float((p[sb] >= p_h).mean())
    out["band"]["rest_share_p_ge_p_h_end"] = float((p[sr] >= p_h).mean())
    if bool(D["gate"]):
        out["gate"] = dict(p_0_Pa=float(D["p_0"]), p_h_Pa=float(D["p_h"]),
                           closed_share_film_end=float((p[film_m] >= p_h).mean()),
                           closed_share_band_end=float((p[sb] >= p_h).mean()),
                           closed_share_rest_end=float((p[sr] >= p_h).mean()),
                           closed_share_per_frame=[float(x) for x in D["gate_closed_share"]],
                           g_mean_film_end=float(np.clip(1 - p[film_m] / p_h, 0, 1).mean()))
    return out


def run_record(tag, D):
    """The json entry of one run: settings, solver counts, wall times, end numbers and z profiles."""
    nutrient, gate, q = RUNS[tag]
    r = dict(nutrient=nutrient, gate=gate, g_over_d_per_mm2=q, stress_every=int(np.diff(D["solved"])[0]),
             stress_mesh=int(D["nm"]) if "nm" in D else NM, amg_iterations_max=int(D["iters"].max()))
    if "cg_iters" in D:
        r["nutrient_solver_iterations_max"] = int(D["cg_iters"].max())
    for k in ("wall_field", "wall_c", "wall_stress"):
        if k in D:
            r[k + "_s"] = float(D[k])
    r.update(stats(D))
    r["profiles"] = profiles(D)
    return r


JSON_PATH = RESDIR / "real_tooth_physics.json"


def save_run_json(tag, D):
    """Merge this run's numbers into real_tooth_physics.json right after the run (a restart loses at most one run)."""
    res = json.loads(JSON_PATH.read_text()) if JSON_PATH.exists() else {}
    res.setdefault("runs", {})[tag] = run_record(tag, D)
    RESDIR.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(res, indent=1) + "\n")
    print(f"[{tag}] numbers merged into {JSON_PATH}", flush=True)


def profiles(D):
    """Band averages by z through the film: von Mises on the stress mesh (phi > 0.3, one element thick), alpha - 1
    and c on the field grid (phi > 0.5)."""
    nm = int(D["nm"]) if "nm" in D else NM
    tm = down(tooth.astype(float), nm) > 0.5
    axn = -L / 2 + (np.arange(nm) + 0.5) * L / nm
    phi = D["phi"][-1].astype(float)
    film_m = (down(phi, nm) > 0.3) & ~tm
    vm = D["vm"][-1].astype(float)
    zv = [float(axn[k]) for k in range(nm) if film_m[:, :, k].any()]
    vmz = [float(vm[:, :, k][film_m[:, :, k]].mean()) for k in range(nm) if film_m[:, :, k].any()]
    film_f = (phi > 0.5) & OPEN
    am1 = D["am1"][-1].astype(float)
    za = [float(ax[k]) for k in range(NF) if film_f[:, :, k].any()]
    amz = [float(am1[:, :, k][film_f[:, :, k]].mean()) for k in range(NF) if film_f[:, :, k].any()]
    out = dict(z_vm_mm=zv, vm_Pa=vmz, z_field_mm=za, am1=amz)
    if "c" in D:
        c = D["c"][-1].astype(float)
        out["c"] = [float(c[:, :, k][film_f[:, :, k]].mean()) for k in range(NF) if film_f[:, :, k].any()]
        out["c_min"] = [float(c[:, :, k][film_f[:, :, k]].min()) for k in range(NF) if film_f[:, :, k].any()]
    return out


# ---------------------------------------------------------------- render (as anim_real_tooth_3d.py)
DEPTH = 1.5                                              # mm, sampling depth of the film colouring
CAM = [(-3.6 * L / 2, -1.9 * L / 2, 1.3 * L / 2), (0, 0, -0.7), (0, 0, 1)]


def fill_nearest(v, bad):
    if bad.any():
        ind = distance_transform_edt(bad, return_distances=False, return_indices=True)
        v = v[tuple(ind)]
    return v


def vm_fine(vm):
    """A stress-mesh field onto the field cells: tooth elements blanked, nearest film value, then linear."""
    nm = vm.shape[0]
    tm = down(tooth.astype(float), nm) > 0.5
    axn = -L / 2 + (np.arange(nm) + 0.5) * L / nm
    v = fill_nearest(np.where(tm, np.nan, vm), tm)
    f = RegularGridInterpolator((axn,) * 3, v, bounds_error=False, fill_value=None)
    return f(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(X.shape)


def film_surface(field, keep, reduce="max"):
    """Envelope of the cells with phi > 0.05 (smoothed), coloured by the maximum (or minimum) of the field through
    the film: sampled along the inward normal over DEPTH mm."""
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
    f = RegularGridInterpolator((ax,) * 3, np.nan_to_num(field), bounds_error=False,
                                fill_value=1.0 if reduce == "min" else 0.0)
    depths = np.linspace(0.0, DEPTH, 10)
    samples = np.stack([f(pts - d * nrm) for d in depths], 0)
    surf.point_data["v"] = {"max": samples.max, "min": samples.min, "mean": samples.mean}[reduce](0)
    return surf


def render(tag, D, png_only=False):
    pv.OFF_SCREEN = True
    pv.global_theme.font.family = "times"
    t, PHI, AM1, VM = D["t"], D["phi"], D["am1"], D["vm"]
    nf = len(t)
    nutrient, gate = bool(D["nutrient"]), bool(D["gate"])
    st = stats(D)
    vtop = st["film"]["vm_p95_Pa"] * 1e3
    atop = float(AM1[-1].astype(float).max()) * 1e3
    tooth_draw = tooth_mesh.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
    tooth_draw.points = np.asarray(tooth_draw.points) - 0.5 * H * np.asarray(tooth_draw["Normals"])
    tooth_draw = tooth_draw.clip_box([-L / 2, L / 2, -L / 2, L / 2, -L / 2, L / 2], invert=False)
    KEEP = lambda i: (PHI[i].astype(float) > 0.05) & OPEN       # noqa: E731
    phi_panel = (lambda i: (PHI[i].astype(float), KEEP(i), "max"), "phase phi", "Greens", (0, 1), "%.1f")
    am1_panel = (lambda i: (AM1[i].astype(float) * 1e3, KEEP(i), "max"), "alpha - 1 (1e-3)", "Oranges",
                 (0, atop), "%.2f")
    vm_panel = (lambda i: (vm_fine(VM[i]) * 1e3, KEEP(i), "max"), "von Mises (mPa), film", "viridis", (0, vtop), "%.1f")
    if nutrient:
        C = D["c"]
        ctop = float(np.ceil(st["film"]["c_max"] * 20) / 20)           # the film's own range (c is 1 at the box faces)
        c_panel = (lambda i: (fill_nearest(C[i].astype(float), tooth), KEEP(i), "mean"), "nutrient c, film mean",
                   "Blues", (0.0, ctop), "%.2f")
        panels = [c_panel, am1_panel, vm_panel]
    else:
        panels = [phi_panel, am1_panel, vm_panel]
    note = (f"molar FDI 46 (grey, E = 1e4 Pa) at real size, cube {L:g} mm; film 0.5 mm, beta = {BETA:.2f} mm^2/T*; "
            f"bottom face clamped; undeformed; film surface = phi > 0.05, colour = maximum of the field through the film")
    if nutrient:
        note += f"; c = 1 on the box faces around the crown, g/d = {float(D['q']):g} /mm^2, c panel = mean through the film"
    if gate:
        note += f"; growth gate 1 - p/p_h, p_h = {float(D['p_h']) * 1e3:.2f} mPa (0.5 x the film mean pressure without the gate)"
    pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
    acts = []
    for j, (get, lab, cm, clim, fmt) in enumerate(panels):
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
    pl.add_text(note, position="lower_right", font_size=8, font="times", color="#555555")
    pl.link_views()

    def frame(i, azim):
        for a, get in acts:
            a.mapper.SetInputData(film_surface(*get(i)))
        title.SetText(2, f"T* = {t[i]:.2f}")
        pl.camera.azimuth = azim

    out = ROOT / "assets" / f"anim_real_tooth_{tag}_3d"
    out.parent.mkdir(exist_ok=True)
    cam0 = pl.camera.azimuth
    t0 = time.time()
    if png_only:
        frame(nf - 1, cam0)
        pl.camera_position = CAM                             # the movie ends 30 degrees past CAM (its last 24 frames
        pl.camera.Azimuth(30.0)                              # turn to +30) and its screenshot keeps that view
        pl.render()
        pl.screenshot(str(out.with_suffix(".png")))
    else:
        pl.open_movie(str(out.with_suffix(".mp4")), framerate=16, quality=9)
        for i in range(nf):
            frame(i, cam0 + 30.0 * i / (nf - 1) - 15.0)
            pl.write_frame()
        for k in range(24):
            frame(nf - 1, cam0 + 15.0 + 15.0 * (k + 1) / 24)
            pl.write_frame()
        frame(nf - 1, cam0)
        pl.screenshot(str(out.with_suffix(".png")))
    pl.close()
    print(f"[{tag}] wrote {out.with_suffix('.png' if png_only else '.mp4')} ({time.time() - t0:.0f} s)", flush=True)


# ---------------------------------------------------------------- summary: json + figure
def summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply(size=9)
    res = {"model": {
        "phi": "d phi/dt = beta lap phi + k_alpha alpha - penalty (max(phi - 1, 0) + min(phi, 0))  [Eq. 34, no front term]",
        "c": "0 = d lap c - g phi c, c = 1 on the five outer box faces except the clamped bottom, zero flux at the "
             "bottom and at the tooth  [Eq. 35 quasi-static, consumption g phi c]",
        "alpha": "d alpha/dt = k_alpha phi c g  [Eq. 36 with the factor c (assumption) and the gate g]",
        "gate": "g = min(1, max(0, 1 - p/p_h)), p = -tr(sigma)/3 at the start of the frame, held over the frame",
        "mechanics": "small strain, hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), tooth 1e4 Pa without growth, nu = 0.49, "
                     "eigenstrain (alpha - 1) I, bottom face clamped",
    }, "parameters": dict(L_mm=L, field_cells=NF, stress_elements=NM, beta_mm2_per_Tstar=BETA, k_alpha=KA, penalty=PEN,
                          T_end=TEND, frames=NFRAMES, g_over_d_per_mm2=10 * Q_TABLE2, g_over_d_table2_per_mm2=Q_TABLE2,
                          penetration_depth_mm=(10 * Q_TABLE2) ** -0.5, penetration_depth_table2_mm=Q_TABLE2 ** -0.5,
                          p_h_fraction_of_p_0=P_H_FRACTION, band_z_mm=[Z_CEJ - BAND_HALF, Z_CEJ + BAND_HALF],
                          film_definition="field grid phi > 0.5; stress mesh phi > 0.5 (p_0, band statistics), "
                                          "phi > 0.3 for the von Mises profile"),
        "runs": {}}
    p0, p_h = p_threshold()
    res["parameters"]["p_0_Pa"], res["parameters"]["p_h_Pa"] = p0, p_h
    caches = {}
    for tag in RUNS:
        cp = cache_path(tag)
        if not cp.exists():
            print("missing", cp)
            continue
        D = dict(np.load(cp))
        if tag == "base":
            D.update(nutrient=False, gate=False, nm=NM, q=np.nan, p_0=np.nan, p_h=np.nan)
        caches[tag] = D
        res["runs"][tag] = run_record(tag, D)
    RESDIR.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(res, indent=1) + "\n")
    # figure
    label = {"base": "no nutrient, no gate", "nutrient": "nutrient", "gate": "gate", "nutrient_gate": "nutrient + gate",
             "nutrient_t2": "nutrient, Table 2 ratio"}
    style = {"base": dict(color="k", ls="-"), "nutrient": dict(color="C0", ls="-"), "gate": dict(color="C3", ls="--"),
             "nutrient_gate": dict(color="C2", ls="-."), "nutrient_t2": dict(color="C0", ls=":")}
    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.3))
    for tag in ("base", "nutrient", "gate", "nutrient_gate", "nutrient_t2"):
        if tag not in res["runs"]:
            continue
        pr = res["runs"][tag]["profiles"]
        axs[0].plot(np.array(pr["vm_Pa"]) * 1e3, pr["z_vm_mm"], marker="o", ms=3, label=label[tag], **style[tag])
        axs[1].plot(np.array(pr["am1"]) * 1e3, pr["z_field_mm"], label=label[tag], **style[tag])
        if "c" in pr:
            axs[2].plot(pr["c"], pr["z_field_mm"], label=label[tag], **style[tag])
    for a in axs:
        a.axhspan(Z_CEJ - BAND_HALF, Z_CEJ + BAND_HALF, color="0.85", zorder=0, lw=0)
        a.set_ylabel("z (mm)")
        a.set_ylim(-3.6, -0.4)
    axs[0].set_xlabel("von Mises (mPa), film mean by z")
    axs[1].set_xlabel(r"$\alpha - 1$ ($10^{-3}$), film mean by z")
    axs[2].set_xlabel("c, film mean by z")
    axs[2].set_xlim(0, 1.05)
    for a, letter in zip(axs, "abc"):
        a.text(0.02, 0.98, f"({letter})", transform=a.transAxes, ha="left", va="top")
    axs[1].legend(loc="lower right", frameon=False, fontsize=6.5)
    fig.tight_layout(w_pad=0.6)
    out = ROOT / "assets" / "fig_real_tooth_physics"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    print("wrote", RESDIR / "real_tooth_physics.json", "and", out.with_suffix(".pdf"))
    return res


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--summary" in sys.argv:
        summary()
        sys.exit(0)
    tag = args[0]
    nutrient, gate, q = RUNS[tag]
    D = run(nutrient, gate, q=q, tag=tag, stress_every=NFRAMES if tag == "nutrient_t2" else STRESS_EVERY)
    if tag == "base":
        D.update(nutrient=False, gate=False, nm=NM, q=np.nan, p_0=np.nan, p_h=np.nan)
    print(json.dumps(stats(D), indent=1))
    save_run_json(tag, D)
    if "--no-video" not in sys.argv:
        render(tag, D, png_only="--png-only" in sys.argv)
