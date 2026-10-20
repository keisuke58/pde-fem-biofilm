#!/usr/bin/env python3
"""Three steps from the idealised film towards a realistic plaque shape on the real molar (demonstration).

The idealised film of anim_real_tooth_3d.py (0.5 mm, one flank, 120 degree sector, phi = 1 inside) is replaced in three
steps, each built synthetically from published descriptions of supragingival plaque and from aggregate statistics of
CLSM/FISH stacks of in-vitro biofilms (no per-stack data, see REAL_TOOTH_PLAQUE.md for the sources):

  1  distribution: a cervical band along the gingival margin all round the crown, thicker on the approximal (mesial,
     distal) surfaces and reaching up to the contact area there, thinner on the buccal and lingual surfaces and
     tapering off coronally, plus deposits in the occlusal groove-fossa system (morphological closing of the voxel
     tooth); thickness map t(theta, z) from the signed distance to the STL and the azimuth about the crown axis;
  2  thickness variation: t multiplied by a lognormal Gaussian random field on the tooth surface (evaluated at the foot
     point of each cell), coefficient of variation and correlation length from the CLSM statistics, fixed seed;
  3  internal structure: water channels and voids inside the film (phi = 0) from a short-correlated Gaussian field
     with a porosity that rises from the tooth to the outer surface (mean porosity from the CLSM statistics), so
     that the outer surface breaks up into colonies; E = 10 Pa (phi^2 + 1e-3) then varies inside the film.

Growth model and mechanics as the demonstration (RK4 field with the masked Laplacian, k_alpha = 1e-3 /T*, penalty 5,
T* = 1.1; beta = 0.02 (t_mean / 0.25 mm)^2 from the mean thickness of the map, the rule of REAL_TOOTH.md), field on
192^3 cells (h_f = 0.083 mm; thin parts need three cells, so the thickness floor is T_MIN = 0.25 mm), stress with the
film submodel of real_tooth_submodel.py (rigid tooth, crop, void phi <= 0.05 dropped, voxel wall) at h = h_f / 2
(six elements across the 0.25 mm band, twelve across the approximal film) and h = h_f (the convergence check).

    python ansys_usermat/apdl/real_tooth_plaque.py --geom                   # signed distance, fissure set, axes (cached)
    python ansys_usermat/apdl/real_tooth_plaque.py --field 1 2 3            # the field runs (cached per step)
    python ansys_usermat/apdl/real_tooth_plaque.py --stress 3 --levels r2 r1  # final-state stress (JSON + npz)
    python ansys_usermat/apdl/real_tooth_plaque.py --video-stress 1 3       # stress of every 6th frame at h_f
    python ansys_usermat/apdl/real_tooth_plaque.py --render 1 3 [--png-only]
    python ansys_usermat/apdl/real_tooth_plaque.py --figure
        -> figdata/real_tooth_plaque_geom.npz, figdata/real_tooth_plaque_<step>.npz,
           figdata/real_tooth_plaque_<step>_stress_<level>.npz, figdata/real_tooth_plaque_<step>_vid_f<k>.npz,
           ansys_usermat/apdl/results/2026-10-3d/real_tooth_plaque.json (merged after every run),
           assets/anim_real_tooth_plaque_<step>_3d.mp4 / .png, assets/fig_real_tooth_plaque.pdf / .png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402
from scipy.special import ndtri  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import real_tooth_convergence as RC  # noqa: E402  (load_tooth, make_lap, constants)
import real_tooth_submodel as SM  # noqa: E402  (solve_submodel, Mesh, crop_box)

L, KA, PEN, TEND = RC.L, RC.KA, RC.PEN, RC.TEND
E_BIO, FLOOR, NU = RC.E_BIO, RC.FLOOR, RC.NU
Z_CEJ = RC.Z_CEJ
NF = 192
H = L / NF
NFRAMES, STORE_EVERY, VID_STRESS_EVERY = 66, 3, 6
FIGDATA = ROOT / "figdata"
OUT_JSON = HERE / "results" / "2026-10-3d" / "real_tooth_plaque.json"
FIG = ROOT / "assets" / "fig_real_tooth_plaque"

# ---------------------------------------------------------------- the plaque shape (step 1): parameters
# Every value below is my choice within the published range unless the comment names a source; the report lists them.
Z_GM = Z_CEJ + 0.5          # gingival margin 0.5 mm coronal to the CEJ (healthy margin lies coronal to the CEJ; choice)
Z_HOC = 2.0                 # height of contour / contact area of the crown (widest at z = 1..3 mm in this STL)
T_IP = 0.5                  # approximal thickness, mm: the idealised film's value, upper end of the literature (choice)
T_BL = 0.25                 # buccal / lingual thickness, mm = T_MIN (literature 0.035..0.22 mm is below the grid floor)
T_MIN = 3 * H               # thickness floor: three field cells (0.25 mm)
T_CUT = T_MIN / 2           # below this map value there is no film (bare tooth); between T_CUT and T_MIN the floor applies
LAM_BL = 1.5                # coronal decay length of the buccal / lingual band, mm: with T_CUT the band ends 1.25 mm above the margin (cervical third; choice)
LAM_IP = 0.75               # coronal decay above the contact area on the approximal surfaces, mm (choice)
SIG_THETA = 25.0            # angular half width of the approximal sectors, degrees (choice)
R_FISSURE = 0.6             # closing ball radius for the groove-fossa deposits, mm (choice)
Z_OCCLUSAL = 2.0            # the groove-fossa system lies above this height
# step 2 (CLSM statistics, relative values only; REAL_TOOTH_PLAQUE.md)
COV_T = 0.25                # coefficient of variation of the local thickness (stacks 0.16..0.26 at day >= 6; Al-Ahmad 0.24)
ELL_REL = 1.0               # correlation length / mean thickness (stacks 0.7..1.1, FOV-limited)
SEED_T = 0
# step 3
POROSITY = 0.40             # mean void fraction inside the envelope (stacks 0.36..0.46 at day >= 6)
P_BASE = 0.03               # porosity at the tooth (lateral coverage of the stacks 0.97..1.0)
ELL_C = 2 * H               # correlation length of the channel field, mm (stacks 0.45 x thickness = 0.13 mm; grid floor 2 cells)
SEED_C = 1
BAND_H = 1.0                # cervical band region: z in [Z_GM, Z_GM + BAND_H]
LAYER = 0.1                 # the 0.1 mm next to the tooth
FILM_PHI = RC.FILM_PHI      # 0.3

# the idealised film of REAL_TOOTH_SUBMODEL.md (converged values; its band = z in [-2.1, -1.1])
IDEAL = {"film vM mean": 3.287e-3, "band vM mean": 3.479e-3, "layer vM mean": 7.76e-3, "vM p95": 10.62e-3,
         "p mean": 2.27e-3, "film volume mm^3": 17.33, "t mean mm": 0.5, "source": "real_tooth_submodel.json level r2"}

STEPS = {1: "distribution", 2: "thickness variation", 3: "internal structure"}
LEVELS = {"r2": 0.5, "r1": 1.0}


# ---------------------------------------------------------------- geometry
def axes():
    return -L / 2 + (np.arange(NF) + 0.5) * H


def geometry(force=False):
    """Signed distance on the 192^3 grid, the tooth and fissure sets, the mesio-distal axis (cached)."""
    path = FIGDATA / "real_tooth_plaque_geom.npz"
    if path.exists() and not force:
        d = np.load(path)
        return {k: d[k] for k in d.files}
    import pyvista as pv
    t0 = time.time()
    ax = axes()
    cen = pv.ImageData(dimensions=(NF,) * 3, spacing=(H,) * 3, origin=(ax[0],) * 3)
    tooth_mesh = RC.load_tooth()
    sdf = np.asarray(cen.compute_implicit_distance(tooth_mesh)["implicit_distance"]).reshape((NF,) * 3, order="F")
    sdf = sdf.astype(np.float32)
    tooth = sdf < 0
    t_sdf = time.time() - t0
    # mesio-distal axis: principal axis of the crown outline between z = 1 and 3 mm, checked against the direction
    # between the two root centroids at z = -6 mm (the lower molar's two roots lie mesially and distally)
    P = np.asarray(tooth_mesh.points)
    crown = P[(P[:, 2] > 1.0) & (P[:, 2] < 3.0), :2]
    c = crown.mean(0)
    w, v = np.linalg.eigh(np.cov((crown - c).T))
    md = v[:, np.argmax(w)]
    if md[1] < 0:
        md = -md
    ext = [float(((crown - c) @ v[:, k]).max() - ((crown - c) @ v[:, k]).min()) for k in (np.argmax(w), np.argmin(w))]
    sl = tooth_mesh.slice(normal=(0, 0, 1), origin=(0, 0, -6.0)).connectivity()
    lab = np.asarray(sl["RegionId"]); Q = np.asarray(sl.points)
    cents = [Q[lab == r, :2].mean(0) for r in np.unique(lab) if (lab == r).sum() > 10]
    root_dir = None
    if len(cents) >= 2:
        d = cents[1] - cents[0]
        root_dir = d / np.linalg.norm(d)
        if root_dir[1] < 0:
            root_dir = -root_dir
    # the groove-fossa system: cells filled by a morphological closing of the voxel tooth (ball R_FISSURE) above
    # Z_OCCLUSAL; the closing fills every concavity narrower than about 2 R
    r = int(round(R_FISSURE / H))
    g = np.arange(-r, r + 1)
    ball = (g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) <= r * r
    closed = ndi.binary_closing(tooth, structure=ball, border_value=0)
    Z = ax[None, None, :]
    fissure = closed & ~tooth & (np.broadcast_to(Z, tooth.shape) > Z_OCCLUSAL)
    # the voxel-surface cells (open, next to a tooth cell), for the contact area
    nb = np.zeros_like(tooth)
    for axn in range(3):
        for sgn in (1, -1):
            s = np.roll(tooth, sgn, axn)
            e = [slice(None)] * 3
            e[axn] = 0 if sgn == 1 else -1
            s[tuple(e)] = False
            nb |= s
    surf = nb & ~tooth
    out = dict(sdf=sdf, tooth=tooth, fissure=fissure, surf=surf, md=md.astype(np.float64),
               md_extent_mm=np.array(ext), root_dir=(root_dir if root_dir is not None else np.full(2, np.nan)),
               crown_centre=c.astype(np.float64), seconds_sdf=np.array(t_sdf))
    FIGDATA.mkdir(exist_ok=True)
    np.savez_compressed(path, **out)
    print(f"geometry: sdf {t_sdf:.0f} s; tooth {tooth.sum()} cells; fissure fill {fissure.sum() * H ** 3:.2f} mm^3; "
          f"mesio-distal axis {np.round(md, 3).tolist()} (crown extent {ext[0]:.1f} x {ext[1]:.1f} mm), "
          f"root axis {None if root_dir is None else np.round(root_dir, 3).tolist()}", flush=True)
    return out


def azimuth_weight(G):
    """w_ip in [0, 1]: 1 at the mesial and distal directions (the mesio-distal axis), falling off with the angular
    distance (Gaussian, SIG_THETA); evaluated per cell from its position about the crown axis."""
    ax = axes()
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    md = G["md"]
    cx, cy = G["crown_centre"]
    ang = np.arctan2(Y - cy, X - cx) - np.arctan2(md[1], md[0])
    dev = np.abs((ang + np.pi / 2) % np.pi - np.pi / 2)        # angular distance to the nearest of the two directions
    w = np.exp(-(np.degrees(dev) / SIG_THETA) ** 2)
    return w[:, :, None]                                         # (x, y, 1)


def thickness_map(G, step):
    """t(x) on the field grid, mm, before the floor: the band (step 1), times the lognormal random field evaluated at
    the surface foot point of each cell (step 2 and 3). Returns t (float32, 0 where no film) and the random field g."""
    ax = axes()
    Z = np.broadcast_to(ax[None, None, :], (NF,) * 3)
    w = azimuth_weight(G)
    f_bl = np.where(Z >= Z_GM, np.exp(-((Z - Z_GM) / LAM_BL) ** 2), 0.0)
    f_ip = np.where(Z >= Z_GM, np.where(Z <= Z_HOC, 1.0, np.exp(-((Z - Z_HOC) / LAM_IP) ** 2)), 0.0)
    t = w * T_IP * f_ip + (1 - w) * T_BL * f_bl
    g = None
    if step >= 2:
        rng = np.random.default_rng(SEED_T)
        t_mean = mean_thickness(G, t)["t mean mm"]
        ell = ELL_REL * t_mean
        white = rng.standard_normal((NF,) * 3)
        g = ndi.gaussian_filter(white, ell / (2 * H), mode="nearest")   # covariance exp(-r^2 / ell^2)
        g /= g.std()
        # evaluate at the foot point x - sdf grad(sdf) so that the field is a surface field (constant across the film)
        sdf = G["sdf"].astype(np.float64)
        gx, gy, gz = np.gradient(sdf, H)
        n = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2) + 1e-12
        X, Y, Zc = np.meshgrid(ax, ax, ax, indexing="ij")
        band = (sdf >= 0) & (sdf < 1.5 * T_IP)
        foot = np.stack([(X - sdf * gx / n)[band], (Y - sdf * gy / n)[band], (Zc - sdf * gz / n)[band]], 1)
        gi = RegularGridInterpolator((ax,) * 3, g, bounds_error=False, fill_value=0.0)(foot)
        gs = np.zeros_like(g)
        gs[band] = gi
        g = gs
        s_g = np.sqrt(np.log(1 + COV_T ** 2))
        t = t * np.exp(s_g * g - 0.5 * s_g ** 2)
    return t.astype(np.float32), (None if g is None else g.astype(np.float32))


def film_set(G, t):
    """Cells of the initial film: open, within the thickness of the map (floor T_MIN applied), above the gingival
    margin; plus the groove-fossa fill."""
    sdf = G["sdf"]
    tt = np.where(t >= T_CUT, np.maximum(t, T_MIN), 0.0)
    film = ~G["tooth"] & (sdf >= 0) & (sdf < tt)
    return film | G["fissure"], tt


def mean_thickness(G, t):
    """Mean thickness = volume / contact area, the area from the voxel surface cells under the film corrected by the
    staircase factor 1.5 (a plane at a general orientation is counted up to that factor on a voxel surface) ... here
    instead the harmonic mean of t over the film cells: V / sum(h^3 / t) = the volume divided by the area it covers."""
    film = ~G["tooth"] & (G["sdf"] >= 0) & (G["sdf"] < np.where(t >= T_CUT, np.maximum(t, T_MIN), 0))
    tt = np.maximum(t[film], T_MIN)
    vol = film.sum() * H ** 3
    area = (H ** 3 / tt).sum()
    return {"t mean mm": float(vol / area) if area > 0 else 0.0, "band volume mm^3": float(vol),
            "contact area mm^2": float(area), "t volume-weighted mm": float(tt.mean()) if film.any() else 0.0}


def porosity_profile(zeta):
    """Void fraction as a function of the relative height zeta = sdf / t in the film: linear from P_BASE at the tooth
    to P_TOP at the envelope with the mean POROSITY."""
    p_top = 2 * POROSITY - P_BASE
    return np.clip(P_BASE + (p_top - P_BASE) * zeta, 0, 0.95)


def build_phi0(G, step):
    """The initial phi on the field grid for a step, with the map and the diagnostics."""
    t, g = thickness_map(G, step)
    film, tt = film_set(G, t)
    phi = film.astype(np.float64)
    info = {"step": step, "name": STEPS[step]}
    info.update(mean_thickness(G, t))
    info["fissure volume mm^3"] = float(G["fissure"].sum() * H ** 3)
    info["film cells"] = int(film.sum())
    info["film volume mm^3"] = float(film.sum() * H ** 3)
    tf = tt[film & ~G["fissure"]]
    info["t on the band: min/mean/max mm"] = [float(tf.min()), float(tf.mean()), float(tf.max())] if tf.size else None
    info["fraction of band cells at the floor"] = float((tf <= T_MIN + 1e-6).mean()) if tf.size else None
    if step >= 2:
        band = (G["sdf"] >= 0) & (G["sdf"] < H) & (t >= T_CUT)
        tb = np.maximum(t[band], T_MIN)
        info["surface t: cov after floor"] = float(tb.std() / tb.mean())
        info["surface t: cov before floor"] = float(t[band].std() / t[band].mean())
        info["random field: cov"] = COV_T
        info["random field: correlation length mm"] = float(ELL_REL * mean_thickness(G, thickness_map(G, 1)[0])["t mean mm"])
    chan = None
    if step == 3:
        rng = np.random.default_rng(SEED_C)
        gc = ndi.gaussian_filter(rng.standard_normal((NF,) * 3), ELL_C / (2 * H), mode="nearest")
        gc /= gc.std()
        zeta = np.clip(G["sdf"] / np.maximum(tt, T_MIN), 0, 1)
        thr = ndtri(porosity_profile(zeta))
        chan = film & (gc < thr)
        chan &= ~((G["sdf"] >= 0) & (G["sdf"] < H) & (gc >= ndtri(P_BASE)))   # the first cell layer keeps P_BASE
        phi[chan] = 0.0
        info["porosity: mean inside the envelope"] = float(chan.sum() / film.sum())
        info["porosity: profile"] = {"base": P_BASE, "top": float(2 * POROSITY - P_BASE), "mean target": POROSITY}
        info["channel field: correlation length mm"] = ELL_C
        first = film & (G["sdf"] < H)
        info["coverage of the first cell layer"] = float(1 - chan[first].mean())
        info["film volume with channels mm^3"] = float(phi.sum() * H ** 3)
    return phi, dict(t=t, tt=tt.astype(np.float32), film=film, g=g, chan=chan), info


# ---------------------------------------------------------------- field run
def beta_of(t_mean):
    return RC.BETA0 * (t_mean / 0.25) ** 2


def cache_path(step):
    return FIGDATA / f"real_tooth_plaque_{step}.npz"


def z_store():
    """Index range of the stored frames in z (the crown and the band; the roots are not needed for the videos)."""
    ax = axes()
    k0 = int(np.searchsorted(ax, Z_GM - 0.6))
    return k0, NF


def run_field(step, force=False):
    if cache_path(step).exists() and not force:
        print(f"step {step}: field cached, skipped", flush=True)
        return
    G = geometry()
    OPEN = ~G["tooth"]
    phi, aux, info = build_phi0(G, step)
    t_mean = build_phi0(G, 1)[2]["t mean mm"] if step > 1 else info["t mean mm"]   # one beta for the three steps
    beta = beta_of(t_mean)
    dt = min(0.15 * H * H / beta, TEND / NFRAMES)
    per = int(np.ceil(TEND / NFRAMES / dt))
    dt = TEND / (NFRAMES * per)
    print(f"step {step} ({STEPS[step]}): t_mean {t_mean:.4f} mm -> beta {beta:.4f} mm^2/T*, dt {dt:.4f}, "
          f"{NFRAMES * per} steps; film {info['film volume mm^3']:.2f} mm^3", flush=True)
    lap = RC.make_lap(OPEN, H)
    alpha = np.ones_like(phi)
    k0, k1 = z_store()
    T, P, A = [], [], []

    def rhs(phi, alpha):
        dphi = beta * lap(phi) + KA * alpha - PEN * (np.maximum(phi - 1, 0) + np.minimum(phi, 0))
        return np.where(OPEN, dphi, 0.0), np.where(OPEN, KA * phi, 0.0)

    t0 = time.time()
    for fr in range(NFRAMES + 1):
        if fr % STORE_EVERY == 0 or fr == NFRAMES:
            T.append(fr * per * dt)
            P.append(phi[:, :, k0:k1].astype(np.float16))
            A.append(np.where(OPEN, alpha - 1, 0)[:, :, k0:k1].astype(np.float16))
        if fr == NFRAMES:
            break
        for _ in range(per):
            k1_ = rhs(phi, alpha)
            k2 = rhs(phi + 0.5 * dt * k1_[0], alpha + 0.5 * dt * k1_[1])
            k3 = rhs(phi + 0.5 * dt * k2[0], alpha + 0.5 * dt * k2[1])
            k4 = rhs(phi + dt * k3[0], alpha + dt * k3[1])
            phi = phi + dt / 6 * (k1_[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            alpha = alpha + dt / 6 * (k1_[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if fr % 11 == 0:
            print(f"  frame {fr}/{NFRAMES}  phi volume {phi.sum() * H ** 3:.3f} mm^3  max alpha-1 "
                  f"{(alpha - 1)[OPEN].max():.2e}  {time.time() - t0:.0f} s", flush=True)
    am1 = np.where(OPEN, alpha - 1.0, 0.0)
    info.update({"beta mm^2/T*": beta, "dt": dt, "steps": NFRAMES * per, "seconds field": time.time() - t0,
                 "phi volume end mm^3": float(phi[OPEN].sum() * H ** 3), "am1 max": float(am1[OPEN].max()),
                 "am1 mean in film": float(am1[aux["film"]].mean()), "phi mean in film": float(phi[aux["film"]].mean())})
    if step == 3:
        info["porosity at the end (phi < 0.3 in the initial envelope)"] = float((phi[aux["film"]] < FILM_PHI).mean())
        info["porosity at the end (phi < 0.5 in the initial envelope)"] = float((phi[aux["film"]] < 0.5).mean())
    np.savez_compressed(cache_path(step), phi=phi.astype(np.float32), am1=am1.astype(np.float32),
                        t=aux["tt"], film0=aux["film"], g=(aux["g"] if aux["g"] is not None else np.zeros(0, np.float32)),
                        chan=(aux["chan"] if aux["chan"] is not None else np.zeros(0, bool)),
                        frames_t=np.array(T), frames_phi=np.array(P), frames_am1=np.array(A), z0=k0,
                        info=json.dumps(info))
    merge_json({"steps": {str(step): {"field": info}}})
    print(f"step {step}: field done in {time.time() - t0:.0f} s; " +
          ", ".join(f"{k} {v:.4g}" for k, v in info.items() if isinstance(v, float)), flush=True)


def load_step(step):
    d = np.load(cache_path(step))
    return {k: d[k] for k in d.files}


# ---------------------------------------------------------------- regions and measures
def regions(G, phi, film0):
    """Element sets on the field grid at the final state: film = phi > 0.3 and open; band = film within BAND_H of the
    gingival margin; approximal / buccal-lingual = film by the azimuth weight (w > 0.5 / <= 0.5) outside the fissures;
    fissure = film in the groove-fossa fill; layer = film within 0.1 mm of the tooth."""
    ax = axes()
    Z = np.broadcast_to(ax[None, None, :], (NF,) * 3)
    w = np.broadcast_to(azimuth_weight(G), (NF,) * 3)
    film = (phi > FILM_PHI) & ~G["tooth"]
    fis = film & G["fissure"]
    return {"film": film, "band": film & (Z >= Z_GM) & (Z < Z_GM + BAND_H), "approximal": film & (w > 0.5) & ~fis,
            "buccal-lingual": film & (w <= 0.5) & ~fis, "fissure": fis, "layer": film & (G["sdf"] < LAYER),
            "core": film & ~((Z >= Z_GM) & (Z < Z_GM + BAND_H))}


def region_measures(res, G, phi, film0):
    """Volume-weighted means of the von Mises stress and the pressure per region on a submodel result."""
    M, vm, p, act, te = res["mesh"], res["vm"], res["p"], res["active"], res["tooth"]
    R = regions(G, phi, film0)
    vmf = np.nan_to_num(np.where(act, vm, 0.0))
    pf = np.nan_to_num(np.where(act, p, 0.0))
    out = {}
    for k, m in R.items():
        wk = M.sample(m.astype(float))
        wk[~act | te] = 0.0
        out[f"{k} vM mean"] = RC.wmean(vmf, wk)
        out[f"{k} p mean"] = RC.wmean(pf, wk)
        out[f"{k} volume mm^3"] = float(wk.sum() * M.h ** 3)
        if k == "film":
            sel = wk > 0
            out["vM p95"] = RC.wpercentile(vmf[sel], wk[sel], 95)
            out["vM p99"] = RC.wpercentile(vmf[sel], wk[sel], 99)
            out["vM median"] = RC.wpercentile(vmf[sel], wk[sel], 50)
            out["vM max"] = float(vmf[sel].max())
    out["band/core"] = out["band vM mean"] / out["core vM mean"] if out.get("core vM mean") else None
    return out


def stress_path(step, level):
    return FIGDATA / f"real_tooth_plaque_{step}_stress_{level}.npz"


def run_stress(step, levels, force=False):
    G = geometry()
    D = load_step(step)
    phi, am1 = D["phi"].astype(np.float64), D["am1"].astype(np.float64)
    tooth, sdf = G["tooth"], G["sdf"].astype(np.float64)
    for lev in levels:
        if stress_path(step, lev).exists() and not force:
            print(f"step {step} level {lev}: cached, skipped", flush=True)
            continue
        ratio = LEVELS[lev]
        print(f"step {step} level {lev}: h = {ratio * H:.4f} mm", flush=True)
        res = SM.solve_submodel(phi, am1, ratio * H, tooth=tooth, sdf=sdf, h_f=H, void="drop", void_thr=0.05,
                                tooth_mode="voxel")
        meas = region_measures(res, G, phi, D["film0"])
        meas.update({"h": ratio * H, "h/h_f": ratio, "solver": res["solver"], "seconds": res["seconds"],
                     "submodel measures": {k: v for k, v in res["measures"].items()
                                           if isinstance(v, (int, float, str, bool)) or v is None}})
        M = res["mesh"]
        np.savez_compressed(stress_path(step, lev), vm=res["vm"].astype(np.float32), p=res["p"].astype(np.float32),
                            active=res["active"], tooth=res["tooth"], origin=M.origin, h=M.h, shape=np.array(M.shape))
        merge_json({"steps": {str(step): {"stress": {lev: meas}}}})
        print(f"  step {step} {lev}: film {meas['film vM mean'] * 1e3:.3f} band {meas['band vM mean'] * 1e3:.3f} "
              f"approx {meas['approximal vM mean'] * 1e3:.3f} b/l {meas['buccal-lingual vM mean'] * 1e3:.3f} "
              f"fissure {(meas['fissure vM mean'] or 0) * 1e3:.3f} layer {meas['layer vM mean'] * 1e3:.3f} "
              f"p95 {meas['vM p95'] * 1e3:.3f} p {meas['film p mean'] * 1e3:.3f} mPa, {res['seconds']:.0f} s",
              flush=True)


# ---------------------------------------------------------------- video stress (every 6th frame, h = h_f)
def vid_path(step, k):
    return FIGDATA / f"real_tooth_plaque_{step}_vid_f{k:02d}.npz"


def vid_frames(nstored):
    """Indices (into the stored frames, every 3rd field frame) of the frames solved for the video."""
    s = list(range(0, nstored, VID_STRESS_EVERY // STORE_EVERY))
    if s[-1] != nstored - 1:
        s.append(nstored - 1)
    return s


def run_video_stress(step, force=False):
    G = geometry()
    D = load_step(step)
    tooth, sdf = G["tooth"], G["sdf"].astype(np.float64)
    k0 = int(D["z0"])
    n = len(D["frames_t"])
    box = SM.crop_box(D["phi"].astype(np.float64), tooth, H)
    for k in vid_frames(n):
        if vid_path(step, k).exists() and not force:
            continue
        phi = np.zeros((NF,) * 3); am1 = np.zeros((NF,) * 3)
        phi[:, :, k0:] = D["frames_phi"][k].astype(np.float64)
        am1[:, :, k0:] = D["frames_am1"][k].astype(np.float64)
        t0 = time.time()
        if am1.max() <= 0:
            M = SM.Mesh(box, H, 1.0)
            te = M.sample(tooth.astype(float)) > 0.5
            act = ~te & (np.clip(M.sample(phi), 0, 1) > 0.05)
            vm = np.where(act, 0.0, np.nan)
            rec = dict(vm=vm, active=act, origin=M.origin, h=M.h, shape=np.array(M.shape), measures="{}")
        else:
            r = SM.solve_submodel(phi, am1, H, tooth=tooth, sdf=sdf, h_f=H, box=box, void="drop", void_thr=0.05,
                                  tooth_mode="voxel", verbose=False)
            M = r["mesh"]
            m = region_measures(r, G, phi, D["film0"])
            rec = dict(vm=r["vm"], active=r["active"], origin=M.origin, h=M.h, shape=np.array(M.shape),
                       measures=json.dumps(_js(m)))
        tmp = vid_path(step, k).with_suffix(".tmp.npz")
        np.savez_compressed(tmp, vm=rec["vm"].astype(np.float32), active=rec["active"], origin=rec["origin"], h=rec["h"],
                            shape=rec["shape"], t=float(D["frames_t"][k]), frame=k, measures=rec["measures"],
                            seconds=time.time() - t0)
        os.replace(tmp, vid_path(step, k))
        m = json.loads(rec["measures"])
        print(f"[step {step}] video frame {k}/{n - 1} t = {D['frames_t'][k]:.3f}  "
              f"{'film %.3f mPa' % (m['film vM mean'] * 1e3) if m else 'no growth'}  {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------- JSON
def _js(o):
    if isinstance(o, dict):
        return {str(k): _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _deep(a, b):
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            _deep(a[k], v)
        else:
            a[k] = v
    return a


def merge_json(part):
    out = json.loads(OUT_JSON.read_text()) if OUT_JSON.exists() else base_json()
    _deep(out, _js(part))
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1) + "\n")


def base_json():
    G = geometry()
    return {"model": {"field grid": NF, "h_f mm": H, "k_alpha 1/T*": KA, "penalty": PEN, "T*": TEND,
                      "beta rule": "0.02 (t_mean / 0.25 mm)^2 mm^2/T*, t_mean of the step-1 map, one value for the three steps",
                      "E_bio Pa": E_BIO, "floor": FLOOR, "nu": NU,
                      "stress": "film submodel (real_tooth_submodel.py): rigid tooth, crop, void phi <= 0.05 dropped, voxel wall",
                      "film": f"cells with phi > {FILM_PHI}, not tooth, at the final state",
                      "regions": {"band": f"film with z in [{Z_GM}, {Z_GM + BAND_H}] mm (gingival margin + {BAND_H} mm)",
                                  "approximal": f"film with azimuth weight > 0.5 (within about {SIG_THETA * 0.83:.0f} deg of the mesio-distal axis), not fissure",
                                  "buccal-lingual": "film with azimuth weight <= 0.5, not fissure",
                                  "fissure": f"film in the groove-fossa fill (closing, ball {R_FISSURE} mm, z > {Z_OCCLUSAL})",
                                  "layer": f"film within {LAYER} mm of the tooth surface", "core": "film outside the band"}},
            "shape": {"Z_GM mm": Z_GM, "Z_HOC mm": Z_HOC, "T_IP mm": T_IP, "T_BL mm": T_BL, "T_MIN mm": T_MIN,
                      "T_CUT mm": T_CUT, "LAM_BL mm": LAM_BL, "LAM_IP mm": LAM_IP, "SIG_THETA deg": SIG_THETA, "R_FISSURE mm": R_FISSURE,
                      "COV_T": COV_T, "ELL_REL": ELL_REL, "SEED_T": SEED_T, "POROSITY": POROSITY, "P_BASE": P_BASE,
                      "ELL_C mm": ELL_C, "SEED_C": SEED_C,
                      "mesio-distal axis (x, y)": G["md"].tolist(), "crown extent along / across mm": G["md_extent_mm"].tolist(),
                      "root separation direction (x, y)": G["root_dir"].tolist(), "crown centre (x, y) mm": G["crown_centre"].tolist(),
                      "fissure fill mm^3": float(G["fissure"].sum() * H ** 3), "seconds sdf": float(G["seconds_sdf"])},
            "idealised film": IDEAL, "steps": {}}


# ---------------------------------------------------------------- render
CAM_DIST = 1.9 * L / 2
DEPTH = 0.8


def _surface(keep):
    import pyvista as pv
    g = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-L / 2,) * 3)
    g.cell_data["keep"] = keep.ravel(order="F").astype(float)
    surf = g.threshold(0.5, scalars="keep").extract_surface(algorithm="dataset_surface")
    if surf.n_points == 0:
        return surf
    surf = surf.smooth_taubin(n_iter=40, pass_band=0.05, feature_angle=150.0)
    return surf.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)


def _colour(surf, f, depth=DEPTH):
    if surf.n_points == 0:
        surf.point_data["v"] = np.zeros(0)
        return surf
    pts = np.asarray(surf.points); nrm = np.asarray(surf["Normals"])
    d = np.linspace(0.0, depth, 8)
    surf.point_data["v"] = np.stack([f(pts - dd * nrm) for dd in d], 0).max(0)
    return surf


def load_vid(step, n):
    VM, ACT, meas, axes_, h = {}, {}, {}, None, None
    for k in vid_frames(n):
        d = np.load(vid_path(step, k))
        vm = d["vm"].astype(np.float64)
        act = np.asarray(d["active"], bool) & np.isfinite(vm)
        if act.any():
            ind = ndi.distance_transform_edt(~act, return_distances=False, return_indices=True)
            vm = np.where(act, vm, 0.0)[tuple(ind)]
        VM[k] = np.nan_to_num(vm).astype(np.float32); ACT[k] = act
        meas[k] = json.loads(str(d["measures"]))
        if axes_ is None:
            o, h, shp = d["origin"], float(d["h"]), d["shape"]
            axes_ = [o[i] + (np.arange(int(shp[i])) + 0.5) * h for i in range(3)]
    return dict(frames=vid_frames(n), vm=VM, act=ACT, measures=meas, axes=axes_, h=h)


def render(step, png_only=False, vtop=None):
    import pyvista as pv
    pv.OFF_SCREEN = True
    pv.global_theme.font.family = "times"
    G = geometry()
    D = load_step(step)
    info = json.loads(str(D["info"]))
    ax = axes()
    k0 = int(D["z0"])
    T, PH, AM = D["frames_t"], D["frames_phi"], D["frames_am1"]
    n = len(T)
    V = load_vid(step, n)
    OPEN = ~G["tooth"]
    full = lambda a: np.concatenate([np.zeros((NF, NF, k0), a.dtype), a], 2)             # noqa: E731

    def phi_at(i):
        """phi, alpha - 1 at a video frame i (3 per stored interval), linear in time."""
        a, w = divmod(i, 3)
        a = min(a, n - 1); b = min(a + 1, n - 1); w = w / 3 if b > a else 0.0
        phi = (1 - w) * full(PH[a].astype(np.float32)) + w * full(PH[b].astype(np.float32))
        am = (1 - w) * full(AM[a].astype(np.float32)) + w * full(AM[b].astype(np.float32))
        return phi, am, (1 - w) * T[a] + w * T[b]

    def vm_at(i):
        fr = V["frames"]
        s = i / 3
        b = next(f for f in fr if f >= s)
        if b == s or b == fr[0]:
            return V["vm"][b]
        a = fr[fr.index(b) - 1]
        w = (s - a) / (b - a)
        return (1 - w) * V["vm"][a] + w * V["vm"][b]

    f_field = lambda a: RegularGridInterpolator((ax,) * 3, a, bounds_error=False, fill_value=0.0)      # noqa: E731
    f_crop = lambda a: RegularGridInterpolator(tuple(V["axes"]), a, bounds_error=False, fill_value=0.0)  # noqa: E731
    atop = float(AM[-1].astype(float).max()) * 1e3
    if vtop is None:
        k = V["frames"][-1]
        vtop = float(V["measures"][k]["vM p95"]) * 1e3
    tooth_draw = RC.load_tooth().compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
    tooth_draw.points = np.asarray(tooth_draw.points) - 0.5 * H * np.asarray(tooth_draw["Normals"])
    tooth_draw = tooth_draw.clip_box([-L / 2, L / 2] * 3, invert=False)
    panels = [("phase phi", "Greens", (0, 1), "%.1f"), ("alpha - 1 (1e-3)", "Oranges", (0, atop), "%.2f"),
              (f"von Mises (mPa), film, h = {V['h']:.3f} mm", "viridis", (0, vtop), "%.1f")]
    cut = step == 3
    ncol = 4 if cut else 3
    W = 800
    pl = pv.Plotter(shape=(1, ncol), off_screen=True, window_size=(W * ncol, 912), border=False)
    acts = []
    phi0, am0, _ = phi_at(0)
    keep0 = (phi0 > 0.05) & OPEN
    s0 = _surface(keep0)
    cam = [(-CAM_DIST * 1.6, -CAM_DIST * 0.9, CAM_DIST * 0.75), (0, 0, 0.3), (0, 0, 1)]
    for j, (lab, cm, clim, fmt) in enumerate(panels):
        pl.subplot(0, j)
        pl.set_background("white")
        a = pl.add_mesh(_colour(s0.copy(), f_field(phi0)), scalars="v", cmap=cm, clim=clim, ambient=0.45, diffuse=0.6,
                        specular=0.1, smooth_shading=True,
                        scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                             width=0.7, height=0.07, label_font_size=22, fmt=fmt, n_labels=5,
                                             font_family="times", color="black"))
        pl.add_mesh(tooth_draw, color="#d8d3c8", ambient=0.35, diffuse=0.7, specular=0.15, smooth_shading=True)
        pl.add_text(lab, position="upper_edge", font_size=13, font="times", color="black")
        pl.add_mesh(pv.Plane(center=(0, 0, -8.01), direction=(0, 0, 1), i_size=24, j_size=24), color="#c8c8c8")
        pl.camera_position = cam
        pl.camera.zoom(1.5)
        acts.append(a)
    if cut:
        # cut-away close-up: the film cells (phi > 0.3) of the approximal sector around the distal contact, clipped by
        # a vertical plane through the crown axis along the mesio-distal direction, coloured by phi
        md = G["md"]; cx, cy = G["crown_centre"]
        pl.subplot(0, 3)
        pl.set_background("white")
        img = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-L / 2,) * 3)
        normal = np.array([md[1], -md[0], 0.0])
        cen_cut = np.array([cx + 4.5 * md[0], cy + 4.5 * md[1], Z_GM + 0.8])
        v_phi = pl.add_mesh(pv.Sphere(radius=0.01, center=cen_cut), opacity=0.0)   # placeholder, replaced per frame
        pl.add_mesh(tooth_draw.clip(normal=normal, origin=cen_cut, invert=False), color="#d8d3c8", ambient=0.35,
                    diffuse=0.7, specular=0.15, smooth_shading=True)
        pl.add_text("cut-away: phi in the film, distal", position="upper_edge", font_size=13, font="times", color="black")
        d_cam = 5.5
        pl.camera_position = [tuple(cen_cut + normal * d_cam + np.array([0, 0, 1.5])), tuple(cen_cut), (0, 0, 1)]
        pl.camera.zoom(1.0)
        cut_state = {"actor": None, "img": img, "normal": normal, "cen": cen_cut}
    pl.subplot(0, 0)
    title = pl.add_text("", position="upper_left", font_size=15, font="times", color="black")
    pl.subplot(0, 2)
    note = (f"molar FDI 46 (grey) at real size, cube {L:g} mm; step {step}: {STEPS[step]}; beta = {info['beta mm^2/T*']:.3f} mm^2/T*; "
            f"field 192^3 (h = {H:.3f} mm); undeformed; film surface = phi > 0.05, colour = maximum of the field through the film (0.8 mm)\n"
            f"von Mises: film submodel, h = {V['h']:.3f} mm (3 elements across the 0.25 mm band, 6 across the approximal film), rigid tooth, "
            f"void dropped, solved every 6th frame, linear in time between; range 0 to {vtop:.1f} mPa = 95th percentile on the last frame")
    pl.add_text(note, position="lower_right", font_size=7, font="times", color="#555555")
    pl.link_views([0, 1, 2])

    def frame(i, azim):
        phi, am, tt = phi_at(i)
        keep = (phi > 0.05) & OPEN
        s = _surface(keep)
        for a, f in zip(acts, (f_field(phi), f_field(am * 1e3), f_crop(vm_at(i) * 1e3))):
            a.mapper.SetInputData(_colour(s.copy(), f))
        if cut:
            pl.subplot(0, 3)
            img = cut_state["img"]
            img.cell_data["phi"] = phi.ravel(order="F")
            sel = (phi > FILM_PHI) & OPEN
            img.cell_data["sel"] = sel.ravel(order="F").astype(float)
            vox = img.threshold(0.5, scalars="sel").clip(normal=cut_state["normal"], origin=cut_state["cen"], invert=False)
            vox = vox.clip_box([cut_state["cen"][0] - 4, cut_state["cen"][0] + 4, cut_state["cen"][1] - 4,
                                cut_state["cen"][1] + 4, Z_GM - 0.3, Z_GM + 4.5], invert=False)
            if cut_state["actor"] is not None:
                pl.remove_actor(cut_state["actor"], render=False)
            cut_state["actor"] = pl.add_mesh(vox, scalars="phi", cmap="Greens", clim=(0, 1), show_scalar_bar=False,
                                             ambient=0.4, diffuse=0.6)
        title.SetText(2, f"T* = {tt:.2f}")
        for j in range(3):
            pl.subplot(0, j)
            pl.camera_position = cam
            pl.camera.zoom(1.5)
            pl.camera.azimuth = azim

    out = ROOT / "assets" / f"anim_real_tooth_plaque_{step}_3d"
    nfr = 3 * (n - 1) + 1
    t0 = time.time()
    if png_only:
        frame(nfr - 1, 35.0)
        pl.screenshot(str(out.with_suffix(".png")))
    else:
        pl.open_movie(str(out.with_suffix(".mp4")), framerate=16, quality=9)
        for i in range(nfr):                                   # growth with a slow turn of 180 degrees
            frame(i, 180.0 * i / (nfr - 1) - 90.0)
            pl.write_frame()
        for k in range(96):                                    # a full turn at the final state
            frame(nfr - 1, 90.0 + 360.0 * (k + 1) / 96)
            pl.write_frame()
        frame(nfr - 1, 35.0)
        pl.screenshot(str(out.with_suffix(".png")))
    pl.close()
    sec = time.time() - t0
    merge_json({"steps": {str(step): {"video": {"file": f"assets/anim_real_tooth_plaque_{step}_3d.mp4",
                                                 "frames": nfr + 96, "fps": 16, "colour range mPa": [0.0, vtop],
                                                 "stress frames": V["frames"], "h mm": V["h"],
                                                 "series film vM mean mPa": [1e3 * (V["measures"][k].get("film vM mean") or 0) for k in V["frames"]],
                                                 "seconds render": sec}}}})
    print(f"[step {step}] wrote {out.with_suffix('.png' if png_only else '.mp4')} ({sec:.0f} s)", flush=True)


# ---------------------------------------------------------------- figure
def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply(size=9)
    G = geometry()
    out = json.loads(OUT_JSON.read_text())
    ax_ = axes()
    fig = plt.figure(figsize=(6.3, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], wspace=0.35, hspace=0.5, left=0.09, right=0.98, top=0.95,
                          bottom=0.08)
    # (a) unrolled thickness map of step 2: azimuth about the crown axis (0 = distal) against z, from the first
    # surface cell of each column
    D2 = load_step(2) if cache_path(2).exists() else load_step(1)
    tt = D2["t"].astype(float)
    sdf = G["sdf"]
    surf = (sdf >= 0) & (sdf < H) & (tt > 0)
    X, Y, Z = np.meshgrid(ax_, ax_, ax_, indexing="ij")
    md = G["md"]; cx, cy = G["crown_centre"]
    ang = np.degrees(np.arctan2(Y[surf] - cy, X[surf] - cx) - np.arctan2(md[1], md[0]))
    ang = (ang + 180) % 360 - 180
    zz = Z[surf]; tv = tt[surf]
    nb_a, nb_z = 72, 40
    ea = np.linspace(-180, 180, nb_a + 1); ez = np.linspace(Z_GM, Z_GM + 4.0, nb_z + 1)
    Hs, _, _ = np.histogram2d(ang, zz, bins=[ea, ez], weights=tv)
    Hn, _, _ = np.histogram2d(ang, zz, bins=[ea, ez])
    tmap = np.where(Hn > 0, Hs / np.maximum(Hn, 1), np.nan)
    axa = fig.add_subplot(gs[0, 0])
    im = axa.pcolormesh(ea, ez, tmap.T, cmap="YlGnBu", vmin=0, vmax=max(0.6, np.nanmax(tmap)), shading="flat")
    for a_, lab in ((0, "distal"), (180, "mesial"), (-180, "mesial"), (90, "buccal /\nlingual"), (-90, "lingual /\nbuccal")):
        axa.axvline(a_, color="#888888", lw=0.5, ls=":")
    axa.set_xticks([-180, -90, 0, 90, 180])
    axa.set_xticklabels(["mesial", "b / l", "distal", "l / b", "mesial"], fontsize=8)
    axa.set_ylabel("$z$ (mm)")
    axa.set_xlabel("azimuth about the crown axis")
    axa.axhline(Z_CEJ, color="#eb6834", lw=0.7, ls="--"); axa.text(-175, Z_CEJ + 0.05, "CEJ", color="#eb6834", fontsize=7)
    axa.set_title("(a) thickness map $t$ (mm), step 2", fontsize=9, loc="left")
    axa.grid(False)
    cb = fig.colorbar(im, ax=axa, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=7.5)
    # (b) mean von Mises by region
    axb = fig.add_subplot(gs[0, 1])
    regs = [("film", "film"), ("band", "cervical band"), ("approximal", "approximal"), ("buccal-lingual", "buccal / lingual"),
            ("fissure", "fissures"), ("layer", "0.1 mm at the tooth")]
    cases = [("idealised", None)] + [(f"step {s}", s) for s in (1, 2, 3)]
    cols = ["#444444", "#2a78d6", "#1baf7a", "#eb6834"]
    nb = len(cases); wdt = 0.8 / nb
    for ci, (lab, s) in enumerate(cases):
        vals = []
        for key, _ in regs:
            if s is None:
                v = IDEAL.get(f"{key} vM mean")
            else:
                st = out["steps"].get(str(s), {}).get("stress", {})
                lev = "r2" if "r2" in st else ("r1" if "r1" in st else None)
                v = st[lev].get(f"{key} vM mean") if lev else None
            vals.append(np.nan if v is None else v * 1e3)
        axb.bar(np.arange(len(regs)) + (ci - (nb - 1) / 2) * wdt, vals, wdt, color=cols[ci], label=lab)
    axb.set_xticks(np.arange(len(regs))); axb.set_xticklabels([r[1] for r in regs], rotation=30, ha="right", fontsize=7.5)
    axb.set_ylabel("mean von Mises (mPa)")
    axb.legend(fontsize=7, frameon=False)
    axb.set_title("(b) mean stress by region, $h$ = 0.042 mm", fontsize=9, loc="left")
    # (c) section through the film of step 3 with the resolved stress: the vertical plane through the crown axis along
    # the mesio-distal direction (mesial and distal films and the fissure deposits in one picture)
    axc = fig.add_subplot(gs[1, :])
    lev = "r2" if stress_path(3, "r2").exists() else "r1"
    d = np.load(stress_path(3, lev))
    vm, act, te = d["vm"] * 1e3, d["active"], d["tooth"]
    org, h, shp = d["origin"], float(d["h"]), d["shape"]
    axes3 = [org[k] + (np.arange(shp[k]) + 0.5) * h for k in range(3)]
    # sample the plane: points s along md (horizontal) and z
    s_ax = np.arange(-7.0, 7.0, h); z_ax = np.arange(Z_GM - 0.5, 5.6, h)
    S, ZZ = np.meshgrid(s_ax, z_ax, indexing="ij")
    px = cx + S * md[0]; py = cy + S * md[1]
    pts = np.stack([px.ravel(), py.ravel(), ZZ.ravel()], 1)
    fi = RegularGridInterpolator(tuple(axes3), np.where(act, vm, np.nan), method="nearest", bounds_error=False, fill_value=np.nan)
    ft = RegularGridInterpolator(tuple(axes3), te.astype(float), method="nearest", bounds_error=False, fill_value=0.0)
    img = fi(pts).reshape(S.shape); tm = ft(pts).reshape(S.shape) > 0.5
    # tooth from the field sdf for the part outside the crop box
    fs = RegularGridInterpolator((ax_,) * 3, G["sdf"].astype(float), bounds_error=False, fill_value=1.0)
    tm |= fs(pts).reshape(S.shape) < 0
    vmax = out["steps"]["3"]["stress"][lev]["vM p99"] * 1e3
    ext = [s_ax[0], s_ax[-1], z_ax[0], z_ax[-1]]
    axc.imshow(np.ma.masked_where(~tm.T, np.ones_like(img.T)), cmap=ListedColormap(["#d0ccc4"]), extent=ext, origin="lower",
               vmin=0, vmax=1, interpolation="nearest", aspect="equal")
    imc = axc.imshow(np.ma.masked_invalid(img.T), cmap="viridis", extent=ext, origin="lower", vmin=0, vmax=vmax,
                     interpolation="nearest", aspect="equal")
    axc.axhline(Z_GM, color="#eb6834", lw=0.7, ls="--"); axc.text(s_ax[0] + 0.1, Z_GM + 0.08, "gingival margin", color="#eb6834", fontsize=7)
    axc.set_xlabel("distance along the mesio-distal axis (mm), mesial left")
    axc.set_ylabel("$z$ (mm)")
    axc.set_title(f"(c) step 3, von Mises in the vertical section along the mesio-distal axis, $h$ = {h:.3f} mm", fontsize=9, loc="left")
    axc.grid(False)
    cb = fig.colorbar(imc, ax=axc, fraction=0.025, pad=0.02)
    cb.set_label("von Mises (mPa), scale: 99th percentile", fontsize=8); cb.ax.tick_params(labelsize=7.5)
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG.with_suffix(".pdf"))
    fig.savefig(FIG.with_suffix(".png"), dpi=300)
    print("wrote", FIG.with_suffix(".pdf"), FIG.with_suffix(".png"))


# ---------------------------------------------------------------- driver
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--geom", action="store_true")
    ap.add_argument("--map", nargs="*", type=int, help="print the map diagnostics of the steps (no run)")
    ap.add_argument("--field", nargs="*", type=int)
    ap.add_argument("--stress", nargs="*", type=int)
    ap.add_argument("--levels", nargs="+", default=["r2"], choices=list(LEVELS))
    ap.add_argument("--video-stress", nargs="*", type=int)
    ap.add_argument("--render", nargs="*", type=int)
    ap.add_argument("--png-only", action="store_true")
    ap.add_argument("--vtop", type=float, help="top of the von Mises colour scale, mPa (default: p95 of the last frame)")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    if a.geom:
        geometry(force=a.force)
        merge_json({})
    if a.map:
        G = geometry()
        for s in a.map:
            _, _, info = build_phi0(G, s)
            print(json.dumps(_js(info), indent=1))
            merge_json({"steps": {str(s): {"map": info}}})
    for s in a.field or []:
        run_field(s, force=a.force)
    for s in a.stress or []:
        run_stress(s, a.levels, force=a.force)
    for s in a.video_stress or []:
        run_video_stress(s, force=a.force)
    for s in a.render or []:
        render(s, png_only=a.png_only, vtop=a.vtop)
    if a.figure:
        figure()


if __name__ == "__main__":
    main()
