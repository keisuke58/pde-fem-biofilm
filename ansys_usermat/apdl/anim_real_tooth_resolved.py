"""anim_real_tooth_resolved.py -- the real-molar videos with the von Mises panel from the film submodel.

The three videos of the molar (anim_real_tooth_3d.py: base; real_tooth_model.py: gate, nutrient + gate) colour their
third panel with the stress of the 32^3 full-box mesh (h = 0.5 mm, the 0.5 mm film one element thick), which
REAL_TOOTH_CONVERGENCE.md shows is not converged (film mean 2.0 mPa against 3.3 mPa resolved, 95th percentile 5.3 against
10.5 mPa). Here the stress of every third frame of each run is solved again with the film submodel of
real_tooth_submodel.py (rigid tooth, crop box around the film, void phi <= 0.05 dropped, wall = the field's own voxel
staircase) at h = h_f / 2 = 0.083 mm on the 96^3 field of the run, six elements across the film, where
REAL_TOOTH_SUBMODEL.md finds the film mean within 1 % of the finest meshes. The growth fields (phi, alpha - 1, c) are
those of the run caches and are not recomputed; the gate runs were driven by the 32^3 pressure during the run, so only
their stress panel changes, not their growth. Between the solved frames the stress is linear in time, as in the
originals. The film surface (phi > 0.05) is coloured with the maximum of the von Mises stress along the inward normal
through the film (as the originals; "--reduce mean" gives the mean over the film elements on that line instead). The
colour range is 0 to the 95th percentile (volume-weighted over the film) of the last frame of the base run, the same
for all three videos so that they can be compared.

    python ansys_usermat/apdl/anim_real_tooth_resolved.py --solve [--tags base gate nutrient_gate] [--worker k --workers n]
        -> figdata/real_tooth_resolved_<tag>_f<k>.npz (one file per solved frame, written as soon as it is solved and
           skipped on a rerun, so a restart loses at most one solve; two workers = two solves at a time)
    python ansys_usermat/apdl/anim_real_tooth_resolved.py --render [--tags ...] [--png-only] [--reduce max|mean]
        -> assets/anim_real_tooth_resolved_<tag>_3d.mp4 / .png,
           ansys_usermat/apdl/results/2026-10-3d/real_tooth_videos.json
    python ansys_usermat/apdl/anim_real_tooth_resolved.py            # both, one solve at a time
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")        # threaded BLAS only slows the CG down; two solves share the four cores
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import real_tooth_submodel as SM  # noqa: E402

TAGS = ["base", "gate", "nutrient_gate"]
STRESS_EVERY = 3                           # as the originals
VOID_THR = 0.05
RATIO = 0.5                                # h = h_f / 2 = 0.083 mm on the 96^3 field
OUT_JSON = HERE / "results" / "2026-10-3d" / "real_tooth_videos.json"


def cache_path(tag):
    return SM.RC.CACHE_DEMO if tag == "base" else ROOT / "figdata" / f"real_tooth_{tag}.npz"


def frame_path(tag, k):
    return ROOT / "figdata" / f"real_tooth_resolved_{tag}_f{k:02d}.npz"


def solved_frames(nf):
    s = list(range(0, nf, STRESS_EVERY))
    if s[-1] != nf - 1:
        s.append(nf - 1)
    return s


def _js(o):
    if isinstance(o, (np.floating, float)):
        return float(o)
    if isinstance(o, (np.integer, int)):
        return int(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, dict):
        return {k: _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    return o


# ---------------------------------------------------------------- solve
class Run:
    """A run cache: frames of phi and alpha - 1 (float16, 96^3), the tooth of the base cache, the crop box of the
    final frame (shared by every frame of the run, so the solved stresses live on one grid)."""

    def __init__(self, tag):
        self.tag = tag
        D = np.load(cache_path(tag))
        self.t = D["t"]
        self.PHI, self.AM1 = D["phi"], D["am1"]
        self.nf = len(self.t)
        self.sdf = np.load(SM.RC.CACHE_DEMO)["sdf"].astype(np.float64)
        self.tooth = self.sdf < 0
        self.h_f = SM.L / self.PHI.shape[1]
        self.box = SM.crop_box(self.PHI[-1].astype(np.float64), self.tooth, self.h_f)

    def fields(self, k):
        return self.PHI[k].astype(np.float64), self.AM1[k].astype(np.float64)


def solve_frame(run, k, verbose=True):
    out = frame_path(run.tag, k)
    if out.exists():
        return None
    phi, am1 = run.fields(k)
    h = RATIO * run.h_f
    t0 = time.time()
    if am1.max() <= 0.0:                                     # the first frame: no growth yet, no stress
        M = SM.Mesh(run.box, run.h_f, RATIO)
        tooth_e = M.sample(run.tooth.astype(float)) > 0.5
        active = ~tooth_e & (np.clip(M.sample(phi), 0, 1) > VOID_THR)
        vm = np.where(active, 0.0, np.nan)
        rec = dict(vm=vm, p=vm.copy(), active=active, origin=M.origin, h=M.h, shape=np.array(M.shape),
                   measures="null", solver="null")
    else:
        r = SM.solve_submodel(phi, am1, h, tooth=run.tooth, sdf=run.sdf, h_f=run.h_f, box=run.box, void="drop",
                              void_thr=VOID_THR, tooth_mode="voxel", verbose=False)
        M = r["mesh"]
        rec = dict(vm=r["vm"], p=r["p"], active=r["active"], origin=M.origin, h=M.h, shape=np.array(M.shape),
                   measures=json.dumps(_js(r["measures"])), solver=json.dumps(_js(r["solver"])))
    sec = time.time() - t0
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, vm=rec["vm"].astype(np.float32), p=rec["p"].astype(np.float32), active=rec["active"],
                        origin=rec["origin"], h=rec["h"], shape=rec["shape"], t=float(run.t[k]), frame=k,
                        measures=rec["measures"], solver=rec["solver"], seconds=sec)
    os.replace(tmp, out)                                     # the file appears complete or not at all
    if verbose:
        m = json.loads(rec["measures"])
        s = (f"film {m['film vM mean'] * 1e3:.3f} band/core {m['band/core']:.3f} p95 {m['vM p95'] * 1e3:.3f} "
             f"max {m['vM max'] * 1e3:.2f} mPa") if m else "no growth, stress 0"
        print(f"[{run.tag}] frame {k:2d}/{run.nf - 1}  t = {run.t[k]:.3f}  {s}  {sec:.0f} s", flush=True)
    return sec


def solve(tags, worker=0, workers=1):
    jobs = []
    runs = {}
    for tag in tags:
        if not cache_path(tag).exists():
            print("missing cache", cache_path(tag))
            continue
        runs[tag] = Run(tag)
        jobs += [(tag, k) for k in solved_frames(runs[tag].nf)]
    mine = jobs[worker::workers]
    todo = [j for j in mine if not frame_path(*j).exists()]
    print(f"worker {worker}/{workers}: {len(mine)} frames, {len(todo)} to solve", flush=True)
    t0 = time.time()
    for tag, k in todo:
        solve_frame(runs[tag], k)
    print(f"worker {worker} done in {time.time() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------- render
def load_resolved(tag, nf):
    """The solved frames of a run: von Mises on the crop grid with the inactive cells (tooth, void) filled from the
    nearest active one, the active mask as float, the measures, the seconds; and the crop axes."""
    import real_tooth_model as RM
    frames = solved_frames(nf)
    VM, ACT, meas, secs, axes, h_el = {}, {}, {}, {}, None, None
    for k in frames:
        p = frame_path(tag, k)
        if not p.exists():
            raise SystemExit(f"{p} missing: run --solve first")
        d = np.load(p)
        vm = d["vm"].astype(np.float64)
        act = np.asarray(d["active"], bool) & np.isfinite(vm)
        VM[k] = RM.fill_nearest(np.where(act, vm, np.nan), ~act).astype(np.float32) if act.any() else np.zeros_like(vm, np.float32)
        ACT[k] = act.astype(np.float32)
        meas[k] = json.loads(str(d["measures"]))
        secs[k] = float(d["seconds"])
        if axes is None:
            o, h, shp = d["origin"], float(d["h"]), d["shape"]
            axes = [o[i] + (np.arange(int(shp[i])) + 0.5) * h for i in range(3)]
            h_el = h
    return dict(frames=frames, vm=VM, act=ACT, measures=meas, seconds=secs, axes=axes, h=h_el)


def render(tag, vtop, png_only=False, reduce="max"):
    import pyvista as pv
    from scipy.interpolate import RegularGridInterpolator
    import real_tooth_model as RM

    pv.OFF_SCREEN = True
    pv.global_theme.font.family = "times"
    L, H, NF, ax, OPEN, tooth = RM.L, RM.H, RM.NF, RM.ax, RM.OPEN, RM.tooth
    D = dict(np.load(cache_path(tag)))
    if tag == "base":
        D.update(nutrient=False, gate=False, q=np.nan, p_h=np.nan)
    t, PHI, AM1 = D["t"], D["phi"], D["am1"]
    nf = len(t)
    nutrient, gate = bool(D["nutrient"]), bool(D["gate"])
    R = load_resolved(tag, nf)
    fr = R["frames"]
    atop = float(AM1[-1].astype(float).max()) * 1e3

    def vm_at(i):
        """Von Mises and active weight on the crop grid at frame i, linear in time between the solved frames."""
        b = next(f for f in fr if f >= i)
        if b == i or b == fr[0]:
            return R["vm"][b], R["act"][b]
        a = fr[fr.index(b) - 1]
        w = (i - a) / (b - a)
        return (1 - w) * R["vm"][a] + w * R["vm"][b], (1 - w) * R["act"][a] + w * R["act"][b]

    def surface(keep):
        """Envelope of the cells with phi > 0.05 (smoothed), with point normals (real_tooth_model.film_surface)."""
        g = pv.ImageData(dimensions=(NF + 1,) * 3, spacing=(H,) * 3, origin=(-L / 2,) * 3)
        g.cell_data["keep"] = keep.ravel(order="F").astype(float)
        surf = g.threshold(0.5, scalars="keep").extract_surface(algorithm="dataset_surface").extract_largest()
        if surf.n_points == 0:
            return surf
        surf = surf.smooth_taubin(n_iter=60, pass_band=0.05, feature_angle=150.0)
        return surf.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)

    def colour(surf, f, red, f_w=None):
        """The field sampled along the inward normal over DEPTH mm: maximum, plain mean, or the mean weighted by
        f_w (the active film elements) where f_w is given."""
        if surf.n_points == 0:
            surf.point_data["v"] = np.zeros(0)
            return surf
        pts = np.asarray(surf.points)
        nrm = np.asarray(surf["Normals"])
        depths = np.linspace(0.0, RM.DEPTH, 10)
        samples = np.stack([f(pts - d * nrm) for d in depths], 0)
        if red == "max":
            v = samples.max(0)
        elif f_w is not None:
            w = np.stack([f_w(pts - d * nrm) for d in depths], 0)
            v = (samples * w).sum(0) / np.maximum(w.sum(0), 1e-9)
        else:
            v = samples.mean(0)
        surf.point_data["v"] = v
        return surf

    def f_field(a, fill=0.0):
        return RegularGridInterpolator((ax,) * 3, np.nan_to_num(a), bounds_error=False, fill_value=fill)

    def f_crop(a):
        return RegularGridInterpolator(tuple(R["axes"]), a, bounds_error=False, fill_value=0.0)

    KEEP = lambda i: (PHI[i].astype(float) > 0.05) & OPEN       # noqa: E731

    def vm_panel(i):
        vm, act = vm_at(i)
        return f_crop(vm * 1e3), reduce, (f_crop(act) if reduce == "mean" else None)

    phi_panel = (lambda i: (f_field(PHI[i].astype(float)), "max", None), "phase phi", "Greens", (0, 1), "%.1f")
    am1_panel = (lambda i: (f_field(AM1[i].astype(float) * 1e3), "max", None), "alpha - 1 (1e-3)", "Oranges",
                 (0, atop), "%.2f")
    vm_lab = "von Mises (mPa), film, h = 0.083 mm"
    vm_p = (vm_panel, vm_lab, "viridis", (0, vtop), "%.1f")
    if nutrient:
        C = D["c"]
        c_max = max(float(C[-1].astype(float)[(PHI[-1].astype(float) > 0.5) & OPEN].max()), 0.05)
        ctop = float(np.ceil(c_max * 20) / 20)
        c_panel = (lambda i: (f_field(RM.fill_nearest(C[i].astype(float), tooth)), "mean", None),
                   "nutrient c, film mean", "Blues", (0.0, ctop), "%.2f")
        panels = [c_panel, am1_panel, vm_p]
    else:
        panels = [phi_panel, am1_panel, vm_p]
    red_txt = {"max": "maximum", "mean": "mean"}[reduce]
    note = (f"molar FDI 46 (grey) at real size, cube {L:g} mm; film 0.5 mm, beta = {RM.BETA:.2f} mm^2/T*; undeformed; "
            f"film surface = phi > 0.05, colour = maximum of the field through the film")
    if nutrient:
        note += f"; c = 1 on the box faces around the crown, g/d = {float(D['q']):g} /mm^2, c panel = mean through the film"
    if gate:
        note += (f"; growth gate 1 - p/p_h, p_h = {float(D['p_h']) * 1e3:.2f} mPa (0.5 x the film mean pressure "
                 f"without the gate, 32^3 mesh)")
    note += (f"\nvon Mises: film submodel, h = {R['h']:.3f} mm (six elements across the film), rigid tooth, void dropped, "
             f"solved every third frame and linear in time between; colour = {red_txt} through the film along the "
             f"inward normal; range 0 to {vtop:.1f} mPa = 95th percentile of the base run's last frame, same in all three videos")

    tooth_draw = _tooth_draw(RM, pv)
    pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 912), border=False)
    acts = []
    for j, (get, lab, cm, clim, fmt) in enumerate(panels):
        pl.subplot(0, j)
        pl.set_background("white")
        a = pl.add_mesh(colour(surface(KEEP(0)), *get(0)), scalars="v", cmap=cm, clim=clim, ambient=0.45, diffuse=0.6,
                        specular=0.1, smooth_shading=True,
                        scalar_bar_args=dict(title=" " * (j + 1), vertical=False, position_x=0.15, position_y=0.03,
                                             width=0.7, height=0.07, label_font_size=24, fmt=fmt, n_labels=5,
                                             font_family="times", color="black"))
        pl.add_mesh(tooth_draw, color="#d8d3c8", ambient=0.35, diffuse=0.7, specular=0.15, smooth_shading=True)
        pl.add_text(lab, position="upper_edge", font_size=14, font="times", color="black")
        pl.add_mesh(pv.Box((-L / 2, L / 2) * 3).extract_feature_edges(), color="#9a9a9a", line_width=1.0)
        pl.add_mesh(pv.Plane(center=(0, 0, -L / 2 - 0.01), direction=(0, 0, 1), i_size=1.3 * L, j_size=1.3 * L),
                    color="#b8b8b8")
        pl.camera_position = RM.CAM
        pl.camera.zoom(1.0)
        acts.append((a, get))
    pl.subplot(0, 0)
    title = pl.add_text("", position="upper_left", font_size=16, font="times", color="black")
    pl.subplot(0, 2)
    pl.add_text(note, position="lower_right", font_size=8, font="times", color="#555555")
    pl.link_views()

    def frame(i, azim):
        keep = KEEP(i)
        for a, get in acts:
            a.mapper.SetInputData(colour(surface(keep), *get(i)))
        title.SetText(2, f"T* = {t[i]:.2f}")
        pl.camera.azimuth = azim

    out = ROOT / "assets" / f"anim_real_tooth_resolved_{tag}_3d"
    out.parent.mkdir(exist_ok=True)
    cam0 = pl.camera.azimuth
    t0 = time.time()
    if png_only:
        frame(nf - 1, cam0)
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
    sec = time.time() - t0
    print(f"[{tag}] wrote {out.with_suffix('.png' if png_only else '.mp4')} ({sec:.0f} s)", flush=True)
    return R, sec


def _tooth_draw(RM, pv):
    """The tooth surface drawn half a field cell inside its voxel set, so that it does not poke through the film."""
    td = RM.tooth_mesh.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
    td.points = np.asarray(td.points) - 0.5 * RM.H * np.asarray(td["Normals"])
    return td.clip_box([-RM.L / 2, RM.L / 2] * 3, invert=False)


def colour_top():
    """0 to the 95th percentile of the last frame of the base run (volume-weighted over the film), in mPa."""
    nf = len(np.load(cache_path("base"))["t"])
    p = frame_path("base", nf - 1)
    if not p.exists():
        raise SystemExit(f"{p} missing: solve the base run first (it sets the colour range)")
    return float(json.loads(str(np.load(p)["measures"]))["vM p95"]) * 1e3


def render_all(tags, png_only=False, reduce="max"):
    vtop = colour_top()
    print(f"colour range 0 to {vtop:.3f} mPa", flush=True)
    res = json.loads(OUT_JSON.read_text()) if OUT_JSON.exists() else {}
    res["colour range mPa"] = [0.0, vtop]
    res["colour"] = f"{reduce} of the von Mises stress along the inward normal through the film"
    res.setdefault("runs", {})
    for tag in tags:
        R, sec = render(tag, vtop, png_only, reduce)
        keys = ("film vM mean", "band vM mean", "core vM mean", "band/core", "layer vM mean", "vM p95", "vM p99",
                "vM max", "p mean")
        series = {k: [(R["measures"][f] or {}).get(k) for f in R["frames"]] for k in keys}
        res["runs"][tag] = dict(h_mm=R["h"], frames_solved=R["frames"], t_solved=[float(np.load(cache_path(tag))["t"][f]) for f in R["frames"]],
                                seconds_solve=sum(R["seconds"].values()), seconds_solve_per_frame=R["seconds"],
                                seconds_render=sec, last=R["measures"][R["frames"][-1]], series=series,
                                video=f"assets/anim_real_tooth_resolved_{tag}_3d.mp4")
        OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
        OUT_JSON.write_text(json.dumps(_js(res), indent=1) + "\n")
    print("wrote", OUT_JSON)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--solve", action="store_true")
    ap.add_argument("--render", action="store_true")
    ap.add_argument("--tags", nargs="+", default=TAGS, choices=TAGS)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--png-only", action="store_true")
    ap.add_argument("--reduce", choices=["max", "mean"], default="max")
    a = ap.parse_args(argv)
    if not a.solve and not a.render:
        a.solve = a.render = True
    if a.solve:
        solve(a.tags, a.worker, a.workers)
    if a.render:
        render_all(a.tags, a.png_only, a.reduce)


if __name__ == "__main__":
    main()
