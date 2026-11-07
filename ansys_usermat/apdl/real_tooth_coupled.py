"""real_tooth_coupled.py -- the agreed coupling (coupled_living_consumption.py, 20 Oct 2026) on the clinical molar
crown tier2b_real/_crown_real_watertight.stl (derived from Open-Full-Jaw, Gholamalizadeh et al. 2022, CC BY-NC-SA;
THIRD_PARTY.md): the numbers and the film stress of the runs

    python ansys_usermat/apdl/coupled_living_consumption.py --nf 80 --front 1 --phi-min 0.001 --v2 --fortran
        --inherit-living --alpha0 100 --gab 1 --gc 1 --tab 0.333333 --tvd --obstacle stl --box 10 --amg --only tooth_ab
    (the same with --only tooth for the run without antibiotic, --nf 40 for the coarse grid)

Set-up (assumptions of this work, the physics and constants are those of the cube examples, all in mm): a 10 mm box,
the crown (6.7 x 6.8 x 5.5 mm) centred in x and y and standing with its capped cervical base on the bottom face
(the root side, zero flux); nutrient c = 1 and, from T* = 1/3, antibiotic a = 1 on the five other faces (the saliva
side); a cervical ring of biofilm 0.25 mm thick and 1 mm high on the crown surface at the start, species A on the half
x < 0, species B on the half x > 0. The crown is a rigid voxel body for the fields (no biofilm, no nutrient inside,
zero flux through its faces) and a rigid wall (u = 0) for the stress.

Stress: real_tooth_submodel.solve_submodel on the field grid (hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), nu = 0.49,
eigenstrain (alpha - 1) I, rigid crown, elements with phi <= 0.05 dropped so the film surface is free), at the
element size h = h_f (and 2 h_f for the mesh comparison). Film = cells with phi_3D > 0.5 that are not tooth, as the
figures of the cube examples; its means are reported with the share, the living fraction, c, a, alpha - 1 and the
von Mises stress at T* = 1/3 (antibiotic onset) and T* = 1.

    python ansys_usermat/apdl/real_tooth_coupled.py [--nf 80] [--box 10] [--cases tooth_ab,tooth]
        -> ansys_usermat/apdl/results/2026-10-3d/real_tooth_coupled.json (merged per run),
           figdata/real_tooth_coupled_stress_<case>_nf<NF>_T<T>_h<ratio>.npz

The cervical window (7 Nov 2026, papers/coupling_paper/REAL_TOOTH_WINDOW.ja.md; the runs of coupled_living_consumption.py
--obstacle stl --window mixed|layered --stl-scale 100 --box 30 --only window_ab): --window adds the roller faces of
the window to the stress (y = +-BOX/2 and z = +BOX/2: the cut faces, u_n = 0; z = -BOX/2: the gingiva side, u_z = 0;
x = -BOX/2 and the film surface free; the crown rigid) and reports the film thickness (volume over the crown area).
    python ansys_usermat/apdl/real_tooth_coupled.py --window --nf 64 --box 30 --cases window_ab --opts <tag>
        (--opts: the run's tag between nf<NF>_ and _<case>, e.g. front1_..._obststl_winmixed_stlx100_box30_amg)
        -> results/2026-10-3d/real_tooth_coupled_<tag after _tab0.333333_>.json
"""
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import real_tooth_convergence as RC  # noqa: E402
import real_tooth_submodel as RS  # noqa: E402


def _opt(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


WINDOW = "--window" in sys.argv
NF = _opt("--nf", 80, int)
BOX = _opt("--box", 30.0 if WINDOW else 10.0)
CASES = _opt("--cases", "window_ab,window" if WINDOW else "tooth_ab,tooth", str).split(",")
# the roller faces of the window (the mesh faces that reach the domain faces; real_tooth_submodel.solve_submodel)
ROLLER = {"y-": 1, "y+": 1, "z+": 2, "z-": 2} if WINDOW else None
RATIOS = [float(x) for x in _opt("--ratios", "1,2", str).split(",")]
# --opts: the run's option tag between nf<NF>_ and _<case> (7 Nov: the runs with dc/dt and superbee are
# ..._tab0.333333_transc1e+08_tvd_superbee_wgmirror_obststl_box10_amg)
OPTS = _opt("--opts", "front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_tvd_wgmirror_obststl_box10_amg", str)
RES = HERE / "results" / "2026-10-3d"
OUT_JSON = RES / ("real_tooth_coupled.json" if "--opts" not in sys.argv else f"real_tooth_coupled_{OPTS.split('_tab0.333333_')[-1]}.json")
FILM_PHI = 0.5
# the submodel's geometry constants are module globals: the box and a "band" of 1 mm above the cervical base
for mod in (RS, RC):
    mod.L = BOX
    mod.Z_CEJ, mod.BAND_HALF = (-BOX / 4, BOX / 4) if WINDOW else (-BOX / 2 + 0.5, 0.5)   # window: band = the lower half
RS.FILM_PHI = RC.FILM_PHI = FILM_PHI


def npz_path(case):
    return ROOT / "figdata" / f"coupled_living_{case}_nf{NF}_{OPTS}_{case}.npz"


def json_path(case):
    return RES / f"coupled_living_consumption_nf{NF}_{OPTS}_{case}.json"


def film_numbers(z, T, obst, seedA, seedB):
    """means over the film (phi_3D > 0.5, not tooth) and over the two halves of the start ring at the snapshot T"""
    k = lambda f: z[f"{f}_{T}"]  # noqa: E731
    phi, L, sh, c, a = k("phi"), k("L"), k("shareA"), k("c"), k("a")
    am1 = k("alpha") - 1.0
    film = (phi > FILM_PHI) & ~obst
    h = BOX / NF

    def means(m):
        if not m.any():
            return dict(cells=0)
        return dict(cells=int(m.sum()), sum_phibar=float(L[m].mean()), share_A=float(np.nanmean(sh[m])),
                    c=float(c[m].mean()), a=float(a[m].mean()), am1=float(am1[m].mean()), phi=float(phi[m].mean()))
    out = dict(T=T, film=means(film), ringA=means(seedA), ringB=means(seedB),
               film_volume_mm3=float(phi[~obst].sum() * h ** 3), colonised_cells=int((~np.isnan(sh)).sum()),
               open_cells=int((~obst).sum()), living_cells=int((film & (L > 0.01)).sum()),
               c_min_film=float(c[film].min()), am1_max=float(am1[~obst].max()))
    if WINDOW:                                   # film thickness: cells of phi > 0.5 along x per (y, z) column, and
        thick = film.sum(0) * h                  # the volume over the crown area (the voxel faces exposed to open cells)
        area = tooth_area(obst, h)
        out["thickness"] = dict(mean_mm=float(thick.mean()), min_mm=float(thick.min()), max_mm=float(thick.max()),
                                volume_over_area_mm=float(film.sum() * h ** 3 / area), crown_area_mm2=float(area),
                                sdf_max_mm=float(sdf_of(z)[film].max()))
    # film by distance to the crown surface (0.5 mm bins): where the living biomass sits
    sdf = sdf_of(z)
    edges = np.arange(0.0, float(np.ceil(sdf[film].max() / 0.5)) * 0.5 + 0.5, 0.5)
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = film & (sdf >= lo) & (sdf < hi)
        if m.any():
            prof.append(dict(d_lo=float(lo), d_hi=float(hi), cells=int(m.sum()), sum_phibar=float(L[m].mean()),
                             share_A=float(np.nanmean(sh[m])), c=float(c[m].mean()), a=float(a[m].mean())))
    out["by_distance"] = prof
    return out


def sdf_of(z):
    return z[[f for f in z.files if f.startswith("sdf_")][0]]


def tooth_area(obst, h):
    """area of the crown surface as the voxel body exposes it: faces between a tooth cell and an open cell"""
    n = 0
    for ax in range(3):
        for sg in (1, -1):
            m = ~obst & np.roll(obst, -sg, ax)
            e = [slice(None)] * 3; e[ax] = -1 if sg == 1 else 0; m[tuple(e)] = False
            n += int(m.sum())
    return n * h * h


def stress(z, T, obst, sdf, case, ratio):
    phi, am1 = z[f"phi_{T}"], z[f"alpha_{T}"] - 1.0
    h_f = BOX / NF
    t0 = time.time()
    r = RS.solve_submodel(phi, am1, h=ratio * h_f, tooth=obst, sdf=sdf, h_f=h_f, void="drop", void_thr=0.05,
                          rigid=True, tol=1e-9, verbose=True, roller=ROLLER)
    M = r["mesh"]
    np.savez_compressed(ROOT / "figdata" / f"real_tooth_coupled_stress_{case}_nf{NF}_T{T}_h{ratio:g}.npz",
                        vm=r["vm"].astype(np.float32), p=r["p"].astype(np.float32), active=r["active"],
                        tooth=r["tooth"], origin=M.origin, h=M.h, box=np.array(M.box))
    m = r["measures"]
    keep = ("h", "h/h_f", "active elements", "tooth elements", "film elements", "film volume mm^3", "film vM mean",
            "band vM mean", "core vM mean", "band/core", "first layer vM mean", "outer vM mean", "layer vM mean",
            "vM p95", "vM p99", "vM median", "vM max", "vM max at", "p mean", "p max", "void", "void thr", "box mm")
    out = {k: m[k] for k in keep if k in m}
    out.update(solver=r["solver"], seconds=time.time() - t0)
    if WINDOW:                                   # means over the film by distance to the crown (FILM_H0 / 2 bins)
        d_el, vm_e, w = r["sdf_el"], np.nan_to_num(r["vm"]), r["w_film"]
        out["roller"] = ROLLER
        out["vM by distance"] = [dict(d_lo=float(lo), d_hi=float(lo + 1.0), vM_mean=RC.wmean(vm_e, w * ((d_el >= lo) & (d_el < lo + 1.0))))
                                 for lo in np.arange(0.0, float(np.nanmax(np.where(w > 0.5, d_el, 0.0))) + 1.0, 1.0)]
    return out


MODEL_WINDOW = {
    "fields and point model": "coupled_living_consumption.py, the agreed coupling (Eqs. 34-36 with the front term "
                              "-r c/(k+c) n_c . grad phi_3D, TVD superbee, dc/dt kept; point model of Klempt et al. 2026 "
                              "with Fritsch et al. 2025 case I; consumption g phi_3D sum_i phibar_i c; antibiotic from "
                              "T* = 1/3 with alpha*_0 = 100; new points inherit the neighbours' living composition)",
    "geometry": "cervical window of edge 300 um (30 code mm; 1 code mm = 10 um) at the -x flank of the clinical molar "
                "crown (tier2b_real/_crown_real_watertight.stl x 100, Open-Full-Jaw, Gholamalizadeh et al. 2022), "
                "window centre on the crown surface 150 um above the cervical base, the bottom face of the window in the "
                "plane of the cervical base; the crown a rigid voxel body (no biofilm, no nutrient, zero flux)",
    "boundary": "c = 1 (and a = 1 from T* = 1/3) on the saliva face x = -W/2; zero flux through the cut faces "
                "y = +-W/2, z = +W/2, the gingiva face z = -W/2 and the crown",
    "start": "film of 20 um (2 code mm) on the crown surface over the whole window; mixed: Fritsch's start shares "
             "(0.25, 0.30); layered: A on the inner 10 um, B on the outer 10 um",
    "stress": "real_tooth_submodel.solve_submodel: hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain "
              "(alpha - 1) I, rigid crown (u = 0), roller on y = +-W/2 (u_y = 0), z = +W/2 and z = -W/2 (u_z = 0), "
              "x = -W/2 and the film surface free (elements with phi <= 0.05 dropped); film = phi_3D > 0.5",
    "constants": "Klempt 2024 Table 2 at physical size in the code's units (1 code mm = 10 um): beta = 0.02 mm^2/T* "
                 "(2 um^2/T*), r = 10 mm/T* (100 um/T*), k = 1, k_alpha = 1e-3 /T*, penalty 5, g/d = 1 /mm^2 "
                 "(depletion length 1 mm = 10 um; nutrient and antibiotic), c*_0 = 100, alpha*_0 = 100, s = 0.1, "
                 "phi_min = 0.001, outer step 0.025",
}


def main():
    res = json.loads(OUT_JSON.read_text()) if OUT_JSON.exists() else {}
    res["model"] = MODEL_WINDOW if WINDOW else {
        "fields and point model": "coupled_living_consumption.py, the agreed coupling (Eqs. 34-36 with the front term "
                                  "-r c/(k+c) n_c . grad phi_3D, TVD; point model of Klempt et al. 2026 with Fritsch "
                                  "et al. 2025 case I; consumption g phi_3D sum_i phibar_i c; antibiotic from T* = 1/3 "
                                  "with alpha*_0 = 100; new points inherit the neighbours' living composition)",
        "geometry": "clinical molar crown (tier2b_real/_crown_real_watertight.stl, Open-Full-Jaw, Gholamalizadeh et "
                    "al. 2022), rigid voxel body, centred in x, y, cervical base on the bottom face of the box",
        "boundary": "c = 1 (and a = 1 from T* = 1/3) on the five outer faces except the bottom; zero flux at the "
                    "bottom and at the crown",
        "start": "ring of biofilm on the crown surface, thickness and height in the run records; A on x < 0, B on "
                 "x > 0, start amounts as the cube examples (Fritsch's total 0.55, the other species at 1e-6)",
        "stress": "real_tooth_submodel.solve_submodel: hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain "
                  "(alpha - 1) I, rigid crown (u = 0), elements with phi <= 0.05 dropped (free film surface), "
                  "crop box with 1 mm margin; film = phi_3D > 0.5",
        "constants": "as the cube examples in mm: beta = 0.02 mm^2/T*, r = 10 mm/T*, k = 1, k_alpha = 1e-3 /T*, "
                     "penalty 5, g/d = 1 /mm^2 (nutrient and antibiotic), c*_0 = 100, alpha*_0 = 100, s = 0.1, "
                     "phi_min = 0.001, outer step 0.025",
    }
    for case in CASES:
        p = npz_path(case)
        if not p.exists():
            print("missing", p)
            continue
        z = np.load(p)
        hist = json.loads(json_path(case).read_text())["runs"][case]
        times = sorted({float(f.split("_", 1)[1]) for f in z.files if f.startswith("phi_")})
        T_end = times[-1]
        obst = z[f"obst_{T_end}"].astype(bool)
        sdf = z[f"sdf_{T_end}"]
        seedA, seedB = z[f"seed_{T_end}"].astype(bool), z[f"edges_{T_end}"].astype(bool)
        key = f"{case}_nf{NF}"
        rec = dict(case=case, nf=NF, h_mm=BOX / NF, box_mm=BOX, antibiotic=case.endswith("_ab"),
                   tooth=hist.get("tooth"), wall_seconds=None, snapshots={}, stress={})
        for T in times:
            rec["snapshots"][str(T)] = film_numbers(z, T, obst, seedA, seedB)
            print(f"[{key}] T* = {T}: film {rec['snapshots'][str(T)]['film']}", flush=True)
        rec["history"] = {k: hist[k] for k in ("t", "living_seed", "share_a_seed", "living_edges", "share_a_edges",
                                                "c_seed", "a_seed", "colonised", "film", "phi_max", "phi_min") if k in hist}
        for T in times:
            for ratio in RATIOS:
                print(f"[{key}] stress at T* = {T}, h = {ratio:g} h_f", flush=True)
                rec["stress"][f"T{T}_h{ratio:g}"] = stress(z, T, obst, sdf, case, ratio)
        res.setdefault("runs", {})[key] = rec
        OUT_JSON.write_text(json.dumps(res, indent=1, default=float) + "\n")
        print(f"[{key}] merged into {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
