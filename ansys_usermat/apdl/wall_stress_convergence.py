#!/usr/bin/env python3
"""Mesh convergence of the stress where a growing seed meets the clamped face (3D demo cases).

The demonstrations of anim_seed_offcentre_3d.py (corner and edge seed, bottom face z = -1 clamped) show the largest
von Mises stress at the clamped edge under the seed, about 20 times the value of the centre seed. The corner of the
elastic problem (growing stiff seed meets a rigid face) need not have a finite stress, so this script separates the
measures that converge with the mesh from those that do not, as the cervical study of p3/uel does (band mean against
band maximum).

Field: d phi/dt = beta lap phi + k_alpha alpha - penalty, d alpha/dt = k_alpha phi (wp3_reference.rhs, zero-flux
walls, beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, RK4) on ONE 96^3 grid per case to T* = 1.1, cached in
figdata/wall_field_<case>_96.npz. 96 is a multiple of 16, 24, 32 and 48, so phi and alpha - 1 are block means of the
same field on every stress mesh and only the mechanics changes between the meshes (the 48^3 field of the demos
does not divide 32^3).
Mechanics: hex8, B-bar (mesh_study_seed.element_matrices), E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain
(alpha - 1) I, small strain; bottom face clamped (corner, edge) or the three corner-node rigid-body constraints of the
deck (centre, the reference of Chapter 4). Direct solve up to 24^3, Jacobi-preconditioned CG above (residual in the
JSON).
Measures per mesh (Pa; p = -tr(sigma)/3, compression positive; biofilm = elements with phi > 0.5):
  seed mean vM (initial seed region); band = biofilm with z in [-1, -0.75] (the quarter mm on the clamped face);
  core = biofilm with z in [-0.75, -0.5]; ratio band/core; wall band = biofilm within 0.25 mm of a free wall the
  seed touches, z > -0.75 (the free re-entrant edges); maximum vM and its element centre; 95th and 99th percentile
  over the biofilm; maximum p over the biofilm. Fits: v = v0 + c h^p for the means (observed order, h -> 0 value),
  max = a h^-s for the maxima and percentiles.

    python ansys_usermat/apdl/wall_stress_convergence.py --cases corner edge centre      # solve (one case per
    python ansys_usermat/apdl/wall_stress_convergence.py --figure                        #   process is fine)
        -> ansys_usermat/apdl/results/2026-10-3d/wall_convergence.json, assets/fig_wall_convergence.pdf/.png
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import mesh_study_seed as S  # noqa: E402
import wp3_reference as W  # noqa: E402

E_BIO, NU, FLOOR = 10.0, 0.49, 1e-3
NF = 96
OUT_JSON = HERE / "results/2026-10-3d/wall_convergence.json"
FIG = ROOT / "assets/fig_wall_convergence"
CASES = ("corner", "edge", "centre")


def seed_mask8(case):
    m = np.zeros((8, 8, 8), bool)
    if case == "corner":
        m[0:4, 0:4, 0:2] = True                     # x < 0, y < 0, z < -0.5: touches three faces
    elif case == "edge":
        m[0:4, 2:6, 0:2] = True                     # x < 0, |y| < 0.5, z < -0.5: touches two faces
    elif case == "centre":
        cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
        m[cells8[:, 0], cells8[:, 1], cells8[:, 2]] = True
    else:
        raise ValueError(case)
    assert m.sum() == 32
    return m


def up(m, n):
    f = n // m.shape[0]
    return np.repeat(np.repeat(np.repeat(m, f, 0), f, 1), f, 2)


def field(case, nf=NF, lam=0.15):
    """phi and alpha - 1 at T* on the nf^3 grid (cached)."""
    cache = ROOT / "figdata" / f"wall_field_{case}_{nf}.npz"
    if cache.exists():
        d = np.load(cache)
        return d["phi"], d["am1"]
    h = W.L / nf
    dt = lam * h * h / W.BETA
    nstep = int(np.ceil(W.T_END / dt))
    dt = W.T_END / nstep
    phi = up(seed_mask8(case), nf).astype(float)
    alpha = np.ones_like(phi)
    lap = np.empty_like(phi)
    t0 = time.time()
    for s in range(nstep):
        k1p, k1a = W.rhs(phi, alpha, h, lap)
        k2p, k2a = W.rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a, h, lap)
        k3p, k3a = W.rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a, h, lap)
        k4p, k4a = W.rhs(phi + dt * k3p, alpha + dt * k3a, h, lap)
        phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        if (s + 1) % max(1, nstep // 5) == 0:
            print(f"  field {case} {nf}^3 step {s + 1}/{nstep} phi in [{phi.min():.4f}, {phi.max():.4f}] "
                  f"{time.time() - t0:.0f} s", flush=True)
    cache.parent.mkdir(exist_ok=True)
    np.savez_compressed(cache, phi=phi, am1=alpha - 1.0, dt=dt, nstep=nstep)
    return phi, alpha - 1.0


def mechanics(n, phi, am1, clamped, rtol=1e-10):
    """von Mises, p = -tr/3 (Pa) per element of the n^3 mesh; CG info."""
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, NU, False)
    idx = np.arange((n + 1) ** 3, dtype=np.int32).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3, dtype=np.int32) for a in range(8)], 1)
    Ee = E_BIO * (np.clip(phi.ravel(), 0.0, 1.0) ** 2 + FLOOR)
    eps0 = am1.ravel()
    ndof = 3 * (n + 1) ** 3
    rows = np.repeat(dofs, 24, 1).ravel()
    cols = np.tile(dofs, (1, 24)).ravel()
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    del rows, cols
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    if clamped:
        fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
    else:
        node = lambda p: idx[tuple(int(round((c + 1) / h)) for c in p)]
        fix = [3 * node(S.CORNERS[0]) + d for d in (0, 1, 2)] + [3 * node(S.CORNERS[1]) + d for d in (1, 2)] \
            + [3 * node(S.CORNERS[2]) + 2]
    free = np.setdiff1d(np.arange(ndof), fix)
    Kff = K[free][:, free].tocsc()
    del K
    u = np.zeros(ndof)
    info = {"solver": "direct"}
    if n <= 24:
        u[free] = spla.spsolve(Kff, F[free])
    else:
        dinv = 1.0 / Kff.diagonal()
        Mop = spla.LinearOperator(Kff.shape, matvec=lambda x: dinv * x)
        it = [0]
        x, flag = spla.cg(Kff, F[free], M=Mop, rtol=rtol, maxiter=200000, callback=lambda _: it.__setitem__(0, it[0] + 1))
        res = float(np.linalg.norm(Kff @ x - F[free]) / np.linalg.norm(F[free]))
        info = {"solver": "cg jacobi", "iterations": it[0], "flag": int(flag), "relative residual": res}
        print(f"    cg: {it[0]} iterations, flag {flag}, relative residual {res:.1e}", flush=True)
        u[free] = x
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = -sig[:, :3].mean(1)
    s = sig[:, :3] + p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape((n,) * 3), p.reshape((n,) * 3), info


def measures(case, n, phi, am1, vm, p):
    h = 2.0 / n
    ax = -1 + (np.arange(n) + 0.5) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    seed = up(seed_mask8(case), n)
    bio = phi > 0.5
    bio3 = phi > 0.3                                   # the colour scale of the demos
    band = bio & (Z < -0.75)
    core = bio & (Z > -0.75) & (Z < -0.5)
    walls = {"corner": (X < -0.75) | (Y < -0.75), "edge": X < -0.75, "centre": np.zeros_like(bio)}[case]
    wall = bio & walls & (Z > -0.75)
    imax = np.unravel_index(np.argmax(vm), vm.shape)
    ipmax = np.unravel_index(np.argmax(np.where(bio, p, -np.inf)), p.shape)
    r = {"h": h, "elements": n ** 3,
         "seed vM mean": float(vm[seed].mean()), "seed p mean": float(p[seed].mean()),
         "biofilm elements": int(bio.sum()), "biofilm vM mean": float(vm[bio].mean()),
         "vM max": float(vm.max()), "vM max at": [float(X[imax]), float(Y[imax]), float(Z[imax])],
         "vM max in biofilm": bool(bio[imax]),
         "vM p95": float(np.percentile(vm[bio], 95)), "vM p99": float(np.percentile(vm[bio], 99)),
         "vM p95 phi>0.3": float(np.percentile(vm[bio3], 95)),
         "p max": float(p[bio].max()), "p max at": [float(X[ipmax]), float(Y[ipmax]), float(Z[ipmax])],
         "p p95": float(np.percentile(p[bio], 95))}
    if case != "centre":
        r.update({"band elements": int(band.sum()), "band vM mean": float(vm[band].mean()),
                  "band vM max": float(vm[band].max()), "band p mean": float(p[band].mean()),
                  "core elements": int(core.sum()), "core vM mean": float(vm[core].mean()),
                  "core vM max": float(vm[core].max()), "band/core": float(vm[band].mean() / vm[core].mean()),
                  "wall elements": int(wall.sum()), "wall vM mean": float(vm[wall].mean()),
                  "wall vM max": float(vm[wall].max())})
    return r


def fit_converging(hs, v):
    """Convergence of a mean: change per refinement, first- and second-order extrapolations (linear least squares
    on all meshes), and a free fit v = v0 + c h^p (order bounded to [0.3, 4]). The free fit is marked reliable only
    when the changes are monotone and the last one is above the solver's resolution (0.2 %); otherwise the first-
    order value is the h -> 0 estimate (the seed stress converges with order about 1, wp3_stress_reference)."""
    from scipy.optimize import least_squares
    hs, v = np.asarray(hs, float), np.asarray(v, float)
    change = [float((v[i + 1] / v[i] - 1) * 100) for i in range(len(v) - 1)]
    ext = {}
    for pp in (1, 2):
        A = np.stack([np.ones_like(hs), hs ** pp], 1)
        ext[pp] = float(np.linalg.lstsq(A, v, rcond=None)[0][0])
    f = lambda q: q[0] + q[1] * hs ** q[2] - v
    best = None
    for p0 in (1.0, 2.0):
        q = least_squares(f, [v[-1], v[0] - v[-1], p0], bounds=([-np.inf, -np.inf, 0.3], [np.inf, np.inf, 4.0]))
        if best is None or q.cost < best.cost:
            best = q
    order = float(best.x[2])
    mono = all(c > 0 for c in change) or all(c < 0 for c in change)
    reliable = bool(mono and abs(change[-1]) > 0.2 and 0.31 < order < 3.99)
    h0 = float(best.x[0]) if reliable else ext[1]
    return {"change % per refinement": change, "h->0 p=1": ext[1], "h->0 p=2": ext[2],
            "free fit order": order, "free fit h->0": float(best.x[0]), "free fit reliable": reliable,
            "h->0": h0, "error % at h->0": [float((x / h0 - 1) * 100) for x in v]}


def fit_power(hs, v):
    """v = a h^-s (log-log least squares) and the slope of each successive pair."""
    hs, v = np.asarray(hs, float), np.asarray(v, float)
    s, loga = np.polyfit(np.log(1 / hs), np.log(v), 1)
    pair = [float(np.log(v[i + 1] / v[i]) / np.log(hs[i] / hs[i + 1])) for i in range(len(v) - 1)]
    return {"s": float(s), "a": float(np.exp(loga)), "pair slopes": pair,
            "change % per refinement": [float((v[i + 1] / v[i] - 1) * 100) for i in range(len(v) - 1)]}


def load_json():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"model": {"beta mm^2/T*": W.BETA, "k_alpha 1/T*": W.K_ALPHA, "penalty": W.PENALTY, "T*": W.T_END,
                      "E_bio Pa": E_BIO, "nu": NU, "floor": FLOOR, "field grid": NF,
                      "biofilm": "phi > 0.5", "band": "biofilm, z in [-1, -0.75]",
                      "core": "biofilm, z in [-0.75, -0.5]",
                      "wall band": "biofilm within 0.25 mm of a free wall the seed touches, z > -0.75",
                      "p": "-tr(sigma)/3, compression positive"},
            "cases": {}}


def run_case(case, meshes):
    print(f"== {case}", flush=True)
    phi_f, am1_f = field(case)
    rows = {}
    for n in meshes:
        t0 = time.time()
        phi, am1 = W.block_mean(phi_f, n), W.block_mean(am1_f, n)
        vm, p, info = mechanics(n, phi, am1, clamped=(case != "centre"))
        r = measures(case, n, phi, am1, vm, p)
        r["solver"] = info
        r["seconds"] = time.time() - t0
        rows[str(n)] = r
        extra = (f"  band {r['band vM mean']:.3e} core {r['core vM mean']:.3e} ratio {r['band/core']:.2f}"
                 f" wall max {r['wall vM max']:.3e}") if case != "centre" else ""
        print(f"  {n:>2}^3 seed {r['seed vM mean']:.4e}{extra}  max {r['vM max']:.3e} at {r['vM max at']}"
              f"  p95 {r['vM p95']:.3e} p99 {r['vM p99']:.3e} pmax {r['p max']:.3e}  ({r['seconds']:.0f} s)",
              flush=True)
        np.savez_compressed(ROOT / "figdata" / f"wall_stress_{case}_{n}.npz", vm=vm.astype(np.float32),
                            p=p.astype(np.float32), phi=phi.astype(np.float32))
    fits = compute_fits(case, rows)
    out = load_json()                                  # merge (the cases may run in parallel processes)
    out["cases"][case] = {"meshes": rows, "fits": fits}
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1))
    print("wrote", OUT_JSON, flush=True)


def compute_fits(case, rows):
    ns = sorted(int(k) for k in rows)
    hs = [2.0 / n for n in ns]
    fits = {}
    conv_keys = ["seed vM mean", "seed p mean", "biofilm vM mean"] + \
        (["band vM mean", "core vM mean", "band/core", "wall vM mean", "band p mean"] if case != "centre" else [])
    for k in conv_keys:
        fits[k] = fit_converging(hs, [rows[str(n)][k] for n in ns])
    for k in ["vM max", "vM p95", "vM p99", "p max"] + (["band vM max", "wall vM max", "core vM max"] if case != "centre" else []):
        fits[k] = fit_power(hs, [rows[str(n)][k] for n in ns])
    print(f"== {case}: fits over {ns}")
    for k, f in fits.items():
        if "h->0" in f:
            tag = f"order {f['free fit order']:.2f}" if f["free fit reliable"] else "order n/a (p = 1 used)"
            print(f"  {k:16s} change {' / '.join(f'{e:+.2f}' for e in f['change % per refinement'])} %  {tag}"
                  f"  h->0 {f['h->0']:.4e}  error % {' / '.join(f'{e:+.1f}' for e in f['error % at h->0'])}")
        else:
            print(f"  {k:16s} ~ h^({-f['s']:+.2f}), pair slopes {' / '.join(f'{e:.2f}' for e in f['pair slopes'])},"
                  f" change {' / '.join(f'{e:+.1f}' for e in f['change % per refinement'])} %")
    return fits


def refit():
    out = load_json()
    for case, d in out["cases"].items():
        d["fits"] = compute_fits(case, d["meshes"])
    OUT_JSON.write_text(json.dumps(out, indent=1))
    print("wrote", OUT_JSON)


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    import figstyle
    figstyle.apply(size=9)
    out = load_json()
    C = out["cases"]
    col = {"max": "#2a78d6", "p99": "#6250d6", "p95": "#1baf7a", "band": "#eb6834", "core": "#8c6a00",
           "seed": "#444444"}
    mk = {"max": "o", "p99": "^", "p95": "v", "band": "s", "core": "D", "seed": "x"}

    def axes_h(ax, hs, ns, yt):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_major_locator(FixedLocator(hs)); ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels([f"{h:.2g}\n${n}^3$" for h, n in zip(hs, ns)])
        ax.yaxis.set_major_locator(FixedLocator(yt)); ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(ScalarFormatter()); ax.set_yticklabels([f"{v:g}" for v in yt])
        ax.set_xlim(hs[0] * 1.15, hs[-1] / 1.15)
        ax.set_xlabel("element size $h$ (mm), mesh")

    fig, axs = plt.subplots(1, 3, figsize=(6.3, 2.9))
    series = (("vM max", "maximum", "max"), ("vM p99", "99th percentile", "p99"),
              ("vM p95", "95th percentile", "p95"), ("band vM mean", "band mean", "band"),
              ("core vM mean", "core mean", "core"), ("seed vM mean", "seed mean", "seed"))
    for k, (ax, case) in enumerate(zip(axs[:2], ("corner", "edge"))):
        rows = C[case]["meshes"]
        ns = sorted(int(x) for x in rows)
        hs = [2.0 / n for n in ns]
        for key, lab, c in series:
            v = np.array([rows[str(n)][key] for n in ns]) * 1e3
            ax.plot(hs, v, marker=mk[c], color=col[c], lw=1.5, ms=5, mfc="white" if c != "seed" else col[c],
                    label=lab)
        f = C[case]["fits"]
        ax.text(0.97, 0.30, f"maximum $\\sim h^{{-{f['vM max']['s']:.2f}}}$\nband mean "
                            f"{f['band vM mean']['change % per refinement'][-1]:+.1f} % at the last step",
                transform=ax.transAxes, va="top", ha="right", fontsize=7.5)
        axes_h(ax, hs, ns, [1, 2, 5, 10, 20, 50])
        ax.set_ylim(1.0, 40)
        ax.set_title(f"({'ab'[k]}) {case} seed, clamped face", fontsize=9)
    axs[0].set_ylabel("von Mises (mPa)")
    axs[1].set_yticklabels([])
    # (c) ratios: solid = corner, dashed = edge
    ax = axs[2]
    cen = C["centre"]["meshes"]
    for case, ls in (("corner", "-"), ("edge", "--")):
        rows = C[case]["meshes"]
        ns = sorted(int(x) for x in rows)
        hs = [2.0 / n for n in ns]
        ref = np.array([cen[str(n)]["seed vM mean"] for n in ns])
        ax.plot(hs, [rows[str(n)]["vM max"] for n in ns] / ref, ls=ls, marker=mk["max"], color=col["max"], lw=1.5,
                ms=5, mfc="white")
        ax.plot(hs, [rows[str(n)]["band vM mean"] for n in ns] / ref, ls=ls, marker=mk["band"], color=col["band"],
                lw=1.5, ms=5, mfc="white")
        ax.plot(hs, [rows[str(n)]["band/core"] for n in ns], ls=ls, marker=mk["core"], color=col["core"], lw=1.5,
                ms=5, mfc="white")
    axes_h(ax, hs, ns, [2, 5, 10, 20, 50])
    ax.set_ylim(1.0, 70)
    ax.set_ylabel("ratio")
    ax.set_title("(c) ratios (solid corner, dashed edge)", fontsize=9)
    hnd = [Line2D([], [], color=col["max"], marker=mk["max"], mfc="white", lw=1.5, ms=5, label="maximum / centre seed mean"),
           Line2D([], [], color=col["band"], marker=mk["band"], mfc="white", lw=1.5, ms=5, label="band mean / centre seed mean"),
           Line2D([], [], color=col["core"], marker=mk["core"], mfc="white", lw=1.5, ms=5, label="band mean / core mean")]
    ax.legend(handles=hnd, loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=6.5, frameon=False,
              handlelength=2.2, labelspacing=0.3)
    fig.tight_layout(w_pad=0.8, rect=[0, 0.09, 1, 1])
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=6, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.8, columnspacing=1.2)
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG.with_suffix(".pdf"))
    fig.savefig(FIG.with_suffix(".png"), dpi=300)
    print("wrote", FIG.with_suffix(".pdf"), FIG.with_suffix(".png"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", default=list(CASES), choices=CASES)
    ap.add_argument("--meshes", type=int, nargs="+", default=[16, 24, 32, 48])
    ap.add_argument("--field-only", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--refit", action="store_true", help="recompute the fits from the stored meshes")
    a = ap.parse_args(argv)
    if a.refit:
        refit()
    if a.figure:
        figure()
    if a.refit or a.figure:
        return
    for case in a.cases:
        if a.field_only:
            field(case)
        else:
            run_case(case, a.meshes)


if __name__ == "__main__":
    main()
