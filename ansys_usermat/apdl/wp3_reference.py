#!/usr/bin/env python3
"""Reference solution of the WP3 seed-growth problem (one species, no front
term) to measure the error of the partner element's NEM field solve.

The element solves, with phi at the integration points and explicit Euler,
    d phi / dt  = beta lap(phi) + k_alpha alpha - P (max(0, phi - 1) + min(0, phi)),
    d alpha / dt = k_alpha phi,
on the 2 mm cube with zero-flux walls, phi(0) = 1 in the 32 seed elements of
the partner's 8^3 mesh and 0 elsewhere, alpha(0) = 1, to T* = 1.1
(RUN_WP3_IKMHIWI03.md, block A/C; deck values beta = 0.02 mm^2/T*,
k_alpha = K_LOCAL1 = 1e-3 /T*, P = PENALTY1 = 5). Apart from the penalty,
which is only active where the source pushes phi above 1 inside the seed,
the system is linear.

This script solves the same problem with second-order finite volumes on a
uniform N^3 grid (the seed is a union of 8^3 cells, so it is resolved exactly
on every N that is a multiple of 8) and classical RK4 in time with a step far
below the stability limit, so the time error is negligible. The grids 48^3,
96^3 and 192^3 give the order of the space error and a Richardson
extrapolation, which is the reference. The same measures as
wp3_alpha_profile.py are reported (seed mean of alpha - 1, domain integral,
alpha - 1 along the x axis through the seed centre), and the NEM runs are
compared with the reference: the error at Delta t = 0.0125 and after
extrapolating the first-order time error to Delta t = 0 from the
Delta t = 0.025 and 0.0125 runs.

    python ansys_usermat/apdl/wp3_reference.py --grids 48 96 192 \
        --nem ansys_usermat/apdl/results/2026-10-wp3_fix \
        --json ansys_usermat/apdl/results/2026-10-wp3_fix/wp3_reference.json
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
BETA = 0.02       # mm^2 / T*
K_ALPHA = 1e-3    # 1 / T*
PENALTY = 5.0
T_END = 1.1
L = 2.0           # cube edge, mm, [-1, 1]^3
OFFSETS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)   # as wp3_alpha_profile.OFFSETS


def seed_cells_8(json_8: Path):
    """Index triples (i, j, k) of the seed cells on the partner's 8^3 mesh."""
    r = json.loads(json_8.read_text())
    a = r["all_stress"]
    elem = np.array(a["elem"], int)
    c = np.stack([a["cx"], a["cy"], a["cz"]], 1)
    sel = np.isin(elem, r["seed_BIOFILM1"])
    idx = np.floor((c[sel] + L / 2) / (L / 8)).astype(int)
    assert len(idx) == 32 and idx.min() >= 0 and idx.max() < 8
    return idx


def seed_mask(n, cells8):
    assert n % 8 == 0
    m = np.zeros((8, 8, 8), bool)
    m[cells8[:, 0], cells8[:, 1], cells8[:, 2]] = True
    f = n // 8
    return np.repeat(np.repeat(np.repeat(m, f, 0), f, 1), f, 2)


def laplacian(u, h, out):
    """7-point Laplacian with zero-flux walls (mirror cells), cell centred."""
    out[...] = -6.0 * u
    out[1:, :, :] += u[:-1, :, :]
    out[:-1, :, :] += u[1:, :, :]
    out[0, :, :] += u[0, :, :]
    out[-1, :, :] += u[-1, :, :]
    out[:, 1:, :] += u[:, :-1, :]
    out[:, :-1, :] += u[:, 1:, :]
    out[:, 0, :] += u[:, 0, :]
    out[:, -1, :] += u[:, -1, :]
    out[:, :, 1:] += u[:, :, :-1]
    out[:, :, :-1] += u[:, :, 1:]
    out[:, :, 0] += u[:, :, 0]
    out[:, :, -1] += u[:, :, -1]
    out /= h * h
    return out


def rhs(phi, alpha, h, lap):
    laplacian(phi, h, lap)
    dphi = BETA * lap + K_ALPHA * alpha - PENALTY * (np.maximum(phi - 1.0, 0.0) + np.minimum(phi, 0.0))
    dalpha = K_ALPHA * phi
    return dphi, dalpha


def solve(n, cells8, lam=0.15, penalty=True, verbose=True):
    """RK4 on an n^3 grid; lam = beta dt / h^2 (RK4 is stable up to 2.78/12 = 0.23
    for the 7-point stencil, whose eigenvalues reach 12 beta / h^2). Returns (axes, alpha - 1 field, phi field)."""
    global PENALTY
    p_saved = PENALTY
    if not penalty:
        PENALTY = 0.0
    h = L / n
    dt = lam * h * h / BETA
    nstep = int(np.ceil(T_END / dt))
    dt = T_END / nstep
    phi = seed_mask(n, cells8).astype(float)
    alpha = np.ones_like(phi)
    lap = np.empty_like(phi)
    t0 = time.time()
    for s in range(nstep):
        k1p, k1a = rhs(phi, alpha, h, lap)
        k2p, k2a = rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a, h, lap)
        k3p, k3a = rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a, h, lap)
        k4p, k4a = rhs(phi + dt * k3p, alpha + dt * k3a, h, lap)
        phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        if verbose and (s + 1) % max(1, nstep // 5) == 0:
            print(f"  n={n} step {s+1}/{nstep} t={dt*(s+1):.3f} "
                  f"phi in [{phi.min():.4f}, {phi.max():.4f}] {time.time()-t0:.0f}s", flush=True)
    PENALTY = p_saved
    ax = -L / 2 + (np.arange(n) + 0.5) * h
    return (ax, ax, ax), alpha - 1.0, phi


def block_mean(field, m):
    """Mean of an n^3 cell field over the cells of the coarser m^3 mesh (n % m == 0)."""
    n = field.shape[0]
    f = n // m
    return field.reshape(m, f, m, f, m, f).mean(axis=(1, 3, 5))


def line_values(axes, field):
    f = RegularGridInterpolator(axes, field, bounds_error=False, fill_value=None)
    return [float(v) for v in f(np.array([[d, 0.0, 0.0] for d in OFFSETS]))]


def mesh_axes(m):
    ax = -L / 2 + (np.arange(m) + 0.5) * (L / m)
    return (ax, ax, ax)


def measures(axes, am1, cells8):
    n = len(axes[0])
    h = axes[0][1] - axes[0][0]
    seed = seed_mask(n, cells8)
    return {"n": n, "seed mean": float(am1[seed].mean()),
            "integral": float(am1.sum() * h ** 3), "line": line_values(axes, am1)}


def richardson(m1, m2, m3):
    """Observed order and extrapolation from three grids with ratio 2."""
    out = {}
    for key in ("seed mean", "integral", "line"):
        a, b, c = (np.atleast_1d(np.array(m[key], float)) for m in (m1, m2, m3))
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.log(np.abs((a - b) / (b - c))) / np.log(2.0)
        ext2 = c + (c - b) / 3.0                       # assuming second order
        out[key] = {"order": p.tolist() if p.size > 1 else float(p[0]),
                    "extrapolated_p2": ext2.tolist() if ext2.size > 1 else float(ext2[0])}
    return out


def nem_fields(folder: Path):
    """Element values alpha - 1 of the NEM runs on their m^3 mesh as arrays,
    from wp3c_n{8,16,24}_dt*.json of the folder (8 Oct, with the surface fix)
    and the dt = 0.025 runs of 8^3 and 16^3 in results/2026-10-07_surface_fix
    (same deck, same build)."""
    surf = HERE / "results/2026-10-07_surface_fix"
    files = sorted(folder.glob("wp3c_n*_dt*.json")) + [surf / "ds8_beta002_dt4.json", surf / "ds16_beta002_dt4.json"]
    runs = {}
    for p in files:
        if not p.exists():
            continue
        rec = json.loads(p.read_text())
        if "all_stress" not in rec:
            print("  (skipped, no all_stress:", p.name, ")")
            continue
        mt = re.search(r"_n(\d+)_dt(\d+)|^ds(\d+)_", p.stem)
        if mt.group(1):
            m, dt = int(mt.group(1)), float("0." + mt.group(2)[1:])
        else:
            m, dt = int(mt.group(3)), 0.025
        a = rec["all_stress"]
        c = np.stack([a["cx"], a["cy"], a["cz"]], 1)
        idx = np.floor((c + L / 2) / (L / m)).astype(int)
        assert idx.min() >= 0 and idx.max() < m and len(c) == m ** 3
        field = np.full((m, m, m), np.nan)
        field[idx[:, 0], idx[:, 1], idx[:, 2]] = a["alpha"]
        assert not np.isnan(field).any()
        runs.setdefault(m, {})[dt] = {"field": field, "file": p.name}
    return runs


def compare(ref_coarse, runs, cells8):
    """Error of the NEM element values against the reference averaged over the
    same elements, at each dt and extrapolated to dt -> 0 (first order,
    2 a(dt) - a(2 dt) from 0.0125 and 0.025). Relative errors of the seed mean,
    the domain integral and the line values; L2 and maximum norm of the
    element-wise error relative to the reference."""
    out = {}
    for m in sorted(runs):
        if m not in ref_coarse:
            continue
        r = ref_coarse[m]
        seed = seed_mask(m, cells8)
        axes = mesh_axes(m)
        r_line = np.array(line_values(axes, r))
        row = {"reference": {"seed mean": float(r[seed].mean()), "integral": float(r.sum() * (L / m) ** 3),
                             "line": r_line.tolist()}}
        series = runs[m]
        cands = {f"dt{dt:g}": s["field"] for dt, s in series.items()}
        if 0.0125 in series and 0.025 in series:
            cands["dt->0"] = 2 * series[0.0125]["field"] - series[0.025]["field"]
        for name, f in cands.items():
            e = f - r
            row[name] = {"seed mean": float(f[seed].mean()), "integral": float(f.sum() * (L / m) ** 3),
                         "line": line_values(axes, f),
                         "err seed mean": float(f[seed].mean() / r[seed].mean() - 1),
                         "err integral": float(f.sum() / r.sum() - 1),
                         "err line": (np.array(line_values(axes, f)) / r_line - 1).tolist(),
                         "err L2": float(np.sqrt((e ** 2).sum() / (r ** 2).sum())),
                         "err max": float(np.abs(e).max() / np.abs(r).max()),
                         "err seed L2": float(np.sqrt((e[seed] ** 2).sum() / (r[seed] ** 2).sum())),
                         "err seed max": float(np.abs(e[seed]).max() / np.abs(r[seed]).max())}
        out[m] = row
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--grids", type=int, nargs="+", default=[48, 96, 192])
    ap.add_argument("--meshes", type=int, nargs="+", default=[8, 16, 24], help="NEM meshes to average onto")
    ap.add_argument("--lam", type=float, default=0.15)
    ap.add_argument("--seed-json", type=Path,
                    default=HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    ap.add_argument("--nem", type=Path, default=HERE / "results/2026-10-wp3_fix")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--no-penalty", action="store_true")
    a = ap.parse_args(argv)

    cells8 = seed_cells_8(a.seed_json)
    res, coarse = [], []
    for n in a.grids:
        t0 = time.time()
        axes, am1, phi = solve(n, cells8, a.lam, penalty=not a.no_penalty)
        m = measures(axes, am1, cells8)
        m["phi min"], m["phi max"] = float(phi.min()), float(phi.max())
        m["seconds"] = time.time() - t0
        res.append(m)
        coarse.append({mm: block_mean(am1, mm) for mm in a.meshes if n % mm == 0})
        print(f"{n}^3: seed mean {m['seed mean']:.6e} integral {m['integral']:.6e} "
              f"line {' / '.join(f'{v:.4e}' for v in m['line'])}  ({m['seconds']:.0f} s)", flush=True)
    out = {"beta": BETA, "k_alpha": K_ALPHA, "penalty": 0.0 if a.no_penalty else PENALTY,
           "T": T_END, "lam_rk4": a.lam, "offsets_mm": OFFSETS, "grids": res}
    if len(res) >= 3:
        out["richardson"] = richardson(*res[-3:])
        r = out["richardson"]
        print("observed order: seed mean %.2f, integral %.2f, line %s" % (
            r["seed mean"]["order"], r["integral"]["order"],
            " / ".join(f"{p:.2f}" for p in r["line"]["order"])))
        print("reference (p = 2 extrapolation): seed mean %.6e integral %.6e line %s" % (
            r["seed mean"]["extrapolated_p2"], r["integral"]["extrapolated_p2"],
            " / ".join(f"{v:.4e}" for v in r["line"]["extrapolated_p2"])))
    if len(res) >= 2:
        # reference on the NEM meshes: p = 2 extrapolation of the element means of the two finest grids
        ref_coarse = {mm: coarse[-1][mm] + (coarse[-1][mm] - coarse[-2][mm]) / 3.0
                      for mm in coarse[-1] if mm in coarse[-2]}
        out["reference_on_mesh"] = {str(mm): f.tolist() for mm, f in ref_coarse.items()}
        if a.nem and a.nem.exists():
            runs = nem_fields(a.nem)
            cmp_ = compare(ref_coarse, runs, cells8)
            out["nem_error"] = {str(mm): v for mm, v in cmp_.items()}
            print("\nNEM against the reference averaged over the same elements (relative errors):")
            for mm, row in cmp_.items():
                for name in sorted(k for k in row if k != "reference"):
                    e = row[name]
                    print(f"  {mm:2d}^3 {name:9s} seed mean {e['err seed mean']*100:+6.2f} % integral "
                          f"{e['err integral']*100:+6.2f} % L2 {e['err L2']*100:5.2f} % max {e['err max']*100:5.2f} % "
                          f"seed L2 {e['err seed L2']*100:5.2f} %  line "
                          + " / ".join(f"{v*100:+5.1f}" for v in e["err line"]) + " %")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
        print("wrote", a.json)


if __name__ == "__main__":
    main()
