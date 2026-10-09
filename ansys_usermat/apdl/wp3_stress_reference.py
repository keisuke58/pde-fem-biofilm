#!/usr/bin/env python3
"""Mesh error of the seed stress with the reference growth field.

The error of the seed stress on the NEM meshes has two parts: the error of
the field solve (alpha - 1 against the finite-volume reference,
wp3_reference.py) and the error of the mechanical solve on the same mesh.
This script measures the second part. It takes phi and alpha - 1 at T* from
the reference solution, averaged over the elements of the 8^3, 16^3, 24^3 and
48^3 meshes (reference_phi_on_mesh / reference_on_mesh in
wp3_reference.json), and solves the mechanics of the element on each mesh
with the hex8 solver of mesh_study_seed.py: E(phi) = (phi^2 + f) E_bio
(Eq. 20 of Klempt et al. 2024 with the void floor f = 1e-3, as in the runs),
nu = 0.49, B-bar, eigenstrain (alpha - 1) I, small strain, rigid-body
constraints at the three corner nodes of the deck. The seed mean of the von
Mises and the mean stress are reported per mesh with the Richardson
extrapolation of 24^3 and 48^3, which gives the mechanical mesh error of
each NEM mesh for this field. The seed von Mises stresses of the element's
runs (dt -> 0) are listed next to them.

    python ansys_usermat/apdl/wp3_stress_reference.py \
        [--json ansys_usermat/apdl/results/2026-10-wp3_fix/wp3_reference.json]
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
sys.path.insert(0, str(HERE))
import mesh_study_seed as S  # noqa: E402

E_BIO, NU, FLOOR = 10.0, 0.49, 1e-3


def solve_fields(n, phi, am1, nu=NU):
    """Hex8 solve on the n^3 mesh with element-wise phi and alpha - 1 (arrays n x n x n).
    Returns von Mises and mean stress per element (Pa) in the array layout of phi."""
    h = 2.0 / n
    Ke0, fe0, B0, D, m = S.element_matrices(h, nu, False)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    conn = np.stack([idx[i:i + n, j:j + n, k:k + n].ravel() for i, j, k in off], 1)
    dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
    Ee = E_BIO * (np.clip(phi.ravel(), 0.0, 1.0) ** 2 + FLOOR)
    eps0 = am1.ravel()
    rows = np.repeat(dofs, 24, 1).ravel()
    cols = np.tile(dofs, (1, 24)).ravel()
    ndof = 3 * (n + 1) ** 3
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    node = lambda p: idx[tuple(int(round((c + 1) / h)) for c in p)]
    fix = [3 * node(S.CORNERS[0]) + d for d in (0, 1, 2)] + [3 * node(S.CORNERS[1]) + d for d in (1, 2)] \
        + [3 * node(S.CORNERS[2]) + 2]
    free = np.setdiff1d(np.arange(ndof), fix)
    u = np.zeros(ndof)
    Kff = K[free][:, free].tocsc()
    u[free] = spla.spsolve(Kff, F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    p = sig[:, :3].mean(1)
    s = sig[:, :3] - p[:, None]
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape(phi.shape), p.reshape(phi.shape)


def seed_mask(n, cells8):
    m = np.zeros((8, 8, 8), bool)
    m[cells8[:, 0], cells8[:, 1], cells8[:, 2]] = True
    f = n // 8
    return np.repeat(np.repeat(np.repeat(m, f, 0), f, 1), f, 2)


def nem_seed_vm():
    """Seed mean von Mises of the element's runs, dt -> 0 (as wp3_figs.fig_convergence)."""
    sys.path.insert(0, str(HERE))
    from wp3_figs import measures, FIX, SURF
    runs = {8: (SURF / "ds8_beta002_dt4", FIX / "wp3c_n8_dt00125"),
            16: (SURF / "ds16_beta002_dt4", FIX / "wp3c_n16_dt00125"),
            24: (FIX / "wp3c_n24_dt0025", FIX / "wp3c_n24_dt00125")}
    out = {}
    for n, pair in runs.items():
        v = []
        for p in pair:
            if not p.with_suffix(".json").exists():
                break
            m, _ = measures(json.loads(p.with_suffix(".json").read_text()))
            v.append((m["seed vM"], m["seed p"]))
        if len(v) == 2:
            out[n] = {"vM dt0.0125": v[1][0], "vM dt->0": 2 * v[1][0] - v[0][0],
                      "p dt0.0125": v[1][1], "p dt->0": 2 * v[1][1] - v[0][1]}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=HERE / "results/2026-10-wp3_fix/wp3_reference.json")
    ap.add_argument("--out", type=Path, default=HERE / "results/2026-10-wp3_fix/wp3_stress_reference.json")
    ap.add_argument("--meshes", type=int, nargs="+", default=[8, 16, 24, 48])
    ap.add_argument("--seed-json", type=Path, default=HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    a = ap.parse_args(argv)
    ref = json.loads(a.json.read_text())
    sys.path.insert(0, str(HERE))
    from wp3_reference import seed_cells_8
    cells8 = seed_cells_8(a.seed_json)
    rows = {}
    for n in a.meshes:
        if str(n) not in ref["reference_on_mesh"]:
            print("no reference on", n); continue
        phi = np.array(ref["reference_phi_on_mesh"][str(n)])
        am1 = np.array(ref["reference_on_mesh"][str(n)])
        t0 = time.time()
        vm, p = solve_fields(n, phi, am1)
        seed = seed_mask(n, cells8)
        rows[n] = {"seed vM": float(vm[seed].mean()), "seed vM max": float(vm[seed].max()),
                   "seed p": float(p[seed].mean()), "outside vM max": float(vm[~seed].max()),
                   "seconds": time.time() - t0}
        r = rows[n]
        print(f"{n:>3}^3  seed vM {r['seed vM']:.4e}  seed p {r['seed p']:+.4e}  "
              f"outside vM max {r['outside vM max']:.3e}  ({r['seconds']:.0f} s)", flush=True)
    out = {"E_bio": E_BIO, "nu": NU, "floor": FLOOR, "python": {str(n): v for n, v in rows.items()}}
    ns = sorted(rows)
    if len(ns) >= 3:
        for key in ("seed vM", "seed p"):
            v = np.array([rows[n][key] for n in ns])
            hs = np.array([2.0 / n for n in ns])
            # order from the last three meshes (unequal ratios): fit v = v0 + c h^p
            from scipy.optimize import least_squares
            fun = lambda q: q[0] + q[1] * hs[-3:] ** q[2] - v[-3:]
            q = least_squares(fun, [v[-1], v[0] - v[-1], 2.0]).x
            ext = {"order": float(q[2]), "h->0": float(q[0]),
                   "error %": {str(n): float((rows[n][key] / q[0] - 1) * 100) for n in ns}}
            out[key + " extrapolation"] = ext
            print(f"{key}: order {q[2]:.2f}, h -> 0: {q[0]:.4e}; error " +
                  ", ".join(f"{n}^3 {e:+.1f} %" for n, e in ext["error %"].items()))
    # the element's own alpha - 1 (dt -> 0) with the reference phi: effect of the field error on the stress
    from wp3_reference import nem_fields
    runs = nem_fields(HERE / "results/2026-10-wp3_fix")
    out["python with NEM alpha"] = {}
    for n in a.meshes:
        if n in runs and 0.0125 in runs[n] and 0.025 in runs[n] and str(n) in ref["reference_on_mesh"]:
            am1 = 2 * runs[n][0.0125]["field"] - runs[n][0.025]["field"]
            phi = np.array(ref["reference_phi_on_mesh"][str(n)])
            vm, p = solve_fields(n, phi, am1)
            seed = seed_mask(n, cells8)
            out["python with NEM alpha"][str(n)] = {"seed vM": float(vm[seed].mean()), "seed p": float(p[seed].mean())}
            print(f"{n:>3}^3  with the element's alpha - 1 (dt -> 0): seed vM {vm[seed].mean():.4e} "
                  f"({(vm[seed].mean() / rows[n]['seed vM'] - 1) * 100:+.1f} % against the reference field), "
                  f"seed p {p[seed].mean():+.4e} ({(p[seed].mean() / rows[n]['seed p'] - 1) * 100:+.1f} %)", flush=True)
    nem = nem_seed_vm()
    if nem:
        out["nem"] = {str(n): v for n, v in nem.items()}
        print("element runs, seed vM dt -> 0:", ", ".join(f"{n}^3 {v['vM dt->0']:.4e}" for n, v in nem.items()))
        if "seed vM extrapolation" in out:
            v0 = out["seed vM extrapolation"]["h->0"]
            print("  against the Python h -> 0 value:", ", ".join(f"{n}^3 {(v['vM dt->0']/v0-1)*100:+.1f} %" for n, v in nem.items()))
    a.out.write_text(json.dumps(out, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
