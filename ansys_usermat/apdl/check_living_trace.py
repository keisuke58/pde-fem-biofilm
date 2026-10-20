#!/usr/bin/env python3
"""check_living_trace.py -- replay check of an ANSYS run of the coupling agreed with Oliver (A-19,
CLOUD_TO_IKMHIWI03_2026-10-20.ja.md): every traced step of the point model is recomputed in Python from the step's
start state and its inputs, and compared with what ANSYS stored.

Per traced Gauss point and step (the converged, i.e. last, row of each (ldstep, isubst)):
    start state gold_1..12, end state gnew_1..12, n_sub, dt_couple (= s dt), and the inputs c0, c1, a0, a1
    (nutrient and antibiotic at the start and end of the increment), L (the living biomass sum phi_i psi_i passed back).
Python: n_sub steps of dt_couple / n_sub with c* = c*_0 ((1 - w) c0 + w c1), alpha* = alpha*_0 ((1 - w) a0 + w a1),
w = k / n_sub, k = 1..n_sub (the rule of coupled_living_consumption.py), constants from the point-model file
(default abaqus_composition/eco_fritsch_2sp.txt). Also checks L = sum phi_i psi_i of gnew, and that
the start state of each step is the end state of the step before (the state is carried, not rescaled).

    python ansys_usermat/apdl/check_living_trace.py JOB [--results DIR] [--eco FILE] [--max-steps N]
Pass: max |gnew(ANSYS) - gnew(Python)| below ~1e-10 (Newton tolerance), carried-state gaps 0, L error ~1e-15.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT / "ansys_usermat/coupling", ROOT):
    sys.path.insert(0, str(p))


def read_eco(path):
    tok = [ln.split() for ln in Path(path).read_text().splitlines() if ln.strip()]
    n = int(tok[0][0]); cs, als = map(float, tok[1]); eta = list(map(float, tok[2])); th = list(map(float, tok[3]))
    nw = [(int(t[1]), float(t[2])) for t in tok if t[0] == "newton"]
    return n, cs, als, eta, th, (nw[0] if nw else (6, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job")
    ap.add_argument("--results", default=str(HERE / "results" / "2026-10-living"))
    ap.add_argument("--eco", default=str(ROOT / "abaqus_composition" / "eco_fritsch_2sp.txt"))
    ap.add_argument("--max-steps", type=int, default=0)
    a = ap.parse_args()
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    import ecology_jax as eco
    import jax_hamilton_0d_5species_demo as demo

    n_act, cs0, al0, eta, th, (nw_max, nw_tol) = read_eco(a.eco)
    A, B = eco.theta_to_matrices(jnp.asarray(th, dtype=jnp.float64))
    res_dir = Path(a.results)
    path = res_dir / f"{a.job}_pm.json"
    if not path.exists():
        path = res_dir / f"trace_{a.job}.json"          # A-19 export (fields_living.py): key "rows", column "hit"
    pm = json.load(open(path))
    if "seed_trace" not in pm:
        pm["seed_trace"] = pm["rows"]
    C = {k: i for i, k in enumerate(pm["columns"])}
    need = ["c0", "c1", "a0", "a1", "L"]
    miss = [k for k in need if k not in C]
    if miss:
        raise SystemExit(f"trace has no columns {miss}: add them in the material routine (A-19, check 1)")
    P0 = eco.default_hparams(1.0)
    P0.update(A=A, b_diag=B, active_mask=eco.active_mask(n_act),
              Eta=jnp.asarray(eta, dtype=jnp.float64), EtaPhi=jnp.asarray(eta, dtype=jnp.float64))
    P0 = jax.tree_util.tree_map(jnp.asarray, P0)
    mask = P0["active_mask"]

    def newton(g_prev, q):
        """the Fortran point model's Newton (ecology_native.f): up to nw_max iterations, stop once |Q|^2 < nw_tol.
        The first step from the start state (psi = 0.999, gamma = 0) needs ~40 iterations (A-19), so a fixed 6 or 12
        leaves it unconverged and the result then depends on the iteration budget."""
        F = lambda gg: demo.residual(gg, g_prev, q)

        def body(s):
            it, g, _ = s
            g = demo.clip_state(g, mask)
            Q = F(g)
            done = jnp.sum(Q * Q) < nw_tol
            step = demo.clip_state(g + jnp.linalg.solve(jax.jacfwd(F)(g), -Q), mask)
            return it + 1, jnp.where(done, g, step), done

        return jax.lax.while_loop(lambda s: (s[0] < nw_max) & ~s[2], body, (0, g_prev, False))[1]

    @jax.jit
    def _scan(g, dt, ws, c0, c1, a0, a1):
        def body(gg, w):
            q = dict(P0)
            q["dt_h"] = dt
            q["c"] = cs0 * ((1 - w) * c0 + w * c1)
            q["alpha"] = al0 * ((1 - w) * a0 + w * a1)
            return newton(gg, q), None
        return jax.lax.scan(body, g, ws)[0]

    def replay(g, dt, n_sub, c0, c1, a0, a1):
        ws = jnp.arange(1, n_sub + 1, dtype=jnp.float64) / n_sub
        return _scan(g, jnp.float64(dt), ws, c0, c1, a0, a1)

    rows = pm["seed_trace"]
    if "hit" in C:                                      # 0 = computed, 1 = cached repeat, 2 = not started
        rows = [r for r in rows if r[C["hit"]].strip() == "0"]
    last = {}
    for r in rows:
        last[(r[C["elem"]], r[C["ip"]], int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
    points = sorted({k[:2] for k in last})
    worst, worst_L, worst_gap, n_checked, n_missing = 0.0, 0.0, 0.0, 0, 0
    for pt in points:
        steps = sorted(k[2:] for k in last if k[:2] == pt)
        if a.max_steps:
            steps = steps[:a.max_steps]
        prev, prev_st = None, None
        for st in steps:
            r = last[pt + st]
            g0 = np.array([float(r[C[f"gold_{i}"]]) for i in range(1, 13)])
            g1 = np.array([float(r[C[f"gnew_{i}"]]) for i in range(1, 13)])
            # the trace records a call only while phi >= the trace threshold, so a step can be missing: compare
            # only consecutive substeps (the first substep of a new load step follows any substep of the last)
            if prev is not None:
                if st == (prev_st[0], prev_st[1] + 1) or (st[0] == prev_st[0] + 1 and st[1] == 1):
                    worst_gap = max(worst_gap, float(np.abs(g0 - prev).max()))
                else:
                    n_missing += 1
            prev, prev_st = g1, st
            if g0[:n_act].sum() == 0:            # point model not started yet
                continue
            n_sub = int(float(r[C["n_sub"]])); dt = float(r[C["dt_couple"]]) / n_sub
            c0, c1, a0, a1 = (float(r[C[k]]) for k in ("c0", "c1", "a0", "a1"))
            g = np.asarray(replay(jnp.asarray(g0), dt, n_sub, c0, c1, a0, a1))
            worst = max(worst, float(np.abs(g - g1).max()))
            phi, psi = g1[:n_act], g1[6:6 + n_act]
            worst_L = max(worst_L, abs(float(r[C["L"]]) - float((phi * psi).sum())))
            n_checked += 1
    print(f"{a.job}: {len(points)} traced points, {n_checked} steps replayed")
    print(f"  max |gnew(ANSYS) - gnew(Python)| = {worst:.3e}")
    print(f"  max |L(ANSYS) - sum(phi psi)| = {worst_L:.3e}")
    print(f"  max gap between a step's start and the previous step's end (carried state) = {worst_gap:.3e}"
          f" ({n_missing} places with substeps missing from the trace)")
    print(f"  Newton as in the point-model file: up to {nw_max} iterations, |Q|^2 < {nw_tol:g}")
    ok = worst < 1e-8 and worst_L < 1e-12 and worst_gap < 1e-14
    print("PASS" if ok else "CHECK: differences above the expected round-off")


if __name__ == "__main__":
    main()
