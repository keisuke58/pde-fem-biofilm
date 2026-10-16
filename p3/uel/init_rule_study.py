"""init_rule_study.py -- WP4 (KEIO_PLAN.ja.md): what a newly colonised Gauss point starts its composition
from, compared on the column of driver.py check 2 and on the 3D seed of make_seed_inp.py through the
Abaqus entry point (emulate_abaqus.py).

    python init_rule_study.py [--column 20,40,80] [--seed 8,12] [--rules 0,1,2,3] [--cref 1] [--out results/init_rules.json]

--cref c_ref > 0: the point model sees the local nutrient c / c_ref (PROPS(16)), so the composition differs
from point to point and the rules differ; c_ref = 0 (the decks' default): every point sees c* and the
composition depends only on the time since the point started.

Rules (unified_uel.f90, init_rule): 0 the file's Day-1 state (default), 1 copy of the nearest colonised
Gauss point, 2 phi-weighted mean of the colonised Gauss points within 1.5 element sizes, 3 the point
model's steady state under the local nutrient (c* here: c_ref = 0 in these decks).
Reported (thesis conventions): living biomass sum_i phi_i psi_i (psibar = phi psi), the share of species
1 phi_1 / sum phi_i, growth alpha - 1, von Mises / E.  "late" = the Gauss points that started after the
first growth increment (the newly colonised points, the same set for every rule: taken from the rule-1
run of the same mesh), "front" = the Gauss points of the elements whose mean nodal phi is between 0.2
and 0.5 at the end.
"""
import ctypes as ct
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "JAXFEM")); sys.path.insert(0, str(HERE.parents[1]))
if not os.environ.get("BIOFILM_ECO_CASE"):
    import tempfile
    from make_column_inp import eco_file
    _eco = Path(tempfile.mkdtemp()) / "uel_eco.txt"
    eco_file(_eco)
    os.environ["BIOFILM_ECO_CASE"] = str(_eco)
import emulate_abaqus as em  # noqa: E402
from driver import NDN, NGP, NSV, column, props_base  # noqa: E402

lib = em.lib
lib.uel_set_init_rule.argtypes = [ct.c_int, ct.c_double]
F = em.F


def gp_table(SV):
    """per element: composition arrays (NGP x NSV) -> living, share of species 1, alpha - 1, started, start code
    (sv(14): 1 the default state, 1 + rule when rule 1-3 was applied)"""
    sv = np.array([s[:NSV * NGP].reshape(NGP, NSV) for s in SV])
    liv = (sv[:, :, 0:5] * sv[:, :, 6:11]).sum(-1)
    share1 = sv[:, :, 0] / np.maximum(sv[:, :, 0:5].sum(-1), 1e-30)
    return liv, share1, sv[:, :, 12], sv[:, :, 13] > 0.5, np.rint(sv[:, :, 13]).astype(int)


def stats(code):
    """how many Gauss points started from the default state (the seed and fallbacks) and by rule 1-3"""
    return dict(started_default=int((code == 1).sum()), started_by_rule=int((code >= 2).sum()))


def run_column(n, rule, cref):
    lib.uel_set_init_rule(rule, 0.0)
    nodes, elems, xs = column(n); props = F(props_base(s=1.0, c_ref=cref)); xn = nodes[:, 0]
    SV = [F(np.zeros(em.NSVT)) for _ in elems]
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy()
    U[3::NDN] = np.where(xn <= 0.25 + 1e-12, 0.9, 1e-4); U[4::NDN] = 1.0
    U, SV, _ = em.increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, 1e-8, 1)
    fixed = np.zeros(len(U), bool); fixed[0::NDN] = fixed[1::NDN] = fixed[2::NDN] = True
    fixed[4::NDN] = np.isclose(xn, 1.0)
    for k in range(60):
        U, SV, _ = em.increment(U.copy(), U, fixed, nodes, elems, SV, props, 0.2375 / 60, k + 1, ttot=1e-8 + k * 0.2375 / 60)
    phi = U[3::NDN].reshape(n + 1, 4)[:, 0]
    liv, share1, xi, started, code = gp_table(SV)
    pe = 0.5 * (phi[:-1] + phi[1:])                       # element mean phi
    # the front: x where phi crosses 0.5 (linear between nodes)
    i = int(np.where(phi >= 0.5)[0].max())
    xf = xs[i] + (phi[i] - 0.5) / max(phi[i] - phi[i + 1], 1e-30) * (xs[i + 1] - xs[i])
    front = started & ((pe > 0.2) & (pe < 0.5))[:, None]
    return dict(n=n, rule=rule, c_ref=cref, front_x=float(xf), alpha_seed=float(xi[0].mean()), alpha_mean_started=float(xi[started].mean()),
                living_seed=float(liv[0].mean()), living_front=float(liv[front].mean()) if front.any() else None,
                share1_seed=float(share1[0].mean()), share1_front=float(share1[front].mean()) if front.any() else None,
                started_gps=int(started.sum()), **stats(code)), (liv, share1, xi, code)


def run_seed(n, rule, cref):
    lib.uel_set_init_rule(rule, 0.0)
    em.props_base = lambda **kw: props_base(**{**kw, "c_ref": cref})     # emulate_abaqus.seed has no c_ref argument
    t0 = time.time()
    r = em.seed(n, md=False, nu=0.3)
    liv, share1, xi, started, code = gp_table(r["SV"])
    phi = r["U"][3::NDN]
    pe = np.array([phi[el].mean() for el in r["elems"]])
    front = started & ((pe > 0.2) & (pe < 0.5))[:, None]
    return dict(n=n, rule=rule, c_ref=cref, alpha_seed=float(r["alpha_seed"]), alpha_mean_started=float(xi[started].mean()),
                vm_mean=float(r["vm_mean"]), vm_max=float(r["vm_max"]), edge=float(r["edge"]),
                living_front=float(liv[front].mean()) if front.any() else None,
                share1_front=float(share1[front].mean()) if front.any() else None,
                living_seed=float(liv[(pe > 0.5)[:, None] & started].mean()),
                started_gps=int(started.sum()), seconds=round(time.time() - t0), **stats(code)), (liv, share1, xi, code)


def add_late(rows, arrays):
    """the late-started Gauss points (code >= 2 in a rule-1/2/3 run of this mesh): living, share of species 1
    and alpha - 1 averaged over them, for every rule of the mesh"""
    late = None
    for r, (liv, share1, xi, code) in zip(rows, arrays):
        if r["rule"] >= 1:
            late = code >= 2
            break
    for r, (liv, share1, xi, code) in zip(rows, arrays):
        if late is None or not late.any():
            r.update(late_gps=0)
        else:
            r.update(late_gps=int(late.sum()), living_late=float(liv[late].mean()), share1_late=float(share1[late].mean()),
                     alpha_late=float(xi[late].mean()))


if __name__ == "__main__":
    get = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
    cols = [int(v) for v in get("--column", "20,40,80").split(",") if v]
    seeds = [int(v) for v in get("--seed", "8,12").split(",") if v]
    rules = [int(v) for v in get("--rules", "0,1,2,3").split(",")]
    cref = float(get("--cref", "1.0"))
    out = Path(get("--out", HERE / "results" / "init_rules.json"))
    res = dict(column=[], seed=[])
    for kind, ns, fn in (("column", cols, run_column), ("seed", seeds, run_seed)):
        for n in ns:
            rows, arrays = [], []
            for rule in rules:
                r, a = fn(n, rule, cref); rows.append(r); arrays.append(a)
            add_late(rows, arrays)
            for r in rows:
                res[kind].append(r); print(kind, r, flush=True)
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    print(f"wrote {out}")
