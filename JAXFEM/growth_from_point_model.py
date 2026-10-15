"""Growth without kappa: the volume increase implied by the point model's own production.

The point model of Klempt et al. 2026 lives in a fixed volume: biomass grows only into its
void phi_0, so a packed biofilm (phi_0 ~ 0, the calibration's Day-1 state) cannot grow.
In a growing body the produced biomass can instead push the material apart.  Reading:
at every step the point model is put to the local amount phi_3D (phi_0 = 1 - phi_3D,
as in scheme a), advanced at fixed volume, and the biomass it produced,
    R = sum_i phi_i(after) / sum_i phi_i(before),
is taken as a volume increase J -> J R (the material expands so that phi_0 returns to
1 - phi_3D).  Then ln J = sum ln R, alpha = J^(1/3), with no rate parameter: the rate is
the calibrated point model's.  Scheme (a) already does this rescaling and discards R.

Assumptions: the reading above (produced biomass expands the volume instead of filling
the void); phi_3D fixed at the seed's value; time link T* = 1 = Day 21 (s = 0.2375);
point model as calibrated (CH MAP of the element runs, c* = 25, alpha* = 0).
Run:  python JAXFEM/growth_from_point_model.py"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat/coupling"))
sys.path.insert(0, str(ROOT))
import ecology_jax as ej  # noqa: E402

R_DIR = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"


def day1_state_and_theta(cond="CH"):
    pm = json.load(open(R_DIR / f"w8_{cond}_cal_k0882_pm.json")); C = {k: i for i, k in enumerate(pm["columns"])}
    theta = np.array([float(x) for x in pm["theta"]])
    rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
    d = {int(r[C["isubst"]]): r for r in rows}
    r0 = d[min(d)]
    g0 = np.array([float(r0[C[f"gold_{i}"]]) for i in range(1, 13)])
    chi = np.median([float(r[C["phi_used"]]) for r in d.values()])
    return g0, theta, chi


def run(cond="CH", phi3d=None, steps=95, s=0.2375, n_sub=25):
    g, theta, chi = day1_state_and_theta(cond)
    phi3d = chi if phi3d is None else phi3d
    if g[0:5].sum() == 0:                      # point model not started at T* = 0: Day-1 shares
        raise SystemExit("no Day-1 state in the trace")
    lnJ, rigid = [0.0], [g[0:5].sum()]
    gr = g.copy()
    for k in range(steps):
        tot = g[0:5].sum()
        g[0:5] *= phi3d / tot; g[5] = 1 - phi3d
        g, _ = ej.ecology_substeps(g, theta, s / steps, n_sub)
        g = np.asarray(g, float).copy()
        lnJ.append(lnJ[-1] + np.log(g[0:5].sum() / phi3d))
        gr, _ = ej.ecology_substeps(gr, theta, s / steps, n_sub); gr = np.asarray(gr, float).copy()
        rigid.append(gr[0:5].sum())
    t = np.linspace(0, 1, steps + 1)
    return t, np.array(lnJ), np.array(rigid), phi3d


def run_with_phi(phi_steps, g, theta, s=0.2375, n_sub=25, phi_cap=0.9, weighted=True):
    """Same reading with a given phi_3D per step (e.g. phi_used of an element trace), capped
    at phi_cap as the element's scaled scheme does.  Returns ln J after every step.
    weighted (15 Oct, A-10 option 3, the element's mode 10): ln J += ln(1 + phi_3D (R - 1)) with
    the uncapped phi_3D in [0, 1], so that a point with little biomass grows little (at the front,
    phi_3D just above phi_min, R was about 10 in one step); False: ln J += ln R."""
    lnJ = [0.0]; g = np.array(g, float)
    for p3 in phi_steps:
        p = min(p3, phi_cap); tot = g[0:5].sum()
        g[0:5] *= p / tot; g[5] = 1 - p
        g, _ = ej.ecology_substeps(g, theta, s / len(phi_steps), n_sub)
        g = np.asarray(g, float).copy()
        R = g[0:5].sum() / p
        w = min(max(p3, 0.0), 1.0) if weighted else 1.0
        lnJ.append(lnJ[-1] + np.log(1.0 + w * (R - 1.0)))
    return np.array(lnJ)


if __name__ == "__main__":
    for cond in ("CH", "CS"):
        t, lnJ, rigid, chi = run(cond)
        J = np.exp(lnJ)
        i3 = int(round(0.1 * 95))
        print(f"{cond} (phi_3D = {chi:.3f}): V(Day 21)/V(Day 1) = {J[-1]:.3f}, "
              f"V(21)/V(3) = {J[-1] / J[i3]:.3f}, alpha-1 at T*=1 = {J[-1] ** (1 / 3) - 1:.4f}; "
              f"rigid point model sum phi: {rigid[0]:.4f} -> {rigid[-1]:.4f}")
    for p in (0.5, 0.7, 0.9):
        t, lnJ, rigid, _ = run("CH", phi3d=p)
        print(f"CH with phi_3D = {p}: V(21)/V(1) = {np.exp(lnJ[-1]):.3f}")
