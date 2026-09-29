#!/usr/bin/env python3
"""klempt2026_reproduction.py -- run this repo's 0D Hamilton integrator on the
numerical examples of Klempt et al. (arXiv:2509.01274, cases in
klempt2026_cases.py) and compare with the paper's figures.

What is under test is the repo's own `newton_step` / `residual`
(jax_hamilton_0d_5species_demo.py, the code behind ecology_jax.py and the
ANSYS ecology bridge), unmodified: n<5 species run with the unused slots
masked out, the Hill gate (a repo extension acting on slot 5 only) is off,
and eta, b, c*, alpha* are set per case. The equations there are the paper's
Eq. 16-18 divided through by eta_i. Result: all 14 cases reproduce the
paper's figures, given the two printing inconsistencies documented in
klempt2026_cases.py (no +gamma in Eq. 17; Eq. 19's A without the 1/2).
The variants below keep the rejected alternatives runnable as evidence.

    python JAXFEM/klempt2026_reproduction.py            # all cases
    python JAXFEM/klempt2026_reproduction.py 2sp_case1 4sp_case1
    python JAXFEM/klempt2026_reproduction.py 4sp_case1 --as-printed   # Eq. 19 literally

Writes JAXFEM/klempt2026_results/<case>.png (paper-style panels) and
summary.json, and prints ours vs. the values read off the paper's figures.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
sys.path.insert(0, str(_HERE))

from jax_hamilton_0d_5species_demo import clip_state, newton_step, residual  # noqa: E402
from klempt2026_cases import CASES, DT, PSI0  # noqa: E402


def residual_psi_gamma(g_new, g_prev, params):
    """The repo residual plus gamma/eta_i in each active psi-equation, i.e.
    Eq. 17 exactly as printed (the repo's psi-residual is Eq. 17 / eta_i
    without the +gamma). Used only to test which form the paper's figures
    were computed with; the repo's own code is not changed."""
    Q = residual(g_new, g_prev, params)
    add = params["active_mask"].astype(jnp.float64) * g_new[11] / params["Eta"]
    return Q.at[6:11].add(add)


def _newton(res_fn, iters):
    """The repo's newton_step loop (clip, full Newton update, clip) on res_fn,
    with a configurable iteration count; iters=6 reproduces newton_step."""
    def step(g_prev, params):
        mask = params["active_mask"]

        def body(g, _):
            g = clip_state(g, mask)
            F = lambda gg: res_fn(gg, g_prev, params)  # noqa: E731
            g = g + jnp.linalg.solve(jax.jacfwd(F)(g), -F(g))
            return clip_state(g, mask), None

        g, _ = jax.lax.scan(body, clip_state(g_prev, mask), jnp.arange(iters))
        return g
    return step


newton_step_psi_gamma = _newton(residual_psi_gamma, 6)

# repo       : the production integrator, unmodified (6 Newton iterations/step)
# newton30   : same residual, 30 iterations -- tests whether 6 is converged
# psi_gamma  : Eq. 17 as printed (+gamma); rejected on the 2-species cases
VARIANTS = {"repo": newton_step, "newton30": _newton(residual, 30),
            "psi_gamma": newton_step_psi_gamma}

OUT = _HERE / "klempt2026_results"
NSLOT = 5


def _pad(v, fill, n=NSLOT):
    v = list(v)
    return v + [fill] * (n - len(v))


def _schedule(value, steps):
    if callable(value):
        return np.array([value(k) for k in range(1, steps + 1)], dtype=float)
    return np.full(steps, float(value))


def build(case, psi0=PSI0, as_printed=False):
    n = case["n"]
    A = np.zeros((NSLOT, NSLOT))
    scale = 1.0 if as_printed else case.get("A_figure_scale", 1.0)
    A[:n, :n] = np.asarray(case["A"], dtype=float) * scale
    eta = jnp.asarray(_pad(case["eta"], 1.0), dtype=jnp.float64)
    base = {
        "dt_h": jnp.float64(DT),
        "Kp1": jnp.float64(1e-4),        # Kp_i = eta_i*1e-4 in the paper; the residual is divided by eta_i
        "Eta": eta,
        "EtaPhi": eta,
        "c": jnp.float64(0.0),           # set per step
        "alpha": jnp.float64(0.0),       # alpha* (antibiotics), set per step
        "K_hill": jnp.float64(0.0),      # Hill gate off: not part of the paper's model
        "n_hill": jnp.float64(4.0),
        "A": jnp.asarray(A),
        "b_diag": jnp.asarray(_pad(case["b"], 0.0), dtype=jnp.float64),
        "active_mask": jnp.asarray([1] * n + [0] * (NSLOT - n), dtype=jnp.int64),
    }
    phi = np.asarray(_pad(case["phi_init"], 0.0))
    g0 = np.zeros(12)
    g0[0:5] = phi
    g0[5] = 1.0 - phi.sum()
    g0[6:11] = [psi0] * n + [0.0] * (NSLOT - n)
    return jnp.asarray(g0), base


def _make_integrator(stepper):
    @jax.jit
    def integrate(g0, base, c_arr, a_arr):
        def step(g, xs):
            c, a = xs
            p = dict(base)
            p["c"] = c
            p["alpha"] = a
            g = stepper(g, p)
            return g, g
        _, traj = jax.lax.scan(step, g0, (c_arr, a_arr))
        return traj
    return integrate


_INTEGRATORS = {k: _make_integrator(v) for k, v in VARIANTS.items()}


def simulate(case, psi0=PSI0, variant="repo", as_printed=False):
    g0, base = build(case, psi0, as_printed)
    steps = case["steps"]
    traj = _INTEGRATORS[variant](g0, base, jnp.asarray(_schedule(case["c_star"], steps)),
                                 jnp.asarray(_schedule(case["alpha_star"], steps)))
    traj = np.vstack([np.asarray(g0)[None, :], np.asarray(traj)])
    n = case["n"]
    return {"phi": traj[:, 0:n], "phi0": traj[:, 5], "psi": traj[:, 6:6 + n],
            "gamma": traj[:, 11]}


def observables(res):
    phibar = res["phi"] * res["psi"]
    hit = np.nonzero(res["phi0"] < 0.02)[0]     # the barrier holds phi0 just above 0; 0.02 is "reached ~0" on the paper's plots
    return {
        "phi_end": res["phi"][-1].tolist(),
        "psi_end": res["psi"][-1].tolist(),
        "phibar_end": phibar[-1].tolist(),
        "phi0_zero_step": int(hit[0]) if hit.size else None,
        "finite": bool(np.all(np.isfinite(res["phi"])) and np.all(np.isfinite(res["psi"]))),
        "sum_end": float(res["phi"][-1].sum() + res["phi0"][-1]),
    }


# Categorical slots 1-4 of the validated reference palette (dataviz skill);
# species colour follows the species, phi0 in primary ink, sum in muted ink.
_SPECIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
_INK, _MUTED, _GRID, _SURF = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"


def plot(name, case, res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = case["n"]
    k = np.arange(res["phi"].shape[0])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.6), facecolor=_SURF)
    for ax in (ax1, ax2):
        ax.set_facecolor(_SURF)
        ax.grid(True, color=_GRID, linewidth=0.6)
        for s in ax.spines.values():
            s.set_color(_MUTED)
        ax.tick_params(colors=_MUTED, labelsize=8)
        ax.set_xlabel("time steps", color=_INK, fontsize=9)
    for i in range(n):
        ax1.plot(k, res["phi"][:, i], color=_SPECIES[i], lw=2, label=f"phi{i+1}")
        ax1.plot(k, res["psi"][:, i], color=_SPECIES[i], lw=1.5, ls="--", label=f"psi{i+1}")
        ax2.plot(k, res["phi"][:, i] * res["psi"][:, i], color=_SPECIES[i], lw=2,
                 label=f"phi{i+1}*psi{i+1}")
    ax1.plot(k, res["phi0"], color=_INK, lw=2, label="phi0")
    ax1.plot(k, res["phi"].sum(axis=1) + res["phi0"], color=_MUTED, lw=1, label="sum")
    ax1.set_ylabel("value", color=_INK, fontsize=9)
    ax1.set_title(f"{name}: phi_i (solid), psi_i (dashed)  [cf. {case['figure']}a]",
                  color=_INK, fontsize=9, loc="left")
    ax2.set_title(f"living fraction phi_i*psi_i  [cf. {case['figure']}b]",
                  color=_INK, fontsize=9, loc="left")
    ax1.legend(fontsize=7, frameon=False, ncol=2, labelcolor=_INK)
    ax2.legend(fontsize=7, frameon=False, labelcolor=_INK)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=_SURF)
    plt.close(fig)


def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, list):
        return "[" + ", ".join(f"{x:.3f}" for x in v) + "]"
    return str(v)


def main(names, as_printed=False):
    OUT.mkdir(exist_ok=True)
    summary = {}
    for name in names:
        case = CASES[name]
        entry = {"figure": case["figure"], "paper_figure_reading": case["paper"],
                 "A_as_printed": as_printed}
        obs = {}
        for variant in VARIANTS:
            res = simulate(case, variant=variant, as_printed=as_printed)
            obs[variant] = observables(res)
            plot(f"{name} [{variant}]", case, res, OUT / f"{name}__{variant}.png")
            entry[variant] = obs[variant]
        summary[name] = entry
        print(f"\n== {name}  ({case['figure']}, {case['steps']} steps)  finite: "
              + ", ".join(f"{v}={obs[v]['finite']}" for v in VARIANTS))
        for key, ref in case["paper"].items():
            cols = "  ".join(f"{v} {_fmt(obs[v][key]):30s}" for v in VARIANTS)
            print(f"   {key:15s} {cols} paper~ {_fmt(ref)}")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    args = sys.argv[1:]
    printed = "--as-printed" in args      # use Eq. 19's A literally (with the 1/2) instead of what the figures used
    args = [a for a in args if a != "--as-printed"]
    main(args or list(CASES), as_printed=printed)
