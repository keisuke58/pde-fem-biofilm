#!/usr/bin/env python3
"""What the zero-flux boundary fix does to the reported 2D condition spread.

`PDE_VERIFICATION_FINDINGS.md` measured the committed Neumann wall against
exact solutions: one order lost in 1D, inconsistent in 2D, and ~5 % of the
species mass leaking through a wall declared zero-flux. That leak sits on the
path `condition_spread_2d.run` takes --
`diffusion_step_species_2d` (core_hamilton_2d_nutrient.py:426) calls
`laplacian_2d_neumann` at line 434 -- so the question this answers is whether
the reported condition spread is a property of the model or of the wall.

Both operators run **in one process, at one horizon, on one set of
constants**, which is the discipline `CLAIMS_AND_EVIDENCE.md` §4 exists to
enforce: a number carried over from another script or another run is not a
baseline.

It also carries its own internal control, for free. The `control` arm of the
study sets `D_eff = 0`, so the Laplacian is multiplied by zero and the
boundary treatment cannot reach it. **The two control spreads must therefore
come out bit-identical**, and if they do not, this comparison is measuring
something other than the boundary.

    python JAXFEM/boundary_fix_impact.py --nx 8 --macro 300 --dt-h 1e-4
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from JAXFEM import core_hamilton_2d_nutrient as H          # noqa: E402
from JAXFEM import condition_spread_2d as cs               # noqa: E402
from JAXFEM import pde_code_verification as v              # noqa: E402
from jax_hamilton_0d_5species_demo import THETA_DEMO       # noqa: E402


def _assert_patch_bites(module, op):
    """Fail loudly if the operator now installed is not the one being used.

    Calls the module's own species-diffusion step on a probe field and checks
    that it agrees with the same step computed from `op` directly. Without
    this, patching the wrong module object reports "no effect" rather than an
    error -- see the comment in main().
    """
    rng = np.random.default_rng(12345)
    phi = jnp.asarray(rng.uniform(0.05, 0.15, (8, 8, 5)))
    D = jnp.asarray([1.0e-3, 1.0e-3, 8.0e-4, 5.0e-4, 2.0e-4])
    dt, h = 2.0e-3, 1.0 / 7.0
    got = np.asarray(module.diffusion_step_species_2d(phi, D, dt, h, h))
    want = np.asarray(
        jnp.clip(jnp.stack([phi[:, :, i] + dt * D[i] * op(phi[:, :, i], h, h)
                            for i in range(5)], axis=-1), 0.0, 1.0))
    s = jnp.sum(jnp.asarray(want), axis=-1)
    want = np.asarray(jnp.asarray(want)
                      * jnp.where(s > 1.0, 1.0 / s, 1.0)[:, :, None])
    if not np.allclose(got, want, rtol=0.0, atol=0.0):
        raise AssertionError(
            "the installed Laplacian is not the one diffusion_step_species_2d "
            f"is using (max difference {np.max(np.abs(got - want)):.3e}) -- "
            "the patch went to a different module object")


def one_study(cfg, comps):
    """Both arms, returning {tag: {cond: alpha}} plus the spread per arm."""
    out, spreads = {}, {}
    for spatial in (False, True):
        tag = "full" if spatial else "control"
        rec = {}
        for name, _, _ in cs.CONDITIONS:
            al, phi2, c = cs.run(THETA_DEMO, cfg, comps[name], spatial=spatial)
            rec[name] = {"alpha_mean": al,
                         "separation": cs.separation(phi2),
                         "c_mean": float(c.mean())}
        al = np.array([rec[n]["alpha_mean"] for n, _, _ in cs.CONDITIONS])
        spreads[tag] = float((al.max() - al.min()) / al.mean())
        out[tag] = rec
    return out, spreads


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=8)
    ap.add_argument("--macro", type=int, default=300)
    ap.add_argument("--react-sub", type=int, default=20)
    ap.add_argument("--dt-h", type=float, default=1.0e-4)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("JAXFEM/klempt2024_results/boundary_fix_impact.json"))
    a = ap.parse_args(argv)

    cfg = H.Config2D(Nx=a.nx, Ny=a.nx, n_macro=a.macro, dt_h=a.dt_h,
                     n_react_sub=a.react_sub, save_every=10 ** 9)
    print(f"grid {cfg.Nx}x{cfg.Ny}  n_macro={cfg.n_macro}  "
          f"dt_macro={cfg.dt_macro:.2e}  total t={cfg.dt_macro * cfg.n_macro:.2e}\n")

    comps = cs.day1_compositions()
    # Patch the module object `condition_spread_2d` itself holds, NOT the one
    # this file imported. `condition_spread_2d` does a top-level
    # `import core_hamilton_2d_nutrient as H` while this file does
    # `from JAXFEM import core_hamilton_2d_nutrient`, and because both the
    # repository root and JAXFEM/ are on sys.path those are TWO DISTINCT
    # module objects. Patching the wrong one is a silent no-op: the first
    # version of this script did exactly that and reported a bit-identical
    # "no effect" for both arms, which is what a comparison of one operator
    # against itself looks like. `_assert_patch_bites` below turns that
    # failure mode into an error instead of a finding.
    target = cs.H
    committed = target.laplacian_2d_neumann
    results = {}
    # Named explicitly rather than by capture order: after the fix landed, the
    # committed operator IS the mirrored one, so "whatever the module had at
    # import time" stopped being a meaningful label for either arm.
    try:
        for label, op in (("pre-fix (ghost = boundary node)",
                           v.laplacian_2d_neumann_half),
                          ("committed (ghost = node 1, mirrored)", committed)):
            target.laplacian_2d_neumann = op
            _assert_patch_bites(target, op)
            out, spreads = one_study(cfg, comps)
            ratio = (spreads["full"] / spreads["control"]
                     if spreads["control"] else float("inf"))
            sep = max(out["full"][n]["separation"] for n, _, _ in cs.CONDITIONS)
            results[label] = {"spreads": spreads, "ratio": ratio,
                              "separation": sep, "conditions": out}
            print(f"{label}")
            print(f"   control {spreads['control']:.4%}  ->  "
                  f"full {spreads['full']:.4%}   ({ratio:.2f}x)")
            print(f"   largest species separation: {sep:.4f}\n")
    finally:
        target.laplacian_2d_neumann = committed

    keys = list(results)
    c0 = results[keys[0]]["spreads"]["control"]
    c1 = results[keys[1]]["spreads"]["control"]
    print(f"internal control (D_eff = 0, so the wall cannot reach it):")
    print(f"   {c0:.10%} vs {c1:.10%}   "
          f"{'identical' if c0 == c1 else 'DIFFERENT -- see the docstring'}")
    f0 = results[keys[0]]["spreads"]["full"]
    f1 = results[keys[1]]["spreads"]["full"]
    print(f"\nfull arm moves {f0:.4%} -> {f1:.4%} "
          f"({(f1 / f0 - 1) * 100:+.2f} % relative) when the wall is fixed")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(
        {"grid": cfg.Nx, "n_macro": cfg.n_macro, "dt_macro": cfg.dt_macro,
         "dt_h": cfg.dt_h, "results": results}, indent=1))
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
