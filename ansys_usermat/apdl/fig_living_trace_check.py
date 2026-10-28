"""fig_living_trace_check.py -- what the ANSYS run of the agreed coupling passes to the point model, and how well the
point model in ANSYS agrees with Python, step by step (A-19 run lv8_sep, 8^3, A in the seed, B in the edge columns).

(a) nutrient at the traced Gauss point of a seed element and of an edge element, the pair (c_n, c_{n+1}) of each
    increment (c_{n+1} of one increment is c_n of the next); (b) max |state(ANSYS) - state(Python)| per increment,
    Python fed with the same start state and the same pair.

    python ansys_usermat/apdl/fig_living_trace_check.py [JOB [TRACE_JSON]]  -> assets/fig_living_trace_check.pdf/.png
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT / "ansys_usermat", ROOT / "ansys_usermat/coupling", ROOT):
    sys.path.insert(0, str(p))
import figstyle  # noqa: E402

JOB = sys.argv[1] if len(sys.argv) > 1 else "lv8_sep"
RES = HERE / "results" / "2026-10-living" / "archive"   # lv8_ and trace_ files archived (22 Oct)


def main():
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    import ecology_jax as eco
    path = Path(sys.argv[2]) if len(sys.argv) > 2 else RES / f"{JOB}_pm.json"
    if not path.exists():
        path = RES / f"trace_{JOB}.json"
    d = json.load(open(path))
    if "rows" not in d:
        d["rows"] = d["seed_trace"]
    C = {k: i for i, k in enumerate(d["columns"])}
    rows = [r for r in d["rows"] if r[C["hit"]].strip() == "0"]
    last = {}
    for r in rows:
        last[(r[C["elem"]], int(r[C["ldstep"]]), int(r[C["isubst"]]))] = r
    th = [1.0, 0.1, 1.0, 1.0, 2.0] + [0.0] * 9 + [1.0, 1.0] + [0.0] * 4
    A, B = eco.theta_to_matrices(jnp.asarray(th))
    P0 = eco.default_hparams(1.0)
    P0.update(A=A, b_diag=B, active_mask=eco.active_mask(2), Eta=jnp.asarray([1.0, 2.0, 1, 1, 1]),
              EtaPhi=jnp.asarray([1.0, 2.0, 1, 1, 1]))
    P0 = jax.tree_util.tree_map(jnp.asarray, P0)

    @jax.jit
    def sc(g, dt, ws, c0, c1, a0, a1):
        def body(gg, w):
            q = dict(P0); q["dt_h"] = dt
            q["c"] = 100 * ((1 - w) * c0 + w * c1); q["alpha"] = 10 * ((1 - w) * a0 + w * a1)
            return eco.newton_step(gg, q), None
        return jax.lax.scan(body, g, ws)[0]

    out = {}
    for e, lab in ((str(d["picked"]["seed"]), "seed"), (str(d["picked"]["edge"]), "edge")):
        T, c0s, c1s, err = [], [], [], []
        for k in sorted(k for k in last if k[0] == e):
            r = last[k]; f = lambda c: float(r[C[c]])  # noqa: E731
            g0 = np.array([f(f"gold_{i}") for i in range(1, 13)]); g1 = np.array([f(f"gnew_{i}") for i in range(1, 13)])
            n = int(f("n_sub")); dt = f("dt_couple") / n
            g = np.asarray(sc(jnp.asarray(g0), dt, jnp.arange(1, n + 1, dtype=float) / n, f("c0"), f("c1"), f("a0"), f("a1")))
            T.append(k[2] * 0.025); c0s.append(f("c0")); c1s.append(f("c1")); err.append(np.abs(g - g1)[:11].max())
        out[lab] = (np.array(T), np.array(c0s), np.array(c1s), np.array(err))
        print(lab, "element", e, "steps", len(T), "max err after step 1", np.max(out[lab][3][1:]), "step 1", out[lab][3][0])

    figstyle.apply(size=9)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.4))
    for lab, col in (("seed", "C0"), ("edge", "C3")):
        T, c0, c1, err = out[lab]
        ax[0].step(T, c1, where="post", color=col, lw=1.2, label=f"{lab}: $c_{{n+1}}$")
        ax[0].plot(T, c0, ".", color=col, ms=3, label=f"{lab}: $c_n$")
        ax[1].semilogy(T, np.maximum(err, 1e-17), "o-", ms=2.5, lw=0.8, color=col, label=lab)
    ax[0].set(xlabel="$T^*$", ylabel="nutrient at the Gauss point", xlim=(0, 1), ylim=(-0.05, 1.1))
    ax[0].legend(fontsize=6.5, frameon=False, ncol=2, loc="lower right")
    ax[1].set(xlabel="$T^*$", ylabel=r"max $|$ANSYS $-$ Python$|$", xlim=(0, 1))
    ax[1].legend(fontsize=7, frameon=False)
    for a, k in zip(ax, "ab"):
        a.text(0.5, 0.97, f"({k})", transform=a.transAxes, va="top", ha="center")
    fig.tight_layout(w_pad=1.0)
    o = ROOT / "assets" / "fig_living_trace_check"
    fig.savefig(o.with_suffix(".pdf")); fig.savefig(o.with_suffix(".png"), dpi=300)
    print("wrote", o.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
