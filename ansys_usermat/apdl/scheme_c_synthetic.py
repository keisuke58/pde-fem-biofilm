"""scheme_c_synthetic.py -- coupling (c), growth from the living biomass, on the Chapter 4 seed problem with the
synthetic four-species point model of Fritsch et al. 2025 (true parameters; no measured data). Not part of the thesis.

Field (as wp3_reference.py, 32^3 cells of the 2 mm cube, RK4): d phi/dt = beta lap phi + k_alpha alpha - penalty,
beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, the 32 seed elements of the partner's 8^3 mesh, zero-flux walls,
T* = 0 .. 1.

Growth:
  (a)  d alpha/dt = k_alpha phi                                    (Eq. 36, the thesis)
  (c)  d alpha/dt = kappa phi sum_i phi_i psi_i                    (living biomass of the point model)
       the point model in every cell keeps its own void (coupling (b)), starts from the synthetic state
       (phi_i = 0.02, psi_i = 0.999), sees c* = 25 c(x) and the antibiotic alpha*(T*), and advances s dT* per coupling
       step dT* = 0.025 with s = 0.15: T* = 1 spans the 1500 steps of 1e-4 of Fritsch's validation run.
       kappa: the seed mean of alpha - 1 at T* = 1 equals that of (a) for uniform c = 1 without antibiotic (the choice
       of the element runs), so (a) and (c) differ only in how the growth is spread over time and space.
Cases of (c): uniform c = 1; antibiotic alpha* = 50 from T* = 0.5 (Fritsch Fig. 18, b = 0.1 .. 0.4); local nutrient,
quasi-static d lap c = g phi c with c = 1 on the face y = -1 mm, Thiele number Lambda^2 = g L^2/d = 16 (L = 1 mm).
Stress: wp3_stress_reference.solve_fields on 16^3 (hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain
(alpha - 1) I, rigid-body constraints only), at T* = 0.5 and 1.

    python ansys_usermat/apdl/scheme_c_synthetic.py
        -> ansys_usermat/apdl/results/2026-10-3d/scheme_c_synthetic.json, assets/fig_scheme_c_synthetic.pdf/.png
"""
import json
import math
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

jax.config.update("jax_enable_x64", True)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (HERE, ROOT / "ansys_usermat", ROOT / "ansys_usermat/coupling", ROOT, ROOT / "JAXFEM"):
    sys.path.insert(0, str(p))
import ecology_jax as eco  # noqa: E402
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402
from fritsch2025_cases import CASES  # noqa: E402

NF, NM, T_END, DT_C, S_LINK, LAM = 32, 16, 1.0, 0.025, 0.15, 4.0
CASE = CASES["fritsch_4sp_val_fig18"]
N_SP = 4
H = W.L / NF
OUT_JSON = HERE / "results" / "2026-10-3d" / "scheme_c_synthetic.json"
OUT_FIG = ROOT / "assets" / "fig_scheme_c_synthetic"


def theta20():
    A, b = CASE["A"], CASE["b"]
    th = [0.0] * 20
    th[0], th[1], th[2], th[3], th[4] = A[0][0], A[0][1], A[1][1], b[0], b[1]
    th[5], th[6], th[7], th[8], th[9] = A[2][2], A[2][3], A[3][3], b[2], b[3]
    th[10], th[11], th[12], th[13] = A[0][2], A[0][3], A[1][2], A[1][3]
    th[14], th[15] = 1.0, 1.0
    return jnp.asarray(th, dtype=jnp.float64)


TH, MASK = theta20(), eco.active_mask(N_SP)
N_SUB = max(1, math.ceil(S_LINK * DT_C / 1e-4 - 1e-9))
STEP_PM = jax.jit(jax.vmap(lambda g, c, a: eco._substep_scan(g, TH, S_LINK * DT_C / N_SUB, jnp.arange(N_SUB), MASK,
                                                               c, a, jnp.ones(5))[0], in_axes=(0, 0, None)))


def g_init(n):
    G = np.zeros((n, 12))
    G[:, :N_SP] = CASE["phi_init"]
    G[:, 5] = 1.0 - sum(CASE["phi_init"])
    G[:, 6:11] = 0.999
    return G


def living(G):
    return (G[:, :N_SP] * G[:, 6:6 + N_SP]).sum(1)


class Nutrient:
    """-lap c + q phi c = 0, c = 1 on the face y = -1 (ghost cell), zero flux elsewhere; q = Lambda^2 /mm^2."""

    def __init__(self, q):
        n3 = NF ** 3
        idx = np.arange(n3).reshape((NF,) * 3)
        rows, cols, diag = [], [], np.zeros((NF,) * 3)
        for ax in range(3):
            for sgn in (1, -1):
                ok = np.ones((NF,) * 3, bool)
                e = [slice(None)] * 3; e[ax] = -1 if sgn == 1 else 0; ok[tuple(e)] = False
                rows.append(idx[ok]); cols.append(np.roll(idx, -sgn, ax)[ok]); diag += ok
        diag[:, 0, :] += 2.0
        self.A0 = sp.csr_matrix((-np.ones(sum(r.size for r in rows)) / H ** 2,
                                 (np.concatenate(rows), np.concatenate(cols))), shape=(n3, n3))
        self.d0 = diag.ravel() / H ** 2
        b = np.zeros((NF,) * 3); b[:, 0, :] = 2.0 / H ** 2
        self.b, self.q, self.c = b.ravel(), q, np.ones(n3)

    def solve(self, phi):
        A = (self.A0 + sp.diags(self.d0 + self.q * np.clip(phi, 0, 1).ravel())).tocsc()
        self.c = spla.spsolve(A, self.b)
        return self.c.reshape((NF,) * 3)


def to_cells(v, np_):
    """Values on np_^3 blocks (or one value) -> the NF^3 cells."""
    if np_ == 1:
        return np.full((NF,) * 3, float(v[0]))
    f = NF // np_
    a = np.asarray(v).reshape((np_,) * 3)
    return np.repeat(np.repeat(np.repeat(a, f, 0), f, 1), f, 2)


def run(case, kappa=None):
    """case: "a", "c", "c_ab" (antibiotic), "c_nut" (local nutrient). Returns time series and end fields."""
    cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")
    seed = W.seed_mask(NF, cells8)
    phi = seed.astype(float)
    am1 = np.zeros_like(phi)
    lap = np.empty_like(phi)
    # the point model does not see phi (own void): with a uniform c every cell is the same point model, so one is
    # solved; with the local nutrient it is solved on the 16^3 blocks (block-mean c) and its living biomass is
    # copied to the 32^3 cells of each block
    NP = NM if case == "c_nut" else 1
    G = g_init(NP ** 3)
    c = np.ones((NF,) * 3)
    nut = Nutrient(LAM ** 2) if case == "c_nut" else None
    dt = 0.15 * H * H / W.BETA
    per = int(np.ceil(DT_C / dt)); dt = DT_C / per
    ts, liv_seed, am1_seed, snaps = [0.0], [float(living(G).mean())], [0.0], {}
    t0 = time.time()
    for k in range(int(round(T_END / DT_C))):
        T = k * DT_C
        if case != "a":
            if nut is not None:
                c = nut.solve(phi)
            a_star = 50.0 if (case == "c_ab" and T >= 0.5 - 1e-12) else 0.0
            L0 = living(G)
            cb = W.block_mean(c, NM).ravel() if NP > 1 else np.ones(1)
            G = np.asarray(STEP_PM(jnp.asarray(G), jnp.asarray(CASE["c_star"] * cb), jnp.float64(a_star)))
            Lv = to_cells(0.5 * (L0 + living(G)), NP)                # trapezoid over the coupling step
        for _ in range(per):
            # growth held over the coupling step: (a) k_alpha phi, (c) kappa phi L
            rate = W.K_ALPHA * np.ones_like(phi) if case == "a" else kappa * Lv
            def f(ph, am):
                W.laplacian(ph, H, lap)
                dph = W.BETA * lap + W.K_ALPHA * (1 + am) - W.PENALTY * (np.maximum(ph - 1, 0) + np.minimum(ph, 0))
                return dph, rate * ph
            k1 = f(phi, am1); k2 = f(phi + 0.5 * dt * k1[0], am1 + 0.5 * dt * k1[1])
            k3 = f(phi + 0.5 * dt * k2[0], am1 + 0.5 * dt * k2[1]); k4 = f(phi + dt * k3[0], am1 + dt * k3[1])
            phi = phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            am1 = am1 + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        ts.append((k + 1) * DT_C)
        liv_seed.append(float(to_cells(living(G), NP)[seed].mean()) if case != "a" else np.nan)
        am1_seed.append(float(am1[seed].mean()))
        if abs(ts[-1] - 0.5) < 1e-9 or abs(ts[-1] - 1.0) < 1e-9:
            snaps[round(ts[-1], 2)] = (phi.copy(), am1.copy(), c.copy(), to_cells(living(G), NP))
    print(f"[{case}] field done, {time.time() - t0:.0f} s; seed alpha-1 at T*=1: {am1_seed[-1]:.4e}", flush=True)
    seed16 = S.seed_mask(NM, cells8)
    out = dict(t=ts, living_seed=liv_seed, am1_seed=am1_seed, stress={})
    for T, (ph, am, cc, lv) in snaps.items():
        vm, p = S.solve_fields(NM, W.block_mean(ph, NM), W.block_mean(am, NM))
        out["stress"][str(T)] = dict(seed_vm_Pa=float(vm[seed16].mean()), seed_p_Pa=float(p[seed16].mean()),
                                     seed_am1=float(am[seed].mean()), max_vm_Pa=float(vm.max()),
                                     seed_c_min=float(cc[seed].min()), seed_c_max=float(cc[seed].max()),
                                     seed_living_min=float(lv[seed].min()), seed_living_max=float(lv[seed].max()))
        if T == 1.0:
            k = NF // 2
            out["mid"] = dict(am1=(0.5 * (am[:, :, k - 1] + am[:, :, k])).tolist(),
                              c=(0.5 * (cc[:, :, k - 1] + cc[:, :, k])).tolist(),
                              vm=(0.5 * (vm[:, :, NM // 2 - 1] + vm[:, :, NM // 2])).tolist())
    return out


def main():
    R = {"a": run("a")}
    # kappa: equal seed growth at T* = 1 for c = 1 without antibiotic, from a unit-kappa run
    probe = run("c", kappa=1.0)
    kappa = R["a"]["am1_seed"][-1] / probe["am1_seed"][-1]
    print(f"kappa = {kappa:.4e} /T*", flush=True)
    for case in ("c", "c_ab", "c_nut"):
        R[case] = run(case, kappa=kappa)
    res = {"model": "Fritsch et al. 2025 Case II (true parameters, validation antibiotic at T* = 0.5), seed problem "
                    "of Chapter 4, 32^3 field, 16^3 stress", "s": S_LINK, "kappa": kappa, "Lambda": LAM,
           "runs": {k: {kk: vv for kk, vv in v.items() if kk != "mid"} for k, v in R.items()}}
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1) + "\n")
    print("\n| run | seed alpha-1 T*=0.5 | T*=1 | seed vM T*=1 (Pa) | seed p T*=1 (Pa) | living seed T*=1 |")
    for k, v in R.items():
        s5, s1 = v["stress"]["0.5"], v["stress"]["1.0"]
        print(f"| {k} | {s5['seed_am1']:.4e} | {s1['seed_am1']:.4e} | {s1['seed_vm_Pa']:.4e} | {s1['seed_p_Pa']:.4e} "
              f"| {v['living_seed'][-1]:.3f} |")
    for k in ("c_nut",):
        s1 = R[k]["stress"]["1.0"]
        print(f"{k}: c in seed {s1['seed_c_min']:.2f}-{s1['seed_c_max']:.2f}, living {s1['seed_living_min']:.3f}-"
              f"{s1['seed_living_max']:.3f}")

    figstyle.apply(size=9)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lab = {"a": "(a) Eq. 36", "c": "(c) uniform $c$", "c_ab": "(c) antibiotic from $T^*=0.5$",
           "c_nut": r"(c) local nutrient, $\Lambda=4$"}
    sty = {"a": dict(color="0.2", ls="--"), "c": dict(color="C0"), "c_ab": dict(color="C3"), "c_nut": dict(color="C2")}
    fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.35))
    for k, v in R.items():
        if k != "a":
            ax[0].plot(v["t"], v["living_seed"], lw=1.2, label=lab[k], **sty[k])
        ax[1].plot(v["t"], np.array(v["am1_seed"]) * 1e3, lw=1.2, label=lab[k], **sty[k])
    ax[0].set(xlabel="$T^*$", ylabel=r"seed mean $\Sigma_i\phi_i\psi_i$", xlim=(0, 1), ylim=(0, None))
    ax[1].set(xlabel="$T^*$", ylabel=r"seed mean $\alpha-1$ ($10^{-3}$)", xlim=(0, 1), ylim=(0, None))
    e = np.linspace(-1, 1, NF + 1)
    im = ax[2].pcolormesh(e, e, np.array(R["c_nut"]["mid"]["am1"]).T * 1e3, cmap="Oranges", shading="flat",
                          rasterized=True)
    cb = fig.colorbar(im, ax=ax[2], fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    cb.set_label(r"$\alpha-1$ ($10^{-3}$)", fontsize=8)
    ax[2].set_aspect("equal"); ax[2].set(xlabel="$x$ (mm)", ylabel="$y$ (mm)", xticks=[-1, 0, 1], yticks=[-1, 0, 1])
    ax[2].grid(False)
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=7, handlelength=1.8, columnspacing=1.0)
    for a, k in zip(ax, "abc"):
        a.text(0.04, 0.96, f"({k})", transform=a.transAxes, va="top", ha="left",
               color="black" if k != "c" else "0.1")
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=0.6)
    fig.savefig(OUT_FIG.with_suffix(".pdf")); fig.savefig(OUT_FIG.with_suffix(".png"), dpi=300)
    print("wrote", OUT_JSON, OUT_FIG.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
