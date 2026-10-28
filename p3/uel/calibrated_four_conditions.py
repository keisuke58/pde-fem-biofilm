"""calibrated_four_conditions.py -- WP5 / item (g) of KEIO_PLAN.ja.md: the calibrated five-species point model
(the final MAPs of Chapter 3: CS, CH, DS at stage (iv), DH at stage (iii)) inside the user element
unified_uel.f90, through the Abaqus entry point UEL (emulate_abaqus.py), for all four conditions.

    python calibrated_four_conditions.py [--import-maps /path/to/tmcmc202601 [--ref origin/claude/gate-off-map-check]]
    python calibrated_four_conditions.py [--conditions CS,CH,DS,DH] [--n 8] [--no-seed] [--figure]

For each condition (each in its own process, because ecology_native.f reads the case file once per process):
  (a) uniform check: one element, phi = 0.9 and c prescribed everywhere, c_ref = 0 (the point model sees c*
      as in the calibration), the INIT increment (1e-8, the point model not yet started, as in the seed) and 95
      increments of T/95 with 25 inner steps each, i.e.
      the point model's step 1e-4 of the calibration (2375 steps; the 60 increments of the earlier CH seed runs
      without inner steps give a different trajectory: the first step from the Day-1 state is stiff in psi,
      and T/60 is 40 times the calibration's step), T = 0.2375 (the end time of the CH seed runs, the
      calibration's Day 1 to Day 21: t = 0.2375 d / 21). The composition at the 8 Gauss points is compared with
      the point model run alone on the same inner grid:
        - the calibration's own integrator, Tmcmc202601 hamilton_ode_jax_paper.newton_step (clone at TMCMC_DIR,
          default /home/user/tmcmc202601) -> the WP5 criterion, 1e-6;
        - ecology_jax.ecology_substeps (the material server's fixed six Newton iterations), for information;
        - JAXFEM/unified_1d_monolithic's damped Newton solved to |Q| < 1e-11 every step (core_hamilton_1d.residual,
          no clip, from the clipped start) and the same with 4 times finer steps: the Newton and the time-step
          error of the 1e-4 grid (information).
  (b) the 8^3 seed of make_seed_inp.py (finite strain, nu = 0.49, mean dilatation, L0 = 2, D = 0.2, lam = 0),
      new Gauss points started with init rule 1 (nearest colonised Gauss point, WP4), 95 increments x 25 inner steps to T;
      recorded at every increment: the seed-mean (Gauss points of the elements with mean nodal phi > 0.5 at
      the start) shares phi_i / sum_j phi_j and living fractions phi_i psi_i, and at the end the same for the
      late-started points, alpha - 1 and von Mises / E of the seed.
The point model's constants come from the case file written by abaqus_composition/write_eco_cfg.py
(c* and alpha* from ecology_constants.py, eta = 1, b = 0 as in the MAP files, phi_init = the phi_init_state of
the calibration as make_initial_state builds it (see day1_phi), newton 12 1e-20 = the calibration's own Newton control
(hamilton_ode_jax_paper.newton_step; the CH file of make_column_inp.eco_file has 60 1e-26, which differs from it by
5.6e-4 at Day 21 for CH because the clipped Newton of the first steps, phi_0 at its bound 1e-10, has no root
below 1e-20 and the two controls stop at different iterates).

MAPs: results/calibrated_maps.json, imported with --import-maps from the Tmcmc202601 clone
(docs/revision/generated/final_theta_MAP/<C>.json and phi_init.json at the given ref; the source commit of
each file is recorded). Results: results/calibrated_<C>.json per condition and results/calibrated_four_conditions.json;
figure results/calibrated_four_conditions.pdf (--figure).
"""
import argparse
import ctypes as ct
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = HERE / "results"
MAPS = RES / "calibrated_maps.json"
COND = ("CS", "CH", "DS", "DH")
SPECIES = ("So", "An", "Vei", "Fn", "Pg")
T, NINC, NSUB, T0 = 0.2375, 95, 25, 1e-8      # 95 x 25 inner steps of 1e-4: the calibration's grid (2375 steps)
MAP_DIR = "docs/revision/generated/final_theta_MAP"


def import_maps(tmcmc, ref):
    """theta, phi_init_state and the source commit of each final MAP file, from the Tmcmc202601 clone"""
    def show(path):
        return subprocess.run(["git", "show", f"{ref}:{path}"], cwd=tmcmc, check=True, capture_output=True, text=True).stdout

    def last_commit(path):
        return subprocess.run(["git", "log", "-1", "--format=%H", ref, "--", path], cwd=tmcmc, check=True,
                              capture_output=True, text=True).stdout.strip()
    sha = subprocess.run(["git", "rev-parse", ref], cwd=tmcmc, check=True, capture_output=True, text=True).stdout.strip()
    pi = json.loads(show(f"{MAP_DIR}/phi_init.json"))
    out = {"source": {"repository": "keisuke58/tmcmc202601", "ref": ref, "commit": sha, "directory": MAP_DIR,
                      "phi_init_commit": last_commit(f"{MAP_DIR}/phi_init.json")}}
    for c in COND:
        d = json.loads(show(f"{MAP_DIR}/{c}.json"))
        theta = [float(d[k]) for k in sorted(d, key=int)]
        assert len(theta) == 20 and [int(k) for k in sorted(d, key=int)] == list(range(20))
        out[c] = dict(theta=theta, phi_init_state=[float(v) for v in pi[c]["phi_init_state"]], run=pi[c]["run"],
                      data_source=pi[c]["data_source"], two_channel=pi[c]["two_channel"],
                      theta_file=f"{MAP_DIR}/{c}.json", theta_commit=last_commit(f"{MAP_DIR}/{c}.json"))
    MAPS.write_text(json.dumps(out, indent=1))
    print(f"wrote {MAPS} from {tmcmc} at {ref} ({sha[:7]})")


def write_eco(cond, path):
    """the case file of the condition, through write_eco_cfg.py (the official writer)"""
    sys.path.insert(0, str(ROOT / "abaqus_composition"))
    import write_eco_cfg as wec
    m = json.loads(MAPS.read_text())[cond]
    tj = path.with_suffix(".theta.json")
    tj.write_text(json.dumps({str(i): v for i, v in enumerate(m["theta"])}))
    wec.main_theta(str(tj), str(path), phi_init=day1_phi(m).tolist(), newton=(12, 1e-20))
    tj.unlink()
    return m


def day1_phi(m):
    """the forward model's start as the calibration builds it (hamilton_ode_jax_paper.make_initial_state: the
    five phi_i scaled up to a sum of 0.999999 when their sum is below it; a sum above it is left as it is, so CH
    (sum 1.001) and DH (1.0004; the clip of Pg / An at 0.001) start with phi_0 < 0, which the first clip_state /
    eco_clip of the integrator moves to 1e-10). run_calibrated.ps1 passes the same values."""
    p = np.array(m["phi_init_state"], float)
    return p * (0.999999 / min(p.sum(), 0.999999))


def start_state(m):
    g = np.zeros(12)
    g[0:5] = day1_phi(m); g[5] = 1 - g[0:5].sum(); g[6:11] = 0.999
    return g


def composition(sv_gp):
    """(..., 12) states -> shares phi_i / sum phi_j, living phi_i psi_i, sum of the living"""
    phi, psi = sv_gp[..., 0:5], sv_gp[..., 6:11]
    share = phi / np.maximum(phi.sum(-1, keepdims=True), 1e-300)
    liv = phi * psi
    return share, liv, liv.sum(-1)


TMCMC = Path(os.environ.get("TMCMC_DIR", "/home/user/tmcmc202601"))


def calibration_stepper(theta, dt):
    """the calibration's own integrator (Tmcmc202601 data_5species/main/hamilton_ode_jax_paper.newton_step: implicit
    Euler, up to 12 Newton iterations to |Q|^2 < 1e-20, clip_state), with the constants of ecology_constants.py
    (c* = 25, alpha* = 0, Hill gate off, Kp1 = 1e-4): (step, converged) functions, or None when the clone is not
    at TMCMC_DIR. converged(g_new, g_prev) says whether the step's Newton reached its tolerance."""
    main = TMCMC / "data_5species" / "main"
    if not (main / "hamilton_ode_jax_paper.py").exists():
        return None
    import jax
    import jax.numpy as jnp
    import ecology_constants as ec
    sys.path.insert(0, str(main))
    import hamilton_ode_jax_paper as hp
    A, b = hp.theta_to_matrices(jnp.asarray(theta, dtype=jnp.float64))
    params = {"dt_h": dt, "Kp1": 1e-4, "Eta": jnp.ones(5), "EtaPhi": jnp.ones(5), "c": ec.C_STAR, "alpha": ec.ALPHA_STAR,
              "K_hill": jnp.array(ec.K_HILL), "n_hill": jnp.array(ec.N_HILL), "A": A, "b_diag": b,
              "active_mask": jnp.ones(5, dtype=jnp.int64)}
    step = jax.jit(lambda gg: hp.newton_step(gg, params))
    res = jax.jit(lambda gn, gp: jnp.sum(hp.residual(gn, gp, params) ** 2))
    return (lambda g: np.asarray(step(jnp.asarray(g, dtype=jnp.float64))),
            lambda gn, gp: bool(res(jnp.asarray(gn), jnp.asarray(gp)) < 1e-20))


def calibration_trajectory(stepper, g, steps):
    """states after every step and whether each step's Newton converged"""
    step, conv = stepper
    out, ok = [], []
    for _ in range(steps):
        gn = step(g); out.append(gn); ok.append(conv(gn, g)); g = gn
    return np.array(out), np.array(ok)


def clip_state(g, e=1e-10):
    """ecology_native.eco_clip / the calibration's clip_state: phi, phi_0, psi into [e, 1 - e]"""
    g = np.array(g, float); g[0:11] = np.clip(g[0:11], e, 1 - e); return g


def converged_trajectory(g, theta, dt, steps):
    """unified_1d_monolithic.point_model_converged (core_hamilton_1d.residual, damped Newton to |Q| < 1e-11
    every step), compiled once, returning the state after every step"""
    import jax
    import jax.numpy as jnp
    import core_hamilton_1d as ch
    import ecology_jax as ej
    from unified_1d_monolithic import theta_to_matrices
    A, b = theta_to_matrices(jnp.asarray(theta))
    p = jax.tree_util.tree_map(jnp.asarray, ej.default_hparams(dt)); p["A"] = A; p["b_diag"] = b
    F = jax.jit(lambda gn, gp: ch.residual(gn, gp, p)); J = jax.jit(jax.jacfwd(lambda gn, gp: ch.residual(gn, gp, p)))
    g = jnp.asarray(g); out = []
    for _ in range(steps):
        gp, x = g, g
        for _ in range(40):
            r = F(x, gp)
            if float(jnp.max(jnp.abs(r))) < 1e-11:
                break
            d = jnp.linalg.solve(J(x, gp), -r)
            v, dv = x[0:11], d[0:11]
            lim = jnp.where(dv < 0, -0.9 * v / jnp.minimum(dv, -1e-300), jnp.where(dv > 0, 0.9 * (1 - v) / jnp.maximum(dv, 1e-300), jnp.inf))
            x = x + min(1.0, float(jnp.min(lim))) * d
        g = x; out.append(np.asarray(g))
    return np.array(out)


def uniform_check(em, m, theta):
    """(a): one element, everything prescribed, against the point model alone on the same time grid"""
    sys.path[:0] = [str(ROOT / "JAXFEM"), str(ROOT / "ansys_usermat/coupling"), str(ROOT / "ansys_usermat"), str(ROOT)]
    import unified_1d_monolithic as m1
    import ecology_jax as ej
    import driver as drv
    nodes, elems, _ = drv.column(1)
    props = em.F(drv.props_base(s=1.0, D=0.2, L0=2.0, nu=0.49, lam=0.0, c_ref=0.0, n_sub=NSUB))
    SV = [em.F(np.zeros(em.NSVT))]
    U0 = np.zeros(em.NDN * len(nodes)); U = U0.copy(); U[3::em.NDN] = 0.9; U[4::em.NDN] = 1.0
    allfix = np.ones(len(U), bool)
    U, SV, _ = em.increment(U0, U, allfix, nodes, elems, SV, props, T0, 1, finite=True, md=True)
    traj = []
    for k in range(NINC):
        U, SV, _ = em.increment(U.copy(), U, allfix, nodes, elems, SV, props, T / NINC, k + 1, ttot=T0 + k * T / NINC,
                                finite=True, md=True)
        traj.append(SV[0][:em.NSV * em.NGP].reshape(em.NGP, em.NSV)[:, 0:12].copy())
    traj = np.array(traj)                         # (NINC, NGP, 12)
    gp = traj[-1]
    # point model alone on the element's inner grid (NINC x NSUB steps of 1e-4, the calibration's own step), the
    # calibration's own integrator, in two ways:
    #  - per increment, from the element's start-of-increment state (the WP5 criterion: does the element take
    #    the calibration's steps? judged on the increments whose inner Newton steps all converged, because an
    #    unconverged step's iterate is not reproducible between implementations: the Day-1 state puts psi of
    #    Pg at the clip bound within a few steps and the Jacobian is singular there);
    #  - from the Day-1 state alone (cumulative; the figure's dashed lines).
    # (the point model starts in the first growth increment, as in the seed: the INIT increment's start state
    # has phi = 0 at the Gauss points, so the element does not advance it there)
    g0 = start_state(m)
    st = calibration_stepper(theta, T / (NINC * NSUB))
    d_inc, ok_inc = [], []
    if st is not None:
        prev = g0
        for k in range(NINC):
            ck, kk = calibration_trajectory(st, prev, NSUB)
            d_inc.append(float(np.abs(traj[k][:, 0:11] - ck[-1][0:11]).max())); ok_inc.append(bool(kk.all()))
            prev = traj[k][0]
        ref, ok_all = calibration_trajectory(st, g0, NINC * NSUB)
        ref = ref[NSUB - 1::NSUB]
        unconv_steps = [int(i) for i in np.where(~ok_all)[0]]
    # the server's scheme (six fixed iterations; ecology_native.f is its Fortran copy), same grid
    gs = np.asarray(ej.ecology_substeps(g0, theta, T, NINC * NSUB)[0])
    if st is None:
        ref = np.array([gs]); unconv_steps = None
    # an independent Newton solved to convergence every step (core_hamilton_1d.residual, no clip, so it starts
    # from the clipped state), and the same with 4 times finer steps: the time-step error of the 1e-4 grid
    gc = clip_state(g0)
    conv = converged_trajectory(gc, theta, T / (NINC * NSUB), NINC * NSUB)[-1]
    gf = converged_trajectory(gc, theta, T / (4 * NINC * NSUB), 4 * NINC * NSUB)[-1]
    d_inc, ok_inc = np.array(d_inc), np.array(ok_inc, bool)
    d_crit = float(d_inc[ok_inc].max()) if ok_inc.any() else float("nan")
    d_cal = np.abs(gp[:, 0:11] - ref[-1][0:11]).max()
    d_spread = np.abs(gp[:, 0:11] - gp[:1, 0:11]).max()
    d_server = np.abs(gp[:, 0:11] - gs[0:11]).max()
    d_conv = np.abs(gp[:, 0:11] - conv[0:11]).max()
    d_fine = np.abs(conv[0:11] - gf[0:11]).max()
    sh, lv, lt = composition(ref)
    return dict(max_diff_per_increment_converged=d_crit, passes_1e6=bool(d_crit < 1e-6),
                increments_with_unconverged_newton=[int(i) for i in np.where(~ok_inc)[0]],
                max_diff_per_increment_all=float(d_inc.max()) if d_inc.size else float("nan"),
                per_increment_diff=d_inc.tolist(), per_increment_newton_converged=ok_inc.tolist(),
                max_diff_vs_point_model_cumulative=float(d_cal), unconverged_inner_steps_point_model_alone=unconv_steps,
                gauss_point_spread=float(d_spread),
                reference="calibration newton_step (Tmcmc202601 hamilton_ode_jax_paper)" if st is not None else "ecology_jax.ecology_substeps (clone not found)",
                max_diff_vs_server_scheme=float(d_server), max_diff_vs_converged_newton=float(d_conv),
                time_step_error_1em4_vs_4x_finer=float(d_fine),
                element_end_state=gp[0].tolist(), point_model_end_state=ref[-1].tolist(),
                point_model_trajectory=dict(t=(T0 + np.arange(1, NINC + 1) * T / NINC).tolist(), share=sh.tolist(),
                                            living=lv.tolist(), living_total=lt.tolist()))


def seed_run(em, n, rule=1, nu=0.49, md=True, L0=2.0):
    """(b): the 8^3 seed as emulate_abaqus.seed (finite strain, nu 0.49, mean dilatation), init rule 1,
    the seed-mean composition recorded at every increment"""
    from make_seed_inp import cube, seed_phi
    sys.path.insert(0, str(ROOT / "JAXFEM")); sys.path.insert(0, str(ROOT))
    from fuel_front_3d_species import edge_disp
    import driver as drv
    em.lib.uel_set_init_rule.argtypes = [ct.c_int, ct.c_double]
    em.lib.uel_set_init_rule(int(rule), 0.0)
    nodes, elems = cube(n)
    P = drv.props_base(s=1.0, D=0.2, L0=L0, lam=0.0, psiref=0.0, nu=nu, Y=1.0, mu_b=np.log(0.1), n_sub=NSUB)
    props = em.F(P)
    SV = [em.F(np.zeros(em.NSVT)) for _ in elems]
    U0 = np.zeros(em.NDN * len(nodes)); U = U0.copy(); U[3::em.NDN] = seed_phi(nodes); U[4::em.NDN] = 0.1
    U[em.NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = 1.0
    pe = np.array([U[3::em.NDN][el].mean() for el in elems])
    seed_el = pe > 0.5
    U, SV, _ = em.increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, T0, 1, finite=True, md=md)
    fixed = np.zeros(len(U), bool)
    base = np.isclose(nodes[:, 2], 0.0)
    for d in range(3):
        fixed[em.NDN * np.where(base)[0] + d] = True
    fixed[em.NDN * np.where(np.isclose(nodes[:, 0], 1.0))[0] + 4] = True
    t_hist, share_hist, liv_hist, livtot_hist, alpha_hist, ngp_hist = [], [], [], [], [], []
    t_start = time.time()
    for k in range(NINC):
        U, SV, it = em.increment(U.copy(), U, fixed, nodes, elems, SV, props, T / NINC, k + 1, ttot=T0 + k * T / NINC,
                                 finite=True, md=md)
        sv = np.array([s[:em.NSV * em.NGP].reshape(em.NGP, em.NSV) for s in SV])     # (nel, NGP, NSV)
        started = sv[:, :, 13] > 0.5
        sel = seed_el[:, None] & started
        sh, lv, lt = composition(sv[sel][:, 0:12])
        t_hist.append(T0 + (k + 1) * T / NINC); share_hist.append(sh.mean(0)); liv_hist.append(lv.mean(0))
        livtot_hist.append(lt.mean()); alpha_hist.append(sv[sel][:, 12].mean()); ngp_hist.append(int(started.sum()))
        if k % 10 == 9:
            print(f"  increment {k + 1}: {it} Newton iterations, {time.time() - t_start:.0f} s, "
                  f"started Gauss points {int(started.sum())}", flush=True)
    sv = np.array([s[:em.NSV * em.NGP].reshape(em.NGP, em.NSV) for s in SV])
    st = np.array([s[em.NSV * em.NGP:].reshape(em.NGP, 2) for s in SV])
    code = np.rint(sv[:, :, 13]).astype(int)
    late = code >= 2
    sel = seed_el[:, None] & (code > 0)
    sh_s, lv_s, lt_s = composition(sv[sel][:, 0:12])
    out = dict(n=n, init_rule=rule, nu=nu, mean_dilatation=md, L0=L0, T=T, increments=NINC, seed_elements=int(seed_el.sum()),
               seed_gauss_points=int(sel.sum()), started_gauss_points=int((code > 0).sum()), late_gauss_points=int(late.sum()),
               seed_share=sh_s.mean(0).tolist(), seed_living=lv_s.mean(0).tolist(), seed_living_total=float(lt_s.mean()),
               seed_alpha_minus_1=float(sv[sel][:, 12].mean()),
               seed_vm_mean=float(st[seed_el][:, :, 0].mean()), seed_vm_max=float(st[seed_el][:, :, 0].max()),
               seconds=round(time.time() - t_start),
               history=dict(t=t_hist, share=np.array(share_hist).tolist(), living=np.array(liv_hist).tolist(),
                            living_total=livtot_hist, alpha_minus_1=alpha_hist, started_gauss_points=ngp_hist))
    if late.any():
        sh_l, lv_l, lt_l = composition(sv[late][:, 0:12])
        out.update(late_share=sh_l.mean(0).tolist(), late_living=lv_l.mean(0).tolist(), late_living_total=float(lt_l.mean()),
                   late_alpha_minus_1=float(sv[late][:, 12].mean()))
    nn = n + 1
    u = U.reshape(-1, em.NDN)[:, 0:3].reshape(nn, nn, nn, 3)
    out["edge_30deg"] = float(edge_disp(np.linspace(0, 1, nn), u, 30, 0.0))
    return out


def run_condition(cond, n, do_seed):
    """child process: the case file is read by ecology_native.f at the first element call"""
    eco = RES / f"calibrated_{cond}_eco.txt"
    m = write_eco(cond, eco)
    os.environ["BIOFILM_ECO_CASE"] = str(eco)
    sys.path.insert(0, str(HERE))
    import emulate_abaqus as em
    theta = np.array(m["theta"], float)
    res = dict(condition=cond, run=m["run"], theta_file=m["theta_file"], theta_commit=m["theta_commit"], eco_file=eco.name,
               species=list(SPECIES), theta=m["theta"], phi_init_state_file=m["phi_init_state"],
               phi_init_sum_file=float(np.sum(m["phi_init_state"])), phi_init_used=day1_phi(m).tolist(),
               c_star=None, alpha_star=None)
    import ecology_constants as ec
    res.update(c_star=ec.C_STAR, alpha_star=ec.ALPHA_STAR)
    res["uniform"] = uniform_check(em, m, theta)
    u = res["uniform"]
    print(f"{cond} uniform: max |element - point model alone| per increment = {u['max_diff_per_increment_converged']:.2e} "
          f"({'PASS' if u['passes_1e6'] else 'FAIL'} at 1e-6) on the increments with converged Newton; unconverged increments "
          f"{u['increments_with_unconverged_newton']} (max difference there {u['max_diff_per_increment_all']:.1e}); "
          f"cumulative from Day 1 {u['max_diff_vs_point_model_cumulative']:.1e}, unconverged inner steps of the point model alone "
          f"{u['unconverged_inner_steps_point_model_alone']}; Gauss-point spread {u['gauss_point_spread']:.1e}; "
          f"vs the server's six-iteration scheme {u['max_diff_vs_server_scheme']:.1e}; vs a converged Newton "
          f"{u['max_diff_vs_converged_newton']:.1e}; time-step error of dt 1e-4 "
          f"{u['time_step_error_1em4_vs_4x_finer']:.1e}", flush=True)
    if do_seed:
        res["seed"] = seed_run(em, n)
        s = res["seed"]
        print(f"{cond} seed {n}^3: shares " + " ".join(f"{a} {b:.3f}" for a, b in zip(SPECIES, s["seed_share"])) +
              f"; living " + " ".join(f"{b:.3f}" for b in s["seed_living"]) + f" (sum {s['seed_living_total']:.3f}); "
              f"alpha - 1 {s['seed_alpha_minus_1']:.4f}; vM/E {s['seed_vm_mean']:.4e} / {s['seed_vm_max']:.4e}; "
              f"late points {s['late_gauss_points']}; {s['seconds']} s", flush=True)
    (RES / f"calibrated_{cond}.json").write_text(json.dumps(res, indent=1))


def figure(summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "ansys_usermat"))
    import figstyle
    figstyle.apply()
    colors = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]      # Okabe-Ito, fixed order
    fig, axs = plt.subplots(2, 2, figsize=(6.3, 5.0), sharex=True, sharey=True)
    for ax, cond, letter in zip(axs.ravel(), COND, "abcd"):
        r = summary.get(cond)
        if r is None or "seed" not in r:
            ax.text(0.5, 0.5, f"{cond}: no run", ha="center", transform=ax.transAxes)
            continue
        h, pm = r["seed"]["history"], r["uniform"]["point_model_trajectory"]
        days = 21 * np.array(h["t"]) / T
        for i, (name, col) in enumerate(zip(SPECIES, colors)):
            ax.plot(days, np.array(h["share"])[:, i], color=col, lw=1.4, label=name)
            ax.plot(21 * np.array(pm["t"]) / T, np.array(pm["share"])[:, i], color=col, lw=1.0, ls="--")
        ax.text(0.03, 0.93, f"({letter}) {cond}", transform=ax.transAxes, va="top")
        ax.set_ylim(0, 1); ax.set_xlim(0, 21)
        ax.grid(True, lw=0.3, alpha=0.4)
    for ax in axs[1]:
        ax.set_xlabel("$21\\,t/T$ (days after Day 1, $T = 0.2375$)")
    for ax in axs[:, 0]:
        ax.set_ylabel("share $\\phi_i/\\sum_j\\phi_j$")
    h1, l1 = axs[0, 0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    h1 += [Line2D([], [], color="k", lw=1.4), Line2D([], [], color="k", lw=1.0, ls="--")]
    l1 += ["seed in the element", "point model alone"]
    fig.legend(h1, l1, loc="lower center", ncol=7, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = RES / "calibrated_four_conditions.pdf"
    fig.savefig(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--import-maps", help="path of the Tmcmc202601 clone")
    ap.add_argument("--ref", default="origin/claude/gate-off-map-check")
    ap.add_argument("--conditions", default=",".join(COND))
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--no-seed", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--run", help=argparse.SUPPRESS)       # child: one condition
    a = ap.parse_args()
    RES.mkdir(exist_ok=True)
    if a.import_maps:
        import_maps(a.import_maps, a.ref)
    if a.run:
        run_condition(a.run, a.n, not a.no_seed)
        sys.exit()
    conds = [c for c in a.conditions.split(",") if c]
    if conds:
        maps = json.loads(MAPS.read_text())
        for c in conds:
            if c not in maps:
                print(f"{c}: no MAP in {MAPS}, skipped"); continue
            cmd = [sys.executable, str(HERE / "calibrated_four_conditions.py"), "--run", c, "--n", str(a.n)] + (["--no-seed"] if a.no_seed else [])
            subprocess.run(cmd, check=True, cwd=HERE)
        summary = {"source": maps["source"], "T": T, "increments": NINC}
        for c in COND:
            f = RES / f"calibrated_{c}.json"
            if f.exists():
                summary[c] = json.loads(f.read_text())
        (RES / "calibrated_four_conditions.json").write_text(json.dumps(summary, indent=1))
        print("\n| condition | run (MAP) | uniform: max diff vs point model alone | shares So An Vei Fn Pg (seed, end) | living total | alpha - 1 |")
        print("|---|---|---|---|---|---|")
        for c in COND:
            r = summary.get(c)
            if not r:
                continue
            u = r["uniform"]; s = r.get("seed")
            print(f"| {c} | {r['run']} | {u['max_diff_per_increment_converged']:.1e} ({'pass' if u['passes_1e6'] else 'FAIL'}; cumulative {u['max_diff_vs_point_model_cumulative']:.1e}) | "
                  + (" / ".join(f"{v:.3f}" for v in s["seed_share"]) + f" | {s['seed_living_total']:.3f} | {s['seed_alpha_minus_1']:.4f} |" if s else "- | - | - |"))
        if a.figure:
            figure(summary)
