"""validation_0d_check.py -- does the point model as the coupled code runs it (ecology_native.f through ctypes, the
file of write_eco_cfg.py) reproduce the calibration pipeline's 0D composition at the data days?

Pipeline (Tmcmc202601, data_5species/main/hamilton_ode_jax_paper.simulate_0d_full): 2500 steps of 1e-4, 12 Newton
iterations per step (stop at |Q|^2 < 1e-20), c* = 25, alpha* = 0, gate off, start make_initial_state(phi_init_state);
day d at step round(d * 0.95 * 0.25 / 21 / 1e-4) (estimate_commensal_static.convert_days_to_model_time): 1 -> 113,
3 -> 339, 6 -> 679, 10 -> 1131, 15 -> 1696, 21 -> 2375.

Checked here, for each condition's file abaqus_composition/eco_cal5_<C>.txt:
  (i)  the same steps (dt 1e-4, 2375 steps, c = 1): the Fortran point model vs the pipeline;
  (ii) the coupled code's stepping (coupled_living_consumption.py --eco-file --dtc 0.0025): 113 pre-steps of 1e-4, then
       400 outer steps of dT* = 0.0025 with s = 0.2262 and N_SUB = ceil(0.2262 * 0.0025 / 1e-4) = 6 inner steps of
       9.425e-5 (c = 1 throughout, as the calibration has no nutrient field), sampled at T* = (d - 1)/20.
Reported: max |phibar_i| difference (phibar_i = phi_i psi_i) over the five species at days 3, 6, 10, 15, 21.

    python ansys_usermat/apdl/validation_0d_check.py --tmcmc /path/to/tmcmc202601 [--cond CH]
        -> ansys_usermat/apdl/results/2026-10-validation/check_0d.json (merged over conditions)
Without --cond, one subprocess per condition: the Fortran module reads its file once per process (eco_inited), so a
second FN.init in the same process keeps the first condition's constants (found 30 Oct: CS/DH/DS then differed by 0.4).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat/coupling"))
import ecology_native_ctypes as FN  # noqa: E402

TMCMC = Path(sys.argv[sys.argv.index("--tmcmc") + 1]) if "--tmcmc" in sys.argv else Path("/home/user/tmcmc202601")
sys.path.insert(0, str(TMCMC / "data_5species/main"))
import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402
import hamilton_ode_jax_paper as HP  # noqa: E402

DAYS = [1, 3, 6, 10, 15, 21]
N_STEPS, DT = 2500, 1e-4
IDX = {d: int(round(d * 0.95 * N_STEPS * DT / 21 / DT)) for d in DAYS}
DT_C, S_LINK, PRE = 0.0025, (IDX[21] - IDX[1]) * DT, IDX[1]
N_SUB = max(1, math.ceil(S_LINK * DT_C / DT - 1e-9))


def read_file(p):
    lines = [l for l in Path(p).read_text().splitlines() if l.strip()]
    th = [float(v) for v in lines[3].split()]
    pi = [float(v) for v in next(l for l in lines[4:] if l.startswith("phi_init")).split()[1:6]]
    return th, pi


def start(pi):
    g = np.zeros(12); g[:5] = pi; g[5] = 1.0 - sum(pi); g[6:11] = 0.999
    return g


def phibar(g):
    return np.asarray(g[..., :5]) * np.asarray(g[..., 6:11])


RES = HERE / "results" / "2026-10-validation"
RES.mkdir(exist_ok=True)
OUT = RES / "check_0d.json"
COND = sys.argv[sys.argv.index("--cond") + 1] if "--cond" in sys.argv else None
if COND is None:
    import subprocess
    for C in ("CH", "CS", "DH", "DS"):
        subprocess.run([sys.executable, __file__, "--tmcmc", str(TMCMC), "--cond", C], check=True)
    print("wrote", OUT)
    sys.exit(0)
out = json.loads(OUT.read_text()) if OUT.exists() else dict(idx=IDX, dt_outer=DT_C, s=S_LINK, presteps=PRE,
                                                           n_sub=N_SUB, conditions={})
for C in (COND,):
    f = ROOT / "abaqus_composition" / f"eco_cal5_{C}.txt"
    th, pi = read_file(f)
    cst = FN.init(f)
    assert cst["nact"] == 5 and cst["c_star"] == 25.0 and cst["newton_max"] == 12
    ref = np.asarray(HP.simulate_0d_full(jnp.asarray(th), n_steps=N_STEPS, dt=DT, phi_init=jnp.asarray(pi),
                                         K_hill=0.0, n_hill=4.0, c_const=25.0))      # (2501, 12)
    # (i) same steps
    g = start(pi)[None]
    same = {}
    k = 0
    for d in DAYS:
        n = IDX[d] - k
        g, ok = FN.step_lin(g, n * DT, n, 1.0, 1.0, 0.0, 0.0); assert ok.all()
        k = IDX[d]
        same[d] = float(np.abs(phibar(g[0]) - phibar(ref[IDX[d]])).max())
    # (ii) the coupled code's stepping
    g, ok = FN.step_lin(start(pi)[None], PRE * DT, PRE, 1.0, 1.0, 0.0, 0.0); assert ok.all()
    coup = {1: float(np.abs(phibar(g[0]) - phibar(ref[IDX[1]])).max())}
    pb = {1: phibar(g[0]).tolist()}
    for kk in range(1, int(round(1.0 / DT_C)) + 1):
        g, ok = FN.step_lin(g, S_LINK * DT_C, N_SUB, 1.0, 1.0, 0.0, 0.0); assert ok.all()
        T = kk * DT_C
        for d in DAYS[1:]:
            if abs(T - (d - 1) / 20.0) < 1e-9:
                coup[d] = float(np.abs(phibar(g[0]) - phibar(ref[IDX[d]])).max())
                pb[d] = phibar(g[0]).tolist()
    out["conditions"][C] = dict(same_steps_max_abs_diff=same, coupled_stepping_max_abs_diff=coup,
                                pipeline_phibar={d: phibar(ref[IDX[d]]).tolist() for d in DAYS},
                                coupled_phibar=pb)
    print(f"{C}: (i) same steps max|dphibar| {max(same[d] for d in DAYS[1:]):.2e}; "
          f"(ii) coupled stepping {max(coup[d] for d in DAYS[1:]):.2e}  "
          f"(per day: {', '.join(f'{d}: {coup[d]:.1e}' for d in DAYS[1:])})")
OUT.write_text(json.dumps(out, indent=1) + "\n")
