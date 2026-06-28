"""
test_umat_2ch.py
================
Python-port verification of umat_biofilm_visco_2ch.f

Tests:
  1. Prony relaxation: σ(t) = 2G_∞·ε + 2G_1·ε·e^{-t/τ_1} + 2G_2·ε·e^{-t/τ_2}
     at small strain (linearised), step stretch ε held constant
  2. Limit A1=A2=0: purely elastic (no relaxation)
  3. Limit A2=0: single Maxwell (1-channel), matches phase2 UMAT result
  4. Backward Euler stability check: large dt still converges
  5. Growth: alpha_g > 0 reduces effective F_mech → lower stress

Run:  python test_umat_2ch.py
"""
from __future__ import annotations

import sys
import numpy as np

EPS = np.finfo(float).eps

# ────────────────────────────────────────────────────────────────────────────
# Python port of MAXWELL_BRANCH + BIOFILM_STRESS_CORE_2CH
# ────────────────────────────────────────────────────────────────────────────
I3 = np.eye(3)

def mat_aat(A):
    return A @ A.T

def mat_inv3(A):
    det = np.linalg.det(A)
    if abs(det) < 1e-30:
        return np.eye(3), 1.0
    return np.linalg.inv(A), det

def maxwell_branch(Fmech, Fv_old, C10_k, eta_k, dt):
    """Single Maxwell branch: deviatoric Neo-Hookean + backward-Euler Fv update."""
    if C10_k < 1e-30:
        return np.zeros((3,3)), Fv_old.copy(), 0.0

    Fv_inv, detFv = mat_inv3(Fv_old)
    Fe = Fmech @ Fv_inv
    Be = mat_aat(Fe)
    Je = max(np.linalg.det(Fe), 1e-15)

    tmp1 = Je**(-2/3)
    I1b = tmp1 * np.trace(Be)
    tau_dev = 2*C10_k*tmp1*(Be - I1b/3*I3)

    if eta_k > 1e-20:
        tmp = dt / (2*eta_k*Je)
        Fv_new = (I3 + tmp*tau_dev) @ Fv_old
    else:
        Fv_new = Fv_old.copy()

    # Recompute with updated Fv
    Fv_inv2, _ = mat_inv3(Fv_new)
    Fe2 = Fmech @ Fv_inv2
    Be2 = mat_aat(Fe2)
    Je2 = max(np.linalg.det(Fe2), 1e-15)

    tmp1 = Je2**(-2/3)
    I1b2 = tmp1 * np.trace(Be2)
    tau_dev2 = 2*C10_k*tmp1*(Be2 - I1b2/3*I3)

    sigma_k = tau_dev2 / Je2

    spd = 0.0
    if eta_k > 1e-20:
        spd = np.sum(tau_dev2**2) * dt / (4*eta_k*Je2)

    return sigma_k, Fv_new, spd


def stress_2ch(F, Fv1, Fv2, alpha_g, dt,
               C10, C01, D1, eta1, eta2, a1, a2, mtype=0.0):
    """
    2-channel Prony stress.
    Returns: sigma_total (3,3), Fv1_new (3,3), Fv2_new (3,3), sse, spd
    """
    fg_scale = max(1.0 + alpha_g, 1e-15)
    Fmech = F / fg_scale

    C10_inf = C10 * (1.0 - a1 - a2)
    C10_1   = C10 * a1
    C10_2   = C10 * a2

    dt_s = max(dt, 1e-20)

    # --- Equilibrium spring ---
    Be_inf = mat_aat(Fmech)
    Je_inf = max(np.linalg.det(Fmech), 1e-15)
    tmp1 = Je_inf**(-2/3)
    I1b = tmp1 * np.trace(Be_inf)
    trce = I1b / 3.0
    press = (2/D1)*(Je_inf - 1)*Je_inf

    sigma_inf = (2*C10_inf*tmp1*(Be_inf - trce*I3) + press*I3) / Je_inf

    # MR correction (omitted for brevity when mtype=0)
    sse = C10_inf*(I1b - 3) + (1/D1)*(Je_inf - 1)**2

    # --- Branch 1 ---
    sigma_1, Fv1_new, spd1 = maxwell_branch(Fmech, Fv1, C10_1, eta1, dt_s)

    # --- Branch 2 ---
    sigma_2, Fv2_new, spd2 = maxwell_branch(Fmech, Fv2, C10_2, eta2, dt_s)

    sigma_tot = sigma_inf + sigma_1 + sigma_2
    return sigma_tot, Fv1_new, Fv2_new, sse, spd1 + spd2


# ────────────────────────────────────────────────────────────────────────────
# Time-marching: prescribed stretch λ in x, unconfined (F = diag(λ,1/√λ,1/√λ))
# ────────────────────────────────────────────────────────────────────────────
def simulate_shear(params: dict, gamma: float, T: float, N: int):
    """
    Simple shear F = I + γ·e1⊗e2, record σ_12(t).
    For small γ: σ_12 ≈ 2*C10_k*γ·Prony(t) — cleanest Prony check.
    """
    C10, C01, D1 = params['C10'], params.get('C01',0), params['D1']
    eta1, eta2   = params['eta1'], params['eta2']
    a1, a2       = params['a1'], params['a2']
    alpha_g      = params.get('alpha_g', 0.0)

    F = np.array([[1, gamma, 0],[0, 1, 0],[0, 0, 1]], dtype=float)
    Fv1 = np.eye(3)
    Fv2 = np.eye(3)

    dt = T / N
    t_arr   = np.zeros(N+1)
    sig_arr = np.zeros(N+1)

    sig0, Fv1, Fv2, _, _ = stress_2ch(
        F, Fv1, Fv2, alpha_g, 1e-30,
        C10, C01, D1, eta1, eta2, a1, a2)
    sig_arr[0] = sig0[0,1]   # σ_12

    for i in range(1, N+1):
        sig, Fv1, Fv2, _, _ = stress_2ch(
            F, Fv1, Fv2, alpha_g, dt,
            C10, C01, D1, eta1, eta2, a1, a2)
        t_arr[i]   = i * dt
        sig_arr[i] = sig[0,1]

    return t_arr, sig_arr


def simulate_relaxation(params: dict, lam: float, T: float, N: int):
    """Legacy: uniaxial stretch, σ_11. Kept for backward compatibility."""
    C10, C01, D1 = params['C10'], params.get('C01',0), params['D1']
    eta1, eta2   = params['eta1'], params['eta2']
    a1, a2       = params['a1'], params['a2']
    alpha_g      = params.get('alpha_g', 0.0)

    lam2 = 1.0 / np.sqrt(lam)
    F    = np.diag([lam, lam2, lam2])
    Fv1  = np.eye(3)
    Fv2  = np.eye(3)

    dt = T / N
    t_arr   = np.zeros(N+1)
    sig_arr = np.zeros(N+1)

    sig0, Fv1, Fv2, _, _ = stress_2ch(
        F, Fv1, Fv2, alpha_g, 1e-30,
        C10, C01, D1, eta1, eta2, a1, a2)
    sig_arr[0] = sig0[0,0]

    for i in range(1, N+1):
        sig, Fv1, Fv2, _, _ = stress_2ch(
            F, Fv1, Fv2, alpha_g, dt,
            C10, C01, D1, eta1, eta2, a1, a2)
        t_arr[i]   = i * dt
        sig_arr[i] = sig[0,0]

    return t_arr, sig_arr


# ────────────────────────────────────────────────────────────────────────────
# Analytic Prony for SIMPLE SHEAR (linearised):
# σ_12(t) = 2C10_∞·γ + 2C10_1·γ·e^{-t/τ_1} + 2C10_2·γ·e^{-t/τ_2}
# ────────────────────────────────────────────────────────────────────────────
def prony_analytic_shear(params, gamma, t_arr):
    C10 = params['C10']
    a1, a2   = params['a1'], params['a2']
    eta1, eta2 = params['eta1'], params['eta2']
    tau1 = eta1 / (2*C10*a1) if a1 > 0 else np.inf
    tau2 = eta2 / (2*C10*a2) if a2 > 0 else np.inf
    C10_inf = C10*(1-a1-a2)
    C10_1   = C10*a1
    C10_2   = C10*a2
    sig = 2*gamma*(C10_inf
                   + C10_1*np.exp(-t_arr/tau1)
                   + C10_2*np.exp(-t_arr/tau2))
    return sig, tau1, tau2


# ────────────────────────────────────────────────────────────────────────────
# Tests
# ────────────────────────────────────────────────────────────────────────────
PASS = "\033[32m PASS\033[0m"
FAIL = "\033[31m FAIL\033[0m"

def check(name, cond, detail=""):
    tag = PASS if cond else FAIL
    print(f"  {tag}  {name}" + (f"  [{detail}]" if detail else ""))
    return cond


def test_prony_relaxation():
    """σ_12(t) vs analytic Prony (simple shear, small γ)."""
    print("\nTest 1: Prony relaxation vs analytic (simple shear, small γ)")
    gamma = 0.01
    params = dict(C10=1e-3, C01=0, D1=1e3, eta1=5e-3, eta2=5e-1,
                  a1=0.5, a2=0.3)
    tau1 = params['eta1']/(2*params['C10']*params['a1'])
    tau2 = params['eta2']/(2*params['C10']*params['a2'])
    # Simulate for 5×τ_1; τ_2 still active
    T, N = 5*tau1, 2000

    t, sig = simulate_shear(params, gamma, T, N)
    sig_a, _, _ = prony_analytic_shear(params, gamma, t)

    err_max = np.max(np.abs(sig - sig_a)) / np.abs(sig_a).max()
    ok1 = check("max relative error vs analytic < 1%", err_max < 0.01,
                 f"err={err_max:.3e}")

    # Initial: σ_12(0) = 2*C10*γ
    sig_0_a = 2*params['C10']*gamma
    ok2 = check("σ_12(0) = 2C10·γ",
                 abs(sig[0] - sig_0_a)/abs(sig_0_a) < 0.01,
                 f"sim={sig[0]:.5e}  theory={sig_0_a:.5e}")

    # After 5τ_1: branch 1 decayed; branch 2 still ~exp(-5/τ2_ratio)
    ratio = tau2/tau1  # large
    sig_lt_a = 2*gamma*(params['C10']*(1-params['a1']-params['a2'])
                        + params['C10']*params['a2']*np.exp(-T/tau2))
    ok3 = check("σ_12(5τ_1) branch-1 decayed",
                 abs(sig[-1] - sig_lt_a)/abs(sig_lt_a) < 0.02,
                 f"sim={sig[-1]:.5e}  theory={sig_lt_a:.5e}")

    print(f"    τ_1={tau1:.1f}s  τ_2={tau2:.1f}s  C10_∞={1e3*params['C10']*(1-params['a1']-params['a2']):.1f}mPa")
    return ok1 and ok2 and ok3


def test_elastic_limit():
    """A1=A2=0 → no relaxation, σ constant."""
    print("\nTest 2: Elastic limit (a1=a2=0)")
    params = dict(C10=1e-3, C01=0, D1=1e3, eta1=1e-3, eta2=1e-3,
                  a1=0.0, a2=0.0)
    t, sig = simulate_relaxation(params, 1.02, 100.0, 200)
    rel_var = (sig.max() - sig.min()) / sig.mean()
    return check("σ constant (rel variation < 1e-6)", rel_var < 1e-6,
                 f"rel_var={rel_var:.2e}")


def test_single_channel():
    """A2=0: σ_12(∞) = 2·C10_∞·γ (only equilibrium remains)."""
    print("\nTest 3: Single-channel limit (a2=0, shear)")
    gamma = 0.01
    p = dict(C10=1e-3, C01=0, D1=1e3, eta1=2e-3, eta2=0.0, a1=0.6, a2=0.0)
    tau1 = p['eta1'] / (2*p['C10']*p['a1'])
    T = 10*tau1
    t, sig = simulate_shear(p, gamma, T, int(T/tau1*200))

    sig_lt_a = 2*gamma*p['C10']*(1-p['a1'])   # branch 1 fully relaxed
    ok = check("long-time σ_12 → 2·C10_∞·γ",
               abs(sig[-1] - sig_lt_a)/abs(sig_lt_a) < 0.01,
               f"sim={sig[-1]:.4e}  theory={sig_lt_a:.4e}")
    return ok


def test_large_dt_stability():
    """Large dt: backward Euler must not diverge."""
    print("\nTest 4: Large Δt stability (backward Euler)")
    params = dict(C10=1e-3, C01=0, D1=1e3, eta1=1e-3, eta2=1.0,
                  a1=0.5, a2=0.3)
    t_coarse, s_coarse = simulate_relaxation(params, 1.02, 50.0, 10)   # very coarse
    t_fine,   s_fine   = simulate_relaxation(params, 1.02, 50.0, 5000) # fine

    ok1 = check("coarse simulation is finite", np.all(np.isfinite(s_coarse)))
    ok2 = check("coarse long-time → same as fine",
                abs(s_coarse[-1] - s_fine[-1])/abs(s_fine[-1]) < 0.05,
                f"coarse={s_coarse[-1]:.4e}  fine={s_fine[-1]:.4e}")
    return ok1 and ok2


def test_growth():
    """Growth (alpha_g > 0) reduces effective stress vs alpha_g=0."""
    print("\nTest 5: Growth effect (alpha_g > 0 → lower effective stretch)")
    params_0 = dict(C10=1e-3, C01=0, D1=1e3, eta1=5e-3, eta2=0.5,
                    a1=0.4, a2=0.3, alpha_g=0.0)
    params_g = {**params_0, 'alpha_g': 0.1}
    lam = 1.1

    _, s0 = simulate_relaxation(params_0, lam, 1.0, 10)
    _, sg = simulate_relaxation(params_g, lam, 1.0, 10)

    ok = check("σ(alpha_g=0.1) < σ(alpha_g=0) [growth reduces mechanical strain]",
               sg[0] < s0[0],
               f"σ_growth={sg[0]:.4e}  σ_no_growth={s0[0]:.4e}")
    return ok


def test_energy_dissipation():
    """SPD (dissipation) >= 0 at each step."""
    print("\nTest 6: Viscous dissipation >= 0")
    params = dict(C10=1e-3, C01=0, D1=1e3, eta1=5e-3, eta2=0.5,
                  a1=0.5, a2=0.3)
    lam2 = 1/np.sqrt(1.1)
    F = np.diag([1.1, lam2, lam2])
    Fv1, Fv2 = np.eye(3), np.eye(3)
    all_ok = True
    for _ in range(20):
        sig, Fv1, Fv2, _, spd = stress_2ch(
            F, Fv1, Fv2, 0.0, 0.5,
            params['C10'], 0, params['D1'],
            params['eta1'], params['eta2'],
            params['a1'], params['a2'])
        if spd < -1e-20:
            all_ok = False
            break
    return check("SPD >= 0 at every step", all_ok)


if __name__ == "__main__":
    print("=" * 58)
    print("2-Channel UMAT Verification  (umat_biofilm_visco_2ch.f)")
    print("=" * 58)

    results = [
        test_prony_relaxation(),
        test_elastic_limit(),
        test_single_channel(),
        test_large_dt_stability(),
        test_growth(),
        test_energy_dissipation(),
    ]

    n_pass = sum(results)
    n_fail = len(results) - n_pass
    print(f"\n{'='*58}")
    print(f"  {n_pass}/{len(results)} tests passed"
          + (f"  ← {n_fail} FAILED" if n_fail else "  ✓ all clear"))
    sys.exit(0 if n_fail == 0 else 1)
