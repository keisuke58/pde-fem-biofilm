"""
afm_prony_fit.py
================
Fit 2-term Prony viscoelastic parameters from AFM stress-relaxation data.

AFM protocol (constant indentation δ₀, record force F(t)):
  F(t) = F₀ · [e_∞ + e₁·exp(-t/τ₁) + e₂·exp(-t/τ₂)]
  e_∞ + e₁ + e₂ = 1

Maps to UMAT parameters:
  C10   = G₀ / 2          (G₀ = instantaneous shear modulus, from Hertz fit)
  alpha_1 = e₁            (PROPS(6))
  alpha_2 = e₂            (PROPS(7))
  eta_1 = G₀·e₁·τ₁ / 2   (PROPS(4))
  eta_2 = G₀·e₂·τ₂ / 2   (PROPS(5))

Usage
-----
  python afm_prony_fit.py                    # synthetic data demo
  python afm_prony_fit.py --data data.csv    # real data (columns: time[s], force[nN])
  python afm_prony_fit.py --data data.csv --f0 1.2  # specify F(0) manually

Input CSV format:
  time_s,force_nN      (header row)
  0.0,1.234
  0.1,1.200
  ...

Literature reference parameters (used for synthetic data):
  Stoodley et al. 2002 (P. aeruginosa):   G₀≈1kPa, τ₁≈3s, τ₂≈120s
  Fabbri et al. 2014   (dental biofilm):  G₀≈2kPa, e_∞≈0.15, e₁≈0.45, τ₁≈5s
  Gloag et al. 2020    (P. gingivalis):   G₀≈0.5kPa, τ_fast≈2s
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
from scipy.optimize import curve_fit, differential_evolution
from scipy.stats import t as t_dist

# ─────────────────────────────────────────────────────────────────────────────
# Prony model
# ─────────────────────────────────────────────────────────────────────────────
def prony2(t: np.ndarray, e_inf: float, e1: float,
           tau1: float, tau2: float) -> np.ndarray:
    """
    Normalised relaxation: F(t)/F(0) = e_inf + e1·exp(-t/τ1) + e2·exp(-t/τ2)
    with e2 = 1 - e_inf - e1.
    """
    e2 = 1.0 - e_inf - e1
    return e_inf + e1 * np.exp(-t / tau1) + e2 * np.exp(-t / tau2)


def prony2_abs(t: np.ndarray, F0: float, e_inf: float, e1: float,
               tau1: float, tau2: float) -> np.ndarray:
    """Absolute force (for fitting raw F vs t)."""
    return F0 * prony2(t, e_inf, e1, tau1, tau2)


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic data generator
# ─────────────────────────────────────────────────────────────────────────────
def make_synthetic(true_params: dict, T: float = 200.0, N: int = 500,
                   noise_frac: float = 0.02) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic AFM relaxation curve with Gaussian noise.

    true_params: dict with keys G0_Pa, e_inf, e1, tau1_s, tau2_s
    noise_frac : RMS noise as fraction of F(0)
    """
    t = np.linspace(0, T, N)
    e_inf = true_params['e_inf']
    e1    = true_params['e1']
    tau1  = true_params['tau1_s']
    tau2  = true_params['tau2_s']
    F0    = true_params['G0_Pa'] * 1e-9   # convert to nN (arbitrary scale)
    F     = prony2_abs(t, F0, e_inf, e1, tau1, tau2)
    rng   = np.random.default_rng(42)
    F    += rng.normal(0, noise_frac * F0, size=F.shape)
    return t, F


# ─────────────────────────────────────────────────────────────────────────────
# Fitting
# ─────────────────────────────────────────────────────────────────────────────
def fit_prony(t: np.ndarray, F: np.ndarray,
              F0: Optional[float] = None) -> dict:
    """
    Fit 2-term Prony to AFM force relaxation.

    Parameters
    ----------
    t  : time array [s]
    F  : force array [nN] (or any consistent unit)
    F0 : initial force (if None, use F[0])

    Returns
    -------
    dict with fitted parameters, confidence intervals, goodness of fit.
    """
    if F0 is None:
        F0 = float(F[0])

    # Normalise
    Fn = F / F0

    # ── Global search first (differential evolution) ──────────────────────
    # Parameters: [e_inf, e1, tau1, tau2]
    # Constraints: e_inf > 0, e1 > 0, e_inf + e1 < 1, 0 < tau1 < tau2
    T_total = float(t[-1])

    def residuals_sq(p):
        e_inf, e1, tau1, tau2 = p
        if e_inf <= 0 or e1 <= 0 or e_inf + e1 >= 1 or tau1 <= 0 or tau2 <= tau1:
            return 1e10
        pred = prony2(t, e_inf, e1, tau1, tau2)
        return np.mean((pred - Fn)**2)

    bounds_de = [
        (0.01, 0.95),           # e_inf
        (0.01, 0.95),           # e1
        (t[1], T_total * 0.3),  # tau1  (fast)
        (t[1], T_total * 2.0),  # tau2  (slow)
    ]
    res_de = differential_evolution(residuals_sq, bounds_de, seed=0,
                                    maxiter=1000, tol=1e-8, polish=True)

    p0 = res_de.x

    # ── Local refinement (Levenberg-Marquardt) ────────────────────────────
    bounds_lm = ([0, 0, 0, 0],
                 [1, 1, T_total, T_total * 5])

    try:
        popt, pcov = curve_fit(
            prony2, t, Fn,
            p0=p0,
            bounds=bounds_lm,
            maxfev=10000,
            method='trf',
        )
    except RuntimeError:
        popt = p0
        pcov = np.diag([np.nan]*4)

    e_inf_f, e1_f, tau1_f, tau2_f = popt
    e2_f = 1.0 - e_inf_f - e1_f

    # Confidence intervals (95%)
    dof   = len(t) - len(popt)
    t_val = t_dist.ppf(0.975, dof)
    perr  = np.sqrt(np.diag(pcov)) * t_val if not np.any(np.isnan(pcov)) else np.full(4, np.nan)

    # Goodness of fit
    pred  = prony2(t, *popt)
    resid = Fn - pred
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((Fn - Fn.mean())**2)
    R2    = 1.0 - ss_res / ss_tot
    rmse  = np.sqrt(np.mean(resid**2))

    return {
        'F0': F0,
        'e_inf': e_inf_f, 'e_inf_ci': perr[0],
        'e1':    e1_f,    'e1_ci':    perr[1],
        'e2':    e2_f,
        'tau1':  tau1_f,  'tau1_ci':  perr[2],
        'tau2':  tau2_f,  'tau2_ci':  perr[3],
        'R2': R2, 'rmse': rmse,
        'popt': popt, 'pcov': pcov,
        't': t, 'Fn': Fn, 'pred': pred,
    }


# ─────────────────────────────────────────────────────────────────────────────
# UMAT parameter translation
# ─────────────────────────────────────────────────────────────────────────────
def to_umat_props(fit: dict, G0_Pa: Optional[float] = None) -> dict:
    """
    Convert fitted Prony fractions to UMAT PROPS array.

    G0_Pa: instantaneous shear modulus [Pa] from Hertz fit (required for UMAT).
           If None, a warning is printed and relative values are given.
    """
    e_inf = fit['e_inf']
    e1    = fit['e1']
    e2    = fit['e2']
    tau1  = fit['tau1']
    tau2  = fit['tau2']

    if G0_Pa is not None:
        C10   = G0_Pa / 2.0        # Neo-Hookean: G = 2*C10
        eta1  = C10 * e1 * tau1 * 2.0   # η = G_branch * τ = 2*C10*e1 * τ
        eta2  = C10 * e2 * tau2 * 2.0
        D1    = 1.0 / (0.495 * G0_Pa)   # nearly incompressible κ ≈ 0.495E
    else:
        C10 = eta1 = eta2 = D1 = None

    return {
        'C10_MPa':    C10 * 1e-6 if C10 else None,
        'C01_MPa':    0.0,
        'D1_1_MPa':  D1  * 1e-6 if D1  else None,
        'eta1_MPa_s': eta1 * 1e-6 if eta1 else None,
        'eta2_MPa_s': eta2 * 1e-6 if eta2 else None,
        'alpha1': e1,
        'alpha2': e2,
        'mtype':  0.0,
        'G0_Pa': G0_Pa,
        'tau1_s': tau1,
        'tau2_s': tau2,
        'G_inf_Pa': G0_Pa * e_inf if G0_Pa else None,
        'G1_Pa':    G0_Pa * e1    if G0_Pa else None,
        'G2_Pa':    G0_Pa * e2    if G0_Pa else None,
    }


def print_report(fit: dict, props: dict, true_params: Optional[dict] = None):
    print("\n" + "=" * 58)
    print("  Prony 2-term Fit Results")
    print("=" * 58)
    print(f"  R²   = {fit['R2']:.6f}")
    print(f"  RMSE = {fit['rmse']:.4e}  (normalised force)")
    print()
    print(f"  {'Parameter':12s}  {'Fitted':>10s}  {'±95%CI':>10s}", end="")
    if true_params:
        print(f"  {'True':>10s}", end="")
    print()
    print(f"  {'-'*12}  {'-'*10}  {'-'*10}", end="")
    if true_params:
        print(f"  {'-'*10}", end="")
    print()

    rows = [
        ('e_∞',  fit['e_inf'], fit['e_inf_ci'],
         true_params.get('e_inf') if true_params else None),
        ('e₁',   fit['e1'],    fit['e1_ci'],
         true_params.get('e1') if true_params else None),
        ('e₂=1-e∞-e₁', fit['e2'], np.nan,
         1-true_params['e_inf']-true_params['e1'] if true_params else None),
        ('τ₁ [s]', fit['tau1'], fit['tau1_ci'],
         true_params.get('tau1_s') if true_params else None),
        ('τ₂ [s]', fit['tau2'], fit['tau2_ci'],
         true_params.get('tau2_s') if true_params else None),
    ]
    for name, val, ci, truth in rows:
        ci_str = f"{ci:>10.4f}" if not np.isnan(ci) else "        —"
        row = f"  {name:12s}  {val:>10.4f}  {ci_str}"
        if true_params and truth is not None:
            row += f"  {truth:>10.4f}"
        print(row)

    print()
    print("  UMAT PROPS (for umat_biofilm_visco_2ch.f):")
    print(f"  {'PROPS':8s}  {'Value':>14s}  {'Meaning'}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*30}")
    g0 = props.get('G0_Pa')
    rows_u = [
        ('PROPS(1)', props['C10_MPa'],    'C10   [MPa]  = G₀/2'),
        ('PROPS(2)', props['C01_MPa'],    'C01   [MPa]  = 0 (Neo-Hookean)'),
        ('PROPS(3)', props['D1_1_MPa'],   'D1    [1/MPa]'),
        ('PROPS(4)', props['eta1_MPa_s'], 'η₁    [MPa·s]'),
        ('PROPS(5)', props['eta2_MPa_s'], 'η₂    [MPa·s]'),
        ('PROPS(6)', props['alpha1'],     'α₁    (Prony branch 1)'),
        ('PROPS(7)', props['alpha2'],     'α₂    (Prony branch 2)'),
        ('PROPS(8)', props['mtype'],      'mtype (0=NH)'),
    ]
    for name, val, meaning in rows_u:
        val_str = f"{val:>14.6g}" if val is not None else "  needs G₀[Pa]"
        print(f"  {name:8s}  {val_str}  {meaning}")

    if g0:
        print()
        print(f"  G₀   = {g0:.1f} Pa   (input)")
        print(f"  G_∞  = {props['G_inf_Pa']:.1f} Pa  (long-time, = G₀·e_∞)")
        print(f"  G₁   = {props['G1_Pa']:.1f} Pa  (branch 1)")
        print(f"  G₂   = {props['G2_Pa']:.1f} Pa  (branch 2)")
        print(f"  τ₁   = {props['tau1_s']:.2f} s   (fast)")
        print(f"  τ₂   = {props['tau2_s']:.1f} s  (slow)")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data",  type=str, default=None,
                        help="CSV file (columns: time_s, force_nN)")
    parser.add_argument("--f0",    type=float, default=None,
                        help="Initial force F(0) [nN]; default = data[0]")
    parser.add_argument("--g0",    type=float, default=None,
                        help="Instantaneous shear modulus G₀ [Pa] from Hertz fit")
    parser.add_argument("--plot",  action="store_true")
    parser.add_argument("--save",  type=str, default="afm_prony_fit.pdf")
    args = parser.parse_args()

    # ── Literature values for synthetic demo ──────────────────────────────
    # Based on: Fabbri et al. 2014 (subgingival biofilm AFM relaxation)
    # G₀ ≈ 2 kPa, τ_1 ≈ 5 s (EPS), τ_2 ≈ 80 s (cell/matrix), e_∞ ≈ 0.20
    TRUE_PARAMS = dict(G0_Pa=2000.0, e_inf=0.20, e1=0.45, tau1_s=5.0, tau2_s=80.0)

    if args.data is None:
        print("No --data provided. Using synthetic data (Fabbri 2014 literature values).")
        t, F = make_synthetic(TRUE_PARAMS, T=300.0, N=600, noise_frac=0.02)
        true_ref = TRUE_PARAMS
        G0_fit = TRUE_PARAMS['G0_Pa']
    else:
        data = np.loadtxt(args.data, delimiter=',', skiprows=1)
        t, F = data[:, 0], data[:, 1]
        true_ref = None
        G0_fit = args.g0   # may be None

    fit   = fit_prony(t, F, F0=args.f0)
    props = to_umat_props(fit, G0_Pa=G0_fit)
    print_report(fit, props, true_params=true_ref)

    if not args.plot:
        return

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.subplots_adjust(wspace=0.35)

    # --- Left: force vs time ---
    ax = axes[0]
    t_dense = np.linspace(t[0], t[-1], 1000)
    F_fit   = fit['F0'] * prony2(t_dense, *fit['popt'])
    ax.scatter(t, F, s=4, alpha=0.5, label="Data", color="steelblue")
    ax.plot(t_dense, F_fit, 'r-', lw=2, label="Prony 2-term fit")
    if true_ref:
        F_true = fit['F0'] * prony2(t_dense,
                                     true_ref['e_inf'], true_ref['e1'],
                                     true_ref['tau1_s'], true_ref['tau2_s'])
        ax.plot(t_dense, F_true, 'g--', lw=1.5, label="True (synthetic)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Force [nN]")
    ax.set_title(f"AFM Stress Relaxation  R²={fit['R2']:.4f}")
    ax.legend(fontsize=8)

    # --- Right: residual ---
    ax2 = axes[1]
    ax2.plot(t, (fit['Fn'] - fit['pred']) * 100, 'k.', ms=2, alpha=0.5)
    ax2.axhline(0, color='r', lw=1)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Residual [% of F₀]")
    ax2.set_title(f"RMSE = {fit['rmse']*100:.3f}% F₀")

    # Annotate fitted params
    e_inf, e1, tau1, tau2 = fit['popt']
    txt = (f"e_∞={e_inf:.3f}  e₁={e1:.3f}  e₂={fit['e2']:.3f}\n"
           f"τ₁={tau1:.1f}s   τ₂={tau2:.1f}s")
    if G0_fit:
        txt = f"G₀={G0_fit:.0f}Pa  " + txt
    ax.text(0.98, 0.97, txt, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', alpha=0.8))

    plt.suptitle("AFM Prony 2-Term Fit → UMAT Parameters", fontsize=11)
    plt.savefig(args.save, bbox_inches='tight', dpi=150)
    print(f"\nSaved: {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
