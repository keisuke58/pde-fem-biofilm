"""
growth_jacobian.py
==================
Composition sensitivity:  dσ / dφ_i  via JAX forward-mode AD.

Pipeline  (from thesis MAP):
  φ = [φ_Pg, φ_Fn, φ_other]   composition fractions, Σφ_i = 1
  k_eff^b(φ) = Σ_i φ_i K_i^b          volumetric growth coupling
  σ̂(φ) = C · (k_eff^b)^b              surrogate (LOO 2.3%)

Question: which species contributes most to σ when its fraction increases by 1%?

  ∂σ̂/∂φ_i = C · b · (k_eff^b)^(b-1) · K_i^b

This is analytic, but JAX jacfwd verifies it and generalises
to any differentiable function of φ.

Parameters (from thesis Table 5.1 / posterior CI analysis)
----------------------------------------------------------
  K^b : volumetric growth coupling coefficients [1/(Pa·s)] (relative, from Klempt)
        K_Pg^b = 2.0  (Pg drives growth: high k_eff^b in CH-dominant region)
        K_Fn^b = 0.8  (Fn moderate)
        K_Ot^b = 0.2  (other, background)

  C   = 1.0  (scale, absorbed into stress units)
  b   = 2.68 (exponent from surrogate fit, LOO 2.3%)

MAP compositions (from thesis Table 4.2):
  CH: φ_Pg=0.42, φ_Fn=0.33, φ_Ot=0.25
  DH: φ_Pg=0.21, φ_Fn=0.48, φ_Ot=0.31

Usage
-----
  python growth_jacobian.py
"""
from __future__ import annotations

import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
K_B = jnp.array([2.0, 0.8, 0.2])    # K^b for [Pg, Fn, Other]
SPECIES = ["Pg", "Fn", "Other"]
C_SCALE = 1.0                         # absorbed into stress normalisation
B_EXP   = 2.68                        # surrogate exponent

# MAP compositions from thesis
PHI_CH = jnp.array([0.42, 0.33, 0.25])
PHI_DH = jnp.array([0.21, 0.48, 0.31])


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
def k_eff_b(phi: jnp.ndarray) -> jnp.ndarray:
    """Effective volumetric growth coupling."""
    return jnp.dot(phi, K_B)


def sigma_hat(phi: jnp.ndarray) -> jnp.ndarray:
    """Surrogate stress (scalar)."""
    return C_SCALE * k_eff_b(phi) ** B_EXP


# ---------------------------------------------------------------------------
# Jacobian  dσ/dφ  via JAX
# ---------------------------------------------------------------------------
grad_sigma = jax.grad(sigma_hat)   # R^n → R, grad is shape (n,)


def sensitivity_report(phi: np.ndarray, label: str) -> None:
    phi_j = jnp.array(phi)
    k     = float(k_eff_b(phi_j))
    sig   = float(sigma_hat(phi_j))
    dsig  = np.array(grad_sigma(phi_j))

    # Analytic check: ∂σ/∂φ_i = C·b·k^(b-1)·K_i^b
    dsig_analytic = C_SCALE * B_EXP * k**(B_EXP - 1) * np.array(K_B)

    # Normalised sensitivity: (φ_i / σ) · ∂σ/∂φ_i  (elasticity)
    elasticity = phi_j * dsig / sig

    print(f"\n{'='*52}")
    print(f"  Condition: {label}")
    print(f"  φ = {dict(zip(SPECIES, [float(x) for x in phi_j]))}")
    print(f"  k_eff^b  = {k:.4f}")
    print(f"  σ̂        = {sig:.4f}  (relative units)")
    print()
    print(f"  {'Species':8s}  {'∂σ/∂φ (JAX)':>14s}  {'∂σ/∂φ (analytic)':>18s}  {'elasticity ε_i':>14s}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*18}  {'-'*14}")
    for i, sp in enumerate(SPECIES):
        print(f"  {sp:8s}  {dsig[i]:>14.4f}  {dsig_analytic[i]:>18.4f}  {float(elasticity[i]):>14.4f}")

    max_i = int(np.argmax(np.abs(dsig)))
    print(f"\n  → 1% increase in φ_{SPECIES[max_i]} changes σ by {dsig[max_i]*0.01:.4f} ({dsig[max_i]*0.01/sig*100:.2f}%)")


def ratio_jacobian() -> None:
    """Jacobian of the CH/DH stress ratio w.r.t. both composition vectors."""
    phi_ch = jnp.array(PHI_CH)
    phi_dh = jnp.array(PHI_DH)
    sig_ch = sigma_hat(phi_ch)
    sig_dh = sigma_hat(phi_dh)
    ratio  = sig_ch / sig_dh

    # ∂ratio/∂φ_CH = (1/σ_DH) ∂σ_CH/∂φ_CH
    dr_dphi_ch = jnp.array(grad_sigma(phi_ch)) / sig_dh
    # ∂ratio/∂φ_DH = -(σ_CH/σ_DH²) ∂σ_DH/∂φ_DH
    dr_dphi_dh = -sig_ch / sig_dh**2 * jnp.array(grad_sigma(phi_dh))

    print(f"\n{'='*52}")
    print(f"  Stress ratio σ_CH/σ_DH = {float(ratio):.3f}  (MAP)")
    print()
    print("  ∂ratio/∂φ_CH  (sensitivity of ratio to CH composition):")
    for sp, v in zip(SPECIES, dr_dphi_ch):
        print(f"    {sp:8s}: {float(v):+.3f}")
    print("  ∂ratio/∂φ_DH  (sensitivity of ratio to DH composition):")
    for sp, v in zip(SPECIES, dr_dphi_dh):
        print(f"    {sp:8s}: {float(v):+.3f}")

    print()
    max_ch = SPECIES[int(jnp.argmax(jnp.abs(dr_dphi_ch)))]
    max_dh = SPECIES[int(jnp.argmax(jnp.abs(dr_dphi_dh)))]
    print(f"  → ratio most sensitive to φ_{max_ch} in CH, φ_{max_dh} in DH")


if __name__ == "__main__":
    sensitivity_report(np.array(PHI_CH), "CH (commensal)")
    sensitivity_report(np.array(PHI_DH), "DH (dysbiotic)")
    ratio_jacobian()

    # Sanity: JAX vs analytic
    phi_t = jnp.array([0.33, 0.33, 0.34])
    jax_g = np.array(grad_sigma(phi_t))
    k_t   = float(k_eff_b(phi_t))
    ana_g = C_SCALE * B_EXP * k_t**(B_EXP - 1) * np.array(K_B)
    max_err = np.max(np.abs(jax_g - ana_g)) / np.max(np.abs(ana_g))
    print(f"\n  JAX vs analytic max relative error: {max_err:.2e}  (expect <1e-14)")
