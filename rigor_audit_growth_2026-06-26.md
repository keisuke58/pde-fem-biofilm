# Rigor audit — growth kinematics & numerics (2026-06-26)

Targeted rigor review of 4 implementation areas. **Headline conclusion: the
core thesis result (σ_CH/σ_DH MAP = 6.44×) is rigorously sound — it is built
from the Klempt-faithful UMAT path.** One real inconsistency found in a
*secondary* path (eigenstrain), plus minor notes.

## #1 — Growth parameter α: definition & cross-path consistency

**Source of truth.** Klempt 2024 §2.1 "Kinematics of growth":
> `F = Fe·Fg`, growth part `Fg = α·I`, α = *local expansion parameter*,
> isotropic growth (α = 1 at no growth).

So α is the **per-direction stretch**. The accumulation variable used in the
code is `α_acc = α − 1` (per-direction engineering growth strain), giving
`Fg = (1+α_acc)·I`, `J_g = (1+α_acc)³`.

| Path | α handling | Klempt-faithful? |
|---|---|---|
| **UMAT** `umat_klempt_voigt.f` / `umat_klempt2025.f` (PRODUCTION + headline) | `s_iso = 1+α_acc`, `Fe = F/s` → `Fg=(1+α_acc)I` | ✅ exact |
| **eigenstrain** `compute_alpha_eigenstrain.py` (secondary biofilm-mode) | `eps = α_acc/3` (volumetric split) | ❌ inconsistent (÷3 per dir) |

- **Headline is safe.** `posterior_klempt_stress_ci.py` fits a surrogate to the
  4 UMAT MAP stresses (large-α, Fg=(1+α)I) and propagates posterior φ. No
  eigenstrain involved. → the 6.44× result is Klempt-faithful.
- **Real finding (secondary).** The eigenstrain helper divides the
  per-direction α by 3 (treating it as a volumetric strain), under-applying
  growth ≈3× per direction vs Klempt/UMAT. Klempt-consistent value is
  `eps = α_acc` (per direction); via Abaqus `*Expansion` (alpha_T=1) this gives
  stretch 1+α_acc = **exact finite-strain match to Fg, no linearization**.
  → Annotated in `compute_alpha_eigenstrain.py` (added `eps_growth_klempt`);
  not silently changed because it alters secondary comparison figures and is
  entangled with the helper's k_alpha tuning. **Decision needed** before
  regenerating biofilm-mode comparison runs.

## #2 — Per-species α split ("Nishioka rate-weighted")

`α_s = α_total · (K_α,s·φ_s)/k_α,eff`, `k_α,eff = Σ_s K_α,s·φ_s`.
The condition-level φ_s are **fixed in time** (TMCMC MAP fractions per
condition), so `α_s/α_total = K_α,s φ_s / Σ K_α,s φ_s` holds **exactly** for the
model as posed. The "approx" label refers only to approximating the full
spatiotemporal per-species PDE. → **Exact under the model premise. No fix.**

## #3 — Mesh convergence

- 2D biofilm PDE (`fem_convergence.py` → `convergence_report.md`):
  domain-averaged φ_i converge **< 0.03 %**; P.g spatial pattern L2 **< 1.5 %**. ✅
- Tooth-FEM stress (`p4_mesh_convergence.py`): through-thickness layer sweep
  {4,8,16}, criterion median MISES Δ < 2 %. Study exists; no red flag.
→ **PDE side rigorously converged; FEM-stress convergence study in place.**

## #4 — Viscoelastic time integration (`umat_biofilm_visco.f`)

Viscous flow uses a single-step **backward-Euler** update (1st-order in Δt,
semi-implicit). The exact consistent tangent (verified, 2.9e-8 vs FD) is
consistent *with this update*, so Newton convergence is unaffected; only the
**time-accuracy** is 1st-order. Acceptable for quasi-static growth; an
exponential-map / 2nd-order update would improve transient accuracy if needed.
This UMAT is a beyond-Klempt extension **not used in production**. → low priority.

## Net
Core result rigorous (Klempt-faithful, converged). Actionable items: (a) decide
on the eigenstrain ÷3 fix for secondary biofilm-mode runs (annotated); (b)
optional 2nd-order viscous update if the visco UMAT is ever run transiently.
