# Gauss-point ecology ODE vs the TMCMC calibration ODE (2026-10-01)

Question: is `ecology_jax` (the Gauss-point point model, dt = 1e−4 inner
steps) the same model as `BiofilmNewtonSolver5S` in tmcmc202601, which the
calibration used with dt = 0.01? Script: `tmcmc_ode_crosscheck.py`.

**Same equations, with one difference, which is a bug here.** The residuals
(φ, φ₀, ψ, γ rows) are written identically. The difference is the Hill gate:

- TMCMC: `if K_hill > 0: Interaction[4] *= gate`. With the gate off, species
  5 keeps its interaction.
- `jax_hamilton_0d_5species_demo.residual`: `factor = (...) * hill_mask`,
  where `hill_mask = 0` when `K_hill = 0`. With the gate off, species 5's
  interaction is **multiplied by 0**.

Measured to T* = 2 with THETA_DEMO, c* = 25:
- as is, the codes disagree by up to 0.44, and ψ₅ is stuck at 0.5;
- with the gate-off branch set to factor 1, ours equals TMCMC to
  **7.8e−15**.

**Scope.** Only species 5 (P. gingivalis) in five-species runs. The one- and
two-species runs (stages 0–5) are unaffected. **Not fixed yet**: the fix moves
every five-species reference value verified on ANSYS, for example
`t_growth_ecology_substep.dat`'s expected block, so it needs a decision.

**Newton loop.** `newton_step` runs a fixed 6 iterations without a
convergence check. On the first step from ψ = 0.999 at dt = 1e−4 it ends at
residual 1.1e3; 20 iterations reach 1e−7. The transient dies out: over the
stage-5 runs, 50 iterations change α by 2e−6 (0.5074717 → 0.5074696), and
every φ, ψ, γ is unchanged to four digits.

**Time step.**

| TMCMC dt | φ at T* = 2 | min φ₀ |
|---|---|---|
| 0.01 | [0.021 0.375 0.703 0.112 0.192] | **−0.40** (unphysical) |
| 1e−3 | [0.0173 0.0353 0.8942 0.0204 0.0192] | 0.0135 |
| 1e−4 | identical to 1e−3 | 0.0135 |

The calibration's own dt = 0.01 is **unstable for THETA_DEMO at c* = 25**.
Whether that also holds for the calibrated θ actually used is not checked
here. **dt = 1e−3 is converged**, which makes the Gauss-point inner step 10×
cheaper, not 100×. The integrator here is implicit and gives the same answer
at dt = 1e−2, 1e−3 and 1e−4.
