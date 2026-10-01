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
two-species runs (stages 0–5) and the Klempt 2026 reproductions (species 5
inactive) are unaffected.

**Fixed 2026-10-01** (decision: fix, re-run the ANSYS baselines at Keio).
With the gate off, the interaction is now left unchanged. The same pattern
was in 9 places: the 0D demo; the 1D and 2D demos/pipelines; the
`JAXFEM/core_hamilton_{1d,1d_nutrient,2d_nutrient}` steppers (two in 2D);
and `hamilton_ode_jax_nsp`. All were fixed together, so 0D, 1D and 2D stay
mutually consistent (`test_pde_uniform_consistency`).

What moved:

| | before | after |
|---|---|---|
| `t_growth_ecology_clsm_phi.dat` (CLSM seed, φ₅ = 0.022) | α 4.972502801e−4, SX −3.00046169e−1, γ 148.88 | α 4.972516806e−4, SX −3.00047014e−1, γ 159.89 (α, SX rel. 2.8e−6) |
| `t_growth_ecology_substep.dat` (default seed, φ₅ = 0) | ψ₅ 0.815483278 | ψ₅ 0.815483971; α, SX unchanged at printed precision |
| other ecology decks | — | unchanged at printed precision |
| ψ-spread sensitivity, α spread at t ≈ 1e−2 | 2.4 % | 1.2 % |

The ANSYS-confirmed pre-fix blocks are kept in the deck headers, with the
new expected values above them. **To be re-run on ANSYS at Keio.**

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
