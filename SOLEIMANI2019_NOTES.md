# Soleimani (2019), visco-elastic growth — the supervisor took the other road

> Soleimani M. **Finite strain visco-elastic growth driven by nutrient
> diffusion: theory, FEM implementation and an application to the biofilm
> growth.** *Computational Mechanics* **64**(5):1289–1301 (2019).
> doi:[10.1007/s00466-019-01708-0](https://doi.org/10.1007/s00466-019-01708-0)

Read 2026-09-30. This is the precedent for this repository's viscoelastic UMAT
(`F = Fe·Fv·Fg`, continuation item D). It turns out to take the opposite
approach, and to carry the integrator that removes a defect we have documented.

---

## 1. He deliberately avoided the three-way split we use

From §2, on viscoelasticity:

> "there are two possibilities ... either to use a linear evolution equation for
> the so-called non-equilibrium part of the stress or a nonlinear one. The
> former is termed *finite linear viscoelasticity* ... The latter is referred as
> the *finite viscoelasticity* ... The second approach entails **extending the
> multiplicative decomposition of the gradient deformation to account for the
> viscous-related part and consequently introducing an additional intermediate
> configuration**. **That is why we adopt the first approach**"

So his kinematics is **`F = Fe·Fg`** (Eq. 1), two-way, and viscosity lives in
the free energy as an internal variable:

```
Ψ(Ĉ, F̂) = Ψ^∞(Ĉ) + γ(Ĉ, F̂)                                    (12)
Ψ^∞(Ĉ) = (λ/2)(Je − 1)² + (μ/2)(tr(Ĉ) − 3 + 2 log Je)          (13)
```

`F̂` is a strain-like internal variable, `Q̂` its work-conjugate stress-like
measure, and `γ` never needs an explicit form — only its derivative is used.

**This repository does the three-way split he chose not to do.** Worth knowing
before proposing it as the shared path.

## 2. His time integrator is unconditionally stable; ours is not

The internal variable follows a 1D Maxwell spring–dashpot ODE

```
Q̂̇ + Q̂/τ = d/dt (∂Ψ^∞/∂Ĉ)                                       (31)
```

solved by a mid-point scheme giving the **recursive exponential update**

```
Q̂_{n+1} = exp(−Δt/τ)·Q̂_n + exp(−Δt/2τ)·( ∂Ψ^∞/∂Ĉ_{n+1} − ∂Ψ^∞/∂Ĉ_n )   (32)
```

`exp(−Δt/τ) ∈ (0, 1]` for every Δt, so this **cannot go unstable however coarse
the step**.

Our viscous update is explicit, `Fv_{n+1} = (I + dt/(2ηJe)·τ)·Fv_n` with the
stress at the old state, and it **flips the sign of the stress at
dt/τ ≈ 0.5** — measured, documented, and the reason the delivered routine has
to inspect the step and refuse. Eq. 32 is the standard fix for exactly that,
and it is in the supervisor's own paper.

**This is the most actionable thing in the three papers read today.**

## 3. Smaller confirmations, each useful

- **`Fg = Jg^{1/3} I = α I`** (Eqs. 2–3), with `α := Jg^{1/3}`. Confirms
  `CITATION_AUDIT.md` F1c: this lineage writes `Fg = αI`, not `Fg = (1+α)I`.
  The two agree with `α_ours = α_paper − 1`.
- **`α̇ = Y·K₁C/(K₂+C)`** (Eq. 4), backward Euler in Eq. 5 — Monod drives the
  growth variable directly.
- **The nutrient is quasi-static**: "the time dependent part is absent in
  Eq. (6) due to the time scale argument". `klempt2024_quantitative.py` drops
  `ċ` for the same reason, so that simplification is the lineage's own, not a
  liberty we took.
- **Eshelby stress, not Mandel**, is the driving force conjugate to the growth
  velocity gradient, "not Mandel stress M̂ **which some authors come up with**"
  — growth is not isochoric, so `J̇g ≠ 0`. A pointed remark worth having in the
  theory chapter.
- Thermodynamic consistency under nutrient-driven growth needs an extra entropy
  source `θS₀` (Eqs. 27–29), because the system is open. Without it
  non-negative dissipation is not guaranteed.

## 4. A caution about our premise, in his own words

> "biofilm growth ... is expected to reach a zero stress-state in large time
> scales compared to the relaxation time. **In such slow growth process, all
> residual stresses are released (dissipated).** That is why biofilms can also
> be modeled as viscous fluid or potential flow, rather than a visco-elastic
> solid."

Our reported runs use **η = 0**, i.e. the unrelaxed elastic stress. That is
defensible for a short horizon, but it is exactly the assumption this sentence
questions at long times, and a reader from this group will know the sentence.

## 5. The Young's modulus problem — and it qualifies the Chu et al. bridge

Table 2 gives the biofilm material properties used here and, via the same
citation, in Klempt 2024:

| | Soleimani 2019 (Table 2) | this repository (`espec_sensitivity.py`) |
|---|---|---|
| Young's modulus | **E = 10 Pa** | per species 10 – 1000 Pa; `E_voigt` ≈ **960 Pa** (CH), **500 Pa** (DH) |
| Poisson ratio | 0.45 | — |
| relaxation time | τ = 0.01 s | — |

So our commensal conditions assume a biofilm about **100× stiffer** than the
lineage we build on, and the dysbiotic ones about 50×. Neither is obviously
wrong — measured biofilm moduli span orders of magnitude, which is what
Billings et al. (2015) reviews — but it is a two-order-of-magnitude open
question sitting under every absolute stress we quote, and `σ ∝ E` is verified
here (2× E gives 2.00× σ).

**Consequence for [`CHU2018_NOTES.md`](CHU2018_NOTES.md).** That note observes
our computed maxima (2.3–13.9 kPa) straddling Chu et al.'s ~5 kPa threshold for
the `rpoH` stress response. At E = 10 Pa the same stresses would be
**0.02–0.3 kPa**, two orders below it. The straddle is therefore a property of
our `E_SPEC` choice as much as of the biology, and must be quoted with that
attached — it is a motivating coincidence of scale, not a result.
