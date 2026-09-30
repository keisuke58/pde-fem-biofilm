# What this work may and may not claim, and what to do next

State of play, 2026-09-30. Several headline claims changed in the last two days,
and two of the changes reverse earlier corrections. This is the one page to read
before writing a results section or walking into a meeting; the detail is in the
linked notes.

---

## 1. Claims that changed — do not state these the old way

### 1.1 "The clinical condition drives the stress" — **withdrawn**

The two cylinder layers were bonded with `NUMMRG,NODE`, which merges only
coincident nodes; in the 4-condition ecology deck that came to **four corner
nodes**. `VGLUE` bonds all 651 and the mesh then converges at second order.

| | before | bonded re-run |
|---|---|---|
| stress spread across the four conditions | 9.1 % | **0.4 %** |
| growth `α` spread | 1.106× | 1.006× |
| growth layer / substrate von Mises | 15.2× | **6.3×** |

The "two-lobe buckling mode" was the layer lifting between bond lines, and the
finer meshes failing at their first substep had the same cause. This is the
**second** correction to the same claim: it went "position" → "condition" →
neither. See `ansys_usermat/apdl/README.md`.

### 1.2 The reason is the **simplex**, not multi-species modelling

The conditions cannot differ much because `φ_tot = Σφᵢψᵢ` is a single total and
CLSM compositions are fractions summing to 1, so the differences cancel *inside
the sum*. Soleimani et al. (2023) — the supervisor's own multi-species paper —
uses **independent Allen–Cahn phase fields** with no such constraint, and there
the species can both be present at one point.

> **State it as: the degeneracy belongs to the simplex formulation, not to
> multi-species modelling.** That is both more precise and more defensible.

[`SOLEIMANI2023_NOTES.md`](SOLEIMANI2023_NOTES.md) §1.

### 1.3 It survives spatial coupling — which is the stronger result

The 0D collapse is a statement about the pointwise setting the thesis exists to
leave behind, so it was tested in 2D with per-species diffusivities and
per-species nutrient consumption — the concrete way `φ_tot` could stop standing
for the composition.

| arm | condition spread in `α` |
|---|---|
| control (no coupling between nodes) | 0.0338 % |
| full (species diffusion + nutrient consumption) | **0.0450 %** |

and this time the mechanism genuinely ran: mean `c` fell to 0.304 and the
species separated (0.1265, against exactly 0 in earlier attempts). Spatial
coupling widens the gap 1.33×, and nothing like the withdrawn 9.1 % returns.

`JAXFEM/condition_spread_2d.py`.

### 1.4 The Klempt 2024 "agreement" must not be quoted

`klempt2024_quantitative.py` reported Fig. 4 "within 0.17 / 0.10". **That is a
variant the paper does not use.** Reading the paper settles it: Eq. 34's growth
is a signed dot product with no absolute value, and Eq. 35 with Eq. 24 fixes
consumption as `g φ`, zeroth order in `c`. The two variants that fit best are
the two the paper does not use, and the paper's own combination is the worst fit
— 0.030 against 1.000 at `T* = 0.20`.

[`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md).

### 1.5 The ~5 kPa "straddle" is half an `E_SPEC` artefact

Chu et al. (2018) measure that confined colonies generate ~10 kPa and that the
`rpoH` stress response turns on near **~5 kPa**, leading to biofilm formation
and antibiotic tolerance. Our computed maxima are 2.3–13.9 kPa, which straddles
it — **but** Soleimani (2019) Table 2, the source Klempt 2024 cites for its own
properties, uses **E = 10 Pa**, while we assume 10–1000 Pa per species
(`E_voigt` ≈ 960 Pa commensal). Since `σ ∝ E` is verified here, at 10 Pa the
same stresses are 0.02–0.3 kPa — two orders *below* the threshold.

> **Quote the straddle as a motivating coincidence of scale, never as a result.**

[`CHU2018_NOTES.md`](CHU2018_NOTES.md), [`SOLEIMANI2019_NOTES.md`](SOLEIMANI2019_NOTES.md) §5.

---

## 2. Actionable now, in order of value

1. **Replace the viscous integrator.** Soleimani (2019) Eq. 32 is a recursive
   exponential update whose factor `exp(−Δt/τ) ∈ (0,1]` for any step, so it
   cannot destabilise. Our explicit `Fv` update **flips the stress sign at
   `dt/τ ≈ 0.5`** — measured, and the reason the delivered routine has to
   inspect the step and refuse. The standard fix is in the supervisor's own
   paper. → [`SOLEIMANI2019_NOTES.md`](SOLEIMANI2019_NOTES.md) §2.
2. **Run `tier2b_real/tie_coverage_check.py` on `tier2b_real.inp`.** The
   tooth/implant `*TIE` coverage has never been checked, and `ADJUST=NO` with
   hand-set tolerances (0.5 / 0.6 / 1.0 / 2.8 mm) drops far slave nodes with
   only a warning — the same silent-partial-bond failure as §1.1. Needs the
   machine that has the `.inp`; takes seconds.
3. **Cite Soleimani et al. (2023).** The nearest precedent to this thesis and it
   appeared nowhere. → [`READING_GAPS.md`](READING_GAPS.md).
4. **Close the UserElement question.** `THESIS_ASSIGNMENT.md` §4.1 says not to
   implement until it is settled. Soleimani 2023 settles it: a multi-field ANSYS
   user element with **four scalar DOFs per node**.
5. **Ask the experimental side for per-species LIVE/DEAD.** That group already
   does SYTO9/PI CLSM, which is a bounded viability measure, and
   `JAXFEM/psi_spread_sensitivity.py` says the species would need a spread of
   about **17 %** before composition reached the stress. That turns "no data for
   ψ" into a specific request.

---

## 3. Open, with a plan

| item | status |
|---|---|
| **Klempt 2024 Fig. 7** | the paper's own equations are the worst fit. Seed and nutrient cleared by measurement, constants cleared by reading. `JAXFEM/klempt2024_resolution.py` is testing resolution — supported by the same group stating phase-field needs 0.5 µm while Klempt 2024 runs 1 µm in 3D. The residue is a question for the authors, and **Meisam is a co-author**. |
| **`E_SPEC`** | uncertain across two orders (10 Pa in the lineage against ~960 Pa here). Carries about half the headline ratio, and all of §1.5. |
| **`k_alpha` time unit** | the real blocker on the partner's deck; `α` reaches 4.7 (a 5.7× stretch) at `TIME INC = 0.1`, so the distortion is the growth itself. Question 4 on the 10/1 agenda. |
| **Partner routine's call site** | still refuses `dt > DT_ECO_MAX` before the hook is reached; not in this repository. `ansys_usermat/apdl/NSUB_WIRING.md`. |
| **Tooth/implant geometry** | the roadmap's open item, gated on item 2 above. |

---

## 4. How the errors were caught, because it keeps mattering

Three wrong claims were found today, and **none of them by reading code**.

- "Spatial structure widens the condition gap 2.4×" — false. It compared a 2D
  number against a 0D figure quoted from a different script with a different
  horizon. An internal control gave 1.4736 % against 1.4736 %, identical.
- "The 2D study measures spatial effects" — it did not. It mirrored
  `run_simulation()`, which passes a *constant* nutrient into the reaction and
  never feeds its own field back, so the grid was 0D replicated 64 times.
  Control equalling full to every digit was the tell.
- "This needs a bigger machine" — it did not. A Python macro loop re-entered XLA
  compilation each iteration until LLVM could not map memory, which reads as an
  environment limit with 15 GB free.

> **Put the control in the same script, with the same horizon and constants. A
> number carried from another file is not a baseline.**

The same rule caught the bonding bug: it surfaced from a mesh-convergence study,
i.e. from a check that *could fail*. Checks that only pass do not raise rigor.
