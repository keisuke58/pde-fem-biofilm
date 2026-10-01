# Point model inside the 3D biofilm field

> **Current decision (1 Oct, later the same day): option (i), section 0.**
> The partner's 3D field is left exactly as it is (it *is* Klempt 2024
> Eq. 34–36). The point model decides only the composition. The
> transport/reaction split in sections 1–6 is **superseded**. It stays in
> the fragments (`prop(28) = 5, 6`) and is pre-flighted, but it is not the
> plan.

## 0. Option (i): amount from the 3D field, composition from the point model (`prop(28) = 7`)

- **Amount:** the 3D field (Klempt 2024, USSFin unchanged) gives
  `φ_3D = bio1 + bio2` at the Gauss point.
- **Composition:** the point model (Klempt et al. 2026, a Table 1 case)
  gives only `χᵢ = φᵢ/Σφ` and ψ. Its carried state is rescaled to
  `Σφᵢ = φ_3D`, with `φ₀ = 1 − φ_3D`, before and after each inner
  integration.
- **Growth:** Eq. 36 with `φ_3D`, and `k_α` from Klempt 2024 Table 2
  (1e−3 T*⁻¹).
- **Below `φ_min`** (`prop(29)`) the state is held and the server is not
  called.
- **Initial composition:** `prop(30) = χ₁`.
- **Constants:** two-species constants come from the server,
  `material_server.py --case 2sp_case3`. The case supplies `c* = 100`,
  `α* = 10` and `ηᵢ`; `prop(8:27)` must carry the case's A and b, and the
  server refuses a mismatch. Through the server, case 3 and case 6 reproduce
  `JAXFEM/klempt2026_reproduction` exactly (difference 0.0; case 3
  φ = 0.5945 / 0.3926, against the paper's 0.60 / 0.40).

Pre-flight `tests/test_composition_fragment.py` (all pass). The test has
three points: one growing, one in the void that is filled later, and one at
the seed (φ_3D = 1, clamped).

- **Trace checks:** once per increment, state carried, `Σφᵢ = φ_3D` to
  4e−16, one server call per point and substep, Eq. 36 exact.
- **Replay:** every server call is reproduced bit for bit.
- **Stand-alone scheme:** every point equals
  `composition_reference.reference` (the scheme run on its own for the same
  `φ_3D(t)`) **bit for bit**. This is the exact target for a gradient-free
  region.

### What the rescaling does (measured, Python, macro dt = 0.01, to T = 0.15)

χ₁ at T = 0.05 / 0.15:

| φ_3D held at | case 3 | case 6 |
|---|---|---|
| stand-alone point model (its own amount) | 0.602 / 0.614 | 0.307 / 0.014 |
| 1e−4 (near void) | **0.474 / 0.474, frozen** | 0.473 / 0.473, frozen |
| 0.4 | 0.647 / 0.669 | 0.308 / 0.055 |
| 0.9 | 0.564 / 0.633 | 0.328 / 0.014 |
| 1 − 1e−6 (seed) | 0.490 / 0.473 | 0.048 / 0.014 |

Every case stays finite, but:
1. **The composition depends on the local amount.** At φ_3D = 1 (the seed
   region), case 3's coexistence ends at 0.47 / 0.53, not the paper's
   0.60 / 0.40. The paper's ratio is set during the filling transient, which
   rescaling to a fixed amount removes. Case 6 (one species wins) is robust.
2. **Near the void the composition freezes** at about its starting value.
3. **Exact-solution check (b) needs restating.** `φᵢ = sinh(k_α t)·χᵢ(t)`
   with χ from the *stand-alone* point model is not what option (i)
   computes: with `φ_3D = sinh(1e−3 t)` ≈ 1e−4 the composition stays frozen
   at about 0.47, while the stand-alone χ goes to 0.61. The exact target is
   the stand-alone *rescaled* scheme driven by `φ_3D(t)`, which the
   pre-flight now matches bit for bit. Check (a) (`φ = sinh`,
   `α_K = cosh`) is unaffected; it concerns the 3D field alone.

### `φ_min` (`prop(29)`), defined

`φ_min` gates the point model on the amount. Where `φ_3D < φ_min` at a Gauss
point, the point model is not run that substep: there is no server call and
no rescaling, and composition, ψ and γ are held at their last values (or at
`prop(30)` if the point was never active). Growth (Eq. 36) and the 3D field
are not affected.

The partner's field raises `bio` in the void by `dt·K_LOCAL·locbio` per
substep (locbio ≈ 1), so the void background at the end of the run is about
`K_LOCAL·T_end`. With `K_LOCAL = k_α = 1e−3` (Klempt 2024 Table 2; Eq. 34
and Eq. 36 share k_α) and `T* = 1`, the background is `1e−3`, so **`φ_min =
1e−2` is 10× above it. Confirmed.** Rule: `φ_min ≥ 10·K_LOCAL·T_end`. With the
bring-up value `K_LOCAL = 0.01` and `T = 1.1` the background reaches 0.011,
which would cross 1e−2. In the sweep below, the void stays at χ₁ = 0.5000
(held) for every s.

### `s` (`prop(31)`), and the sensitivity sweep

Option (B) is decided. The point model advances `s·dt` per substep. **s is an
explicit assumption**, because no paper links Klempt 2024's T* to Klempt
2026's time unit. Pre-flighted: the mock with s = 0.5 matches the reference
bit for bit, and α is identical for every s.

**Stress is not swept:** in this mode α depends only on φ_3D, so the stress
is identical for every s by construction, and the test pins it. Only the
composition moves. Script: `ansys_usermat/composition_s_sweep.py`.

χ₁ at T* = 1, deck dt = 0.1, φ_min = 0.01:

| amount φ_3D(t) | case 3, s = 0.05 / 0.15 / 0.5 | case 6, s = 0.05 / 0.15 / 0.5 |
|---|---|---|
| seed (1) | 0.478 / 0.486 / 0.513 | 0.014 / 0.014 / 0.014 |
| interior 0.9 | 0.602 / 0.609 / 0.635 | **0.329** / 0.014 / 0.014 |
| interior 0.4 | 0.646 / 0.669 / 0.676 | **0.307** / 0.049 / 0.014 |
| front 0.05 → 1 | 0.638 / 0.672 / 0.701 | **0.230** / 0.014 / 0.014 |
| void (K_LOCAL·t) | 0.500 (held) | 0.500 (held) |

- **Case 6 (one species wins):** the outcome is set by s. At s = 0.05 the
  point model has not finished: the winner has not emerged, except in the
  seed. From s = 0.15 up, the paper's end state is reached everywhere.
- **Case 3 (coexistence):** χ₁ moves by 0.03–0.06 over the sweep, and depends
  more on the local amount (0.48–0.70) than on s.

### Coupling interval, and the seed region (open)

The deck dt is the coupling interval: rescaling happens once per substep. It
must converge on its own. For case 3, s = 0.15:

| deck dt | seed (φ_3D = 1) | interior 0.9 | interior 0.4 | front |
|---|---|---|---|---|
| 0.1 | 0.486 | 0.609 | 0.669 | 0.672 |
| 0.05 | 0.459 | 0.647 | 0.668 | 0.671 |
| 0.025 | 0.395 | 0.667 | 0.668 | 0.670 |
| 0.0125 | **0.201** | 0.668 | 0.668 | 0.669 |

**Interior and front converge at dt ≤ 0.025. Use deck dt = 0.025:** 40
substeps to T* = 1, with n_sub ≈ 38 at s = 0.15, which is cheap.

**The seed does not converge.** There `φ₀` is pinned at 1e−6, below the
point model's own packing limit (stand-alone case 3 settles at φ₀ ≈ 0.013).
Each inner integration pushes φ₀ up off its barrier and each rescale pushes
it back, so the finer the coupling, the more that fight dominates. Capping
the amount fed to the point model does not cure it above φ₀ ≈ 0.05:

| cap on φ_3D | dt 0.1 / 0.05 / 0.025 / 0.0125 |
|---|---|
| 1 − 1e−6 | 0.486 / 0.459 / 0.395 / 0.201 |
| 0.99 | 0.504 / 0.499 / 0.488 / 0.465 |
| 0.95 | 0.562 / 0.598 / 0.638 / 0.664 |
| 0.9 (the interior row) | converged, 0.668 |

**A decision is needed** on how composition is defined where the 3D field
is fully packed:
- (a) evaluate it at a capped amount (φ_cap ≈ 0.9; converged). This is one
  more stated assumption.
- (b) in packed regions, take χ from the free-running point model, with no
  rescaling.

Either way, the seed numbers above (and the ANSYS smoke test on element
220, which is in the seed) are coupling-interval artefacts, not results.

### ANSYS smoke test, IKMHIWI03, 1 Oct (prop(28) = 7)

One element (220, in the seed), v222, `-np 1`, 30 constants (before s
existed, so s = 1), server `--case 2sp_case3`, `DELTIM 0.01`, end time 0.16.

- 0 errors.
- `check_comp_trace`: once / carried / amount / one_call / eq36 all True;
  replay True (worst 0.0).
- The whole trajectory equals `composition_reference.reference` (max
  difference **0.0**).
- χ₁ ends at 0.4708, ψ at 0.9835 / 0.9800.

So the Fortran path is exact on real ANSYS. The seed caveat above applies
to the value 0.4708 itself.

### Time: matching the papers

- **Klempt 2024** normalises its own time, `T* = t/t_ref ∈ [0, 1]`, with
  10³ nominal substeps.
- **Klempt 2026** runs 500–1500 steps of 1e−4, which is 0.05–0.15 in *its*
  time unit.
- **PAMM 2023** (Klempt, Soleimani, Junker) plots 0–0.2.

Nothing links the 2024 and 2026 units: no `t_ref` or day count is given in
either paper. Two ways to proceed:

- **(A) Deck end = 0.15**, the point model at its own pace. The composition
  reaches the paper's horizon, but the 3D growth covers only 15 % of Klempt
  2024's horizon: with `k_α = 1e−3`, α ≤ 1.5e−4, so there is hardly any
  stress. Cheap: 15 substeps of 0.01, `n_sub = 100`.
- **(B) Deck end T* = 1** (Klempt 2024's own horizon), with the point model
  run on a scaled clock, `dt_pm = s·dt` and `s = 0.15`. Each paper keeps
  its own horizon, and `s` is the single, explicit modelling assumption
  linking them. It is not built yet: it would be one more property,
  `prop(31) = s`.

**Recommendation: (B).** It keeps both published parameter sets unchanged
and puts the unknown into one named factor, the same unknown as the
T* ↔ days question. Either way, `k_α = 1e−3` with T* ≤ 1 gives
α ≤ 1e−3. The stresses will be small: Klempt 2024's own growth is slow, and
that is a property of the paper's parameters, not of the coupling.

### For the thesis: composition does not feed back into the stress

Suggested wording:

> The amount of biofilm and its growth follow Klempt et al. (2024); the
> species composition is computed from the multi-species point model of
> Klempt et al. (2026) and is a one-way output: it does not alter the
> mechanical response. This is a property of the 2024 formulation, in
> which no material parameter depends on the species. In Klempt, Soleimani
> and Junker (PAMM 2023), species identity enters only through the
> front-growth factor R_s and the consumption g, both stated to depend on
> the type of microorganism; with these taken equal for all species, as
> here, the composition cannot reach the stress. A composition-dependent
> R_s = Σχᵢ R_s,i would be the natural two-way extension.

The PAMM paper also supports two existing choices:
- the nutrient diffuses "instantly" compared with growth, which supports a
  quasi-static or constant c;
- the hydrostatic pressure is compressive inside the biofilm, with a ring of
  tension at the 0 < φ < 1 edge. This is a qualitative check on our whole-
  model stress, where the ring SEQV is the largest.

### Run sheet (IKMHIWI03), option (i)

1. `git pull`. Re-paste both fragments; keep `split_rates.f` in the build
   (its cache is used) and `USE biofilm_split`. **USSFin is not changed.**
2. Server: `material_server.py --case 2sp_case3`. This sets two species,
   c* = 100, α* = 10 and η = (1, 2).
3. `TB,USER` with 31 constants (`prop(31) = s`, e.g. 0.15; deck
   `DELTIM 0.025` to T* = 1):
   - `prop(6) = 0`, `prop(28) = 7`;
   - `prop(7) = 1e−3` (k_α, Klempt 2024 Table 2);
   - `prop(8) = 1`, `prop(9) = 1`, `prop(10) = 1`, `prop(11) = 0`,
     `prop(12) = 0` (case 3's A and b); the rest of `prop(8:27)` = 0;
   - `prop(29) = φ_min` (e.g. 0.01);
   - `prop(30) = 0.5` (case 3 starts 0.2 / 0.2).
4. Gradient-free check first: in a region where `φ_3D = sinh(k_α t)`
   (check (a)), the composition must equal
   `composition_reference.reference(φ_3D series, ...)`.
5. Judge with
   `composition_reference.check_comp_trace(read_comp_trace('comp_trace.csv'), 1e-3, theta, hp)`,
   where `theta`/`hp` are `material_server.ECOLOGY_CASE` after
   `set_case('2sp_case3')`.

---

Decision (1 Oct): for one and two species, connect the calibrated point model
(Klempt et al. 2026, Hamilton ecology ODE) to the partner's 3D field by
**operator splitting**:
- **transport** in the partner's field (USSFin);
- **reaction** from the point model at each Gauss point;
- **growth and mechanics** by Klempt 2024: `α̇ = k_α φ` with `φ = Σφᵢ`.

The identification is Klempt's `φ = Σᵢ φᵢ = 1 − φ₀`: the point model's void
is Klempt's void `1 − φ`.

Pre-flight: `tests/test_split_coupling_fragment.py`. The fragments run in the
mock usermat, linked against the real bridge and a live `material_server`,
together with a stand-in for USSFin (explicit zero-flux diffusion + `dt·Rᵢ`).
All checks pass for one and two species.

## 1. What Klempt 2024 says about each term

Eq. 34 (Eq. 30 divided by η, with the mechanical term dropped) is

    φ̇ = β∇²φ + k_α α − r c/(k+c) · ‖∇φ‖ n∇φ·n∇c

Taking the terms one at a time:

- **`+ k_α α` is local**, with no gradient. This is the partner's
  `+K_LOCAL·bioloc` (`bioloc = α_K`, starting at 1). It acts at every point,
  void included, which explains bio1 = 0, 0.001, 0.002 … in empty elements.
  It is the paper's own term, not a bug. **This is the term the reaction
  replaces.**
- **The Monod term is transport, not reaction.**
  `‖∇φ‖ n∇φ·n∇c = ∇φ·∇c/|∇c|`, so the term is `v·∇φ` with
  `v = r c/(k+c) · ∇c/|∇c|`. That is advection of the biofilm front up the
  nutrient gradient. The paper says so itself: it "promotes growth in the
  interface of biofilm and void … no non-local growth may happen in a region
  with homogeneous biofilm". **Recommendation: keep `−Ori·Growth` in the
  transport step.** Two notes:
  - the partner weights it with `|∇²bio|` where the paper has `‖∇φ‖`;
  - with a constant nutrient (`c* = 25`, point 4 of the decision) `∇c = 0`,
    so the paper's form of the term is zero. It survives only while the
    partner's own nutrient field is solved.
- **`β∇²φ`** stays transport, together with the penalty term.
- **Growth law.** Eq. 36, `α̇ = k_α φ`, is used as decided. Table 1 of the
  paper writes the local equation differently,
  `(α_{n+1}−α_n)/((1+α_{n+1})Δt) = k_α/η_α · (φ_{n+1}−φ_n)/Δt`, which is not
  a discretisation of Eq. 36. Recorded here, not resolved.

## 2. The scheme, per substep (Lie splitting, explicit like USSFin)

1. **Reaction, in the usermat**, at each Gauss point, from the converged
   field `bioᵢ` of the previous substep:
   - set the point-model state to `φᵢ := bioᵢ` and `φ₀ := 1 − Σφᵢ`;
   - carry ψ and γ in `ustatev(72:83)`;
   - advance by `dt` with `n_sub = ⌈dt / 1e−4⌉`;
   - `Rᵢ = (φᵢ_after − φᵢ_before) / dt`.

   It is computed once per substep: ANSYS hands back the same converged
   state on every iteration, and a cache makes it one server call per point
   and substep. Rᵢ goes to USSFin through the module `biofilm_split`
   (`callsite/split_rates.f`), indexed by the NEM point `ID = (elem−1)·8 +
   kDomIntPt`. No new pool offsets are needed.
2. **Transport, in USSFin**:

       bioᵢ_n = bioᵢ_{n−1} + dt · ( Rᵢ + Betaᵢ·∇²bioᵢ + Penᵢ − Oriᵢ·Growthᵢ )

   This is the existing update with `+K_LOCALᵢ·biolocᵢ` replaced by `Rᵢ`.
   Drop `−Ori·Growth` as well only if the decision is to treat it as
   reaction, which section 1 argues against.
3. **Growth**: `α_n = α_{n−1} + k_α · Σbioᵢ · dt` (Eq. 36, explicit), in
   `ustatev(84)`, as in the one- and two-species φ modes. There is no ψ
   weighting.

Time: deck time is T*. Nutrient: constant `c* = 25` inside the point model.
A local `c` is a later option; it needs a fourth argument to the bridge.

## 3. The void: the point model makes biomass from nothing

Measured (one species, `dt = 0.1`): started from `φ = 0`, the point model
gives `φ = 0.117` after one step, and from `φ = 1e−8` it gives `0.121`.
The barrier term makes `φ = 0` a non-equilibrium. Applying the reaction at
every point would fill the void everywhere within a substep, a much
stronger form of the `K_LOCAL` effect above.

So **`prop(29) = φ_min`**: below it there is no reaction, the point is not
sent to the server, and the void is reached by transport only. The value is
a modelling decision, still open. 0 means the reaction runs everywhere.

## 4. Checks (pre-flight, all pass)

| | check | result |
|---|---|---|
| a | spatially uniform field, n = 1 and 2 | stays uniform; equals the stand-alone point model to < 1e−12 |
| b | mass | each substep, `Σbio_n − Σbio_{n−1} = dt·ΣRᵢ` to < 1e−13 (transport conserves) |
| c | reaction off (`φ_min` > 1) | transport only: mass conserved to 1e−14, no server call, α follows Eq. 36 exactly (`check_trace`) |
| d | below `φ_min` | `Rᵢ = 0` exactly; reaction above it |
| e | cost | one server call per point and substep (three equilibrium iterations: miss, hit, hit) |

"One species reduces to Klempt 2024" holds for growth and mechanics (c):
Eq. 36 is exact. The φ-dynamics deliberately do not reduce to it, because
the local source `k_α α` is replaced by the point model's reaction.

## 5. Run sheet (IKMHIWI03)

1. `git pull`. Copy `callsite/split_rates.f` into the wired folder.
2. In `Usermat_P21-V21_v222.F`, add `USE biofilm_split` next to
   `USE biofilm_py_bridge`, and re-paste both fragments.
3. In `Ussfin_P21-V21_Conection_Test.F`:
   - add `USE biofilm_split` to the routine that does the bio update;
   - at each NEM point,
     `CALL SPLIT_GET(ID, R1, R2, LD, IS)`;
   - replace `+K_LOCAL1*bioloc1` with `+R1` (and the same for species 2);
   - keep `−Ori·Growth`, `Beta·Lap` and `Pen`.

   `LD = −1` means no rate was written for that point (not a split
   material); R is then 0. Leave the `bioloc` update in place; it is unused
   by the split growth.
4. Build order: `usermat_py_hook.f`, then `split_rates.f`, then the rest.
   Delete stale `.mod`, `.obj` and `ANSYS.exe` first.
5. Server: `material_server.py --active-species 1` (or `2`).
6. Material: `TB,USER` with **29** constants:
   - `prop(6) = 0`;
   - `prop(28) = 5` (one species) or `6` (two);
   - `prop(7) = k_α`;
   - θ in `prop(8:27)` (a₁₁ = `prop(8)`, b₁ = `prop(11)`; two species add
     a₁₂ = `prop(9)`, a₂₂ = `prop(10)`, b₂ = `prop(12)`);
   - `prop(29) = φ_min`.
7. Bring-up order:
   1. **reaction off** (`φ_min = 2`): the field must match the partner's
      run with `K_LOCAL = 0`;
   2. **uniform** `bio1 = 0.3` everywhere, `Beta = 0`: every point must
      follow the stand-alone point model (`ecology_substeps`);
   3. the real initial condition.

   `-np 1`.
8. Cost: 4096 points × about 0.07 s per server call (1000 inner steps) is
   about 5 min per substep. Reaching T* ≈ 15–25 (where the calibration ODE
   settles) at `dt = 0.1` is 150–250 substeps.
   `ODE_TMCMC_CROSSCHECK.md`: an inner step of 1e−3 is converged, which
   would make that 10× cheaper. The step is the `PM_DTMAX` parameter in the
   fragment.

Check after each run: `split_trace.csv` has `hit = 0` exactly once per point
and substep. Mass per substep: `Σbio_n − Σbio_{n−1} = dt·ΣR` over the NEM
points, which needs a pool dump.

## 6. For the partner (no internal names)

> We would like to replace the local source term in the biofilm update
> (`K_LOCAL·bioloc`, Klempt 2024's `k_α α` in Eq. 34) by a reaction rate
> computed at each Gauss point from our calibrated multi-species point
> model, and keep diffusion, penalty and the Monod front term as they are.
> Growth stays Eq. 36 with φ the total biofilm fraction. Two questions:
> (1) we read the Monod term `‖∇φ‖ n∇φ·n∇c` as advection of the front
> along the nutrient gradient (transport), not as a local reaction — do you
> agree, and is `|∇²bio|` in place of `‖∇φ‖` intended? (2) The point model
> does not keep φ = 0 at zero (it creates biomass from nothing), so we
> switch the reaction off below a threshold φ_min and let transport carry
> the biofilm into the void — is that acceptable to you, or would you
> prefer a different treatment?
