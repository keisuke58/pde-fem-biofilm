# Point model inside the 3D biofilm field — transport/reaction split

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
