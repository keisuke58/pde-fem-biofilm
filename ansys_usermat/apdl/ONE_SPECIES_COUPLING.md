# One-species coupling into the partner's element

Decided 2026-10-01 in the meeting with Oliver: couple **one species fully**
first, then extend to two. Main results follow Klempt et al. (2024); nothing
is added to their formulation.

## Scope decision (2026-10-01)

**The thesis completes the coupling for one and two species, and stops
there.** Five species move to the Keio continuation.

- Oliver's element is a two-species code (`Bio1`, `Bio2`, `Interaction12/21`),
  so one and two species run inside what he gave us, with a shared variable
  to verify against (`φ`, or `φ₁ + φ₂`).
- From three species on, the code itself has to be generalised: his
  interaction term is multiplicative, the paper's is additive, so there is no
  shared variable space to check one against the other, and the five-species
  trial already showed a coupling ceiling between 0.0005 and 0.0006. That is
  the "one to two months" that made a full merge future work.
- "Complete" for each of the two means: stages 0–2 of the bring-up pass on
  real ANSYS, the alpha history matches `one_species_reference`, and the
  stresses are reported in the Klempt setting.

## Why one species first

Every obstacle that made the five-species merge "future work" goes away:

| obstacle with five species | with one species |
|---|---|
| no shared variable space (their 2 species / 2 nutrients, our 5 species) | one `φ` on both sides |
| who computes `φ` — two ecology models would both evolve it | their NEM solves `φ` and `c`; we only grow and stress |
| nothing to verify against | the chain reduces to closed form, and `JAXFEM/klempt2024_quantitative.py` is single-species Klempt 2024 |
| a Python server per Gauss point | none: Eq. 36 is one line of Fortran |

## The law

Klempt et al. (2024) Eq. 36, read off the bundled PDF:

    α̇ − k_α φ = 0

Their `Fg = α_K I` with `α_K(0) = 1`; ours `Fg = (1+α) I` with `α(0) = 0`
(Soleimani, Haverich & Wriggers 2021, Eq. 15). `α = α_K − 1`, a constant
shift, so it is the same equation.

## What is in this repository (verified, no ANSYS needed)

- `ansys_usermat/growth_from_phi.f` — `BIOFILM_ALPHA_FROM_PHI`, explicit:
  `α_{n+1} = α_n + k_α φ Δt`. No cap, no clamp, no gate.
- `ansys_usermat/one_species_reference.py` — the same, in Python, written
  independently; also predicts `SVAR(84)` from a printed `φ` history.
- `ansys_usermat/crosscheck/one_species_driver.f` — the chain at one Gauss
  point: `φ → BIOFILM_ALPHA_FROM_PHI → BIOFILM_GROWTH_VISCO_V01 → σ`, with
  `η = 0`, `C01 = 0`, neo-Hookean.
- `tests/test_one_species_coupling.py` — 7 tests: Eq. 36 for constant `φ`;
  bit-identical to the Python reference on a varying `φ`; no biofilm, no
  growth; nothing added to Eq. 36; free growth gives zero stress; constrained
  growth matches closed form (plus the documented spherical term,
  `DEVIATOR_SCALING_FINDING.md`).

## The call-site change (apply on IKMHIWI03, `F:\biofilm_upf_wired`)

The partner's file stays off this public repository; this is the recipe. It is
superseded by the run sheet below (switch `prop(28)`, target
`Sbio_GrowthConst` in v222.F); kept for the reasoning:

```fortran
      else if (prop(6) .gt. 1.5d0) then
C       one-species coupling: the partner's local phi drives growth,
C       Klempt 2024 Eq. 36. ustatev(84) is the CONVERGED alpha -- ANSYS
C       hands it back unchanged on every equilibrium iteration, so growth
C       is counted once per increment, not once per iteration.
        PHI_LOC = <local biofilm phi from the pool -- see question 1>
        call BIOFILM_ALPHA_FROM_PHI(ustatev(84), PHI_LOC, prop(7),
     &                              dTime, ALPHA_NEW)
        sGrowth     = ALPHA_NEW
        ustatev(84) = ALPHA_NEW
```

Declare `double precision PHI_LOC, ALPHA_NEW` explicitly. Material constants
for the Klempt setting: `prop(2) = 0` (η), `prop(3) = 0` (C01 ratio),
`prop(4) = 0` (neo-Hookean). Add `growth_from_phi.f` to the compile list; it
has no module dependency, so build order does not matter for it. No material
server is needed in this mode.

## Questions that must be answered first (one answered)

1. **Which pool variable is the local `φ` at the Gauss point.** Oliver (1 Oct):
   the macroscopic `φ` is defined at each quadrature point, not at the nodes,
   and is the fraction of the maximum local biofilm density, in `[0, 1]`.
   Still not said which of
   `Sdp_bio1_n` or `Sdp_locbio1_n`? Their file uses both (`Sdp_sumBio` from
   the first, `Sdp_sumLocal` from the second), and the code does not say which
   is the Gauss-point value.
2. ~~**`Sdp_sumBio = Sdp_bio1_n + Sdp_bio1_n`.**~~ **Discussed 2026-10-01
   (meeting transcript): Oliver assumes it is a typo** for
   `Sdp_bio1_n + Sdp_bio2_n` **but has not verified it**, and asked for it to
   go in the meeting summary. He also pointed out that in the minimal working
   example the starting values of biofilm 2 are zero ("a leftover" from
   Felix's investigations), so today it has no influence. **The decision to
   fix it in our working copy is ours**, not his; tell him when it is in.
   `apply_partner_patches.py` makes the one-line change on IKMHIWI03 and
   refuses unless the line occurs exactly once:

   ```bat
   python ansys_usermat\apdl\apply_partner_patches.py ^
       F:\biofilm_upf_wired\Usermat_P21-V21_Conection_Test.F
   ```

   With biofilm 2 at zero, `Sdp_sumBio` is then exactly `φ`. The fix also
   changes the stiffness blend of the original AceGen path (`prop(1) = 0`),
   so note it beside any earlier run of that path.
3. **The physical time unit of the deck**, which sets `k_α` (`prop(7)`).
   **Oliver does not know it either** ("0.1 second or 0.1 minute, 0.1 hour, I
   don't know"). Until it is fixed, choose `k_α` by its result instead of by a
   rate: `one_species_reference.k_alpha_for_target(α_target, φ_max,
   TIME_total)` keeps the total growth near a small, safe `α_target`, and the
   thesis reports it that way. `k_α = 50` (chosen for millisecond decks) takes
   `α` to about 4.7 per increment at `TIME INC = 0.1`.

## What Oliver said about the bridge (meeting, 1 Oct)

- **Target ANSYS 2024** — the version Felix gave him, "likely the most
  up-to-date regarding the biofilm implementation".
- **Keep to Klempt 2024; add nothing** — no viscosity, no Mooney-Rivlin. "It is
  always good when you have a paper that you can reference." For a start he
  even suggested small deformations / plain elasticity.
- **The two models stay separate and are bridged.** The 3D model has one
  homogeneous `φ` and deliberately ignores what it is made of; the point model
  resolves species. At each quadrature point, run the point model with many
  small inner steps until one outer step is covered, and carry its end state
  into the next outer step so there are no jumps. That is exactly the `n_sub`
  sub-stepping already verified in `usermat_biofilm.f`.
- **Open on his side too — `φ₀`.** The point model's `φ₀` is an empty volume
  inside a simplex; the outer `φ` is a fraction of the available volume. How
  the two relate is not settled. With **one** inner species the unexplained
  `φ₀` becomes larger, so he suggested considering **two** species. Start
  simple, get it running, then improve.

## Staged bring-up (answers incomplete: confirm each stage by running)

Each stage has a pass criterion; do not move on until it holds. Stages 1–3
add a one-row-per-call trace (`elem, ip, ldstep, isubst, dtime, bio1,
locbio1, alpha_n, alpha_new`) written from the call site for two or three
elements at integration point 1, run with `-np 1`;
`one_species_reference.check_trace` reads it and returns the checks below.

| stage | what | passes when |
|---|---|---|
| 0 | typo patch only, `prop(1) = 0` | the minimal working example runs to the end with 0 errors, as before |
| 1 | `prop(28) = 1` with `k_α = 0` | 0 errors, `α` stays 0 everywhere, and the trace shows which of `bio1` / `locbio1` differs between points and changes in time — that one is the Gauss-point `φ` |
| 2 | small `k_α` (`α_target` ~ 0.02) | `check_trace`: `once_per_increment`, `eq36` and `carried` all true |
| 3 | raise `α_target` step by step | runs to the end; record where it stops being stable |
| 4 | two species: `φ = φ₁ + φ₂` | same checks, with biofilm 2 seeded non-zero |

## Run sheet — one and two species (IKMHIWI03, working copy)

**Target file: `Usermat_P21-V21_v222.F`** — the one the 9/7 `ANSYS.exe` was
built from (it has the `.obj`, the `prop(6)` ecology branch and
`ustatev(84)`). `Usermat_P21-V21_Conection_Test.F` is not in the build; the
typo fix already applied to it on 2026-10-01 is harmless and has a backup.
Found by the IKMHIWI03 session on first contact with the real folder: the
first version of this sheet assumed a `sGrowth` variable and a free
`prop(6)`, neither of which exists in v222.F.

The fragments are pre-flighted here as strict 72-column fixed form under
`IMPLICIT NONE`, with bounds checking, inside a mock usermat driven the way
ANSYS drives one (`tests/test_partner_callsite_fragment.py`). What remains
untested is only their position in the real file.

1. `git pull` this branch on IKMHIWI03.
2. Typo fix on the target file (refuses unless the line occurs exactly once
   and `Sdp_bio2_n` is defined elsewhere in the file):
   `python ansys_usermat\apdl\apply_partner_patches.py F:\biofilm_upf_wired\Usermat_P21-V21_v222.F`
3. Paste `apdl/callsite/phi_mode_decl.inc` with the local declarations, and
   `apdl/callsite/phi_mode_exec.inc` **after the ecology block has finished
   setting `Sbio_GrowthConst`, immediately before `CALL
   BIOFILM_GROWTH_VISCO_V01`**. The switch is `prop(28)`; `prop(6)` is left
   to the ecology mode. Check these names exist in scope first: `nProp`,
   `prop`, `ustatev`, `dTime`, `elemId`, `kDomIntPt`, `ldstep`, `isubst`,
   `Sdp_bio1_n`, `Sdp_bio2_n`, `Sdp_locbio1_n`, `Sdp_locbio2_n`,
   `Sbio_GrowthConst`.
4. Build. Copy `ansys_usermat\growth_from_phi.f` into the folder, then
   (`-Filter` takes one pattern, and two files define `usermat`, so the list
   is built explicitly):

   ```powershell
   Copy-Item ansys_usermat\growth_from_phi.f F:\biofilm_upf_wired\
   $skip = 'Usermat_P21-V21_Conection_Test.F'
   $hook = 'usermat_py_hook.f'
   $rest = Get-ChildItem F:\biofilm_upf_wired -File |
       Where-Object { $_.Extension -in '.f','.F' -and
                      $_.Name -ne $skip -and $_.Name -ne $hook } |
       ForEach-Object Name
   .\ansys_usermat\apdl\link_v222.ps1 -WorkDir F:\biofilm_upf_wired `
       -Sources (@($hook) + $rest)
   ```

   `usermat_py_hook.f` goes **first** (v222.F `use`s its module; a stale
   `.mod` is the trap `CLAUDE.md` records). **Use the copy already in
   `F:\biofilm_upf_wired`, not the repository's current one**: the
   repository's hook gained `n_sub` and `phi_int` in `bdddf8a`, while
   v222.F still makes the old five-argument ecology call, so mixing them
   fails to compile. The `n_sub` call-site patch (NSUB_WIRING.md) is a
   separate step and not needed for the phi mode.
5. Material constants — `TB,USER,<mat>,1,28` (28 constants), Klempt setting:

   | prop | value | meaning |
   |---|---|---|
   | 1 | 1 | use the growth law |
   | 2 | 0 | η — viscosity off |
   | 3 | 0 | C01 ratio — neo-Hookean |
   | 4 | 0 | mtype — neo-Hookean |
   | 5 | 0 | constant growth (overwritten in phi mode) |
   | 6 | **0** | ecology mode **off** |
   | 7 | `k_α` | stage 1: `0`; later: `k_alpha_from_trace(<stage-1 trace>, α_target)` |
   | 8–27 | 0 | ecology θ, unused |
   | 28 | 1, 2 or 3 | phi mode: one species (`φ = bio1`), two (`φ = bio1 + bio2`), or 3 = the partner's own α (`Sdp_sumLocal − 1`) |

6. **Two species only:** seed biofilm 2 — `sGdp_Bio2start` is zero in the
   minimal working example. Oliver noted biofilm 1 cannot grow into biofilm-2
   nodes unless their boundary value is released there.
7. Run with `-np 1` (the trace is one file).
8. Judge it: `python -c "import sys; sys.path.insert(0,'ansys_usermat'); import one_species_reference as r; print(r.check_trace(r.read_trace('phi_trace.csv'), K_ALPHA))"`
   with the `k_α` used. Pass = `once_per_increment`, `eq36`, `carried` all
   `True`, plus `NUMBER OF ERROR MESSAGES = 0` and the load step complete.
   The `phi_bio1` / `phi_locbio1` entries say which variable varies between
   points and in time — that is the Gauss-point `φ`.

## The point model with one or two species (ecology bridge)

Separate from the phi mode above. When the inner point model (the Hamilton
ecology ODE) runs at the Gauss point, a one- or two-species run uses the
**unchanged five-species interface** — `g(12)`, `theta(20)`,
`ustatev(72:83)`, `prop(8:27)` — with species masked off in the server:

    python ansys_usermat/coupling/material_server.py --active-species 2

**Masked, not zeroed.** The model clips every `φ` up to `1e-10` each step,
so a species merely started at zero is revived and grows (measured: to
0.036 in 2000 steps, `γ` off by 95). Masked, it stays exactly zero, and the
masked five-species model reproduces `JAXFEM/hamilton_ode_jax_nsp.py` at
**n = 1 to 1.7e-12 and n = 2 to 3.7e-12** — it *is* the n-species model
(`tests/test_ecology_active_species.py`). The five-species path is
bit-identical to before. The five-species default the Fortran seeds is made
consistent on the first call (switched-off `φ`, `ψ` zeroed;
`φ₀ = 1 − Σ active φ`).

Parameters keep their five-species positions (`prop(8:27) = theta(1:20)`):

| | `A` | `b` |
|---|---|---|
| n = 1 | `prop(8) = a₁₁` | `prop(23) = b₁` |
| n = 2 | `prop(8) = a₁₁`, `prop(9) = a₁₂`, `prop(10) = a₂₂` | `prop(23) = b₁`, `prop(24) = b₂` |

All other `prop(8:27)` entries are ignored by the mask; set them to 0.

Still needed for this path in v222.F: the `n_sub` call-site patch
(`NSUB_WIRING.md`), because the macro step 0.1 is 1000 times the ODE's
limit. At about 0.07 s per Gauss-point call for 1000 sub-steps, the full
minimal working example (about 150,000 points) costs hours per equilibrium
iteration, so bring it up on one element first.

## Results on IKMHIWI03, 2026-10-01 — stages 0–2 PASS

New `ANSYS.exe` built from `Usermat_P21-V21_v222.F` with the typo fix and the
phi-mode fragments, run with `-np 1`. All three runs: 0 errors, completed to
the same 11 sub-steps (`TIME = 1.1`) as the 9/7 run. The 9/7 executable is
kept as `F:\biofilm_upf_wired_ANSYS_20260907.exe`.

| stage | deck | result |
|---|---|---|
| 0 | `ds_oliver_wired_baseline.dat`, `prop(1) = 0` | 0 errors, 2 warnings (same as 9/7), no trace — as expected |
| 1 | `ds_phi_stage1_k0.dat`, 28 constants, `prop(28) = 1`, `k_α = 0` | 0 errors, 46 trace rows, α stays 0 |
| 2 | `ds_phi_stage2_k002.dat`, `k_α = 0.0181818…` | 0 errors; `check_trace`: `once_per_increment`, `eq36`, `carried` all True; `eq36_worst_abs_error = 0.0` |

**The re-entrancy question is settled on real hardware**: 46 calls for 11
sub-steps, so usermat is called several times per sub-step, and α still
rises exactly once per sub-step — never once per Newton iteration.

Found by these runs and fixed afterwards:

- **Sampling.** `MOD(elemId, 997) = 1` hit only element 1 in the 512-element
  deck, so `varies_between_points` could not be judged. Now
  `TRACE_STRIDE = 37` (about 14 elements there).
- **Which φ.** In time alone, `bio1` (0 → 0.0099990) and `bio2`
  (0 → 0.0060601) rise inside `[0, 1]`, while `locbio1/2` sit at
  1.0000–1.000045, slightly above 1. **`bio` looks like the Gauss-point φ**;
  the wider trace in stage 3 confirms or refutes it between points.
- **`k_α`.** φ peaks near 0.01, not 1, so the `φ_max = 1` rule left the final
  α at 9.9994e-05, 1/200 of its 0.02 target. `k_alpha_from_trace` now takes
  the φ the run actually produced from the stage-1 trace; for this deck it
  gives about 1.8 for α ≈ 0.02.

## Stage 3 on IKMHIWI03, 2026-10-01

Rebuilt with `TRACE_STRIDE = 37`, 0 build errors, `-np 1`.

- **`k_alpha_from_trace` over-estimated `k_α` about 100 times.** The deck
  seeds `MY_BIOSTART1 = 1.0`, so φ = 1 in the seed region, but no traced point
  lay there (traced φ_max ≈ 0.015). `k = 3.64` for α 0.02 distorted an element
  at sub-step 1. Ran instead with `k = α_target / (1.0 × 1.1)`. Fixed since:
  the fragment also traces every point with φ ≥ `TRACE_PHI_MIN` (0.5).
- **α_target 0.05 / 0.1 / 0.2 / 0.3 / 0.4: 0 errors, 11/11 sub-steps.**
  `k ≈ 0.55` (α_target 0.5): element 109 "turning inside out" at sub-step 9.
  **Break point about 0.4–0.45**, assuming φ ≈ 1 at the seed-region points —
  not yet traced, so unconfirmed.
- At α_target 0.3, `check_trace`: `once_per_increment`, `eq36`, `carried` all
  True, and `phi_bio1.varies_between_points` True. **`bio` is the
  Gauss-point φ.**
- **`locbio` is not a φ candidate — it is the partner's own `α_K`.** USSFin
  integrates Eq. 36 into it and their AceGen law uses
  `(locbio1 + locbio2)/2` as `Fg`. With `k_α = K_LOCAL1 = 0.01`, our α equals
  `locbio1 − 1` to all printed digits, one sub-step apart. See
  `OLIVER_MODEL_NOTES.md`, Correction 2.
- **`prop(28) = 3`** added (first locally on IKMHIWI03, now in the fragment
  and pre-flighted): `Sbio_GrowthConst = Sdp_sumLocal − 1`, i.e. the
  partner's own growth driving our material (η = 0). 0 errors, 11/11.
  `one_species_reference.partner_alpha_gap(rows)` measures our α against
  theirs. Note `Sdp_sumLocal` averages the two species' `α_K`, so with
  biofilm 2 absent it gives half of `locbio1 − 1` — a question for Oliver.

## Stage 3 confirmed with the seed region traced, 2026-10-01

Rebuilt at `62c5766` (threshold tracing on). The seed region now appears in
the trace: 32 elements, φ = 0.991–1.002.

| deck | result | from the trace |
|---|---|---|
| α_target 0.4 | 0 errors, 11/11 | α_max = 0.4007; all three checks True |
| `k = 0.5` | element 109 inverts at sub-step 10 | converged through sub-step 9 (α ≈ 0.45), fails on the way to α ≈ 0.50; all three checks True up to there |
| `prop(28) = 3` (partner's α) | 0 errors, 11/11 | `partner_alpha_gap = 0.0` |

**One-species coupling: stages 0–3 PASS. Usable range α ≤ 0.45 on this mesh;
the first failure is at α ≈ 0.50, by element inversion.** For scale, Klempt's
own `k_α = 1e-3` keeps α near 1e-4, far inside it.

**The halving is real, measured.** In the seed region (φ ≈ 1,
`K_LOCAL1 = 0.01`, `T = 1.1`) `locbio1 − 1` should reach about 0.011; the
partner's α (`Sdp_sumLocal − 1`) peaks at 0.005 — half, as the average with
an inactive `locbio2 = 1` predicts. This is the evidence for the question to
Oliver.

## Stage 4 — two species (next)

- `prop(28) = 2` (φ = bio1 + bio2); seed biofilm 2, e.g. `MY_BIOSTART2 = 1.0`
  in a region apart from biofilm 1 (Oliver: biofilm 1 cannot grow into
  biofilm-2 nodes unless their boundary value is released).
- `k = α_target / (φ_max × T)` with α_target well inside the window
  (0.2); confirm φ_max in the trace.
- Pass: 0 errors, 11/11, the three checks True, and `phi_bio2`
  `varies_between_points` True.
- Averaging check: run once with `k_α = K_LOCAL1`. (The prediction first
  written here — "about 2 where one species is present, closer to 1 where
  both are" — was wrong; see the stage 4 results below.)

## Stage 4 on IKMHIWI03, 2026-10-01 — PASS. Both cases complete.

`prop(28) = 2`, `MY_BIOSTART2 = 1.0` (the deck's own BIOFILM2: 8 core
elements inside BIOFILM1's 32), `k = 0.2/(2 × 1.1)`. 0 errors, 11/11, all
three checks True, `phi_bio2.varies_between_points` True, φ_max = 2.002
(both species at 1 in the core).

**One- and two-species coupling now both pass on real ANSYS — the thesis
scope (one and two species) is complete.**

- **The averaging, exactly.** own(s − 1) / Oliver(s) = **2.000000 at every
  traced point, whether biofilm 2 is present or not** (n = 46). With
  `K_LOCAL1 = K_LOCAL2 = k`, ours is `k∫(φ₁ + φ₂) = (α_K1 − 1) + (α_K2 − 1)`
  and theirs is `(α_K1 + α_K2)/2 − 1`, half of that sum everywhere. So
  `sAlpha` always gives half the summed growth; "halved when biofilm 2 is
  absent" was a special case, and the prediction above was wrong.
- **The partner's biofilm-2 transport is unstable at the deck's
  `MY_BETA2 = 0.05`** (500 × BETA1) with `dt = 0.1`: core `bio2` oscillates
  sub-step to sub-step (1.0, 0.05, 0.40, 0.10, 0.31, 0.08, …) and goes
  negative (min −0.021), so `phi_bio2.in_0_1` is False and α_max only
  reaches 0.147 of the 0.2 target. With `MY_BETA2 = 0.0001`, `bio2` decays
  smoothly 1.0 → 0.991, no negatives, α_max 0.2003, 0 errors, 11/11, checks
  True. Not our code — an explicit diffusion step past its stability limit,
  the same failure `PDE_VERIFICATION_FINDINGS.md` §6d describes — but every
  two-species result needs a smaller `BETA2` or `dt`, and it is a question
  for Oliver whether 500 × BETA1 is intended.

## Stage 5 — the point model at the Gauss point (`prop(28) = 4`)

The ecology ODE runs inside the element, with `n_sub` inner steps and one or
two species masked in the server. α grows from the ODE's own `φ`, integrated
over the inner steps (`α ← α + k_α · phi_int`), not from `bio1`. This mode
uses the fragments' own call to the hook, so v222.F's ecology branch
(`prop(6)`) and its `dt > DT_ECO_MAX` refusal stay switched off and untouched.

Pre-flighted here (`tests/test_point_model_fragment.py`): the mock usermat,
linked against the repository's hook and C shim, talks to a real
`material_server --active-species 1` at `dt = 0.1` (1000 inner steps). The
trace passes `once_per_increment`, `carried`, `nsub` and `inactive_zero`, and
a replay of every increment with `ecology_substeps` matches **bit for bit**.

1. Re-paste both fragments; they gained the mode-4 block.
2. **Use the repository's hook this time** (it has `n_sub` and `phi_int`).
   v222.F's own ecology call is the old five-argument form and will no longer
   compile against it. Change only that one call to the new form, which with
   `n_sub = 1` returns exactly what the old one did:

   ```fortran
         CALL biofilm_ecology_hook(<g>, <theta>, <dt>, 1, <g_new>,
        &     PM_PHIINT, <ok>)
   ```

   Keep v222.F's own argument names in the `<…>` positions. `PM_PHIINT` is
   declared by `phi_mode_decl.inc` and is only a scratch output there.
3. Build as in step 4 of the run sheet, with three changes. Copy
   `ansys_usermat\coupling\usermat_py_hook.f` from the repository over the
   folder's copy. Compile the C shim with `cl /c /O2
   ansys_usermat\coupling\biofilm_py_eval.c` and put the `.obj` in the
   folder. Delete the stale `ANSYS.exe` / `.lib` / `.exp` / `.map` and
   `biofilm_py_bridge.mod` first.
4. Start the server in the run directory's shell:
   `python ansys_usermat\coupling\material_server.py --active-species 1`.
5. Material constants: the stage table with these changes:

   | prop | n = 1 | n = 2 |
   |---|---|---|
   | 6 | 0 | 0 |
   | 7 | `k_α` (start with 0.5) | same |
   | 8 | a₁₁ = 1.34 | a₁₁ |
   | 9, 10 | 0 | a₁₂, a₂₂ |
   | 23 | b₁ = THETA_DEMO[15] | b₁ |
   | 24 | 0 | b₂ |
   | other 8–27 | 0 | 0 |
   | 28 | **4** | **4** |

   For n = 2, start the server with `--active-species 2`.
6. **One element first**, `-np 1`. At about 0.07 s per Gauss-point call for
   1000 inner steps, the minimal working example would cost hours per
   iteration.
7. Judge it:

   ```
   python -c "import sys; sys.path.insert(0,'ansys_usermat'); sys.path.insert(0,'ansys_usermat/coupling'); import one_species_reference as r; print(r.check_pm_trace(r.read_pm_trace('pm_trace.csv'), K_ALPHA, 1, THETA))"
   ```

   `THETA` is the list `prop(8:27)`. Pass means:
   - `once_per_increment`, `carried`, `nsub`, `inactive_zero` and `replay` are all `True`;
   - `NUMBER OF ERROR MESSAGES = 0`;
   - no `keycut` cut-backs in the `.out`.

   A `keycut` means the returned state failed the sanity check (NaN, `Σφ > 1.5`, or `|γ| > 1e5`). α is then held, and the trace still records the attempt.

## Stage 5 on IKMHIWI03, 2026-10-01 — PASS for n = 1 and n = 2

Build: the repository's `usermat_py_hook.f` and `biofilm_py_eval.c`
(`cl /c /O2 /MD`), v222.F's old ecology call changed to
`(g, theta, dTime, 1, g_new, PM_PHIINT, ok)`, both fragments from `b002a41`.
0 build errors.

The NEM deck cannot be shrunk to one element, so the point model runs on one
element of it. Element 220 gets mat 2 (`prop(28) = 4`) via `MPCOPY,,1,2` and
`MPCHG,2,220`. The other 511 elements stay mat 1 with `prop(28) = 0`. The run
used `-np 1`, `n_sub = 1000`, and took about 200 s per run.

| | n = 1 | n = 2 |
|---|---|---|
| server | `--active-species 1` | `--active-species 2` |
| θ (k_α = 0.5) | a₁₁ = 1.34, b₁ = 0.32 | + a₁₂ = −0.18, a₂₂ = 1.79, b₂ = 1.49 |
| errors / substeps / keycut | 0 / 11 of 11 / none | 0 / 11 of 11 / none |
| `check_pm_trace` | all five `True` | all five `True` |
| `replay_worst_g`, `replay_worst_alpha` | 0.0, 0.0 | 0.0, 0.0 |
| α at t = 1.1 | 0.507 | 0.519 |
| φ at the end | φ₁ = 0.976 (fixed point from substep 2) | φ₁ = 0.018, φ₂ = 0.965, φ₃₋₅ = 0 |

The only extra warning is "elapsed time > CPU time", which is ANSYS waiting on
the server.

**Reproduced off-machine.** `ecology_substeps` from the same seed and θ,
11 × (dt = 0.1, 1000 inner steps), gives α = 0.5075 / 0.5187 and the same φ.
The Gauss point and the standalone point model agree. The one- and
two-species bridge is therefore complete on real ANSYS:
- the phi modes (stages 0–4);
- the point model (stage 5).

## Verifying the first run

1. Record `NUMBER OF ERROR MESSAGES`, whether the load step reaches its end,
   and `SVAR(84)` plus the local `φ` at a few Gauss points over time.
2. Feed that `φ` history to `one_species_reference.alpha_history`: it must
   reproduce `SVAR(84)` to printed precision. A mismatch means the call site
   reads a different `φ` than it prints, or counts growth per iteration.
3. Then, on a matching problem, compare the `α` field against
   `JAXFEM/klempt2024_quantitative.py`. That is the first cross-check between
   the two codes that has a shared variable to compare.

## Then: two species

The same routine with `PHI_LOC = φ₁ + φ₂` — total biofilm drives growth. That
needs question 2 settled, because the sum is exactly what `Sdp_sumBio` was
meant to be.
