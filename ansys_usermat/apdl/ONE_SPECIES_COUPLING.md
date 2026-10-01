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
| 2 | small `k_α` from `k_alpha_for_target` (`α_target` ~ 0.02) | `check_trace`: `once_per_increment`, `eq36` and `carried` all true |
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
   | 7 | `k_α` | stage 1: `0`; stage 2: `k_alpha_for_target(0.02, 1.0, TIME_total)` |
   | 8–27 | 0 | ecology θ, unused |
   | 28 | 1 or 2 | phi mode: one species (`φ = bio1`) or two (`φ = bio1 + bio2`) |

6. **Two species only:** seed biofilm 2 — `sGdp_Bio2start` is zero in the
   minimal working example. Oliver noted biofilm 1 cannot grow into biofilm-2
   nodes unless their boundary value is released there.
7. Run with `-np 1` (the trace is one file).
8. Judge it: `python -c "import sys; sys.path.insert(0,'ansys_usermat'); import one_species_reference as r; print(r.check_trace(r.read_trace('phi_trace.csv'), K_ALPHA))"`
   with the `k_α` used. Pass = `once_per_increment`, `eq36`, `carried` all
   `True`, plus `NUMBER OF ERROR MESSAGES = 0` and the load step complete.
   The `phi_bio1` / `phi_locbio1` entries say which variable varies between
   points and in time — that is the Gauss-point `φ`.

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
