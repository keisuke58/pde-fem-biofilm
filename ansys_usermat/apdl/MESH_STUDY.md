# Mesh study — the light version

Purpose: show that **the ratio Ch5 reports is mesh-stable**. Not a convergence
study of absolute stress; that is a bigger job and a different question (§3).

## Why the ratio and not the stress

§5.6 compares conditions — σ_CH/σ_DH and the like. Both runs use the same mesh,
so discretisation error largely cancels in the ratio. That is a much weaker
requirement than converged absolute stress, and it is what the claim actually
rests on.

The same logic is why §5.4 needs no mesh study at all: wrapper against core is
a build-to-build comparison on one identical mesh, so the error cancels exactly.

## The run

Three mesh levels, `ESIZE` halved each time — the current deck is **2 elements
through each layer's thickness**, which is why level 0 is a floor, not a
baseline:

| level | substrate / growth `ESIZE` | elements through thickness |
|---|---|---|
| 0 (current) | 0.15 / 0.05 | 2 / 2 |
| 1 | 0.075 / 0.025 | 4 / 4 |
| 2 | 0.0375 / 0.0125 | 8 / 8 |

`python ansys_usermat/apdl/make_mesh_levels.py <deck>.dat` writes the variants.
Run each condition at each level. Extract max and mean SEQV per run.

## Acceptance

**The ratio between conditions changes by less than 5 % from level 1 to level 2.**

If it does, report the level-2 ratio and state that it moved less than 5 % under
a halved element size. If it does not, the ratio is not yet mesh-independent and
the honest options are to refine further or to report the comparison
qualitatively — not to pick the level that looks best.

Absolute stresses will still be moving at these levels. Say so rather than
implying otherwise.

## Known hazard — expect level 2 to be the hard one

This is not hypothetical. The deck's own header records `element highly
distorted` failures at **`ESIZE = 0.033` and `ESIZE = 0.05`**, with different
elements failing each time (1216/1461/2, then 6061/1298/6232/6182), and
independent of α magnitude (0.02 and 0.05 both failed) and of substep count
(up to `NSUBST` 200 with `AUTOTS,ON`). The levels above put the growth layer at
0.025 and 0.0125 — inside and below that range.

So plan for level 2 to fail rather than being surprised by it. If it does:

- a level that will not solve is **a result about the mesh**, to be reported,
  not dropped quietly;
- level 0→1 alone still says something — if the ratio barely moves across one
  halving, that is weaker evidence than two, and should be described as such
  rather than presented as convergence;
- the fix that helped last time was a boundary condition, not a finer mesh
  (a symmetry-style `UZ=0` at the axial ends took the error count from 12 to
  3). Distortion here has been a BC problem as much as a resolution one, which
  is worth knowing before spending a session on element sizes.

## Not in scope

A convergence study of absolute stress on the curved-shell deck. Its own header
says it is a smoke test and that such a study "should precede using this for
anything quantitative" — whether that deck becomes the quantitative vehicle at
all is undecided, and settling it is not this study's job.

## 2026-09-03 — absolute-SEQV convergence check (single condition, α=0.01)

Not the ratio-between-conditions study above (still blocked on per-condition
alpha, see `RUN_PREP.md`) — this is the "not in scope" absolute-stress check,
done anyway because it's cheap and honest to know before anything downstream
leans on this deck's numbers. `make_mesh_levels.py t_growth_cylinder_shell.dat
--levels 2`, run on `F:\biofilm_upf\ANSYS.exe` (the plain `usermat_biofilm.f`
build), `ETABLE,SEQV,S,EQV` parsed with a small Python script (element table →
count/min/max/mean).

| level | ESIZE (substrate/growth) | elements | SEQV min | SEQV max | SEQV mean |
|---|---|---|---|---|---|
| 0 (original) | 0.15 / 0.05 | 12,240 | 4.352e-09 | 1.1862e-05 | 9.144e-06 |
| 1 | 0.075 / 0.025 | 58,284 | 8.763e-12 | 1.0324e-05 | 8.356e-06 |
| 2 | 0.0375 / 0.0125 | 775,680 | — | — | — |

Level 0→1: **max −13.0%, mean −8.6%**. Errors 0 both levels (4 benign warnings
at level 1 — NUMMRG node-association note, USERCK stock-material note,
Newton-Raphson reference-force note — none new or concerning). Confirms what
this file already warned: **absolute SEQV is not converged at level 1**, only
level-0→1 stability was ever claimed as weak evidence, and it doesn't hold up
— an 8–13% swing on one halving is not small.

Level 2 did **not** fail from element distortion, unlike the history in
`t_growth_cylinder_shell.dat`'s own header — a different failure mode:
element count jumped to 775,680 (63× level 0, far more than the ~5x from
0→1), and the run hit `*** FATAL *** This model requires more scratch space
than available` at ANSYS's default memory allocation. Retried with explicit
`-m`/`-db`; still in progress as of this writing (a 775k-element nonlinear
solve is genuinely slow, not stuck — converging substeps observed). Whether
it finishes and what mean/max SEQV level 2 gives is open; update this table
once it lands. Note for next time: the jump from ~58k to ~776k elements
between one ESIZE halving suggests the free mesher is refining more than
just the thickness direction once the growth-layer ESIZE gets small enough
to interact with the circumferential/axial sizing — worth checking with
`NLIST`/`ESIZE` diagnostics before assuming "elements through thickness"
alone predicts run cost at finer levels.

## The 512-element model of the 5 Oct slides: seed stress against mesh (3 Oct, Python)

`mesh_study_seed.py` solves the mechanics of Check 3 alone (the 32-element
seed of the partner's deck grows by alpha - 1 = 1.1e-3, stiffness E (phi^2 + f),
E = 10 Pa, nu = 0.49, f = 1e-3, three corner nodes constrained) with a
linear-elastic hex8 solver on three meshes. Stresses in Pa:

| mesh | elements | seed von Mises, mean | seed von Mises, max | seed mean stress | largest von Mises outside |
|---|---|---|---|---|---|
| 8^3 (as the ANSYS model) | 512 | 5.70e-5 | 7.51e-5 | -3.77e-5 | 1.63e-5 |
| 16^3 | 4096 | 3.59e-5 | 8.97e-5 | -2.52e-5 | 2.13e-5 |
| 32^3 | 32768 | 3.26e-5 | 1.09e-4 | -2.26e-5 | 2.97e-5 |

- Averages over the seed converge: the 8^3 mesh is about 75 % (von Mises) and
  67 % (mean stress) above the 32^3 values, the 16^3 mesh within 10 %.
- Maxima do not converge. The seed is a staircase of cubes, and its re-entrant
  corners are stress singularities: the peaks grow with every refinement.
- So the 512-element ANSYS stresses are right in pattern and sign but too
  large in magnitude by up to a factor of about 1.7; quote seed averages from
  a 16^3 mesh, never peaks. Next: the same comparison in ANSYS itself (a 16^3
  deck), on IKMHIWI03.
- The ratio "neighbours about 25x the seeded element" of Check 3 is not this
  table's last column: there the seeded element is one interior element of
  the seed, here the comparison is the seed average against the void.

### The same in ANSYS: a 16^3 deck (prepared 3 Oct, run on IKMHIWI03)

`refine_deck.py` writes the partner's deck on an n^3 mesh, every component
mapped by its place in space (seed 32 -> 256 elements, nutrient layer 64 -> 512,
the three constrained corner nodes unchanged; `--selftest` reproduces the base
at n = 8; all 4096 elements keep the base orientation). Steps:

```powershell
git pull
python ansys_usermat\apdl\refine_deck.py <stage-1 base>.dat F:\biofilm_upf_wired\ds16_base.dat --n 16
#   prints: old element 220 -> new elements [2151, 2152, 2167, 2168, 2407, 2408, 2423, 2424]
python ansys_usermat\apdl\make_wired_deck.py F:\biofilm_upf_wired\ds16_base.dat F:\biofilm_upf_wired\ds16_pv_eq36.dat `
    --set K_LOCAL1=1e-3 --set K_LOCAL2=0 --set MY_BIOSTART2=0.0 `
    --set YOUNG_BIO=1e-5 --set POISSON_BIO=0.49 --set YOUNG_VOID=-1e-3 `
    --props 7=1e-3,28=1 --post both --post-elem 2151
.\ansys_usermat\apdl\run_apdl.ps1 -Deck ds16_pv_eq36.dat -WorkDir F:\biofilm_upf_wired
```
Then compare the seed average von Mises and mean stress of the 8^3 run
(`all_stress_ds_pv_eq36.csv`) with the 16^3 run: Python predicts the 8^3 values
about 1.6-1.75 times the 16^3 ones. Use the same `--set`/`--props` as the 8^3
paper-value runs (see the header of `figs_1005.py`); only the mesh may differ.

Then, in one command (3 Oct):
```powershell
python ansys_usermat\apdl\compare_mesh.py "$R\all_stress_ds_pv_eq36.csv" `
    F:\biofilm_upf_wired\all_stress_ds16_pv_eq36.csv --grid8 8 --track 220
```
It prints the seed averages on both meshes, the 16^3/8^3 ratios next to the
Python ones (0.63 von Mises, 0.67 mean stress), and element 220 against the
eight 16^3 elements inside it. The seed must come out as 32 and 256 elements
with the same mean alpha; otherwise the deck mapping is wrong.

## Volumetric locking (nu = 0.49), 3 Oct

`locking_check.py` solves the same seed problem with B-bar (volumetric part
at the centre point, what SOLID185 does with its default KEYOPT(2) = 0; the
partner's deck sets no KEYOPT) and with full 2x2x2 integration, for several
nu. Seed averages in Pa:

| mesh | nu | B-bar vM | B-bar p | full vM | full p |
|---|---|---|---|---|---|
| 8^3 | 0.30 | 2.48e-5 | -2.56e-5 | 2.30e-5 | -2.74e-5 |
| 8^3 | 0.49 | 5.70e-5 | -3.77e-5 | 5.96e-5 | -1.41e-4 |
| 8^3 | 0.499 | 6.44e-5 | -4.08e-5 | 1.01e-4 | -1.13e-3 |
| 16^3 | 0.30 | 2.19e-5 | -2.19e-5 | 2.12e-5 | -2.26e-5 |
| 16^3 | 0.49 | 3.59e-5 | -2.52e-5 | 3.41e-5 | -5.69e-5 |
| 16^3 | 0.499 | 3.79e-5 | -2.58e-5 | 4.81e-5 | -3.20e-4 |

- With B-bar the mean stress hardly moves from nu = 0.49 to 0.499 (8^3: -3.8e-5
  to -4.1e-5): no locking.
- Full integration locks: at nu = 0.49 its mean stress is 3.7x (8^3) and
  2.3x (16^3) the B-bar value, at 0.499 it is 28x and 12x. The von Mises
  stress is affected much less, as expected (locking is in the pressure).
- So the ANSYS results are free of locking only if SOLID185 really runs with
  B-bar together with the user material. To check on IKMHIWI03: `ETLIST`
  shows KEYOPT(2) = 0 for type 1. The constrained-cube check (homogeneous
  strain) cannot see locking, so it does not answer this.

### Independent check in CalculiX (4 Oct, cloud)

The same problem in CalculiX 2.21 with C3D8, C3D8I and C3D20R up to 32³:
`ansys_usermat/calculix/README.md`. The Python solver with full integration
matches C3D8 to four digits; at 32³ all locking-free elements give a seed von
Mises mean of 3.05–3.27e-5 Pa and a seed mean stress of −2.3e-5 Pa, so the factor
of about 1.7 for the 8³ mesh stands. Without B-bar the 8³ mean stress is six times
too large, which makes the ETLIST check on IKMHIWI03 more than a formality.
