# Review of the implementation, 20 Nov 2026

Question: is anything overlooked in the implementation? Three independent readings (the main path of
`coupled_living_consumption.py` against the agreed equations; the options and scripts added since 25 Oct; the failing tests of
the user element) and the full test suite. Result: no defect in the configuration that produced the thesis and paper numbers;
several defects in side options, scripts and tests, fixed or listed below.

## 1. Main configuration: nothing found that changes results

`--fast --fortran --front 1 --gc 1 --gab 1 --alpha0 100 --tab 0.333333 --phi-min 0.001 --tvd [--limiter superbee]
--wall-grad mirror --v2 --ansys-bc --inherit-living --transient-c 1e8 [--dtc 0.0025]`.

Checked by running:
- quarter domain against the full cube (16^3): seed means agree to 1e-13;
- a manufactured solution with the MUSCL front term $-r\,\frac{c}{k+c}\,\mathbf n_c\cdot\nabla\phi_{3D}$: L2 error of
  $\phi_{3D}$ 1.86e-3, 5.34e-4, 1.50e-4 on 16^3, 32^3, 64^3 (order about 1.8), nutrient second order, minmod and superbee;
- the Fortran point model replayed over every substep of the separate and the mixed example (1.92 million steps): residual
  of the Python twin below 1e-20 everywhere, also at the antibiotic onset;
- one 32^3 field on a 16^3 and a 32^3 stress mesh: seed von Mises 5.360e-4 and 5.366e-4 Pa.

Checked by reading: Eqs. 34 to 36 term by term, the nutrient and antibiotic solves and their held cells, the time levels
(consumption with $\sum_i\bar\phi_i$ of the previous outer step; $c^*$, $\alpha^*$ interpolated between the two ends of the
step), the residuals and the Jacobian of `ecology_native.f` against Klempt 2026 Eqs. 16 to 18, the seed mask, the stress
solve.

## 2. Fixed on 20 Nov

| what | where | effect |
|---|---|---|
| `--tend T` divided the time link by T (s = 0.1/T) | `coupled_living_consumption.py`, `S_LINK` | the runs tagged `_tend0.05`, `_tend0.25`, `_tend0.5`, `_tend2` ran the point model 20, 4, 2, 0.5 times too fast; they are used only in `fig_contact_zone.py` (drafts), not in the thesis or the paper. Now s does not depend on `--tend`; those files are not comparable with T_END = 1 runs |
| thickness of the window film plotted one outer step early | `fig_tooth_window_thickness.py` | T* at 50 um: 0.99, 0.95, 0.99 instead of 0.96, 0.92, 0.97; slopes within 0.5 um/T*; "one T* is five to six days" unchanged (`TOOTH_WINDOW.md` corrected) |
| stress file of the window without the run's tag | `real_tooth_coupled.py`, `fig_tooth_window_section.py` | mixed, layered and saliva variants overwrote each other; the section figures of 128^3 mixed and of 64^3 may show the stress of another variant (fields differ by 0.1 to 1 %); table numbers come from tagged JSON files and are right. New runs write tagged files; the figure script warns when it falls back to an untagged one |
| reading of the step dependence at K = 1e5 | `results/2026-10-3d/MECH_K_SWEEP.md`, section 2b | mostly the inner step of the base scheme (nutrient held over the Runge-Kutta stages), not the term |
| a checkpoint was resumed by any run with the same tag | `coupled_living_consumption.py` | the checkpoint stores the command line and is ignored when it differs |
| no provenance in the result files | same | `setting.argv` and `setting.commit` are stored |
| `--tooth-seed`, `--tooth-seed-h` missing from the tag | same | in the tag when given |
| `--alpha0 0` became the case's value in the JAX path | same | fixed |
| 12 failing tests of the user element | `tests/conftest.py`, `p3/uel/homeostatic_test.py` | test isolation and $\phi=1$ prescribed at the barrier; no result affected; suite: 528 passed |

Results of the main configuration are unchanged by these edits (16^3 before and after: identical to 1e-14, the run-to-run
noise of the threaded solver).

## 3. Not changed: deliberate, documented, or for the December read-through of the thesis

- At outer step 0.0025 the point model takes 3 substeps of 8.33e-5 per outer step (1200 per T*), not 1000 of 1e-4 as the
  thesis text says; the difference is below 1e-4 in the state.
- The antibiotic is switched on over the first outer step that ends after $T^*=1/3$, so at step 0.025 it starts to act
  0.008 early (seed stress within 1.5 to 2 % on 16^3); negligible at 0.0025.
- With `--transient-c` the nutrient at the end of a step comes from the last inner solve; with d = 1e8 mm^2/T* this is
  immaterial.
- $\phi_{3D}$ crosses $\phi_{\min}=0.001$ far ahead of the visible front (the source $k_\alpha\alpha$ alone reaches it at
  $T^*=1$), so the count of "colonised" cells is not a measure of the front.

## 4. Open hazards (no change made)

1. `ecology_native.f` reports success when Newton's method ends above its tolerance (only a singular Jacobian is flagged).
   It did not happen in the runs checked. Returning the iteration count would make it visible; the same file is linked
   into the ANSYS and Abaqus builds, so the change needs a rebuild there.
2. No test runs the solver end to end, and none covers the mechanical term, growth from the living biomass, the window,
   the obstacles, the saliva option, snapshots, checkpoints or the tag. `biofilm_ecology_hook_lin` with antibiotic (the path
   of the main runs) has no test either; the replay described above would be a cheap one.
3. `fig_reference_style.py` selects the limiter runs by a glob that would also take new variants of the same grid; the
   minmod values are set by hand for that reason.
4. Side paths with known gaps: `--lw` with `--mech-exp` or `--phi-clip`; `--transport-state` and `--budget` without the newer
   terms; `--mech-feedback` or `--homeostatic` with `--obstacle stl` (the tooth is treated as void in that stress solve);
   `--stiff-species` does not check the number of moduli; a 32^3 window run cannot be post-processed (no `nf` in its name).
5. The crown example (`--obstacle stl --box 10`): the signed distance has the wrong sign above the cap, so the start ring
   can also place biofilm inside the crown. The paper uses the cervical window instead.
6. Stale checkpoints of killed runs may lie in `figdata/ckpt/` on each machine; with the check above they are ignored and
   removed when a run with other options meets them.
