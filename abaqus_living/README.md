# The agreed coupling as an Abaqus user element (A-46, 22 Oct 2026)

The coupling agreed with Oliver (`MEETING_2026-10-20.md`; Klempt et al. 2024 Eqs. 34-36 with the front term, the point
model of Klempt et al. 2026 as published at the Gauss points, consumption by the living biomass, the antibiotic field)
in a standard Galerkin finite element, as a third implementation next to the ANSYS element (NEM, `F:\biofilm_upf_living8`)
and the Python finite-volume reference (`ansys_usermat/apdl/coupled_living_consumption.py`). The equations are the same;
only the discretisation differs. In a Galerkin element the phi Laplacian is symmetric and conserves phi and the walls
need no treatment (A-38..A-45: both had to be fixed in the other two codes).

Abaqus is not part of the LUH thesis (CLAUDE.md); this is the Keio continuation and a cross-check.

| file | what |
|---|---|
| `living_uel.f90` | the element (free form): 8-node hexahedron, dofs u (1-3), phi (11), c (12), a (13); equations in its header |
| `make_abaqus_source.py` | `living_uel_abaqus.for` = `ecology_native.f` + the element in fixed form (the p3 converter) |
| `test_tangent.f90` | closed-form tangent against central differences of the residual (22 Oct: 1e-10; `gate`/`igate` arguments: with the homeostatic law, below) |
| `homeostatic_driver.f90` | one element under the homeostatic-pressure law without Abaqus (free or fully constrained cube; below) |
| `make_living_inp.py` | deck: n^3 cube, `sep` or `mixeq`, increment, antibiotic start; the set-up of the Python reference `--v2 --ansys-bc --start-position --front 1 --artdiff` |
| `summary_living.py` | means per time from the Gauss-point file (seed: wp3 seed cells) |
| `eco_u3_a100*.txt` | the point-model files of the ANSYS runs (Fritsch et al. 2025 case I, alpha*_0 = 100) |

Run on IKMHIWI03 (Abaqus 2024; not from a shell with `dev-env.ps1` dot-sourced):

```
python abaqus_living/make_abaqus_source.py
python abaqus_living/make_living_inp.py 8 sep abaqus_living/lv8_sep_ab.inp --dt 0.005
$env:BIOFILM_ECO_CASE = "<repo>\abaqus_living\eco_u3_a100.txt"; $env:BIOFILM_UEL_GP = "F:\abaqus_work\lv8_sep_ab\living_gp.txt"
.\run_abaqus.ps1 -InpFile abaqus_living/lv8_sep_ab.inp -UserFile abaqus_living/living_uel_abaqus.for -LibsDll
python abaqus_living/summary_living.py F:\abaqus_work\lv8_sep_ab\living_gp.txt 8
```

Procedure as p3/uel (A-11): a 1e-8 step prescribes phi, c, a at every node, then the growth step; transient
`*COUPLED TEMPERATURE-DISPLACEMENT` with fixed increments, `UNSYMM`, dummy C3D8T elements with 1e-9 material.
All `*BOUNDARY` of a step must use the same OP (Abaqus refuses OP=NEW mixed with OP=MOD).

## Homeostatic-pressure growth law (27 Oct 2026; a Keio extension, off by default)

The growth law of the P3 element (`p3/uel/unified_uel.f90`, `p3/uel/README.md` "Homeostatic-pressure growth law";
an assumption of this work, not from a paper) as an option of this element, so that the agreed coupling (living
consumption, antibiotic, inherit rule) and the law can run together. Eq. 36 becomes

    alpha_dot = k_alpha phi g,   g = min(1, max(0, 1 - p / p_h)),   p = -tr(sigma)/3 (positive in compression).

Properties (PROPS(19:21) are taken by mms, SUPG and the rule, so the slots differ from p3's 19, 20; the convention is
the same): **PROPS(22) = p_ref** (absent or <= 0: off), **PROPS(23) = form**: 0 p_h = p_ref, 1 p_h = p_ref (phi_n^2 + f),
2 g = 1/(1 + max(p, 0)/p_ref) (the ANSYS cal5gp form); form + 10 (10, 11, 12) is the implicit gate. Deck:
`make_living_inp.py ... --p-h-e 0.05 [--gate-form 0|1|2] [--implicit-gate]` (p_ref = p_h/E_bio times E_bio, so
`--p-h-e 0.05` with E_bio = 10 Pa is p_h = 0.5 Pa; PROPERTIES=23). No new state variable: with small strain the
start-of-increment pressure is recomputed at the Gauss point from U_n, phi_n and alpha_n (VARIABLES stays 128).

- Explicit gate (forms 0-2, as p3): p of the start of the increment, g constant within it; the tangent keeps its
  form with dt k_alpha replaced by dt k_alpha g. p can overshoot p_h by at most one increment of growth.
- Implicit gate (forms 10-12): p of the end of the increment. With small strain p = p_0 + 3 E' Kb xi is linear in
  xi = alpha - 1 (E' = E_bio (phi_c^2 + f), Kb = 1/(3 (1 - 2 nu)), p_0 = -E' Kb tr eps with the B-bar strain), so
  xi = xi_n + dt k_alpha phi g(p) is a monotone scalar equation per Gauss point (`solve_xi`, safeguarded Newton to
  round-off, no coupling between Gauss points, unlike p3's finite-strain mean dilatation). The tangent takes
  d xi/d phi and d xi/d u from the implicit function theorem, with one-sided dg/dp at the kinks. p never exceeds p_h.

Default unchanged: a deterministic dump of R, K, sigma and the state of one distorted element with started Gauss
points (with and without SUPG) is bit-for-bit identical before and after the change (gfortran -O2, 27 Oct); with the
law on and p = 0 (free growth) alpha - 1 differs from the law-off run by 1e-16 (the round-off of p, 1e-17 Pa).

Verification without Abaqus (`homeostatic_driver.f90`: one unit cube, phi = c = 1, no point model, E = 10 Pa,
nu = 0.49, f = 1e-3, k_alpha = 1e-3; `tests/test_living_homeostatic.py`, about 5 s, in the style of
`tests/test_uel_homeostatic.py`):

- free growth (3-2-1 support): p = 0, alpha - 1 = k_alpha t = 1.1e-3 at T = 1.1 for every form (Eq. 36);
- fully constrained (all u = 0, so p = E (1 + f) xi / (1 - 2 nu) = 500.5 xi Pa), p_h = 0.25 Pa (xi_h = 4.995e-4,
  below the free growth): the element reproduces the scalar forward-Euler sequence of the ODE alpha_dot = k_alpha
  g(p(alpha)) to 1e-19 (explicit gate) and the backward-Euler sequence to 2e-19 (implicit gate); against the ODE
  solved by scipy (alpha - 1 = 4.44276e-4 at T = 1.1) the error is 1.24e-5, 6.16e-6, 3.06e-6, 1.53e-6 for 11, 22,
  44, 88 increments (ratios 2.02, 2.01, 2.01; implicit: 1.19e-5 ... 1.52e-6, ratios 1.98, 1.99, 1.99): first order;
  forms 1 and 2 reproduce their Euler sequences as well;
- homeostatic limit (T = 4, 40 increments, free growth would give 4e-3): alpha - 1 = 4.99435e-4 against
  xi_h = 4.99501e-4, p = 0.24997 Pa, last increment of alpha - 1 1.6e-8 (explicit); implicit: 4.99163e-4, p = 0.24983
  from below, largest p/p_h 0.9993; form 2 keeps growing (alpha - 1 = 1.58e-3, p = 0.79 Pa);
- tangent (`test_tangent.f90 gate|gate1|gate2|igate|igate1|igate2`, random distorted element, g at the Gauss points
  inside (0, 1)): explicit gate 5e-11 to 1e-10 against central differences (the law changes the tangent by 1e-6
  relative at k_alpha = 1e-3); implicit gate 8e-10 to 1.2e-9 with a 1e-7 step of the displacement dofs (the residual
  is linear in u between the kinks; with a 1e-9 step the round-off of the local solve gives 1e-7), the implicit terms
  change the tangent by 1 to 2 % relative.

The shape-gradient bug of p3 (17 Oct 2026, `p3/uel/README.md` "Shape gradient") is not in this element: `gauss`
has dN = J^-1 dN/dr; on a rotated (45 and 30 degrees) and on a skewed element the gradient of x is I to 4e-16 and
free growth u = xi x with alpha - 1 = xi leaves a mechanics residual of 2e-14 (checked 27 Oct).

Not yet run in Abaqus (IKMHIWI03 or fifa): `living_uel_abaqus.for` is regenerated and compiles in fixed form with
gfortran; a deck with `--p-h-e` has not been run.
