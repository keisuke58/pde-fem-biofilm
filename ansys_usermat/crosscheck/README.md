# Abaqus UMAT ↔ ANSYS USERMAT cross-check

Proves the ANSYS `USERMAT` constitutive core reproduces the **verified Abaqus
UMAT** law to **machine precision** — the confidence needed before handing the
port to Felix and for the thesis verification chapter.

It compiles the two **actual Fortran cores** (not re-implementations) and drives
them over a battery of deformation states, reconstructing the full 3×3 Cauchy
stress tensor from each (accounting for the different Voigt orderings) and
comparing stress, viscous state `Fv`, and `Je`.

| file | role |
|---|---|
| `xcheck_driver_abq.f` | reads `F, Fv, (α,C10,C01,D1,η,mtype,dt)` from stdin, calls `BIOFILM_STRESS_CORE` from `umat_biofilm_visco.f` (Abaqus order 11,22,33,12,13,23) |
| `xcheck_driver_ans.f` | same, calls the core from `usermat_biofilm.f` (ANSYS order 11,22,33,12,23,13) |
| `crosscheck.py` | builds both (gfortran), runs 20 cases, reconstructs tensors, asserts equal |
| `fuzz_driver_abq.f` / `fuzz_driver_ans.f` | batch drivers: loop over stdin until EOF, one output line per case (so thousands of cases stream through a single process) |
| `adversarial.py` | **tries to break the equivalence** — wide random fuzz + pathological battery + a physical frame-indifference probe |

## Run

```bash
python ansys_usermat/crosscheck/crosscheck.py         # 20-case report
python ansys_usermat/crosscheck/adversarial.py        # adversarial hunt
python ansys_usermat/crosscheck/adversarial.py -n 20000   # bigger sweep
python -m pytest ansys_usermat/crosscheck/            # both, as tests
```

Requires `gfortran` + `numpy`. No Abaqus/ANSYS needed.

## Adversarial hunt (`adversarial.py`)

`crosscheck.py` confirms 20 nice states. `adversarial.py` instead *hunts* for a
state where the two cores disagree, so "identical" becomes a claim over thousands
of cases — including the nasty ones — not twenty. Three passes:

1. **Pathological battery** (17 named cases) — each pins one core branch to an
   extreme: near-incompressible (`D1→0`), near-singular `Fe` (the `detFe` clamp),
   shrinkage (`α→−1`) and extreme growth, the Mooney-Rivlin `C01` branch, the
   elastic (`η=0`) and frozen (`η→∞`) limits, tiny/huge `dt`, large simple shear in
   each of the 12/13/23 planes (Voigt 5↔6 map integrity), pure rotation, and a
   sheared prior `Fv`.
2. **Wide random fuzz** (8000 cases by default) over finite-strain `F`, prior `Fv`,
   and every material parameter swept across decades, with a non-finite guard.
3. **Frame indifference** (physical, per core) — in the elastic limit the Cauchy
   stress must co-rotate, `σ(QF)=Qσ(F)Qᵀ`; any sign flip or index transposition in
   the stress assembly (including the Voigt map) would break it.

**Result: PASS — equivalence `|Δ| = 0.00e+00` (0 ULP) across all 8017 cases**, zero
non-finite outputs, and frame-indifference residual `4.9e-17` (machine precision,
*identical* for both cores). "Bit-identical" now rests on a hunt that actively
tried, and failed, to find a counterexample — not on a fixed list.

## Result

**20/20 cases bit-identical (`|Δσ|=|ΔFv|=|ΔJe|=0`)** — identity+growth, uniaxial,
shear, elastic (`η=0`), frozen (`η→∞`), large growth, Neo-Hookean, Mooney-Rivlin,
and 12 random finite-strain states.

> This harness caught a real bug: the initial ANSYS port applied the
> Mooney-Rivlin `C01` term only to the final Cauchy stress, not to the viscous
> flow driver (as the Abaqus reference does), so `Fv` evolved differently for
> `mtype=1, η>0`. Fixed — now bit-identical. (Since the top-level tangent is a
> finite-difference of this same core, core equivalence implies tangent
> equivalence.)
