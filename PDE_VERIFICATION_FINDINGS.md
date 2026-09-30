# Code verification of the PDE solvers, and what it found

Done 2026-09-30. Code: [`JAXFEM/pde_code_verification.py`](JAXFEM/pde_code_verification.py),
pinned by `tests/test_pde_code_verification.py` (10 tests).

This is **verification**, not validation: it asks whether the discretisation
solves the equations it claims to, and at what order. It needs no experimental
data, no calibration and no ANSYS — which is why it was the first thing worth
doing. Every validation claim in this repository rests on the solver being
right, and until now nothing tested that.

**Headline: both nutrient/species solvers had their zero-flux wall wrong — one
order lost in 1D, no convergence at all at the 2D wall. Fixed in all three
places it appeared. It turns out not to move the reported 2D condition spread,
which is measured rather than assumed (§6), and an earlier over-claim to the
contrary is corrected in §5.**

---

## 1. Why the existing test could not see this

`tests/test_pde_uniform_consistency.py` compares the 1D and 2D codes to the 0D
integrator on a **uniform** field. A uniform field makes every Laplacian vanish
identically, so that test cannot see the spatial operator or either boundary
condition. It is a real test of the reaction coupling and the nutrient scaling —
it caught the `c* = 1` bug — but it is blind to everything below.

## 2. The exact solutions used

No implementation is imported to build a reference; all three targets are
derived from the continuum statement, the same way
`ansys_usermat/apdl/closed_form_reference.py` avoids circularity.

For `dc/dt = D c'' − g φ c/(k+c)`, `c'(0) = 0`, `c(L) = 1`, the Monod factor has
no closed-form steady state, but both of its asymptotic limits do and both are
reachable exactly by choosing `k_monod`:

| limit | steady solution | why this one |
|---|---|---|
| first order, `c ≪ k` | `c(x) = cosh(θx/L)/cosh(θ)`, `θ = L√(gφ/(kD))` | the classical Thiele-modulus slab profile; `core_hamilton_1d_nutrient`'s own docstring already quotes a Thiele number for its defaults without ever having checked the profile |
| zeroth order, `c ≫ k` | `c(x) = 1 − (gφ/2D)(L² − x²)` | **this is Klempt et al. (2024) Eq. 35**, fixed by their Eq. 24 as "the simplest possible functional dependency, a linear relation" |

For the 2D operator, the zero-flux eigenfunctions on `[0,L]²` are
`cos(nπx/L)cos(mπy/L)` with eigenvalue `−((nπ/L)² + (mπ/L)²)`.

## 3. 1D — the wall was first order, not second

`nutrient_step` wrote the Neumann row as `lap[0] = (c₁ − c₀)/dx²`. A ghost node
reflecting about node 0 (`c₋₁ = c₁`, the zero-flux wall) gives
`2(c₁ − c₀)/dx²`, and the finite-volume reading — a half-width control volume at
node 0, flux `D(c₁−c₀)/dx` in across a cell of width `dx/2` — gives the same
factor 2. The old row was half of both, while its comment said "ghost node
approach".

Sup-norm error against the exact solution, steady state reached to zero residual:

| | pre-fix | committed (reflected about the node) |
|---|---|---|
| cosh, `θ = 2` | 1.20e−2 → 1.59e−3, **order 0.95, 0.98, 0.99** | 2.14e−4 → 3.44e−6, **order 2.00, 1.99, 1.97** |
| parabola (Klempt Eq. 35) | 2.50e−2 → 3.13e−3, **order 1.00, 1.00, 1.00** | **8.8e−9, flat** |

The sup norm is attained at `x = 0` on every grid in every case, which localises
it to the wall.

**The flat 8.8e−9 is not a failure to converge — it is exactness.** A
second-order central difference reproduces a quadratic with no error at all, and
the parabola is a quadratic. The residual 8.8e−9 is the `c ≫ k_monod`
asymptotic, not the discretisation: it measures `0.878 × k_monod` across four
decades (`k_monod` = 1e−6 … 1e−10). So with the wall reflected correctly, the
solver reproduces the exact steady solution of Klempt 2024's own nutrient
equation to the limit of the asymptotic approximation.

## 4. 2D — the same treatment on four walls, and there it was inconsistent

`core_hamilton_2d_nutrient.laplacian_2d_neumann` wrote the identical rows on all
four walls. Its comment — "ghost `u[-1] = u[0]`" — was exactly what the
arithmetic did, so this was a deliberate choice, not a slip. It *is* zero flux,
but for a wall sitting half a cell outside the node.

Against the exact eigenfunction (`−2(π/L)² = −19.74`):

| | sup error | interior only |
|---|---|---|
| pre-fix | **9.890 → 9.870, order 0.00** | 3.96e−2 → 6.34e−4, order 2.00 |
| committed (reflected) | 4.06e−2 → 6.34e−4, **order 2.00** | same |

9.87 is exactly half of 19.74: the wall rows return half the correct second
difference, and since the exact value there does not shrink with `h`, neither
does the error. **The operator was inconsistent at the boundary — not merely
first order.** The interior was, and is, clean.

## 5. The consequence in the abstract, and why it does not transfer

Zero-flux diffusion must conserve mass. Pure diffusion of a blob, 41×41, 2000
explicit steps at `dt·D/h² = 0.2`, trapezoidal (finite-volume) measure:

| | total species mass |
|---|---|
| pre-fix | 0.0314159 → 0.0299024, **−4.82 %** |
| reflected about the node | 0.0314159 → 0.0314159, **0.0000 %** (to 10 digits) |

**An earlier version of this section carried that 4.82 % over to the reported
2D condition spread and concluded the leak was "two orders of magnitude larger
than the signal". That was wrong, and it was wrong in the exact way §4 of
`CLAIMS_AND_EVIDENCE.md` warns about: a number measured in one configuration
was quoted as though it applied to another.** The corrected measurement is
below.

The leak scales with how much diffusion actually happens. That blob test runs
a total diffusion number of `2000 × 0.2 = 400`. `condition_spread_2d` runs
`dt_macro·D/dx² ≈ 9.8e−5` per macro step for the fastest species, so even at
900 macro steps the total is `≈ 0.09` — species diffusion is very nearly inert
there, and a boundary error on a nearly-inert operator is nearly nothing.

## 6. The fix, and its measured effect on the reported number

Fixed 2026-09-30 in all three places the treatment appeared:

- `core_hamilton_1d_nutrient.nutrient_step` — 1D nutrient, `x = 0`
- `core_hamilton_2d_nutrient.laplacian_2d_neumann` — 2D species, four walls
- `core_hamilton_2d_nutrient._make_nutrient_step_mixed` — 2D nutrient, three
  Neumann walls (found last; its padding sets `ghost = boundary`, the same
  choice in a different idiom)

`JAXFEM/boundary_fix_impact.py` runs the pre-fix and committed operators
**in one process, at one horizon, on one set of constants**, and carries its
own internal control: the study's `control` arm sets `D_eff = 0`, so the wall
cannot reach it and the two control spreads must come out bit-identical. They
do, at every horizon below.

| horizon | condition spread, pre-fix → committed | relative change |
|---|---|---|
| 60 macro steps | 1.1072260468 % → 1.1072260278 % | **−1.7e−8** |
| 300 | 0.2416998005 % → 0.2416998843 % | **+3.5e−7** |
| 900 | 0.0871953 % → 0.0871955 % | **+1.7e−6** |

So the wall is worth about **one part in 10⁶** at 900 steps, against a
condition effect of 1.16× there (and 1.33× at the 2000 steps the committed
result used). It grows with horizon at roughly the 1.5th power of the step
count, not explosively.

**Conclusion: the reported condition separation is not an artifact of the
boundary treatment.** The number stands. The 2000-step horizon itself was not
re-measured — that run exhausts memory in this container — so this is three
measured horizons plus a tame trend, not a measurement at the reported one.

**A side observation worth its own test.** Because species diffusion is nearly
inert at these settings, whatever separates `full` from `control` is most
likely the nutrient coupling rather than the species diffusion — the `full` arm
switches on both at once. A one-at-a-time run would settle it and has not been
done.

## 6a. A silent no-op, and the guard now standing over it

The first version of `boundary_fix_impact.py` reported the two arms as
**bit-identical at every digit** — i.e. no effect at all. That was not a
finding; it was a comparison of one operator against itself.
`condition_spread_2d` does a top-level `import core_hamilton_2d_nutrient as H`
while the driver did `from JAXFEM import core_hamilton_2d_nutrient`, and
because both the repository root and `JAXFEM/` sit on `sys.path` **those are
two distinct module objects**. The patch went to the one nobody was using.

What caught it was not suspicion of the plumbing but a physical check that
refused to add up: the reaction step was measured to *amplify* a perturbation
(input 2.9e−12 → output 1.4e−10, ~50× per macro step), so a difference that
existed at step 1 could not possibly be absent after twenty. The driver now
patches the module object `condition_spread_2d` itself holds and calls
`_assert_patch_bites` before each arm, which re-computes the species-diffusion
step from the operator it just installed and raises if the module disagrees.
A wrong patch is now an error rather than a result.

## 7. Bearing on the Klempt 2024 gap

`SOLEIMANI2021_NOTES.md` §3 quotes the group's own Remark that the advection
term needs "particular numerical remedies" in the advection-dominant regime
Klempt 2024 runs (`r = 100`), and `KLEMPT2024_REPRODUCTION.md` records the
scheme as a suspicion. This is the first *measured* support for a
discretisation explanation: the nutrient solver these reproductions sit on is
one order below its design order, and its 2D sibling is inconsistent at the
wall. It does not explain the factor of ten on its own — the root cause found
there (the best-fitting variants are not the paper's) stands — but "our scheme
is crude compared to their Galerkin FEM" is no longer a guess.

## 8. What this does and does not buy

It raises **verification** rigor, which needs no data — the question the thesis
faces is how to be rigorous without experimental targets, and this is the
standard answer. It says nothing about **validation**: whether these are the
right equations for oral biofilm remains open and still needs the per-species
LIVE/DEAD data requested from the experimental side.
