# Code verification of the PDE solvers, and what it found

Done 2026-09-30. Code: [`JAXFEM/pde_code_verification.py`](JAXFEM/pde_code_verification.py),
pinned by `tests/test_pde_code_verification.py` (10 tests).

This is **verification**, not validation: it asks whether the discretisation
solves the equations it claims to, and at what order. It needs no experimental
data, no calibration and no ANSYS — which is why it was the first thing worth
doing. Every validation claim in this repository rests on the solver being
right, and until now nothing tested that.

**Headline: both nutrient/species solvers get their zero-flux wall wrong, and
in 2D the operator does not converge there at all.**

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

## 3. 1D — the committed solver is first order, not second

`nutrient_step` writes the Neumann row as `lap[0] = (c₁ − c₀)/dx²`. A ghost node
reflecting about node 0 (`c₋₁ = c₁`, the zero-flux wall) gives
`2(c₁ − c₀)/dx²`, and the finite-volume reading — a half-width control volume at
node 0, flux `D(c₁−c₀)/dx` in across a cell of width `dx/2` — gives the same
factor 2. The committed row is half of both, while its comment says "ghost node
approach".

Sup-norm error against the exact solution, steady state reached to zero residual:

| | as committed | reflected about the node |
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

## 4. 2D — the same treatment on four walls, and there it is inconsistent

`core_hamilton_2d_nutrient.laplacian_2d_neumann` writes the identical rows on all
four walls. Its comment — "ghost `u[-1] = u[0]`" — is exactly what the
arithmetic does, so this is a deliberate choice, not a slip. It *is* zero flux,
but for a wall sitting half a cell outside the node.

Against the exact eigenfunction (`−2(π/L)² = −19.74`):

| | sup error | interior only |
|---|---|---|
| as committed | **9.890 → 9.870, order 0.00** | 3.96e−2 → 6.34e−4, order 2.00 |
| reflected | 4.06e−2 → 6.34e−4, **order 2.00** | same |

9.87 is exactly half of 19.74: the wall rows return half the correct second
difference, and since the exact value there does not shrink with `h`, neither
does the error. **The operator is inconsistent at the boundary — not merely
first order.** The interior is clean.

## 5. The consequence in the units the 2D study reports

Zero-flux diffusion must conserve mass. Pure diffusion of a blob, 41×41, 2000
explicit steps, trapezoidal (finite-volume) measure:

| | total species mass |
|---|---|
| as committed | 0.0314159 → 0.0299024, **−4.82 %** |
| reflected about the node | 0.0314159 → 0.0314159, **0.0000 %** (to 10 digits) |

**This is on the reported path.** `condition_spread_2d.py:173` calls
`diffusion_step_species_2d`, which at `core_hamilton_2d_nutrient.py:434` calls
`laplacian_2d_neumann`. The committed condition-separation result — control
0.0338 % → full 0.0450 %, a 1.33× effect — was computed with a species
diffusion operator that leaks about 5 % of the species mass through walls
declared zero-flux. **The leak is two orders of magnitude larger than the
signal being reported.** That does not mean the 1.33× is wrong; it means it has
not been shown to be a property of the model rather than of the boundary
treatment, and the two arms are affected similarly so some of it may cancel.
Re-running both arms with the wall reflected is the way to find out, and it is
the first thing to do next.

## 6. What is deliberately NOT changed here

Neither solver is modified by this branch. The fix is one coefficient in each
place, it is verified above, and it is the obvious thing to do — but it moves a
number that has already been reported, so it is the user's call and not a
drive-by edit. `pde_code_verification.py` carries the corrected operators as
named variants (`nutrient_step_neumann2`, `laplacian_2d_neumann_mirror`), and a
test asserts that the 1D variant at coefficient 1.0 reproduces `nutrient_step`
**bit for bit** — so the comparison is one routine with one coefficient changed,
not a second opinion about the scheme.

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
