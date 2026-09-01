# A PINN for the biofilm field — design, for Keio

Written 2026-09-01. **Scoped deliberately as a plan, not an implementation**:
per [`ROADMAP_2026.md`](ROADMAP_2026.md), the weeks to November belong to the
ANSYS contribution, and this lands in the JAXFEM material that is being held
back for the continuation at Keio. Nothing here is started.

The point of writing it now is that the design decisions are clearest while the
PDE, the calibration and the solver are fresh — and they are all in this
repository to check against.

## 1. What it would be for

Not "solve the PDE with a network instead". `JAXFEM/` already solves the
Klempt φ–c–α system, and it is fast, verified and differentiable. Replacing a
working solver with a slower, less accurate one is not a contribution.

Three uses that are not that:

| Use | Why a PINN and not the solver |
|---|---|
| **Surrogate for parameter sweeps** | TMCMC needs the forward model thousands of times. A network trained once amortises that; the solver pays full price per sample. This is the strongest case. |
| **Inverse problems on sparse data** | CLSM gives a few timepoints at a few depths. A PINN takes the residual and the data in one loss, rather than wrapping a solver in an outer optimisation. |
| **Fields the mesh makes awkward** | Moving biofilm/void interfaces, where remeshing costs more than a collocation-point resolve. |

Pick **one** before writing code. Trying to serve all three produces something
that does none well, and the surrogate case is the one with a clear customer
(the calibration loop) already in the repository.

## 2. What already exists to build on

This is why the design is cheap now and expensive later.

- **The PDE, stated and solved** — `JAXFEM/hamilton_pde_jaxfem.py`, backward-Euler
  in time, verified against Klempt Eq. 34–36. That is the residual to embed.
- **A reference solution generator.** The same solver gives supervised data,
  which turns "does the PINN work" from an argument into a measurement.
- **A calibrated posterior.** The TMCMC work gives realistic parameter ranges,
  so the surrogate can be trained over the region that is actually sampled
  rather than an arbitrary box.
- **JAX throughout.** `jax.grad` for the residual, one framework, no bridge.

## 3. The design decisions, and what makes each go wrong

**Loss weighting is the whole problem.** Residual, initial, boundary and data
terms have different scales, and fixed weights are the usual reason a PINN
looks trained and is wrong. Decide up front: adaptive weighting, or
non-dimensionalise the system so the terms are comparable by construction. The
second is more work and fails less often.

**Multi-species stiffness.** Five species with interaction terms — the ecology
has fast and slow modes, and a network fits the slow ones and quietly ignores
the fast. Check per-species residuals, never a total.

**Hard-constrain what you can.** Non-negativity of concentrations and the
composition summing to one should be built into the output layer, not asked of
the loss. A PINN that must *learn* `φ ≥ 0` will spend its capacity there.

**Decide the validation metric before training.** Against the solver, on
parameters *not* in the training set. A PINN that reproduces its training
trajectories is not a surrogate.

## 4. What would make it publishable rather than a re-implementation

PINNs for reaction–diffusion are not new. Two things here are less common:

- **A calibrated posterior to sample over**, so the surrogate can be trained and
  assessed over the region the data actually supports, rather than a guessed
  domain. That is a statement about accuracy where it matters.
- **A verified mechanical consumer downstream.** The α field feeds a
  constitutive law that is 0 ULP against an independent implementation. Error
  in the surrogate can be propagated to a stress, which is a far more useful
  accuracy statement than an L2 norm on the field.

Neither is available to someone starting from the equations alone. Both come
from work that is already done here.

## 5. First week at Keio, if it goes ahead

1. Pick the use case from §1 and write down the validation metric (§3).
2. Generate a reference set with `hamilton_pde_jaxfem.py` over posterior
   samples, held out properly.
3. Non-dimensionalise; get the 1-D single-species case to match the solver
   before adding species or dimensions.
4. Only then the five-species system, with per-species residuals reported.

If step 3 does not match the solver, stop — that is the honest exit, and it is
much cheaper to find in week one than in month three.

## 6. Relationship to option (D)

[`OLIVER_MODEL_NOTES.md`](ansys_usermat/OLIVER_MODEL_NOTES.md) records a
separate continuation idea: borrowing the partner group's NEM derivative
operator to strengthen `JAXFEM/`. These are alternatives for the same slot, not
complements — NEM is a better discretisation of the same solve, a PINN is a
different object with different uses. Decide between them at Keio with §1 in
hand; do not start both.
