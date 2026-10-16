/-
  The implicit homeostatic gate cannot overshoot at a single material point.

  Growth law (p3/uel/unified_uel.f90, implicit gate, PROPS(20) = form + 10):
      alpha_{n+1} = alpha_n + c * g(p(alpha_{n+1})),   c = dt * k_alpha * phi >= 0,
      g(p) = min(1, max(0, 1 - p/p_h)),  so g(p) = 0 whenever p > p_h,
  with p the pressure (positive in compression) at the point.

  Claim: if p(alpha_0) <= p_h, then p(alpha_n) <= p_h for every n, and alpha never decreases.
  The proof needs only that g vanishes above p_h and that the pressure of the point is a function of
  its own growth alone. It does not need p to be monotone or any property of the hyperelastic law.

  In the finite element model the pressure at a Gauss point also depends on the growth of its
  neighbours, so the second assumption fails there: a point whose own growth has stopped (g = 0) is
  still loaded by its neighbours. This is the overshoot above p_h seen in the implant sector and the
  cervical study with either gate (p3/uel/README.md, "Implicit gate").

  The statements are written for any type R with the few order and arithmetic facts used, given as
  hypotheses, so that the file checks with core Lean alone (no Mathlib). The real numbers satisfy them.

  Check:  lean p3/uel/lean/HomeostaticGate.lean      (Lean 4, core only)
-/

namespace HomeostaticGate

variable {R : Type} [LE R] [LT R] [Add R] [Mul R] [Zero R]

/-- One implicit step does not take the pressure above `ph` if it started at or below `ph`. -/
theorem step_no_overshoot
    (p g : R → R) (ph c a a' : R)
    (step : a' = a + c * g (p a'))
    (gate_zero : ∀ x, ph < x → g x = 0)
    (mul_zero : ∀ x : R, x * 0 = 0)
    (add_zero : ∀ x : R, x + 0 = x)
    (lt_of_not_le : ∀ x y : R, ¬ x ≤ y → y < x)
    (h0 : p a ≤ ph) :
    p a' ≤ ph := by
  apply Classical.byContradiction
  intro hneg
  have hzero : g (p a') = 0 := gate_zero _ (lt_of_not_le _ _ hneg)
  have hsame : a' = a := by
    rw [step, hzero, mul_zero, add_zero]
  apply hneg
  rw [hsame]
  exact h0

/-- Every implicit step keeps the pressure at or below `ph`. -/
theorem no_overshoot
    (p g : R → R) (ph c : R) (alpha : Nat → R)
    (step : ∀ n, alpha (n + 1) = alpha n + c * g (p (alpha (n + 1))))
    (gate_zero : ∀ x, ph < x → g x = 0)
    (mul_zero : ∀ x : R, x * 0 = 0)
    (add_zero : ∀ x : R, x + 0 = x)
    (lt_of_not_le : ∀ x y : R, ¬ x ≤ y → y < x)
    (h0 : p (alpha 0) ≤ ph) :
    ∀ n, p (alpha n) ≤ ph := by
  intro n
  induction n with
  | zero => exact h0
  | succ n ih =>
    exact step_no_overshoot p g ph c (alpha n) (alpha (n + 1)) (step n)
      gate_zero mul_zero add_zero lt_of_not_le ih

omit [LT R] in
/-- Growth never decreases: the gate is not negative and the rate factor `c` is not negative. -/
theorem growth_monotone
    (p g : R → R) (c : R) (alpha : Nat → R)
    (step : ∀ n, alpha (n + 1) = alpha n + c * g (p (alpha (n + 1))))
    (increment_nonneg : ∀ x, 0 ≤ c * g x)
    (le_add : ∀ x d : R, 0 ≤ d → x ≤ x + d) :
    ∀ n, alpha n ≤ alpha (n + 1) := by
  intro n
  rw [step n]
  exact le_add _ _ (increment_nonneg _)

/-!
## Two coupled points: the implicit gate can overshoot

Minimal model of a Gauss point loaded by a neighbour. Point 1's pressure depends on its own growth and on its
neighbour's, p1 = a1 + a2; point 2's only on its own, p2 = a2. Both points use the implicit gate,
a_i' = a_i + c * g(p_i'), with the step gate g(p) = 0 above p_h and 1 at or below (a special case of the gate class
above; the real gate min(1, max(0, 1 - p/p_h)) behaves the same way above p_h). Integer values for a check by
evaluation: p_h = 2, c = 1, start a1 = a2 = 1, so p1 = 2 <= p_h and p2 = 1 <= p_h.

The state a1' = 1, a2' = 2 solves both implicit equations (point 1 has stopped, g(p1') = g(3) = 0; point 2 still
grows, g(p2') = g(2) = 1), and p1' = 3 > p_h. So the single-point theorem fails once the pressure of a point depends
on a neighbour: its own growth stops, but the neighbour's growth keeps raising its pressure.
-/

/-- The step gate on integers: 0 above `ph`, 1 at or below. -/
def stepGate (ph p : Int) : Int := if ph < p then 0 else 1

theorem coupled_implicit_can_overshoot :
    let ph : Int := 2
    let c : Int := 1
    let a1 : Int := 1
    let a2 : Int := 1
    let a1' : Int := 1
    let a2' : Int := 2
    let p1 := fun (x y : Int) => x + y      -- point 1: own growth plus the neighbour's
    let p2 := fun (y : Int) => y            -- point 2: own growth only
    -- both implicit equations hold
    a1' = a1 + c * stepGate ph (p1 a1' a2') ∧
    a2' = a2 + c * stepGate ph (p2 a2') ∧
    -- the start is at or below p_h ...
    p1 a1 a2 ≤ ph ∧ p2 a2 ≤ ph ∧
    -- ... and point 1 ends above p_h, although its own growth has stopped
    ph < p1 a1' a2' ∧ a1' = a1 := by
  decide

end HomeostaticGate
