/- Negative test. Run:  lake env lean ShouldFail.lean
   Expected: Lean reports an ERROR for `explicit_no_overshoot_false`. If this file checks without error, the
   setup is not checking proofs (stop and report). The statement claims that the EXPLICIT gate never
   overshoots, which is false (see `explicit_can_overshoot`); the proof below copies the implicit one. -/
import HomeostaticGateReal
open HomeostaticGateReal

theorem explicit_no_overshoot_false (p : ℝ → ℝ) (ph c : ℝ) (hph : 0 < ph) (alpha : ℕ → ℝ)
    (step : ∀ n, alpha (n + 1) = alpha n + c * gate ph (p (alpha n)))
    (h0 : p (alpha 0) ≤ ph) :
    ∀ n, p (alpha n) ≤ ph := by
  intro n
  induction n with
  | zero => exact h0
  | succ n ih =>
    by_contra hneg
    push_neg at hneg
    have hz : gate ph (p (alpha (n + 1))) = 0 := gate_zero ph _ hph hneg
    have hsame : alpha (n + 1) = alpha n := by
      rw [step n, hz, mul_zero, add_zero]
    rw [hsame] at hneg
    linarith
