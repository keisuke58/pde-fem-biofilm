/-
  The homeostatic gate over the real numbers, with Mathlib.
  Same statements as ../HomeostaticGate.lean (core Lean, facts of R passed as hypotheses), here for ℝ itself,
  plus a counterexample showing that the explicit gate can overshoot at a single point.

  gate p_h p = min 1 (max 0 (1 - p / p_h))           (p3/uel/unified_uel.f90, form 0)
  implicit:  alpha_{n+1} = alpha_n + c * gate p_h (p alpha_{n+1})
  explicit:  alpha_{n+1} = alpha_n + c * gate p_h (p alpha_n)

  Build (fifa, see CLOUD_TO_FIFA.ja.md):  cd p3/uel/lean/mathlib && lake exe cache get && lake build
  Not yet checked: written in the cloud session, where Mathlib does not fit on the disk.
  Also: CheckAxioms.lean (no sorry, standard axioms only) and ShouldFail.lean (must give an error).
-/
import Mathlib

namespace HomeostaticGateReal

/-- The gate of the homeostatic-pressure growth law. -/
noncomputable def gate (ph p : ℝ) : ℝ := min 1 (max 0 (1 - p / ph))

lemma gate_nonneg (ph p : ℝ) : 0 ≤ gate ph p :=
  le_min zero_le_one (le_max_left _ _)

lemma gate_zero (ph p : ℝ) (hph : 0 < ph) (h : ph < p) : gate ph p = 0 := by
  have h1 : 1 < p / ph := (one_lt_div hph).mpr h
  have h2 : 1 - p / ph ≤ 0 := by linarith
  unfold gate
  rw [max_eq_left h2, min_eq_right zero_le_one]

/-- Implicit gate, one material point: the pressure never goes above `ph`. -/
theorem implicit_no_overshoot (p : ℝ → ℝ) (ph c : ℝ) (hph : 0 < ph) (alpha : ℕ → ℝ)
    (step : ∀ n, alpha (n + 1) = alpha n + c * gate ph (p (alpha (n + 1))))
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

/-- Growth never decreases when the rate factor is not negative. -/
theorem growth_monotone (p : ℝ → ℝ) (ph c : ℝ) (hc : 0 ≤ c) (alpha : ℕ → ℝ)
    (step : ∀ n, alpha (n + 1) = alpha n + c * gate ph (p (alpha (n + 1)))) :
    ∀ n, alpha n ≤ alpha (n + 1) := by
  intro n
  have h := mul_nonneg hc (gate_nonneg ph (p (alpha (n + 1))))
  rw [step n]
  linarith

/-- The explicit gate can overshoot: p(α) = α, p_h = 1, α_0 = 0, c = 2 gives α_1 = 2 > p_h. -/
theorem explicit_can_overshoot :
    ∃ (p : ℝ → ℝ) (ph c a0 : ℝ), 0 < ph ∧ 0 ≤ c ∧ p a0 ≤ ph ∧ ph < p (a0 + c * gate ph (p a0)) := by
  refine ⟨fun x => x, 1, 2, 0, by norm_num, by norm_num, by norm_num, ?_⟩
  norm_num [gate]

end HomeostaticGateReal
