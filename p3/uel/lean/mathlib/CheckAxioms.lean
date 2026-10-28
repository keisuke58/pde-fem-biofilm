/- Run:  lake env lean CheckAxioms.lean
   Expected: each theorem depends only on propext, Classical.choice and Quot.sound (the standard axioms);
   no `sorryAx`. -/
import HomeostaticGateReal
open HomeostaticGateReal
#print axioms implicit_no_overshoot
#print axioms growth_monotone
#print axioms explicit_can_overshoot
#print axioms gate_zero
