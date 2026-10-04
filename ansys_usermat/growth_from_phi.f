C=======================================================================
C  growth_from_phi.f -- one-species growth at a Gauss point, Klempt 2024
C=======================================================================
C
C  The whole growth law for the one-species coupling (2026-10-01): the
C  partner's NEM solves the biofilm field phi and the nutrient c; at each
C  Gauss point this routine turns the local phi into the growth variable,
C  and BIOFILM_GROWTH_VISCO_V01 turns that into stress through Fg.
C
C  Klempt et al. (2024), Biomech Model Mechanobiol 23:2091, Eq. 36:
C
C      alpha_dot - k_alpha * phi = 0
C
C  with their Fg = alpha_K * I, alpha_K(0) = 1. This repository writes
C  Fg = (1 + alpha) I with alpha(0) = 0 (Soleimani, Haverich & Wriggers
C  2021, Eq. 15), i.e. alpha = alpha_K - 1. The shift is constant, so the
C  evolution equation is the same equation and nothing else changes.
C
C  Discretisation: explicit in phi,
C
C      alpha_{n+1} = alpha_n + k_alpha * phi * dt
C
C  where phi is the value the partner's pool holds when the material is
C  called -- the field their USSFin solved at the end of the previous
C  sub-step. Nothing is added to Eq. 36: no cap, no clamp, no gate. That is
C  deliberate (decision of 2026-10-01: follow Klempt for the main results).
C
C  Re-entrancy: ANSYS calls usermat several times per increment, once per
C  equilibrium iteration, and hands back the state of the last CONVERGED
C  increment each time. The caller must therefore pass ALPHA_N from
C  ustatev (the converged value), never the value this routine returned on
C  the previous iteration -- otherwise growth would be counted once per
C  Newton iteration. This routine is a pure function of its inputs, so as
C  long as ALPHA_N comes from ustatev, repeated calls give the same answer.
C=======================================================================
      subroutine BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI, KALPHA, DTIME,
     &                                  ALPHA_N1)
      implicit none
      double precision, intent(in)  :: ALPHA_N, PHI, KALPHA, DTIME
      double precision, intent(out) :: ALPHA_N1
      ALPHA_N1 = ALPHA_N + KALPHA * PHI * DTIME
      end subroutine BIOFILM_ALPHA_FROM_PHI
