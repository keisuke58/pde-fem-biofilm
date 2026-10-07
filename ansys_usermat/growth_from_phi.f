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

C=======================================================================
C  BIOFILM_HOMEOSTATIC_FACTOR -- the homeostatic-pressure growth law of
C  the Keio WP2 study (keio_wp2/README.ja.md sec. 4-7; an assumption of
C  this work, not from a paper):
C
C      alpha_dot = k_alpha * phi * max(0, 1 - p / p_h)
C
C  p = -tr(sigma)/3 (compression positive) at the start of the increment.
C  This routine returns the factor G = max(0, 1 - p/p_h) that the caller
C  applies to the increment of Eq. 36. Two forms of p_h:
C      IFORM = 0:  p_h = PREF                    (P_h E k_alpha T*)
C      IFORM = 1:  p_h = PREF * (phi^2 + F)      (local stiffness E(phi^2+f))
C  with phi clipped to [0, 1]. PREF <= 0 (or p_h <= 0) gives G = 1, i.e.
C  Eq. 36 unchanged.
C=======================================================================
      subroutine BIOFILM_HOMEOSTATIC_FACTOR(P, PREF, IFORM, PHI, F, G)
      implicit none
      double precision, intent(in)  :: P, PREF, PHI, F
      integer,          intent(in)  :: IFORM
      double precision, intent(out) :: G
      double precision PH, PC
      G = 1.0d0
      if (PREF .le. 0.0d0) return
      PH = PREF
      if (IFORM .eq. 1) then
        PC = min(max(PHI, 0.0d0), 1.0d0)
        PH = PREF * (PC * PC + F)
      end if
      if (PH .le. 0.0d0) return
      G = max(0.0d0, 1.0d0 - P / PH)
      end subroutine BIOFILM_HOMEOSTATIC_FACTOR
