C     --- phi mode, declarations (pde-fem-biofilm, 2026-10-01) ---------
C     Paste with the other local declarations of the partner's usermat.
      DOUBLE PRECISION PHI_USED, ALPHA_N, ALPHA_NEW
C     trace integration point 1 of every element with
C     MOD(elemId, TRACE_STRIDE) = 1. 37 gives about 14 elements in the
C     512-element test deck and about 500 in the 18,750-element MWE.
C     (997 left only element 1 in the test deck, so phi could not be
C     told apart between points -- found in the stage 1 run, 1 Oct.)
      INTEGER TRACE_STRIDE
      PARAMETER (TRACE_STRIDE = 37)
C     ...and every integration point 1 where phi reaches TRACE_PHI_MIN:
C     stride sampling missed the seeded region (phi = 1) in the stage 3
C     run, so the trace never saw the points that grow fastest.
      DOUBLE PRECISION TRACE_PHI_MIN
      PARAMETER (TRACE_PHI_MIN = 0.5D0)
      LOGICAL TRACE_OPEN
      SAVE TRACE_OPEN
      DATA TRACE_OPEN /.FALSE./
C     --- point-model mode (prop(28) = 4) -----------------------------
C     PM_ prefix: no clash with the partner's own local names.
      DOUBLE PRECISION PM_GOLD(12), PM_GNEW(12), PM_TH(20)
      DOUBLE PRECISION PM_PHIINT, PM_PHITOT, PM_SUM
      INTEGER PM_NSUB, PM_I
      LOGICAL PM_OK, PM_TRACE_OPEN
      SAVE PM_TRACE_OPEN
      DATA PM_TRACE_OPEN /.FALSE./
      DOUBLE PRECISION PM_DTMAX
      PARAMETER (PM_DTMAX = 1.0D-4)
C     --- transport/reaction split (prop(28) = 5, 6) -------------------
C     Needs USE biofilm_split next to USE biofilm_py_bridge.
      DOUBLE PRECISION SP_PHI(2), SP_SUM, SP_MIN, SP_R(2)
      INTEGER SP_N, SP_ID, SP_HITI
      LOGICAL SP_HIT, SP_TRACE_OPEN
      SAVE SP_TRACE_OPEN
      DATA SP_TRACE_OPEN /.FALSE./
      DOUBLE PRECISION SP_PHIMAX
      PARAMETER (SP_PHIMAX = 1.0D0 - 1.0D-6)
C     --- composition from the point model (prop(28) = 7) --------------
      DOUBLE PRECISION CM_CHI1, CM_S, CM_F
      LOGICAL CM_TRACE_OPEN
      SAVE CM_TRACE_OPEN
      DATA CM_TRACE_OPEN /.FALSE./
