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
      LOGICAL TRACE_OPEN
      SAVE TRACE_OPEN
      DATA TRACE_OPEN /.FALSE./
