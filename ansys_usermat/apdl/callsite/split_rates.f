C=======================================================================
C  split_rates.f -- transport/reaction split, hand-over of the reaction
C  rate from the usermat (Gauss points) to USSFin (NEM points).
C  (pde-fem-biofilm, 2026-10-01)
C
C  The usermat's prop(28) = 5 / 6 branch runs the point model at every
C  Gauss point once per substep and stores R_i = (phi_i_after -
C  phi_i_before) / dt here. USSFin reads R_i at the same point and adds
C  it to its explicit update, in place of its local source term.
C  NEM point ID = (elem - 1) * SPLIT_GPPE + kDomIntPt: NEM points are
C  the Gauss points (4096 = 512 elements x 8 in the minimal working
C  example). Shared memory only: -np 1 or -smp, never DMP (each process
C  would see its own copy).
C
C  The usermat is called several times per substep (equilibrium
C  iterations, then the output call) with the same converged inputs, so
C  the point model gives the same answer each time. The cache below
C  makes that one server call per point and substep instead of one per
C  call, and every call still writes the same R.
C=======================================================================
      MODULE biofilm_split
      IMPLICIT NONE
      INTEGER SPLIT_NMAX, SPLIT_GPPE
      PARAMETER (SPLIT_NMAX = 65536, SPLIT_GPPE = 8)
      DOUBLE PRECISION SPLIT_R(2, SPLIT_NMAX)
      INTEGER SPLIT_STAMP(2, SPLIT_NMAX)
      DOUBLE PRECISION SPLIT_CKEY(12, SPLIT_NMAX)
      DOUBLE PRECISION SPLIT_CVAL(13, SPLIT_NMAX)
      INTEGER SPLIT_CSTAMP(2, SPLIT_NMAX)
      SAVE
      DATA SPLIT_R /131072*0.0D0/
      DATA SPLIT_STAMP /131072*-1/
      DATA SPLIT_CSTAMP /131072*-1/
      CONTAINS

      INTEGER FUNCTION SPLIT_ID(ELEM, IP)
      INTEGER ELEM, IP
      SPLIT_ID = (ELEM - 1) * SPLIT_GPPE + IP
      END FUNCTION

      LOGICAL FUNCTION SPLIT_IN_RANGE(ID)
      INTEGER ID
      SPLIT_IN_RANGE = (ID .GE. 1 .AND. ID .LE. SPLIT_NMAX)
      END FUNCTION

C     usermat side: store this substep's reaction rates
      SUBROUTINE SPLIT_PUT(ID, R1, R2, LD, IS)
      INTEGER ID, LD, IS
      DOUBLE PRECISION R1, R2
      IF (.NOT. SPLIT_IN_RANGE(ID)) RETURN
      SPLIT_R(1, ID) = R1
      SPLIT_R(2, ID) = R2
      SPLIT_STAMP(1, ID) = LD
      SPLIT_STAMP(2, ID) = IS
      END SUBROUTINE

C     USSFin side: read them; LD = IS = -1 if never written
      SUBROUTINE SPLIT_GET(ID, R1, R2, LD, IS)
      INTEGER ID, LD, IS
      DOUBLE PRECISION R1, R2
      R1 = 0.0D0
      R2 = 0.0D0
      LD = -1
      IS = -1
      IF (.NOT. SPLIT_IN_RANGE(ID)) RETURN
      R1 = SPLIT_R(1, ID)
      R2 = SPLIT_R(2, ID)
      LD = SPLIT_STAMP(1, ID)
      IS = SPLIT_STAMP(2, ID)
      END SUBROUTINE

C     cache: same substep and same inner state in -> same answer out
      SUBROUTINE SPLIT_CACHE_GET(ID, LD, IS, GOLD, GNEW, PHIINT, HIT)
      INTEGER ID, LD, IS, I
      DOUBLE PRECISION GOLD(12), GNEW(12), PHIINT
      LOGICAL HIT
      HIT = .FALSE.
      IF (.NOT. SPLIT_IN_RANGE(ID)) RETURN
      IF (SPLIT_CSTAMP(1, ID) .NE. LD) RETURN
      IF (SPLIT_CSTAMP(2, ID) .NE. IS) RETURN
      DO I = 1, 12
        IF (SPLIT_CKEY(I, ID) .NE. GOLD(I)) RETURN
      END DO
      DO I = 1, 12
        GNEW(I) = SPLIT_CVAL(I, ID)
      END DO
      PHIINT = SPLIT_CVAL(13, ID)
      HIT = .TRUE.
      END SUBROUTINE

      SUBROUTINE SPLIT_CACHE_PUT(ID, LD, IS, GOLD, GNEW, PHIINT)
      INTEGER ID, LD, IS, I
      DOUBLE PRECISION GOLD(12), GNEW(12), PHIINT
      IF (.NOT. SPLIT_IN_RANGE(ID)) RETURN
      DO I = 1, 12
        SPLIT_CKEY(I, ID) = GOLD(I)
        SPLIT_CVAL(I, ID) = GNEW(I)
      END DO
      SPLIT_CVAL(13, ID) = PHIINT
      SPLIT_CSTAMP(1, ID) = LD
      SPLIT_CSTAMP(2, ID) = IS
      END SUBROUTINE

      END MODULE biofilm_split
