C     Batch driver for the ANSYS USERMAT constitutive core
C     (BIOFILM_STRESS_CORE with the VI/VJ Voigt map, ANSYS order 11,22,33,12,23,13).
C     Same per-case stdin format as fuzz_driver_abq.f; loops until EOF and prints
C     all 16 outputs (SV(6), Fv_new(9), detFe) on ONE line per case.
C     Links usermat_biofilm.f.
      PROGRAM FANS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MTYPE
      DIMENSION F(3,3), FGI(3,3), FV(3,3), FVN(3,3), SV(6)
      INTEGER VI(6), VJ(6)
      DATA VI /1,2,3,1,2,1/
      DATA VJ /1,2,3,2,3,3/
 10   CONTINUE
      READ(*,*,END=99) ((F(I,J),J=1,3),I=1,3)
      READ(*,*) ((FV(I,J),J=1,3),I=1,3)
      READ(*,*) ALPHA, C10, C01, D1, ETA, MTYPE, DT
      FGSC = MAX(1.0D0+ALPHA, 1.0D-15)
      DO I=1,3
        DO J=1,3
          FGI(I,J)=0.0D0
        END DO
        FGI(I,I)=1.0D0/FGSC
      END DO
      CALL BIOFILM_STRESS_CORE(F, FGI, FV, C10, C01, D1, ETA, MTYPE,
     1     DT, VI, VJ, SV, FVN, SSE, SPD, DETFE)
      WRITE(*,'(16E26.17)') (SV(K),K=1,6),
     1     ((FVN(I,J),J=1,3),I=1,3), DETFE
      GOTO 10
 99   CONTINUE
      END
