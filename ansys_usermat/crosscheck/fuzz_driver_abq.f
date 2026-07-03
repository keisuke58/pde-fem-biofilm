C     Batch driver for the Abaqus UMAT constitutive core.
C     Same per-case stdin format as xcheck_driver_abq.f, but loops until EOF
C     and prints all 16 outputs (SV(6), Fv_new(9), detFe) on ONE line per case,
C     so an adversarial fuzzer can stream thousands of cases through one process.
C     Links umat_biofilm_visco.f.
      PROGRAM FABQ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MTYPE
      DIMENSION F(3,3), FGI(3,3), FV(3,3), FVN(3,3), SV(6)
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
     1     DT, SV, FVN, SSE, SPD, DETFE)
      WRITE(*,'(16E26.17)') (SV(K),K=1,6),
     1     ((FVN(I,J),J=1,3),I=1,3), DETFE
      GOTO 10
 99   CONTINUE
      END
