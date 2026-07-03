C     Driver: exercise the ANSYS USERMAT constitutive core
C     (BIOFILM_STRESS_CORE with VI/VJ, stress in ANSYS order 11,22,33,12,23,13).
C     Same stdin format as the Abaqus driver.  Links usermat_biofilm.f.
      PROGRAM XANS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MTYPE
      DIMENSION F(3,3), FGI(3,3), FV(3,3), FVN(3,3), SV(6)
      INTEGER VI(6), VJ(6)
      DATA VI /1,2,3,1,2,1/
      DATA VJ /1,2,3,2,3,3/
      READ(*,*) ((F(I,J),J=1,3),I=1,3)
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
      WRITE(*,'(6E26.17)') (SV(K),K=1,6)
      WRITE(*,'(9E26.17)') ((FVN(I,J),J=1,3),I=1,3)
      WRITE(*,'(E26.17)') DETFE
      END
