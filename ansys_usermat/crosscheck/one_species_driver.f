C     one_species_driver.f -- the one-species chain at one Gauss point:
C     local phi -> BIOFILM_ALPHA_FROM_PHI -> BIOFILM_GROWTH_VISCO_V01.
C
C     It mimics what the partner-side call site will do on each converged
C     increment, so a test can hold the chain against closed form without
C     ANSYS. Viscosity off and neo-Hookean (eta = 0, C01 = 0, mtype = 0):
C     the Klempt-2024 setting.
C
C     stdin : F(3x3)
C             Young, YoungL, Nu, NuL, Biofilm, Kalpha, Dt, Nsteps
C             phi_1 ... phi_Nsteps   (one per line)
C     stdout: alpha after every step (Nsteps lines), then
C             vCauchy(6) [ANSYS order] at the final alpha, then keycut
      PROGRAM ONE_SPECIES_DRV
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), FVN(3,3), FVN1(3,3)
      DOUBLE PRECISION VCAU(6), TANG(6,6)
      DOUBLE PRECISION YOUNG, YOUNGL, XNU, XNUL, BIOF, XKA, DT
      DOUBLE PRECISION ALPHA, ALPHA1, PHI, EWORK
      INTEGER NSTEP, K, I, J, KEYCUT, ID

      READ(*,*) ((F(I,J),J=1,3),I=1,3)
      READ(*,*) YOUNG, YOUNGL, XNU, XNUL, BIOF, XKA, DT, NSTEP

      ALPHA = 0.0D0
      DO K = 1, NSTEP
        READ(*,*) PHI
        CALL BIOFILM_ALPHA_FROM_PHI(ALPHA, PHI, XKA, DT, ALPHA1)
        ALPHA = ALPHA1
        WRITE(*,'(E27.17E3)') ALPHA
      END DO

      DO I = 1, 3
        DO J = 1, 3
          FVN(I,J) = 0.0D0
        END DO
        FVN(I,I) = 1.0D0
      END DO
      ID = 1
      EWORK = 0.0D0
      KEYCUT = 0
      CALL BIOFILM_GROWTH_VISCO_V01(F, VCAU, TANG,
     &     YOUNG, YOUNGL, XNU, XNUL, BIOF,
     &     ALPHA, FVN, FVN1,
     &     0.0D0, DT, 0.0D0, 0.0D0,
     &     EWORK, KEYCUT, ID)
      WRITE(*,'(6E27.17E3)') (VCAU(I),I=1,6)
      WRITE(*,'(I4)') KEYCUT
      END
