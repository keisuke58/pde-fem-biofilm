C     wrapper_driver.f — drives BIOFILM_GROWTH_VISCO_V01 (the routine we
C     hand to Oliver) so it can be held against the verified core.
C
C     The wrapper takes (E, nu, biofilm fraction) the way their framework
C     passes material constants; the core takes (C10, C01, D1). This driver
C     exposes the wrapper's output so a test can check the adapter does not
C     change the physics — i.e. that the 0-ULP Abaqus equivalence survives
C     the repackaging.
C
C     stdin : F(3x3), Fv_n(3x3),
C             Young, YoungL, Nu, NuL, Biofilm, Growth, Eta, Dt, C01Ratio, Mtype
C     stdout: vCauchy(6) [ANSYS order], Fv_n1(3x3), keycut, elasticWork,
C             mTangCC(6x6) row-major
      PROGRAM WRAPPER_DRV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F(3,3), FVN(3,3), FVN1(3,3)
      DOUBLE PRECISION VCAU(6), TANG(6,6)
      DOUBLE PRECISION YOUNG, YOUNGL, XNU, XNUL, BIOF, GROW
      DOUBLE PRECISION ETA, DT, C01R, XMTYPE, EWORK
      INTEGER KEYCUT, ID, I, J

      READ(*,*) ((F(I,J),J=1,3),I=1,3)
      READ(*,*) ((FVN(I,J),J=1,3),I=1,3)
      READ(*,*) YOUNG, YOUNGL, XNU, XNUL, BIOF, GROW, ETA, DT,
     &          C01R, XMTYPE

      ID = 1
      EWORK = 0.0D0
      KEYCUT = 0

      CALL BIOFILM_GROWTH_VISCO_V01(F, VCAU, TANG,
     &     YOUNG, YOUNGL, XNU, XNUL, BIOF,
     &     GROW, FVN, FVN1,
     &     ETA, DT, C01R, XMTYPE,
     &     EWORK, KEYCUT, ID)

      WRITE(*,'(6E27.17E3)') (VCAU(I),I=1,6)
      WRITE(*,'(9E27.17E3)') ((FVN1(I,J),J=1,3),I=1,3)
      WRITE(*,'(I4,1X,E27.17E3)') KEYCUT, EWORK
      DO I=1,6
        WRITE(*,'(6E27.17E3)') (TANG(I,J),J=1,6)
      END DO
      END
