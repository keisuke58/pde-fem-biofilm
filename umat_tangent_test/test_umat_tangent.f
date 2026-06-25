C======================================================================
C     test_umat_tangent.f
C     Standalone verification of the consistent tangent in
C     umat_biofilm_visco.f -- mirrors phase2_patch_test.py.
C
C     For each test case:
C       (a) call UMAT  -> DDSDDE (forward-perturbation tangent)
C       (b) build a central-difference reference tangent by calling
C           BIOFILM_STRESS_CORE directly (independent of the UMAT path)
C       (c) report symmetry error and max relative error vs FD.
C
C     PASS criteria (same as Python oracle):
C       symmetry error < 5e-3   (Jaumann asymmetry tolerance)
C       max rel err vs FD < 1e-2
C======================================================================
      PROGRAM TEST_UMAT_TANGENT
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION F1(3,3), F3(3,3), FVI(3,3), FV3(3,3), PROPS(5)
      DOUBLE PRECISION E_REF, NU, C10, C01, D1V, DT, AG, ETA0, ETAV

      E_REF = 300.0D0
      NU    = 0.45D0
      C10   = E_REF / (4.0D0*(1.0D0+NU))
      C01   = 0.0D0
      D1V   = 2.0D0*(1.0D0 - 2.0D0*NU) / E_REF
      DT    = 0.1D0
      AG    = 0.05D0
      ETA0  = 0.0D0
      ETAV  = 1.0D3

C     --- F for Case 1 (small elastic) ---
      CALL ZERO33(F1)
      F1(1,1)=1.02D0
      F1(1,2)=0.01D0
      F1(2,1)=0.005D0
      F1(2,2)=1.01D0
      F1(3,3)=1.0D0

C     --- identity Fv ---
      CALL EYE33(FVI)

C     --- Case 3 (viscoelastic) ---
      CALL ZERO33(F3)
      F3(1,1)=1.05D0
      F3(1,2)=0.02D0
      F3(2,1)=0.01D0
      F3(2,2)=1.03D0
      F3(3,3)=1.0D0
      CALL ZERO33(FV3)
      FV3(1,1)=1.01D0
      FV3(1,2)=0.005D0
      FV3(2,1)=0.002D0
      FV3(2,2)=0.99D0
      FV3(3,3)=1.0D0

      PROPS(1)=C10
      PROPS(2)=C01
      PROPS(3)=D1V
      PROPS(5)=0.0D0

      WRITE(*,'(A)') '================================================'
      WRITE(*,'(A)') ' Fortran UMAT consistent-tangent verification'
      WRITE(*,'(A)') '================================================'

C     Case 1: elastic
      PROPS(4)=ETA0
      CALL CHECK_CASE('Case 1 elastic (eta=0)   ',
     1                F1, FVI, AG, DT, PROPS)

C     Case 3: viscoelastic
      PROPS(4)=ETAV
      CALL CHECK_CASE('Case 3 viscoelastic eta=1e3',
     1                F3, FV3, AG, DT, PROPS)

      END

C======================================================================
      SUBROUTINE CHECK_CASE(LABEL, FTOT, FVOLD, AG, DT, PROPS)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*(*) LABEL
      DIMENSION FTOT(3,3), FVOLD(3,3), PROPS(5)
      DIMENSION DDS(6,6), DFD(6,6), FGINV(3,3)
      DIMENSION STRESS(6), STATEV(9), DDSDDE(6,6)
      DIMENSION DDSDDT(6), DRPLDE(6), STRAN(6), DSTRAN(6)
      DIMENSION TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      DIMENSION DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      CHARACTER*80 CMNAME
      DOUBLE PRECISION SERR, FERR, FGSC, AG, DT

C     --- Build DDSDDE via UMAT ---
      DO I=1,6
        STRESS(I)=0.0D0
        DDSDDT(I)=0.0D0
        DRPLDE(I)=0.0D0
        STRAN(I)=0.0D0
        DSTRAN(I)=0.0D0
      END DO
      K=0
      DO I=1,3
        DO J=1,3
          K=K+1
          STATEV(K)=FVOLD(I,J)
        END DO
      END DO
      CALL EYE33(DROT)
      CALL EYE33(DFGRD0)
      DO I=1,3
        DO J=1,3
          DFGRD1(I,J)=FTOT(I,J)
        END DO
      END DO
      TIME(1)=0.0D0
      TIME(2)=0.0D0
      SSE=0.0D0
      SPD=0.0D0
      SCD=0.0D0
      RPL=0.0D0
      DRPLDT=0.0D0
      DTIME=DT
      TEMPV=0.0D0
      DTEMPV=AG
      PNEWDT=1.0D0
      CELENT=1.0D0
      CMNAME='BIOFILM'

      CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMPV, DTEMPV, PREDEF, DPRED,
     3 CMNAME, 3, 3, 6, 9, PROPS, 5, COORDS,
     4 DROT, PNEWDT, CELENT, DFGRD1, DFGRD1, 1, 1, 1,
     5 1, 1, 1)

      DO I=1,6
        DO J=1,6
          DDS(I,J)=DDSDDE(I,J)
        END DO
      END DO

C     --- Build central-difference reference via core directly ---
      FGSC = MAX(1.0D0 + AG, 1.0D-15)
      CALL EYE33(FGINV)
      DO I=1,3
        DO J=1,3
          FGINV(I,J)=FGINV(I,J)/FGSC
        END DO
      END DO
      CALL FD_REF(FTOT, FGINV, FVOLD, PROPS, DT, DFD)

C     --- Metrics ---
      CALL MAXABS66(DDS, DMAX)
      SERR=0.0D0
      DO I=1,6
        DO J=1,6
          IF (DABS(DDS(I,J)-DDS(J,I)) .GT. SERR)
     1        SERR = DABS(DDS(I,J)-DDS(J,I))
        END DO
      END DO
      SERR = SERR / (DMAX + 1.0D-20)

      CALL MAXABS66(DFD, FMAX)
      FERR=0.0D0
      DO I=1,6
        DO J=1,6
          IF (DABS(DDS(I,J)-DFD(I,J)) .GT. FERR)
     1        FERR = DABS(DDS(I,J)-DFD(I,J))
        END DO
      END DO
      FERR = FERR / (FMAX + 1.0D-20)

      WRITE(*,'(1X,A)') LABEL
      WRITE(*,'(4X,A,E12.4,A)') 'symmetry error    = ', SERR,
     1   '   (OK if <5e-3)'
      WRITE(*,'(4X,A,E12.4,A)') 'max rel err vs FD = ', FERR,
     1   '   (OK if <1e-2)'
      IF (SERR .LT. 5.0D-3 .AND. FERR .LT. 1.0D-2) THEN
        WRITE(*,'(4X,A)') 'RESULT: PASS'
      ELSE
        WRITE(*,'(4X,A)') 'RESULT: FAIL'
      END IF
      WRITE(*,*)
      RETURN
      END

C======================================================================
C     Central-difference reference tangent (pert=1e-5), mirrors
C     compute_ddsdde_fd() in phase2_patch_test.py.
C======================================================================
      SUBROUTINE FD_REF(FTOT, FGINV, FVOLD, PROPS, DT, DFD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FTOT(3,3), FGINV(3,3), FVOLD(3,3), PROPS(5), DFD(6,6)
      DIMENSION FP(3,3), FM(3,3), SVP(6), SVM(6)
      DIMENSION FVD(3,3)
      INTEGER VI(6), VJ(6)
      DOUBLE PRECISION PERT, DT, SDUM, PDUM
      DATA VI /1,2,3,1,1,2/
      DATA VJ /1,2,3,2,3,3/

      PERT=1.0D-5
      DO IP=1,6
        II=VI(IP)
        JJ=VJ(IP)
        DO I=1,3
          DO J=1,3
            FP(I,J)=FTOT(I,J)
            FM(I,J)=FTOT(I,J)
          END DO
        END DO
        IF (II .EQ. JJ) THEN
          DO M=1,3
            FP(II,M)=FP(II,M)+0.5D0*PERT*FTOT(JJ,M)
            FM(II,M)=FM(II,M)-0.5D0*PERT*FTOT(JJ,M)
          END DO
        ELSE
          DO M=1,3
            FP(II,M)=FP(II,M)+0.25D0*PERT*FTOT(JJ,M)
            FP(JJ,M)=FP(JJ,M)+0.25D0*PERT*FTOT(II,M)
            FM(II,M)=FM(II,M)-0.25D0*PERT*FTOT(JJ,M)
            FM(JJ,M)=FM(JJ,M)-0.25D0*PERT*FTOT(II,M)
          END DO
        END IF
        CALL BIOFILM_STRESS_CORE(FP, FGINV, FVOLD, PROPS(1), PROPS(2),
     1       PROPS(3), PROPS(4), PROPS(5), DT, SVP, FVD, SDUM, PDUM)
        CALL BIOFILM_STRESS_CORE(FM, FGINV, FVOLD, PROPS(1), PROPS(2),
     1       PROPS(3), PROPS(4), PROPS(5), DT, SVM, FVD, SDUM, PDUM)
        DO IQ=1,6
          DFD(IQ,IP)=(SVP(IQ)-SVM(IQ))/PERT
        END DO
      END DO
      RETURN
      END

C======================================================================
      SUBROUTINE ZERO33(A)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3)
      DO I=1,3
        DO J=1,3
          A(I,J)=0.0D0
        END DO
      END DO
      RETURN
      END

      SUBROUTINE EYE33(A)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3)
      DO I=1,3
        DO J=1,3
          A(I,J)=0.0D0
        END DO
        A(I,I)=1.0D0
      END DO
      RETURN
      END

      SUBROUTINE MAXABS66(A, AMX)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(6,6)
      DOUBLE PRECISION AMX
      AMX=0.0D0
      DO I=1,6
        DO J=1,6
          IF (DABS(A(I,J)) .GT. AMX) AMX=DABS(A(I,J))
        END DO
      END DO
      RETURN
      END
