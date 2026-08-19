C     TEMPLATE — driver for cross-checking a THIRD constitutive core
C     (e.g. Felix's UserElement/UserMat implementation, once it arrives —
C     see THESIS_ASSIGNMENT.md SS4.2 and ansys_usermat/USERELEM_NOTES.md).
C
C     Copy this file to xcheck_driver_<name>.f and fill in the TODOs. Model
C     it on xcheck_driver_abq.f / xcheck_driver_ans.f, which are the two
C     working examples for this same stdin/stdout contract.
C
C     stdin  (same for every driver, so crosscheck.py can drive all of them
C            identically):
C       line 1: F(i,j), i=1,3 j=1,3   (row-major, 9 values)
C       line 2: Fv(i,j), i=1,3 j=1,3  (prior viscous deformation, 9 values)
C       line 3: alpha, C10, C01, D1, eta, mtype, dt
C
C     stdout (same for every driver):
C       line 1: 6 stress components, in THIS core's own Voigt shear order
C               (tell crosscheck.py via --right-voigt/--left-voigt so it can
C               reconstruct the tensor correctly regardless of convention)
C       line 2: 9 components of the updated viscous deformation Fvn
C       line 3: 1 value, Je (or whatever this core's determinant-of-elastic-
C               part quantity is called)
C
C     TODO before this compiles:
C       1. Replace CALL BIOFILM_STRESS_CORE(...) below with a call to
C          Felix's actual entry point — name and argument list will differ.
C          If his top-level routine is only reachable via the full
C          UserElem/UserMat interface (no separable "core" subroutine), this
C          template is not enough — a from-scratch driver stubbing out the
C          rest of the UserElem contract (see USERELEM_NOTES.md) is needed
C          instead, and that's a bigger job worth flagging back rather than
C          quietly doing.
C       2. Fix VI/VJ (or whatever his shear-index arrays are called) to
C          match his Voigt convention, and pass --left-voigt/--right-voigt
C          on the crosscheck.py command line to match.
C       3. Confirm the growth law used to turn ALPHA into Fg^-1 matches ours
C          (FGSC = MAX(1+ALPHA, 1e-15) below) — if his uses a different
C          growth kinematics convention, that's a finding, not a driver bug.
C
      PROGRAM XTEMPLATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MTYPE
      DIMENSION F(3,3), FGI(3,3), FV(3,3), FVN(3,3), SV(6)
      INTEGER VI(6), VJ(6)
C     TODO: confirm this matches Felix's shear-component ordering
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
C     TODO: replace with Felix's actual entry point and argument list
      CALL BIOFILM_STRESS_CORE(F, FGI, FV, C10, C01, D1, ETA, MTYPE,
     1     DT, VI, VJ, SV, FVN, SSE, SPD, DETFE)
      WRITE(*,'(6E26.17)') (SV(K),K=1,6)
      WRITE(*,'(9E26.17)') ((FVN(I,J),J=1,3),I=1,3)
      WRITE(*,'(E26.17)') DETFE
      END
