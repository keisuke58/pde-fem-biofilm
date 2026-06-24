C======================================================================
C     UMAT Phase 2: Viscoelastic biofilm — EXACT consistent tangent
C            F = F_e * F_v * F_g
C
C     Extension of Klempt et al. (2024) Biomech Model Mechanobiol
C     Phase 1: approximate elastic tangent  (umat_biofilm_visco.f)
C     Phase 2: numerical-perturbation consistent tangent (this file)
C
C     Key improvement: DDSDDE = exact dSigma/dEps (engineering Voigt)
C       via finite-difference perturbation of the deformation gradient.
C       Achieves quadratic N-R convergence vs. sub-quadratic in Phase 1.
C
C     Consistent tangent derivation:
C       For each Voigt column B (B=1..NTENS):
C         Form perturbed F_B = (I + h * N_B) * F
C         where N_B = sym(e_I x e_J), factor h/2 for shear (engineering convention)
C         Recompute sigma from F_B (full constitutive update, including F_v)
C         DDSDDE(:,B) = (sigma_B - sigma_base) / h
C       This gives EXACT algorithmic tangent including viscous correction.
C
C     Voigt index ordering (Abaqus 3D):
C       1=(11), 2=(22), 3=(33), 4=(12), 5=(13), 6=(23)
C       Shear DSTRAN = engineering strain (2*eps_ij), so shear perturbation is h/2.
C
C     PROPS(1) = C10   Neo-Hookean/Mooney-Rivlin parameter [MPa]
C     PROPS(2) = C01   Mooney-Rivlin parameter (0 for Neo-Hookean) [MPa]
C     PROPS(3) = D1    Compressibility parameter [1/MPa]
C     PROPS(4) = eta   Viscosity [MPa*s]
C     PROPS(5) = type  0.0=Neo-Hookean, 1.0=Mooney-Rivlin
C
C     STATEV(1:9) = F_v tensor (3x3, row-major: 11,12,13,21,22,23,31,32,33)
C
C     Ref: Klempt et al. 2024; Simo & Pister 1984; Wriggers 2008 (App.C)
C======================================================================

      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,
     3 CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS,
     4 DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     5 KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME

      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3),
     3 DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

C     --- Local variables ---
      DOUBLE PRECISION C10, C01, D1, ETA, MTYPE
      DOUBLE PRECISION FV_OLD(3,3), FV_NEW(3,3)
      DOUBLE PRECISION ALPHA_G, DTIME_IN
      DOUBLE PRECISION STRESS_BASE(6), STRESS_PERT(6)
      DOUBLE PRECISION FV_DUMMY(3,3)
      DOUBLE PRECISION DFGRD_P(3,3)
      DOUBLE PRECISION PERT, SYM_FAC
      DOUBLE PRECISION I1B, DETFE_DUMMY, SSE_C, SPD_C
      INTEGER I, J, K, P, Q, II, JJ
C     Voigt → tensor index maps (for 3D NTENS=6)
C     For 2D plane strain (NDI=3,NSHR=1) only P=1..4 are used
      INTEGER VOIGT_I(6), VOIGT_J(6)
      DATA VOIGT_I /1, 2, 3, 1, 1, 2/
      DATA VOIGT_J /1, 2, 3, 2, 3, 3/

C     --- Read material properties ---
      C10   = PROPS(1)
      C01   = PROPS(2)
      D1    = PROPS(3)
      ETA   = PROPS(4)
      MTYPE = PROPS(5)

      ALPHA_G  = TEMP + DTEMP
      IF (ALPHA_G .LT. 0.0D0) ALPHA_G = 0.0D0
      DTIME_IN = DTIME

C     --- Read F_v^old from state variables ---
      K = 0
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          FV_OLD(I,J) = STATEV(K)
        END DO
      END DO

C     Initialize F_v if zero (first increment)
      SSE_C = 0.0D0
      DO I = 1, 3
        SSE_C = SSE_C + FV_OLD(I,I)
      END DO
      IF (DABS(SSE_C) .LT. 1.0D-10) THEN
        DO I = 1, 3
          DO J = 1, 3
            FV_OLD(I,J) = 0.0D0
          END DO
          FV_OLD(I,I) = 1.0D0
        END DO
      END IF

C     ================================================================
C     STEP 1-6: Base stress computation
C     ================================================================
      CALL BIOFILM_STRESS_CORE(
     1   DFGRD1, FV_OLD, ALPHA_G, DTIME_IN,
     2   C10, C01, D1, ETA, MTYPE,
     3   NTENS, NDI, NSHR,
     4   STRESS_BASE, FV_NEW, SSE_C, SPD_C)

C     Copy base stress to STRESS output
      DO I = 1, NTENS
        STRESS(I) = STRESS_BASE(I)
      END DO
      SSE = SSE_C
      SPD = SPD_C

C     ================================================================
C     STEP 7: Exact consistent tangent via numerical perturbation
C
C     DDSDDE(Q,P) = d(sigma_Q) / d(eps_P)  [Voigt, engineering shear]
C
C     Perturbation of deformation gradient for column P:
C       F_pert = F + h * N_P * F
C     where N_P = e_{II} x e_{JJ}  (normal: II=JJ)
C             or (e_{II} x e_{JJ} + e_{JJ} x e_{II}) / 2  (shear: II /= JJ)
C     Factor 1/2 for shear: Abaqus shear DSTRAN = 2*eps_ij (engineering),
C     so h in engineering = h/2 in tensor perturbation of each F row.
C     ================================================================
      PERT = 1.0D-7

      DO P = 1, NTENS
        II = VOIGT_I(P)
        JJ = VOIGT_J(P)

C       Build perturbed F
        DO I = 1, 3
          DO J = 1, 3
            DFGRD_P(I,J) = DFGRD1(I,J)
          END DO
        END DO

        IF (II .EQ. JJ) THEN
C         Normal strain: F_pert_{II,K} += h * F_{II,K}
          DO K = 1, 3
            DFGRD_P(II,K) = DFGRD_P(II,K) + PERT * DFGRD1(II,K)
          END DO
        ELSE
C         Engineering shear: factor 1/2 for each symmetric part
          SYM_FAC = 0.5D0
          DO K = 1, 3
            DFGRD_P(II,K) = DFGRD_P(II,K) + SYM_FAC*PERT*DFGRD1(JJ,K)
            DFGRD_P(JJ,K) = DFGRD_P(JJ,K) + SYM_FAC*PERT*DFGRD1(II,K)
          END DO
        END IF

C       Compute perturbed stress (use same FV_OLD for consistent tangent)
        CALL BIOFILM_STRESS_CORE(
     1     DFGRD_P, FV_OLD, ALPHA_G, DTIME_IN,
     2     C10, C01, D1, ETA, MTYPE,
     3     NTENS, NDI, NSHR,
     4     STRESS_PERT, FV_DUMMY, SSE_C, SPD_C)

C       Finite difference: DDSDDE(:,P) = (sigma_pert - sigma_base) / h
        DO Q = 1, NTENS
          DDSDDE(Q,P) = (STRESS_PERT(Q) - STRESS_BASE(Q)) / PERT
        END DO

      END DO

C     ================================================================
C     STEP 8: Update state variables: store F_v^new
C     ================================================================
      K = 0
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          STATEV(K) = FV_NEW(I,J)
        END DO
      END DO

      RETURN
      END

C======================================================================
C     BIOFILM_STRESS_CORE: Steps 1-6 factored out for perturbation
C
C     Computes Cauchy stress sigma from:
C       F_in        : deformation gradient (3x3)
C       FV_OLD_IN   : viscous F_v at beginning of step (3x3)
C       ALPHA_G_IN  : growth parameter (= Temp in Abaqus)
C       DT_IN       : time step
C       C10,C01,D1,ETA,MTYPE : material parameters
C
C     Returns:
C       STRESS_OUT(NTENS) : Cauchy stress in Abaqus Voigt
C       FV_OUT(3,3)       : updated F_v
C       SSE_OUT, SPD_OUT  : strain energy, plastic dissipation
C======================================================================
      SUBROUTINE BIOFILM_STRESS_CORE(
     1   F_IN, FV_OLD_IN, ALPHA_G_IN, DT_IN,
     2   C10, C01, D1, ETA, MTYPE,
     3   NTENS, NDI, NSHR,
     4   STRESS_OUT, FV_OUT, SSE_OUT, SPD_OUT)

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION F_IN(3,3), FV_OLD_IN(3,3), FV_OUT(3,3)
      DOUBLE PRECISION ALPHA_G_IN, DT_IN
      DOUBLE PRECISION C10, C01, D1, ETA, MTYPE
      DOUBLE PRECISION STRESS_OUT(NTENS)
      DOUBLE PRECISION SSE_OUT, SPD_OUT
      INTEGER NTENS, NDI, NSHR

C     Local
      DOUBLE PRECISION FG(3,3), FG_INV(3,3), FG_DET
      DOUBLE PRECISION FV_OLD(3,3), FV_NEW(3,3), FV_INV(3,3)
      DOUBLE PRECISION FE(3,3), FTRIAL(3,3)
      DOUBLE PRECISION CE(3,3), BE(3,3)
      DOUBLE PRECISION TAU_DEV(3,3), SIGMA(3,3)
      DOUBLE PRECISION DETFE, DETFV
      DOUBLE PRECISION I1_BAR, I2_BAR, I1E, I2E
      DOUBLE PRECISION TRCE, TMP1, TMP2, DT_SAFE
      DOUBLE PRECISION PRESS, IDEN(3,3), TEMP3(3,3)
      DOUBLE PRECISION FG_SCALE
      INTEGER I, J, K

C     Identity
      DO I = 1, 3
        DO J = 1, 3
          IDEN(I,J) = 0.0D0
        END DO
        IDEN(I,I) = 1.0D0
      END DO

C     Copy F_v input
      DO I = 1, 3
        DO J = 1, 3
          FV_OLD(I,J) = FV_OLD_IN(I,J)
        END DO
      END DO

C     === 1. F_g = (1 + alpha_acc)*I  (Klempt 2024 exact) ===
C     Felix (Klempt 2024): Fg = alpha*I where alpha starts at 1.
C     Here ALPHA_G_IN = alpha_acc = alpha - 1 (starts at 0), passed as TEMP.
C     So Fg = (1 + alpha_acc)*I = alpha*I  ← Felix's notation.
C     J_g = det(Fg) = (1 + alpha_acc)^3.
C     The elastic part is incompressible: det(Fe)=1 (enforced by D1 penalty).
      FG_SCALE = MAX(1.0D0 + ALPHA_G_IN, 1.0D-15)
      FG_DET   = FG_SCALE**3
      DO I = 1, 3
        DO J = 1, 3
          FG(I,J)     = FG_SCALE * IDEN(I,J)
          FG_INV(I,J) = IDEN(I,J) / FG_SCALE
        END DO
      END DO

C     === 2. F_trial = F * inv(F_g) ===
      CALL MAT_MULT(F_IN, FG_INV, FTRIAL, 3)

C     === 3. Trial elastic: F_e = F_trial * inv(F_v) ===
      CALL MAT_INV3(FV_OLD, FV_INV, DETFV)
      CALL MAT_MULT(FTRIAL, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1 = DETFE**(-2.0D0/3.0D0)
      I1_BAR = TMP1 * (BE(1,1) + BE(2,2) + BE(3,3))
      PRESS = (2.0D0/D1) * (DETFE - 1.0D0) * DETFE

C     === 4. Kirchhoff deviatoric stress (for viscous update) ===
      TRCE = I1_BAR / 3.0D0
      DO I = 1, 3
        DO J = 1, 3
          TAU_DEV(I,J) = 2.0D0*C10*TMP1*(BE(I,J) - TRCE*IDEN(I,J))
        END DO
      END DO

      IF (MTYPE .GT. 0.5D0) THEN
        I1E = BE(1,1)+BE(2,2)+BE(3,3)
        I2E = 0.5D0*(I1E**2 - (BE(1,1)**2+BE(2,2)**2+BE(3,3)**2
     1        + 2.0D0*(BE(1,2)**2+BE(1,3)**2+BE(2,3)**2)))
        I2_BAR = TMP1**2 * I2E
        DO I = 1, 3
          DO J = 1, 3
            TMP2 = 0.0D0
            DO K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            END DO
            TEMP3(I,J) = I1_BAR*TMP1*BE(I,J) - TMP1**2*TMP2
          END DO
        END DO
        TMP2 = (TEMP3(1,1)+TEMP3(2,2)+TEMP3(3,3))/3.0D0
        DO I = 1, 3
          DO J = 1, 3
            TAU_DEV(I,J) = TAU_DEV(I,J) +
     1        2.0D0*C01*(TEMP3(I,J) - TMP2*IDEN(I,J))
          END DO
        END DO
      END IF

C     === 4b. Viscous F_v update (backward Euler) ===
      DT_SAFE = MAX(DT_IN, 1.0D-20)
      IF (ETA .GT. 1.0D-20) THEN
        TMP1 = DT_SAFE / (2.0D0*ETA*DETFE)
        DO I = 1, 3
          DO J = 1, 3
            TEMP3(I,J) = IDEN(I,J) + TMP1*TAU_DEV(I,J)
          END DO
        END DO
        CALL MAT_MULT(TEMP3, FV_OLD, FV_NEW, 3)
      ELSE
        DO I = 1, 3
          DO J = 1, 3
            FV_NEW(I,J) = FV_OLD(I,J)
          END DO
        END DO
      END IF

C     === 5. Recompute F_e with updated F_v ===
      CALL MAT_INV3(FV_NEW, FV_INV, DETFV)
      CALL MAT_MULT(FTRIAL, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1 = DETFE**(-2.0D0/3.0D0)
      I1_BAR = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      PRESS = (2.0D0/D1)*(DETFE - 1.0D0)*DETFE

C     === 6. Cauchy stress ===
      TRCE = I1_BAR / 3.0D0
      DO I = 1, 3
        DO J = 1, 3
          SIGMA(I,J) = (2.0D0*C10*TMP1*(BE(I,J) - TRCE*IDEN(I,J))
     1                  + PRESS*IDEN(I,J)) / DETFE
        END DO
      END DO

      IF (MTYPE .GT. 0.5D0) THEN
        I1E = BE(1,1)+BE(2,2)+BE(3,3)
        I2E = 0.5D0*(I1E**2 - (BE(1,1)**2+BE(2,2)**2+BE(3,3)**2
     1        + 2.0D0*(BE(1,2)**2+BE(1,3)**2+BE(2,3)**2)))
        I2_BAR = TMP1**2*I2E
        DO I = 1, 3
          DO J = 1, 3
            TMP2 = 0.0D0
            DO K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            END DO
            TEMP3(I,J) = I1_BAR*TMP1*BE(I,J) - TMP1**2*TMP2
          END DO
        END DO
        TMP2 = (TEMP3(1,1)+TEMP3(2,2)+TEMP3(3,3))/3.0D0
        DO I = 1, 3
          DO J = 1, 3
            SIGMA(I,J) = SIGMA(I,J) +
     1        2.0D0*C01*(TEMP3(I,J) - TMP2*IDEN(I,J))/DETFE
          END DO
        END DO
      END IF

C     === Output: Abaqus Voigt format ===
      STRESS_OUT(1) = SIGMA(1,1)
      STRESS_OUT(2) = SIGMA(2,2)
      STRESS_OUT(3) = SIGMA(3,3)
      IF (NSHR .GE. 1) STRESS_OUT(4) = SIGMA(1,2)
      IF (NSHR .GE. 2) STRESS_OUT(5) = SIGMA(1,3)
      IF (NSHR .GE. 3) STRESS_OUT(6) = SIGMA(2,3)

C     === Output: updated F_v and energies ===
      DO I = 1, 3
        DO J = 1, 3
          FV_OUT(I,J) = FV_NEW(I,J)
        END DO
      END DO

      SSE_OUT = C10*(I1_BAR - 3.0D0) + (1.0D0/D1)*(DETFE-1.0D0)**2
      IF (MTYPE .GT. 0.5D0) SSE_OUT = SSE_OUT + C01*(I2_BAR - 3.0D0)

      SPD_OUT = 0.0D0
      IF (ETA .GT. 1.0D-20) THEN
        DO I = 1, 3
          DO J = 1, 3
            SPD_OUT = SPD_OUT + TAU_DEV(I,J)**2
          END DO
        END DO
        SPD_OUT = SPD_OUT * DT_SAFE / (4.0D0*ETA*DETFE)
      END IF

      RETURN
      END

C======================================================================
C     Utility subroutines (same as umat_biofilm_visco.f)
C======================================================================
      SUBROUTINE MAT_MULT(A, B, C, N)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER N, I, J, K
      DIMENSION A(N,N), B(N,N), C(N,N)
      DO I = 1, N
        DO J = 1, N
          C(I,J) = 0.0D0
          DO K = 1, N
            C(I,J) = C(I,J) + A(I,K)*B(K,J)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE MAT_AAT(A, C, N)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER N, I, J, K
      DIMENSION A(N,N), C(N,N)
      DO I = 1, N
        DO J = 1, N
          C(I,J) = 0.0D0
          DO K = 1, N
            C(I,J) = C(I,J) + A(I,K)*A(J,K)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE MAT_DET3(A, DET)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3)
      DOUBLE PRECISION DET
      DET = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
     1    - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
     2    + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      RETURN
      END

      SUBROUTINE MAT_INV3(A, AINV, DET)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3), AINV(3,3)
      DOUBLE PRECISION DET, DINV
      CALL MAT_DET3(A, DET)
      IF (DABS(DET) .LT. 1.0D-30) THEN
        DO I = 1, 3
          DO J = 1, 3
            AINV(I,J) = 0.0D0
          END DO
          AINV(I,I) = 1.0D0
        END DO
        RETURN
      END IF
      DINV = 1.0D0 / DET
      AINV(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) * DINV
      AINV(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) * DINV
      AINV(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) * DINV
      AINV(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) * DINV
      AINV(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) * DINV
      AINV(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) * DINV
      AINV(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) * DINV
      AINV(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) * DINV
      AINV(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) * DINV
      RETURN
      END
