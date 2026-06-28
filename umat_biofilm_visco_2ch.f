C======================================================================
C     UMAT 2-channel viscoelastic biofilm — Prony 2-term + growth
C
C     Kinematics:  F = Fe_∞ · Fg   (equilibrium spring)
C                  F = Fe_k · Fv_k · Fg   k=1,2  (Maxwell branches)
C
C     All branches share the same F (total) and Fg (growth).
C     F_mech = F · inv(Fg)  is computed once.
C     Fe_∞  = F_mech          (equilibrium, purely elastic)
C     Fe_k  = F_mech · inv(Fv_k)  (branch k elastic part)
C
C     Stress (additive, Cauchy):
C       sigma = sigma_inf + sigma_1 + sigma_2
C       sigma_inf = Neo-Hookean(Fe_inf) with C10_inf = C10*(1-a1-a2), includes D1 pressure
C       sigma_k   = deviatoric Neo-Hookean(Fe_k) with C10_k = C10*alpha_k
C
C     Viscous update (backward Euler, same as phase2):
C       Fv_k_new = (I + dt/(2*eta_k*Je_k) * tau_k_dev) · Fv_k_old
C
C     Tangent: FD perturbation of F (exact consistent, same as phase2)
C
C     PROPS:
C       PROPS(1) = C10    instantaneous shear modulus [MPa]  (at t=0)
C       PROPS(2) = C01    Mooney-Rivlin C01 [MPa]  (0 = Neo-Hookean)
C       PROPS(3) = D1     compressibility [1/MPa]
C       PROPS(4) = eta_1  fast viscosity [MPa·s]   (τ_1 = eta_1/(2*C10*α_1))
C       PROPS(5) = eta_2  slow viscosity [MPa·s]   (τ_2 = eta_2/(2*C10*α_2))
C       PROPS(6) = alpha_1  Prony fraction branch 1  (0 < α_1 < 1)
C       PROPS(7) = alpha_2  Prony fraction branch 2  (α_1+α_2 < 1)
C       PROPS(8) = mtype  0=Neo-Hookean, 1=Mooney-Rivlin
C
C     STATEV(1:9)   = Fv_1 (row-major 3×3)
C     STATEV(10:18) = Fv_2 (row-major 3×3)
C     → NSTATV = 18
C
C     Growth: ALPHA_G = TEMP (accumulated volumetric growth; Fg=(1+α)·I)
C
C     Relation to 1-channel UMAT (umat_biofilm_visco_phase2.f):
C       Set alpha_2=0, eta_2 very large → reduces to 1-ch with C10_eff=C10*(1-α_1)
C       Set alpha_1=alpha_2=0 (no Prony fractions) → purely elastic
C
C     AFM calibration guide (2-stage creep):
C       Observe two relaxation times τ_fast (1–10 s) and τ_slow (30–300 s).
C       G_0   = C10 (instantaneous)
C       G_∞   = C10*(1-α_1-α_2) (long-time)
C       G_1   = C10*α_1 (branch 1 amplitude)
C       G_2   = C10*α_2 (branch 2 amplitude)
C       η_1   = G_1 * τ_1  →  eta_1 = C10*α_1 * τ_fast
C       η_2   = G_2 * τ_2  →  eta_2 = C10*α_2 * τ_slow
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

C     --- Local ---
      DOUBLE PRECISION C10, C01, D1, ETA1, ETA2, A1, A2, MTYPE
      DOUBLE PRECISION FV1_OLD(3,3), FV2_OLD(3,3)
      DOUBLE PRECISION FV1_NEW(3,3), FV2_NEW(3,3)
      DOUBLE PRECISION FV1_DUM(3,3), FV2_DUM(3,3)
      DOUBLE PRECISION ALPHA_G, DTIME_IN
      DOUBLE PRECISION STRESS_BASE(6), STRESS_PERT(6)
      DOUBLE PRECISION DFGRD_P(3,3)
      DOUBLE PRECISION PERT, SYM_FAC, SSE_C, SPD_C
      INTEGER I, J, K, P, Q, II, JJ
      INTEGER VOIGT_I(6), VOIGT_J(6)
      DATA VOIGT_I /1, 2, 3, 1, 1, 2/
      DATA VOIGT_J /1, 2, 3, 2, 3, 3/

C     --- Read material properties ---
      C10   = PROPS(1)
      C01   = PROPS(2)
      D1    = PROPS(3)
      ETA1  = PROPS(4)
      ETA2  = PROPS(5)
      A1    = PROPS(6)
      A2    = PROPS(7)
      MTYPE = PROPS(8)

      ALPHA_G  = TEMP + DTEMP
      IF (ALPHA_G .LT. 0.0D0) ALPHA_G = 0.0D0
      DTIME_IN = DTIME

C     --- Read Fv1, Fv2 from STATEV ---
      K = 0
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          FV1_OLD(I,J) = STATEV(K)
        END DO
      END DO
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          FV2_OLD(I,J) = STATEV(K)
        END DO
      END DO

C     Initialize if zero (first increment)
      CALL INIT_FV_IF_ZERO(FV1_OLD)
      CALL INIT_FV_IF_ZERO(FV2_OLD)

C     ================================================================
C     STEP 1-6: Base stress computation
C     ================================================================
      CALL BIOFILM_STRESS_CORE_2CH(
     1   DFGRD1, FV1_OLD, FV2_OLD, ALPHA_G, DTIME_IN,
     2   C10, C01, D1, ETA1, ETA2, A1, A2, MTYPE,
     3   NTENS, NDI, NSHR,
     4   STRESS_BASE, FV1_NEW, FV2_NEW, SSE_C, SPD_C)

      DO I = 1, NTENS
        STRESS(I) = STRESS_BASE(I)
      END DO
      SSE = SSE_C
      SPD = SPD_C

C     ================================================================
C     STEP 7: Exact consistent tangent via FD perturbation of F
C     ================================================================
      PERT = 1.0D-7

      DO P = 1, NTENS
        II = VOIGT_I(P)
        JJ = VOIGT_J(P)

        DO I = 1, 3
          DO J = 1, 3
            DFGRD_P(I,J) = DFGRD1(I,J)
          END DO
        END DO

        IF (II .EQ. JJ) THEN
          DO K = 1, 3
            DFGRD_P(II,K) = DFGRD_P(II,K) + PERT*DFGRD1(II,K)
          END DO
        ELSE
          SYM_FAC = 0.5D0
          DO K = 1, 3
            DFGRD_P(II,K) = DFGRD_P(II,K) + SYM_FAC*PERT*DFGRD1(JJ,K)
            DFGRD_P(JJ,K) = DFGRD_P(JJ,K) + SYM_FAC*PERT*DFGRD1(II,K)
          END DO
        END IF

C       Perturb with SAME Fv_old (consistent tangent)
        CALL BIOFILM_STRESS_CORE_2CH(
     1     DFGRD_P, FV1_OLD, FV2_OLD, ALPHA_G, DTIME_IN,
     2     C10, C01, D1, ETA1, ETA2, A1, A2, MTYPE,
     3     NTENS, NDI, NSHR,
     4     STRESS_PERT, FV1_DUM, FV2_DUM, SSE_C, SPD_C)

        DO Q = 1, NTENS
          DDSDDE(Q,P) = (STRESS_PERT(Q) - STRESS_BASE(Q)) / PERT
        END DO

      END DO

C     ================================================================
C     STEP 8: Store updated state variables
C     ================================================================
      K = 0
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          STATEV(K) = FV1_NEW(I,J)
        END DO
      END DO
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          STATEV(K) = FV2_NEW(I,J)
        END DO
      END DO

      RETURN
      END


C======================================================================
C     BIOFILM_STRESS_CORE_2CH
C     2-channel Prony viscoelastic stress routine.
C
C     σ = σ_∞(Fe_∞) + σ_1(Fe1) + σ_2(Fe2)
C
C     Fe_∞ = F_mech          (equilibrium spring; includes pressure)
C     Fe_k = F_mech·Fvk⁻¹  (Maxwell branch k; purely deviatoric)
C
C     C10_inf = C10*(1-α1-α2)  equilibrium stiffness
C     C10_k   = C10*αk         branch k stiffness
C======================================================================
      SUBROUTINE BIOFILM_STRESS_CORE_2CH(
     1   F_IN, FV1_OLD, FV2_OLD, ALPHA_G_IN, DT_IN,
     2   C10, C01, D1, ETA1, ETA2, A1, A2, MTYPE,
     3   NTENS, NDI, NSHR,
     4   STRESS_OUT, FV1_OUT, FV2_OUT, SSE_OUT, SPD_OUT)

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION F_IN(3,3)
      DOUBLE PRECISION FV1_OLD(3,3), FV2_OLD(3,3)
      DOUBLE PRECISION FV1_OUT(3,3), FV2_OUT(3,3)
      DOUBLE PRECISION ALPHA_G_IN, DT_IN
      DOUBLE PRECISION C10, C01, D1, ETA1, ETA2, A1, A2, MTYPE
      DOUBLE PRECISION STRESS_OUT(NTENS)
      DOUBLE PRECISION SSE_OUT, SPD_OUT
      INTEGER NTENS, NDI, NSHR

      DOUBLE PRECISION FG_SCALE, FG_INV(3,3)
      DOUBLE PRECISION FMECH(3,3)
      DOUBLE PRECISION FV_INV(3,3), DETFV
      DOUBLE PRECISION FE(3,3), BE(3,3), DETFE
      DOUBLE PRECISION IDEN(3,3), TEMP3(3,3)
      DOUBLE PRECISION TAU_DEV(3,3)
      DOUBLE PRECISION SIGMA_INF(3,3), SIGMA_K(3,3), SIGMA_TOT(3,3)
      DOUBLE PRECISION TMP1, TMP2, TRCE, PRESS, DT_SAFE
      DOUBLE PRECISION I1B, I2B, I1E, I2E
      DOUBLE PRECISION C10_INF, C10_1, C10_2
      DOUBLE PRECISION JE_INF, SPD_TMP
      INTEGER I, J, K

C     Identity
      DO I = 1, 3
        DO J = 1, 3
          IDEN(I,J) = 0.0D0
        END DO
        IDEN(I,I) = 1.0D0
      END DO

C     Prony stiffnesses
      C10_INF = C10 * (1.0D0 - A1 - A2)
      C10_1   = C10 * A1
      C10_2   = C10 * A2

C     DT safe
      DT_SAFE = MAX(DT_IN, 1.0D-20)

C     === 1. Growth: Fg = (1+α)·I ===
      FG_SCALE = MAX(1.0D0 + ALPHA_G_IN, 1.0D-15)
      DO I = 1, 3
        DO J = 1, 3
          FG_INV(I,J) = IDEN(I,J) / FG_SCALE
        END DO
      END DO

C     === 2. F_mech = F · inv(Fg) = F / FG_SCALE (isotropic Fg) ===
      DO I = 1, 3
        DO J = 1, 3
          FMECH(I,J) = F_IN(I,J) / FG_SCALE
        END DO
      END DO

C     ================================================================
C     EQUILIBRIUM BRANCH: Fe_inf = F_mech
C     Includes compressibility (D1); deviatoric + pressure.
C     ================================================================
      CALL MAT_AAT(FMECH, BE, 3)
      CALL MAT_DET3(FMECH, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15
      JE_INF = DETFE

      TMP1  = DETFE**(-2.0D0/3.0D0)
      I1B   = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      TRCE  = I1B / 3.0D0
      PRESS = (2.0D0/D1)*(DETFE-1.0D0)*DETFE

      DO I = 1, 3
        DO J = 1, 3
          SIGMA_INF(I,J) = (2.0D0*C10_INF*TMP1*(BE(I,J)-TRCE*IDEN(I,J))
     1                      + PRESS*IDEN(I,J)) / DETFE
        END DO
      END DO

C     Mooney-Rivlin correction for infinity branch
      IF (MTYPE .GT. 0.5D0) THEN
        I1E = BE(1,1)+BE(2,2)+BE(3,3)
        I2E = 0.5D0*(I1E**2-(BE(1,1)**2+BE(2,2)**2+BE(3,3)**2
     1        +2.0D0*(BE(1,2)**2+BE(1,3)**2+BE(2,3)**2)))
        I2B = TMP1**2 * I2E
        DO I = 1, 3
          DO J = 1, 3
            TMP2 = 0.0D0
            DO K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            END DO
            TEMP3(I,J) = I1B*TMP1*BE(I,J) - TMP1**2*TMP2
          END DO
        END DO
        TMP2 = (TEMP3(1,1)+TEMP3(2,2)+TEMP3(3,3))/3.0D0
        DO I = 1, 3
          DO J = 1, 3
            SIGMA_INF(I,J) = SIGMA_INF(I,J) +
     1        2.0D0*C01*(TEMP3(I,J)-TMP2*IDEN(I,J))/DETFE
          END DO
        END DO
      END IF

C     SSE: equilibrium contribution
      SSE_OUT = C10_INF*(I1B-3.0D0) + (1.0D0/D1)*(DETFE-1.0D0)**2
      IF (MTYPE .GT. 0.5D0) SSE_OUT = SSE_OUT + C01*(I2B-3.0D0)
      SPD_OUT = 0.0D0

C     ================================================================
C     BRANCH 1: Maxwell element (C10_1, eta_1)
C     Fe1 = F_mech · inv(Fv1)
C     ================================================================
      CALL MAXWELL_BRANCH(FMECH, FV1_OLD, C10_1, ETA1, DT_SAFE,
     1   IDEN, NTENS, NDI, NSHR,
     2   SIGMA_K, FV1_OUT, SPD_TMP)

      DO I = 1, 3
        DO J = 1, 3
          SIGMA_TOT(I,J) = SIGMA_INF(I,J) + SIGMA_K(I,J)
        END DO
      END DO
      SPD_OUT = SPD_OUT + SPD_TMP

C     ================================================================
C     BRANCH 2: Maxwell element (C10_2, eta_2)
C     Fe2 = F_mech · inv(Fv2)
C     ================================================================
      CALL MAXWELL_BRANCH(FMECH, FV2_OLD, C10_2, ETA2, DT_SAFE,
     1   IDEN, NTENS, NDI, NSHR,
     2   SIGMA_K, FV2_OUT, SPD_TMP)

      DO I = 1, 3
        DO J = 1, 3
          SIGMA_TOT(I,J) = SIGMA_TOT(I,J) + SIGMA_K(I,J)
        END DO
      END DO
      SPD_OUT = SPD_OUT + SPD_TMP

C     ================================================================
C     Pack Cauchy stress into Abaqus Voigt
C     ================================================================
      STRESS_OUT(1) = SIGMA_TOT(1,1)
      STRESS_OUT(2) = SIGMA_TOT(2,2)
      STRESS_OUT(3) = SIGMA_TOT(3,3)
      IF (NSHR .GE. 1) STRESS_OUT(4) = SIGMA_TOT(1,2)
      IF (NSHR .GE. 2) STRESS_OUT(5) = SIGMA_TOT(1,3)
      IF (NSHR .GE. 3) STRESS_OUT(6) = SIGMA_TOT(2,3)

      RETURN
      END


C======================================================================
C     MAXWELL_BRANCH
C     Single Maxwell branch: deviatoric Neo-Hookean stress + Fv update.
C     Purely deviatoric (no pressure in branch stress).
C
C     Given:
C       FMECH : mechanical deformation gradient (3×3)
C       FV_OLD: old viscous Fv for this branch (3×3)
C       C10_K : Neo-Hookean modulus for this branch
C       ETA_K : viscosity for this branch
C       DT    : time step
C
C     Returns:
C       SIGMA_K(3,3) : Cauchy stress contribution (deviatoric)
C       FV_OUT(3,3)  : updated Fv
C       SPD_K        : viscous dissipation
C======================================================================
      SUBROUTINE MAXWELL_BRANCH(FMECH, FV_OLD, C10_K, ETA_K, DT,
     1   IDEN, NTENS, NDI, NSHR,
     2   SIGMA_K, FV_OUT, SPD_K)

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION FMECH(3,3), FV_OLD(3,3), FV_OUT(3,3)
      DOUBLE PRECISION IDEN(3,3)
      DOUBLE PRECISION SIGMA_K(3,3)
      DOUBLE PRECISION C10_K, ETA_K, DT, SPD_K
      INTEGER NTENS, NDI, NSHR

      DOUBLE PRECISION FV_INV(3,3), DETFV
      DOUBLE PRECISION FE(3,3), BE(3,3), DETFE
      DOUBLE PRECISION TAU_DEV(3,3), TEMP3(3,3)
      DOUBLE PRECISION TMP1, TRCE, I1B
      INTEGER I, J, K

C     --- If branch is inactive (C10_K ~ 0 or ETA_K huge → no stress) ---
      IF (C10_K .LT. 1.0D-30) THEN
        DO I = 1, 3
          DO J = 1, 3
            SIGMA_K(I,J) = 0.0D0
            FV_OUT(I,J)  = FV_OLD(I,J)
          END DO
        END DO
        SPD_K = 0.0D0
        RETURN
      END IF

C     === Trial elastic: Fe = F_mech · inv(Fv) ===
      CALL MAT_INV3(FV_OLD, FV_INV, DETFV)
      CALL MAT_MULT(FMECH, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1 = DETFE**(-2.0D0/3.0D0)
      I1B  = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      TRCE = I1B / 3.0D0

C     === Kirchhoff deviatoric ===
      DO I = 1, 3
        DO J = 1, 3
          TAU_DEV(I,J) = 2.0D0*C10_K*TMP1*(BE(I,J)-TRCE*IDEN(I,J))
        END DO
      END DO

C     === Viscous Fv update (backward Euler) ===
      IF (ETA_K .GT. 1.0D-20) THEN
        TMP1 = DT / (2.0D0*ETA_K*DETFE)
        DO I = 1, 3
          DO J = 1, 3
            TEMP3(I,J) = IDEN(I,J) + TMP1*TAU_DEV(I,J)
          END DO
        END DO
        CALL MAT_MULT(TEMP3, FV_OLD, FV_OUT, 3)
      ELSE
C       Frozen dashpot (elastic limit for this branch)
        DO I = 1, 3
          DO J = 1, 3
            FV_OUT(I,J) = FV_OLD(I,J)
          END DO
        END DO
      END IF

C     === Recompute Fe with updated Fv ===
      CALL MAT_INV3(FV_OUT, FV_INV, DETFV)
      CALL MAT_MULT(FMECH, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1 = DETFE**(-2.0D0/3.0D0)
      I1B  = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      TRCE = I1B / 3.0D0

      DO I = 1, 3
        DO J = 1, 3
          TAU_DEV(I,J) = 2.0D0*C10_K*TMP1*(BE(I,J)-TRCE*IDEN(I,J))
        END DO
      END DO

C     === Cauchy stress (deviatoric, divide by Je) ===
      DO I = 1, 3
        DO J = 1, 3
          SIGMA_K(I,J) = TAU_DEV(I,J) / DETFE
        END DO
      END DO

C     === Dissipation ===
      SPD_K = 0.0D0
      IF (ETA_K .GT. 1.0D-20) THEN
        DO I = 1, 3
          DO J = 1, 3
            SPD_K = SPD_K + TAU_DEV(I,J)**2
          END DO
        END DO
        SPD_K = SPD_K * DT / (4.0D0*ETA_K*DETFE)
      END IF

      RETURN
      END


C======================================================================
C     INIT_FV_IF_ZERO: initialise Fv=I if all-zero (first increment)
C======================================================================
      SUBROUTINE INIT_FV_IF_ZERO(FV)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION FV(3,3)
      DOUBLE PRECISION TRACE
      INTEGER I, J
      TRACE = FV(1,1) + FV(2,2) + FV(3,3)
      IF (DABS(TRACE) .LT. 1.0D-10) THEN
        DO I = 1, 3
          DO J = 1, 3
            FV(I,J) = 0.0D0
          END DO
          FV(I,I) = 1.0D0
        END DO
      END IF
      RETURN
      END


C======================================================================
C     Utility subroutines (shared with umat_biofilm_visco_phase2.f)
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
      INTEGER I, J
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
