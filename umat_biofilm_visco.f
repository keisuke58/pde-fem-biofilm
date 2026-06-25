C======================================================================
C     UMAT: Viscoelastic biofilm with multiplicative decomposition
C            F = F_e * F_v * F_g
C
C     Extension of Klempt et al. (2024) Biomech Model Mechanobiol
C     https://pmc.ncbi.nlm.nih.gov/articles/PMC11554842/
C     Klempt: F = Fe * Fg  (elastic + growth, single species)
C     Here  : F = Fe * Fv * Fg  (+ viscous dashpot, 5-species DI-driven)
C
C     Hyperelastic base: Neo-Hookean or Mooney-Rivlin
C     Viscous flow: deviatoric Newtonian dashpot (backward Euler)
C     Growth: isotropic via thermal analogy (F_g = (1 + alpha_T*T) I)
C       TEMP field = alpha_Monod(x) from multiscale_coupling_1d.py
C       This corresponds to Klempt's alpha_growth (local expansion)
C
C     DI -> material params pipeline (set by gen_abaqus_eigenstrain.py):
C       E(DI) = E_max*(1-r)^2 + E_min*r,  r = DI/DI_scale
C       C10   = E(DI)/(4*(1+nu))  [Neo-Hookean: E=6*C10 small strain]
C       D1    = 2*(1-2*nu) / E(DI)   [compressibility]
C
C     PROPS(1) = C10   Neo-Hookean/Mooney-Rivlin parameter [MPa]
C     PROPS(2) = C01   Mooney-Rivlin parameter (0 for Neo-Hookean) [MPa]
C     PROPS(3) = D1    Compressibility parameter [1/MPa]
C     PROPS(4) = eta   Viscosity [MPa*s]
C     PROPS(5) = type  0.0=Neo-Hookean, 1.0=Mooney-Rivlin
C
C     STATEV(1:9) = F_v (3x3 row-major: 11,12,13,21,22,23,31,32,33)
C
C     CONSISTENT TANGENT (Section 7) -- EXACT (2026-06-26):
C     DDSDDE is now the algorithmic consistent tangent, computed by
C     numerical perturbation of the deformation gradient (Sun, Chaikof &
C     Levenston 2008, IJNME 76:1233; Miehe 1996).  Each perturbed stress
C     evaluation re-runs the full backward-Euler viscous update with the
C     OLD viscous state F_v^old frozen, so the viscous-update Jacobian
C     d(F_v^new)/d(eps) is captured automatically -> quadratic N-R.
C     Verified vs finite-difference reference in phase2_patch_test.py
C     (Python replica of the same core): 10/10, rel.err vs FD ~ 2.6e-8.
C     Replaces the previous approximate initial-elastic tangent.
C
C     Ref: Klempt et al. 2024; Junker & Balzani 2021; Soleimani 2019;
C          Sun, Chaikof & Levenston 2008
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
      DOUBLE PRECISION FG_INV(3,3), FV_OLD(3,3), FV_NEW(3,3)
      DOUBLE PRECISION FV_DUM(3,3)
      DOUBLE PRECISION IDEN(3,3), DFG_P(3,3)
      DOUBLE PRECISION SV0(6), SVP(6)
      DOUBLE PRECISION ALPHA_GROWTH, FG_SCALE
      DOUBLE PRECISION SSE0, SPD0, SSE_D, SPD_D, EPS, TMP1
      DOUBLE PRECISION DETFE0, DETFE_D
      INTEGER I, J, K, M, IP, JP, IQ, IPERT
      INTEGER VOIGT_I(6), VOIGT_J(6)

C     Voigt (Abaqus) component -> (i,j): 11,22,33,12,13,23
      DATA VOIGT_I /1,2,3,1,1,2/
      DATA VOIGT_J /1,2,3,2,3,3/

C     --- Read material properties ---
      C10   = PROPS(1)
      C01   = PROPS(2)
      D1    = PROPS(3)
      ETA   = PROPS(4)
      MTYPE = PROPS(5)

C     --- Identity tensor ---
      DO I = 1, 3
        DO J = 1, 3
          IDEN(I,J) = 0.0D0
        END DO
        IDEN(I,I) = 1.0D0
      END DO

C     ================================================================
C     1. Compute inv(F_g) (Klempt 2024 exact): Fg = (1+alpha_acc)*I
C        TEMP = alpha_acc (accumulated growth, starts at 0).
C        alpha_g is a growth field (NOT strain) -> NOT perturbed below.
C     ================================================================
      ALPHA_GROWTH = TEMP + DTEMP
      IF (ALPHA_GROWTH .LT. 0.0D0) ALPHA_GROWTH = 0.0D0
      FG_SCALE = MAX(1.0D0 + ALPHA_GROWTH, 1.0D-15)
      DO I = 1, 3
        DO J = 1, 3
          FG_INV(I,J) = IDEN(I,J) / FG_SCALE
        END DO
      END DO

C     ================================================================
C     2. Read F_v^old from state variables; initialize to I on first inc
C     ================================================================
      K = 0
      DO I = 1, 3
        DO J = 1, 3
          K = K + 1
          FV_OLD(I,J) = STATEV(K)
        END DO
      END DO
      TMP1 = 0.0D0
      DO I = 1, 3
        TMP1 = TMP1 + FV_OLD(I,I)
      END DO
      IF (DABS(TMP1) .LT. 1.0D-10) THEN
        DO I = 1, 3
          DO J = 1, 3
            FV_OLD(I,J) = IDEN(I,J)
          END DO
        END DO
      END IF

C     ================================================================
C     3. Base stress: full algorithmic update at the actual DFGRD1.
C        Returns Cauchy stress (Voigt), updated F_v, SSE, SPD.
C     ================================================================
      CALL BIOFILM_STRESS_CORE(DFGRD1, FG_INV, FV_OLD,
     1     C10, C01, D1, ETA, MTYPE, DTIME, SV0, FV_NEW, SSE0, SPD0,
     2     DETFE0)

C     Robustness: if the elastic volume ratio J_e is non-physical,
C     request a time-increment cutback instead of returning garbage.
      IF (DETFE0 .LT. 0.05D0 .OR. DETFE0 .GT. 20.0D0) THEN
        PNEWDT = 0.5D0
        RETURN
      END IF

      DO I = 1, NTENS
        STRESS(I) = SV0(I)
      END DO
      SSE = SSE0
      SPD = SPD0

C     ================================================================
C     7. Consistent tangent by deformation-gradient perturbation.
C        F-perturbation (Sun et al. 2008), engineering-shear convention:
C          normal P (i=j):  dF = eps * (e_i (x) e_j) F
C          shear  P (i/=j): dF = (eps/2)*((e_i(x)e_j)+(e_j(x)e_i)) F
C        DDSDDE(Q,P) = ( sigma_Q(F+dF^P) - sigma_Q(F) ) / eps
C        Each eval uses the SAME F_v^old (frozen) -> algo tangent.
C     ================================================================
      EPS = 1.0D-7
      DO IPERT = 1, NTENS
        IP = VOIGT_I(IPERT)
        JP = VOIGT_J(IPERT)
        DO I = 1, 3
          DO J = 1, 3
            DFG_P(I,J) = DFGRD1(I,J)
          END DO
        END DO
        IF (IP .EQ. JP) THEN
          DO M = 1, 3
            DFG_P(IP,M) = DFG_P(IP,M) + EPS * DFGRD1(JP,M)
          END DO
        ELSE
          DO M = 1, 3
            DFG_P(IP,M) = DFG_P(IP,M) + 0.5D0 * EPS * DFGRD1(JP,M)
            DFG_P(JP,M) = DFG_P(JP,M) + 0.5D0 * EPS * DFGRD1(IP,M)
          END DO
        END IF

        CALL BIOFILM_STRESS_CORE(DFG_P, FG_INV, FV_OLD,
     1       C10, C01, D1, ETA, MTYPE, DTIME, SVP, FV_DUM, SSE_D, SPD_D,
     2       DETFE_D)

        DO IQ = 1, NTENS
          DDSDDE(IQ,IPERT) = (SVP(IQ) - SV0(IQ)) / EPS
        END DO
      END DO

C     ================================================================
C     8. Update state variables: store F_v^new (from the base eval)
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
C     BIOFILM_STRESS_CORE
C     Given the total deformation gradient DFG, inv(F_g), and the OLD
C     viscous state FV_OLD, perform one backward-Euler viscous update
C     and return:
C        SV(6)   Cauchy stress, Voigt order (11,22,33,12,13,23)
C        FV_NEW   updated viscous deformation gradient
C        SSE_OUT  specific elastic strain energy
C        SPD_OUT  specific viscous dissipation
C     Mirrors the Python replica biofilm_stress_core() in
C     phase2_patch_test.py exactly (same algebra, same conventions).
C======================================================================
      SUBROUTINE BIOFILM_STRESS_CORE(DFG, FG_INV, FV_OLD,
     1     C10, C01, D1, ETA, MTYPE, DTIME, SV, FV_NEW,
     2     SSE_OUT, SPD_OUT, DETFE_OUT)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION DFG(3,3), FG_INV(3,3), FV_OLD(3,3), FV_NEW(3,3), SV(6)
      DOUBLE PRECISION C10, C01, D1, ETA, MTYPE, DTIME, SSE_OUT, SPD_OUT
      DOUBLE PRECISION DETFE_OUT

      DOUBLE PRECISION FTRIAL(3,3), FV_INV(3,3), FE(3,3), BE(3,3)
      DOUBLE PRECISION TAU_DEV(3,3), SIGMA(3,3), IDEN(3,3), TEMP3(3,3)
      DOUBLE PRECISION DETFV, DETFE, TMP1, TMP2, TRCE
      DOUBLE PRECISION I1_BAR, I2_BAR, I1E, I2E, PRESS, DTIME_SAFE
      INTEGER I, J, K

      DO I = 1, 3
        DO J = 1, 3
          IDEN(I,J) = 0.0D0
        END DO
        IDEN(I,I) = 1.0D0
      END DO

C     --- F_trial = F * inv(F_g) ---
      CALL MAT_MULT(DFG, FG_INV, FTRIAL, 3)

C     --- Trial elastic: Fe = Ftrial * inv(Fv_old) ---
      CALL MAT_INV3(FV_OLD, FV_INV, DETFV)
      CALL MAT_MULT(FTRIAL, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1   = DETFE**(-2.0D0/3.0D0)
      I1_BAR = TMP1 * (BE(1,1) + BE(2,2) + BE(3,3))

C     --- Trial deviatoric Kirchhoff stress (viscous flow rule) ---
      TRCE = I1_BAR / 3.0D0
      DO I = 1, 3
        DO J = 1, 3
          TAU_DEV(I,J) = 2.0D0 * C10 * TMP1 *
     1                   (BE(I,J) - TRCE * IDEN(I,J))
        END DO
      END DO
      IF (MTYPE .GT. 0.5D0) THEN
        DO I = 1, 3
          DO J = 1, 3
            TMP2 = 0.0D0
            DO K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            END DO
            TEMP3(I,J) = I1_BAR * TMP1 * BE(I,J) - TMP1**2 * TMP2
          END DO
        END DO
        TMP2 = (TEMP3(1,1) + TEMP3(2,2) + TEMP3(3,3)) / 3.0D0
        DO I = 1, 3
          DO J = 1, 3
            TAU_DEV(I,J) = TAU_DEV(I,J) +
     1        2.0D0 * C01 * (TEMP3(I,J) - TMP2 * IDEN(I,J))
          END DO
        END DO
      END IF

C     --- Viscous update (backward Euler): Fv_new = (I + dDv) Fv_old ---
      DTIME_SAFE = MAX(DTIME, 1.0D-20)
      IF (ETA .GT. 1.0D-20) THEN
        TMP1 = DTIME_SAFE / (2.0D0 * ETA * DETFE)
        DO I = 1, 3
          DO J = 1, 3
            TEMP3(I,J) = IDEN(I,J) + TMP1 * TAU_DEV(I,J)
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

C     --- Recompute Fe with updated Fv, then Cauchy stress ---
      CALL MAT_INV3(FV_NEW, FV_INV, DETFV)
      CALL MAT_MULT(FTRIAL, FV_INV, FE, 3)
      CALL MAT_AAT(FE, BE, 3)
      CALL MAT_DET3(FE, DETFE)
      IF (DETFE .LT. 1.0D-15) DETFE = 1.0D-15

      TMP1   = DETFE**(-2.0D0/3.0D0)
      I1_BAR = TMP1 * (BE(1,1) + BE(2,2) + BE(3,3))
      PRESS  = (2.0D0/D1) * (DETFE - 1.0D0) * DETFE

      TRCE = I1_BAR / 3.0D0
      DO I = 1, 3
        DO J = 1, 3
          SIGMA(I,J) = (2.0D0*C10*TMP1*(BE(I,J) - TRCE*IDEN(I,J))
     1                  + PRESS * IDEN(I,J)) / DETFE
        END DO
      END DO

      I2_BAR = 0.0D0
      IF (MTYPE .GT. 0.5D0) THEN
        I1E = BE(1,1) + BE(2,2) + BE(3,3)
        I2E = 0.5D0*(I1E**2 - (BE(1,1)**2 + BE(2,2)**2 + BE(3,3)**2
     1        + 2.0D0*(BE(1,2)**2 + BE(1,3)**2 + BE(2,3)**2)))
        I2_BAR = TMP1**2 * I2E
        DO I = 1, 3
          DO J = 1, 3
            TMP2 = 0.0D0
            DO K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            END DO
            TEMP3(I,J) = I1_BAR * TMP1*BE(I,J) - TMP1**2 * TMP2
          END DO
        END DO
        TMP2 = (TEMP3(1,1) + TEMP3(2,2) + TEMP3(3,3)) / 3.0D0
        DO I = 1, 3
          DO J = 1, 3
            SIGMA(I,J) = SIGMA(I,J) +
     1        2.0D0*C01*(TEMP3(I,J) - TMP2*IDEN(I,J)) / DETFE
          END DO
        END DO
      END IF

C     --- Cauchy stress in Abaqus Voigt order: 11,22,33,12,13,23 ---
      SV(1) = SIGMA(1,1)
      SV(2) = SIGMA(2,2)
      SV(3) = SIGMA(3,3)
      SV(4) = SIGMA(1,2)
      SV(5) = SIGMA(1,3)
      SV(6) = SIGMA(2,3)

      DETFE_OUT = DETFE

C     --- Energies ---
      SSE_OUT = C10*(I1_BAR - 3.0D0) + (1.0D0/D1)*(DETFE - 1.0D0)**2
      IF (MTYPE .GT. 0.5D0) SSE_OUT = SSE_OUT + C01*(I2_BAR - 3.0D0)

      SPD_OUT = 0.0D0
      IF (ETA .GT. 1.0D-20) THEN
        DO I = 1, 3
          DO J = 1, 3
            SPD_OUT = SPD_OUT + TAU_DEV(I,J)**2
          END DO
        END DO
        SPD_OUT = SPD_OUT * DTIME_SAFE / (4.0D0 * ETA * DETFE)
      END IF

      RETURN
      END

C======================================================================
C     Utility: 3x3 matrix multiply C = A * B
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

C======================================================================
C     Utility: 3x3 matrix A*A^T  ->  C
C======================================================================
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

C======================================================================
C     Utility: Determinant of 3x3 matrix
C======================================================================
      SUBROUTINE MAT_DET3(A, DET)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3)
      DOUBLE PRECISION DET
      DET = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
     1    - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
     2    + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      RETURN
      END

C======================================================================
C     Utility: Inverse of 3x3 matrix (Cramer's rule)
C======================================================================
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
