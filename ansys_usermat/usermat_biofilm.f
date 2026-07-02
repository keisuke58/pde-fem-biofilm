C=======================================================================
C  usermat_biofilm.f  —  ANSYS USERMAT port of the biofilm growth /
C  viscoelastic constitutive law (F = Fe . Fv . Fg, Fg = (1+alpha) I).
C
C  PURPOSE (master-thesis starting point)
C  --------------------------------------
C  Port of the *verified* Abaqus UMAT (umat_biofilm_visco.f /
C  umat_biofilm_visco_phase2.f) to the ANSYS Mechanical APDL USERMAT
C  interface, so Felix/IKM's existing ANSYS FE model can call the same
C  material law at each Gauss point in place of its phenomenological
C  constitutive model.  The constitutive algebra below MIRRORS the
C  Abaqus core exactly (same Neo-Hookean deviator + D1 pressure, same
C  backward-Euler viscous update, same F-perturbation consistent
C  tangent) — only the *interface* (argument list, storage order,
C  cut-back signalling, tangent convention) is ANSYS-specific.
C
C  A clearly-marked hook (see "PYTHON MATERIAL HOOK") shows where the
C  per-Gauss-point call to the Python material model would replace the
C  inline Fortran law (via ISO_C_BINDING / local socket).
C
C  STATUS: syntax-checked with gfortran; NOT yet run inside ANSYS.
C  The constitutive core it reproduces IS verified (tangent vs FD
C  ~2.4e-8; patch tests 13/13 in phase2_patch_test.py).
C
C  ABAQUS UMAT  <->  ANSYS USERMAT  mapping (key porting knowledge)
C  ---------------------------------------------------------------
C   Abaqus                         ANSYS USERMAT
C   DFGRD1 / DFGRD0                defGrad / defGrad_t      (3x3 F)
C   STRESS(NTENS)                 stress(ncomp)            (Cauchy)
C   DDSDDE(NTENS,NTENS)           dsdePl(ncomp,ncomp)      (material Jacobian)
C   STATEV(NSTATV)                ustatev(nStatev)
C   PROPS(NPROPS)                 prop(nProp)
C   DTIME                         dTime
C   PNEWDT < 1  (cut-back)        keycut = 1 (+ cutFactor)
C   SSE / SPD                     sedEl / sedPl
C   comp. order 11,22,33,12,13,23 ANSYS order 11,22,33,12,23,13  (<-- 5<->6 swap!)
C
C  PROPS (prop):
C    prop(1) = C10   instantaneous shear modulus [stress]
C    prop(2) = C01   Mooney-Rivlin C01 (0 = Neo-Hookean)
C    prop(3) = D1    compressibility [1/stress]
C    prop(4) = eta   viscosity [stress*time]  (0 = elastic)
C    prop(5) = mtype 0=Neo-Hookean, 1=Mooney-Rivlin
C    prop(6) = kUsePy 0=inline Fortran law, 1=call Python material hook
C
C  STATE (ustatev), NSTATV = 10:
C    ustatev(1:9) = Fv  (row-major 3x3)
C    ustatev(10)  = alpha  (accumulated volumetric growth; growth driver,
C                   set from the JAXFEM alpha-field mapped to this IP, or
C                   evolved via a user field / TB,STATE table)
C=======================================================================
      subroutine usermat(
     &   matId, elemId, kDomIntPt, kLayer, kSectPt,
     &   ldstep, isubst, keycut,
     &   nDirect, nShear, ncomp, nStatev, nProp,
     &   Time, dTime, Temp, dTemp,
     &   stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &   Strain, dStrain, epsPl, prop, coords,
     &   defGrad_t, defGrad, tsstif, epsZZ, cutFactor)

      implicit none
C     --- ANSYS USERMAT argument list (3-D / plane-strain solid) ---
      integer          matId, elemId, kDomIntPt, kLayer, kSectPt,
     &                 ldstep, isubst, keycut,
     &                 nDirect, nShear, ncomp, nStatev, nProp
      double precision Time, dTime, Temp, dTemp, sedEl, sedPl, epseq,
     &                 epsZZ, cutFactor
      double precision stress(ncomp), ustatev(nStatev),
     &                 dsdePl(ncomp,ncomp), Strain(ncomp),
     &                 dStrain(ncomp), epsPl(ncomp), prop(nProp),
     &                 coords(3), defGrad_t(3,3), defGrad(3,3),
     &                 tsstif(2)

C     --- locals ---
      double precision C10, C01, D1, ETA, MTYPE, KUSEPY
      double precision ALPHA, FGSC, FG_INV(3,3)
      double precision FV_OLD(3,3), FV_NEW(3,3), FV_DUM(3,3)
      double precision SV0(6), SVP(6), DFP(3,3)
      double precision SSE_C, SPD_C, DETFE, PERT, SYMF
      integer          I, J, K, P, Q, IP, JP
C     ANSYS Voigt map (11,22,33,12,23,13)
      integer          VI(6), VJ(6)
      data VI /1, 2, 3, 1, 2, 1/
      data VJ /1, 2, 3, 2, 3, 3/

C     --- material properties ---
      C10    = prop(1)
      C01    = prop(2)
      D1     = prop(3)
      ETA    = prop(4)
      MTYPE  = prop(5)
      KUSEPY = 0.0d0
      if (nProp .ge. 6) KUSEPY = prop(6)

C     --- growth driver + viscous state from ustatev ---
      ALPHA = 0.0d0
      if (nStatev .ge. 10) ALPHA = ustatev(10)
      if (ALPHA .lt. 0.0d0) ALPHA = 0.0d0
      K = 0
      do I = 1, 3
        do J = 1, 3
          K = K + 1
          FV_OLD(I,J) = ustatev(K)
        end do
      end do
      call INIT_FV_IF_ZERO(FV_OLD)

C     --- Fg = (1+alpha) I ;  inv(Fg) ---
      FGSC = max(1.0d0 + ALPHA, 1.0d-15)
      do I = 1, 3
        do J = 1, 3
          FG_INV(I,J) = 0.0d0
        end do
        FG_INV(I,I) = 1.0d0 / FGSC
      end do

C=======================================================================
C  PYTHON MATERIAL HOOK  (per Gauss point)
C  ---------------------------------------
C  When kUsePy=1 the constitutive response would be delegated to the
C  Python material model (the paper's calibrated law) instead of the
C  inline Fortran core below.  Intended mechanism: an ISO_C_BINDING /
C  local-socket bridge that sends (defGrad, FV_OLD, alpha, dTime, prop)
C  and receives (Cauchy stress, FV_NEW, dsdePl).  Left as an explicit
C  extension point for the thesis; the inline core is the reference /
C  fallback used for verification.
C
C      if (KUSEPY .gt. 0.5d0) then
C        call BIOFILM_PY_MATERIAL(defGrad, FV_OLD, ALPHA, dTime, prop,
C     &                           nProp, SV0, FV_NEW, dsdePl, ncomp,
C     &                           VI, VJ, keycut)
C        goto 900
C      end if
C=======================================================================

C     --- base stress (ANSYS Voigt order) ---
      call BIOFILM_STRESS_CORE(defGrad, FG_INV, FV_OLD,
     &     C10, C01, D1, ETA, MTYPE, dTime, VI, VJ,
     &     SV0, FV_NEW, SSE_C, SPD_C, DETFE)

C     cut-back if the elastic Jacobian collapsed (Je -> 0)
      if (DETFE .le. 1.0d-12) then
        keycut   = 1
        cutFactor = 0.5d0
        return
      end if

      do I = 1, ncomp
        stress(I) = SV0(I)
      end do
      sedEl = SSE_C
      sedPl = SPD_C

C     --- consistent tangent dsdePl by F-perturbation (Sun et al. 2008) ---
      PERT = 1.0d-7
      do P = 1, ncomp
        IP = VI(P)
        JP = VJ(P)
        do I = 1, 3
          do J = 1, 3
            DFP(I,J) = defGrad(I,J)
          end do
        end do
        if (IP .eq. JP) then
          do K = 1, 3
            DFP(IP,K) = DFP(IP,K) + PERT*defGrad(IP,K)
          end do
        else
          SYMF = 0.5d0
          do K = 1, 3
            DFP(IP,K) = DFP(IP,K) + SYMF*PERT*defGrad(JP,K)
            DFP(JP,K) = DFP(JP,K) + SYMF*PERT*defGrad(IP,K)
          end do
        end if
        call BIOFILM_STRESS_CORE(DFP, FG_INV, FV_OLD,
     &       C10, C01, D1, ETA, MTYPE, dTime, VI, VJ,
     &       SVP, FV_DUM, SSE_C, SPD_C, DETFE)
        do Q = 1, ncomp
          dsdePl(Q,P) = (SVP(Q) - SV0(Q)) / PERT
        end do
      end do

C     --- store updated viscous state ---
      K = 0
      do I = 1, 3
        do J = 1, 3
          K = K + 1
          ustatev(K) = FV_NEW(I,J)
        end do
      end do

      return
      end

C=======================================================================
C  BIOFILM_STRESS_CORE — growth + (Mooney-Rivlin) Neo-Hookean + 1 viscous
C  branch; returns Cauchy stress in ANSYS Voigt order via (VI,VJ).
C  Mirrors the verified Abaqus core (umat_biofilm_visco*.f).
C=======================================================================
      subroutine BIOFILM_STRESS_CORE(DFG, FG_INV, FV_OLD,
     &     C10, C01, D1, ETA, MTYPE, DT, VI, VJ,
     &     SV, FV_NEW, SSE_OUT, SPD_OUT, DETFE_OUT)
      implicit none
      double precision DFG(3,3), FG_INV(3,3), FV_OLD(3,3), FV_NEW(3,3)
      double precision C10, C01, D1, ETA, MTYPE, DT
      integer          VI(6), VJ(6)
      double precision SV(6), SSE_OUT, SPD_OUT, DETFE_OUT

      double precision FMECH(3,3), FV_INV(3,3), FE(3,3), BE(3,3)
      double precision TAU(3,3), SIG(3,3), IDEN(3,3), T3(3,3)
      double precision DETFV, DETFE, TMP1, TMP2, TRCE, I1B, PRESS, DTS
      integer          I, J, K

      do I = 1, 3
        do J = 1, 3
          IDEN(I,J) = 0.0d0
        end do
        IDEN(I,I) = 1.0d0
      end do

C     F_mech = F . inv(Fg)
      call MAT_MULT(DFG, FG_INV, FMECH, 3)

C     trial Fe = F_mech . inv(Fv_old)
      call MAT_INV3(FV_OLD, FV_INV, DETFV)
      call MAT_MULT(FMECH, FV_INV, FE, 3)
      call MAT_AAT(FE, BE, 3)
      call MAT_DET3(FE, DETFE)
      if (DETFE .lt. 1.0d-15) DETFE = 1.0d-15
      TMP1 = DETFE**(-2.0d0/3.0d0)
      I1B  = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      TRCE = I1B/3.0d0
      do I = 1, 3
        do J = 1, 3
          TAU(I,J) = 2.0d0*C10*TMP1*(BE(I,J)-TRCE*IDEN(I,J))
        end do
      end do

C     viscous update (backward Euler)
      DTS = max(DT, 1.0d-20)
      if (ETA .gt. 1.0d-20) then
        TMP1 = DTS/(2.0d0*ETA*DETFE)
        do I = 1, 3
          do J = 1, 3
            T3(I,J) = IDEN(I,J) + TMP1*TAU(I,J)
          end do
        end do
        call MAT_MULT(T3, FV_OLD, FV_NEW, 3)
      else
        do I = 1, 3
          do J = 1, 3
            FV_NEW(I,J) = FV_OLD(I,J)
          end do
        end do
      end if

C     recompute Fe with updated Fv -> Cauchy stress
      call MAT_INV3(FV_NEW, FV_INV, DETFV)
      call MAT_MULT(FMECH, FV_INV, FE, 3)
      call MAT_AAT(FE, BE, 3)
      call MAT_DET3(FE, DETFE)
      if (DETFE .lt. 1.0d-15) DETFE = 1.0d-15
      DETFE_OUT = DETFE
      TMP1  = DETFE**(-2.0d0/3.0d0)
      I1B   = TMP1*(BE(1,1)+BE(2,2)+BE(3,3))
      TRCE  = I1B/3.0d0
      PRESS = (2.0d0/D1)*(DETFE-1.0d0)*DETFE
      do I = 1, 3
        do J = 1, 3
          SIG(I,J) = (2.0d0*C10*TMP1*(BE(I,J)-TRCE*IDEN(I,J))
     &               + PRESS*IDEN(I,J)) / DETFE
        end do
      end do

C     Mooney-Rivlin correction
      if (MTYPE .gt. 0.5d0) then
        do I = 1, 3
          do J = 1, 3
            TMP2 = 0.0d0
            do K = 1, 3
              TMP2 = TMP2 + BE(I,K)*BE(K,J)
            end do
            T3(I,J) = I1B*TMP1*BE(I,J) - TMP1**2*TMP2
          end do
        end do
        TMP2 = (T3(1,1)+T3(2,2)+T3(3,3))/3.0d0
        do I = 1, 3
          do J = 1, 3
            SIG(I,J) = SIG(I,J)
     &        + 2.0d0*C01*(T3(I,J)-TMP2*IDEN(I,J))/DETFE
          end do
        end do
      end if

C     pack Cauchy stress into ANSYS Voigt order (VI,VJ)
      do I = 1, 6
        SV(I) = SIG(VI(I), VJ(I))
      end do

C     energies
      SSE_OUT = C10*(I1B-3.0d0) + (1.0d0/D1)*(DETFE-1.0d0)**2
      SPD_OUT = 0.0d0
      if (ETA .gt. 1.0d-20) then
        TMP2 = 0.0d0
        do I = 1, 3
          do J = 1, 3
            TMP2 = TMP2 + TAU(I,J)**2
          end do
        end do
        SPD_OUT = TMP2*DTS/(4.0d0*ETA*DETFE)
      end if
      return
      end

C=======================================================================
C  utilities (self-contained; ANSYS build normally adds impcom.inc)
C=======================================================================
      subroutine INIT_FV_IF_ZERO(FV)
      implicit none
      double precision FV(3,3), TR
      integer I, J
      TR = FV(1,1)+FV(2,2)+FV(3,3)
      if (abs(TR) .lt. 1.0d-10) then
        do I = 1, 3
          do J = 1, 3
            FV(I,J) = 0.0d0
          end do
          FV(I,I) = 1.0d0
        end do
      end if
      return
      end

      subroutine MAT_MULT(A, B, C, N)
      implicit none
      integer N, I, J, K
      double precision A(N,N), B(N,N), C(N,N)
      do I = 1, N
        do J = 1, N
          C(I,J) = 0.0d0
          do K = 1, N
            C(I,J) = C(I,J) + A(I,K)*B(K,J)
          end do
        end do
      end do
      return
      end

      subroutine MAT_AAT(A, C, N)
      implicit none
      integer N, I, J, K
      double precision A(N,N), C(N,N)
      do I = 1, N
        do J = 1, N
          C(I,J) = 0.0d0
          do K = 1, N
            C(I,J) = C(I,J) + A(I,K)*A(J,K)
          end do
        end do
      end do
      return
      end

      subroutine MAT_DET3(A, DET)
      implicit none
      double precision A(3,3), DET
      DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     &    - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     &    + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      return
      end

      subroutine MAT_INV3(A, AINV, DET)
      implicit none
      double precision A(3,3), AINV(3,3), DET, DINV
      integer I, J
      call MAT_DET3(A, DET)
      if (abs(DET) .lt. 1.0d-30) then
        do I = 1, 3
          do J = 1, 3
            AINV(I,J) = 0.0d0
          end do
          AINV(I,I) = 1.0d0
        end do
        return
      end if
      DINV = 1.0d0/DET
      AINV(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))*DINV
      AINV(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3))*DINV
      AINV(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))*DINV
      AINV(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))*DINV
      AINV(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))*DINV
      AINV(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3))*DINV
      AINV(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))*DINV
      AINV(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2))*DINV
      AINV(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))*DINV
      return
      end
