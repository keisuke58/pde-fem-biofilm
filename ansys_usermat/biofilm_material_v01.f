C=======================================================================
C  biofilm_material_v01.f
C
C  The growth + viscoelastic biofilm law, packaged in the shape Oliver's
C  framework calls its constitutive routines with, so it can be swapped in
C  at the AceGenNeoHookV04 call site in Usermat_P21-V21_*.F.
C
C  Their routine, for comparison:
C
C    SUBROUTINE AceGenNeoHookV04(v, mDefGrad, vCauchy, mTangCC,
C   &   sYoung, sYoungL, sNu, sNuL, sBiofilm, sAlpha, sElasticWork, sID)
C
C  Like theirs this is self-contained: no ANSYS includes, no common block,
C  no UPF calls. Deformation gradient and material constants in, Cauchy
C  stress and consistent tangent out. It therefore does not care which
C  ANSYS release it is linked into — only the surrounding usermat glue is
C  release-specific.
C
C  What it adds over AceGenNeoHookV04, which is purely elastic:
C
C    * growth kinematics  F = Fe.Fv.Fg  with  Fg = (1+alpha) I
C    * a single-step viscous update carrying Fv (see note 3 -- the repo
C      calls this backward Euler, but the flow increment is evaluated at
C      the old state, which is what bounds the usable time step)
C
C  The physics is NOT reimplemented here. This is a thin adapter around
C  BIOFILM_STRESS_CORE in usermat_biofilm.f, which is verified 0 ULP
C  against the Abaqus UMAT over 8017 states (crosscheck/, adversarial.py).
C  Keeping it a pure adapter is the point: the verification travels with
C  the routine instead of having to be redone.
C
C  ---------------------------------------------------------------------
C  Three interface notes for whoever wires this in
C  ---------------------------------------------------------------------
C
C  1. sGrowth is the GROWTH variable, i.e. the alpha of Fg=(1+alpha)I.
C     It is deliberately NOT called sAlpha: in AceGenNeoHookV04 the
C     argument named sAlpha is fed Sdp_sumLocal, a local biofilm average,
C     and is not a growth variable at all. Reusing the name would invite
C     exactly the wrong wiring.
C
C  2. The viscous state is Fv(3,3), NOT a 6-component Cauchy-Green.
C     Fv does not stay symmetric under this update -- measured at ~6e-5
C     relative asymmetry after 20 steps over 200 random states -- so it
C     cannot be reconstructed from Cv(6) without losing information.
C     That means 9 state slots where Vdp_Cv_n currently reserves 6.
C
C  3. sDt must resolve the viscous relaxation time eta/(2*C10), and this
C     routine ENFORCES it rather than only documenting it. The core
C     advances Fv with the flow increment evaluated at the old state, so
C     the step is only accurate while dt is small against that time; past
C     dt/tau ~ 0.5 the elastic strain overshoots through zero and the
C     stress changes sign, and past dt/tau ~ 1 it diverges.
C
C     Because the caller sets dt, leaving this to a comment would mean a
C     wrong step returns a plausible-looking wrong stress. Instead
C     dt/tau > DTMAX_RATIO sets sKeyCut = 1 and returns without computing
C     a stress, so the solver cuts the increment and retries. Note that
C     growth shrinks tau as C10 rises: a step that is fine at the start
C     of a solve can stop being fine later, so it is checked every call.
C
C     sEta = 0 selects the purely elastic path and removes the limit
C     entirely. Behaviour inside the valid range is unchanged.
C     Pinned in tests/test_material_wrapper.py.
C
C  On a cut-back (sKeyCut = 1) every output is still defined: stress and
C  tangent are zeroed and mFvN1 is returned unchanged, so a caller that
C  reads them before checking the flag gets zeros rather than whatever
C  was in its work array.
C
C  Voigt order is ANSYS (11,22,33,12,23,13), matching what usermat passes
C  straight through to the constitutive call.
C
C  Build: plain Fortran, no includes. Compile alongside usermat_biofilm.f
C  (which supplies BIOFILM_STRESS_CORE).
C=======================================================================
      subroutine BIOFILM_GROWTH_VISCO_V01(
     &   mDefGrad, vCauchy, mTangCC,
     &   sYoung, sYoungL, sNu, sNuL, sBiofilm,
     &   sGrowth, mFvN, mFvN1,
     &   sEta, sDt, sC01Ratio, sMtype,
     &   sElasticWork, sKeyCut, sID)

      implicit none

C     --- arguments -------------------------------------------------
C     in
      double precision mDefGrad(3,3)      ! total deformation gradient F
      double precision sYoung, sYoungL    ! E of biofilm / of void ("leer")
      double precision sNu, sNuL          ! nu of biofilm / of void
      double precision sBiofilm           ! biofilm fraction, blends the two
      double precision sGrowth            ! growth alpha, Fg = (1+alpha) I
      double precision mFvN(3,3)          ! viscous state at t_n
      double precision sEta               ! viscosity (0 = elastic)
      double precision sDt                ! time increment
      double precision sC01Ratio          ! C01/C10 (0 => neo-Hookean)
      double precision sMtype             ! 0 neo-Hookean, 1 Mooney-Rivlin
      integer          sID                ! element/point id, for debug only
C     out
      double precision vCauchy(6)         ! Cauchy stress, ANSYS Voigt
      double precision mTangCC(6,6)       ! consistent tangent
      double precision mFvN1(3,3)         ! viscous state at t_n+1
      double precision sElasticWork
      integer          sKeyCut            ! 1 => ask the solver to cut back

C     --- locals ----------------------------------------------------
      double precision C10, C01, D1, MU, BULK, EBLEND, NUBLEND
      double precision FG_INV(3,3), FGSC
      double precision SV0(6), SVP(6), DFP(3,3), FV_DUM(3,3)
      double precision SSE_C, SPD_C, DETFE, PERT, SYMF
      double precision TRELAX, DTRAT
      integer          I, J, K, P, Q, IP, JP

C     Largest dt/tau_relax this routine will answer for. Past this the
C     explicit flow increment overshoots through zero and the stress comes
C     back with the wrong sign (see note 3 in the header), so returning a
C     number would be worse than asking the solver to cut back.
C
C     0.5 is where the sign flip sets in, measured. It is deliberately not
C     tighter: below it the answer degrades in accuracy but stays
C     qualitatively right, and how much accuracy to buy with step size is
C     the caller's engineering choice, not ours to force. Cutting back on
C     merely-inaccurate steps would also risk stalling a solve. Lower it
C     here if you want the routine to insist on more resolution.
      double precision DTMAX_RATIO
      parameter       (DTMAX_RATIO = 0.5d0)
      integer          VI(6), VJ(6)
      data VI /1, 2, 3, 1, 2, 1/
      data VJ /1, 2, 3, 2, 3, 3/

      sKeyCut = 0

C     --- define every output up front. A cut-back returns early, and the
C         caller is entitled to read vCauchy/mTangCC even then: their
C         framework passes work arrays that are not zeroed between calls,
C         so leaving them untouched hands back whatever was in memory.
      do I = 1, 6
        vCauchy(I) = 0.0d0
        do J = 1, 6
          mTangCC(I,J) = 0.0d0
        end do
      end do
      do I = 1, 3
        do J = 1, 3
          mFvN1(I,J) = mFvN(I,J)
        end do
      end do
      sElasticWork = 0.0d0

C     --- blend biofilm against void, the same convention their deck
C         carries as YOUNG_BIO / YOUNG_VOID and POISSON_BIO / POISSON_VOID.
C         A linear blend in the biofilm fraction; if their AceGen routine
C         uses a different interpolation this is the one line to match.
      EBLEND  = sBiofilm * sYoung + (1.0d0 - sBiofilm) * sYoungL
      NUBLEND = sBiofilm * sNu    + (1.0d0 - sBiofilm) * sNuL
      if (EBLEND .lt. 1.0d-12) EBLEND = 1.0d-12
      if (NUBLEND .gt.  0.49999d0) NUBLEND =  0.49999d0
      if (NUBLEND .lt. -0.99999d0) NUBLEND = -0.99999d0

C     --- (E, nu) -> (C10, C01, D1), the small-strain-consistent map
C         mu = 2(C10+C01),  K = 2/D1
      MU   = EBLEND / (2.0d0 * (1.0d0 + NUBLEND))
      BULK = EBLEND / (3.0d0 * (1.0d0 - 2.0d0 * NUBLEND))
      if (BULK .lt. 1.0d-20) BULK = 1.0d-20
      C10 = 0.5d0 * MU / (1.0d0 + sC01Ratio)
      C01 = C10 * sC01Ratio
      D1  = 2.0d0 / BULK

C     --- Fg = (1+alpha) I ; inv(Fg) ---
      FGSC = max(1.0d0 + sGrowth, 1.0d-15)
      do I = 1, 3
        do J = 1, 3
          FG_INV(I,J) = 0.0d0
        end do
        FG_INV(I,I) = 1.0d0 / FGSC
      end do

C     --- refuse a step that does not resolve the viscous relaxation time.
C         The caller sets dt, so this cannot be left to documentation: the
C         failure is silent and severe -- the stress changes sign rather
C         than erroring. sKeyCut is exactly the channel for this; the
C         solver cuts the increment and retries at a step that works.
C
C         Elastic runs (sEta = 0) have no relaxation time and are never
C         restricted. Growth shrinks tau as C10 rises, so a step that was
C         fine early in a solve can stop being fine later -- which is why
C         this is checked every call rather than once.
      if (sEta .gt. 1.0d-20 .and. C10 .gt. 1.0d-20) then
        TRELAX = sEta / (2.0d0 * C10)
        DTRAT  = sDt / TRELAX
        if (DTRAT .gt. DTMAX_RATIO) then
          sKeyCut = 1
          do I = 1, 3
            do J = 1, 3
              mFvN1(I,J) = mFvN(I,J)
            end do
          end do
          return
        end if
      end if

C     --- base stress and the updated viscous state ---
      call BIOFILM_STRESS_CORE(mDefGrad, FG_INV, mFvN,
     &     C10, C01, D1, sEta, sMtype, sDt, VI, VJ,
     &     SV0, mFvN1, SSE_C, SPD_C, DETFE)

C     cut-back if the elastic Jacobian collapsed. Restore mFvN1: the core
C     wrote its own update into it above, and advancing the viscous state
C     off a collapsed configuration would corrupt the state the solver
C     restarts the cut increment from.
      if (DETFE .le. 1.0d-12) then
        sKeyCut = 1
        do I = 1, 3
          do J = 1, 3
            mFvN1(I,J) = mFvN(I,J)
          end do
        end do
        return
      end if

      do I = 1, 6
        vCauchy(I) = SV0(I)
      end do
      sElasticWork = SSE_C

C     --- consistent tangent by F-perturbation (Sun et al. 2008), the
C         same scheme and step usermat_biofilm.f uses, so the tangent
C         matches the one already confirmed to converge under
C         SOLID185/NLGEOM,ON.
      PERT = 1.0d-7
      do P = 1, 6
        IP = VI(P)
        JP = VJ(P)
        do I = 1, 3
          do J = 1, 3
            DFP(I,J) = mDefGrad(I,J)
          end do
        end do
        if (IP .eq. JP) then
          do K = 1, 3
            DFP(IP,K) = DFP(IP,K) + PERT*mDefGrad(IP,K)
          end do
        else
          SYMF = 0.5d0
          do K = 1, 3
            DFP(IP,K) = DFP(IP,K) + SYMF*PERT*mDefGrad(JP,K)
            DFP(JP,K) = DFP(JP,K) + SYMF*PERT*mDefGrad(IP,K)
          end do
        end if
        call BIOFILM_STRESS_CORE(DFP, FG_INV, mFvN,
     &       C10, C01, D1, sEta, sMtype, sDt, VI, VJ,
     &       SVP, FV_DUM, SSE_C, SPD_C, DETFE)
        do Q = 1, 6
          mTangCC(Q,P) = (SVP(Q) - SV0(Q)) / PERT
        end do
      end do

      return
      end
