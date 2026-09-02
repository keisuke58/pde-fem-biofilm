C=======================================================================
C  usermat_wrapper_v01_smoketest.f
C
C  A minimal ANSYS USERMAT entry point whose ONLY job is to call
C  BIOFILM_GROWTH_VISCO_V01 (biofilm_material_v01.f) -- the actual routine
C  meant for handover to Oliver's framework -- so it can be exercised in a
C  real ANSYS solve on v222 without touching any of Oliver's pool files.
C
C  This is a verification harness, not part of the deliverable: Oliver's
C  own framework will call BIOFILM_GROWTH_VISCO_V01 from their
C  Usermat_P21-V21_*.F at the AceGenNeoHookV04 call site (step 4 of
C  INTEGRATION_PLAN.md, their job). This file exists only to prove the
C  routine produces correct growth-driven stress inside a real v222 solve
C  ahead of that, per ROADMAP_2026.md Week 1 ("a working local v222 build
C  lets us run our own ANSYS jobs with the wrapper... without waiting for
C  step 4 at all").
C
C  State/property layout (deliberately matching t_growth_free.dat /
C  t_growth_baseclamped.dat's existing convention where possible, so those
C  decks need only the property block changed, not rewritten):
C    ustatev(1:9)  = Fv(3,3), row-major, prior viscous state
C    ustatev(10)   = alpha (growth driver)
C    prop(1) = sYoung (E, biofilm)     prop(2) = sNu
C    prop(3) = sEta (0 = elastic)      prop(4) = sC01Ratio (0 = neo-Hookean)
C    prop(5) = sMtype (unused by the core's neo-Hookean path, passed through)
C  sBiofilm is fixed at 1.0 (pure biofilm, no void blend) -- this harness
C  has no second material to blend against, unlike the real framework.
C=======================================================================
      subroutine usermat(
     &   matId, elemId, kDomIntPt, kLayer, kSectPt,
     &   ldstep, isubst, keycut,
     &   nDirect, nShear, ncomp, nStatev, nProp,
     &   Time, dTime, Temp, dTemp,
     &   stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &   Strain, dStrain, epsPl, prop, coords,
     &   var0, defGrad_t, defGrad, tsstif, epsZZ, cutFactor,
     &   var1, var2, var3, var4, var5, var6, var7, var8)

      implicit none
C     --- ANSYS USERMAT argument list (v222, same as usermat_biofilm.f --
C     see that file's header for the per-release trailing-argument note) ---
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
      double precision var0, var1, var2, var3, var4, var5, var6,
     &                 var7, var8

C     --- locals ---
      double precision sYoung, sNu, sEta, sC01Ratio, sMtype
      double precision sGrowth, sElasticWork
      double precision mFvN(3,3), mFvN1(3,3)
      double precision vCauchy(6), mTangCC(6,6)
      integer          sKeyCut, I, J, K

      sYoung    = prop(1)
      sNu       = prop(2)
      sEta      = prop(3)
      sC01Ratio = prop(4)
      sMtype    = prop(5)

      sGrowth = 0.0d0
      if (nStatev .ge. 10) sGrowth = ustatev(10)
      if (sGrowth .lt. 0.0d0) sGrowth = 0.0d0

      K = 0
      do I = 1, 3
        do J = 1, 3
          K = K + 1
          mFvN(I,J) = ustatev(K)
        end do
      end do
C     identity if never initialised (ustatev starts at 0 on the first call)
      if (abs(mFvN(1,1)) .lt. 1.0d-30 .and.
     &    abs(mFvN(2,2)) .lt. 1.0d-30 .and.
     &    abs(mFvN(3,3)) .lt. 1.0d-30) then
        do I = 1, 3
          do J = 1, 3
            mFvN(I,J) = 0.0d0
          end do
          mFvN(I,I) = 1.0d0
        end do
      end if

      call BIOFILM_GROWTH_VISCO_V01(
     &   defGrad, vCauchy, mTangCC,
     &   sYoung, sYoung, sNu, sNu, 1.0d0,
     &   sGrowth, mFvN, mFvN1,
     &   sEta, dTime, sC01Ratio, sMtype,
     &   sElasticWork, sKeyCut, elemId)

      if (sKeyCut .ne. 0) then
        keycut = 1
        cutFactor = 0.5d0
        return
      end if

      do I = 1, ncomp
        stress(I) = vCauchy(I)
        do J = 1, ncomp
          dsdePl(I,J) = mTangCC(I,J)
        end do
      end do

      K = 0
      do I = 1, 3
        do J = 1, 3
          K = K + 1
          ustatev(K) = mFvN1(I,J)
        end do
      end do
      if (nStatev .ge. 10) ustatev(10) = sGrowth

      sedEl = sElasticWork

      return
      end
