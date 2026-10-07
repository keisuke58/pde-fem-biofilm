C     --- phi mode, executable part (pde-fem-biofilm, 2026-10-01) -----
C     Target: Usermat_P21-V21_v222.F (the file the 9/7 ANSYS.exe was
C     built from). Paste AFTER the ecology block has finished setting
C     Sbio_GrowthConst, immediately BEFORE the call to
C     BIOFILM_GROWTH_VISCO_V01.
C
C     Klempt et al. 2024, Eq. 36: alpha_dot = k_alpha * phi, nothing
C     added. The switch is prop(28), NOT prop(6): in v222.F any
C     prop(6) /= 0 already means the ecology mode. Keep prop(6) = 0.
C       prop(28) = 1: one species,  phi = bio1
C       prop(28) = 2: two species,  phi = bio1 + bio2
C       prop(28) = 4: the POINT MODEL (Hamilton ecology ODE, up to five
C                     species; set the number with the server's
C                     --active-species) at this Gauss point. One macro
C                     step is split into ceiling(dTime/1e-4) inner steps
C                     inside one call; the 12-component inner state is
C                     carried across macro steps in ustatev(72:83), and
C                     alpha += k_alpha * phi_int, phi_int being the
C                     integral of sum(phi_i psi_i) over the inner steps.
C                     theta in prop(8:27). The partner's prop(6) ecology
C                     block is left alone: keep prop(6) = 0.
C       prop(28) = 8: COMPOSITION from an independent point model that
C                     starts when biofilm arrives (no scaling to the
C                     field); see the branch below.
C       prop(28) = 7: COMPOSITION from the point model, amount from the
C                     3D field (unchanged): see the branch below.
C       prop(28) = 5, 6: REACTION STEP of a transport/reaction split
C                     for one (5) or two (6) species: phi_i from the
C                     3D field, the point model gives the reaction
C                     rate R_i, handed to USSFin through biofilm_split
C                     (split_rates.f); alpha from Eq. 36 with
C                     phi = bio1 (+ bio2). prop(29) = phi_min.
C       prop(28) = 3: the partner's OWN growth. USSFin already integrates
C                     Eq. 36 as locbio = alpha_K (Klempt convention,
C                     starting at 1; rate K_LOCAL1, explicit in bio), and
C                     Sdp_sumLocal = (locbio1 + locbio2)/2 is what their
C                     AceGen law uses as Fg = alpha_K I. Our alpha is
C                     alpha_K - 1. Found on IKMHIWI03, 1 Oct (stage 3).
C       prop(7)  = k_alpha (same slot the ecology mode uses for it)
C     Needs TB,USER with 28 constants; the nProp guard below keeps a
C     27-constant deck from reading past the end of prop.
C     ustatev(84) holds the CONVERGED alpha: ANSYS hands it back
C     unchanged on every equilibrium iteration, so growth is counted
C     once per increment.
C     If stage 1 shows locbio is the Gauss-point value, change bio to
C     locbio in the two PHI_USED lines -- nothing else.
      IF (nProp .GE. 28) THEN
      IF (prop(28) .GT. 0.5D0) THEN
        IF ((prop(28) .GT. 1.5D0 .AND. prop(28) .LT. 2.5D0) .OR.
     &      prop(28) .GT. 5.5D0) THEN
          PHI_USED = Sdp_bio1_n + Sdp_bio2_n
        ELSE
          PHI_USED = Sdp_bio1_n
        END IF
        ALPHA_N = ustatev(84)
        IF (prop(28) .GT. 7.5D0) THEN
C         --- composition from an independent point model -----------
C         Klempt et al. 2026 formulate a material point model and call
C         phi0 a numerical auxiliary, not a physical void. Here the
C         point model runs on its own at each Gauss point, NOT scaled
C         to the field, from the substep in which the field first
C         reaches phi_min (prop(29)). It supplies only the composition
C         chi_i = phi_i / sum(phi); the amount stays the field's
C         (Klempt 2024), and growth is Eq. 36 with phi_3D.
C         prop(30), prop(32) = phi_1, phi_2 at the start (default
C         0.2, 0.2); prop(31) = s, the point model's clock.
          SP_SUM = MIN(MAX(PHI_USED, 0.0D0), SP_PHIMAX)
          DO PM_I = 1, 12
            PM_GOLD(PM_I) = ustatev(71 + PM_I)
          END DO
          SP_MIN = 0.0D0
          IF (nProp .GE. 29) SP_MIN = prop(29)
          DO PM_I = 1, 20
            PM_TH(PM_I) = prop(7 + PM_I)
          END DO
          CM_TS = 1.0D0
          IF (nProp .GE. 31) THEN
            IF (prop(31) .GT. 0.0D0) CM_TS = prop(31)
          END IF
          CM_DT = CM_TS * dTime
          PM_NSUB = 1
          IF (CM_DT .GT. PM_DTMAX) PM_NSUB = CEILING(CM_DT / PM_DTMAX)
          SP_ID = SPLIT_ID(elemId, kDomIntPt)
          PM_OK = .TRUE.
          PM_PHIINT = 0.0D0
          IF (PM_GOLD(7) .EQ. 0.0D0 .AND. SP_SUM .LT. SP_MIN) THEN
C           no biofilm here yet: the point model has not started
            SP_HITI = 2
            DO PM_I = 1, 12
              PM_GNEW(PM_I) = PM_GOLD(PM_I)
            END DO
          ELSE
            IF (PM_GOLD(7) .EQ. 0.0D0) THEN
C             biofilm has arrived: start from the case's initial state
              DO PM_I = 1, 12
                PM_GOLD(PM_I) = 0.0D0
              END DO
              PM_GOLD(1) = 0.2D0
              PM_GOLD(2) = 0.2D0
              IF (nProp .GE. 30) THEN
                IF (prop(30) .GT. 0.0D0) PM_GOLD(1) = prop(30)
              END IF
              IF (nProp .GE. 32) THEN
                IF (prop(32) .GT. 0.0D0) PM_GOLD(2) = prop(32)
              END IF
              PM_GOLD(6) = 1.0D0 - (PM_GOLD(1) + PM_GOLD(2))
              PM_GOLD(7) = 0.999D0
              PM_GOLD(8) = 0.999D0
            END IF
            CALL SPLIT_CACHE_GET(SP_ID, ldstep, isubst, PM_GOLD,
     &                           PM_GNEW, PM_PHIINT, SP_HIT)
            IF (SP_HIT) THEN
              SP_HITI = 1
            ELSE
              SP_HITI = 0
              CALL biofilm_ecology_hook(PM_GOLD, PM_TH, CM_DT,
     &             PM_NSUB, PM_GNEW, PM_PHIINT, PM_OK)
              IF (PM_OK) THEN
                DO PM_I = 1, 12
                  IF (PM_GNEW(PM_I) .NE. PM_GNEW(PM_I)) PM_OK = .FALSE.
                END DO
                IF (PM_GNEW(1) + PM_GNEW(2) .LE. 0.0D0) PM_OK = .FALSE.
              END IF
              IF (PM_OK) CALL SPLIT_CACHE_PUT(SP_ID, ldstep, isubst,
     &                        PM_GOLD, PM_GNEW, PM_PHIINT)
            END IF
          END IF
          IF (PM_OK) THEN
            DO PM_I = 1, 12
              ustatev(71 + PM_I) = PM_GNEW(PM_I)
            END DO
          ELSE
            keycut = 1
          END IF
          CALL BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI_USED, prop(7),
     &                                dTime, ALPHA_NEW)
          IF (kDomIntPt .EQ. 1 .AND.
     &        (MOD(elemId, TRACE_STRIDE) .EQ. 1 .OR.
     &         PHI_USED .GE. TRACE_PHI_MIN)) THEN
            IF (.NOT. AG_TRACE_OPEN) THEN
              OPEN(95, FILE='age_trace.csv', STATUS='UNKNOWN',
     &             POSITION='APPEND')
              AG_TRACE_OPEN = .TRUE.
            END IF
            WRITE(95,'(I0,5(",",I0),6(",",ES24.16E3),'//
     &            '24(",",ES24.16E3))')
     &        elemId, kDomIntPt, ldstep, isubst, PM_NSUB, SP_HITI,
     &        dTime, CM_DT, PHI_USED, SP_SUM, ALPHA_N, ALPHA_NEW,
     &        (PM_GOLD(PM_I), PM_I = 1, 12),
     &        (PM_GNEW(PM_I), PM_I = 1, 12)
            FLUSH(95)
          END IF
        ELSE IF (prop(28) .GT. 6.5D0) THEN
C         --- composition from the point model (Klempt 2026) --------
C         The 3D field (Klempt 2024) decides the AMOUNT phi_3D = bio1
C         + bio2 and is not touched. The point model decides only the
C         COMPOSITION chi_i = phi_i / sum(phi) and psi: its state is
C         rescaled to sum(phi_i) = phi_3D before and after each inner
C         integration (phi0 = 1 - phi_3D). Growth: Eq. 36 with phi_3D.
C         prop(29) = phi_min: where phi_3D < phi_min the point model is
C         not run (no server call, no rescaling) and composition, psi
C         and gamma are held; growth (Eq. 36) is not affected.
C         prop(30) = chi_1 at the first call (default 0.5).
C         prop(31) = s, the point model's clock: it advances s * dTime
C         per substep (default 1). No paper links Klempt 2024's T* to
C         Klempt 2026's time unit; s is that assumption, made explicit.
C         prop(32) = phi_cap: the amount the point model sees is
C         min(phi_3D, phi_cap) (default 1 - 1e-6). In fully packed
C         regions phi0 would otherwise be pinned below the point model's
C         own packing limit and the composition would depend on the
C         coupling interval; at phi_cap = 0.9 it converges.
C         prop(33) = c_ref > 0: local nutrient in the point model
C         (two-way step 1, ROADMAP_TWO_WAY.md). The server uses
C         c* = c*_0 * c_rel with c_rel = min(max(CM_NUT/c_ref, 0), 1).
C         CM_NUT is the field's nutrient at this Gauss point; the next
C         executable line is its source, and paste_fragments.py
C         --nut-var NAME replaces it with CM_NUT = NAME. As shipped,
C         CM_NUT = -1D30 (not available) and prop(33) has no effect.
C         The field's c itself can be < 0 (zero-order consumption), so
C         "not available" is -1D30, not -1: a negative c gives c_rel = 0
C         (5 Oct: with -1 as the mark, c < 0 fell back to c_rel = 1).
C         prop(37) = n, the number of species (2..5, default 2; 5 Oct):
C         the point model's state holds phi_1..phi_n, psi_1..psi_n. Two
C         species start from chi_1 = prop(30) / prop(35); more species
C         start from equal shares 1/n. The server must run the same case
C         (--case, e.g. 4sp_case1, sets n there).
C         prop(35) > 0: chi_1 at the first call for points the biofilm
C         reaches later (phi_3D below phi_min at their first call); the
C         points that hold biofilm from the start keep prop(30). Without
C         it prop(30) sets both, so the composition behind the front and
C         the seed's start cannot be varied apart (5 Oct runs).
          CM_NUT = -1.0D30
          CM_NSP = 2
          IF (nProp .GE. 37) THEN
            IF (prop(37) .GT. 0.0D0) CM_NSP = MIN(MAX(NINT(prop(37)),
     &                                          2), 5)
          END IF
          CM_CAP = SP_PHIMAX
          IF (nProp .GE. 32) THEN
            IF (prop(32) .GT. 0.0D0) CM_CAP = MIN(prop(32), SP_PHIMAX)
          END IF
          SP_SUM = MIN(MAX(PHI_USED, 0.0D0), CM_CAP)
          DO PM_I = 1, 12
            PM_GOLD(PM_I) = ustatev(71 + PM_I)
          END DO
          SP_MIN = 0.0D0
          IF (nProp .GE. 29) SP_MIN = prop(29)
          IF (PM_GOLD(7) .EQ. 0.0D0) THEN
            CM_CHI1 = 0.5D0
            IF (nProp .GE. 30) CM_CHI1 = prop(30)
            IF (nProp .GE. 35) THEN
              IF (prop(35) .GT. 0.0D0 .AND. (SP_SUM .LT. SP_MIN
     &            .OR. SP_SUM .LE. 0.0D0)) CM_CHI1 = prop(35)
            END IF
            DO PM_I = 1, 12
              PM_GOLD(PM_I) = 0.0D0
            END DO
            IF (CM_NSP .EQ. 2) THEN
              PM_GOLD(1) = CM_CHI1
              PM_GOLD(2) = 1.0D0 - CM_CHI1
            ELSE
              DO PM_I = 1, CM_NSP
                PM_GOLD(PM_I) = 1.0D0 / DBLE(CM_NSP)
              END DO
            END IF
            DO PM_I = 1, CM_NSP
              PM_GOLD(6 + PM_I) = 0.999D0
            END DO
          END IF
          DO PM_I = 1, 20
            PM_TH(PM_I) = prop(7 + PM_I)
          END DO
          CM_TS = 1.0D0
          IF (nProp .GE. 31) THEN
            IF (prop(31) .GT. 0.0D0) CM_TS = prop(31)
          END IF
          CM_DT = CM_TS * dTime
          PM_NSUB = 1
          IF (CM_DT .GT. PM_DTMAX) PM_NSUB = CEILING(CM_DT / PM_DTMAX)
          SP_ID = SPLIT_ID(elemId, kDomIntPt)
          PM_OK = .TRUE.
          PM_PHIINT = 0.0D0
          IF (SP_SUM .LT. SP_MIN .OR. SP_SUM .LE. 0.0D0) THEN
            SP_HITI = 2
            DO PM_I = 1, 12
              PM_GNEW(PM_I) = PM_GOLD(PM_I)
            END DO
          ELSE
            CM_S = PM_GOLD(1)
            DO PM_I = 2, CM_NSP
              CM_S = CM_S + PM_GOLD(PM_I)
            END DO
            CM_F = SP_SUM / CM_S
            DO PM_I = 1, CM_NSP
              PM_GOLD(PM_I) = PM_GOLD(PM_I) * CM_F
            END DO
            PM_GOLD(6) = 1.0D0 - SP_SUM
            CALL SPLIT_CACHE_GET(SP_ID, ldstep, isubst, PM_GOLD,
     &                           PM_GNEW, PM_PHIINT, SP_HIT)
            IF (SP_HIT) THEN
              SP_HITI = 1
            ELSE
              SP_HITI = 0
              CM_CREF = 0.0D0
              IF (nProp .GE. 33) CM_CREF = prop(33)
              IF (CM_CREF .GT. 0.0D0 .AND. CM_NUT .GT. -1.0D29) THEN
                CM_CREL = MIN(MAX(CM_NUT / CM_CREF, 0.0D0), 1.0D0)
                CALL biofilm_ecology_hook_c(PM_GOLD, PM_TH, CM_DT,
     &               PM_NSUB, CM_CREL, PM_GNEW, PM_PHIINT, PM_OK)
              ELSE
                CALL biofilm_ecology_hook(PM_GOLD, PM_TH, CM_DT,
     &               PM_NSUB, PM_GNEW, PM_PHIINT, PM_OK)
              END IF
              IF (PM_OK) THEN
                DO PM_I = 1, 12
                  IF (PM_GNEW(PM_I) .NE. PM_GNEW(PM_I)) PM_OK = .FALSE.
                END DO
                CM_S = PM_GNEW(1)
                DO PM_I = 2, CM_NSP
                  CM_S = CM_S + PM_GNEW(PM_I)
                END DO
                IF (CM_S .LE. 0.0D0) PM_OK = .FALSE.
              END IF
              IF (PM_OK) THEN
                CM_S = PM_GNEW(1)
                DO PM_I = 2, CM_NSP
                  CM_S = CM_S + PM_GNEW(PM_I)
                END DO
                CM_F = SP_SUM / CM_S
                DO PM_I = 1, CM_NSP
                  PM_GNEW(PM_I) = PM_GNEW(PM_I) * CM_F
                END DO
                PM_GNEW(6) = 1.0D0 - SP_SUM
                CALL SPLIT_CACHE_PUT(SP_ID, ldstep, isubst, PM_GOLD,
     &                               PM_GNEW, PM_PHIINT)
              END IF
            END IF
          END IF
          IF (PM_OK) THEN
            DO PM_I = 1, 12
              ustatev(71 + PM_I) = PM_GNEW(PM_I)
            END DO
          ELSE
            keycut = 1
          END IF
          CALL BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI_USED, prop(7),
     &                                dTime, ALPHA_NEW)
C         prop(36) = d (not 0): species-weighted growth law, the
C         composition into the growth (5 Oct, two_way_growth_rate.py):
C         alpha_dot = k_alpha phi sum_i f_i chi_i, chi from the point
C         model's state after this step. Two species: f = 1 +- d. More
C         species (prop(37) = n > 2): f_i = prop(37 + i), i = 1..n.
C         Not from a paper; d = 0 (default) is Eq. 36 unchanged.
          IF (nProp .GE. 36) THEN
            IF (prop(36) .NE. 0.0D0 .AND. PM_OK) THEN
              CM_S = PM_GNEW(1)
              DO PM_I = 2, CM_NSP
                CM_S = CM_S + PM_GNEW(PM_I)
              END DO
              IF (CM_S .GT. 0.0D0) THEN
                IF (CM_NSP .EQ. 2) THEN
                  CM_GF = 1.0D0 + prop(36) * (2.0D0 * PM_GNEW(1)
     &                    / CM_S - 1.0D0)
                ELSE
                  CM_GF = 0.0D0
                  DO PM_I = 1, CM_NSP
                    IF (nProp .GE. 37 + PM_I) CM_GF = CM_GF
     &                  + prop(37 + PM_I) * PM_GNEW(PM_I) / CM_S
                  END DO
                END IF
                ALPHA_NEW = ALPHA_N + CM_GF * (ALPHA_NEW - ALPHA_N)
              END IF
            END IF
          END IF
          IF (kDomIntPt .EQ. 1 .AND.
     &        (MOD(elemId, TRACE_STRIDE) .EQ. 1 .OR.
     &         PHI_USED .GE. CM_TRACE_PHI_MIN)) THEN
            IF (.NOT. CM_TRACE_OPEN) THEN
              OPEN(96, FILE='comp_trace.csv', STATUS='UNKNOWN',
     &             POSITION='APPEND')
              CM_TRACE_OPEN = .TRUE.
            END IF
            WRITE(96,'(I0,5(",",I0),6(",",ES24.16E3),'//
     &            '24(",",ES24.16E3))')
     &        elemId, kDomIntPt, ldstep, isubst, PM_NSUB, SP_HITI,
     &        dTime, CM_DT, PHI_USED, SP_SUM, ALPHA_N, ALPHA_NEW,
     &        (PM_GOLD(PM_I), PM_I = 1, 12),
     &        (PM_GNEW(PM_I), PM_I = 1, 12)
            FLUSH(96)
          END IF
        ELSE IF (prop(28) .GT. 4.5D0) THEN
C         --- reaction step of the transport/reaction split ---------
C         The 3D field (bio1, bio2 at this Gauss point) gives the
C         species' volume fractions; psi and gamma are carried here.
C         One point-model step over dTime gives the reaction rate
C         R_i = (phi_i after - phi_i before) / dTime, handed to USSFin
C         (biofilm_split). Growth stays Klempt 2024 Eq. 36 with
C         phi = sum of the species: the PHI_USED call below.
          SP_N = 1
          IF (prop(28) .GT. 5.5D0) SP_N = 2
          SP_PHI(1) = MIN(MAX(Sdp_bio1_n, 0.0D0), 1.0D0)
          SP_PHI(2) = 0.0D0
          IF (SP_N .EQ. 2) SP_PHI(2) = MIN(MAX(Sdp_bio2_n, 0.0D0),
     &                                     1.0D0)
          SP_SUM = SP_PHI(1) + SP_PHI(2)
          IF (SP_SUM .GT. SP_PHIMAX) THEN
            SP_PHI(1) = SP_PHI(1) * (SP_PHIMAX / SP_SUM)
            SP_PHI(2) = SP_PHI(2) * (SP_PHIMAX / SP_SUM)
            SP_SUM = SP_PHI(1) + SP_PHI(2)
          END IF
          DO PM_I = 1, 12
            PM_GOLD(PM_I) = ustatev(71 + PM_I)
          END DO
C         first call at this point: psi = 0.999 for the active species
          IF (PM_GOLD(7) .EQ. 0.0D0) THEN
            DO PM_I = 7, 6 + SP_N
              PM_GOLD(PM_I) = 0.999D0
            END DO
            PM_GOLD(12) = 0.0D0
          END IF
          DO PM_I = 1, 5
            PM_GOLD(PM_I) = 0.0D0
          END DO
          DO PM_I = 7 + SP_N, 11
            PM_GOLD(PM_I) = 0.0D0
          END DO
          PM_GOLD(1) = SP_PHI(1)
          IF (SP_N .EQ. 2) PM_GOLD(2) = SP_PHI(2)
          PM_GOLD(6) = 1.0D0 - SP_SUM
C         below phi_min (prop(29)) no reaction: the point model makes
C         biomass out of nothing at phi = 0, so the void is reached by
C         transport only
          SP_MIN = 0.0D0
          IF (nProp .GE. 29) SP_MIN = prop(29)
          DO PM_I = 1, 20
            PM_TH(PM_I) = prop(7 + PM_I)
          END DO
          PM_NSUB = 1
          IF (dTime .GT. PM_DTMAX) PM_NSUB = CEILING(dTime / PM_DTMAX)
          SP_ID = SPLIT_ID(elemId, kDomIntPt)
          SP_HITI = 0
          PM_OK = .TRUE.
          IF (SP_SUM .LT. SP_MIN) THEN
            SP_HITI = 2
            DO PM_I = 1, 12
              PM_GNEW(PM_I) = PM_GOLD(PM_I)
            END DO
            PM_PHIINT = 0.0D0
          ELSE
            CALL SPLIT_CACHE_GET(SP_ID, ldstep, isubst, PM_GOLD,
     &                           PM_GNEW, PM_PHIINT, SP_HIT)
            IF (SP_HIT) THEN
              SP_HITI = 1
            ELSE
              CALL biofilm_ecology_hook(PM_GOLD, PM_TH, dTime,
     &             PM_NSUB, PM_GNEW, PM_PHIINT, PM_OK)
              IF (PM_OK) THEN
                DO PM_I = 1, 12
                  IF (PM_GNEW(PM_I) .NE. PM_GNEW(PM_I)) PM_OK = .FALSE.
                END DO
                IF (ABS(PM_GNEW(12)) .GT. 1.0D5) PM_OK = .FALSE.
              END IF
              IF (PM_OK) CALL SPLIT_CACHE_PUT(SP_ID, ldstep, isubst,
     &                        PM_GOLD, PM_GNEW, PM_PHIINT)
            END IF
          END IF
          IF (PM_OK) THEN
            SP_R(1) = (PM_GNEW(1) - PM_GOLD(1)) / dTime
            SP_R(2) = 0.0D0
            IF (SP_N .EQ. 2) SP_R(2) = (PM_GNEW(2) - PM_GOLD(2)) / dTime
            DO PM_I = 1, 12
              ustatev(71 + PM_I) = PM_GNEW(PM_I)
            END DO
          ELSE
            SP_R(1) = 0.0D0
            SP_R(2) = 0.0D0
            keycut = 1
          END IF
          CALL SPLIT_PUT(SP_ID, SP_R(1), SP_R(2), ldstep, isubst)
          CALL BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI_USED, prop(7),
     &                                dTime, ALPHA_NEW)
          IF (kDomIntPt .EQ. 1 .AND.
     &        (MOD(elemId, TRACE_STRIDE) .EQ. 1 .OR.
     &         PHI_USED .GE. TRACE_PHI_MIN)) THEN
            IF (.NOT. SP_TRACE_OPEN) THEN
              OPEN(99, FILE='split_trace.csv', STATUS='UNKNOWN',
     &             POSITION='APPEND')
              SP_TRACE_OPEN = .TRUE.
            END IF
            WRITE(99,'(I0,5(",",I0),7(",",ES24.16E3),'//
     &            '24(",",ES24.16E3))')
     &        elemId, kDomIntPt, ldstep, isubst, PM_NSUB, SP_HITI,
     &        dTime, Sdp_bio1_n, Sdp_bio2_n, SP_R(1), SP_R(2),
     &        ALPHA_N, ALPHA_NEW, (PM_GOLD(PM_I), PM_I = 1, 12),
     &        (PM_GNEW(PM_I), PM_I = 1, 12)
            FLUSH(99)
          END IF
        ELSE IF (prop(28) .GT. 3.5D0) THEN
C         --- point model, inner steps inside this one call ---
          PM_SUM = 0.0D0
          DO PM_I = 1, 12
            PM_GOLD(PM_I) = ustatev(71 + PM_I)
            PM_SUM = PM_SUM + ABS(PM_GOLD(PM_I))
          END DO
          IF (PM_SUM .EQ. 0.0D0) THEN
C           first call at this point: the five-species default seed
C           (same numbers as INIT_ECO_IF_ZERO); the server masks it
C           down to --active-species and restores phi0 itself.
            PM_GOLD(1) = 0.12D0
            PM_GOLD(2) = 0.12D0
            PM_GOLD(3) = 0.08D0
            PM_GOLD(4) = 0.05D0
            PM_GOLD(5) = 0.0D0
            PM_SUM = 0.37D0
            DO PM_I = 1, 5
              PM_GOLD(PM_I) = PM_GOLD(PM_I) * (0.999999D0 / PM_SUM)
            END DO
            PM_GOLD(6) = 1.0D0 - (PM_GOLD(1) + PM_GOLD(2) + PM_GOLD(3)
     &                   + PM_GOLD(4) + PM_GOLD(5))
            DO PM_I = 7, 11
              PM_GOLD(PM_I) = 0.999D0
            END DO
            PM_GOLD(12) = 0.0D0
          END IF
          DO PM_I = 1, 20
            PM_TH(PM_I) = prop(7 + PM_I)
          END DO
          PM_NSUB = 1
          IF (dTime .GT. PM_DTMAX) PM_NSUB = CEILING(dTime / PM_DTMAX)
          CALL biofilm_ecology_hook(PM_GOLD, PM_TH, dTime, PM_NSUB,
     &                              PM_GNEW, PM_PHIINT, PM_OK)
C         the 9/7 sanity check, kept: finite, phi_tot plausible,
C         gamma well inside the solver's failure clamp
          IF (PM_OK) THEN
            PM_PHITOT = 0.0D0
            DO PM_I = 1, 5
              PM_PHITOT = PM_PHITOT + PM_GNEW(PM_I) * PM_GNEW(5 + PM_I)
            END DO
            DO PM_I = 1, 12
              IF (PM_GNEW(PM_I) .NE. PM_GNEW(PM_I)) PM_OK = .FALSE.
            END DO
            IF (PM_PHITOT .GT. 1.5D0) PM_OK = .FALSE.
            IF (ABS(PM_GNEW(12)) .GT. 1.0D5) PM_OK = .FALSE.
          END IF
          IF (PM_OK) THEN
            DO PM_I = 1, 12
              ustatev(71 + PM_I) = PM_GNEW(PM_I)
            END DO
            ALPHA_NEW = ALPHA_N + prop(7) * PM_PHIINT
          ELSE
C           hold the converged state and ask ANSYS to cut the step
            ALPHA_NEW = ALPHA_N
            keycut = 1
          END IF
          IF (kDomIntPt .EQ. 1 .AND.
     &        (MOD(elemId, TRACE_STRIDE) .EQ. 1 .OR.
     &         PHI_USED .GE. TRACE_PHI_MIN)) THEN
            IF (.NOT. PM_TRACE_OPEN) THEN
              OPEN(98, FILE='pm_trace.csv', STATUS='UNKNOWN',
     &             POSITION='APPEND')
              PM_TRACE_OPEN = .TRUE.
            END IF
            WRITE(98,'(I0,4(",",I0),3(",",ES24.16E3),'//
     &            '24(",",ES24.16E3))')
     &        elemId, kDomIntPt, ldstep, isubst, PM_NSUB, dTime,
     &        ALPHA_N, ALPHA_NEW, (PM_GOLD(PM_I), PM_I = 1, 12),
     &        (PM_GNEW(PM_I), PM_I = 1, 12)
            FLUSH(98)
          END IF
        ELSE IF (prop(28) .GT. 2.5D0) THEN
          ALPHA_NEW = Sdp_sumLocal - 1.0D0
        ELSE
          CALL BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI_USED, prop(7),
     &                                dTime, ALPHA_NEW)
        END IF
C       prop(49) = p_ref > 0: homeostatic-pressure growth law (Keio WP2,
C       keio_wp2/README.ja.md; an assumption, not from a paper). The
C       increment of alpha above is multiplied by g = max(0, 1 - p/p_h),
C       p = -(s11 + s22 + s33)/3 of the stress the host hands in, i.e. at
C       the start of the increment (compression positive). p_h = prop(49)
C       (= P_h E k_alpha T*, in stress units), or with prop(50) = 1
C       p_h = prop(49) (phi^2 + prop(51)) (local stiffness). 0 = off
C       (default): nothing changes. Not in mode 3 (the partner's alpha).
C       ustatev(95) = p, ustatev(96) = g.
        IF (nProp .GE. 49) THEN
        IF (prop(49) .GT. 0.0D0 .AND. (prop(28) .LT. 2.5D0 .OR.
     &      prop(28) .GT. 3.5D0)) THEN
          HP_P = -(stress(1) + stress(2) + stress(3)) / 3.0D0
          HP_FORM = 0
          HP_F = 0.0D0
          IF (nProp .GE. 50) HP_FORM = NINT(prop(50))
          IF (nProp .GE. 51) HP_F = prop(51)
          CALL BIOFILM_HOMEOSTATIC_FACTOR(HP_P, prop(49), HP_FORM,
     &                                    PHI_USED, HP_F, HP_G)
          ALPHA_NEW = ALPHA_N + HP_G * (ALPHA_NEW - ALPHA_N)
          ustatev(95) = HP_P
          ustatev(96) = HP_G
        END IF
        END IF
        Sbio_GrowthConst = ALPHA_NEW
        ustatev(84)      = ALPHA_NEW
C       trace: integration point 1 of one element in TRACE_STRIDE,
C       and of every element where phi has reached TRACE_PHI_MIN
        IF (kDomIntPt .EQ. 1 .AND.
     &      (MOD(elemId, TRACE_STRIDE) .EQ. 1 .OR.
     &       PHI_USED .GE. TRACE_PHI_MIN)) THEN
          IF (.NOT. TRACE_OPEN) THEN
            OPEN(97, FILE='phi_trace.csv', STATUS='UNKNOWN',
     &           POSITION='APPEND')
            TRACE_OPEN = .TRUE.
          END IF
          WRITE(97,'(I0,3(",",I0),8(",",ES24.16E3))')
     &      elemId, kDomIntPt, ldstep, isubst, dTime,
     &      Sdp_bio1_n, Sdp_locbio1_n, Sdp_bio2_n, Sdp_locbio2_n,
     &      PHI_USED, ALPHA_N, ALPHA_NEW
          FLUSH(97)
        END IF
      END IF
      END IF
