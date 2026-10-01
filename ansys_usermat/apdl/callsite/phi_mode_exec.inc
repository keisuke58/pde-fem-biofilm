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
        IF (prop(28) .GT. 1.5D0 .AND. prop(28) .LT. 2.5D0) THEN
          PHI_USED = Sdp_bio1_n + Sdp_bio2_n
        ELSE
          PHI_USED = Sdp_bio1_n
        END IF
        ALPHA_N = ustatev(84)
        IF (prop(28) .GT. 3.5D0) THEN
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
