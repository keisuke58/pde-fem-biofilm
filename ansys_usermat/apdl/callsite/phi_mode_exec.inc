C     --- phi mode, executable part (pde-fem-biofilm, 2026-10-01) -----
C     Paste inside the prop(1) /= 0 branch, AFTER sGrowth has been set
C     from prop(5) or the ecology mode, BEFORE the call to
C     BIOFILM_GROWTH_VISCO_V01.
C
C     Klempt et al. 2024, Eq. 36: alpha_dot = k_alpha * phi, nothing
C     added. prop(6) = 2: one species, phi = bio1.
C              prop(6) = 3: two species, phi = bio1 + bio2.
C     prop(7) = k_alpha. ustatev(84) holds the CONVERGED alpha: ANSYS
C     hands it back unchanged on every equilibrium iteration, so growth
C     is counted once per increment.
C     If stage 1 shows locbio is the Gauss-point value, change bio to
C     locbio in the two PHI_USED lines -- nothing else.
      IF (prop(6) .GT. 1.5D0) THEN
        IF (prop(6) .GT. 2.5D0) THEN
          PHI_USED = Sdp_bio1_n + Sdp_bio2_n
        ELSE
          PHI_USED = Sdp_bio1_n
        END IF
        ALPHA_N = ustatev(84)
        CALL BIOFILM_ALPHA_FROM_PHI(ALPHA_N, PHI_USED, prop(7),
     &                              dTime, ALPHA_NEW)
        sGrowth     = ALPHA_NEW
        ustatev(84) = ALPHA_NEW
C       trace: integration point 1 of about one element in 997
        IF (kDomIntPt .EQ. 1 .AND. MOD(elemId, 997) .EQ. 1) THEN
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
