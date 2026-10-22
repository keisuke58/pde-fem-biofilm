      LOGICAL FUNCTION BIOFILM_TRACE_WANT(IEL)
C     comp_trace.csv / phi_trace.csv: every element (default), or only
C     the elements listed in the environment variable
C     BIOFILM_TRACE_ELEMS (ids separated by commas). A-53: the traces
C     of a 16^3 run reached 11 GB; fields_living.py reads only two.
      INTEGER IEL, NEL, EL(16), I, J, K, ST
      CHARACTER*256 S
      LOGICAL INIT
      SAVE INIT, NEL, EL
      DATA INIT /.FALSE./
      IF (.NOT. INIT) THEN
        NEL = 0
        S = ' '
        CALL GET_ENVIRONMENT_VARIABLE('BIOFILM_TRACE_ELEMS', S,
     &       STATUS=ST)
        IF (ST .EQ. 0 .AND. LEN_TRIM(S) .GT. 0) THEN
          I = 1
          DO WHILE (I .LE. LEN_TRIM(S) .AND. NEL .LT. 16)
            J = INDEX(S(I:), ',')
            IF (J .EQ. 0) THEN
              K = LEN_TRIM(S)
            ELSE
              K = I + J - 2
            END IF
            IF (K .GE. I) THEN
              NEL = NEL + 1
              READ(S(I:K), *) EL(NEL)
            END IF
            I = K + 2
          END DO
        END IF
        INIT = .TRUE.
      END IF
      BIOFILM_TRACE_WANT = (NEL .EQ. 0)
      DO I = 1, NEL
        IF (IEL .EQ. EL(I)) BIOFILM_TRACE_WANT = .TRUE.
      END DO
      END
