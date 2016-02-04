      SUBROUTINE CAL3INTG(Y,SGY,INTG)

        IMPLICIT NONE

        REAL Y,SGY
        REAL LIM,T,ARG0
        REAL X0,X1,X2,X3,X4,X5
        DOUBLE PRECISION INTG

        PARAMETER(X0=0.3989,X1=0.3194,X2=0.3566)
        PARAMETER(X3=1.7815,X4=1.8213,X5=1.3303)

        LIM=ABS(Y/SGY)

        IF (LIM.GT.5.) THEN
          INTG=0.
        ELSE
          ARG0=LIM**2/(-2.)
          T=1./(1.+0.23164*LIM)
          INTG=X0*EXP(ARG0)*(X1*T-X2*T**2+X3*T**3-X4*T**4+X5*T**5)
        END IF

        RETURN

      END SUBROUTINE
