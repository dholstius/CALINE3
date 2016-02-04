      SUBROUTINE CAL3CONC(U,W,W2,H,D,Z,UWL,DWL,
     +                    MIXH,BASE,PHI,Q1,
     +                    PY1,PY2,PZ1,PZ2,
     +                    DSTR,
     +                    CONC)

      USE CAL3DATA
      IMPLICIT NONE

      REAL CONC,INC

C     Arguments passed to this subroutine
      DOUBLE PRECISION D
      REAL U,W,W2,H,Z,UWL,DWL
      REAL MIXH,BASE,PHI,Q1
      REAL PY1,PY2,PZ1,PZ2,DSTR

C     Local integers (loop counters)
      INTEGER I

C     Local scalars
      DOUBLE PRECISION FAC2,HYP,SIDE,PD,FET
      REAL FAC1,FAC3,FAC4,FAC5,FACT
      REAL ARG1,ARG2,ARG3
      REAL EXP1,EXP2
      REAL ED1,ED2,EL2,ELL2,EM2,EN2,ECLD,EFRC
      REAL EXLS
      REAL CNT
      REAL STP
      REAL CSL2
      REAL NE
      REAL YE
      REAL QE
      REAL SGZ,SGY,KZ
      REAL SGN,FINI
      REAL Y
      REAL T
      REAL VS,V1
      REAL HDS

C     Local arrays
      DOUBLE PRECISION INTG
      DIMENSION INTG(6),Y(6)

#ifdef debug
      PRINT *, "U,W,W2,H,D,Z,UWL,DWL are: ", U,W,W2,H,D,Z,UWL,DWL
      PRINT *, "MIXH,BASE,PHI are: ", MIXH,BASE,PHI,Q1
      PRINT *, "PY1,PY2,PZ1,PZ2 are: ", PY1,PY2,PZ1,PZ2
      PRINT *, "DSTR is: ", DSTR
#endif

 3050 SGN=1.

#ifdef debug
      PRINT *, "CONC is: ", CONC
#endif

 3060 NE=0.
      STP=1.
      FINI=1.
      IF (SGN.EQ.1. .AND. UWL.LE.0. .AND. DWL.LT.0.) SGN=-1.
 3080 IF (SGN.EQ.-1. .AND. UWL.GT.0. .AND. DWL.GE.0.) GO TO 6000

C *****  ELEMENT LOOP  *****
      ED1=0.
      ED2=SGN*W

 3110 IF (SGN.EQ.-1.) GO TO 3160
      IF (ED1.LE.DWL .AND. ED2.LE.DWL) GO TO 3770
      IF (ED1.GT.DWL .AND. ED2.LT.UWL) GO TO 3250
      IF (ED1.LE.DWL) ED1=DWL
      IF (ED2.LT.UWL) GO TO 3250
      ED2=UWL
      SGN=-1.
      NE=-1.
      GO TO 3250

 3160 IF (ED1.GE.UWL .AND. ED2.GE.UWL) GO TO 3770
      IF (ED1.LT.UWL .AND. ED2.GT.DWL) GO TO 3250
      IF (ED1.GE.UWL) ED1=UWL
      IF (ED2.GT.DWL) GO TO 3250
      ED2=DWL
      FINI=0.

 3250 EL2=ABS(ED2-ED1)/2.
      ECLD=(ED1+ED2)/2.
      ELL2=W2/COS(PHI)+(EL2-W2*TAN(PHI))*SIN(PHI)
      IF (PHI.GE.ATAN(W2/EL2)) THEN
        CSL2=W2/SIN(PHI)
      ELSE
        CSL2=EL2/COS(PHI)
      END IF
      EM2=ABS((EL2-W2/TAN(PHI))*SIN(PHI))
      EN2=(ELL2-EM2)/2.

C *****  RECEPTOR DISTANCE LOOP  *****
      QE=Q1*CSL2/W2
      FET=(ECLD+D*TAN(PHI))*COS(PHI)
      HYP=ECLD**2+D**2
      SIDE=FET**2
      IF (SIDE.GT.HYP) YE=0.
      IF (SIDE.LE.HYP) YE=SQRT(HYP-SIDE)

#ifdef debug
      PRINT *, ""
      PRINT *, "QE,FET,HYP,SIDE are: ", Q1,PY1,PY2,DSTR
#endif

C *****  DETERMINE SIGMA Y AND SIGMA Z  *****
      IF (FET.LE.-CSL2) GO TO 3830


      IF (FET.GE.CSL2) GO TO 3320

      QE=QE*(FET+CSL2)/(2.*CSL2)
      FET=(CSL2+FET)/2.

 3320 CALL CAL3SGMA(FET,PZ1,PZ2,SGZ)
      CALL CAL3SGMA(FET,PY1,PY2,SGY)

#ifdef debug
      PRINT *, "SGZ,SGY are: ", SGZ,SGY
#endif

      FAC1=0.399/(SGZ*U)

C *****  ADJUSTMENT FOR ELEMENT END EFFECT  *****
C           (POLYNOMIAL APPROXIMATION)
      Y(1)=YE+ELL2
      Y(2)=Y(1)-EN2
      Y(3)=Y(2)-EN2
      Y(4)=Y(3)-2*EM2
      Y(5)=Y(4)-EN2
      Y(6)=Y(5)-EN2

C ***  SUB-ELEMENT SOURCE STRENGTH LOOP
      DO 3480 I=1,6
        CALL CAL3INTG(Y(I),SGY,INTG(I))
 3480 CONTINUE

#ifdef debug
      PRINT *, ""
      PRINT *, "INTG are: ", INTG
      PRINT *, "Y are: ", Y
#endif

      FAC2=0.
      DO 3530 I=1,5
        IF ((SIGN(1.,Y(I))).EQ.(SIGN(1.,Y(I+1)))) THEN
          PD=ABS(INTG(I+1)-INTG(I))
        ELSE
          PD=1.-INTG(I)-INTG(I+1)
        END IF
        FAC2=FAC2+PD*QE*WT(I)
 3530 CONTINUE

      FACT=FAC1*FAC2

#ifdef debug
      PRINT *, "FAC1,FAC2,FACT are: ", FAC1,FAC2,FACT
#endif

C *****  DEPRESSED SECTION  *****

      IF (HDS.LT.-1.5 .AND. ABS(D).LT.(W2-3.*HDS)) THEN
 3560   IF (ABS(D).LE.W2) THEN
          FACT=FACT*DSTR
        ELSE
          FACT=FACT*(DSTR-(DSTR-1.)*(ABS(D)-W2)/(-3.*HDS))
        END IF
      END IF
C
C *****  DEPOSITION CORRECTION  *****
C

 3580 FAC3=0.

      IF (V1.NE.0.) THEN
        KZ=SGZ**2*U/(2.*FET)
        ARG3=V1*SGZ/(KZ*SQRT(2.))+(Z+H)/(SGZ*SQRT(2.))
        IF (ARG3.GT.5.) GO TO 3770
        T=1./(1.+0.47047*ARG3)
        EFRC=(.3480242*T-.0958798*T**2+.7478556*T**3)*EXP(-1.*ARG3**2)
        FAC3=(SQRT(2.*PI)*V1*SGZ*EXP(V1*(Z+H)/KZ+.5*(V1*SGZ/KZ)**2)
     *    *EFRC)/KZ
        IF (FAC3.GT.2.) FAC3=2.
      END IF

C *****  SETTLING CORRECTION  *****
 3670 IF (VS.NE.0.) THEN
        FAC4=EXP(-VS*(Z-H)/(2.*KZ)-(VS*SGZ/KZ)**2/8.)
        FACT=FACT*FAC4
      END IF


C *****  INCREMENTAL CONCENTRATION  *****
 3710 FAC5=0.
      CNT=0.
 3720 EXLS=0.
 3730 ARG1=-0.5*((Z+H+2.*CNT*MIXH)/SGZ)**2
      IF (ARG1.LT.-44.) THEN
        EXP1=0.
      ELSE
        EXP1=EXP(ARG1)
      END IF
      ARG2=-0.5*((Z-H+2.*CNT*MIXH)/SGZ)**2
      IF (ARG2.LT.-44.) THEN
        EXP2=0.
      ELSE
        EXP2=EXP(ARG2)
      END IF
      FAC5=FAC5+EXP1+EXP2
      IF (MIXH.GE.1000.) GO TO 3760
      IF ((EXP1+EXP2+EXLS).EQ.0. .AND. CNT.LE.0.) GO TO 3760
 3740 IF (CNT.GT.0.) GO TO 3750
      CNT=ABS(CNT)+1.
      GO TO 3720
 3750 CNT=-1.*CNT
      EXLS=EXP1+EXP2
      GO TO 3730

 3760 INC=FACT*(FAC5-FAC3)
      CONC=CONC+INC

#ifdef debug
      PRINT *, "FACT,FAC5,FAC3 are: ", FACT,FAC5,FAC3
      PRINT *, "INC is: ", INC
      PRINT *, "CONC is: ", CONC
#endif

 3770 IF (FINI.EQ.0.) GO TO 6000
      NE=NE+1.
      STP=BASE**NE
      IF (NE.EQ.0.) GO TO 3080
      ED1=ED2
      ED2=ED2+SGN*STP*W
      GO TO 3110

 3830 IF (SGN.EQ.1.) GO TO 3770

 6000 RETURN

      END SUBROUTINE
