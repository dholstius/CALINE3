      SUBROUTINE CAL3RXL(
     +                   NR,      ! number of receptors
     +                   XR,      ! x-coordinates of receptors
     +                   YR,      ! y-coordinates of receptors
     +                   ZR,      ! z-coordinates (heights) of receptors
     +                   NL,      ! number of links
     +                   XL1,     ! x-coordinates of link vertices
     +                   YL1,     ! y-coordinates of link vertices
     +                   XL2,     ! x-coordinates of link vertices
     +                   YL2,     ! y-coordinates of link vertices
     +                   WL,      ! link widths
     +                   HL,      ! link heights
     +                   NTYP,    ! link classifications *as integers*
     +                   VPHL,    ! traffic volume (per link)
     +                   EFL,     ! emission factor (per link)
     +                   U,       ! wind speed
     +                   BRG,     ! wind bearing
     +                   CLAS,    ! atmospheric stability class
     +                   MIXH,    ! mixing height
     +                   ATIM,    ! averaging time
     +                   Z0,      ! surface roughness
     +                   VS,      ! settling velocity
     +                   VD,      ! deposition velocity
     +                   LXR,
     +                   C)       ! resulting concentrations (per receptor)

      USE CAL3DATA

      INTEGER NR,NL,NM
      DIMENSION C(NR,NL)

      LOGICAL LXR
      DIMENSION LXR(NR,NL)

      DOUBLE PRECISION HYP,SIDE,FAC2,PD,A,B,L,D,
     +    XPRI,YPRI,APRI,BPRI,LPRI,DPRI,XD,YD,D1,D2,
     +    LL(NL),INTG(6)

      INTEGER CLAS

      REAL CONC
      REAL NE,LIM,KZ,LB,INC
      REAL V1,YE,Z,EXP1,EXP2,DVIR,CSL2

      REAL XR(NR),YR(NR),ZR(NR)
      REAL XL1(NL),YL1(NL),XL2(NL),YL2(NL),WL(NL),HL(NL)
      REAL VPHL(NL),EFL(NL)
      REAL U,BRG,MIXH

      INTEGER NTYP(NL)
      PARAMETER(NTYP_AG=1,NTYP_BR=2,NTYP_FL=3,NTYP_DP=4)

      V1=VD-VS/2.

      DO 1050 I=1,NL
        LL(I)=SQRT((XL1(I)-XL2(I))**2+(YL1(I)-YL2(I))**2)
 1050 CONTINUE

      BRG=BRG+180.
      IF (BRG.GE.360.) BRG=BRG-360.

      AFAC=(ATIM/3.0)**.2
      IF (CLAS.GT.6) CLAS=6
      SY1=ALOG(AY1(CLAS)*((Z0/3.)**.2)*AFAC)
      SY10=ALOG(AY2(CLAS)*((Z0/3.)**.07)*AFAC)
      PY1=EXP(SY1)
      PY2=(SY10-SY1)/DREF
      SZ10=ALOG(AZ(CLAS)*((Z0/10.)**.07)*AFAC)

C *****  LINK LOOP  *****
      DO 8000 IL=1,NL

      IF (NTYP(IL).EQ.NTYP_DP .OR. NTYP(IL).EQ.NTYP_FL) THEN
        H=0.
      ELSE
        H=HL(IL)
      END IF
      W=WL(IL)
      W2=W/2.

      Q1=0.1726*VPHL(IL)*EFL(IL)

      XD=XL2(IL)-XL1(IL)
      YD=YL2(IL)-YL1(IL)
C      ABS(XD)=ABS(XD)
      IF(ABS(XD).GT.LL(IL)) LL(IL)=ABS(XD)
      LB=DEG*(ACOS(ABS(XD)/LL(IL)))
      IF (XD.GT.0. .AND. YD.GE.0.) LB=90.-LB
      IF (XD.GE.0. .AND. YD.LT.0.) LB=90.+LB
      IF (XD.LT.0. .AND. YD.LE.0.) LB=270.-LB
      IF (XD.LE.0. .AND. YD.GT.0.) LB=270.+LB

      PHI=ABS(BRG-LB)
      IF (PHI.LE.90.) GO TO 7600
      IF (PHI.GE.270.) GO TO 5000
      PHI=ABS(PHI-180.)
      GO TO 7600
 5000 PHI=ABS(PHI-360.)
 7600 IF (PHI.LT.20.) GO TO 7630
      IF (PHI.LT.50.) GO TO 7620
      IF (PHI.LT.70.) GO TO 7610
      BASE=4.
      GO TO 7650
 7610 BASE=2.
      GO TO 7650
 7620 BASE=1.5
      GO TO 7650
 7630 BASE=1.1

 7650 PHI=RAD*(PHI)
      IF (PHI.GT.1.5706) PHI=1.5706
      IF (PHI.LT.0.00017) PHI=0.00017

C *****  DEPRESSED SECTION  *****
      IF (HL(IL).LT.-1.5) GO TO 7700
      DSTR=1.
      HDS=1.
      GO TO 7800
 7700 HDS=HL(IL)
      DSTR=0.72*ABS(HDS)**0.83
C     ! RESIDENCE TIME FACTOR

C *****  SIGMA Z POWER CURVE  *****
 7800 TR=DSTR*W2/U
C     ! RESIDENCE TIME
      SGZ1=ALOG((1.8+0.11*TR)*(ATIM/30.)**0.2)
C     ! ALOG(SIGMA Z) AT W2
      PZ2=(SZ10-SGZ1)/(DREF-ALOG(W2))
      PZ1=EXP((SZ10+SGZ1-PZ2*(DREF+ALOG(W2)))/2.)

C *****  END OF LINK ROUTINE  *****

C *****  RECEPTOR LOOP  *****
      DO 6000 IR=1,NR

      IF (.NOT.(LXR(IR,IL))) GO TO 6000

      A=(XR(IR)-XL1(IL))**2+(YR(IR)-YL1(IL))**2
      B=(XR(IR)-XL2(IL))**2+(YR(IR)-YL2(IL))**2

      L=(B-A-LL(IL)**2)/(2.*LL(IL))
C     ! OFFSET LENGTH
      IF (A.GT.L**2) D=DSQRT(A-L**2)
      IF (A.LE.L**2) D=0.
C     ! RECEPTOR DISTANCE
      UWL=LL(IL)+L
C     ! UPWIND LENGTH
      DWL=L
C     ! DOWNWIND LENGTH
      IF(D.EQ.0.D0)DVIR=1.D0
      IF(D.NE.0.D0)DVIR=D
      XPRI=XR(IR)+DVIR*COS(RAD*(450.-BRG))
      YPRI=YR(IR)+DVIR*SIN(RAD*(450.-BRG))
      APRI=(XPRI-XL1(IL))**2+(YPRI-YL1(IL))**2
      BPRI=(XPRI-XL2(IL))**2+(YPRI-YL2(IL))**2
      LPRI=(BPRI-APRI-LL(IL)**2)/(2.*LL(IL))
      IF (APRI.GT.LPRI**2) DPRI=DSQRT(APRI-LPRI**2)
      IF (APRI.LE.LPRI**2) DPRI=0.
      IF (DPRI.LT.D) D=-D

      IF (LPRI-L) 5725,5735,5735
 5725 TEMP=UWL
      UWL=-DWL
      DWL=-TEMP
 5735 IF (NTYP(IL).EQ.NTYP_AG .OR. NTYP(IL).EQ.NTYP_BR) GO TO 5750
C
      D1=W2+2.*ABS(HL(IL))
      D2=W2
C     ! SINGLE PRECISION TO DOUBLE PRECISION FOR LOGICAL 'IF'
      IF (DABS(D).GE.D1) GO TO 5750
C     ! 2:1 SLOPE ASSUMED
      IF (DABS(D).LE.D2) THEN
        Z=ZR(IR)-HL(IL)
      ELSE
        Z=ZR(IR)-HL(IL)*(1.-(DABS(D)-W2)/(2.*ABS(HL(IL))))
      END IF
      GO TO 3059
 5750 Z=ZR(IR)

 3059 CONC = 0.
      CALL CAL3CONC(U,W,W2,H,D,Z,UWL,DWL,
     +              MIXH,BASE,PHI,Q1,
     +              PY1,PY2,PZ1,PZ2,DSTR,
     +              CONC)
      C(IR,IL)=C(IR,IL)+CONC

 6000 CONTINUE

 8000 CONTINUE

 9000 END SUBROUTINE


