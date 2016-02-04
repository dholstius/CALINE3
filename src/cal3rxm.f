      SUBROUTINE CAL3RXM(
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
     +                   NM,      ! number of met conditions
     +                   UM,      ! wind speeds
     +                   BRGM,    ! wind bearings
     +                   CLASM,   ! atmospheric stability classes
     +                   MIXHM,   ! mixing heights
     +                   ATIM,    ! averaging time
     +                   Z0,      ! surface roughness
     +                   VS,      ! settling velocity
     +                   VD,      ! deposition velocity
     +                   LXR,
     +                   C)       ! resulting concentrations (per receptor)
C
C    Computes incremental concentrations at a number of receptors,
C    given an array of link geometries + traffic volumes, an array of
C    meteorological conditions, and parameters describing the site.
C
C    The returned value will be an array of dimension NM x NR,
C    where NM is the number of meteorological conditions and NR is the
C    number of receptors.
C
C    Peak concentrations and average concentrations can be computed as
C    margined values from the returned array.
C

      USE CAL3DATA

      INTEGER NR,NL,NM
      DIMENSION C(NR,NM)

      LOGICAL LXR
      DIMENSION LXR(NR,NL)

      DOUBLE PRECISION HYP,SIDE,FAC2,PD,A,B,L,D,
     +    XPRI,YPRI,APRI,BPRI,LPRI,DPRI,XD,YD,D1,D2,
     +    LL(NL),INTG(6)

      INTEGER CLAS

      REAL CONC
      REAL U,W,W2,H,Z,UWL,DWL
      REAL MIXH,BASE,PHI,Q1
      REAL PY1,PY2,PZ1,PZ2,DSTR

      REAL MOWT,NE,LIM,KZ,LB,INC
      REAL V1,YE,EXP1,EXP2,DVIR,CSL2

      REAL XR(NR),YR(NR),ZR(NR)
      REAL XL1(NL),YL1(NL),XL2(NL),YL2(NL),WL(NL),HL(NL)
      REAL VPHL(NL),EFL(NL)
      REAL UM(NM),BRGM(NM),MIXHM(NM)
      INTEGER CLASM(NM)

      INTEGER NTYP(NL)
      PARAMETER(NTYP_AG=0,NTYP_BR=1,NTYP_FL=2,NTYP_DP=3)

C
C
C *****  INITIALIZATION OF CONSTANTS AND COUNTERS  *****
C
      MOWT=28.
C     ! MOLECULAR WEIGHT OF CO
      FPPM=0.0245/MOWT

      V1=VD-VS/2.

C     MET LOOP BEGINS
      DO 8500 IM=1,NM

      U=UM(IM)
      BRG=BRGM(IM)
      CLAS=CLASM(IM)
      MIXH=MIXHM(IM)
C      PRINT *, 'U, BRG, CLAS, MIXH are: ', U,BRG,CLAS,MIXH

      DO 1050 I=1,NL
      LL(I)=SQRT((XL1(I)-XL2(I))**2+(YL1(I)-YL2(I))**2)
C     ! LINK LENGTH
 1050 CONTINUE

C        U = WIND SPEED (M/S)
C      BRG = WIND DIRECTION (DEGREES)
C     CLAS = STABILITY CLASS (A-F)
C     MIXH = MIXING HEIGHT (M)
C      AMB = AMBIENT CONCENTRATION (PPM)
C
      BRG1=BRG
C     ! WIND ANGLE FOR OUTPUT
C
      BRG=BRG+180.
      IF (BRG.GE.360.) BRG=BRG-360.
C     ! CONVERSION TO VECTOR ORIENTATION
C
C ***  VIRTUAL DISPLACEMENT VECTORS
C
      XVEC=COS(RAD*(450.-BRG))
      YVEC=SIN(RAD*(450.-BRG))
C
C *****  CORRECTIONS FOR AVERAGING TIME AND SURFACE ROUGHNESS
C
      IF (CLAS.GT.6) CLAS=6

      AFAC=(ATIM/3.0)**.2
      SY1=ALOG(AY1(CLAS)*((Z0/3.)**.2)*AFAC)
C     ! ALOG(SIGMA Y) AT 1 M
      SY10=ALOG(AY2(CLAS)*((Z0/3.)**.07)*AFAC)
C     ! ALOG(SIGMA Y) AT 10 KM
      PY1=EXP(SY1)
      PY2=(SY10-SY1)/DREF
      SZ10=ALOG(AZ(CLAS)*((Z0/10.)**.07)*AFAC)
C     ! ALOG(SIGMA Z) AT 10 KM
C
C *****  LINK LOOP  *****
C
      !PRINT *, 'CLAS is: ', CLAS
      !PRINT *, 'AY1(CLAS), AY2(CLAS) are: ', AY1(CLAS), AY2(CLAS)
      !PRINT *, 'AFAC, ATIM are: ', AFAC, ATIM
      !PRINT *, 'AFAC, ATIM are: ', AFAC, ATIM
      !PRINT *, 'SY1, SY10 are: ', SY1, SY10
      !PRINT *, 'PY1, PY2 are: ', PY1, PY2

      DO 8000 IL=1,NL

      VPH=VPHL(IL)
      EF=EFL(IL)
      IF (NTYP(IL).EQ.NTYP_DP .OR. NTYP(IL).EQ.NTYP_FL) GO TO 870
      H=HL(IL)
      GO TO 880
  870 H=0.
  880 W=WL(IL)
C
C *****  LINK ROUTINE  *****
C **************************
C
      W2=W/2.
      Q1=0.1726*VPH*EF
C     ! LINEAL SOURCE STRENGTH PARALLEL TO HIGHWAY IN MICRO-GRAMS/
C                                                   (METER*SEC)
      XD=XL2(IL)-XL1(IL)
      YD=YL2(IL)-YL1(IL)
      IF(ABS(XD).GT.LL(IL)) LL(IL)=ABS(XD)

      LB=DEG*(ACOS(ABS(XD)/LL(IL)))
C     ! LINK BEARING

      IF (XD.GT.0. .AND. YD.GE.0.) LB=90.-LB
      IF (XD.GE.0. .AND. YD.LT.0.) LB=90.+LB
      IF (XD.LT.0. .AND. YD.LE.0.) LB=270.-LB
      IF (XD.LE.0. .AND. YD.GT.0.) LB=270.+LB

C     LBRG(IL)=LB

      !PRINT *, 'LL(IL), LB are:', LL(IL), LB

C     ! LINK BEARING MATRIX FOR OUTPUT
      PHI=ABS(BRG-LB)
C     ! WIND ANGLE WITH RESPECT TO LINK
      IF (PHI.LE.90.) GO TO 7600
      IF (PHI.GE.270.) GO TO 5000
      PHI=ABS(PHI-180.)
      GO TO 7600
 5000 PHI=ABS(PHI-360.)
C     ! SET ELEMENT GROWTH BASE
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
C     ! CONVERSION OF PHI FROM DEGREES TO RADIANS
      IF (PHI.GT.1.5706) PHI=1.5706
      IF (PHI.LT.0.00017) PHI=0.00017
C
C *****  DEPRESSED SECTION  *****
C
      IF (HL(IL).LT.-1.5) GO TO 7700
      DSTR=1.
      HDS=1.
      GO TO 7800
 7700 HDS=HL(IL)
      DSTR=0.72*ABS(HDS)**0.83
C     ! RESIDENCE TIME FACTOR
C
C *****  SIGMA Z POWER CURVE  *****
C
 7800 TR=DSTR*W2/U
C     ! RESIDENCE TIME
      SGZ1=ALOG((1.8+0.11*TR)*(ATIM/30.)**0.2)
C     ! ALOG(SIGMA Z) AT W2
      PZ2=(SZ10-SGZ1)/(DREF-ALOG(W2))
      PZ1=EXP((SZ10+SGZ1-PZ2*(DREF+ALOG(W2)))/2.)

      !PRINT *, 'PZ1, PZ2 are: ', PZ1, PZ2
C
C *****  END OF LINK ROUTINE  *****
C
C
C *****  RECEPTOR LOOP  *****
C
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
      XPRI=XR(IR)+DVIR*XVEC
      YPRI=YR(IR)+DVIR*YVEC
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
      IF (DABS(D).LE.D2) Z=ZR(IR)-HL(IL)
      IF (DABS(D).GT.D2)
     *    Z=ZR(IR)-HL(IL)*(1.-(DABS(D)-W2)/(2.*ABS(HL(IL))))
      GO TO 3059
 5750 Z=ZR(IR)

 3059 CONC = 0.
      CALL CAL3CONC(U,W,W2,H,D,Z,UWL,DWL,
     +              MIXH,BASE,PHI,Q1,
     +              PY1,PY2,PZ1,PZ2,DSTR,
     +              CONC)
      C(IR,IM)=C(IR,IM)+CONC

 6000 CONTINUE

C     PRINT *, "C is:", C

 8000 CONTINUE

 8500 CONTINUE

 9000 END SUBROUTINE

