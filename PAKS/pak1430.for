C=======================================================================
C
C  MODEL 30 : KULONOV OTPOR ZA TRENJE VLAKNA (FIBER), SRBIN MODEL 30
C             ZA STAP
C
C               D1M30
C               TAU130
C
C======================================================================
      SUBROUTINE D1M30(DEF,TAU,NLM,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA POZIVANJE PODPROGRAMA ODREDJIVANJE VISKOZNE SILE
C     (NAPONA) I KRUTOSTI
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
C
      LFUN=MREPER(1)
C
      NLM1=(NLM-1)*IDVA
      NLL=NE*IDVA
      NLL10 = 10*NLM1
      NLL5 = 5*NLM1 
      NE10 = NE*6*IDVA
      NE5 = NE*16*IDVA
C
      LEL=LPLAST+NLM1
      LP0 = LEL + NLL
      LSH = LP0 + NLL
      LSIG = LSH + NLL
      LEPS = LSIG + NLL
      LFT = LEPS + NLL
      LL10 = LPLAST + NE10 + NLL10
      LCYC = LPLAST + NE5 + NLL5  
C
C
      LEL1=LPLAS1+NLM1
      LP01 = LEL1 + NLL
      LSH1 = LP01 + NLL
      LSIG1 = LSH1 + NLL
      LEPS1 = LSIG1 + NLL
      LFT1 = LEPS1 + NLL
      LL101 = LPLAS1 + NE10 + NLL10
      LCYC1 = LPLAS1 + NE5 + NLL5  
C
      CALL TAU130(A(LEL),A(LP0),A(LSH),A(LSIG),A(LEPS),A(LFT),
     1            A(LL10),A(LCYC),A(LEL1),A(LP01),
     2            A(LSH1),A(LSIG1),A(LEPS1),A(LFT1),
     3            A(LL101),A(LCYC1),A(LFUN),TAU,DEF,IRAC)
      RETURN
      END
C======================================================================
      SUBROUTINE TAU130(ELST,P0T,PSH,SIGV,EPSV,FT,EL10,CYC,
     1                  ELST1,P0T1,PSH1,SIGV1,EPSV1,FT1,EL101,
     2                  CYC1,FUN,TAU,DEF,IRAC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTIC MATERIAL WITH INTERNAL FRICTION
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /VISK1/ APR
C
      DIMENSION EL10(*),EL101(*),CYC(*),CYC1(*)
      DIMENSION FUN(14,*)
      DIMENSION C2(10),SS(10),SST(10),RAT(10),RATT(10),DELTA(10),
     1          DELTAT(10),EL10W(10)
      DATA EL10W/19.D0,17.D0,15.D0,13.D0,11.D0,9.D0,7.D0,5.D0,3.D0,1.D0/
      DATA C2/35.D0,32.D0,28.D0,24.D0,20.D0,16.D0,12.D0,8.D0,4.D0,0.D0/
C
CE  INITIAL DATA
C
      TOL = 1.D-8
      TOLK = 1.D-5
      TOLE = 1.D-10 
      BETA = 0.D0
C      
      ISH=PSH
      ISH1=PSH1
C
      EE=FUN(1,MAT)
      P0=FUN(10,MAT)
      AMI=FUN(11,MAT)
      CNI=FUN(12,MAT)
      EL=FUN(13,MAT)
      AT=FUN(14,MAT)
C
      IF (KOR.LE.1) THEN
          P0T = P0
          DO 2 I = 1,5
          CYC1(I) = 0.D0
    2     CYC(I) = 0.D0 
          DO 4 I = 1,10
          RAT(I) = 0.D0
    4     EL10(I) = EL10W(I)
      ENDIF
      DO 7 I = 1,5
    7 CYC1(I) = CYC(I) 
      DEPST = CYC(1)
      SIGPL = CYC(2)
      SIGMN = CYC(3)
      RBARP = CYC(4)
      RBARM = CYC(5)
      DO 3 I = 1,10
      IF (EL10(I).GT.0.D0) THEN
         RAT(I) = EL10(I) - EL10W(I)
         DELTA(I) = 1.D0
      ELSE
         RAT(I) = -(EL10W(I)+EL10(I))
         DELTA(I) = -1.D0
      ENDIF   
      RATT(I) = RAT(I)
    3 DELTAT(I) = DELTA(I)
      EPST = DEF
      EPSA = DABS(EPST)
      SIG = SIGV
      ELSTT = ELST
      TEL = AT/EL
      AMIP0 = AMI*P0T
      AKE = 0.5*EE*TEL/AMIP0
      AKE4 = 4.*AKE
      AKESQ = SQRT(1.+AKE4*EPSA) 
      ET = EE/AKESQ
      ELAST(1,1) = ET
      IF (ISH.EQ.-1) ELAST(1,1) = TOLK*ELAST(1,1)
      P0TT = P0T
      ELSTT = ELST
      IF(IRAC.EQ.2) RETURN
C     ZA ITER=0 UZIMA TAU=SIGV I VRACA SE
      IF (ITER.EQ.0) THEN
         TAU = SIGV
         RETURN
      ENDIF
      IT = 0
      ITMAX = 30
      P0IN = P0TT 
      EEPS = EE*DEF
      AMITEL = AMI/TEL
      AMI5 = 0.2*AMITEL
      AMI50 = 0.1*AMI5
      DEPS = DEF - EPSV
      DEPSA = DABS(DEPS)
C
C     ITERACIJE PO P0TT
C     NEMA ITERACIJA PO P0TT ZA STAP 
      P0TT = P0 
      IZNAK = 1
      IF (DEPS.LT.0.D0) IZNAK = -1
      IF (P0TT.LT.TOL) ISH = -1
      IF  (ISH.EQ.-1) THEN
          TAU = 0.D0
          ELAST(1,1) = TOLK*ELAST(1,1)
          GO TO 200 
      ENDIF
      IF (DEPSA.LT.TOLE) GO TO 30
      IF (ISH.EQ.1) GO TO 120
C
      IPL = 0
      IMIN = 0 
      DO 25 I = 1,10
      IF (EL10(I).GT.0.D0) IPL = 1
      IF (EL10(I).LT.0.D0) IMIN = 1
   25 CONTINUE
      IF (DEPS.GT.TOLE .AND. IMIN.EQ.1) GO TO 120
      IF (DEPS.LT.TOLE. AND. IPL.EQ.1) GO TO 120
C
C     CELA DUZINA PROKLIZAVANJA IMA JEDAN SMER NAPONA KLIZANJA
C 
   30 ISH = 0
      AMIP0 = AMI*P0TT
      AKE = 0.5*EE*TEL/AMIP0
      AKE4 = 4.*AKE
      AKESQ = SQRT(1.+AKE4*EPSA) 
      ET = EE/AKESQ
      TAU = AMIP0/TEL*(-1.+AKESQ)
      IF(DEF.LT.0.D0) TAU = - TAU
      GO TO 190
C
C     PETLJA PROMENJLIVOG OPTERECENJA
C
  120 SIGOR = SIGPL
      RBARS = RBARP
      IF (IZNAK.LT.0) THEN
         SIGOR = SIGMN
         RBARS = RBARM
      ENDIF 
      CON5 = 2.*AMI5*P0TT*ELSTT
      CON50 = AMI50*P0TT*ELSTT*ELSTT
C
C     BROJ SEGMENATA SA PROMENJENIM SMEROM SMICANJA
C
      ISH = 1
      A11T = 0.D0
      DFSH = CON5*RBARS
      DO 125 I = 1,10
      DELTAT(I) = DELTA(I)
  125 RATT(I) = RAT(I)
      DO 126 I = 2,9
  126 SS(I) = EL10W(I) - (C2(I)+2.*RAT(I))*RAT(I)
      SS(1) = 19.0 - (35.+3.*RAT(1))*RAT(1)        
      SS(10) = 1. - 2.*RAT(10)*RAT(10)
      DO 127 I = 1,10
      SST(I) = SS(I)
  127 A11T = A11T + SS(I)*DELTA(I)
      SIG0 = EE*EPST 
      SIGNN = SIG0 - CON50*A11T
      DSIG = SIGNN - SIGOR
      NN = 0
  130 NN = NN + 1
      NN10 = 10-NN+1
      DSIG1 = DSIG
      DFSH1 = DFSH
      A11T1 = A11T
      SIGNN1 = SIGNN 
      RATIO = RAT(NN10) 
      RATT(NN10) = 0.D0    
      IF (IZNAK.GT.0) THEN
         IF (DELTA(NN10).GT.0.D0) THEN
            IF(RATIO.LT.TOL) RATIO = 0.D0
            IF (RATIO.LT.TOL.AND.NN.LT.10) GO TO 130   
            RBAR = RATIO
         ELSE
            RBAR = 1. - RATIO 
            DELTAT(NN10) = -DELTA(NN10)
         ENDIF
      ELSE
         IF (DELTA(NN10).LT.0.D0) THEN
            IF (RATIO.LT.TOL) RATIO = 0.D0
            IF (RATIO.LT.TOL.AND.NN.LT.10) GO TO 130
            RBAR = RATIO
         ELSE
            RBAR = 1. - RATIO
            DELTAT(NN10) = -DELTA(NN10)
         ENDIF
      ENDIF 
      SST(NN10) = EL10W(NN10)
      A11T = A11T - SS(NN10)*DELTA(NN10) + SST(NN10)*DELTAT(NN10)
      SIGNN = SIG0 - CON50*A11T
      DSIG = SIGNN - SIGOR
      DSIGA = DABS(DSIG)
      DFSH = DFSH + CON5*RBAR
      RBARS = RBARS + RBAR
      SIGNNA = DABS(SIGNN)
      DFSHR = DFSH - DSIGA
      IF (SIGNNA.GT.TOL) DFSHR = DFSHR/SIGNNA
      IF (DABS(DFSHR).LE.TOL) GO TO 175 
      IF(DFSHR.GT.0.D0) GO TO 170
      IF (NN.EQ.10) GO TO 30
      GO TO 130
C     RATIO ZA SEGMENT NN
  170 DELTAT(NN10) = -DELTAT(NN10)
      RAT1 = RATIO
      IF (IZNAK.LT.0) THEN
          IF (DELTA(NN10).LT.0.D0) RAT1 = 0.D0
          DELTAT(NN10) = 1.D0
          IF (NN.EQ.10) THEN
            AA = -3.*CON50
            BB = -CON5 - CON50*(35.+6.*RAT1)
            CC = -DSIG1 - CON50*(SS(NN10)-19.+(35.+3.*RAT1)*RAT1)-DFSH1
          ELSE
            AA = -2.*CON50
            BB = -CON5 - CON50*(C2(NN10)+4.*RAT1)
            CC = -DSIG1 - CON50*(SS(NN10)-EL10W(NN10)+(C2(NN10)+
     1           2.*RAT1)*RAT1) - DFSH1
          ENDIF     
      ELSE
         IF (DELTA(NN10).GT.0.D0) RAT1 = 0.D0
         DELTAT(NN10) = -1.D0
         IF (NN.EQ.10) THEN
            AA = 3.*CON50
            BB = CON5 + CON50*(35.+6.*RAT1)
            CC = -DSIG1 + CON50*(SS(NN10)-19.+(35.+3.*RAT1)*RAT1)+DFSH1
         ELSE
            AA = 2.*CON50
            BB = CON5+CON50*(C2(NN10)+4.*RAT1)
            CC = -DSIG1 + CON50*(SS(NN10)-EL10W(NN10)+(C2(NN10)+
     1            2.*RAT1)*RAT1) + DFSH1
         ENDIF
      ENDIF
      DET = BB*BB - 4.*AA*CC
      IF (DET.LT.0.D0) THEN
          RATNN = 0.D0
          GO TO 177
      ENDIF 
      AA2 = 2.*AA
      DETSQ = SQRT(DET)
      RATN1 = (-BB+DETSQ)/AA2
      RATN2 = (-BB-DETSQ)/AA2
      I1 = 0
      I2 = 0
      IF (RATN1.GE.0.D0 .AND. RATN1.LE.1.D0) I1 = 1
      IF (RATN2.GE.0.D0 .AND. RATN2.LE.1.D0) I2 = 1
      IF (I1.EQ.0 .AND.I2.EQ.0) THEN
         IF (NN.EQ.10. AND. (RATN1.GT.1.D0.OR.RATN2.GT.1.D0)) GO TO 30   
         RATNN = 0.D0
         GO TO 177 
      ENDIF
      IF (I1.EQ.1 .AND. I2.EQ.0) RATNN = RATN1
      IF (I1.EQ.0 .AND. I2.EQ.1) RATNN = RATN2 
      IF (I1.EQ.1 .AND. I2.EQ.1) THEN
         RATNN = RATN1
         IF (RATN1.GT.RATN2) RATNN = RATN2
      ENDIF  
  177 RATA = RATIO + RATNN
      IF (IZNAK.GT.0) THEN
         IF (DELTA(NN10).GT.0.D0) THEN
            RATA = RATIO - RATNN
            IF (RATA.LT.TOL) RATA = 0.D0
            DELTAT(NN10) = 1.D0
C             RATA = 1.-RATIO+RATNN 
C             DELTAT(NN10) = -1.D0
C             IF (RATA.GT.1.D0) THEN
C                RATA = 0.D0
C                DELTAT(NN10) = 1.D0
C             ENDIF
         ENDIF
      ELSE
         IF (DELTA(NN10).LT.0.D0) THEN
            RATA = RATIO - RATNN
            IF (RATA.LT.TOL) RATA = 0.D0
            DELTAT(NN10) = -1.D0
C             DELTAT(NN10) = 1.D0
C             RATA = 1.-RATIO+RATNN
C             IF (RATA.GT.1.D0) THEN
C                RATA = 0.D0
C                DELTAT(NN10) = -1.D0
C             ENDIF
         ENDIF
      ENDIF
      IF (RATA.GE.1. .OR. (1.-RATA).LE.TOL) THEN
         RATA = 0.D0
         DELTAT(NN10) = 1.D0
         IF (IZNAK.LT.0.D0) DELTAT(NN10) = -1.D0
      ENDIF
      RATT(NN10) = RATA 
      RBARS = RBARS - RBAR + RATNN
      IF (NN.EQ.10) THEN
         SST(1) = 19.0 - (35.+3.*RATT(1))*RATT(1) 
      ELSE IF (NN.EQ.1) THEN
         SST(10) = 1. - 2.*RATT(10)*RATT(10)
      ELSE
         SST(NN10) = EL10W(NN10) - (C2(NN10)+2.*RATT(NN10))*RATT(NN10)
      ENDIF
      A11T = A11T1 - SS(NN10)*DELTA(NN10) + SST(NN10)*DELTAT(NN10)
      SIGNN = SIG0 - CON50*A11T
  175 TAU = SIGNN
      DSIGA = DABS(SIGNN-SIG)
      A50 = AMI50*P0TT*A11T
      A5 = CON5/ELSTT*RBARS
      ET = EE/(1.+2.*A50*DSIGA/(A5*A5)) 
C*** zasada uzeto et=ee
      AMIP0 = AMI*P0TT
      ET = EE
      ELAST(1,1) = ET 
C
C     RACUNANJE ELST = ELDIM/EL  (EL-STAR)
  190 IF (ISH.EQ.1) GO TO 195
      ELSTT = 0.5*TEL*DABS(TAU)/AMIP0
  195 IF ((ELSTT-0.5).GT.TOL) THEN
         ISH = -1
         TAU = 0.D0 
         ELAST(1,1) = TOLK*ELAST(1,1)
      ENDIF
C
C     AZURIRANJE VELECINA ZA SLEDECI KORAK 
C
  200 ISH1 = ISH
      PSH1 = ISH1
      P0T1 = P0TT
      ELST1 = ELSTT
      SIG1 = TAU
      FT1 = TAU*APR
      SIGV1 = TAU
      EPSV1 = DEF
      CYC1(1) = DEPS
      IF (ISH.EQ.1) GO TO 210
C     KOREKCIJA ZBOG VRACANJA SA LABELE 120
      CYC1(4) = 0.D0
      CYC1(5) = 0.D0
      IF (IZNAK.GT.0) THEN
         CYC1(3) = SIG1
         DO 202 I = 1,10
  202    EL101(I) = EL10W(I)
      ELSE
         CYC1(2) = SIG1
         DO 203 I = 1,10
  203    EL101(I) = -EL10W(I)
      ENDIF
  210 IF (ISH.NE.1) RETURN  
      IF (DEPS.GT.0.D0) THEN
         CYC1(3) = SIG1
         CYC1(5) = 0.D0
         CYC1(4) = RBARS 
         IF (DEPST.LT.0.D0) CYC1(2) = SIG
      ELSE
         CYC1(2) = SIG1
         CYC1(4) = 0.D0
         CYC1(5) = RBARS
         IF (DEPST.GT.0.D0) CYC1(3) = SIG
      ENDIF 
      DO 205 I = 1,10
      IF (DELTAT(I).GT.0.D0) THEN
         EL101(I) = EL10W(I) + RATT(I)
      ELSE
         EL101(I) = -EL10W(I) - RATT(I)
      ENDIF
  205 CONTINUE
      RETURN
C
      END
