C=======================================================================
C
C   PLASTICNOST LJUSKE  -  MESOVITO OJACANJE    (17.10.1992)
C
C=======================================================================
      SUBROUTINE D4M6(TAU,DEF,IRAC,LPOCG,LPOC1,IND3D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE  SPACE IN WORKING VECTOR 
C
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFPT=LDEFT  + 6*IDVA
      LALFAT=LDEFPT + 6*IDVA
      LTEQT =LALFAT + 6*IDVA
      LDQPT =LTEQT  + 1*IDVA
      LIPL  =LDQPT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6*IDVA
      LDEFP1=LDEFT1 + 6*IDVA
      LALFA1=LDEFP1 + 6*IDVA
      LTEQT1=LALFA1 + 6*IDVA
      LDQPT1=LTEQT1 + 1*IDVA
      LIPL1 =LDQPT1 + 1*IDVA
C
      CALL TAUI46(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC,IND3D)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUI46( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC,IND3D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC MATERIAL , MIXED HARDENING
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /COEFSM/ COEF(3),ICOEF
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /ITERBR/ ITER
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(6),GLITL(6),COEFE(2)
      DIMENSION FUN(2,MATE,*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C      COEFE(1)=0.833333333D0
      COEFE(1)=COEF(1)
      COEFE(2)=COEF(2)
C
CE  INITIAL DATA
C
      IPL =PL
C
CE  E,V,TEQY0,CY,AN,EM
C
      E     =FUN(1,MAT,1)
      V   =FUN(2,MAT,1)
      TEQY0 =FUN(1,MAT,2)
      CY    =FUN(2,MAT,2)
      AN    =FUN(1,MAT,3)
      EM    =FUN(2,MAT,3)
      HS    =FUN(1,MAT,4)
      ISIMO=0
      IF(DABS(HS).GT.1.D-10) ISIMO=1
C
CE    AUXILIARY CONSTANTS
C
      DVT =2.D0/3.
      C1    =(2.-V)/3./(1.-V)
      CNI   =1.-C1
      CYAN23=DVT*CY*AN
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-2.
      AE    =(1.+V)/E
      AEINV =1./AE
      CM    =E/(1.-2.*V)
      C3    =DVT*(1.+V)
C
      AE5   =AE/COEFE(2)
      AE6   =AE/COEFE(1)
      BA    =1.5/E
C
CE    YIELD STRESS
C     NAPON TECENJA (2.3.38)
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      CALL MEL96(COEFE,ETP,1,IND3D)
      IF(IRAC.EQ.2) RETURN
      CALL CLEAR(TAU,6)
C
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
    5 DEFPT(I)=0.5*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, ESEKUNDUM
C     ELASTICNE DEFORMACIJE, (PROBNA RESENJA ZA DDEFQP=0)
      IF(IATYP.LT.4)THEN
C       SREDNJA ELASTICNA DEFORMACIJA (2.3.69) (ZA DDEFQP=0)
        EMT=(1.-2.*V)*(DEF(3)-DEFPT(3))/3.
C       PROBNE DEVIJATORSKE ELASTICNE DEFORMACIJE (2.3.71)
        DEFDS(3)=C3*(DEF(3)-DEFPT(3))
        DO 10 I=5,6
   10   DEFDS(I)=0.5*DEF(I)-DEFPT(I)
      ELSE
        EMT = (1.-2.*V)*(DEF(3))/3.
        DEFDS(3)= DEF(3)
        DO 11 I=5,6
   11   DEFDS(I)=0.5*DEF(I)
      ENDIF
        DEFDS(1)=-0.5*DEFDS(3)
        DEFDS(2)=-0.5*DEFDS(3)
        DEFDS(4)= 0.
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD) , (GLITL), (SHET)
C
      DO 40 I=1,6
C     SHET/2G - ZA GREDU
      GLITL(I)=DEFDS(I)-AE*ALFAT(I)
C     PROBNI ELASTICNI DEVIJATORSKI NAPONI (2.3.54)
      TAUD(I) =DEFDS(I)*AEINV
C     PROBNI ELASTICNI RADIJUS NAPONI (2.3.74)
   40 SHET(I) =TAUD(I)-ALFAT(I)
C
CE   2)  CHECK FOR YIELDING
C     PROBNI EFEKTIVNI ELASTICNI RADIJUS NAPON (2.3.50)
      TEQ=DSQRT(1.5*TENDOT(SHET))
C      WRITE(3,*) 'ITER,JG,EMT',ITER,JG,EMT
C      WRITE(3,*) 'TEQY,TEQ',TEQY,TEQ
C      CALL WRR6(TAU,6,'TAU0')
C      CALL WRR6(DEF,6,'DEFU')
C      CALL WRR6(DEFPT,6,'DFPT')
C      CALL WRR6(ALFAT,6,'ALFT')
C      CALL WRR6(DEFDS,6,'DFDS')
C      CALL WRR6(GLITL,6,'GTL0')
C      CALL WRR6(TAUD,6,'TAD0')
C      CALL WRR6(SHET,6,'SHE0')
C     PROVERA TECENJA (2.3.49)
      IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
C       ELASTICNA OBLAST
        TEQ  =TEQY
        DEFQP=DEFQPT
        DO 600 I=1,6
  600   DEFP(I)=DEFPT(I)
        GO TO 500
      ENDIF
C
CE   3)   SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C     ELASTO-PLASTICNA OBLAST 
      PL1=1.0D0
C
CE    BISECTION
C     ODREDJIVANJE GRANICA ZA DDEFQP
      DEFQP=DEFQPT
      IF(DEFQP.LE.1.D-4) DEFQP=1.D-4
      IF(ISIMO.EQ.0)THEN
        EP2=AN*CY*DEFQP**AN1      
      ELSE
        EP2=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
C
      IB = 0
      IT = 0
      AF    = 5.D0
      DDD   = 0.1*(TEQ-TEQY)/EP2
      FP    = TEQ - TEQY
      TAUY  = TEQY
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
C     PROBNO RESENJE ZA DDEFQP A DOBICE SE RESAVANJEM FUNKCIJE TECENJA
      DDEFQP= DDD
  100 IT=IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      ENDIF
C
C     EFEKTIVNA PLASTICNA DEFORMACIJA (2.3.57)
      DEFQP=DEFQPT+DDEFQP
      IF(ISIMO.EQ.0)THEN
C       MODUL (2.3.29)
        CC=CYAN23*DEFQP**AN1
C       MODUL (2.3.46) ILI (2.3.30)
        EM1CC =EM1*CC
      ELSE
        CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
        EM1CC =EM1*CC
      ENDIF
C
        TAUY=TEQY0+CY*(EM*DEFQP)**AN
        DLAM=1.5*DDEFQP/TAUY
        SHET(3)=GLITL(3)/(AE+(C3+EM1CC*AE)*DLAM)
        SHET(5)=GLITL(5)/(AE+(1.+EM1CC*AE)*DLAM)
        SHET(6)=GLITL(6)/(AE+(1.+EM1CC*AE)*DLAM)
        SHET(1)=-0.5*SHET(3)
        SHET(2)=-0.5*SHET(3)
        SHET(4)= 0.
        TEQ =DSQRT(1.5*TENDOT(SHET))
        FB = TEQ-TAUY
C
C      WRITE(3,*) 'IT,FP,FB',IT,FP,FB
C      WRITE(3,*) 'FM,AF',FM,AF
C      WRITE(3,*) 'DEPBM,DEPBP',DEPBM,DEPBP
C      WRITE(3,*) 'DDD,DEFQPT',DDD,DEFQPT
C      WRITE(3,*) 'DEFQP,DDEFQP',DEFQP,DDEFQP
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
      DUM =1.+EM1CC*DLAM
      DO 165 I=1,6
C     DEVIJATORSKI NAPON (2.3.43)
  165 TAUD(I)=ALFAT(I)+DUM*SHET(I)
C
CE   4)  DETERMINE SOLUTION 
C
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
      DO 170 I=1,6
C       PRIRASTAJ PLASTICNE DEFORMACIJE U KORAKU (2.3.19)
        DDEFPS  =DLAM*SHET(I)
C       UKUPNE PLASTICNE DEFORMACIJE (2.3.39)
        DEFP(I) =DEFPT(I)+DDEFPS
C       POLOZAJNI NAPONI (2.3.42)
        ALFA1(I)=ALFAT(I)+EM1CC*DDEFPS
  170 CONTINUE
C
      IF(IATYP.LT.4)THEN
C       SREDNJA ELASTICNA DEFORMACIJA (2.3.69)
        EMT=(1.-2.*V)*(DEF(3)-DEFP(3))/3.
      ELSE
        EMT=(1.-2.*V)*(DEF(3))/3.
      ENDIF
C      WRITE(3,*) ' DLAM,TEQ',DLAM,TEQ
C      CALL WRR6(SHET,6,'SHET')
C      CALL WRR6(TAUD,6,'TAUD')
C      CALL WRR6(DEFP,6,'DEFP')
C      CALL WRR6(ALFA1,6,'ALF1')
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2) THEN
C        GLITL ODGOVARA 3D ELEMENTU A NE GREDI
         IF(IATYP.LT.4)THEN
            GLITL(3)=DEF(3)-EMT-DEFPT(3)-AE*ALFAT(3)
            GLITL(1)=-0.5*GLITL(3)
            GLITL(2)=GLITL(1)
            GLITL(5)=0.5*DEF(5)-DEFPT(5)-AE*ALFAT(5)
            GLITL(6)=0.5*DEF(6)-DEFPT(6)-AE*ALFAT(6)
            GLITL(4)=0.
         ELSE
            GLITL(3)=DEF(3)-EMT-AE*ALFAT(3)
            GLITL(1)=-0.5*GLITL(3)
            GLITL(2)=GLITL(1)
            GLITL(5)=0.5*DEF(5)-AE*ALFAT(5)
            GLITL(6)=0.5*DEF(6)-AE*ALFAT(6)
            GLITL(4)=0.
         ENDIF
         IF(INDCEP.EQ.0)
     1      CALL CEP96(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &                 HS,ISIMO)
C         IF(INDCEP.NE.0)
C     1      CALL CEP46(SHET,DLAM,TEQ,DEFQP,E,EM,AE5,AE6,BA,CC,CY,AN)
C         CALL WRR6(ETP,36,'ETP ')
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      IF(IATYP.LT.4)THEN
C       SREDNJI NAPON (2.3.6)
        TAUM=CM*EMT
      ELSE
C       TAUM=
      ENDIF
C     NAPONI NA KRAJU KORAKA (2.3.56)
C      TAU(1)=0.D0
C      TAU(2)=0.D0
C      TAU(3)=1.5*TAUD(3)
      DO 45 I=1,3
   45 TAU(I)=TAUD(I)+TAUM
      TAU(4)=0.D0
      DO 46 I=5,6
      TAU(I)=TAUD(I)
   46 DEFP(I)=2.*DEFP(I)
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF(IATYP.LT.4)THEN
C       UKUPNE DEFORMACIJE U RAVNI PRESEKA (2.3.75)
        DEF(1)=-V*(DEF(3)-DEFP(3))+DEFP(1)
        DEF(2)=-V*(DEF(3)-DEFP(3))+DEFP(2)
      ELSE
        DEF(1)=-V*DEF(3)
        DEF(2)=-V*DEF(3)
      ENDIF
      DO 290 I=1,6
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
C      WRITE(3,*) 'TAUM,EMT',TAUM,EMT
C      CALL WRR6(DEF1,6,'DEF1')
C      CALL WRR6(TAU1,6,'TAU1')
      RETURN
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI46')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TAUI46')
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE CEP96(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &                 HS,ISIMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ETP )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      DIMENSION GLITL(*),CP(3,3),EL(6,6) 
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      DVT   =DVA*TR
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-DVA
      CHET  =EM1*CC
C
      EMEPH = 0.D0
      IF(ISIMO.EQ.0)THEN
        EPPRIM=AN*AN1*CY*DEFQP**AN2      
        IF(EM.GT.1.D-10) EMEPH = EM * AN*CY*(EM*DEFQP)**AN1      
      ELSE
        EPPRIM=AN*AN*CY*DEXP(-AN*DEFQP)
        IF(EM.GT.1.D-10) EMEPH = EM * (AN*CY*DEXP(-AN*DEFQP)+HS)
      ENDIF
      EPHETP=EMEPH+EM1*EPPRIM*DDEFQP
C
      AP=(AE*(DVT*EPHETP+CHET)+1.)*DSQRT(TRIPO*TENDOT(GLITL))
      BP=TRIPO/AP/TEQ*(1.-EMEPH*DDEFQP/TEQ)
      DP=AE+(1.+AE*CHET)*DLAM
      DP=(BP-DVT*EM1*DLAM*DLAM*EPPRIM/AP)/DP/DP
      AA=(1.+CHET*DLAM)/(AE+(1.+AE*CHET)*DLAM)
C
      DO 25 I=1,3
       DO 20 J=I,3
        CP(I,J)   =-DP*GLITL(I)*GLITL(J)
   20  CONTINUE
       CP(I,I)   =CP(I,I)+AA
   25 CONTINUE
      DO 30 I=1,3
      DO 30 J=4,6
        EL(I,J)=-DP*GLITL(I)*GLITL(J)
   30 CONTINUE
      DO 40 I=4,6
      DO 40 J=I,6
        EL(I,J)=-DP*GLITL(I)*GLITL(J)
   40 CONTINUE
      DO 45 I=4,6
        EL(I,I)=EL(I,I)+0.5*AA
   45 CONTINUE
C
      EL(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      EL(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      EL(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      EL(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      EL(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      EL(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
C      
      DO 50 I=1,6
      DO 50 J=I,6
        EL(J,I)=EL(I,J)
   50 CONTINUE
C      
      CALL CLEAR(ETP,36)
      E12=EL(1,1)*EL(2,2)
      A12=1.-EL(2,1)*EL(1,2)/E12
      A13=(EL(1,2)*EL(2,3)/E12-EL(1,3)/EL(1,1))/A12
      A23=(EL(2,1)*EL(1,3)/E12-EL(2,3)/EL(2,2))/A12
      A15=(EL(1,2)*EL(2,5)/E12-EL(1,5)/EL(1,1))/A12
      A25=(EL(2,1)*EL(1,5)/E12-EL(2,5)/EL(2,2))/A12
      A16=(EL(1,2)*EL(2,6)/E12-EL(1,6)/EL(1,1))/A12
      A26=(EL(2,1)*EL(1,6)/E12-EL(2,6)/EL(2,2))/A12
      ETP(3,3)=EL(3,1)*A13+EL(3,2)*A23+EL(3,3)
      ETP(3,5)=EL(3,1)*A15+EL(3,2)*A25+EL(3,5)
      ETP(3,6)=EL(3,1)*A16+EL(3,2)*A26+EL(3,6)
      ETP(5,5)=EL(5,1)*A15+EL(5,2)*A25+EL(5,5)
      ETP(5,6)=EL(5,1)*A16+EL(5,2)*A26+EL(5,6)
      ETP(6,6)=EL(6,1)*A16+EL(6,2)*A26+EL(6,6)
      ETP(5,3)=ETP(3,5)
      ETP(6,3)=ETP(3,6)
      ETP(6,5)=ETP(5,6)
C
      RETURN
      END
C======================================================================
      SUBROUTINE CEP46(SHET,DLAM,TEQ,DEFQP,E,EM,AE5,AE6,BA,CC,CY,AN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ETP )    LJUSKA
CE     ELASTO-PLASTIC  CEP MATRIX        SHELL
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      DIMENSION SHET(*)
C
      DVA   =2.D0
      TR    =1.D0/3.
      DVT   =DVA*TR
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-DVA
      CHET  =EM1*CC
C
      EPPRIM=AN*AN1*CY*DEFQP**AN2
      EPSEK=EPPRIM*DVT*EM1      
      EMEPH = 0.D0
      IF(EM.GT.1.D-10) EMEPH = EM * AN*CY*(EM*DEFQP)**AN1      
      ELP=(DVT-DLAM*EMEPH)/TEQ
C
      DB=BA+(1.+BA*CHET)*DLAM
      DA5=AE5+(1.+AE5*CHET)*DLAM
      DA6=AE6+(1.+AE6*CHET)*DLAM
      A3=1./DB
      A5=1./DA5
      A6=1./DA6
      B3=-SHET(3)*((1.+BA*CHET)*ELP+BA*EPSEK*DLAM)/DB
      B5=-SHET(5)*((1.+AE5*CHET)*ELP+AE5*EPSEK*DLAM)/DA5
      B6=-SHET(6)*((1.+AE6*CHET)*ELP+AE6*EPSEK*DLAM)/DA6
      W3=9./4.*SHET(3)*A3
      W5=3.*SHET(5)*A5
      W6=3.*SHET(6)*A6
      W0=9./4.*SHET(3)*B3+3.*SHET(5)*B5+3.*SHET(6)*B6-EMEPH*TEQ
C
      DEPQ3=-W3/W0
      DEPQ5=-W5/W0
      DEPQ6=-W6/W0
      DEP33=(SHET(3)*ELP+DLAM*B3)*DEPQ3+DLAM*A3
      DEP35=(SHET(3)*ELP+DLAM*B3)*DEPQ5
      DEP36=(SHET(3)*ELP+DLAM*B3)*DEPQ6
      DEP55=(SHET(5)*ELP+DLAM*B5)*DEPQ5+DLAM*A5
      DEP56=(SHET(5)*ELP+DLAM*B5)*DEPQ6
      DEP66=(SHET(6)*ELP+DLAM*B6)*DEPQ6+DLAM*A6
C
      CALL CLEAR(ETP,36)
      ETP(3,3)=E*(1.-DEP33)
      ETP(3,5)=-0.5*E*DEP35
      ETP(5,3)=ETP(3,5)
      ETP(3,6)=-0.5*E*DEP36
      ETP(6,3)=ETP(3,6)
      ETP(5,5)=0.5*(1.-DEP55)/AE5
C PROVERI
C      ETP(5,6)=-0.5*DEP56/AE
      ETP(5,6)=-DEP56/(AE5+AE6)
C
      ETP(6,5)=ETP(5,6)
      ETP(6,6)=0.5*(1.-DEP66)/AE6
      RETURN
      END
C======================================================================
      SUBROUTINE MEL96(COEFE,ETP,IOPT,IND3D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE ELAST
CE     ELAST MATRIX
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /MATIZO/ E,V
      COMMON /PODTIP/ IPODT
C
      DIMENSION ETP(6,*),COEFE(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MEL96 '
C
C      WRITE(3,*) 'E,V,IND3D,IPODT',E,V,IND3D,IPODT
C      WRITE(3,*) 'COEFE,IOPT',COEFE(1),COEFE(2),IOPT
C     NULOVANJE ETP
      CALL CLEAR(ETP,36)
C     MATRICA ETP ZA 3D
      IF(IPODT.EQ.2) THEN
         IF(IND3D.EQ.0) ETP(3,3)=E
         IF(IND3D.EQ.1) THEN
            ETP(3,3)=E
            ETP(5,5)=E/(2.*(1.+V))
            ETP(6,6)=ETP(5,5)
         ENDIF
         IF(IND3D.EQ.2) THEN
            ETP(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
            ETP(2,2)=ETP(1,1)
            ETP(3,3)=ETP(1,1)
            ETP(1,2)=ETP(1,1)*V/(1.-V)
            ETP(2,1)=ETP(1,2)
            ETP(1,3)=ETP(1,2)
            ETP(3,1)=ETP(1,2)
            ETP(2,3)=ETP(1,2)
            ETP(3,2)=ETP(1,2)
            ETP(4,4)=ETP(1,1)*(1.-2.*V)/(1.-V)/2.
            ETP(5,5)=ETP(4,4)
            ETP(6,6)=ETP(4,4)
         ENDIF
      ENDIF
C     MATRICA ETP ZA LJUSKU
      IF(IPODT.EQ.1) THEN
         EA=E/(1.D0-V*V)
         E1=0.5D0*(1.D0-V)
         IF(IOPT.EQ.0) THEN
C           NULOVANJE ELAST
            CALL CLEAR(ELAST,36)
            ELAST(2,2)=EA
            ELAST(3,3)=ELAST(2,2)
            ELAST(4,4)=ELAST(2,2)*E1
            ELAST(6,6)=ELAST(4,4)
            RETURN
         ENDIF
C
         IF(IND3D.EQ.-2) ETP(5,5)=E/(2.*(1.+V))
         IF(IND3D.EQ.-1) ETP(2,2)=E
         IF(IND3D.EQ.0) ETP(3,3)=E
         IF(IND3D.EQ.1) THEN
             ETP(3,3)=E
             ETP(5,5)=E/(2.*(1.+V))
             IF(NMODM.LT.5) THEN
                ETP(6,6)=COEFE(2)*ETP(5,5)
             ELSE
                ETP(6,6)=ETP(5,5)
             ENDIF
         ENDIF
         IF(IND3D.EQ.2) THEN
            ETP(2,2)=EA
            ETP(3,3)=EA
            ETP(2,3)=EA*V
            ETP(3,2)=ETP(2,3)
            ETP(5,5)=EA*E1
            IF(NMODM.LT.5) THEN
               ETP(4,4)=COEFE(1)*ETP(5,5)
               ETP(6,6)=COEFE(2)*ETP(5,5)
            ELSE
               ETP(4,4)=ETP(5,5)
               ETP(6,6)=ETP(5,5)
            ENDIF
         ENDIF
      ENDIF
C     MATRICA ETP ZA GREDU
      IF(IPODT.EQ.0) THEN
         IF(IND3D.EQ.0) ETP(3,3)=E
         IF(IND3D.EQ.1) THEN
            ETP(3,3)=E
            IF(NMODM.LT.5) THEN
               ETP(5,5)=COEFE(1)*E/(2.*(1.+V))
               ETP(6,6)=COEFE(2)*E/(2.*(1.+V))
            ELSE
               ETP(5,5)=E/(2.*(1.+V))
               ETP(6,6)=E/(2.*(1.+V))
            ENDIF
         ENDIF
         IF(IND3D.EQ.2) THEN
            ETP(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
            ETP(2,2)=ETP(1,1)
            ETP(3,3)=ETP(1,1)
            ETP(1,2)=ETP(1,1)*V/(1.-V)
            ETP(2,1)=ETP(1,2)
            ETP(1,3)=ETP(1,2)
            ETP(3,1)=ETP(1,2)
            ETP(2,3)=ETP(1,2)
            ETP(3,2)=ETP(1,2)
            ETP(4,4)=ETP(1,1)*(1.-2.*V)/(1.-V)/2.
            ETP(5,5)=ETP(4,4)
            ETP(6,6)=ETP(4,4)
         ENDIF
      ENDIF
      RETURN
      END
