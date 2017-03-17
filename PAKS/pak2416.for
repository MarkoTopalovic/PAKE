C=======================================================================
C
C   TERMO-PLASTICNOST 2/D ELEMENT  -  MESOVITO OJACANJE    (05.01.1994)
C
C    SUBROUTINE D2M14
C               TI214
C               MEL2T
C               CEP214
C
C=======================================================================
      SUBROUTINE D2M14(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
      LTAU  =LPOCG
      LDEFT =LTAU   + 4*IDVA
      LDEFPT=LDEFT  + 4*IDVA
      LALFAT=LDEFPT + 4*IDVA
      LTEQT =LALFAT + 4*IDVA
      LDQPT =LTEQT  + 1*IDVA
      LIPL  =LDQPT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LALFA1=LDEFP1 + 4*IDVA
      LTEQT1=LALFA1 + 4*IDVA
      LDQPT1=LTEQT1 + 1*IDVA
      LIPL1 =LDQPT1 + 1*IDVA
C
      CALL TI214(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI214( PL ,DEFPT,ALFAT,TEMT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEM, DEFQP,
     1                   FUN,MATE,TAU,DEF,TGT,TREF,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTOPLASTICAN MATERIJAL SA MESOVITIM OJACANJEM
C
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
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(4),SEHET(4),ETHERM(4)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR ZA KONTROLNU STAMPU
C
      KSS=4
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) KSS=3
C
CE  INITIAL DATA
C
      IPL =PL
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      DVT   =2.D0/3.
      CZZ   =1.-V
      CNI   =CZZ-V
      CM    =E/CNI
      CNI   =CNI/3./CZZ
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2    =E/(1.+V)
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
      G2C1  =G2*C1
      G2CNI =G2*CNI
C
CE    YIELD STRESS
C
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      CALL MEL2T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      TEM=TGT
      IF(KOR.GT.1.AND.ITER.EQ.0) TGT=TGT-TEMT+TEMP0
      CALL STERM2(ETHERM,TGT)
      ETH = ETHERM(1)
C     ZA ITER=0, TERMICKE SILE ZA POMERANJA OD TERMICKIH DEFORMACIJA
      IF(ITER.EQ.0) GO TO 900
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFPT(3)=0.5*DEFPT(3)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2)-2*ETH)+ETH
          DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNI*(DEF(2)-DEFPT(2))-CTH*ETH
          DEFDS(2)=-CNI*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))-CTH*ETH
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT-DEFPT(1)
          DEFDS(2)=DEF(2)-EMT-DEFPT(2)
          DEFDS(4)=DEF(4)-EMT-DEFPT(4)
        END IF
        DEFDS(3)=0.5*DEF(3)-DEFPT(3)
      ELSE
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2))
          DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
          DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(4)=DEF(4)-EMT
        END IF
        DEFDS(3)=0.5*DEF(3)
      END IF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      DO 40 I=1,KSS
        TAUD(I) =G2*DEFDS(I)
        SEHET(I)=TAUD(I)-ALFAT(I)
   40   SHET(I) =SEHET(I)
      IF (KSS.EQ.3) THEN
        TAUD(4)  =-TAUD(1)-TAUD(2)
        SEHET(4) =-SEHET(1)-SEHET(2)
        SHET(4)  =SEHET(4)
      END IF
      TEQ=DSQRT(1.5*TENDO2(SHET))
      TEQSEH=TEQ
C
CE   2)  CHECK FOR YIELDING
C
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5) THEN
        DEFQP=DEFQPT
        DO 600 I=1,4
  600     DEFP(I)=DEFPT(I)
        GO TO 500
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
C
CE    BISECTION
C
      DEFQP=DEFQPT
      IF (DEFQP.LE.1.D-8) DEFQP=1.D-8
      EALFA=ANCY*DEFQP**AN1
      IF (EALFA.LE.1.D-8) EALFA=1.D-8
C
      IT=0
      IB = 0
      AF = 3.D0
      FP = TEQ-TEQY
      DDD = 0.1*FP/EALFA
      FM     = 0.D0
      DEPBM  = 0.D0
      DEPBP  = 0.D0
      DDEFQP = DDD
C
  100 IT=IT+1
      IB1 = IB
C
      IF (IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS
        STOP
      END IF
C
      DEFQP=DEFQPT+DDEFQP
      TEQ=TEQY0+CY*(EM*DEFQP)**AN
      EALFA=ANCY*DEFQP**AN1
      DTEQ=EALFA*EM**AN
      CHET=DVT*EM1*EALFA
      CHETP=AN1*CHET/DEFQP
      CG2=G2+CHET
C
      IF (KSS.EQ.3) THEN
        TAUY=TEQ
        DLAM=1.5*DDEFQP/TAUY
        DEL=1.+CG2*DLAM
        P1 =1.+(G2C1+CHET)*DLAM
        P2 =G2CNI*DLAM
        P12=P1*P1-P2*P2
        SHET(1)=(P1*SEHET(1)+P2*SEHET(2))/P12
        SHET(2)=(P2*SEHET(1)+P1*SEHET(2))/P12
        SHET(4)=-SHET(1)-SHET(2)
        SHET(3)=SEHET(3)/DEL
        TEQ=DSQRT(1.5*TENDO2(SHET))
        FB=TEQ-TAUY
      ELSE
        FB=-(TEQ+1.5*CG2*DDEFQP-TEQSEH)
      END IF
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE   4)   DEVIATORIC STRESS, PLASTIC STRAIN, BACK STRESS
C
      IF (KSS.NE.3) THEN
        FHETP=DTEQ+1.5*CG2+1.5*CHETP*DDEFQP
        DLAM=1.5*DDEFQP/TEQ
        DEL=1.+CG2*DLAM
        DO 160 I=1,4
  160     SHET(I)=SEHET(I)/DEL
      END IF
C
      DO 170 I=1,4
        DDEFPS  =DLAM*SHET(I)
        DEFP(I) =DEFPT(I)+DDEFPS
        ALFA1(I)=ALFAT(I)+CHET*DDEFPS
  170   TAUD(I) =SHET(I)+ALFA1(I)
C
      IF (KSS.EQ.3) THEN
        IF (IATYP.LT.4) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-2*ETH)+ETH
        ELSE
          EMT = -TAUD(4)/CM
        END IF
      END IF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF (ISKNP.NE.2) THEN
        IF (KSS.EQ.3) THEN
          IF (IATYP.LT.4) THEN
            SEHET(1)=G2*(DEF(1)-EMT-DEFPT(1))-ALFAT(1)
            SEHET(2)=G2*(DEF(2)-EMT-DEFPT(2))-ALFAT(2)
          ELSE
            SEHET(1)=G2*(DEF(1)-EMT)-ALFAT(1)
            SEHET(2)=G2*(DEF(2)-EMT)-ALFAT(2)
          END IF
          SEHET(4)=-SEHET(1)-SEHET(2)
          TEQSEH=DSQRT(1.5*TENDO2(SEHET))
        END IF
        FHETP=DTEQ+1.5*CG2+1.5*CHETP*DDEFQP
        IF (FHETP.LE.1.D-8) FHETP=1.D-8
      IF(INDCEP.EQ.0)
     1  CALL CEP214(SEHET,DLAM,TEQ,G2,CM,FHETP,TEQSEH,CHET,CHETP,DEL,
     &              DTEQ)
      END IF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)/3.D0
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=TAUD(3)
      IF (KSS.EQ.3) THEN
        TAU(4)=0.D0
      ELSE
        TAU(4)=TAUD(4)+TAUM
      END IF
      DEFP(3) =2.*DEFP(3)
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-2*ETH)+DEFP(4)+ETH
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      ELSE
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2))
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      END IF
  900 IF(ITER.EQ.0) THEN
         ETH=ETH*(ELAST(1,1)+ELAST(1,2)+ELAST(1,4))
         TAU(1)=TAU(1)-ETH
         TAU(2)=TAU(2)-ETH
         IF(IETYP.EQ.1) TAU(4)=TAU(4)-ETH
      ENDIF
      DO 290 I=1,4
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI214')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI214'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI214')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI214'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE CEP214(SEHET,DLAM,TEQ,G2,CM,FHETP,TEQSEH,CHET,CHETP,
     &                  DEL,DTEQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION SEHET(*),CP(4,4)
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      AMULT =G2/DEL
C
      AP=TRIPO*G2/FHETP/TEQSEH
      BP=(TRIPO-DLAM*DTEQ)/TEQ
      DP=AMULT*(BP-DLAM*DLAM*CHETP)/DEL*AP
      AA=AMULT*(1.+CHET*DLAM)
C
      DO 25 I=1,4
        DO 20 J=I,4
          CP(I,J)=-DP*SEHET(I)*SEHET(J)
   20   CONTINUE
        CP(I,I)=CP(I,I)+AA
   25 CONTINUE
      ELAST(1,3)=CP(1,3)
      ELAST(2,3)=CP(2,3)
      ELAST(3,3)=CP(3,3)-AA/DVA
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,4)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,4)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,4)+CM)
C
      IF(IETYP.NE.2) THEN
        ELAST(1,4)=TR*(DVA*CP(1,4)-CP(1,1)-CP(1,2)+CM)
        ELAST(2,4)=TR*(DVA*CP(2,4)-CP(1,2)-CP(2,2)+CM)
        ELAST(3,4)=CP(3,4)
        ELAST(4,4)=TR*(DVA*CP(4,4)-CP(1,4)-CP(2,4)+CM)
      ENDIF
C
CE  STATIC CONDESATION
C
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        DO 60 I=1,3
          DO 60 J=I,3
            ELAST(I,J)=ELAST(I,J)-ELAST(I,4)*ELAST(J,4)/ELAST(4,4)
  60    CONTINUE
        DO 61 I=1,4
  61      ELAST(I,4)=0.D0
      END IF
C
      DO 50 I=1,4
        DO 50 J=I,4
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MEL2T
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
CE     NULL ( ELAST )
C
      DO 15 I=1,4
        DO 15 J=I,4
   15     ELAST(I,J)=0.0D0
C
CS     RAVANSKO STANJE NAPONA
CE     PLANE STRESS
C
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        ELAST(1,1)=E/(1.-V*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V
        ELAST(3,3)=ELAST(1,1)*(1.-V)/2.
      ELSE
C
CS     RAVANSKA DEFORMACIJA
CE     PLANE STRAIN
C
        ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V/(1.-V)
        ELAST(3,3)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
c        IF (IETYP.NE.1) GO TO 40
        ELAST(4,4)=ELAST(1,1)
        ELAST(1,4)=ELAST(1,2)
        ELAST(2,4)=ELAST(1,2)
      END IF
C
c   40 DO 50 I=1,4
      DO 50 I=1,4
        DO 50 J=I,4
   50     ELAST(J,I)=ELAST(I,J)
C
      RETURN
      END
C=======================================================================
C
C   TERMO-ELASTICNOST SA PUZANJEM 2/D ELEMENT              (18.04.1994)
C
C    SUBROUTINE D2M15
C               TI215
C               CEC215
C
C=======================================================================
      SUBROUTINE D2M15(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 4*IDVA
      LDEFCT=LDEFT  + 4*IDVA
      LOPLUT=LDEFCT + 4*IDVA
      LOMINT=LOPLUT + 4*IDVA
      LTEFT =LOMINT + 4*IDVA
      LDQCT =LTEFT  + 1*IDVA
      LIOR  =LDQCT  + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
      LTEQT =LA2CT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 4*IDVA
      LDEFC1=LDEF1  + 4*IDVA
      LOPLU1=LDEFC1 + 4*IDVA
      LOMIN1=LOPLU1 + 4*IDVA
      LTEF1 =LOMIN1 + 4*IDVA
      LDQC1 =LTEF1  + 1*IDVA
      LIOR1 =LDQC1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
      LTEQT1=LA2C1  + 1*IDVA
C
      CALL TI215(A(LIOR),A(LDEFCT),A(LOPLUT),A(LOMINT),A(LTEFT),
     &            A(LDQCT),A(LTEQT),A(LTEQT1),
     &            A(LIOR1),A(LTAU1),A(LDEF1),A(LDEFC1),A(LOPLU1),
     &            A(LOMIN1),A(LTEF1),A(LDQC1),
     &            A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),
     &            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
      RETURN
      END
C=======================================================================
      SUBROUTINE TI215(ORI ,DEFCT,OPLUT,OMINT,TEFT ,DEFQCT,TEMT ,TEM  ,
     1                 ORI1,TAU1 ,DEF1 ,DEFC ,OPLUS,OMINS ,TEF  ,DEFQC,
     &                 A2CT,A2C  ,PTAUT,PTAU ,
     1                 FUN ,MATE ,TAU  ,DEF  ,TGT  ,TREF  ,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTICAN MATERIJAL SA PUZANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      DIMENSION DEFCT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFC(*),SE(4),
     1          OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),ETHERM(4)
      DIMENSION FUN(2,MATE*3,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR ZA KONTROLNU STAMPU
C
      KSS=4
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) KSS=3
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
C
CE  E,V,ALFA,TEMP0
C
      MAT3=(MAT-1)*3
      MATE3=MATE*3
      DO 6 J=1,3
        NFE=MAT3+J
        CALL BTAB(FUN,NTFUN,NFE,MATE3,TGT,NL,IND,2)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
C
CE    AUXILIARY CONSTANTS
C
      CZZ   =1.-V
      CNI   =CZZ-V
      CM    =E/CNI
      CNI   =CNI/3./CZZ
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2INV =(1.+V)/E
      G2C1  =C1/G2INV
      G2CNI =CNI/G2INV
      DVT   =2.D0/3.
      TRIPO =3.D0/2.
C
CE    ELASTIC CONSTUTIVE MATRIX
C
      CALL MEL2T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM2(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFCT(3)=0.5*DEFCT(3)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFCT(1)-DEFCT(2)-2*ETH)+ETH
          DEFDS(1)= C1 *(DEF(1)-DEFCT(1))-CNI*(DEF(2)-DEFCT(2))-CTH*ETH
          DEFDS(2)=-CNI*(DEF(1)-DEFCT(1))+C1 *(DEF(2)-DEFCT(2))-CTH*ETH
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT-DEFCT(1)
          DEFDS(2)=DEF(2)-EMT-DEFCT(2)
          DEFDS(4)=DEF(4)-EMT-DEFCT(4)
        END IF
        DEFDS(3)=0.5*DEF(3)-DEFCT(3)
      ELSE
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2))
          DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
          DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(4)=DEF(4)-EMT
        END IF
        DEFDS(3)=0.5*DEF(3)
      END IF
      IF (KSS.NE.3) D2=1.5*TENDO2(DEFDS)
C
CE    ELASTIC DEVIATORIC STRESS SOLUTION  (SE)
C
      DO 40 I=1,KSS
   40   SE(I) =DEFDS(I)/G2INV
      IF (KSS.EQ.3) SE(4) =-SE(1)-SE(2)
      TEF=DSQRT(1.5*TENDO2(SE))
C
      LPU=0
      IF (TEF.LT.1.D-12) THEN
        LPU=-1
        DEFQC=DEFQCT
        DO 600 I=1,4
          TAUD(I)=SE(I)
  600     DEFC(I)=DEFCT(I)
        GO TO 500
      END IF
C
CE   1)    OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      A2C    = A2CT
C
      IF (KSS.EQ.3) THEN
        P1      = 1.+G2C1*GAMDT
        P2      = G2CNI*GAMDT
        P12     = P1*P1-P2*P2
        S1      = P1*SE(1)+P2*SE(2)
        S2      = P2*SE(1)+P1*SE(2)
        TAUD(1) = S1/P12
        TAUD(2) = S2/P12
        TAUD(4) = -TAUD(1)-TAUD(2)
        TAUD(3) = DEFDS(3)/AGAMA
        ETEF    = DSQRT(1.5*TENDO2(TAUD))
        FL      = ETEF-TEFL
      ELSE
        FL      = D2-(AGAMA*TEFL)**2
      END IF
C
      IF (ITER.EQ.0) GO TO 50
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
      IF (KSS.EQ.3) THEN
        P1      = 1.+G2C1*GAMDT
        P2      = G2CNI*GAMDT
        P12     = P1*P1-P2*P2
        S1      = P1*SE(1)+P2*SE(2)
        S2      = P2*SE(1)+P1*SE(2)
        TAUD(1) = S1/P12
        TAUD(2) = S2/P12
        TAUD(4) = -TAUD(1)-TAUD(2)
        TAUD(3) = DEFDS(3)/AGAMA
        ETEF    = DSQRT(1.5*TENDO2(TAUD))
        F       = ETEF-TEF
      ELSE
        F       = D2-(AGAMA*TEF)**2
      END IF
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   2)   DEVIATORIC STRESS AND CREEP STRAIN
C
   50 IF (KSS.NE.3) THEN
        DO 165 I=1,4
  165     TAUD(I) = DEFDS(I)/AGAMA
      END IF
      DO 170 I=1,4
  170   DEFC(I) = DEFCT(I)+GAMDT*TAUD(I)
C
      IF (KSS.EQ.3) THEN
        IF (IATYP.LT.4) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+ETH
        ELSE
          EMT = -TAUD(4)/CM
        END IF
      END IF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF (ISKNP.NE.2) THEN
        IF (IATYP.LT.4) THEN
          DEFDS(1)=DEF(1)-EMT-DEFCT(1)
          DEFDS(2)=DEF(2)-EMT-DEFCT(2)
          DEFDS(4)=CZZ*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+DEFC(4)+
     &             ETH-EMT-DEFCT(4)
        ELSE
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(4)=CZZ*(DEF(1)+DEF(2))-EMT
        END IF
        CALL CEC215(AGAMA,A2C,TEF,CM)
      END IF
C
CE   3)    CALCULATE STRESS
C
  500 CONTINUE
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)/3.D0
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=TAUD(3)
      IF (KSS.EQ.3) THEN
        TAU(4)=0.D0
      ELSE
        TAU(4)=TAUD(4)+TAUM
      END IF
      DEFC(3) =2.*DEFC(3)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,4
          OPLUS(I)=OPLUT(I)
  410     OMINS(I)=OMINT(I)
        CALL ORNL2(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1=IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+DEFC(4)+ETH
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      ELSE
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2))
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      END IF
      DO 290 I=1,4
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI215')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI215'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI215 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI215')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI215'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI215 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE CEC215(AGAMA,A2C,TEF,CM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEC ( ELAST )
CE     ELASTO-CREEP  CEC MATRIX
C
C      ZA PLANE STRESS "CEC" KORISTENA STATICKA KONDENZACIJA
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION CP(4,4)
C
      DVA   =2.D0
      TR    =1.D0/3.D0
      A2CDT =A2C*DT
C
      AA=1.D0/AGAMA
      IF (TEF.LT.1.D-8) TEF=1.D-8
      DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
      DO 25 I=1,4
        DO 20 J=I,4
          CP(I,J)=-DP*DEFDS(I)*DEFDS(J)
   20   CONTINUE
        CP(I,I)=CP(I,I)+AA
   25 CONTINUE
      ELAST(1,3)=CP(1,3)
      ELAST(2,3)=CP(2,3)
      ELAST(3,3)=CP(3,3)-AA/DVA
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,4)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,4)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,4)+CM)
C
      IF(IETYP.NE.2) THEN
        ELAST(1,4)=TR*(DVA*CP(1,4)-CP(1,1)-CP(1,2)+CM)
        ELAST(2,4)=TR*(DVA*CP(2,4)-CP(1,2)-CP(2,2)+CM)
        ELAST(3,4)=CP(3,4)
        ELAST(4,4)=TR*(DVA*CP(4,4)-CP(1,4)-CP(2,4)+CM)
      END IF
C
CE  STATIC CONDESATION
C
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        DO 60 I=1,3
          DO 60 J=I,3
            ELAST(I,J)=ELAST(I,J)-ELAST(I,4)*ELAST(J,4)/ELAST(4,4)
  60    CONTINUE
        DO 61 I=1,4
  61      ELAST(I,4)=0.D0
      END IF
C
      DO 50 I=1,4
        DO 50 J=I,4
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END
C=======================================================================
C
C   TERMO-ELASTO-PLASTICNOST SA PUZANJEM -  2/D ELEMENT     (21.04.1994)
C               (IZOTROPAN SA MESOVITIM OJACANJEM)
C
C    SUBROUTINE D2M16
C               TI216
C               EPC216
C
C=======================================================================
      SUBROUTINE D2M16(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 4*IDVA
      LDEFPT=LDEFT  + 4*IDVA
      LALFAT=LDEFPT + 4*IDVA
      LDEFCT=LALFAT + 4*IDVA
      LOPLUT=LDEFCT + 4*IDVA
      LOMINT=LOPLUT + 4*IDVA
      LTEQT =LOMINT + 4*IDVA
      LTEFT =LTEQT  + 1*IDVA
      LDQPT =LTEFT  + 1*IDVA
      LDQCT =LDQPT  + 1*IDVA
      LIPL  =LDQCT  + 1*IDVA
      LIOR  =LIPL   + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 4*IDVA
      LDEFP1=LDEF1  + 4*IDVA
      LALFA1=LDEFP1 + 4*IDVA
      LDEFC1=LALFA1 + 4*IDVA
      LOPLU1=LDEFC1 + 4*IDVA
      LOMIN1=LOPLU1 + 4*IDVA
      LTEQ1 =LOMIN1 + 4*IDVA
      LTEF1 =LTEQ1  + 1*IDVA
      LDQP1 =LTEF1  + 1*IDVA
      LDQC1 =LDQP1  + 1*IDVA
      LIPL1 =LDQC1  + 1*IDVA
      LIOR1 =LIPL1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
C
      CALL TI216(A(LIOR) ,A(LDEFCT),A(LOPLUT),A(LOMINT),A(LDQCT),
     &            A(LIPL) ,A(LDEFPT),A(LALFAT),A(LDQPT) ,A(LTEQT),
     &            A(LIOR1),A(LDEFC1),A(LOPLU1),A(LOMIN1),A(LDQC1),
     &            A(LIPL1),A(LDEFP1),A(LALFA1),A(LDQP1) ,A(LTEQ1),
     &            A(LTAU1),A(LDEF1) ,A(LTEFT) ,A(LTEF1) ,
     &            A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),
     &            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI216( ORI ,DEFCT,OPLUT,OMINT ,DEFQCT,
     &                  PL  ,DEFPT,ALFAT,DEFQPT,TEMT  ,
     &                  ORI1,DEFC ,OPLUS,OMINS ,DEFQC ,
     &                  PL1 ,DEFP ,ALFA1,DEFQP ,TEM   ,
     &                  TAU1,DEF1 ,TEFT ,TEF   ,
     &                  A2CT,A2C  ,PTAUT,PTAU  ,
     &                  FUN ,MATE ,TAU  ,DEF   ,TGT   ,TREF,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTICAN MATERIJAL SA PUZANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      DIMENSION TAU(*),DEF(*),TAU1(*),DEF1(*),ETHERM(4),
     &          DEFCT(*),DEFC(*),OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),
     &          DEFPT(*),DEFP(*),ALFAT(*),ALFA1(*),SHET(4) ,GLITL(4),
     &          AHET(4)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR ZA KONTROLNU STAMPU
C
      KSS=4
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) KSS=3
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
      IPL =PL
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      CZZ   =1.-V
      CNI   =CZZ-V
      CM    =E/CNI
      CNI   =CNI/3./CZZ
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2INV =(1.+V)/E
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
      DVT   =2.D0/3.
      TRIPO =3.D0/2.
C
CE    YIELD STRESS
C
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      CALL MEL2T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM2(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFPT(3)=0.5*DEFPT(3)
      DEFCT(3)=0.5*DEFCT(3)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2)-DEFCT(1)-DEFCT(2)-
     &          2*ETH)+ETH
          DEFDS(1)=  C1*(DEF(1)-DEFPT(1)-DEFCT(1))
     &             -CNI*(DEF(2)-DEFPT(2)-DEFCT(2))-CTH*ETH
          DEFDS(2)=-CNI*(DEF(1)-DEFPT(1)-DEFCT(1))
     &             + C1*(DEF(2)-DEFPT(2)-DEFCT(2))-CTH*ETH
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT-DEFPT(1)-DEFCT(1)
          DEFDS(2)=DEF(2)-EMT-DEFPT(2)-DEFCT(2)
          DEFDS(4)=DEF(4)-EMT-DEFPT(4)-DEFCT(4)
        END IF
        DEFDS(3)=0.5*DEF(3)-DEFPT(3)-DEFCT(3)
      ELSE
        IF (KSS.EQ.3) THEN
          EMT = CNI*(DEF(1)+DEF(2))
          DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
          DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(4)=DEF(4)-EMT
        END IF
        DEFDS(3)=0.5*DEF(3)
      END IF
      D0=0.D0
      IF (KSS.NE.3) D0=1.5*TENDO2(DEFDS)
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      AGAMA = G2INV
      DO 40 I=1,KSS
        TAUD(I) =DEFDS(I)/AGAMA
        SHET(I) =TAUD(I)-ALFAT(I)
   40   GLITL(I)=SHET(I)*AGAMA
      IF (KSS.EQ.3) THEN
        TAUD(4) =-TAUD(1)-TAUD(2)
        SHET(4) =-SHET(1)-SHET(2)
        GLITL(4)=SHET(4)*AGAMA
      END IF
      TEQ=DSQRT(TRIPO*TENDO2(SHET))
      TEQE=TEQ
      TEF=DSQRT(TRIPO*TENDO2(TAUD))
C
CE   2)  CHECK FOR YIELDING
C
      LPU=0
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5.OR.ITER.EQ.0) THEN
        IPL1=0
        DEFQP=DEFQPT
        IF (TEF.LT.1.D-10) THEN
          LPU=-1
          DEFQC=DEFQCT
          DO 601 I=1,4
  601       DEFC(I)=DEFCT(I)
          GO TO 500
        END IF
      ELSE
        PL1 =1.0D0
        IPL1=1
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC WITH CREEP.
C
CE         OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      DQMIN = 0.D0
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
      D2     = D0
      A2C    = A2CT
C
      IF (KSS.EQ.3) THEN
        QGAM    = G2INV+C1*GAMDT
        CGAM    = CNI*GAMDT
        QC2     = QGAM*QGAM-CGAM*CGAM
        S1      = QGAM*DEFDS(1)+CGAM*DEFDS(2)
        S2      = CGAM*DEFDS(1)+QGAM*DEFDS(2)
        TAUD(1) = S1/QC2
        TAUD(2) = S2/QC2
        TAUD(4) = -TAUD(1)-TAUD(2)
        TAUD(3) = DEFDS(3)/AGAMA
        ETEF    = DSQRT(1.5*TENDO2(TAUD))
        FL      = ETEF-TEFL
      ELSE
        FL      = D2-(AGLHET*TEFL)**2
      END IF
C
      IF (ITER.EQ.0) GO TO 50
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
CE    FINDING THE RADIUS OF THE YIELD SURFACE (TEQ)
C
        DLAM  = 0.D0
        DEFQP = DEFQPT
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        IF (DEFQP.LT.1.D-8) DEFQP=1.D-8
        EALFA = ANCY*DEFQP**AN1
        IF (EALFA.LT.1.D-8) EALFA=1.D-8
C
        KB    = 0
        KT    = 0
        DEPL  = 0.D0
        DDEP  = 0.1*(TEQE-TEQY)/EALFA
C
        IF (KSS.EQ.3) THEN
          QGAM    = G2INV+C1*GAMDT
          CGAM    = CNI*GAMDT
          QC2     = QGAM*QGAM-CGAM*CGAM
          GLITL(1)= DEFDS(1)-QGAM*ALFAT(1)+CGAM*ALFAT(2)
          GLITL(2)= DEFDS(2)+CGAM*ALFAT(1)-QGAM*ALFAT(2)
          GLITL(3)= DEFDS(3)-ALFAT(3)*AGAMA
          S1H     = QGAM*GLITL(1)+CGAM*GLITL(2)
          S2H     = CGAM*GLITL(1)+QGAM*GLITL(2)
          SHET(1) = S1H/QC2
          SHET(2) = S2H/QC2
          SHET(4) = -SHET(1)-SHET(2)
          SHET(3) = GLITL(3)/AGLHET
          TEQ     = DSQRT(1.5*TENDO2(SHET))
          FHETL   = TEQ-TEQY
        ELSE
          DO 240 I=1,4
  240       GLITL(I) = DEFDS(I)-ALFAT(I)*AGAMA
          GNORM = DSQRT(TRIPO*TENDO2(GLITL))
          FHETL   = GNORM-AGAMA*TEQY
        END IF
C
        IF (FHETL.LE.0.D0) THEN
          IPL1=0
          GO TO 260
        END IF
        IPL1=1
        DDEFQP= DDEP
C
  250   KT    = KT+1
        KB1   = KB
C
        IF(KT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2015)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6015)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        DEFQP = DEFQPT+DDEFQP
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        EALFA = ANCY*DEFQP**AN1
        CHET  = DVT*EM1*EALFA
        AGPOT = 1.+AGAMA*CHET
C
        IF (KSS.EQ.3) THEN
          DLAM    = TRIPO*DDEFQP/TEQY
          AGLHET  = AGAMA+AGPOT*DLAM
          P1      = QGAM+(C1+QGAM*CHET)*DLAM
          P2      = CGAM+(CNI+CGAM*CHET)*DLAM
          P12     = P1*P1-P2*P2
          S1H     = P1*GLITL(1)+P2*GLITL(2)
          S2H     = P2*GLITL(1)+P1*GLITL(2)
          SHET(1) = S1H/P12
          SHET(2) = S2H/P12
          SHET(4) = -SHET(1)-SHET(2)
          SHET(3) = GLITL(3)/AGLHET
          TEQ     = DSQRT(1.5*TENDO2(SHET))
          FHET    = TEQ-TEQY
        ELSE
          FHET    = GNORM-AGAMA*TEQY-TRIPO*AGPOT*DDEFQP
        END IF
C
        CALL BISECB(DDEFQP,DEPL,DEPD,DDEP,
     &              FHET,FHETL,FHETD,AF,KB,DQMIN,DQTOL)
        IF (KB.EQ.-1) GO TO 255
        IF (KB1.EQ.0) GO TO 250
      IF (DABS(DDEP).GT.EPSIL.AND.
     1      (DABS(DDEP)/(DEPL+DEPD)).GT.EPSIL) GO TO 250
C
  255   IF (KSS.NE.3) THEN
          DLAM = TRIPO*DDEFQP/TEQY
          DO 251 I=1,4
  251       AHET(I) = ALFAT(I)+CHET*DEFDS(I)
          AGLHET = AGAMA+AGPOT*DLAM
          D2 = D0+3*DLAM*TDO2AB(DEFDS,AHET)+TRIPO*DLAM*DLAM*TENDO2(AHET)
        END IF
C
  260   IF (KSS.EQ.3) THEN
          DX      = DEFDS(1)-C1*DLAM*SHET(1)+CNI*DLAM*SHET(2)
          DY      = DEFDS(2)+CNI*DLAM*SHET(1)-C1*DLAM*SHET(2)
          S1      = QGAM*DX+CGAM*DY
          S2      = CGAM*DX+QGAM*DY
          TAUD(1) = S1/QC2
          TAUD(2) = S2/QC2
          TAUD(4) = -TAUD(1)-TAUD(2)
          TAUD(3) = (DEFDS(3)-DLAM*SHET(3))/AGAMA
          ETEF    = DSQRT(1.5*TENDO2(TAUD))
          F       = ETEF-TEF
        ELSE
          F       = D2-(AGLHET*TEF)**2
        END IF
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   4)   DEVIATORIC STRESS,PLASTIC STRAIN,BACK STRESS AND CREEP STRAIN
C
   50 DO 165 I=1,4
        IF (KSS.NE.3) TAUD(I) = DEFDS(I)/AGAMA
        IF ((IPL1.EQ.1).AND.(ITER.NE.0)) THEN
          IF (KSS.NE.3) SHET(I) = GLITL(I)/AGLHET
          DDEFPS   = DLAM*SHET(I)
          DEFP(I)  = DEFPT(I)+DDEFPS
          ALFA1(I) = ALFAT(I)+CHET*DDEFPS
          IF (KSS.NE.3) TAUD(I) = TAUD(I)-DDEFPS/AGAMA
        ELSE
          DEFP(I)  = DEFPT(I)
        END IF
        DEFC(I) = DEFCT(I)+GAMDT*TAUD(I)
  165 CONTINUE
C
      IF (KSS.EQ.3) THEN
        IF (IATYP.LT.4) THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &          2*ETH)+ETH
        ELSE
          EMT = -TAUD(4)/CM
        END IF
      END IF
C
CE     E L A S T I C  -  P L A S T I C  -  C R E E P   M A T R I X  CEPC
C
      IF (ISKNP.NE.2) THEN
        IF (KSS.EQ.3) THEN
          IF (IATYP.LT.4) THEN
            DEFDS(1)=DEF(1)-EMT-DEFPT(1)-DEFCT(1)
            DEFDS(2)=DEF(2)-EMT-DEFPT(2)-DEFCT(2)
            DEFDS(4)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &               2*ETH)+DEFP(4)+DEFC(4)+ETH-EMT-DEFPT(4)-DEFCT(4)
          ELSE
            DEFDS(1)=DEF(1)-EMT
            DEFDS(2)=DEF(2)-EMT
            DEFDS(4)=CZZ*(DEF(1)+DEF(2))-EMT
          END IF
          IF ((IPL1.EQ.1).AND.(ITER.NE.0)) THEN
            DO 245 I=1,4
              AHET(I)  = ALFAT(I)+CHET*DEFDS(I)
  245         GLITL(I) = DEFDS(I)-ALFAT(I)*AGAMA
            GNORM = DSQRT(TRIPO*TENDO2(GLITL))
          END IF
        END IF
        CALL EPC216(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &              AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQY,
     &              ALFAT)
      END IF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      IF (IPL1.NE.1) THEN
        DO 600 I=1,4
  600     DEFP(I)=DEFPT(I)
        TEQ=TEQE
      END IF
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)/3.D0
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=TAUD(3)
      IF (KSS.EQ.3) THEN
        TAU(4)=0.D0
      ELSE
        TAU(4)=TAUD(4)+TAUM
      END IF
      DEFP(3) =2.*DEFP(3)
      DEFC(3) =2.*DEFC(3)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,4
          OPLUS(I)=OPLUT(I)
  410     OMINS(I)=OMINT(I)
        CALL ORNL2(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1=IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF (IATYP.LT.4) THEN
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &           2*ETH)+DEFP(4)+DEFC(4)+ETH
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      ELSE
        IF (KSS.EQ.3) THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2))
        ELSE IF (IETYP.EQ.2) THEN
          DEF(4)=0.D0
        END IF
      END IF
      DO 290 I=1,4
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI216')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI216'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI216 ',
     &        '( PSEUDO-TIME )')
 2015 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI216 ',
     &        '( RADIJUS POVRSI TECENJA )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI216')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI216'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI216 ',
     &        '( PSEUDO-TIME )')
 6015 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI216 ',
     &        '( THE RADIUS OF YIELD SURFACE )')
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE EPC216(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &                  AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQ,
     &                  ALFAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE EPC ( ELAST )
CE     ELASTO-PLASTIC-CREEP  EPC MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      DIMENSION CP(4,4),AHET(*),EAHET(4),GLITL(*),ALFAT(*)
C
      DVA   =2.D0
      TR    =1.D0/3.
      A2CDT =A2C*DT
      IF (TEF.LT.1.D-8) TEF=1.D-8
C
      IF ((IPL1.NE.1).OR.(ITER.EQ.0)) THEN
        AA=1./AGAMA
        DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
C
        DO 25 I=1,4
          DO 20 J=I,4
            CP(I,J)=-DP*DEFDS(I)*DEFDS(J)
   20     CONTINUE
          CP(I,I)=CP(I,I)+AA
   25   CONTINUE
      ELSE
C
        EMEAL=(EM**AN)*EALFA
        CHETP=CHET*AN1/DEFQP
        AMGAM=AGAMA*(EMEAL+1.5*CHETP*DDEFQP)+1.5*AGPOT
        ASPOT=1.5/GNORM/AMGAM
        BGPOT=A2CDT*(TEQ+1.5*CHET*DDEFQP+1.5*TDO2AB(GLITL,ALFAT)/GNORM)
     &              /AMGAM
        DLFAK=(1.5-DLAM*EMEAL)/TEQ
        ASIG =ASPOT*DLFAK
        BGAM =BGPOT*DLFAK
        CHDL1=1.+CHET*DLAM
        BIFAK=CHDL1*A2CDT
        DIFAK=CHETP*DLAM
        AA   =CHDL1/AGLHET
        ATEF2=AGLHET*TEF*TEF
        DO 100 I=1,4
  100     EAHET(I)=DEFDS(I)+DLAM*AHET(I)
C
        ASIGMA=DVA*BIFAK*ATEF2
        BSIGMA=DVA*AGPOT*ATEF2-3.*TDO2AB(EAHET,AHET)
        CSIGMA=DLAM*(DVA*AGAMA*ATEF2-3.*TDO2AB(EAHET,DEFDS))
        DSIGMA=ASIGMA-BGAM*BSIGMA-BGPOT*CHETP*CSIGMA+2*AGLHET*AGLHET*TEF
        IF (DSIGMA.LT.1.D-8) DSIGMA=1.D-8
        HJFAK =ASIG*BSIGMA+ASPOT*CHETP*CSIGMA
C
        DO 110 I=1,4
          BI=BIFAK*TAUD(I)
          CI=AGPOT*TAUD(I)-AHET(I)
          DI=DIFAK*(AGAMA*TAUD(I)-DEFDS(I))
          DO 120 J=I,4
            HJ=3*CHDL1*EAHET(J)-HJFAK*GLITL(J)
            SJOT=HJ/DSIGMA
            CP(I,J)=(-BI*SJOT-CI*(ASIG*GLITL(J)-BGAM*SJOT)
     &              -DI*(ASPOT*GLITL(J)-BGPOT*SJOT))/AGLHET
  120     CONTINUE
          CP(I,I)=CP(I,I)+AA
  110   CONTINUE
      END IF
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,4)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,4)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,4)+CM)
      ELAST(1,3)=CP(1,3)
      ELAST(2,3)=CP(2,3)
      ELAST(3,3)=CP(3,3)-AA/DVA
C
      IF(IETYP.NE.2) THEN
        ELAST(1,4)=TR*(DVA*CP(1,4)-CP(1,1)-CP(1,2)+CM)
        ELAST(2,4)=TR*(DVA*CP(2,4)-CP(1,2)-CP(2,2)+CM)
        ELAST(3,4)=CP(3,4)
        ELAST(4,4)=TR*(DVA*CP(4,4)-CP(1,4)-CP(2,4)+CM)
      END IF
C
CE  STATIC CONDESATION
C
      IF (IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        DO 60 I=1,3
          DO 60 J=I,3
            ELAST(I,J)=ELAST(I,J)-ELAST(I,4)*ELAST(J,4)/ELAST(4,4)
  60    CONTINUE
        DO 61 I=1,4
  61      ELAST(I,4)=0.D0
      END IF
C
      DO 50 I=1,4
        DO 50 J=I,4
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END

