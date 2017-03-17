C=======================================================================
C
C   VISKOPLASTICNOST LJUSKE  -  ORTOTROPNI MATERIJAL
C                               MESOVITO OJACANJE     (24.09.1994)
C
C=======================================================================
      SUBROUTINE D8M21(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TI821 (A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI821 ( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC ANISOTROPIC MATERIAL , MIXED HARDENING
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /UGAOV6/ TE(6,6)
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),FUN(*)
      DIMENSION COEFE(2)
      DIMENSION CE(6,6),CPE(3,3),CM(3),CE2(3)
      DIMENSION SS(6),SSH(6),EVPT(6),DEVP(6),
     1          CP(6,6),AMAT(3,3),AV(3),BV(6),Y0(6),CY(6),
     2          AN(6),AN1(6),E(6),DEFL(6),WJ(6),CPHM(3,3),
     2          TAUL(6),EN(3,3),B(3,3),CPH(6),
     2          SE(6),BVP(6),EOLD(3),DE(3),
     3          DEN(3,3),ENOLD(3,3),BVOLD(3),DBV(3)
      DATA MAXBIS/100/,EPSIL/1.0D-10/
      COEFE(1)=0.83333333D0
      COEFE(2)=COEFE(1)
CE  INITIAL DATA
C
C.. CONSTANTS
      KDIM=6
      IPL =PL
      IPL1=PL1
      DVA =2.D0
      DVT =DVA/3.
      SQRT3 = SQRT(3.D0)
      ZER = 0.D0
      IF (IRAC.EQ.2) RETURN
      TOL0 = 1.D-9
      TOLDX = 1.D-10
      TOLF = 1.D-10
      FACTVP = 1.D-5
C
C     FORM ELASTIC MATRICES AND YIELD CURVES
      CALL       VP810(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,
     &                  COEFE,ETP,BETA,TE,ETA)
C
      CALL  ENTY(TAUY,E,Y0,CY,AN,EM,DEFQPT,DEFQPT,KDIM)
C
      CE2(1) = 2.*CE(4,4)
      CE2(2) = 2.*CE(5,5)
      CE2(3) = 2.*CE(6,6)
      TEQY=TEQY0
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
      ANQ1 = ANQ - 1.
      EM1   =1.-EM
      EM1DVT = EM1*DVT
      IMIX  = 1
      IF(DABS(EM1).LT.EPSIL) IMIX=0
C
C... TRANSFORM STRAIN INTO MATERIAL AXES DIRECTION
C
      CALL CLEAR(DEFL,6)
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL MNOZI1(DEFL,TE,DEF,6,6)
      ELSE
         CALL JEDNA1(DEFL,DEF,6)
      ENDIF
      IF(IRAC.EQ.2)THEN
C...................OVO ZAMENITI JER OVAKO JE UVEK U ITERACIJI 0 ELASTICNA
C...................TREBA RACUNATI C ELASTOPLASTIC.
C****** OVO NEMA VEZE !
C******        CALL JEDNA1(ETP,CE,36)
        RETURN
      ENDIF
C
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=1,3
      I3 = I + 3
      EVPT(I) = DEFPT(I)
   10 EVPT(I3)=0.5*DEFPT(I3)
      ETADT = ETA/DT
      ETAD15 = 1.5*ETADT
      IVP = 0
      DDEFQP = ZER
      DO 8 I = 1,6
    8 CPH(I) = ZER
      DO 9 I = 1,3
      DO 9 J = 1,3
    9 CPHM(I,J) = ZER
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.NE.4)THEN
        DO 20 I=1,3
   20   DEFDS(I)=DEFL(I)-EVPT(I)
        DO 25 I=4,6
   25   DEFDS(I)=0.5*DEFL(I)-EVPT(I)
      ELSE
        DO 21 I=1,3
   21   DEFDS(I)=DEFL(I)
        DO 26 I=4,6
   26   DEFDS(I)=0.5*DEFL(I)
      ENDIF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      DO 40 I=1,2
      DUM=0.D0
      DO 30 J=1,2
   30 DUM=DUM+CPE(I,J)*DEFDS(J)
      TAUD(I) =DUM
   40 SE(I) =DUM-ALFAT(I)
      TAUD(3) = - TAUD(1) - TAUD(2)
      SE(3)=-SE(1)-SE(2)
      DO 45 I=1,3
      I3 = I+3
      TAUD(I3) = CE2(I)*DEFDS(I3)
   45 SE(I3)  =TAUD(I3)-ALFAT(I3)
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5*TDOTAN(SE,E))
      DEPMX = (TEQ-TEQY)/ETADT
      IF ((TEQ-TEQY).LT.0.D0. OR. DEPMX.LT.TOL0) THEN
        DEFQP=DEFQPT
        TAUY = TEQY
        CALL JEDNA1(DEFP,DEFPT,6)
        GO TO 500
      ENDIF
C
CE   3)  SOLUTION IS VICOPLASTIC.  OBTAIN ZERO OF THE GOVER. FN.
C
      IVP = 1
      PL1=1.0D0
      ETADT = ETA/DT
      ETAD15 = 1.5*ETADT
C
CE    BISECTION ON DDEFQP   
C     
      DEFQP=DEFQPT
      DUM=DEFQP
      IF(DEFQP.LE.1.D-4) DUM=1.D-4
C***************  OVO TREBA DA BUDE EKVIVALENTNI MODUL
      EP2=ANQ*CYQ*DUM**ANQ1 
C
      IB = 0
      IT = 0
      AF = 3.D0
      DDD = FACTVP*(TEQ-TEQY)/ETADT
      FP    = TEQ - TEQT
      TAUY  = TEQT
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      DDEFQP= DDD
      TOLD= TAUY
      DEFQP=DEFQPT+DDEFQP
      DOLD= DEFQPT
      EP  = EP2
      DO 24 I = 1,3
      DBV(I) = 0.D0
      DE(I) = 0.D0
      DO 24 J = 1,3
      DEN(I,J) = 0.D0
   24 CONTINUE
      DBH11 = 0.D0
      DBH12 = 0.D0
      DBH21 = 0.D0
      DBH22 = 0.D0
C
      IDER = 0
  100 IT=IT+1
      IB1 = IB
C
      DEFQP=DEFQPT+DDEFQP
C
      IF (IT.GT.1) CALL  ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      AK=ETAD15*DDEFQP/TAUY
      EN(1,1) = E(1) + E(2)
      EN(1,2) = - E(1)
      EN(1,3) = - E(2)
      EN(2,2) = E(1) + E(3)
      EN(2,3) = - E(3)
      EN(3,3) = E(2) + E(3)
      EN(2,1) = EN(1,2)
      EN(3,1) = EN(1,3)
      EN(3,2) = EN(2,3)
      IF (IMIX.EQ.1)THEN
         CALL XMXMAT(CPH,DEFQP,EM1DVT,CY,AN,AN1,KDIM)
         DO 27 I = 1,3
         DO 27 J = 1,3
   27    CPHM(I,J) = CPH(J)*EN(I,J)
         CPH(4) = SQRT3*E(4)*CPH(4)
         CPH(5) = SQRT3*E(5)*CPH(5)
         CPH(6) = SQRT3*E(6)*CPH(6) 
      ENDIF
      DO 70 I = 1,3
      DO 70 J = 1,3
      DUM = 0.D0
      DO 68 K = 1,3
   68 DUM = DUM + (CPHM(I,K)+CPE(I,K))*EN(K,J)
      B(I,J) = EN(I,J) + DUM/ETADT
   70 CONTINUE
      BH11 = B(2,2) - B(2,3)
      BH12 = B(1,2) - B(1,3)
      BH21 = B(2,1) - B(2,3)
      BH22 = B(1,1) - B(1,3)
      SH11 = 1. + AK*BH11
      SH12 = AK*BH12
      SH21 = AK*BH21
      SH22 = 1. + AK*BH22
      DK = SH11*SH22 - SH12*SH21
      SSH(1) = (SH11*SE(1) - SH12*SE(2))/DK
      SSH(2) = (-SH21*SE(1) + SH22*SE(2))/DK
      SSH(3) = - SSH(1) - SSH(2)
      DO 22 I = 1,3
      I3 = I + 3
      BV(I3) = (1. + (CPH(I3)+CE2(I))/ETADT)*E(I3)
      SSH(I3) = SE(I3)/(1.+AK*BV(I3))
   22 CONTINUE
C
      TEQ=DSQRT(1.5*TDOTAN(SSH,E))
      FB = TEQ-TAUY
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
C
      DDF   = DEFQP-DOLD
      DOLD = DEFQP
      IF (IB.EQ.0)  GO TO 33
      DDFA = ABS(DDF)
      IF  (DDFA.GE.1.D-10) THEN
          IF (IDER.EQ.1) GO TO 35
   32     DO 34 I = 1,3
          I3 = I + 3
          BVOLD(I) = BV(I3)
          EOLD(I) = E(I3)
          DO 34 J = 1,3
   34     ENOLD(I,J) = EN(I,J)
          BH11OL = BH11
          BH12OL = BH12
          BH21OL = BH21
          BH22OL = BH22
          TOLD= TAUY
          IDER = 1
          GO TO 33
   35     EP  = (TAUY-TOLD)/DDF
          TOLD = TAUY
          DO 36 I = 1,3
          I3 = I + 3
          DBV(I) = (BV(I3)-BVOLD(I))/DDF
          DE(I) = (E(I3)-EOLD(I))/DDF
          DO 36 J = I,3
          DEN(I,J) = (EN(I,J)-ENOLD(I,J))/DDF
          DEN(J,I) = DEN(I,J) 
   36     CONTINUE
          DBH11 = (BH11-BH11OL)/DDF
          DBH12 = (BH12-BH12OL)/DDF
          DBH21 = (BH21-BH21OL)/DDF
          DBH22 = (BH22-BH22OL)/DDF
          GO TO 32    
      ENDIF
C
   33 IF  (I.GT.MAXBIS) THEN                                                  
C                                                                               
C         NUMBER OF TRIALS FOR DDDEFQP EXCEEDS MAXBIS - STOP                      
C                                                                               
          WRITE (IZLAZ,1030) I,DEPBM,DEPBP,DDEFQP,FB,FM,FP
          STOP                                                                  
      ENDIF
      IF (IB1.EQ.0) GO TO 100
      IF (DEPBM.LT.TOL0) GO TO 162
      IF (ABS(DDD).LT.TOLDX .AND. ABS(FB).LT.TOLF) GO TO 162
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE   4)  DETERMINE SOLUTION 
C
  162 COEF = AK/ETADT
      DEFQPT = DEFQP
      TEQT = TAUY
      TEQY = TAUY
      DO 200 I = 1,2
      DUM = 0.D0
      DO 198 K = 1,3
  198 DUM = DUM + EN(I,K)*SSH(K)
      DEVP(I) = COEF*DUM
      DEFP(I) = DEFPT(I) + DEVP(I)
  200 CONTINUE
      DEVP(3) = - DEVP(1) - DEVP(2)
      DEFP(3) = - DEFP(1) - DEFP(2)
      DO 201 I = 4,6
      DEVP(I) = COEF*E(I)*SSH(I)
      DEFP(I) = DEFPT(I) + DEVP(I) + DEVP(I) 
  201 CONTINUE
      IF  (IMIX.EQ.1) THEN
          DO 205 I = 1,2
          DUM = 0.D0
          DO 204 J = 1,3
  204     DUM = DUM + CPHM(I,J)*DEVP(J)
          ALFA1(I) = ALFAT(I) + DUM
  205     CONTINUE
          ALFA1(3) = - ALFA1(1) - ALFA1(2)
          DO 206 I = 4,6
          ALFA1(I) = ALFAT(I) + CPH(I)*DEVP(I)
  206     CONTINUE
      ENDIF
      DO 207 I = 1,6
      SS(I) = ALFA1(I) + SSH(I)
  207 TAUD(I) = SS(I) + DEVP(I)*ETADT
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      TAUM = CM(1)*(DEFL(1)-DEFPT(1)-DEVP(1)) +
     1       CM(2)*(DEFL(2)-DEFPT(2)-DEVP(2))
      DO 505 I=1,2
  505 TAUL(I)=TAUD(I)+TAUM
      TAUL(3) = 0.D0
      TAUL(4)=TAUD(4)
      TAUL(5)=TAUD(5)
      TAUL(6)=TAUD(6)
C
C...  TRANSFORM STRESS INTO GLOBAL AXES DIRECTION
C
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL CLEAR(TAU,6)
         CALL MNOZI2(TAU,TE,TAUL,6,6)
      ELSE
         CALL JEDNA1(TAU,TAUL,6)
      ENDIF
C
      IF(ISKNP.EQ.2) GO TO 700
C
C     CONSTITUTIVE MATRIX
C
      IF (IVP.EQ.1) GO TO 550
C
C     ELASTIC MATRIX
C
      CALL JEDNA1(ELAST,CE,36)
      GO TO 570
C
C     ELASTIC - VISCOPLASTIC MATRIX
C
  550 C11 = SH11/DK
      C12 = SH12/DK
      C21 = SH21/DK
      C22 = SH22/DK
      DO 510 J = 1,3
      AMAT(1,J) = C11*CPE(1,J) - C12*CPE(2,J)
      AMAT(2,J) = -C21*CPE(1,J) + C22*CPE(2,J)
      AMAT(3,J) = -AMAT(1,J) - AMAT(2,J)
  510 CONTINUE
      A1 = (1.5*ETADT - AK*EP)/TAUY
      A3 = BH11 + BH22 + 2.*AK*(BH11*BH22-BH12*BH21)
      A4 = AK*((1.+AK*BH22)*DBH11 + (1.+AK*BH11)*DBH22 -
     1     AK*(BH12*DBH21+DBH12*BH21))
      A1DK = A1/DK
      BV(1) = A1DK*(BH11*SE(1) - BH12*SE(2) - A3*SSH(1))
      BV(1) = BV(1) + (-A4*SSH(1) + AK*(SE(1)*DBH11-SE(2)*DBH12))/DK
      BV(2) = A1DK*(-BH21*SE(1) + BH22*SE(2) - A3*SSH(2))
      BV(2) = BV(2) + (-A4*SSH(2) + AK*(-SE(1)*DBH21+SE(2)*DBH22))/DK
      BV(3) = -BV(1) - BV(2)
      DO 512 I = 1,3
      I3 = I + 3
      DENA = 1. + AK*BV(I3)
      AV(I) = CE2(I)/DENA
      BV(I3) = -(A1*BV(I3)+AK*DBV(I))*SSH(I3)/DENA 
  512 CONTINUE
      DO 515 I = 1,3
      I3 = I + 3
      CP(1,I3) = SSH(I3)*E(I3)
      CP(2,I3) = DE(I)*SSH(I3)
      CP(1,I) = 0.D0
      CP(2,I) = 0.D0
      DO 514 J = 1,3
      CP(2,I) = CP(2,I) + DEN(J,I)*SSH(J)
  514 CP(1,I) = CP(1,I) + EN(J,I)*SSH(J)
      CP(2,I) = 0.5*CP(2,I)
  515 CONTINUE
C
      W0 = 0.D0
      DO 520 I = 1,3
      I3 = I + 3
      WJ(I) = 0.D0
      DO 518 K = 1,3
  518 WJ(I) = WJ(I) + CP(1,K)*AMAT(K,I)
      W0 = W0 + CP(1,I)*BV(I) + 2.*CP(1,I3)*BV(I3) +
     1          CP(2,I)*SSH(I) + CP(2,I3)*SSH(I3)
      WJ(I3) = 2.*CP(1,I3)*AV(I)
  520 CONTINUE
      W0 = W0 + DVT*TAUY*EP
      DO 522 I = 1,6
  522 WJ(I) = -WJ(I)/W0
      DO 530 I = 1,3
      I3 = I + 3
      BVP(I) = 0.D0
      BVP(I3) = AK*SSH(I3)*DE(I)
      DO 529 J = 1,3
      J3 = J + 3
      BVP(I) = BVP(I) + DEN(I,J)*SSH(J)
      ELAST(I,J) = AMAT(I,J) + BV(I)*WJ(J)
      ELAST(I,J3) = BV(I)*WJ(J3)
      ELAST(I3,J) = BV(I3)*WJ(J)
      ELAST(I3,J3) = BV(I3)*WJ(J3)
  529 CONTINUE
      ELAST(I3,I3) = ELAST(I3,I3) + AV(I)
      BVP(I) = AK*BVP(I)
  530 CONTINUE
      DO 535 I = 1,6
      DO 535 J = 1,6
      ELAST(I,J) = A1*SSH(I)*WJ(J) + AK*ELAST(I,J)
  535 CONTINUE
C
      DO 540 I = 1,3
      I3 = I + 3
      DO 540 J = 1,3
      J3 = J + 3
      CP(I,J) = 0.D0
      CP(I,J3) = 0.D0
      DO 536 K = 1,3
      CP(I,J) = CP(I,J) + EN(I,K)*ELAST(K,J)
  536 CP(I,J3) = CP(I,J3) + EN(I,K)*ELAST(K,J3)
      CP(I,J) = CP(I,J) + BVP(I)*WJ(J)
      CP(I,J3) = CP(I,J3) + BVP(I)*WJ(J3)
      CP(I3,J) = E(I3)*ELAST(I3,J) + BVP(I3)*WJ(J)
      CP(I3,J3) = E(I3)*ELAST(I3,J3) + BVP(I3)*WJ(J3)
  540 CONTINUE
      DO 545 I = 1,3
      I3 = I + 3
      DO 547 J = 1,3
      J3 = J + 3
      DUM = 0.D0
      DUM3 = 0.D0
      DO 544 K = 1,3
      DUM1 = CPE(I,K) + CM(K)
      DUM = DUM + DUM1*CP(K,J)
      DUM3 = DUM3 + DUM1*CP(K,J3)
  544 CONTINUE
      ELAST(I,J) = CPE(I,J) + CM(J) - DUM/ETADT
      ELAST(I,J3) = -0.5*DUM3/ETADT
      ELAST(I3,J) = - CE2(I)*CP(I3,J)/ETADT
      ELAST(I3,J3) = -0.5*CE2(I)*CP(I3,J3)/ETADT
  547 CONTINUE
      ELAST(I3,I3) = ELAST(I3,I3) + 0.5*CE2(I)
  545 CONTINUE
      DO 560 I = 1,6
      DO 560 J = I,6
      ELAST(I,J) = 0.5*(ELAST(I,J) + ELAST(J,I))
      ELAST(J,I) = ELAST(I,J)
  560 CONTINUE
C
  570 IF(DABS(BETA).GT.1.D-6) THEN
        CALL TRAETP(ELAST,ETP,TE)
      ELSE
        CALL JEDNA1(ETP,ELAST,36)
      ENDIF
C
CE  UPDATE FROM PREVIOUS STEP
C
  700 DEFL(3)=-VZX/EZ*TAUL(1)-VYZ/EY*TAUL(2)+DEFP(3)
      DO 290 I=1,6
      DEF1(I)=DEFL(I)
  290 TAU1(I)=TAU(I)
C
      RETURN
C
 1030 FORMAT (/ ' SUBROUTINE TI821'/
     1 ' NO SOLUTION IN BISECTING FOR F=0. '/
     1 ' IBIS =',I5/
     3 ' X-MINUS =',D14.6/' X-PLUS=',D14.6/' DDEFQP=',D14.6/
     4        ' F =',D14.6/' F-MINUS=',D14.6/' F-PLUS=',D14.6)
      END
C=======================================================================
      SUBROUTINE VP810(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,
     &                  COEFE,ETP,BETA,TE,ETA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORM ELASTICITY CONSTANTS - ANISOTROPIC VISCOPLASTICITY
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      DIMENSION FUN(30,*),CE(6,*),CPE(3,*),CM(*),COEFE(*),ETP(6,*)
      DIMENSION Y0(*),CY(*),AN(*),AN1(*)
C
      D13 =1.D0/3.
      DVA =2.D0
      ZER =0.D0
      EX=FUN(1,MAT)
      EY=FUN(2,MAT)
      EZ=FUN(3,MAT)
      VXY=FUN(4,MAT)
      VYZ=FUN(5,MAT)
      VZX=FUN(6,MAT)
      GXY=FUN(7,MAT)
      GYZ=FUN(8,MAT)
      GZX=FUN(9,MAT)
      DO 10 I=1,6
       Y0(I) =FUN( 9+I,MAT)
       CY(I) =FUN(15+I,MAT)
       AN(I) =FUN(21+I,MAT)
       AN1(I)=AN(I)-1.
   10 CONTINUE
      TEQY0 =FUN(28,MAT)
      EM    =FUN(29,MAT)
      ETA   =FUN(30,MAT)
      CYQ   =CY(1)
      ANQ   =AN(1)
      DO 15 I=1,6
      DO 15 J=1,6
   15 CE(I,J)=0.D0
C     MATRICA CE
      POM=EX-EY*VXY*VXY
      CE(1,1)=EX*EX/POM
      CE(2,2)=EX*EY/POM
      CE(1,2)=EX*EY*VXY/POM
      CE(4,4)=GXY
      CE(5,5)=COEFE(1)*GYZ
      CE(6,6)=COEFE(2)*GZX
      DO 50 I=1,6
      DO 50 J=I,6
   50 CE(J,I)=CE(I,J)
C...   MATRIX  C'E
      CPE(1,1)=D13*(DVA*CE(1,1)-CE(1,2))
      CPE(1,2)=D13*(DVA*CE(1,2)-CE(2,2))
      CPE(2,1)=D13*(DVA*CE(1,2)-CE(1,1))
      CPE(2,2)=D13*(DVA*CE(2,2)-CE(1,2))
        DO 60 I=1,3
        CPE(I,3)=ZER
   60   CPE(3,I)=ZER
C...   VECTOR   CM
      CM(1)=D13*(CE(1,1)+CE(1,2))
      CM(2)=D13*(CE(1,2)+CE(2,2))
      CM(3)=ZER
C
      IF(DABS(BETA).GT.1.D-6) THEN
        CALL TRAETP(CE,ETP,TE)
      ELSE
        CALL JEDNA1(ETP,CE,36)
      ENDIF
      RETURN
      END
