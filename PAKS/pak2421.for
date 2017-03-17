C=======================================================================
C
C   VISKO-PLASTICNOST 2/D ELEMENT  -  ANIZOTROPNI MATERIJAL  
C                               MESOVITO OJACANJE  2D  (4.09.1994)
C
C=======================================================================
      SUBROUTINE D2M21(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
C      MATE=MREPER(4)
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
      CALL TI221 (A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI221 ( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
C     VISKO-PLASTICAN ANIZOTROPAN MATERIJAL SA MESOVITIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /UGAOVL/ TE(4,4)
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFP(*),
     1          ALFAT(*),ALFA1(*),FUN(*)
      DIMENSION SE(4),DEFL(4),TAUL(4)
      DIMENSION DEVP(4),SS(4),SSH(4)
      DIMENSION CE(4,4),CPE(3,3),CM(3),E(6),
     1          Y0(6),CY(6),AN(6),AN1(6)
      DIMENSION 
     1          CP(4,4),AMAT(3,3),BV(3),WJ(4),
     2          EN(3,3),B(3,3),CPH(6),CPHM(3,3),
     2          BVP(3),DEN(3,3),ENOLD(3,3),IPERM(4)
C
      DATA MAXBIS/100/,EPSIL/1.0D-10/
      DATA IPERM/1,2,4,3/
CE  INITIAL DATA
C
      IPL =PL
      IPL1=PL1
      KDIM=4
      DVA =2.D0
      DVT =DVA/3.
      SQRT3 = SQRT(3.D0)
      ZER = 0.D0
      TOL0 = 1.D-9
      TOLDX = 1.D-10
      TOLF = 1.D-10
      FACTVP = 1.D-5
C
C     FORM ELASTIC MATRICES AND YIELD CURVES
      IF (IRAC.EQ.2) RETURN
C
      CALL  VP210(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,ETA)
      TAUY   =0.D0
      CALL  ENTY(TAUY,E,Y0,CY,AN,EM,DEFQPT,DEFQPT,KDIM)
C
      CEXY = 2.*CE(3,3)
      TEQY=TEQY0
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
C
C... TRANSFORM STRAIN INTO MATERIAL AXES DIRECTION
C
      CALL CLEAR(DEFL,4)
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL MNOZI1(DEFL,TE,DEF,4,4)
      ELSE
         CALL JEDNA1(DEFL,DEF,4)
      ENDIF
      IF(IRAC.EQ.2)THEN
C...................OVO ZAMENITI JER OVAKO JE UVEK U ITERACIJI 0 ELASTICNA
C...................TREBA RACUNATI C ELASTOPLASTIC.
        DO 10 I=1,4
        DO 10 J=1,4
   10   ELAST(I,J)=CE(I,J)
CS.... TRANSFORMACIJA MATRICE  ELAST()
CE     TRANSFORM ELAST MATRIX
      IF(DABS(BETA).GT.1.0D-10.AND.IETYP.NE.1) 
     1CALL TRAEL(ELAST,TE,4,3,3,ELAST)
      IF(DABS(BETA).GT.1.0D-10.AND.IETYP.EQ.1) 
     1CALL TRAEL(ELAST,TE,4,4,4,ELAST)
        RETURN
      ENDIF
C
      ANQ1 = ANQ - 1.
      EM1   =1.-EM
      EM1DVT = EM1*DVT
      IMIX  = 1
      IF(DABS(EM1).LT.EPSIL) IMIX=0
C
      IVP = 0
      DDEFQP = ZER
      ETADT = ETA/DT
      DO 8 I = 1,4
    8 CPH(I) = ZER
      DO 9 I = 1,3
      DO 9 J = 1,3
    9 CPHM(I,J) = ZER
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.NE.4)THEN
        DO 20 I=1,4
   20   DEFDS(I)=DEFL(I)-DEFPT(I)
      ELSE
        DO 21 I=1,4
   21   DEFDS(I)=DEFL(I)
      ENDIF
      DEFDS(3)=0.5*DEFDS(3)
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      TAUD(1) =CPE(1,1)*DEFDS(1)+CPE(1,2)*DEFDS(2)+CPE(1,3)*DEFDS(4)
      TAUD(2) =CPE(2,1)*DEFDS(1)+CPE(2,2)*DEFDS(2)+CPE(2,3)*DEFDS(4)
      TAUD(3) = CEXY*DEFDS(3)
      TAUD(4) =-TAUD(1)-TAUD(2)
      SE(1)=TAUD(1)-ALFAT(1)
      SE(2)=TAUD(2)-ALFAT(2)
      SE(3)=TAUD(3)-ALFAT(3)
      SE(4)   =-SE(1)-SE(2)
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5*TDOTA2(SE,E))
      DEPMX = (TEQ-TEQY)/ETADT
      IF ((TEQ-TEQY).LT.0.D0. OR. DEPMX.LT.TOL0) THEN
        DEFQP=DEFQPT
        TAUY = TEQY
        CALL JEDNA1(DEFP,DEFPT,4)
        GO TO 500
      ENDIF
C
CE   3)  SOLUTION IS VICOPLASTIC.  OBTAIN ZERO OF THE GOVER. FN.
C
      IVP = 1
      PL1=1.0D0
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
      DBV2 = 0.D0
      DO 24 I = 1,3
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
         CALL     XMXMAT(CPH,DEFQP,EM1DVT,CY,AN,AN1,KDIM)
C
         DO 27 I = 1,3
         DO 27 J = 1,3
   27    CPHM(I,J) = CPH(J)*EN(I,J)
         CPH(4) = SQRT3*E(4)*CPH(4)
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
      SSH(4) = - SSH(1) - SSH(2)
      BV2 = (1. + (CPH(4)+CEXY)/ETADT)*E(4)
      SSH(3) = SE(3)/(1.+AK*BV2)
C
      TEQ=DSQRT(1.5*TDOTA2(SSH,E))
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
   32     BVOLD = BV2
          EOLD = E(4)
          DO 34 I = 1,3
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
          DBV2 = (BV2-BVOLD)/DDF
          DE2 = (E(4)-EOLD)/DDF
          DO 36 I = 1,3
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
C         NUMBER OF TRIALS FOR DEP EXCEEDS MAXBIS - STOP                      
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
C
  162 COEF = AK/ETADT
      DO 200 I = 1,2
      DUM = 0.D0
      DO 198 K = 1,3
      KA = IPERM(K)
  198 DUM = DUM + EN(I,K)*SSH(KA)
      DEVP(I) = COEF*DUM
  200 CONTINUE
      DEVP(3) = COEF*E(4)*SSH(3)
      DEVP(4) = -DEVP(1) - DEVP(2)
      IF  (IMIX.EQ.1) THEN
          DO 205 I = 1,2
          DUM = 0.D0
          DO 204 J = 1,3
          JJ = IPERM(J)
  204     DUM = DUM + CPHM(I,J)*DEVP(JJ)
          ALFA1(I) = ALFAT(I) + DUM
  205     CONTINUE
          ALFA1(3) = ALFAT(3) + CPH(4)*DEVP(3)
          ALFA1(4) = - ALFA1(1) - ALFA1(2)
      ENDIF
      DO 206 I = 1,4
      SS(I) = ALFA1(I) + SSH(I)
  206 TAUD(I) = SS(I) + DEVP(I)*ETADT
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      IF  (IETYP.EQ.0 .OR. IETYP.EQ.3) THEN
          TAUM = CM(1)*(DEFL(1)-DEFPT(1)-DEVP(1)) +
     1           CM(2)*(DEFL(2)-DEFPT(2)-DEVP(2))
      ELSE
          TAUM=CM(1)*(DEFL(1)-DEFPT(1))+CM(2)*(DEFL(2)-DEFPT(2))+
     &    CM(3)*(DEFL(4)-DEFPT(4))
      ENDIF
      TAUL(1)=TAUD(1)+TAUM
      TAUL(2)=TAUD(2)+TAUM
      TAUL(4)=TAUD(4)+TAUM
      IF (IETYP.EQ.0 .OR. IETYP.EQ.3) TAUL(4) = 0.D0
      TAUL(3)=TAUD(3)
C
C... TRANSFORM STRESS INTO GLOBAL AXES DIRECTION
C
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL CLEAR(TAU,4)
         CALL MNOZI2(TAU,TE,TAUL,4,4)
      ELSE
         CALL JEDNA1(TAU,TAUL,4)
      ENDIF
      IF(ISKNP.EQ.2) GO TO 700
C
C     CONSTITUTIVE MATRIX
C
      IF (IVP.EQ.1) GO TO 550
C
C     ELASTIC MATRIX
C
      DO 110 I=1,4
      DO 110 J=1,4
  110 ELAST(I,J)=CE(I,J)
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
      DENA = 1. + AK*BV2
      AV2 = CEXY/DENA
      BV2 = -(A1*BV2+AK*DBV2)*SSH(3)/DENA 
      DO 515 I = 1,3
      CP(1,I) = 0.D0
      CP(2,I) = 0.D0
      DO 514 J = 1,3
      JJ = IPERM(J)
      CP(2,I) = CP(2,I) + DEN(J,I)*SSH(JJ)
  514 CP(1,I) = CP(1,I) + EN(J,I)*SSH(JJ)
      CP(2,I) = 0.5*CP(2,I)
  515 CONTINUE
C
      W0 = 0.D0
      DO 520 I = 1,3
      I3 = I + 3
      II = IPERM(I)
      WJ(II) = 0.D0
      DO 518 K = 1,3
  518 WJ(II) = WJ(II) + CP(1,K)*AMAT(K,I)
      W0 = W0 + CP(1,I)*BV(I) + CP(2,I)*SSH(II)
  520 CONTINUE
      W0 = W0 + DVT*TAUY*EP + 2.*SSH(3)*E(4)*BV2 + DE2*SSH(3)*SSH(3)
      WJ(3) = 2.*SSH(3)*E(4)*AV2
      DO 522 I = 1,4
  522 WJ(I) = -WJ(I)/W0
      BVP2 = AK*SSH(3)*DE2
      DO 530 I = 1,3
      II = IPERM(I)
      BVP(I) = 0.D0
      DO 529 J = 1,3
      JJ = IPERM(J)
      BVP(I) = BVP(I) + DEN(I,J)*SSH(JJ)
      ELAST(II,JJ) = AMAT(I,J) + BV(I)*WJ(JJ)
  529 CONTINUE
      ELAST(II,3) = BV(I)*WJ(3)
      ELAST(3,II) = BV2*WJ(II)
      BVP(I) = AK*BVP(I)
  530 CONTINUE
      ELAST(3,3) = BV2*WJ(3) + AV2
      DO 535 I = 1,4
      DO 535 J = 1,4
      ELAST(I,J) = A1*SSH(I)*WJ(J) + AK*ELAST(I,J)
  535 CONTINUE
C
      DO 540 I = 1,3
      II = IPERM(I)
      DO 540 J = 1,3
      JJ = IPERM(J)
      CP(II,JJ) = 0.D0
      CP(II,3) = 0.D0
      DO 536 K = 1,3
      KA = IPERM(K)
      CP(II,JJ) = CP(II,JJ) + EN(I,K)*ELAST(KA,JJ)
  536 CP(II,3) = CP(II,3) + EN(I,K)*ELAST(KA,3)
      CP(II,JJ) = CP(II,JJ) + BVP(I)*WJ(JJ)
      CP(II,3) = CP(II,3) + BVP(I)*WJ(3)
      CP(3,JJ) = E(4)*ELAST(3,JJ) + BVP2*WJ(JJ)
  540 CONTINUE
      CP(3,3) = E(4)*ELAST(3,3) + BVP2*WJ(3)
      DO 545 I = 1,3
      II = IPERM(I)
      DO 545 J = 1,3
      JJ = IPERM(J)
      DUM = 0.D0
      DUM3 = 0.D0
      DO 544 K = 1,3
      KA = IPERM(K)
      DUM1 = CPE(I,K) + CM(K)
      DUM = DUM + DUM1*CP(KA,JJ)
      DUM3 = DUM3 + DUM1*CP(KA,3)
  544 CONTINUE
      ELAST(II,JJ) = CPE(I,J) + CM(J) - DUM/ETADT
      ELAST(II,3) = -0.5*DUM3/ETADT
      ELAST(3,JJ) = - CEXY*CP(3,JJ)/ETADT
  545 CONTINUE
      ELAST(3,3) = 0.5*CEXY*(1. - CP(3,3)/ETADT)
C
      IF  (IETYP.NE.1) THEN
          DO 549 I = 1,4
          ELAST(I,4) = 0.D0
  549     ELAST(4,I) = 0.D0
          ELAST(4,4) = 1.D0
      ENDIF
      DO 560 I = 1,4
      DO 560 J = I,4
      ELAST(I,J) = 0.5*(ELAST(I,J) + ELAST(J,I))
      ELAST(J,I) = ELAST(I,J)
  560 CONTINUE
C
CS....TRANSFORMACIJA MATRICE  ELAST()
C
  570 IF(DABS(BETA).GT.1.0D-10.AND.IETYP.NE.1) 
     1CALL TRAEL(ELAST,TE,4,3,3,ELAST)
      IF(DABS(BETA).GT.1.0D-10.AND.IETYP.EQ.1) 
     1CALL TRAEL(ELAST,TE,4,4,4,ELAST)
C
CE  UPDATE FROM PREVIOUS STEP
C
  700 IF(IETYP.EQ.0.OR.IETYP.EQ.3)
     1 DEFL(4)=-VZX/EZ*TAUL(1)-VYZ/EY*TAUL(2)+DEFPT(4)+DEVP(4)
      IF (IETYP.EQ.2) DEFL(4) = 0.D0
      DO 290 I=1,4
      DEFP(I) = DEFPT(I) + DEVP(I)
      DEF1(I)=DEFL(I)
  290 TAU1(I)=TAU(I)
      DEFP(3) = DEFP(3) + DEVP(3)
      RETURN
 1030 FORMAT (/ ' SUBROUTINE TI221'/
     1 ' NO SOLUTION IN BISECTING FOR F=0. '/
     1 ' IBIS =',I5/
     3 ' X-MINUS =',D14.6/' X-PLUS=',D14.6/' DDEFQP=',D14.6/
     4        ' F =',D14.6/' F-MINUS=',D14.6/' F-PLUS=',D14.6)
C
      END
C=======================================================================
      SUBROUTINE VP210(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,ETA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORM ELASTICITY CONSTANTS
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      DIMENSION FUN(30,*),CE(4,*),CPE(3,*),CM(*)
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
       DO 15 J=1,4
        CE(J,3)=ZER
   15   CE(3,J)=ZER
C     MATRICA CE
      POM=(1.-2.*VXY*VYZ*VZX-EX/EZ*VZX*VZX-EY/EX*VXY*VXY
     1-EZ/EY*VYZ*VYZ)/(EX*EY*EZ)
      CE(1,1)=(1./EZ-VYZ*VYZ/EY)/(EY*POM)
      CE(2,2)=(1./EX-VZX*VZX/EZ)/(EZ*POM)
      CE(4,4)=(1./EY-VXY*VXY/EX)/(EX*POM)
      CE(1,2)=(VZX*VYZ/EY+VXY/EX)/(EZ*POM)
      CE(1,4)=(VXY*VYZ/EX+VZX/EZ)/(EY*POM)
      CE(2,4)=(VXY*VZX/EZ+VYZ/EY)/(EX*POM)
      CE(3,3)=GXY
C  PLANE STRESS
      IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
        CE(1,1)=CE(1,1)-CE(1,4)*CE(1,4)/CE(4,4)
        CE(1,2)=CE(1,2)-CE(2,4)*CE(1,4)/CE(4,4)
        CE(2,2)=CE(2,2)-CE(2,4)*CE(2,4)/CE(4,4)
        CE(1,4)=ZER
        CE(2,4)=ZER
        CE(4,4)=ZER
      ENDIF
      DO 50 I=1,4
      DO 50 J=I,4
   50 CE(J,I)=CE(I,J)
C...   MATRIX  C'E
        CPE(1,1)=D13*(DVA*CE(1,1)-CE(1,2)-CE(1,4))
        CPE(1,2)=D13*(DVA*CE(1,2)-CE(2,2)-CE(2,4))
        CPE(2,1)=D13*(DVA*CE(1,2)-CE(1,1)-CE(1,4))
        CPE(2,2)=D13*(DVA*CE(2,2)-CE(1,2)-CE(2,4))
      IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
        DO 60 I=1,3
        CPE(I,3)=ZER
   60   CPE(3,I)=ZER
      ELSE
        CPE(1,3)=D13*(DVA*CE(1,4)-CE(2,4)-CE(4,4))
        CPE(2,3)=D13*(DVA*CE(2,4)-CE(1,4)-CE(4,4))
        CPE(3,1)=D13*(DVA*CE(1,4)-CE(1,1)-CE(1,2))
        CPE(3,2)=D13*(DVA*CE(2,4)-CE(1,2)-CE(2,2))
        CPE(3,3)=D13*(DVA*CE(4,4)-CE(1,4)-CE(2,4))
      ENDIF
C...   VECTOR   CM
      CM(1)=D13*(CE(1,1)+CE(1,2)+CE(1,4))
      CM(2)=D13*(CE(1,2)+CE(2,2)+CE(2,4))
      CM(3)=D13*(CE(1,4)+CE(2,4)+CE(4,4))
      RETURN
      END
