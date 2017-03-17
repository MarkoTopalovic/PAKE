C=======================================================================
C
C   VISKO-PLASTICNOST 3/D ELEMENT  -  MESOVITO OJACANJE    (20.08.1994)
C
C=======================================================================
      SUBROUTINE D3M20(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
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
      MATE=MREPER(4)
C      
      LTAU  =LPOCG
      LDEFT =LTAU   + 6
      LDEFPT=LDEFT  + 6
      LALFAT=LDEFPT + 6
      LTEQT =LALFAT + 6
      LDQPT =LTEQT  + 1
      LIPL  =LDQPT  + 1
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6
      LDEFP1=LDEFT1 + 6
      LALFA1=LDEFP1 + 6
      LTEQT1=LALFA1 + 6
      LDQPT1=LTEQT1 + 1
      LIPL1 =LDQPT1 + 1
C
      CALL TAU320(PLAST(LIPL),PLAST(LDEFPT),
     1            PLAST(LALFAT),PLAST(LTEQT),PLAST(LDQPT),
     1            PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LALFA1),PLAS1(LTEQT1),PLAS1(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAU320( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
C     VISKOPLASTICAN MATERIJAL SA MESOVITIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*)
      DIMENSION FUN(2,MATE,*)
      DIMENSION S(6),DEVP(6),SHATE(6),EN(6),CP(6,6),
     4          DDVP(6),SHBAR(6)
      DATA MAXBIS/100/
C
CE  INITIAL DATA
C
      IPL =PL
      DVA   =2.D0
      TR    =1.D0/3.
      DVT =2.D0/3.
C
      EE    =FUN(1,MAT,1)
      ANY   =FUN(2,MAT,1)
      SIGYV =FUN(1,MAT,2)
      TEQY0 =SIGYV
      CY    =FUN(2,MAT,2)
      ENEXP =FUN(1,MAT,3)
      EM    =FUN(2,MAT,3)
      ETA   =FUN(1,MAT,4)
      EM1   =1.-EM
      ENEXP1   =ENEXP-1.
      ENEXP2   =ENEXP-2.
C
      TEQY=TEQY0
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
      CALL MEL36(FUN,MATE)
      IF(IRAC.EQ.2) RETURN
C
      AF = 3.D0
      TOL = 1.D-10
      TOL0 = 1.D-10
      TOLDX = 1.D-12
      TOLF = 1.D-12
      SQRT15 = SQRT(1.5D0)
      SQRT23 = SQRT(2.D0/3.D0)
      AE = (1.+ANY)/EE
      COEFCH = 2./3.*(1.-EM)
C
      DO 3 I = 1,6
      DO 3 J = 1,6
      CP(I,J) = 0.D0
    3 CONTINUE
      CM = EE/(1.-2.*ANY)
      SIGY = TEQY
      FACTVP = 1.D-5
      IPR = 0
      DTETA = DT/ETA
      ETADT = ETA/DT
      AVP = 1. + DTETA/AE
      ETADTA = ETADT*AVP
C
      IVP = 0
      EPSM = (DEF(1)+DEF(2)+DEF(3))/3.
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.NE.4)THEN
        DO 18 I=1,3
   18   DEFDS(I)=DEF(I)-EPSM-DEFPT(I)
        DO 19 I=4,6
   19   DEFDS(I)=0.5*DEF(I)-DEFPT(I)
      ELSE
        DO 21 I=1,3
   21   DEFDS(I)=DEF(I)-EPSM
        DO 26 I=4,6
   26   DEFDS(I)=0.5*DEF(I)
      ENDIF
C
C     ELASTIC DEVIATORIC STRESS
C
      DEPVP = 0.D0
      DO 20 I = 1,6
      DEVP(I) = 0.D0
      S(I) = DEFDS(I)/AE
      SHATE(I) = S(I) - ALFAT(I)
   20 CONTINUE
      SIGYE = SQRT(1.5*TENDOT(SHATE)) 
      IF ((SIGYE-SIGY).LT.0.D0) GO TO 200
C
C     VISCOPLASTIC FLOW
C
      DEPMX = DTETA*(SIGYE-SIGY)
      IF (DEPMX.LE.TOL0) GO TO 200
      FM = 0.D0
      FP = -SIGY + SIGYE
      XM = 0.D0
      XP = 0.D0
      IVP = 1
      DEPVP = FACTVP*DEPMX
      DX = DEPVP
      IB = 0
      I = 0 
C
C     BISECTION LOOP
C
   10 I = I + 1
      DEFQP = DEFQPT + DEPVP
      SIGYA = SIGYV + CY*(EM*DEFQP)**ENEXP
      CHAT = COEFCH*ENEXP*CY*(DEFQP)**ENEXP1
      COEF = 1.5*(ETADTA + CHAT)
      F = -SIGYA - COEF*DEPVP + SIGYE
      CALL BISEC (DEPVP,XM,XP,DX,F,FM,FP,AF,IB)
      IF  (I.GT.MAXBIS) THEN                                                  
C                                                                               
C         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
C                                                                               
          WRITE (IZLAZ,1010) I,XM,XP,DEPVP,F,FM,FP
          STOP                                                                  
      ENDIF                                                                     
C                                                                               
      IF (IB.EQ.0) GO TO 10                                                   
C                                                                               
      IF (XM.LT.TOL0) GO TO 162
      IF (ABS(DX).LT.TOLDX .AND. ABS(F).LT.TOLF) GO TO 162
      IF ((ABS(DX)/(XP+XM)).GT.TOL) GO TO 10
C
  162 SINT = SQRT23*SIGYE
      TEQ  = SIGYA
      DEVPIN = SQRT15*DEPVP
      DO 25 I = 1,6
      EN(I) = SHATE(I)/SINT
      DEVP(I) = DEVPIN*EN(I)
   25 S(I) = (DEFDS(I) - DEVP(I))/AE
C
  200 CONTINUE
C
      SIGM = CM*EPSM
      DO 201 I = 1,6
      TAU(I) = S(I)
      IF (I.LE.3) TAU(I) = TAU(I) + SIGM
  201 CONTINUE
      IF  (IVP.EQ.0) THEN 
C
C         ELASTIC SOLUTION
C
          TEQ  =TEQY
          DEFQP=DEFQPT
          DO 600 I=1,6
  600     DEFP(I)=DEFPT(I)
          GO TO 700
      ENDIF
C
      IF(ISKNP.EQ.2) GO TO 230
C
C     ELASTIC-VISCOPLASTIC MATRIX
C
      EPB = ENEXP*CY*(EM*DEFQP)**ENEXP1
      DEP2 = ENEXP*ENEXP1*CY*DEFQP**ENEXP2
      AAVP = SQRT23*AE*SINT*(EM*EPB+COEF+EM1*DEPVP*DEP2)
      DO 221 I = 1,6
      SHBAR(I) = SHATE(I)
      IF (I.GT.3) SHBAR(I) = SHBAR(I) + SHBAR(I)
      DDVP(I) = SHBAR(I)/AAVP
  221 CONTINUE
      CDIAG = (1.-SQRT15*DEPVP/(AE*SINT))/AE
      COEF1 = SQRT15*DEPVP/(AE*AE*SINT*SINT*SINT)
      COEF2 = SQRT15/AE
      DO 223 I = 1,6
      DO 222 J = I,6
  222 CP(I,J) = COEF1*SHATE(I)*SHBAR(J) - COEF2*EN(I)*DDVP(J)
      CP(I,I) = CP(I,I) + CDIAG
  223 CONTINUE
C
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ELAST(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ELAST(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ELAST(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
      DO 226 I = 4,6
      DO 224 J = 1,3
      ELAST(J,I) = 0.5*CP(J,I)
  224 CONTINUE
      DO 226 J = I,6
      ELAST(I,J) = 0.5*CP(I,J)
  226 CONTINUE
C
      DO 227 I = 1,6
      DO 227 J = I,6
  227 ELAST(J,I) = ELAST(I,J)
C
  230 DO 320 I = 1,6
      DEFP(I) = DEFPT(I) + DEVP(I)
      IF (I.GT.3) DEFP(I) = DEFP(I) + DEVP(I)
  320 ALFA1(I) = ALFAT(I) + CHAT*DEVP(I)
      PL1 = 1.D0
C
CE  UPDATE FROM PREVIOUS STEP
C
  700 DO 290 I=1,6
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
C
 1010 FORMAT (/' SUBROUTINE TAU3420'//' NO SOLUTION FOR F = 0.'/
     1 ' NUMBER OF BISECTIONS =',I5/
     2 ' DEPVP-MINUS =',D15.8/' DEPVP-PLUS=',D15.8/' DEPVP =',D15.8/
     3 ' F =',D15.8/' F-MINUS=',D15.8/' F-PLUS ',D15.8/
     4 ' S T O P ')
C
      RETURN
      END
