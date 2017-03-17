C=======================================================================
C
C   VISKO-PLASTICNOST 2/D ELEMENT  -  MESOVITO OJACANJE    (17.2.1994)
C
C=======================================================================
      SUBROUTINE D2M20(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TAU220(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAU220( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   VISCOPLASTIC MATERIAL , MIXED HARDENING
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*)
      DIMENSION FUN(2,MATE,*)
      DIMENSION S(4),DEVP(4),SHATE(4),SHAT(4),EN(4),CP(4,4),
     4          DDVP(4),SHBAR(4)
      DATA MAXBIS/100/
C
C
CE  INITIAL DATA
C
      IAX = 1
      IF (IETYP.EQ.2) IAX = 0
C     IAX = 0 - PLANE STRAIN
C     IAX = 1 - AXIAL SYMMETRY  OR PLANE STRESS
C      
      IPL =PL
      IPL1=PL1
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
      CALL MEL26(FUN,MATE)
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
      IF  (IETYP.EQ.0 .OR. IETYP.EQ.3) THEN
          CNY = (1.-2.*ANY)/(3.*(1.-ANY))
          C1 = 1. - CNY
          CZZ   =-ANY/(1.-ANY)
          AEM1 = 1./AE
          CE11 = AEM1*C1
          CE12 = -AEM1*CNY
      ENDIF
C
      DO 3 I = 1,4
      DO 3 J = 1,4
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
      ETAD15 = 1.5*ETADT 
C
      IVP = 0
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.LT.4)THEN
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNY*(DEF(2)-DEFPT(2))
          DEFDS(2)=-CNY*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))
          EPSM = CNY*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2))
        ELSE
          EPSM = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EPSM-DEFPT(1)
          DEFDS(2)=DEF(2)-EPSM-DEFPT(2)
          DEFDS(4)=DEF(4)-EPSM-DEFPT(4)
        ENDIF
        DEFDS(3)=0.5*DEF(3)-DEFPT(3)
      ELSE
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          DEFDS(1)= C1 *DEF(1)-CNY*DEF(2)
          DEFDS(2)=-CNY*DEF(1)+C1 *DEF(2)
          EPSM = CNY*(DEF(1)+DEF(2))
        ELSE
          EPSM = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EPSM
          DEFDS(2)=DEF(2)-EPSM
          DEFDS(4)=DEF(4)-EPSM
        ENDIF
        DEFDS(3)=0.5*DEF(3)
      ENDIF
C
C     ELASTIC DEVIATORIC STRESS
C
      DEPVP = 0.D0
      DO 20 I = 1,3
      DEVP(I) = 0.D0
      S(I) = DEFDS(I)/AE
      SHATE(I) = S(I) - ALFAT(I)
   20 CONTINUE
      S(4) = -S(1) - S(2)
      SHATE(4) = -SHATE(1) - SHATE(2)
      SIGYE = SQRT(1.5*TENDO2(SHATE))
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
      IF (IETYP.EQ.0 .OR.IETYP.EQ.3) GO TO 11
C
C     PLANE STRAIN OR AXIAL SYMMETRY
C
      COEF = 1.5*(ETADTA + CHAT)
      F = -SIGYA - COEF*DEPVP + SIGYE
      GO TO 12
C
C     PLANE STRESS
C
   11 AK = ETAD15*DEPVP/SIGYA
      COEF = ETADT + CHAT
      AMAT11 = CE11 + COEF
      AMAT12 = CE12
      DET = AMAT11*AMAT11 - AMAT12*AMAT12
      DMAT11 = AMAT11/DET
      DMAT12 = -AMAT12/DET
      BE1 = DMAT11*SHATE(1) + DMAT12*SHATE(2)
      BE2 = DMAT12*SHATE(1) + DMAT11*SHATE(2)
      AETA = AEM1 + COEF
      AKD = AK*DTETA
      DET = (DMAT11+AKD)*(DMAT11+AKD) - DMAT12*DMAT12
      BMAT11 = (DMAT11+AKD)/DET
      BMAT12 = -DMAT12/DET
      SHAT(1) = BMAT11*BE1 + BMAT12*BE2
      SHAT(2) = BMAT12*BE1 + BMAT11*BE2
      SHAT(4) = - SHAT(1) - SHAT(2)
      SHAT(3) = SHATE(3)/(1.+AETA*AKD)
      F = SQRT(1.5*TENDO2(SHAT)) - SIGYA
C
   12 CALL BISEC (DEPVP,XM,XP,DX,F,FM,FP,AF,IB)
      IF  (I.GT.MAXBIS) THEN                                                  
C                                                                               
C         NUMBER OF TRIALS FOR DEP EXCEEDS MAXBIS - STOP                      
C                                                                               
          WRITE (IZLAZ,1010) I,XM,XP,DEPVP,F,FM,FP
          STOP                                                                  
      ENDIF
C                                                                               
      IF (IB.EQ.0) GO TO 10                                                   
C     
      IF (XM.LT.TOL0) GO TO 162
      IF ((ABS(DX)/(XP+XM)).GT.TOL) GO TO 10
C
  162 IF (IETYP.EQ.0 .OR. IETYP.EQ.3) GO TO 163
C
C     PLANE STRAIN OR AXIAL SYMMETRY
C
      SINT = SQRT23*SIGYE
      TEQ  = SIGYA
      DEVPIN = SQRT15*DEPVP
      DO 25 I = 1,4
      EN(I) = SHATE(I)/SINT
      DEVP(I) = DEVPIN*EN(I)
   25 S(I) = (DEFDS(I) - DEVP(I))/AE
      GO TO 200
C
C     PLANE STRESS
C
  163 DO 26 I = 1,4
   26 DEVP(I) = AKD*SHAT(I)
      S(1) = S(1) - CE11*DEVP(1) - CE12*DEVP(2)
      S(2) = S(2) - CE12*DEVP(1) - CE11*DEVP(2)
      S(3) = S(3) - AEM1*DEVP(3)
      S(4) = - S(1) - S(2)
      IF(IATYP.LT.4)THEN
        EPSEXX = DEF(1) - DEFPT(1) - DEVP(1)
        EPSEYY = DEF(2) - DEFPT(2) - DEVP(2)
        EPSEZZ = CZZ*(EPSEXX+EPSEYY)
        DEF(4) = EPSEZZ + DEFPT(4) + DEVP(4)
        EPSM = CNY*(EPSEXX+EPSEYY)
      ELSE
        DEF(4)=CZZ*(DEF(1)+DEF(2))
        EPSM = -S(4)/CM
      ENDIF
C
  200 CONTINUE
C
       IF(IATYP.GE.4) THEN
          SIGM = CM*(AJOT-1.)/3.
       ELSE
         SIGM = CM*EPSM
       ENDIF
      DO 201 I = 1,4
      TAU(I) = S(I)
      IF (I.NE.3) TAU(I) = TAU(I) + SIGM
  201 CONTINUE
      IF (IETYP.EQ.0 .OR. IETYP.EQ.3) TAU(4) = 0.D0
      IF  (IVP.EQ.0) THEN 
C
C         ELASTIC SOLUTION
C
          TEQ  =TEQY
          DEFQP=DEFQPT
          DO 600 I=1,4
  600     DEFP(I)=DEFPT(I)
          GO TO 700
      ENDIF
C
      IF(ISKNP.EQ.2) GO TO 230
C
C     ELASTIC-VISCOPLASTIC MATRIX
C
      IF (IETYP.NE.0 .AND. IETYP.NE.3) GO TO 93
C
      COEF = 1.5*(ETADTA + CHAT)
      EPSM = (DEF(1)+DEF(2)+DEF(4))/3.
      DEFDS(1)=DEF(1)-EPSM-DEFPT(1)
      DEFDS(2)=DEF(2)-EPSM-DEFPT(2)
      DEFDS(4)=DEF(4)-EPSM-DEFPT(4)
      DO 90 I = 1,4
   90 SHATE(I) = DEFDS(I)/AE - ALFAT(I)
      SINT = SQRT(TENDO2(SHATE))
      DO 91 I = 1,4
   91 EN(I) = SHATE(I)/SINT
C
   93 EPB = ENEXP*CY*(EM*DEFQP)**ENEXP1
      DEP2 = ENEXP*ENEXP1*CY*(DEFQP)**ENEXP2
      AAVP = SQRT23*AE*SINT*(EM*EPB+COEF+EM1*DEPVP*DEP2)
      DO 221 I = 1,4
      SHBAR(I) = SHATE(I)
      IF (I.EQ.3) SHBAR(I) = SHBAR(I) + SHBAR(I)
      DDVP(I) = SHBAR(I)/AAVP
  221 CONTINUE
      CDIAG = (1.-SQRT15*DEPVP/(AE*SINT))/AE
      COEF1 = SQRT15*DEPVP/(AE*AE*SINT*SINT*SINT)
      COEF2 = SQRT15/AE
      DO 223 I = 1,4
      DO 222 J = 1,4
  222 CP(I,J) = COEF1*SHATE(I)*SHBAR(J) - COEF2*EN(I)*DDVP(J)
      CP(I,I) = CP(I,I) + CDIAG
  223 CONTINUE
C
      IF  (IAX.EQ.0) THEN
          ELAST(1,1) = (2.*CP(1,1)-CP(1,2)-CP(1,4)+CM)/3.
          ELAST(1,2) = (2.*CP(1,2)-CP(1,1)-CP(1,4)+CM)/3.
          ELAST(2,2) = (2.*CP(2,2)-CP(2,1)-CP(2,4)+CM)/3.
          ELAST(2,1) = ELAST(1,2)
      ELSE
          ELAST(1,1) = (2.*CP(1,1)-CP(1,2)-CP(1,4)+CM)/3.
          ELAST(1,2) = (2.*CP(1,2)-CP(1,1)-CP(1,4)+CM)/3.
          ELAST(2,2) = (2.*CP(2,2)-CP(2,1)-CP(2,4)+CM)/3.
          ELAST(1,4) = (2.*CP(1,4)-CP(1,1)-CP(1,2)+CM)/3.
          ELAST(2,4) = (2.*CP(2,4)-CP(2,1)-CP(2,2)+CM)/3.
          ELAST(4,4) = (2.*CP(4,4)-CP(4,1)-CP(4,2)+CM)/3.
          ELAST(2,1) = ELAST(1,2)
          ELAST(4,1) = ELAST(1,4)
          ELAST(4,2) = ELAST(2,4)
      ENDIF
      DO 224 J = 1,4
      ELAST(J,3) = 0.5*CP(J,3)
      ELAST(3,J) = ELAST(J,3)
  224 CONTINUE
C
      IF (IETYP.NE.0 .AND. IETYP.NE.3) GO TO 230
C
C     STATIC CONDENSATION FOR PLANE STRESS
C
      DO 227 I = 1,3
      COEFR = ELAST(I,4)/ELAST(4,4)
      DO 227 J = I,3
      ELAST(I,J) = ELAST(I,J) - COEFR*ELAST(J,4)
      ELAST(J,I) = ELAST(I,J)
  227 CONTINUE
C
  230 DO 320 I = 1,4
      DEFP(I) = DEFPT(I) + DEVP(I)
      IF (I.EQ.3) DEFP(I) = DEFP(I) + DEVP(I)
  320 ALFA1(I) = ALFAT(I) + CHAT*DEVP(I)
      PL1 = 1.D0
C
CE  UPDATE FROM PREVIOUS STEP
C
  700 IF(IATYP.LT.4)THEN
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2))-DEFP(1)-DEFP(2)
        ELSEIF(IETYP.EQ.2)THEN
          DEF(4)=0.D0
        ENDIF
      ELSE
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          DEF(4)=CZZ*(DEF(1)+DEF(2))
        ELSEIF(IETYP.EQ.2)THEN
          DEF(4)=0.D0
        ENDIF
      ENDIF
      DO 290 I=1,4
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
C
 1010 FORMAT (/' SUBROUTINE TAU2420'//' NO SOLUTION FOR F = 0.'/
     1 ' NUMBER OF BISECTIONS =',I5/
     2 ' DEPVP-MINUS =',D15.8/' DEPVP-PLUS=',D15.8/' DEPVP =',D15.8/
     3 ' F =',D15.8/' F-MINUS=',D15.8/' F-PLUS ',D15.8/
     4 ' S T O P ')
C
      RETURN
      END
