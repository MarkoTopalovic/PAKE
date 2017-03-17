C=======================================================================
C
C  VISKO-PLASTICNOST LJUSKE  -  MESOVITO OJACANJE    (19.8.1994)
C
C=======================================================================
      SUBROUTINE D8M20(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TAU820(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAU820( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC MATERIAL , MIXED HARDENING
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
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
     1          DEFP(*),ALFAT(*),ALFA1(*),EPR(6),COEFE(2)
      DIMENSION FUN(2,MATE,*)
      DIMENSION S(6),DEVP(6),SHATE(6),SHAT(6),EN(6),CP(6,6),
     4          DDVP(6),SHBAR(6)
      DATA MAXBIS/100/
      COEFE(1)=0.83333333D0
      COEFE(2)=COEFE(1)
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
      DUM=0.
      CALL MEL81(FUN,COEFE,ETP,1,DUM)
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
      CNY = (1.-2.*ANY)/(3.*(1.-ANY))
      C1 = 1. - CNY
      CZZ   =-ANY/(1.-ANY)
      AEM1 = 1./AE
      CE11 = AEM1*C1
      CE12 = -AEM1*CNY
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
      ETAD15 = 1.5*ETADT 
C
      IVP = 0
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.LT.4)THEN
          EPR(1)= C1 *(DEF(1)-DEFPT(1))-CNY*(DEF(2)-DEFPT(2))
          EPR(2)=-CNY*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))
          EPR(4)=0.5*DEF(4)-DEFPT(4)
          EPR(5)=0.5*DEF(5)-DEFPT(5)
          EPR(6)=0.5*DEF(6)-DEFPT(6)
          EPSM = CNY*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2))
      ELSE
          EPR(1)= C1 *DEF(1)-CNY*DEF(2)
          EPR(2)=-CNY*DEF(1)+C1 *DEF(2)
          EPR(4)=0.5*DEF(4)
          EPR(5)=0.5*DEF(5)
          EPR(6)=0.5*DEF(6)
          EPSM = CNY*(DEF(1)+DEF(2))
      ENDIF
C
C     ELASTIC DEVIATORIC STRESS
C
      DEPVP = 0.D0
      DO 20 I = 1,6
      DEVP(I) = 0.D0
      S(I) = EPR(I)/AE
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
C
      AK = ETAD15*DEPVP/SIGYA
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
      SHAT(3) = - SHAT(1) - SHAT(2)
      DENSH = 1. + AETA*AKD
      DO 11 I = 4,6
   11 SHAT(I) = SHATE(I)/DENSH 
      F = SQRT(1.5*TENDOT(SHAT)) - SIGYA 
C
      CALL BISEC (DEPVP,XM,XP,DX,F,FM,FP,AF,IB)
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
  162 DO 26 I = 1,6
   26 DEVP(I) = AKD*SHAT(I)
      S(1) = S(1) - CE11*DEVP(1) - CE12*DEVP(2)
      S(2) = S(2) - CE12*DEVP(1) - CE11*DEVP(2)
      S(3) = - S(1) - S(2)
      DO 163 I = 4,6
  163 S(I) = S(I) - AEM1*DEVP(I)
      IF(IATYP.LT.4)THEN
        EPSEXX = DEF(1) - DEFPT(1) - DEVP(1)
        EPSEYY = DEF(2) - DEFPT(2) - DEVP(2)
        EPSEZZ = CZZ*(EPSEXX+EPSEYY)
        DEF(3) = EPSEZZ + DEFPT(3) + DEVP(3)
        EPSM = CNY*(EPSEXX+EPSEYY)
      ELSE
        DEF(3)=CZZ*(DEF(1)+DEF(2))
        EPSM = -S(3)/CM
      ENDIF
C
  200 CONTINUE
C
       IF(IATYP.GE.4) THEN
          SIGM = CM*(AJOT-1.)/3.
       ELSE
         SIGM = CM*EPSM
       ENDIF
      DO 201 I = 1,6
      TAU(I) = S(I)
      IF (I.LT.3) TAU(I) = TAU(I) + SIGM
  201 CONTINUE
      TAU(3) = 0.D0
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
C
      COEF = 1.5*(ETADTA + CHAT)
      EPSM = (DEF(1)+DEF(2)+DEF(3))/3.
      EPR(1)=DEF(1)-EPSM-DEFPT(1)
      EPR(2)=DEF(2)-EPSM-DEFPT(2)
      EPR(3)=DEF(3)-EPSM-DEFPT(3)
      DO 90 I = 1,3
   90 SHATE(I) = EPR(I)/AE - ALFAT(I)
      SINT = SQRT(TENDOT(SHATE))
      DO 91 I = 1,6
   91 EN(I) = SHATE(I)/SINT
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
      ETP(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ETP(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ETP(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ETP(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ETP(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ETP(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
      ETP(3,1) = ETP(1,3)
      ETP(3,2) = ETP(2,3)
      DO 226 I = 4,6
      DO 224 J = 1,3
      ETP(J,I) = 0.5*CP(J,I)
  224 CONTINUE
      DO 226 J = I,6
      ETP(I,J) = 0.5*CP(I,J)
  226 CONTINUE
C
CE  STATIC CONDESATION
C
      DO 228 I=1,6
      IF (I.EQ.3) GO TO 228
      COEFR = ETP(3,I)/ETP(3,3)
      DO 227 J=I,6
      ETP(I,J)=ETP(I,J)-COEFR*ETP(3,J)
  227 ETP(J,I) = ETP(I,J)
  228 CONTINUE
C      
      DO 229 I=1,6
      ETP(3,I)=0.D0
      ETP(I,3)=0.D0
  229 CONTINUE
C
  230 DO 320 I = 1,6
      DEFP(I) = DEFPT(I) + DEVP(I)
      IF (I.GT.3) DEFP(I) = DEFP(I) + DEVP(I)
  320 ALFA1(I) = ALFAT(I) + CHAT*DEVP(I)
      PL1 = 1.D0
C
CE  UPDATE FROM PREVIOUS STEP
C
  700 IF(IATYP.LT.4)THEN
        DEF(3)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2))-DEFP(1)-DEFP(2)
      ELSE
        DEF(3)=CZZ*(DEF(1)+DEF(2))
      ENDIF
      DO 290 I=1,6
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
C
 1010 FORMAT (/' SUBROUTINE TAU8420'//' NO SOLUTION FOR F = 0.'/
     1 ' NUMBER OF BISECTIONS =',I5/
     2 ' DEPVP-MINUS =',D15.8/' DEPVP-PLUS=',D15.8/' DEPVP =',D15.8/
     3 ' F =',D15.8/' F-MINUS=',D15.8/' F-PLUS ',D15.8/
     4 ' S T O P ')
C
      RETURN
      END
