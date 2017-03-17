C=======================================================================
C
C   PLASTICNOST 3/D ELEMENT  -  ANIZOTROPNI MATERIJAL  (17.04.1993)
C                               MESOVITO OJACANJE      (31.05.1993)
C                               GENERALNA ANIZOTROPIJA (25.01.1994)
C
C=======================================================================
      SUBROUTINE D3M10(TAU,DEF,IRAC,LPOCG,LPOC1)
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
C      MATE=MREPER(4)
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
      CALL TI310 (PLAST(LIPL),PLAST(LDEFPT),
     1            PLAST(LALFAT),PLAST(LTEQT),PLAST(LDQPT),
     1            PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LALFA1),PLAS1(LTEQT1),PLAS1(LDQPT1),
     1            A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI310 ( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC ANISOTROPIC MATERIAL , MIXED HARDENING  3D
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFP(*),
     1          ALFAT(*),ALFA1(*),FUN(*)
      DIMENSION SE(6),DEFL(6),TAUL(6)
      DIMENSION CE(6,6),CPE(3,3),CM(3),E(6),GB(3),XMX(6),
     1          Y0(6),CY(6),AN(6),AN1(6)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      KDIM=6
      IPL =PL
      IPL1=PL1
      DVA =2.D0
      DVT =DVA/3.
C.. CONSTANTS
      CALL       PEL310(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM)
      ANQ1  =ANQ-1.
      EM1   =1.-EM
      EM1DVT=EM1*DVT
      IMIX  = 1
      IF(DABS(EM1).LT.EPSIL) IMIX=0
      DUM   =0.D0
      CALL       ENTY(DUM,E,Y0,CY,AN,EM,DEFQPT,DEFQPT,KDIM)
C
      TEQY=TEQY0
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
C
C... TRANSFORM STRAIN INTO MATERIAL AXES DIRECTION
C
      IF(IRAC.EQ.2)THEN
C...................OVO ZAMENITI JER OVAKO JE UVEK U ITERACIJI 0 ELASTICNA
C...................TREBA RACUNATI C ELASTOPLASTIC.
        CALL JEDNA1(ELAST,CE,36)
C      IF(DABS(BETA).GT.1.0D-10) CALL TRAEL(ELAST,TSG,6,6,6,ELAST)
        RETURN
      ENDIF
      CALL CLEAR(DEFL,6)
      CALL MNOZI1(DEFL,TSG,DEF,6,6)
C
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
   10 DEFPT(I)=0.5*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
C
      IF(IATYP.NE.4)THEN
        DO 20 I=1,3
   20   DEFDS(I)=DEFL(I)-DEFPT(I)
        DO 25 I=4,6
   25   DEFDS(I)=0.5*DEFL(I)-DEFPT(I)
      ELSE
        DO 21 I=1,3
   21   DEFDS(I)=DEFL(I)
        DO 26 I=4,6
   26   DEFDS(I)=0.5*DEFL(I)
      ENDIF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      DO 40 I=1,3
       I3=I+3
       DUM=0.D0
       DO 30 J=1,3
   30  DUM=DUM+CPE(I,J)*DEFDS(J)
       TAUD(I) =DUM
       SE(I)   =DUM-ALFAT(I)
       DUM     =2.*DEFDS(I3)*CE(I3,I3)
       TAUD(I3)=DUM
       SE(I3)  =DUM-ALFAT(I3)
   40 CONTINUE
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5*TDOTAN(SE,E))
      IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
        TEQ  =TEQY
        DEFQP=DEFQPT
        CALL JEDNA1(DEFP,DEFPT,6)
        GO TO 500
      ENDIF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
C
CE       ITERATIONS
C     
      DEFQP=DEFQPT
      DUM=DEFQP
      IF(DEFQP.LE.1.D-4) DUM=1.D-4
C***************  OVO TRBA DA BUDE EKVIVALENTNI MODUL
      EP2=ANQ*CYQ*DUM**ANQ1      
C
      IB = 0
      IT = 0
      AF = 5.D0
      DDD= 0.1*(TEQ-TEQT)/EP2
      FP    = TEQ - TEQT
      TAUY  = TEQT
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      DDEFQP= DDD
C
      CALL       XMXMAT(XMX,DUM,EM1DVT,CY,AN,AN1,KDIM)
      CALL       ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      CALL       INI310(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
C
      TOLD= TAUY
      DOLD= DEFQPT
      EP  = EP2
  100 IT=IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        WRITE(IZLAZ,2000)
        STOP
      ENDIF
C
      DEFQP=DEFQPT+DDEFQP
C
      CALL       ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      IF(IMIX.EQ.1)THEN
        CALL     XMXMAT(XMX,DEFQP,EM1DVT,CY,AN,AN1,KDIM)
      ENDIF
      CALL     INI310(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
C
      DLAM=1.5*DDEFQP/TAUY
      CALL       DEV310(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
C
      TEQ=DSQRT(1.5*TDOTAN(TAUD,E))
      FB = TEQ-TAUY
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
C
      DDF   = DEFQP-DOLD
      IF((DABS(DDF).GT.1.D-8 .OR. IT.EQ.1).AND.IB1.EQ.1)THEN
        EP  = DABS((TAUY-TOLD)/DDF)
        TOLD= TAUY
        DOLD= DEFQP
      ENDIF
C
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
 2000 FORMAT(' ','DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI36')
C
CE   4)  DETERMINE SOLUTION 
C
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2)THEN
C**** BRISI OVU LINIJU
C      EP  = AN(1)*CY(1)*DEFQP**AN1(1)
C**** 
      CALL       CEP310(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                  DL,DXX,DYY,DLAM,TAUY,XMX,EP)
      ENDIF
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
        DDEFP1 =DLAM*( (E(1)+E(2))*TAUD(1)-E(1)*TAUD(2)-E(2)*TAUD(3))
        DDEFP2 =DLAM*(-E(1)*TAUD(1)+(E(1)+E(3))*TAUD(2)-E(3)*TAUD(3))
        DDEFP4 =DLAM*E(4)*TAUD(4)
        DDEFP5 =DLAM*E(5)*TAUD(5)
        DDEFP6 =DLAM*E(6)*TAUD(6)
        DEFP(1)=DDEFP1+DEFPT(1)
        DEFP(2)=DDEFP2+DEFPT(2)
CC        DEFP(3)=DLAM*(-E(2)*TAUD(1)-E(3)*TAUD(2)+(E(2)+E(3))*TAUD(3))+
CC     &          DEFPT(3)
        DEFP(3)=-DDEFP1-DDEFP2+DEFPT(3)
        DEFP(4)=DDEFP4+DEFPT(4)
        DEFP(5)=DDEFP5+DEFPT(5)
        DEFP(6)=DDEFP6+DEFPT(6)
C
      IF(IMIX.EQ.1)THEN
        ALFA1(1)=XMX(1)*DDEFP1+ALFAT(1)
        ALFA1(2)=XMX(2)*DDEFP2+ALFAT(2)
        ALFA1(3)=-ALFA1(1)-ALFA1(2)
        ALFA1(4)=XMX(4)*DDEFP4+ALFAT(4)
        ALFA1(5)=XMX(5)*DDEFP5+ALFAT(5)
        ALFA1(6)=XMX(6)*DDEFP6+ALFAT(6)
        DO 160 I=1,6
  160   TAUD(I)=TAUD(I)+ALFA1(I)
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      TAUM=CM(1)*(DEFL(1)-DEFP(1))+CM(2)*(DEFL(2)-DEFP(2))+
     &     CM(3)*(DEFL(3)-DEFP(3))
      DO 200 I=1,3
  200 TAUL(I)=TAUD(I)+TAUM
      DO 205 I=4,6
      TAUL(I)=TAUD(I)
  205 DEFP(I)=2.*DEFP(I)
C
C... TRANSFORM STRESS INTO GLOBAL AXES DIRECTION
C
      CALL CLEAR(TAU,6)
      CALL MNOZI2(TAU,TSG,TAUL,6,6)
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
      DEF1(I)=DEFL(I)
  290 TAU1(I)=TAU(I)
      RETURN
      END
C=======================================================================
      SUBROUTINE PEL310(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORM ELASTICITY CONSTANTS
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      DIMENSION FUN(29,*),CE(6,*),CPE(3,*),CM(*)
      DIMENSION Y0(*),CY(*),AN(*),AN1(*)
C
      D13 =1.D0/3.
      ONE =1.D0
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
      CYQ   =CY(1)
      ANQ   =AN(1)
      DO 15 I=1,6
      DO 15 J=1,6
   15 CE(I,J)=ZER
C     MATRICA CE
      POM=(ONE-DVA*VXY*VYZ*VZX-EX/EZ*VZX*VZX-EY/EX*VXY*VXY
     1-EZ/EY*VYZ*VYZ)/(EX*EY*EZ)
      CE(1,1)=(ONE/EZ-VYZ*VYZ/EY)/(EY*POM)
      CE(2,2)=(ONE/EX-VZX*VZX/EZ)/(EZ*POM)
      CE(3,3)=(ONE/EY-VXY*VXY/EX)/(EX*POM)
      CE(1,2)=(VZX*VYZ/EY+VXY/EX)/(EZ*POM)
      CE(1,3)=(VXY*VYZ/EX+VZX/EZ)/(EY*POM)
      CE(2,3)=(VXY*VZX/EZ+VYZ/EY)/(EX*POM)
      CE(4,4)=GXY
      CE(5,5)=GYZ
      CE(6,6)=GZX
      DO 50 I=1,6
      DO 50 J=I,6
   50 CE(J,I)=CE(I,J)
C...   MATRIX  C'E
      CPE(1,1)=D13*(DVA*CE(1,1)-CE(1,2)-CE(1,3))
      CPE(1,2)=D13*(DVA*CE(1,2)-CE(2,2)-CE(2,3))
      CPE(1,3)=D13*(DVA*CE(1,3)-CE(2,3)-CE(3,3))
      CPE(2,1)=D13*(DVA*CE(1,2)-CE(1,1)-CE(1,3))
      CPE(2,2)=D13*(DVA*CE(2,2)-CE(1,2)-CE(2,3))
      CPE(2,3)=D13*(DVA*CE(2,3)-CE(1,3)-CE(3,3))
      CPE(3,1)=D13*(DVA*CE(1,3)-CE(1,1)-CE(1,2))
      CPE(3,2)=D13*(DVA*CE(2,3)-CE(1,2)-CE(2,2))
      CPE(3,3)=D13*(DVA*CE(3,3)-CE(1,3)-CE(2,3))
C...   VECTOR   CM
      CM(1)=D13*(CE(1,1)+CE(1,2)+CE(1,3))
      CM(2)=D13*(CE(1,2)+CE(2,2)+CE(2,3))
      CM(3)=D13*(CE(1,3)+CE(2,3)+CE(3,3))
      RETURN
      END
C======================================================================
      SUBROUTINE DEV310(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TAUD(*),E(*),SE(*),GB(*),XMX(*)
        DL = 1. + DLAM*(CXX+CYY + DLAM*CED)
        TAUD(1) = (SE(1) + DLAM*DXX)/DL
        TAUD(2) = (SE(2) + DLAM*DYY)/DL
        TAUD(3) = - TAUD(1) - TAUD(2)
        TAUD(4) = SE(4)/(1.+DLAM*(GB(1)+XMX(4)*E(4)))
        TAUD(5) = SE(5)/(1.+DLAM*(GB(2)+XMX(5)*E(5)))
        TAUD(6) = SE(6)/(1.+DLAM*(GB(3)+XMX(6)*E(6)))
      RETURN
      END      
C======================================================================
      SUBROUTINE INI310(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SE(*),E(*),CE(6,*),CPE(3,*),GB(*),XMX(*)
C
      DVA=2.D0
      XX1 = E(1)+DVA*E(2)
      XX2 = E(1)-E(3)
      XX3 = DVA*E(2)+E(3)
      CXX =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2) - XX3*CPE(1,3)
      CYX = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2)) + XX3*CPE(2,3)
      XX1 = E(1)-E(2)
      XX2 = E(1)+DVA*E(3)
      XX3 = E(2)+DVA*E(3)
      CYY = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2)) - XX3*CPE(2,3)
      CXY =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2) + XX3*CPE(1,3)
      CED = CXX*CYY - CXY*CYX
      DXX = CYY*SE(1)+CXY*SE(2)
      DYY = CYX*SE(1)+CXX*SE(2)
      GB(1)=CE(4,4)*DVA*E(4)
      GB(2)=CE(5,5)*DVA*E(5)
      GB(3)=CE(6,6)*DVA*E(6)
      RETURN
      END
C======================================================================
      SUBROUTINE CEP310(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                  DL,DXX,DYY,DLAM,TAUY,XMX,EP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION E(*),CE(6,*),CPE(3,*),TAUD(*),CM(*),GB(*),XMX(*)
      DIMENSION P(3,3),Q(3),W(6),DDL(6),T(6),R(6),CMB(6),XEN(6,6),
     &          CP(6,6)
      EQUIVALENCE (P(1,1),XEN(1,1))
C
      ZERO  =0.D0
      ONE   =1.D0
      DVA   =2.D0
      PO    =0.5D0
      DVT   =DVA/3.
C
      ELP = (1.5 - EP*DLAM)/TAUY
      CALL CLEAR(ELAST,36)
C
      DS12 = TAUD(1) - TAUD(2)
      DS13 = TAUD(1) - TAUD(3)
      DS23 = TAUD(2) - TAUD(3)
      AXX = E(1)*DS12 + DVA*E(2)*DS13 + E(3)*DS23
      AYY = DVA*E(3)*DS23 - E(1)*DS12 + E(2)*DS13
      DO 25 I=1,3
        I3=I+3
        P(1,I) = (CPE(1,I) + DLAM*(CYY*CPE(1,I) + CXY*CPE(2,I)))/DL
        P(2,I) = (CPE(2,I) + DLAM*(CXX*CPE(2,I) + CYX*CPE(1,I)))/DL
        P(3,I) = -P(1,I)-P(2,I)
        DUM  = 1. + DLAM*GB(I)
        R(I3) = DVA*CE(I3,I3)/DUM
        T(I3) = GB(I)*TAUD(I3)/DUM
C
        W(I) = AXX*P(1,I) + AYY*P(2,I)
        W(I3)= DVA*E(I3)*TAUD(I3)*R(I3)
   25 CONTINUE
      DLL = CXX + CYY + DVA*DLAM*CED
      Q(1) = (DXX - DLL*TAUD(1))/DL
      Q(2) = (DYY - DLL*TAUD(2))/DL
      Q(3) = -Q(1)-Q(2)
      W0 = ELP*(AXX*Q(1) + AYY*Q(2) - 
     &   DVA*(E(4)*TAUD(4)*T(4)+E(5)*TAUD(5)*T(5)+E(6)*TAUD(6)*T(6)))
     &   - DVT*EP*TAUY
      DUM  =  -ELP/W0
      DO 15 I=1,6 
   15 DDL(I) = DUM*W(I)
C
      DO 40 I=1,3
      DO 20 J=1,3
        ELAST(I,J)=P(I,J)+Q(I)*DDL(J)
   20 CONTINUE
      DO 35 J=4,6
        ELAST(I,J)=Q(I)*DDL(J)
   35 CONTINUE
   40 CONTINUE
      DO 60 I=4,6
        DO 50 J=1,6
   50   ELAST(I,J)=-T(I)*DDL(J)
        ELAST(I,I)=ELAST(I,I)+R(I)
   60 CONTINUE
C      
      AXX =  (E(1)+E(2))*CM(1)-E(1)*CM(2)-E(2)*CM(3)
      AYY = -E(1)*CM(1)+(E(1)+E(3))*CM(2)-E(3)*CM(3)
      AZZ = -E(2)*CM(1)-E(3)*CM(2)+(E(2)+E(3))*CM(3)
      DUM  = AXX*TAUD(1)+AYY*TAUD(2)+AZZ*TAUD(3)
      DO 70 I=1,3 
        PM = AXX*ELAST(1,I)+AYY*ELAST(2,I)+AZZ*ELAST(3,I)
        QM = DUM*DDL(I)
        CMB(I)= CM(I)-DLAM*PM-QM
   70 CONTINUE
      DO 71 I=4,6 
        PM = AXX*ELAST(1,I)+AYY*ELAST(2,I)+AZZ*ELAST(3,I)
        QM = DUM*DDL(I)
        CMB(I)= -DLAM*PM-QM
   71 CONTINUE
C
      CALL CLEAR(XEN,36)
      XEN(1,1)= XMX(1)*(E(1)+E(2))
      XEN(1,2)=-XMX(1)*E(1)
      XEN(1,3)=-XMX(1)*E(2)
      XEN(2,1)=-XMX(2)*E(1)
      XEN(2,2)= XMX(2)*(E(1)+E(3))
      XEN(2,3)=-XMX(2)*E(3)
      XEN(3,1)=-XMX(3)*E(2)
      XEN(3,2)=-XMX(3)*E(3)
      XEN(3,3)= XMX(3)*(E(2)+E(3))
      XEN(4,4)= XMX(4)*E(4)
      XEN(5,5)= XMX(5)*E(5)
      XEN(6,6)= XMX(6)*E(6)
C      CALL WRR(XEN,36,'XEN ')
C... (A)  ,  T  DOBIJA NOVE VREDNOSTI
      T(1)=XEN(1,1)*TAUD(1)+XEN(1,2)*TAUD(2)+XEN(1,3)*TAUD(3)
      T(2)=XEN(2,1)*TAUD(1)+XEN(2,2)*TAUD(2)+XEN(2,3)*TAUD(3)
      T(3)=XEN(3,1)*TAUD(1)+XEN(3,2)*TAUD(2)+XEN(3,3)*TAUD(3)
      T(4)=XEN(4,4)*TAUD(4)
      T(5)=XEN(5,5)*TAUD(5)
      T(6)=XEN(6,6)*TAUD(6)
C... (B)
      DO 78 I=1,3
      DO 77 J=1,3
   77 XEN(I,J)=DLAM*XEN(I,J)
      XEN(I,I)=XEN(I,I)+ONE
      I3=I+3
      XEN(I3,I3)=DLAM*XEN(I3,I3)+ONE
   78 CONTINUE
C      CALL WRR(XEN,36,'XEN ')
C      CALL WRR(T   ,6,'T   ')
C
      DO 85 I=1,6
      DO 85 J=I,6
       DUM=ZERO
        DO 83 K=1,6
   83   DUM=DUM+XEN(I,K)*ELAST(K,J)
       CP(I,J)=DUM+T(I)*DDL(J)
   85 CONTINUE
      DO 87 I=4,6
      DO 87 J=I,6
        CP(I,J)=PO*CP(I,J)
   87 CONTINUE
      DO 88 I=1,6
      DO 88 J=I,6
        ELAST(I,J)=CP(I,J)
   88 CONTINUE
C
      DO 82 I=1,3
      DO 80 J=I,3
   80   ELAST(I,J)=ELAST(I,J)+CMB(J)
      DO 81 J=4,6
   81   ELAST(I,J)=PO*(ELAST(I,J)+CMB(J))
   82 CONTINUE
C
      DO 90 I=1,6
      DO 90 J=I,6
        ELAST(J,I)=ELAST(I,J)
   90 CONTINUE
      RETURN
      END
