C=======================================================================
C
C   PLASTICNOST 2/D ELEMENT  -  ANIZOTROPNI MATERIJAL  (17.04.1993)
C                               MESOVITO OJACANJE  2D  (08.06.1993)
C                               GENERALNA ANIZOTROPIJA (24.09.1993)
C
C=======================================================================
      SUBROUTINE D2M10(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TI210 (A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI210 ( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA ELASTOPLASTIC
C     ELASTOPLASTICAN MATERIJAL SA IZOTROPNIM OJACANJEM
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
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /UGAOVL/ TE(4,4)
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFP(*),
     1          ALFAT(*),ALFA1(*),FUN(*)
      DIMENSION SE(4),DEFL(4),TAUL(4)
      DIMENSION CE(4,4),CPE(3,3),CM(3),E(6),GB(1),XMX(4),
     1          Y0(6),CY(6),AN(6),AN1(6)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      KDIM=4
      IPL =PL
      IPL1=PL1
      DVA =2.D0
      DVT =DVA/3.
C.. CONSTANTS
      CALL       PEL210(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM)
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
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFPT(3)=0.5*DEFPT(3)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
C
      IF(IATYP.NE.4)THEN
        DEFDS(1)=DEFL(1)-DEFPT(1)
        DEFDS(2)=DEFL(2)-DEFPT(2)
        DEFDS(4)=DEFL(4)-DEFPT(4)
        DEFDS(3)=0.5*DEFL(3)-DEFPT(3)
      ELSE
        DEFDS(1)=DEFL(1)
        DEFDS(2)=DEFL(2)
        DEFDS(4)=DEFL(4)
        DEFDS(3)=0.5*DEFL(3)
      ENDIF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      TAUD(1) =CPE(1,1)*DEFDS(1)+CPE(1,2)*DEFDS(2)+CPE(1,3)*DEFDS(4)
      TAUD(2) =CPE(2,1)*DEFDS(1)+CPE(2,2)*DEFDS(2)+CPE(2,3)*DEFDS(4)
      TAUD(3) =2.*DEFDS(3)*CE(3,3)-ALFAT(3)
      TAUD(4) =-TAUD(1)-TAUD(2)
      SE(1)=TAUD(1)-ALFAT(1)
      SE(2)=TAUD(2)-ALFAT(2)
      SE(3)=TAUD(3)-ALFAT(3)
      SE(4)   =-SE(1)-SE(2)
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5*TDOTA2(SE,E))
C      IF((TEQ-TEQY)/TEQY.LT.1.D-5.OR.ITER.EQ.0)THEN
      IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
        TEQ  =TEQY
        DEFQP=DEFQPT
        CALL JEDNA1(DEFP,DEFPT,4)
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
      CALL       INI210(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX,
     &                  IETYP)
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
      CALL     INI210(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX,
     &                IETYP)
C
      DLAM=1.5*DDEFQP/TAUY
      CALL       DEV210(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
C
      TEQ=DSQRT(1.5*TDOTA2(TAUD,E))
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
C      CALL WRR(TAUD,4,'shet')
C
CE      ...   ( DEVIATORIC STRESS )
C
 2000 FORMAT(' ','DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI210')
C
CE   4)  DETERMINE SOLUTION 
C
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2) THEN
C      EP  =ANQ*CYQ*DEFQP**ANQ1      
         CALL  CEP210(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                DL,DXX,DYY,DLAM,TAUY,XMX,EP)
CS.... TRANSFORMACIJA MATRICE  ELAST()
CE     TRANSFORM ELAST MATRIX
      IF(DABS(BETA).GT.1.0D-10.AND.IETYP.NE.1) 
     1CALL TRAEL(ELAST,TE,4,3,3,ELAST)
      IF(DABS(BETA).GT.1.0D-10.AND.IETYP.EQ.1) 
     1CALL TRAEL(ELAST,TE,4,4,4,ELAST)
      ENDIF
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
        DDEFP1 =DLAM*( (E(1)+E(2))*TAUD(1)-E(1)*TAUD(2)-E(2)*TAUD(4))
        DDEFP2 =DLAM*(-E(1)*TAUD(1)+(E(1)+E(3))*TAUD(2)-E(3)*TAUD(4))
        DDEFP3 =DLAM*E(4)*TAUD(3)
        DEFP(1)=DDEFP1+DEFPT(1)
        DEFP(2)=DDEFP2+DEFPT(2)
CC        DEFP(4)=DLAM*(-E(2)*TAUD(1)-E(3)*TAUD(2)+(E(2)+E(3))*
CC     &          TAUD(4))+DEFPT(4)
        DEFP(4)=-DDEFP1-DDEFP2+DEFPT(4)
        DEFP(3)=DDEFP3+DEFPT(3)
C
      IF(IMIX.EQ.1)THEN
        ALFA1(1)=XMX(1)*DDEFP1+ALFAT(1)
        ALFA1(2)=XMX(2)*DDEFP2+ALFAT(2)
        ALFA1(4)=-ALFA1(1)-ALFA1(2)
CC        ALFA1(4)=XMX(4)*DLAM*(-E(2)*TAUD(1)-E(3)*TAUD(2)+(E(2)+E(3))*
CC     &           TAUD(4))+ALFAT(4)
        ALFA1(3)=XMX(4)*DDEFP3+ALFAT(3)
        DO 160 I=1,4
  160   TAUD(I)=TAUD(I)+ALFA1(I)
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      TAUM=CM(1)*(DEFL(1)-DEFP(1))+CM(2)*(DEFL(2)-DEFP(2))+
     &     CM(3)*(DEFL(4)-DEFP(4))
      TAUL(1)=TAUD(1)+TAUM
      TAUL(2)=TAUD(2)+TAUM
      TAUL(4)=TAUD(4)+TAUM
      TAUL(3)=TAUD(3)
      DEFP(3)=2.*DEFP(3)
C
C... TRANSFORM STRESS INTO GLOBAL AXES DIRECTION
C
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL CLEAR(TAU,4)
         CALL MNOZI2(TAU,TE,TAUL,4,4)
      ELSE
         CALL JEDNA1(TAU,TAUL,4)
      ENDIF
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
        DEFL(4)=-VZX/EZ*TAUL(1)-VYZ/EY*TAUL(2)+DEFP(4)
      ELSE
        DEFL(4)=0.D0
      ENDIF
      DO 290 I=1,4
      DEF1(I)=DEFL(I)
  290 TAU1(I)=TAU(I)
      RETURN
      END
C=======================================================================
      SUBROUTINE PEL210(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORM ELASTICITY CONSTANTS
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      DIMENSION FUN(29,*),CE(4,*),CPE(3,*),CM(*)
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
C======================================================================
      SUBROUTINE DEV210(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TAUD(*),E(*),SE(*),GB(*),XMX(*)
        DL = 1. + DLAM*(CXX+CYY + DLAM*CED)
        TAUD(1) = (SE(1) + DLAM*DXX)/DL
        TAUD(2) = (SE(2) + DLAM*DYY)/DL
        TAUD(4) = - TAUD(1) - TAUD(2)
        TAUD(3) = SE(3)/(1.+DLAM*(GB(1)+XMX(4)*E(4)))
      RETURN
      END      
C======================================================================
      SUBROUTINE INI210(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX,
     &                  IETYP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SE(*),E(*),CE(4,*),CPE(3,*),GB(*),XMX(*)
C
      DVA=2.D0
      XX1 = E(1)+DVA*E(2)
      XX2 = E(1)-E(3)
      XX3 = 0.D0
      IF(IETYP.NE.0.AND.IETYP.NE.3) XX3 = DVA*E(2)+E(3)
      CXX =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2) - XX3*CPE(1,3)
      CYX = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2)) + XX3*CPE(2,3)
      XX1 = E(1)-E(2)
      XX2 = E(1)+DVA*E(3)
      IF(IETYP.NE.0.AND.IETYP.NE.3) XX3 = E(2)+DVA*E(3)
      CYY = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2)) - XX3*CPE(2,3)
      CXY =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2) + XX3*CPE(1,3)
      CED = CXX*CYY - CXY*CYX
      DXX = CYY*SE(1)+CXY*SE(2)
      DYY = CYX*SE(1)+CXX*SE(2)
      GB(1)=CE(3,3)*DVA*E(4)
      RETURN
      END
C======================================================================
      SUBROUTINE CEP210(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                  DL,DXX,DYY,DLAM,TAUY,XMX,EP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C*****  NIJE PRILAGODJENO ZA MESOVITO OJACANJE
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      DIMENSION E(*),CE(4,*),CPE(3,*),TAUD(*),CM(*),GB(*),XMX(*)
      DIMENSION T(4),DDL(4),XEN(4,4),ETP(4,4)
C
      ZERO=0.D0
      ONE =1.D0
      PO  =0.5D0
      DVA =2.D0
      DVT =DVA/3.
      ELP = (1.5 - EP*DLAM)/TAUY
C
      DS12 = TAUD(1) - TAUD(2)
      DS13 = TAUD(1) - TAUD(4)
      DS23 = TAUD(2) - TAUD(4)
      AXX = E(1)*DS12 + DVA*E(2)*DS13 + E(3)*DS23
      AYY = DVA*E(3)*DS23 - E(1)*DS12 + E(2)*DS13
C
      P11 = (CPE(1,1) + DLAM*(CYY*CPE(1,1) + CXY*CPE(2,1)))/DL
      P12 = (CPE(1,2) + DLAM*(CYY*CPE(1,2) + CXY*CPE(2,2)))/DL
      P21 = (CPE(2,1) + DLAM*(CXX*CPE(2,1) + CYX*CPE(1,1)))/DL
      P22 = (CPE(2,2) + DLAM*(CXX*CPE(2,2) + CYX*CPE(1,2)))/DL
      P31 = -P11-P21
      P32 = -P12-P22
      DLL = CXX + CYY + DVA*DLAM*CED
      Q1 = (DXX - DLL*TAUD(1))/DL
      Q2 = (DYY - DLL*TAUD(2))/DL
      Q3 = -Q1-Q2
C
      DUM  = 1. + DLAM*GB(1)
      R4= DVA*CE(3,3)/DUM
      T4= GB(1)*TAUD(3)/DUM
      W1 = AXX*P11 + AYY*P21
      W2 = AXX*P12 + AYY*P22
      W3 = DVA*E(4)*TAUD(3)*R4
      W0 = ELP*(AXX*Q1 + AYY*Q2 - DVA*E(4)*TAUD(3)*T4) - DVT*EP*TAUY
      DUM  =  -ELP/W0
      DDL(1) = DUM*W1
      DDL(2) = DUM*W2
      DDL(3) = DUM*W3
C
      ELAST(1,1) = P11 + Q1*DDL(1)
      ELAST(1,2) = P12 + Q1*DDL(2)
      ELAST(2,1) = P21 + Q2*DDL(1)
      ELAST(2,2) = P22 + Q2*DDL(2)
      ELAST(1,3) = Q1*DDL(3)
      ELAST(2,3) = Q2*DDL(3)
      ELAST(3,1) = - T4*DDL(1)
      ELAST(3,2) = - T4*DDL(2)
      ELAST(3,3) = R4 - T4*DDL(3)
C
      IF(IETYP.NE.0.AND.IETYP.NE.3)THEN
        AXX =  (E(1)+E(2))*CM(1)-E(1)*CM(2)-E(2)*CM(3)
        AYY = -E(1)*CM(1)+(E(1)+E(3))*CM(2)-E(3)*CM(3)
        AZZ = -E(2)*CM(1)-E(3)*CM(2)+(E(2)+E(3))*CM(3)
      ELSE
        AXX = (E(1)+DVA*E(2))*CM(1)-(E(1)-E(3))*CM(2)
        AYY = (E(2)-E(1))*CM(1)+(E(1)+DVA*E(3))*CM(2)
        AZZ = 0.D0
      ENDIF
      DUM  = AXX*TAUD(1)+AYY*TAUD(2)+AZZ*TAUD(4)
      PM1 = AXX*ELAST(1,1) + AYY*ELAST(2,1)
      PM2 = AXX*ELAST(1,2) + AYY*ELAST(2,2)
      PM3 = AXX*ELAST(1,3) + AYY*ELAST(2,3)
      QM1 = DUM*DDL(1)
      QM2 = DUM*DDL(2)
      QM3 = DUM*DDL(3)
C
      CM1 = CM(1) - DLAM*PM1 - QM1
      CM2 = CM(2) - DLAM*PM2 - QM2
      CM3 =       - DLAM*PM3 - QM3
      CALL CLEAR(XEN,16)
      XEN(1,1)= XMX(1)*(E(1)+E(2))
      XEN(1,2)=-XMX(1)*E(1)
      XEN(1,4)=-XMX(1)*E(2)
      XEN(2,1)=-XMX(2)*E(1)
      XEN(2,2)= XMX(2)*(E(1)+E(3))
      XEN(2,4)=-XMX(2)*E(3)
      XEN(4,1)=-XMX(4)*E(2)
      XEN(4,2)=-XMX(4)*E(3)
      XEN(4,4)= XMX(4)*(E(2)+E(3))
      XEN(3,3)= XMX(3)*E(4)
C      CALL WRR(XEN,16,'XEN ')
C... (A) 
      T(1)=XEN(1,1)*TAUD(1)+XEN(1,2)*TAUD(2)+XEN(1,4)*TAUD(4)
      T(2)=XEN(2,1)*TAUD(1)+XEN(2,2)*TAUD(2)+XEN(2,4)*TAUD(4)
      T(3)=XEN(3,3)*TAUD(3)
C      T(4)=XEN(4,1)*TAUD(4)+XEN(4,2)*TAUD(2)+XEN(4,4)*TAUD(4)
C... (B)
      DO 78 I=1,4
      DO 77 J=1,4
   77 XEN(I,J)=DLAM*XEN(I,J)
      XEN(I,I)=XEN(I,I)+ONE
   78 CONTINUE
C
C*** ZA PLANE STRAIN TREBA DO 4 !
      DO 85 I=1,3
      DO 85 J=I,3
       DUM=ZERO
        DO 83 K=1,3
   83   DUM=DUM+XEN(I,K)*ELAST(K,J)
       ETP(I,J)=DUM+T(I)*DDL(J)
   85 CONTINUE
C
      ELAST(1,1) = ETP(1,1) + CM1
      ELAST(1,2) = ETP(1,2) + CM2
      ELAST(1,3) = PO*(ETP(1,3) + CM3)
      ELAST(2,2) = ETP(2,2) + CM2
      ELAST(2,3) = PO*(ETP(2,3) + CM3)
      ELAST(3,3) = PO*ETP(3,3)
C
      DO 90 I=1,4
      DO 90 J=I,4
        ELAST(J,I)=ELAST(I,J)
   90 CONTINUE
      RETURN
      END
