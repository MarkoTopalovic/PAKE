C=======================================================================
C
C   PLASTICNOST 2/D ELEMENT  -  MESOVITO OJACANJE    (01.09.1992)
C
C=======================================================================
      SUBROUTINE D2M6G(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE  SPACE IN WORKING VECTOR 
C
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /DDUPLAP/ LPLAS,LPLA1
C
      LPLAS=LPOCG
      LPLA1=LPOC1
      LFUN=MREPER(1)
      MATE=MREPER(4)
      LTAU  =LPOCG
      LDEFT =LTAU   + 4*IDVA
      LDEFPT=LDEFT  + 4*IDVA
      LFPO=LDEFPT + 4*IDVA
      LTEQT =LFPO + 4*IDVA
      LDQPT =LTEQT  + 1*IDVA
      LIPL  =LDQPT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LFPO1=LDEFP1 + 4*IDVA
      LTEQT1=Lfpo1 + 4*IDVA
      LDQPT1=LTEQT1 + 1*IDVA
      LIPL1 =LDQPT1 + 1*IDVA
C
      CALL TAUI26G(A(LIPL),A(LDEFPT),A(LFPO),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LFPO1),
     1            A(LTEQT1),A(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUI26G( PL ,DEFPT,DFPO,TEQT,DEFQPT,
     1                  PL1,TAU1,DEF1,DEFP, DFPO1, TEQ,DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
      
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA ELASTOPLASTIC
C     ELASTOPLASTICAN MATERIJAL SA IZOTROPNIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT,TKOR
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR

      COMMON /FIER/ FIER(3,3)
      COMMON /CEPMAT/ INDCEP
      COMMON /DDUPLAP/ LPLAS,LPLA1
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),SHET(4),POM(3,3)
      DIMENSION FUN(2,MATE,*),SII(4),GLITL(4)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      KSS=4
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) KSS=3
      DO 387 I=5,6
         TAU(I)=0.
  387 DEF(I)=0.
      IPL =PL
      IPL1=PL1
      PL11=0.
      DVT =2.D0/3.
C
      E     =FUN(1,MAT,1)
      ANI   =FUN(2,MAT,1)
      TEQY0 =FUN(1,MAT,2)
      CY    =FUN(2,MAT,2)
      AN    =FUN(1,MAT,3)
      EM    =FUN(2,MAT,3)
      DPOR  =FUN(2,MAT,5)
      Dff   =FUN(2,MAT,6)
      Dq1   =FUN(2,MAT,7)
      Dq2   =FUN(2,MAT,8)
      Dq3   =FUN(2,MAT,9)
      Dfc   =FUN(2,MAT,10)
      IF(kor.Eq.1)DFPO=dpor
c      IF(VREME.LE.1)DFPO=dpor
      CZZ   =1.-ANI
      C1    =(2.-ANI)/3./CZZ
      CZZ   =-ANI/CZZ
      CNI   =1.-C1
      CYAN23=DVT*CY*AN
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-2.
      AE    =(1.+ANI)/E
      AEINV =1./AE
      CM    =E/(1.-2.*ANI)
      DGT=E/(2.*(1.+ANI))
      dk=(1/dq1-dfc)/(dff-dfc)
C
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
      CALL MEL26G(FUN,MATE)
      IF(IRAC.EQ.2)RETURN
C      IF(ITER.EQ.0)RETURN
      DO I=1,3
      DO J=1,3
         POM(I,J)=0.
      ENDDO
      ENDDO
C
C...  TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFPT(3)=0.5*DEFPT(3)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
C
      EMPT = (DEFPT(1)+DEFPT(2)+DEFPT(4))/3.
      IF(IATYP.LT.4)THEN
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2))
          DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNI*(DEF(2)-DEFPT(2))
          DEFDS(2)=-CNI*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))
        ELSE
          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT-DEFPT(1)
          DEFDS(2)=DEF(2)-EMT-DEFPT(2)
          DEFDS(4)=DEF(4)-EMT-DEFPT(4)
        ENDIF
       DEFDS(3)=0.5*DEF(3)-DEFPT(3)
      ELSE
        IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
          EMT = CNI*(DEF(1)+DEF(2))
          DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
          DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        ELSE
          EMT=(detg-1.D0)/3.D0
C          EMT = (DEF(1)+DEF(2)+DEF(4))/3.
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(4)=DEF(4)-EMT
        DEFDS(3)=0.5*DEF(3)
        ENDIF
        DEFDS(3)=0.5*DEF(3)
      ENDIF
      if(IATYP.EQ.5) THEN
      EMT=DLOG(detg)/3.D0
	      DEFDS(1)=DEF(1)
            DEFDS(2)=DEF(2)
            DEFDS(4)=DEF(4)
        DEFDS(3)=0.5*DEF(3)

C      EMT=DABS(EMT)
       ENDIF
       TAUM=CM*EMT
c       IF(TAUM.LT.1.D-07)then
c       CALL D2M6(TAU,DEF,IRAC,LPLAS,LPLA1)
C	dfpo=0.
c       go to 339
c       endif
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      DO 40 I=1,KSS
      GLITL(I)=DEFDS(I)
      TAUD(I) =DEFDS(I)*AEINV
   40  SHET(I) =TAUD(I)
      IF(IETYP.EQ.0.OR.IETYP.EQ.3)THEN
        GLITL(4)=-GLITL(1)-GLITL(2)
        TAUD(4) =-TAUD(1)-TAUD(2)
        SHET(4) =-SHET(1)-SHET(2)
      ENDIF
      TEQ=DSQRT(1.5*TENDO2(SHET))
C			
CE   2)  CHECK FOR YIELDING
C
      IF((TEQ-(TEQY))/(TEQY).LT.1.D-5)THEN
        TEQ  =TEQY
        DEFQP=DEFQPT
         DO 600 I=1,4
  600    DEFP(I)=DEFPT(I)
	DFPOO=DFPO
	  TAUM=CM*EMT
        GO TO 500
      ENDIF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
      PL11=1.0D0
C
CE       ITERATIONS
C
      DEFQP=DEFQPT
      IF(DEFQP.LE.1.D-4) DEFQP=1.D-4
        EP2=AN*CY*DEFQP**AN1      
      IB = 0
      IT = 0
      IT1 = 0
      IBB=0 
      AF = 1.5
      AF1 =1.5
CE   POCETNI SREDNJI PRIRASTAJ PLASTI. DEFORM (DDEMMP)
        DDD=0.0001*(TEQ-TEQY)/EP2
	DDEMMP=0.1*DDD*DFPO*DSINH((CM*EMT/TEQY))/(3*TEQY)
      FM1    = 0.D0
      DEPBM1 = 0.D0
      DEPBP1 = 0.D0
      DDEMP= DDEMMP
        TAUY  = TEQT
  109    IT1=IT1+1
       IBB1=IBB
      IF(IT1.GT.ITMAX) THEN
        WRITE(IZLAZ,2000)
        STOP
      ENDIF
     	 SHITM=CM*(EMT-DDEMP)
 	DDFPO=(3*(1.-DFPO)*DABS(DDEMP)/(1.+3*DABS(DDEMP)))
	DFPOO=DFPO+DDFPO
	IF(DFPOO.GT.DFC)DFPOO=DFC+DK*(DFPOO-DFC)
CE  SREDNJI NAPON (SIGMAM)
CE       POROZNOST
       SH=DSINH((1.5*SHITM/TEQY))
       DL=3*DDEMP/(DQ1*DQ2*DFPOO*SH)/TEQY
       DCH=DCOSH(1.5*SHITM/TEQY)
C ODREDIVANJE KONSTANTE BISEKCIJE (FP)
       DO 345 I=1,4
  345     SII(I)=TAUD(I)/(1+2*DGT*DL)
       SQQ=TENDO2(SII)
       FP=DL*SQQ+3*SHITM*DDEMP
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      IF(DDD.LE.1.D-4) DDD=1.D-4
      DDEFQP=DDD
       IB=0
	 IT=0
  100 IT=IT+1
      IB1 = IB
C_____________________________________________________________
C ODREDIVANJE PROMENLJIVE BISEKCIJE (FB)
C_____________________________________________________________
      IF(IT.GT.ITMAX) THEN
        WRITE(IZLAZ,2000)
	write(izlaz,*)fb,fb1,ddemp,ddefqp,ne
        STOP
      ENDIF
C
      DEFQP=DEFQPT+DDEFQP
        TEQY=TEQY0+CY*(EM*DEFQP)**AN
       SH=DSINH((1.5*SHITM/TEQY))
       DL=3*DDEMP/(DQ1*DQ2*DFPOO*SH)/TEQY
       DCH=DCOSH(1.5*SHITM/TEQY)
       DO 346 I=1,4
  346     SII(I)=TAUD(I)/(1+2*DGT*DL)
       SQQ=TENDO2(SII)
       FB=(DL*SQQ+3*SHITM*DDEMP-(1-DFPOO)*TEQY*DDEFQP)
C____________________________________________________________
C UNUTRASNJA BISEKCIJA
C_____________________________________________________________
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF ((DABS(DDD)/(DEPBM+DEPBP)).GT.1.D-10)GO TO 100
       DDD=DDEFQP
C______________________________________________________________________
C_____________________________________________________________
C ODREDIVANJE KONSTANTE BISEKCIJE (FP1)
C_____________________________________________________________
C_____________________________________________________________
C_____________________________________________________________
C ODREDIVANJE PROMENLJIVE BISEKCIJE (FB1)
C_____________________________________________________________
      DEFQP=DEFQPT+DDEFQP
        TEQY=TEQY0+CY*(EM*DEFQP)**AN
       IF(IT1.EQ.1)THEN
     	 SHITM=CM*(EMT-EMPT)
       DCH=DCOSH((1.5*SHITM/TEQY))
       DO 347 I=1,4
  347     SII(I)=TAUD(I)
       SQQQ=TENDO2(SII)
       FP1=+0.5*SQQQ+(2*DQ1*DQ2*DFPO*DCH-1-(DQ3*DFPO)**2)*(TEQY**2)/3
       ENDIF
     	 SHITM=CM*(EMT-DDEMP)
       SH=DSINH((1.5*SHITM/TEQY))
       DCH=DCOSH((1.5*SHITM/TEQY))
 	DDFPO=(3*(1.-DFPO)*DABS(DDEMP)/(1.+3*DABS(DDEMP)))
	DFPOO=DFPO+DDFPO
       DL=3*DDEMP/(DQ1*DQ2*DFPOO*SH)/TEQY
       DO 349 I=1,4
  349     SII(I)=TAUD(I)/(1+2*DGT*DL)
       SQQ=TENDO2(SII)
       FB1=0.5*SQQ+(2*DQ1*DQ2*DFPOO*DCH-1-(DQ3*DFPOO)**2)*(TEQY**2)/3.
      CALL BISEC(DDEmp,DEPBM1,DEPBP1,DDemmp,FB1,FM1,FP1,AF1,IBB)
      IF (IBB1.EQ.0) GO TO 109
      IF ((DABS(DDEMMP)/(DEPBM1+DEPBP1)).GT.1.D-10) GO TO 109
C____________________________________________________________
C UNUTRASNJA BISEKCIJA
C_____________________________________________________________
C	 DDEMMP=DDEMP
      IF (DABS(FB).GT.1.D-7)GO TO 109
C
CE      ...   ( DEVIATORIC STRESS )
        TEQ =TEQY0+CY*(EM*DEFQP)**AN
        DO 165 I=1,4
        SHET(I)=SII(I)
  165   TAUD(I)=SHET(I)
 2000 FORMAT(' ','DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI36')
C
CE   4)  DETERMINE SOLUTION 
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
	IF(IATYP.EQ.5) THEN
	  CALL DIJAD(POM,FIER,FIER,SII(1),SII(2),SII(4))
        SII(1)=POM(1,1)
        SII(2)=POM(2,2)
        SII(4)=POM(3,3)
        SII(3)=POM(1,2)
       ENDIF
	DEFP(1)=DEFPT(1)+DL*SII(1)+DDEMP
	DEFP(2)=DEFPT(2)+DL*SII(2)+DDEMP
	DEFP(4)=DEFPT(4)+DL*SII(4)+DDEMP
	DEFP(3)=DEFPT(3)+DL*SII(3)
C
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
	DFPO1=DFPOO
      IF(INDCEP.EQ.0)
     1   CALL CEP26G(TEQ,DEFQP,DDEFQP,EM,AE,CM,CY,AN,
     1       DL,DFPO,DGT,DDEMP,SHET,SHITM,DQ1,DQ2,DQ3,DFPOO,GLITL,AEINV)
C     1       DL,DFPO,DGT,DDEMP,TAUD,SHITM,DQ1,DQ2,DQ3,DFPOO,GLITL,AEINV)
C
CE   5)    CALCULATE STRESS
C
C	  TAUM=CM*EMT
		TAUM=SHITM
	IF(TAUM.GT.1D-5)THEN
	ENDIF
  500 CONTINUE
	DFPO1=DFPOO
      IF(IATYP.lt.4)THEN
      ELSEIF(IATYP.EQ.4) THEN
        DO 210 I=1,2
  210   DEF(I)=TAUD(I)*AE+TAUM/CM
        DEF(4)=TAUD(4)*AE+TAUM/CM
        DEF(3)=TAUD(3)*AE*2.D0
      ELSEIF(IATYP.EQ.5) THEN
C SREDNJI NAPON
C NORMIRANA ELASTICNA DEFORMACIJA
	aee=ae
C        DO 230 I=1,2
        DEF(2)=TAUD(2)*AEE
        DEF(1)=TAUD(1)*AEE
        DEF(4)=TAUD(4)*AEE
	def(3)=2.*taud(3)*ae
      ENDIF
      DETGG=1.
      taum=taum
      TAU(1)=(TAUD(1)+TAUM)*DETGG
      TAU(2)=(TAUD(2)+TAUM)*DETGG
      TAU(4)=(TAUD(4)+TAUM)*DETGG
      TAU(3)=TAUD(3)*DETGG
      DEFP(3)=2.*DEFP(3)
C
C
CE  UPDATE FROM PREVIOUS STEP
C
      CALL JEDNA1(SHET,DEF,4)
      IF(IATYP.EQ.5) THEN
	  CALL DIJAD(POM,FIER,FIER,TAU(1),TAU(2),TAU(4))
        TAU(1)=POM(1,1)
        TAU(2)=POM(2,2)
        TAU(4)=POM(3,3)
        TAU(3)=POM(1,2)
C	  PRINT*,POM(1,1),POM(2,2),POM(3,3),POM(1,2)
         CALL DIJAD(POM,FIER,FIER,DEF(1),DEF(2),DEF(4))
        DEF(1)=POM(1,1)
        DEF(2)=POM(2,2)
        DEF(4)=POM(3,3)
        DEF(3)=POM(1,2)
      ENDIF
C	  PRINT*,DEF(1),DEF(2),DEF(3),DEF(4)
      DO 290 I=1,4
        DEF1(I)=DEF(I)
  290      TAU1(I)=TAU(I)
      IF(IATYP.EQ.4)THEN
C KORIGOVANJE NORMIRANOG ELASTICNOG DEFORMACIONOG TENZORA NA KRAJU KORAKA
        DO 300 I=1,2
  300   DEF(I)=2.D0*DEF(I)+1.D0
        DEF(4)=2.D0*DEF(4)+1.D0
        POM(1,1)=DEF(1)
        POM(2,2)=DEF(2)
        POM(3,3)=DEF(4)
        POM(1,2)=DEF(3)
      ELSEIF(IATYP.EQ.5) THEN
        DO 310 I=1,2
  310   DEF(I)=DEXP(2.D0*SHET(I))
        DEF(4)=DEXP(2.D0*SHET(4))
        CALL DIJAD(POM,fier,fier,DEF(1),DEF(2),DEF(4))
        DEF(1)=POM(1,1)
        DEF(2)=POM(2,2)
        DEF(3)=POM(3,3)
        DEF(4)=POM(1,2)
C        CALL DETER3(POM,DETB)
	endif
 339     RETURN
      END
C=======================================================================
      SUBROUTINE MEL26G(FUN,MATE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION FUN(2,MATE,*)
C
      E=FUN(1,MAT,1)
      V=FUN(2,MAT,1)
      DO 15 I=1,4
      DO 15 J=I,4
   15 ELAST(I,J)=0.0D0
CS     RAVANSKO STANJE NAPONA
CE     PLANE STRESS
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        ELAST(1,1)=E/(1.-V*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V
        ELAST(3,3)=ELAST(1,1)*(1.-V)/2.
      ELSE
CS     RAVANSKA DEFORMACIJA
CE     PLANE STRAIN
        ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V/(1.-V)
        ELAST(3,3)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
        IF(IETYP.NE.1) GO TO 40
        ELAST(4,4)=ELAST(1,1)
        ELAST(1,4)=ELAST(1,2)
        ELAST(2,4)=ELAST(1,2)
      ENDIF
C
   40 DO 50 I=1,4
      DO 50 J=I,4
   50 ELAST(J,I)=ELAST(I,J)
      RETURN
      END
C======================================================================
      SUBROUTINE CEP26G(TEQ,DEFQP,DDEFQP,EM,AE,CM,CY,AN,
     &    DL,DFPO,DGT,DDEMP,SI,SHITM,DQ1,DQ2,DQ3,DFPOO,GLITL,AEINV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )  2D
CE     ELASTO-PLASTIC  CEP MATRIX        2D
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION SI(4),GLITL(4) 
      DIMENSION CBC(4),BBB(4),CPC(4,4),AA(4,4),AK(16),LJA(16),MJA(16) 
      DIMENSION CPCC(4,4) 
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      DVT   =DVA*TR
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-DVA
C
      CHET  =0.
	 CBC(1)=0.
	 CBC(2)=0.
	 CBC(3)=0.
	 CBC(4)=0.
C
C-----------------------------------------------------------------
      GGG=1/(1+2*DGT*DL)**2
      GGl=1/(1+2*DGT*DL)
      GG1=(1+2*DGT*DL)**2
      GG2=(1+2*DGT*DL)
      DA0=-CM/TEQ
       DAP=(CY*AN*(DEFQP**(AN-1)))
       DA1=-SHITM*DAP/TEQ**2
      DFSM=3*(1-DFPO)/(1+3*dabs(DDEMP))**2
      DFYM=DQ1*DQ2*DFPOO*DSINH((1.5*SHITM/TEQ))*TEQ/3.
      DFYMP=DQ1*DQ2*DFPOO*DSINH((1.5*SHITM/TEQ))/3.
      AAA=DAP*DFYMP
      DFSSYM0=(DFYM*DFSM/DFPOO+1.5*TEQ*DQ1*DQ2*DFPOO*DCOSH(1.5*SHITM/
     1 TEQ)*DA0)
      DFSSYM1=AAA+1.5*TEQ*DQ1*DQ2*DFPOO*DCOSH(1.5*SHITM/TEQ)*DA1
      DL00=3*(1-DL*DFSSYM0/3)/DFYM
      DL11=-DFSSYM1*DL/DFYM
      DEPM0=(2*DQ1*DQ2*DFSM*DCOSH(1.5*SHITM/TEQ)+3*DFYMP*DA0-2*DQ3
     1 **2*DFPOO*DFSM)*(TEQ**2)/3
      DEPM1=TEQ**2*DFYMP*DA1+2*TEQ*DAP*(2*DQ1*DQ2*DFPOO*DCOSH(1.5*
     1 SHITM/TEQ)-1-(DQ3*DFPOO)**2)/3
       XX=2*DGT
C--------------------------------------------------------------
C-----------------------------------------------------------------
      DFQPM=-(-XX*TENDO2(SI)*DL00+GG2*DEPM0)/(-XX*TENDO2(SI)*DL11+GG2*
     1 DEPM1)
      DLM=DL00+DL11*DFQPM
C--------------------------------------------------------------------------
	AA(1,1)=2./3.
	AA(2,2)=2./3.
	AA(4,4)=2./3.
	AA(3,3)=1./2.
	AA(1,2)=-1./3.
	AA(1,4)=-1./3.
	AA(2,1)=-1./3.
	AA(4,1)=-1./3.
	AA(4,2)=-1./3.
	AA(2,4)=-1./3.
	AA(1,3)=0.
	AA(3,1)=0.
	AA(2,3)=0.
	AA(3,2)=0.
	AA(4,3)=0.
	AA(3,4)=0.
	AB=2*DGT/GG2
	DO 567 I=1,4
	BB=0.
	DO 577 J=1,4
 577	BB=BB+4*DGT*GGL*DL*SI(J)*AA(J,I)
 567	BBB(I)=-BB+CM*DDEMP
	AMMM=DA0+DA1*DFQPM
      DLLL=DLM*TENDO2(SI)+3*SHITM-(1-DFPOO)*DAP*DFQPM*DDEFQP-TEQ*
     1 ((1-DFPOO)*DFQPM-DDEFQP*DFSM)-3*CM*DDEMP
	BA=3./(GG2*DFYM)
	CAB=1-DQ1*DQ2*DDEMP*(DAP*DFQPM*DFPOO*DSINH((1.5*SHITM/TEQ))
     1  +TEQ*(DFSM*DSINH((1.5*SHITM/TEQ))+1.5*DCOSH(1.5*SHITM/TEQ)
     1   *AMMM))/DFYM
	CBB=CAB*BA
	WWW=-AB*DLM
      DO 28 I=1,4
        DO 29 J=I,4
        IF(I.NE.J) CPC(I,J)=-2*DL*AB*AA(I,J)*CBB*TENDO2(SI)
   29   CONTINUE
        CPC(I,I)=DLLL+2*DL*TENDO2(SI)*(WWW-AB*CBB*AA(I,I))
   28 CONTINUE
      DO 57 I=1,4
      DO 57 J=I,4
        CPC(J,I)=CPC(I,J)
   57 CONTINUE
	 K=1
	DO 347 I=1,4
	DO 348 J=1,4
 	AK(K)=CPC(J,I)
 348	K=K+1
 347	CONTINUE
	CALL MINV(AK,4,DD,LJA,MJA)
	 K=0
	DO 297 I=1,4
	DO 247 J=1,4
	 K=K+1
 247	CPCC(J,I)=AK(K)
 297	CONTINUE
	CALL MNOZI1(CBC,CPCC,BBB,4,4)
	ELAST(1,1)=AA(1,1)*AB*(1-CBB*SI(1)*CBC(1))-AB*AA(1,2)*CBB*SI(1)*
     1CBC(2)-AB*AA(1,3)*CBB*SI(1)*CBC(3)-AB*AA(1,4)*CBB*SI(1)*CBC(4)-AB
     1*SI(1)*DLM*CBC(1)+CM/3.-CM*CBC(1)
	ELAST(2,2)=-AA(1,2)*AB*CBB*SI(2)*CBC(1)+AB*AA(2,2)*(1-CBB*SI(2)*
     1CBC(2))-AB*AA(3,2)*CBB*SI(2)*CBC(3)-AB*AA(4,2)*CBB*SI(2)*CBC(4)-AB
     1*SI(2)*DLM*CBC(2)+CM/3.-CM*CBC(2)
	ELAST(3,3)=-AA(1,3)*AB*CBB*SI(3)*CBC(1)-AB*AA(2,3)*CBB*SI(3)*
     1CBC(2)+AB*AA(3,3)*(1-CBB*SI(3)*CBC(3))-AB*AA(4,3)*CBB*SI(3)*CBC
     1(4)-AB*SI(3)*DLM*CBC(3)
	ELAST(4,4)=-AA(1,4)*AB*CBB*SI(4)*CBC(1)-AB*AA(2,4)*CBB*SI(4)*
     1CBC(2)-AB*AA(3,4)*CBB*SI(4)*CBC(3)+AB*AA(4,4)*(1-CBB*SI(4)*CBC
     1 (4))-AB*SI(4)*DLM*CBC(4)+CM/3.-CM*CBC(4)
	ELAST(1,2)=AA(1,2)*AB*(1-CBB*SI(1)*CBC(1))-AB*AA(2,2)*CBB*SI(1)*
     1CBC(2)-AB*AA(3,2)*CBB*SI(1)*CBC(3)-AB*AA(4,2)*CBB*SI(1)*CBC(4)-AB
     1*SI(1)*DLM*CBC(2)+CM/3.-CM*CBC(2)
	ELAST(2,1)=-AA(1,1)*AB*CBB*SI(2)*CBC(1)+AB*AA(2,1)*(1-CBB*SI(2)*
     1CBC(2))-AB*AA(3,1)*CBB*SI(2)*CBC(3)-AB*AA(4,1)*CBB*SI(2)*CBC(4)-AB
     1*SI(2)*DLM*CBC(1)+CM/3.-CM*CBC(1)
	ELAST(1,2)=(ELAST(1,2)+ELAST(2,1))/2
	ELAST(1,3)=AA(1,3)*AB*(1-CBB*SI(1)*CBC(1))-AB*AA(2,3)*CBB*SI(1)*
     1CBC(2)-AB*AA(3,3)*CBB*SI(1)*CBC(3)-AB*AA(4,3)*CBB*SI(1)*CBC(4)-AB
     1*SI(1)*DLM*CBC(3)-CM*CBC(3)
	ELAST(1,4)=AA(1,4)*AB*(1-CBB*SI(1)*CBC(1))-AB*AA(2,4)*CBB*SI(1)*
     1CBC(2)-AB*AA(3,4)*CBB*SI(1)*CBC(3)-AB*AA(4,4)*CBB*SI(1)*CBC(4)-AB
     1*SI(1)*DLM*CBC(4)+CM/3.-CM*CBC(4)
	ELAST(4,1)=-AA(1,1)*AB*CBB*SI(4)*CBC(1)-AB*AA(2,1)*CBB*SI(4)*
     1CBC(2)-AB*AA(3,1)*CBB*SI(4)*CBC(3)+AB*AA(4,1)*(1-CBB*SI(4)*CBC
     1 (4))-AB*SI(4)*DLM*CBC(1)+CM/3.-CM*CBC(1)
	ELAST(1,4)=(ELAST(1,4)+ELAST(4,1))/2
	ELAST(2,3)=-AA(1,3)*AB*CBB*SI(2)*CBC(1)+AB*AA(2,3)*(1-CBB*SI(2)*
     1CBC(2))-AB*AA(3,3)*CBB*SI(2)*CBC(3)-AB*AA(4,3)*CBB*SI(2)*CBC(4)-AB
     1*SI(2)*DLM*CBC(3)-CM*CBC(3)
	ELAST(3,2)=-AA(1,2)*AB*CBB*SI(3)*CBC(1)-AB*AA(2,2)*CBB*SI(3)*
     1CBC(2)+AB*AA(3,2)*(1-CBB*SI(3)*CBC(3))-AB*AA(4,2)*CBB*SI(3)*CBC
     1(4)-AB*SI(3)*DLM*CBC(2)-CM*CBC(2)
	ELAST(2,3)=(ELAST(3,2)+ELAST(2,3))/2.
	ELAST(2,4)=-AA(1,4)*AB*CBB*SI(2)*CBC(1)+AB*AA(2,4)*(1-CBB*SI(2)*
     1CBC(2))-AB*AA(3,4)*CBB*SI(2)*CBC(3)-AB*AA(4,4)*CBB*SI(2)*CBC(4)-AB
     1*SI(2)*DLM*CBC(4)+CM/3.-CM*CBC(4)
	ELAST(4,2)=-AA(1,2)*AB*CBB*SI(4)*CBC(1)-AB*AA(2,2)*CBB*SI(4)*
     1CBC(2)-AB*AA(3,2)*CBB*SI(4)*CBC(3)+AB*AA(4,2)*(1-CBB*SI(4)*CBC
     1 (4))-AB*SI(4)*DLM*CBC(2)+CM/3.-CM*CBC(2)
	ELAST(2,4)=(ELAST(2,4)+ELAST(4,2))/2.
	ELAST(3,4)=-AA(1,4)*AB*CBB*SI(3)*CBC(1)-AB*AA(2,4)*CBB*SI(3)*
     1CBC(2)+AB*AA(3,4)*(1-CBB*SI(3)*CBC(3))-AB*AA(4,4)*CBB*SI(3)*CBC
     1(4)-AB*SI(3)*DLM*CBC(4)-CM*CBC(4)
	ELAST(4,3)=-AA(1,3)*AB*CBB*SI(4)*CBC(1)-AB*AA(2,3)*CBB*SI(4)*
     1CBC(2)-AB*AA(3,3)*CBB*SI(4)*CBC(3)+AB*AA(4,3)*(1-CBB*SI(4)*CBC
     1 (4))-AB*SI(4)*DLM*CBC(3)-CM*CBC(3)
	ELAST(3,4)=(ELAST(3,4)+ELAST(4,3))/2.
	ELAST(3,1)=-AA(1,1)*AB*CBB*SI(3)*CBC(1)-AB*AA(2,1)*CBB*SI(3)*
     1CBC(2)+AB*AA(3,1)*(1-CBB*SI(3)*CBC(3))-AB*AA(4,1)*CBB*SI(3)*CBC
     1(4)-AB*SI(3)*DLM*CBC(1)-CM*CBC(1)
	ELAST(1,3)=(ELAST(1,3)+ELAST(3,1))/2.
      DO 51 I=1,4
      DO 51 J=I,4
        ELAST(J,I)=ELAST(I,J)
  51  CONTINUE
      RETURN
      END

