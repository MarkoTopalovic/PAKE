C=======================================================================
C
C   SHAPE-MEMORY ALLOY 3/D ELEMENT    (14.03.2006) VLADA RANKOVIC
C
C=======================================================================
      SUBROUTINE D3M54(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
C      COMMON A(17000)
C      REAL A
      include 'paka.inc'
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
C     LMAT=MREPER(1)+(MAT-1)*20*IDVA
C     MATE=MREPER(4)
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C      
      LTAU  =LPOCG
      LDEFT =LTAU   + 6
      LDEFPT =LDEFT + 6
      LKSIT=LDEFPT  + 6
	LFAST = LKSIT + 1
	LFSAT = LFAST + 1
	LELASTT=LFSAT + 1
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1   + 6
      LDEFPT1=LDEFT1 + 6
      LKSIT1=LDEFPT1 + 6
	LFAST1 = LKSIT1 + 1
	LFSAT1 = LFAST1 + 1
	LELAST1=LFSAT1 + 1
C
      CALL TI3454(PLAST(LKSIT),PLAST(LFAST),PLAST(LFSAT),PLAST(LELASTT),
     1           PLAST(LTAU),PLAS1(LTAU1),PLAST(LDEFT),
     1           PLAS1(LDEFT1),PLAS1(LDEFPT1),
     1			 PLAS1(LKSIT1),PLAS1(LFAST1),PLAS1(LFSAT1),PLAS1(LELAST1),
     1			 A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI3454( dKsiT, dFasT, dFsaT,ELASTT, TAUT,
     1                   TAU1, DEFT, DEF1, DEFP1,
     1				   dKsi1, dFas1, dFsa1,ELAST1,
     1                   dMatProperties,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA ELASTOPLASTIC
C     ELASTOPLASTICAN MATERIJAL SA IZOTROPNIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /POCNAP/ IPOCET
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM      

	DIMENSION TAU(*),TAU0(6),TAUT(*),TAU1(6)
	DIMENSION DEF(*),DEFT(6),DEF1(*),DEF0(6),DEFP0(6),DEFP1(*),DEFP(6)
      DIMENSION dMatProperties(*),dLambdaElastic(3),dLambdaTotal(3)
	DIMENSION dLambdaTotalT(3), dYieldStretch(3)
      DIMENSION dLambdaElDev(3),dLambdaTotDev(3), dLambdaTotDevT(3)
	DIMENSION dTau(3), dNa(3), dDefDev(3)
      DIMENSION ELASTT(6,6),ELAST1(6,6)
      
	DIMENSION VDEF(3,3),DDEFPS(6),CT(3,3,3,3),TSG(6,6),TSS(6,6)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C	Indikator tipa (0 Nelinearno, 1 linearno)
      nType=0
C
CE  E,ANI,TEQY0,CY,AN,EM
C
C	Material properties	
      dE = dMatProperties(1)
      dNi = dMatProperties(2)
      dEpsilonL = dMatProperties(3)
      dT = dMatProperties(4)
      dTsas = dMatProperties(5)
      dTfas = dMatProperties(6)
      dTssa = dMatProperties(7)
      dTfsa = dMatProperties(8)
      dCas = dMatProperties(9)
      dCsa = dMatProperties(10)
      dSsas = dMatProperties(11)
      dSfas = dMatProperties(12)
      dSssa = dMatProperties(13)
      dSfsa = dMatProperties(14)
      dBas = dMatProperties(15)
      dBsa = dMatProperties(16)
	dAlpha = dMatProperties(17)
C
C	Elastic matrix
C	================================
      IF(IRAC.EQ.2) THEN
C	 IF(KOR.GT.1) THEN
C	  CALL JEDNA1(ELAST,ELASTT,36)
C	 ELSE
C	  CALL ELAST_TANGENT_MATRIX(dE, dNi)
	  CALL MEL3EL(ELAST,dE,dNi)
C	 ENDIF
	 RETURN
	ENDIF
C	================================
C	First iteration
C	================================
      IF(ITER.EQ.0) THEN
	 IF(KOR.GT.1) THEN
	  CALL JEDNA1(ELAST,ELASTT,36)
	 ELSE
C	  CALL ELAST_TANGENT_MATRIX(dE, dNi)
	  CALL MEL3EL(ELAST,dE,dNi)
	 ENDIF
	 CALL JEDNA1(TAU,TAUT,6)
	 RETURN
	ENDIF
C	================================
	dG = dE / (2.D0 * (1.D0 + dNi))
	dK = dE / (3.D0 * (1.D0 - 2.D0 * dNi))
	
C  999	continue

C	Total stretches
	DO I=1,3
	 dLambdaTotal(I) = DEXP(DEF(I))
      ENDDO

	dJ = dLambdaTotal(1)*dLambdaTotal(2)*dLambdaTotal(3)
	dJ3 = dJ**(-1.D0/3)

C	Total deviatoric stretches
	DO I=1,3
	 dLambdaTotDev(I) = dLambdaTotal(I)*dJ3
      ENDDO
C	
	DO I=1,3
	 dDefDev(I) = DLOG(dLambdaTotDev(I))
      ENDDO
      
	dDefDevNorm=dDefDev(1)*dDefDev(1)+dDefDev(2)*dDefDev(2)
	1           +dDefDev(3)*dDefDev(3)
	dDefDevNorm=DSQRT(dDefDevNorm)

	dDefNormTol = 1.D-8
	IF(dDefDevNorm.LT.dDefNormTol) THEN
	 CALL CLEAR(TAU,6)
	 CALL CLEAR(TAU1,6)
	 RETURN
	ENDIF
	DO I=1,3
	 dNa(I) = dDefDev(I)/dDefDevNorm
      ENDDO

C	Elastic deviatoric stretches
	DO I=1,3
	 dLambdaElDev(I)=dLambdaTotDev(I)/DEXP(dEpsilonL*dKsiT*dNa(I))
      ENDDO
	
	dJe = dJ/DEXP(3.D0*dAlpha*dEpsilonL*dKsiT)

C	Pressure
	dP = dK * DLOG(dJe)
C	Deviatoric stresses
	DO I=1,3
	 dTau(I) = 2.D0 * dG * DLOG(dLambdaElDev(I))
      ENDDO

C	Calculate functions Fas and Fsa from previous step
	CALL CALCULATE_FAS(dTau,dP,dMatProperties,dFas)
	CALL CALCULATE_FSA(dTau,dP,dMatProperties,dFsa)

C	Check phase transformation
	CALL CALCULATE_RSAS(dMatProperties,dRsas)
	CALL CALCULATE_RFAS(dMatProperties,dRfas)
	CALL CALCULATE_RSSA(dMatProperties,dRssa)
	CALL CALCULATE_RFSA(dMatProperties,dRfsa)
	
	dFsas = dFas - dRsas
	dFfas = dFas - dRfas
	dFssa = dFsa - dRssa
	dFfsa = dFsa - dRfsa

	dDeltaFas = dFas - dFasT
	dDeltaFsa = dFsa - dFsaT

	nHas = 0.0
	nHsa = 0.0
C	Check phase transformation
	IF(dFsas.GT.0.D0.AND.dFfas.LT.0.D0.AND.dDeltaFas.GT.0.D0) THEN
	 nHas = 1.0
	ENDIF

	IF(dFssa.LT.0.D0.AND.dFfsa.GT.0.D0.AND.dDeltaFsa.LT.0.D0) THEN
	 nHsa = 1.0
	ENDIF
	
	dLambdaAS = 0.D0
	dLambdaSA = 0.D0
	dKsiNew = dKsiT
CCC	dRASTol = 1.E-5
CCC	dRSATol = 1.E-5
	dRASTol = 1.E+3
	dRSATol = 1.E+3
	dKsiTol = 1.E-8
	IF(nHas.GT.0.D0.OR.nHsa.GT.0.D0) THEN

C	Korekcija granicnih vrednosti
	IF(DABS(dKsiT).LT.dKsiTol) dFasT=dRsas
	IF(DABS(dKsiT).GT.(1.0-dKsiTol)) dFsaT=dRssa
	
C	  Calculate functions Ras and Rsa in previous iteration
	 IF(nType.EQ.0) THEN
	   dRas=dFfas*dFfas*dLambdaAS
	1       -nHas*dBas*(1-dKsiT)*(dFas-dFasT)
	   dRsa=dFfsa*dFfsa*dLambdaSA
	1      -nHsa*dBsa*dKsiT*(dFsa-dFsaT)
	 ELSE
	   dRas=dFfas*dLambdaAS
	1       +nHas*(1-dKsiT)*(dFas-dFasT)
	   dRsa=dFfsa*dLambdaSA
	1      -nHsa*dKsiT*(dFsa-dFsaT)
	 ENDIF


C	 Iterations

C      Brojac
	 nCounter = 0
	 nMaxCount = 10
	 DO WHILE((DABS(dRas).GT.dRASTol.OR.
	1	      DABS(dRsa).GT.dRSATol).AND.nCounter.LT.nMaxCount)
	
	  dKsiNew = dKsiT + (dLambdaAS + dLambdaSA)

C	  Elastic stretches
	  DO I=1,3
	   dLambdaElDev(I)=dLambdaTotDev(I)
	1				    /DEXP(dEpsilonL*dKsiNew*dNa(I))
	  ENDDO

	  dJe = dJ/DEXP(3*dAlpha*dEpsilonL*dKsiNew)

C	  Pressure
	  dP = dK * DLOG(dJe)
C	  Deviatoric stresses
	  DO I=1,3
	   dTau(I) = 2.D0 * dG * DLOG(dLambdaElDev(I))
	  ENDDO

C	  Calculate functions Fas and Fsa from previous step
	  CALL CALCULATE_FAS(dTau,dP,dMatProperties,dFas)
	  CALL CALCULATE_FSA(dTau,dP,dMatProperties,dFsa)
C	  dFfas = dFas - dRfas
C	  dFfsa = dFsa - dRfsa

	  dG1 = (2.D0*dG + 9.D0*dK*dAlpha*dAlpha)*dEpsilonL

C	  Elements of matrix R
	  IF(nType.EQ.0) THEN
	   dA11=nHas*dBas*(dFas-dFasT)+dG1*(nHas*dBas*(1-dKsiNew)
	1        -2*dLambdaAS*dFfas)+dFfas*dFfas
	   dA12=nHas*dBas*(dFas-dFasT)+dG1*(nHas*dBas*(1-dKsiNew)
	1        -2*dLambdaAS*dFfas)
	   dA21=-nHsa*dBsa*(dFsa-dFsaT)
	1        +dG1*(nHsa*dBsa*dKsiNew-2*dLambdaSA*dFfsa)
	   dA22=-nHsa*dBsa*(dFsa-dFsaT)
	1        +dG1*(nHsa*dBsa*dKsiNew-2*dLambdaSA*dFfsa)+dFfsa*dFfsa
c	   dA11=dBas*(dFas-dFasT)+dG1*(dBas*(1-dKsiNew)
c	1        -2*dLambdaAS*dFfas)+dFfas*dFfas
c	   dA12=dBas*(dFas-dFasT)+dG1*(dBas*(1-dKsiNew)
c	1        -2*dLambdaAS*dFfas)
c	   dA21=-dBsa*(dFsa-dFsaT)
c	1        +dG1*(dBsa*dKsiNew-2*dLambdaSA*dFfsa)
c	   dA22=-dBsa*(dFsa-dFsaT)
c	1        +dG1*(dBsa*dKsiNew-2*dLambdaSA*dFfsa)+dFfsa*dFfsa
	  ELSE
	   dA11=-dG1*dLambdaAS-nHas*((dFas-dFasT)+(1-dKsiNew)*dG1)+dFfas
	   dA12=-dG1*dLambdaAS-nHas*((dFas-dFasT)+(1-dKsiNew)*dG1)
	   dA21=-dG1*dLambdaSA-nHsa*((dFsa-dFsaT)-dKsiNew*dG1)
	   dA22=-dG1*dLambdaSA-nHsa*((dFsa-dFsaT)-dKsiNew*dG1)+dFfsa
	  ENDIF


C	  Matrix R determinant	  
	  dDetRInv = dA11*dA22-dA21*dA12
C	  Matrix R exp(-1)
	  dDetRInv = 1.D0/dDetRInv

C	  Inverse matrix R elements
	  dA11Inv = nHas*dDetRInv * dA22
	  dBB = dA11Inv
	  dA12Inv = nHas*nHsa*dDetRInv * (-dA12)
	  dCC = dA12Inv
	  dA21Inv = nHsa*nHas*dDetRInv * (-dA21)
	  dDD = dA21Inv
	  dA22Inv = nHsa*dDetRInv * dA11
	  dEE = dA22Inv
c	  dA11Inv = dDetRInv * dA22
c	  dBB = dA11Inv
c	  dA12Inv = dDetRInv * (-dA12)
c	  dCC = dA12Inv
c	  dA21Inv = dDetRInv * (-dA21)
c	  dDD = dA21Inv
c	  dA22Inv = dDetRInv * dA11
c	  dEE = dA22Inv

C	  Save lambda from previous iteration
	  dLambdaASPrevious = dLambdaAS
	  dLambdaSAPrevious = dLambdaSA

C	  Calculate functions Ras and Rsa in previous iteration
	  IF(nType.EQ.0) THEN
	   dRas=dFfas*dFfas*dLambdaASPrevious
	1       -nHas*dBas*(1-dKsiNew)*(dFas-dFasT)
	   dRsa=dFfsa*dFfsa*dLambdaSAPrevious
	1       -nHsa*dBsa*dKsiNew*(dFsa-dFsaT)
	  ELSE
	   dRas=dFfas*dLambdaASPrevious
	1       +nHas*(1-dKsiNew)*(dFas-dFasT)
	   dRsa=dFfsa*dLambdaSAPrevious
	1       -nHsa*dKsiNew*(dFsa-dFsaT)
	  ENDIF


C	  New lambda
	  dLambdaAS = dLambdaASPrevious
	1			- (dA11Inv * dRas + dA12Inv * dRsa)
	  dLambdaSA = dLambdaSAPrevious
	1			- (dA21Inv * dRas + dA22Inv * dRsa)

C	  Ako nema fazne transformacije onda ne bi trebalo da se menja lambda te fazne transformacije		  
c	  dLambdaAS = nHas * dLambdaAS
c	  dLambdaSA = nHsa * dLambdaSA
	  
	  nCounter = nCounter + 1

	  ENDDO

C	 Tangent matrix
c	 dAas=-2.D0*dLambdaAS*dFfas+dBas*(1.D0-dKsiNew)
c	 dAsa=-2.D0*dLambdaSA*dFfsa+dBsa*dKsiNew
	 dAas=-2.D0*dLambdaAS*dFfas+nHas*dBas*(1.D0-dKsiNew)
	 dAsa=-2.D0*dLambdaSA*dFfsa+nHsa*dBsa*dKsiNew

	 dT1as=2.D0*dG*(dBB*dAas+dCC*dAsa)
	 dT2as=3.D0*dK*dAlpha*(dBB*dAas+dCC*dAsa)
	 dT1sa=2.D0*dG*(dDD*dAas+dEE*dAsa)
	 dT2sa=3.D0*dK*dAlpha*(dDD*dAas+dEE*dAsa)

c       CALL TANGENTMATRIX(dT1as,dT2as,dT1sa,dT2sa,dK,dG,dEpsilonL,
c	1                         dAlpha,dKsiNew,dDefDevNorm,dNa,dE,dNi)
       CALL TANGENTMATRIX_N(dT1as,dT2as,dT1sa,dT2sa,dK,dG,dEpsilonL,
     1                  dAlpha,dKsiNew,dDefDevNorm,dNa,dE,dNi)
	ENDIF



C	Stress Napon      
	DO I=1,3
	 TAU0(I) = dTau(I) + dP
	ENDDO
	DO I=4,6
	 TAU0(I) = 0.0
	ENDDO	

	IF(NLM.EQ.1.OR.NLM.EQ.75) THEN
C	WRITE(3,*) "Element=", NLM
C	WRITE(3,*) "TAU(1)=", TAU0(1)
C	WRITE(3,*) "TAU(2)=", TAU0(2)
C	WRITE(3,*) "TAU(3)=", TAU0(3)
C	WRITE(3,*) "DEF1=", DEF(1)
C	WRITE(3,*) "DEF2=", DEF(2)
C	WRITE(3,*) "DEF3=", DEF(3)
C	WRITE(3,*) "KsiT=", dKsiT
C	WRITE(3,*) "KsiNew=", dKsiNew
C	WRITE(3,*) "LambdaAS=", dLambdaAS
C	WRITE(3,*) "LambdaSA=", dLambdaSA
	ENDIF

	CALL TRANSE(TSG,QP)
	CALL CLEAR(TAU1,6)
      CALL MNOZI2(TAU1,TSG,TAU0,6,6)

	DO I=1,6
	 TAU(I) = TAU1(I)
	ENDDO

C	pert = 1.D-6
C	if(IPROL.EQ.0.AND.nHas.GT.0.D0) then
C	 stress0 = TAU(2)
C	 DEF(1) = DEF(1)+pert
C	 IPROL=1
C	 goto 999
C	endif
C	deriv = (TAU(2)-stress0)/pert


	
C	Total deformation - Potpuna deformacija
	DO I=1,3
	 DEF0(I) = DLOG(dLambdaTotal(I))
	ENDDO
	DO I=4,6
	 DEF0(I) = 0.0
	ENDDO

	CALL TRANSS(TSS,QP)
	CALL CLEAR(DEF1,6)
      CALL MNOZI2(DEF1,TSS,DEF0,6,6)

c	DO I=1,6
c	 DEF(I) = DEF1(I)
c	ENDDO
         
         call clear(DEF0,6)
         DO I=1,3
	   DEF0(I) = DLOG(dLambdaElDev(I))
	   ENDDO
	   DO I=4,6
	   DEF0(I) = 0.0
	   ENDDO
         DO 310 I=1,3
  310    DEF(I)=DEXP(2.D0*DEF0(I))
         CALL DIJADS(QP,DEF) 
      
C	Plastic deformation - Plasticna deformacija
	DO I=1,3
	 DEFP0(I) = dEpsilonL*dKsiNew*dNa(I)
	ENDDO
	DO I=4,6
	 DEFP0(I) = 0.0
	ENDDO

	CALL CLEAR(DEFP1,6)
      CALL MNOZI2(DEFP1,TSS,DEFP0,6,6)

C	Transform ELAST from principal to global
      IF(nHas.GT.0.OR.nHsa.GT.0) CALL TRAETP(ELAST,ELAST,TSG)

      dKsi1 = dKsiNew
	dFas1 = dFas
	dFsa1 = dFsa
	CALL JEDNA1(ELAST1, ELAST,36)

CC	WRITE(3,*) "TAU(1)=", TAU1(1)
CC	WRITE(3,*) "TAU(2)=", TAU1(2)
CC	WRITE(3,*) "TAU(3)=", TAU1(3)
CC	WRITE(3,*) "KsiT=", dKsiT
CC	WRITE(3,*) "El1=", ELAST(1,1)
CC	WRITE(3,*) "El12=", ELAST(1,2)
CC	WRITE(3,*) "El13=", ELAST(1,3)
CC	WRITE(3,*) "DEF1=", DEF(1)
CC	WRITE(3,*) "DEF2=", DEF(2)
CC	WRITE(3,*) "DEF3=", DEF(3)
CC	WRITE(3,*) "LambdaAS=", dLambdaAS
CC	WRITE(3,*) "LambdaSA=", dLambdaSA

C	CALL WRR3(QP,9,"QP")

      END
C=======================================================================
C=======================================================================

C=============================================================================================
C=============================================================================================


      SUBROUTINE ELAST_TANGENT_MATRIX(dE, dNi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C      DIMENSION FUN(2,MATE,*)
C
      E=dE
      V=dNi
CE  NULL ( ELAST )
      DO 15 I=1,6
      DO 15 J=I,6
   15 ELAST(I,J)=0.0D0
C
      ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
      ELAST(2,2)=ELAST(1,1)
      ELAST(3,3)=ELAST(1,1)
      ELAST(1,2)=ELAST(1,1)*V/(1.-V)
      ELAST(1,3)=ELAST(1,2)
      ELAST(2,3)=ELAST(1,2)
      ELAST(4,4)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
      ELAST(5,5)=ELAST(4,4)
      ELAST(6,6)=ELAST(4,4)
      DO 50 I=1,6
      DO 50 J=I,6
   50 ELAST(J,I)=ELAST(I,J)
      RETURN
      END

C=======================================================================
C	SHAPE-MEMORY ALLOY MODEL FUNCTIONS
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE TANGENTMATRIX(dT1as,dT2as,dT1sa,dT2sa,dK,dG,dEpsilonL,
	1                         dAlpha,dKsi,dDefNorm,dNa,dE,dNi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION dNa(3), dNaDyad(3,3), dJedinicna(3,3), dJedDev(3,3)
	DIMENSION dNaDyad1(3,3), d1DyadNa(3,3), dWWW(3,3)

	dKi=dK*(1-3.D0*dEpsilonL*dAlpha*(dT2as+dT2sa))
	dGi=dG*(1.D0-dEpsilonL*dKsi/dDefNorm)
	dM1i=2.D0*dG*dEpsilonL*(dKsi/dDefNorm-(dT1as+dT1sa))
	dM2i=2.D0*dG*dEpsilonL*(dT2as+dT2sa)
	
	DO I=1,3
	 DO J=1,3
	  dNaDyad(I,J)=dNa(I)*dNa(J)
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  IF(I.EQ.J) THEN
	   dJedinicna(I,J)=1.D0
	  ELSE
	   dJedinicna(I,J)=0.D0
	  ENDIF
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  dJedDev(I,J)=dJedinicna(I,J)-1.D0/3.D0
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  dNaDyad1(I,J)=dNa(I)
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  d1DyadNa(I,J)=dNa(J)
	 ENDDO
	ENDDO

C	Elasticna tangentna matrica
C	dEplasticno = dE*0.05
C	CALL ELAST_TANGENT_MATRIX(dEplasticno,dNi)
	CALL ELAST_TANGENT_MATRIX(dE,dNi)

	DO I=1,3
	 DO J=1,3
	  ELAST(I,J)=dKi+2.D0*dGi*dJedDev(I,J)
	1            +dM1i*dNaDyad(I,J)-dM2i*(dNaDyad1(I,J)+d1DyadNa(I,J))
	 ENDDO
	ENDDO

      RETURN
      END

C=======================================================================
C=======================================================================
C	Nasa
C=======================================================================
      SUBROUTINE TANGENTMATRIX_N(dT1as,dT2as,dT1sa,dT2sa,dK,dG,
	1				dEpsilonL,dAlpha,dKsi,dDefNorm,dNa,dE,dNi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION dNa(3), dNaDyad(3,3), dJedinicna(3,3), dJedDev(3,3)
	DIMENSION dNaDyad1(3,3), d1DyadNa(3,3), dWWW(3,3)

	dKi=dK*(1.D0-3.D0*dEpsilonL*dAlpha*(dT2as+dT2sa))
	dPi=3.D0*dK*dAlpha*dEpsilonL*(dT1as+dT1sa)
	dGi=2.D0*dG*(1.D0-dEpsilonL*dKsi/dDefNorm)
	dRi=2.D0*dG*dEpsilonL*dKsi/dDefNorm
	dM1i=2.D0*dG*dEpsilonL*(dT1as+dT1sa)
	dM2i=2.D0*dG*dEpsilonL*(dT2as+dT2sa)
	
	DO I=1,3
	 DO J=1,3
	  dNaDyad(I,J)=dNa(I)*dNa(J)
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  IF(I.EQ.J) THEN
	   dJedinicna(I,J)=1.D0
	  ELSE
	   dJedinicna(I,J)=0.D0
	  ENDIF
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  dJedDev(I,J)=dJedinicna(I,J)-1.D0/3.D0
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  dNaDyad1(I,J)=dNa(I)
	 ENDDO
	ENDDO

	DO I=1,3
	 DO J=1,3
	  d1DyadNa(I,J)=dNa(J)
	 ENDDO
	ENDDO

	CALL MNOZM1(dWWW,dNaDyad,dJedDev,3,3,3)
		
C	Elasticna tangentna matrica
	dEplasticno = dE*0.05
C	CALL ELAST_TANGENT_MATRIX(dEplasticno,dNi)
C	CALL ELAST_TANGENT_MATRIX(dE,dNi)
	CALL CLEAR(ELAST,36)

	DO I=1,3
	 DO J=1,3
	  ELAST(I,J)=dKi-dPi*dNaDyad1(I,J)+dGi*dJedDev(I,J)
	1            +dRi*dWWW(I,J)
	1            -dM1i*dNaDyad(I,J)-dM2i*d1DyadNa(I,J)
	 ENDDO
	ENDDO

C	Modul smicanja izracunat preko ekvivalentnog modula elasticnosti dobijenog iz plasticne kons. matrice
	DO I=4,6
	  ELAST(I,I)=ELAST(I-3,I-3)*(1.D0-2.D0*dNi)/(2.D0*(1.D0-dNi))
	ENDDO


      RETURN
      END

C=======================================================================
C=======================================================================
C	Calculate Drucker-Prager-type loading function
C	transformation austenit to martensit A->S (Auricchio, Taylor, 1996)
C=======================================================================
      SUBROUTINE CALCULATE_FAS(dTau,dP,dMatProperties,dFas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
      DIMENSION dTau(*)
C
	dAlpha = dMatProperties(17)
	dCas = dMatProperties(9)
	dT = dMatProperties(4)

      dTauNorm2=dTau(1)*dTau(1)+dTau(2)*dTau(2)+dTau(3)*dTau(3)
	dTauNorm=DSQRT(dTauNorm2)

	dFas = dTauNorm + 3 * dAlpha * dP - dCas * dT

      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE CALCULATE_RSAS(dMatProperties,dRsas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
C
	dSsas = dMatProperties(11)
	dCas = dMatProperties(9)
	dTsas = dMatProperties(5)
	dAlpha = dMatProperties(17)

	dRsas = dSsas * ( DSQRT(2.D0/3.D0) + dAlpha) - dCas * dTsas

      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE CALCULATE_RFAS(dMatProperties,dRfas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
C
	dSfas = dMatProperties(12)
	dCas = dMatProperties(9)
	dTfas = dMatProperties(6)
	dAlpha = dMatProperties(17)

	dRfas = dSfas * ( DSQRT(2.D0/3.D0) + dAlpha) - dCas * dTfas

      RETURN
      END
C======================================================================
C=======================================================================
C	Calculate Drucker-Prager-type loading function
C	transformation martensit to austenit S->A (Auricchio, Taylor, 1996)
C=======================================================================
      SUBROUTINE CALCULATE_FSA(dTau,dP,dMatProperties,dFsa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
      DIMENSION dTau(*)
C
	dAlpha = dMatProperties(17)
	dCsa = dMatProperties(10)
	dT = dMatProperties(4)

      dTauNorm2=dTau(1)*dTau(1)+dTau(2)*dTau(2)+dTau(3)*dTau(3)
	dTauNorm=DSQRT(dTauNorm2)

	dFsa = dTauNorm + 3 * dAlpha * dP - dCsa * dT

      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE CALCULATE_RSSA(dMatProperties,dRssa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
C
	dSssa = dMatProperties(13)
	dCsa = dMatProperties(10)
	dTssa = dMatProperties(7)
	dAlpha = dMatProperties(17)

	dRssa = dSssa * ( DSQRT(2.D0/3.D0) + dAlpha) - dCsa * dTssa

      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE CALCULATE_RFSA(dMatProperties,dRfsa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION dMatProperties(*)
C
	dSfsa = dMatProperties(14)
	dCsa = dMatProperties(10)
	dTfsa = dMatProperties(8)
	dAlpha = dMatProperties(17)

	dRfsa = dSfsa * ( DSQRT(2.D0/3.D0) + dAlpha) - dCsa * dTfsa

      RETURN
      END
C======================================================================