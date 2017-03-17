C=======================================================================
C=======================================================================
CS    GENERALISANI MODEL SA KAPOM  3/D ELEMENT
CE    GENERALIZED CAP MODEL 3D ELEMENT
C=======================================================================
C=======================================================================
CE    SUBROUTINE D3M41
CE               TI3441
C
      SUBROUTINE D3M41(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
C
      include 'paka.inc'
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M41'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LEMP=LDEFPP + 6
      LXT=LEMP + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LEMP1=LDEFP1 + 6
      LXTDT=LEMP1 + 1
C
      CALL TI3441(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3441(TAUT,DEFT,DEFPP,EMP,XT,
     1                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     1                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    DRUCKER-PRAGER MODEL
CE    PROGRAM FOR STRESS INTEGRATION FOR DRUCKER-PRAGER MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1               DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8 
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ITERBR/ ITER
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(6),DEFT(6),DEFPP(6),TAU(6),DEF(6),TAU1(6),DEF1(6), 
     +          DEFP1(6)
      DIMENSION FUN(11,*),NTA(*),DSIG(6),DEPS(6),DFDS(6),DGDS(6),ALAM(6)
     +         ,CP(6,6),POM(6), CEP(6,6),DDEFE(6),TAUDE(6)
     +         ,DFdDS(6),DFcDS(6),DGdDS(6),DGcDS(6),ALAMd(6),ALAMc(6)
     +         ,al1p(6),al2p(6),ddefpd(6),ddefpc(6),taue(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3441'
C
      IF(IRAC.EQ.2) RETURN
CE    BASIC KONSTANTS
      TOL =   1.D-8
	  TOLL=  -1.D-5
      MAXT=   500.D0
C====================================================================
C
CE.   MATERIAL CONSTANTS 
      E       = FUN(1,MAT) 
      ANI     = FUN(2,MAT) 
C
      AK      = FUN(3,MAT)
      ALF     = FUN(6,MAT)
      T       = FUN(7,MAT)
      X0      =-FUN(8,MAT)
C
      W       =-FUN(9,MAT)
      D       =-FUN(10,MAT)
C====================================================================
C
      CM      = E/(1-2.D0*ANI)
      G       = 0.5D0*E/(1.D0+ ANI)
      DUM     = 0.01D0
      XTMAX   =-DLOG(DUM)/D+X0
C
      IF(KOR.EQ.1) THEN
	   XT=X0
	ENDIF
C====================================================================
C
CE    ELASTICITY MATRIX FOR 3D
      CALL MEL3EL(ELAST,E,ANI)
      IF(IRAC.EQ.2) RETURN
C
CE    MEAN STRAIN emE
      EMT=(DEF(1)+DEF(2)+DEF(3))/3.0D0
C
CE    TOTAL DEVIATORIC STRAINS e'=e-em
      DO I=1,6
         IF(I.LE.3) THEN
            DEFDPR(I)=DEF(I)-EMT
         ELSE
            DEFDPR(I)=0.5D0*DEF(I)
         ENDIF 
      ENDDO
C
CE    TRIAL ELASTIC DEVIATORIC STRAINS e''=e'-emp'
      DO I=1,6
         IF(I.LE.3) THEN
            DEFDS(I)=DEFDPR(I)-DEFPP(I)+EMP
         ELSE
            DEFDS(I)=DEFDPR(I)-DEFPP(I)
         ENDIF
      ENDDO
C
CE    TRIAL ELASTIC MEAN STRAIN em''=em-emp
      EMS=EMT-EMP
C
CE    ELASTIC SOLUTION SE=2G*e''
      DO I=1,6
          TAUD(I)=2.0D0*G*DEFDS(I)
      ENDDO
CD
      DO I=1,6
          TAUDE(I)=TAUD(I)
      ENDDO
CE    TRIAL ELASTIC MEAN STRESS sigmaE=cm*em''
      SMTE=EMS*CM
C=================================================================
C
CE    SECOND INVARIANT OF DEVIATORIC STRESSES J2d
      AJ2DE=0.5*(TAUDE(1)*TAUDE(1)+TAUDE(2)*TAUDE(2)+TAUDE(3)*TAUDE(3))+
     1           TAUDE(4)*TAUDE(4)+TAUDE(5)*TAUDE(5)+TAUDE(6)*TAUDE(6)
c	CALL WRR6 (TAUDE,6,'TAUDE')
      AJ2DEQ  =DSQRT(AJ2DE)
C
CE    FIRST STRESS INVARIANT I1E=3*sigmamE
      AI1E    =3.0D0*SMTE
C
C ================================================================     
CE    TENSION CUTOFF (FT>0)
      IF(ALF.GT.TOL.AND.T.GT.(AK/ALF))THEN
          T=0.99*(AK/ALF)
     	ENDIF
C
     	IF(AI1E.GT.T) THEN
          DO I=1,6
     	        IF(I.LE.3)THEN
     	            TAU(I)=T/3.0D0
     	        ELSE
     	            TAU(I)=0
     	        ENDIF
          ENDDO
          GO TO 500
      ENDIF 
C
C     DRUCKER-PRAGER LINE
      FDP=ALF*AI1E+AJ2DEQ-AK
C      
C     CAP LINE
      FC=XT-AI1E
C
     	SMTDT=SMTE
      DELEMP=0.
CE    ELASTIC SOLUTION
     	IF(FDP.LE.TOL) GO TO 400
C==================================================================
C
CE    PLASTIC STRAIN
C
CE    DRUCKER-PRAGER YIELDING FDP>0 & FC<0
      IF(FDP.GT.TOL)GOTO 100
C
CE    CUP YIELDING FC>0 & FDP<0
c      IF(FDP.LT.TOL.AND.FC.GT.TOL)GOTO 200
C
CE    VERTEX YIELDING FC>0 & FDP>0
c     	IF(FDP.GE.TOL.AND.FC.GE.TOL)GOTO 300
C      
C      
      STOP ' NONDEFINED AREA - DRUCKER-PRAGER '
C===================================================================
C===================================================================
C
CE    DRUCKER-PRAGER YIELDING FDP>0 & FC<0
  100 CONTINUE
      IYIELD=1
c      WRITE(3,*)'OBLAST 1'
      IF(DABS(ALF).LT.TOL) THEN
CE       INCREMENT OF MEAN PLASTIC STRAIN Demp
         AJ2DQ=AK
C        sqrt(J2D)=sqrt(J2DE)-Dlam*G
CE       Dlambda
         DLAM=(AJ2DEQ-AJ2DQ)/G
	   		 DELEMP=0.D0
      	 ELSE
CE       INCREMENT OF MEAN PLASTIC STRAIN Demp
         DELEMP=(ALF*AI1E+AJ2DEQ-AK)/(3.0D0*ALF*CM+G/ALF)
CE       Dlambda.
         DLAM=DELEMP/ALF
C        sqrt(J2D)=sqrt(J2DE)-Dlam*G (6.2.29)
         AJ2DQ=AJ2DEQ-DLAM*G
      ENDIF
C
CE    MEAN STRESS sigmam=sigmamE-cm*Demp
			SMTDT=SMTE-CM*DELEMP
C
CE    FIRST STRESS INVARIANT
      AI1=3.0D0*SMTDT
C
CE    MEAN PLASTIC STRAIN emp=emp+Demp
      EMP1=EMP+DELEMP
C
C====================================================================
CE    DEVIATORIC STRESS S
      DEN=1.0D0+DLAM*G/AJ2DQ
      DO I=1,6
         TAUD(I)=TAUDE(I)/DEN
      ENDDO
C
CE    ADD INCREMENT OF DEVIATORIC PLASTIC STRAIN ep'=ep'+De'p 
      COEF=DLAM/AJ2DQ/2.0D0
      DO I=1,6
         DEFP1(I)=DEFPP(I)+COEF*TAUD(I)
      ENDDO
C
			AJ2D=AJ2DQ*AJ2DQ
      GO TO 400
C====================================================================
C====================================================================
C
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
c      IF(XTDT.GT.T) STOP '(XTDT.GT.T) KRAJ'
CE    ADD INCREMENT OF MEAN PLASTIC STRAIN ep=ep+Demp
      DO I=1,6
         IF(I.LE.3) THEN
            DEFP1(I)=DEFP1(I)+DELEMP
         ELSE
            DEFP1(I)=DEFP1(I)
         ENDIF
      ENDDO
C
CE    STRESS CALCULATION sigma=S+sigmam
      DO I=1,6
         IF(I.LE.3) THEN
            TAU(I)=TAUD(I)+SMTDT
         ELSE
            TAU(I)=TAUD(I)
         ENDIF
      ENDDO
cd      CALL CEPDP3(TAUD,TAUDE,ALF,G,CM,AJ2DEQ,AJ2D,
cd    +            AJ2DQ,DLAM,DEFDPR,IYIELD,X0,XTDT,D,W)
c      CALL CEPDP3(TAU,6,'TAU ')
c      CALL WRR6(DEFP1,6,'DEFP')
C========================================================================
C
      IF(ITER.EQ.0) GO TO 500
C
  500 CONTINUE
      DO 290 I=1,6
      DEF1(I)=DEF(I) 
  290 TAU1(I)=TAU(I)
      RETURN
      END