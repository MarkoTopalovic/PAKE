C=======================================================================
C
C   VISCOPLASTICNOST ZA OPSTI MODEL SA KAPOM  - 3/D ELEMENT
C        GENERALIZED CAP MODEL
C
C    SUBROUTINE D3M23
C               TAUI23
C
      SUBROUTINE D3M23(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
      include 'paka.inc'
      
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M23'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
CC      MATE=1
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LSIGS=LDEFPP + 6
      LEMP=LSIGS + 6
      LXT=LEMP + 1
      LELT=LXT + 1
      LEMP0=LELT + 1
      LIPL=LEMP0 + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LSIGS1=LDEFP1 + 6
      LEMP1=LSIGS1 + 6
      LXTDT=LEMP1 + 1
      LELTDT=LXTDT + 1
      LEMP01=LELTDT + 1
      LIPL1=LEMP01 + 1
C
      CALL TI3423(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),PLAST(LSIGS),
     1            PLAST(LEMP),PLAST(LXT),PLAST(LELT),
     1            PLAST(LEMP0),PLAST(LIPL),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LSIGS1),PLAS1(LEMP1),PLAS1(LXTDT),
     1            PLAS1(LELTDT),PLAS1(LEMP01),PLAS1(LIPL1), 
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C  =====================================================================
C
C
C  =====================================================================
C
      SUBROUTINE TI3423(TAUT,DEFT,DEFPP,SIGS,EMP,XT,ELT,EMP0,PL,
     1           TAU1,DEF1,DEFP1,SIGS1,EMP1,XTDT,ELTDT,EMP01,PL1,
     1           FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
C     VISKOPLASTICNI GENERALISANI MATERIJALNI MODEL SA KAPOM
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1              DETAU(6),DDEF(6)
      COMMON/MAT2D/E,ANI,ET,TEQY0
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON/PLASTI/LPLAST,LPLAS1,LSIGMA
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
C
      COMMON/ITERBR/ITER
C
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /CDEBUG/ IDEBUG
C      COMMON /ICT/ IBTC
C
      DIMENSION TAUT(*),DEFT(*),DEFPP(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP1(*),SIGS(*),SIGS1(*),SE(6),SS(6),DEPVP(6)
      DIMENSION FUN(34,*),NTA(*),FUN2(34)
      DIMENSION CMP(6),WJ(6),CP(6,6),A33(3,3)
      IF(IDEBUG.EQ.1) PRINT *, 'TI3423'
C
C     OSNOVNE KONSTANTE
C
      TOLFY=1.0D-8
C
C
      DVT=2.0D0/3.0D0
      SQ2=DSQRT(2.0D0)
      DJP=DSQRT(1.5D0)
      MAXIT = 999
      TOL = 1.D-9
      TOLL = 1.D-8
      TOLYLD = 1.D-8
      SQRT2 = SQRT(2.D0)
      AF =1.0D0
      AF3 = AF*AF*AF
C
C     MATERIAL CONSTANTS
      A     =  FUN( 1,MAT)
      B1    = -FUN( 2,MAT)
      C     =  FUN( 3,MAT)
      TETA  = -FUN( 4,MAT)
      AT  =   -FUN( 5,MAT)
      AI1A0  = FUN( 6,MAT)
      IEL    = INT(FUN( 7,MAT))
C
      W     = -FUN( 8,MAT)
      D     = -FUN( 9,MAT)
      ALFAM =  FUN( 10,MAT)
      W1    = -FUN( 11,MAT)
      D1    = -FUN(12,MAT)
      D2    =  FUN(13,MAT)
      D3    = -FUN(14,MAT)
C
      R0    =  FUN(15,MAT)
      R1    =  FUN(16,MAT)
      R2    = -FUN(17,MAT)
      R3    =  FUN(18,MAT)
      R4    =  FUN(19,MAT)
      R5    =  FUN(20,MAT)
C
      AKEI  =  FUN(21,MAT)
      AKS   =  FUN(22,MAT)
C     PROMENJENE OZNAKE ZA AK1 I AK2 U AK1PL I AK2PL
      AK1PL = -FUN(23,MAT)
      AK2PL = -FUN(24,MAT)
C     PROMENJENA OZNAKA ZA BETA U BETAPL
      BETAPL  =  FUN(25,MAT)
      DELTA =  FUN(26,MAT)
C
      GEI   =  FUN(27,MAT)
      GS    =  FUN(28,MAT)
      G1    =  FUN(29,MAT)
      G2    =  FUN(30,MAT)
      G3    = -FUN(31,MAT)
C     PROMENJENA OZNAKA ZA ETA U ETAPL
      ETAPL   =  FUN(32,MAT)
      GAMA  =  FUN(33,MAT)
C**   ZA VISKOPLASTICNOST
      ETA = FUN(34,MAT)
C
      AT=0.D0
C      IF (DABS(TETA).GT.TOL) THEN
C      AT=(C-A)/TETA
C      ENDIF
      B1C = B1*C
      R0R1R2 = R0*R1*R2/(1.+R1)
      R3R4 = 2.*R3*R4
C
      E     = FUN(21,MAT)
      ANI   = FUN(27,MAT)
      R     = FUN(15,MAT)   
      AE = (1.+ANI)/E	
CC     VISCOPLASTICITY CONSTANTS
      ETADT = ETA/DT
      DTETA = DT/ETA
      BETA =1./(1.+DTETA/AE)
      BETA1 = 1.+BETA*DTETA/AE
      CETA = AE*ETADT
C      TETA =0.D0
      B1C = B1*C
      R0R1R2 = R0*R1*R2/(1.+R1)
      R3R4 = 2.*R3*R4
C
      IF(IEL.EQ.0) THEN
       E     = FUN(21,MAT)
       ANI   = FUN(27,MAT)
       R     = FUN(15,MAT)	
      ENDIF 
C
C     INICIJALIZACIJA OSNOVNIH VELICINA
C     INITIAL VALUES OF VARIABLES
C
      IF(KOR.EQ.1) THEN
       R = R0 + R3*DEXP(-R4*R5*R5)
       R=1.0D0
      ELSE
       EL=ELT        
       ELR5 =  EL + R5
       ELR5 = ELR5*ELR5
C       R = R0/(1.+R1)*(1.+R1*DEXP(-R2*EL)) + R3*DEXP(-R4*ELR5)
      ENDIF
      R=1.0D0
      IF(KOR.EQ.1.OR.ITER.GT.1) GO TO 482
      IF(KOR.NE.1) GO TO 483
      EPVT = 0.D0
      EPPV = 0.D0
      AI1 = 0.D0
      AJ2D = 0.D0
      EL = 0.D0
      ELT = 0.D0
C
       AI1A = AI1A0
       DELEMP=0.0D0
C
      XP = 0.D0
      AI1AS = AI1A
      ELS = ELT
CC      EL =ELT
      EL =0.0D0
      XT=AI1A0
      I = 0 
C
C     BISECTION LOOP
C
  481 I = I + 1
C     SOLUTION FOR EL (NEWTON ITERATIONS)
      RP=0.0D0
C
      F = EL-AI1A-R*(A-C*DEXP(-B1*EL)+TETA*EL) 
C
      FP = 1.- R*(C*B1*DEXP(-B1*EL)+TETA) 
C
      DEL = F/FP
      EL = EL - DEL
C
      IF  (I.GT.MAXIT) THEN                                                  
          WRITE (6,1009) I,X,F,FP
          STOP 'DP RESENJE X'                                                                  
1009  FORMAT (' CAPVP-INITIAL VALUE FOR X'/
     1 ' NO SOLUTION IN BISECTING FOR F=0. '/
     1 ' IBIS =',I5/
     3 ' X=',D14.6/
     4 ' F =',D14.6/ ' F-PLUS=',D14.6)
      ENDIF                                                                     
C                                                                               
      IF (DABS(DEL).GT.TOLL) GO TO 481
C
      AI1A0=AI1A
      XT = AI1A0
      XTDT = AI1A0
      ELT = EL
      ELTDT = EL
C
 482  AI1A = AI1A0
      XT = AI1A0
      XTDT = XT
      ELTDT = ELT      
C
      XSQ = AI1A*AI1A
      DAI1A = D*AI1A
CC      DEPV =-W*(1.-DEXP(-D*AI1A))-EPVT
      DEPV =W*(1.-DEXP(-D*AI1A))
C      DEPV = -W*(1.-DEXP(-DAI1A)-ALFAM*DAI1A*DEXP(-D1*AI1A-D2*XSQ)) -
C     1        W1*XSQ*DEXP(-D3*AI1A) - EPVT
      EMP0=DEPV/3.
C
 483  CM = E/(1-2.*ANI)
      G = 0.5*E/(1.+ ANI)
       AE2 = 1./G
       AE = 0.5*AE2
       DDL=1.0D-7
       AM =1.0/CM
       CM6 = 6.*CM*DTETA
C
      IF(KOR.GT.1) SM0=0.0D0
      FUN2(1)=E
      FUN2(2)=ANI
C
      CALL MEL31(FUN2)
      IF(IRAC.EQ.2) RETURN
C
      AE=(1.0D0+ANI)/E
      DO 16 I=1,6
   16 DDEF(I)=DEF(I)-DEFPP(I)
C
C     STRESS ESTIMATION
C
      IF (ITER.EQ.0) THEN
         IF (KOR.EQ.1) THEN
            DO 166 I=1,6
  166       TAU(I)=TAUT(I)
            RETURN
         ENDIF 
        DO 167 I = 1,6  
  167   DEPVP(I)=(TAU(I)-SIGS(I))/ETADT
         DO 169 I=1,6
         TAU(I) = 0.D0
         DO 169 J = 1,6
         TAU(I)=ELAST(I,J)*(DEF(J)-DEFPP(J)-DEPVP(J))
  169    CONTINUE
         RETURN
      ENDIF
C 
C     ODREDIVANJE DEVIJATORA UKUPNE DEFORMACIJE
C
      EMT = (DEF(1)+DEF(2)+DEF(3))/3.0
      DO 30 I=1,6
      IF(I.LE.3) THEN
       DEFDPR(I)=DEF(I)-EMT
        DEFDS(I)=DEFDPR(I)-DEFPP(I)+EMP
      ELSE
       DEFDPR(I)=0.5D0*DEF(I)
        DEFDS(I)=DEFDPR(I)-0.5*DEFPP(I)
      ENDIF
   30 CONTINUE
C
C     PROVERA ELASTICNOG RESENJA
C     ODREDIVANJE NAPONA KOJI ODGOVARA ELASTICNOM RESENJU
C
      DO 32 I=1,6
      TAUD(I)=DEFDS(I)/AE
   32 SE(I) = TAUD(I)
C
      AJ2DE=0.5*(TAUD(1)*TAUD(1)+TAUD(2)*TAUD(2)+TAUD(3)*TAUD(3))+
     1   TAUD(4)*TAUD(4)+TAUD(5)*TAUD(5)+TAUD(6)*TAUD(6)
      AJ2DQ=DSQRT(AJ2DE)
      AJ2DQE = AJ2DQ
C
      EMS = EMT - EMP
      SMTE =EMS*CM
      SMTDT=SMTE
      IF(3.0*SMTE.GT.AT.OR.PL1.LT.0.0D0) THEN
        DO 117 I=1,6
  117   TAU(I)=0.0D0
        PL1=-1.0D0
C       PL1=1.0D0
        GO TO 500
      ENDIF 
C
      AI1A=XT
      AI1 = 3.0*SMTE
      AI1L = AI1-ELT
      XML = AI1A - ELT
      RSQ = R*R
C
      FC = AI1L*AI1L + RSQ*AJ2DE - XML*XML
      FDP = AJ2DQ - (A - C*DEXP(-B1*AI1))-TETA*AI1
C
      IF((FC.LE.TOLYLD .AND. FDP.LE.TOLYLD).OR.
     * (FDP.LE.TOLYLD.AND.AI1.GT.ELT)) THEN
          IELAST = 0
          GO TO 500
      ENDIF
C
C     PLASTIC DEFORMATION
C
      IELAST = 1
      IF (FDP.GT.TOLYLD .AND. AI1L.GT.TOLYLD) GO TO 100
      IF (FC.GT.TOLYLD .AND. AI1L.GT.TOLYLD) THEN
          IELAST = 0
          GO TO 400
      ENDIF
      IF (FC.GT.TOLYLD .AND.AI1L.LT.TOLYLD) GO TO 200
C
C     YIELDING ON THE FAILURE SURFACE
C
  100 IYIELD=1
      CETA = AE*ETADT
      CETA1 = (1.+CETA)/CETA
      CETA1H = 0.5*CETA1
      CETA12 = 2.*CETA1
      D2N = CM + ETADT
      D23N = 3.*D2N
      B13 = 3.*B1   
      PL1=1.0D0
      X = 0.D0
      I = 0 
C
C     NEWTON ITERATION LOOP
C
   10 I = I + 1
C
      SMTDT = CM*(EMS-X)
      SIGMS = SMTDT - ETADT*X
      AI1 = 3.*SIGMS
      D1N = C*DEXP(-B1*AI1)
      D1BN = B1*D1N
      AK = (-SIGMS+SMTDT)/(B1*D1N+TETA)
      C14 = B13*AK*D1BN
      DEN = D1BN + TETA
      C9 = (ETADT+C14*D2N)/DEN
      AJ2DQ = AJ2DQE - CETA1H*AK
      F = AJ2DQ - A + D1N - TETA*AI1
      F1P = C9/CETA12 + D23N*DEN
      DX = -F/F1P
      X=X+DX
      IF  (I.GT.MAXIT) THEN                                                  
C                                                                               
C         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
C                                                                               
          WRITE(3,2416) NLM,IBTC
          WRITE (3,1010) I,X,DX,F,F1P
          STOP                                                                  
      ENDIF                                                                     
 2416 FORMAT(' NLM,IBTC',2I5)
1010  FORMAT (' CAPVP-YIELDING ON FAILURE SURFACE'/
     1 ' NO SOLUTION IN NEWTON ITERATIONS FOR F=0. '/
     1 ' I =',I5/
     3 '  DEPVP=',D14.6/ ' DX =',D14.6/
     4 ' F =',D14.6/' F1P =',D14.6)
C                                                                               
C                                                                               
      IF (DABS(DX).GT.TOLL) GO TO 10
C
      XTDT = XT
      DELEMP = X
      EMP1 = EMP + DELEMP
C      IF(IBTC.EQ.1) 
C     * WRITE(3,7234)AI1A,EL
C 7234 FORMAT(' AI1A',2D16.6)
      AI1A = -DLOG(1.-3.*EMP1/W)/D + AI1A0
C      IF(IBTC.EQ.1) 
C     * WRITE(3,7234)AI1A
C
      AI1AS = AI1A
      ELS = ELT
CC      EL =ELT
      EL =0.0D0
      IF(KOR.EQ.1) THEN
       EPVT=3.*EMP0
CC       EMS=EMT-EMP0
      ENDIF 
C
C    NEWTON ITERATIONS LOOP FOR EL
C
  281 I = I + 1
C
C     SOLUTION FOR EL (NEWTON ITERATIONS)
C      IF(IBTC.EQ.1) 
C     * WRITE(3,1998) R,RP,EL
C 1998 FORMAT(' **R,RP,EL',3D12.4)
      RP=0.0D0
C
      F = EL-AI1A-R*(A-C*DEXP(-B1*EL)+TETA*EL) 
C
      FP = 1.- R*(C*B1*DEXP(-B1*EL)+TETA) 
C
      DEL = F/FP
      EL = EL - DEL
C
      IF  (I.GT.MAXIT) THEN                                                  
      WRITE(3,2416) NLM,IBTC
      WRITE (3,1012) I,EL,DEL,F,FP
1012  FORMAT (' CAPVP-CALCULATION OF EL AFTER YIELDING ON THE ',
     1 ' FAILURE SURFACE'/
     1 ' I =',I5/
     3 '  EL=',D14.6/ ' DX =',D14.6/
     4 ' F =',D14.6/' FP =',D14.6)
      STOP
      ENDIF                                                                     
C                                                                               
      IF (DABS(DEL).GT.TOLL) GO TO 281
C                                                                               
      ELTDT = EL
      XTDT = AI1A
      AKAJ = AK/(2.*AJ2DQ)
      COEF = CETA/(CETA+(1.+CETA)*AKAJ)
      COEFS = 1. + AKAJ
      DO 615 I=1,6
      SS(I)=COEF*SE(I)
      TAUD(I) = COEFS*SS(I)
      DEPVP(I) = AE*(SE(I)-TAUD(I))
      IF  (I.LE.3) THEN
          DEPVP(I) = DEPVP(I) + EMP1
      ELSE
          DEPVP(I) = DEPVP(I) + DEPVP(I)
      ENDIF
  615 CONTINUE
      GO TO 400
C
C     CAP YIELDING
C
C     CAP YIELDING - NEWTON ITERATIONS
C 
200   ETADT2 = 0.5*ETADT
       IF (KOR.GT.27) THEN
       PERILO=3.
       ENDIF
C      CETA = AE*DTETA
      IYIELD=2
      ALFAD = ALFAM*D
      D22 = 2.*D2
      AI1AS = AI1A
      ELS = ELT
      EL =ELT
C      DX = 0.02*AI1A
      X = 0.D0
      I = 0 
      EPVT=3.*EMP
      IF(KOR.EQ.1) THEN
       EPVT=3.*EMP0
CC       EMS=EMT-EMP0
      ENDIF 
C
C     NEWTON ITERATIONS LOOP
      X = -TOL
      IF (KOR.GT.27) X=-1.0D-3
C
  210 I = I + 1
C
C      AI1A = AI1AS + X
C      IF (KOR.EQ.8.AND.I.GT.10) THEN
C      IF (ITER.GT.0) AF=1.0D0
C      ENDIF
C
      EL = ELS + X
      ELR5 =  EL + R5
      ELR5 = ELR5*ELR5
      RP=0.0D0
      EXP1 = C*DEXP(-B1*EL)
      ACTETA = A-EXP1+TETA*EL
      IF (ACTETA.LE.0) THEN 
      WRITE(3,*) 'GRESKA '
      ENDIF
      AI1A = EL - R*ACTETA
      AI1AP = 1. - RP*ACTETA - R*(B1*EXP1+TETA)
      XSQ = AI1A*AI1A
      DAI1A = ALFAD*AI1A
C  KOJIC    D1D2 = D1-D22*AI1A
      D1D2 = D1+D22*AI1A
      EXP1 = DEXP(-D*AI1A)
      EXP2 = DEXP(-D1*AI1A-D2*XSQ)
      WEXP3 = W1*DEXP(-D3*AI1A)
      DEPV = W*(1.-EXP1-DAI1A*EXP2) + XSQ*WEXP3
      DELEMP = DEPV/3.-EMP
C
      DEPVPR = W*(D*EXP1-(ALFAD-DAI1A*D1D2)*EXP2) + (2.*AI1A-D3*XSQ)*
     1        WEXP3
C
      B4 = DEPVPR*AI1AP
      B5 = B4*(CM+ETADT)
      B4 = B4/3.
      B4C = B4
      SMTDT = CM*(EMS-DELEMP)
      AI1 = 3.*SMTDT
C      XL = DELEMP/(2.*(AI1-EL))
      SIGMS = SMTDT - ETADT*DELEMP
      AI1S = 3.*SIGMS
      AI1L = AI1S-EL
      AK = ETADT2/AI1L*DELEMP
      RSQ = R*R
C      BKR = BETA/(BETA1+AK*RSQ)
      BKR = CETA/(2.D0+CETA+(1.D0+CETA)*AK*RSQ)
      AJ2D = BKR*BKR*AJ2DE
      AJ2DQ = DSQRT(AJ2D)
      XML = AI1A - EL
      B1N = BETA1 + AK*RSQ
      AK2 = 2.*AK
CC      AK1 = (B4*(1.+CM6*AK)+AK2*DTETA)*AK/DELEMP
C      AK1 = 0.5*AK/(DELEMP*AI1L)*(B4-(2*AI1S-1.)*DELEMP)
      AK1=-AK/AI1L*(-B5-1.)+AK/DELEMP*B4
CC      AK12R = (R*AK1+AK2*RP)*R
CC      B2N = BKR*AK12R/B1N
C      BKRP=-BKR*AK1/(BETA1+AK*RSQ)
       BKRP=-BKR/(2.D0+CETA+(1.D0+CETA)*AK*RSQ)*(1.+CETA)*RSQ*AK1
C
C      CALL BISECR (X,XM,XP,DX,F,FM,FP,AF,IB)
CCC      FP = -AI1L*(B5+1.) - R*(BKR*B2N*R*AJ2DE-RP*AJ2D) - XML*(AI1AP-1.)
       F = AI1L*AI1L + RSQ*AJ2D - XML*XML
       FP = -AI1L*(B5+1.)+RSQ*AJ2DE*BKR*BKRP-XML*(AI1AP-1.)
       DX = 0.5* F/FP
       X = X - DX
C       IF (X.GT.0.D0) X=0.D0
C
C       WRITE (3,1211) F,X,I
C 1211  FORMAT('F=',D14.6,'   X=',D14.6,'  I=',I3)
      IF  (I.GT.MAXIT) THEN                                                  
C                                                                               
C         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
C                                                                               
          WRITE(3,2416) NLM,IBTC
          WRITE (3,1011) I,X,DX,F,FP
1011  FORMAT (' CAPVP-YIELDING ON THE CAP'/
     1 ' NO SOLUTION IN NEWTON ITERATIONS FOR F=0. '/
     1 ' I =',I5/
     3 '  X =',D14.6/ ' DX =',D14.6/
     4 ' F =',D14.6/' FP =',D14.6)
          STOP                                                                  
      ENDIF                                                                     
C
C      IF (DABS(DEL).GT.TOLL) GO TO 210
      IF (DABS(DX).GT.TOLL) GO TO 210
C
      ELTDT = EL
      XTDT = AI1A
      EMP1 = EMP + DELEMP
      COEF = 1. + AK*RSQ
      DO 625 I=1,6
      SS(I) = BKR*SE(I)
      TAUD(I)=COEF*SS(I)
      DEPVP(I) = DTETA*(TAUD(I)-SS(I))
      IF  (I.LE.3) THEN
          DEPVP(I) = DEPVP(I) + DELEMP
      ELSE
          DEPVP(I) = DEPVP(I) + DEPVP(I)
      ENDIF
  625 CONTINUE
C
C     ODREDIVANJE VISKOPLASTICNIH DEFORMACIJA I NAPONA
C
  400 DO 710 I=1,6
      DEFP1(I)=DEFPP(I)+DEPVP(I)
  710 CONTINUE
      DO 260 I=1,6
       IF(I.LE.3) THEN
        TAU(I)=TAUD(I)+SMTDT
         ELSE
        TAU(I)=TAUD(I)
       ENDIF
  260 CONTINUE
C
      IF(ITER.EQ.0) GO TO 500
      IF((METOD.NE.4).AND.(METOD.NE.3).AND.(METOD.NE.63).AND.
     *(METOD.NE.73)) GO TO 500
C
      IF (IYIELD.EQ.2) GO TO 350
C
C     YIELDING ON THE FAILURE SURFACE
C
      DEN = CETA+(1.+CETA)*AK/(2.*AJ2DQ)
      C1=CETA/DEN
      CETA1=1.+CETA
      DEN2=DEN*DEN
      C2=0.25*CETA*G/(DEN2*AJ2DQE*AJ2DQ*AJ2DQ*AJ2DQ)
      CETA2=0.5/CETA
      C3=CETA*CETA1/DEN2*(0.5/AJ2DQ+CETA2)
      C4=C9
      D1B3=3.*(D1BN+TETA)
      DEN=C4*CETA1*CETA2+D1B3*D2N
      COEF=G/(AJ2DQE*DEN)
      DO 310 I=1,6
  310 WJ(I)=-COEF*SE(I)
      COEFM=D1B3*CM/(3.*DEN)
      WJ(1)=WJ(1)+COEFM
      WJ(2)=WJ(2)+COEFM
      WJ(3)=WJ(3)+COEFM
      DO 311 I = 1,6
      DO 311 J = 1,6
      CP(I,J)=0.D0
  311 ELAST(I,J)=0.D0
      C33=DVT/3.
      DO 312 I=1,3
      DO 312 J=1,3
      A33(I,J)=-C33
      IF (I.EQ.J) A33(I,J)=DVT
  312 CONTINUE
      COEF=(CETA+C1)/CETA1*2.*G
      COEF1=0.5*COEF
      DO 314 I=1,3
      DO 313 J=1,3
  313 CP(I,J)=COEF*A33(I,J)
      I3=I+3
      CP(I3,I3)=COEF1
  314 CONTINUE
      COEF=C2/CETA1
      COEF1=C3*C4/CETA1
      CM3=CM/3.
      DO 315 I=1,6
      CMP(I)=-CM*WJ(I)
      IF (I.LE.3) CMP(I)=CMP(I)+CM3
  315 CONTINUE
      C2=C2/CETA1
      C3C4=C3*C4/CETA1
      DO 320 I=1,6
      C2I=C2*SE(I)
      C3I=C3C4*SE(I)
      DO 320 J=1,6
      ELAST(I,J)=CP(I,J)+C2I*SE(J)+C3I*WJ(J)
      IF (I.LE.3) ELAST(I,J)=ELAST(I,J)+CMP(J)
  320 CONTINUE
      GO TO 500
C
C        CAP YIELDING
C
  350 CONTINUE
      DENK=2.+CETA+AK*RSQ
      BMBAR=(B4*(ETADT2+3.*(CM+ETADT))+1.)/AI1L
      DBAR=RSQ*BKR*BMBAR/DENK
      B4=AJ2DE*DBAR
      COEFB=AK*CM*RSQ*BKR/(DENK*AI1L)
      COEF=AJ2DE*COEFB
      COEF1=2.*G*BKR
      DO 351 I=1,3
      CMP(I)=COEF+COEF1*SE(I)
  351 CMP(I+3)=COEF1*SE(I+3)
      BP=C*B1*DEXP(-B1*EL)+TETA
      B=(EL-AI1A)/R
      DEN=-2.*AI1L*((CM+ETADT)*DEPVPR*AI1AP+1.)-RSQ*(B4+B*BP)
      COEF=-2.*AI1L/DEN*CM
      COEF1=-RSQ/DEN
      DO 352 I=1,3
      WJ(I)=COEF+COEF1*CMP(I)
  352 WJ(I+3)=COEF1*CMP(I+3)
      C33=DVT/3.*2.*G
      C33D=2.*C33
      DO 353 I=1,3
      DO 353 J=1,3
      A33(I,J)=-C33
      IF (I.EQ.J) A33(I,J)=C33D
  353 CONTINUE
      COEFG=BKR*G
      DO 355 I=1,3
      I3=I+3
      DO 355 J=1,3
      CP(I,J)=-SE(I)*COEFG+BKR*A33(I,J)-DBAR*SE(I)*WJ(J)
      J3=J+3
      CP(I3,J3)=-DBAR*SE(I3)*WJ(J3)
      IF (I.EQ.J) CP(I3,J3)=CP(I3,J3)+COEFG
  355 CONTINUE
      AJBAR=-AK*CM/AI1L
      DO 358 I=1,3
      CMP(I)=AJBAR+BMBAR*WJ(I)
  358 CMP(I+3)=BMBAR*WJ(I+3)
      COEF1=2.*G*DTETA*RSQ      
      COEF2=COEF1*AK
      DO 360 I=1,3
      I3=I+3
      DO 360 J=1,3
      J3=J+3
      CP(I,J)=A33(I,J)-COEF1*SS(I)*CMP(J)-COEF2*CP(I,J)
      CP(I3,J3)=-COEF1*SS(I3)*CMP(J3)-COEF2*CP(I3,J3)
      IF (I.EQ.J) CP(I3,J3)=CP(I3,J3)+COEFG
  360 CONTINUE
      CM3=CM/3.
      CMB=CM*B4C
      DO 362 I=1,6
      CMP(I)=-CMB*WJ(I)
      IF (I.LE.3) CMP(I)=CMP(I)+CM3
  362 CONTINUE
      DO 365 I=1,6
      DO 365 J=1,6
      ELAST(I,J)=CP(I,J)
      IF (I.LE.3) ELAST(I,J)=ELAST(I,J)+CMP(J)
  365 CONTINUE
C
C     KORIGOVANJE VELICINA IZ PRETHODNOG KORAKA KAD JE POSTIGNUTA
C     KONVERGENCIJA
C
  500 CONTINUE
      DO 290 I=1,6
      SIGS1(I)=TAU(I)-ETADT*DEPVP(I)
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
      RETURN
      END
