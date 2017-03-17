C=======================================================================
C
CE   CONCRETE FOR 3/D ELEMENTS
CS   BETON 3/D ELEMENT     (22.12.1994)
C
C=======================================================================
      SUBROUTINE D3M27(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    LOKATIONS ON GAUSS LEVEL POINT
CS    LOKACIJE NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /CONFAC/ FACTOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' D3M27'
C
      LFUN=MREPER(1)
C
      LTAUT =LPOCG
      LDEFT =LTAUT  + 6
      LDEFPT=LDEFT  + 6
      LSIGOT=LDEFPT + 6
      LEPSCT=LSIGOT + 3
      LPRINT=LEPSCT + 3
      LVT   =LPRINT + 3
      LEANIT=LVT    + 9
      LTEQT =LEANIT + 4
      LIPLT =LTEQT  + 1
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 6
      LDEFP1=LDEF1  + 6
      LSIGO1=LDEFP1 + 6
      LEPSC1=LSIGO1 + 3
      LPRIN1=LEPSC1 + 3
      LV1   =LPRIN1 + 3
      LEANI1=LV1    + 9
      LTEQ1 =LEANI1 + 4
      LIPL1 =LTEQ1  + 1
C
      IF(FACTOR.GT.1.D-12) THEN
         CALL JEDNAK(PLAS1(LTAU1),PLAS1(LTAU1),1.D0/FACTOR,6)
         CALL JEDNAK(TAU,TAU,1.D0/FACTOR,6)
         CALL JEDNAK(ELAST,ELAST,1.D0/FACTOR,36)
	ENDIF
      CALL TUI327(PLAST(LIPLT),PLAST(LTAUT),PLAST(LDEFT),
     1            PLAST(LDEFPT),PLAST(LSIGOT),PLAST(LEPSCT),
     1            PLAST(LPRINT),PLAST(LVT),PLAST(LEANIT),PLAST(LTEQT),
     1            PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEF1),
     1            PLAS1(LDEFP1),PLAS1(LSIGO1),PLAS1(LEPSC1),
     1            PLAS1(LPRIN1),PLAS1(LV1),PLAS1(LEANI1),PLAS1(LTEQ1),
     1            A(LFUN),TAU,DEF,IRAC)
C
      IF(FACTOR.GT.1.D-12) THEN
         CALL JEDNAK(A(LTAU1),A(LTAU1),FACTOR,6)
         CALL JEDNAK(TAU,TAU,FACTOR,6)
         CALL JEDNAK(ELAST,ELAST,FACTOR,36)
	ENDIF
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TUI327( PLT,TAUT,DEFT,DEFPT,SIGOT,
     1                   EPSCT,PRINT,VT,EANIT,TEQT,
     1                   PL1,TAU1,DEF1,DEFP1,SIGO1,
     1                   EPSC1,PRIN1,V1,EANI1,TEQ,
     1                   FUN,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    INTEGRACIJA KONSTITUTIVNIH RELACIJA ZA BETON
CE    INTEGRATION OF CONSTITUTIVE RELATION FOR CONCRETE
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
      COMMON /CRACKC/ CBTC(10,2),NBTC(10,3),KBTC,ICRACK,NBRCR
      COMMON /PRINCI/ PRINC(3)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION TAUT(*),DEFT(*),DEFPT(*),SIGOT(*),
     1          EPSCT(*),PRINT(*),VT(3,*),EANIT(*),
     1          TAU1(*),DEF1(*),DEFP1(*),SIGO1(*),
     1          EPSC1(*),PRIN1(*),V1(3,*),EANI1(*),
     1          TAU(*),DEF(*)
      DIMENSION FUN(6,*),DL(6,6),TSG(6,6),TPRE(6),IEXTC(3),TN(100)
      DATA ITMAX/100/,EPSIL/1.0D-06/
      IF(IDEBUG.GT.0) PRINT *, ' TAUI327'
      IF(IRAC.EQ.2)RETURN
C
CC     ALF=-0.5D0
      ALF=0.D0
C
CE    INITIAL DATA
C
      IPLT=PLT
      IPL1=PL1
      IF(IPL1.GT.10003) STOP 'FATAL ERROR 1'
C
C
CE    CRUSHED POINTS: ZERO STRESS 
      IF(IPL1.EQ.10003)THEN
        CALL CLEAR(TAU,6)
        CALL CLEAR(ELAST,36)
        CALL JEDNA1(DEF1,DEF,6)
        CALL JEDNA1(TAU1,TAU,6)
        RETURN
      ENDIF
CE    ALL CRACKS CLOSED
      IF(IPL1.LT.10)THEN
        PL1=0.D0
        IPL1=0
      ENDIF
C
CE    DECODE INFORMATION ABOUT TYPES OF CRACKS
CE    IEXTC = 0   NO TENSILE CRACK   J, J=1,3
CE    IEXTC = 1   TENSILE CRACK 
      DO 5 I=1,3
        IEX=10**(4-I)
        II=IPL1/IEX
        IF(II.GT.0) IPL1=IPL1-II*IEX
        IEXTC(4-I)=II
5     CONTINUE
C      WRITE(3,*)'IPLT,IPL1,ICRACK,NBRCR,KBTC',
C     1           IPLT,IPL1,ICRACK,NBRCR,KBTC
C      WRITE(3,*)'PL1,IEXTC',PL1,IEXTC
      TEQY=TEQT
      E   =EANI1(1)
      ANI =EANI1(2)
      ES  =EANI1(3)
      ANIS=EANI1(4)
      FC  =FUN(1,MAT)
      SRF =FUN(2,MAT)
      E0  =FUN(3,MAT)
      ANI0=FUN(4,MAT)
      G0  =E0/2./(1.+ANI0)
C-----------------------------------------------------------------------
C.............................................
C.    A)   UNCRACKED, NONLINEAR CONCRETE
C.         OR FIRST CRACKING
C.............................................
C
      IF(IPL1.EQ.0.AND.ICRACK.GE.0)THEN
C
CE      CALCULATE ELASTIC SOLUTION
C
        CALL MEL27(ELAST,E0,ANI0,SRF,0,IEXTC)
        CALL ODUZ2B(DDEF,DEF,DEFPT,6)
        CALL CLEAR(TAU,6)
        CALL MNOZI1(TAU,ELAST,DDEF,6,6)
C
CE      CHECK FOR YIELDING
C
        TEQ=CONYLD(TAU,ALF)
        IF(TEQY.GT.1.D-20)THEN
          IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
            E   =E0
            ANI =ANI0
            ES  =E0
            ANIS=ANI0
            TEQ =TEQY
            CALL JEDNA1(DEFP1,DEFPT,6)
            CALL OKTAS(TAU,SIGO1)
            GO TO 500
          ENDIF
        ENDIF
C
CE      ITERATIONS
C     
        IT = 0
  100   IT=IT+1
C
        IF(IT.GT.ITMAX) THEN
          WRITE(IZLAZ,2000)TNORM0,TNORMI/TNORM0,EPSIL
          CALL WRR(TN,100,'NORM')
          STOP
        ENDIF
C
CE      OKTAL STRESS CORESPONDING TO PREVIOUS ITERATION STRESS
C
        CALL OKTAS(TAU,SIGO1)
C
CE      CONSTANTS AND CONSTITUTIVE MATRIX CORESPONDING TO PREVIOUS ITERATION
C
        CALL CONCON(FC,SIGO1,AK,G,AKS,GS,E,ANI,ES,ANIS)
        CALL MEL27(ELAST,E,ANI,SRF,0,IEXTC)
C
CE      TOTAL STRAINS CORESPONDING TO PREVIOUS ITERATION STRESS
C
        CALL GLODEF(TAU,ES,ANIS,GS,DEFDPR)
C
CE      NEW STRESS ESTIMATION
C
        CALL ODUZ2B(DDEF,DEF,DEFDPR,6)
        CALL CLEAR(DETAU,6)
        CALL MNOZI1(DETAU,ELAST,DDEF,6,6)
        TNORMI=RNORM(DETAU,6)
        CALL ZBIR2B(TAU,TAU,DETAU,6)
C
        IF(IT.EQ.1)THEN
           TNORM0=TNORMI
           IF(TNORM0.LT.1.D-15) TNORM0=1.D0
        ENDIF
        TN(IT)=TNORMI/TNORM0
        IF (TN(IT).GT.EPSIL) GO TO 100
C
 2000   FORMAT(' ','MAXIMUM NUMBER OF ITERATION EXIDED IN : TUI327'
     1        /'TNORM0        =',1PD12.5
     2        /'TNORMI/TNORM0 =',1PD12.5
     3        /'EPSIL         =',1PD12.5)
500     CONTINUE
        EANI1(1)=E
        EANI1(2)=ANI
        EANI1(3)=ES
        EANI1(4)=ANIS
C
CE      FIND PLASTIC STRAINS
C
        CALL GLODEF(TAU,E0,ANI0,G0,DEFDPR)
        CALL ODUZ2B(DEFP1,DEF,DEFDPR,6)
C
CE      DETERMINE IF CRACKING
C
        CALL GLAVN(TAU)
        IF(ICRACK.GE.0) CALL CRCCHK(FC,SIGO1,ICRK)
C
CE      UPDATE FROM PREVIOUS STEP
C
        CALL JEDNA1(DEF1,DEF,6)
        CALL JEDNA1(TAU1,TAU,6)
C
      ENDIF
C
      IF(IPL1.EQ.0.AND.ICRACK.LT.0) CALL MEL27(ELAST,E,ANI,SRF,0,IEXTC)
C-----------------------------------------------------------------------
C
C.............................................
C.    B)   CRACKED CONCRETE
C.............................................
C
CE    ALREDY CRACKED AND NO CHANGE : LINEAR SOLUTION
      IF(IPL1.GT.0)THEN
CE....  S T A R T     TO BLOCK DEACTIVATION OF CRACKS
CE      DEACTIVATE CRACK
        ICLOSE=0
        IF(PLT.GE.PL1.AND.ICRACK.GE.0)THEN
          CALL TRANSE(TSG,V1)
          CALL CLEAR(DEFDS,6)
          CALL MNOZI1(DEFDS,TSG,DEF,6,6)
          DO 700 I=1,IPL1
            TREL=EPSC1(I)-DEFDS(I)
            IF(IEXTC(I).EQ.1.AND.TREL.GT.0.D0)THEN
               ICLOSE=1
               WRITE(*,*)'* MARKED TO CLOSE CRACK '
               CALL CRACON(TREL,I,2)
            ENDIF
700       CONTINUE
        ENDIF
CE....  E N D    BLOCK DEACTIVATION OF CRACKS
CX       IF(IPL1.EQ.3)THEN
        IF(IPL1.EQ.3.AND.IEXTC(1).NE.0.AND.IEXTC(2).NE.0
     1                                .AND.IEXTC(3).NE.0)THEN
          CALL MEL27(ELAST,ES,ANIS,SRF,3,IEXTC)
          CALL CLEAR(TAU,6)
          CALL MNOZI1(TAU,ELAST,DEF,6,6)
CE        UPDATE FROM PREVIOUS STEP
          CALL JEDNA1(DEF1,DEF,6)
          CALL JEDNA1(TAU1,TAU,6)
          RETURN
        ELSE
          CALL MEL27(DL,E,ANI,SRF,IPL1,IEXTC)
          CALL TRANSE(TSG,V1)
          CALL TRAETP(DL,ELAST,TSG)
          IF(ICRACK.GE.0)THEN
CA           CALL CLEAR(TAU,6)
CA           CALL MNOZI1(TAU,ELAST,DEF,6,6)
            CALL ODUZ2B(DDEF,DEF,DEF1,6)
            CALL CLEAR(DETAU,6)
            CALL MNOZI1(DETAU,ELAST,DDEF,6,6)
            CALL ZBIR2B(TAU,TAU1,DETAU,6)
CA
            CALL OKTAS(TAU,SIGO1)
            IF(PLT.LE.PL1.AND.ICLOSE.EQ.0) CALL CRCCHK(FC,SIGO1,ICRK)
CE          UPDATE FROM PREVIOUS STEP
            CALL JEDNA1(DEF1,DEF,6)
            CALL JEDNA1(TAU1,TAU,6)
          ENDIF
        ENDIF
      ENDIF
C
CE    NEW CRACKING
C
      IF(ICRACK.LT.0)THEN
        IF(ICRACK.EQ.-1)THEN
CE        CRUSH
          PL1=10003.
          IPL1=3
          CALL CLEAR(DETAU,6)
          CALL CLEAR(DL,36)
          WRITE(*,*)'*** CRUSH ***'
        ELSEIF(ICRACK.EQ.-2)THEN
CE        NEW CRACK
          IF(IPL1.EQ.3)THEN
            IF(IEXTC(3).EQ.0)THEN
               IEXTC(3)=1
               IEPS=3
            ELSEIF(IEXTC(2).EQ.0)THEN
               IEXTC(2)=1
               IEPS=2
            ELSEIF(IEXTC(1).EQ.0)THEN
               IEXTC(1)=1
               IEPS=1
            ENDIF
            PL1=PL1+10.**IEPS
          ELSE
            IPL1=IPL1+1
            PL1 =PL1+10.**IPL1+1.
            IEXTC(IPL1)=1
            IEPS=IPL1
          ENDIF
          IF(IPL1.EQ.3)THEN
             WRITE(*,*)'*** TOTAL CRACKING ***'
          ELSE
             WRITE(*,*)'*** CRACK ***'
          ENDIF
          CALL GLAVN(TAU1)
          CALL GLAPR3(TAU1,V1)
          CALL MEL27(DL,E,ANI,SRF,IPL1,IEXTC)
          CALL TRANSE(TSG,V1)
          CALL TRAETP(DL,DL,TSG)
CE...     SET STRAINS IN PRINCIPAL STRESS DIRECTIONS; SAVE FOR UNCRACKING
          CALL CLEAR(DEFDS,6)
CC         CALL MNOZI1(DEFDS,TSG,DEF,6,6)
          CALL MNOZI1(DEFDS,TSG,DEFP1,6,6)
          EPSC1(IEPS)=DEFDS(IEPS)
C
          CALL CLEAR(TPRE,6)
          CALL JEDNA1(TPRE,PRINC,3)
CCCC       CALL CLEAR(TPRE,IPL1)
          DO 212 IZ=1,3
212       IF(IEXTC(IZ).EQ.1) TPRE(IZ)=0.D0
          CALL CLEAR(DETAU,6)
          CALL MNOZI2(DETAU,TSG,TPRE,6,6)
        ELSE
CE        CLOSE CRACK
          II=IABS(ICRACK)-2
          PL1=PL1-10.**II
          IEXTC(II)=0
          EPSC1(II)=0.D0
          CALL MEL27(DL,E,ANI,SRF,IPL1,IEXTC)
          CALL TRANSE(TSG,V1)
          CALL TRAETP(DL,DL,TSG)
          CALL ODUZ2B(DDEF,DEF,DEF1,6)
          CALL CLEAR(DETAU,6)
          CALL MNOZI1(DETAU,DL,DDEF,6,6)
          CALL ZBIR2B(DETAU,TAU1,DETAU,6)
          WRITE(*,*)'*** CLOSE CRACK ***'
        ENDIF
C
CC       CALL CLEAR(DETAU,6)
CC       CALL MNOZI1(DETAU,DL,DEF,6,6)
        CALL ODUZ2B(TAU,DETAU,TAU1,6)
        CALL JEDNA1(TAU1,DETAU,6)
        CALL ODUZ2B(ELAST,DL,ELAST,36)
      ENDIF
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MEL27(D,E,ANI,SRF,IPL,IEXTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    FORM ( D ) MATRIX
C
      DIMENSION D(6,*),IEXTC(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' MEL27'
      CALL CLEAR(D,36)
C
      I1=IPL+1
      D11=E*(1.-ANI)/(1.+ANI)/(1.-2.*ANI)
      D12=D11*ANI/(1.-ANI)
      D44=D11*(1.-2.*ANI)/(1.-ANI)/2.
      GO TO(1,2,3,4),I1
1     D(1,1)=D11
      D(1,2)=D12
      D(1,3)=D12
2     D(2,2)=D11
      D(2,3)=D12
3     D(3,3)=D11
C
4     D(4,4)=D44
      D(5,5)=D(4,4)
      D(6,6)=D(4,4)
CE  SHEAR RETENTION FACTOR
      GO TO(7,5,6,6),I1
5     D(5,5)=D(5,5)*SRF
      D(6,6)=D(6,6)*SRF
      GO TO 7
6     D(5,5)=D(5,5)*SRF
      D(6,6)=D(6,6)*SRF
      D(4,4)=D(4,4)*SRF
C
7     CONTINUE
      IF(IPL.GT.0)THEN
        IF(IEXTC(1).EQ.0)THEN
          D(1,1)=D11
          IF(IEXTC(2).EQ.0)THEN
             D(1,2)=D12
             D(4,4)=D44
          ENDIF
          IF(IEXTC(3).EQ.0)THEN
             D(1,3)=D12
             D(6,6)=D44
          ENDIF
        ENDIF
        IF(IEXTC(2).EQ.0)THEN
          D(2,2)=D11
          IF(IEXTC(3).EQ.0)THEN
             D(2,3)=D12
             D(5,5)=D44
          ENDIF
        ENDIF          
        IF(IEXTC(3).EQ.0) D(3,3)=D11
      ENDIF
C
      CALL SIMETR(D,6)
      RETURN
      END
C======================================================================
      SUBROUTINE OKTAS(TAU,SIGO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C--------------------------------------------------------------------
CE    OCTAHEDRAL STRSSES 
CE       NOTE: COMPRESS IS POSITIVE
C
C--------------------------------------------------------------------
C
      COMMON /PRINCI/ PRINC(3)
      DIMENSION TAU(*),SIGO(*),DEV(3)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OKTAS'
C
      ONE=0.9999999999999999D0
      AI1=TAU(1)+TAU(2)+TAU(3)
      AI2=TAU(1)*TAU(2)+TAU(2)*TAU(3)+TAU(3)*TAU(1) 
     1   -TAU(4)*TAU(4)-TAU(5)*TAU(5)-TAU(6)*TAU(6) 
      SIGO(1)=-AI1/3.
      DO 10 I=1,3
10    DEV(I)=TAU(I)+SIGO(1)
      AI3=DEV(1)*DEV(2)*DEV(3)+TAU(4)*TAU(5)*TAU(6)*2.
     1   -DEV(1)*TAU(5)*TAU(5)-DEV(2)*TAU(6)*TAU(6)-DEV(3)*TAU(4)*TAU(4)
      SIGO(2)=DSQRT(2.*SIGO(1)*SIGO(1)-2.*AI2/3.)
      IF(DABS(SIGO(2)).LT.1.D-15)THEN
        SIGO(3)=one
c        SIGO(3)=0.D0
      ELSE
        TETA3=DSQRT(2.D0)*AI3/(SIGO(2)*SIGO(2)*SIGO(2))
        IF(DABS(TETA3).GT.ONE) TETA3=DSIGN(ONE,TETA3)
        SIGO(3)=DCOS(DACOS(TETA3)/3.)
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE OKTAE(SIGO,EPSO,AKS,GS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C--------------------------------------------------------------------
CE    OCTAHEDRAL STRAINS
C
C--------------------------------------------------------------------
C
      DIMENSION SIGO(*),EPSO(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OKTAE'
C
      EPSO(1)=SIGO(1)/AKS/3.
      EPSO(2)=SIGO(2)/GS/2.
C NIJE DEFINISANO (A.7c)
      EPSO(3)=0.
      RETURN
      END
C=======================================================================
      SUBROUTINE CONCON(FC,SIGO,AK,G,AKS,GS,E,ANI,ES,ANIS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C--------------------------------------------------------------------
CE    FC: UNIAXIAL CONCRETE CYLINDER STRENGTH
CE                  20MPA  < FC < 70MPA
C--------------------------------------------------------------------
C
      DIMENSION SIGO(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' CONCON'
      RSO=SIGO(1)/FC
      RTO=SIGO(2)/FC
C
CE       COMPUTE K0 AND G0
C
         AK0=11000.+3.2*FC*FC
         G0 =9224.+136.*FC+3296.D-15*FC**8.273
C
CE       COMPUTE PARAMETERS A,B,C,D
C
         IF(FC.LE.31.7D0) THEN
            A=.516D0
            B=2.+1.81D-8*FC**4.461
            C=3.573D0
            D=2.12+.0183*FC
         ELSE
            FCC=FC-31.7
            A=.516/(1.+.0027*FCC**2.397)
            B=2.+1.81D-8*FC**4.461
            C=3.573/(1.+.0134*FCC**1.414)
            D=2.7D0
         ENDIF
C
CE    COMPUTE TANGENT MODULI AND SECANT MODULI
C
      ARSO=DABS(RSO)
      B1=B-1.
      D1=D-1.
      IF(ARSO.LT.2.D0) THEN
         AK=AK0/(1.+A*B*ARSO**B1)
      ELSE
         AK=AK0/(1.+(2.D0**B1)*A*B)
      ENDIF
      G=G0/(1.+D*C*RTO**D1)
      IF(ARSO.LT.2.D0) THEN
         AKS=AK0/(1.+A*ARSO**B1)
      ELSE
         AKS=AK0/(1.+(2.**B1)*A*B-(2.**B)*A*B1*(1./ARSO))
      ENDIF
      GS=G0/(1.+C*RTO**D1)
C
CE    YOUNG'S AND POISSON'S RATIO
C
      E=(9.*AK*G)/(3.*AK+G)
      ANI=(3.*AK-2.*G)/(6.*AK+2.*G)
      ES=(9.*AKS*GS)/(3.*AKS+GS)
      ANIS=(3.*AKS-2.*GS)/(6.*AKS+2.*GS)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE CRCCHK(FC,SIGO,ICRK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    STRENGTH ENVELOPE: CHECK FOR CRACKING
C
      COMMON /CRACKC/ CBTC(10,2),NBTC(10,3),KBTC,ICRACK,NBRCR
      COMMON /PRINCI/ PRINC(3)
      DIMENSION SIGO(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' CRCCHK'
C
      RSO=SIGO(1)/FC
      CT =SIGO(3)
      RSO5=RSO+0.05
      ICRK=0            
      IF(RSO5.LT.0.D0)THEN
        ICRK=-1
        TREL=1.D20
      ELSEIF(PRINC(3).LE.-FC.AND.
     1       PRINC(1).LE.0.D0.AND.PRINC(2).LE.0.D0)THEN
        ICRK=1
        TREL=2.D20
      ENDIF
C
      IF(ICRK.EQ.0)THEN
        TOC=FC*.944*RSO5**.724
        TOE=FC*.633*RSO5**.857
        T2M=2.*(TOC*TOC-TOE*TOE)*CT
        TE2=2.*TOE-TOC
        DUM1=2.*T2M*CT
        DUM2=TOC*(T2M+TE2*DSQRT(DUM1+5*TOE*TOE-4.*TOC*TOE))
        DUM3=DUM1+TE2*TE2
        TOU=DUM2/DUM3 
        IF(TOU-SIGO(2).LT.0.D0) ICRK=-1
        TREL=SIGO(2)/TOU      
      ENDIF
C
      IF(ICRK.NE.0)THEN
        IF(NBRCR.GT.0)THEN
           IF(ICRK.EQ.1)  WRITE(*,*)'* MARKED CRUSH'
           IF(ICRK.EQ.-1) WRITE(*,*)'* MARKED TENSILE CRACK'
        ENDIF
        CALL CRACON(TREL,ICRK,1)
      ENDIF
C
      RETURN
      END
C=======================================================================
      SUBROUTINE GLODEF(TAU,ES,ANIS,GS,DEFDPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C--------------------------------------------------------------------
CE    COMPUTE GLOBAL STRAINS
C--------------------------------------------------------------------
C
      DIMENSION TAU(*),DEFDPR(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' GLODEF'
C
      DEFDPR(1)=(TAU(1)-ANIS*(TAU(2)+TAU(3)))/ES
      DEFDPR(2)=(TAU(2)-ANIS*(TAU(1)+TAU(3)))/ES
      DEFDPR(3)=(TAU(3)-ANIS*(TAU(1)+TAU(2)))/ES
      DUM=1./GS
      DO 20 I=4,6
20    DEFDPR(I)=DUM*TAU(I)
      RETURN
      END
C=======================================================================
      SUBROUTINE CRACON(TREL,ICRK,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C--------------------------------------------------------------------
CE    ARANGE CRACKING POINTS
C--------------------------------------------------------------------
C
      COMMON /CRACKC/ CBTC(10,2),NBTC(10,3),KBTC,ICRACK,NBRCR
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' CRACON'
C
      ICRACK=ICRACK+1
      DD=0.D0
      DO 10 I=1,NBRCR
        DUM=TREL-CBTC(I,IND)
        IF(DUM.GT.DD)THEN
          DD=DUM
          II=I
        ENDIF
10    CONTINUE
      IF(DD.GT.0.D0)THEN
        CBTC(II,IND)=TREL
        IF(ICRK.GT.0)  NBTC(II,IND)= KBTC
        IF(IND.EQ.2)   NBTC(II,3)  = ICRK
	IF(ICRK.EQ.-1) NBTC(II,IND)=-KBTC
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION CONYLD(TAU,ALF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   LOADING-UNLOADING (YIELDING) FUNCTION
C
      DIMENSION TAU(*)
C
      SM3=TAU(1)+TAU(2)+TAU(3)
      SM =SM3/3.
      CONYLD=0.D0
      DO 10 I=1,3
      DUM=TAU(I)-SM
10    CONYLD=CONYLD+DUM*DUM
      CONYLD=CONYLD + 2.*(TAU(4)*TAU(4)+TAU(5)*TAU(5)+TAU(6)*TAU(6))
      CONYLD=DSQRT(CONYLD)+ALF*SM3
      RETURN
      END
