C=======================================================================
C
C   TERMO-PLASTICNOST LJUSKE  -  MESOVITO OJACANJE    (20.01.1994)
C
C    SUBROUTINE D8M14
C               TI814
C               MEL8T
C               CEP814
C
C=======================================================================
      SUBROUTINE D8M14(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE  SPACE IN WORKING VECTOR
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFPT=LDEFT  + 6*IDVA
      LALFAT=LDEFPT + 6*IDVA
      LTEQT =LALFAT + 6*IDVA
      LDQPT =LTEQT  + 1*IDVA
      LIPL  =LDQPT  + 1*IDVA
      LTHI  =LIPL   + 1*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6*IDVA
      LDEFP1=LDEFT1 + 6*IDVA
      LALFA1=LDEFP1 + 6*IDVA
      LTEQT1=LALFA1 + 6*IDVA
      LDQPT1=LTEQT1 + 1*IDVA
      LIPL1 =LDQPT1 + 1*IDVA
      LTHI1 =LIPL1  + 1*IDVA
C
      CALL TI814(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1           A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1           A(LTEQT1),A(LDQPT1),A(LTHI1),
     1           A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC,A(LDEFT))
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI814( PL ,DEFPT,ALFAT,TEMT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEM, DEFQP, THI1,
     1                   FUN,MATE,TAU,DEF,TGT,TREF,NTFUN,IRAC,DEFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   THERMO-ELASTOPLASTIC MATERIAL , MIXED HARDENING
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
C
      COMMON /INDNAP/ NAPON
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
C
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFT(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(6),SEHET(6),ETHERM(3),
     2          THI1(*),ALFATG(6),ALFATL(6)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*),COEFE(2)
C
      DIMENSION DDEFPS(6),POM(3,3),VDEF(3,3),TAU0(6),CT(3,3,3,3)
C
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C      IF(NLM.EQ.NE) IST=1
C     THICKNESS IN INTEGRATED POINT
C      IF(IATYP.GE.4.AND.NNSL.EQ.1) THEN
      IF(IATYP.GE.4) THEN
         IF(ITER.EQ.0) THEN
            IF(KOR.EQ.1) THI1(1)=THIC0
            THI1(2)=THI1(1)
         ENDIF
C        DEBLJINA IZ PRETHODNE ITERACIJE
         THICK=THI1(1)
C        DEBLJINA IZ PRETHODNOG KORAKA
         THICT=THI1(2)
         IF(IRAC.NE.2.AND.IST.EQ.1) THEN
            WRITE(3,*)'NLM,NX,NY,ITER,KOR',NLM,NGAUSX,NGAUSY,ITER,KOR
            WRITE(3,*)'T1,TT,T0',THI1(1),THI1(2),THIC0
         ENDIF
      ENDIF
C
      COEFE(1)=0.83333333D0
      COEFE(2)=COEFE(1)
C
CE  INITIAL DATA
C
      IPL =PL
C
C ?????????????????????????
      HS=0.D0
C      HS    =FUN(1,MAT,4)  
      ISIMO=0
      IF(DABS(HS).GT.1.D-10) ISIMO=1
C     PRIVREMENO
C     INDIKATOR ZA NAPONE (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
      NAPGL=1
C     INDIKATOR ZA TRANSFORMACIJU NA KOSIJEVE NAPONE (0-NE,1-DA)
      NAPKO=1
C     INDIKATOR ZA PLASTICNOST (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
C     VAZNO! PLASTICNE DEFORMACIJE MORAJU U LOKALNOM SISTEMU ZBOG DEBLJINE
      LOKAL=1
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 441
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      DVT   =2.D0/3.
      CZZ   =1.-V
      CNI   =CZZ-V
C?????
      AE    =(1.+V)/E
      CYAN23=DVT*CY*AN
C?????
      CM    =E/CNI
      CNI   =CNI/3./CZZ
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2    =E/(1.+V)
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
      G2C1  =G2*C1
      G2CNI =G2*CNI
C       write(3,*)'E, V',E,V
C
CE    YIELD STRESS
C
C      IBRAH=1
C      IF(IBRAH.EQ.1.AND.KOR.EQ.1) THEN
C         DEFQPT=0.000769D0
C      ENDIF
      IF(ISIMO.EQ.0) THEN
         TEQY =TEQY0+CY*(EM*DEFQPT)**AN
      ELSE
         TEQY =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT
      ENDIF
C      WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C      CALL WRR6(DEF,6,'DEFU')
C      CALL WRR6(TAU,6,'TAUU')
C
CE    ELASTIC CONSTUTIVE MATRIX
C
      CALL MEL8T(COEFE,ETP)
      IF (IRAC.EQ.2) RETURN
C      WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C      CALL WRR6(DEF,6,'DEFU')
C      CALL WRR6(TAU,6,'TAUU')
C
CE    THERMAL STRAIN
C
      TEM=TGT
      IF(KOR.GT.1.AND.ITER.EQ.0) TGT=TGT-TEMT+TEMP0
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C
C     ZA ITER=0, TERMICKE SILE ZA POMERANJA OD TERMICKIH DEFORMACIJA
      IF(ITER.EQ.0) GO TO 900
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
    5   DEFPT(I)=0.5*DEFPT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      DEFDS(3)=0.D0
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2)-2*ETH)+ETH
        DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNI*(DEF(2)-DEFPT(2))-CTH*ETH
        DEFDS(2)=-CNI*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))-CTH*ETH
        DO 10 I=4,6
   10     DEFDS(I)=0.5*DEF(I)-DEFPT(I)
      ELSE
        EMT = CNI*(DEF(1)+DEF(2))
        DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
        DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        DO 11 I=4,6
   11     DEFDS(I)=0.5*DEF(I)
      END IF
C      DEFDS(3)  =-DEFDS(1)-DEFDS(2)  ????????????????
C
C       WRITE(3,*) 'ETH=',ETH
C       WRITE(3,*) 'EMT=',EMT
C       WRITE(3,*) 'ISIMO=',ISIMO
C       WRITE(3,*) 'IGLPR=',IGLPR
C       WRITE(3,*) 'DEFDS'
C       WRITE(3,145) (DEFDS(I),I=1,6)
C 145   FORMAT(6(1PD18.9))
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      IF(IATYP.GE.4) THEN
         IF(LOKAL.EQ.0) THEN
C           Al=Qs*Ad
            CALL CLEAR(ALFATL,6)
            CALL MNOZI1(ALFATL,TSGN,ALFAT,6,6)
         ELSE
            CALL JEDNA1(ALFATL,ALFAT,6)
         ENDIF
      ELSE
         CALL JEDNA1(ALFATL,ALFAT,6)
      ENDIF
C
      DO 40 I=1,6
        TAUD(I) =G2*DEFDS(I)
        SEHET(I)=TAUD(I)-ALFATL(I)
   40   SHET(I) =SEHET(I)
      TAUD(3)  =-TAUD(1)-TAUD(2)
      SEHET(3) =-SEHET(1)-SEHET(2)
      SHET(3)  =SEHET(3)
      TEQ=DSQRT(1.5*TENDOT(SHET))
      TEQSEH=TEQ
C??
C      EMT0=EMT
C      EMT=-TAUD(3)/CM
C      IF(DABS(EMT0-EMT).GT.1.D-8) STOP 'EMT'
C      IF(IST.EQ.1) call wrr6(ALFAT,6,'ALFT')
C      IF(IST.EQ.1) call wrr6(ALFATL,6,'ALFL')
C      IF(IST.EQ.1) call wrr6(defds,6,'defs')
C      IF(IST.EQ.1) call wrr6(glitl,6,'glit')
C      IF(IST.EQ.1) call wrr6(taud,6,'taud')
C??
C
CE   2)  CHECK FOR YIELDING
C
cccc      IF(NGE.EQ.4.AND.NLM.GT.31.AND.NLM.LT.38) THEN
C      WRITE(3,*) 'NLM,JG,TEQ,TEQY',NLM,JG,TEQ,TEQY
C      CALL WRR6(DEF,6,'DEF ')
C      CALL WRR6(TAU,6,'TAU ')
C      CALL WRR6(DEFPT,6,'DEPT')
C      CALL WRR6(DEFDS,6,'DEFS')
C      CALL WRR6(SEHET,6,'SEHT')
C      CALL WRR6(SHET,6,'SHET')
cccc      ENDIF
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5) THEN
         TEQ  =TEQY
         DEFQP=DEFQPT
         CALL JEDNA1(DEFP,DEFPT,6)
         CALL CLEAR(DDEFPS,6)
         GO TO 500
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
C
CE    BISECTION
C
      DEFQP=DEFQPT
      IF (DEFQP.LE.1.D-8) DEFQP=1.D-8
      IF(ISIMO.EQ.0) THEN
         EALFA=ANCY*DEFQP**AN1
      ELSE
         EALFA=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
      IF (EALFA.LE.1.D-8) EALFA=1.D-8
C
      IB     = 0
      IT     = 0
      AF     = 3.D0
      FP     = TEQ-TEQY
      DDD    = 0.1*FP/EALFA
      FM     = 0.D0
      DEPBM  = 0.D0
      DEPBP  = 0.D0
      DDEFQP = DDD
C
  100 IT=IT+1
      IB1 = IB
C
      IF (IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
      DEFQP=DEFQPT+DDEFQP
      EALFA=ANCY*DEFQP**AN1
      DTEQ=EALFA*EM**AN
      CHET=DVT*EM1*EALFA
      CHETP=AN1*CHET/DEFQP
      CG2=G2+CHET
C
      TEQY=TEQY0+CY*(EM*DEFQP)**AN
      DLAM=1.5*DDEFQP/TEQY
      DEL=1.+CG2*DLAM
      P1 =1.+(G2C1+CHET)*DLAM
      P2 =G2CNI*DLAM
      P12=P1*P1-P2*P2
      SHET(1)=(P1*SEHET(1)+P2*SEHET(2))/P12
      SHET(2)=(P2*SEHET(1)+P1*SEHET(2))/P12
      SHET(3)=-SHET(1)-SHET(2)
      DO 30 I=4,6
   30   SHET(I)=SEHET(I)/DEL
      TEQ=DSQRT(1.5*TENDOT(SHET))
      FB=TEQ-TEQY
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE   4)   DEVIATORIC STRESS, PLASTIC STRAIN, BACK STRESS
C
C??????
      IF(ISIMO.EQ.0) THEN
C         IF(IST.EQ.1)WRITE(3,*)'DEFQP,DEFQPT,DDEFQP',DEFQP,DEFQPT,DDEFQP
         CC=CYAN23*DEFQP**AN1
         EM1CC =EM1*CC
      ELSE
         CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
         EM1CC =EM1*CC
      ENDIF
C??????
C
C ????????????
         DUM =1.D0+EM1CC*DLAM
         DO 160 I=1,6
  160    TAUD(I)=ALFATL(I)+DUM*SHET(I)
C ????????????
C
      DO 170 I=1,6
        DDEFPS(I) =DLAM*SHET(I)
  170 CONTINUE
C      DEFZ=DDEFPS(3)
      IF(IST.EQ.1) call wrr6(DDEFPS,6,'DEFP')
      IF(IATYP.GE.4) THEN
         CALL JEDNA1(ALFATG,DDEFPS,6)
C        Pd=QeT*Pl
         CALL CLEAR(DDEFPS,6)
         CALL MNOZI2(DDEFPS,TSGD,ALFATG,6,6)
         IF(LOKAL.EQ.0) CALL JEDNA1(ALFATG,DDEFPS,6)
      ELSE
         CALL JEDNA1(ALFATG,DDEFPS,6)
      ENDIF
      DO 171 I=1,6
C??        DEFP(I) =DEFPT(I)+DDEFPS(I)
         DEFP(I) =DEFPT(I)+ALFATG(I)
C??        ALFA1(I)=ALFAT(I)+CHET*DDEFPS
         ALFA1(I)=ALFAT(I)+EM1CC*ALFATG(I)
C         TAUD(I) =SHET(I)+ALFA1(I)
  171 CONTINUE
C
C       WRITE(3,*) 'TAUD'
C       WRITE(3,145) (TAUD(I),I=1,6)
C
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-2*ETH)+ETH
      ELSE
        EMT = -TAUD(3)/CM
      END IF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF (ISKNP.NE.2) THEN
        IF (IATYP.LT.4) THEN
          SEHET(1)=G2*(DEF(1)-EMT-DEFPT(1))-ALFAT(1)
          SEHET(2)=G2*(DEF(2)-EMT-DEFPT(2))-ALFAT(2)
        ELSE
          SEHET(1)=G2*(DEF(1)-EMT)-ALFATL(1)
          SEHET(2)=G2*(DEF(2)-EMT)-ALFATL(2)
        END IF
        SEHET(3)=-SEHET(1)-SEHET(2)
        TEQSEH=DSQRT(1.5*TENDOT(SEHET))
        FHETP=DTEQ+1.5*CG2+1.5*CHETP*DDEFQP
        IF (FHETP.LE.1.D-8) FHETP=1.D-8
      IF(INDCEP.EQ.0)
     1  CALL CEP814(SEHET,DLAM,TEQ,G2,CM,FHETP,TEQSEH,CHET,CHETP,DEL,
     &             DTEQ)
      END IF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
C
C       WRITE(3,*) 'ETH=',ETH
C       WRITE(3,*) 'EMT=',EMT
C       WRITE(3,*) 'CM=',CM
C       WRITE(3,*) 'TAUM=',TAUM
C       WRITE(3,*) 'TAUD'
C       WRITE(3,145) (TAUD(I),I=1,6)
C
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=0.D0
      DO 46 I=4,6
        TAU(I)=TAUD(I)
   46   DEFP(I) =2.*DEFP(I)
C
C NORMIRANA ELASTICNA DEFORMACIJA
C     OVO ZA MALE DEFORMACIJE PROVERI DA LI TREBA
      IF (IATYP.LT.4) THEN
        DEF(3)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-2*ETH)+DEFP(3)+ETH
      END IF
C
      IF(IATYP.GE.4) THEN
         DO 210 I=1,2
  210    DEF(I)=TAUD(I)*AE+EMT
         DEF(3)=-(TAU(1)+TAU(2))*V/E
         DO 211 I=4,6
  211    DEF(I)=TAUD(I)*AE*2.D0
      ENDIF
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
      IF(IST.EQ.1) call wrr6(tau,6,'tau ')
      IF(IST.EQ.1) call wrr6(defp,6,'defp')
      IF(IST.EQ.1) write(3,*) 'emt,e,V',emt,e,V
C
CE    UPDATE FROM PREVIOUS STEP
C
C OVO PROVERITI ZA MALE DEFORMACIJE
      CALL JEDNA1(SHET,DEF,6)
      IF(IST.EQ.1) call wrr6(SHET,6,'SHET')
      IF(IGLPR.EQ.1) THEN
C        NAPONI U LOKALNOM DEKARTOVOM
         CALL JEDNA1(ALFATG,TAU,6)
         TAU(3)=0.D0
         CALL JEDNA1(TAU0,TAU,6)
C        NAPONI U GLOBALNOM DEKARTOVOM SISTEMU
C        Sd=QeT*Sl
         CALL CLEAR(TAU,6)
         CALL MNOZI2(TAU,TSGD,TAU0,6,6)
C         call wrr(tau0,6,'tau0')
C         call wrr(tau,6,'tau ')
C         call wrr6(tsg,36,'tsg ')
         IF(NAPKO.EQ.1) THEN
CV            IF(IATYP.EQ.4.AND.ILEDE.EQ.0) THEN
C              GLAVNE VREDNOSTI
C              LAMBDA
CV               P1=DSQRT(PRINC(1))
CV               P2=DSQRT(PRINC(2))
CV               P3=DSQRT(PRINC(3))
C              LEVI KOSI-GRINOV DEFORMACIONI TENZOR (V)
CV               CALL DIJAD(VDEF,QP,QP,P1,P2,P3)
CS             TRANSF. ROTIRANI PIOLA KIRKOFOV - KOSIJEV NAPON 
CE             TRANSFORM. ROTATED PIOLA KIRCKOF - CAUCHY STRESS
C              s = V * S * V
CV               CALL PIOKOS(VDEF,TAU0)
CV               CALL JEDNA1(TAU,TAU0,6)
CV               CALL CEPMT(ELAST,CT,0)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  WRITE(3,*) 'NLM,JG,ITER',NLM,JG,ITER
C                  CALL WRR3(VDEF,9,'VDEF')
C                  CALL WRR6(ELAST,36,'ELAP')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTDP')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTGP')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTDP')
C               ENDIF
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
CV               CALL RRRRC(ELAST,CT,VDEF,1)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  CALL WRR6(ELAST,36,'ELAS')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTD ')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTG ')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTD ')
C               ENDIF
CV            ENDIF
            IF(ILEDE.EQ.1) THEN
C              GLAVNE VREDNOSTI
C              INVERZNO LAMBDA
               P1=1.D0/DSQRT(PRINC(1))
               P2=1.D0/DSQRT(PRINC(2))
               P3=1.D0/DSQRT(PRINC(3))
C              INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
               CALL DIJAD(POM,QP,QP,P1,P2,P3)
C              TENZOR ROTACIJE R
C              R = Fe * Ue**-1 
               CALL MNOZM1(VDEF,XJ,POM,3,3,3)
CS             TRANSF. UNAZAD ROTIRANI KOSIJEV - KOSIJEV NAPON 
CE             TRANSFORM. BACK ROTATED COUCHY - CAUCHY STRESS
C              s = R * S * RT
               CALL PIOKOS(VDEF,TAU)
CS             TRANSFORM. MATRICE ELAST. - ETP (LOKALNI - GLOBALNI DEKART.)
CE             TRANSFORM ELASTICITY MATRIX - ETP (LOCAL - GLOBAL CARTESIAN)
               CALL TRAETP(ETP,ELAST,TSGD)
               CALL CEPMT(ELAST,CT,0)
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
               CALL RRRRC(ELAST,CT,VDEF,1)
CS             TRANSFORM. MATRICE ELAST. - ETP (GLOBALN - LOKALNI DEKART.)
CE             TRANSFORM ELASTICITY MATRIX - ETP (GLOBAL - LOCAL CARTESIAN)
               CALL TRLETP(ELAST,ETP,TSGN)
            ENDIF
         ENDIF
         IF(NAPGL.EQ.0) CALL JEDNA1(TAU0,TAU,6)
C        NAPONI U LOKALNOM DEKARTOVOM SISTEMU
         CALL JEDNA1(TAU,ALFATG,6)
      ENDIF
C
C     NAPON I DEFORMACIJA ZA STAMPANJE
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
         DEF1(I)=DEF(I)
         IF(IGLPR.EQ.1) THEN
            TAU1(I)=TAU0(I)
         ELSE
            TAU1(I)=TAU(I)
         ENDIF
  290 CONTINUE
C
  900 IF(ITER.EQ.0) THEN
         ETH=ETH*(ETP(2,1)+ETP(2,2)+ETP(2,3))
         TAU(1)=TAU(1)-ETH
         TAU(2)=TAU(2)-ETH
C         TAU(3)=TAU(3)-ETH
      ENDIF
C
C      IF(NGE.EQ.4.AND.NLM.GT.31.AND.NLM.LT.38) THEN
C      WRITE(3,*) 'TAUM,ETH,CM,EMT)',TAUM,ETH,CM,EMT
C      CALL WRR6(DEF,6,'DEFI')
C      CALL WRR6(TAU,6,'TAUI')
C      CALL WRR6(DEFP,6,'DEPI')
C      ENDIF
C
C???????????????????????????????????????????????????
C
      IF(IST.EQ.1) call wrr6(def1,6,'def1')
      IF(IST.EQ.1) call wrr6(tau1,6,'tau1')
      IF(IATYP.GE.4) THEN
C         DEFN=DEF1(3)-DEFT(3)+DEFZ
         DEFN=DEF1(3)+DEFP(3)-DEFT(3)-DEFPT(3)
         THI1(1)=THI1(2)*DEXP(DEFN)
         IF(THI1(1).LT.0.1D0*THIC0) THEN
            WRITE(3,*) 'NLM,NGAUSX,NGAUSY,NGAUSZ,THI1',
     +                  NLM,NGAUSX,NGAUSY,NGAUSZ,THI1(1),THI1(2)
            THI1(1)=0.1D0*THIC0
         ENDIF
         IF(IST.EQ.1) WRITE(3,*) 'DEFN,THI1',DEFN,THI1(1),THI1(2)
      ENDIF
C
      IF(NAPON.EQ.1.AND.IATYP.GE.4) THEN
         IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN
C           PRIRASTAJ PLASTICNE U GLOBALNOM DEKARTOVOM
            CALL GLAVN(DDEFPS)
            CALL GLAPR3(DDEFPS,QP)
            CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
            IF(IST.EQ.1) CALL WRR3(PRINC,3,'PP  ') 
            IF(IST.EQ.1) CALL WRR3(QP,9,'QP  ') 
            CALL CLEAR(DEF,6)
            CALL JEDNA1(DEF,PRINC,3)
C            CALL JEDNA1(DEF,DDEFPS,6)
            RETURN
         ENDIF
C
CS       TRANSFORMACIJA DEFORMACIJA NA GLOBALNI SISTEM (SMIC.INZ.)
CE       TRANSFORM STRAIN TO GLOBAL COORDS.
C        ed=QsT*el
         CALL CLEAR(DEF,6) 
         CALL MNOZI2(DEF,TSGN,SHET,6,6)
         IF(IST.EQ.1) CALL WRR3(DEF,6,'DEFD')
C         call wrr(def,6,'def ')
C         call wrr6(tsgn,36,'tss ')
C
C        KORIGOVANJE NORMIRANOG ELAST. DEF. TENZORA Be NA KRAJU KORAKA
C
         IF(IATYP.EQ.4) THEN
            DO 300 I=1,3
  300       DEF(I)=2.D0*DEF(I)+1.D0
         ELSEIF(IATYP.EQ.5) THEN
            NESAD=0
            IF(NESAD.EQ.1) THEN
               DO 301 I=4,6
  301          DEF(I)=0.5D0*DEF(I)
               IF(IST.EQ.1) CALL WRR3(DEF,6,'DEFD')
               CALL GLAVN(DEF)
               CALL GLAPR3(DEF,QP)
               CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
               P1=PRINC(1)
               P2=PRINC(2)
               P3=PRINC(3)
               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
               IF(IST.EQ.1) CALL WRR3(POM,9,'QPN ') 
            ENDIF
            IF(NESAD.EQ.0) THEN
C               WRITE(3,*) 'NLM,NGAUSX,NGAUSY,KOR,ITER',
C     +                     NLM,NGAUSX,NGAUSY,KOR,ITER
C               CALL WRR3(PRINC,3,'PPNS') 
C               CALL WRR3(QP,9,'QPNS') 
CS             TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM 
CE             TRANSFORM STRAIN TO GLOBAL COORDS.
C              Sd=QeT*Sl
               CALL CLEAR(DEF,6) 
               CALL MNOZI2(DEF,TSGD,TAU,6,6)
               IF(IST.EQ.1) CALL WRR3(DEF,6,'STRG')
               CALL GLAVN(DEF)
               CALL GLAPR3(DEF,QP)
               CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
               IF(IST.EQ.1) CALL WRR3(QP,9,'QPN ') 
C               CALL WRR3(PRINC,3,'PPN ') 
C               CALL WRR3(QP,9,'QPN ') 
CE             ELASTIC LOGARITHMIC STRAIN IN PRINCIPAL DIRECTIONS
               P1=(PRINC(1)-V*(PRINC(2)+PRINC(3)))/E
               P2=(PRINC(2)-V*(PRINC(1)+PRINC(3)))/E
               P3=(PRINC(3)-V*(PRINC(2)+PRINC(1)))/E
               IF(IST.EQ.1) WRITE(3,*) 'E,V',E,V
            ENDIF
            IF(IST.EQ.1) WRITE(3,*) 'LOGDEF-P1,P2,P3',P1,P2,P3
CE          NEW LEFT CAUCHY-GREEN DEFORMATION TENSOR Be
            DEF(1)=DEXP(2.D0*P1)
            DEF(2)=DEXP(2.D0*P2)
            DEF(3)=DEXP(2.D0*P3)
            CALL DIJADS(QP,DEF)
            IF(IST.EQ.1) CALL WRR3(DEF,6,'BENG')
         ENDIF
      ENDIF
C
C???????????????????????????????????????????????????
C
      RETURN
  441 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI814')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI814'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI814')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI814'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE CEP814(SEHET,DLAM,TEQ,G2,CM,FHETP,TEQSEH,CHET,CHETP,
     &                  DEL,DTEQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ETP )    LJUSKA
CE     ELASTO-PLASTIC  CEP MATRIX        SHELL
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      DIMENSION SEHET(*),CP(6,6)
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      AMULT =G2/DEL
C
      AP=TRIPO*G2/FHETP/TEQSEH
      BP=(TRIPO-DLAM*DTEQ)/TEQ
      DP=AMULT*(BP-DLAM*DLAM*CHETP)/DEL*AP
      AA=AMULT*(1.+CHET*DLAM)
C
      DO 25 I=1,3
        DO 20 J=I,3
          ETP(I,J)=-DP*SEHET(I)*SEHET(J)
   20   CONTINUE
        ETP(I,I)=ETP(I,I)+AA
   25 CONTINUE
      DO 30 I=1,3
        DO 30 J=4,6
          CP(I,J)=-DP*SEHET(I)*SEHET(J)
   30 CONTINUE
      DO 40 I=4,6
        DO 45 J=I,6
   45     CP(I,J)=-DP*SEHET(I)*SEHET(J)
        CP(I,I)=CP(I,I)+0.5*AA
   40 CONTINUE
C
      CP(1,1)=TR*(DVA*ETP(1,1)-ETP(1,2)-ETP(1,3)+CM)
      CP(1,2)=TR*(DVA*ETP(1,2)-ETP(1,1)-ETP(1,3)+CM)
      CP(1,3)=TR*(DVA*ETP(1,3)-ETP(1,1)-ETP(1,2)+CM)
      CP(2,2)=TR*(DVA*ETP(2,2)-ETP(1,2)-ETP(2,3)+CM)
      CP(2,3)=TR*(DVA*ETP(2,3)-ETP(1,2)-ETP(2,2)+CM)
      CP(3,3)=TR*(DVA*ETP(3,3)-ETP(1,3)-ETP(2,3)+CM)
      CP(3,1)=CP(1,3)
      CP(3,2)=CP(2,3)
C
CE  STATIC CONDESATION
C
       DO 60 I=1,6
         ETP(I,3)=0.D0
         IF (I.EQ.3) GO TO 60
         COEFR = CP(3,I)/CP(3,3)
         DO 65 J=I,6
           IF (J.NE.3) ETP(I,J)=CP(I,J)-COEFR*CP(3,J)
           ETP(J,I)=ETP(I,J)
  65     CONTINUE
  60   CONTINUE
C      CALL WRR6(ETP,36,'ETP ')
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MEL8T(COEFE,ETP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /MATIZO/ E,V
      COMMON /PODTIP/ IPODT
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION ETP(6,*),COEFE(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MEL8T '
C
CE     NULL ( ETP )
C
      DO 13 I=1,6
        DO 13 J=1,6
   13     ETP(I,J)=0.D0
C
      EA=E/(1.D0-V*V)
      E1=0.5D0*(1.D0-V)
      ETP(1,1)=EA
      ETP(2,2)=EA
      ETP(2,1)=EA*V
      ETP(1,2)=ETP(2,1)
      ETP(4,4)=EA*E1
C      ETP(5,5)=COEFE(1)*ETP(4,4)
C      ETP(6,6)=COEFE(2)*ETP(4,4)
      ETP(5,5)=ETP(4,4)
      ETP(6,6)=ETP(4,4)
C      IF(IPODT.EQ.2) ETP(3,3)=ETP(4,4)
C
C      CALL WRR6(ETP,36,'ELAS')
      RETURN
      END
C=======================================================================
C
C   TERMO-ELASTICNOST SA PUZANJEM LJUSKE                   (26.04.1994)
C
C    SUBROUTINE D8M15
C               TI815
C               CEC815
C
C=======================================================================
      SUBROUTINE D8M15(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE  SPACE IN WORKING VECTOR
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFCT=LDEFT  + 6*IDVA
      LOPLUT=LDEFCT + 6*IDVA
      LOMINT=LOPLUT + 6*IDVA
      LTEFT =LOMINT + 6*IDVA
      LDQCT =LTEFT  + 1*IDVA
      LIOR  =LDQCT  + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
      LTEQT =LA2CT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 6*IDVA
      LDEFC1=LDEF1  + 6*IDVA
      LOPLU1=LDEFC1 + 6*IDVA
      LOMIN1=LOPLU1 + 6*IDVA
      LTEF1 =LOMIN1 + 6*IDVA
      LDQC1 =LTEF1  + 1*IDVA
      LIOR1 =LDQC1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
      LTEQT1=LA2C1  + 1*IDVA
C
      CALL TI815(A(LIOR),A(LDEFCT),A(LOPLUT),A(LOMINT),A(LTEFT),
     &            A(LDQCT),A(LTEQT),A(LTEQT1),
     &            A(LIOR1),A(LTAU1),A(LDEF1),A(LDEFC1),A(LOPLU1),
     &            A(LOMIN1),A(LTEF1),A(LDQC1),
     &            A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),
     &            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI815(ORI ,DEFCT,OPLUT,OMINT,TEFT ,DEFQCT,TEMT ,TEM  ,
     1                 ORI1,TAU1 ,DEF1 ,DEFC ,OPLUS,OMINS ,TEF  ,DEFQC,
     &                 A2CT,A2C  ,PTAUT,PTAU ,
     1                 FUN ,MATE ,TAU  ,DEF  ,TGT  ,TREF  ,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   THERMO-ELASTIC MATERIAL , WITH CREEP
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      DIMENSION DEFCT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFC(*),
     1          OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),ETHERM(3),SE(6)
      DIMENSION FUN(2,MATE*3,*),TREF(*),NTFUN(*),COEFE(2)
      DATA ITMAX/100/,EPSIL/1.0D-10/
      COEFE(1)=0.83333333D0
      COEFE(2)=COEFE(1)
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
C
CE  E,V,ALFA,TEMP0
C
      MAT3=(MAT-1)*3
      MATE3=MATE*3
      DO 6 J=1,3
        NFE=MAT3+J
        CALL BTAB(FUN,NTFUN,NFE,MATE3,TGT,NL,IND,2)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
C
CE    AUXILIARY CONSTANTS
C
      CZZ   =1.-V
      CNI   =CZZ-V
      CM    =E/CNI
      CNI   =CNI/3./CZZ
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2INV =(1.+V)/E
      G2C1  =C1/G2INV
      G2CNI =CNI/G2INV
      DVT   =2.D0/3.
      TRIPO =3.D0/2.
C
CE    ELASTIC CONSTUTIVE MATRIX
C
      CALL MEL8T(COEFE,ETP)
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
    5   DEFCT(I)=0.5*DEFCT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFCT(1)-DEFCT(2)-2*ETH)+ETH
        DEFDS(1)= C1 *(DEF(1)-DEFCT(1))-CNI*(DEF(2)-DEFCT(2))-CTH*ETH
        DEFDS(2)=-CNI*(DEF(1)-DEFCT(1))+C1 *(DEF(2)-DEFCT(2))-CTH*ETH
        DO 10 I=4,6
   10     DEFDS(I)=0.5*DEF(I)-DEFCT(I)
      ELSE
        EMT = CNI*(DEF(1)+DEF(2))
        DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
        DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        DO 11 I=4,6
   11     DEFDS(I)=0.5*DEF(I)
      END IF
      DEFDS(3)  =-DEFDS(1)-DEFDS(2)
C
CE    ELASTIC DEVIATORIC STRESS SOLUTION  (SE)
C
      DO 40 I=1,6
   40   SE(I) =DEFDS(I)/G2INV
c      SE(3) =-SE(1)-SE(2)
      TEF=DSQRT(1.5*TENDOT(SE))
C
      LPU=0
      IF (TEF.LT.1.D-12) THEN
        LPU=-1
        DEFQC=DEFQCT
        DO 600 I=1,6
          TAUD(I)=SE(I)
  600     DEFC(I)=DEFCT(I)
        GO TO 500
      END IF
C
CE   1)    OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      A2C    = A2CT
C
      P1      = 1.+G2C1*GAMDT
      P2      = G2CNI*GAMDT
      P12     = P1*P1-P2*P2
      S1      = P1*SE(1)+P2*SE(2)
      S2      = P2*SE(1)+P1*SE(2)
      TAUD(1) = S1/P12
      TAUD(2) = S2/P12
      TAUD(3) = -TAUD(1)-TAUD(2)
      DO 80 I=4,6
   80   TAUD(I) = DEFDS(I)/AGAMA
      ETEF    = DSQRT(1.5*TENDOT(TAUD))
      FL      = ETEF-TEFL
C
      IF (ITER.EQ.0) GO TO 50
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
      P1      = 1.+G2C1*GAMDT
      P2      = G2CNI*GAMDT
      P12     = P1*P1-P2*P2
      S1      = P1*SE(1)+P2*SE(2)
      S2      = P2*SE(1)+P1*SE(2)
      TAUD(1) = S1/P12
      TAUD(2) = S2/P12
      TAUD(3) = -TAUD(1)-TAUD(2)
      DO 85 I=4,6
   85   TAUD(I) = DEFDS(I)/AGAMA
      ETEF    = DSQRT(1.5*TENDOT(TAUD))
      F       = ETEF-TEF
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   2)   CREEP STRAIN
C
   50 DO 170 I=1,6
  170   DEFC(I) =DEFCT(I)+GAMDT*TAUD(I)
C
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+ETH
      ELSE
        EMT = -TAUD(3)/CM
      END IF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF (ISKNP.NE.2) THEN
        IF (IATYP.LT.4) THEN
          DEFDS(1)=DEF(1)-EMT-DEFCT(1)
          DEFDS(2)=DEF(2)-EMT-DEFCT(2)
          DEFDS(3)=CZZ*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+DEFC(3)+
     &             ETH-EMT-DEFCT(3)
        ELSE
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(3)=CZZ*(DEF(1)+DEF(2))-EMT
        END IF
        CALL CEC815(AGAMA,A2C,TEF,CM)
      END IF
C
CE   3)    CALCULATE STRESS
C
  500 CONTINUE
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=0.D0
      DO 46 I=4,6
        TAU(I)   = TAUD(I)
        DEFCT(I) = 2.*DEFCT(I)
   46   DEFC(I)  = 2.*DEFC(I)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,6
          OPLUS(I)=OPLUT(I)
  410     OMINS(I)=OMINT(I)
        CALL ORNL3(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1=IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF (IATYP.LT.4) THEN
        DEF(3)=CZZ*(DEF(1)+DEF(2)-DEFC(1)-DEFC(2)-2*ETH)+DEFC(3)+ETH
      ELSE
        DEF(3)=CZZ*(DEF(1)+DEF(2))
      END IF
      DO 290 I=1,6
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI815')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI815'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI815 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI815')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI815'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI815 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
       END
C=======================================================================
      SUBROUTINE CEC815(AGAMA,A2C,TEF,CM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEC ( ETP )    LJUSKA
CE     ELASTO-PLASTIC  CEC MATRIX        SHELL
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      DIMENSION CP(6,6)
C
      DVA   =2.D0
      TR    =1.D0/3.
      A2CDT =A2C*DT
C
      AA=1./AGAMA
      IF (TEF.LT.1.D-8) TEF=1.D-8
      DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
C
      DO 25 I=1,3
        DO 20 J=I,3
          ETP(I,J)=-DP*DEFDS(I)*DEFDS(J)
   20   CONTINUE
        ETP(I,I)=ETP(I,I)+AA
   25 CONTINUE
      DO 30 I=1,3
        DO 30 J=4,6
          CP(I,J)=-DP*DEFDS(I)*DEFDS(J)
   30 CONTINUE
      DO 40 I=4,6
        DO 45 J=I,6
   45     CP(I,J)=-DP*DEFDS(I)*DEFDS(J)
        CP(I,I)=CP(I,I)+0.5*AA
   40 CONTINUE
C
      CP(1,1)=TR*(DVA*ETP(1,1)-ETP(1,2)-ETP(1,3)+CM)
      CP(1,2)=TR*(DVA*ETP(1,2)-ETP(1,1)-ETP(1,3)+CM)
      CP(1,3)=TR*(DVA*ETP(1,3)-ETP(1,1)-ETP(1,2)+CM)
      CP(2,2)=TR*(DVA*ETP(2,2)-ETP(1,2)-ETP(2,3)+CM)
      CP(2,3)=TR*(DVA*ETP(2,3)-ETP(1,2)-ETP(2,2)+CM)
      CP(3,3)=TR*(DVA*ETP(3,3)-ETP(1,3)-ETP(2,3)+CM)
      CP(3,1)=CP(1,3)
      CP(3,2)=CP(2,3)
C
CE  STATIC CONDESATION
C
       DO 60 I=1,6
         ETP(I,3)=0.D0
         IF (I.EQ.3) GO TO 60
         COEFR = CP(3,I)/CP(3,3)
         DO 65 J=I,6
           IF (J.NE.3) ETP(I,J)=CP(I,J)-COEFR*CP(3,J)
           ETP(J,I)=ETP(I,J)
  65     CONTINUE
  60   CONTINUE
C
      RETURN
      END
C=======================================================================
C
C   TERMO-ELASTO-PLASTICNOST SA PUZANJEM -  LJUSKE     (21.04.1994)
C               (IZOTROPAN SA MESOVITIM OJACANJEM)
C
C    SUBROUTINE D8M16
C               TI816
C               EPC816
C
C=======================================================================
      SUBROUTINE D8M16(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFPT=LDEFT  + 6*IDVA
      LALFAT=LDEFPT + 6*IDVA
      LDEFCT=LALFAT + 6*IDVA
      LOPLUT=LDEFCT + 6*IDVA
      LOMINT=LOPLUT + 6*IDVA
      LTEQT =LOMINT + 6*IDVA
      LTEFT =LTEQT  + 1*IDVA
      LDQPT =LTEFT  + 1*IDVA
      LDQCT =LDQPT  + 1*IDVA
      LIPL  =LDQCT  + 1*IDVA
      LIOR  =LIPL   + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
      LTHI  =LA2CT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 6*IDVA
      LDEFP1=LDEF1  + 6*IDVA
      LALFA1=LDEFP1 + 6*IDVA
      LDEFC1=LALFA1 + 6*IDVA
      LOPLU1=LDEFC1 + 6*IDVA
      LOMIN1=LOPLU1 + 6*IDVA
      LTEQ1 =LOMIN1 + 6*IDVA
      LTEF1 =LTEQ1  + 1*IDVA
      LDQP1 =LTEF1  + 1*IDVA
      LDQC1 =LDQP1  + 1*IDVA
      LIPL1 =LDQC1  + 1*IDVA
      LIOR1 =LIPL1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
      LTHI1 =LA2C1  + 1*IDVA
C
      CALL TI816(A(LIOR) ,A(LDEFCT),A(LOPLUT),A(LOMINT),A(LDQCT),
     &           A(LIPL) ,A(LDEFPT),A(LALFAT),A(LDQPT) ,A(LTEQT),
     &           A(LIOR1),A(LDEFC1),A(LOPLU1),A(LOMIN1),A(LDQC1),
     &           A(LIPL1),A(LDEFP1),A(LALFA1),A(LDQP1) ,A(LTEQ1),
     &           A(LTAU1),A(LDEF1) ,A(LTEFT) ,A(LTEF1) ,
     &           A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),A(LTHI1),
     &           A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC,A(LDEFT))
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI816(ORI  ,DEFCT,OPLUT,OMINT ,DEFQCT,
     &                  PL  ,DEFPT,ALFAT,DEFQPT,TEMT  ,
     &                  ORI1,DEFC ,OPLUS,OMINS ,DEFQC ,
     &                  PL1 ,DEFP ,ALFA1,DEFQP ,TEM   ,
     &                  TAU1,DEF1 ,TEFT ,TEF   ,
     &                  A2CT,A2C  ,PTAUT, PTAU  , THI1,
     &                  FUN ,MATE ,TAU ,DEF, TGT ,TREF,NTFUN,IRAC,DEFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTICAN MATERIJAL SA PUZANJEM
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
C
      COMMON /INDNAP/ NAPON
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
C
      DIMENSION TAU(*),DEF(*),TAU1(*),DEF1(*),ETHERM(3),DEFT(*),
     &          DEFCT(*),DEFC(*),OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),
     &          DEFPT(*),DEFP(*),ALFAT(*),ALFA1(*),SHET(6) ,GLITL(6),
     &          AHET(6),THI1(*),ALFATG(6),ALFATL(6)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*),COEFE(2)
C
      DIMENSION DDEFPS(6),POM(3,3),VDEF(3,3),TAU0(6),CT(3,3,3,3)
C
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C      IF(NLM.EQ.NE) IST=1
C     THICKNESS IN INTEGRATED POINT
C      IF(IATYP.GE.4.AND.NNSL.EQ.1) THEN
      IF(IATYP.GE.4) THEN
         IF(ITER.EQ.0) THEN
            IF(KOR.EQ.1) THI1(1)=THIC0
            THI1(2)=THI1(1)
         ENDIF
C        DEBLJINA IZ PRETHODNE ITERACIJE
         THICK=THI1(1)
C        DEBLJINA IZ PRETHODNOG KORAKA
         THICT=THI1(2)
         IF(IRAC.NE.2.AND.IST.EQ.1) THEN
            WRITE(3,*)'NLM,NX,NY,ITER,KOR',NLM,NGAUSX,NGAUSY,ITER,KOR
            WRITE(3,*)'T1,TT,T0',THI1(1),THI1(2),THIC0
         ENDIF
      ENDIF
C
      COEFE(1)=0.83333333D0
      COEFE(2)=COEFE(1)
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
      IPL =PL
C
C ?????????????????????????
      HS=0.D0
C      HS    =FUN(1,MAT,4)  
      ISIMO=0
      IF(DABS(HS).GT.1.D-10) ISIMO=1
C     PRIVREMENO
C     INDIKATOR ZA NAPONE (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
      NAPGL=1
C     INDIKATOR ZA TRANSFORMACIJU NA KOSIJEVE NAPONE (0-NE,1-DA)
      NAPKO=1
C     INDIKATOR ZA PLASTICNOST (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
C     VAZNO! PLASTICNE DEFORMACIJE MORAJU U LOKALNOM SISTEMU ZBOG DEBLJINE
      LOKAL=1
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 331
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      CZZ   =1.-V
      CNI   =CZZ-V
      CM    =E/CNI
      CNI   =CNI/3./CZZ
C?????
      AE    =(1.+V)/E
      CYAN23=DVT*CY*AN
C?????
      CZZ   =-V/CZZ
      C1    =1.-CNI
      CTH   =C1-CNI
      G2INV =(1.+V)/E
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
      DVT   =2.D0/3.
      TRIPO =3.D0/2.
C
CE    YIELD STRESS
C
C      IBRAH=1
C      IF(IBRAH.EQ.1.AND.KOR.EQ.1) THEN
C         DEFQPT=0.000769D0
C      ENDIF
      IF(ISIMO.EQ.0) THEN
         TEQY =TEQY0+CY*(EM*DEFQPT)**AN
      ELSE
         TEQY =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT
      ENDIF
C      WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C      CALL WRR6(DEF,6,'DEFU')
C      CALL WRR6(TAU,6,'TAUU')
C
CE    ELASTIC CONSTUTIVE MATRIX
C
      CALL MEL8T(COEFE,ETP)
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
        DEFPT(I)=0.5*DEFPT(I)
    5   DEFCT(I)=0.5*DEFCT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2)-DEFCT(1)-DEFCT(2)-
     &        2*ETH)+ETH
        DEFDS(1)=  C1*(DEF(1)-DEFPT(1)-DEFCT(1))
     &           -CNI*(DEF(2)-DEFPT(2)-DEFCT(2))-CTH*ETH
        DEFDS(2)=-CNI*(DEF(1)-DEFPT(1)-DEFCT(1))
     &           + C1*(DEF(2)-DEFPT(2)-DEFCT(2))-CTH*ETH
        DO 10 I=4,6
   10     DEFDS(I)=0.5*DEF(I)-DEFPT(I)-DEFCT(I)
      ELSE
        EMT = CNI*(DEF(1)+DEF(2))
        DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
        DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
        DO 11 I=4,6
   11     DEFDS(I)=0.5*DEF(I)
      END IF
      DEFDS(3)=0.D0
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      IF(IATYP.GE.4) THEN
         IF(LOKAL.EQ.0) THEN
C           Al=Qs*Ad
            CALL CLEAR(ALFATL,6)
            CALL MNOZI1(ALFATL,TSGN,ALFAT,6,6)
         ELSE
            CALL JEDNA1(ALFATL,ALFAT,6)
         ENDIF
      ELSE
         CALL JEDNA1(ALFATL,ALFAT,6)
      ENDIF
C
      AGAMA = G2INV
      DO 40 I=1,6
        TAUD(I) =DEFDS(I)/AGAMA
        SHET(I) =TAUD(I)-ALFATL(I)
   40   GLITL(I)=SHET(I)*AGAMA
      TAUD(3) =-TAUD(1)-TAUD(2)
      SHET(3) =-SHET(1)-SHET(2)
      GLITL(3)=SHET(3)*AGAMA
      TEQ=DSQRT(TRIPO*TENDOT(SHET))
      TEQE=TEQ
      TEF=DSQRT(TRIPO*TENDOT(TAUD))
C??
C      EMT0=EMT
C      EMT=-TAUD(3)/CM
C      IF(DABS(EMT0-EMT).GT.1.D-8) STOP 'EMT'
C
CE   2)  CHECK FOR YIELDING
C
      LPU=0
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5.OR.ITER.EQ.0) THEN
        IPL1=0
        DEFQP=DEFQPT
C
        TEQ  =TEQY
        CALL JEDNA1(DEFP,DEFPT,6)
        CALL CLEAR(DDEFPS,6)
C
        IF (TEF.LT.1.D-10) THEN
          LPU=-1
          DEFQC=DEFQCT
          DO 601 I=1,6
  601       DEFC(I)=DEFCT(I)
          GO TO 500
        END IF
      ELSE
        PL1 =1.0D0
        IPL1=1
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC WITH CREEP.
C
CE         OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      DQMIN = 0.D0
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
      A2C    = A2CT
C
      QGAM    = G2INV+C1*GAMDT
      CGAM    = CNI*GAMDT
      QC2     = QGAM*QGAM-CGAM*CGAM
      S1      = QGAM*DEFDS(1)+CGAM*DEFDS(2)
      S2      = CGAM*DEFDS(1)+QGAM*DEFDS(2)
      TAUD(1) = S1/QC2
      TAUD(2) = S2/QC2
      TAUD(3) = -TAUD(1)-TAUD(2)
      DO 15 I=4,6
   15   TAUD(I) = DEFDS(I)/AGAMA
      ETEF    = DSQRT(1.5*TENDOT(TAUD))
      FL      = ETEF-TEFL
C
      IF (ITER.EQ.0) GO TO 50
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
CE    FINDING THE RADIUS OF THE YIELD SURFACE (TEQ)
C
        DLAM  = 0.D0
        DEFQP = DEFQPT
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        IF (DEFQP.LT.1.D-8) DEFQP=1.D-8
C
      IF(ISIMO.EQ.0) THEN
         EALFA=ANCY*DEFQP**AN1
      ELSE
         EALFA=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
        IF (EALFA.LE.1.D-8) EALFA=1.D-8
C
        KB    = 0
        KT    = 0
        DEPL  = 0.D0
        DDEP  = 0.1*(TEQE-TEQY)/EALFA
C
        QGAM    = G2INV+C1*GAMDT
        CGAM    = CNI*GAMDT
        QC2     = QGAM*QGAM-CGAM*CGAM
        GLITL(1)= DEFDS(1)-QGAM*ALFAT(1)+CGAM*ALFAT(2)
        GLITL(2)= DEFDS(2)+CGAM*ALFAT(1)-QGAM*ALFAT(2)
        DO 33 I=4,6
   33     GLITL(I)= DEFDS(I)-ALFAT(I)*AGAMA
        S1H     = QGAM*GLITL(1)+CGAM*GLITL(2)
        S2H     = CGAM*GLITL(1)+QGAM*GLITL(2)
        SHET(1) = S1H/QC2
        SHET(2) = S2H/QC2
        SHET(3) = -SHET(1)-SHET(2)
        DO 22 I=4,6
   22     SHET(I) = GLITL(I)/AGLHET
        TEQ    = DSQRT(1.5*TENDOT(SHET))
        FHETL   = TEQ-TEQY
        IF (FHETL.LE.0.D0) THEN
          IPL1=0
          GO TO 260
        END IF
        IPL1=1
        DDEFQP= DDEP
C
  250   KT    = KT+1
        KB1   = KB
C
        IF(KT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2015)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6015)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        DEFQP = DEFQPT+DDEFQP
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        EALFA = ANCY*DEFQP**AN1
        CHET  = DVT*EM1*EALFA
        AGPOT = 1.+AGAMA*CHET
C
        DLAM    = TRIPO*DDEFQP/TEQY
        AGLHET  = AGAMA+AGPOT*DLAM
        P1      = QGAM+(C1+QGAM*CHET)*DLAM
        P2      = CGAM+(CNI+CGAM*CHET)*DLAM
        P12     = P1*P1-P2*P2
        S1H     = P1*GLITL(1)+P2*GLITL(2)
        S2H     = P2*GLITL(1)+P1*GLITL(2)
        SHET(1) = S1H/P12
        SHET(2) = S2H/P12
        SHET(3) = -SHET(1)-SHET(2)
        DO 44 I=4,6
   44     SHET(I) = GLITL(I)/AGLHET
        TEQ    = DSQRT(1.5*TENDOT(SHET))
        FHET    = TEQ-TEQY
C
        CALL BISECB(DDEFQP,DEPL,DEPD,DDEP,
     &              FHET,FHETL,FHETD,AF,KB,DQMIN,DQTOL)
        IF (KB.EQ.-1) GO TO 260
        IF (KB1.EQ.0) GO TO 250
      IF (DABS(DDEP).GT.EPSIL.AND.
     1      (DABS(DDEP)/(DEPL+DEPD)).GT.EPSIL) GO TO 250
C
  260   DX      = DEFDS(1)-C1*DLAM*SHET(1)+CNI*DLAM*SHET(2)
        DY      = DEFDS(2)+CNI*DLAM*SHET(1)-C1*DLAM*SHET(2)
        S1      = QGAM*DX+CGAM*DY
        S2      = CGAM*DX+QGAM*DY
C??????
      IF(ISIMO.EQ.0) THEN
         CC=CYAN23*DEFQP**AN1
         EM1CC =EM1*CC
      ELSE
         CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
         EM1CC =EM1*CC
      ENDIF
C??????
C
C ????????????
C         DUM =1.D0+EM1CC*DLAM
C         DO 160 I=1,6
C  160    TAUD(I)=ALFATL(I)+DUM*SHET(I)
C ????????????
C
        TAUD(1) = S1/QC2
        TAUD(2) = S2/QC2
        TAUD(3) = -TAUD(1)-TAUD(2)
        DO 55 I=4,6
   55     TAUD(I) = (DEFDS(I)-DLAM*SHET(I))/AGAMA
        ETEF    = DSQRT(1.5*TENDOT(TAUD))
        F       = ETEF-TEF
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   4)   DEVIATORIC STRESS,PLASTIC STRAIN,BACK STRESS AND CREEP STRAIN
C
C
   50 DO 161 I=1,6
        DDEFPS(I) =DLAM*SHET(I)
  161 CONTINUE
C
      IF(IST.EQ.1) call wrr6(DDEFPS,6,'DEFP')
      IF(IATYP.GE.4) THEN
         CALL JEDNA1(ALFATG,DDEFPS,6)
C        Pd=QeT*Pl
         CALL CLEAR(DDEFPS,6)
         CALL MNOZI2(DDEFPS,TSGD,ALFATG,6,6)
         IF(LOKAL.EQ.0) CALL JEDNA1(ALFATG,DDEFPS,6)
      ELSE
         CALL JEDNA1(ALFATG,DDEFPS,6)
      ENDIF
C
C   50 DO 165 I=1,6
      DO 165 I=1,6
        IF ((IPL1.EQ.1).AND.(ITER.NE.0)) THEN
C          DDEFPS(I)   = DLAM*SHET(I)
C          DEFP(I)  = DEFPT(I)+DDEFPS(I)
          DEFP(I) =DEFPT(I)+ALFATG(I)
C          ALFA1(I) = ALFAT(I)+CHET*DDEFPS(I)
          ALFA1(I)=ALFAT(I)+EM1CC*ALFATG(I)
        ELSE
          DEFP(I)  = DEFPT(I)
        END IF
        DEFC(I) = DEFCT(I)+GAMDT*TAUD(I)
  165 CONTINUE
C
      IF (IATYP.LT.4) THEN
        EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &        2*ETH)+ETH
      ELSE
        EMT = -TAUD(3)/CM
      END IF
C
CE     E L A S T I C  -  P L A S T I C  -  C R E E P   M A T R I X  CEPC
C
      IF (ISKNP.NE.2) THEN
        IF (IATYP.LT.4) THEN
          DEFDS(1)=DEF(1)-EMT-DEFPT(1)-DEFCT(1)
          DEFDS(2)=DEF(2)-EMT-DEFPT(2)-DEFCT(2)
          DEFDS(3)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &             2*ETH)+DEFP(3)+DEFC(3)+ETH-EMT-DEFPT(3)-DEFCT(3)
        ELSE
          DEFDS(1)=DEF(1)-EMT
          DEFDS(2)=DEF(2)-EMT
          DEFDS(3)=CZZ*(DEF(1)+DEF(2))-EMT
        END IF
        IF ((IPL1.EQ.1).AND.(ITER.NE.0)) THEN
          DO 245 I=1,6
            AHET(I)  = ALFAT(I)+CHET*DEFDS(I)
  245       GLITL(I) = DEFDS(I)-ALFAT(I)*AGAMA
          GNORM = DSQRT(TRIPO*TENDOT(GLITL))
        END IF
        CALL EPC816(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &              AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQY,
     &              ALFAT)
      END IF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      IF (IPL1.NE.1) THEN
        DO 600 I=1,6
  600     DEFP(I)=DEFPT(I)
        TEQ=TEQE
      END IF
C RADOVAN
C       IF (IATYP.GE.4) THEN
C          TAUM=CM*(AJOT-1.D0)/3.D0
C       ELSE
         TAUM=CM*(EMT-ETH)
C       END IF
C RADOVAN
      TAU(1)=TAUD(1)+TAUM
      TAU(2)=TAUD(2)+TAUM
      TAU(3)=0.D0
      DO 66 I=4,6
        TAU(I)  = TAUD(I)
        DEFP(I) = 2.*DEFP(I)
   66   DEFC(I) = 2.*DEFC(I)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,6
          OPLUS(I)=OPLUT(I)
  410     OMINS(I)=OMINT(I)
        CALL ORNL3(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1=IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      IF (IATYP.LT.4) THEN
        DEF(3)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2)-DEFC(1)-DEFC(2)-
     &         2*ETH)+DEFP(3)+DEFC(3)+ETH
C      ELSE
C        DEF(3)=CZZ*(DEF(1)+DEF(2))
      END IF
C
      IF(IATYP.GE.4) THEN
         DO 210 I=1,2
  210    DEF(I)=TAUD(I)*AE+EMT
         DEF(3)=-(TAU(1)+TAU(2))*V/E
         DO 211 I=4,6
  211    DEF(I)=TAUD(I)*AE*2.D0
      ENDIF
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
      IF(IST.EQ.1) call wrr6(tau,6,'tau ')
      IF(IST.EQ.1) call wrr6(defp,6,'defp')
      IF(IST.EQ.1) write(3,*) 'emt,e,V',emt,e,V
C
CE    UPDATE FROM PREVIOUS STEP
C
C OVO PROVERITI ZA MALE DEFORMACIJE
      CALL JEDNA1(SHET,DEF,6)
      IF(IST.EQ.1) call wrr6(SHET,6,'SHET')
      IF(IGLPR.EQ.1) THEN
C        NAPONI U LOKALNOM DEKARTOVOM
         CALL JEDNA1(ALFATG,TAU,6)
         TAU(3)=0.D0
         CALL JEDNA1(TAU0,TAU,6)
C        NAPONI U GLOBALNOM DEKARTOVOM SISTEMU
C        Sd=QeT*Sl
         CALL CLEAR(TAU,6)
         CALL MNOZI2(TAU,TSGD,TAU0,6,6)
C         call wrr(tau0,6,'tau0')
C         call wrr(tau,6,'tau ')
C         call wrr6(tsg,36,'tsg ')
         IF(NAPKO.EQ.1) THEN
CV            IF(IATYP.EQ.4.AND.ILEDE.EQ.0) THEN
C              GLAVNE VREDNOSTI
C              LAMBDA
CV               P1=DSQRT(PRINC(1))
CV               P2=DSQRT(PRINC(2))
CV               P3=DSQRT(PRINC(3))
C              LEVI KOSI-GRINOV DEFORMACIONI TENZOR (V)
CV               CALL DIJAD(VDEF,QP,QP,P1,P2,P3)
CS             TRANSF. ROTIRANI PIOLA KIRKOFOV - KOSIJEV NAPON 
CE             TRANSFORM. ROTATED PIOLA KIRCKOF - CAUCHY STRESS
C              s = V * S * V
CV               CALL PIOKOS(VDEF,TAU0)
CV               CALL JEDNA1(TAU,TAU0,6)
CV               CALL CEPMT(ELAST,CT,0)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  WRITE(3,*) 'NLM,JG,ITER',NLM,JG,ITER
C                  CALL WRR3(VDEF,9,'VDEF')
C                  CALL WRR6(ELAST,36,'ELAP')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTDP')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTGP')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTDP')
C               ENDIF
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
CV               CALL RRRRC(ELAST,CT,VDEF,1)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  CALL WRR6(ELAST,36,'ELAS')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTD ')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTG ')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTD ')
C               ENDIF
CV            ENDIF
            IF(ILEDE.EQ.1) THEN
C              GLAVNE VREDNOSTI
C              INVERZNO LAMBDA
               P1=1.D0/DSQRT(PRINC(1))
               P2=1.D0/DSQRT(PRINC(2))
               P3=1.D0/DSQRT(PRINC(3))
C              INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
               CALL DIJAD(POM,QP,QP,P1,P2,P3)
C              TENZOR ROTACIJE R
C              R = Fe * Ue**-1 
               CALL MNOZM1(VDEF,XJ,POM,3,3,3)
CS             TRANSF. UNAZAD ROTIRANI KOSIJEV - KOSIJEV NAPON 
CE             TRANSFORM. BACK ROTATED COUCHY - CAUCHY STRESS
C              s = R * S * RT
               CALL PIOKOS(VDEF,TAU)
CS             TRANSFORM. MATRICE ELAST. - ETP (LOKALNI - GLOBALNI DEKART.)
CE             TRANSFORM ELASTICITY MATRIX - ETP (LOCAL - GLOBAL CARTESIAN)
               CALL TRAETP(ETP,ELAST,TSGD)
               CALL CEPMT(ELAST,CT,0)
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
               CALL RRRRC(ELAST,CT,VDEF,1)
CS             TRANSFORM. MATRICE ELAST. - ETP (GLOBALN - LOKALNI DEKART.)
CE             TRANSFORM ELASTICITY MATRIX - ETP (GLOBAL - LOCAL CARTESIAN)
               CALL TRLETP(ELAST,ETP,TSGN)
            ENDIF
         ENDIF
         IF(NAPGL.EQ.0) CALL JEDNA1(TAU0,TAU,6)
C        NAPONI U LOKALNOM DEKARTOVOM SISTEMU
         CALL JEDNA1(TAU,ALFATG,6)
      ENDIF
C
C     NAPON I DEFORMACIJA ZA STAMPANJE
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
         DEF1(I)=DEF(I)
         IF(IGLPR.EQ.1) THEN
            TAU1(I)=TAU0(I)
         ELSE
            TAU1(I)=TAU(I)
         ENDIF
  290 CONTINUE
C      IF(NGE.EQ.4.AND.NLM.GT.31.AND.NLM.LT.38) THEN
C      WRITE(3,*) 'TAUM,ETH,CM,EMT)',TAUM,ETH,CM,EMT
C      CALL WRR6(DEF,6,'DEFI')
C      CALL WRR6(TAU,6,'TAUI')
C      CALL WRR6(DEFP,6,'DEPI')
C      ENDIF
C
C???????????????????????????????????????????????????
C
      IF(IST.EQ.1) call wrr6(def1,6,'def1')
      IF(IST.EQ.1) call wrr6(tau1,6,'tau1')
      IF(IATYP.GE.4) THEN
C         DEFN=DEF1(3)-DEFT(3)+DEFZ
         DEFN=DEF1(3)+DEFP(3)-DEFT(3)-DEFPT(3)
         THI1(1)=THI1(2)*DEXP(DEFN)
         IF(THI1(1).LT.0.1D0*THIC0) THEN
            WRITE(3,*) 'NLM,NGAUSX,NGAUSY,NGAUSZ,THI1',
     +                  NLM,NGAUSX,NGAUSY,NGAUSZ,THI1(1),THI1(2)
            THI1(1)=0.1D0*THIC0
         ENDIF
         IF(IST.EQ.1) WRITE(3,*) 'DEFN,THI1',DEFN,THI1(1),THI1(2)
      ENDIF
C
      IF(NAPON.EQ.1.AND.IATYP.GE.4) THEN
         IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN
C           PRIRASTAJ PLASTICNE U GLOBALNOM DEKARTOVOM
            CALL GLAVN(DDEFPS)
            CALL GLAPR3(DDEFPS,QP)
            CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
            IF(IST.EQ.1) CALL WRR3(PRINC,3,'PP  ') 
            IF(IST.EQ.1) CALL WRR3(QP,9,'QP  ') 
            CALL CLEAR(DEF,6)
            CALL JEDNA1(DEF,PRINC,3)
C            CALL JEDNA1(DEF,DDEFPS,6)
            RETURN
         ENDIF
C
CS       TRANSFORMACIJA DEFORMACIJA NA GLOBALNI SISTEM (SMIC.INZ.)
CE       TRANSFORM STRAIN TO GLOBAL COORDS.
C        ed=QsT*el
         CALL CLEAR(DEF,6) 
         CALL MNOZI2(DEF,TSGN,SHET,6,6)
         IF(IST.EQ.1) CALL WRR3(DEF,6,'DEFD')
C         call wrr(def,6,'def ')
C         call wrr6(tsgn,36,'tss ')
C
C        KORIGOVANJE NORMIRANOG ELAST. DEF. TENZORA Be NA KRAJU KORAKA
C
         IF(IATYP.EQ.4) THEN
            DO 300 I=1,3
  300       DEF(I)=2.D0*DEF(I)+1.D0
         ELSEIF(IATYP.EQ.5) THEN
            NESAD=0
            IF(NESAD.EQ.1) THEN
               DO 301 I=4,6
  301          DEF(I)=0.5D0*DEF(I)
               IF(IST.EQ.1) CALL WRR3(DEF,6,'DEFD')
               CALL GLAVN(DEF)
               CALL GLAPR3(DEF,QP)
               CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
               P1=PRINC(1)
               P2=PRINC(2)
               P3=PRINC(3)
               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
               IF(IST.EQ.1) CALL WRR3(POM,9,'QPN ') 
            ENDIF
            IF(NESAD.EQ.0) THEN
C               WRITE(3,*) 'NLM,NGAUSX,NGAUSY,KOR,ITER',
C     +                     NLM,NGAUSX,NGAUSY,KOR,ITER
C               CALL WRR3(PRINC,3,'PPNS') 
C               CALL WRR3(QP,9,'QPNS') 
CS             TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM 
CE             TRANSFORM STRAIN TO GLOBAL COORDS.
C              Sd=QeT*Sl
               CALL CLEAR(DEF,6) 
               CALL MNOZI2(DEF,TSGD,TAU,6,6)
               IF(IST.EQ.1) CALL WRR3(DEF,6,'STRG')
               CALL GLAVN(DEF)
               CALL GLAPR3(DEF,QP)
               CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
               IF(IST.EQ.1) CALL WRR3(QP,9,'QPN ') 
C               CALL WRR3(PRINC,3,'PPN ') 
C               CALL WRR3(QP,9,'QPN ') 
CE             ELASTIC LOGARITHMIC STRAIN IN PRINCIPAL DIRECTIONS
               P1=(PRINC(1)-V*(PRINC(2)+PRINC(3)))/E
               P2=(PRINC(2)-V*(PRINC(1)+PRINC(3)))/E
               P3=(PRINC(3)-V*(PRINC(2)+PRINC(1)))/E
               IF(IST.EQ.1) WRITE(3,*) 'E,V',E,V
            ENDIF
            IF(IST.EQ.1) WRITE(3,*) 'LOGDEF-P1,P2,P3',P1,P2,P3
CE          NEW LEFT CAUCHY-GREEN DEFORMATION TENSOR Be
            DEF(1)=DEXP(2.D0*P1)
            DEF(2)=DEXP(2.D0*P2)
            DEF(3)=DEXP(2.D0*P3)
            CALL DIJADS(QP,DEF)
            IF(IST.EQ.1) CALL WRR3(DEF,6,'BENG')
         ENDIF
      ENDIF
C
C???????????????????????????????????????????????????
C
      RETURN
  331 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI816')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI816'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI816 ',
     &        '( PSEUDO-TIME )')
 2015 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI816 ',
     &        '( RADIJUS POVRSI TECENJA )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI816')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI816'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI816 ',
     &        '( PSEUDO-TIME )')
 6015 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI816 ',
     &        '( THE RADIUS OF YIELD SURFACE )')
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE EPC816(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &                  AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQ,
     &                  ALFAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE EPC ( ETP )
CE     ELASTO-PLASTIC-CREEP  EPC MATRIX
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      DIMENSION CP(6,6),AHET(*),EAHET(6),GLITL(*),ALFAT(*)
C
      DVA   =2.D0
      TR    =1.D0/3.
      A2CDT =A2C*DT
      IF (TEF.LT.1.D-8) TEF=1.D-8
C
      IF ((IPL1.NE.1).OR.(ITER.EQ.0)) THEN
        AA=1./AGAMA
        DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
C
        DO 25 I=1,6
          DO 20 J=I,6
            EIJ=-DP*DEFDS(I)*DEFDS(J)
            IF (J.LT.4) THEN
              ETP(I,J)=EIJ
            ELSE
              CP(I,J)=EIJ
            END IF
   20     CONTINUE
          IF (I.LT.4) THEN
            ETP(I,I)=ETP(I,I)+AA
          ELSE
            CP(I,I)=CP(I,I)+.5*AA
          END IF
   25   CONTINUE
      ELSE
C
        EMEAL=(EM**AN)*EALFA
        CHETP=CHET*AN1/DEFQP
        AMGAM=AGAMA*(EMEAL+1.5*CHETP*DDEFQP)+1.5*AGPOT
        ASPOT=1.5/GNORM/AMGAM
        BGPOT=A2CDT*(TEQ+1.5*CHET*DDEFQP+1.5*TDOTAB(GLITL,ALFAT)/GNORM)
     &              /AMGAM
        DLFAK=(1.5-DLAM*EMEAL)/TEQ
        ASIG =ASPOT*DLFAK
        BGAM =BGPOT*DLFAK
        CHDL1=1.+CHET*DLAM
        BIFAK=CHDL1*A2CDT
        DIFAK=CHETP*DLAM
        AA   =CHDL1/AGLHET
        ATEF2=AGLHET*TEF*TEF
        DO 100 I=1,6
  100     EAHET(I)=DEFDS(I)+DLAM*AHET(I)
C
        ASIGMA=DVA*BIFAK*ATEF2
        BSIGMA=DVA*AGPOT*ATEF2-3.*TDOTAB(EAHET,AHET)
        CSIGMA=DLAM*(DVA*AGAMA*ATEF2-3.*TDOTAB(EAHET,DEFDS))
        DSIGMA=ASIGMA-BGAM*BSIGMA-BGPOT*CHETP*CSIGMA+2*AGLHET*AGLHET*TEF
        IF (DSIGMA.LT.1.D-8) DSIGMA=1.D-8
        HJFAK =ASIG*BSIGMA+ASPOT*CHETP*CSIGMA
C
        DO 110 I=1,6
          BI=BIFAK*TAUD(I)
          CI=AGPOT*TAUD(I)-AHET(I)
          DI=DIFAK*(AGAMA*TAUD(I)-DEFDS(I))
          DO 120 J=I,6
            HJ=3*CHDL1*EAHET(J)-HJFAK*GLITL(J)
            SJOT=HJ/DSIGMA
            EIJ=(-BI*SJOT-CI*(ASIG*GLITL(J)-BGAM*SJOT)
     &           -DI*(ASPOT*GLITL(J)-BGPOT*SJOT))/AGLHET
            IF (J.LT.4) THEN
              ETP(I,J)=EIJ
            ELSE
              CP(I,J)=EIJ
            END IF
  120     CONTINUE
          IF (I.LT.4) THEN
            ETP(I,I)=ETP(I,I)+AA
          ELSE
            CP(I,I)=CP(I,I)+.5*AA
          END IF
  110   CONTINUE
      END IF
C
      CP(1,1)=TR*(DVA*ETP(1,1)-ETP(1,2)-ETP(1,3)+CM)
      CP(1,2)=TR*(DVA*ETP(1,2)-ETP(1,1)-ETP(1,3)+CM)
      CP(1,3)=TR*(DVA*ETP(1,3)-ETP(1,1)-ETP(1,2)+CM)
      CP(2,2)=TR*(DVA*ETP(2,2)-ETP(1,2)-ETP(2,3)+CM)
      CP(2,3)=TR*(DVA*ETP(2,3)-ETP(1,2)-ETP(2,2)+CM)
      CP(3,3)=TR*(DVA*ETP(3,3)-ETP(1,3)-ETP(2,3)+CM)
      CP(3,1)=CP(1,3)
      CP(3,2)=CP(2,3)
C
CE  STATIC CONDESATION
C
       DO 60 I=1,6
         ETP(I,3)=0.D0
         IF (I.EQ.3) GO TO 60
         COEFR = CP(3,I)/CP(3,3)
         DO 65 J=I,6
           IF (J.NE.3) ETP(I,J)=CP(I,J)-COEFR*CP(3,J)
           ETP(J,I)=ETP(I,J)
  65     CONTINUE
  60   CONTINUE
C
      RETURN
      END

