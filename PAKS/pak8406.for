C=======================================================================
C
C   PLASTICNOST LJUSKE  -  MESOVITO OJACANJE    (17.10.1992)
C
C=======================================================================
      SUBROUTINE D8M6(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TAUI86(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),A(LTHI1),
     1            A(LFUN),MATE,TAU,DEF,IRAC,A(LDEFT))
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUI86( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,THI1,
     1                   FUN,MATE,TAU,DEF,IRAC,DEFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC MATERIAL , MIXED HARDENING
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
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /INDNAP/ NAPON
      COMMON /MATIZO/ E,V
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
C
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFT(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(6),GLITL(6),COEFE(2),
     1          THI1(*),ALFATG(6),ALFATL(6)
      DIMENSION FUN(2,MATE,*),DDEFPS(6),POM(3,3),
     1          VDEF(3,3),TAU0(6),CT(3,3,3,3)
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
      IPL =PL
C
CE  E,V,TEQY0,CY,AN,EM
C
      E     =FUN(1,MAT,1)
      V   =FUN(2,MAT,1)
      TEQY0 =FUN(1,MAT,2)
      CY    =FUN(2,MAT,2)
      AN    =FUN(1,MAT,3)
      EM    =FUN(2,MAT,3)
      HS    =FUN(1,MAT,4)
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
CE    AUXILIARY CONSTANTS
C
      DVT   =2.D0/3.
      CZZ   =1.-V
      C1    =(2.-V)/3./CZZ
      CZZ   =-V/CZZ
      CNI   =1.-C1
      CYAN23=DVT*CY*AN
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-2.
      AE    =(1.+V)/E
      AEINV =1./AE
      CM    =E/(1.-2.*V)
C
      CALL MEL81(FUN,COEFE,ETP,1,THICK)
C
      IF(IRAC.EQ.2) RETURN
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
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
    5 DEFPT(I)=0.5*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      DEFDS(3)=0.D0
      IF(IATYP.LT.4) THEN
         EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2))
         DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNI*(DEF(2)-DEFPT(2))
         DEFDS(2)=-CNI*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))
         DO 10 I=4,6
   10    DEFDS(I)=0.5D0*DEF(I)-DEFPT(I)
      ELSE
         EMT = CNI*(DEF(1)+DEF(2))
         DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
         DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
         DO 11 I=4,6
   11    DEFDS(I)=0.5D0*DEF(I)
      ENDIF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD) , (GLITL), (SHET)
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
      DO 40 I=1,6
         GLITL(I)=DEFDS(I)-AE*ALFATL(I)
         TAUD(I) =DEFDS(I)*AEINV
   40    SHET(I) =TAUD(I)-ALFATL(I)
      GLITL(3)=-GLITL(1)-GLITL(2)
      TAUD(3) =-TAUD(1)-TAUD(2)
      SHET(3) =-SHET(1)-SHET(2)
      EMT0=EMT
      EMT=-TAUD(3)/CM
      IF(DABS(EMT0-EMT).GT.1.D-8) STOP 'EMT'
      IF(IST.EQ.1) call wrr6(ALFAT,6,'ALFT')
      IF(IST.EQ.1) call wrr6(ALFATL,6,'ALFL')
      IF(IST.EQ.1) call wrr6(defds,6,'defs')
      IF(IST.EQ.1) call wrr6(glitl,6,'glit')
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5D0*TENDOT(SHET))
      IF(IST.EQ.1) WRITE(3,*)'emt,teq,teqy',emt,teq,teqy
      IF(DABS(TEQY).GT.1.D-5) THEN
      IF((TEQ-TEQY)/TEQY.LT.1.D-5) THEN
         TEQ  =TEQY
         DEFQP=DEFQPT
         CALL JEDNA1(DEFP,DEFPT,6)
         CALL CLEAR(DDEFPS,6)
         GO TO 500
      ENDIF
      ENDIF
C
CE   3)   SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
C
CE       ITERATIONS
C     
      DEFQP=DEFQPT
      IF(DEFQP.LE.1.D-4) DEFQP=1.D-4
      IF(ISIMO.EQ.0) THEN
         EP2=AN*CY*DEFQP**AN1      
      ELSE
         EP2=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
C
      IB = 0
      IT = 0
      AF    = 5.D0
      DDD   = 0.1D0*(TEQ-TEQY)/EP2
      FP    = TEQ - TEQY
      TAUY  = TEQY
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      DDEFQP= DDD
  100 IT=IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
         IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
         IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
         WRITE(IZLAZ,2001) NLM,NGR,NGS,NGT
         STOP
      ENDIF
C
      DEFQP=DEFQPT+DDEFQP
      IF(ISIMO.EQ.0) THEN
c         IF(IST.EQ.1)WRITE(3,*)'DEFQP,DEFQPT,DDEFQP',DEFQP,DEFQPT,DDEFQP
         CC=CYAN23*DEFQP**AN1
         EM1CC =EM1*CC
      ELSE
         CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
         EM1CC =EM1*CC
      ENDIF
C
      TAUY=TEQY0+CY*(EM*DEFQP)**AN
      DLAM=1.5D0*DDEFQP/TAUY
      P1 =(C1+EM1CC*AE)*DLAM+AE
      P2 =CNI*DLAM
      P12=P1*P1-P2*P2
      SHET(1)=(P1*GLITL(1)+P2*GLITL(2))/P12
      SHET(2)=(P2*GLITL(1)+P1*GLITL(2))/P12
      SHET(3)=-SHET(1)-SHET(2)
      DO 30 I=4,6
   30 SHET(I)=1.D0/(AE+DLAM+EM1CC*DLAM*AE)*GLITL(I)
      TEQ =DSQRT(1.5D0*TENDOT(SHET))
C
      FB = TEQ-TAUY
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
         DUM =1.D0+EM1CC*DLAM
         DO 160 I=1,6
  160    TAUD(I)=ALFATL(I)+DUM*SHET(I)
C
CE   4)  DETERMINE SOLUTION 
C
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
      DO 170 I=1,6
         DDEFPS(I)=DLAM*SHET(I)
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
         DEFP(I) =DEFPT(I)+ALFATG(I)
         ALFA1(I)=ALFAT(I)+EM1CC*ALFATG(I)
  171 CONTINUE
C
      IF(IATYP.LT.4) THEN
         EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2))
      ELSE
         EMT = -TAUD(3)/CM
      ENDIF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2) THEN
         IF(IST.EQ.1) call wrr6(GLITL,6,'GLI0')
         IF(IATYP.LT.4) THEN
            GLITL(1)=DEF(1)-EMT-DEFPT(1)-AE*ALFAT(1)
            GLITL(2)=DEF(2)-EMT-DEFPT(2)-AE*ALFAT(2)
            GLITL(3)=-GLITL(1)-GLITL(2)
            DO 35 I=4,6
   35       GLITL(I)=0.5D0*DEF(I)-DEFPT(I)-AE*ALFAT(I)
         ELSE
            GLITL(1)=DEF(1)-EMT-AE*ALFATL(1)
            GLITL(2)=DEF(2)-EMT-AE*ALFATL(2)
            GLITL(3)=-GLITL(1)-GLITL(2)
            DO 36 I=4,6
   36       GLITL(I)=0.5D0*DEF(I)-AE*ALFATL(I)
         ENDIF
         IF(IST.EQ.1) call wrr6(GLITL,6,'GLIT')
         IF(INDCEP.EQ.0)
     1   CALL CEP86(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &              HS,ISIMO)
         IF(IST.EQ.1.AND.INDCEP.EQ.0) CALL WRR6(ETP,36,'ETP ')
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C
C     SREDNJI NAPON
      TAUM=CM*EMT
      IF(IST.EQ.1) write(3,*) 'taum,emt,cm',taum,emt,cm
C
      DO 200 I=1,2
  200 TAU(I)=TAUD(I)+TAUM
      TAU(3)=0.D0
      DO 46 I=4,6
      TAU(I)=TAUD(I)
   46 DEFP(I)=2.D0*DEFP(I)
C NORMIRANA ELASTICNA DEFORMACIJA
C     OVO ZA MALE DEFORMACIJE PROVERI DA LI TREBA
      IF(IATYP.LT.4) THEN
         DEF(3)=CZZ*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2))-DEFP(1)-DEFP(2)
      ENDIF
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
      DO 290 I=1,6
         DEF1(I)=DEF(I)
         IF(IGLPR.EQ.1) THEN
            TAU1(I)=TAU0(I)
         ELSE
            TAU1(I)=TAU(I)
         ENDIF
  290 CONTINUE
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
CE             TRANSFORM STRESS TO GLOBAL COORDS.
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
      RETURN
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI86')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TAUI86')
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE CEP86(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &                 HS,ISIMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ETP )    LJUSKA
CE     ELASTO-PLASTIC  CEP MATRIX        SHELL
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      DIMENSION GLITL(*),CP(6,6) 
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      DVT   =DVA*TR
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-DVA
      CHET  =EM1*CC
C
      EMEPH = 0.D0
      IF(ISIMO.EQ.0) THEN
         EPPRIM=AN*AN1*CY*DEFQP**AN2      
         IF(EM.GT.1.D-10) EMEPH = EM * AN*CY*(EM*DEFQP)**AN1      
      ELSE
         EPPRIM=AN*AN*CY*DEXP(-AN*DEFQP)
         IF(EM.GT.1.D-10) EMEPH = EM * (AN*CY*DEXP(-AN*DEFQP)+HS)
      ENDIF
      EPHETP=EMEPH+EM1*EPPRIM*DDEFQP
C
      AP=(AE*(DVT*EPHETP+CHET)+1.)*DSQRT(TRIPO*TENDOT(GLITL))
      BP=TRIPO/AP/TEQ*(1.-EMEPH*DDEFQP/TEQ)
      DP=AE+(1.+AE*CHET)*DLAM
      DP=(BP-DVT*EM1*DLAM*DLAM*EPPRIM/AP)/DP/DP
      AA=(1.+CHET*DLAM)/(AE+(1.+AE*CHET)*DLAM)
C
      DO 22 I=1,3
        DO 20 J=I,3
   20   CP(I,J)=-DP*GLITL(I)*GLITL(J)
        CP(I,I)=CP(I,I)+AA
   22 CONTINUE
      DO 30 I=1,3
      DO 30 J=4,6
        ETP(I,J)=-DP*GLITL(I)*GLITL(J)
   30 CONTINUE
      DO 40 I=4,6
      DO 40 J=I,6
        ETP(I,J)=-DP*GLITL(I)*GLITL(J)
   40 CONTINUE
      DO 42 I=4,6
        ETP(I,I)=ETP(I,I)+0.5*AA
   42 CONTINUE
      ETP(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ETP(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ETP(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ETP(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ETP(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ETP(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
C
CE  STATIC CONDESATION
C
      ETP(3,1)=ETP(1,3)
      ETP(3,2)=ETP(2,3)
      DO 60 I=1,6
        IF(I.EQ.3) GO TO 60
        COEFR = ETP(3,I)/ETP(3,3)
        DO 55 J=I,6
          IF (J.EQ.3) GO TO 55
          ETP(I,J)=ETP(I,J)-COEFR*ETP(3,J)
          ETP(J,I)=ETP(I,J)
  55    CONTINUE
  60  CONTINUE
C      
      DO 45 I=1,6
        ETP(3,I)=0.D0
        ETP(I,3)=0.D0
   45 CONTINUE
C      
      RETURN
      END
