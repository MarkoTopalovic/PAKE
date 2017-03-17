C=======================================================================
C
C   PLASTICNOST LJUSKE  -  GENERALNA ANIZOTROPIJA (04.12.1993)
C                          CEP                    (14.01.1994)
C
C=======================================================================
      SUBROUTINE D8M10(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      CALL TI810 (A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),A(LTHI1),
     1            A(LFUN),TAU,DEF,IRAC,A(LDEFT))
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI810 ( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,THI1,
     1                   FUN,TAU,DEF,IRAC,DEFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC ANISOTROPIC MATERIAL , MIXED HARDENING
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
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /UGAOV6/ TE(6,6)
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFT(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),FUN(*),
     1          THI1(*),ALFATG(6),ALFATL(6)
      DIMENSION SE(6),DEFL(6),TAUL(6),COEFE(2)
      DIMENSION CE(6,6),CPE(3,3),CM(3),E(6),GB(3),XMX(6),
     1          Y0(6),CY(6),AN(6),AN1(6)
      DIMENSION DDEFPS(6),POM(3,3),
     1          VDEF(3,3),TAU0(6),CT(3,3,3,3),FLEX(6,6)
C      EQUIVALENCE (DEFDPR(1),DEFL(1)),(DDEF(1),TAUL(1)),(DEFDS(1),SE(1))
C     & ,(DDEFP(1),XMX(1))
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
      KDIM=6
C     PRIVREMENO
C     INDIKATOR ZA NAPONE (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
      NAPGL=1
C     INDIKATOR ZA TRANSFORMACIJU NA KOSIJEVE NAPONE (0-NE,1-DA)
      NAPKO=1
C     INDIKATOR ZA PLASTICNOST (0-U DEKARTOVOM SISTEMU, 1-U LOKALNOM SISTEMU)
C     VAZNO! PLASTICNE DEFORMACIJE MORAJU U LOKALNOM SISTEMU ZBOG DEBLJINE
      LOKAL=1
C
CE  INITIAL DATA
C
      IPL =PL
      IPL1=PL1
      DVA =2.D0
      DVT =DVA/3.D0
C.. CONSTANTS
      CALL       PEL810(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,
     &                  COEFE,ETP,BETA,TE)
      ANQ1  =ANQ-1.D0
      EM1   =1.D0-EM
      EM1DVT=EM1*DVT
      IMIX  = 1
      IF(DABS(EM1).LT.EPSIL) IMIX=0
      CALL       ENTY(DUM,E,Y0,CY,AN,EM,DEFQPT,DEFQPT,KDIM)
C
      TEQY=TEQY0
      IF(IPL.EQ.0) TEQT=TEQY0
      IF(IPL.EQ.1) TEQY=TEQT
C
      IF(IRAC.EQ.2) RETURN
C
C... TRANSFORM STRAIN INTO MATERIAL AXES DIRECTION
C
C      CALL WRR(DEF,6,'DEF ')
      CALL CLEAR(DEFL,6)
      IF(DABS(BETA).GT.1.0D-10) THEN
         CALL MNOZI1(DEFL,TE,DEF,6,6)
      ELSE
         CALL JEDNA1(DEFL,DEF,6)
      ENDIF
C
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 5 I=4,6
    5 DEFPT(I)=0.5*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
C
      IF(IATYP.LT.4) THEN
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
      TAUD(1) =CPE(1,1)*DEFDS(1)+CPE(1,2)*DEFDS(2)
      TAUD(2) =CPE(2,1)*DEFDS(1)+CPE(2,2)*DEFDS(2)
      SE(1)=TAUD(1)-ALFAT(1)
      SE(2)=TAUD(2)-ALFAT(2)
      SE(3)  =-SE(1)-SE(2)
      TAUD(3)=-TAUD(1)-TAUD(2)
      DO 45 I=4,6
        DUM=2.*DEFDS(I)*CE(I,I)
        TAUD(I) =DUM
        SE(I)=DUM-ALFAT(I)
   45 CONTINUE
      IF(IST.EQ.1) call wrr6(ALFAT,6,'ALFT')
      IF(IST.EQ.1) call wrr6(defds,6,'defs')
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5D0*TDOTAN(SE,E))
      IF(IST.EQ.1) WRITE(3,*)'teq,teqy',teq,teqy
      IF((TEQ-TEQY)/TEQY.LT.1.D-5) THEN
        TEQ  =TEQY
        DEFQP=DEFQPT
        CALL JEDNA1(DEFP,DEFPT,6)
        CALL CLEAR(DDEFPS,6)
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
      AF    = 5.D0
      DDD   = 0.1D0*(TEQ-TEQT)/EP2
      FP    = TEQ - TEQT
      TAUY  = TEQT
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      DDEFQP= DDD
C
      CALL       XMXMAT(XMX,DUM,EM1DVT,CY,AN,AN1,KDIM)
      CALL       ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      CALL       INI810(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
C
      TOLD= TAUY
      DOLD= DEFQPT
      EP  = EP2
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
C
      CALL       ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      IF(IMIX.EQ.1)THEN
        CALL     XMXMAT(XMX,DEFQP,EM1DVT,CY,AN,AN1,KDIM)
      ENDIF
      CALL     INI810(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
C
      DLAM=1.5D0*DDEFQP/TAUY
      CALL       DEV810(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
C
      TEQ=DSQRT(1.5D0*TDOTAN(TAUD,E))
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
         IF(INDCEP.EQ.0)
     1   CALL  CEP810(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                DL,DXX,DYY,DLAM,TAUY,XMX,BETA,TE,EP)
         IF(IST.EQ.1.AND.INDCEP.EQ.0) CALL WRR6(ETP,36,'ETP ')
      ENDIF
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
        DDEFPS(1) =DLAM*( (E(1)+E(2))*TAUD(1)-E(1)*TAUD(2)-E(2)*TAUD(3))
        DDEFPS(2) =DLAM*(-E(1)*TAUD(1)+(E(1)+E(3))*TAUD(2)-E(3)*TAUD(3))
        DDEFPS(3) =-DDEFPS(1)-DDEFPS(2)
        DDEFPS(4) =DLAM*E(4)*TAUD(4)
        DDEFPS(5) =DLAM*E(5)*TAUD(5)
        DDEFPS(6) =DLAM*E(6)*TAUD(6)
      IF(IST.EQ.1) call wrr6(DDEFPS,6,'DEFP')
        DEFP(1)=DDEFPS(1)+DEFPT(1)
        DEFP(2)=DDEFPS(2)+DEFPT(2)
CC        DEFP(3)=DLAM*(-E(2)*TAUD(1)-E(3)*TAUD(2)+(E(2)+E(3))*TAUD(3))+
CC     &          DEFPT(3)
        DEFP(3)=DDEFPS(3)+DEFPT(3)
        DEFP(4)=DDEFPS(4)+DEFPT(4)
        DEFP(5)=DDEFPS(5)+DEFPT(5)
        DEFP(6)=DDEFPS(6)+DEFPT(6)
C
      IF(IMIX.EQ.1)THEN
        ALFA1(1)=XMX(1)*DDEFPS(1)+ALFAT(1)
        ALFA1(2)=XMX(2)*DDEFPS(2)+ALFAT(2)
        ALFA1(3)=-ALFA1(1)-ALFA1(2)
        ALFA1(4)=XMX(4)*DDEFPS(4)+ALFAT(4)
        ALFA1(5)=XMX(5)*DDEFPS(5)+ALFAT(5)
        ALFA1(6)=XMX(6)*DDEFPS(6)+ALFAT(6)
        DO 160 I=1,6
  160   TAUD(I)=TAUD(I)+ALFA1(I)
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C
C     SREDNJI NAPON
      IF(IATYP.GE.4) THEN
         TAUM=CM(1)*(DEFL(1)-DDEFPS(1))+CM(2)*(DEFL(2)-DDEFPS(2))
      ELSE
         TAUM=CM(1)*(DEFL(1)-DEFP(1))+CM(2)*(DEFL(2)-DEFP(2))
      ENDIF
      IF(IST.EQ.1) write(3,*) 'taum,cm',taum,cm(1),cm(2)
      IF(IST.EQ.1) call wrr6(ddefps,6,'dfps')
      IF(IST.EQ.1) call wrr6(defl,6,'defl')
C
      DO 200 I=1,3
  200 TAUL(I)=TAUD(I)+TAUM
      DO 205 I=4,6
      TAUL(I)=TAUD(I)
  205 DEFP(I)=2.*DEFP(I)
C
      IF(IATYP.GE.4) THEN
         CALL CLEAR(FLEX,36)
         FLEX(1,1)=1.D0/EX
         FLEX(2,1)=-VXY/EX
         FLEX(3,1)=-VZX/EZ
         FLEX(1,2)=FLEX(2,1)
         FLEX(2,2)=1.D0/EY
         FLEX(3,2)=-VYZ/EY
         FLEX(1,3)=FLEX(3,1)
         FLEX(2,3)=FLEX(3,2)
         FLEX(3,3)=1.D0/EZ
         FLEX(4,4)=1.D0/GXY
         FLEX(5,5)=1.D0/GYZ
         FLEX(6,6)=1.D0/GZX
         CALL CLEAR(DEFL,6)
         CALL MNOZI1(DEFL,FLEX,TAUL,6,6)
      ELSE
         DEFL(3)=-VZX/EZ*TAUL(1)-VYZ/EY*TAUL(2)+DEFP(3)
      ENDIF
C
C... TRANSFORM STRESS INTO GLOBAL AXES DIRECTION - OVO NIJE TACNO
C... TRANSFORM STRESS INTO LOCAL AXES DIRECTION
C
      IF(DABS(BETA).GT.1.0D-10) THEN
        CALL CLEAR(TAU,6)
        CALL MNOZI2(TAU,TE,TAUL,6,6)
        CALL CLEAR(DEF,6)
        DEFL(4)=.5D0*DEFL(4)
        DEFL(5)=.5D0*DEFL(5)
        DEFL(6)=.5D0*DEFL(6)
        CALL MNOZI2(DEF,TE,DEFL,6,6)
        DEF(4)=2.D0*DEF(4)
        DEF(5)=2.D0*DEF(5)
        DEF(6)=2.D0*DEF(6)
      ELSE
         CALL JEDNA1(DEF,DEFL,6)
         CALL JEDNA1(TAU,TAUL,6)
      ENDIF
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
      IF(IST.EQ.1) call wrr6(taul,6,'taul')
      IF(IST.EQ.1) call wrr6(tau,6,'tau ')
      IF(IST.EQ.1) call wrr6(defp,6,'defp')
      IF(IST.EQ.1) write(3,*) 'EX,VXY',EX,VXY
      IF(IST.EQ.1) call wrr6(DEFL,6,'DEFL')
      IF(IST.EQ.1) call wrr6(DEF,6,'DEF')
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
      IF(IST.EQ.1) call wrr6(def1,6,'def1')
      IF(IST.EQ.1) call wrr6(tau1,6,'tau1')
C
      IF(IATYP.GE.4) THEN
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
         CALL MNOZI2(DEF,TSGN,DEF1,6,6)
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
            NESAD=1
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
               P1=(PRINC(1)-VXY*(PRINC(2)+PRINC(3)))/EX
               P2=(PRINC(2)-VXY*(PRINC(1)+PRINC(3)))/EX
               P3=(PRINC(3)-VXY*(PRINC(2)+PRINC(1)))/EX
               IF(IST.EQ.1) WRITE(3,*) 'E,ANI',EX,VXY
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
C=======================================================================
      SUBROUTINE PEL810(FUN,CE,TEQY0,CYQ,ANQ,Y0,CY,AN,AN1,CPE,CM,EM,
     &                  COEFE,ETP,BETA,TE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORM ELASTICITY CONSTANTS
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      DIMENSION FUN(29,*),CE(6,*),CPE(3,*),CM(*),COEFE(*),ETP(6,*)
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
   15 CE(I,J)=0.D0
C     MATRICA CE
      POM=EX-EY*VXY*VXY
      CE(1,1)=EX*EX/POM
      CE(2,2)=EX*EY/POM
      CE(1,2)=EX*EY*VXY/POM
      CE(4,4)=GXY
C      CE(5,5)=COEFE(1)*GYZ
C      CE(6,6)=COEFE(2)*GZX
      CE(5,5)=GYZ
      CE(6,6)=GZX
      DO 50 I=1,6
      DO 50 J=I,6
   50 CE(J,I)=CE(I,J)
C...   MATRIX  C'E
      CPE(1,1)=D13*(DVA*CE(1,1)-CE(1,2))
      CPE(1,2)=D13*(DVA*CE(1,2)-CE(2,2))
      CPE(2,1)=D13*(DVA*CE(1,2)-CE(1,1))
      CPE(2,2)=D13*(DVA*CE(2,2)-CE(1,2))
        DO 60 I=1,3
        CPE(I,3)=ZER
   60   CPE(3,I)=ZER
C...   VECTOR   CM
      CM(1)=D13*(CE(1,1)+CE(1,2))
      CM(2)=D13*(CE(1,2)+CE(2,2))
      CM(3)=ZER
C
      IF(DABS(BETA).GT.1.D-6) THEN
        CALL TRAETP(CE,ETP,TE)
      ELSE
        CALL JEDNA1(ETP,CE,36)
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE DEV810(TAUD,E,SE,GB,CXX,CYY,DXX,DYY,CED,DLAM,DL,XMX)
C************  ISTO KAO ZA 3D
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
      SUBROUTINE INI810(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,XMX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SE(*),E(*),CE(6,*),CPE(3,*),GB(*),XMX(*)
C
      DVA=2.D0
      XX1 = E(1)+DVA*E(2)
      XX2 = E(1)-E(3)
      CXX =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2)
      CYX = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2))
      XX1 = E(1)-E(2)
      XX2 = E(1)+DVA*E(3)
      CYY = -XX1*CPE(2,1) + XX2*(CPE(2,2)+XMX(2))
      CXY =  XX1*(CPE(1,1)+XMX(1)) - XX2*CPE(1,2)
      CED = CXX*CYY - CXY*CYX
      DXX = CYY*SE(1)+CXY*SE(2)
      DYY = CYX*SE(1)+CXX*SE(2)
      GB(1)=CE(4,4)*DVA*E(4)
      GB(2)=CE(5,5)*DVA*E(5)
      GB(3)=CE(6,6)*DVA*E(6)
      RETURN
      END
C======================================================================
      SUBROUTINE CEP810(TAUD,E,CE,CPE,CM,GB,CED,CXX,CYY,CXY,CYX,
     &                  DL,DXX,DYY,DLAM,TAUY,XMX,BETA,TE,EP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION E(*),CE(6,*),CPE(3,*),TAUD(*),CM(*),GB(*),XMX(*)
      DIMENSION P(3,3),Q(3),W(6),DDL(6),T(6),R(6),CMB(6),XEN(6,6)
      EQUIVALENCE (P(1,1),XEN(1,1))
C
      ZERO  =0.D0
      ONE   =1.D0
      DVA   =2.D0
      PO    =0.5D0
      DVT   =DVA/3.
C
      ELP = (1.5 - EP*DLAM)/TAUY
C
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
CC        P(3,I) = -P(1,I)-P(2,I)
        P(3,I) = ZERO
        DUM  = ONE + DLAM*GB(I)
        R(I3) = DVA*CE(I3,I3)/DUM
        T(I3) = GB(I)*TAUD(I3)/DUM
C
        W(I) = AXX*P(1,I) + AYY*P(2,I)
        W(I3)= DVA*E(I3)*TAUD(I3)*R(I3)
   25 CONTINUE
C
      W(3)=ZERO
C
      DLL = CXX + CYY + DVA*DLAM*CED
      Q(1) = (DXX - DLL*TAUD(1))/DL
      Q(2) = (DYY - DLL*TAUD(2))/DL
CC      Q(3) = -Q(1)-Q(2)
      Q(3) = ZERO
      W0 = ELP*(AXX*Q(1) + AYY*Q(2) - 
     &   DVA*(E(4)*TAUD(4)*T(4)+E(5)*TAUD(5)*T(5)+E(6)*TAUD(6)*T(6)))
     &   - DVT*EP*TAUY
      DUM  =  -ELP/W0
      DO 15 I=1,6 
   15 DDL(I) = DUM*W(I)
C
      DO 40 I=1,2
        DO 20 J=1,2
   20   ELAST(I,J)=P(I,J)+Q(I)*DDL(J)
        DO 35 J=4,6
   35   ELAST(I,J)=Q(I)*DDL(J)
   40 CONTINUE
      DO 60 I=4,6
        DO 50 J=1,6
   50   ELAST(I,J)=-T(I)*DDL(J)
        ELAST(I,I)=ELAST(I,I)+R(I)
   60 CONTINUE
C      
CC      AXX =  (E(1)+E(2))*CM(1)-E(1)*CM(2)-E(2)*CM(3)
CC      AYY = -E(1)*CM(1)+(E(1)+E(3))*CM(2)-E(3)*CM(3)
CC      AZZ = -E(2)*CM(1)-E(3)*CM(2)+(E(2)+E(3))*CM(3)
        AXX = (E(1)+DVA*E(2))*CM(1)-(E(1)-E(3))*CM(2)
        AYY = (E(2)-E(1))*CM(1)+(E(1)+DVA*E(3))*CM(2)
CC        AZZ = ZERO
      DUM  = AXX*TAUD(1)+AYY*TAUD(2)
      DO 70 I=1,2 
        PM = AXX*ELAST(1,I)+AYY*ELAST(2,I)
        QM = DUM*DDL(I)
C        WRITE(3,*)'I,PM,QM',I,PM,QM
        CMB(I) = CM(I)-DLAM*PM-QM
   70 CONTINUE
      DO 71 I=4,6 
        PM = AXX*ELAST(1,I)+AYY*ELAST(2,I)
        QM = DUM*DDL(I)
C        WRITE(3,*)'ELAST(1,I),ELAST(2,I)',ELAST(1,I),ELAST(2,I)
C        WRITE(3,*)'I,PM,QM',I,PM,QM
        CMB(I)= -DLAM*PM-QM
   71 CONTINUE
      CMB(3)=ZERO
C      DO 75 I=1,6
C      DO 75 J=I,6
C        ELAST(J,I)=ELAST(I,J)
C   75 CONTINUE
C        WRITE(3,*)'P11,P12,P21,P22,Q1,Q2',P(1,1),P(1,2),
C     &             P(2,1),P(2,2),Q(1),Q(2)
C        WRITE(3,*)'W1,W2,W0,XL1,XL2,R4,T4',W(1),W(2),W0,DDL(1),
C     &             DDL(2),R(4),T(4)
C        WRITE(3,*)'AXX,AYY,DUM',AXX,AYY,DUM
C        WRITE(3,*)'CM(1),CM(2),DLAM',CM(1),CM(2),DLAM
C        WRITE(3,*)'PM1,PM2,PM3,QM1,QM2,QM3',PM1,PM2,PM3,QM1,QM2,QM3
C      CALL WRR(CMB,6,'CMB ')
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
C      T(3)=XEN(3,1)*TAUD(1)+XEN(3,2)*TAUD(2)+XEN(3,3)*TAUD(3)
      T(3)=ZERO
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
       ETP(I,J)=DUM+T(I)*DDL(J)
   85 CONTINUE
      DO 87 I=4,6
      DO 87 J=I,6
        ETP(I,J)=PO*ETP(I,J)
   87 CONTINUE
      DO 88 I=1,6
      DO 88 J=I,6
        ELAST(I,J)=ETP(I,J)
   88 CONTINUE
C
      DO 82 I=1,2
      DO 80 J=I,2
   80   ELAST(I,J)=ELAST(I,J)+CMB(J)
      DO 81 J=4,6
   81   ELAST(I,J)=PO*(ELAST(I,J)+CMB(J))
   82 CONTINUE
C
      DO 90 I=1,6
      DO 90 J=I,6
        ELAST(J,I)=ELAST(I,J)
   90 CONTINUE
C
      IF(DABS(BETA).GT.1.D-6) THEN
        CALL TRAETP(ELAST,ETP,TE)
      ELSE
        CALL JEDNA1(ETP,ELAST,36)
      ENDIF
      RETURN
      END
