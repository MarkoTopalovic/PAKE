C=======================================================================
C
C   PLASTICNOST 2/D ELEMENT  -  MESOVITO OJACANJE    (01.09.1992)
C
C=======================================================================
      SUBROUTINE D2M6(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      LDEFT =LTAU   + 4*IDVA
      LDEFPT=LDEFT  + 4*IDVA
      LALFAT=LDEFPT + 4*IDVA
      LTEQT =LALFAT + 4*IDVA
      LDQPT =LTEQT  + 1*IDVA
      LIPL  =LDQPT  + 1*IDVA
      LTHI  =LIPL   + 1*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LALFA1=LDEFP1 + 4*IDVA
      LTEQT1=LALFA1 + 4*IDVA
      LDQPT1=LTEQT1 + 1*IDVA
      LIPL1 =LDQPT1 + 1*IDVA
      LTHI1 =LIPL1  + 1*IDVA
C
      CALL TAUI26(A(LIPL),A(LDEFPT),A(LALFAT),A(LTEQT),A(LDQPT),
     1            A(LIPL1),A(LTAU1),A(LDEFT1),A(LDEFP1),A(LALFA1),
     1            A(LTEQT1),A(LDQPT1),A(LTHI1),
     1            A(LFUN),MATE,TAU,DEF,IRAC,A(LDEFT))
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUI26( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,THI1,
     1                   FUN,MATE,TAU,DEF,IRAC,DEFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTOPLASTIC MATERIAL , MIXED HARDENING
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
      COMMON /KOREKJ/ AJOT
      COMMON /SRPSKI/ ISRPS
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /INDNAP/ NAPON
      COMMON /MATIZO/ E,V
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /TRA2DN/ TSGD(4,4),TSGN(4,4),INTGL
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
C
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFT(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(4),GLITL(4),THI1(*),
     1          ALFATG(6),ALFATL(6)
      DIMENSION FUN(2,MATE,*),DDEFPS(6),POM(3,3),
     1          VDEF(3,3),TAU0(6),CT(3,3,3,3)
      DATA ITMAX/100/,EPSIL/1.0D-10/,TOLD/1.0D-4/
C 
      IF(IRAC.EQ.2) RETURN
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C     THICKNESS IN INTEGRATED POINT
C     URADI ZA VISESLOJNE
      IF(IATYP.GE.4.AND.IETYP.NE.1.AND.IETYP.NE.2) THEN
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
CE  INITIAL DATA
C
      KSS=4
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) KSS=3
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
      itunel=0
      if(itunel.eq.1) then
      if(korz(6,mat,1).gt.0.and.korz(6,mat,1).le.kor) then
         e=e*(1.-float(kor-korz(6,mat,1)+1)/10.)
      endif
      if(korz(6,mat,2).gt.0.and.korz(6,mat,2).le.kor) then
         e=e/1.e9
      endif
c     odumiranje i zamena materijala
      if(korz(6,mat,3).gt.0.and.korz(6,mat,3).le.kor) then
         e=evg(6,mat,1)
         V=evg(6,mat,2)
         do i=1,4
            def(i)=def(i)-deft(i)
         enddo
      endif
      endif
C     PRIVREMENO
C     INDIKATOR ZA NAPONE (0-U DEKARTOVOM SISTEMU, 1-U GLAVNIM PRAVCIMA)
      NAPGL=0
C     INDIKATOR ZA TRANSFORMACIJU NA KOSIJEVE NAPONE (0-NE,1-DA)
      NAPKO=1
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
      CALL MEL26(FUN,MATE)
C
      IF(IRAC.EQ.2) RETURN
C
CE    YIELD STRESS
C
C      IBRAH=0
C      IF(IBRAH.EQ.1.AND.KOR.EQ.1) THEN
C         DEFQPT=0.000769D0
C      ENDIF
      IF(ISIMO.EQ.0) THEN
         TEQY =TEQY0+CY*(EM*DEFQPT)**AN
      ELSE
         TEQY =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT
      ENDIF
C
C... TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DEFPT(3)=0.5*DEFPT(3)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM
C
      IF(IATYP.LT.4) THEN
         IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
            EMT = CNI*(DEF(1)+DEF(2)-DEFPT(1)-DEFPT(2))
            DEFDS(1)= C1 *(DEF(1)-DEFPT(1))-CNI*(DEF(2)-DEFPT(2))
            DEFDS(2)=-CNI*(DEF(1)-DEFPT(1))+C1 *(DEF(2)-DEFPT(2))
         ELSE
            EMT = (DEF(1)+DEF(2)+DEF(4))/3.
            DEFDS(1)=DEF(1)-EMT-DEFPT(1)
            DEFDS(2)=DEF(2)-EMT-DEFPT(2)
            DEFDS(4)=DEF(4)-EMT-DEFPT(4)
         ENDIF
         DEFDS(3)=0.5D0*DEF(3)-DEFPT(3)
      ELSE
         IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
            EMT = CNI*(DEF(1)+DEF(2))
C            WRITE(3,*) 'emt',emt
            DEFDS(1)= C1 *DEF(1)-CNI*DEF(2)
            DEFDS(2)=-CNI*DEF(1)+C1 *DEF(2)
         ELSE
C            EMT = (DEF(1)+DEF(2)+DEF(4))/3.
            IF(IATYP.EQ.4) EMT=(DETG-1.D0)/3.D0
C            IF(IATYP.EQ.5) EMT=0.
            IF(IATYP.EQ.5) EMT=DLOG(DETG)/3.D0
            DEFDS(1)=DEF(1)-EMT
            DEFDS(2)=DEF(2)-EMT
            DEFDS(4)=DEF(4)-EMT
         ENDIF
         DEFDS(3)=0.5D0*DEF(3)
      ENDIF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD) , (GLITL), (SHET)
C
      IF(IATYP.GE.4.AND.IGLPR.EQ.1.AND.INTGL.EQ.0) THEN
C        Ag=Qs*Ad
         CALL CLEAR(ALFATL,4)
         CALL MNOZI1(ALFATL,TSGN,ALFAT,4,4)
      ELSE
         CALL JEDNA1(ALFATL,ALFAT,4)
      ENDIF
      DO 40 I=1,KSS
      GLITL(I)=DEFDS(I)-AE*ALFATL(I)
      TAUD(I) =DEFDS(I)*AEINV
   40 SHET(I) =TAUD(I)-ALFATL(I)
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
C PROVERI DA LI OVO VAZI I U DEKARTOVOM I U GLAVNIM PRAVCIMA
         GLITL(4)=-GLITL(1)-GLITL(2)
         TAUD(4) =-TAUD(1)-TAUD(2)
         SHET(4) =-SHET(1)-SHET(2)
         EMT0=EMT
         EMT=-TAUD(4)/CM
         IF(DABS(EMT0-EMT).GT.1.D-8) STOP 'EMT'
      ENDIF
      IF(IST.EQ.1) call wrr6(ALFAT,4,'ALFT')
      IF(IST.EQ.1) call wrr6(ALFATL,4,'ALFL')
      IF(IST.EQ.1) call wrr6(defds,4,'defs')
      IF(IST.EQ.1) call wrr6(glitl,4,'glit')
      IF(IST.EQ.1) call wrr6(taud,4,'taud')
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5D0*TENDO2(SHET))
      IF(IST.EQ.1) WRITE(3,*)'emt,teq,teqy',emt,teq,teqy
      IF(DABS(TEQY).GT.1.D-5) THEN
      IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
         TEQ  =TEQY
         DEFQP=DEFQPT
         CALL JEDNA1(DEFP,DEFPT,4)
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
      IF(DEFQP.LE.TOLD) DEFQP=TOLD
      IF(ISIMO.EQ.0) THEN
         EP2=AN*CY*DEFQP**AN1      
      ELSE
         EP2=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
C
      IF(CY.LT.1.D-6.AND.(IETYP.EQ.1.OR.IETYP.EQ.2)) THEN
         DUM=DSQRT(1.5D0*TENDO2(GLITL))
         DDEFQP=DVT*(DUM-AE*TEQY)
         DEFQP=DEFQPT+DDEFQP
C         IF(DEFQP.LE.TOLD) DEFQP=TOLD
         CC=0.D0
         EM1CC=0.D0
         GO TO 199
      ENDIF
      IB = 0
      IT = 0
      AF = 5.D0
      DDD= 0.1D0*(TEQ-TEQY)/EP2
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
         FP    = TEQ - TEQY
      ELSE
         GG =AE*TEQY0-DSQRT(1.5D0*TENDO2(GLITL))
         IF(ISIMO.EQ.0) THEN
            AA =AE*CY*(EM**AN)
            FP=-AA*DEFQPT**AN-GG
         ELSE
            FP = -AE*(CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT)-GG
         ENDIF
      ENDIF
      TAUY  = TEQT
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
         WRITE(IZLAZ,2001) NLM,NGR,NGS
         STOP
      ENDIF
C
      DEFQP=DEFQPT+DDEFQP
      IF(ISIMO.EQ.0) THEN
         IF(IST.EQ.1)WRITE(3,*)'DEFQP,DEFQPT,DDEFQP',DEFQP,DEFQPT,DDEFQP
         CC=CYAN23*DEFQP**AN1
         EM1CC =EM1*CC
      ELSE
         CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
         EM1CC =EM1*CC
      ENDIF
C
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
         TAUY=TEQY0+CY*(EM*DEFQP)**AN
         DLAM=1.5D0*DDEFQP/TAUY
         P1 =(C1+EM1CC*AE)*DLAM+AE
         P2 =CNI*DLAM
         P12=P1*P1-P2*P2
         SHET(1)=(P1*GLITL(1)+P2*GLITL(2))/P12
         SHET(2)=(P2*GLITL(1)+P1*GLITL(2))/P12
         SHET(4)=-SHET(1)-SHET(2)
         SHET(3)=1.D0/(AE+DLAM+EM1CC*DLAM*AE)*GLITL(3)
         TEQ =DSQRT(1.5D0*TENDO2(SHET))
C
         FB = TEQ-TAUY
      ELSE
        IF(ISIMO.EQ.0) THEN
          BB=1.5D0+1.5D0*EM1CC*AE
          FB=-AA*DEFQP**AN-BB*DDEFQP-GG
        ELSE
          BB=1.5D0+1.5D0*EM1CC*AE
          FB=-AE*(CY*(1.D0-DEXP(-AN*EM*DEFQP))-HS*EM*DEFQP)-BB*DDEFQP-GG
        ENDIF
      ENDIF
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
  199 IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
         DUM =1.D0+EM1CC*DLAM
         DO 160 I=1,4
  160    TAUD(I)=ALFATL(I)+DUM*SHET(I)
      ELSE
         IF(ISIMO.EQ.0) THEN
            TEQ =TEQY0+CY*(EM*DEFQP)**AN
         ELSE
            TEQ =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQP))+HS*EM*DEFQP
         ENDIF
         DLAM=1.5D0*DDEFQP/TEQ
         DUM =EM1CC*DLAM
         DUH =1.D0/(AE+DLAM+DUM*AE)
         DUM =1.D0+DUM
         DO 165 I=1,4
            SHET(I)=DUH*GLITL(I)
  165    TAUD(I)=ALFATL(I)+DUM*SHET(I)
      ENDIF
C
CE   4)  DETERMINE SOLUTION 
C
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
      DO 170 I=1,4
         DDEFPS(I)=DLAM*SHET(I)
  170 CONTINUE
      IF(IST.EQ.1) call wrr6(DDEFPS,4,'DEFP')
      IF(IATYP.GE.4.AND.IGLPR.EQ.1.AND.INTGL.EQ.0) THEN
C        Pd=QeT*Pg
         CALL CLEAR(ALFATG,4)
         CALL MNOZI2(ALFATG,TSGD,DDEFPS,4,4)
      ELSE
         CALL JEDNA1(ALFATG,DDEFPS,4)
      ENDIF
C      DEFZ=ALFATG(4)
      DO 171 I=1,4
         DEFP(I) =DEFPT(I)+ALFATG(I)
         ALFA1(I)=ALFAT(I)+EM1CC*ALFATG(I)
  171 CONTINUE
C
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
         IF(IATYP.LT.4) THEN
            EMT = CNI*(DEF(1)+DEF(2)-DEFP(1)-DEFP(2))
         ELSE
            EMT = -TAUD(4)/CM
         ENDIF
      ENDIF
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2) THEN
         IF(IST.EQ.1) call wrr6(GLITL,4,'GLI0')
         IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
            IF(IATYP.LT.4) THEN
               GLITL(1)=DEF(1)-EMT-DEFPT(1)-AE*ALFAT(1)
               GLITL(2)=DEF(2)-EMT-DEFPT(2)-AE*ALFAT(2)
               GLITL(3)=0.5D0*DEF(3)-DEFPT(3)-AE*ALFAT(3)
               GLITL(4)=-GLITL(1)-GLITL(2)
            ELSE
               GLITL(1)=DEF(1)-EMT-AE*ALFATL(1)
               GLITL(2)=DEF(2)-EMT-AE*ALFATL(2)
               GLITL(3)=0.5D0*DEF(3)-AE*ALFATL(3)
               GLITL(4)=-GLITL(1)-GLITL(2)
            ENDIF
         ENDIF
         IF(IST.EQ.1) call wrr6(GLITL,4,'GLIT')
         IF(IATYP.GE.4.AND.IGLPR.EQ.1.AND.INTGL.EQ.0) THEN
            CALL JEDNA1(ALFATG,GLITL,4)
C           Ed=QeT*Eg
            CALL CLEAR(GLITL,4)
            CALL MNOZI2(GLITL,TSGD,ALFATG,4,4)
         ENDIF
         IF(IST.EQ.1) call wrr6(GLITL,4,'GLID')
         IF(INDCEP.EQ.0)
     1   CALL CEP26(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &              HS,ISIMO)
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C
C     SREDNJI NAPON
      TAUM=CM*EMT
C
      DO 200 I=1,2
  200 TAU(I)=TAUD(I)+TAUM
      TAU(4)=TAUD(4)+TAUM
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) TAU(4)=0.D0
      TAU(3)=TAUD(3)
      DEFP(3)=2.D0*DEFP(3)
C NORMIRANA ELASTICNA DEFORMACIJA
      IF(IATYP.GE.4) THEN
         DO 210 I=1,2
  210    DEF(I)=TAUD(I)*AE+EMT
         IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
            DEF(4)=-(TAU(1)+TAU(2))*V/E
C            WRITE(3,*) 'EMT,DEF(4)',EMT,DEF(4)
         ELSE
            DEF(4)=TAUD(4)*AE+EMT
         ENDIF
         DEF(3)=TAUD(3)*AE*2.D0
      ENDIF
      IF(IST.EQ.1) call wrr(taud,4,'taud')
      IF(IST.EQ.1) call wrr(tau,4,'tau ')
      IF(IST.EQ.1) call wrr(defp,4,'defp')
      IF(IST.EQ.1) write(3,*) 'emt,e,V',emt,e,V
C
CE  UPDATE FROM PREVIOUS STEP
C
C OVO PROVERITI ZA MALE DEFORMACIJE
      CALL JEDNA1(SHET,DEF,4)
      IF(IST.EQ.1) call wrr(SHET,4,'SHET')
      IF(IGLPR.EQ.1) THEN
C        NAPONI U GLAVNIM PRAVCIMA ILI U DEKARTOVOM SISTEMU
         CALL JEDNA1(TAU0,TAU,4)
C        NAPONI U DEKARTOVOM SISTEMU
         IF(INTGL.EQ.0) THEN
C           Sd=QeT*Sg
            CALL CLEAR(TAU,6)
            CALL MNOZI2(TAU,TSGD,TAU0,4,4)
         ENDIF
C         call wrr(tau0,4,'tau0')
C         call wrr(tau,4,'tau ')
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
               DUM=TAU(3)
               TAU(3)=TAU(4)
               TAU(4)=DUM
               CALL PIOKOS(VDEF,TAU)
               DUM=TAU(3)
               TAU(3)=TAU(4)
               TAU(4)=DUM
               CALL CEPMT(ELAST,CT,0)
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
               CALL RRRRC(ELAST,CT,VDEF,1)
            ENDIF
         ENDIF
         IF(NAPGL.EQ.0) CALL JEDNA1(TAU0,TAU,4)
CS       TRANSFORMACIJA DEFORMACIJA NA GLOBALNI SISTEM (SMIC.INZ.)
CE       TRANSFORM STRAIN TO GLOBAL COORDS.
         IF(INTGL.EQ.0) THEN
C           ed=QsT*eg
            CALL CLEAR(DEF,6)
            CALL MNOZI2(DEF,TSGN,SHET,4,4)
         ENDIF
C         call wrr(def,4,'def ')
C         call wrr(tsgn,16,'tss ')
      ENDIF
C
CE    UPDATE FROM PREVIOUS STEP
C     NAPON I DEFORMACIJA ZA STAMPANJE
C
      DO 290 I=1,4
         DEF1(I)=DEF(I)
         IF(IGLPR.EQ.1) THEN
            TAU1(I)=TAU0(I)
         ELSE
            TAU1(I)=TAU(I)
         ENDIF
  290 CONTINUE
      IF(IST.EQ.1) call wrr(def1,4,'def1')
      IF(IST.EQ.1) call wrr(tau1,4,'tau1')
      IF(IATYP.GE.4.AND.(IETYP.EQ.0.OR.IETYP.EQ.3)) THEN
C         DEFN=DEF1(4)-DEFT(4)+DEFZ
         DEFN=DEF1(4)+DEFP(4)-DEFT(4)-DEFPT(4)
         THI1(1)=THI1(2)*DEXP(DEFN)
         IF(THI1(1).LT.0.1D0*THIC0) THEN
            WRITE(3,*) 'NLM,NGAUSX,NGAUSY,THI1',
     +                  NLM,NGAUSX,NGAUSY,THI1(1),THI1(2)
            THI1(1)=0.1D0*THIC0
         ENDIF
         IF(IST.EQ.1) WRITE(3,*) 'DEFN,THI1',DEFN,THI1(1),THI1(2)
      ENDIF
C
      IF(NAPON.EQ.1.AND.IATYP.GE.4) THEN
      IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN
         IF(INTGL.EQ.1) THEN
C           Pg=Qs*Pd
            CALL CLEAR(DEF,4)
            CALL MNOZI1(DEF,TSGN,DDEFPS,4,4)
         ELSE
            CALL JEDNA1(DEF,DDEFPS,4)
         ENDIF
         DUM=DEF(3)
         DEF(3)=DEF(4)
         DEF(4)=DUM
         RETURN
      ENDIF
C
C     KORIGOVANJE NORMIRANOG ELAST. DEF. TENZORA Be NA KRAJU KORAKA
C
      IF(IATYP.EQ.4) THEN
         DO 300 I=1,2
  300    DEF(I)=2.D0*DEF(I)+1.D0
         DUM=DEF(3)
         DEF(3)=2.D0*DEF(4)+1.D0
         DEF(4)=DUM
      ELSEIF(IATYP.EQ.5) THEN
         NESAD=0
         IF(NESAD.EQ.1) THEN
            DUM=DEF(4)
            DEF(4)=0.5D0*DEF(3)
            DEF(3)=DUM
            IF(IST.EQ.1) CALL WRR3(DEF,4,'DEFD')
            CALL GLAVN(DEF)
            CALL GLAPR3(DEF,QP)
            CALL GLASOR(PRINC,QP,0.D0,0.D0,1.D0,IMAX)
            P1=PRINC(1)
            P2=PRINC(2)
            P3=PRINC(3)
            IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
            IF(IST.EQ.1) CALL WRR3(POM,9,'QPN ') 
         ENDIF
         IF(NESAD.EQ.0) THEN
C            WRITE(3,*) 'NLM,NGAUSX,NGAUSY,KOR,ITER',
C     +                  NLM,NGAUSX,NGAUSY,KOR,ITER
C            CALL WRR3(PRINC,3,'PPNS') 
C            CALL WRR3(QP,9,'QPNS') 
            IF(IST.EQ.1) CALL WRR3(TAU,4,'STRG')
            DUM=TAU(3)
            TAU(3)=TAU(4)
            TAU(4)=DUM
            IF(IETYP.EQ.0.OR.IETYP.EQ.3) TAU(3)=0.D0
            CALL GLAVN(TAU)
            CALL GLAPR3(TAU,QP)
            CALL GLASOR(PRINC,QP,0.D0,0.D0,1.D0,IMAX)
            IF(IST.EQ.1) CALL WRR3(PRINC,3,'PPN ') 
            IF(IST.EQ.1) CALL WRR3(QP,9,'QPN ') 
C            CALL WRR3(PRINC,3,'PPN ') 
C            CALL WRR3(QP,9,'QPN ') 
CE          ELASTIC LOGARITHMIC STRAIN IN PRINCIPAL DIRECTIONS
            P1=(PRINC(1)-V*(PRINC(2)+PRINC(3)))/E
            P2=(PRINC(2)-V*(PRINC(1)+PRINC(3)))/E
            P3=(PRINC(3)-V*(PRINC(2)+PRINC(1)))/E
            IF(IST.EQ.1) WRITE(3,*) 'E,V',E,V
         ENDIF
         IF(IST.EQ.1) WRITE(3,*) 'LOGDEF-P1,P2,P3',P1,P2,P3
CE       NEW LEFT CAUCHY-GREEN DEFORMATION TENSOR Be
         DEF(1)=DEXP(2.D0*P1)
         DEF(2)=DEXP(2.D0*P2)
         DEF(3)=DEXP(2.D0*P3)
         CALL DIJADS(QP,DEF)
         IF(IST.EQ.1) CALL WRR3(DEF,6,'BENG')
C         IF(INTGL.EQ.1) THEN
C           Eg=Qe*Ed
C            CALL CLEAR(SHET,4)
C            CALL MNOZI1(SHET,TSGD,DEF,4,4)
C         ENDIF
C         DO 310 I=1,2
C  310    DEF(I)=DEXP(2.D0*SHET(I))
C         DEF(3)=DEXP(2.D0*SHET(4))
C         CALL DIJADS(QP,DEF)
      ENDIF
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI26')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TAUI26')
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE MEL26(FUN,MATE)
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
      CALL CLEAR(ELAST,36)
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
CS      RAVANSKO STANJE NAPONA
CE      PLANE STRESS
        ELAST(1,1)=E/(1.-V*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V
        ELAST(3,3)=ELAST(1,1)*(1.-V)/2.
      ELSE
CS      RAVANSKA DEFORMACIJA
CE      PLANE STRAIN
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
      SUBROUTINE CEP26(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &                 HS,ISIMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )  2D
CE     ELASTO-PLASTIC  CEP MATRIX        2D
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION GLITL(*),CP(4,4) 
C
      CALL CLEAR(ELAST,36)
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
      AP=(AE*(DVT*EPHETP+CHET)+1.)*DSQRT(TRIPO*TENDO2(GLITL))
      BP=TRIPO/AP/TEQ*(1.-EMEPH*DDEFQP/TEQ)
      DP=AE+(1.+AE*CHET)*DLAM
      DP=(BP-DVT*EM1*DLAM*DLAM*EPPRIM/AP)/DP/DP
      AA=(1.+CHET*DLAM)/(AE+(1.+AE*CHET)*DLAM)
C
      DO 22 I=1,4
      IF(I.NE.3)THEN
        DO 20 J=I,4
        IF(J.NE.3) CP(I,J)=-DP*GLITL(I)*GLITL(J)
   20   CONTINUE
        CP(I,I)=CP(I,I)+AA
      ENDIF
   22 CONTINUE
      ELAST(1,3)=-DP*GLITL(1)*GLITL(3)
      ELAST(2,3)=-DP*GLITL(2)*GLITL(3)
      ELAST(3,3)=0.5*AA-DP*GLITL(3)*GLITL(3)
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,4)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,4)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,4)+CM)
C
      IF(IETYP.NE.2) THEN
        ELAST(1,4)=TR*(DVA*CP(1,4)-CP(1,1)-CP(1,2)+CM)
        ELAST(2,4)=TR*(DVA*CP(2,4)-CP(1,2)-CP(2,2)+CM)
        ELAST(3,4)=-DP*GLITL(3)*GLITL(4)
        ELAST(4,4)=TR*(DVA*CP(4,4)-CP(1,4)-CP(2,4)+CM)
      ENDIF
C
CE  STATIC CONDESATION
C
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
       DO 60 I=1,3
       DO 60 J=I,3
        ELAST(I,J)=ELAST(I,J)-ELAST(I,4)*ELAST(J,4)/ELAST(4,4)
  60   CONTINUE
      ENDIF
C      
      DO 50 I=1,4
      DO 50 J=I,4
        ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      IF(IETYP.NE.1.AND.IETYP.NE.2) THEN
         DO 45 I=1,6
            ELAST(4,I)=0.D0
            ELAST(I,4)=0.D0
   45    CONTINUE
      ENDIF
C      
      RETURN
      END



