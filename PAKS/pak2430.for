C=======================================================================
C
C  ELASTICAN MATERIJAL SA UNUTRASNJIM TRENJEM I VISKOZNOSCU
C  (OVDE SAM UBACIO I VISKOZNO TRENJE,RACA)
C   ZADNJA VERZIJA 15.07.1996 !!!!!!!!!!!!!
C
C=======================================================================
      SUBROUTINE D2M30(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE  SPACE IN WORKING VECTOR 
C
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      DIMENSION TAU(*),DEF(*)
C
      LFUN=MREPER(1)
      LEL  =LPOCG
      LP0 = LEL + IDVA
      LSH= LP0  + IDVA
      LSIG = LSH + IDVA
      LEPS = LSIG + 4*IDVA
      LL10 = LEPS + 4*IDVA
      LCYC = LL10 + 10*IDVA
C
      LEL1  =LPOC1
      LP01 = LEL1 + IDVA
      LSH1= LP01  + IDVA
      LSIG1 = LSH1 + IDVA
      LEPS1 = LSIG1 + 4*IDVA
      LL101 = LEPS1 + 4*IDVA
      LCYC1 = LL101 + 10*IDVA
C
      CALL TAU230(A(LEL),A(LP0),A(LSH),A(LSIG),A(LEPS),
     1            A(LL10),A(LCYC),A(LEL1),A(LP01),
     2            A(LSH1),A(LSIG1),A(LEPS1),
     3            A(LL101),A(LCYC1),A(LFUN),TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAU230(ELST,P0T,PSH,SIGV,EPSV,DELTAT,CYC,
     1                  ELST1,P0T1,PSH1,SIGV1,EPSV1,DELTA1,CYC1,
     2                  FUN,TAU,DEF,IRAC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE   ELASTIC MATERIAL WITH INTERNAL FRICTION
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KOREKJ/ AJOT
      DIMENSION TAU(*),DEF(*),CYC(*),CYC1(*),SIGV(*),
     1          EPSV(*),SIGV1(*),EPSV1(*)
      DIMENSION FUN(14,*),FUNT(9)
      DIMENSION A(11),AN(11),W(11),WN(11),U(11),UN(11),DELTAT(10),
     1          DELTA1(10),BROZ(9),SUMD(9),DPROV(10),DELTA(10)
C
C
CE  INITIAL DATA
C
      TOL = 1.D-8
      TOLK = 1.D-5
      TOLE = 1.D-10 
      BETA = 0.D0
      MATSV = MAT
      IAX = 1
      IF (IETYP.EQ.2) IAX = 0
C     IAX = 0 - PLANE STRAIN
C     IAX = 1 - AXIAL SYMMETRY  OR PLANE STRESS
C      
      ISH=PSH
      ISH1=PSH1
C
      DO 1 I = 1,9
    1 FUNT(I)=FUN(I,MAT)
      EE=FUN(1,MAT)
      GG=FUN(7,MAT) 
      P0=FUN(10,MAT)
      AMI=FUN(11,MAT)
      CNI=FUN(12,MAT)
      EL=FUN(13,MAT)
      AT=FUN(14,MAT)
C
C
      IF (KOR.LE.1) THEN
          P0T = P0
          DO 2 I = 1,5
          CYC1(I) = 0.D0
    2     CYC(I) = 0.D0 
          DO 4 I = 1,10
    4     DELTAT(I) = 1.D0
      ENDIF
C
      DO 7 I = 1,5
    7 CYC1(I) = CYC(I) 
      DEPST = CYC(1)
      SIGPL = CYC(2)
      SIGMN = CYC(3)
      RBARP = CYC(4)
      RBARM = CYC(5)
C
      DO 3 I = 1,10
      IF (DELTAT(I).GT.0.D0) THEN
         DELTA(I) = 1.D0
      ELSE
         DELTA(I) = -1.D0
      ENDIF   
    3 DELTAT(I) = DELTA(I)
C
      EPST = DEF(1)
      GAMAT = DEF(3)
      EPSA = DABS(EPST)
      GAMAA = DABS(GAMAT)
      SIG = SIGV(1)
      ELSTT = ELST
      TEL=AT/EL
      AMIP0 = AMI*P0T
      AKE = 0.5*EE*TEL/AMIP0
      AKE4 = 4.*AKE
      AKESQ = SQRT(1.+AKE4*EPSA) 
      ET = EE/AKESQ
      FUNT(1)=ET
      AKG = 0.5*GG/AMIP0
      AKG4 = 4.*AKG
      GTXY = GG/SQRT(1.+AKG4*GAMAA)
      GTXZ = GTXY
      FUNT(7)=GTXY
      FUNT(9)=GTXZ
      CALL ANICHK(FUNT,IZLAZ,ISRPS)
      MAT = 1
      CALL MEL02 (FUNT)
      MAT = MATSV
      IF (ISH.EQ.-1) THEN
         DO 5 I = 1,4
         DO 5 J = 1,4
    5    ELAST(I,J) = TOLK*ELAST(I,J)
      ENDIF
      P0TT = P0T
      ELSTT = ELST
      IF(IRAC.EQ.2) RETURN
C     ZA ITER=0 UZIMA TAU=SIGV I VRACA SE
      IF (ITER.EQ.0) THEN
         DO 11 I = 1,4
   11    TAU(I) = SIGV(I)
         RETURN
      ENDIF
      IT = 0
      ITMAX = 30
      P0IN = P0TT 
      EEPS = EE*DEF(1)
      AV1 = CNI*EL*EL/(100.D0*EE*AT*DT)
      BK1 = AMI*EL*EL/(100.D0*EE*AT)
      CL1 = 0.1D0*EL/EE
      DEPS = DEF(1) - EPSV(1)
      DEPSA = DABS(DEPS)
C
C     ITERACIJE PO P0TT
    6 IT = IT + 1
      DO 10 I = 1,4
      TAU(I) = 0.D0
      DO 10 J = 1,4
   10 TAU(I) = TAU(I) + ELAST(I,J)*DEF(J)
      P0TT = P0 - TAU(2)
      IZNAK = 1
      IF (DEPS.LT.0.D0) IZNAK = -1
      IF (P0TT.LT.TOL) ISH = -1
      IF  (ISH.EQ.-1) THEN
          DO 15 I = 1,4
          TAU(I) = 0.D0
          DO 15 J = 1,4
   15     ELAST(I,J) = TOLK*ELAST(I,J)
          GO TO 200 
      ENDIF
C      IF (DEPSA.LT.TOLE) GO TO 30
      IF (ISH.EQ.1) GO TO 120
C
      IPL = 0
      IMIN = 0 
      DO 25 I = 1,10
      IF (DELTAT(I).GT.0.D0) IPL = 1
      IF (DELTAT(I).LT.0.D0) IMIN = 1
   25 CONTINUE
C      IF (DEPS.GT.TOLE .AND. IMIN.EQ.1) GO TO 120
C      IF (DEPS.LT.TOLE. AND. IPL.EQ.1) GO TO 120
C
C     CELA DUZINA PROKLIZAVANJA IMA JEDAN SMER NAPONA KLIZANJA
C 
      IF (KOR.EQ.1) THEN
      ISH = 0
      AMIP0 = AMI*P0TT
      AKE = 0.5*EE*TEL/AMIP0
      AKE4 = 4.*AKE
      AKESQ = SQRT(1.+AKE4*EPSA) 
      ET = EE/AKESQ
      FUNT(1)=ET
      AKG = 0.5*GG/AMIP0
      AKG4 = 4.*AKG
      AKGSQ = SQRT(1.+AKG4*GAMAA)
      GTXY = GG/AKGSQ
      GTXZ = GTXY
      FUNT(7)=GTXY
      FUNT(9)=GTXZ
      CALL ANICHK(FUNT,IZLAZ,ISRPS)
      MAT = 1
      CALL MEL02 (FUNT)
      MAT = MATSV
      IF (EPSA.LT.TOL) GO TO 50
C     NORMALNI NAPON U PRAVCU VLAKNA
      TAU(1) = AMIP0/TEL*(-1.+AKESQ)
      IF(DEF(1).LT.0.D0) TAU(1) = - TAU(1)
C *** ?     GO TO 200
   50 IF (GAMAA.GE.TOL) THEN
         TAU(3) = AMIP0*(-1.+AKGSQ)
         IF (DEF(3).LT.0.D0) TAU(3) = - TAU(3)
      ENDIF
      GO TO 180
      ENDIF
C
C     PETLJA PROMENJLIVOG OPTERECENJA
C
  120 IF (P0TT.LT.P0T) ELSTT = ELST*P0T/P0TT
      SIGOR = SIGPL
      IF (IZNAK.LT.0) THEN
         SIGOR = SIGMN
      ENDIF 
C
      ISH = 1
  123 CONTINUE     
      AV = AV1*ELSTT*ELSTT
      BK = BK1*P0TT*ELSTT*ELSTT
      CL = CL1*ELSTT
C      A1 = 1.D0/EE*(0.5D0-EL*ELSTT)
      A1 = 1.D0/EE*(0.5D0*EL-EL*ELSTT)
C
C     BROJ SEGMENATA SA PROMENJENIM SMEROM SMICANJA
C
      DO 125 I = 1,10
  125 DELTAT(I) = DELTA(I)
C
      AN(1)=0.D0
      A(1)=0.D0
      AN(2)=A1+CL
      A(2)=(-AV*DELTAT(1)*AN(2)+CL+A1)/(1.D0-AV*DELTAT(1))
      WN(2)=-BK*DELTAT(1)/(1.D0-2.D0*AV*DELTAT(1))
      W(2)=(BK*DELTAT(1)-AV*DELTAT(1)*WN(2))/(1.D0-AV*DELTAT(1))
      UN(2)=TAU(1)*AN(2)+WN(2)
      U(2)=TAU(1)*A(2)+W(2)
C
      call clear(BROZ,9)
      call clear(SUMD,9)
c      BROZ(0)=0.D0
c      SUMD(0)=0.D0
C
      DO 126 I = 3,11
      BROZ(I-2)=BROZ(I-3)+DELTAT(I-2)*(U(I-1)-UN(I-1))
      SUMD(I-2)=SUMD(I-3)+2.D0*DELTAT(I-2)
C
      AN(I) = CL+((1.D0-AV*DELTAT(I-1))*AN(I-1)-AV*DELTAT(I-1)*A(I-1))
     1 /(1.D0-2.D0*AV*DELTAT(I-1))
      A(I) = (A(I-1)-AV*DELTAT(I-1)*AN(I)+CL)/(1.D0-AV*DELTAT(I-1))
C
      WN(I)=((1.D0-AV*DELTAT(I-1))*WN(I-1)-AV*DELTAT(I-1)*W(I-1)-BK*
     1 (SUMD(I-2)+DELTAT(I-1))-2.D0*AV*BROZ(I-2))/(1.D0-2.D0*AV
     2 *DELTAT(I-1))
      W(I)=(-AV*DELTAT(I-1)*WN(I)+BK*(SUMD(I-2)+DELTAT(I-1))+2.D0*AV
     1 *BROZ(I-2)+W(I-1))/(1.D0-AV*DELTAT(I-1))
C  POMERANJA DONJEG I GORNJEG VLAKNA
      UN(I)=TAU(1)*AN(I)+WN(I)
  126 U(I)= TAU(1)*A(I)+W(I)
C
      SIG0 = 0.5D0*EL/A(11)*EPST
C      SIGNN = (SIG0 - W(11)/A(11))*100.D0
      SIGNN = (SIG0 - W(11)/A(11))
      DSIG = SIGNN - SIGOR
      DFSH1 = 0.D0
      TRENJE=0.D0
      VISC=0.D0
C
      NN = 0
  130 NN = NN + 1
      N0 = 10-NN+1
      RATIO=0.D0
      IF (IZNAK.GT.0) THEN
         IF (DELTA(N0).GT.0.D0) THEN
            RBAR = RATIO
         ELSE
            RBAR = 1. - RATIO 
            DELTAT(N0) = -DELTA(N0)
         ENDIF
      ELSE
         IF (DELTA(N0).LT.0.D0) THEN
            RBAR = RATIO
         ELSE
            RBAR = 1. - RATIO
           DELTAT(N0) = -DELTA(N0)
         ENDIF
      ENDIF 
C
      RAZL = U(N0+1)-UN(N0+1)
      VISC = VISC+4.D0*AT*CNI/DT*(U(N0+1)-UN(N0+1))*ELSTT*EL*0.1D0
      TRENJE = TRENJE+4.D0*AT*AMI*P0TT*ELSTT*EL*0.1D0
C
      DFSH1 = DFSH1 + 4.D0*AT*(AMI*P0TT + CNI/DT*
     1 (DABS(U(N0+1)-UN(N0+1)))) * ELSTT*EL*0.1D0
      DFSH = DFSH1/AT/AT
C
C      SIGNN = SIG0 - W(11)/A(11)
      DSIG = SIGNN - SIGOR
      DSIGA = DABS(DSIG)
      SIGNNA = DABS(SIGNN)
      DFSHR = DFSH - DSIGA
      DPROV(NN) = DABS(DFSHR)
      IF (SIGNNA.GT.TOL) DFSHR = DFSHR/SIGNNA
      IF (DABS(DFSHR).LE.TOL) GO TO 170 
      IF(DFSHR.GT.0.D0) GO TO 170
      IF (NN.EQ.10) THEN
      ISH = 0
      ELSTT = ELSTT + 0.1D0*ELSTT
      GO TO 123
      ENDIF
      GO TO 130
C
C
  170 IF (NN.EQ.1) THEN
      DELTAT(N0)=-DELTAT(N0)
      ELSE
      IF(DPROV(NN).GT.DPROV(NN-1)) THEN
      DELTAT(N0)=-DELTAT(N0)
      ENDIF
      ENDIF
C  170 IF(DPROV(NN).GT.DPROV(NN-1)) THEN
C      DELTAT(N0)=-DELTAT(N0)
C      ENDIF
C
      TAU(1) = SIGNN
C      DSIGA = DABS(SIGNN-SIG)
C      A50 = AMI50*P0TT*A11T
C      A5 = CON5/ELSTT*RBARS
C      ET = EE/(1.+2.*A50*DSIGA/(A5*A5)) 
C*** ZASADA UZETO ET=EE
      AMIP0 = AMI*P0TT
      ET = EE
      FUNT(1)=ET
      AKG = 0.5*GG/AMIP0
      AKG4 = 4.*AKG
      AKGSQ = SQRT(1.+AKG4*GAMAA)
      GTXY = GG/AKGSQ
      GTXZ = GTXY
      FUNT(7)=GTXY
      FUNT(9)=GTXZ
      CALL ANICHK(FUNT,IZLAZ,ISRPS)
      MAT = 1
      CALL MEL02 (FUNT)
      MAT = MATSV
C
  180 K = 2
      DO 60 I = 1,2
      TAU(K) = 0.D0
      DO 59 J = 1,4
   59 TAU(K) = TAU(K) + ELAST(K,J)*DEF(J)
      K = 4
   60 CONTINUE 
      P0TT = P0 - TAU(2)
      IF ((DABS(P0TT-P0IN)/P0IN).LT.TOL) GO TO 190
      IF  (IT.GT.ITMAX) THEN
          WRITE(3,2000) IT,P0TT,TAU(2),TAU(1),TAU(3)
 2000     FORMAT(' STOP U TAU230: NEMA RESENJA ZA P0T'/
     1    ' BROJ ITERACIJA =',I5/
     2    ' P0T =   ',D15.6/
     3    ' TAU(2) =',D15.6/
     4    ' TAU(3) =',D15.6)
          STOP
      ENDIF
      P0IN = P0TT
      GO TO 6
C
C     RACUNANJE ELST = ELDIM/EL  (EL-STAR)
  190 IF (ISH.EQ.1) GO TO 195
      IF (KOR.EQ.1) THEN
      IF (EPSA.GT.TOL) THEN
         ELSTT = 0.5*TEL*DABS(TAU(1))/AMIP0
      ELSE
         IF (GAMAA.GT.TOL) ELSTT = 0.5*DABS(TAU(3))/AMIP0
      ENDIF
      IF (CNI.GT.0.D0) ELSTT = 0.5D0*ELSTT
      ENDIF
C
  195 IF ((ELSTT-0.5).GT.TOL) THEN
         ISH = -1
         DO 196 I = 1,4
         TAU(I) = 0.D0 
         DO 196 J = 1,4
  196    ELAST(I,J) = TOLK*ELAST(I,J)
      ENDIF
C
C     AZURIRANJE VELICINA ZA SLEDECI KORAK 
C
  200 ISH1 = ISH
      PSH1 = ISH1
      P0T1 = P0TT
      ELST1 = ELSTT
      SIG1 = TAU(1)
      DO 208 I = 1,4
      SIGV1(I) = TAU(I)
  208 EPSV1(I) = DEF(I)
      CYC1(1) = DEPS
      IF (ISH.EQ.1) GO TO 210
C     KOREKCIJA ZBOG VRACANJA SA LABELE 120
      IF (IZNAK.GT.0) THEN
         CYC1(3) = SIG1
         DO 202 I = 1,10
  202    DELTA1(I) = 1.D0
      ELSE
         CYC1(2) = SIG1
         DO 203 I = 1,10
  203    DELTA1(I) = -1.D0
      ENDIF
           WRITE(3,*) ISH
           WRITE(3,456) (DELTA1(I),I=1,10)
  456      FORMAT(5F10.5) 
  210 IF (ISH.NE.1) RETURN  
      IF (DEPS.GT.0.D0) THEN
         CYC1(3) = SIG1
         IF (DEPST.LT.0.D0) CYC1(2) = SIG
      ELSE
         CYC1(2) = SIG1
         IF (DEPST.GT.0.D0) CYC1(3) = SIG
      ENDIF 
      DO 205 I = 1,10
  205 DELTA1(I) = DELTAT(I)    
C
           WRITE(3,*) ISH
           WRITE(3,486) (DELTA1(I),I=1,10)
  486      FORMAT(5F10.5) 
      RETURN
      END
