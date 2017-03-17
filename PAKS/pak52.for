C=======================================================================
C
C     SUBROUTINE PLGL3M 
C                SISTEP (AE,AU)
C                INT3M  (ID,IMAT,NEL,NAPV,CORD,DEB, SKEM, F,EMS, SKEB,
C                                      DD,SK,FG,SKBA,SKAA,SKBB,TRC,CC,
C                                      KOSIC,ELK,KDEB,KIMAT,KNAP)
C                JAKOBP (R,S,X2,X3,Y3,A,B1,C,D)
C                READP  (AU)
C                MEL51
C=======================================================================
      SUBROUTINE PLGL3M
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C***********************************************************************
C GLAVNI UPRAVLJACKI PROGRAM ZA FORMIRANJE MATRICA I VEKTORA ELEMENATA *
C                       TROUGAONIH LJUSKI                              *
C***********************************************************************
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERM/ MREPER(4)
      COMMON/PL3CEL/ NPP,NGAUS,INULAZ,KDT,
     1IPP,IPMK,IPP2,MAXSK,INDPF,MATP,INDNA,NAPONG(3,5),LM(18),NDEB
      COMMON /PL3REP/ LSKEB,LSKEM,LEMB,LSIGB,LEN,LFUE
      COMMON /CDEBUG/ IDEBUG
      IF (IDEBUG.GT.0) PRINT *,'PLGL3M'
C
C   CITANJE SA DISKA ----  READP
C
      IPP   = 3 * NPP
      IPP2  = 2 * NPP
      IPMK  = IPP * (IPP+1)/2
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LAU=LMAX
      CALL READP(A(LAU))
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LAE=LMAX
      MXAE=0
C      MXAE=66*NE*IDVA
      LMAX = LAE + MXAE
      IF(LMAX.GT.MTOT) CALL ERROR(1)
C
C     FORMIRANJE MATRICE KRUTOSTI ELEMENATA I PAKOVANJE U SISTEM
C
      CALL SISTEP(A(LAE),A(LAU))
      RETURN
      END
C=======================================================================
      SUBROUTINE READP (AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
      include 'paka.inc'
      
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON/PL3CEL/ NPP,NGAUS,INULAZ,KDT,
     1IPP,IPMK,IPP2,MAXSK,INDPF,MATP,INDNA,NAPONG(3,5),LM(18),NDEB
      COMMON /PL3REP/ LSKEB,LSKEM,LEMB,LSIGB,LEN,LFUE
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
C
      DIMENSION AU(*)
      REAL AU
C....   ULAZNI PODACI IZ DATOTEKE  IELEM
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     1 NPP,INULAZ,INDNA,KDT,MATP,NE,INDPF,NDEB,
     2 LNEL,LTHID,LNMAT,LIPGC,LIPRC,LISNA,MXAU,
     3 MSLOJ,MXS,MSET,
     4 LNSLOJ,LMATSL,LDSLOJ,LBBET,LBET0,(CPP(J,1),J=1,3),IBB0
C
      LSTAZA(2)=LMAX8+1
      CALL READDD(AU(LNSLOJ),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAU+MXAU
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
C     PROVERI
      IF(IATYP.EQ.0) GO TO 10
      IF(NMODM.LE.4) GO TO 10
C     LPLAST=LMAX
C     NPROS=NE*NGAUS *NGAUS *15*IDVA
C     LPLAS1=LPLAST+NPROS
C     LMAX=LPLAS1+NPROS
C     IF(LMAX.GT.MTOT) CALL ERROR(1)
C     LSTAZA(3)=LMAX8+1
C     CALL READDD(A(LPLAST),NPROS/IDVA,IELEM,LMAX8,LDUZI)
C     LSTAZA(4)=LMAX8+1
C     CALL READDD(A(LPLAS1),NPROS/IDVA,IELEM,LMAX8,LDUZI)
C     RETURN
   10 LSIGB=LMAX
      LEN  =LSIGB+NE* 3* MXS* IDVA
      LFUE =LEN  +NE* 3* MXS* IDVA
      LEMB =LFUE +NE* 18*MXS* IDVA
      LMAX=LEMB +NE* 9* MXS* IDVA 
      NPROS=LMAX-LSIGB
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      LSTAZA(5)=LMAX8+1
      CALL READDD(A(LSIGB),NPROS/IDVA,IELEM,LMAX8,LDUZI)
      LSKEM=LMAX
      LSKEB=LSKEM+NE*21*IDVA
      LMAX =LSKEB+NE*45*IDVA
      NPROS=LMAX-LSKEM
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      LSTAZA(4)=LMAX8+1
      CALL READDD(A(LSKEM),NPROS/IDVA,IELEM,LMAX8,LDUZI)
      RETURN
      END
C=======================================================================
      SUBROUTINE SISTEP(AE,AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA(KCAL=1)
C     RACUNANJE NAPONA (KCAL=2)
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERM/ MREPER(4)
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /GREDE4/ IPRIR,INDGR4
      COMMON /PL3REP/ LSKEB,LSKEM,LEMB,LSIGB,LEN,LFUE
      COMMON /PL3CEL/ NPP,NGAUS,INULAZ,KDT,
     1IPP,IPMK,IPP2,MAXSK,INDPF,MATP,INDNA,NAPONG(3,5),LM(18),NDEB
      COMMON /UPRIRI/ LUPRI
C                     -----
      DIMENSION AE(*),AU(*)
      REAL AE,AU
C
C     REPERI U VEKTORU ELEMENATA AE
      LFUN=MREPER(1)
C
C     PAKOVANJE MATRICA ELEMENATA U MATRICU SISTEMA
C
      KORD=LCORD
      KPRIR=LRTDT
      IF(IATYP.GT.1) THEN
        KORD=LCORUL
        KPRIR=LUPRI
      ENDIF
C
      CALL INT3M(A(LID),AU(LNMAT),AU(LNEL),AU(LISNA),A(KORD),AU(LTHID)
     1,A(LSKEM),A(LFTDT),A(LSKEB),A(KPRIR),A(LSK),A(LFUE),
     2 A(LSIGB),A(LEMB),A(LEN),A(LMAXA),A(LFUN),
     3 AU(LNSLOJ),AU(LMATSL),AU(LBBET),AU(LDSLOJ),AU(LBET0))
C
C     IF(IATYP.EQ.0) GO TO 110
C     IF(NMODM.LE.4) GO TO 110
C     NPROS=NE*NGAUS *NGAUS *15*IDVA
C     LMA8=LSTAZA(4)-1
C     CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
C     RETURN
      NPROS=33*NE*MXS
      LMA8=LSTAZA(5)-1
      CALL WRITDD(A(LSIGB),NPROS,IELEM,LMA8,LDUZI)
      IF(ISKNP.EQ.2) RETURN
      NPROS=66*NE
      LMA8=LSTAZA(4)-1
      CALL WRITDD(A(LSKEM),NPROS,IELEM,LMA8,LDUZI)
      RETURN
      END
C=======================================================================
      SUBROUTINE INT3M(ID,IMAT,NEL,NAPV,CORD,DEB,SKEM,F,SKEB,DD,SK,FUE,
     1                 SIGB,EMB,EN,MAXA,ELK,
     1                 NSLOJ,MATSL,BBET,DSLOJ,BET0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C***********************************************************************
C   ********   GLAVNI PODPROGRAM ZA INTEGRALJENJE    *******************
C***********************************************************************
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERM/ MREPER(4)
      COMMON /ITERBR/ ITER
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
C
      COMMON/PL3CEL/ NPP,NGAUS,INULAZ,KDT,
     1IPP,IPMK,IPP2,MAXSK,INDPF,MATP,INDNA,NAPONG(3,5),LM(18),NDEB
      COMMON/PL3INT/DET,B(3,9),EEP(5,5)
      COMMON /PL3REP/ LSKEB,LSKEM,LEMB,LSIGB,LEN,LFUE
C
C      COMMON /OUT/ FUE1(20,18,1),SKEB1(20,45)  
C
      DIMENSION EMB(NE,9,*),EN(NE,3,*),SIGB(NE,3,*),FUE(NE,18,*)
C
      DIMENSION SK(*),CORD(NP,*),F(*),DEB(*),DD(*),SKEB(NE,*),SKEM(NE,*)
      DIMENSION ID(NP,*),IMAT(*),NEL(NE,*),ELK(2,*),NAPV(*),MAXA(*)
      DIMENSION SKE(171),XGT(7),XDL(3),YDL(3)
     1,TR(3,3),SP(9,9),SPP(18,9),TRG(9,18),S18(18,18)
     2,BM(3,6),FUM(9),FUMEM(3),BNL(6,9),SPS(9,6),UB(9),UM(6)
     3,ZD(3),AMB(9),MBI(9),SIGL(3),EME(171),EME18(18,18)
     4,XD(3),YD(3),DUZ(3),AA(3),B1(3),C(3),D(3),COEFE(2)
      DIMENSION NSLOJ(*),MATSL(MSLOJ,*),BBET(MSLOJ,*),DSLOJ(MSLOJ,*),
     1          BET0(*)
      DATA XGT/ .5D0, .5D0, 0.D0, 0.D0, .5D0, .5D0, .000000001D0/
C
      COEFE(1)=0.833333333D0
      COEFE(2)=COEFE(1)
CS.....   REPERI ZA MODEL MATERIJALA  
CE.....   MATERIAL MODEL POINTER
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
CS.....   MATRICA  CFS()  KOREKCIJA SMICANJA
CE.....   SHEAR FACTORS
CCC      IF(MSET.GT.0.AND.(NMODM.EQ.1.OR.NMODM.EQ.2)) THEN
CCC        DO 45 MST=1,MSET  
CCC        NNSL=NSLOJ(MST)
CCC   45   CALL   SMICC (DSLOJ,MATSL,CFS(1,MST),TGT,MST,NNSL,
CCC     1                NMODM,LFUN,LNTA,LTEM,MATE)
CCC      ENDIF
CS.....   MATRICA TRANSFORMISANJA KONSTITUITIVNIH RELACIJA LAMINATA
CE.....   FIBER ORIENTATION TRANSFORMATION MATRIX
      IF(MSET.EQ.0.AND.(NMODM.EQ.2.OR.NMODM.EQ.4.OR.
     1   NMODM.EQ.6.OR.NMODM.EQ.8).AND.IBB0.EQ.0)CALL STBETA(TBETA,BETA)
C
CS            P E T L J A    P O    E L E M N T I M A
CE            L O O P    O V E R    E L E M E N T S
C
      DO 11 NLM=1,NE
      Z0 = 0.D0
      DEB0=DEB(NLM)
      DEBNLM=DEB0
C............................................................VEKTOR LM
      DO 125 I=1,3
      M1=NEL(NLM,I)
      KK=6*(I-1)
        DO 125 J=1,6
  125   LM(KK+J)=ID(M1,J)
      DO 4 J=1,3
      K=J+1
      IF (J.EQ.3) K=1
      M1=NEL(NLM,J)
      M2=NEL(NLM,K)
      XD(J)=CORD(M1,1)-CORD(M2,1)
      YD(J)=CORD(M1,2)-CORD(M2,2)
      ZD(J)=CORD(M1,3)-CORD(M2,3)
   4  DUZ(J)= XD(J)*XD(J)+YD(J)*YD(J)+ZD(J)*ZD(J)
C...............................................MATRICA TRANSFORMACIJE
      TR(3,1)=-YD(1)*ZD(3)+ZD(1)*YD(3)
      TR(3,2)=-ZD(1)*XD(3)+XD(1)*ZD(3)
      TR(3,3)=-XD(1)*YD(3)+YD(1)*XD(3)
      TRI=DSQRT(TR(3,1)*TR(3,1)+TR(3,2)*TR(3,2)+TR(3,3)*TR(3,3))
      DO 103 I=1,3
  103 TR(3,I)= TR(3,I)/TRI
      TRI=DSQRT(DUZ(1))
      CA=-XD(1)/TRI
      SA=-YD(1)/TRI
      TR(1,1)=CA
      TR(1,2)=SA
      TR(1,3)=-ZD(1)/TRI
      TR(2,1)=TR(3,2)*TR(1,3)-TR(3,3)*TR(1,2)
      TR(2,2)=TR(3,3)*TR(1,1)-TR(3,1)*TR(1,3)
      TR(2,3)=TR(3,1)*TR(1,2)-TR(3,2)*TR(1,1)
C.................................................RELATIVNE KOORDINATE
      Y3=XD(3)*TR(2,1)+YD(3)*TR(2,2)+ZD(3)*TR(2,3)
      X3=XD(3)*TR(1,1)+YD(3)*TR(1,2)+ZD(3)*TR(1,3)
      X2=DSQRT(DUZ(1))
      XDL(1)=-X2
      XDL(2)=X2-X3
      XDL(3)=X3
      YDL(1)=0.D0
      YDL(2)=-Y3
      YDL(3)=Y3
      DO 44 J=1,3
      AA(J)=(-XDL(J)/DUZ(J))*6.
      B1(J)=(3.*XDL(J)*YDL(J)/4./DUZ(J))*4.
      D(J)=(-YDL(J)/DUZ(J))*6.
      C(J)=3*YDL(J)*YDL(J)/DUZ(J)
   44 CONTINUE
C
CS         P E T L J A    P O    S L O J E V I M A
CE         L O O P    O V E R    L A Y E R S
C
      MAT=IMAT(NLM)
      MST=MAT
      NNSL=1
      IF(MSET.GT.0) THEN
        NNSL=NSLOJ(MST)
CCC        COEFE(1)=CFS(1,MST)
CCC        COEFE(2)=CFS(2,MST)
        TTT=-1.D0
      ENDIF
C .......................................................................
      DO 25 MSL=1,NNSL
        DO 8 I=1,IPMK
        EME(I)=0.D0
    8   SKE(I)=0.D0
C
C
        IF(MSET.GT.0) THEN
           BETA=BBET (MSL,MST)
           IF(IBB0.EQ.1) BETA=BETA+BET0(NLM)
           IF(NMODM.EQ.2.OR.NMODM.EQ.4.OR.NMODM.EQ.6.OR.NMODM.EQ.8)
     1     CALL STBETA(TBETA,BETA)
           DDD =DSLOJ(MSL,MST)
           TTT =TTT+2.D0*DDD
           Z0  =(TTT-DDD)*0.5D0*DEB0
           MAT =MATSL(MSL,MST)
           DEBNLM=DEB0*DDD
        ENDIF   
C
      GO TO (1,2,14,14,14,14,14,14),NMODM
    1 CALL MEL51(A(LFUN),COEFE,EEP,1)
      GO TO 14
    2 CALL MEL52(A(LFUN),COEFE,EEP,TBETA,BETA,1)
C......................FORMIRANJE MATRICE BL ZA SAVIJANJE - PRVI KORAK
C...........................................PETLJA PO GAUSOVIM TACKAMA
C
   14 CONTINUE
C. ...................................................LINEARNI DEO SKE
      POVR=X2*Y3
      DO 113 I=1,3
      DO 113 J=1,6
  113 BM(I,J)=0.D0
C..........................FORMIRANJE MATRICE BM - LINEARNA MEMBRANSKA
      BM(1,1)=-Y3/POVR
      BM(1,3)=-BM(1,1)
      BM(2,2)=(X3-X2)/POVR
      BM(2,4)=-X3/POVR
      BM(2,6)=X2/POVR
      BM(3,1)=BM(2,2)
      BM(3,2)=BM(1,1)
      BM(3,3)=BM(2,4)
      BM(3,4)=BM(1,3)
      BM(3,5)=BM(2,6)
C.....................................MNOZENJE MATRICOM TRANSFORMACIJE
      DO 107 I=1,9
      DO 107 J=1,18
  107 TRG(I,J)=0.D0
C.............................MATRICA TRANSFORMACIJE (TB) ZA SAVIJANJE
      DO 109 J=1,3
      DO 109 K=1,3
      JK=(J-1)*6+K
      JK3=(J-1)*6+3+K
      TRG(J,JK)   =TR(3,K)
      TRG(J+3,JK3)=TR(1,K)
  109 TRG(J+6,JK3)=TR(2,K)
C...............................MATRICA KRUTOSTI SAVIJANJA
      IF(ISKNP.EQ.2) GO TO 123 
      DO 5 NG=1,3
      R=XGT(NG)
      S=XGT(NG+3)
      DO 7 I=1,3
      DO 7 J=1,9
    7 B(I,J)=0.D0
      CALL JAKOBP(R,S,X2,X3,Y3,AA,B1,C,D)
      ZPH = Z0+0.5D0*DEBNLM
      ZMH = ZPH - DEBNLM
      CONST = (ZPH*ZPH*ZPH - ZMH*ZMH*ZMH)/3.D0
      CALL MATKEP(SKE,CONST)
    5 CONTINUE
C
      DZ=0.333333333333333*DET
      DO 31 I=1,IPMK
   31 SKE(I)=SKE(I)*DZ
       DO 112 I=1,45
  112 SKEB(NLM,I)=SKE(I)
C..............................MATRICA KL - SAVIJANJA ZA OSTALE KORAKE
      KK=0
      DO 106 I=1,9
      DO 106 J=I,9
      KK=KK+1
      SP(I,J)=SKEB(NLM,KK)
  106 SP(J,I)=SP(I,J)
      DO 110 I=1,18
      DO 110 J=1,9
      SPP(I,J)=0.D0
      DO 110 K=1,9
  110 SPP(I,J)=SPP(I,J)+TRG(K,I)*SP(K,J)
      KK=0
      DO 111 I=1,18
      DO 111 J=I,18
      KK=KK+1
      SKE(KK)=0.D0
      DO 111 K=1,9
  111 SKE(KK)=SKE(KK)+SPP(I,K)*TRG(K,J)
C
C.....................VEKTOR UNUTRASNJIH SILA ZBOG SAVIJANJA PRIRASTAJ
C***       GO TO 119
C................................VEKTOR LOKALNOG POMERANJA - SAVIJANJE
  123 IF(ISKNP.EQ.1) GO TO 133
      DO 126 I=1,3
      IG=NEL(NLM,I)
      JM=(I-1)*6
      DO 126 J=1,6
      IP = ID(IG,J)
      IF (IP.EQ.0) GO TO 20
      SPP(JM+J,1)=DD(IP)
      GO TO 126
   20 SPP(JM+J,1)=0.D0
  126 CONTINUE
      DO 128 I=1,9
      UB(I)=0.D0
      DO 128 J=1,18
  128 UB(I)=UB(I)+TRG(I,J)*SPP(J,1)
C     MAX. NAPONI U TEZISTU ELEMENTA USLED SAVIJANJA
      DO 36 I=1,3
      DO 36 J=1,9
   36 B(I,J)=0.D0
      R=1.D0/3.D0
      S=R
      CALL JAKOBP(R,S,X2,X3,Y3,AA,B1,C,D)
      DO 37 I=1,3
      SIGL(I)=0.D0
      DO 38 J=1,9
   38 SIGL(I)=SIGL(I)+B(I,J)*UB(J)
   37 SIGL(I)=0.5*DEBNLM*SIGL(I)
      DO 39 I=1,3
      DO 39 J=1,3
      SIGB(NLM,I,MSL)=SIGB(NLM,I,MSL)+EEP(I,J)*SIGL(J)
   39 CONTINUE
      KK=0
      DO 121 I=1,9
      DO 121 J=I,9
      KK=KK+1
      SP(I,J)=SKEB(NLM,KK)
  121 SP(J,I)=SP(I,J)
      IF(NGENL.GT.0) THEN
      DO 124 I=1,9
      FUM(I)=0.D0
      DO 124 J=1,9
  124 FUM(I)=FUM(I)+SP(I,J)*UB(J)
C.....................................RASPOREDJIVANJE FUM U EMB(NLM,9)
      DO 129 I=1,9
  129 EMB(NLM,I,MSL)=EMB(NLM,I,MSL)+FUM(I)
C
C... VEKTOR UNUTRASNJIH SILA PO ELEMENTU ZA TRENUTAK (T)
C
      DO 134 I=1,18
        K=LM(I)
        IF(K.EQ.0) GO TO 134
        FGLS=0.D0
         DO 135 J=1,9
  135    FGLS=FGLS+TRG(J,I)*EMB(NLM,J,MSL)
         F (K)=F (K) + FGLS
  134 CONTINUE
      ENDIF
C..LINEARNI DEO MEMBRANSKE MATRICE KRUTOSTI TRANSFORMACIJOM SPP(6,9)
C....................IZ SKEM(NLM,LOKALNO) U GLOBALNI OBLIK SP(SIM,9)
  133 DO 117 I=1,6
      DO 117 J=1,9
  117 SPP(I,J)=0.D0
      DO 118 I=1,3
      SPP(1,I)=TR(1,I)
      SPP(2,I)=TR(2,I)
      SPP(3,3+I)=TR(1,I)
      SPP(4,3+I)=TR(2,I)
      SPP(5,6+I)=TR(1,I)
  118 SPP(6,6+I)=TR(2,I)
C
       IF(ISKNP.EQ.1) GO TO 158
C........................VEKTOR UNUTRASNJIH,MEMBRANSKIH,SILA,EN(NLM,3)
      IF(NGENL.GT.0) THEN
      KK=0
      DO 160 I=1,3
      IG=NEL(NLM,I)
      DO 160 J=1,3
      KK=KK+1
      IP = ID(IG,J)
      IF (IP.EQ.0) GO TO 1004
      SP(KK,1)=DD(IP)
      GO TO 160
 1004 SP(KK,1)=0.D0
  160 CONTINUE
      DO 147 I=1,6
      UM(I)=0.D0
      DO 147K=1,9
  147 UM(I)=UM(I)+SPP(I,K)*SP(K,1)
      DO 148 I=1,3
      SPS(I,1)=0.D0
      DO 148 K=1,6
  148 SPS(I,1)=SPS(I,1)+BM(I,K)*UM(K)
C.....RASPOREDJIVANJE PRIRASTAJA FUMEM(3) U EN(NLM,3)
      DO 150 I=1,3
      FUMEM(I)=0.D0
      DO 149 K=1,3
  149 FUMEM(I)=DEBNLM*EEP(I,K)*SPS(K,1)+FUMEM(I)
  150 EN(NLM,I,MSL)=EN(NLM,I,MSL)+FUMEM(I)
C
C.........VEKTOR UNUTRASNJIH MEMBRANSKIH SILA NA ELEMENTU 
C
      POVR2=POVR/2.D0
      DO 151 I=1,6
      UM(I)=0.D0
      DO 151 J=1,3
  151 UM(I)=UM(I)+POVR2*BM(J,I)*EN(NLM,J,MSL)
C
      KK=0
      DO 154 I=1,3
      DO 154 J=1,3
        I3=(I-1)*6+J
        KK=KK+1
        JM=LM(I3)
        IF(JM.EQ.0) GO TO 154
        FGLS=0.D0
          DO 155 K=1,6
  155     FGLS=FGLS+SPP(K,KK)*UM(K)
         F (JM)=F (JM) + FGLS
  154 CONTINUE
      ENDIF
C....LINEARNI DEO MEMBRANSKE MATRICE KRUTOSTI U SKEM(NLM,LOKALNI)
  158 IF(ISKNP.EQ.2) GO TO 25 
      DO 114 I=1,6
      DO 114 J=1,3
      SP(I,J)=0.D0
      DO 114 K=1,3
  114 SP(I,J)=SP(I,J)+BM(K,I)*EEP(K,J)
      KK=0
      DO 115 I=1,6
      DO 115 J=I,6
      KK=KK+1
      SKEM(NLM,KK)=0.D0
      DO 116 K=1,3
  116 SKEM(NLM,KK)=SKEM(NLM,KK)+SP(I,K)*BM(K,J)
      SKEM(NLM,KK)=SKEM(NLM,KK)*DEBNLM*POVR/2.
  115 CONTINUE
C............................LINEARNI DEO IZ SKEM(NLM,LOK.) USP(SIM.9)
      KK=0
      DO 131 I=1,6
      DO 131 J=I,6
      KK=KK+1
      SP(I,J)=SKEM(NLM,KK)
  131 SP(J,I)=SP(I,J)
      DO 132 I=1,9
      DO 132 J=1,6
      TRG(I,J)=0.D0
      DO 132 K=1,6
  132 TRG(I,J)=TRG(I,J)+SPP(K,I)*SP(K,J)
      DO 136 I=1,9
      DO 136 J=1,9
      SP(I,J)=0.D0
      DO 136 K=1,6
  136 SP(I,J)=SP(I,J)+TRG(I,K)*SPP(K,J)
      IF(IATYP.LE.1) GO TO 140
C...........................NELINEARNI DEO MEMBRANSKE MATRICE KRUTOSTI
C....................................SLAZE SE SA LINEARNIM U SP(SIM,9)
      DO 137 I=1,6
      DO 137 J=1,9
  137 BNL(I,J)=0.D0
      BNL(1,1)=-Y3/POVR
      BNL(1,3)=-BNL(1,1)
      BNL(2,1)=(X3-X2)/POVR
      BNL(2,3)=-X3/POVR
      BNL(2,5)=X2/POVR
      BNL(3,2)=BNL(1,1)
      BNL(3,4)=BNL(1,3)
      BNL(4,2)=BNL(2,1)
      BNL(4,4)=BNL(2,3)
      BNL(4,6)=BNL(2,5)
      DO 138 I=1,6
      DO 138 J=1,6
  138 TRG(I,J)=0.D0
      TRG(1,1)=EN( NLM,1,MSL)
      TRG(1,3)=EN( NLM,3,MSL)
      TRG(2,2)=EN( NLM,1,MSL)
      TRG(2,4)=TRG(1,3)
      TRG(3,1)=TRG(1,3)
      TRG(3,3)=EN( NLM,2,MSL)
      TRG(4,2)=TRG(1,3)
      TRG(4,4)=TRG(3,3)
      DO 139 I=1,6
      DO 139 J=1,4
      SPS(I,J)=0.D0
      DO 139 K=1,4
  139 SPS(I,J)=SPS(I,J)+BNL(K,I)*TRG(K,J)
      DO 143 I=1,6
      DO 143 J=1,6
      TRG(I,J)=0.D0
      DO 143 K=1,4
  143 TRG(I,J)=TRG(I,J)+POVR*SPS(I,K)*BNL(K,J)/2.
      DO 144 I=1,9
      DO 144 J=1,6
      SPS(I,J)=0.D0
      DO 144 K=1,6
  144 SPS(I,J)=SPS(I,J)+SPP(K,I)*TRG(K,J)
      DO 145 I=1,9
      DO 145 J=1,9
      DO 145 K=1,6
  145 SP(I,J)=SP(I,J)+SPS(I,K)*SPP(K,J)
C .....................SLAGANJE MEMBRANSKE MATRICE KRUTOSTI U SKE
  140 DO 141 I=1,18
      DO 141 J=1,18
      EME18(I,J)=0.D0
  141 S18(I,J)=0.D0
      DO 142 I=1,3
      DO 142 J=1,3
      S18(I,J)=SP(I,J)
      S18(I,J+6)=SP(I,J+3)
      S18(I,J+12)=SP(I,J+6)
      S18(I+6,J+6)=SP(I+3,J+3)
      S18(I+6,J+12)=SP(I+3,J+6)
      S18(I+12,J+12)=SP(I+6,J+6)
  142 CONTINUE
      IF (NDIN.EQ.0) GO TO 1061
      IF (MATP.GT.1) GO TO 1062
      RO=ELK(1,3)
      GO TO 1063
 1062 I=IMAT(NLM)
      RO=ELK(I,3)
 1063 CM=DEBNLM*DET*RO/3.
      DO 1060 I=1,3
      DO 1060 J=1,3
      KK=J+(I-1)*3
      EME18(KK,KK)=CM
      EME18(9+KK,9+KK)=CM*DEBNLM*DEBNLM*(TR(1,J)*TR(1,J)+TR(2,J)*
     *TR(2,J))/4.
 1060 CONTINUE
 1061 KK=0
      DO 146 I=1,18
      DO 146 J=I,18
      KK=KK+1
      IF (NDIN.EQ.0) GO TO 987
      EME(KK)=EME(KK)+EME18(I,J)
  987 SKE(KK)=SKE(KK)+S18(I,J)
  146 CONTINUE
C....
      JL=0
      KLL=0
      KL=0
      KK=-18
      DO 170 I=1,3
      DO 170 J=1,6
      KL=KL+1
      KK=KK+20-KL
      IF (J.LE.3) GO TO 170
      JL=JL+1
      IF(SKE(KK).GT.XGT(7)) GO TO 171
      KLL=1
      MBI(JL)=KK
      AMB(JL)=0.D0
      GO TO 170
  171 MBI(JL)=KK
      AMB(JL)=SKE(KK)
  170 CONTINUE
      IF (KLL.EQ.0) GO TO 15
      KLL=0
      DO 172 I=1,9
      IF(DABS(AMB(I)).LT.1.0D-20) GO TO 172
      KLL=KLL+1
      IF(KLL.GT.1) GO TO 173
      AMIN=AMB(I)
      GO TO 172
  173 IF (AMB(I).GT.AMIN) GO TO 172
      AMIN=AMB(I)
  172 CONTINUE
      AMIN=0.0001*AMIN
      DO 174 I=1,9
      IF (AMB(I).GT.0) GO TO 174
      KK=MBI(I)
      SKE(KK)=AMIN
  174 CONTINUE
   15 NDFE = NPP * 6
C
C....  RAZMESTANJE U MATRICU KRUTOSTI  S I S T E M A
      CALL SPAKUJ(SK, MAXA,SKE,LM,NDFE)
C
C....  KRAJ PETLJE PO SLOJEVIMA
C
   25 CONTINUE
C
   11 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE JAKOBP(R,S,X2,X3,Y3,A,B1,C,D)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PL3INT/DET,B(3,9),EEP(5,5)
      DIMENSION HYR(9),HYS(9),HXR(9),HXS(9),A(*),B1(*),C(*),D(*)
      HXR(1)=A(1)*(1.-2.*R)+(A(3)-A(1))*S
      HXR(4)=B1(1)*(1.-2.*R)-(B1(3)+B1(1))*S
      HXR(7)=                  -4.+6.*(R+S)+C(1)*(1.-2.*R)-S*(C(3)+C(1))
      HXR(2)=-A(1)*(1.-2.*R)+(A(2)+A(1))*S
      HXR(5)=B1(1)*(1.-2.*R)-(B1(1)-B1(2))*S
      HXR(8)=                       6.*R-2.+C(1)*(1.-2.*R)+S*(C(2)-C(1))
      HXR(3)=-S*(A(3)+A(2))
      HXR(6)=+S*(B1(2)-B1(3))
      HXR(9)=-S*(C(3)-C(2))
      HYR(1)=D(1)*(1.-2.*R)+S*(D(3)-D(1))
      HYR(4)=                      C(1)*(1.-2.*R)-S*(C(3)+C(1))      +1.
      HYR(7)=-B1(1)*(1.-2.*R)+S*(B1(3)+B1(1))
      HYR(2)=-D(1)*(1.-2.*R)+S*(D(2)+D(1))
      HYR(5)=                      C(1)*(1.-2.*R)+S*(C(2)-C(1))      -1.
      HYR(8)=-B1(1)*(1.-2.*R)-S*(B1(2)-B1(1))
      HYR(3)=-S*(D(2)+D(3))
      HYR(6)=S*(C(2)-C(3))
      HYR(9)=-S*(B1(2)-B1(3))
      HXS(1)=-A(3)*(1.-2.*S)-R*(A(1)-A(3))
      HXS(4)=B1(3)*(1.-2.*S)-R*(B1(1)+B1(3))
      HXS(7)=            -4.+6.*(R+S)+C(3)*(1.-2.*S)-R*(C(3)+C(1))
      HXS(2)=R*(A(2)+A(1))
      HXS(5)=R*(B1(2)-B1(1))
      HXS(8)=-R*(C(1)-C(2))
      HXS(3)=A(3)*(1.-2.*S)-R*(A(3)+A(2))
      HXS(6)=B1(3)*(1.-2.*S)+R*(B1(2)-B1(3))
      HXS(9)=              -2.+6.*S+C(3)*(1.-2.*S)+R*(C(2)-C(3))
      HYS(1)=-D(3)*(1.-2.*S)-R*(D(1)-D(3))
      HYS(4)=C(3)*(1.-2.*S)-R*(C(3)+C(1))      +1.
      HYS(7)=-B1(3)*(1.-2.*S)+R*(B1(3)+B1(1))
      HYS(2)=R*(D(2)+D(1))
      HYS(5)=R*(C(2)-C(1))
      HYS(8)=-R*(B1(2)-B1(1))
      HYS(3)=D(3)*(1.-2.*S)-R*(D(2)+D(3))
      HYS(6)=C(3)*(1.-2.*S)+R*(C(2)-C(3))       -1.
      HYS(9)=-B1(3)*(1.-2.*S)-R*(B1(2)-B1(3))
      POV2=X2*Y3
      DET=POV2/2.
      DO 206 I=1,9
      B(1,I)=(Y3*HXR(I))/POV2
      B(2,I)=(-X3*HYR(I)+X2*HYS(I))/POV2
  206 B(3,I)=(-X3*HXR(I)+X2*HXS(I)+Y3*HYR(I))/POV2
      RETURN
      END
C=======================================================================
      SUBROUTINE MEL51(FUN,COEFE,ETP,IOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE ELAST
CE     ELAST MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      DIMENSION FUN(2,*),ETP(5,*),COEFE(*)
C
      EE=FUN(1,MAT)
      AM=FUN(2,MAT)
C
      IF(IOPT.EQ.0) THEN
        DO 12 I=1,6
        DO 12 J=1,6
   12   ELAST(I,J)=0.D0
        ELAST(1,1)=EE/(1.D0-AM*AM)
        ELAST(2,2)=ELAST(1,1)
        ELAST(5,5)=ELAST(1,1)*0.5D0*(1.D0-AM)
        ELAST(6,6)=ELAST(5,5)
      RETURN
      ENDIF
C
      DO 13 I=1,5
      DO 13 J=1,5
   13 ETP(I,J)=0.D0
C
      EA=EE/(1.D0-AM*AM)
      E1=0.5D0*(1.D0-AM)
      ETP(1,1)=EA
      ETP(2,2)=EA
      ETP(2,1)=EA*AM
      ETP(1,2)=ETP(2,1)
      ETP(3,3)=EA*E1
      ETP(4,4)=COEFE(1)*ETP(4,4)
      ETP(5,5)=COEFE(2)*ETP(4,4)
      RETURN
      END
C======================================================================
      SUBROUTINE MEL52(FUN,COEFE,ETP,TBETA,BETA,IOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE ELAST
CE     ELAST MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      DIMENSION FUN(9,*),COEFE(*)
C
      EX=FUN(1,MAT)
      EY=FUN(2,MAT)
      VXY=FUN(4,MAT)
      GXY=FUN(7,MAT)
      GYZ=FUN(8,MAT)
      GZX=FUN(9,MAT)
C
      IF(IOPT.EQ.0) THEN
        POM=EX-EY*VXY*VXY
        ELAST(1,1)=EX*EX/POM
        ELAST(2,2)=EX*EY/POM
        ELAST(5,5)=GYZ
        ELAST(6,6)=GZX
      RETURN
      ENDIF
C
      DO 15 I=1,6
      DO 15 J=I,6
   15 ELAST(I,J)=0.D0
C 
      POM=EX-EY*VXY*VXY
      ELAST(1,1)=EX*EX/POM
      ELAST(2,2)=EX*EY/POM
      ELAST(1,2)=EX*EY*VXY/POM
      ELAST(4,4)=GXY
      ELAST(5,5)=COEFE(1)*GYZ
      ELAST(6,6)=COEFE(2)*GZX
      DO 50 I=1,6
      DO 50 J=I,6
   50 ELAST(J,I)=ELAST(I,J)
C
      IF(DABS(BETA).GT.1.D-6) THEN
CCCCCC*****        CALL TRAETP(ELAST,ETP,TBETA)
      ELSE
CCCCCC*****        CALL JEDNA1(ETP,ELAST,36)
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE MATKEP(SKE,CONST)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C***********************************************************************
C    PODPROGRAM ZA FORMIRANJE MATRICE KRUTOSTI ELEMENTA                *
C***********************************************************************
      COMMON/PL3INT/DET,B(3,9),EEP(5,5)
      DIMENSION SKE(*),BT(9,3),BP(9,3)
      KK=0
      DO 207 I=1,9
      DO 207 J=1,3
  207 BT(I,J)=B(J,I)
      DO 78 I=1,9
      DO 78 J=1,3
      BP(I,J)=0.D0
      DO 79 K=1,3
   79 BP(I,J)=BP(I,J)+BT(I,K)*EEP(K,J)
   78 CONTINUE
      DO 710 I=1,9
      DO 710 J=I,9
      KK=KK+1
      DO 712 K=1,3
  712 SKE(KK)=SKE(KK)+BP(I,K)*B(K,J)*CONST
  710 CONTINUE
      RETURN
      END
