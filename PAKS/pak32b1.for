C=======================================================================
C
C=======================================================================
      SUBROUTINE ELTE3B(BLT,SKE,UEL,LM,NOP,NMAT,TSGE,HE,CORD,UPRI,RTDT,
     1                FTDT,TAU,ISNA,IPGC,LMEL,ALFE,HAEM,HINV,GEEK,IALFA,
     1                COR0,TEMGT,CORGT,AU,ZAPS,NPRZ,INDZS,GUSM,LA,CEGE,
     1                ESILA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       FOR STRAIN AND STRESS CALCULATION AND INTEGRATION OF ELEMENT 
CE.       STIFFNESS MATRIX AND INTERNAL FORCES
C .
CE.    V A R I A B L E S
CE.       NGAUSX- NUMBER OF INTEGRATION POINTS IN R-DIRECTION, /13-3/(A)
CE.       NGAUSY- NUMBER OF INTEGRATION POINTS IN S-DIRECTION, /13-3/(A)
CE.       NGAUSZ- NUMBER OF INTEGRATION POINTS IN T-DIRECTION, /13-3/(A)
C .
CE.    V E C T O R S
CE.    STRAIN(6)- VECTOR OF STRAIN COMPONENTS, E=B*U  
CE.        TA(6)- VECTOR OF STRESS COMPONENTS, S=C*E  
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /REPERM/ MREPER(4)
      COMMON /TEMPCV/ LTECV,ITEMP
      COMMON /restap/ irestp,lplas0,lstaz0
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /INTGRA/ INCOTX,INCOTY,INCOTZ
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /GAUSVR/ LTEMGT,LCORGT,ICORGT
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /DUPLAP/ IDVA
      COMMON /ITERBR/ ITER
      COMMON /INDKON/ IKONV
      COMMON /INDNAP/ NAPON
      COMMON /FILTER/ TOLNAP,TOLPOM
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /UGAOV6/ TE(6,6)
      COMMON /SRPSKI/ ISRPS
      COMMON /CRACKC/ CBTC(10,2),NBTC(10,3),KBTC,ICRACK,NBRCR
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /AKUSTI/ ROAKUS,IAKUS
      COMMON /RMISIC/ VMS(3,3),VMS1(3,3),RACGR(3,3)
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /CRKLI/ MRTVI,MRTVIU
      COMMON /IKOVAR/ INDKOV
      COMMON /BOBAN3/ IBOB,IC69
C
      DIMENSION AU(*)
      REAL AU
C
      DIMENSION BLT(6,*),SKE(*),UEL(*),LM(*),NOP(NE,*),NMAT(*),ISNA(*),
     1          CORD(NP,*),HE(NCVE,*),IPGC(*),TSGE(6,6,*),UPRI(*),
     1          LMEL(ND,*),RTDT(*),FTDT(*),TAU(NLD,NGS12,*),
     1          TEMGT(NGS12,*),CORGT(3,NGS12,*),ZAPS(*),NPRZ(*),
     1          GUSM(50,*),AMASC(21),TSS(6,6),TBETA(6,6),ESILA(ND,*),
     1          ALFE(LA,*),HAEM(LA,*),HINV(LA,LA,*),GEEK(LA,24,*)
      DIMENSION STRAIN(6),STRESS(6),TA(6)
      DIMENSION XG(55),WGT(55),NREF(11),XNC(15),WNC(15)
      DIMENSION XG9(9),YG9(9),ZG9(9),WG9(9)
      DIMENSION COR(21,3),CORT(21,3),CON(21,3),COR0(NE,3,*),
     1          FOT(6,6),EKOR(6,30),GERS(6,30),LJA(30),MJA(30),
     1          CEGE(LA,*),HAEML(30),XJT(3,3),GRAN(3,3),XJ0(3,3),
     1          CT(3,3,3,3),FE(63),DUEL(63)
C     BOBAN: VEKTORI R,S,T,S JEDINICNO,H-OVI
      DIMENSION RST(8,3),SJED(8),HOVI(8,4)
      DIMENSION CM(3,3),BOVI(8,7,3),GAME(8,4)
      DIMENSION BZHU(9,24),BNAD(9,24,7),BNL(9,24,7),BLT9(9,24)
      DIMENSION XG4(4),YG4(4),ZG4(4)
      DIMENSION ELAST9(9,9),STRAIN9(9,7),SIGMA9(9,7)
C
      DATA NREF/0,1,3,6,10,15,21,28,36,45,55/
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0,
     5       .171324492379170D0, .360761573048139D0, .467913934572691D0,
     6       .467913934572691D0, .360761573048139D0, .171324492379170D0,
     7       .129484966168870D0, .279705391489277D0, .381830050505119D0,
     8       .417959183673469D0, .381830050505119D0, .279705391489277D0,
     9       .129484966168870D0, .101228536290376D0, .222381034453374D0,
     9       .313706645877887D0, .362683783378362D0, .362683783378362D0,
     1       .313706645877887D0, .222381034453374D0, .101228536290376D0,
     2       .081274388361574D0, .180648160694857D0, .260610696402935D0,
     3       .312347077040003D0, .330239355001260D0, .312347077040003D0,
     4       .260610696402935D0, .180648160694857D0, .081274388361574D0,
     5       .066671344308688D0, .149451349150581D0, .219086362515982D0,
     6       .269266719309996D0, .295524224714753D0, .295524224714753D0,
     7       .269266719309996D0, .219086362515982D0, .149451349150581D0,
     8       .066671344308688D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0,
     5      -.932469514203152D0,-.661209386466265D0,-.238619186083197D0,
     6       .238619186083197D0, .661209386466265D0, .932469514203152D0,
     7      -.949107912342759D0,-.741531185599394D0,-.405845151377397D0,
     8                     0.D0, .405845151377397D0, .741531185599394D0,
     9       .949107912342759D0,-.960289856497536D0,-.796666477413627D0,
     9      -.525532409916329D0,-.183434642495650D0, .183434642495650D0,
     1       .525532409916329D0, .796666477413627D0, .960289856497536D0,
     2      -.968160239507626D0,-.836031107326636D0,-.613371432700590D0,
     3      -.324253423403809D0,               0.D0, .324253423403809D0,
     4       .613371432700590D0, .836031107326636D0, .968160239507626D0,
     5      -.973906528517172D0,-.865063366688985D0,-.679409568299024D0,
     6      -.433395394129247D0,-.148874338981631D0, .148874338981631D0,
     7       .433395394129247D0, .679409568299024D0, .865063366688985D0,
     8       .973906528517172D0/
C
      DATA WNC/            2.D0,               1.D0,               1.D0,
     1       .333333333333333D0,1.333333333333333D0, .333333333333333D0,
     2       .250000000000000D0, .750000000000000D0, .750000000000000D0,
     3       .250000000000000D0, .155555555555556D0, .711111111111111D0,
     4       .266666666666667D0, .711111111111111D0, .155555555555556D0/
C
      DATA XNC/            0.D0,              -1.D0,               1.D0,
     1                    -1.D0,               0.D0,               1.D0,
     2                    -1.D0,-.333333333333333D0, .333333333333333D0,
     3                     1.D0,              -1.D0,             -0.5D0,
     4                     0.D0,              0.5D0,               1.D0/
C
      DATA WG9/3.55555555555556D0,.555555555555556D0,.555555555555556D0,
     1         .555555555555556D0,.555555555555556D0,.555555555555556D0,
     1         .555555555555556D0,.555555555555556D0,.555555555555556D0/
      DATA XG9/            0.D0,-.774596669241483D0,-.774596669241483D0,
     1                          -.774596669241483D0,-.774596669241483D0,
     1                           .774596669241483D0, .774596669241483D0,
     1                           .774596669241483D0, .774596669241483D0/
      DATA YG9/            0.D0,-.774596669241483D0,-.774596669241483D0,
     1                           .774596669241483D0, .774596669241483D0,
     1                          -.774596669241483D0,-.774596669241483D0,
     1                           .774596669241483D0, .774596669241483D0/
      DATA ZG9/            0.D0,-.774596669241483D0, .774596669241483D0,
     1                          -.774596669241483D0, .774596669241483D0,
     1                          -.774596669241483D0, .774596669241483D0,
     1                          -.774596669241483D0, .774596669241483D0/
      DATA RST/ 1.D0,-1.D0,-1.D0, 1.D0, 1.D0,-1.D0,-1.D0, 1.D0,
     1          1.D0, 1.D0,-1.D0,-1.D0, 1.D0, 1.D0,-1.D0,-1.D0,
     1          1.D0, 1.D0, 1.D0, 1.D0,-1.D0,-1.D0,-1.D0,-1.D0/
      DATA HOVI/   1.D0,-1.D0, 1.D0,-1.D0, 1.D0,-1.D0, 1.D0,-1.D0,
     1             1.D0,-1.D0,-1.D0, 1.D0,-1.D0, 1.D0, 1.D0,-1.D0,
     1             1.D0, 1.D0,-1.D0,-1.D0,-1.D0,-1.D0, 1.D0, 1.D0,
     1             1.D0,-1.D0, 1.D0,-1.D0,-1.D0, 1.D0,-1.D0, 1.D0/
      DATA SJED/             1.D0,                1.D0,
     1                       1.D0,                1.D0,
     1                       1.D0,                1.D0,
     1                       1.D0,                1.D0/
      DATA XG4/ -.577350269189626D0,-.577350269189626D0,
     1           .577350269189626D0, .577350269189626D0/
      DATA YG4/ -.577350269189626D0, .577350269189626D0,
     1          -.577350269189626D0, .577350269189626D0/
      DATA ZG4/ -.577350269189626D0,-.577350269189626D0,
     1           .577350269189626D0, .577350269189626D0/
C
C     INDIKATOR DA LI JE BOBANOV ELEMENT I ZA DIMENZIJU MATRICE C
C      IBOB=1
C      IC69=1
C
      NWE=ND*(ND+1)/2
C
CS            P E T L J A    P O    E L E M N T I M A
CE            L O O P    O V E R    E L E M E N T S
C
C     INDIKATOR ZA INTEGRACIJU: (0-U GL.PRAVCIMA; 1-U DEKARTOVOM) 
      INTGL=0
C     INDIKATOR ZA GLAVNE PRAVCE
      IGLPR=0
      IF(IATYP.EQ.5.OR.ILEDE.EQ.1.OR.ICPM1.EQ.2) IGLPR=1
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C
      GRZAP=0.D0
      GRMAS=0.D0
      IF(IBB0.EQ.0) CALL STBETA(TBETA,BETA)
      IF(NBRCR.GT.0)THEN
        CALL CLEAR(CBTC,20)
        CALL ICLEAR(NBTC,30)
      ENDIF
      ICRACK=0
  900 IF(ICRACK.GT.0) ICRACK=-1
      KBTC=0
      DETGM=0.D0
      DETGP=0.D0
      MRTVI=0
      IF(KOR.EQ.1.AND.ITER.EQ.0) MRTVIU=0
      IF(MRTVIU.EQ.144.AND.ISKNP.NE.2) STOP 'PROGRAM STOP 144'
      DO 10 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 10
C
         NNXYZ=(NLM-1)*NGS12-1
         IPG=IPGC(NLM)
         ISN=ISNA(NLM)
         IF(ISN.GE.10)THEN
            ISN=ISN-10
         ENDIF
C
CS       VEKTOR  LM I VEKTOR LOKALNIH POMERANJA
CE       VECTOR  LM  AND  LOCALL DISPLACEMENTS
C
         CALL IJEDN1(LM,LMEL(1,NLM),ND)
C
      DO 1222 I=1,NCVE
         II=NOP(NLM,I)
         DO 1222 J=1,3
 1222    COR(I,J)=CORD(II,J)
C
      DO 242 I=1,ND
         UEL(I)=0.D0
         DUEL(I)=0.D0
         IF(LM(I).EQ.0) GO TO 242
         J=LM(I)
         if(j.lt.0) THEN
            uel(i)=condof(rtdt,a(lcmpc),a(lmpc),j)
            duel(i)=condof(upri,a(lcmpc),a(lmpc),j)
         ELSE
            UEL(I)=RTDT(J)
            DUEL(I)=UPRI(J)
         ENDIF
  242 CONTINUE
C
      IF(IATYP.GE.3) THEN
C
CS          KOORDINATE U TRENUTKU - 0
CE          COORDINATE IN TIME - 0
C
            II=0
            DO 11 I=1,NCVE
            DO 11 K=1,3
               II=II+1
               CON(I,K)=COR(I,K)-UEL(II)
   11       CONTINUE
C         WRITE(3,9997) ((CON(I,J),J=1,3),I=1,NCVE)
C 9997 FORMAT(' CON',3(1PD12.5))
C         WRITE(3,9998) ((COR(I,J),J=1,3),I=1,NCVE)
C 9998 FORMAT(' COR',3(1PD12.5))
C         WRITE(3,9999) (UEL(I),I=1,3*NCVE)
C 9999 FORMAT(' UEL',3(1PD12.5))
      ENDIF
C
      IF(IATYP.GE.4) THEN
C
CS          KOORDINATE U TRENUTKU - T
CE          COORDINATE IN TIME - T
C
         DO 1224 I=1,NCVE
         DO 1224 J=1,3
            IF(ITER.EQ.0) THEN
               COR0(NLM,J,I)=COR(I,J)
            ENDIF
            CORT(I,J)=COR0(NLM,J,I)
 1224    CONTINUE
C         WRITE(3,9996) ((CORT(I,J),J=1,3),I=1,NCVE)
C 9996 FORMAT(' CORT',3(1PD12.5))
      ENDIF
C
C
      IF(ISKNP.NE.2) CALL CLEAR(SKE,NWE)
      CALL CLEAR(FE,ND)
      CALL CLEAR(AMASC,21)
      CALL CLEAR(BLT,6*ND)
      CALL CLEAR(HE,4*NCVE)
C
      R=0.0D0
      S=0.0D0
      T=0.0D0
      IF(IALFA.GE.0) THEN
         IF(ISKNP.NE.1) THEN
C
CS          RACUNANJE IZRAZA - F * U + h
CE          CALCULATE EXPRESSION - F * U + h
C
            DO 180 I=1,LA
               HAEML(I)=0.D0
               DO 190 J=1,ND
                  LL=LM(J)
                  IF(LL.EQ.0) GO TO 190
                  IF(NGENL.GT.0) THEN
                     HAEML(I)=HAEML(I)+GEEK(I,J,NLM)*DUEL(J)
C                     HAEML(I)=HAEML(I)+GEEK(I,J,NLM)*UPRI(LL)
                  ELSE
                     HAEML(I)=HAEML(I)+GEEK(I,J,NLM)*UEL(J)
C                     HAEML(I)=HAEML(I)+GEEK(I,J,NLM)*RTDT(LL)
                  ENDIF
  190          CONTINUE
               IF(NGENL.GT.0) HAEML(I)=HAEML(I)+HAEM(I,NLM)
  180       CONTINUE
C
CS          RACUNANJE PARAMETARA - ALFA
CE          CALCULATE PARAMETERS - ALFA
C           ALFA = ALFA - (H**-1) * (F * U + h )
C
            IF(NGENL.EQ.0) CALL CLEAR(ALFE(1,NLM),LA)
            CALL INTEV1(ALFE(1,NLM),HINV(1,1,NLM),HAEML,-1.D0,LA,LA)
C
         ENDIF
C
CS       BRISANJE PROSTORA ZA MATRICE H, F, h, E
CS       CLEAR SPACE FOR MATRIX H, F, h, E
C
         IF(ISKNP.NE.2) CALL CLEAR(HINV(1,1,NLM),LA*LA)
         IF(ISKNP.NE.2) CALL CLEAR(GEEK(1,1,NLM),LA*ND)
         IF(NGENL.GT.0) CALL CLEAR(HAEM(1,NLM),LA)
         CALL CLEAR(EKOR,6*LA)
C
CS       INTERPOLACIJSKE FUNKCIJE I JAKOBIJAN U TACKI R, S, T
CE       INTERPOLATED FUNCTIONS AND JACOBIAN MATRIX IN R,S,T POINT
C
         CALL JACTE3(NOP,CORD,HE,R,S,T,0)
C
CS       INVERTOVAN JAKOBIJAN - XJ
CE       INVERSE JACOBIAN - XJ
C
         DET0=DET
C
CS       MATRICA TRANSFORMACIJE - FOT (KOVARIJANTNI - GLOBALNI DEKARTOV) 
CE       TRANSFORMATION MATRIX - FOT (COVARIANT - GLOBAL CARTESIAN)
C
         CALL TRANSE(FOT,XJ)
C
      ENDIF 
C
C
CS    P E T L J A    P O    S L O J E V I M A
CE    L O O P    O V E R    L A Y E R S
C
C
      IBTC=0
      IRAC=2
      MAT=NMAT(NLM)
      GUST=GUSM(NMODM,MAT)
         GO TO (  1,  2,999, 40,999,999,999,999,999, 40,
     1          999,999,999,999,999,999, 40, 40, 40,999,
     2           40,999,999,999,999,999,999,999,999,999,
     3          999,999,999,999,999,999,999,999,999,999,
     4          999,999,999,999,999,999,999,999,999,999,
     5          999,999,999,999,999,999,999,999,999,999,
     6          999,999,999,999,999,999,999,999,999,999,
     7          999,999,999,999,999,999,999,999,999,999,
     8          999,999,999,999,999,999,999,999,999,999,
     9          999,999,999,999,999,999,999,999,999,999),NMODM
CE    POINTER FOR MATERIAL DATA
    1 LFUN=MREPER(1)
      CALL MEL31(A(LFUN))
      GO TO 999
    2 LFUN=MREPER(1)
      CALL MEL32(A(LFUN))
CS    MATRICA TRANSFORMACIJE
      IF(IBB0.EQ.1) THEN
         IF(IATYP.GT.1) THEN
            CALL JACTE3(NOP,CORD,HE,R,S,T,0)
            CALL TRANAL(XJJ,TSS,0)
            CALL MNOZM1(TSG,TSGE(1,1,NLM),TSS,6,6,6)
         ENDIF
      ELSE
         IF(IATYP.LE.1) THEN
            CALL JEDNA1(TSG,TSGE(1,1,NLM),36)
         ELSE
            CALL JACTE3(NOP,CORD,HE,R,S,T,0)
            CALL TRANAL(XJJ,TSS,0)
            CALL MNOZM1(TSG,TBETA,TSS,6,6,6)
         ENDIF
      ENDIF
      CALL TRAETP(ELAST,ELAST,TSG)
      GO TO 999
CS    MATRICA TRANSFORMACIJE
   40 IF(IBB0.EQ.1) THEN
         IF(IATYP.GT.1) THEN
            CALL JACTE3(NOP,CORD,HE,R,S,T,0)
            CALL TRANAL(XJJ,TSS,0)
            CALL MNOZM1(TSG,TSGE(1,1,NLM),TSS,6,6,6)
         ENDIF
      ELSE
         IF(IATYP.LE.1) THEN
            CALL JEDNA1(TSG,TSGE(1,1,NLM),36)
         ELSE
            CALL JACTE3(NOP,CORD,HE,R,S,T,0)
            CALL TRANAL(XJJ,TSS,0)
            CALL MNOZM1(TSG,TBETA,TSS,6,6,6)
         ENDIF
      ENDIF
C     BOBAN: MATRICA C
 999  IF(IBOB.EQ.1) THEN
C      (5A), (4)
       CALL CMAT(RST,COR,CM)
C      (7)
       CALL NAPUNIB(CM,RST,HOVI,COR,BOVI,GAME)
      ENDIF
C
CS       PETLJA PO GAUSOVIM TACKAMA 
CE       LOOP OVER GAUSS POINTS
C
         JG=0
         JGCR=0
         ZAPRE=0.D0
         AMASA=0.D0
         SNAP=0.D0
         DO 20 NGR=1,NGAUSX
            JGR=NREF(NGAUSX)+NGR
            IF(INCOTX.EQ.0)THEN
              IF(NGAUSX.EQ.9)THEN
                R=XG9(NGR)
                S=YG9(NGR)
                T=ZG9(NGR)
                WR=WG9(NGR)
                WS=1.D0
                WT=1.D0
                NGAUSY=1
              ELSE
C     ZA ZHU-OVU INTEGRACIJU U 4 TACKE;   BOBAN
               IF(NGAUSX.EQ.4)THEN
                R=XG4(NGR)
                S=YG4(NGR)
                T=ZG4(NGR)
                WR=2.D0
                WS=1.D0
                WT=1.D0
               ELSE
                R=XG(JGR)
                WR=WGT(JGR)
               ENDIF
              ENDIF
            ELSE
               R =XNC(JGR)
               WR=WNC(JGR)
            ENDIF
C
         DO 20 NGS=1,NGAUSY
            JGR=NREF(NGAUSY)+NGS
            IF(INCOTY.EQ.0)THEN
              IF(NGAUSX.NE.9 .AND. NGAUSX.NE.4)THEN
                S=XG(JGR)
                WS=WGT(JGR)
              ENDIF
            ELSE
               S =XNC(JGR)
               WS=WNC(JGR)
            ENDIF
C
         DO 20 NGT=1,NGAUSZ
            JGR=NREF(NGAUSZ)+NGT
            IF(INCOTZ.EQ.0)THEN
              IF(NGAUSX.NE.9 .AND. NGAUSX.NE.4)THEN
                T=XG(JGR)
                WT=WGT(JGR)
              ENDIF
            ELSE
               T =XNC(JGR)
               WT=WNC(JGR)
            ENDIF
C
            JG=JG+1
       WRITE(3,*) 'NLM,NGR,NGS,NGT,R,S,T',NLM,NGR,NGS,NGT,R,S,T
C
CE          INTERPOLATION FUNCTIONS AND JACOBIAN MATRIX IN CURRENT
CE          INTEGRATION POINT (R,S,T - ARE NATURAL COORDINATES)
CS          INTERPOLACIJSKE FUNKCIJE I JAKOBIJAN U TACKI R, S, T
C
            IF(IBOB.EQ.0) THEN
             CALL JACTE3(NOP,CORD,HE,R,S,T,0)
            ELSE
C            (11), (12)
             CALL JACBOB(COR,RST,HOVI,R,S,T)
C            (16)
             CALL DOPUNIB(BOVI,GAME)
C            (6)(15)
             CALL NAPUNIHE(HE,HOVI,BOVI,SJED,RST,R,S,T)
            ENDIF
            
C            CALL WRR3(XJJ,9,'XJJ ')
C
CS          TEMPERATURA, GLOBALNE KOORDINATE GAUSOVE TACKE
CE          TEMPHERATURE AND GLOBAL COORDINATES OF GAUSS POINT
C
            IF(ICORGT.EQ.1.OR.ITERME.EQ.1)
     1      CALL JACGA3(NOP,CORD,A(LTECV),HE,CORGT(1,JG,NLM),TGT)
            IF(ITERME.EQ.1) TEMGT(JG,NLM)=TGT
C
            DETT=1.D0
            IF(IALFA.GE.0) THEN
               DETT=DET0/DET
C
CS             FORMIRANJE MATRICE - E (KOVARIJANTNA)
CE             FORM MATRIX - E (COVARIANT)
C
               IF(IALFA.EQ.1.OR.IALFA.EQ.2) THEN
                  EKOR(1,1)=R*DETT
                  EKOR(2,2)=S*DETT
                  EKOR(3,3)=T*DETT
                  EKOR(4,4)=R*DETT
                  EKOR(5,5)=S*DETT
                  EKOR(6,6)=T*DETT
                  EKOR(6,7)=R*DETT
                  EKOR(4,8)=S*DETT
                  EKOR(5,9)=T*DETT
               ENDIF
               IF(IALFA.EQ.2) THEN
                  EKOR(1,10)=R*S*DETT
                  EKOR(2,10)=-R*S*DETT
                  EKOR(4,10)=(R*R-S*S)*DETT
                  EKOR(5,11)=R*T*DETT
                  EKOR(6,11)=S*T*DETT
                  EKOR(4,12)=R*T*DETT
                  EKOR(6,12)=R*S*DETT
                  EKOR(2,13)=S*T*DETT
                  EKOR(3,13)=-S*T*DETT
                  EKOR(5,13)=(S*S-T*T)*DETT
                  EKOR(4,14)=S*T*DETT
                  EKOR(5,14)=R*S*DETT
                  EKOR(1,15)=-R*T*DETT
                  EKOR(3,15)=R*T*DETT
                  EKOR(6,15)=(T*T-R*R)*DETT
                  EKOR(1,16)=R*S*T*DETT
                  EKOR(4,16)=R*R*T*DETT
                  EKOR(6,16)=R*R*S*DETT
                  EKOR(2,17)=R*S*T*DETT
                  EKOR(4,17)=S*S*T*DETT
                  EKOR(5,17)=R*S*S*DETT
                  EKOR(3,18)=R*S*T*DETT
                  EKOR(5,18)=R*T*T*DETT
                  EKOR(6,18)=S*T*T*DETT
               ENDIF
               IF(IALFA.EQ.3) THEN
                  EKOR(1,1)=-R*XJJ(1,1)*DETT
                  EKOR(1,2)=-R*XJJ(1,2)*DETT
                  EKOR(1,3)=-R*XJJ(1,3)*DETT
                  EKOR(2,4)=-S*XJJ(2,1)*DETT
                  EKOR(2,5)=-S*XJJ(2,2)*DETT
                  EKOR(2,6)=-S*XJJ(2,3)*DETT
                  EKOR(3,7)=-T*XJJ(3,1)*DETT
                  EKOR(3,8)=-T*XJJ(3,2)*DETT
                  EKOR(3,9)=-T*XJJ(3,3)*DETT
                  EKOR(4,1)=-R*XJJ(2,1)*DETT
                  EKOR(4,2)=-R*XJJ(2,2)*DETT
                  EKOR(4,3)=-R*XJJ(2,3)*DETT
                  EKOR(4,4)=-S*XJJ(1,1)*DETT
                  EKOR(4,5)=-S*XJJ(1,2)*DETT
                  EKOR(4,6)=-S*XJJ(1,3)*DETT
                  EKOR(5,4)=-S*XJJ(3,1)*DETT
                  EKOR(5,5)=-S*XJJ(3,2)*DETT
                  EKOR(5,6)=-S*XJJ(3,3)*DETT
                  EKOR(5,7)=-T*XJJ(2,1)*DETT
                  EKOR(5,8)=-T*XJJ(2,2)*DETT
                  EKOR(5,9)=-T*XJJ(2,3)*DETT
                  EKOR(6,1)=-R*XJJ(3,1)*DETT
                  EKOR(6,2)=-R*XJJ(3,2)*DETT
                  EKOR(6,3)=-R*XJJ(3,3)*DETT
                  EKOR(6,7)=-T*XJJ(1,1)*DETT
                  EKOR(6,8)=-T*XJJ(1,2)*DETT
                  EKOR(6,9)=-T*XJJ(1,3)*DETT
                  EKOR(1,10)=S*T*XJJ(1,1)*DETT
                  EKOR(1,11)=S*T*XJJ(1,2)*DETT
                  EKOR(1,12)=S*T*XJJ(1,3)*DETT
                  EKOR(2,10)=R*T*XJJ(2,1)*DETT
                  EKOR(2,11)=R*T*XJJ(2,2)*DETT
                  EKOR(2,12)=R*T*XJJ(2,3)*DETT
                  EKOR(3,10)=S*R*XJJ(3,1)*DETT
                  EKOR(3,11)=S*R*XJJ(3,2)*DETT
                  EKOR(3,12)=S*R*XJJ(3,3)*DETT
                  EKOR(4,10)=(S*XJJ(2,1)+R*XJJ(1,1))*T*DETT
                  EKOR(4,11)=(S*XJJ(2,2)+R*XJJ(1,2))*T*DETT
                  EKOR(4,12)=(S*XJJ(2,3)+R*XJJ(1,3))*T*DETT
                  EKOR(5,10)=(T*XJJ(3,1)+S*XJJ(2,1))*R*DETT
                  EKOR(5,11)=(T*XJJ(3,2)+S*XJJ(2,2))*R*DETT
                  EKOR(5,12)=(T*XJJ(3,3)+S*XJJ(2,3))*R*DETT
                  EKOR(6,10)=(T*XJJ(3,1)+R*XJJ(1,1))*S*DETT
                  EKOR(6,11)=(T*XJJ(3,2)+R*XJJ(1,2))*S*DETT
                  EKOR(6,12)=(T*XJJ(3,3)+R*XJJ(1,3))*S*DETT
               ENDIF
               IF(IALFA.EQ.6) THEN
                  EKOR(1,1)=R*DETT
                  EKOR(2,2)=S*DETT
                  EKOR(3,3)=T*DETT
                  EKOR(4,4)=R*DETT
                  EKOR(4,5)=S*DETT
                  EKOR(5,6)=R*DETT
                  EKOR(5,7)=T*DETT
                  EKOR(6,8)=S*DETT
                  EKOR(6,9)=T*DETT
                  EKOR(4,10)=R*T*DETT
                  EKOR(4,11)=S*T*DETT
                  EKOR(5,12)=R*S*DETT
                  EKOR(5,13)=S*T*DETT
                  EKOR(6,14)=R*S*DETT
                  EKOR(6,15)=R*T*DETT
                  EKOR(1,16)=R*S*DETT
                  EKOR(1,17)=R*T*DETT
                  EKOR(2,18)=R*S*DETT
                  EKOR(2,19)=S*T*DETT
                  EKOR(3,20)=R*T*DETT
                  EKOR(3,21)=S*T*DETT
                  EKOR(4,22)=R*S*DETT
                  EKOR(5,23)=R*T*DETT
                  EKOR(6,24)=S*T*DETT
                  EKOR(1,25)=R*S*T*DETT
                  EKOR(2,26)=R*S*T*DETT
                  EKOR(3,27)=R*S*T*DETT
                  EKOR(4,28)=R*S*T*DETT
                  EKOR(5,29)=R*S*T*DETT
                  EKOR(6,30)=R*S*T*DETT
               ENDIF
C
CS             FORMIRANJE MATRICE - G (GLOBALNI DEKARTOV)
CE             FORM MATRIX - G (GLOBAL CARTESIAN)
C
               CALL MNOZM1(GERS,FOT,EKOR,6,LA,6)
            ENDIF
C
CS          FORMIRANJE KONSTITUITIVNE MATRICE
CE          CONSTITUITIVE MATRIX
C
            GO TO ( 70, 70,  3,  4,  5,  5,  5,  5,  5,  5,
     1               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     2               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     3               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     4               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     5               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     6               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     7               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     8               5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     9               5,  5,  5,  5,  5,  5,  5,  5,  5,  5),NMODM
    3       LFUN=MREPER(1)
            LNTA=MREPER(2)
            LTEM=MREPER(3)
            MATE=MREPER(4)
            CALL MEL33(A(LFUN),A(LNTA),A(LTEM),MATE,TGT)
            GO TO 70
    4       LFUN=MREPER(1)
            LNTA=MREPER(2)
            LTEM=MREPER(3)
            MATE=MREPER(4)
            CALL MEL34(A(LFUN),A(LNTA),A(LTEM),MATE,TGT)
CS          MATRICA TRANSFORMACIJE
            CALL TRAETP(ELAST,ELAST,TSG)
            GO TO 70
    5       IBTC=IBTC+1
            KBTC=KBTC+1
            IRAC=2
            NPROS=(NNXYZ+IBTC)*MODPRO( NMODM )*IDVA
!            LPLAS=LPLAST+NPROS     
!            LPLA1=LPLAS1+NPROS
!            LPLA0=LPLAS0+NPROS
            LPLAS=1+NPROS/IDVA     
            LPLA1=1+NPROS/IDVA
            LPLA0=1+NPROS/IDVA
            IF(ICRACK.LT.0) ICRACK=-1
            IF(ICRACK.EQ.-1) THEN
               DO 810 JGI=1,NBRCR
                  IF(NBTC(JGI,1).LT.0) ICRACK=-2
                  IF(KBTC.EQ.IABS(NBTC(JGI,1)).OR.KBTC.EQ.NBTC(JGI,2))
     1              GO TO 820 
  810          CONTINUE
               GO TO 20
  820          JGCR=1
               IF(NBTC(JGI,2).EQ.KBTC) ICRACK=-NBTC(JGI,3)-2
            ENDIF
C  RAD     
            IF(IRAC.EQ.2) THEN
              CALL CLEAR(TA,6)
              CALL CLEAR(STRAIN,6)
            ENDIF
C  RAD     
      CALL MDMAT3(STRAIN,TA,NMODM,IRAC,LPLAS,LPLA1,IBTC,TGT,INTGL,lpla0)
CS          MATRICA TRANSFORMACIJE *** PROBA10
            IF(MODORT(NMODM).EQ.1) THEN
               CALL TRAETP(ELAST,ELAST,TSG)
            ENDIF
   70       WTU=WR*WS*WT*DET
C           FORMIRANJE VEKTORA H ZA MATRICU MASA
            IF(ITER.EQ.0.AND.INDZS.GT.0)THEN
               WD=WTU*GUST
               ZAPRE=ZAPRE+WTU
               AMASA=AMASA+WD
               DO 149 I=1,NCVE
                  AMASC(I)=AMASC(I)+HE(I,1)*WD
  149          CONTINUE
            ELSE
               ZAPRE=ZAPRE+WTU
            ENDIF
C
CS          FORMIRANJE MATRICE BL LINEARNO
CE          FORM LINEAR STRAIN-DISPLACEMENT MATRIX - B
C
            IF(IBOB.EQ.0) THEN
             CALL BETBE3(HE,BLT,NOP)
            ELSE
C            (20), (21)
             IF(NGAUSX.EQ.1) THEN
              CALL FORMBNAD(BOVI,BNAD)
C              CALL FLATMAT(BNAD,BZHU,R,S,T)
C            KONVERTOVANJE IZ ZHU-OVOG 9X24 U PAK 9X24
              CALL BNAD2BNL(BNAD,BNL)
             ELSE
C             (20) BEZ POBOLJSANJA
C              CALL FORMBZHU(BOVI,R,S,T,BZHU)              
C             (44) SA POBOLJSANJEM
              CALL FORMBZH2(BOVI,R,S,T,BZHU)
C            KONVERTOVANJE IZ ZHU-OVOG 9X24 U PAK 6X24
              IF(IC69.EQ.0) CALL BZHU2BLT(BZHU,BLT)
C            KONVERTOVANJE IZ ZHU-OVOG 9X24 U PAK 9X24
              IF(IC69.EQ.1) CALL BZHU2B9(BZHU,BLT9)
             ENDIF
            ENDIF
c            CALL WRR6(BLT,6*24,'BLT ')
C
CS          FORMIRANJE MATRICE B1 LINEARNO
CE          FORM LINEAR B1 MATRIX
C
            IF(IATYP.EQ.2) CALL BETL13(BLT,NOP,LM,UEL,STRAIN,0)
C            IF(IATYP.EQ.2) CALL BETL13(BLT,NOP,LM,RTDT,STRAIN,0)
C
CS    MATRICA TRANSFORM. DEFORMACIJA - TSG (GLOBALNI - LOKALNI DEKARTOV)
CE    STRAIN TRANSFORMATION MATRIX - TSG (GLOBAL - LOCAL CARTESIAN) 
C
C           CALL WRR3(XJJ,9,'XJJ ')
            CALL JEDNA1(XJ0,XJJ,9)
            IF(MODORT(NMODM).NE.1) CALL TRANAL(XJ0,TSG,0)
C           CALL WRR3(XJJ,9,'XJ1 ')
C           CALL WRR6(TSG,36,'TSG ')
C          IF(IST.EQ.0) CALL WRR(TSG,36,'TSG ')
C
C     KONVERTOVANJE MATRICE C 6X6 U MATRICU C 9X9
C      IF(IBOB.GT.0.AND.NGAUSX.EQ.1) CALL C62C9(ELAST,ELAST9)
      IF(IBOB.GT.0) CALL C62C9(ELAST,ELAST9)
C
      IF(ISKNP.NE.1) THEN
C
CS       RACUNANJE DEFORMACIJA U GLOBALNOM DEKARTOVOM
CE       CALCULATE STRAINS IN GLOBAL CARTESIAN COORDINATE SYSTEM
C
C
         CALL CLEAR(STRAIN,6) 
C
CS       LINEARNOST I M.N.O.
CE       LINEAR PART AND M.N.O., E=B*U
C       
        IF(IBOB.GT.0.AND.NGAUSX.EQ.1) THEN
         CALL FORMSTR9(STRAIN9,BNL,UEL)
        ELSE 
         IF(IATYP.EQ.0.OR.IATYP.EQ.1)
     1   CALL MNOZI1(STRAIN,BLT,UEL,6,ND)
        ENDIF
C         call wrr6(blt,nd*6,'blt ')
C         call wrr6(uel,nd,'uel ')
C         call wrr6(strain,6,'stra')
C
CS       GEOMETRIJSKA NELINEARNOST - UKUPNE DEFORMACIJE ZA T.L. I U.L.
CE       GEOMETRICAL NONLINEARITY - TOTAL STRAIN FOR T.L. AND U.L.
C
         IF(IATYP.EQ.2) CALL BETL13(BLT,NOP,LM,UEL,STRAIN,1)
         IF(IATYP.EQ.3) CALL STRUL3(BLT,NOP,LM,UEL,STRAIN)
C         IF(IATYP.EQ.2) CALL BETL13(BLT,NOP,LM,RTDT,STRAIN,1)
C         IF(IATYP.EQ.3) CALL STRUL3(BLT,NOP,LM,RTDT,STRAIN)
C
C     RACA - !!!!!!!!!!!!!!!!!!!!!!!
         IF (NMODM.EQ.31) THEN
C
CS          KOORDINATE U TRENUTKU - 0
CE          COORDINATE IN TIME - 0
C
            II=0
            DO 14 I=1,NCVE
            DO 14 K=1,3
               II=II+1
               IF (IATYP.EQ.1) THEN 
                 CON(I,K)=COR(I,K)+UEL(II)
               ELSE             
                 CON(I,K)=COR(I,K)-UEL(II)
               ENDIF
   14       CONTINUE
C
C           CALL WRR(COR,18,'COR ')
C           CALL WRR(CON,18,'CON ')
C           CALL WRR(UEL,ND,'UEL ')
C
C           ZBOG RACUNANJA VEKTORA GR U T+DT
CS          JAKOBIJEVA MATRICA U GAUSOVOJ TACKI U TRENUTKU - T+DT
CE          JACOBIAN MATRIX IN GAUSS POINT IN TIME - T+DT
C
               IF (IATYP.EQ.1) THEN 
                 CALL GRGSG3(CON,XJJ,HE,NCVE)
               ELSE             
                 CALL GRGSG3(COR,XJJ,HE,NCVE)
               ENDIF
C
              CALL JEDNA1(VMS,XJJ,9) 
C
CS          JAKOBIJEVA MATRICA U GAUSOVOJ TACKI U TRENUTKU - 0
CE          JACOBIAN MATRIX IN GAUSS POINT IN TIME - 0
C
               IF (IATYP.EQ.1) THEN 
                 CALL GRGSG3(COR,XJJ,HE,NCVE)
               ELSE             
                 CALL GRGSG3(CON,XJJ,HE,NCVE)
               ENDIF
C
              CALL JEDNA1(VMS1,XJJ,9) 
C
CS          INVERTOVAN JAKOBIJAN - XJJ
CE          INVERSE JACOBIAN - XJJ
C
            CALL MINV3(XJJ,DET)
C
CS          RACUNANJE GRADIJENTA DEFORMACIJE OD 0 DO T+DT
CE          CALCULATE DEFORMATION GRADIENT FROM 0 TO T+DT
C
               IF (IATYP.EQ.1) THEN 
                  CALL GRAD3T(XJJ,HE,CON,COR,GRAD,NCVE)
               ELSE             
                  CALL GRAD3T(XJJ,HE,COR,CON,GRAD,NCVE)
               ENDIF
C
            CALL JEDNA1(GRAD1R,GRAD,9)
            DLAMZ=GRAD(3,3)
C
CS          INVERTOVAN GRADIJENT - GRAD
CE          INVERSE GRADIENT - GRAD
C
            CALL MINV3(GRAD,DET)
            CALL JEDNA1(RACGR,GRAD,9)
            RACGR(3,3)=DLAMZ 
C
         ENDIF
C     RACA - !!!!!!!!!!!!!!!!!!!!!!!
C
CS       VELIKE DEFORMACIJE
CE       LARGE STRAIN
C
         IF(IATYP.GE.4) THEN
C
CS          STARO - Be, Cp**-1, Fp**-1
CE          OLD - B
C
            IF(ILEDE.EQ.0.AND.ICPM1.LE.1) THEN
               CALL JEDNA1(STRESS,TAU(1,JG,NLM),NLD)
               IF(IST.EQ.1) CALL WRR(STRESS,NLD,'BOLD')
            ENDIF
            IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN 
               CALL JEDNA1(GRAP,TAU(1,JG,NLM),NLD)
               CALL MINV3(GRAP,DUM)
               IF(IST.EQ.1) CALL WRR(GRAP,NLD,'BOLD')
            ENDIF
C
CS          JAKOBIJEVA MATRICA U TACKI (R,S,T) U TRENUTKU - T,0
CE          JACOBIAN MATRIX IN POINT (R,S,T) IN TIME - T,0
C
            CALL GRGSG3(CORT,XJT,HE,NCVE)
            CALL GRGSG3(CON,XJ0,HE,NCVE)
C
CS          INVERTOVAN JAKOBIJAN - XJT,XJ0
CE          INVERSE JACOBIAN - XJT,XJ0
C
            CALL MINV3(XJT,DUM)
            CALL MINV3(XJ0,DUM)
CS          INKOMPATIBILNE DEFORMACIJE
CE          INCOMPATIBILE STRAIN
C
            IF(IALFA.GE.0) THEN
               GRAD(1,1)=-R*ALFE(1,NLM)
               GRAD(2,1)=-R*ALFE(2,NLM)
               GRAD(3,1)=-R*ALFE(3,NLM)
               GRAD(1,2)=-S*ALFE(4,NLM)
               GRAD(2,2)=-S*ALFE(5,NLM)
               GRAD(3,2)=-S*ALFE(6,NLM)
               GRAD(1,3)=-T*ALFE(7,NLM)
               GRAD(2,3)=-T*ALFE(8,NLM)
               GRAD(3,3)=-T*ALFE(9,NLM)
C OVO PROVERITI
               CALL MNOZM3(GRAN,GRAD,XJT,3,3,3)
C               CALL WRR(GRAN,9,'GRAA')
            ENDIF
C
CS          RACUNANJE GRADIJENTA DEFORMACIJE OD T DO T+DT
CE          CALCULATE DEFORMATION GRADIENT FROM T TO T+DT
C
            IF(ILEDE.EQ.0) THEN
               IF(ICPM1.EQ.0) CALL MNOZM4(GRAD,XJJ,XJT,3,3,3)
               IF(ICPM1.GE.1) CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
            ENDIF
            IF(ILEDE.EQ.1) CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
            IF(IST.EQ.1) CALL WRR3(XJJ,9,'XJJ ')
            IF(IST.EQ.1) CALL WRR3(XJT,9,'XJT ')
            IF(IST.EQ.1) CALL WRR3(GRAD,9,'GRAD')
            IF(IALFA.GE.0) CALL ZBIRM1(GRAD,GRAN,9)
         KOREKG=0
         IF(KOREKG.EQ.1) THEN
CS          RACUNANJE NORMIRANOG GRADIJENTA DEFORMACIJE OD T DO T+DT
            CALL DETER3(GRAD,DETG)
            IF(DETG.GT.DETGM) DETGM=DETG
C SA KOREKCIJOM
            DETG=DEXP(-1.D0/3.D0*DLOG(DETG))
C BEZ KOREKCIJE
C            DETG=1.D0
            CALL JEDNAK(GRAN,GRAD,DETG,9)
         ELSE
            CALL JEDNA1(GRAN,GRAD,9)
         ENDIF
C            CALL WRR(GRAN,9,'GRAN')
CS          RACUNANJE GRADIJENTA DEFORMACIJE OD 0 DO T+DT
CE          CALCULATE DEFORMATION GRADIENT FROM 0 TO T+DT
            CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
            CALL DETER3(GRAD,DETG)
            IF(DETG.LT.1.D-15) STOP 'DETG=0, PAK32.FOR'
            IF(IST.EQ.1) CALL WRR3(XJJ,9,'XJJ ')
            IF(IST.EQ.1) CALL WRR3(XJ0,9,'XJ0 ')
            IF(IST.EQ.1) CALL WRR3(GRAD,9,'GRAD')
            IF(IST.EQ.1) WRITE(3,*) 'DETG',DETG
C
CS          TRANSF. ELASTICNOG LEVOG KOSI - GRINOVOG TENZORA DEFORMACIJE
CE          TRANSFORM. ELASTIC LEFT CAUCHY - GREEN DEFORMATION TENSOR
C           Be* = Ft * Be * FtT  (ICPM1=0)
C           Be* = F * (Cp**-1) * FT  (ICPM1=1)
            IF(ILEDE.EQ.0.AND.ICPM1.LE.1) CALL PIOKOS(GRAN,STRESS)
CS          RACUNANJE ELASTICNOG GRADIJENTA DEFORMACIJE
CE          ELASTIC DEFORMATION GRADIENT
C           Fe = F * Fp**-1         Be = Fe * FeT
            IF(ILEDE.EQ.0.AND.ICPM1.EQ.2) THEN
C  GRAD ILI GRAN
               CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
               CALL MNOZM3(GRAP,GRAE,GRAE,3,3,3)
               CALL TENVEK(GRAP,STRESS,1)
               CALL JEDNA1(XJ,GRAE,9)
            ENDIF
C           Fe = F * Fp**-1         Ce = FeT * Fe
            IF(ILEDE.EQ.1) THEN
C  GRAD ILI GRAN
               CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
               CALL MNOZM2(GRAP,GRAE,GRAE,3,3,3)
               CALL TENVEK(GRAP,STRESS,1)
               CALL JEDNA1(XJ,GRAE,9)
            ENDIF
C
C           GLAVNE VREDNOSTI I PRAVCI
C
            IF(IST.EQ.1) CALL WRR6(STRESS,6,'BECE')
            IF(IGLPR.EQ.1) THEN
               CALL GLAVN(STRESS)
               CALL GLAPR3(STRESS,QP)
               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PP  ') 
               IF(IST.EQ.1) CALL WRR3(QP,9,'QP  ') 
            ENDIF
C
CS          PROBNE ELASTICNE DEFORMACIJE
CE          TRIAL ELASTIC STRAIN
C
            IF(IATYP.EQ.4) THEN
               IF(IGLPR.EQ.1) THEN
                  CALL CLEAR(STRESS,6)
                  STRESS(1)=PRINC(1)
                  STRESS(2)=PRINC(2)
                  STRESS(3)=PRINC(3)
               ENDIF
               CALL UKDEFV(STRAIN,STRESS)
            ENDIF
            IF(IATYP.EQ.5) THEN
C              GLAVNE VREDNOSTI
C              LAMBDA
               P1=DSQRT(PRINC(1))
               P2=DSQRT(PRINC(2))
               P3=DSQRT(PRINC(3))
               IF(IST.EQ.1) WRITE (3,*) 'L1,L2,L3',P1,P2,P3
               STRAIN(1)=DLOG(P1)
               STRAIN(2)=DLOG(P2)
               STRAIN(3)=DLOG(P3)
            ENDIF
CD          ZA INTEGRACIJU U DEKARTOVOM SISTEMU (INTGL.EQ.1)
C           TRANSFORMACIJA DEFORMACIJA: GLAVNI PRAVCI - DEKARTOV SISTEM
            IF(IGLPR.EQ.1) THEN
               IF(INTGL.EQ.1) CALL DIJADS(QP,STRAIN)
               STRAIN(4)=2.D0*STRAIN(4)
               STRAIN(5)=2.D0*STRAIN(5)
               STRAIN(6)=2.D0*STRAIN(6)
            ENDIF
            IF(IST.EQ.1) CALL WRR3(STRAIN,6,'STDE')
         ENDIF
C
CS       INKOMPATIBILNE DEFORMACIJE
CE       INCOMPATIBILE STRAIN
C                          OVAJ USLOV PROVERITI
         IF(IALFA.GE.0.AND.IATYP.LT.4)
     1   CALL MNOZI1(STRAIN,GERS,ALFE(1,NLM),6,LA)
C         call wrr6(gers,6*la,'gers')
C         call wrr6(alfe(1,nlm),la,'alfa')
C         CALL WRR6(STRAIN,6,'DEFO')
C
CS       TERMICKE DEFORMACIJE   ETH=ALFA*(T-T0)
CE       THERMAL STRAINS
C
         IF (NMODM.EQ.3.OR.NMODM.EQ.4) CALL MAMOD3(STRAIN,TGT)
C
C
CS       RACUNANJE NAPONA
CE       CALCULATE STRESSES 
C
C
         CALL CLEAR(TA,6) 
C
         IF(NMODM.LT.5) THEN
C
CS          MATERIJALNA LINEARNOST
CE          STRESS CALCULATION FOR LINEAR MATERIAL MODELS, S=C*E 
C
C            CALL WRR6(STRAIN,6,'STRA')
C            CALL WRR6(ELAST,36,'ELAS')
         IF(IBOB.GT.0.AND.NGAUSX.EQ.1) THEN
          CALL FORMSIG9(SIGMA9,ELAST9,STRAIN9)
         ELSE
            CALL MNOZI1(TA,ELAST,STRAIN,6,6)
         ENDIF
C            CALL WRR6(TA,6,'TA  ')
C
CS          RACUNANJE UNUTRASNJIH SILA
CE          CALCULATE INTERNAL FORCES
C           r = BT * S 
C

         IF(IBOB.GT.0.AND.NGAUSX.EQ.1) THEN
          CALL FORMFE(FE,BNAD,SIGMA9)
         ELSE
            IF(NZADP.GT.0.AND.ISKNP.EQ.2.AND.NGENL.EQ.0)
     1      CALL INTEV2(FE,BLT,TA,WTU,ND,6)
C     1      CALL INTEGF(FTDT,BLT,TA,LM,WTU,ND,6)
         ENDIF
C
CS          CISCENJE NUMERICKIH GRESAKA ZA NAPONE
CE          CLEANING NUMERICAL ERRORS FOR STRESS
C
            CALL CISTIN(TA,6)
C
            DO 21 I=1,6
CS             KOSIJEVI NAPONI U GLOBAL. ILI LOKAL. DEKART. SISTEMU
CE             CAUCHY STRESS IN GLOBAL OR LOCAL CARTESIAN SYSTEM 
               CALL GLLOKN(TA,TSG,TAU(I,JG,NLM),ISN,I)
   21       CONTINUE
            IF(SNAP.LT.TAU(2,JG,NLM)) SNAP=TAU(2,JG,NLM)
C            SNAP=SNAP+TAU(2,JG,NLM)
         ELSE
C
CS          NAPONI ZA PLASTICAN MODEL
CE          STRESS FOR MATERIAL NONLINEARITY 
C
C            CALL WRR6(STRAIN,6,'STRU')
            IRAC=1
      CALL MDMAT3(STRAIN,TA,NMODM,IRAC,LPLAS,LPLA1,IBTC,TGT,INTGL,lpla0)
            IF(SNAP.LT.TA(2)) SNAP=TA(2)
C            SNAP=SNAP+TA(2)
            IF(IST.EQ.1) CALL WRR3(STRAIN,6,'DEFL')
            IF(IST.EQ.1) CALL WRR3(TA,6,'STRL')
C            CALL WRR(ELAST,36,'ELAS')
CS          TRANSFORMACIJA ELAST
            IF(ISKNP.NE.2.AND.MODORT(NMODM).EQ.1) THEN
C OVO PROVERITI ZA ANIZOTROPNE KOD 3D I GREDE
               CALL TRAETP(ELAST,ELAST,TSG)
            ENDIF
C
CS          PROMENA ELASTICNOG LEVOG KOSI - GRINOVOG TENZORA DEFORMACIJE
CE          UPDATE ELASTIC LEFT CAUCHY - GREEN DEFORMATION TENSOR
C
C  OVO PROVERITI ????????????????????????????
            IF(NAPON.EQ.1.AND.IATYP.GE.4) THEN
C               CALL WRR6(TAU(1,JG,NLM),NLD,'BST ')
               IF(ILEDE.EQ.0) THEN
                  IF(ICPM1.EQ.0) CALL JEDNA1(TAU(1,JG,NLM),STRAIN,NLD)
                  IF(ICPM1.EQ.1) THEN
CS                   TRANSF. ELAS. LEVOG KOSI - GRINOVOG TENZORA DEFOR.
CS                   U DESNI PLASTICNI INVERZNI
CE                   TRANSFORM. ELASTIC LEFT CAUCHY - GREEN DEFOR. TENS.
C                    CP**-1 = F**-1 * Be * F**-T
                     CALL MINV3(GRAN,DUM)
C                     CALL WRR(GRAN,6,'GR-1')
                     CALL PIOKOS(GRAN,STRAIN)
C                     CALL WRR6(STRAIN,6,'CP-1')
C                     CALL WRR6(TAU(1,JG,NLM),6,'BST ')
                     CALL JEDNA1(TAU(1,JG,NLM),STRAIN,6)
C                     CALL WRR6(TAU(1,JG,NLM),6,'BNOV')
C                    PROVERA PLASTICNE ROTACIJE
C                     IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                        WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C                        CALL WRR3(PRINC,3,'PRI0')
C                        CALL WRR3(QP,9,'QP0 ')
C                        CALL GLAVN(STRAIN)
C                        CALL GLAPR3(STRAIN,QP)
C                        CALL WRR3(PRINC,3,'PRIN')
C                        CALL WRR3(QP,9,'QP  ')
C                     ENDIF
                  ENDIF
               ENDIF
               IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN 
                  CALL JEDNA1(GRAD,TAU(1,JG,NLM),NLD)
                  P1=DEXP(STRAIN(1))
                  P2=DEXP(STRAIN(2))
                  P3=DEXP(STRAIN(3))
                  CALL DIJAD(GRAN,QP,QP,P1,P2,P3)
                  CALL MNOZM1(GRAP,GRAN,GRAD,3,3,3)
                  CALL JEDNA1(TAU(1,JG,NLM),GRAP,NLD)
C                 PROVERA PLASTICNE ROTACIJE
C                  IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                     WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C                     CALL WRR3(GRAD,9,'GRAD')
C                     CALL WRR3(PRINC,3,'PRI0')
C                     CALL WRR3(QP,9,'QP0 ')
C                     CALL WRR3(GRAP,9,'GRAP')
C                     CALL MNOZM2(GRAN,GRAP,GRAP,3,3,3)
C                     CALL WRR3(GRAN,9,'GRAN')
C                     CALL TENVEK(GRAN,STRESS,1)
C                     CALL GLAVN(STRESS)
C                     CALL GLAPR3(STRESS,QP)
C                     CALL WRR3(PRINC,3,'PRIN')
C                     CALL WRR3(QP,9,'QP  ')
C                     CALL MNOZM3(GRAN,GRAP,GRAP,3,3,3)
C                     CALL WRR3(GRAN,9,'GRA1')
C                     CALL TENVEK(GRAN,STRESS,1)
C                     CALL GLAVN(STRESS)
C                     CALL GLAPR3(STRESS,QP)
C                     CALL WRR3(PRINC,3,'PRI1')
C                     CALL WRR3(QP,9,'QP1 ')
C                  ENDIF
                  NAPKO=0
                  IF(NAPKO.EQ.1) THEN
                     IF(IATYP.EQ.4.AND.ILEDE.EQ.1) THEN
                        CALL MINV3(GRAP,DETP)
                        IF(DETP.GT.DETGP) DETGP=DETP
                        CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
CS                      TRANSF. PIOLA KIRKOFOV - KOSIJEV NAPON 
CE                      TRANSFORM. PIOLA KIRCKOF - CAUCHY STRESS
C                       s = F * S * FT
                        CALL PIOKOS(GRAE,TA)
                        CALL CEPMT(ELAST,CT,0)
C                       Cmnop = Fmi Fnj Fok Fpl Cijkl
                        CALL RRRRC(ELAST,CT,GRAE,1)
                     ENDIF
                  ENDIF
               ENDIF
               IF(IST.EQ.1) CALL WRR6(TAU(1,JG,NLM),NLD,'BNOV')
            ENDIF
C
         ENDIF
C
C
CS       RACUNANJE UNUTRASNJIH SILA
CE       CALCULATE INTERNAL FORCES 
C
C
         IF(NGENL.GT.0) THEN
C
CS          INTEGRACIJA VEKTORA - h
CE          INTEGRATE VEKTOR - h
C           h = GT * S
C
            IF(IALFA.GE.0)
     1      CALL INTEV2(HAEM(1,NLM),GERS,TA,WTU,LA,6)
C
CS          RACUNANJE UNUTRASNJIH SILA
CE          CALCULATE INTERNAL FORCES
C           r = BT * S 
C
         IF(IBOB.GT.0.AND.NGAUSX.EQ.1) THEN
          CALL FORMFE(FE,BNAD,SIGMA9)
         ELSE
            CALL INTEV2(FE,BLT,TA,WTU,ND,6)
C            CALL INTEGF(FTDT,BLT,TA,LM,WTU,ND,6)
         ENDIF
C
         ENDIF
C
         IF(IPG.GT.0) THEN
            CALL JEDNAK(STRAIN,TA,WTU,6)
            CALL MNOZI2(ESILA(1,IPG),BLT,STRAIN,ND,6)
         ENDIF
C
      ENDIF
C
C
CS    INTEGRACIJA MATRICE KRUTOSTI ELEMENTA - SKE
CE    INTEGRATE ELEMENT STIFFNESS MATRIX - SKE
C     SKE = BT * C * B
C
C
      IF(ISKNP.EQ.2) GO TO 20
CC      WRITE(3,'(6G10.4)')((ELAST(I,J),J=1,6),I=1,6)
      IF(IBOB.GT.0.AND.NGAUSX.EQ.1) THEN
       CALL FORMK(SKE,BNL,ELAST9,LM)
      ELSE
       IF(IC69.EQ.0) CALL INTEGK(SKE,BLT,ELAST,LM,WTU,ND,6)
       IF(IC69.EQ.1) CALL INTEGK(SKE,BLT9,ELAST9,LM,WTU,ND,9)
      ENDIF
C
CS    GEOMETRIJSKI NELINEARAN DEO MATRICE KRUTOSTI
CE    GEOMETRIC NONLINEAR PART OF STIFFNESS MATRIX   
C
      IF(IATYP.GT.1.AND.INDKOV.NE.2) CALL KNL3(SKE,HE,NOP,LM,TA,WTU,ND)
C      IF(IATYP.GT.1) CALL KNL3(SKE,HE,NOP,LM,TA,WTU,ND)
C
      IF(IALFA.GE.0) THEN
C
CS       PROIZVOD MATRICA - GT * C
CE       MULTIPLY MATRIX - GT * C
C
         CALL MNOZM2(CEGE,GERS,ELAST,LA,6,6)
C
CS       INTEGRACIJA MATRICE - H
CE       INTEGRATE MATRIX - H
C        H = GT * C * G
C
         CALL INTEM1(HINV(1,1,NLM),CEGE,GERS,WTU,LA,LA,6)
C
CS       INTEGRACIJA MATRICE - GAMA
CE       INTEGRATE MATRIX - GAMA
C        GAMA = GT * C * B
C
         CALL INTEM1(GEEK(1,1,NLM),CEGE,BLT,WTU,LA,ND,6)
      ENDIF
C
   20    CONTINUE
C KOSOVO
         IF(NAPON.EQ.1.AND.INDDTH.EQ.1.AND.IPG.GT.0) 
     +   CALL KOSOVO(SNAP,1,NLM,AU(LTDTH),VREME,BETA)
C     +   CALL KOSOVO(SNAP,NGS12,NLM,AU(LTDTH),VREME,BETA)
C KOSOVO
            IF(ICRACK.LT.0.AND.JGCR.EQ.0) GO TO 10 
C
CS-------------------------- KRAJ PETLJE PO GAUSOVIM TACKAMA --------
CE-------------------------- END  LOOP  OVER  GAUSS  POINTS  --------
C
         IF(IALFA.GE.0) THEN
C
CS          INVERTOVANJE MATRICE - H
CE          INVERSE MATRIX - H
C
            IF(ISKNP.NE.2) CALL MINV(HINV(1,1,NLM),LA,DET,LJA,MJA)
C
            IF(NGENL.GT.0) THEN
C
CS             RACUNANJE IZRAZA - (H**-1) * h
CE             CALCULATE EXPRESSION - (H**-1) * h
C
               CALL CLEAR(HAEML,LA)
               CALL MNOZI1(HAEML,HINV(1,1,NLM),HAEM(1,NLM),LA,LA)
C
CS             KOREKCIJA UNUTRASNJIH SILA
CE             CORECTION INTERNAL FORCES
C              r = r - FT * (H**-1) * h
C
               CALL MNOZ2I(FE,GEEK(1,1,NLM),HAEML,ND,LA)
C               CALL INTEGF(FTDT,GEEK(1,1,NLM),HAEML,LM,-1.D0,ND,LA)
C
               IF(IPG.GT.0) 
     1         CALL MNOZ2I(ESILA(1,IPG),GEEK(1,1,NLM),HAEML,ND,LA)
C
            ENDIF
C
CS          KOREKCIJA MATRICE KRUTOSTI
CE          CORECTION STIFFNESS MATRIX
C           K = K - FT * (H**-1) * F
C
            IF(ISKNP.NE.2)
     1    CALL INTEGK(SKE,GEEK(1,1,NLM),HINV(1,1,NLM),LM,-1.D0,ND,LA)
C
         ENDIF
C
CS       RASPOREDJIVANJE MATRICE KRUTOSTI (SKE)
CE       ASSEMBLE STIFFNESS MATRIX
C
         IF(ISKNP.NE.2.AND.IAKUS.EQ.1) CALL AKUST(SKE,NCVE,NLM,NWE)
         IF(ISKNP.NE.2) CALL SPAKUJ(A(LSK),A(LMAXA),SKE,LM,ND)
         CALL SILPAK(FTDT,FE,LM,ND,A(LCMPC),A(LMPC))
C        RASPOREDJIVANJE KONCENTRISANIH MATRICE MASA AMASC U VEKTOR ZAPS
         IF(ITER.EQ.0.AND.INDZS.GT.0)THEN
C        WRITE(3,*) 'AMASA,NCVE,ND,NLM',AMASA,NCVE,ND,NLM
C        CALL WRR(AMASC,NCVE,'AMAS')
            IK=0
            DO 500 I=1,NCVE
            DO 500 K=1,3
               IK=IK+1
               NJ=LM(IK)
               IF(NJ.GT.0) THEN
                  NPRNJ=NPRZ(NJ)
                  IF(NPRNJ.GT.0.AND.NPRNJ.NE.K) STOP 'NPRNJ'
                  NPRZ(NJ)=K
                  ZAPS(NJ)=ZAPS(NJ)+AMASC(I)
               ENDIF
  500       CONTINUE
         ENDIF
         IF(ITER.EQ.0.AND.INDZS.GT.0) THEN
            GRZAP=GRZAP+ZAPRE
            GRMAS=GRMAS+AMASA
         ELSE
            GRZAP=GRZAP+ZAPRE
         ENDIF
   10 CONTINUE
      IF(ICRACK.GT.0) GO TO 900
      IF(ITER.EQ.0.AND.INDZS.GT.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NGE,GRZAP,NGE,GRMAS
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NGE,GRZAP,NGE,GRMAS
      ELSE
C      WRITE(3,*) 'DETGM,DETGP',DETGM,DETGP
      WRITE(3,*) 'NGE,GRZAP,ITER',NGE,GRZAP,ITER
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(///
     111X,'GRUPA ELEMENATA',I5,' ZAPREMINA =',1PD12.5/
     111X,'GRUPA ELEMENATA',I5,'      MASA =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(///
     111X,'ELEMENT GROUP',I5,' VOLUME =',1PD12.5/
     111X,'ELEMENT GROUP',I5,'   MASS =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C     RUTINA ZA FORMIRANJE MATRICE C(3,3) (4) ZHU
C=======================================================================
      SUBROUTINE CMAT(RST,COR,CM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RST(8,3),COR(21,3),CM(3,3)
      DO I=1,3
       DO J=1,3
        CM(I,J)=DOT(RST(1,J),COR(1,I),8)/8
       ENDDO
      ENDDO
      CALL MINV3(CM,D)
      RETURN
      END
C=======================================================================
C     RUTINA ZA FORMIRANJE MATRICE KOJA SADRZI B-OVE BOVI(8,7,3)
C     BOVI(1-8=CVOROVI ; 1=BEZ INDEKSA,2-4=R-S-T,5-7=RS-RT-ST ; 1-3=X-Y-Z)
C=======================================================================
      SUBROUTINE NAPUNIB(CM,RST,HOVI,COR,BOVI,GAME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CM(3,3),RST(8,3),HOVI(8,4),COR(21,3),
     +          BOVI(8,7,3),GAME(8,4)
      DIMENSION POM(8,3)
C
C     BX,BY,BZ
      ONE8=1./8
      DO I=1,3
       CALL ZBIR3(BOVI(1,1,I),RST(1,1),RST(1,2),RST(1,3),
     +          CM(1,I)/8,CM(2,I)/8,CM(3,I)/8,8)
      ENDDO
C
C     GAME
      DO I=1,4
       DO J=1,3
        CALL JEDNAK( POM(1,J), BOVI(1,1,J), DOT(HOVI(1,I),
     +                  COR(1,J),8), 8)
       ENDDO
       CALL ZBIR4(GAME(1,I),HOVI(1,I),POM(1,1),POM(1,2),POM(1,3),
     +          ONE8,-ONE8,-ONE8,-ONE8,8)       
      ENDDO
C       CALL ZBIR4(DELTA,SJED,BOVI(1,1,1),BOVI(1,1,2),BOVI(1,1,3),
C     1      ONE8,-DOT(SJED,COR(1,1),8)/8,-DOT(SJED,COR(1,2),8)/8,
C     1      -DOT(SJED,COR(1,3),8)/8,8)
      RETURN
      END
C=======================================================================
C     RUTINA ZA DOPUNJAVANJE MATRICE KOJA SADRZI B-OVE BOVI(8,7,3)
C     BOVI(1-8=CVOROVI ; 1=BEZ INDEKSA,2-4=R-S-T,5-7=RS-RT-ST ; 1-3=X-Y-Z)
C=======================================================================
      SUBROUTINE DOPUNIB(BOVI,GAME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION BOVI(8,7,3),GAME(8,4)
      DIMENSION PG1(8),PG2(8)
C
C     B(X-Y-Z)(R-S-T)
      DO I=1,3
       DO J=1,3
        IF(I.EQ.1) THEN
         PJ1=XJ(J,2)
         PJ2=XJ(J,3)
         CALL JEDNA1(PG1,GAME(1,1),8)
         CALL JEDNA1(PG2,GAME(1,2),8)
        ELSE 
         IF(I.EQ.2) THEN
          PJ1=XJ(J,1)
          PJ2=XJ(J,3)
          CALL JEDNA1(PG1,GAME(1,1),8)
          CALL JEDNA1(PG2,GAME(1,3),8)
         ELSE
          PJ1=XJ(J,1)
          PJ2=XJ(J,2)
          CALL JEDNA1(PG1,GAME(1,2),8)
          CALL JEDNA1(PG2,GAME(1,3),8)
         ENDIF
        ENDIF
        CALL ZBIR2(BOVI(1,1+I,J),PG1,PG2,PJ1,PJ2,8)
       ENDDO
      ENDDO
C
C     B(X-Y-Z)(RS-RT-ST)
      DO I=1,3
       DO J=1,3
        CALL JEDNAK(BOVI(1,4+I,J),GAME(1,4),XJ(J,4-I),8)
       ENDDO
      ENDDO
      RETURN
      END
C======================================================================
      SUBROUTINE JACBOB(COR,RST,HOVI,R,S,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       FOR BOBAN JACOBIAN MATRIX IN CURRENT
CE.       INTEGRATION POINT (R,S,T - ARE NATURAL COORDINATES)
C .
C ......................................................................
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      DIMENSION COR(21,3),RST(8,3),HOVI(8,4)
      DIMENSION POM(8)
C
      ONE=1.
      DO I=1,3
       IF(I.EQ.1) THEN
        CALL ZBIR4(POM,RST(1,1),HOVI(1,1),HOVI(1,2),
     +                  HOVI(1,4),ONE,S,T,S*T,8)
       ELSE 
        IF(I.EQ.2) THEN
         CALL ZBIR4(POM,RST(1,2),HOVI(1,1),HOVI(1,3),
     +                  HOVI(1,4),ONE,R,T,R*T,8)
        ELSE
         CALL ZBIR4(POM,RST(1,3),HOVI(1,2),HOVI(1,3),
     +                  HOVI(1,4),ONE,R,S,R*S,8)
        ENDIF
       ENDIF
       DO J=1,3
        XJ(I,J)=DOT(POM,COR(1,J),8)/8
       ENDDO
      ENDDO
      CALL JEDNA1(XJJ,XJ,9)
      CALL MINV3(XJ,DET)
      RETURN
      END
C======================================================================
C=======================================================================
C     RUTINA ZA FORMIRANJE ZHU-OVE MATRICE B BEZ POBOLJSANJA
C=======================================================================
      SUBROUTINE FORMBZHU(BOVI,R,S,T,BZHU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BOVI(8,7,3),BZHU(9,24)
C
      CALL CLEAR(BZHU,9*24)
      DO I=1,3
       DO J=1,8
        BZHU(I,J)=BOVI(J,1,I)+R*BOVI(J,2,I)+S*BOVI(J,3,I)+T*BOVI(J,4,I)+
     1            R*S*BOVI(J,5,I)+R*T*BOVI(J,6,I)+S*T*BOVI(J,7,I)
        BZHU(I+3,J+8)=BZHU(I,J)
        BZHU(I+6,J+16)=BZHU(I,J)
       ENDDO
      ENDDO
      RETURN
      END
C======================================================================
C=======================================================================
C     RUTINA ZA FORMIRANJE ZHU-OVIH MATRICA BNAD
C=======================================================================
      SUBROUTINE FORMBNAD(BOVI,BNAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      DIMENSION BOVI(8,7,3),BNAD(9,24,7)
C
C     **********
      E1=2.D0/3
      E2=-1.D0/3
      E3=E2
c      BETA=0.2
C     **********
C
      CALL CLEAR(BNAD,9*24*7)
C     FORMIRANJE B0
      DO I=1,3
       DO J=1,8
        BNAD(I,J,1)=BOVI(J,1,I)
        BNAD(I+3,J+8,1)=BNAD(I,J,1)
        BNAD(I+6,J+16,1)=BNAD(I,J,1)
       ENDDO
      ENDDO
C
      DO I=2,7
       CALL DPODM(BOVI,1,E1,BNAD,1,1,I)
       CALL DPODM(BOVI,2,E3,BNAD,1,9,I)
       CALL DPODM(BOVI,3,E2,BNAD,1,17,I)
       CALL DPODM(BOVI,2,BETA,BNAD,2,1,I)
       CALL DPODM(BOVI,3,BETA,BNAD,3,1,I)
C                    
       CALL DPODM(BOVI,1,E2,BNAD,5,1,I)
       CALL DPODM(BOVI,2,E1,BNAD,5,9,I)
       CALL DPODM(BOVI,3,E3,BNAD,5,17,I)
       CALL DPODM(BOVI,1,BETA,BNAD,4,9,I)
       CALL DPODM(BOVI,3,BETA,BNAD,6,9,I)
C              
       CALL DPODM(BOVI,1,E3,BNAD,9,1,I)
       CALL DPODM(BOVI,2,E2,BNAD,9,9,I)
       CALL DPODM(BOVI,3,E1,BNAD,9,17,I)
       CALL DPODM(BOVI,1,BETA,BNAD,7,17,I)
       CALL DPODM(BOVI,2,BETA,BNAD,8,17,I)
      ENDDO
C      
      RETURN
      END
C=======================================================================
C     RUTINA ZA FORMIRANJE ZHU-OVE MATRICE B SA POBOLJSANJEM
C=======================================================================
      SUBROUTINE FORMBZH2(BOVI,R,S,T,BZHU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      DIMENSION BOVI(8,7,3),BZHU(9,24)
C
C     **********
      E1=2.D0/3
      E2=-1.D0/3
      E3=E2
c      BETA=0.2
C     **********
C
      CALL CLEAR(BZHU,9*24)
      DO I=1,3
       DO J=1,8
        BZHU(I,J)=BOVI(J,1,I)
        BZHU(I+3,J+8)=BZHU(I,J)
        BZHU(I+6,J+16)=BZHU(I,J)
       ENDDO
      ENDDO
C
      CALL DPODM2(BOVI,1,E1,BZHU,1,1,R,S,T)
      CALL DPODM2(BOVI,2,E3,BZHU,1,9,R,S,T)
      CALL DPODM2(BOVI,3,E2,BZHU,1,17,R,S,T)
      CALL DPODM2(BOVI,2,BETA,BZHU,2,1,R,S,T)
      CALL DPODM2(BOVI,3,BETA,BZHU,3,1,R,S,T)
C              
      CALL DPODM2(BOVI,1,E2,BZHU,5,1,R,S,T)
      CALL DPODM2(BOVI,2,E1,BZHU,5,9,R,S,T)
      CALL DPODM2(BOVI,3,E3,BZHU,5,17,R,S,T)
      CALL DPODM2(BOVI,1,BETA,BZHU,4,9,R,S,T)
      CALL DPODM2(BOVI,3,BETA,BZHU,6,9,R,S,T)
C              
      CALL DPODM2(BOVI,1,E3,BZHU,9,1,R,S,T)
      CALL DPODM2(BOVI,2,E2,BZHU,9,9,R,S,T)
      CALL DPODM2(BOVI,3,E1,BZHU,9,17,R,S,T)
      CALL DPODM2(BOVI,1,BETA,BZHU,7,17,R,S,T)
      CALL DPODM2(BOVI,2,BETA,BZHU,8,17,R,S,T)
C      
      RETURN
      END
C=======================================================================
C     RUTINA ZA DODAVANJE PODMATRICE U MATRICU BNAD
C=======================================================================
      SUBROUTINE DPODM(BOVI,IXYZ,COEF,BNAD,IROW,ICOL,INDEX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BOVI(8,7,3),BNAD(9,24,7)
C
      DO I=1,8
       BNAD(IROW,ICOL+I-1,INDEX)=COEF*BOVI(I,INDEX,IXYZ)
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA DODAVANJE PODMATRICE U MATRICU B
C=======================================================================
      SUBROUTINE DPODM2(BOVI,IXYZ,COEF,BZHU,IROW,ICOL,R,S,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BOVI(8,7,3),BZHU(9,24)
C
      DO I=1,8
       BZHU(IROW,ICOL+I-1)=BZHU(IROW,ICOL+I-1)+COEF*
     1        (R*BOVI(I,2,IXYZ)+S*BOVI(I,3,IXYZ)+T*BOVI(I,4,IXYZ)+
     1         R*S*BOVI(I,5,IXYZ)+R*T*BOVI(I,6,IXYZ)+S*T*BOVI(I,7,IXYZ))
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA FLAT-OVANJE MATRICE BNAD U MATRICU BZHU
C=======================================================================
      SUBROUTINE FLATMAT(BNAD,BZHU,R,S,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BNAD(9,24,7),BZHU(9,24)
C
      DO I=1,9
       DO J=1,24
        BZHU(I,J)=BNAD(I,J,1)+R*BNAD(I,J,2)+S*BNAD(I,J,3)+T*BNAD(I,J,4)+
     1            R*S*BNAD(I,J,5)+R*T*BNAD(I,J,6)+S*T*BNAD(I,J,7)
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA KONVERTOVANJE ZHU-OVE MATRICE 9X24 B U PAK-OVU 6X24
C=======================================================================
      SUBROUTINE BZHU2BLT(BZHU,BLT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BZHU(9,24),BLT(6,*)
C
      CALL COPYBZ(BZHU,BLT,1,1)
      CALL COPYBZ(BZHU,BLT,5,2)
      CALL COPYBZ(BZHU,BLT,9,3)
      CALL COPYBZ2(BZHU,BLT,2,4,4)
      CALL COPYBZ2(BZHU,BLT,6,8,5)
      CALL COPYBZ2(BZHU,BLT,3,7,6)
      RETURN
      END
C=======================================================================
C     RUTINA ZA KONVERTOVANJE ZHU-OVE MATRICE 9X24 B U PAK-OVU 9X24
C=======================================================================
      SUBROUTINE BZHU2B9(BZHU,BLT9)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BZHU(9,24),BLT9(9,24)
C
      DO I=1,9
       DO J=1,8
        DO K=1,3
         BLT9(I,(J-1)*3+K)=BZHU(I,(K-1)*8+J)
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA KOPIRANJE ZHU-OVE VRSTE U PAK-OVU
C=======================================================================
      SUBROUTINE COPYBZ(BZHU,BLT,IROWZHU,IROWPAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BZHU(9,24),BLT(6,*)
C
      DO I=1,8
       DO J=1,3
        BLT(IROWPAK,(I-1)*3+J)=BZHU(IROWZHU,(J-1)*8+I)
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA KOPIRANJE ZBIRA DVE ZHU-OVE VRSTE U PAK-OVU
C=======================================================================
      SUBROUTINE COPYBZ2(BZHU,BLT,IROWZHU1,IROWZHU2,IROWPAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BZHU(9,24),BLT(6,*)
C
      DO I=1,8
       DO J=1,3
        BLT(IROWPAK,(I-1)*3+J)=BZHU(IROWZHU1,(J-1)*8+I)+
     +                  BZHU(IROWZHU2,(J-1)*8+I)
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA PUNJENJE INTERPOLACIONE MATRICE (N U TEORIJI)
C=======================================================================
      SUBROUTINE NAPUNIHE(HE,HOVI,BOVI,SJED,RST,R,S,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      DIMENSION HE(NCVE,*),HOVI(8,4),BOVI(8,7,3),SJED(8),RST(8,3)
C
      DO I=1,8
       HE(I,1)=SJED(I)+R*RST(I,1)+S*RST(I,2)+T*RST(I,3)+
     1         HOVI(I,1)*R*S+HOVI(I,2)*R*T+HOVI(I,3)*S*T+
     1         HOVI(I,4)*R*S*T
       HE(I,1)=HE(I,1)/8
      ENDDO
C
      DO J=1,3
       DO I=1,8
        HE(I,J+1)=BOVI(I,1,J)+R*BOVI(I,2,J)+S*BOVI(I,3,J)+T*BOVI(I,4,J)+
     1          BOVI(I,5,J)*R*S+BOVI(I,6,J)*R*T+BOVI(I,7,J)*S*T
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA KONVERTOVANJE C 6X6 U C 9X9
C=======================================================================
      SUBROUTINE C62C9(ELAST,ELAST9)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ELAST(6,6),ELAST9(9,9),IA(9)
      DATA IA/1,4,6,4,2,5,6,5,3/
C
      DO I=1,9
       DO J=1,9
        ELAST9(I,J)=ELAST(IA(I),IA(J))
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA KONVERTOVANJE ZHU-OVE MATRICE 9X24 B U PAK-OVU 9X24
C=======================================================================
      SUBROUTINE BNAD2BNL(BNAD,BNL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BNAD(9,24,7),BNL(9,24,7)
C
      DO M=1,7
       DO I=1,9
        DO J=1,8
         DO K=1,3
          BNL(I,(J-1)*3+K,M)=BNAD(I,(K-1)*8+J,M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA RACUNANJE DEFORMACIJA
C=======================================================================
      SUBROUTINE FORMSTR9(STRAIN9,BNL,UEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION STRAIN9(9,7),BNL(9,24,7),UEL(24)
C
      CALL CLEAR(STRAIN9,9*7)
      DO I=1,7
       CALL MNOZI1(STRAIN9(1,I),BNL(1,1,I),UEL,9,24)
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA RACUNANJE NAPONA
C=======================================================================
      SUBROUTINE FORMSIG9(SIGMA9,ELAST9,STRAIN9)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SIGMA9(9,7),ELAST9(9,9),STRAIN9(9,7)
C
      CALL CLEAR(SIGMA9,9*7)
      DO I=1,7
       CALL MNOZI1(SIGMA9(1,I),ELAST9,STRAIN9(1,I),9,9)
      ENDDO
      RETURN
      END
C=======================================================================
C     RUTINA ZA RACUNANJE SILA
C=======================================================================
      SUBROUTINE FORMFE(FE,BNAD,SIGMA9)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION FE(24),BNAD(9,24,7),SIGMA9(9,7),S1(24),S2(24)
C
      V=DET*8
      CALL MNOZM2(FE,BNAD(1,1,1),SIGMA9(1,1),24,1,9)
      DO I=1,24
       FE(I)=FE(I)*V
      ENDDO
C
      CALL CLEAR(S1,24)
      CALL CLEAR(S2,24)
      DO I=2,4
       CALL MNOZI2(S1,BNAD(1,1,I),SIGMA9(1,I),24,9)
      ENDDO
      DO I=5,7
       CALL MNOZI2(S2,BNAD(1,1,I),SIGMA9(1,I),24,9)
      ENDDO
      CALL ZBIRM2(FE,S1,S2,V/3,V/9,24)
      RETURN
      END
C=======================================================================
C     RUTINA ZA RACUNANJE MATRICE KRUTOSTI
C=======================================================================
      SUBROUTINE FORMK(SKE,BNL,ELAST9,LM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION SKE(*),BNL(9,24,7),ELAST9(9,9),LM(*),S1(300),S2(300)
C
      V=DET*8
C      CALL CLEAR(SKE,300)
      CALL INTEGK(SKE,BNL(1,1,1),ELAST9,LM,V,24,9)
C
      CALL CLEAR(S1,300)
      CALL CLEAR(S2,300)
      DO I=2,4
       CALL INTEGK(S1,BNL(1,1,I),ELAST9,LM,V/3,24,9)
      ENDDO
      DO I=5,7
       CALL INTEGK(S2,BNL(1,1,I),ELAST9,LM,V/9,24,9)
      ENDDO
      CALL ZBIRM2(SKE,S1,S2,1.D0,1.D0,300)
      RETURN
      END
C=======================================================================
C     RUTINA ZA FORMIRANJE C+S
C=======================================================================
      SUBROUTINE DODS(ELAST9,SIGMA9,R,S,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ELAST9(9,9),SIGMA9(9,7)
C
      DO I=1,3
       DO J=1,3
        S=SIGMA9((J-1)*3+I,1)+
     1              R*SIGMA9((J-1)*3+I,2)+S*SIGMA9((J-1)*3+I,3)+
     1              T*SIGMA9((J-1)*3+I,4)+SIGMA9((J-1)*3+I,5)*R*S+
     1              SIGMA9((J-1)*3+I,6)*R*T+SIGMA9((J-1)*3+I,7)*S*T
        ELAST9(I,J)=ELAST9(I,J)+S
        ELAST9(I+3,J+3)=ELAST9(I+3,J+3)+S
        ELAST9(I+6,J+6)=ELAST9(I+6,J+6)+S
       ENDDO
      ENDDO
      RETURN
      END
