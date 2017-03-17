C=======================================================================
C
CE        PRINTOUT STRESS 3/D ELEMENTS
CS        STAMPANJE NAPONA 3/D ELEMENATA
C
C   SUBROUTINE K21NAP
C              STANP3
C              STANG3
C              STAN31
C              STAN33
C              STAN35
C              STAN39
C              STAG31
C
C=======================================================================
      SUBROUTINE K21NAP
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3D ELEMENTS
CS.   P R O G R A M
CS.       ZA POZIVANJE PROGRAMA ZA STAMPANJE NAPONA 3/D ELEMENATA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /restap/ irestp,lplas0,lstaz0
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /OBNOVA/ IRUSI
      COMMON /RESTAR/ TSTART,IREST
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' K21NAP'
      LAU=LMAX
CS    UCITAVANJE PODATAKA O ELEMENTIMA I NAPONA SA DISKA
CE    READING DATA FROM RECORDS OF THE FILE (ZIELEM)
      CALL READE3(A(LAU))
CS    PODPROGRAM ZA STAMPANJE NAPONA
CE    PRINT STRESSES TO OUTPUT FILE (*.LST)
      IF(INDPR.EQ.0.OR.INDPR.EQ.2) CALL STANP3(A(LAU))
CS    PODPROGRAM ZA STAMPANJE NAPONA ZA POSTPROCESIRANJE
CE    PRINT STRESS FOR POSTPROCESSING TO THE FILE (*.UNV)
      IF(INDGR.EQ.0.OR.INDGR.EQ.2) CALL STANG3(A(LAU))
      IF(NMODM.GT.4) THEN
         NPROS=NE*NGS12*MODPRO( NMODM )*IDVA
C         CALL WRR(A(LPLAS1),NPROS,'LP23')
C
C        PRIPREMA ZA RESTART CAM-CLY (SOPSTVENA TEZINA)
C
c        ponistavanje samo posle prvog pustanja
         IF(IRUSI.EQ.1.AND.IRESTP.EQ.0.AND.KOR.EQ.NDT) then
            LMA8=LSTAZ0-1
            CALL WRITDD(PLAS1,NPROS/IDVA,IELEM,LMA8,LDUZI)
!            CALL WRITDD(A(PLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
            CALL NULL39(PLAS1,A(LAU))
!            CALL NULL39(A(LPLAS1),A(LAU))
         endif
C
         LMA8=LSTAZA(3)-1
         CALL WRITDD(PLAS1,NPROS/IDVA,IELEM,LMA8,LDUZI)
         CALL WRITDD(PLAS1,NPROS/IDVA,IELEM,LMA8,LDUZI)
!         CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
!         CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
      ENDIF
      IF ((NMODM.GE.5) .AND. (IATYP.GT.0)) THEN
      DEALLOCATE(PLAST)
      DEALLOCATE(PLAS1)
      IF(IRUSI.EQ.1) DEALLOCATE(PLAS0)
      ENDIF
!      IF(IATYP.EQ.0) RETURN
!      IF(NMODM.LE.4) RETURN
!      NPROS=NE*NGS12*MODPRO( NMODM )*IDVA
!      LMA8=LSTAZA(3)-1
!C      CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
      RETURN
      END
C======================================================================
      SUBROUTINE STANP3(AU)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO CALL ROUTINES FOR PRINTOUT STRESSES FOR 3D ELEMENTS,(*.LST)
CS.   P R O G R A M
CS.       ZA POZIVANJE PROGRAMA ZA STAMPANJE NAPONA 3/D ELEMENATA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /GAUSVR/ LTEMGT,LCORGT,ICORGT
      COMMON /CVSILE/ NSILA,LESILA
      COMMON /DUPLAP/ IDVA
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STANP3'
C
      LDEFOR=LSIGMA
      IF(NMODM.LE.4) LDEFOR=LSIGMA+NE*NGS12*6*IDVA
C
      GO TO (     1,  1,  3,  3,  5,  5,  9,  9,  9,  5,
     1          999,999,  5,  5,  5,  5,  5,  5,  5,  5,
     2            5,  9, 23,  5,  9,  5,  5,  9,  5, 30,
     3           31,999,999,999,999,999,999,999,999,999,
     4            9,  9,  9,  9,  9,999,999,999,999,999,
     5          999,  5,  5,  5,999,999,999,999,999,999,
     6            5,999,999,999,999,999,999,999,999,999,
     7          999,999,999,999,999,999,999,999,999,999,
     8          999,999,999,999,999,999,999,999,999,999,
     9          999,999,999,999,999,999,999,999,999,999),NMODM
C
CE    PRINT STRESSES FOR MATERIAL MODEL 1 AND 2
    1 CALL STAN31(A(LSIGMA),A(LTEMGT),A(LCORGT),A(LAU),AU(LISNA),
     1            AU(LIPRC),AU(LCEL),AU(LNEL),A(LCVEL),ICVEL,  0,
     1            A(LDEFOR))
      GO TO 100
C
    3 CALL STAN31(A(LSIGMA),A(LTEMGT),A(LCORGT),A(LAU),AU(LISNA),
     1            AU(LIPRC),AU(LCEL),AU(LNEL),A(LCVEL),ICVEL,  1,
     1            A(LDEFOR))
      GO TO 100
C
    5 CALL STAN35(PLAS1,A(LCORGT),A(LAU),AU(LISNA),AU(LIPRC),
     1            AU(LCEL),ICVEL)
      GO TO 100
C
    9 CALL STAN39(PLAS1,A(LCORGT),A(LAU),AU(LISNA),AU(LIPRC),
     1            AU(LCEL),ICVEL,
     1            AU(LNMAT),AU(LNSLOJ),AU(LMATSL),AU(LBBET),AU(LDSLOJ))
      GO TO 100
C
   23 CALL STA323(PLAS1,A(LCORGT),A(LAU),AU(LISNA),AU(LIPRC),
     1            AU(LCEL),ICVEL)
      GO TO 100
C
   30 CALL STA330(PLAS1,A(LCORGT),A(LAU),AU(LISNA),AU(LIPRC),
     1            AU(LCEL),ICVEL)
      GO TO 100
C
   31 CALL STA331(PLAS1,A(LCORGT),A(LAU),AU(LISNA),AU(LIPRC),
     1            AU(LCEL),ICVEL)
C
  100 IF(NSILA.GT.0) CALL STSIL3(A(LAU),AU(LESILA),AU(LIPGC),
     1                           AU(LCEL),AU(LNEL),A(LCVEL),ICVEL)
C
  999 RETURN
      END
C======================================================================
      SUBROUTINE STANG3(AU)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO CALL ROUTINES FOR PRINTOUT STRESSES FOR 3D ELEMENTS
CE.       FOR GRAPHICS TO FILE (*.UNV)
CS.   P R O G R A M
CS.       ZA POZIVANJE PROGRAMA ZA STAMPANJE NAPONA 3/D ELEMENATA
CS.       ZA GRAFIKU
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /NIDEAS/ IDEAS
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STANG3'
      GO TO (     1,  1,  3,  3,  5,  5,  5,  5,  5,  5,
     1          999,999,  5,  5,  5,  5,  5,  5,  5,  5,
     2            5,999,999,  5, 25,  5,  5,  5,  5,999,
     3            5,999,999,999,999,999,999,999,999,999,
     4            5,  5,  5,  5,  5,999,999,999,999,999,
     5          999,999, 53,999,999,999,999,999,999,999,
     6            5,999,999,999,999,999,999,999,999,999,
     7          999,999,999,999,999,999,999,999,999,999,
     8          999,999,999,999,999,999,999,999,999,999,
     9          999,999,999,999,999,999,999,999,999,999),NMODM
CE    PRINT STRESSES FOR MATERIAL MODEL 1 AND 2
    1 IF(IDEAS.GE.0)
     1CALL STAG31(A(LSIGMA),A(LAU),AU(LISNA),AU(LCEL),ICVEL)
      CALL STAU31(A(LSIGMA),A(LAU),AU(LISNA),AU(LCEL),ICVEL,A(LSIGMA),0)
      RETURN
    3 IF(IDEAS.GE.0)
     1CALL STAG31(A(LSIGMA),A(LAU),AU(LISNA),AU(LCEL),ICVEL)
      CALL STAU31(A(LSIGMA),A(LAU),AU(LISNA),AU(LCEL),ICVEL,A(LSIGMA),0)
      RETURN
    5 ISUME=ISUMEL
      IF(IDEAS.GE.0)
     1CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,2)
      CALL STAU35(A(LSIGMA),AU(LISNA),AU(LCEL),ICVEL,PLAS1,A(LAU),1,
     1            AU(LNMAT))
      RETURN
   53 ISUME=ISUMEL
      IF(IDEAS.GE.0)
     1CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,2)
      CALL STAU353(A(LSIGMA),AU(LISNA),AU(LCEL),ICVEL,PLAS1,A(LAU)
     1            ,1,AU(LNMAT))
      IF(ISTDE.NE.-1) THEN
         ISUMEL=ISUME 
         IF(IDEAS.GE.0)
     1   CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,3)
      ENDIF
      IF(NMODM.EQ.26) RETURN
      IF(NMODM.NE.15.AND.NMODM.NE.18.AND.NMODM.NE.61) THEN
         ISUMEL=ISUME 
         IF(IDEAS.GE.0)
     1   CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,0)
      ENDIF
      IF(NMODM.EQ.9) THEN
         ISUMEL=ISUME 
         IF(IDEAS.GE.0)
     1   CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,99)
      ENDIF
      IF(NMODM.EQ.15.OR.NMODM.EQ.16.OR.NMODM.EQ.18.OR.NMODM.EQ.19) THEN
         ISUMEL=ISUME 
         IF(IDEAS.GE.0)
     1   CALL STAG35(PLAS1,A(LAU),AU(LISNA),AU(LCEL),ICVEL,1)
      ENDIF
      RETURN
   25 ISUME=ISUMEL
      IF(IDEAS.GE.0)
     1CALL STG325(PLAS1,A(LAU),AU(LNMAT),AU(LNSLOJ),AU(LISNA),
     &            AU(LCEL),ICVEL,2)
      CALL STAU35(A(LSIGMA),AU(LISNA),AU(LCEL),ICVEL,PLAS1,A(LAU),1,
     1            AU(LNMAT))
      ISUMEL=ISUME 
      IF(IDEAS.GE.0)
     1CALL STG325(PLAS1,A(LAU),AU(LNMAT),AU(LNSLOJ),AU(LISNA),
     &            AU(LCEL),ICVEL,3)
      ISUMEL=ISUME 
      IF(IDEAS.GE.0)
     1CALL STG325(PLAS1,A(LAU),AU(LNMAT),AU(LNSLOJ),AU(LISNA),
     &            AU(LCEL),ICVEL,0)
  999 RETURN
      END
C======================================================================
      SUBROUTINE STAN31(TAU,TEMGT,CORGT,AU,ISNA,IPRC,MCVEL,NEL,NCVEL,
     &                  ICVEL,INDCAL,DEF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 1 AND 2
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 1 I 2
C .
C ......................................................................
C
      CHARACTER*11 IMD
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /PRINCI/ PRINC(3)
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /CDEBUG/ IDEBUG
      DIMENSION AU(*)
      REAL AU
      DIMENSION TAU(6,NGS12,*),TEMGT(NGS12,*),CORGT(3,NGS12,*),NEL(NE,*)
     1          ,DEF(6,NGS12,*)
      DIMENSION ISNA(*),IPRC(*),MCVEL(*),NCVEL(*),STRESL(6),STRESS(7),
     1          STRAIN(7),SRED(6),SRED2(6),V(3,3),princs(3),vs(3,3)
      DIMENSION XC(3),LR(21),LS(21),LT(21)
      DATA LR/1,2,2,1,1,2,2,1, 3,2,3,1,3,2,3,1, 1,2,2,1, 3/
      DATA LS/1,1,2,2,1,1,2,2, 1,3,2,3,1,3,2,3, 1,1,2,2, 3/
      DATA LT/1,1,1,1,2,2,2,2, 1,1,1,1,2,2,2,2, 3,3,3,3, 3/
      DATA XC/1.D0,-1.D0,0.D0/
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN31'
C PRIVREMENO Rakic za Micuna
cr         IDIS=78
cr         IMD='DEFORMACIJE'
cr         OPEN (IDIS,FILE=IMD,STATUS='NEW',FORM='FORMATTED',
cr     1         ACCESS='SEQUENTIAL')
C PRIVREMENO     
      NGR1=1
      NGS1=1
      NGT1=1
      NNLM=0
      INDSN=0
      INNOLD=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         ISN=ISNA(NLM)
         INN=0
         IF(ISN.GE.10)THEN
            ISN=ISN-10
            INN=1
         ENDIF
         IF(INN.NE.INNOLD) INDSN=0
         IF(ISN.EQ.1) GO TO 20
         IF(INDSN.EQ.0.AND.INN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NGE
            ENDIF
         ENDIF
         IF(INDSN.EQ.0.AND.INN.EQ.1) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2013) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6013) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2012) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6012) NGE
            ENDIF
         ENDIF
C
         INDSN=1
         INNOLD=INN
         IPC=IPRC(NLM)
         NMM=NLM
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
C
C  G A U S O V E    T A C K E
C
      IF(INN.EQ.0) THEN
         JG=0
         C=1./NGS12
         SEFE1=0.D0
         SEFE2=0.D0
         CALL CLEAR(SRED,6)
         CALL CLEAR(SRED2,6)
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            JG=JG+1
            CALL JEDNA1(STRESS,TAU(1,JG,NLM),6)
            IF(ISTNA.EQ.1) CALL ZBIRM1(SRED,STRESS,6)
            CALL GLAVN3(STRESS)
            CALL GLAPR3(STRESS,Vs)
            call jedna1(princs,princ,3)
               IF(ISTNA.EQ.2) THEN
                  IF(STRESS(7).GT.SEFE1) THEN
                     NGR1=NGR
                     NGS1=NGS
                     NGT1=NGT
                     SEFE1=STRESS(7)
                     CALL JEDNA1(SRED,STRESS,6)
                  ENDIF
               ENDIF
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NMM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               NNCV=0
               SEFE=STRESS(7)
            ENDIF
            IF(ISTNA.EQ.0) THEN
               WRITE(IZLAZ,5001) NGR,NGS,NGT,(STRESS(K),K=1,6)
            ENDIF
            IF(ISTDE.NE.-1) THEN
               CALL JEDNA1(STRAIN,DEF(1,JG,NLM),6)
               IF(ISTNA.EQ.1) CALL ZBIRM1(SRED2,STRAIN,6)
               CALL GLAVN3(STRAIN)
               CALL GLAPR3(STRain,V)
               IF(ISTNA.EQ.2) THEN
               CALL GLAVN3(STRAIN)
                  IF(STRAIN(7).GT.SEFE2) THEN
                     SEFE2=STRAIN(7)
                     CALL JEDNA1(SRED2,STRAIN,6)
                  ENDIF
               ENDIF
               IF(ISTNA.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (STRAIN(K),K=1,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (STRAIN(K),K=1,6)
cr      IF(NGR.EQ.2.AND.NGS.EQ.2.AND.NGT.EQ.2) THEN
cr      WRITE(IDIS,8002) NLM,(STRAIN(K),K=1,6)
cr 8002 FORMAT(I5,6(1PE11.3))
cr      ENDIF
               ENDIF
            ENDIF
            IF(ISTNA.EQ.0) THEN
               IF(INDCAL.EQ.0)THEN
                  WRITE(IZLAZ,5002) STRESS(7),PRINCs
               ELSEIF(INDCAL.EQ.1)THEN
                  WRITE(IZLAZ,5002) STRESS(7),PRINCs,TEMGT(JG,NLM)
               ENDIF
               WRITE(IZLAZ,5004) 1,(Vs(I,1),I=1,3)
               WRITE(IZLAZ,5004) 2,(Vs(I,2),I=1,3)
               WRITE(IZLAZ,5004) 3,(Vs(I,3),I=1,3)
               
            IF(ISTDE.NE.-1) THEN
               WRITE(IZLAZ,5002) STRain(7),PRINC
               WRITE(IZLAZ,5004) 1,(V(I,1),I=1,3)
               WRITE(IZLAZ,5004) 2,(V(I,2),I=1,3)
               WRITE(IZLAZ,5004) 3,(V(I,3),I=1,3)
            endif
               IF(IPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) (CORGT(K,JG,NLM),K=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) (CORGT(K,JG,NLM),K=1,3)
               ENDIF
            ENDIF
   10    CONTINUE
         IF(ISTNA.EQ.1) THEN
            CALL JEDNAK(STRESS,SRED,C,6)
            WRITE(IZLAZ,1010) (STRESS(I),I=1,6) 
         ENDIF
         IF(ISTNA.EQ.2) THEN
            WRITE(IZLAZ,5001) NGR1,NGS1,NGT1,(SRED(I),I=1,6) 
         ENDIF
         IF(ISTDE.NE.-1) THEN
            IF(ISTNA.EQ.1) THEN 
               CALL JEDNAK(STRAIN,SRED2,C,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (STRAIN(K),K=1,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (STRAIN(K),K=1,6)
            ENDIF
            IF(ISTNA.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (SRED2(K),K=1,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (SRED2(K),K=1,6)
            ENDIF
         ENDIF
      ENDIF
C
C  C V O R O V I
C
      IF(INN.EQ.1) THEN
         DO 12 N=1,NCVE
C ...TRANSFORMACIJA ZA ANIZOTROPIJU ??
            NN=NEL(NLM,N)
            IF(NN.EQ.0) GO TO 12
            IF(ICVEL.EQ.1) NN=NCVEL(NN)
            CALL SMOOTH(TAU(1,1,NLM),STRESL,6,NGAUSX,NGAUSY,NGAUSZ,
     1                  XC(LR(N)),XC(LS(N)),XC(LT(N)))
            CALL JEDNA1(STRESS,STRESL,6)
            CALL GLAVN3(STRESS)
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NMM
               NNGR=0
               NNGS=0
               NNGT=0
               NNCV=NN
               SEFE=STRESS(7)
            ENDIF
            WRITE(IZLAZ,5003) NN,(STRESS(K),K=1,6)
            IF(ISTDE.NE.-1) THEN
               CALL SMOOTH(DEF(1,1,NLM),STRAIN,6,NGAUSX,NGAUSY,NGAUSZ,
     1                     XC(LR(N)),XC(LS(N)),XC(LT(N)))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (STRAIN(K),K=1,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (STRAIN(K),K=1,6)
            ENDIF
            WRITE(IZLAZ,5002) STRESS(7),PRINC
   12    CONTINUE
      ENDIF
C
   20 CONTINUE
      IF(NNLM.GT.0)THEN
         IF(NNCV.EQ.0)THEN
      IF(ISRPS.EQ.0.AND.INDSN.EQ.1)
     1WRITE(IZLAZ,2030) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1.AND.INDSN.EQ.1)
     1WRITE(IZLAZ,6030) NNLM,NNGR,NNGS,NNGT,SEFE
         ELSE
      IF(ISRPS.EQ.0.AND.INDSN.EQ.1)
     1WRITE(IZLAZ,2035) NNLM,NNCV,SEFE
      IF(ISRPS.EQ.1.AND.INDSN.EQ.1)
     1WRITE(IZLAZ,6035) NNLM,NNCV,SEFE
         ENDIF
         IF(SMKOR.LT.SEFE) THEN
            SMKOR=SEFE
            NGRKOR=NGE
            NSMKOR=NNLM
         ENDIF
         IF(SMALL.LT.SEFE) THEN
            SMALL=SEFE
            NGRALL=NGE
            NSMALL=NNLM
            KSMALL=KOR
         ENDIF
      ENDIF
      RETURN
 5000 FORMAT(/I8)
 5001 FORMAT(3I3,6(1PE11.3))
 5002 FORMAT(9X,6(1PE11.3))
 5003 FORMAT(I8,1X,6(1PE11.3))
 5004 FORMAT(9X,'V',I1,3(1PE11.3))
 1010 FORMAT(9X,6(1PE11.3))
C-----------------------------------------------------------------------
 2002 FORMAT(' UK.DEFO.',6(1PE11.3))
 2010 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX'/
     112X,'EF.NAPON',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)',
     14X,'TEMPER.')
 2011 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-RR',2X,' NAPON-SS',2X,
     1' NAPON-TT',2X,' NAPON-RS',2X,' NAPON-ST',2X,' NAPON-TR'/
     112X,'EF.NAPON',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)',
     14X,'TEMPER.')
 2012 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' C V O R    NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX'/
     112X,'EF.NAPON',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)')
 2013 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' C V O R    NAPON-RR',2X,' NAPON-SS',2X,
     1' NAPON-TT',2X,' NAPON-RS',2X,' NAPON-ST',2X,' NAPON-TR'/
     112X,'EF.NAPON',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)')
 2020 FORMAT(
     1' KOORD.-X=',1PE10.3,2X,' KOORD.-Y=',1PE10.3,2X,' KOORD.-Z=',
     11PE10.3)
 2030 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2035 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,
     1'   CVOR =',I6,'  MAX.EFE.NAP. =',1PE11.3//)
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6002 FORMAT(' TO.STRA.',6(1PE11.3))
 6010 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX'/
     111X,'EF.STRESS',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)',
     14X,'TEMPER.')
 6011 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR'/
     111X,'EF.STRESS',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)',
     14X,'TEMPER.')
 6012 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' N O D E   STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX'/
     111X,'EF.STRESS',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)')
 6013 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' N O D E   STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR'/
     111X,'EF.STRESS',3X,'PRINC(1)',3X,'PRINC(2)',3X,'PRINC(3)')
 6020 FORMAT(
     1' COORD.-X=',1PE10.3,2X,' COORD.-Y=',1PE10.3,2X,' COORD.-Z=',
     11PE10.3)
 6030 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  IT =',I3,'  MAX EFF STR  =',1PE11.3//)
 6035 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,
     1'   NODE =',I6,'  MAX EFF STR  =',1PE11.3//)
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE STAN35(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 5
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 5
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*)
      DIMENSION STRESS(7),STRES(7),SRED(6),SRED2(6),SRED3(6),SRED4(6)
      DIMENSION CORGT(3,NGS12,*)
      DIMENSION AU(*)
      REAL AU
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN35'
      NGR1=1
      NGS1=1
      NGT1=1
      NNLM=0
      NGXYZ=NGS12
      NPR56=MODPRO( NMODM )
      INDSN=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         NPROS=(NLM-1)*NGXYZ-1
         IBTC=0
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
C        IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10.AND.IT.EQ.0) THEN
C           IT=1
C           CALL MATRTE
C        ELSE
C           GLOBAL STRESS
            ISN=2
C        ENDIF
         IF(INDSN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) NGE
            ENDIF
            IF(ISN.EQ.2) THEN

      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NGE
            ENDIF
         ENDIF
         INDSN=1
         IPC=IPRC(NLM)
         NMM=NLM
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
         C=1./NGS12
         SEFE1=0.D0
         SEFE2=0.D0
         SEFE3=0.D0
         SEFE4=0.D0
         CALL CLEAR(SRED,6)
         CALL CLEAR(SRED2,6)
         CALL CLEAR(SRED3,6)
         CALL CLEAR(SRED4,6)
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NPROS+IBTC)*NPR56
            IPL=0
            IF(NMODM.NE.31) IPL=PLAST(LL+NPR56-1)
C       IF(IBTC.EQ.1) WRITE(44,*)-PLAST(LL+8),-PLAST(LL+2)
            CALL JEDNA1(STRES,PLAST(LL),6)
C           IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10) THEN
C              DO 33 II=1,6
C                 STRESS(II)=0.D0
C              DO 33 K=1,6
C                STRESS(II)=STRESS(II)+TE(II,K)*STRES(K)
C  33          CONTINUE
C           ENDIF
C            IF(ISN.EQ.2) CALL JEDNA1(STRESS,STRES,6)
            CALL JEDNA1(STRESS,STRES,6)
            IF(ISTNA.EQ.1) CALL ZBIRM1(SRED,STRESS,6)
            CALL GLAVN3(STRESS)
            IF(ISTNA.EQ.2) THEN
               IF(STRESS(7).GT.SEFE1) THEN
                  NGR1=NGR
                  NGS1=NGS
                  NGT1=NGT
                  SEFE1=STRESS(7)
                  CALL JEDNA1(SRED,STRESS,6)
               ENDIF
            ENDIF
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NLM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               SEFE=STRESS(7)
            ENDIF
            IF(ISTNA.EQ.0) 
     +      WRITE(IZLAZ,5001) NGR,NGS,NGT,(STRESS(K),K=1,6)
            IF(ISTDE.NE.-1) THEN
               LL6=LL+6
               CALL JEDNA1(STRES,PLAST(LL6),6)
               IF(ISTNA.EQ.1) CALL ZBIRM1(SRED2,STRES,6)
               CALL GLAVN3(STRES)
               IF(ISTNA.EQ.2) THEN
                  IF(STRES(7).GT.SEFE2) THEN
                     SEFE2=STRES(7)
                     CALL JEDNA1(SRED2,STRES,6)
                  ENDIF
               ENDIF
               IF(ISTNA.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (PLAST(II),II=LL+6,LL+11)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (PLAST(II),II=LL+6,LL+11)
               ENDIF
            ENDIF
            IPRPL=0
C              Isotropic Damage Model (Oliver 1996); PRINT DAMAGE
               IF(NMODM.EQ.61) THEN
      IF(ISRPS.EQ.0)
     1          WRITE(IZLAZ,2008) PLAST(LL+12),PLAST(LL+13),PLAST(LL+14)
      IF(ISRPS.EQ.1)
     1          WRITE(IZLAZ,6008) PLAST(LL+12),PLAST(LL+13),PLAST(LL+14)
C
      GO TO 136
               ENDIF
C
            IF(NMODM.EQ.26) GO TO 31
            DO 30 II=LL+12,LL+17
               IF(PLAST(II).GT.1.0D-10) IPRPL=1
   30       CONTINUE
            IF(IPRPL.EQ.1) THEN
               LL6=LL+12
               CALL JEDNA1(STRES,PLAST(LL6),6)
               IF(ISTNA.EQ.1) CALL ZBIRM1(SRED3,STRES,6)
               CALL GLAVN3(STRES)
               IF(ISTNA.EQ.2) THEN
                  IF(STRES(7).GT.SEFE3) THEN
                     SEFE3=STRES(7)
                     CALL JEDNA1(SRED3,STRES,6)
                  ENDIF
               ENDIF
               IF(ISTNA.EQ.0) THEN
                  IF(NMODM.NE.15.AND.NMODM.NE.18) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) (PLAST(II),II=LL+12,LL+17)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) (PLAST(II),II=LL+12,LL+17)
                  ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (PLAST(II),II=LL+12,LL+17)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (PLAST(II),II=LL+12,LL+17)
                  ENDIF
               ENDIF
            ENDIF
  136       IPRCR=0
            IF(NMODM.EQ.16.OR.NMODM.EQ.19) THEN
               DO 40 II=LL+24,LL+29
                  IF(PLAST(II).GT.1.0D-10) IPRCR=1
   40          CONTINUE
               IF(IPRCR.EQ.1) THEN
                  LL6=LL+24
                  CALL JEDNA1(STRES,PLAST(LL6),6)
                  IF(ISTNA.EQ.1) CALL ZBIRM1(SRED4,STRES,6)
                  CALL GLAVN3(STRES)
                  IF(ISTNA.EQ.2) THEN
                     IF(STRES(7).GT.SEFE4) THEN
                        SEFE4=STRES(7)
                        CALL JEDNA1(SRED4,STRES,6)
                     ENDIF
                  ENDIF
                  IF(ISTNA.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (PLAST(II),II=LL+24,LL+29)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (PLAST(II),II=LL+24,LL+29)
                  ENDIF
               ENDIF
            ENDIF
   31       IF(ISTNA.EQ.0) THEN
               IF(IPC.EQ.0.OR.IPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2015) STRESS(7),IPL
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6015) STRESS(7),IPL
               ENDIF
            ENDIF
   10    CONTINUE
         IF(ISTNA.EQ.1) THEN
            CALL JEDNAK(STRESS,SRED,C,6)
            WRITE(IZLAZ,1010) (STRESS(I),I=1,6) 
            IF(ISTDE.NE.-1) THEN
               CALL JEDNAK(STRESS,SRED2,C,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (STRESS(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (STRESS(I),I=1,6) 
            ENDIF
            IF(NMODM.EQ.26) GO TO 20
            IF(IPRPL.EQ.1) THEN
               CALL JEDNAK(STRESS,SRED3,C,6)
               IF(NMODM.NE.15.AND.NMODM.NE.18) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) (STRESS(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) (STRESS(I),I=1,6) 
               ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (STRESS(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (STRESS(I),I=1,6) 
               ENDIF
            ENDIF
            IF(IPRCR.EQ.1) THEN
               CALL JEDNAK(STRESS,SRED4,C,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (STRESS(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (STRESS(I),I=1,6) 
            ENDIF
         ENDIF
         IF(ISTNA.EQ.2) THEN
            WRITE(IZLAZ,5001) NGR1,NGS1,NGT1,(SRED(I),I=1,6) 
            IF(ISTDE.NE.-1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (SRED2(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (SRED2(I),I=1,6) 
            ENDIF
            IF(NMODM.EQ.26) GO TO 20
            IF(IPRPL.EQ.1) THEN
               IF(NMODM.NE.15.AND.NMODM.NE.18) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) (SRED3(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) (SRED3(I),I=1,6)
               ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (SRED3(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (SRED3(I),I=1,6) 
               ENDIF
            ENDIF
            IF(IPRCR.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (SRED4(I),I=1,6) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (SRED4(I),I=1,6) 
            ENDIF
         ENDIF
   20 CONTINUE
      IF(NNLM.GT.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NNLM,NNGR,NNGS,NNGT,SEFE
         IF(SMKOR.LT.SEFE) THEN
            SMKOR=SEFE
            NGRKOR=NGE
            NSMKOR=NNLM
         ENDIF
         IF(SMALL.LT.SEFE) THEN
            SMALL=SEFE
            NGRALL=NGE
            NSMALL=NNLM
            KSMALL=KOR
         ENDIF
      ENDIF
      RETURN
 5000 FORMAT(/I8)
 5001 FORMAT(3I3,6(1PE11.3))
 1010 FORMAT(9X,6(1PE11.3))
C-----------------------------------------------------------------------
 2010 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2011 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-RR',2X,' NAPON-SS',2X,
     1' NAPON-TT',2X,' NAPON-RS',2X,' NAPON-ST',2X,' NAPON-TR')
 2015 FORMAT(' EF.NAPON=',1PE10.3,'   STATUS=',I10)
 2030 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2002 FORMAT(' UK.DEFO.',6(1PE11.3))
 2003 FORMAT(' PL.DEFO.',6(1PE11.3))
 2004 FORMAT(' PU.DEFO.',6(1PE11.3))
 2008 FORMAT(' OSTECENJE DN=',1PE11.4,'  RN=',1PE11.4,'  QN=',1PE11.4)
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6010 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX')
 6011 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR')
 6015 FORMAT(' EF.STRE.=',1PE10.3,'   STATUS=',I10)
 6030 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  IT =',I3,'  MAX EFF STR  =',1PE11.3//)
 6002 FORMAT(' TO.STRA.',6(1PE11.3))
 6003 FORMAT(' PL.STRA.',6(1PE11.3))
 6004 FORMAT(' CR.STRA.',6(1PE11.3))
 6008 FORMAT(' DAMAGE DN= ',1PE11.4,'  RN= ',1PE11.4,'  QN= ',1PE11.4)
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE STAN39(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL,
     1                  NMAT,NSLOJ,MATSL,BBET,DSLOJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 5
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 5
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /PRINCI/ PRINC(3)
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NSLOJ(*),MATSL(MSLOJ,*),BBET(MSLOJ,*),DSLOJ(MSLOJ,*)
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*),STRESS(7),STRES(7)
      DIMENSION CORGT(3,NGS12,*),NMAT(*)
      DIMENSION AU(*)
      REAL AU
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN39'
c
	imi=91
      if(kor.eq.1) then
         OPEN(imi,file='dijagram',STATUS='unknown',FORM='FORMATTED',
     1        ACCESS='SEQUENTIAL')
      if(nmodm.eq.9)
     1write(imi,*)'        Sxx,       Exx,       Po,       Sm,       Ev
     1         J2d,     Eyy,     Ezz '
      if(nmodm.eq.7.OR.nmodm.eq.41.OR.nmodm.eq.42.OR.nmodm.eq.43
     1  .OR.nmodm.eq.44)
     1write(imi,*)'        Sxx,       Exx,        X,       I1,       Ev
     1         J2d,       L,      Syy'
      endif
      NNLM=0
      INDSN=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
        MAT=NMAT(NLM)
        NNSL=1
        IF(MSET.GT.0) NNSL=NSLOJ(MAT)
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         IBTC=0
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
C        IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10.AND.IT.EQ.0) THEN
C           IT=1
C           CALL MATRTE
C        ELSE
            ISN=2
C        ENDIF
         IF(INDSN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NGE
            ENDIF
         ENDIF
         INDSN=1
         IPC=IPRC(NLM)
         NMM=NLM
         NMMPR=MODPRO(NMODM)
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NLM-1)*NGS12*NMMPR+(IBTC-1)*NMMPR
            CALL JEDNA1(STRES,PLAST(LL),6)
C           IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10) THEN
C              DO 33 II=1,6
C                 STRESS(II)=0.D0
C              DO 33 K=1,6
C                STRESS(II)=STRESS(II)+TE(II,K)*STRES(K)
C  33          CONTINUE
C           ENDIF
            IF(ISN.EQ.2) CALL JEDNA1(STRESS,STRES,6)
            CALL GLAVN3(STRESS)
C           RACUNANJE GLAVNIH NAPONA I DSQRT(3*J2D)
            SM=-(PRINC(3)+PRINC(2)+PRINC(1))/3.
            S1D=-PRINC(3)-SM
            S2D=-PRINC(2)-SM
            S3D=-PRINC(1)-SM
            AJ2D=DSQRT(1.5*(S1D*S1D+S2D*S2D+S3D*S3D))
            if(nmodm.eq.7.or.nmodm.eq.41.or.nmodm.eq.42.or.nmodm.eq.43
     1         .or.nmodm.eq.44.or.nmodm.eq.45)
     1      AJ2D=DSQRT(0.5*(S1D*S1D+S2D*S2D+S3D*S3D))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) (PRINC(K),K=1,3),AJ2D
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) (PRINC(K),K=1,3),AJ2D
C
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NLM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               SEFE=STRESS(7)
            ENDIF
            IF(NMODM.EQ.28)THEN
               CALL KOSI3(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL,LL)
               WRITE(IZLAZ,5001) NGR,NGS,NGT,(PLAST(K),K=LL,LL+5)
               CALL JEDNA1(STRESS,PLAST(LL),6)
               CALL GLAVN3(STRESS)
               IF(STRESS(7).GT.SEFE) THEN
                  NNLM=NLM
                  NNGR=NGR
                  NNGS=NGS
                  NNGT=NGT
                  SEFE=STRESS(7)
               ENDIF
               GO TO 10
            ENDIF
            WRITE(IZLAZ,5001) NGR,NGS,NGT,(PLAST(K),K=LL,LL+5)
            DO 300 MSL=1,NNSL
               MSLOF =MSL*6
               MSLOFF=MSLOF+5
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (PLAST(II),II=LL+MSLOF,LL+MSLOFF)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (PLAST(II),II=LL+MSLOF,LL+MSLOFF)
               IPRPL=0
               MSLOF =MXS*6+MSL*6
               MSLOFF=MSLOF+5
               DO 30 II=LL+MSLOF,LL+MSLOFF
                  IF(DABS(PLAST(II)).GT.1.0D-10) IPRPL=1
   30          CONTINUE
               IF(IPRPL.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) (PLAST(II),II=LL+MSLOF,LL+MSLOFF)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) (PLAST(II),II=LL+MSLOF,LL+MSLOFF)
               ENDIF
               SMT=(PLAST(LL)+PLAST(LL+1)+PLAST(LL+2))/3.0D0
C
               IF(NMODM.EQ.9) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) PLAST(LL+18),SMT,PLAST(LL+19),PLAST(LL+20)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) PLAST(LL+18),SMT,PLAST(LL+19),PLAST(LL+20)
               ENDIF
C
               IF(NMODM.EQ.25) THEN
                  MSLOF=MXS*12+5+MSL
      IF(ISRPS.EQ.0)THEN
         WRITE(IZLAZ,2004) PLAST(LL+MSLOF),SMT,PLAST(LL+MSLOF+MXS)/3.,
     +                     PLAST(LL+MSLOF+2*MXS)
         WRITE(IZLAZ,2254) PLAST(LL+MSLOF+3*MXS),PLAST(LL+MSLOF+4*MXS)
      ELSE
         WRITE(IZLAZ,2004) PLAST(LL+MSLOF),SMT,PLAST(LL+MSLOF+MXS)/3.,
     +                     PLAST(LL+MSLOF+2*MXS)
         WRITE(IZLAZ,2254) PLAST(LL+MSLOF+3*MXS),PLAST(LL+MSLOF+4*MXS)
      ENDIF
               ENDIF
C
               IF(NMODM.EQ.7) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18),PLAST(LL+20)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18),PLAST(LL+20)
               ENDIF
c               (rakic)               
               IF(NMODM.EQ.41) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18)
               ENDIF
c               (rakic)               
               IF(NMODM.EQ.42) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18)
               ENDIF
c               (rakic)               
               IF(NMODM.EQ.43) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18)
               ENDIF
c               (rakic)               
               IF(NMODM.EQ.44) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18)
               ENDIF
c               (rakic)               
               IF(NMODM.EQ.45) THEN
                  SMT3=3.0D0*SMT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) PLAST(LL+19),SMT3,PLAST(LL+18)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) PLAST(LL+19),SMT3,PLAST(LL+18)
               ENDIF
C
               IF(NMODM.EQ.22) THEN
      DIST=(PLAST(LL)-PLAST(LL+18))**2+(PLAST(LL+1)-PLAST(LL+19))**2
     1 +(PLAST(LL+2)-PLAST(LL+20))**2+2.*(PLAST(LL+3)-PLAST(LL+21))**2
     2 +2.*(PLAST(LL+4)-PLAST(LL+22))**2+2.*(PLAST(LL+5)-PLAST(LL+23))
     3 **2
                  DIST=SQRT(DIST)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2103) (PLAST(II),II=LL+18,LL+23)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) PLAST(LL+24),SMT,PLAST(LL+26),DIST
      ENDIF
C
      if(nlm.eq.1.and.ibtc.eq.1) then
      if(NMODM.EQ.9)
     1write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+18),SMT,PLAST(LL+19)
     1               ,aj2d,PLAST(LL+7),PLAST(LL+8)
      if(NMODM.EQ.7)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+20),plast(ll+1)
      endif
      if(NMODM.EQ.41)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+19),plast(ll+1)
      endif
      if(NMODM.EQ.42)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+19),plast(ll+1)
      endif
      if(NMODM.EQ.43)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+19),plast(ll+1)
      endif
C     stampanje tabele za model 44
      if(NMODM.EQ.44)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+19),plast(ll+1)
      endif
C     stampanje tabele za model 45
      if(NMODM.EQ.45)then
         ev=3.*PLAST(LL+18)
      write(imi,7000)PLAST(LL),PLAST(LL+6),PLAST(LL+19),SMT3,ev
     1               ,aj2d,PLAST(LL+19),plast(ll+1)
      endif
c      endif
 7000 format(8(1PE12.4))
c	endif
      endif
  300 CONTINUE
   10 CONTINUE
   20 CONTINUE
      if(kor.eq.ndt) close(imi,status='keep')
      IF(NNLM.GT.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NNLM,NNGR,NNGS,NNGT,SEFE
         IF(SMKOR.LT.SEFE) THEN
            SMKOR=SEFE
            NGRKOR=NGE
            NSMKOR=NNLM
         ENDIF
         IF(SMALL.LT.SEFE) THEN
            SMALL=SEFE
            NGRALL=NGE
            NSMALL=NNLM
            KSMALL=KOR
         ENDIF
      ENDIF
      RETURN
 5000 FORMAT(/I8)
 5001 FORMAT(3I3,6(1PE11.3))
C-----------------------------------------------------------------------
 2010 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2011 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2030 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2002 FORMAT(' UK.DEFO.',6(1PE11.3))
 2003 FORMAT(' PL.DEFO.',6(1PE11.3))
 2100 FORMAT(' GLAVNI NAPONI',3(1PE11.3),'     DSQRT(3*J2D)=',1PE11.3)
 2103 FORMAT(' NAPON-S.',6(1PE11.3))
 2004 FORMAT(' PK.P. ',1PE11.3,'  SMT ',1PE11.3,'  SR.PL.D. ',
     11PE11.3,'  OCR ', 1PE11.3)
 2254 FORMAT(' EPSPD ',1PE11.3,'  DAMAGE ',1PE11.3)
 2005 FORMAT(' PREK. PR  ',1PE11.3,'SMT ',1PE11.3,' SR.PL.DEF. ',
     11PE11.3,' DIST ',1PE11.3/)
 2014 FORMAT(' XTDT  ',1PE11.3,'  SMT ',1PE11.3,'  SR.PL.DEF.',
     11PE11.3,'   EL ', 1PE11.3)
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6010 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX')
 6011 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR')
 6030 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  IT =',I3,'  MAX EFF STR  =',1PE11.3//)
 6002 FORMAT(' TO.STRA.',6(1PE11.3))
 6003 FORMAT(' PL.STRA.',6(1PE11.3))
 6100 FORMAT(' PRINC. STRESSES',3(1PE11.3),'     DSQRT(3*J2D)=',1PE11.3)
 6004 FORMAT(' PC.P.',1PE11.3,' SMT ',1PE11.3,' VOL.P.D.',
     11PE11.3,' OCR ', 1PE11.3)
 6014 FORMAT(' AI1A  ',1PE11.3,'  SI1 ',1PE11.3,'  VOL.P.DEF.',
     11PE11.3,'  EL ', 1PE11.3)
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C-----------------------------------------------------------------------
      SUBROUTINE STA323(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 23
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 23
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*),STRESS(7),STRES(6)
      DIMENSION CORGT(3,NGS12,*)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN23'
      NNLM=0
      INDSN=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2140) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6140) NLM
         GO TO 20
         ENDIF
C
         IBTC=0
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
C        IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10.AND.IT.EQ.0) THEN
C           IT=1
C           CALL MATRTE
C        ELSE
            ISN=2
C        ENDIF
         IF(INDSN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2111) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6111) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2110) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6110) NGE
            ENDIF
         ENDIF
         INDSN=1
         IPC=IPRC(NLM)
         NMM=NLM
         NMMPR=MODPRO(NMODM)
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5100) NMM
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NLM-1)*NGS12*NMMPR+(IBTC-1)*NMMPR
            CALL JEDNA1(STRES,PLAST(LL),6)
            IF(ISN.EQ.2) CALL JEDNA1(STRESS,STRES,6)
            CALL GLAVN3(STRESS)
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NLM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               SEFE=STRESS(7)
            ENDIF
            WRITE(IZLAZ,5101) NGR,NGS,NGT,(PLAST(K),K=LL,LL+5)
      DIST=(PLAST(LL)-PLAST(LL+18))**2+(PLAST(LL+1)-PLAST(LL+19))**2
     1 +(PLAST(LL+2)-PLAST(LL+20))**2+2.*(PLAST(LL+3)-PLAST(LL+21))**2
     2 +2.*(PLAST(LL+4)-PLAST(LL+22))**2+2.*(PLAST(LL+5)-PLAST(LL+23))
     3 **2
      DIST=SQRT(DIST)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2102) (PLAST(II),II=LL+6,LL+11)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6102) (PLAST(II),II=LL+6,LL+11)
            IPRPL=0
            DO 30 II=LL+12,LL+17
C               IF(PLAST(II).GT.1.0D-10) IPRPL=1
   30       CONTINUE
C            IF(IPRPL.EQ.1) THEN
      IF(ISRPS.EQ.0) 
     1WRITE(IZLAZ,2103) (PLAST(II),II=LL+12,LL+17)
      IF(ISRPS.EQ.1) 
     1WRITE(IZLAZ,6103) (PLAST(II),II=LL+12,LL+17)
      SMT=(PLAST(LL)+PLAST(LL+1)+PLAST(LL+2))/3.0D0
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2104) (PLAST(II),II=LL+18,LL+23)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2105) PLAST(LL+24),PLAST(LL+25),PLAST(LL+26),DIST
C
   10    CONTINUE
   20 CONTINUE
      IF(NNLM.GT.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2130) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6130) NNLM,NNGR,NNGS,NNGT,SEFE
      ENDIF
      RETURN
 5100 FORMAT(/I8)
 5101 FORMAT(3I3,6(1PE11.3))
C-----------------------------------------------------------------------
 2110 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2111 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2130 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2102 FORMAT(' UK.DEFO.',6(1PE11.3))
 2103 FORMAT(' PL.DEFO.',6(1PE11.3))
 2104 FORMAT(' NAPON-S.',6(1PE11.3))
 2105 FORMAT(' SR.VP.DEF. ',1PE11.3,' X-CAP ',1PE11.3,' EL-CAP ',
     11PE11.3,' DIST ',1PE11.3/)
 2140 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6110 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX')
 6111 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR')
 6130 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  IT =',I3,'  MAX EFF STR  =',1PE11.3//)
 6102 FORMAT(' TO.STRA.',6(1PE11.3))
 6103 FORMAT(' PL.STRA.',6(1PE11.3))
 6140 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE STA330(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 30
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 30
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*),STRESS(7),STRES(6)
      DIMENSION CORGT(3,NGS12,*)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN30'
      NNLM=0
      INDSN=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         IBTC=0
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
C        IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10.AND.IT.EQ.0) THEN
C           IT=1
C           CALL MATRTE
C        ELSE
            ISN=2
C        ENDIF
         IF(INDSN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NGE
            ENDIF
         ENDIF
         INDSN=1
         IPC=IPRC(NLM)
         NMM=NLM
         NMMPR=MODPRO(NMODM)
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NLM-1)*NGS12*NMMPR+(IBTC-1)*NMMPR
            CALL JEDNA1(STRES,PLAST(LL),6)
            IF(ISN.EQ.2) CALL JEDNA1(STRESS,STRES,6)
            CALL GLAVN3(STRESS)
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NLM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               SEFE=STRESS(7)
            ENDIF
            WRITE(IZLAZ,5001) NGR,NGS,NGT,(PLAST(K),K=LL+3,LL+8)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (PLAST(II),II=LL+9,LL+14)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (PLAST(II),II=LL+9,LL+14)
            IPRPL=0
      INDSH = PLAST(LL+2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) PLAST(LL),PLAST(LL+1),PLAST(LL+25),INDSH
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) PLAST(LL),PLAST(LL+1),PLAST(LL+25),INDSH
C
   10    CONTINUE
   20 CONTINUE
      IF(NNLM.GT.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NNLM,NNGR,NNGS,NNGT,SEFE
      ENDIF
      RETURN
 5000 FORMAT(/I8)
 5001 FORMAT(3I3,6(1PE11.3))
C-----------------------------------------------------------------------
 2010 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2011 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2030 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2002 FORMAT(' UK.DEFO.',6(1PE11.3))
 2003 FORMAT(' DUZ PROKL.',1PE11.3,' ATH NAPON',1PE11.3,
     1 ' BRZ KLIZ.',1PE11.3,' IND KLIZ.',I2/)
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6010 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 2/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS     STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-XY',2X,'STRESS-ZZ',2X,'STRESS-1',2X,'STRESS-2  ALFA'/
     2'                     TOTAL STRAINS',15X,'EFF STRESS' /
     3' REL SLIP LENGTH  ADH STRESS  STRAIN VELOC',3X,'IND OF SLIP')
 6011 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 2/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS     STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-RS',2X,'STRESS-ZZ',2X,'STRESS-1',2X,'STRESS-2  ALFA'/
     2'                     TOTAL STRAINS',15X,'EFF STRESS' /
     3' REL SLIP LENGTH  ADH STRESS  STRAIN VELOC',3X,'IND OF SLIP')
 6030 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  LAYER =',I3,'  MAX EFF STR  =',1PE11.3//)
 6002 FORMAT(' TO.STRA.',6(1PE11.3))
 6003 FORMAT(' SLIP LENGTH ',1PE11.3,' ADH STRESS ',1PE11.3,
     1 ' STRAIN VELOC ',1PE11.3,' IND OF SLIP ',I2/)
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE STA331(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS IN 3/D ELEMENTS FOR MATERIAL MODEL 31
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA ZA MATERIJALNI MODEL 31
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*),STRESS(7),STRES(6)
      DIMENSION CORGT(3,NGS12,*)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STA331'
      NNLM=0
      INDSN=0
      IT=0
      SEFE=0.D0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         IBTC=0
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
C        IF(ISN.EQ.0.AND.DABS(BETA).GT.1.0D-10.AND.IT.EQ.0) THEN
C           IT=1
C           CALL MATRTE
C        ELSE
            ISN=2
C        ENDIF
         IF(INDSN.EQ.0) THEN
            IF(ISN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) NGE
            ENDIF
            IF(ISN.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NGE
            ENDIF
         ENDIF
         INDSN=1
         IPC=IPRC(NLM)
         NMM=NLM
         NMMPR=MODPRO(NMODM)
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NLM-1)*NGS12*NMMPR+(IBTC-1)*NMMPR
            CALL JEDNA1(STRES,PLAST(LL),6)
            IF(ISN.EQ.2) CALL JEDNA1(STRESS,STRES,6)
            CALL GLAVN3(STRESS)
            IF(STRESS(7).GT.SEFE) THEN
               NNLM=NLM
               NNGR=NGR
               NNGS=NGS
               NNGT=NGT
               SEFE=STRESS(7)
            ENDIF
            WRITE(IZLAZ,5001) NGR,NGS,NGT,(PLAST(K),K=LL,LL+5)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (PLAST(II),II=LL+6,LL+11)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (PLAST(II),II=LL+6,LL+11)
            IPRPL=0
            DO 30 II=LL+12,LL+17
C               IF(PLAST(II).GT.1.0D-10) IPRPL=1
   30       CONTINUE
C            IF(IPRPL.EQ.1) THEN
      SMT=(PLAST(LL)+PLAST(LL+1)+PLAST(LL+2))/3.0D0
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) PLAST(LL+12),PLAST(LL+13),PLAST(LL+14)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) PLAST(LL+12),PLAST(LL+13),PLAST(LL+14)
C
   10    CONTINUE
   20 CONTINUE
      IF(NNLM.GT.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NNLM,NNGR,NNGS,NNGT,SEFE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NNLM,NNGR,NNGS,NNGT,SEFE
      ENDIF
      RETURN
C
 5000 FORMAT(/I8)
 5001 FORMAT(3I3,6(1PE11.3))
C-----------------------------------------------------------------------
 2010 FORMAT(/'1'///' GLOBALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRU
     1PE ELEMENATA',I6/1X,68('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2011 FORMAT(/'1'///' LOKALNI NAPONI IZOPARAMETARSKIH 3/D ELEMENATA GRUP
     1E ELEMENATA',I6/1X,67('-')//
     1' ELEMENT /'/' IR IS IT   NAPON-XX',2X,' NAPON-YY',2X,
     1' NAPON-ZZ',2X,' NAPON-XY',2X,' NAPON-YZ',2X,' NAPON-ZX')
 2030 FORMAT(//' MAKSIMALNI EFEKTIVNI NAPON:'/' ELEMENT =',I6,'  IR =',
     1I3,'  IS =',I3,'  IT =',I3,'  MAX.EFE.NAP. =',1PE11.3//)
 2002 FORMAT(' UK.DEFO.',6(1PE11.3))
 2004 FORMAT(' NAPON SIGS',1PE11.3,' LAMDA-SARKOMER',1PE11.3,
     1' LAMDA-TENDON',1PE11.3)
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6010 FORMAT(/'1'///' GLOBAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEM
     1ENTS FOR GROUP',I6/1X,70('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-XX',2X,'STRESS-YY',2X,
     1'STRESS-ZZ',2X,'STRESS-XY',2X,'STRESS-YZ',2X,'STRESS-ZX')
 6011 FORMAT(/'1'///' LOCAL STRESS COMPONENTS OF ISOPARAMETRIC 3/D ELEME
     1NTS FOR GROUP',I6/1X,69('-')//
     1' ELEMENT /'/' IR IS IT  STRESS-RR',2X,'STRESS-SS',2X,
     1'STRESS-TT',2X,'STRESS-RS',2X,'STRESS-ST',2X,'STRESS-TR')
 6030 FORMAT(//' MAXIMUM EFFECTIVE STRESS:'/' ELEMENT =',I6,'  IR =',I3,
     1'  IS =',I3,'  IT =',I3,'  MAX EFF STR  =',1PE11.3//)
 6002 FORMAT(' TO.STRA.',6(1PE11.3))
 6004 FORMAT(' STRESS SIGS',1PE11.3,' STRECH-SARCOMERE',1PE11.3,
     1' STRECH-TENDON',1PE11.3)
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE STAG31(TAU,AU,ISNA,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS IN UNIVERSAL FILE - ELASTI
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA U UNIVERZALNI FILE - ELASTIC
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /NASLOV/ NASLOV
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
C
      DIMENSION MCVEL(*),ISNA(*),SR(6)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAG31'
      NNCVE=NCVE
      IF(NCVE.LT.20) NCVE=8
      IF(NCVE.EQ.21) NCVE=20
C     STRUKTURNA ANALIZA = 1
      IMOTY=1
C     STACIONARAN = 1; NESTACIONARAN = 4   
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN 
        IANTY=4
        IFAT1=2
        FATY8=VREME
      ENDIF
C     SIMETRICAN TENZOR = 4
      IDACH=4
C     NAPONI = 2
      ISDTY=2
C     PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4 
      IDATY=2
C     BROJ PODATAKA = 6
      NDVPN=6
      IEXP=1
      NNODS=NCVE
      IF(ISTNA.GT.0) NNODS=-NNODS
      NVPN=6
      IND1=-1
      INA1=57
      INDSN=0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 20
C
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
         IF(INDSN.EQ.0) THEN
      WRITE(IGRAF,5000) IND1
      WRITE(IGRAF,5000) INA1
      WRITE(IGRAF,5003) NASLOV
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2004) KOR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6004) KOR
      WRITE(IGRAF,1000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR,KOR
      WRITE(IGRAF,5001) FATY8
         ENDIF
         INDSN=1
         NUM=NLM+ISUMEL
         IF(ICVEL.EQ.1) NUM=MCVEL(NLM)
         WRITE(IGRAF,1000) NUM,IEXP,NNODS,NVPN
         C=1./NCVE
         SE=0.D0
         CALL CLEAR(SR,6)
         IF(NCVE.EQ.8.AND.NGAUSX.EQ.2) THEN
            CALL WRIT31(TAU,NLM,2,2,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,2,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,2,SR,SE)
C
            CALL WRIT31(TAU,NLM,2,2,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,1,SR,SE)
         ENDIF
         IF(NCVE.EQ.8.AND.NGAUSX.EQ.3) THEN
            CALL WRIT31(TAU,NLM,3,3,3,SR,SE)
            CALL WRIT31(TAU,NLM,1,3,3,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,3,SR,SE)
            CALL WRIT31(TAU,NLM,3,1,3,SR,SE)
C
            CALL WRIT31(TAU,NLM,3,3,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,3,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,3,1,1,SR,SE)
         ENDIF
         IF(NCVE.EQ.20.AND.NGAUSX.EQ.3) THEN
            CALL WRIT31(TAU,NLM,3,3,3,SR,SE)
            CALL WRIT31(TAU,NLM,2,3,3,SR,SE)
            CALL WRIT31(TAU,NLM,1,3,3,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,3,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,3,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,3,SR,SE)
            CALL WRIT31(TAU,NLM,3,1,3,SR,SE)
            CALL WRIT31(TAU,NLM,3,2,3,SR,SE)
C
            CALL WRIT31(TAU,NLM,3,3,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,3,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,2,SR,SE)
            CALL WRIT31(TAU,NLM,3,1,2,SR,SE)
C
            CALL WRIT31(TAU,NLM,3,3,1,SR,SE)
            CALL WRIT31(TAU,NLM,2,3,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,3,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,3,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,3,2,1,SR,SE)
         ENDIF
         IF(NCVE.EQ.20.AND.NGAUSX.EQ.2) THEN
            CALL WRIT31(TAU,NLM,2,2,2,SR,SE)
            CALL WRIT32(TAU,NLM,2,2,2,1,2,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,2,SR,SE)
            CALL WRIT32(TAU,NLM,1,2,2,1,1,2,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,2,SR,SE)
            CALL WRIT32(TAU,NLM,1,1,2,2,1,2,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,2,SR,SE)
            CALL WRIT32(TAU,NLM,2,1,2,2,2,2,SR,SE)
C
            CALL WRIT32(TAU,NLM,2,2,2,2,2,1,SR,SE)
            CALL WRIT32(TAU,NLM,1,2,2,1,2,1,SR,SE)
            CALL WRIT32(TAU,NLM,1,1,2,1,1,1,SR,SE)
            CALL WRIT32(TAU,NLM,2,1,2,2,1,1,SR,SE)
C
            CALL WRIT31(TAU,NLM,2,2,1,SR,SE)
            CALL WRIT32(TAU,NLM,2,2,1,1,2,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,2,1,SR,SE)
            CALL WRIT32(TAU,NLM,1,2,1,1,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,1,1,1,SR,SE)
            CALL WRIT32(TAU,NLM,1,1,1,2,1,1,SR,SE)
            CALL WRIT31(TAU,NLM,2,1,1,SR,SE)
            CALL WRIT32(TAU,NLM,2,1,1,2,2,1,SR,SE)
         ENDIF
         IF(ISTNA.GT.0) CALL STGRNA(SR,C,ISTNA,IGRAF,6)
   20 CONTINUE
      IF(INDSN.EQ.1) WRITE(IGRAF,5000) IND1
      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
 1000 FORMAT(8I10)
 5000 FORMAT(I6)
 5001 FORMAT(6(1PE13.5))
 5003 FORMAT(A80)
C-----------------------------------------------------------------------
 2004 FORMAT('NAPONI 3/D ELEMENATA'/
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6004 FORMAT('STRESS OF 3/D ELEMENTS'/
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRIT31(TAU,NLM,I,J,K,SRED,SEFE1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS IN UNIVERSAL FILE - ELASTI
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA U UNIVERZALNI FILE - ELASTIC
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      DIMENSION TAU(6,NGS12,*),SIGM(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WRIT31'
      IJ=(I-1)*NGAUSY*NGAUSZ+(J-1)*NGAUSZ+K
      SRBA=1.
C      SRBA=1.E7
      SIGM(1)=TAU(1,IJ,NLM)*SRBA
      SIGM(2)=TAU(4,IJ,NLM)*SRBA
      SIGM(3)=TAU(2,IJ,NLM)*SRBA
      SIGM(4)=TAU(6,IJ,NLM)*SRBA
      SIGM(5)=TAU(5,IJ,NLM)*SRBA
      SIGM(6)=TAU(3,IJ,NLM)*SRBA
      IF(ISTNA.EQ.0) THEN
         WRITE(IGRAF,5000) SIGM
      ELSE
         CALL SREDNA(SRED,SIGM,SEFE1,ISTNA,1)
      ENDIF
      RETURN
 5000 FORMAT(6(1PE13.5))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRIT32(TAU,NLM,I,J,K,L,M,N,SRED,SEFE1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS WITH 20 NODAL POINT
CE.       IN UNIVERSAL FILE, ORDER GAUSS INTEGRATION 2 X 2 X 2 - ELASTIC
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA SA 20 COROVA U UNIVERZALNI
CS.       FILE, RED GAUSOVE INTEGRACIJE 2 X 2 X 2 - ELASTICNOST
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      DIMENSION TAU(6,NGS12,*),SIGM(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WRIT32'
      IJ=(I-1)*NGAUSY*NGAUSZ+(J-1)*NGAUSZ+K
      MN=(L-1)*NGAUSY*NGAUSZ+(M-1)*NGAUSZ+N
      SIGM(1)=(TAU(1,IJ,NLM)+TAU(1,MN,NLM))/2.D0
      SIGM(2)=(TAU(4,IJ,NLM)+TAU(4,MN,NLM))/2.D0
      SIGM(3)=(TAU(2,IJ,NLM)+TAU(2,MN,NLM))/2.D0
      SIGM(4)=(TAU(6,IJ,NLM)+TAU(6,MN,NLM))/2.D0
      SIGM(5)=(TAU(5,IJ,NLM)+TAU(5,MN,NLM))/2.D0
      SIGM(6)=(TAU(3,IJ,NLM)+TAU(3,MN,NLM))/2.D0
      IF(ISTNA.EQ.0) THEN
         WRITE(IGRAF,5000) SIGM
      ELSE
         CALL SREDNA(SRED,SIGM,SEFE1,ISTNA,1)
      ENDIF
      RETURN
 5000 FORMAT(6(1PE13.5))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAG35(PLAST,AU,ISNA,MCVEL,ICVEL,ISDTY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS IN UNIVERSAL FILE - PLASTI
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA U UNIVERZALNI FILE - PLASTIC
C . 
C ......................................................................
C 
      CHARACTER*250 NASLOV    
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /NASLOV/ NASLOV
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
C
      DIMENSION MCVEL(*),PLAST(*),ISNA(*),SR(6)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAG35'
      NNCVE=NCVE
      IF(NCVE.LT.20) NCVE=8
      IF(NCVE.EQ.21) NCVE=20
C     STRUKTURNA ANALIZA = 1
      IMOTY=1
C     STACIONARAN = 1; NESTACIONARAN = 4   
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN 
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
C     SIMETRICAN TENZOR = 4
      IDACH=4
C     NAPONI = 2
C     ISDTY=2
C     UKUPNE DEFORMACIJE = 3
C     ISDTY=3
C     PLASTICNE DEFORMACIJE = 0
C     ISDTY=0
C     PARAMETRI CAM-CLAY MODELA = 99
C     ISDTY=99
C     DEFORMACIJE PUZANJA = 1
C     ISDTY=1
C     PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4 
      IDATY=2
C     BROJ PODATAKA = 6
      NDVPN=6
      IEXP=1
      NNODS=NCVE
      IF(ISTNA.GT.0) NNODS=-NNODS
C     BROJ PODATAKA PO CVORU
      NVPN=6
      IND1=-1
      INA1=57
C
      NPROS=MODPRO( NMODM )
      NGXYZ=NGS12*NPROS
      IF(ISDTY.EQ.2) THEN
         KK=1
      ELSEIF(ISDTY.EQ.3) THEN
         KK=7
      ELSEIF(ISDTY.EQ.0) THEN
         KK=13
      ELSEIF(ISDTY.EQ.1) THEN
         KK=13
         IF(NMODM.EQ.16.OR.NMODM.EQ.19) KK=25
      ELSEIF(ISDTY.EQ.99) THEN
         KK=26
         IF(NMODM.EQ.15.OR.NMODM.EQ.18) KK=32
         IF(NMODM.EQ.16.OR.NMODM.EQ.19) KK=45
         IF(NMODM.EQ.61) KK=13
      ENDIF
      INDSN=0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 20
C
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
         IF(INDSN.EQ.0) THEN
      WRITE(IGRAF,5000) IND1
      WRITE(IGRAF,5000) INA1
      WRITE(IGRAF,5003) NASLOV
      IF(ISDTY.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6001)
      ELSEIF(ISDTY.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2002)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6002)
      ELSEIF(ISDTY.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2003)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6003)
      ELSEIF(ISDTY.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6005)
      ELSEIF(ISDTY.EQ.99) THEN
         if(nmodm.eq.61) then
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2007)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6007)
         else
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2006)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6006)
         endif
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2004) KOR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6004) KOR
      WRITE(IGRAF,1000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR,KOR
      WRITE(IGRAF,5001) FATY8
         ENDIF
         INDSN=1
         LL=KK+(NLM-1)*NGXYZ
         NUM=NLM+ISUMEL
         IF(ICVEL.EQ.1) NUM=MCVEL(NLM)
         WRITE(IGRAF,1000) NUM,IEXP,NNODS,NVPN
         C=1./NCVE
         SE=0.D0
         CALL CLEAR(SR,6)
         IF(NCVE.EQ.8.AND.NGAUSX.EQ.2) THEN
            CALL WRIT35(PLAST(LL+7*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+3*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+  NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+5*NPROS),SR,SE)
C
            CALL WRIT35(PLAST(LL+6*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+2*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL        ),SR,SE)
            CALL WRIT35(PLAST(LL+4*NPROS),SR,SE)
         ENDIF
         IF(NCVE.EQ.8.AND.NGAUSX.EQ.3) THEN
            CALL WRIT35(PLAST(LL+26*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 8*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 2*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+20*NPROS),SR,SE)
C
            CALL WRIT35(PLAST(LL+24*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 6*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL         ),SR,SE)
            CALL WRIT35(PLAST(LL+18*NPROS),SR,SE)
         ENDIF
         IF(NCVE.EQ.20.AND.NGAUSX.EQ.3) THEN
            CALL WRIT35(PLAST(LL+26*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+17*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 8*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 5*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 2*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+11*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+20*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+23*NPROS),SR,SE)
C
            CALL WRIT35(PLAST(LL+25*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 7*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+   NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+19*NPROS),SR,SE)
C
            CALL WRIT35(PLAST(LL+24*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+15*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 6*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+ 3*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL         ),SR,SE)
            CALL WRIT35(PLAST(LL+ 9*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+18*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+21*NPROS),SR,SE)
         ENDIF
         IF(NCVE.EQ.20.AND.NGAUSX.EQ.2) THEN
            CALL WRIT35(PLAST(LL+7*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+7*NPROS),PLAST(LL+3*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+3*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+3*NPROS),PLAST(LL+  NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+  NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+  NPROS),PLAST(LL+5*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+5*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+5*NPROS),PLAST(LL+7*NPROS),SR,SE)
C
            CALL WRIT36(PLAST(LL+7*NPROS),PLAST(LL+6*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+3*NPROS),PLAST(LL+2*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+  NPROS),PLAST(LL        ),SR,SE)
            CALL WRIT36(PLAST(LL+5*NPROS),PLAST(LL+4*NPROS),SR,SE)
C
            CALL WRIT35(PLAST(LL+6*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+6*NPROS),PLAST(LL+2*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+2*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+2*NPROS),PLAST(LL        ),SR,SE)
            CALL WRIT35(PLAST(LL        ),SR,SE)
            CALL WRIT36(PLAST(LL        ),PLAST(LL+4*NPROS),SR,SE)
            CALL WRIT35(PLAST(LL+4*NPROS),SR,SE)
            CALL WRIT36(PLAST(LL+4*NPROS),PLAST(LL+6*NPROS),SR,SE)
         ENDIF
         IF(ISTNA.GT.0) CALL STGRNA(SR,C,ISTNA,IGRAF,6)
   20 CONTINUE
      IF(INDSN.EQ.1) WRITE(IGRAF,5000) IND1
      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
 1000 FORMAT(8I10)
 5000 FORMAT(I6)
 5001 FORMAT(6(1PE13.5))
 5003 FORMAT(A80)
C-----------------------------------------------------------------------
 2001 FORMAT('NAPONI 3/D ELEMENATA')
 2002 FORMAT('UKUPNE DEFORMACIJE 3/D ELEMENATA')
 2003 FORMAT('PLASTICNE DEFORMACIJE 3/D ELEMENATA')
 2005 FORMAT('DEFORMACIJE PUZANJA 3/D ELEMENATA')
 2006 FORMAT('EFEKTIVNE DEFORMACIJE PLASTICNE,PUZANJA')
 2007 FORMAT('OSTECENJE')
 2004 FORMAT(
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6001 FORMAT('STRESS OF 3/D ELEMENTS')
 6002 FORMAT('TOTAL STRAIN OF 3/D ELEMENTS')
 6003 FORMAT('PLASTIC STRAIN OF 3/D ELEMENTS')
 6005 FORMAT('CREEP STRAIN OF 3/D ELEMENTS')
 6006 FORMAT('EFFECTIVE STRAIN PLASTIC,CREEP OF SHELL ELEMENT')
 6007 FORMAT('DAMAGE')
 6004 FORMAT(
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRIT35(PLAST1,SRED,SEFE1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................

C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS IN UNIVERSAL FILE - PLASTI
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA U UNIVERZALNI FILE - PLASTIC
C . 
C ......................................................................
C 
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      DIMENSION PLAST1(*),SIGM(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WRIT35'
      SIGM(1)=PLAST1(1)
      SIGM(2)=PLAST1(4)
      SIGM(3)=PLAST1(2)
      SIGM(4)=PLAST1(6)
      SIGM(5)=PLAST1(5)
      SIGM(6)=PLAST1(3)
      IF(ISTNA.EQ.0) THEN
         WRITE(IGRAF,5000) SIGM
      ELSE
         CALL SREDNA(SRED,SIGM,SEFE1,ISTNA,1)
      ENDIF
      RETURN
 5000 FORMAT(6(1PE13.5))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRIT36(PLAST1,PLAST2,SRED,SEFE1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS WITH 20 NODAL POINT
CE.       IN UNIVERSAL FILE, ORDER GAUSS INTEGRATION 2 X 2 X 2 - PLASTIC
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA SA 20 COROVA U UNIVERZALNI
CS.       FILE, RED GAUSOVE INTEGRACIJE 2 X 2 X 2 - PLASTICNOST
C . 
C ......................................................................
C 
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      DIMENSION PLAST1(*),PLAST2(*),SIGM(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WRIT36'
      SIGM(1)=(PLAST1(1)+PLAST2(1))/2.D0
      SIGM(2)=(PLAST1(4)+PLAST2(4))/2.D0
      SIGM(3)=(PLAST1(2)+PLAST2(2))/2.D0
      SIGM(4)=(PLAST1(6)+PLAST2(6))/2.D0
      SIGM(5)=(PLAST1(5)+PLAST2(5))/2.D0
      SIGM(6)=(PLAST1(3)+PLAST2(3))/2.D0
      IF(ISTNA.EQ.0) THEN
         WRITE(IGRAF,5000) SIGM
      ELSE
         CALL SREDNA(SRED,SIGM,SEFE1,ISTNA,1)
      ENDIF
      RETURN
 5000 FORMAT(6(1PE13.5))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STSIL3(AU,ESILA,IPGC,MCVEL,NEL,NCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRINTOUT FORCES IN 3D ELEMENTS NODES
CS.   P R O G R A M
CS.      ZA STAMPANJE SILA U CVOROVIMA 3D ELEMENATA 
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION AU(*)
      REAL AU
C
      DIMENSION ESILA(ND,*),IPGC(*),MCVEL(*),NEL(NE,*),NCVEL(*),SILA(3)
C
      IF(IDEBUG.GT.0) PRINT *, ' STSIL3'
C
      INDSN=0
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLM
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLM
         GO TO 20
         ENDIF
C
         IPG=IPGC(NLM)
         IF(IPG.EQ.0) GO TO 20
         IF(INDSN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2012) NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6012) NGE
         ENDIF
C
         INDSN=1
         NMM=NLM
         IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
         WRITE(IZLAZ,5000) NMM
C
C        C V O R O V I
C
         DO 12 N=1,NCVE
            NN=NEL(NLM,N)
            IF(NN.EQ.0) GO TO 12
            IF(ICVEL.EQ.1) NN=NCVEL(NN)
            N1=(N-1)*3
            DO 13 I=1,3
               SILA(I)=ESILA(N1+I,IPG)
   13       CONTINUE
            WRITE(IZLAZ,5003) NN,(SILA(K),K=1,3)
   12    CONTINUE
C
   20 CONTINUE
      RETURN
C
 5000 FORMAT(/I8)
 5003 FORMAT(I6,3X,3(1PE11.3))
C-----------------------------------------------------------------------
 2012 FORMAT(/'1'///' GLOBALNE CVORNE SILE 3-D ELEMENATA GRUPE ELEMENATA
     1',I6/1X,56('-')//
     1' ELEMENT /'/' C V O R     SILA-X ',2X,'  SILA-Y ',2X,'  SILA-Z ')
 2040 FORMAT(/' NESTAO JE ELEMENT -',I5)
C-----------------------------------------------------------------------
 6012 FORMAT(/'1'///' GLOBAL NODAL FORCES COMPONENTS OF 3-D ELEMENTS FOR
     1 GROUP',I6/1X,62('-')//
     1' ELEMENT /'/' N O D E    FORCE-X ',2X,' FORCE-Y ',2X,' FORCE-Z ')
 6040 FORMAT(/' DEATH ELEMENT -',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STG325(PLAST,AU,NMAT,NSLOJ,ISNA,MCVEL,ICVEL,ISDTY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT STRESS FOR 3/D ELEMENTS IN UNIVERSAL FILE - PLASTI
CS.   P R O G R A M
CS.       ZA STAMPANJE NAPONA 3/D ELEMENATA U UNIVERZALNI FILE - PLASTIC
C . 
C ......................................................................
C 
      CHARACTER*250 NASLOV    
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /NASLOV/ NASLOV
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      DIMENSION MCVEL(*),PLAST(*),NMAT(*),NSLOJ(*),ISNA(*),SR(6)
      DIMENSION AU(*)
      DIMENSION VEK(6),ISHIF8(8)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
      DATA ISHIF8/7,3,1,5,6,2,0,4/
C
      IF(IDEBUG.GT.0) PRINT *, ' STAG35'
      NNCVE=NCVE
      IF(NCVE.LT.20) NCVE=8
      IF(NCVE.EQ.21) NCVE=20
C     STRUKTURNA ANALIZA = 1
      IMOTY=1
C     STACIONARAN = 1; NESTACIONARAN = 4   
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN 
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
C     SIMETRICAN TENZOR = 4
      IDACH=4
C     NAPONI = 2
C     ISDTY=2
C     UKUPNE DEFORMACIJE = 3
C     ISDTY=3
C     PLASTICNE DEFORMACIJE = 0
C     ISDTY=0
C     PARAMETRI CAM-CLAY MODELA = 99
C     ISDTY=99
C     DEFORMACIJE PUZANJA = 1
C     ISDTY=1
C     PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4 
      IDATY=2
C     BROJ PODATAKA = 6
      IEXP=1
      NDVPN=6
      NNODS=NCVE
      IF(ISTNA.GT.0) NNODS=-NNODS
C     BROJ PODATAKA PO CVORU
      NVPN=6
      IND1=-1
      INA1=57
C
      NPROS=MODPRO( NMODM )
      NGXYZ=NGS12*NPROS
C
C      IF(ISDTY.EQ.2) THEN
C         KK=1
C      ELSEIF(ISDTY.EQ.3) THEN
C         KK=7
C      ELSEIF(ISDTY.EQ.0) THEN
C         KK=13
C      ELSEIF(ISDTY.EQ.99) THEN
C         KK=19
C      ELSEIF(ISDTY.EQ.1) THEN
C         KK=25
C      ENDIF
      INDSN=0
      DO 20 NLM=1,NE
C
      MST =NMAT(NLM)
      NNSL=1
      IF(MSET.GT.0) NNSL=NSLOJ(MST)
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 20
C
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
         IF(INDSN.EQ.0) THEN
      WRITE(IGRAF,5000) IND1
      WRITE(IGRAF,5000) INA1
      WRITE(IGRAF,5003) NASLOV
      IF(ISDTY.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6001)
      ELSEIF(ISDTY.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2002)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6002)
      ELSEIF(ISDTY.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2003)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6003)
      ELSEIF(ISDTY.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6005)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2004) KOR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6004) KOR
      WRITE(IGRAF,1000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR,KOR
      WRITE(IGRAF,5001) FATY8
         ENDIF
         INDSN=1
C         LL=KK+(NLM-1)*NGXYZ
         LL=1+(NLM-1)*NGXYZ
         NUM=NLM+ISUMEL
         IF(ICVEL.EQ.1) NUM=MCVEL(NLM)
         WRITE(IGRAF,1000) NUM,IEXP,NNODS,NVPN
         C=1./NCVE
         SE=0.D0
         CALL CLEAR(SR,6)
         IF(NNODS.EQ.8.AND.NGAUSX.EQ.2) THEN
            IF(ISDTY.EQ.2)THEN
              CALL WRIT35(PLAST(LL+7*NPROS),SR,SE)
              CALL WRIT35(PLAST(LL+3*NPROS),SR,SE)
              CALL WRIT35(PLAST(LL+  NPROS),SR,SE)
              CALL WRIT35(PLAST(LL+5*NPROS),SR,SE)
C
              CALL WRIT35(PLAST(LL+6*NPROS),SR,SE)
              CALL WRIT35(PLAST(LL+2*NPROS),SR,SE)
              CALL WRIT35(PLAST(LL        ),SR,SE)
              CALL WRIT35(PLAST(LL+4*NPROS),SR,SE)
            ENDIF
            IF(ISDTY.EQ.3)THEN
              DO 200 NC=1,8
                CALL CLEAR(VEK,6)
                DO 150 MSL=1,NNSL
                  LS=LL+ISHIF8(NC)*NPROS+6*MXS+6*(MSL-1)+5
                  VEK(1)=VEK(1)+PLAST(LS+1)
                  VEK(2)=VEK(2)+PLAST(LS+4)
                  VEK(3)=VEK(3)+PLAST(LS+2)
                  VEK(4)=VEK(4)+PLAST(LS+6)
                  VEK(5)=VEK(5)+PLAST(LS+5)
                  VEK(6)=VEK(6)+PLAST(LS+3)
150             CONTINUE
                WRITE(IGRAF,1500) VEK
200           CONTINUE
            ENDIF
            IF(ISDTY.EQ.0)THEN
              DO 100 NC=1,8
                CALL CLEAR(VEK,6)
                DO 250 MSL=2,NNSL
C  ZAPREMINSKA DEF
C                  LS=LL+ISHIF8(NC)*NPROS+13*MXS+MSL+5
C  KOMPONENTE PLAST.DEF.
                  LS=LL+ISHIF8(NC)*NPROS+6*MXS+(MSL-1)*6+6
                  DEFQP=EQDEF3(PLAST(LS))
                  VEK(1)=VEK(1)+DEFQP
C++++ KRPEZ ZA 5 FAMILIJA PRSLINA SAMO ZA DIVCA!!!
                  IF(MSL.GT.1.AND.MSL.LT.7) VEK(MSL)=DEFQP
C++++
CC                  LS=LL+ISHIF8(NC)*NPROS+17*MXS+MSL+5
CC                  VEK(1)=VEK(1)+PLAST(LS)/3.
CCC++++ KRPEZ ZA 5 FAMILIJA PRSLINA SAMO ZA DIVCA!!!
CC                  IF(MSL.GT.1.AND.MSL.LT.7) VEK(MSL)=PLAST(LS)/3.
CCC++++
250             CONTINUE
                WRITE(IGRAF,1500) VEK
100           CONTINUE
            ENDIF
         ENDIF
C         IF(ISTNA.GT.0) CALL STGRNA(SR,C,ISTNA,IGRAF,6)
   20 CONTINUE
      IF(INDSN.EQ.1) WRITE(IGRAF,5000) IND1
C      WRITE(IGRAF,5000) IND1
      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
 1000 FORMAT(8I10)
 1500 FORMAT(6(1PE13.5))
 5000 FORMAT(I6)
 5001 FORMAT(6(1PE13.5))
 5003 FORMAT(A80)
C-----------------------------------------------------------------------
 2001 FORMAT('NAPONI 3/D ELEMENATA')
 2002 FORMAT('UKUPNE DEFORMACIJE 3/D ELEMENATA')
 2003 FORMAT('PLASTICNE DEFORMACIJE 3/D ELEMENATA')
 2005 FORMAT('DEFORMACIJE PUZANJA 3/D ELEMENATA')
 2004 FORMAT(
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6001 FORMAT('STRESS OF 3/D ELEMENTS')
 6002 FORMAT('TOTAL STRAIN OF 3/D ELEMENTS')
 6003 FORMAT('PLASTIC STRAIN OF 3/D ELEMENTS')
 6005 FORMAT('CREEP STRAIN OF 3/D ELEMENTS')
 6004 FORMAT(
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      FUNCTION EQDEF3(P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  EKVIVALENTNA PLASTICNA DEFORMACIJA
C
      DIMENSION P(*)
      EQDEF3=DSQRT(2.D0/3.D0*((P(1)-P(2))**2+(P(2)-P(3))**2+
     &       (P(1)-P(3))**2+2.D0*(P(4)**2+P(5)**2+P(6)**2)))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SREDNA(SRED,SIGM,SEFE1,ISTNA,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     
      DIMENSION SRED(*),SIGM(*),STRESS(7)
C
      IF(ISTNA.EQ.1) CALL ZBIRM1(SRED,SIGM,6)
      IF(IND.EQ.0) THEN
         CALL JEDNA1(STRESS,SIGM,6)
      ELSE
         STRESS(1)=SIGM(1)
         STRESS(2)=SIGM(3)
         STRESS(3)=SIGM(6)
         STRESS(4)=SIGM(2)
         STRESS(5)=SIGM(5)
         STRESS(6)=SIGM(4)
      ENDIF
      CALL GLAVN3(STRESS)
      IF(ISTNA.EQ.2) THEN
         IF(STRESS(7).GT.SEFE1) THEN
            SEFE1=STRESS(7)
            CALL JEDNA1(SRED,SIGM,6)
         ENDIF
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE NULL39(PLAST,AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO NULL AND PREPEAR FOR RESTART
CS.   P R O G R A M
CS.       ZA NULOVANJE I PRIPREMU ZA RESTART
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION PLAST(*)
      DIMENSION AU(*)
      REAL AU
C
      IF(IDEBUG.GT.0) PRINT *, ' NULL39'
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1)GO TO 20
C
         IBTC=0
         NNPR=MODPRO(NMODM)
         DO 10 NGR=1,NGAUSX
         DO 10 NGS=1,NGAUSY
         DO 10 NGT=1,NGAUSZ
            IBTC=IBTC+1
            LL=1+(NLM-1)*NGS12*NNPR+(IBTC-1)*NNPR
C
            CALL CLEAR(PLAST(LL),18)
C  PREK.PR    P0  PLAST(LL+18)
            PLAST(LL+19)=0.D0
            PLAST(LL+20)=0.D0
   10    CONTINUE
   20 CONTINUE
      END
