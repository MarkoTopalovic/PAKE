C=======================================================================
C
CE       PRINT DISPLACEMENTS AND CALL ROUTINE FOR PRINT STRESSES
CS       STAMPANJE POMERANJA I POZIVANJE PROGRAMA ZA STAMPANJE NAPONA
C
C   SUBROUTINE STAMPA
C              DALIPR
C              STAPO1
C              STANAP
C              ERROR
C              STAGTE
C              STAGP1
C
C=======================================================================
      SUBROUTINE STAMPA(NPODS)
      USE FSIDENT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT RESULTS FOR DISPLACEMENTS AND STRESSES
CS.   PROGRAM
CS.      ZA STAMPANJE REZULTATA
C .
C ......................................................................
C
      include 'paka.inc'
      
      CHARACTER*6 IMD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GRUPER/ LIGRUP
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /MINMAX/ KMINP,KMAXP,KMING,KMAXG
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /DUPLAP/ IDVA
      COMMON /NELINE/ NGENN
      COMMON /SRPSKI/ ISRPS
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /UPRIRI/ LUPRI
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /AUTST2/ PARAM(5)
      COMMON /NIDEAS/ IDEAS
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
      COMMON /CDEBUG/ IDEBUG
      COMMON /ODUZPOM/ KODUZ
      COMMON /ISTRES/ ISTRESS
      COMMON /INTERA/ IINTER,NPTI
C
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' STAMPA'
      LMAX=LRAD
      IF(NBLPR.LE.0) INDPR=NBLPR
      IF(NBLGR.LE.0) INDGR=NBLGR
      IF(NBLPR.GT.0)
     1CALL DALIPR(A(LNDTPR),KOR,NDT,KMINP,KMAXP,NBLPR,INDPR)
      IF(NBLGR.GT.0)
     1CALL DALIPR(A(LNDTGR),KOR,NDT,KMING,KMAXG,NBLGR,INDGR)
C
C     TEZINSKE FUNKCIJE I ELEMENTI PO PRSTENOVIMA ZA EDI INTEGRAL
      IF(NCRACK.GT.0) THEN
         CALL STAQS(A(LQ),NP,MAXRIN)
      ENDIF
C
      KORD=LCORD
      IF(IATYP.GE.3) KORD=LCORUL
CE ZAPISIVANJE NOVIH KOORDINATA ZA FLUID-STRUKTURA INTERAKCIJA
      IF(NPTI.GT.0) THEN
      ALLOCATE (IDENT(2,NPTI))
      REWIND IINTER
      CALL READD(IDENT,NPTI,IINTER)
      CALL NOVEK(A(KORD),IDENT,NPTI,IINTER,NP)
      DEALLOCATE (IDENT)
      ENDIF
CR
CR    UBACENO ISPRED ZBOG ODUZIMANJA 
      CALL CLEAR(A(LUPRI),JEDN)
C
      IF(NZADP.NE.0.AND.(INDPR.GE.0.OR.INDGR.GE.0)) THEN
C           LFTDT=NPODS(JPBR,44)
C           LMAX13=NPODS(JPBR,45)-1
C           CALL READDD(A(LFTDT),JEDN,IPODS,LMAX13,LDUZI)
            LMAX13=NPODS(JPBR,36)-1
            LZADFM = LRAD
            LNZADF = LZADFM + NZADP*IDVA
            LNZADJ = LNZADF + NZADP
            LMAX   = LNZADJ + NZADP
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            NPRO=LMAX-LRAD
            CALL DELJIV(NPRO,2,INDL)
            IF(INDL.EQ.1) NPRO=NPRO+1
            CALL READDD(A(LZADFM),NPRO/IDVA,IPODS,LMAX13,LDUZI)
C            CALL CLEAR(A(LUPRI),JEDN)
            CALL REACTP(A(LFTDT),A(LNZADJ),NZADP,IZLAZ,ISRPS,
     +                  A(LRTDT),A(LUPRI),INDPR)
            IF(INDGR.GE.0.AND.IDEAS.GE.0)
     1      CALL STAGP1(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,IGRAF,4,
     +                  A(LNCVP),NCVPR)
            IF(INDGR.GE.0)
     1      CALL STAU09(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,49,51,
     +                  A(LNCVP),NCVPR)
            ISRBA=1
            IF(ISRBA.EQ.1) THEN
C              ZAPISIVANJE SILA NA ZADATIM POMERANJIMA ZA SRBU
               IDUM=37
C               IMD='ZFORCE'
C               IF(KOR.EQ.1)
C     1         OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',FORM='FORMATTED',
C     1               ACCESS='SEQUENTIAL')
               CALL SRBAP1(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,IDUM,4,
     +                     A(LNCVP),NCVPR,A(KORD))
C               CLOSE(IDUM,STATUS='KEEP')
            ENDIF
            CALL CLEAR(A(LUPRI),JEDN)
C            LMAX13=NPODS(JPBR,47)-1
C            CALL READDD(A(LUPRI),JEDN,IPODS,LMAX13,LDUZI)
      ENDIF
      IF(NDIN.NE.0) GO TO 10
C
CZ      
CZ    STAMPANJE TABELE ZA BOCAC
CZ      
      IF(NP.GT.400000) THEN
         CALL BOCAC1(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,0,
     +               A(LNCVP),NCVPR,A(kord),A(LELCV))
         CALL BOCAC2(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,0,
     +               A(LNCVP),NCVPR,A(kord),A(LELCV))
         CALL BOCAC3(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,0,
     +               A(LNCVP),NCVPR,A(kord),A(LELCV))
      ENDIF      
CR      
CR    ZAPISIVANJE POMERANJA ZA ODUZIMANJE (treba uraditi za dinamiku!)
CR      write(*,*)'kor, koduz',kor, koduz
      IF (KOR.GT.KODUZ.AND.KODUZ.GT.0) THEN
        CALL RSTAZ(A(LIPODS),LUPRI,92)
CR        write(3,*)'kor', kor
CR        call wrr6(a(lupri), jedn, 'upri')
      ENDIF
CR
      CALL ODUZBA(A(LUPRI),A(LRTDT),JEDN)
CE    PRINT DISPLACEMENTS TO THE OUTPUT FILE (*.LST)
c      IF(INDPR.EQ.0.OR.INDPR.EQ.1) THEN
         CALL STAPO1(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,0,
     +               A(LNCVP),NCVPR,A(kord),A(LELCV))
c
      if(istress.ne.0)
     +   CALL NOVIK(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,IND,
     +              A(kord),IZLAZ)
CE       TIME FUNCTIONS
         CALL STAVFN(VREME)
c      ELSE
c      IF(ISRPS.EQ.0)
c     1WRITE(IZLAZ,2000) KOR,NDT
c      IF(ISRPS.EQ.1)
c     1WRITE(IZLAZ,6000) KOR,NDT
c      ENDIF
CE    PRINT DISPLACEMENTS TO THE GRAPHIC FILE (*.UNV)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.IDEAS.GE.0)
     1CALL STAGP1(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,IGRAF,0,
     +            A(LNCVP),NCVPR)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.NCXFEM.EQ.0)
     1CALL STAU09(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,49,1,
     +            A(LNCVP),NCVPR)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.NCXFEM.GT.0)
     1CALL STAU09(A(LUPRI),A(LID),A(LCVEL),ICVEL,NP,49,1,
     +            A(LNCVP),NCVPR)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.METOD.GT.5)
     1 CALL STAUAL(A(LCVEL),ICVEL,NP,49,32)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.METOD.GT.5)
     1 CALL STAUAL(A(LCVEL),ICVEL,NP,49,33)
C     OVO NIJE REGULARNO PO UPUTSTVU ZA IDEAS (TREBA SAMO VN 748)     
C      IF(IDIREK.LE.0) THEN
      IF(IDIREK.LE.0.AND.NCVPR.EQ.0) THEN
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.ISTVN.NE.-1)
     1CALL STAGVN(A(LDRV),A(LCVEL),ICVEL,NP,IGRAF,1)
      IF((INDGR.EQ.0.OR.INDGR.EQ.1).AND.ISTVN.NE.-1)
     1CALL STAGVN(A(LDRV1),A(LCVEL),ICVEL,NP,IGRAF,2)
      ENDIF
CR
      CALL CLEAR(A(LUPRI),JEDN)
      GO TO 20
C
CE    CALCULATE OF INCREMENT OF DISPLACEMENTS, NEW VELOCITIES
CE    AND ACCELERATIONS IN (T+DT)
CS    RACUNANJE PRIRASTAJA POMERANJA, NOVIH BRZINA I UBRZANJA 
CS    U TRENUTKU (T+DT)
C
   10 CALL RPRUVA
      IF(INDPR.EQ.0.OR.INDPR.EQ.1) THEN
         CALL STAMP
      ELSE
      IF(NBLPR.GE.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) KOR,NDT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) KOR,NDT
      ENDIF
      ENDIF
      IF(INDGR.EQ.0.OR.INDGR.EQ.1)
     1CALL STAGP
      IF(ITEST.EQ.1) RETURN
   20 IF(JPS.GT.1.AND.JPBR.EQ.JPS1) RETURN
C     IF(INDPR.EQ.0.OR.INDPR.EQ.2.OR.INDGR.EQ.0.OR.INDGR.EQ.2)
CE    PRINT STRESSES
      CALL STANAP(A(LIGRUP))
C
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(' *** K O R A K    B R O J ',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
 6000 FORMAT(' *** S T E P    N U M B E R ',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE DALIPR(NDTPG,KOR,NDT,KMIN,KMAX,NBL,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      ZA ODLUCIVANJE DA LI SE VRSI STAMPANJE U KORAKU
CS.   PROGRAM
CS.      ZA ODLUCIVANJE DA LI SE VRSI STAMPANJE U KORAKU
C .
C ......................................................................
C
      DIMENSION NDTPG(NBL,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' DALIPR'
      KMIN=100000
      KMAX=0
      DO 20 I=1,NBL
      DO 20 J=1,3,2
         KM=NDTPG(I,J)
         IF(KM.LT.KMIN) KMIN=KM
         IF(KM.GT.KMAX) KMAX=KM
   20 CONTINUE
      IF(KOR.LT.KMIN) GO TO 100
      DO 30 IBL=1,NBL
      IF(KOR.GE.NDTPG(IBL,1).AND.KOR.LE.NDTPG(IBL,3).
     1OR.KOR.GE.NDTPG(IBL,3).AND.KOR.LE.NDTPG(IBL,1)) GO TO 40
   30 CONTINUE
      IF(KOR.GT.KMAX) GO TO 100
   40 IF(KOR.EQ.NDTPG(IBL,1).OR.KOR.EQ.NDTPG(IBL,3)) GO TO 200
      IAUT= (NDTPG(IBL,3)-NDTPG(IBL,1))/NDTPG(IBL,2)-1
      KN=NDTPG(IBL,1)
      DO 50 J=1,IAUT
      KN=KN+NDTPG(IBL,2)
      IF(KOR.EQ.KN) GO TO 200
   50 CONTINUE
C
  100 IND=-1
      RETURN
C
  200 IND=NDTPG(IBL,4)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAPO1(RTH,ID,NCVEL,ICVEL,NP,KOR,VREME,IND,NCVP,NCVPR,
     1                  CORD,NELCV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERBR/ ITER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /SRPSKI/ ISRPS
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /AUTST2/ PARAM(5)
      COMMON /AUTSTK/ DELUM,DELFM,KRAJP
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /MESLESS/ IGBM,ndif,idif(50),NKI,IKI(10)
      COMMON /MATIZO/ E,V
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION RTH(*),ID(NP,*),NCVEL(*),NCVP(*),FSP(6),cord(np,*),
     1          NELCV(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' STAPO1'
      IF(IND.EQ.0) THEN
         IF(JPS.GT.1) THEN
            IF(JPBR.EQ.JPS1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060)
            ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) JPBR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) JPBR
            ENDIF
         ENDIF
         IF(ICCGG.NE.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) AMAXK,AMINK,AMAXF,AMINF
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) AMAXK,AMINK,AMAXF,AMINF
         ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002)
      ENDIF
      IF(METOD.LE.5) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) KOR,VREME
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) KOR,VREME
      WRITE(78,6020) KOR,VREME
      ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2025) KOR,VREME
      IF(ISRPS.EQ.0)
     1WRITE(*,2025) KOR,VREME
      IF(ISRPS.EQ.0)
     1WRITE(*,2026) KPNOD,KPDIR,RTH(KEQ),PARAM(1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6025) KOR,VREME
      IF(ISRPS.EQ.1)
     1WRITE(*,6025) KOR,VREME
      IF(ISRPS.EQ.1)
     1WRITE(*,6026) KPNOD,KPDIR,RTH(KEQ),PARAM(1)
C     KRITERIJUMI ZA ZAUSTAVLJANJE PRORACUNA ZA AUTOMATSKO INKREMENTIRANJE
      IF(DABS(RTH(KEQ)).GT.DABS(DELUM).AND.DABS(DELUM).GT.1.D-9) KRAJP=1
      IF(DABS(VREME).GT.DABS(DELFM).AND.DABS(DELFM).GT.1.D-9) KRAJP=1
      ENDIF
c
      IF(NGENL.GT.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) ITER
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) ITER
      ENDIF
      IF(INDPR.NE.0.AND.INDPR.NE.1) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      PMAX=0.D0
      PK1=0.
      PK2=0.
      AK1=0.
      DO 20 I=1,NP
         DO 10 J=1,6
            FSP(J) = 0.D0
            K=ID(I,J)
            IF(K.EQ.0) GO TO 10
            IF(K.GT.0)THEN
               FSP(J)=RTH(K)
            ELSE
               FSP(J)=CONDOF(RTH,A(LCMPC),A(LMPC),K)
            ENDIF
   10    CONTINUE
         IF(ICVEL.EQ.0) THEN
            NI=I
         ELSE
            NI=NCVEL(I)
         ENDIF
         POM=FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3)
         IF(POM.GE.PMAX) THEN
            NNI=NI
            PMAX=POM
         ENDIF
         IF(NCVPR.GT.0) THEN
            IJ=NCVP(I)
            IF(IJ.EQ.0) GO TO 20
         ENDIF
C
C        KOEFICIJENT INTEZIVNOSTI NAPONA
C        --------------------------------
         IF(IGBM.GT.0.and.NKI.GT.1) THEN
            IF(NKI.LT.1)THEN
               WRITE(IZLAZ,*) 'PRSLINA NIJE DOBRO DEFINISANA'
	         STOP
            ENDIF
C
            IF(ICVEL.EQ.0) THEN
               NUL=IKI(1)
            ELSE
               NUL=NELCV(IKI(1))
            ENDIF
C
            DO 671 II=2,nKI
               IF(NI.eq.IKI(II))then
c            WRITE(IZLAZ,*)'NUL,NI,IKI(II),I',NUL,NI,IKI(II),I  
c            WRITE(IZLAZ,*)'CORD',CORD(I,1),CORD(NUL,1)  
c            WRITE(IZLAZ,*)'CORD',CORD(I,2),CORD(NUL,2)  
                  PI=3.14159265
                  AKA=(3.-V)/(1.+V)
                  AMI=E/(2.*(1.+V))
C
                  DX=CORD(I,1)-CORD(NUL,1)
                  DY=CORD(I,2)-CORD(NUL,2)
                  RAD=DSQRT(DX*DX+DY*DY)
                  WRITE(IZLAZ,*) 'RAD=',RAD
C
                  IF(RAD.LT.1.D-8)THEN
                     WRITE(IZLAZ,*)'RASTOJANJE IZMEDJU VRHA PRSLINE I',  
     1                             ' CVORA PREVISE MALO'
	               STOP
                  ENDIF
C
                  TETA=ATAN(DY/DX)
c	      TETA=0.
c	      RAD=1.E-4
                  SKOR=DSQRT(RAD/(2.*PI))/(4.*AMI)
                  AK=(2.*AKA-1.)*DCOS(TETA/2.)-DCOS(3.*TETA/2.)
                  BK=(2.*AKA+3.)*DSIN(TETA/2.)+DSIN(3.*TETA/2.)
                  CK=(2.*AKA+1.)*DSIN(TETA/2.)-DSIN(3.*TETA/2.)
                  DK=(2.*AKA-3.)*DCOS(TETA/2.)+DCOS(3.*TETA/2.)
C
                  SKOR1=DSQRT(5./(2.*PI))/(2.*AMI)
                  IF(II.EQ.2)AK1=AK1+4.*FSP(1)/SKOR1
                  IF(II.EQ.3)AK1=AK1-FSP(1)/SKOR1
                  WRITE(IZLAZ,*)' AKI=',AK1
                  PK11=(FSP(1)*DK+FSP(2)*BK)/(SKOR*(AK*DK+CK*BK))
                  PK22=(FSP(1)*CK-FSP(2)*AK)/(SKOR*(AK*DK+CK*BK))
                  WRITE(IZLAZ,*)' KI=',PK11
                  WRITE(IZLAZ,*)' KII=',PK22
                  PK1=PK1+PK11
                  PK2=PK2+PK22
               ENDIF
  671	      continue
         ENDIF
C        ----------------------------------------------------
         WRITE(IZLAZ,5100) NI,(FSP(J),J=1,6)
   20 CONTINUE
      IF(IGBM.GT.0.AND.NKI.GT.0)then
        WRITE(IZLAZ,*)'KOEFICIJENT INTEZIVNOSTI NAPONA  KI=',PK1/(NKI-1)
        WRITE(IZLAZ,*)'KOEFICIJENT INTEZIVNOSTI NAPONA KII=',PK2/(NKI-1)
        WRITE(IZLAZ,*)' AKI=',AK1
      endif
      PMAX=DSQRT(PMAX)
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NNI,PMAX
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NNI,PMAX
      WRITE(78,6040) NNI,PMAX
         IF(PMALL.LT.PMAX) THEN
            PMALL=PMAX
            NPMALL=NNI
            KPMALL=KOR
         ENDIF
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2140) NNI,PMAX
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6140) NNI,PMAX
         IF(BMALL.LT.PMAX) THEN
            BMALL=PMAX
            NBMALL=NNI
            KBMALL=KOR
         ENDIF
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2240) NNI,PMAX
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6240) NNI,PMAX
         IF(AMALL.LT.PMAX) THEN
            AMALL=PMAX
            NAMALL=NNI
            KAMALL=KOR
         ENDIF
      ENDIF
      RETURN
C
 5100 FORMAT(' ',I10,6(1PE12.4))
C-----------------------------------------------------------------------
 2060 FORMAT(//////6X,'K O N T U R N I    C V O R O V I')
 2070 FORMAT(//////6X,
     1'P O D S T R U K T U R A    B R O J ..... JPBR =',I5)
 2100 FORMAT(///
     111X,'NAJVECI  DIJAGONALNI CLAN MATRICE KRUTOSTI      =',1PE12.5/
     111X,'NAJMANJI DIJAGONALNI CLAN MATRICE KRUTOSTI      =',1PE12.5//
     111X,'NAJVECI  DIJAGONALNI CLAN FAKTORIZOVANE MATRICE =',1PE12.5/
     111X,'NAJMANJI DIJAGONALNI CLAN FAKTORIZOVANE MATRICE =',1PE12.5)
 2000 FORMAT(///'1'/' K O M P O N E N T E      P O M E R A N J A     I 
     1    R O T A C I J A'/' ',69('-'))
 2001 FORMAT(///'1'/' K O M P O N E N T E     B R Z I N E     I     U G
     1A O N E     B R Z I N E'/' ',73('-'))
 2002 FORMAT(///'1'/' K O M P O N E N T E   U B R Z A N J A   I   U G A 
     1O N O G   U B R Z A N J A'/' ',75('-'))
 2020 FORMAT(//' K O R A K    B R O J =',I5,20X,'V R E M E =',1PE12.5)
 2025 FORMAT(//' K O R A K    B R O J =',I5,20X,'L A M D A =',1PE12.5)
 2026 FORMAT(' CVOR =',I5,5X,'PRAVAC =',I5,5X,'POMERANJE =',1PE12.5,
     +       5X,'LUK =',2(1PE12.5)//)
 2030 FORMAT(//' B R O J    I T E R A C I J A =',I5)
 2010 FORMAT(//
     1' ',' CVOR',5X,'X',11X,'Y',11X,'Z',11X,'XT',10X,'YT',10X,'ZT')
 2040 FORMAT(//' MAKSIMALNO POMERANJE:'/' CVOR =',I6,'   MAX.POM. =',
     11PE12.4//)
 2140 FORMAT(//' MAKSIMALNA BRZINA:'/' CVOR =',I6,'   MAX.BRZ. =',
     11PE12.4//)
 2240 FORMAT(//' MAKSIMALNO UBRZANJE:'/' CVOR =',I6,'   MAX.UBR. =',
     11PE12.4//)
C-----------------------------------------------------------------------
 6060 FORMAT(//////6X,'F O R    C O N T O U R    N O D E S')
 6070 FORMAT(//////6X,
     1'S U B S T R U C T U R E   N U M B E R ..... JPBR =',I5)
 6100 FORMAT(///
     *11X,'LARGEST DIAGONAL ELEMENT OF STIFFNESS MATRIX ..... =',
     *1PE12.5/
     *11X,'SMALLEST DIAGONAL ELEMENT OF STIFFNESS MATRIX .... =',
     *1PE12.5//
     *11X,'LARGEST DIAGONAL ELEMENT OF THE FACTORIZED MATRIX  =',
     *1PE12.5/
     *11X,'SMALLEST DIAGONAL ELEMENT OF THE FACTORIZED MATRIX =',
     *1PE12.5)
 6000 FORMAT(///'1'/' COMPONENTS  OF  DISPLACEMENTS  AND  ROTATIONS'/
     1' ',45('-'))
 6001 FORMAT(///'1'/' COMPONENTS OF VELOCITIES AND ANGULAR VELOCITIES'/
     1' ',47('-'))
 6002 FORMAT(///'1'/' COMPONENTS OF ACCELERATIONS AND ANGULAR ACCELERATI
     1ONS'/' ',53('-'))
 6020 FORMAT(//' S T E P  N U M B E R =',I5,20X,'  T I M E =',1PE12.5)
 6025 FORMAT(//' S T E P  N U M B E R =',I5,20X,'L A M D A =',1PE12.5)
 6026 FORMAT(' NODE=',I5,4X,'DIRECTION=',I5,4X,'DISPLACEMENT=',1PE12.5,
     +       4X,'ARCH=',1PE12.5//)
 6030 FORMAT(//' NUMBER OF ITERATIONS =',I5)
 6010 FORMAT(//
     1' ',' NODE',5X,'X',11X,'Y',11X,'Z',11X,'XT',10X,'YT',10X,'ZT')
 6040 FORMAT(//' MAXIMUM DISPLACEMENT:'/' NODE =',I6,'   MAX DISPL  =',
     11PE12.4//)
 6140 FORMAT(//' MAXIMUM VELOCITY:'/' NODE =',I6,'   MAX VELOC  =',
     11PE12.4//)
 6240 FORMAT(//' MAXIMUM ACCELERATION:'/' NODE =',I6,'   MAX ACCEL  =',
     11PE12.4//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STANAP(IGRUP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   LOOP TO PRINT STRESS AT INTEGRATION POINT OVER ELEMENT GROUPS 
CS.   PETLJA PO GRUPAMA ELEMENATA RADI STAMPANJA NAPONA U GAUS TACKAMA
C .
C ......................................................................
C
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SUMELE/ ISUMEL,ISUMGR
      DIMENSION IGRUP(NGELEM,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STANAP'
      ISUMEL=0
CE    LOOP OVER ELEMENT GROUPS
      DO 100 NGE = 1,NGELEM
         NETIP = IGRUP(NGE,1)
         NE    = IGRUP(NGE,2)
         IATYP = IGRUP(NGE,3)
         NMODM = IGRUP(NGE,4)
         LMAX8 = IGRUP(NGE,5)
         LMAX=LRAD
C
CE       PRINT STRESSES
         CALL ELEME(NETIP,4)
C
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ERROR(I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT ERRORS
CS.   PROGRAM
CS.      ZA STAMPANJE GRESAKA I UPOZORENJA
C .
C ......................................................................
C
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ERROR '
C
      GO TO (1,2,3),I
C
    1 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) LMAX,MTOT
      STOP
C
    2 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) LMAX,MTOT
      STOP
C
    3 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) LMAX,MTOT
      STOP
C-----------------------------------------------------------------------
 2000 FORMAT(///' NEDOVOLJNA DIMENZIJA U OSNOVNOM RADNOM VEKTORU A ZA UC
     1ITAVANJE ULAZNIH PODATAKA'/' POTREBNA DIMENZIJA, LMAX=',I10/
     2' RASPOLOZIVA DIMENZIJA, MTOT=',I10)
C-----------------------------------------------------------------------
 6000 FORMAT(///'  ISUFFICIENT DIMENSION OF WORKING VECTOR FOR READING I
     1NPUT DATA'/
     1' REQUIRED STORAGE  .......  LMAX =',I10/
     2' MAXIMUM STORAGE AVAILABLE  MTOT =',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGTE(TEMP,NCVEL,ICVEL,NP,II,NCVP,NCVPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      TO PRINT TEMPERATURES IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE TEMPERATURA U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /NASLOV/ NASLOV
      COMMON /SRPSKI/ ISRPS
      DIMENSION TEMP(*),NCVEL(*),NCVP(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGTE'
      IF(NCVPR.GT.0) RETURN
      NPP=NP
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
CE    HEAT TRANSFER = 2
CS    PROVODJENJE TOPLOTE = 2
      IMOTY=2
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
CE    SCALAR = 1
CS    SKALAR = 1
      IDACH=1
CE    TEMPERATURES = 5
CS    TEMPERATURE = 5
      ISDTY=5
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=2
CE    NUMBER DATA = 1
CS    BROJ PODATAKA = 1
      NDVPN=1
      IND=-1
      ITYP=55
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 10 I=1,NPP
C         IF(DABS(TEMP(I)).LT.1.D-10) GO TO 10
         IF(ICVEL.EQ.0) THEN
            WRITE(II,5000) I
         ELSE
            WRITE(II,5000) NCVEL(I)
         ENDIF
         WRITE(II,5200) TEMP(I)
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(1PE13.5)
C-----------------------------------------------------------------------
 2000 FORMAT('CVORNE TEMPERATURE'/
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6000 FORMAT('NODAL TEMPERATURES'/
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGP1(RTH,ID,NCVEL,ICVEL,NODOVI,II,IND,NCVP,NCVPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS, VELOSITIES AND ACCELERATIONS
CE.      IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /NASLOV/ NASLOV
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /ITERBR/ ITER
      COMMON /SRPSKI/ ISRPS
      COMMON /NIDEAS/ IDEAS
      COMMON /KOJISR/ KOJSET
CSTOS
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /STOSZP/ STOSZ
CSTOS
      DIMENSION RTH(*),ID(NP,*),NCVEL(*),NCVP(*),FSP(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGP1'
      IF(ideas.eq.-1) return
      JEDAN=1
      NULA=0
      ZERO=0.
      M10=-10
      M11=-11
      IDV=2
      I8=8
C      
      NPP=NODOVI
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NODOVI-NPK
C      
C     MODEL TYPE
C     0-UNKNOWN; 1-STRUCTURAL; 2-HEAT TRANSFER; 3-FLUID FLOW
CE    STRUCTURAL ANALYSIS = 1
CS    STRUKTURNA ANALIZA = 1
      IMOTY=1
C
C     ANALYSIS TYPE
C     0-UNKNOWN; 1-STATIC; 2-NORMAL MODE; 4-TRANSIENT; 6-BUCKLING; 
C     9-STATIC NONLINEAR
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      IANTY=1
      IF(NDIN.GT.0.AND.ISOPS.GT.0) IANTY=2
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1.OR.NZADP.GT.0) THEN
         IANTY=4
         IF(NGENL.GT.0) IANTY=9
         IF(NDIN.EQ.0.AND.ISOPS.GT.0.AND.NGENL.GT.0) IANTY=6
         IFAT1=2
         FATY8=VREME
CSTOS
C	BOBAN: Iskljucio stampanje STOSZ umesto VREME
c      IF(NZADP.GT.0) FATY8=STOSZ
CSTOS
      ENDIF
C
C     DATA CHARACTERISTIC
C     0-UNKNOWN; 1-SCALAR; 2-3 DOF; 3-6 DOF; 4-SYMMETRIC TENSOR
C     6 DOF = 3
      IDACH=3
C
C     RESULT TYPE
C     2-STRESS; 3-STRAIN; 4-ELEMENT FORCE; 5-TEMPERATURE; 8-DISPLACEMENT
C     9-REACTION FORCE; 11-VELOCITY; 12-ACCELERATION; 32-CONSTRAINT FORCE
C     34-PLASTIC STRAIN; 35-CREEP STRAIN; 49-CONTACT FORCES;
C     94-UN SCALAR; 95-UN 3 DOF; 96-UN 6 DOF; 97-UN SYMMETRIC TENSOR;
C     >1000-USER DEFINED RESULT TYPE
CE    DISPLACEMENTS = 8
CS    POMERANJA = 8
      ISDTY=8
CE    FORCE = 4
CS    SILE = 4 
C      IF(IDEAS.GT.6) THEN
C         IF(IND.EQ.4) ISDTY=9
C      ELSE
C         IF(IND.EQ.4) ISDTY=4
C      ENDIF
      IF(IND.EQ.4) ISDTY=4
      IF(IND.EQ.32) ISDTY=32
      IF(IND.EQ.44) ISDTY=44
CE    VELOCITIES = 11
CS    BRZINE = 11
      IF(IND.EQ.1) ISDTY=11
CE    ACCELERATIONS = 12
CS    UBRZANJA = 12
      IF(IND.EQ.2) ISDTY=12
C
C     DATA TYPE
C     1-INTEGER; 2-SINGLE PRECISION; 4-DOUBLE PRECISION FLOATING POINT
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=2
C
C     NUMBER OF DATA VALUES FOR THE DATA COMPONENT (NVALDC)
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=6
C
C     DATASET LOCATION
C     1-DATA AT NODES; 2-DATA ON ELEMENTS; 3-DATA AT NODES ON ELEMENTS
      IELS=1
C
      IND1=-1
      IF(IDEAS.GT.6) THEN
         ITYP=2414
      ELSE
         ITYP=55
      ENDIF
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      IF(IDEAS.GT.6) THEN
         IF(II.EQ.18) THEN
            KOJSET=KOJSET+1
            WRITE(II,5000) KOJSET
         ELSE
            WRITE(II,5000) ITER
         ENDIF
         WRITE(II,5003) NASLOV
         WRITE(II,5000) IELS
      ENDIF
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(IND.EQ.4.OR.IND.EQ.44) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2004)
      IF(ISRPS.EQ.1)
     1WRITE(II,6004)
      ENDIF
      IF(IND.EQ.32) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2032)
      IF(ISRPS.EQ.1)
     1WRITE(II,6032)
      ENDIF
      IF(IND.EQ.44) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2044)
      IF(ISRPS.EQ.1)
     1WRITE(II,6044)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(IDEAS.GT.6) THEN
         IF(IANTY.GT.2) THEN
            WRITE(II,5000) M11,NULA,IDV,NULA,JEDAN,JEDAN,KOR,KOR
            WRITE(II,5000) I8,IDV
            WRITE(II,5200) FATY8,FATY8,ZERO,ZERO,ZERO,ZERO
            WRITE(II,5200) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO
         ELSE
            WRITE(II,5000) M10,NULA,JEDAN,JEDAN,JEDAN,NULA,NULA,NULA
            WRITE(II,5000) IDV,NULA
            WRITE(II,5200) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO
            WRITE(II,5200) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO
         ENDIF
      ELSE
         IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
         IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
         WRITE(II,5200) FATY8
      ENDIF
      DO 10 I=1,NPP
         IF(NCVPR.GT.0) THEN
            IJ=NCVP(I)
            IF(IJ.EQ.0) GO TO 10
         ENDIF
         IMA=0
         DO 20 J=1,6
            FSP(J) = 0.0D0
            IF(ID(I,J).EQ.0) GO TO 20
            K = ID(I,J)
            IF(K.GT.0) THEN
               FSP(J)=RTH(K)
            ELSE
               FSP(J)=CONDOF(RTH,A(LCMPC),A(LMPC),K)
            ENDIF
            IF(DABS(FSP(J)).GT.1.D-10) IMA=1
   20    CONTINUE
         IF(IMA.EQ.0.AND.(IND.EQ.4.OR.IND.EQ.32.OR.IND.EQ.44)) GO TO 10
         IF(ICVEL.EQ.0) THEN
            WRITE(II,5000) I
         ELSE
            WRITE(II,5000) NCVEL(I)
         ENDIF
         WRITE(II,5200) (FSP(J),J=1,6)
   10 CONTINUE
      WRITE(II,5100) IND1
      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(8I10)
 5200 FORMAT(6(1PE13.5))
C-----------------------------------------------------------------------
 2004 FORMAT('SILE U CVOROVIMA')
 2032 FORMAT('SILE NA ZADATIM POMERANJIMA')
 2044 FORMAT('FILTRACIONE SILE')
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6004 FORMAT('NODAL FORCES')
 6032 FORMAT('CONSTRAINT FORCES')
 6044 FORMAT('FILTRATION FORCES')
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOSITIES AND ANGLE VELOSITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGKS(RTH,ID,NCVEL,ICVEL,NP,II,NCVP,NCVPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT CONCENTRATED LOAD SET IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE SETA OPTERECENJA OD KONCENTRISANIH SILA
C .
C ......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /NIDEAS/ IDEAS
      COMMON /SRPSKI/ ISRPS
      DIMENSION RTH(*),ID(NP,*),NCVEL(*),NCVP(*)
      DIMENSION FSP(6),ISP(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGKS'
      IF(ideas.eq.-1) return
      IF(NCVPR.GT.0) RETURN
      JEDAN=1
      NULA=0
      ZERO=0.
      NPP=NP
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
CS    KONCENTRISANE SILE = 1
      IANTY=1
CS    BOJA CRVENA = 11
      ICOL=11
      IND1=-1
      IF(IDEAS.GT.6) THEN
         ITYP=790
      ELSE
         ITYP=782
      ENDIF   
      IMAS=0
      DO 10 I=1,NPP
         IF(ICVEL.EQ.0) THEN
            NC=I
         ELSE
            NC=NCVEL(I)
         ENDIF
         IMA=0
         DO 20 J=1,6
            FSP(J) = 0.0D0
            ISP(J) = 0
            IF(ID(I,J).EQ.0) GO TO 20
            K = ID(I,J)
            IF(DABS(RTH(K)).LT.1.0D-10) GO TO 20
            FSP(J)=RTH(K)
            ISP(J) = 1
            IMA=1
            IF(IMA.EQ.1.AND.IMAS.EQ.0) IMAS=1
            IF(IMAS.EQ.1) THEN
               IMAS=-1
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5000) KOR,IANTY
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
            ENDIF
   20    CONTINUE
         IF(IMA.EQ.1) THEN
            WRITE(II,5000) NC,ICOL,(ISP(J),J=1,6)
            IF(IDEAS.GT.6) THEN
               WRITE(II,5300) (FSP(J),J=1,6),
     1                        NULA,NULA,NULA,NULA,NULA,NULA
            ELSE
               WRITE(II,5200) (FSP(J),J=1,6)
            ENDIF   
         ENDIF
   10 CONTINUE
      IF(IMAS.EQ.-1) WRITE(II,5100) IND1
      RETURN
C
 5100 FORMAT(I6)
 5000 FORMAT(2I10,6I2)
 5200 FORMAT(6(1PE13.5))
 5300 FORMAT(3(1PE25.16)/3(1PE25.16)/6I10)
C-----------------------------------------------------------------------
 2000 FORMAT('KONCENTRISANE SILE SET :',I10)
C-----------------------------------------------------------------------
 6000 FORMAT('CONCENTRATED FORCES SET :',I10)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE MELGR(MAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     M A T E R I J A L I
C
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /REPERM/ MREPER(4)
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /NIDEAS/ IDEAS
      DIMENSION AMAT(33)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MELGR '
      CALL CLEAR(AMAT,33)
      GO TO (  1,  2,  3,  4,999,999,999,999,999,999,
     1       999,999,999,999,999,999,999,999,999,999,
     1       999,999,999,999,999,999,999,999,999,999),NMODM
C
    1 LFUN=MREPER(1)
      LNTA=MREPER(1)
      LTEM=MREPER(1)
      MATE=MREPER(4)
      IM1=2
      IM2=MATE
      ICL=27
      IPH=1
      NPR=27
      GO TO 11
C
    2 LFUN=MREPER(1)
      LNTA=MREPER(1)
      LTEM=MREPER(1)
      MATE=MREPER(4)
      IM1=9
      IM2=MATE
      ICL=33
      IPH=2
      NPR=33
      GO TO 11
C
    3 LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
      IM1=2
      IM2=MATE*3
      ICL=27
      IPH=1
      NPR=27
      GO TO 11
C
    4 LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
      IM1=2
      IM2=MATE*12
      ICL=33
      IPH=2
      NPR=33
C
   11 CALL MELG2(AMAT,A(LFUN),A(LNTA),A(LTEM),IM1,IM2,MAT,A(LGUSM))
      IF(IDEAS.EQ.8) THEN
         CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
      ELSEIF(IDEAS.EQ.7) THEN
         CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
      ELSE
         IF(ideas.eq.-1) return
         IND=-1
         ITYP=773
         WRITE(IGRAF,1100) IND
         WRITE(IGRAF,1100) ITYP
         WRITE(IGRAF,1000) ISUMGR,IPH,NPR
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2000) ISUMGR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6000) ISUMGR
         WRITE(IGRAF,1200) (AMAT(I),I=1,ICL)
         WRITE(IGRAF,1100) IND
	ENDIF
C
  999 RETURN 
C
 1000 FORMAT(8I10)
 1100 FORMAT(I6)
 1200 FORMAT(6(1PE13.5))
C-----------------------------------------------------------------------
 2000 FORMAT('MATERIJALNI SET :',I10)
C-----------------------------------------------------------------------
 6000 FORMAT('MATERIAL PROPERTIES SET :',I10)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE MELG2(AMAT,FUN,NTA,TEM,IM1,IM2,MAT,GUSM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     M A T E R I J A L I
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION AMAT(*),FUN(IM1,IM2,*),NTA(*),TEM(*),GUSM(50,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MELG2 '
      GO TO (  1,  2,  3,  4,999,999,999,999,999,999,
     1       999,999,999,999,999,999,999,999,999,999,
     1       999,999,999,999,999,999,999,999,999,999),NMODM
C
C    1 IF(NDIN.GT.0) AMAT(3)=GUSM(NMODM,MAT)
    1 AMAT(3)=GUSM(NMODM,MAT)
      AMAT(1)=FUN(1,MAT,1)
      AMAT(2)=FUN(2,MAT,1)
      AMAT(4)=AMAT(1)/2.D0/(1.D0+AMAT(2))
      RETURN
C
C    2 IF(NDIN.GT.0) AMAT(7)=GUSM(NMODM,MAT)
    2 AMAT(7)=GUSM(NMODM,MAT)
      AMAT(1)=FUN(1,MAT,1)
      AMAT(2)=FUN(2,MAT,1)
      AMAT(3)=FUN(3,MAT,1)
      AMAT(4)=FUN(4,MAT,1)
      AMAT(5)=FUN(5,MAT,1)
      AMAT(6)=FUN(6,MAT,1)
      AMAT(8)=FUN(7,MAT,1)
      AMAT(9)=FUN(8,MAT,1)
      AMAT(10)=FUN(9,MAT,1)
      RETURN
C
    3 TEMP0=TEM(MAT)
C      IF(NDIN.GT.0) AMAT(3)=GUSM(NMODM,MAT)
      AMAT(3)=GUSM(NMODM,MAT)
      AMAT(6)=TEM(MAT)
      DO 15 J=1,3
      NFE=(MAT-1)*3+J
      CALL TABF(FUN,NTA,NFE,IM2,TEMP0,EVA,NTMX,IND)
      IF(IND.GT.0) GO TO 120
      IF(J.EQ.1) AMAT(1)=EVA
      IF(J.EQ.2) AMAT(2)=EVA
   15 CONTINUE
      AMAT(5)=EVA
      AMAT(4)=AMAT(1)/2.D0/(1.D0+AMAT(2))
      RETURN
C
    4 TEMP0=TEM(MAT)
C      IF(NDIN.GT.0) AMAT(7)=GUSM(NMODM,MAT)
      AMAT(7)=GUSM(NMODM,MAT)
      AMAT(14)=TEM(MAT)
      DO 16 J=1,12
      NFE=(MAT-1)*12+J
      CALL TABF(FUN,NTA,NFE,IM2,TEMP0,EVA,NTMX,IND)
      IF(IND.GT.0) GO TO 120
      IF(J.EQ.1) AMAT(1)=EVA
      IF(J.EQ.2) AMAT(2)=EVA
      IF(J.EQ.3) AMAT(3)=EVA
      IF(J.EQ.4) AMAT(4)=EVA
      IF(J.EQ.5) AMAT(5)=EVA
      IF(J.EQ.6) AMAT(6)=EVA
      IF(J.EQ.7) AMAT(8)=EVA
      IF(J.EQ.8) AMAT(9)=EVA
      IF(J.EQ.9) AMAT(10)=EVA
      IF(J.EQ.10) AMAT(11)=EVA
      IF(J.EQ.11) AMAT(12)=EVA
      IF(J.EQ.12) AMAT(13)=EVA
   16 CONTINUE
C
  999 RETURN
C
  120 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NFE,TEMP0
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NFE,TEMP0
      STOP
C-----------------------------------------------------------------------
 2000 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U MELG2'/
     1' FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PE12.4)
C-----------------------------------------------------------------------
 6000 FORMAT(///' ARGUMENT IS OUT OF RANGE IN MELG2'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',D12.4)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE REACTP(FTDT,NZADJ,NZADP,IZLAZ,ISRPS,RTDT,UPRI,INDPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C .
CE.   PROGRAM
CE.      TO PRINT LOADS CORESPONDING TO PRESCRIBED DISPLACEMENTS
C
      DIMENSION FTDT(*),NZADJ(*),RTDT(*),UPRI(*)
CSTOS
      COMMON /STOSZP/ STOSZ
      STOSZ=0.D0
C PRIVREMENO ZA SIMETRICNE PROFILE KOD BRAZIEROVOG PROBLEMA
C      BRAZS=2.D0
      BRAZS=1.D0
CSTOS
      UREAK=0.D0
      DREAK=0.D0
CP      IF(ISRPS.EQ.0.AND.INDPR.GE.0)
CP     1WRITE(IZLAZ,2000)
CP      IF(ISRPS.EQ.1.AND.INDPR.GE.0)
CP     1WRITE(IZLAZ,6000)
         DO 10 I=1,NZADP
            NJ=NZADJ(I)
            IF(NJ.GT.0) THEN
               UPRI(NJ)=FTDT(NJ)
CP               IF(INDPR.GE.0) WRITE(IZLAZ,5100) NJ,RTDT(NJ),FTDT(NJ)
C               IF(NZADP.LT.20) WRITE(*,5100) KOR,VREME,RTDT(NJ),FTDT(NJ)
CSTOS
CK    IF(NJ.LT.19) THEN
         UREAK=UREAK+FTDT(NJ)*BRAZS
CK     ELSEIF(NJ.GT.26) THEN
CK       DREAK=DREAK+FTDT(NJ)
CK     ENDIF
CSTOS
            ENDIF
   10    CONTINUE
CSTOS
      STOSZ=UREAK
CKD    STOSZ=DREAK
CK    WRITE(*,*) UREAK,DREAK
C      IF(INDPR.GE.0) WRITE(IZLAZ,5001) UREAK,DREAK
C      WRITE(*,5001) UREAK,DREAK
C 5001 FORMAT( 'PRIVREMENO - UKUPNA REAKCIJA JE',2(1PE12.4,2X))
CSTOS
      RETURN
C
 5100 FORMAT(10X,I5,6(1PE12.4))
C-----------------------------------------------------------------------
 2000 FORMAT(///'1'/' S I L E   N A   Z A D A T I M    P O M E R A N J I
     1 M A'/' ',54('-')// 10X,' JEDN.',3X,'POMER.',5X,'SILA')
C-----------------------------------------------------------------------
 6000 FORMAT(///'1'/' LOADS CORESPONDING TO PRESCRIBED DISPLACEMENTS'/
     1' ',46('-')// 10X,' DOF ',3X,'DISPL.',5X,'LOAD')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGVN(RTH,NCVEL,ICVEL,NP,II,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT NODAL VECTORS VN I V1 IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE VEKTORA VN I V1 U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /NASLOV/ NASLOV
      COMMON /SRPSKI/ ISRPS
      COMMON /NIDEAS/ IDEAS
      DIMENSION RTH(NP,*),NCVEL(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGVN'
      IF(ideas.eq.-1) return
      NPP=NP
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
CE    STRUCTURAL ANALYSIS = 1
CS    STRUKTURNA ANALIZA = 1
      IMOTY=1
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
C     6 DOF = 3
      IDACH=3
CE    DISPLACEMENTS = 8
CS    POMERANJA = 8
C      IF(IND.EQ.0) ISDTY=8
CE    NORMAL VECTOR = 11
CS    VEKTOR NORMALE = 11
      IF(IND.EQ.1) ISDTY=11
CE    V1 VECTOR = 12
CS    VEKTOR V1 = 12
      IF(IND.EQ.2) ISDTY=12
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=2
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=6
      IND1=-1
      ITYP=55
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 10 I=1,NPP
         IF(ICVEL.EQ.0) THEN
            WRITE(II,5000) I
         ELSE
            WRITE(II,5000) NCVEL(I)
         ENDIF
         WRITE(II,5200) (RTH(I,J),J=1,3)
   10 CONTINUE
      WRITE(II,5100) IND1
      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
C-----------------------------------------------------------------------
 2006 FORMAT('NORMALE VN U CVOROVIMA ')
 2007 FORMAT('VEKTOR V1 U CVOROVIMA')
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6006 FORMAT('NODAL VECTOR VN')
 6007 FORMAT('NODAL VECTOR V1')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAVFN(VREME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  VREDNOST VREMENSKE FUNKCIJE
C
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /FVREME/ NTABFT,MAXTFT,LNTFT,LTABFT
      COMMON /SRPSKI/ ISRPS
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      DO 10 NF=1,NTABFT
         CALL TABF(A(LTABFT),A(LNTFT),NF,NTABFT,VREME,FT,NTMX,IND)
         WRITE(IZLAZ,1000) NF, FT
10    CONTINUE
1000  FORMAT(5X,I5,5X,1PE12.5)
2000  FORMAT(///' V R E M E N S K E     F U N K C I J E'/' ',38('-')
     1/' ',' FUNKCIJA',5X,'   VREDNOST')
6000  FORMAT(///' T I M E    F U N C T I O N S'/' ',28('-')/
     1' ',' FUNCTION',9X,'VALUE')
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STGRNA(SR,C,ISTNA,IGRFF,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  VREDNOST NAPONA U TEZISTU ELEMENTA ILI MAKSIMALNOG U ELEMENTU
C
      DIMENSION SR(*)
C
        IF(ISTNA.EQ.1) THEN
           CALL JEDNAK(SR,SR,C,KK)
           WRITE(IGRFF,2000) (SR(I),I=1,KK) 
        ENDIF
        IF(ISTNA.EQ.2) THEN
           WRITE(IGRFF,2000) (SR(I),I=1,KK) 
        ENDIF
      RETURN
 2000 FORMAT(6(1PE13.5))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGDC(ITEMP,NCVEL,ICVEL,NP,II,ITER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      TO PRINT DIVERGENT NODES IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE DIVERGENTINIH CVOROVA U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /NASLOV/ NASLOV
      COMMON /SRPSKI/ ISRPS
      DIMENSION ITEMP(*),NCVEL(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGDC'
      KOR=ITER+1
      VREME=FLOAT(ITER+1)
      NPP=NP
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
CE    HEAT TRANSFER = 2
CS    PROVODJENJE TOPLOTE = 2
      IMOTY=2
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      IFAT2=1
      IANTY=4
      IFAT1=2
      FATY8=VREME
CE    SCALAR = 1
CS    SKALAR = 1
      IDACH=1
CE    TEMPERATURES = 5
CS    TEMPERATURE = 5
      ISDTY=5
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=2
CE    NUMBER DATA = 1
CS    BROJ PODATAKA = 1
      NDVPN=1
      IND=-1
      ITYP=55
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 10 I=1,NPP
         IF(ITEMP(I).EQ.0) GO TO 10
         IF(ICVEL.EQ.0) THEN
            WRITE(II,5000) I
         ELSE
            WRITE(II,5000) NCVEL(I)
         ENDIF
         WRITE(II,5200) FLOAT(ITEMP(I))
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(1PE13.5)
C-----------------------------------------------------------------------
 2000 FORMAT('DIVERGENTNI CVOROVI'/
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6000 FORMAT('DIVERGENT NODES'/
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SRBAP1(RTH,ID,NCVEL,ICVEL,NP,II,IND,NCVP,NCVPR,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS, VELOSITIES AND ACCELERATIONS
CE.      IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      include 'paka.inc'
      
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /NASLOV/ NASLOV
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SRPSKI/ ISRPS
CSTOS
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /STOSZP/ STOSZ
CSTOS
      DIMENSION RTH(*),ID(NP,*),NCVEL(*),NCVP(*),FSP(6),CORD(NP,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' SRBAP1'
C
      IF(ISRPS.EQ.0.AND.INDPR.GE.0)
     1WRITE(IZLAZ,3000)
      IF(ISRPS.EQ.1.AND.INDPR.GE.0)
     1WRITE(IZLAZ,7000)
C
      NPP=NP
      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
CE    STRUCTURAL ANALYSIS = 1
CS    STRUKTURNA ANALIZA = 1
      IMOTY=1
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1.OR.NZADP.GT.0) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
CSTOS
      IF(NZADP.GT.0) FATY8=STOSZ
CSTOS
      ENDIF
C     6 DOF = 3
      IDACH=3
CE    DISPLACEMENTS = 8
CS    POMERANJA = 8
      ISDTY=8
CE    FORCE = 4
CS    SILE = 4 
      IF(IND.EQ.4) ISDTY=4
CE    VELOCITIES = 11
CS    BRZINE = 11
      IF(IND.EQ.1) ISDTY=11
CE    ACCELERATIONS = 12
CS    UBRZANJA = 12
      IF(IND.EQ.2) ISDTY=12
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=2
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=6
      IND1=-1
      ITYP=55
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(IND.EQ.4) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2004)
      IF(ISRPS.EQ.1)
     1WRITE(II,6004)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 10 I=1,NPP
         IF(NCVPR.GT.0) THEN
            IJ=NCVP(I)
            IF(IJ.EQ.0) GO TO 10
         ENDIF
         IMA=0
         DO 20 J=1,6
            FSP(J) = 0.0D0
            IF(ID(I,J).EQ.0) GO TO 20
            K = ID(I,J)
            IF(K.GT.0) THEN
               FSP(J)=RTH(K)
            ELSE
               FSP(J)=CONDOF(RTH,A(LCMPC),A(LMPC),K)
            ENDIF
            IF(DABS(FSP(J)).GT.1.D-10) IMA=1
   20    CONTINUE
         IF(IMA.EQ.0.AND.IND.EQ.4) GO TO 10
         IF(ICVEL.EQ.0) THEN
            WRITE(II,5201) I,(CORD(I,J),J=1,3),(FSP(J),J=1,6)
            WRITE(IZLAZ,7201) I,(FSP(J),J=1,6)
         ELSE
            WRITE(II,5201) NCVEL(I),(CORD(I,J),J=1,3),(FSP(J),J=1,6)
            WRITE(IZLAZ,7201) NCVEL(I),(FSP(J),J=1,6)
         ENDIF
   10 CONTINUE
      WRITE(II,5100) IND1
      RETURN
 7100 FORMAT(10X,I10,6(1PE12.4))
 7201 FORMAT(I10,9(1PE13.5))
C-----------------------------------------------------------------------
 3000 FORMAT(///'1'/' S I L E   N A   Z A D A T I M    P O M E R A N J I
     1 M A'/' ',54('-')// 1X,' CVOR',13X,'SILE FX,FY,FZ,MX,MY,MZ')
C-----------------------------------------------------------------------
 7000 FORMAT(///'1'/' LOADS CORESPONDING TO PRESCRIBED DISPLACEMENTS'/
     1' ',46('-')// 1X,' NODE',13X,'FORCE FX,FY,FZ,MX,MY,MZ')
C-----------------------------------------------------------------------
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
 5201 FORMAT(I5,9(1PE13.5))
C-----------------------------------------------------------------------
 2004 FORMAT('CVOR N, KOORDINATE X,Y,Z I SILE FX,FY,FZ,MX,MY,MZ')
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6004 FORMAT('NODE N, COORDINATES X,Y,Z AND FORCES FX,FY,FZ,MX,MY,MZ')
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOSITIES AND ANGLE VELOSITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      SUBROUTINE NOVEK(CORD,IDENT,NPTI,IINTER,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(NP,*),IDENT(2,NPTI)
      DO I=1,NPTI
C PROVERITI ZA SLOBODNU NUMERACIJU
         WRITE(IINTER) IDENT(2,I),(CORD(IDENT(1,I),J),J=1,3)
      ENDDO      
      RETURN
      END
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      SUBROUTINE NOVIK(RTH,ID,NCVEL,ICVEL,NP,KOR,VREME,IND,CORD,IZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      include 'paka.inc'
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION CORD(NP,*),NCVEL(*),ID(NP,*),ID6(6),FSP(6),RTH(*)

         WRITE(IZLAZ,1011) 
         DO 2 I=1,NP
         DO 10 J=1,3
            FSP(J) = 0.D0 
            K=ID(I,J)
            IF(K.EQ.0) GO TO 9
            IF(K.GT.0)THEN
               FSP(J)=RTH(K)
            ELSE
               FSP(J)=CONDOF(RTH,A(LCMPC),A(LMPC),K)
            ENDIF
    9       FSP(J)=FSP(J)+CORD(I,J)
   10    CONTINUE
            II=I
            IF(ICVEL.EQ.1) II=NCVEL(I)
            DO 3 J=1,6
               IF(ID(I,J).EQ.0) THEN
                  ID6(J)=1
               ELSE
                  ID6(J)=0
               ENDIF
    3       CONTINUE
         WRITE(IZLAZ,5000) II,(ID6(J),J=1,6),(FSP(J),J=1,3),0,1
    2    CONTINUE
 1011 FORMAT(///' NOVE KOORDINATE NA KRAJU PRORACUNA'/)
 5000 FORMAT(I5,1X,6I2,2X,3F10.4,2I5)
      RETURN
      END
