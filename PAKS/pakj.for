C=======================================================================
C
C=======================================================================
      SUBROUTINE UCRACK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ CRACKS
CS.    P R O G R A M
CS.        ZA UCITAVANJE PRSLINA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UCRACK'
C     MAKSIMALAN BROJ ELEMENATA PO JEDNOM PRSTENU
      N100=1000
      IF(NCRACK.EQ.0) RETURN
      CALL CLEAR(FK123,30)
      DO 10 J=1,NCRACK
         CALL ISPITA(ACOZ)
         IF(J.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (NODCR(J,I),I=1,14)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (NODCR(J,I),I=1,14)
         IF(NODCR(J,1).NE.2.AND.NODCR(J,1).NE.3) THEN
            STOP 'NODCR(J,1).NE.2.OR.NODCR(J,1).NE.3-UCRACK'
         ENDIF
C        PROCENAT INKREMENTA RASTA PRSLINE U ODNOSU NA POCETNU DUZINU PRVOG SEGMENTA PRSLINE	  
         IF(NODCR(J,2).EQ.0) NODCR(J,2)=100
C        POLUPRECNICI PRSTENOVA PO KOJIMA SE VRSI INTEGRACIJA U PROCENTIMA OD DUZINE PRSLINE	  
         IF(NODCR(J,6).EQ.0) NODCR(J,6)=20
         IF(NODCR(J,7).EQ.0) NODCR(J,7)=30
         IF(NODCR(J,8).EQ.0) NODCR(J,8)=40
         IF(NODCR(J,9).EQ.0) NODCR(J,9)=50
C        KADA JE NODCR(J,10)=2 RADI SE POLA MODELA I J-INTEGRAL TREBA DA SE MNOZI SA 2 
C        KADA JE NODCR(J,10)=1 RADI SE CEO  MODEL  I J-INTEGRAL TREBA DA SE MNOZI SA 1 
C        PROVERI DA LI ZA XFEM MOZE BILO KADA DA BUDE POLOVINA MODELA, AKO NE ONDA UVEK =1         
c         IF(NODCR(J,10).EQ.0.AND.NCXFEM.NE.0) NODCR(J,10)=1
         IF(NODCR(J,10).EQ.0) NODCR(J,10)=1
c         IF(NODCR(J,10).EQ.0) NODCR(J,10)=2
   10 CONTINUE   
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      DO 20 J=1,NCRACK
         WRITE(IZLAZ,5010) (NODCR(J,I),I=1,12)
   20 CONTINUE
      RETURN
C
 1000 FORMAT(14I5)
 5010 FORMAT(11X,14I6)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'P O D A C I   O   P R S L I N A M A'/
     16X,43('-')///11X,
     2' CRTYP  NPDA NODE1 NODE2 NODE3  NDP1  NDP2  NDP3  NDP4  ELTY  NSE
     1G IRING ')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'D A T A   F O R   C R A C K S'/
     16X,29('-')///11X,
     2' CRTYP  NPDA NODE1 NODE2 NODE3  NDP1  NDP2  NDP3  NDP4  ELTY  NSE
     1G IRING')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VEZACC(NELCV,NELC,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS NODES ELEMENTS IN FREE NUMERATION
CS.   P R O G R A M
CS.      ZA POVEZIVANJE BROJEVA CVOROVA ELEMENATA U SLOBODNOJ NUMERACIJI
C .
C ......................................................................
C
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NELCV(*),NELC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' VEZACC'
      DO 20 I=1,NCRACK
         IF(NODCR(I,1).NE.IND) GO TO 20
         DO 10 J=3,5
            N=NODCR(I,J)
            IF(N.LE.0) GO TO 10
            NN=N
            IF(N.GE.NPI.AND.N.LE.NPA) THEN
               N=NN-NPI+1
            ELSE
               N=NN-NGI+2+NPA-NPI
            ENDIF
            K=NELCV(N)
            IF(K.EQ.0) THEN
               N=NN-NGI+2+NPA-NPI
               K=NELCV(N)
            ENDIF
            IF(K.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NN
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NN
      STOP ' PROGRAM STOP - PAK21 - VEZACC'
            ENDIF
            NODCR(I,J)=K
   10    CONTINUE
C 
C         N=NODCR(I,2)
C         IF(N.LE.0) GO TO 20
C         NN=N
C         IF(N.GE.NMI.AND.N.LE.NMA) N=NN-NMI+1
C         K=NELC(N)
C         IF(K.EQ.0) THEN
C      IF(ISRPS.EQ.0)
C     1WRITE(IZLAZ,2001) NN
C      IF(ISRPS.EQ.1)
C     1WRITE(IZLAZ,6001) NN
C      STOP ' PROGRAM STOP - PAK21 - VEZACC'
C         ENDIF
C         NODCR(I,2)=K
C
C         N=NODCR(I,10)
C         IF(N.LE.0) GO TO 20
C         NN=N
C         IF(N.GE.NMI.AND.N.LE.NMA) N=NN-NMI+1
C         K=NELC(N)
C         IF(K.EQ.0) THEN
C      IF(ISRPS.EQ.0)
C     1WRITE(IZLAZ,2001) NN
C      IF(ISRPS.EQ.1)
C     1WRITE(IZLAZ,6001) NN
C      STOP ' PROGRAM STOP - PAK21 - VEZACC'
C         ENDIF
C         NODCR(I,10)=K
   20 CONTINUE
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' GRESKA U UCITAVANJU ULAZNIH PODATAKA O PRSLINAMA'/
     1' CVOR',I5,' NE POSTOJI'//' PROGRAM STOP - PAK21 - VEZACC'//)
 2001 FORMAT(//' GRESKA U UCITAVANJU ULAZNIH PODATAKA O PRSLINAMA'/
     1' ELEMENT',I5,' NE POSTOJI'//' PROGRAM STOP - PAK21 - VEZACC'//)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ERROR IN INPUT DATA ABOUT CRACKS'/
     1' NODE',I5,' NOT EXIST'//' PROGRAM STOP - PAK21 - VEZACC'//)
 6001 FORMAT(//' ERROR IN INPUT DATA ABOUT CRACKS'/
     1' ELEMENT',I5,' NOT EXIST'//' PROGRAM STOP - PAK21 - VEZACC'//)
C-----------------------------------------------------------------------
      END
C=======================================================================

C============================================================
      SUBROUTINE STAQS(Q,NP,MAXRIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
C    **************************************************
C    *   CRTANJE TEZINSKE FUNKCIJE ZA EDI INTEGRAL    *
C    **************************************************
      DIMENSION Q(NP,MAXRIN)
C
      DO I=1,MAXRIN
         CALL STAUML(Q(1,I),NP,49,35+I-1)
      ENDDO
      return
      end
C============================================================

C============================================================
      SUBROUTINE STAUML(RTH,NP,II,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      include 'paka.inc'
      
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION RTH(*)
C      DIMENSION NCVEL(*),NCVP(*)
C
      NPP=NP
C      IF(JPS.GT.1.AND.JPBR.LT.JPS1) NPP=NP-NPK
C      
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
C KOJI JE PROLAZ
      KOR=1
      WRITE(II,1000) KOR,IND,JEDAN
      IF(IND.EQ.35) WRITE(II,3001) 1 
      IF(IND.EQ.36) WRITE(II,3001) 2 
      IF(IND.EQ.37) WRITE(II,3001) 3 
      IF(IND.EQ.38) WRITE(II,3001) 4 
      IF(IND.EQ.39) WRITE(II,3001) 5 
      IF(IND.EQ.40) WRITE(II,3001) 6 
      WRITE(II,1010) ZERO,ZERO,ZERO
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,M6,M7
      WRITE(II,1000) JEDAN,NULA,JEDAN
      DO 10 I=1,NPP
C         IF(NCVPR.GT.0) THEN
C            IJ=NCVP(I)
C            IF(IJ.EQ.0) GO TO 10
C         ENDIF
C         IF(ICVEL.EQ.0) THEN
            WRITE(II,5000) I,RTH(I)
C         ELSE
C            WRITE(II,5000) NCVEL(I),RTH(I)
C         ENDIF
   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      RETURN
 1000 FORMAT(10(I5,','))
 1001 FORMAT(I6,F12.4)
 3035 FORMAT('PSI')
 3036 FORMAT('DISTANT')
 3037 FORMAT('TIME1')
 3038 FORMAT('VELOCITY')
 3001 FORMAT('INTEGRATION AREA',I3)
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME1',1PD12.4)
 1010 FORMAT(3(1PD12.4,','))
 5000 FORMAT(I5,',',3(1PD12.4,','))
      END
