C======================================================================
C
CE    DEFINE INPUT-OUTPUT UNITS
CS    DEFINISANJE ULAZNO - IZLAZNIH JEDINICA
C
C     ZATVOR
C     OTVORI
C     BRIS
C     BLOKOP
C     BFGSOP
C     VREM
C     EOFF
C======================================================================
C
C======================================================================
      SUBROUTINE PAKPAR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO SET PARAMETERS FOR PAK
CS.   P R O G R A M
CS.      ZA SETOVANJE PARAMETARA ZA PAK
C .
C ......................................................................
C
      LOGICAL OLDNEW
      CHARACTER*250    ACOZ
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /DUPLAP/ IDVA
      COMMON /PRIMER/ IPRBR,INDIZL,INDGRA,INDNEU
      COMMON /TEMPCV/ LTECV,ITEMP
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /IMPERF/ NMODS,LIDIM,LSCIM,MODES
      COMMON /BFGSUN/ IFILE
      COMMON /AXISYM/ INDAX
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /RADIZA/ INDBG
      COMMON /SRPSKI/ ISRPS
      COMMON /GREDE4/ IPRIR,INDGR4
      COMMON /REDJA / IREDJA
      COMMON /SCRATC/ ISCRC
      COMMON /PROBAS/ IILS
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
C
      NWK=0
      ISKDSK=0
C
CE    CODE TO SET ANALYSIS TYPE
CS    DEFINISANJE VRSTE ANALIZE
C
CE    LINEAR = 0 - LINEAR ANALYSIS
CE    LINEAR = 1 - NONLINEAR ANALYSIS
CS    LINEAR = 0 - LINEARNA
CS    LINEAR = 1 - NELINEARNA
C
      LINEAR=1
C
CE    INDDIN = 0 - STATIC
CE    INDDIN = 1 - DYNAMIC
CS    INDDIN = 0 - STATICKA
CS    INDDIN = 1 - DINAMICKA
C
      INDDIN=1
C
CE    INDIZL - INDICATOR FOR OPEN OUTPUT FILE (EQ.1, OPEN)
CE    INDGRA - INDICATOR FOR OPEN GRAPHICS FILE (EQ.1, OPEN)
CE    IPRBR  - INDICATOR FOR FIRST CASE
CE    KARTIC - NUMBER CARD IN INPUT DATA
CE    INDAX  - INDICATOR FOR AXISYMMETRY (INDAX=1, GROUP IS
CE             AXISYMMETRY)
CE    LUL    - INDICATOR FOR U.L. FORMULATION
CS    INDIZL - INDIKATOR OTVARANJA IZLAZNE DATOTEKE (=1, OTVORENA)
CE    ICONT  - INDICATOR FOR CONTACT ANALYSIS
CS    INDGRA - INDIKATOR OTVARANJA DATOTEKE ZA GRAFIKU (=1, OTVORENA)
CS    IPRBR  - INDIKATOR PRVOG PRIMERA
CS    KARTIC - BROJ TEKUCE KARTICE U ULAZNOJ DATOTECI
CS    INDAX  - INDIKATOR OSNE SIMETRIJE (ZAVISI OD TIPA ELEMENATA)
CS             AKO JE JEDNA GRUPA OSNOSIMETRICNA,  INDAX=1
CS    LUL    - INDIKATOR U.L. FORMULACIJE
CS    ICONT  - INDIKATOR KONTAKTNOG PROBLEMA
C
      INDIZL=0
      INDGRA=0
      INDNEU=0
      IPRBR =0
      KARTIC=0
      INDAX =0
      LUL   =0
      ICONT =0
      IILS  =0
      INDGR4=0
C
CE    CONFIG.PAK FILE
CS    CONFIG.PAK DATOTEKA
C
C
CE    T H E   F O L L O W I N G   U N I T S   A R E   U S E D
CS    D E F I N I S A N J E   K O R I S C E N I H   J E D I N I C A
C
CE    IULAZ - 'PAK.DAT' UNIT FOR INPUT DATA
CS    IULAZ - 'PAK.DAT' ULAZNA JEDINICA 'PAK.DAT'
      IULAZ=1
C
      IRTDT=2
CE    IZLAZ - 'PAK.LST' UNIT FOR OUTPUT DATA
CS    IZLAZ - 'PAK.LST' IZLAZNA JEDINICA
      IZLAZ=3
C
      ILISK=4
C DRAKCE KORISTI 7 U PAK62.FOR PRIVREMENO ZA STAMPANJE ?
      ILISE=7
CE    IELEM - 'ZIELEM' UNIT STORING ELEMENT DATA
CS    IELEM - 'ZIELEM' JEDINICA ZA PODATKE O ELEMENTIMA
      IELEM=8
C
C DRAKCE KORISTI 9 U PAK62.FOR PRIVREMENO ZA STAMPANJE ?
      IFTDT=9
      ILDLT=10
      ILIMC=11
      IDINA=12
      IPRIT=14
      ISILE=15
      IPOME=16
      IPRIR=17
CE    IGRAF - 'PAK.UNV' UNIT FOR GRAPHICS
CS    IGRAF - 'PAK.UNV' JEDINICA ZA GRAFIKU
      IGRAF=18
CE    ITEMP - UNIT FOR NODAL TEMPERATURES    
CS    ITEMP - JEDINICA ZA TEMPERATURE CVOROVA
      ITEMP=19
CE    IPODS - 'ZIPODS' UNIT STORING SUBSTRUCTURE DATA
CS    IPODS - 'ZIPODS' JEDINICA ZA PODATKE O PODSTRUKTURAMA
      IPODS=20
CE    MODES - 'ZMODES' UNIT FOR EIGEN VECTORS
CS    MODES - 'ZMODES' JEDINICA ZA SOPSTVENE VEKTORE - IMPERFEKCIJU
      MODES=21
C     IDUM - 'ZDIVER' POMERANJA U ITERACIJAMA 
      IDUM=27
C     IDUM - 'ZBELEM' BACKUP ZA 'ZIELEM'
      IDUM=28
C     ISUPB - FILE ZA SUPERGREDNI ELEMENT
      ISUPB=33
C     IFILT - FILE SA FILTRACIONIM SILAMA *.ZS 
      IFILT=35
C     IDUM - 'ZFORCE' ZAPISIVANJE SILA NA ZADATIM POMERANJIMA ZA SRBU
      IDUM=37
C     IDRAKCE - 'FDRAK' PODACI ZA ICCG=1
      IDRAKCE=39
CE    IBACK - 'ZIBACK' UNIT FOR RESTART
CS    IBACK - 'ZIBACK' JEDINICA ZA RESTART
      IBACK=40
C     IODUZ - 'POMER1' POMERANJA PRE RESTAROVANJA ZA MODELE ZEMLJE
      IODUZ=41
C     ISCRC - 'ZSKLIN', 'ZSKNEL', 'ZSCRCH' ZA BLOKOVE 
      ISCRC=45
CE    DETERMINE DATA FOR THE BLOCKS
CS    DEFINISANJE PODATAKA O BLOKOVIMA
CE    IBLK  - UNIT FOR BLOCKS
CS    IBLK  - JEDINICA ZA BLOKOVE
      IBLK =46
CE    IFILE - UNIT STORING BFGS VECTORS
CS    IFILE - JEDINICA ZA BFGS VEKTORE
      IFILE=47
C     IPAKS - 'ZIPAKS' PODACI ZA FLUIDE
      IPAKS=48
C     IZNEU - 'PAK.NEU' FILE ZA GRAFIKU 
      IZNEU=49
C DJORDJE KORISTI 50 U PAK61.FOR ?
C
C     IPAKF - 'ZIPAKF' PODACI ZA FLUIDE
      IPAKF=51
C     IINTER - 'ZIPAKI' ZA FLUIDE
      IINTER=62
C DRAKCE KORISTI 99 U PAK62.FOR PRIVREMENO ZA STAMPANJE ?
C
      INDFOR=2
CE    P A K   C O N F I G U R A T I O N
CS    K O N F I G U R A C I J A    P A K A
CE    IDEBUG - PRINT NAME SUBROUTINES (=0, NO; =1, YES)
CS    IDEBUG - STAMPANJE IMENA POTPROGRAMA (=0, NE; =1, DA)
      IDEBUG=0
CE    INDPC  - INDICATOR COMPUTER
CS    INDPC  - INDIKATOR VRSTE KOMPJUTERA
C              (=1, PC; =2, IBM; =3, VAX; =4, HP; =5, VAX)
      INDPC=1
CE    INDBG  - INDICATOR BACK GROUND ANALYSIS (=0, NO; =1, YES)
CS    INDBG  - INDIKATOR IZVODJ. PROGRAMA U POZADINI (=0, NE; =1, DA)
      INDBG=0
CE    ISRPS  - INDICATOR PROGRAM VERSION
CE             (=0, SERBIAN; =1, ENGLICH)
CS    ISRPS  - INDIKATOR VERZIJE PROGRAMA
CS             (=0, SRPSKA; =1, ENGLESKA)
      ISRPS=1
CE    CODE TO SET ARITHMETIC PRECISION
CE    IDVA = 1 - SINGLE PRECISION
CE    IDVA = 2 - DOUBLE PRECISION
CS    DEFINISANJE ARITMETICKE PRECIZNOSTI
CS    IDVA = 1 - JEDNOSTRUKA PRECIZNOST
CS    IDVA = 2 - DUPLA PRECIZNOST
      IDVA=2
CE    NBMAX - MAXIMUM NUMBER OF BLOCKS
CS    NBMAX - MAKSIMALAN BROJ BLOKOVA
      NBMAX=100
CE    LDUZI - RECORD LENGTH OF A DIRECT ACCESS FILE
CS    LDUZI - DUZINA STAZE KOD JEDINICE SA DIREKTNIM PRISTUPOM
      LDUZI=5120
CE    IREDJA=0, NUMBERS OF EQUATIONS ARE IN READING ORDER
CE    IREDJA=1, NUMBERS OF EQUATIONS ARE IN INCREASING ORDER
CS    IREDJA=0, BROJEVI JEDNACINA SU U UCITANOM REDOSLEDU
CS    IREDJA=1, BROJEVI JEDNACINA SU U RASTUCEM REDOSLEDU
      IREDJA=0
C
C      FIPNEU='PAK.NEU             '
      FIULAZ='PAK.DAT             '
      FIZLAZ='PAK.LST             '
      FIGRAF='PAK.UNV             '
      FIELEM='ZIELEM              '
      FIPODS='ZIPODS              '
      FIFILE='ZIFILE              '
      FMODES='ZMODES              '
      FITEMP='ZITEMP              '
      FIBLOK='ZIBLOK              '
      FIBACK='ZIBACK              '
      IF(INDPC.EQ.2) THEN
         OPEN (IULAZ,STATUS='OLD',FORM='FORMATTED',
     1         ACCESS='SEQUENTIAL')
      ELSE
         INQUIRE(FILE='CONFIG.PAK',EXIST=OLDNEW)
         IF(.NOT.OLDNEW) THEN
            RETURN
         ELSE
            OPEN (IULAZ,FILE='CONFIG.PAK',STATUS='OLD',FORM='FORMATTED',
     1            ACCESS='SEQUENTIAL')
         ENDIF
      ENDIF
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IDEBUG
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) INDPC 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) INDBG 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) ISRPS 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IDVA  
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IULA ,FIULAZ
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IZLAZ,FIZLAZ
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IGRAF,FIGRAF
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IELEM,FIELEM
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IPODS,FIPODS
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IFILE,FIFILE
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) MODES,FMODES
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) ITEMP,FITEMP
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IBLK ,FIBLOK 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IDINA,FIBACK
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) NBMAX 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) LDUZI 
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IREDJA
      CLOSE (IULAZ,STATUS='KEEP')
      IULAZ=IULA
      KARTIC=0
      IF(IDEBUG.GT.0) PRINT *, ' PAKPAR'
      RETURN
 1000 FORMAT (5X,I5,5X,A20)
      END
C=======================================================================
C
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZAGLAV
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO PRINT HEADER - PAKNEL
CS.    P R O G R A M
CS.        ZA STAMPANJA ZAGLAVLJA - PAKNEL
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ZAGLAV'
      IF(LINEAR.EQ.0.AND.INDDIN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      ENDIF
      IF(LINEAR.EQ.0.AND.INDDIN.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      ENDIF
      IF(LINEAR.EQ.1.AND.INDDIN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030)
      ENDIF
      IF(LINEAR.EQ.1.AND.INDDIN.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2103)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6103)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2115)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6115)
      RETURN
C----------------------------------------------------------------------
 2010 FORMAT('1'/////////33X,'P R O G R A M'/
     1 20X,'ZA STATICKU LINEARNU ANALIZU KONSTRUKCIJA'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 2020 FORMAT('1'/////////33X,'P R O G R A M'/
     1 14X,'ZA STATICKU I DINAMICKU LINEARNU ANALIZU KONSTRUKCIJA'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 2030 FORMAT('1'/////////33X,'P R O G R A M'/
     1 19X,'ZA STATICKU NELINEARNU ANALIZU KONSTRUKCIJA'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 2040 FORMAT('1'/////////33X,'P R O G R A M'/
     1 13X,'ZA STATICKU I DINAMICKU NELINEARNU ANALIZU KONSTRUKCIJA'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 2100 FORMAT(/////////////
     1 19X,' PPPPPPPPP         AAAA        KK      KK'/
     2 19X,' PPPPPPPPPP       AAAAAA       KK     KK'/
     3 19X,' PP      PP      AA    AA      KK    KK'/
     4 19X,' PP      PP     AA      AA     KK   KK'/
     5 19X,' PPPPPPPPPP     AA      AA     KK  KK'/
     6 19X,' PPPPPPPPP      AA      AA     KKKKK'/
     7 19X,' PP             AAAAAAAAAA     KKKKKKK'/
     8 19X,' PP             AAAAAAAAAA     KK    KK'/
     1 19X,' PP             AA      AA     KK     KK'/
     2 19X,' PP             AA      AA     KK      KK')
 2103 FORMAT(//////5X,'---------------------------  VERZIJA  8.00  ',
     1                '---------------------------')
 2115 FORMAT(///////////
     1 20X,'UNIVERZITET "SVETOZAR MARKOVIC"'/
     1 20X,'MASINSKI FAKULTET'/
     1 20X,'LABORATORIJA ZA INZENJERSKI SOFTVER'/
     1 20X,'34000 KRAGUJEVAC'/
     1 20X,'JUGOSLAVIJA'//////)
C----------------------------------------------------------------------
 6010 FORMAT('1'/////////28X,'FINITE ELEMENT PROGRAM'/
     1 26X,'FOR STATIC LINEAR ANALYSIS'/)
 6020 FORMAT('1'/////////28X,'FINITE ELEMENT PROGRAM'/
     1 20X,'FOR STATIC AND DYNAMIC LINEAR ANALYSIS'/)
 6030 FORMAT('1'/////////28X,'FINITE ELEMENT PROGRAM'/
     1 24X,'FOR STATIC NONLINEAR ANALYSIS'/)
 6040 FORMAT('1'/////////28X,'FINITE ELEMENT PROGRAM'/
     1 18X,'FOR STATIC AND DYNAMIC NONLINEAR ANALYSIS'/)
 6100 FORMAT(/////////////
     1 19X,' PPPPPPPPP         AAAA        KK      KK'/
     2 19X,' PPPPPPPPPP       AAAAAA       KK     KK'/
     3 19X,' PP      PP      AA    AA      KK    KK'/
     4 19X,' PP      PP     AA      AA     KK   KK'/
     5 19X,' PPPPPPPPPP     AA      AA     KK  KK'/
     6 19X,' PPPPPPPPP      AA      AA     KKKKK'/
     7 19X,' PP             AAAAAAAAAA     KKKKKKK'/
     8 19X,' PP             AAAAAAAAAA     KK    KK'/
     1 19X,' PP             AA      AA     KK     KK'/
     2 19X,' PP             AA      AA     KK      KK')
 6103 FORMAT(//////5X,'---------------------------  VERSION  8.00  ',
     1                '---------------------------')
 6115 FORMAT(///////////
     1 20X,'"SVETOZAR MARKOVIC" UNIVERSITY'/
     1 20X,'FAKULTY OF MECHANICAL ENGINEERING'/
     1 20X,'SOFTWARE LABORATORY'/
     1 20X,'34000 KRAGUJEVAC'/
     1 20X,'YUGOSLAVIA'///////)
C----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WIMIDZ
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ INPUT DATA AND PRINT CARD IMAGE
CS.    P R O G R A M
CS.        ZA UCITAVANJE ULAZNE DATOTEKE I STAMPANJE CARD IMAGE
C .
CS.        KARTIC - BROJ KARTICA VEC UCITANIH PRE POZIVA OVOG PROGRAMA
C .
C ......................................................................
C
      CHARACTER*250 INPDAT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WIMIDZ'
      IF (KARTIC.EQ.0) GO TO 10
      REWIND IULAZ
C
   10 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      KARD  = 0
      KARDNO = 0
C
   30 READ(IULAZ,1000,END=40) INPDAT
      KARD=KARD + 1
      KARDNO=KARDNO + 1
      IF (KARDNO.LE.50) GO TO 50
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      KARDNO=1
   50 WRITE(IZLAZ,5020) KARD,INPDAT
      IF(INPDAT(1:4).EQ.'STOP') GO TO 40
      GO TO 30
   40 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030)
      KARD=KARD - KARTIC + 1
C
      REWIND IULAZ
      DO 60 I=1,KARTIC
      READ(IULAZ,1000,END=60) INPDAT
   60 CONTINUE
      RETURN
C
 1000 FORMAT (A80)
 5020 FORMAT(I5,1X,A70)
C-----------------------------------------------------------------------
 2000 FORMAT(///'1'//' S L E D E   K A R T I C E   L I S T I N G A   ',
     1'U L A Z N E   D A T O T E K E'/28X,'K O L O N A    B R O J'/
     1' KART',10X,
     1'1         2         3         4         5         6         7'/
     1' BROJ 1234567890123456789012345678901234567890123456789012345',
     1 '678901234567890'/' ',4('-'),1X,70('-'))
 2010 FORMAT(' ',4('-'),1X,70('-')/
     1' KART',10X,
     1'1         2         3         4         5         6         7'/
     1' BROJ 1234567890123456789012345678901234567890123456789012345',
     1 '678901234567890'/28X,'K O L O N A    B R O J')
 2030 FORMAT(//' ',19('*'),' K R A J   U L A Z N O G   L I S T I N G A '
     1,13('*'))
C-----------------------------------------------------------------------
 6000 FORMAT(///'1'/' A   C A R D   I M A G E   L I S T I N G   O F   ',
     1'T H E   I N P U T   D A T A'/28X,'C O L U M N   N U M B E R'/
     1' CARD',10X,
     1'1         2         3         4         5         6         7'/
     1' NUMB 1234567890123456789012345678901234567890123456789012345',
     1 '678901234567890'/' ',4('-'),1X,70('-'))
 6010 FORMAT(' ',4('-'),1X,70('-')/
     1' CARD',10X,
     1'1         2         3         4         5         6         7'/
     1' NUMB 1234567890123456789012345678901234567890123456789012345',
     1 '678901234567890'/28X,'C O L U M N   N U M B E R')
 6030 FORMAT(//' ',19('*'),' E N D   O F   I N P U T   L I S T I N G ',
     115('*'))
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE ZATVOR(IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CLOSE FILE
CS.   P R O G R A M
CS.      ZA ZATVARANJE DATOTEKA
C .
C ......................................................................
C
      COMMON /PRIMER/ IPRBR,INDIZL,INDGRA,INDNEU
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /BFGSUN/ IFILE
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ZATVOR'
      IF(INDPC.EQ.2) RETURN
      IF(IND.EQ.1) THEN
        CLOSE (IULAZ,STATUS='KEEP')
        RETURN
      ENDIF
      IF(INDIZL.EQ.1) CLOSE (IZLAZ,STATUS='KEEP')
      IF(INDGRA.EQ.1) CLOSE (IGRAF,STATUS='KEEP')
      IF(METOD.EQ.5)  CLOSE (IFILE,STATUS='DELETE')
      IF(NBLOCK.GT.1) CLOSE (IBLK,STATUS='KEEP')
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE OTVORI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN INPUT FILE
CS.   P R O G R A M
CS.      ZA OTVARANJE ULAZNE DATOTEKE
C .
C ......................................................................
C
      CHARACTER*20    IME
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /RADIZA/ INDBG
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' OTVORI'
C
CE    INPUT FILE
CS    ULAZNA DATOTEKA
C
      IF(INDPC.EQ.2) RETURN
      IF(INDBG.EQ.1) THEN
      IME=FIULAZ
      GO TO 10
      ENDIF
CC      IF(INDPC.EQ.1) THEN
CC      IF(INDPC.GE.3.AND.INDPC.LE.5) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   IME ULAZNE DATOTEKE   <pak.dat>'
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   INPUT FILE NAME   <pak.dat>'
      READ(5,1000) IME
      IF(IME.EQ.'                    ') IME = 'pak.dat'
C
      CALL IMENA(IME)
C

   10 OPEN (IULAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
      RETURN
 1000 FORMAT (A20)
      END
C======================================================================
C
C======================================================================
      SUBROUTINE OTVSCR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILES
CS.   P R O G R A M
CS.      ZA OTVARANJE DATOTEKA
C .
C ......................................................................
C
      CHARACTER*6 IMD
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /RESTAR/ TSTART,IREST
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' OTVSCR'
C
CE    WORKING FILES
CS    RADNE DATOTEKE
C
         ISRE=78
         OPEN (ISRE,FILE='CONTROL.SRE',STATUS='UNKNOWN',
     1         FORM='FORMATTED',ACCESS='SEQUENTIAL')
C              ZAPISIVANJE POCETNIH NAPONA
               IDUM=77
               IMD='ZSTRIN'
               OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1               ACCESS='SEQUENTIAL')
      IF(INDPC.EQ.2) THEN
         OPEN (IELEM,STATUS='UNKNOWN',
     1         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
         OPEN (IPODS,STATUS='UNKNOWN',
     2         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
      ELSE
         IF(IREST.EQ.0) THEN
C CERNI
            IDUM=28
            IMD='ZBELEM'
            OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',
     1            FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
            CLOSE (IDUM,STATUS='DELETE')
C CERNI
C              ZAPISIVANJE SILA NA ZADATIM POMERANJIMA ZA SRBU
               IDUM=37
               IMD='ZFORCE'
               OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',FORM='FORMATTED',
     1               ACCESS='SEQUENTIAL')
C
            OPEN (IELEM,FILE=FIELEM,STATUS='UNKNOWN',
     1            FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
            CLOSE (IELEM,STATUS='DELETE')
            OPEN (IPODS,FILE=FIPODS,STATUS='UNKNOWN',
     2            FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
            CLOSE (IPODS,STATUS='DELETE')
         ENDIF
         OPEN (IELEM,FILE=FIELEM,STATUS='UNKNOWN',
     1         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
         OPEN (IPODS,FILE=FIPODS,STATUS='UNKNOWN',
     2         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE OTVIZL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN OUTPUT FILE
CS.   P R O G R A M
CS.      ZA OTVARANJE IZLAZNE DATOTEKE
C .
C ......................................................................
C
      CHARACTER*20    IME
      CHARACTER *24 PAKLST,PAKUNV,PAKNEU
      LOGICAL OLDNEW
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /RESTAR/ TSTART,IREST
      COMMON /IMEULZ/ PAKLST,PAKUNV,PAKNEU
      COMMON /IMEDUZ/ IDUZIN
      COMMON /RADIZA/ INDBG
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' OTVIZL'
C
CE    OUTPUT FILE
CS    IZLAZNA DATOTEKA
C
      IF(INDPC.EQ.2) RETURN
      IF(INDBG.EQ.1) THEN
         IME=FIZLAZ
         GO TO 10
      ENDIF
    5 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   IME IZLAZNE DATOTEKE   <',
     1               PAKLST(1:IDUZIN),'>'
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   OUTPUT FILE NAME   <',
     1               PAKLST(1:IDUZIN),'>'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKLST(1:20)
      ENDIF
      GO TO 10
C
    6 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   UNESI NOVO IME IZLAZNE DATOTEKE ZA RESTART' 
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   INPUT NEW OUTPUT FILE NAME FOR RESTART'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKLST(1:20)
      ENDIF
C
   10 INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(IREST.NE.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IZLAZ,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            OPEN (IZLAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ENDIF
      ENDIF
C     RESTART
      IF(IREST.EQ.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IZLAZ,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            GO TO 6
C           OPEN (IZLAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
C    1                           ACCESS='SEQUENTIAL')
C           CALL EOFF(IZLAZ,1)
         ENDIF
      ENDIF
      IF(INDBG.EQ.1) RETURN
C
      IND=0
      IF(OLDNEW.AND.IREST.NE.1) CALL BRIS (IME,IZLAZ,IND)
      IF(IND.EQ.1) GO TO 10
      IF(IND.EQ.2) GO TO 5
      RETURN
 1000 FORMAT (A20)
      END
C======================================================================
C
C======================================================================
      SUBROUTINE OTVGRA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILE FOR GRAFIC
CS.   P R O G R A M
CS.      ZA OTVARANJE DATOTEKE ZA GRAFIKU
C .
C ......................................................................
C
      CHARACTER*20    IME
      CHARACTER *24 PAKLST,PAKUNV,PAKNEU
      LOGICAL OLDNEW
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /RESTAR/ TSTART,IREST
      COMMON /IMEULZ/ PAKLST,PAKUNV,PAKNEU
      COMMON /IMEDUZ/ IDUZIN
      COMMON /RADIZA/ INDBG
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' OTVGRA'
C
CE    FILE FOR GRAFIC
CS    IZLAZ ZA GRAFIKU
C
      IF(INDPC.EQ.2) RETURN
      IF(INDBG.EQ.1) THEN
         IME=FIGRAF
         GO TO 20
      ENDIF
   15 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   IME GRAFICKE DATOTEKE   <',
     1               PAKUNV(1:IDUZIN),'>'
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   OUTPUT GRAPHICS FILE NAME   <',
     1               PAKUNV(1:IDUZIN),'>'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKUNV(1:20)
      ENDIF
      GO TO 20
C
    6 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   UNESI NOVO IME GRAFICKE DATOTEKE ZA RESTART' 
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   INPUT NEW OUTPUT GRAPHICS FILE NAME FOR RESTART'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKUNV(1:20)
      ENDIF
C
   20 INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(IREST.NE.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IGRAF,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            OPEN (IGRAF,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ENDIF
      ENDIF
C     RESTART
      IF(IREST.EQ.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IGRAF,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            GO TO 6
C           OPEN (IGRAF,FILE=IME,STATUS='OLD',FORM='FORMATTED',
C    1                        ACCESS='SEQUENTIAL')
C           CALL EOFF(IGRAF,1)
         ENDIF
      ENDIF
      IF(INDBG.EQ.1) RETURN
C
      IND=0
      IF(OLDNEW.AND.IREST.NE.1) CALL BRIS (IME,IGRAF,IND)
      IF(IND.EQ.1) GO TO 20
      IF(IND.EQ.2) GO TO 15
      RETURN
 1000 FORMAT (A20)
      END
C======================================================================
C
C======================================================================
      SUBROUTINE OTVNEU
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN NEU FILE FOR GRAFIC
CS.   P R O G R A M
CS.      ZA OTVARANJE NEU DATOTEKE ZA GRAFIKU
C .
C ......................................................................
C
      CHARACTER*20    IME
      CHARACTER*24    PAKLST,PAKUNV,PAKNEU
      LOGICAL         OLDNEW
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK,FIPNEU
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /RESTAR/ TSTART,IREST
      COMMON /IMEULZ/ PAKLST,PAKUNV,PAKNEU
      COMMON /IMEDUZ/ IDUZIN
      COMMON /RADIZA/ INDBG
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' OTVNEU'
C
      IZNEU=49
      FIPNEU='PAK.NEU             '
C
CE    FILE FOR GRAFIC
CS    IZLAZ ZA GRAFIKU
C
      IF(INDPC.EQ.2) RETURN
      IF(INDBG.EQ.1) THEN
         IME='PAK.NEU'
         IME=FIPNEU
         GO TO 20
      ENDIF
   15 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   IME NEU GRAFICKE DATOTEKE   <',
     1               PAKNEU(1:IDUZIN),'>'
      IF(ISRPS.EQ.1)
     1WRITE(6,*) '   OUTPUT NEU GRAPHICS FILE NAME   <',
     1               PAKNEU(1:IDUZIN),'>'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKNEU(1:20)
      ENDIF
      GO TO 20
C
    6 IF(INDPC.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(6,*) '   UNESI NOVO IME NEU GRAFICKE DATOTEKE ZA RESTART' 
      IF(ISRPS.EQ.1)
     1WRITE(6,*)'   INPUT NEW OUTPUT NEU GRAPHICS FILE NAME FOR RESTART'
         READ(5,1000) IME
         IF(IME.EQ.'                    ') IME = PAKNEU(1:20)
      ENDIF
C
   20 INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(IREST.NE.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IZNEU,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            OPEN (IZNEU,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ENDIF
      ENDIF
C     RESTART
      IF(IREST.EQ.1) THEN
         IF(.NOT.OLDNEW) THEN
            OPEN (IZNEU,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1                           ACCESS='SEQUENTIAL')
         ELSE
            GO TO 6
C           OPEN (IZNEU,FILE=IME,STATUS='OLD',FORM='FORMATTED',
C    1                        ACCESS='SEQUENTIAL')
C           CALL EOFF(IZNEU,1)
         ENDIF
      ENDIF
      IF(INDBG.EQ.1) RETURN
C
      IND=0
      IF(OLDNEW.AND.IREST.NE.1) CALL BRIS (IME,IZNEU,IND)
      IF(IND.EQ.1) GO TO 20
      IF(IND.EQ.2) GO TO 15
      RETURN
 1000 FORMAT (A20)
      END
C======================================================================
C
C=======================================================================
      SUBROUTINE BRIS (IME,IUN,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ERASE FILE
CS.   P R O G R A M
CS.      ZA BRISANJE FILE-A
C .
C ......................................................................
C
      CHARACTER*1     CH
      CHARACTER*20    IME 
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BRIS  '
      IND=2
      CH = 'N'
      IF(ISRPS.EQ.0)
     1WRITE(6,2000) IME
      IF(ISRPS.EQ.1)
     1WRITE(6,6000) IME
      READ(5,1000) CH
      IF(CH.EQ.' ') CH = 'D'
      IF(CH.EQ.'D'.OR.CH.EQ.'d') THEN
         CLOSE (IUN,STATUS='DELETE')
         IND=1
      ELSE
         CLOSE (IUN,STATUS='KEEP')
      ENDIF
      RETURN
 1000 FORMAT(A1)
C-----------------------------------------------------------------------
 2000 FORMAT ('   DATOTEKA    ',A20,' VEC POSTOJI'//
     1'  DA LI JE PREBRISATI  ( D / N ) "D" ? ')
C-----------------------------------------------------------------------
 6000 FORMAT ('   FILE    ',A20,' ALREADY EXISTS'//
     1'                 PRESS "ENTER" TO DELETE OR'/
     1'                 KEY   "N"     TO BYPASS')
C-----------------------------------------------------------------------
      END
C======================================================================
C
C======================================================================
      SUBROUTINE BLOKOP(IBLK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILE FOR BLOCKS
CS.   P R O G R A M
CS.      ZA FILE ZA BLOKOVE - OTVARA SE SAMO PO POTREBI
C .
C ......................................................................
C
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BLOKOP'
      IF(INDPC.EQ.2) THEN
         OPEN (IBLK,ACCESS='DIRECT',
     1         STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
      ELSE
         OPEN (IBLK,FILE=FIBLOK,ACCESS='DIRECT',
     1         STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
      ENDIF
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE MODEOP(MODES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILE FOR EIGEN VECTORS
CS.   P R O G R A M
CS.      FILE ZA SOPSTVENE VEKTORE - INICIJALANA IMPERFEKCIJA
C .
C ......................................................................
C
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /RESTAR/ TSTART,IREST
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MODEOP'
      IF(INDPC.EQ.2) RETURN
      OPEN (MODES,FILE=FMODES,STATUS='UNKNOWN',
     1      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE BFGSOP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILE FOR BFGS VECTORS
CS.   P R O G R A M
CS.      ZA FILE ZA BFGS VEKTORE - OTVARA SE SAMO PO POTREBI
C .
C ......................................................................
C
      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
     1                FITEMP,FIBLOK,FIBACK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BFGSUN/ IFILE
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BFGSOP'
      IF(INDPC.EQ.2) THEN
         OPEN(IFILE,ACCESS='DIRECT',
     1        STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
      ELSE
         OPEN(IFILE,FILE=FIFILE,ACCESS='DIRECT',
     1        STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
      ENDIF
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE VREM(KV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRINT TIME OF PROGRAM EXECUTION 
CS.   P R O G R A M
CS.      ZA STAMPANJE VREMENA IZVRSAVANJA PROGRAMA 
C .
C ......................................................................
C
      CHARACTER*26    BUF
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /VREMEN/ IVREM(4,4),IDATE(3,4)
      COMMON /VREMEH/ BUF(4)
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /SRPSKI/ ISRPS
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /CDEBUG/ IDEBUG    
      integer Dtime(8)
C
      IF(IDEBUG.GT.0) PRINT *, ' VREM  '
C
      IF(INDPC.EQ.1) THEN
         CALL DATE_AND_TIME(VALUES=Dtime)
         do I=1,3
            IDATE(4-I,KV)=Dtime(I)
         enddo
         Dtime(8)=Dtime(8)/10
         do I=5,8
            IVREM(I-4,KV)=Dtime(I)
         enddo
C        CALL GETDAT(IDATE(3,KV),IDATE(2,KV),IDATE(1,KV))
C        CALL GETTIM(IVREM(1,KV),IVREM(2,KV),IVREM(3,KV),IVREM(4,KV))
        IF(KV.LT.4) RETURN
C      IF(NULAZ.GE.0.OR.NBLPR.GE.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001) (IVREM(I,1),I=1,4),(IDATE(I,1),I=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001) (IVREM(I,1),I=1,4),(IDATE(I,1),I=1,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) (IVREM(I,2),I=1,4),(IDATE(I,2),I=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) (IVREM(I,2),I=1,4),(IDATE(I,2),I=1,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003) (IVREM(I,3),I=1,4),(IDATE(I,3),I=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003) (IVREM(I,3),I=1,4),(IDATE(I,3),I=1,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2004) (IVREM(I,4),I=1,4),(IDATE(I,4),I=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6004) (IVREM(I,4),I=1,4),(IDATE(I,4),I=1,3)
      WRITE(78,6004) (IVREM(I,4),I=1,4),(IDATE(I,4),I=1,3)
C      ENDIF
      ENDIF
      IF(INDPC.EQ.4) THEN
C        CALL HVREME(BUF(KV))
         IF(KV.LT.4) RETURN
C      IF(NULAZ.GE.0.OR.NBLPR.GE.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) (BUF(I)(11:20),I=1,KV) 
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) (BUF(I)(11:20),I=1,KV)
C      ENDIF
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(///'1'///6X,'I Z V E S T A J    O    V R E M E N U    I Z V
     1 R S A V A N J A'/6X,62('-')///)
 2001 FORMAT(11X,
     2'POCETAK UCITAVANJA ULAZNIH PODATAKA ............ ',
     2 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 2002 FORMAT(11X,
     4'POCETAK IZVRSAVANJA PO PERIODIMA ............... ',
     3 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 2003 FORMAT(11X,
     4'POCETAK IZRACUNAVANJA SOPSTVENIH VREDNOSTI ..... ',
     4 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 2004 FORMAT(11X,
     5'ZAVRSETAK PRORACUNA ............................ ',
     5 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4///)
 2100 FORMAT(///'1'///6X,'I Z V E S T A J    O    V R E M E N U    I Z V
     1 R S A V A N J A'/6X,62('-')///11X,
     2'POCETAK UCITAVANJA ULAZNIH PODATAKA ............ ',
     2 A10//11X,
     4'POCETAK IZVRSAVANJA PO PERIODIMA ............... ',
     3 A10//11X,
     4'POCETAK IZRACUNAVANJA SOPSTVENIH VREDNOSTI ..... ',
     4 A10//11X,
     5'ZAVRSETAK PRORACUNA ............................ ',
     5 A10///)
C-----------------------------------------------------------------------
 6000 FORMAT(///'1'///6X,'E X E C U T I O N    T I M I N G    I N F O R 
     1M A T I O N'/6X,57('-')///)
 6001 FORMAT(5X,
     2'TIME AT THE START OF READING OF INPUT DATA ..... ',
     2 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 6002 FORMAT(5X,
     3'TIME AT THE START OF FORMING THE MATRICES ...... ',
     3 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 6003 FORMAT(5X,
     4'TIME AT THE START OF EIGEN VALUE SOLUTION ...... ',
     4 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4//)
 6004 FORMAT(5X,
     5'END OF JOB ..................................... ',
     5 I2,':',I2,':',I2,'.',I2,3X,I2,'/',I2,'/',I4///)
 6100 FORMAT(///'1'///6X,'E X E C U T I O N    T I M I N G    I N F O R 
     1M A T I O N'/6X,57('-')///11X,
     2'TIME AT THE START OF READING OF INPUT DATA ..... ',
     2 A10//11X,
     3'TIME AT THE START OF FORMING THE MATRICES ...... ',
     3 A10//11X,
     4'TIME AT THE START OF EIGEN VALUE SOLUTION ...... ',
     4 A10//11X,
     5'END OF JOB ..................................... ',
     5 A10///)
C-----------------------------------------------------------------------
      END
C======================================================================
C
C=======================================================================
      SUBROUTINE EOFF(IUNIT,IPOZIV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO REWIND TO END SEQUENTIAL FILE
CS.   P R O G R A M
CS.      ZA DOLAZAK NA KRAJ SEKVENCIJALNE DATOTEKE
C .
C ......................................................................
C
      GO TO (10,20)IPOZIV
CE    FORMATTED
CS    FORMATIRANA
   10 READ(IUNIT,*,END=999)
      GO TO 10
CE    UNFORMATTED
CS    NEFORMATIRANA
   20 READ(IUNIT,END=999)
      GO TO 20
C
  999 BACKSPACE IUNIT
      RETURN
      END
C=======================================================================
      SUBROUTINE IMENA(IME)
      COMMON /IMEULZ/ PAKLST,PAKUNV,PAKNEU
      COMMON /IMEDUZ/ IDUZIN
      COMMON /CERSIL/ ZSILE
      CHARACTER *20 IME
      CHARACTER *3 VLST,MLST
      CHARACTER *3 VUNV,MUNV
      CHARACTER *3 VNEU,MNEU
      CHARACTER *3 LST,UNV,NEU
      CHARACTER *24 PAKLST,PAKUNV,PAKNEU
      CHARACTER *24 ZSILE
      VLST='LST'
      VUNV='UNV'
      VNEU='NEU'
      MLST='lst'
      MUNV='unv'
      MNEU='neu'
      IF(IME(1:1).GE.'A' .AND. IME(1:1).LE.'Z') THEN
         UNV=VUNV
         LST=VLST
         NEU=VNEU
      ELSE
         UNV=MUNV
         LST=MLST
         NEU=MNEU
      ENDIF
      IB=INDEX(IME,'.')
      DO 20 I=1,20
         IF(IME(I:I).EQ.' ') GOTO 30
         IA=I
  20  CONTINUE      
  30  IF(IB.EQ.0) THEN
         PAKLST=IME(1:IA)//'.'//LST
         PAKUNV=IME(1:IA)//'.'//UNV
         PAKNEU=IME(1:IA)//'.'//NEU
         ZSILE =IME(1:IA)//'.zs'
         IDUZIN=IA+4
      ELSE
         PAKLST=IME(1:IB-1)//'.'//LST
         PAKUNV=IME(1:IB-1)//'.'//UNV
         PAKNEU=IME(1:IB-1)//'.'//NEU
         ZSILE =IME(1:IB-1)//'.zs'
         IDUZIN=IB+3
      ENDIF
      END
