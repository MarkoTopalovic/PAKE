C==========================================================================
C==========================================================================
C   SUBROUTINE ULAZF3
C              TGRAFF
C              ISPITF
C              OTVGRF
C              OTVORF
C              OTVIZF
C              BRISF 
C              IMENAF
C              INPU01
C              INPU02
C              INPU03
C              INPU04
C              INPU05
C              INPU09
C              INPU10
C              INPU11
C              MEMORY
C              RESTAF
C              TIMFUN
C              AXISYF
C              WRRF
C              IWRRF
C              FGRAF3
C              TGRBCF
C              TGRAF2
C              ZATVOF
C              ADDSTF
C              PROMEM
C              TGRAFC
C              ULTAFF
C              PERIOF
C              KONVTF
C              MAXATF
C              INHEAD
C              HEADIN
C              REAINT
C              INTREA
C              UACTCF
C              ZFLUID
C              ZFLU
C              ZFLU1
C              IDENSF
C              WALLPS
C              SOLFLU
C              TUBE
C              RTUBE
C              RTUBEC
C              RTUBE1
C              RCAROT
CC             IDENTI
C              OUTFLU
C              OUTPAK
C              OUTPAF
C==========================================================================
C==========================================================================
      SUBROUTINE ULAZF3(NZAD,ZADVRE,NUMZAD,IULAZ,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /ULAZNI/ IULAZ,IIZLAZ

C
CE Subroutine ULAZF3 is used for reading input data for precribed values
C
      CHARACTER*250 ACOZ
      DIMENSION NZAD(3,*),ZADVRE(*)
      DO 28 I=1,NUMZAD
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1012) NZAD(1,I),NZAD(2,I),NZAD(3,I),ZADVRE(I)
  28  CONTINUE
C3000 FORMAT (//
C    111X,'GRESKA U ULAZNIM PODACIMA ???'//
C    111X,'U CVORU ',I5,' JE ZADATO OGRANICENJE U PRAVCU',I5///)
 1012 FORMAT(3I5,F10.3)
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,1101)
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,3101)
 1101 FORMAT(//
     16X,'BROJ CVORA, IND. PROMENLJIVE, VREMEN. FUNK., ZAD. VREDNOST'
     1/)
 3101 FORMAT(
     16X,'NUM. NODE, IND. ABOUT QUANTITIES, TIME FUNC. NUM., VALUE'
     1/)
      DO 32 I=1,NUMZAD
      WRITE(IIZLAZ,1013) NZAD(1,I),NZAD(2,I),NZAD(3,I),ZADVRE(I)
C      ZADVRE(I)=ZADVRE(I)*IFAKT
  32  CONTINUE
 1013 FORMAT (I20,2I10,F10.3)
      END
C==========================================================================
C==========================================================================
C==========================================================================
C=======================================================================
      SUBROUTINE TGRAFF(NEL,NBR2,NET,NGE,II,NETIP,NDIMM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /TIPELM/ NETIP
C      COMMON /TIPEL/ NP2DMX
	DIMENSION NEL(NDIMM,*)
C
CE Subroutine TGRAFF is used for printing basic data for elements
CE in file *.UNV which is maked for graphical postprocessing
C

      IF (NETIP.EQ.2) THEN
  	  CALL TGRAF2(NEL,NBR2,1,NET,NGE,II,ISRPS,NDIMM)  
        CALL TGRAU2F(NEL,NDIMM,NET,69,1) 
      ELSE IF (NETIP.EQ.3) THEN
  	  CALL FGRAF3(NEL,NBR2,NET,NGE,II,ISRPS,NDIMM)
        CALL TGRAU3F(NEL,NBR2,NET,69,1)  
      ENDIF   
      
      END
C=======================================================================
C======================================================================
      SUBROUTINE ISPITF(ACOZ,IULAZ)

C
CE Subroutine ISPITF is used for reading one line from input data
CE If 'C' on first column is found this line (data) is ignored 
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO TEST INPUT CARD
CS.    P R O G R A M
CS.        ZA ISPITIVANJE ULAZNE KARTICE
CS.        ZA BROJANJE UCITANIH KARTICA - KARTIC=KARTIC
C .
CS.        AKO U PRVE DVE KOLONE KARTICE STOJI: 'C ' ILI 'C-' ILI 'C*'
CS.        KARTICA SE IGNORISE
C .
C ......................................................................
C
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
       COMMON /BROJFK/ INDFOR,NULAZ
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C      COMMON /CDEBUG/ IDEBUG
      CHARACTER*250 ACOZ


C      IF(IDEBUG.GT.0) PRINT *, ' ISPITF'

C   10 KARTIC = KARTIC + 1
  10  READ(IULAZ,1000) ACOZ
      IF(ACOZ(1:2).EQ.'C '.OR.ACOZ(1:2).EQ.'C*'.OR.ACOZ(1:2).EQ.'C-')
     1 GO TO 10
      IF(INDFOR.EQ.2) RETURN
      BACKSPACE 1
      RETURN
 1000 FORMAT(A80)
      END

C.......................................................................
C.                                                                     .
C=======================================================================
C======================================================================
      SUBROUTINE OTVGRF(IME,PAKLST,PAKUNV,IDUZIN,IIUNV)

C
CE Subroutine OTVGRF is used for open file for graphic postprocessing
C
      COMMON /SRPSKI/ ISRPS
      CHARACTER*1 IME*20,STAT*3
       CHARACTER *24 PAKLST,PAKUNV
      LOGICAL OLDNEW
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OTVGRA'
C
CS IZLAZ ZA GRAFIKU
CE OUTPUT FOR GRAPHIC
C
   15 CONTINUE
      IF(ISRPS.EQ.0)
     *WRITE(*,*)' UNETI IME IZLAZNE DATOTEKE ZA GRAFIKU /"'
     1,PAKUNV(1:IDUZIN),'"'
      IF(ISRPS.EQ.1)
     *WRITE(*,*)' ENTER NAME OF OUTPUT FILE FOR GRAPHICS /"' 
     1,PAKUNV(1:IDUZIN),'"'
C   15 WRITE(*,820)
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = PAKUNV
  910 FORMAT (A)
C
   20 STAT='NEW'
      INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(OLDNEW) STAT='OLD'
      IF(STAT.EQ.'NEW') THEN
      OPEN (IIUNV,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ELSE
      OPEN (IIUNV,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ENDIF
C
      IND=0
      IF(STAT.EQ.'OLD') CALL BRISF (IME,IIUNV,IND)
      IF(IND.EQ.1)GO TO 20
      IF(IND.EQ.2)GO TO 15
      RETURN
      END
C=======================================================================
C==========================================================================
      SUBROUTINE OTVORF(ISRPS,IULAZ,IIZLAZ,IME,PAKLST,PAKUNV,PAKNEU,
     1           IDUZIN)
C
CE Subroutine OTVORF is used for open input file
C
      CHARACTER*1 IME*20
       CHARACTER *24 PAKLST,PAKUNV,PAKNEU
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OTVORI'
C
CS ULAZNA DATOTEKA
CE INPUT FILE
C
      IF(ISRPS.EQ.0)
     *WRITE(*,2000)
      IF(ISRPS.EQ.1)
     *WRITE(*,6000)
 2000 FORMAT('   UNETI IME ULAZNE DATOTEKE / "pakf.dat" :')
 6000 FORMAT('   ENTER NAME OF INPUT FILE / "pakf.dat" :')
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = 'pakf.dat            '
  910 FORMAT (A)
      CALL IMENAF(IME,PAKLST,PAKUNV,PAKNEU,IDUZIN)
      OPEN (IULAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
C
CS  SCRATCH DATOTEKE
CE  SCRATCH FILES
C
C      OPEN (15,STATUS='SCRATCH',FORM='UNFORMATTED',
C     1 ACCESS='SEQUENTIAL')
C      OPEN (19,FILE='ZITEMP',STATUS='UNKNOWN',FORM='UNFORMATTED',
C     1 ACCESS='SEQUENTIAL')
C      OPEN (20,STATUS='SCRATCH',FORM='UNFORMATTED',
C     1 ACCESS='SEQUENTIAL')
      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE OTVIZF(ISRPS,IIZLAZ,IME,PAKLST,PAKUNV,IDUZIN)
C
CE Subroutine OTVIZF is used for open file *.LST
C
C      COMMON /IMEF/ IME
C      COMMON /SRPSKI/ ISRPS
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C       COMMON /IMEULF/ PAKLST,PAKUNV,IDUZIN
      CHARACTER*1 IME*20,STAT*3
       CHARACTER *24 PAKLST,PAKUNV
      LOGICAL OLDNEW
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OTVIZL'
C
C
CS IZLAZNA DATOTEKA
CE OUTPUT FILE
C
    5 CONTINUE
      IF(ISRPS.EQ.1)
     *WRITE(*,*)'   ENTER NAME OF OUTPUT FILE /"',PAKLST(1:IDUZIN),'"'
      IF(ISRPS.EQ.0)
     *WRITE(*,*)'   UNETI IME IZLAZNE DATOTEKE /"',PAKLST(1:IDUZIN),'"'
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = PAKLST
  910 FORMAT (A)
C
   10 STAT='NEW'
      INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(OLDNEW) STAT='OLD'
      IF(STAT.EQ.'NEW') THEN
      OPEN (IIZLAZ,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ELSE
      OPEN (IIZLAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ENDIF
C
      IND=0
      IF(STAT.EQ.'OLD') CALL BRISF (IME,IIZLAZ,IND)
      IF(IND.EQ.1)GO TO 10
      IF(IND.EQ.2)GO TO 5
C
      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE BRISF (IME,IUN,IND)

C
CE Subroutine BRISF is used for deleting files
C
CS BRISANJE FILE-A
CE DELETE FILES
      CHARACTER*1 IME(20),CH*1
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' BRIS  '
      IND=2
      IF(ISRPS.EQ.0)
     *WRITE(*,2000) (IME(I),I=1,20)
      IF(ISRPS.EQ.1)
     *WRITE(*,6000) (IME(I),I=1,20)
 2000 FORMAT ('   FILE    ',20A1,' VEC POSTOJI'//
     1 '                 <ENTER>     PREBRISATI'/
     1 '                   "N"       OBICI ')
 6000 FORMAT ('   FILE    ',20A1,' ALREADY EXISTS'//
     1 '                 PRESS "ENTER" TO DELETE OR'/
     1 '                 KEY   "N"     TO BYPASS')
      READ(*,910) CH
  910 FORMAT(A)
      IF(CH.EQ.' ') CH = 'D'
      IF(CH.EQ.'D'.OR.CH.EQ.'d') THEN
                    CLOSE (IUN,STATUS='DELETE')
                    IND=1
                        ELSE
                        CLOSE (IUN,STATUS='KEEP')
                    ENDIF
      RETURN
      END
C======================================================================
C=======================================================================
       SUBROUTINE IMENAF(IME,PAKLST,PAKUNV,PAKNEU,IDUZIN)
C       COMMON /IMEULF/ PAKLST,PAKUNV,IDUZIN
C
CE Subroutine IMENAF is used for definition of names for files *.LST and *.UNV
C
       CHARACTER *20 IME
       CHARACTER *24 PAKLST,PAKUNV,PAKNEU
        IB=INDEX(IME,'.')
        DO 20 I=1,20
       IF (IME(I:I).EQ.' ') GOTO 30
        IA=I
  20    CONTINUE      
     
  30     IF (IB.EQ.0) THEN
          PAKLST=IME(1:IA)//'.lst'
          PAKUNV=IME(1:IA)//'.unv'
          PAKNEU=IME(1:IA)//'.neu'
        IDUZIN=IA+4
        ELSE
          PAKLST=IME(1:IB-1)//'.lst'
          PAKUNV=IME(1:IB-1)//'.unv'
          PAKNEU=IME(1:IB-1)//'.neu'
        IDUZIN=IB+3
        ENDIF

       
        END
C=======================================================================
C==========================================================================
      SUBROUTINE INPU01(IULAZ,NASLOV,IIZLAZ,ISRPS)
C
CE Subroutine INPU01 is used for reading head of problem from input data file
C

      CHARACTER*250 ACOZ
      CHARACTER*250 NASLOV
       
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1000) NASLOV
 1000 FORMAT(A80)

      IF (ISRPS.EQ.0) 
     *WRITE(IIZLAZ,1019)
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,1020)

      WRITE(IIZLAZ,3019) NASLOV

 1019 FORMAT(//6X,'N A S L O V   P R O B L E M A'/6X,29('-'))
 1020 FORMAT(//6X,'H E A D I N G  O F  P R O B L E M'/6X,33('-'))
 3019 FORMAT(//78('*')/,A80,/78('*')///)

      END
C=======================================================================
      SUBROUTINE INPU02(IULAZ,INDFOR,IIZLAZ,ISRPS)
C
CE Subroutine INPU02 is used for reading kind of format from input data file
C
      CHARACTER*250 ACOZ
      COMMON /INDNOVF/ INDNOV
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1001) INDFOR,INDNOV
      WRITE(*,*)'INDNOV',INDNOV
 1001 FORMAT(5I5)

      IF(ISRPS.EQ.0) THEN
       WRITE(IIZLAZ,1010) 
       WRITE(IIZLAZ,3020) INDFOR
      ELSEIF(ISRPS.EQ.1) THEN
       WRITE(IIZLAZ,3010) 
       WRITE(IIZLAZ,7020) INDFOR
      ENDIF

 1010 FORMAT(//6X,'F O R M A T  Z A  U L A Z N E  P O D A T K E'/6X,
     &44('-')//)
 3010 FORMAT(//6X,'F O R M A T  F O R  I N P U T  D A T A'/6X,38('-')//)

 3020 FORMAT(
     111X,'NACIN UCITAVANJA ULAZNIH PODATAKA ............. INDFOR =',I5/
     116X,'EQ.0; INDFOR = 1'/
     116X,'EQ.1; U SLOBODNOM FORMATU'/
     116X,'EQ.2; U OPISANOM FORMATU',///)
 7020 FORMAT(
     111X,'FORMAT INDIKATOR .............................. INDFOR =',I5/
     116X,'EQ.0; INDFOR = 1'/
     116X,'EQ.1; FREE FORMAT'/
     116X,'EQ.2; FIXED FORMAT, AS DESCRIBED',///)

      END
C=======================================================================
      SUBROUTINE INPU03 (IULAZ,NPT,NSTAC,NPER,INDFL,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine INPU03 is used for reading basic data for problem
CE from input data file
C
      CHARACTER*250 ACOZ
      COMMON /INDNOVF/ INDNOV
      CALL ISPITF(ACOZ,IULAZ)
      IF(INDNOV.EQ.0)
     1READ(ACOZ,1002) NPT,NGET,NMATT,NSTAC,NPER,NPRINT
      IF(INDNOV.EQ.1)
     1READ(ACOZ,6002) NPT,NGET,NMATT,NSTAC,NPER,NPRINT
 1002 FORMAT(15I5)
 6002 FORMAT(I10,14I5)
      IF(NPER.EQ.0) NPER = 1
      IF(NGET.EQ.0) NGET = 1
      IF(NMATT.EQ.0) NMATT=1
      IF(NPRINT.EQ.0) NPRINT=1

      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2001) NPT,NGET,NMATT,NSTAC,NPER,NPRINT
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6001) NPT,NGET,NMATT,NSTAC,NPER,NPRINT,INDFL
 2001 FORMAT(6X,'O S N O V N I    P O D A C I    O    P R O B L E M U'
     1/6X,51('-')///
     111X,'UKUPAN BROJ CVORNIH TACAKA ....................... NP =',I10/
     116X,'EQ.0; PREKIDA SE IZVRSAVANJE PROGRAMA'///
     211X,'BROJ GRUPA ELEMENATA ............................ NGET =',I5/
     216X,'EQ.0; POSTAJE "1"; (MAX. 10 GRUPA)'///
     311X,'BROJ RAZLICITIH MATERIJALA ..................... NMATT =',I5/
     316X,'EQ.0; POSTAJE "1"'///
     411X,'INDIKATOR STACIONARNOSTI ....................... NSTAC =',I5/
     416X,'EQ.1; STACIONARAN PROBLEM'/
     416X,'EQ.0; NESTACIONARAN PROBLEM'///
     511X,'BROJ PERIODA SA KONSTANTNIM VREMENSKIM KORACIMA . NPER =',I5/
     516X,'EQ.0; POSTAJE "1"'///
     611X,'DEFINISANJE STAMPARSKOG KORAKA ................ NPRINT =',I5/
     616X,'EQ.0; POSTAJE "1"')
 6001 FORMAT(6X,'B A S I C    D A T A    F O R   T H E   P R O B L E M'
     1/6X,53('-')///
     111X,'TOTAL NUMBER OF NODAL POINTS ..................... NP =',I10/
     116X,'EQ.0; PROGRAM STOP'///
     211X,'NUMBER OF ELEMENT GROUPS ........................ NGET =',I5/
     216X,'EQ.0; DEFAULT SET "1"'///
     311X,'NUMBER OF DIFERENT MATERIALS ................... NMATT =',I5/
     316X,'EQ.0; DEFAULT SET "1"'///
     411X,'STEADINES INDICATOR ............................ NSTAC =',I5/
     416X,'EQ.1; STEADY STATE'/
     416X,'EQ.0; TRANSIENT'///
     511X,'NUMBER OF CONSTANT TIME STEP PERIODS ............ NPER =',I5/
     516X,'EQ.0; DEFAULT SET "1";'///
     611X,'OUTPUT PRINTING INTERVAL ...................... NPRINT =',I5/
     616X,'EQ.0; DEFAULT SET "1"'///
     611X,'OUTPUT PRINTING INTERVAL ...................... NPRINT =',I5/
     616X,'EQ.0; DEFAULT SET "1"')



      END
C==========================================================================
C=======================================================================
      SUBROUTINE INPU04(IULAZ,METOD,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP,
     &MBAND,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine INPU04 is used for reading basic data for problem
CE from input data file
C
      CHARACTER*250 ACOZ
      COMMON /INDNOVF/ INDNOV
      CALL ISPITF(ACOZ,IULAZ)
      
      READ(ACOZ,1003) METOD,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP,MBAND
 1003 FORMAT(4I5,2F10.2,2I5)
      IF(METOD.EQ.0) METOD=1
      IF(NJRAP.EQ.0) NJRAP=1
      IF(MAXIT.EQ.0) MAXIT=15
      IF(DABS(EPSTA).LT.1.D-10.AND.DABS(EPSTR).LT.1.D-10) EPSTR=.001
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2002) METOD,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP,MBAND
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6002) METOD,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP,MBAND
 2002 FORMAT(6X,'O S N O V N I    P O D A C I    O    P R O B L E M U'
     1/6X,51('-')///
     111X,'PRIMENJEN METOD VREMENSKE INTEGRACIJE .......... METOD =',I5/
     116X,'EQ.0; POSTAJE "1"'/
     116X,'EQ.1; EULER BACWARD INTEGRACIJA (ALFA = 1)'/
     116X,'EQ.2; NE KORISTI SE'///
     211X,'STAMPANJE REZULTATA U ZELJENIM CVOROVIMA ....... INDSC =',I5/
     216X,'EQ.0; U SVIM CVOROVIMA'/
     216X,'EQ.1; U ZELJENIM CVOROVIMA'///
     311X,'STAMPANJE REZULTATA U ZELJENOM FORMATU ......... IFORM =',I5/
     316X,'EQ.0; U FORMATU  D13.5'/
     316X,'EQ.1; U FORMATU  F10.3'///
     411X,'MAXIMALNI BROJ RAVNOTEZNIH ITERACIJA ........... MAXIT =',I5/
     416X,'EQ.0; POSTAJE "15"'///11X,
     5'APSOLUTNA TACNOST PRI ITERACIJAMA ......... EPSTA =',1PD10.3///
     611X,'RELATIVNA TACNOST PRI ITERACIJAMA ......... EPSTR =',1PD10.3/
     616X,'EQ.0; POSTAJE "1.E-3"'///
     111X,'PRIMENJEN ITERATIVNI METOD ..................... NJRAP =',I5/
     116X,'EQ.0; POSTAJE "1"'/
     116X,'EQ.1; MODIFIKOVAN NJUTN-RAPSONOV METOD'/
     116X,'EQ.2; POTPUN NJUTNOV'///
     111X,'METOD DEFINISANJA PAKOVANJA JEDNACINA........... MBAND =',I5/
     116X,'EQ.0; JEDNACINE SE PAKUJU REDOM '/
     116X,'EQ.1; JEDNACINE PRITISKA SE PAKUJU NA KRAJU SISTEMA'///)
 6002 FORMAT(6X,'B A S I C    D A T A    F O R   T H E   P R O B L E M'
     1/6X,53('-')///
     111X,'TIME INTEGRATION METHOD USED ................... METOD =',I5/
     116X,'EQ.0; DEFAULT SET "1"'/
     116X,'EQ.1; EULER BACKWARD INTEGRATION (ALPFA = 1)'///
     211X,'PRINT OF RESULTS IN PRESCRIBED NODES ........... INDSC =',I5/
     216X,'EQ.0; IN ALL NODES'/
     216X,'EQ.1; IN PRESCRIBED NODES'///
     311X,'PRINT OF RESULTS IN PRESCRIBED FORMAT .......... IFORM =',I5/
     316X,'EQ.0; IN FORMAT  D13.5'/
     316X,'EQ.1; IN FORMAT  F10.3'///
     411X,'MAXIMUM NUMBER OF ITERATIONS ................... MAXIT =',I5/
     416X,'EQ.0; DEFAULT SET "15"'///11X,
     5'APSOLUTLY ACCURANCU AT ITERATIONS ......... EPSTA =',1PD10.3///
     611X,'RELATIVELY ACCURANCU AT ITERATIONS ........ EPSTR =',1PD10.3/
     616X,'EQ.0; DEFAULT SET "1.E-3"'///
     111X,'ITERATION METHODS EMPLOYED ..................... NJRAP =',I5/
     116X,'EQ.0; DEFAULT SET "1"'/
     116X,'EQ.1; MODIFIED NEWTON'/
     116X,'EQ.2; FULL NEWTON'///
     111X,'DEFINITION EQUATION METHOD EMPLOYED............. MBAND =',I5/
     116X,'EQ.0; EQUATIONS ARE DEFINED SUCCESSIVELY '/
     116X,'EQ.1; EQUATIONS OF PRESSURE ARE DEFINED AT THE END OF'/
     116X,'      GLOBAL SYSTEM'///)


      END
C=======================================================================
      SUBROUTINE INPU05 (IULAZ,IREST,IIZLAZ,ISRPS)
      CHARACTER*250 ACOZ
C
CE Subroutine IREST is used for reading data about execution of problem
C

      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1004) IREST
 1004 FORMAT(I5)
      IF (ISRPS.EQ.0) THEN
       WRITE(IIZLAZ,4001)IREST
      ELSE
       WRITE(IIZLAZ,4002)IREST
      ENDIF
 4001 FORMAT(6X,'P O D A C I   Z A   R E S T A R T O V A NJ E'
     1/6X,44('-')///
     111X,'INDIKATOR IZVRSENJA PROGRAMA ................... IREST =',I5/
     116X,'EQ.0; KONTROLA ULAZNIH PODATAKA'/
     116X,'EQ.1; IZVRSENJE PROGRAMA'//)
 4002 FORMAT(6X,'D A T A   F O R   R E S T A R T'
     1/6X,31('-')///
     111X,'INDICATOR FOR JOB EXECUTION ................... IREST =',I5/
     116X,'EQ.0; CHECK OF INPUT DATA'/
     116X,'EQ.1; JOB EXECUTION'//)

      END
C=======================================================================
C==========================================================================
      SUBROUTINE INPU09(IULAZ,POCU,POCV,POCP,POCT,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine INPU09 is used for reading basic data for initial values
C
      CHARACTER*250 ACOZ

      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1014) POCU,POCV,POCP,POCW,POCT
 1014 FORMAT(5F10.3)

      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,3008)POCU,POCV,POCP,POCW,POCT
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6008)POCU,POCV,POCW,POCP,POCT

 3008 FORMAT(6X,'P O D A C I   O   P O C E T N I M   V R E D N O S T I M
     1 A'/6X,58('-')///
     111X,'POCETNA BRZINA U PRAVCU OSE X ................U= ',F10.3//
     111X,'POCETNA BRZINA U PRAVCU OSE Y ................V= ',F10.3//
     111X,'POCETNA BRZINA U PRAVCU OSE Z ................V= ',F10.3//
     111X,'POCETNI PRITISAK .............................P= ',F10.3//
     111X,'POCETNA TEMPERATURA ..........................T= ',F10.3///)

 6008 FORMAT(///6X,'D A T A   A B O U T   I N I T I A L   V A L U E S'/
     16X,49('-')///
     111X,'INITIAL FLUID VELOCITY Uo ................. Uo = ',F10.3//
     111X,'INITIAL FLUID VELOCITY Vo ................. Vo = ',F10.3//
     111X,'INITIAL FLUID VELOCITY Wo ................. Wo = ',F10.3//
     111X,'INITIAL PRESSURE ........................... P = ',F10.3//
     111X,'INITIAL TEMPERATURE ........................ T = ',F10.3///)
   
      END
C=======================================================================
C==========================================================================
      SUBROUTINE INPU10(IULAZ,NTABFT,MAXTFT,IIZLAZ,ISRPS)
      CHARACTER*250 ACOZ
C
CE Subroutine INPU10 is used for reading basic data for time functions
C

      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1002) NTABFT,MAXTFT

      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,2001) NTABFT,MAXTFT
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6001) NTABFT,MAXTFT

 2001 FORMAT(6X,'P O D A C I   O   V R E M E N S K I M   F U N K C I J A
     1 M A'/6X,59('-')///
     111X,'UKUPAN BROJ VREMENSKIH FUNKCIJA ............. NTABFT =',I5/
     111X,'MAKSIM. BR. TACAKA ZA VREM. FUNK. ........... MAXTFT =',I5//)
 6001 FORMAT(6X,'D A T A   A B O U T   T I M E   F U N C T I O N S'/
     16X,49('-')///
     111X,'NUMBER OF DIFFERENT TIME FUNCTIONS .......... NTABFT =',I5//
     111X,'MAX. NUM. OF POINTS FOR TIME FUNCTIONS ...... MAXTFT =',I5//)

 1002 FORMAT(15I5)

      END
C=======================================================================
C==========================================================================
      SUBROUTINE INPU11(IULAZ)
C
CE Subroutine INPU03 is used for reading end of input file
C
      CHARACTER*250 ACOZ
      CHARACTER*80 KRAJ

      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1000) KRAJ
      IF(KRAJ(1:4).NE.'STOP') WRITE(*,*)'Missing end of file'

 1000 FORMAT(A80)

      END
C=======================================================================
C==========================================================================
      SUBROUTINE MEMORY(LSTART,LMAX,NBYTE,ITWO,MAXVEC,IIZLAZ)
C
CE Subroutine MEMORY is used for reservation of memory 
C
      
      IF (MOD(LMAX,2).EQ.0) LMAX=LMAX+1
      LSTART=LMAX
      LMAX=LSTART+NBYTE*ITWO

      CALL PROMEM(LMAX,MAXVEC,IIZLAZ)

      END
C=======================================================================
C==========================================================================
      SUBROUTINE RESTAF (IREST,IIZLAZ,ISRPS)
C
CE Subroutine RESTAF is used for check input data analysis
C

      IF (IREST.EQ.0) THEN 
      IF (ISRPS.EQ.1) THEN 
        WRITE(IIZLAZ,1000)
      ELSE
        WRITE(IIZLAZ,2000)
      ENDIF
 1000 FORMAT(//
     111X,'ZAVRSENA JE KONTROLA ULAZNIH PODATAKA'/
     111X,'ZA PRORACUN JE POTREBNO STAVITI IREST=1'//)
 2000 FORMAT(//
     111X,'The control for input data is finished'/
     111X,'If you want CFD analysis you must input IREST=1'//)
      STOP      
      ENDIF
      
      END
C=======================================================================
C==========================================================================
      SUBROUTINE TIMFUN(TABF,SILAF,UVREME,IBRT,IFUN,NTABFT,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine TIMFUN is used for definition variable in time
C

C      COMMON /VREPER/ NPER,NTABFT
C      COMMON /VREPER/ NPER,NASLOV,IBRT
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C
      DIMENSION TABF(2,NTABFT,*)

      TOL=1.D-20
      IF (UVREME.LT.TOL) UVREME=1.D-20
       VREME=UVREME+UVREME/1.D10

      DO I=1,IBRT-1
      IF((VREME.GT.TABF(1,IFUN,I)).AND.(TABF(1,IFUN,I+1).GT.VREME)) THEN
        SRT=TABF(1,IFUN,I+1)-TABF(1,IFUN,I) 
        SRF=TABF(2,IFUN,I+1)-TABF(2,IFUN,I) 
        DELTAT=UVREME-TABF(1,IFUN,I)
        SILAF=TABF(2,IFUN,I)+(DELTAT/SRT)*SRF
        RETURN
       ENDIF
      ENDDO
C       WRITE(IIZLAZ,*)'GRESKA U ULAZNIM PODACIMA ZA VREMENSKU FUNKCIJU'
C       WRITE(IIZLAZ,*)'VREMENSKI KORAK JE IZVAN OPSEGA VREMENSKE 
C     1FUNKCIJE'
       WRITE(IIZLAZ,*)'WRONG DATA INPUT FOR TIME FUNCTION'
       WRITE(IIZLAZ,*)'TIME STEP OVER SIZE OF TIME FUNCTION'
      END
C==========================================================================
C==========================================================================
      SUBROUTINE AXISYF(INDAX,DEBLJ,X,H,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION X(*),H(*)

C
CE Subroutine AXISYF is used for axi-symmetric analysis
C

      IF (INDAX.EQ.0) THEN
        DEBLJ=1.D0
        RETURN
      ELSE IF (INDAX.EQ.1) THEN
      DEBLJ=0.D0
      DO I=1,NDIM
        DEBLJ=DEBLJ+H(I)*X(I)
      ENDDO
      ENDIF

      END
C=======================================================================
      SUBROUTINE WRRF(A,N,CHAR,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
      COMMON /CDEBUG/ IDEBUG
C
CE Subroutine WRRF is used for printing real data to output file *.LST
C

C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
      CHARACTER*5 CHAR
      DIMENSION A(*)
      IF(IDEBUG.GT.0) PRINT *, ' WRR   '
C
      WRITE(IIZLAZ,5010) CHAR
      WRITE(IIZLAZ,5000) (A(I),I=1,N)
      RETURN
C
 5010 FORMAT(A5)
 5000 FORMAT(4(1PD18.9))
      END
C=======================================================================
      SUBROUTINE IWRRF(M,N,CHAR,IIZLAZ)
C
CE Subroutine IWRRF is used for printing integer data to output file *.LST
C
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE INTEGER VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE CELOBROJNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
C
      CHARACTER*5 CHAR
      DIMENSION M(*)
C      DIMENSION M(N)
C      COMMON /CDEBUG/ IDEBUG
C      IF(IDEBUG.GT.0) PRINT *, ' IWRRF  '
C
      WRITE(IIZLAZ,5010) CHAR
 
C      DO 10 I=1,N
C      WRITE(IIZLAZ,5000) M(I)
C  10  CONTINUE
     
      WRITE(IIZLAZ,5000) (M(I),I=1,N)
      RETURN
C
 5010 FORMAT(A5)
 5000 FORMAT(14I5)
      END
C=======================================================================
C==========================================================================
C=======================================================================
      SUBROUTINE FGRAF3(NEL,NBR2,NET,NGE,II,ISRPS,NDIMM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine FGRAF3 is used for printing data for finite elements 
CE to output file *.UNV
C

      DIMENSION NEL(NDIMM,*)
C
C     E L E M E N T I   3/D
C
C      IF(ISRPS.EQ.0.AND.(NBR2.NE.8.AND.NBR2.NE.20))
C     1WRITE(II,2200) NGE
C      IF(ISRPS.EQ.1.AND.(NBR2.NE.8.AND.NBR2.NE.20))
C     1WRITE(II,6200) NGE
      NBR=NBR2
      IF(NBR2.LT.20) NBR2=8
      IF(NBR2.EQ.21) NBR2=21
C     GRAFICKI OPIS ELEMENTA: SA 8 CVOROVA = 19, SA 20 CVOROVA = 20
      ITYPE=19
C-----------------------------------------------------------------------
C      IF(NBR2.EQ.20) ITYPE=20
C 04.02.2011.
C-----------------------------------------------------------------------
            IF(NBR2.EQ.21) ITYPE=19
C     VRSTA 3/D ELEMENTA: 
      IE1=115
C----------------------------------------------------------------------- 
C     IF(NBR2.EQ.20) IE1=115
C 04.02.2011.
C-----------------------------------------------------------------------
      IF(NBR2.EQ.21) IE1=115
C     TABELA FIZICKIH OSOBINA
      IE2=1
C     TABELA MATERIJALA
      IE3=1
C     BOJA  
      ICOL=8
      IND=-1
      ITYP=71
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 10 I=1,NET
         WRITE(II,1000) I,ITYPE,IE1,IE2,IE3,ICOL,NBR2
C----------------------------------------------------------------------- 
C 04.02.2011.
C-----------------------------------------------------------------------
         IF(NBR2.EQ.8) THEN
            WRITE(II,1000) (NEL(J,I),J=1,8)
         ELSE
           WRITE(II,1000) (NEL(J,I),J=1,8),(NEL(J,I),J=9,16),
     1      (NEL(J,I),J=17,21)
C-----------------------------------------------------------------------
C           WRITE(II,1000) (NEL(J,I),NEL(J+8,I),J=1,4),
C     1                    (NEL(J,I),J=17,21),(NEL(J,I),NEL(J+8,I),J=5,8) 
         ENDIF
   10 CONTINUE
      WRITE(II,1100) IND

      NBR2=NBR
      RETURN
C
 1100 FORMAT(I6)
 1000 FORMAT(8I10)
C-----------------------------------------------------------------------
 2200 FORMAT(//' PROGRAM ZA GRAFICKO PRIKAZIVANJE REZULTATA "IDEAS"'/
     1' ZAHTEVA 3/D ELEMENT SA 8 ILI 20 CVOROVA U GRUPI ELEMENATA NGE ='
     1,I5)
C-----------------------------------------------------------------------
 6200 FORMAT(//' GRAPHIC PACKAGE   "IDEAS"'/
     1' PERMITS ONLY 3/D ELEMENTS WITH 8 OR 20 NODES PER ELEMENT IN',
     1' GROUP   NGE =',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C======================================================================
      SUBROUTINE TGRBCF(CORD,NPT,II,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine TGRBCF is used for printing boundary conditions 
CE data to output file *.UNV
C

      DIMENSION CORD(3,*),ID(5,*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' TGRBFC'
C      REWIND II
      IND=-1
      ITYP=757
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      WRITE(II,1000) 2
      WRITE(II,*)' DOF SET         2'
 1100 FORMAT(I6)
      IT1=0
      IT2=0
      ICOL=8
      DO 10 I=1,NPT
      I1=0
      I2=0
      I3=0
      I4=0
      I5=0
      IF (ID(1,I).EQ.0) I1=1
      IF (ID(2,I).EQ.0) I2=1
      IF (ID(3,I).EQ.0) I3=1
      IF (ID(4,I).EQ.0) I4=1
      IF (ID(5,I).EQ.0) I5=1
      WRITE(II,1000) I,7,I1,I2,I3,I4,I5,0
 1000 FORMAT(2I10,6I2)
   10 CONTINUE
      WRITE(II,1100) IND
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE TGRAF2(NEL,NBR2,ISTART,IEND,NGE,II,ISRPS,NDIMM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine TGRAF2 is used for printing data about 2D finite elements 
CE to graphical file *.UNV
C
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /TIPELM/ NETIP

C      COMMON /TIPEL/ NP2DMX

C      CHARACTER*250 ACOZ

C      COMMON /SRPSKI/ ISRPS

      DIMENSION NEL(NDIMM,*)
C
C     E L E M E N T I   2/D
C
C      IF(ISRPS.EQ.0.AND.(NBR2.NE.4.AND.NBR2.NE.8))
C     1WRITE(IIZLAZ,2200) NGE
C      IF(ISRPS.EQ.1.AND.(NBR2.NE.4.AND.NBR2.NE.8))
C     1WRITE(IIZLAZ,6200) NGE


      NBR=NBR2
      IF(NBR2.LT.8) NBR2=4
      IF(NBR2.EQ.9) NBR2=8
C     GRAFICKI OPIS RAVANSKOG ELEMENTA: SA 4 CVORA = 27, SA 8 CVOROVA = 28
      ITYPE=27 
      IF(NBR2.EQ.8) ITYPE=28
C     VRSTA 2/D ELEMENTA: 
      IE1=44
      IF(NBR2.EQ.8) IE1=45
C     TABELA FIZICKIH OSOBINA
      IE2=1
C     TABELA MATERIJALA
      IE3=1
C     BOJA  
      ICOL=8
      IND=-1
      ITYP=71
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 10 I=ISTART,IEND
         WRITE(II,1000) I,ITYPE,IE1,IE2,IE3,ICOL,NBR2
         IF(NBR2.EQ.4) THEN
            WRITE(II,1000) (NEL(J,I),J=1,4)
         ELSE
            WRITE(II,1000) (NEL(J,I),NEL(J+4,I),J=1,4)
         ENDIF
   10 CONTINUE
      WRITE(II,1100) IND
      NBR2=NBR
      RETURN
C
 1100 FORMAT(I6)
 1000 FORMAT(8I10)
C-----------------------------------------------------------------------
C 2200 FORMAT(//' PROGRAM ZA GRAFICKO PRIKAZIVANJE REZULTATA "IDEAS"'/
C     1' ZAHTEVA 2/D ELEMENT SA 4 ILI 8 CVOROVA U GRUPI ELEMENATA NGE ='
C     1,I5)
C-----------------------------------------------------------------------
C 6200 FORMAT(//' GRAPHIC PACKAGE   "IDEAS"'/
C     1' PERMITS ONLY 2/D ELEMENTS WITH 4 OR 8 NODES PER ELEMENT IN',
C     1' GROUP   NGE =',I5)
C-----------------------------------------------------------------------
C=======================================================================
C=======================================================================
      END
C==========================================================================
      SUBROUTINE ZATVOF(IPRBR,INDIZL,INDGRA,IIZLAZ)
C
CE Subroutine ZATVOF is used for close files
C
CS ZATVARANJE DATOTEKA
CE CLOSE FILES
C      COMMON /PRIMER/ IPRBR,INDIZL,INDGRA
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' ZATVOF'
C
      CLOSE (1,STATUS='KEEP')
      IF(INDIZL.EQ.1) CLOSE (IIZLAZ,STATUS='KEEP')
      IF(INDGRA.EQ.1) THEN
        CLOSE (18,STATUS='KEEP')
      ENDIF
      RETURN
      END
C======================================================================

C
C======================================================================
      SUBROUTINE ADDSTF(A,B,C,S,P,JDIAG,LD,NEL,INDSK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ADDSTF is used for assemble global arrays
C

C
C.... ASSEMBLE GLOBAL ARRAYS
C
      DIMENSION A(*),B(*),JDIAG(*),P(*),S(NEL,*),LD(*),C(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' ADDSTF'
      DO 200 J = 1,NEL
      K = LD(J)
      IF(K.EQ.0) GO TO 200
      B(K) = B(K) + P(J)
      IF(INDSK.EQ.0) GO TO 200
      L = JDIAG(K) - K
      DO 100 I = 1,NEL
      M = LD(I)
      IF(M.GT.K.OR.M.EQ.0) GO TO 100
      M = L + M
      A(M) = A(M) + S(I,J)
      C(M) = C(M) + S(J,I)
  100 CONTINUE
  200 CONTINUE
      RETURN
      END    
C=======================================================================
       SUBROUTINE PROMEM(MEMNEW,MAXVEC,IIZLAZ)
C      COMMON /RADNIV/ MAXVEC,LMAX
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C
CE Subroutine PROMEM is used for checking memory space
C        

C    WRITE(1,*)'NTOTF',NTOTF
C    WRITE(1,*)'MAXVEC',MAXVEC

       IF (MEMNEW.GT.MAXVEC) THEN
       WRITE(*,*) 'NEMA DOVOLJNO MEMORIJE ZA RESAVANJE PROBLEMA!!!'
      WRITE(IIZLAZ,*)'NIJE DOVOLJNO ZA RADNI VEKTOR ',MAXVEC,' MEMORIJE'
       WRITE(IIZLAZ,*) 'POTREBNO JE MEMORIJE VISE OD: ',MEMNEW
       STOP
       ENDIF

       END
C=======================================================================


C======================================================================
      SUBROUTINE TGRAFC(CORD,NPT,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,MESH,NUMALV,DISTAL,DIS1AL,RE,
     &SVEXC,VALST,PERIOD,TOL,NUMST
C
CE Subroutine TGRAFC is used for printing coordinates of nodes
CE to output graphical file *.UNV
C

      DIMENSION CORD(3,*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' TGRAFC'
      REWIND II
      IND=-1
      ITYP=15
C=================================================================     
C FOR AKIRA'S CALCULATIONS
C=================================================================     
      IF (PERIOD.GT.0.D0) THEN
       WRITE(II,2000) RD,RA,ALEN,GAMA,MESH,NUMALV,DISTAL
       WRITE(II,2001) DIS1AL,RE,SVEXC,NUMST,VALST,PERIOD,TOL
      ENDIF
 2000  FORMAT(4F10.7,2I5,F10.7)
 2001  FORMAT(3F10.7,I5,3F10.7)
C=================================================================     
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
 1100 FORMAT(I6)
      IT1=0
      IT2=0
      ICOL=8
      DO 10 I=1,NPT
      WRITE(II,1000) I,IT1,IT2,ICOL,(CORD(J,I),J=1,3)
C SAMO PRIVREMENO ZA PRIMER KONVEKTIVNE DIFUZIJE
C      WRITE(II,1000) I,IT1,IT2,ICOL,CORD(1,I)*10.D0,CORD(2,I)/2.D0,
C     &CORD(3,I)
 1000 FORMAT(4I10,3E13.5)
   10 CONTINUE
      WRITE(II,1100) IND
      RETURN
      END
C=======================================================================
C==========================================================================
      SUBROUTINE ULTAFF (TABF,ITFMAX,NPER,NTABFT,IULAZ,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ULTABF is used for reading input data for time functions
C
C      COMMON /VREPER/ NPER,NASLOV,IBRT
C      COMMON /VREPER/ NPER,NTABFT
C      COMMON /ULAZNI/ IULAZ,IIZLAZ

      CHARACTER*250 ACOZ

      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)



      DO I=1,NTABFT
      CALL ISPITF(ACOZ,IULAZ)
       READ(ACOZ,1000) IBR,IMAX
       IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,1001) IBR,IMAX
       IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000) IBR,IMAX

       ITFMAX(IBR)=IMAX
       DO J=1,IMAX
      CALL ISPITF(ACOZ,IULAZ)
       READ(ACOZ,1004) TABF(1,IBR,J),TABF(2,IBR,J)

       IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,1002) TABF(1,IBR,J),TABF(2,IBR,J)
       IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,4002) TABF(1,IBR,J),TABF(2,IBR,J)
       ENDDO
      ENDDO

 1004 FORMAT (2D10.3)
 1000 FORMAT (2I5)
 1002 FORMAT(//
     111X,'VREME (ARGUMENT).................................t =',D10.3/
     116X,'VREDNOST VREMENSKE FUNKCIJE...................f(t) =',D10.3
     1//)
 4002 FORMAT(//
     111X,'ARGUMENT (TIME) ............................ t =',D10.3/
     116X,'FUNCTION ................................ f(t) =',D10.3//)
 1001 FORMAT(//
     111X,'REDNI BROJ VREMENSKE FUNKCIJE ................... IFSC =',I5/
     116X,'MAKSIM. BR. TACAKA ZA DATU VREM. FUNK. .......... IMAX =',I5
     1//)
 6000 FORMAT(//
     111X,'FUNCTION SEQUENCE NUMBER ........................ IFSC =',I5/
     116X,'TOTAL NUMBER OF PAIRS ARGUMENT-FUNCTION ......... IMAX =',I5
     1//)
C5000 FORMAT(//
C    111X,'UKUPAN BROJ VREMENSKIH FUNKCIJA (NTABFT) ............. =',I5/
C    116X,'MAKSIM. BR. TACAKA ZA VREM. FUNK. (MAXTFT) ........... =',I5
C    1//)


      END
C==========================================================================
C==========================================================================
      SUBROUTINE PERIOF(VREME,NPER,NTABF,IULAZ,IIZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine PERIOF is used for reading data for time steps
C
C      COMMON /VREPER/ NPER,NTABFT
C      COMMON /VDP/ NKORP,DT,NN,NZAV
      CHARACTER*250 ACOZ

C      COMMON /ULAZNI/ IULAZ,IIZLAZ
      DIMENSION VREME(*),NTABF(*)

      IF (ISRPS.EQ.0) THEN
       WRITE(IIZLAZ,102)
      ELSE
       WRITE(IIZLAZ,104)
      ENDIF

 102  FORMAT(6X,'P O D A C I  Z A  V R E M E N S K E  K O R A K E'/6X,
     &48('-')//)
 104  FORMAT(6X,'D A T A  F O R  T I M E  S T E P S'/6X,
     &48('-')//)

      DO I=1,NPER
       NTABF(I)=0
       VREME(I)=0.D0
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1004) NKORP,DT
 1004 FORMAT (I5,D10.3)
      IF (ISRPS.EQ.0) THEN
       WRITE(IIZLAZ,3002)I
       WRITE(IIZLAZ,4002)NKORP,DT
      ELSE
       WRITE(IIZLAZ,3004)I
       WRITE(IIZLAZ,4004)NKORP,DT
      ENDIF
 3002 FORMAT(6X,'PERIOD BROJ..... NPER=',I5/)
 3004 FORMAT(6X,'PERIOD NUMBER..... NPER=',I5/)
 4002 FORMAT(
     111X,'BROJ KORAKA PO PERIODU.......................... NKORP =',I5/
     116X,'VREDNOST KORAKA...................DT  =',D10.3
     1//)
 4004 FORMAT(
     111X,'NUMBER OF TIME STEPS PER PERIOD................. NKORP =',I5/
     116X,'TIME IN STEP...........DT  =',D10.3
     1//)
        NTABF(I)=NKORP
        VREME(I)=DT
      ENDDO



      END
C==========================================================================
C=======================================================================

C======================================================================
      SUBROUTINE KONVTF(TT1,F,KONV,N,ID,ITER,NPT,EPSTR,ISRPS,IDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine KONVTF is used for criterium of convergence
C
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /TACNOS/ EPSTR,MAXIT
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /SRPSKI/ ISRPS

      DIMENSION TT1(*),F(*),ID(IDIM,*)
      KONV=1

      DTMOD=0.
      TMAX=0.
   
      DO 10 J=1,NPT
      I=ID(N,J)
      IF (I.EQ.0) GO TO 10
      TMAX=TMAX+TT1(I)**2
      DTMOD=DTMOD+F(I)**2
 10   CONTINUE       

      IF (DTMOD.GT.1.D-10) DTMOD=DSQRT(DTMOD)
      IF (TMAX.GT.1.D-10) TMAX=DSQRT(TMAX)

      TOLT=TMAX*EPSTR

      IF (DTMOD.GT.TOLT) KONV=0
      IF (TMAX.LT.1.D-10) THEN
      KONV=1
      ELSE 
      ENE=DTMOD/TMAX
      IF(ISRPS.EQ.0)
     *WRITE(*,1000) ITER,DTMOD,ITER,ENE
      IF(ISRPS.EQ.1)
     *WRITE(*,2000) ITER,DTMOD,ITER,ENE
      ENDIF
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / NORMA ( 0 )  = ',1PD12.4)
 2000 FORMAT(/'      ITERATION = ',I5,5X,'NORM   = ',1PD12.4/
     1        '          NORM (',I3,') / NORM ( 0 )  = ',1PD12.4)
      END
C======================================================================
C==========================================================================
      SUBROUTINE MAXATF(MAXA,MHT,ID,NEL,NE,NTE,JEDN,NWK,IID,NDIMM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine MAXATF is used 
CE for determination of column heights vector and maxa
C
C
CS    PODPROGRAM ZA FORMIRANJE VEKTORA VISINA STUBOVA I MAXA
CS    KONACNO SE SMESTAJU U ISTI PROSTOR
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
C
C      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,HT,DETJ,
C     1DETJS,X,Y,FS2,FS3
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN1,NWK1,NET
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C
      DIMENSION MAXA(*),MHT(*),NEL(NDIMM,*),LM(92),ID(IID,*)
      COMMON /CDEBUG/ IDEBUG

      IF(IDEBUG.GT.0) WRITE(*,*) 'MAXATE'
     
      DO I=1,JEDN
       MHT(I)=0
       MAXA(I)=0
      ENDDO
C
CS    PETLJA PO ELEMENTIMA
CE    ELEMENT LOOP
C
      DO 100 NLM=1,NE
         KK=0
         DO 2 I=1,NTE
            IF(NEL(I,NLM).EQ.0) GO TO 2
            N=NEL(I,NLM)
               DO 1 J=1,IID
                  IF(ID(J,N).LE.0) GO TO 1
                  KK=KK+1
                  LM(KK)=ID(J,N)
    1          CONTINUE
C            ENDIF
    2    CONTINUE
C         CALL IWRR(LM,KK,'   LM')
C
         LS=JEDN+1
         DO 10 I=1,KK
            IF (LM(I).LT.LS) LS=LM(I)
   10    CONTINUE

C         WRITE(IIZLAZ,*)'LS=',LS
C
         DO 20 I=1,KK
            II=LM(I)
            ME=II-LS
            IF(ME.GT.MHT(II)) MHT(II)=ME
C         WRITE(IIZLAZ,*)'MHT(II),II',MHT(II),II
   20    CONTINUE
C
  100 CONTINUE
C
CS    VEKTOR MAXA
CE    VECTOR MAXA
C
      MAXA(1)=1
      DO 200 I=2,JEDN
  200 MAXA(I+1)=MAXA(I)+MHT(I)+1
      NWK=MAXA(JEDN+1)-1
         DO 300 I=2,JEDN
            MAXA(I)=MAXA(I-1)+MAXA(I+1)-MAXA(I)
  300    CONTINUE
C      LS = JEDN+1
C      DO 210 I=1,LS
C  210 MHT(I)=MAXA(I)
C
      NWK=MAXA(JEDN)
      RETURN
      END
C==========================================================================
C======================================================================
      FUNCTION NENPER(NTABF,NNPER)
      DIMENSION NTABF(*)

        NENPER=NTABF(NNPER)
      END
C======================================================================
C==========================================================================
      SUBROUTINE INHEAD (IA,NASLOV)
C
CE Subroutine INHEAD is used for re-writing head of problem
C
      DIMENSION IA(*)
      CHARACTER*250 NASLOV

       DO I=1,80
         NASLOV(I:I)=CHAR(IA(I))
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE HEADIN (IA,NASLOV)
      DIMENSION IA(*)
      CHARACTER*250 NASLOV
C
CE Subroutine HEADIN is used for re-writing head of problem
C

       DO I=1,80
         IA(I)=ICHAR(NASLOV(I:I))
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE REAINT (IA,INT)
C
CE Subroutine REAINT is used converting real to integer data vector
C
      DIMENSION IA(*)
       IA(1)=INT
      END
C==========================================================================
C==========================================================================
      SUBROUTINE INTREA (INT,IA)
C
CE Subroutine INTREA is used converting integer to real data vector
C
      DIMENSION IA(*)
       INT= IA(1)
      END
C==========================================================================
C======================================================================
      SUBROUTINE UACTCF(A,C,B,JDIAG,NEQ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine UACTCF is used for equation solver
CE It is unsymmetric, active column equation solver 
C
      DIMENSION A(*),B(*),JDIAG(*),C(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' UACTCL'
C
C.... UNSYMMETRIC,ACTIVE COLUMN PROFILE EQUATION SOLVER
C
      IF(KKK.EQ.1)THEN
        DO 800 J = 1,NEQ
        JD = JDIAG(J)
  800   C(JD)=1.D0
      ENDIF
C
C.... FACTOR A TO UT*D*U REDUCE B TO Y
C
      JR = 0
      DO 300 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD - JR
      IF(JH.LE.1) GO TO 300
      IS = J + 1 - JH
      IE = J - 1
      IF(KKK.EQ.2) GO TO 250
      K = JR + 1
      ID = 0
C
C.... REDUCE ALL EQUATIONS EXCEPT DIAGONAL
C
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID - IR - 1,I - IS)
      IF(IH.EQ.0) GO TO 150
      A(K) = A(K) - DOT(A(K - IH),C(ID - IH),IH)
      C(K) = C(K) - DOT(C(K - IH),A(ID - IH),IH)
  150 IF(DABS(A(ID)).GT.1.D-20) C(K) = C(K)/A(ID)
C 150 IF(DABS(A(ID)).GT.1.D-40) C(K) = C(K)/A(ID)
  200 K = K + 1
C
C.... REDUCE DIAGONAL TERM
C
      A(JD) = A(JD) - DOT(A(JR + 1),C(JR + 1),JH - 1)
C
C.... FORWARD REDUCE THE R.H.S.
C
  250 IF(KKK.EQ.2) B(J) = B(J) - DOT(C(JR + 1),B(IS),JH - 1)
  300 JR = JD
      IF(KKK.EQ.1) RETURN
C
C.... BACKSUBSTITUTION
C
      J = NEQ
      JD = JDIAG(J)
  500 IF(DABS(A(JD)).GT.1.D-20) B(J) = B(J)/A(JD)
C 500 IF(DABS(A(JD)).GT.1.D-40) B(J) = B(J)/A(JD)
      D = B(J)
      J = J - 1
      IF(J.LE.0) RETURN
      JR = JDIAG(J)
      IF(JD - JR.LE.1) GO TO 700
      IS = J - JD + JR + 2
      K = JR - IS + 1
      DO 600 I = IS,J
      B(I) = B(I) - A(I + K)*D
  600 CONTINUE
  700 JD = JR
      GO TO 500
      END                 
C======================================================================
C==========================================================================
      SUBROUTINE ZFLUID(A,SPSIL,NPT,NETIP,LMAX,NTOTF,IZLAZF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ZFLUID is used for solid-fluid interaction analysis
CE It is writing forces from fluid calculation to file
C
C
      COMMON /INTERA/ IINTER,NPTI
C
      DIMENSION SPSIL(NETIP,*)
      DIMENSION A(*)
      REAL A
C


      CALL MEMORY(LIDENT,LMAX,NPTI*2,1,NTOTF,IZLAZF)
      REWIND IINTER
      CALL READD(A(LIDENT),NPTI,IINTER)
      CALL ZFLU(A(LIDENT),NPTI,SPSIL,NETIP)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZFLU(IDENT,NPTI,SPSIL,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SPSIL(NETIP,*)
      DIMENSION IDENT(2,*),FORCE(3)
C
      IFILE=61

      OPEN(IFILE,FILE='ZFLUID')
      WRITE(IFILE,200) NPTI

      DO I=1,NPTI
        CALL CLEAR(FORCE,3)
        NODES=IDENT(1,I)
        NODEF=IDENT(2,I)
        DO J=1,NETIP
         FORCE(J)=SPSIL(J,NODEF)
        ENDDO 
         WRITE(IFILE,300) NODES,(FORCE(J),J=1,3)
      ENDDO

      CLOSE(IFILE)      

 200  FORMAT (I5)
 300  FORMAT (I5,3D13.5)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZFLU1(IDENT,NPTI,SPSIL,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine ZFLU is used for writing forces from fluid calculation to file
C
      DIMENSION SPSIL(NETIP,*)
      DIMENSION IDENT(2,*),FORCE(3),FORCE0(3)
C
      IFILE=61

      OPEN(IFILE,FILE='ZFLUID')
      WRITE(IFILE,200) NPTI

      DO I=1,NPTI
         READ(IFILE,300) NODES,(FORCE0(J),J=1,3)
        CALL CLEAR(FORCE,3)
        NODES=IDENT(1,I)
        NODEF=IDENT(2,I)
        DO J=1,NETIP
C         FORCE(J)=SPSIL(J,NODEF)
         FORCE(J)=0.5D0*(SPSIL(J,NODEF)+FORCE0(J))
        ENDDO 
         WRITE(IFILE,300) NODES,(FORCE(J),J=1,3)
      ENDDO

      CLOSE(IFILE)      

 200  FORMAT (I5)
 300  FORMAT (I5,3D13.5)

      END
C==========================================================================
      SUBROUTINE IDENSF(IDENT,NPTI,CORDS,CORDF,NPTS,NPTF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IDENT(2,*),CORDS(NPTS,*),CORDF(3,*)
C
CE Subroutine IDENSF is used for solid-fluid interaction analysis
CE It is used for identification same nodes 
CE from meshes of fluid and solid domain
C

      TOL=1.D-10
     
	NPTI=0
      DO I=1,NPTS
       DO J=1,NPTF
         IF( DABS(CORDS(I,1)-CORDF(1,J)).LT.TOL.AND.
     &      DABS(CORDS(I,2)-CORDF(2,J)).LT.TOL.AND.
     &      DABS(CORDS(I,3)-CORDF(3,J)).LT.TOL ) THEN
          NPTI=NPTI+1
C SAMO PRIVREMENO ZBOG SLOBODNE NUMERACIJE CVOROVA 
C KOJA SE KORISTI ZA PRIMERE INTERAKCIJE
          IDENT(1,NPTI)=I

          IDENT(2,NPTI)=J
C          WRITE(3,*) I,J
         ENDIF
       ENDDO
      ENDDO


      END
C==========================================================================
C==========================================================================
      SUBROUTINE WALLPS(A,INIZ,LMAX,GNODE,NPODS,NETIPS,NSTAC,TIME,CCORD,
     &NPTF,NEL,NDIM,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1


      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT

      COMMON /INTERA/ IINTER,NPTI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8


      DIMENSION NPODS(JPS1,*),NEL(NDIM,*)
      DIMENSION INIZ(100)
      DIMENSION A(*),GNODE(2,5,*),CCORD(3,*)
      REAL A
C
      IDVA=2
     
      LIDENT=INIZ(44)
      LMAXAF=INIZ(11)
      NTOTF=INIZ(23)
C      LIDS=IS(3)
      LIDF=INIZ(3)
      LCORDF=INIZ(4)
C      LBRZ=IS(23)
C      LBRZ0=IS(24)
C      LTT1S=IS(20)
C      LTT10S=IS(21)
      LCCORD=INIZ(15)
C      NPTI=IF(45)
      IZLAZF=INIZ(27)
     
      CALL MEMORY(LIDENT,LMAX,NPTI*2,1,NTOTF,IZLAZF)
C===========
      CALL MEMORY(LNZAD,LMAX,NPTI*3*3,1,NTOTF,IZLAZF)
      CALL MEMORY(LZADVR,LMAX,NPTI*3,2,NTOTF,IZLAZF)
      CALL MEMORY(LIDSS,LMAX,NPTF*3,1,NTOTF,IZLAZF)
C===========
      CALL MEMORY(LBRZ,LMAX,JEDN,2,NTOTF,IZLAZF)
      CALL MEMORY(LTT1S,LMAX,JEDN,2,NTOTF,IZLAZF)
      CALL CLEAR (A(LBRZ),JEDN)
      CALL CLEAR (A(LTT1S),JEDN)

      REWIND IINTER
      CALL READD(A(LIDENT),NPTI,IINTER)
C
C
C
        LMAX13=NPODS(JPS1,1)-1
        LIDS=LMAX
        LMAX=LIDS+NP*6
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1
        LCORDS=LMAX
        LMAX=LCORDS+NP*3*IDVA
        CALL READDD(A(LCORDS),NP*3,IPODS,LMAX13,LDUZI)
        CALL IREADD(A(LIDS),NP*6,IPODS,LMAX13,LDUZI)

C
      LMAX13=NPODS(JPS1,55)-1
      IF (LMAX13.GT.0)
     &CALL READDD(A(LBRZ),JEDN,IPODS,LMAX13,LDUZI)

      LMAX13=NPODS(JPS1,52)-1
      IF (NSTAC.EQ.0) LMAX13=NPODS(JPS1,59)-1
      IF (IATYP.EQ.0) LMAX13=NPODS(JPS1,87)-1
       
      IF (LMAX13.GT.0)
     &CALL READDD(A(LTT1S),JEDN,IPODS,LMAX13,LDUZI)
      
C      CALL SOLFLU(A(LIDENT),GNODE,A(LMAXAF),A(LIDS),A(LIDF),A(LBRZ),
C     &A(LTT1S),CCORD,A(LCORDF),NPTI,IZLAZF,NP,NETIP,NSTAC,TIME)

      CALL MOVMSH(A,CCORD,A(LNZAD),A(LZADVR),NEL,A(LIDSS),A(LIDENT),
     &A(LIDS),A(LIDF),A(LBRZ),A(LTT1S),A(LCORDF),NDIM,NET,
     &LMAX,NP,NTOTF,NPTF,NETIPS,NPTI)

								 
      END
C==========================================================================
C==========================================================================
      SUBROUTINE SOLFLU(IDENT,GNODE,MAXAF,IDS,IDF,BRZ,TT1S,
     &CCORD,CORDF,NPTI,IZLAZF,NP,NETIP,NSTAC,TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IDENT(2,*),MAXAF(*),IDS(NP,*),IDF(4,*)
      DIMENSION GNODE(2,5,*),BRZ(*),TT1S(*)
      DIMENSION CCORD(3,*),CORDF(3,*)

      DIMENSION V(3),D(3),NQS(3),NQF(3),DD(3),VOLD(3)

      DO 10 NODEI=1,NPTI
        NODES=IDENT(1,NODEI)
        NODEF=IDENT(2,NODEI)

	  DO I=1,NETIP
C WE PUT NODEI INSTEAD NODES BECAUSE FREE NUMERATION NODES
        NQS(I)=IDS(NODES,I)
C        NQS(I)=IDS(NODEI,I)

        NQF(I)=IDF(I,NODEF)
        V(I)=0.D0
        D(I)=0.D0
        
        IF(NQS(I).NE.0) THEN
          IF (NSTAC.EQ.0) V(I)=BRZ(NQS(I))
          D(I)=TT1S(NQS(I))
        ENDIF
          DD(I)=D(I)-(CCORD(I,NODEF)-CORDF(I,NODEF))
          VOLD(I)=GNODE(2,I,NODEF)
          IF(NSTAC.EQ.1) V(I)=VOLD(I)+DD(I)/TIME 
          IF (NSTAC.EQ.0) GNODE(2,I,NODEF)=V(I)
          CCORD(I,NODEF)=CORDF(I,NODEF)+D(I)
        ENDDO


 10    CONTINUE


C AUTOMATIC REZONING COORDINATES OF NODES:
C       CALL TUBE(CCORD,IZLAZF,5,5,50)
       CALL RCAROT(CCORD,IZLAZF,5,5)

      END
C==========================================================================
      SUBROUTINE IDENTI(NPODS,LSTART,IPAKF,KOLKF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1


      COMMON /INTERA/ IINTER,NPTI
      DIMENSION NPODS(JPS1,*)
      CHARACTER*6    FIPAKI
 

        IOUTS=64
        IOUTF=63
        IINTER=62

      
	  IDVA=2  



        LMAX13=NPODS(JPS1,1)-1
        NPP=NP
        LIDS=LSTART
        LMAX=LIDS+NPP*6
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1
        LCORDS=LMAX
        LMAX=LCORDS+NPP*3*IDVA
C        IF(LMAX.GT.MTOT) CALL ERROR(1)
        CALL READDD(A(LCORDS),NPP*3,IPODS,LMAX13,LDUZI)
        CALL IREADD(A(LIDS),NP*6,IPODS,LMAX13,LDUZI)
   
         NPS=NP
         CALL OUTPAK(IOUTS,A(LCORDS),A(LIDS),NPS)


         
      

         REWIND IPAKF
         LPAKF=LMAX
         LADD=LMAX-LSTART
         CALL READD(A(LPAKF),KOLKF,IPAKF)
         CALL INTREA(LIDATF,A(LPAKF+2))
         CALL OUTFLU(A(LIDATF+LADD),LIDF,LCORDF,NPF)

         CALL OUTPAF(IOUTF,A(LCORDF+LADD),A(LIDF+LADD),NPF)

         
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1

         LIDENT=LMAX						 
         CALL IDENSF(A(LIDENT),NPTI,A(LCORDS),A(LCORDF+LADD),NPS,NPF)
        
            FIPAKI='ZIPAKI'
C            OPEN (IINTER,FILE=FIPAKI,STATUS='UNKNOWN',
C     1                  FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      OPEN (62,FILE='FIPAKI',STATUS='UNKNOWN',
     1FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         KOLKI=(NPTI*2)/IDVA
         IF(KOLKI.GT.0) CALL WRITED(A(LIDENT),KOLKI,IINTER)


      END
C==========================================================================
C==========================================================================