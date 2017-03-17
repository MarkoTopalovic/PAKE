C=======================================================================
C
CE        GENERAL SUBROUTINES
CS        OPSTI PODPROGRAMI
C
C   SUBROUTINE ISPITA
C              WBROJK
C              WIMIDZ
C              READTE
C              READTT
C              VISINE
C              ADRESE
C              ULTABF
C              TABF
C              DELJIV
C              WRITDD
C              IWRITD
C              READDD
C              IREADD
C              CLEAR
C              ICLEAR
C              ZBIRM
C              ZBIRM1
C              JEDNAK
C              JEDNA1
C              IJEDN1
C              READD
C              WRITED
C              WRR
C              IWRR
C              MAXAPR
C              MAXAPRI
C              MINV
C              GMPRD
C              MTRA
C              VEZACE
C              VEZAEL
C              BTAB
C              ORNL3
C              TDOTAB
C              BISECB
C              EC
C              ECDOT
C
C=======================================================================
      SUBROUTINE ISPITA(ACOZ)
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
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ISPITA'
   10 KARTIC = KARTIC + 1
      READ(IULAZ,1000) ACOZ
      IF(ACOZ(1:2).EQ.'C '.OR.ACOZ(1:2).EQ.'C*'.OR.ACOZ(1:2).EQ.'C-'.OR.
     1ACOZ(1:2).EQ.'c ') GO TO 10
      IF(INDFOR.EQ.2) RETURN
      BACKSPACE IULAZ
      RETURN
 1000 FORMAT(A250)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WBROJK(KARTI,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO PRINT NUMBER OF READ CARDS
CS.    P R O G R A M
CS.        ZA STAMPANJE BROJA UCITANIH KARTICA
C .
CE.            KARTI  - FIRST NUMBER CARDS IN GROUP
CS.            KARTI  - POCETNI BROJ KARTICE U GRUPI
C .
CE.            KARTIC - LAST NUMBER CARDS IN GROUP
CS.            KARTIC - POSLEDNJI BROJ KARTICE U GRUPI
C .
CE.            IND    - INDICATOR (=0, CONTINUE;
CE.                                =1, NEW PAGE)
CS.            IND    - INDIKATOR (=0,NASTAVLJA;=1,NOVA STRANA)
C .
C ......................................................................
C
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' WBROJK'
      IF(KARTI.EQ.KARTIC) GO TO 10
C
CE    READ MORE THEN ONE CARDS
CS    UCITANO VISE KARTICA OD JEDNE
C
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) KARTI,KARTIC
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) KARTI,KARTIC
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) KARTI,KARTIC
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) KARTI,KARTIC
      ENDIF
      RETURN
C
CE    READ ONLY ONE CARD
CS    UCITANA JEDNA KARTICA
C
   10 IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) KARTIC
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) KARTIC
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2110) KARTIC
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6110) KARTIC
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//////' KARTICE OD BROJA',I5,' DO',I5///)
 2010 FORMAT(///'1'/' KARTICE OD BROJA',I5,' DO',I5///)
 2100 FORMAT(//////' KARTICA BROJ',I5///)
 2110 FORMAT(///'1'/' KARTICA BROJ',I5///)
C-----------------------------------------------------------------------
 6000 FORMAT(//////' CARD FROM NUMBER ',I5,' TO',I5///)
 6010 FORMAT(///'1'/' CARD FROM NUMBER ',I5,' TO',I5///)
 6100 FORMAT(//////' CARD NUMBER ',I5///)
 6110 FORMAT(///'1'/' CARD NUMBER ',I5///)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VISINE(MHT,ND,LM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE COLUMN HEIGHTS (MHT)
CE.   P R O G R A M
CS.      ZA FORMIRANJE VEKTORA MHT - VISINA STUBOVA U MATRICI SK
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /CDEBUG/ IDEBUG
      DIMENSION MHT(*),LM(*),IA(1)
      EQUIVALENCE (A(1),IA(1))
      MPC(I1,J1)=IA(LMPC+I1+(J1-1)*(NEZAV+1))
      IF(IDEBUG.GT.0) PRINT *, ' VISINE'
C
CE    EQUATION NUMBER
CS    REDNI BROJ JEDNACINE NAJDALJE OD DIJAGONALE
C
      LS=JEDN+1
      DO 120 I=1,ND
      II=LM(I)
      IF(II) 90,120,100
   90 IIP=-II
      DO 95 L=1,NEZAV
        II=MPC(L,IIP)
        IF(II.GT.0)THEN
          IF(II-LS)92,95,95
   92     LS=II
        ENDIF
   95 CONTINUE
      GO TO 120
  100 IF(II-LS)110,120,120
  110 LS=II
  120 CONTINUE
CE  COLUMN HEIGHTS
CS  VISINE STUBOVA(MAXIMUM)
      DO 250 I=1,ND
      II=LM(I)
      IF(II) 160,250,200
  160 IIP=-II
      DO 170 L=1,NEZAV
        II=MPC(L,IIP)
        IF(II.GT.0)THEN
          ME=II-LS
          IF(ME-MHT(II))170,170,165
  165     MHT(II)=ME
        ENDIF
  170 CONTINUE
      GO TO 250
  200 ME=II-LS
      IF(ME-MHT(II))250,250,220
  220 MHT(II)=ME
  250 CONTINUE
      RETURN
      END
C======================================================================
C
C=======================================================================
      SUBROUTINE VEZACE(NEL,NELCV,NE,NCVE)
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
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NELCV(*),NEL(NE,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' VEZACE'
      DO 20 I=1,NE
         DO 10 J=1,NCVE
            N=NEL(I,J)
            IF(N.LE.0)GO TO 10
            NN=N
            IF(N.EQ.0) GO TO 10
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
      STOP ' PROGRAM STOP - PAK01 - VEZACE'
            ENDIF
            NEL(I,J)=K
   10    CONTINUE
   20 CONTINUE
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' GRESKA U UCITAVANJU ULAZNIH PODATAKA O ELEMENTIMA'/
     1' CVOR',I5,' NE POSTOJI'//' PROGRAM STOP - PAK01 - VEZACE'//)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ERROR IN INPUT DATA ABOUT ELEMENTS'/
     1' NODE',I5,' NOT EXIST'//' PROGRAM STOP - PAK01 - VEZACE'//)
C-----------------------------------------------------------------------
      END
C======================================================================
C
C=======================================================================
      SUBROUTINE VEZAEL(MCVEL,MELCV,NE,NMI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS ELEMENTS IN FREE NUMERATION
CS.   P R O G R A M
CS.      ZA POVEZIVANJE BROJEVA ELEMENATA U SLOBODNOJ NUMERACIJI
C .
CE.    V E C T O R S                                                   
CE.               MCVEL(NE) - VECTOR NUMBER OF ELEMENTS       
CE.        MELCV(NMA-NMI+1) - VEKCOR CONECT NUMBER OF ELEMENT
CS.    V E K T O R I                                                   
CS.               MCVEL(NE) - VEKTOR BROJEVA ELEMENATA                
CS.        MELCV(NMA-NMI+1) - VEKTOR VEZE BROJEVA ELEMENATA I NJIHOVOG  
CS.                           REDNOG BROJA U VEKTORU BROJEVA ELEMENATA  
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION MCVEL(*),MELCV(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' VEZAEL'
      DO 10 I=1,NE
         J=MCVEL(I)
         J=J-NMI+1
         MELCV(J)=I
   10 CONTINUE
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE DELJIV(IBROJ,IDELI,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO TEST DIVIDE OF NUMBER WITHOUT RESIDUE
CS.    P R O G R A M
CS.        ZA ISPITIVANJE DELJIVOSTI BROJA BEZ OSTATKA
C .
CE.            IBROJ  - BROJ KOJI DELIMO
CS.            IBROJ  - BROJ KOJI DELIMO
C .
CE.            IDELI  - BROJ KOJIM DELIMO
CS.            IDELI  - BROJ KOJIM DELIMO
C .
CE.            IND    - INDIKATOR OF DIVIDE
CE.                     =0; WITHOUT RESIDUE
CE.                     =0; WITH RESIDUE
CS.            IND    - INDIKATOR DELJIVOSTI
CS.                     =0; BEZ OSTATKA
CS.                     =0; SA OSTATKOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' DELJIV'
      BROJ=IBROJ
      DELI=IDELI
      REZ=BROJ/DELI
      IREZ=REZ
      OST=REZ-IREZ
      TOL=0.0D00
      IND=0
      IF(DABS(OST).GT.TOL) IND=1
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRR(A,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRR'
      WRITE(IZLAZ,5010) CHAR
      WRITE(IZLAZ,5000) (A(I),I=1,N)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(4(1PD18.9))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IWRR(M,N,CHAR)
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
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION M(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IWRR'
      WRITE(IZLAZ,5010) CHAR
      WRITE(IZLAZ,5000) (M(I),I=1,N)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(7I10)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READTE(A,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   P R O G R A M
CS.      TO READING REAL VEKTOR
CS.   P R O G R A M
CS.      ZA UCITAVANJE REALNOG VEKTORA
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READTE'
      BACKSPACE IULAZ
      READ(IULAZ,1000) (A(I),I=1,II)
      WRITE(IZLAZ,1000) (A(I),I=1,II)
 1000 FORMAT(7F10.2)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READTT(ID,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   P R O G R A M
CS.      TO READING INTEGER VEKTOR
CS.   P R O G R A M
CS.      ZA UCITAVANJE CELOBROJNOG VEKTORA
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION ID(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READTT'
      BACKSPACE IULAZ
      READ(IULAZ,1000) (ID(I),I=1,II)
      WRITE(IZLAZ,1000) (ID(I),I=1,II)
 1000 FORMAT(14I5)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TABF(FN,NTAV,IBR,NFN,ARG,Y,NTMX,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO CALCULATE FUNCTIONS IN TABLE FORM
CS.    P R O G R A M
CS        ZA ODREDJIVANJE VREDNOSTI TABELARNO ZADATE FUNKCIJE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION FN(2,NFN,*),NTAV(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' TABF'
      IND=0
      NTMX=NTAV(IBR)
      IF(NTMX.GT.1) GO TO 10
      Y=FN(2,IBR,1)
      RETURN
C
CE    MULTILINEAR CURVE
CS    MULTILINEARNA KRIVA
C
   10 DO 50 I=2,NTMX
         I1=I-1
         IF(ARG.GE.FN(1,IBR,I1).AND.ARG.LE.FN(1,IBR,I)) GO TO 20
         GO TO 50
   20    Y=FN(2,IBR,I1)+(FN(2,IBR,I)-FN(2,IBR,I1))/
     &                  (FN(1,IBR,I)-FN(1,IBR,I1))*(ARG-FN(1,IBR,I1))
         GO TO 100
   50 CONTINUE
C
      IND=1
  100 RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE WRITDD(A,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR AT DIRECT ACCESS DISK
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA NA DISK SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITDD'
C      IF(ISKDSK.EQ.0.AND.N.EQ.NWK)RETURN
      IF(N.EQ.0) RETURN
      NK=N*8/LD
      NN=NK*LD/8
      IF(N.GT.NN) NK=NK+1
C      PRINT *, 'II,N,LD,NK',II,N,LD,NK
      DO 10 K=1,NK
         IK=(K-1)*LD/8+1
         NN=K*LD/8
         IF(K.EQ.NK) NN=N
         LS=LS+1
C         PRINT *,'K,LS,IK,NN',K,LS,IK,NN
         WRITE(II,REC=LS) (A(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'WN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL WRR(A,N,'WN  ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IWRITD(IA,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE INTEGER VECTOR AT DIRECT ACCESS DISK
CS.    P R O G R A M
CS        ZA ZAPIS CELOBROJNOG VEKTORA NA DISK SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IWRITD'
      IF(N.EQ.0) RETURN
      NK=N*4/LD
      NN=NK*LD/4
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/4+1
         NN=K*LD/4
         IF(K.EQ.NK) NN=N
         LS=LS+1
         WRITE(II,REC=LS) (IA(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN 
C         WRITE(3,*) 'IWN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL IWRR(IA,N,'IWN ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READDD(A,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ REAL VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE REALNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READDD'
C      IF(ISKDSK.EQ.0.AND.N.EQ.NWK)RETURN
      IF(N.EQ.0) RETURN
      NK=N*8/LD
      NN=NK*LD/8
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/8+1
         NN=K*LD/8
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (A(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'RN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL WRR(A,N,'RN  ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IREADD(IA,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ INTEGER VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE CELOBROJNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IREADD'
      IF(N.EQ.0) RETURN
      NK=N*4/LD
      NN=NK*LD/4
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/4+1
         NN=K*LD/4
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (IA(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'IRN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL IWRR(IA,N,'IRN ')
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CLEAR(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR -  A)
CE.         A(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA -  A)
CS.         A(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEAR'
      DO 10 I=1,N
   10 A(I)=0.0D0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ICLEAR(IA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR INTEGER ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE CELOBROJNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR - IA)
CE.        IA(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA - IA)
CS.        IA(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ICLEAR'
      DO 10 I=1,N
   10 IA(I)=0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIR4(A,B,C,D,E,B1,C1,D1,E1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 4 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 4 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B1*B(I)+C1*C(I)+D1*D(I)+E1*E(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*),D(*),E(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIR4'
      DO 10 I=1,N
   10 A(I)=B1*B(I)+C1*C(I)+D1*D(I)+E1*E(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIR3(A,B,C,D,B1,C1,D1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 3 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 3 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B1*B(I)+C1*C(I)+D1*D(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*),D(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIR3'
      DO 10 I=1,N
   10 A(I)=B1*B(I)+C1*C(I)+D1*D(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIR2(A,B,C,B1,C1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B1*B(I)+C1*C(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIR2'
      DO 10 I=1,N
   10 A(I)=B1*B(I)+C1*C(I)
      RETURN
      END
C=======================================================================
      SUBROUTINE ZBIR2B(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B(I)+C(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIR2B'
      DO 10 I=1,N
   10 A(I)=B(I)+C(I)
      RETURN
      END
C=======================================================================
      SUBROUTINE ODUZ2B(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO SUBTRACT 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA ODUZIMANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B(I)-C(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ODUZ2B'
      DO 10 I=1,N
   10 A(I)=B(I)-C(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRM3(A,B,C,D,B1,C1,D1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 4 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 4 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=A(I)+B1*B(I)+C1*C(I)+D1*D(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*),D(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRM3'
      DO 10 I=1,N
   10 A(I)=A(I)+B1*B(I)+C1*C(I)+D1*D(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRM2(A,B,C,B1,C1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 3 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 3 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=A(I)+B1*B(I)+C1*C(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRM2'
      DO 10 I=1,N
   10 A(I)=A(I)+B1*B(I)+C1*C(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRM(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=A(I)+C*B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRM'
      DO 10 I=1,N
   10 A(I)=A(I)+C*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRAB(A,B,A1,B1,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=A1*A(I)+B1*B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRAB'
      DO 10 I=1,N
   10 A(I)=A1*A(I)+B1*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRM1(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=A(I)+B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRM1'
      DO 10 I=1,N
   10 A(I)=A(I)+B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRMM(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO SUBTRACTION 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA ODUZIMANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=A(I)-B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRMM'
      DO 10 I=1,N
   10 A(I)=A(I)-B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JEDNAK(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO EQUALIZING 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA IZJEDNACAVANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=C*B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDNAK'
      DO 10 I=1,N
   10 A(I)=C*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JEDNA1(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO EQUALIZING 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA IZJEDNACAVANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDNA1'
      DO 10 I=1,N
   10 A(I)=B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JEDNAS(A,B,N)
C 
C          A(I)=B(I)
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDNAS'
      DO 10 I=1,N
   10 A(I)=B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE IJEDN1(IA,IB,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO EQUALIZING 2 INTEGER VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA IZJEDNAC. 2 CELOBROJNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*),IB(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IJEDN1'
      DO 10 I=1,N
   10 IA(I)=IB(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READD(A,N,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ REAL VECTOR ON SEQUENTIAL ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE REALNOG VEKTORA IZ SEKVENCIJALNOG FILE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READD'
      LDUZZ=30000
      N1=1
      N2=N
      IF(N2.GT.LDUZZ) N2=LDUZZ
   10 READ(II) (A(I),I=N1,N2)
      N1=N1+LDUZZ
      IF(N1.GT.N) GO TO 20
      N2=N2+LDUZZ
      IF(N2.GT.N) N2=N
      GO TO 10
   20 CONTINUE
C      WRITE(3,*) 'RNS,II',N,II
C      CALL WRR(A,N,'RNS ')
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRITED(A,N,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN SEQUENTIAL ACCESS FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U SEKVENCIJALNI FILE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITED'
      LDUZZ=30000
      N1=1
      N2=N
      IF(N2.GT.LDUZZ) N2=LDUZZ
   10 WRITE(II) (A(I),I=N1,N2)
      N1=N1+LDUZZ
      IF(N1.GT.N) GO TO 20
      N2=N2+LDUZZ
      IF(N2.GT.N) N2=N
      GO TO 10
   20 CONTINUE
C      WRITE(3,*) 'WNS,II',N,II
C      CALL WRR(A,N,'WNS ')
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MAXAPR(A,V,R,MAXA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.  P R O G R A M
CE.     TO MULTIPLY MATRIX STORED BY COLUMNS AND VECTOR
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE SLOZENE PO STUPCIMA I VEKTORA
C .
CE.       A - MATRIX STORE BY COLUMNS
CE.       V - VECTOR
CE.       R - R=R+A*V
CE.    MAXA - ADDRESSES OF DIAGONAL ELEMENTS IN MATRIX A
CE.       N - NUMBER OF ROW IN A, ALSO NUMBER OF MEMBERS IN VECTORS V
CE.           AND R
CS.       A - MATRICA UREDJENA PO STUPCIMA
CS.       V - VEKTOR KOJI SE MNOZI MATRICOM
CS.       R - PROIZVOD R=R+A*V
CS.    MAXA - VEKTOR ADRESA DIJAGONALNIH CLANOVA MATRICE A (DIMENZ. N+1)
CS.       N - BROJ VRSTA U MATRICI A, TJ. CLANOVA U VEKTORIMA V I R
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),MAXA(*),V(*),R(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MAXAPR'
      DO 10 IV =1,N
C
CE          RIGHT FROM DIAGONAL
CS          DESNO OD DIJAGONALE
C
            DO 1 I=IV,N
               IMH = MAXA(I+1)-MAXA(I)-1
               IF(IMH.GT.I) GO TO 1
               II=I+MAXA(I)-IV
               IF(II.GE.MAXA(I+1)) GO TO 1
               R(IV)=R(IV)+A(II)*V(I)
    1       CONTINUE
C
CE          LEFT FROM DIAGONAL
CS          LEVO OD DIJAGONALE
C
            IF(IV.EQ.1) GO TO 10
            IMH = MAXA(IV+1) - MAXA(IV) - 1
            IF(IMH.EQ.0) GO TO 10
            IP=IV-IMH
            IK=IV-1
            DO 2 I=IP,IK
               II=MAXA(IV)+IMH-(I-IP)
    2          R(IV)=R(IV)+A(II)*V(I)
   10    CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MAXAPRI(A,V,R,N,irn,jcn,nz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.  P R O G R A M
CE.     TO MULTIPLY MATRIX STORED BY COLUMNS AND VECTOR
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE SLOZENE PO STUPCIMA I VEKTORA
C .
CE.       A - MATRIX STORE BY COLUMNS
CE.       V - VECTOR
CE.       R - R=R+A*V
CE.    MAXA - ADDRESSES OF DIAGONAL ELEMENTS IN MATRIX A
CE.       N - NUMBER OF ROW IN A, ALSO NUMBER OF MEMBERS IN VECTORS V
CE.           AND R
CS.       A - MATRICA UREDJENA PO STUPCIMA
CS.       V - VEKTOR KOJI SE MNOZI MATRICOM
CS.       R - PROIZVOD R=R+A*V
CS.    MAXA - VEKTOR ADRESA DIJAGONALNIH CLANOVA MATRICE A (DIMENZ. N+1)
CS.       N - BROJ VRSTA U MATRICI A, TJ. CLANOVA U VEKTORIMA V I R
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),V(*),R(*)
      integer irn(nz),jcn(nz),nz
C
      IF(IDEBUG.GT.0) PRINT *, ' MAXAPRI'
c      
         do I=1,nz
            if (irn(I).eq.jcn(I)) then
              R(irn(I))=R(irn(I))+A(I)*V(irn(I))
            else
              R(irn(I))=R(irn(I))+A(I)*V(irn(I))
              R(jcn(I))=R(jcn(I))+A(I)*V(jcn(I))
            endif
         enddo
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MAXAPM(A,V,R,MAXA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.  P R O G R A M
CE.     TO MULTIPLY MATRIX STORED BY COLUMNS AND VECTOR
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE SLOZENE PO STUPCIMA I VEKTORA
C .
CE.       A - MATRIX STORE BY COLUMNS
CE.       V - VECTOR
CE.       R - R=R-A*V
CE.    MAXA - ADDRESSES OF DIAGONAL ELEMENTS IN MATRIX A
CE.       N - NUMBER OF ROW IN A, ALSO NUMBER OF MEMBERS IN VECTORS V
CE.           AND R
CS.       A - MATRICA UREDJENA PO STUPCIMA
CS.       V - VEKTOR KOJI SE MNOZI MATRICOM
CS.       R - PROIZVOD R=R-A*V
CS.    MAXA - VEKTOR ADRESA DIJAGONALNIH CLANOVA MATRICE A (DIMENZ. N+1)
CS.       N - BROJ VRSTA U MATRICI A, TJ. CLANOVA U VEKTORIMA V I R
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),MAXA(*),V(*),R(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MAXAPM'
      DO 10 IV =1,N
C
CE       RIGHT FROM DIAGONAL
CS       DESNO OD DIJAGONALE
C
         DO 1 I=IV,N
            IMH = MAXA(I+1)-MAXA(I)-1
            IF(IMH.GT.I) GO TO 1
            II=I+MAXA(I)-IV
            IF(II.GE.MAXA(I+1)) GO TO 1
            R(IV)=R(IV)-A(II)*V(I)
    1    CONTINUE
C
CE       LEFT FROM DIAGONAL
CS       LEVO OD DIJAGONALE
C
         IF(IV.EQ.1) GO TO 10
         IMH = MAXA(IV+1) - MAXA(IV) - 1
         IF(IMH.EQ.0) GO TO 10
         IP=IV-IMH
         IK=IV-1
         DO 2 I=IP,IK
            II=MAXA(IV)+IMH-(I-IP)
    2    R(IV)=R(IV)-A(II)*V(I)
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE MINV3(XJJ,DETT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INVERSE MATRIX 3 X 3
CS.   P R O G R A M
CS.      ZA INVERTOVANJE MATRICE 3 X 3
C .
C ......................................................................
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION XJ1(3,3),XJJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MINV3 '
C
      CALL DETER3(XJJ,DETT)
      IF(DETT.GT.1.D-13) GOTO 10
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) DETT,NLM,NGE
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) DETT,NLM,NGE
      PRINT *, ' DET < 0'
C      IF(DABS(DET).LT.1.D-12) STOP 'MINV3'
      IF(DABS(DETT).LT.1.D-12) RETURN
   10 XJ1(1,1)= (XJJ(2,2)*XJJ(3,3)-XJJ(3,2)*XJJ(2,3))/DETT
      XJ1(2,1)=-(XJJ(2,1)*XJJ(3,3)-XJJ(3,1)*XJJ(2,3))/DETT
      XJ1(3,1)= (XJJ(2,1)*XJJ(3,2)-XJJ(2,2)*XJJ(3,1))/DETT
      XJ1(1,2)=-(XJJ(1,2)*XJJ(3,3)-XJJ(3,2)*XJJ(1,3))/DETT
      XJ1(2,2)= (XJJ(1,1)*XJJ(3,3)-XJJ(3,1)*XJJ(1,3))/DETT
      XJ1(3,2)=-(XJJ(1,1)*XJJ(3,2)-XJJ(3,1)*XJJ(1,2))/DETT
      XJ1(1,3)= (XJJ(1,2)*XJJ(2,3)-XJJ(2,2)*XJJ(1,3))/DETT
      XJ1(2,3)=-(XJJ(1,1)*XJJ(2,3)-XJJ(2,1)*XJJ(1,3))/DETT
      XJ1(3,3)= (XJJ(1,1)*XJJ(2,2)-XJJ(2,1)*XJJ(1,2))/DETT
      CALL JEDNA1(XJJ,XJ1,9)
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(/
     1' DETERMINANTA MANJA ILI JEDNAKA NULI, DET = ',1PD10.3/
     1' ELEMENT =',I8,'     GRUPA =',I5)
C-----------------------------------------------------------------------
 6000 FORMAT(/
     1' ZERO OR NEGATIVE DETERMINANTE, DET = ',1PD10.3/
     1' ELEMENT =',I8,'     GROUP =',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE DETER3(GRAD,DETG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE DETERMINANT OF MATRIX 3 X 3
CS.   P R O G R A M
CS.      ZA RACUNANJE DETERMINANTE MATRICE 3 X 3
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION GRAD(3,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' DETER3'
      DETG=GRAD(1,1)*(GRAD(2,2)*GRAD(3,3)-GRAD(2,3)*GRAD(3,2))+
     1     GRAD(1,2)*(GRAD(2,3)*GRAD(3,1)-GRAD(2,1)*GRAD(3,3))+
     2     GRAD(1,3)*(GRAD(2,1)*GRAD(3,2)-GRAD(2,2)*GRAD(3,1))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MINV(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.  P R O G R A M
CE.      TO INVERT A MATRIX STORED BY COLUMNS IN VECTOR
CS.   P R O G R A M
CS.      ZA INVERTOVANJE MATRICE A SLOZENE PO STUPCIMA U VEKTORU
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),L(*),M(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MINV'
CE    SEARCH FOR LARGEST ELEMENT
CS    NACI NAJVECI ELEMENT
      D=1.0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
CE    INTERCHANGE ROWS
CS    IZMENITI VRSTE
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C     INTERCHANGE COLUMNS
C     IZMENITI KOLONE
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
CE    DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
CE    CONTAINED IN BIGA)
CS    PODELITI KOLONU SA NEGATIVNIM STOZEROM-
CS    VREDNOST STOZERNOG ELEMENTA JE U BIGA.
   45 IF(BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
CE    REDUCE MATRIX
CS    REDUKOVATI MATRICU
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
CE    DIVIDE ROW BY PIVOT
CS    PODELITI VRSTU STOZEROM
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
CE    PRODUCT OF PIVOTS
CS    PROIZVOD STOZERA
      D=D*BIGA
CE    REPLACE PIVOT BY RECIPROCAL
CS    ZAMENITI STOZER RECIPROCNOM VREDNOSCU
      A(KK)=1.0/BIGA
   80 CONTINUE
CE    FINAL ROW AND COLUMN INTERCHANGE
CS    POSLEDNJA IZMENA VRSTA I KOLONA
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE GMPRD(A,B,R,N,M,L)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TWO MATRICES
CS.   P R O G R A M
CS.      ZA MNOZENJE DVE MATRICE
C .
CE.       A= NAME OF FIRST INPUT MATRIX
CE.       B= NAME OF SECOND INPUT MATRIX
CE.       R= NAME OF OUTPUT MATRIX
CE.       N= NUMBER OF ROWS IN A
CE.       M= NUMBER OF COLUMNS IN A AND ROWS IN B
CE.       L= NUMBER OF COLUMNS IN B
CS.       A= IME PRVE ULAZNE MATRICE
CS.       B= IME DRUGE ULAZNE MATRICE
CS.       R= IME IZLAZNE MATRICE
CS.       N= BROJ VRSTA U A
CS.       M= BROJ KOLONA U A I VRSTA U B
CS.       L= BROJ KOLONA U B
C .
CE.       NUMBER OF COLUMNS OF MATRIX A MUST BE EQUAL TO NUMBER OF ROW
CE.       OF MATRIX B
CS.       BROJ KOLONA U A = BROJ VRSTA U B
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),R(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' GMPRD'
      IR=0
      IK=-M
      DO 10 K=1,L
      IK=IK+M
      DO 10 J=1,N
      IR=IR+1
      JI=J-N
      IB=IK
      R(IR)=0
      DO 10 I=1,M
      JI=JI+N
      IB=IB+1
   10 R(IR)=R(IR)+A(JI)*B(IB)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MTRA(A,R,N,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSPOSE A MATRIX
CS.   P R O G R A M
CS.      ZA TRANSPONOVANJE MATRICE A
C .
CE.       A= INPUT MATRIX
CE.       R= OUTPUT MATRIX
CE.       N= NUMBER OF ROWS IN A AND COLUMNS IN R
CE.       M= NUMBER OF COLUMNS IN A AND ROWS IN R
CS.       A= ULAZNA MATRICA
CS.       R= IZLAZNA MATRICA
CS.       N= BROJ VRSTA U A I KOLONA U R
CS.       M= BROJ KOLONA U A I VRSTA U R
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),R(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MTRA'
      IR=0
      DO 30 I=1,N
      IJ=I-N
      DO 30 J=1,M
      IJ=IJ+N
      IR=IR+1
   30 R(IR)=A(IJ)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RSTAZ(NPODS,LRTD,NS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FROM FILE - IPODS
CS.    P R O G R A M
CS.        ZA CITANJE PODATAKA IZ FILE - IPODS
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RSTAZ '
C      WRITE(3,*) 'R,NS,LRTD',NS,LRTD
      LMAX13=NPODS(JPBR,NS)-1
      CALL READDD(A(LRTD),JEDN,IPODS,LMAX13,LDUZI)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WSTAZ(NPODS,LRTD,NS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO WRITE DATA IN FILE - IPODS
CS.    P R O G R A M
CS.        ZA PISANJE PODATAKA NA FILE - IPODS
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WSTAZ '
C      WRITE(3,*) 'W,NS,LRTD',NS,LRTD
      LMAX13=NPODS(JPBR,NS)-1
      CALL WRITDD(A(LRTD),JEDN,IPODS,LMAX13,LDUZI)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RSTAZK(NPODS,LRTD,NS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FROM FILE - IPODS
CS.    P R O G R A M
CS.        ZA CITANJE PODATAKA SA FILE - IPODS
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RSTAZK'
C      WRITE(3,*) 'RK,NS,LRTD',NS,LRTD
      LMAX13=NPODS(JPBR,NS)-1
      IF(JPBR.LT.JPS1.AND.NS.EQ.60) THEN
         CALL READDD(A(LRTD),NWP,IPODS,LMAX13,LDUZI)
         LSKG=LRTD+NWP*IDVA
         NWKP=NWK-NWP
         CALL READDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
      ELSE
         IF(ISKDSK.NE.0) THEN
            IF(NBLOCK.EQ.1) THEN
               CALL READDD(A(LRTD),NWK,IPODS,LMAX13,LDUZI)
            ELSE
               CALL READDB(A(LSK),A(LMAXA),A(LMNQ),A(LLREC),
     1                     NBLOCK,LR,IBLK,LMAX13)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WSTAZK(NPODS,LRTD,NS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO WRITE DATA IN FILE - IPODS
CS.    P R O G R A M
CS.        ZA PISANJE PODATAKA NA FILE - IPODS
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WSTAZK'
C      WRITE(3,*) 'WK,NS,LRTD',NS,LRTD
      LMAX13=NPODS(JPBR,NS)-1
      IF(JPBR.LT.JPS1.AND.NS.EQ.60) THEN
         CALL WRITDD(A(LRTD),NWP,IPODS,LMAX13,LDUZI)
         LSKG=LRTD+NWP*IDVA
         NWKP=NWK-NWP
         CALL WRITDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
      ELSE
         IF(ISKDSK.NE.0) THEN
            IF(NBLOCK.EQ.1) THEN
               CALL WRITDD(A(LRTD),NWK,IPODS,LMAX13,LDUZI)
            ELSE
               CALL WRITEB(A(LSK),A(LMAXA),A(LMNQ),A(LLREC),
     1                     NBLOCK,LR,IBLK,LMAX13)
            ENDIF
         ENDIF 
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE V1V2(EF1,EF2,EF3,I12)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C        FORMIRANJE TRIEDRA JEDINICNIH VEKTORA  V1   V2   VN
C 
      COMMON /CDEBUG/ IDEBUG
      DIMENSION EF1(*),EF2(*),EF3(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' V1V2'
C .                 JEDINICNI VEKTOR  V1  --> EF1
C .                 JEDINICNI VEKTOR  V2  --> EF2
      IF(I12.EQ.1) THEN
C V1=jxVN
         DINT=DSQRT(EF3(1)*EF3(1)+EF3(3)*EF3(3))
         IF(DINT.LT.1.D-6)THEN
            EF1(1)=1.
            EF1(2)=0.
            EF1(3)=0.
         ELSE
            EF1(1)= EF3(3)/DINT
            EF1(2)=0.
            EF1(3)=-EF3(1)/DINT
         ENDIF
      ELSE
C V2=VNxV1
         EF2(1)=EF1(3)*EF3(2)-EF1(2)*EF3(3)
         EF2(2)=EF1(1)*EF3(3)-EF1(3)*EF3(1)
         EF2(3)=EF1(2)*EF3(1)-EF1(1)*EF3(2)
      ENDIF
      RETURN
      END
C=======================================================================
C
C=====================================================================
      SUBROUTINE TRAETP(ETP,ELAST,TSG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS      TRANSFORMACIJA MATRICE ETP(LOKALNO)  U  ELAST(GLOBALNO)
CS      TRANSPONOVANO TSG  *  ETP  *  TSG  =  ELAST
CE      TRANSFORM  MATRIX  ETP(LOCAL)  TO  ELAST(GLOBAL)
CE      TRANSPOSED TSG  *  ETP  *  TSG  =  ELAST
C
      DIMENSION ETP(6,*),ELAST(6,*),TSG(6,*),P(6,6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRAETP'
      DO 80 I=1,6
      DO 80 J=1,6
        P(I,J)=0.D0
          DO 75 K=1,6
   75     P(I,J)=P(I,J)+TSG(K,I)*ETP(K,J)
   80 CONTINUE
      DO 84 I=1,6
      DO 84 J=I,6
        ELAST(I,J)=0.D0
          DO 82 K=1,6
   82     ELAST(I,J)=ELAST(I,J)+P(I,K)*TSG(K,J)
        ELAST(J,I)=ELAST(I,J)
   84 CONTINUE
      RETURN
      END
C=======================================================================
C
C=====================================================================
      SUBROUTINE TRLETP(ETP,ELAST,TSG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS      TRANSFORMACIJA MATRICE ETP(GLOBALNO)  U  ELAST(LOKALNO)
CS       TSG  *  ETP  *  (TSG)T  =  ELAST
CE      TRANSFORM  MATRIX  ETP(GLOBAL)  TO  ELAST(LOCAL)
CE       TSG  *  ETP  *  (TSG)T  =  ELAST
C
      DIMENSION ETP(6,*),ELAST(6,*),TSG(6,*),P(6,6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRLETP'
      DO 80 I=1,6
      DO 80 J=1,6
        P(I,J)=0.D0
          DO 75 K=1,6
   75     P(I,J)=P(I,J)+TSG(I,K)*ETP(K,J)
   80 CONTINUE
      DO 84 I=1,6
      DO 84 J=I,6
        ELAST(I,J)=0.D0
          DO 82 K=1,6
   82     ELAST(I,J)=ELAST(I,J)+P(I,K)*TSG(J,K)
        ELAST(J,I)=ELAST(I,J)
   84 CONTINUE
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE TRANAL (XJJ,TSG,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE BASIS UNIT VECTORS AND STRAIN TRANSFORMATION 
CE.      MATRIX BETWEEN GLOBAL AND LOCAL CARTESIAN SYSTEM
CS.   P R O G R A M
CS.      ZA RACUNANJE JEDINICNIH BAZNIH VEKTORA I MATRICE TRANSFORMACIJE
CS.      DEFORMACIJA IZMEDJU GLOBALNOG I LOKALNOG DEKARTOVOG SISTEMA
CS.      IND = 0  -  FORMIRA JEDINICNE VEKTORE I MATRICU TRANSFORMACIJE
CS.      IND = 1  -  FORMIRA JEDINICNE VEKTORE XJJ
CS.      IND = 2  -  FORMIRA MATRICU TRANSFORMACIJE TSG
C .
C ......................................................................
C
      DIMENSION XJJ(3,*),TSG(6,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRANAL'
C
CS    JEDINICNI BAZNI VEKTORI
CE    BASIS UNIT VECTORS
C
      IF(IND.NE.2) THEN
        XB1=DSQRT(XJJ(1,1)*XJJ(1,1)+XJJ(1,2)*XJJ(1,2)+XJJ(1,3)*XJJ(1,3))
        DO 7 I=1,3
    7   XJJ(1,I)=XJJ(1,I)/XB1
        XJJ(3,1)=XJJ(1,2)*XJJ(2,3)-XJJ(2,2)*XJJ(1,3)
        XJJ(3,2)=XJJ(2,1)*XJJ(1,3)-XJJ(1,1)*XJJ(2,3)
        XJJ(3,3)=XJJ(1,1)*XJJ(2,2)-XJJ(2,1)*XJJ(1,2)
        XB1=DSQRT(XJJ(3,1)*XJJ(3,1)+XJJ(3,2)*XJJ(3,2)+XJJ(3,3)*XJJ(3,3))
        DO 8 I=1,3
    8   XJJ(3,I)=XJJ(3,I)/XB1
        XJJ(2,1)=XJJ(3,2)*XJJ(1,3)-XJJ(1,2)*XJJ(3,3)
        XJJ(2,2)=XJJ(1,1)*XJJ(3,3)-XJJ(3,1)*XJJ(1,3)
        XJJ(2,3)=XJJ(3,1)*XJJ(1,2)-XJJ(1,1)*XJJ(3,2)
C       DO 145 I=1,3
C       DO 145 J=1,3
C          IF(DABS(XJJ(I,J)).LT.1.D-6) XJJ(I,J)=0.D0
C 145   CONTINUE
      ENDIF
C
CS    FORMIRANJE MATRICE TRANSFORMACIJE - TSG
CE    FORM TRANSFORMATION MATRIX - TSG
C
      IF(IND.NE.1) CALL TRANSE(TSG,XJJ)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TRANSE(TSG,XJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE STRAIN TRANSFORMATION MATRIX 
CS.   P R O G R A M
CS.      ZA RACUNANJE MATRICE TRANSFORMACIJE DEFORMACIJA
CS.      XJ - JAKOBIJEVA MATRICA, TSG = TE = (TS)**-T
CS.      XJ - INVERZNA JAKOBIJEVA MATRICA, TSG = (TS)T = (TE)**-1
C .
C ......................................................................
C
      DIMENSION TSG(6,*),XJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRANSE'
C
      TSG(1,1)=XJ(1,1)*XJ(1,1)
      TSG(2,1)=XJ(2,1)*XJ(2,1)
      TSG(3,1)=XJ(3,1)*XJ(3,1)
      TSG(4,1)=2.D0*XJ(1,1)*XJ(2,1)
      TSG(5,1)=2.D0*XJ(2,1)*XJ(3,1)
      TSG(6,1)=2.D0*XJ(3,1)*XJ(1,1)
C
      TSG(1,2)=XJ(1,2)*XJ(1,2)
      TSG(2,2)=XJ(2,2)*XJ(2,2)
      TSG(3,2)=XJ(3,2)*XJ(3,2)
      TSG(4,2)=2.D0*XJ(1,2)*XJ(2,2)
      TSG(5,2)=2.D0*XJ(2,2)*XJ(3,2)
      TSG(6,2)=2.D0*XJ(3,2)*XJ(1,2)
C
      TSG(1,3)=XJ(1,3)*XJ(1,3)
      TSG(2,3)=XJ(2,3)*XJ(2,3)
      TSG(3,3)=XJ(3,3)*XJ(3,3)
      TSG(4,3)=2.D0*XJ(1,3)*XJ(2,3)
      TSG(5,3)=2.D0*XJ(2,3)*XJ(3,3)
      TSG(6,3)=2.D0*XJ(3,3)*XJ(1,3)
C
      TSG(1,4)=XJ(1,1)*XJ(1,2)
      TSG(2,4)=XJ(2,1)*XJ(2,2)
      TSG(3,4)=XJ(3,1)*XJ(3,2)
      TSG(4,4)=XJ(1,1)*XJ(2,2)+XJ(1,2)*XJ(2,1)
      TSG(5,4)=XJ(2,1)*XJ(3,2)+XJ(2,2)*XJ(3,1)
      TSG(6,4)=XJ(3,1)*XJ(1,2)+XJ(3,2)*XJ(1,1)
C
      TSG(1,5)=XJ(1,2)*XJ(1,3)
      TSG(2,5)=XJ(2,2)*XJ(2,3)
      TSG(3,5)=XJ(3,2)*XJ(3,3)
      TSG(4,5)=XJ(1,2)*XJ(2,3)+XJ(1,3)*XJ(2,2)
      TSG(5,5)=XJ(2,2)*XJ(3,3)+XJ(2,3)*XJ(3,2)
      TSG(6,5)=XJ(3,2)*XJ(1,3)+XJ(3,3)*XJ(1,2)
C
      TSG(1,6)=XJ(1,1)*XJ(1,3)
      TSG(2,6)=XJ(2,1)*XJ(2,3)
      TSG(3,6)=XJ(3,1)*XJ(3,3)
      TSG(4,6)=XJ(1,1)*XJ(2,3)+XJ(1,3)*XJ(2,1)
      TSG(5,6)=XJ(2,1)*XJ(3,3)+XJ(2,3)*XJ(3,1)
      TSG(6,6)=XJ(3,1)*XJ(1,3)+XJ(3,3)*XJ(1,1)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TRANSS(TSG,XJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE STRAIN TRANSFORMATION MATRIX 
CS.   P R O G R A M
CS.      ZA RACUNANJE MATRICE TRANSFORMACIJE DEFORMACIJA
CS.      XJ - JAKOBIJEVA MATRICA, TSG = TS = (TE)**-T
CS.      XJ - INVERZNA JAKOBIJEVA MATRICA, TSG = (TE)T = (TS)**-1
C .
C ......................................................................
C
      DIMENSION TSG(6,*),XJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRANSE'
C
      TSG(1,1)=XJ(1,1)*XJ(1,1)
      TSG(2,1)=XJ(2,1)*XJ(2,1)
      TSG(3,1)=XJ(3,1)*XJ(3,1)
      TSG(4,1)=XJ(1,1)*XJ(2,1)
      TSG(5,1)=XJ(2,1)*XJ(3,1)
      TSG(6,1)=XJ(3,1)*XJ(1,1)
C
      TSG(1,2)=XJ(1,2)*XJ(1,2)
      TSG(2,2)=XJ(2,2)*XJ(2,2)
      TSG(3,2)=XJ(3,2)*XJ(3,2)
      TSG(4,2)=XJ(1,2)*XJ(2,2)
      TSG(5,2)=XJ(2,2)*XJ(3,2)
      TSG(6,2)=XJ(3,2)*XJ(1,2)
C
      TSG(1,3)=XJ(1,3)*XJ(1,3)
      TSG(2,3)=XJ(2,3)*XJ(2,3)
      TSG(3,3)=XJ(3,3)*XJ(3,3)
      TSG(4,3)=XJ(1,3)*XJ(2,3)
      TSG(5,3)=XJ(2,3)*XJ(3,3)
      TSG(6,3)=XJ(3,3)*XJ(1,3)
C
      TSG(1,4)=2.D0*XJ(1,1)*XJ(1,2)
      TSG(2,4)=2.D0*XJ(2,1)*XJ(2,2)
      TSG(3,4)=2.D0*XJ(3,1)*XJ(3,2)
      TSG(4,4)=XJ(1,1)*XJ(2,2)+XJ(1,2)*XJ(2,1)
      TSG(5,4)=XJ(2,1)*XJ(3,2)+XJ(2,2)*XJ(3,1)
      TSG(6,4)=XJ(3,1)*XJ(1,2)+XJ(3,2)*XJ(1,1)
C
      TSG(1,5)=2.D0*XJ(1,2)*XJ(1,3)
      TSG(2,5)=2.D0*XJ(2,2)*XJ(2,3)
      TSG(3,5)=2.D0*XJ(3,2)*XJ(3,3)
      TSG(4,5)=XJ(1,2)*XJ(2,3)+XJ(1,3)*XJ(2,2)
      TSG(5,5)=XJ(2,2)*XJ(3,3)+XJ(2,3)*XJ(3,2)
      TSG(6,5)=XJ(3,2)*XJ(1,3)+XJ(3,3)*XJ(1,2)
C
      TSG(1,6)=2.D0*XJ(1,1)*XJ(1,3)
      TSG(2,6)=2.D0*XJ(2,1)*XJ(2,3)
      TSG(3,6)=2.D0*XJ(3,1)*XJ(3,3)
      TSG(4,6)=XJ(1,1)*XJ(2,3)+XJ(1,3)*XJ(2,1)
      TSG(5,6)=XJ(2,1)*XJ(3,3)+XJ(2,3)*XJ(3,1)
      TSG(6,6)=XJ(3,1)*XJ(1,3)+XJ(3,3)*XJ(1,1)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TRA2SE(TSG,XJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE STRAIN TRANSFORMATION MATRIX 
CS.   P R O G R A M
CS.      ZA RACUNANJE MATRICE TRANSFORMACIJE DEFORMACIJA
CS.      XJ - JAKOBIJEVA MATRICA, TSG = TE = (TS)**-T
CS.      XJ - INVERZNA JAKOBIJEVA MATRICA, TSG = (TS)T = (TE)**-1
C .
C ......................................................................
C
      DIMENSION TSG(4,*),XJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRA2SE'
C
      TSG(1,1)=XJ(1,1)*XJ(1,1)
      TSG(2,1)=XJ(2,1)*XJ(2,1)
      TSG(4,1)=XJ(3,1)*XJ(3,1)
      TSG(3,1)=2.D0*XJ(1,1)*XJ(2,1)
C
      TSG(1,2)=XJ(1,2)*XJ(1,2)
      TSG(2,2)=XJ(2,2)*XJ(2,2)
      TSG(4,2)=XJ(3,2)*XJ(3,2)
      TSG(3,2)=2.D0*XJ(1,2)*XJ(2,2)
C
      TSG(1,4)=XJ(1,3)*XJ(1,3)
      TSG(2,4)=XJ(2,3)*XJ(2,3)
      TSG(4,4)=XJ(3,3)*XJ(3,3)
      TSG(3,4)=2.D0*XJ(1,3)*XJ(2,3)
C
      TSG(1,3)=XJ(1,1)*XJ(1,2)
      TSG(2,3)=XJ(2,1)*XJ(2,2)
      TSG(4,3)=XJ(3,1)*XJ(3,2)
      TSG(3,3)=XJ(1,1)*XJ(2,2)+XJ(1,2)*XJ(2,1)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TRA2SS(TSG,XJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE STRAIN TRANSFORMATION MATRIX 
CS.   P R O G R A M
CS.      ZA RACUNANJE MATRICE TRANSFORMACIJE DEFORMACIJA
CS.      XJ - JAKOBIJEVA MATRICA, TSG = TS = (TE)**-T
CS.      XJ - INVERZNA JAKOBIJEVA MATRICA, TSG = (TE)T = (TS)**-1
C .
C ......................................................................
C
      DIMENSION TSG(4,*),XJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRA2SS'
C
      TSG(1,1)=XJ(1,1)*XJ(1,1)
      TSG(2,1)=XJ(2,1)*XJ(2,1)
      TSG(4,1)=XJ(3,1)*XJ(3,1)
      TSG(3,1)=XJ(1,1)*XJ(2,1)
C
      TSG(1,2)=XJ(1,2)*XJ(1,2)
      TSG(2,2)=XJ(2,2)*XJ(2,2)
      TSG(4,2)=XJ(3,2)*XJ(3,2)
      TSG(3,2)=XJ(1,2)*XJ(2,2)
C
      TSG(1,4)=XJ(1,3)*XJ(1,3)
      TSG(2,4)=XJ(2,3)*XJ(2,3)
      TSG(4,4)=XJ(3,3)*XJ(3,3)
      TSG(3,4)=XJ(1,3)*XJ(2,3)
C
      TSG(1,3)=2.D0*XJ(1,1)*XJ(1,2)
      TSG(2,3)=2.D0*XJ(2,1)*XJ(2,2)
      TSG(4,3)=2.D0*XJ(3,1)*XJ(3,2)
      TSG(3,3)=XJ(1,1)*XJ(2,2)+XJ(1,2)*XJ(2,1)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE AXBV( A, B, C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS  VEKTORSKI PROIZVOD
C        C = A X B
C
      DIMENSION A(*),B(*),C(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' AXBV  '
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JEDV(A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDV  '
C
CS  IZRACUNAVANJE JEDINICNOG VEKTORA
C        D^2 = A^2 + B^2 + C^2
C        A=A/D, B=B/D, C=C/D
C
      D=DSQRT(A*A+B*B+C*C)
      A=A/D
      B=B/D
      C=C/D
      RETURN
      END
C=======================================================================
C
C=======================================================================
      FUNCTION INCOT(NGT)
      INCOT=0
      IF(NGT.EQ.0) NGT=2
      IF(NGT.GT.10)THEN
         INCOT=1
         IF(NGT.LT.13)NGT=13
         IF(NGT.GT.13)NGT=15
         NGT=NGT-10
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SIMETR(A,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO EQUALS TRIANGULAR PART OF MATRIX
CS.   P R O G R A M
CS.      ZA IZJEDNACAVANJE DONJEG I GORNJEG TROUGAONOG DELA MATRICE 
C .
CE.       A(II,II) - MATRIX
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SIMETR'
      DO 10 I=1,II
      DO 10 J=I,II
   10 A(J,I)=A(I,J)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZI1(A,B,C,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY MATRIX AND VECTOR 
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE I VEKTORA 
C .
CE.       A(II)    - VECTOR
CE.       B(II,JJ) - MATRIX
CE.       C(JJ)    - VECTOR
CS.       A(II)    - VEKTOR
CS.       B(II,JJ) - MATRICA
CS.       C(JJ)    - VEKTOR
C .       A = A + B * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(II,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZI1'
      DO 10 I=1,II
      DO 10 J=1,JJ
   10 A(I)=A(I)+B(I,J)*C(J)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZ1I(A,B,C,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY MATRIX AND VECTOR 
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE I VEKTORA 
C .
CE.       A(II)    - VECTOR
CE.       B(II,JJ) - MATRIX
CE.       C(JJ)    - VECTOR
CS.       A(II)    - VEKTOR
CS.       B(II,JJ) - MATRICA
CS.       C(JJ)    - VEKTOR
C .       A = A - B * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(II,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZI1'
      DO 10 I=1,II
      DO 10 J=1,JJ
   10 A(I)=A(I)-B(I,J)*C(J)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZI2(A,B,C,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TRANSPONSE MATRIX AND VECTOR 
CS.   P R O G R A M
CS.      ZA MNOZENJE TRANSPONOVANE MATRICE I VEKTORA 
C .
CE.       A(II)    - VECTOR
CE.       B(JJ,II) - MATRIX
CE.       C(JJ)    - VECTOR
CS.       A(II)    - VEKTOR
CS.       B(JJ,II) - MATRICA
CS.       C(JJ)    - VEKTOR
C .       A = A + BT * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZI2'
      DO 10 I=1,II
      DO 10 J=1,JJ
   10 A(I)=A(I)+B(J,I)*C(J)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZ2I(A,B,C,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TRANSPONSE MATRIX AND VECTOR 
CS.   P R O G R A M
CS.      ZA MNOZENJE TRANSPONOVANE MATRICE I VEKTORA 
C .
CE.       A(II)    - VECTOR
CE.       B(JJ,II) - MATRIX
CE.       C(JJ)    - VECTOR
CS.       A(II)    - VEKTOR
CS.       B(JJ,II) - MATRICA
CS.       C(JJ)    - VEKTOR
C .       A = A - BT * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZ2I'
      DO 10 I=1,II
      DO 10 J=1,JJ
   10 A(I)=A(I)-B(J,I)*C(J)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZM1(A,B,C,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY MATRIX - B AND MATRIX - C 
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE - B I MATRICE - C 
C .
CE.       A(II,JJ) - MATRIX
CE.       B(II,KK) - MATRIX
CE.       C(KK,JJ) - MATRIX
CS.       A(II,JJ) - MATRICA
CS.       B(II,KK) - MATRICA
CS.       C(KK,JJ) - MATRICA
C .       A = B * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B(II,*),C(KK,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZM1'
      DO 10 I=1,II
      DO 10 J=1,JJ
         X=0.D0
         DO 20 K=1,KK
   20    X=X+B(I,K)*C(K,J)
   10 A(I,J)=X
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZM2(A,B,C,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TRANSPONSE MATRIX - B AND MATRIX - C 
CS.   P R O G R A M
CS.      ZA MNOZENJE TRANSPONOVANE MATRICE - B I MATRICE - C 
C .
CE.       A(II,JJ) - MATRIX
CE.       B(KK,II) - MATRIX
CE.       C(KK,JJ) - MATRIX
CS.       A(II,JJ) - MATRICA
CS.       B(KK,II) - MATRICA
CS.       C(KK,JJ) - MATRICA
C .       A = BT * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B(KK,*),C(KK,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZM2'
      DO 10 I=1,II
      DO 10 J=1,JJ
         X=0.D0
         DO 20 K=1,KK
   20    X=X+B(K,I)*C(K,J)
   10 A(I,J)=X
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZM3(A,B,C,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TRANSPONSE MATRIX - B AND MATRIX - C 
CS.   P R O G R A M
CS.      ZA MNOZENJE TRANSPONOVANE MATRICE - B I MATRICE - C 
C .
CE.       A(II,JJ) - MATRIX
CE.       B(II,KK) - MATRIX
CE.       C(JJ,KK) - MATRIX
CS.       A(II,JJ) - MATRICA
CS.       B(II,KK) - MATRICA
CS.       C(JJ,KK) - MATRICA
C .       A = B * CT
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B(II,*),C(JJ,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZM3'
      DO 10 I=1,II
      DO 10 J=1,JJ
         X=0.D0
         DO 20 K=1,KK
   20    X=X+B(I,K)*C(J,K)
   10 A(I,J)=X
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZM4(A,B,C,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY TRANSPONSE MATRIX - B AND MATRIX - C 
CS.   P R O G R A M
CS.      ZA MNOZENJE TRANSPONOVANE MATRICE - B I MATRICE - C 
C .
CE.       A(II,JJ) - MATRIX
CE.       B(KK,II) - MATRIX
CE.       C(JJ,KK) - MATRIX
CS.       A(II,JJ) - MATRICA
CS.       B(KK,II) - MATRICA
CS.       C(JJ,KK) - MATRICA
C .       A = BT * CT
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B(KK,*),C(JJ,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZM4'
      DO 10 I=1,II
      DO 10 J=1,JJ
         X=0.D0
         DO 20 K=1,KK
   20    X=X+B(K,I)*C(J,K)
   10 A(I,J)=X
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEM1(A,B,C,W,II,JJ,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(II,JJ) - MATRIX
CE.       B(II,KK) - MATRIX
CE.       C(KK,JJ) - MATRIX
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(II,JJ) - MATRICA
CS.       B(II,KK) - MATRICA
CS.       C(KK,JJ) - MATRICA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (B * C) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B(II,*),C(KK,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEM1'
      DO 10 I=1,II
      DO 10 J=1,JJ
         X=0.D0
         DO 20 K=1,KK
   20    X=X+B(I,K)*C(K,J)
   10 A(I,J)=A(I,J)+X*W
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEV1(A,B,C,W,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE VECTOR - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU VEKTORA - A
C .
CE.       A(II)    - VECTOR
CE.       B(II,JJ) - MATRIX
CE.       C(JJ)    - VECTOR
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(II)    - VEKTOR 
CS.       B(II,JJ) - MATRICA
CS.       C(JJ)    - VEKTOR 
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (B * C) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(II,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEV1'
      DO 10 I=1,II
         X=0.D0
         DO 20 J=1,JJ
   20    X=X+B(I,J)*C(J)
   10 A(I)=A(I)+X*W
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEV2(A,B,C,W,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE VECTOR - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU VEKTORA - A
C .
CE.       A(II)    - VECTOR
CE.       B(JJ,II) - MATRIX
CE.       C(JJ)    - VECTOR
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(II)    - VEKTOR 
CS.       B(JJ,II) - MATRICA
CS.       C(JJ)    - VEKTOR 
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (BT * C) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEV2'
      DO 10 I=1,II
         X=0.D0
         DO 20 J=1,JJ
   20    X=X+B(J,I)*C(J)
   10 A(I)=A(I)+X*W
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEGF(A,B,C,LM,W,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE VECTOR - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU VEKTORA - A
C .
CE.       A(JEDN)  - VECTOR
CE.       B(JJ,II) - MATRIX
CE.       C(JJ)    - VECTOR
CE.       LM(II)   - EQUATION NUMBERS
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(JEDN)  - VEKTOR 
CS.       B(JJ,II) - MATRICA
CS.       C(JJ)    - VEKTOR 
CS.       LM(II)   - BROJEVI JEDNACINA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (BT* C) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),C(*),LM(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGF'
      DO 10 I=1,II
         L=LM(I)
         IF(L.EQ.0) GO TO 10
         X=0.D0
         DO 20 J=1,JJ
   20    X=X+B(J,I)*C(J)
         A(L)=A(L)+X*W
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEGK(A,B,C,LM,W,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(NWE)   - MATRIX
CE.       B(JJ,II) - MATRIX
CE.       C(JJ,JJ) - MATRIX
CE.       LM(II)   - EQUATION NUMBERS
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(NWE)   - MATRICA
CS.       B(JJ,II) - MATRICA
CS.       C(JJ,JJ) - MATRICA
CS.       LM(II)   - BROJEVI JEDNACINA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (BT * C * B) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),C(JJ,*),P(30),LM(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGK'
      IJ=0
      DO 10 I=1,II
c         IF(LM(I).GT.0.OR.LM(I).LT.0) THEN
             DO 20 J=1,JJ
               P(J)=0.D0
            DO 20 K=1,JJ
   20       P(J)=P(J)+B(K,I)*C(K,J)
c         ENDIF
      DO 10 J=I,II
         IJ=IJ+1
c         IF(LM(I).EQ.0.OR.LM(J).EQ.0) GO TO 10
         X=0.D0
         DO 30 K=1,JJ
   30    X=X+P(K)*B(K,J)
         A(IJ)=A(IJ)+X*W
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE PIOKOS(GRAD,TAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO TRANSFORM STRESS (PIOLAKIRCHHOFF - CAUCHY)
CS.   P R O G R A M
CS.      ZA TRANSFORMACIJU NAPONA (PIOLAKIRKHOFOV - KOSIJEV)
C .      T = X * S * XT
C .
C ......................................................................
C
      DIMENSION GRAD(3,*),TAK(*),TAU(3,3),S0T(3,3),POM(3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' PIOKOS'
C
      S0T(1,1)=TAK(1)
      S0T(2,2)=TAK(2)
      S0T(3,3)=TAK(3)
      S0T(1,2)=TAK(4)
      S0T(2,3)=TAK(5)
      S0T(1,3)=TAK(6)
      S0T(2,1)=S0T(1,2)
      S0T(3,2)=S0T(2,3)
      S0T(3,1)=S0T(1,3)
      DO 10 J=1,3
         DO 20 I=1,3
   20    POM(I)=S0T(I,1)*GRAD(J,1)+S0T(I,2)*GRAD(J,2)+S0T(I,3)*GRAD(J,3)
         DO 30 I=J,3
   30    TAU(I,J)=GRAD(I,1)*POM(1)+GRAD(I,2)*POM(2)+GRAD(I,3)*POM(3)
   10 CONTINUE
      TAK(1)=TAU(1,1)
      TAK(2)=TAU(2,2)
      TAK(3)=TAU(3,3)
      TAK(4)=TAU(2,1)
      TAK(5)=TAU(3,2)
      TAK(6)=TAU(3,1)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UKDEFV(STRAIN,TAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO FORM TRIAL ELASTIC STRAIN 
CS.   P R O G R A M
CS.      ZA FORMIRANJE PROBNIH ELASTICNIH DEFORMACIJA
C .      E* = .5 * ( B* - I )
C .
C ......................................................................
C
      DIMENSION STRAIN(*),TAK(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UKDEFV'
C
      STRAIN(1)=.5D0*(TAK(1)-1.D0)
      STRAIN(2)=.5D0*(TAK(2)-1.D0)
      STRAIN(3)=.5D0*(TAK(3)-1.D0)
      STRAIN(4)=TAK(4)
      STRAIN(5)=TAK(5)
      STRAIN(6)=TAK(6)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ROTGRL(GRAD,TAK,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO TRANSFORM STRAIN 
CE.      IND=0; (GREEN LAGRANGE - ROTI GREEN LAGRANGE) 
CE.      IND=1; (ROTI GREEN LAGRANGE - GREEN LAGRANGE) 
CS.   P R O G R A M
CS.      ZA TRANSFORMACIJU DEFORMACIJA
CS.      IND=0; (GRIN LAGRANZEOVA - ROTIRANA GRIN LAGRANZEOVA) 
CS.      IND=1; (ROTIRANA GRIN LAGRANZEOVA - GRIN LAGRANZEOVA) 
C .      IND=0; ER = (X**-T) * E * XT
C .      IND=1; E = (X**-1) * ER * X
C .
C ......................................................................
C
      DIMENSION GRAD(3,*),TAK(*),TAU(3,3),S0T(3,3),POM(3),GRADI(3,3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ROTGRL'
C
C      IST=0
      CALL JEDNA1(GRADI,GRAD,9)
      CALL MINV3(GRADI,DET)
C
      S0T(1,1)=2.D0*TAK(1)+1.D0
      S0T(2,2)=2.D0*TAK(2)+1.D0
      S0T(3,3)=2.D0*TAK(3)+1.D0
      S0T(1,2)=TAK(4)
      S0T(2,3)=TAK(5)
      S0T(1,3)=TAK(6)
      S0T(2,1)=S0T(1,2)
      S0T(3,2)=S0T(2,3)
      S0T(3,1)=S0T(1,3)
C      IF(IST.EQ.1) CALL WRR(S0T,9,'S0T ')
C      IF(IST.EQ.1) CALL WRR(GRAD,9,'GRAD')
C      IF(IST.EQ.1) CALL WRR(GRADI,9,'GRAI')
C
      IF(IND.EQ.0) THEN
      DO 10 J=1,3
         DO 20 I=1,3
   20    POM(I)=S0T(I,1)*GRAD(J,1)+S0T(I,2)*GRAD(J,2)+S0T(I,3)*GRAD(J,3)
C      IF(IST.EQ.1) CALL WRR(POM,3,'POM ')
C         DO 30 I=J,3
         DO 30 I=1,3
   30    TAU(I,J)=GRADI(1,I)*POM(1)+GRADI(2,I)*POM(2)+GRADI(3,I)*POM(3)
   10 CONTINUE
C      IF(IST.EQ.1) CALL WRR(TAU,9,'TAU ')
      ENDIF
C
      IF(IND.EQ.1) THEN
      DO 50 J=1,3
         DO 60 I=1,3
   60    POM(I)=S0T(I,1)*GRAD(1,J)+S0T(I,2)*GRAD(2,J)+S0T(I,3)*GRAD(3,J)
         DO 70 I=J,3
   70    TAU(I,J)=GRADI(I,1)*POM(1)+GRADI(I,2)*POM(2)+GRADI(I,3)*POM(3)
   50 CONTINUE
      ENDIF
C
      TAK(1)=.5D0*(TAU(1,1)-1.D0)
      TAK(2)=.5D0*(TAU(2,2)-1.D0)
      TAK(3)=.5D0*(TAU(3,3)-1.D0)
      TAK(4)=TAU(2,1)
      TAK(5)=TAU(3,2)
      TAK(6)=TAU(3,1)
C      IF(IST.EQ.1) CALL WRR(TAK,9,'TAK ')
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ALMGRL(GRADI,TAK,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO TRANSFORM STRAIN 
CE.      IND=0; (GREEN LAGRANGE - ALMANSI) 
CE.      IND=1; (ALMANSI - GREEN LAGRANGE) 
CS.   P R O G R A M
CS.      ZA TRANSFORMACIJU DEFORMACIJA
CS.      IND=0; (GRIN LAGRANZEOVA - ALMANSIJEVA) 
CS.      IND=1; (ALMANSIJEVA - GRIN LAGRANZEOVA) 
C .      IND=0; EC = (X**-T) * E * (X**-1)
C .      IND=1; E = XT * EC * X
C .
C ......................................................................
C
      DIMENSION GRADI(3,*),TAK(*),TAU(3,3),S0T(3,3),POM(3),GRAD(3,3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ALMGRL'
C
      CALL JEDNA1(GRAD,GRADI,9)
      IF(IND.EQ.0) CALL MINV3(GRAD,DET)
C
      S0T(1,1)=TAK(1)
      S0T(2,2)=TAK(2)
      S0T(3,3)=TAK(3)
      S0T(1,2)=.5D0*TAK(4)
      S0T(2,3)=.5D0*TAK(5)
      S0T(1,3)=.5D0*TAK(6)
      S0T(2,1)=S0T(1,2)
      S0T(3,2)=S0T(2,3)
      S0T(3,1)=S0T(1,3)
C
      DO 10 J=1,3
         DO 20 I=1,3
   20    POM(I)=S0T(I,1)*GRAD(1,J)+S0T(I,2)*GRAD(2,J)+S0T(I,3)*GRAD(3,J)
         DO 30 I=J,3
   30    TAU(I,J)=GRAD(1,I)*POM(1)+GRAD(2,I)*POM(2)+GRAD(3,I)*POM(3)
   10 CONTINUE
C
      TAK(1)=TAU(1,1)
      TAK(2)=TAU(2,2)
      TAK(3)=TAU(3,3)
      TAK(4)=2.D0*TAU(2,1)
      TAK(5)=2.D0*TAU(3,2)
      TAK(6)=2.D0*TAU(3,1)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE GLLOKN(TA,TSG,SIG,ISNA,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO TRANSFORM CAUCHY STRESS (GLOBAL - LOCAL CARTESIAN)
CS.   P R O G R A M
CS.      ZA TRANSFORM. KOSIJEVOG NAPONA (GLOBALNI - LOKALNI DEKARTOV)
C .
C ......................................................................
C
      DIMENSION TA(*),TSG(6,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' GLLOKN'
C
      IF(ISNA.EQ.2) THEN
CS       KOSIJEV NAPON U GLOBALNOM DEKARTOVOM SISTEMU
CE       CAUCHY STRESS IN GLOBAL CARTESIAN SYSTEM
         SIG=TA(I)
      ELSE
CS       KOSIJEV NAPON U LOKALNOM DEKARTOVOM SISTEMU
CE       CAUCHY STRESS IN LOCAL CARTESIAN SYSTEM
         SIG=0.D0
         IF(I.LE.3) THEN
            DO 10 J=1,3
   10       SIG=SIG+TSG(I,J)*TA(J)
            DO 20 J=4,6
   20       SIG=SIG+2.0D0*TSG(I,J)*TA(J)
         ELSE
            DO 30 J=1,3
   30       SIG=SIG+0.5D0*TSG(I,J)*TA(J)
            DO 40 J=4,6
   40       SIG=SIG+TSG(I,J)*TA(J)
         ENDIF     
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE BISEC (X,XM,XP,DX,F,FM,FP,AF,IB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BISEC '
C
C     SOLVE A NONLINEAR EQUATION BY BISECTION
C
C     X  = X-ARGUMENT (POSITIVE)
C     XM = X-MINUS
C     XP = X-PLUS
C     F  = FUNCTION
C     FM = F-MINUS
C     FP = F-PLUS
C     DX = ICREMENT OF X
C     IB = BISECTION FLAG
C          0 = SEARCH FOR XM AND XP
C          1 = XM AND XP ARE OBTAINED (SEARCH FOR SOLUTION BY BISECTION)
C     AF = ACCELERATION FACTOR FOR DX
C
      IF (IB.EQ.1) GO TO 100
C
      DX = AF*DX
      IF (F.LE.0.D0) THEN
          FM = F
          XM = X
          X  = 0.5*(XM+XP)
          IB = 1
      ELSE
          FP = F
          XP = X
          X  = X + DX
      ENDIF
      RETURN
C
C     BISECTION
C
  100 IF (F.LT.0.D0) THEN
        XM1 = X
        X = X - F/(FM-F)*(XM-X)
        IF (X.LE.XP) X = 0.5*(XM1+XP)
        XM = XM1
        FM = F
        DX = XM - X
      ELSEIF (F.GT.0.D0) THEN
        XP1 = X
        X = X + F/(FP-F)*(X-XP)
        IF (X.GE.XM) X = 0.5*(XM+XP1)
        XP = XP1
        FP = F
        DX = X - XP
      ELSE
        DX = 0.D0
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION TENDOT(TENV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TENV(*)
      TENDOT=TENV(1)*TENV(1)+TENV(2)*TENV(2)+TENV(3)*TENV(3)+
     1    2.*(TENV(4)*TENV(4)+TENV(5)*TENV(5)+TENV(6)*TENV(6))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE DTHBTH(TBTH,TDTH,VREM,NLM,IBD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ELEMENT BIRTH AND DEATH OPTION
CS.   P R O G R A M
CS.      ZA UTVRDJIVANJE NASTAJANJA I NESTAJANJA ELEMENATA
C .
C ......................................................................
C
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
C
      DIMENSION TBTH(*),TDTH(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' DTHBTH'
C
      BTH=0.D0
      DTH=0.D0
      IF(INDBTH.EQ.1) BTH=TBTH(NLM)
      IF(INDDTH.EQ.1) DTH=TDTH(NLM)
C
CS    IGNORISE SE NESTAJANJE KADA JE UCITANO  DTH=0
C
      INDTH=INDDTH
      IF(DABS(DTH).LT.1.D-15.AND.INDDTH.EQ.1) INDTH=0
C
CS    NEMA NASTAJANJA I NESTAJANJA ELEMENATA
C
      IF(INDBTH.EQ.0.AND.INDTH.EQ.0) RETURN
C
CS    POMERANJE GRANICA ZA MALI BROJ U LEVO
C
      BT=BTH-1.D-15
      DT=DTH-1.D-15
C
CS    NASTAJANJE ELEMENATA
C
      IF(INDBTH.EQ.1.AND.INDTH.EQ.0.AND.VREM.LT.BT) IBD=1
C
CS    NESTAJANJE ELEMENATA
C
      IF(INDBTH.EQ.0.AND.INDTH.EQ.1.AND.VREM.GT.DT) IBD=1
C
CS    NASTAJANJE PA NESTAJANJE ELEMENATA I OBRNUTO
C
      IF(INDBTH.EQ.1.AND.INDTH.EQ.1) THEN
CS       KADA SU VREMENA ISTA ILI KADA SU OBA VREMENA JEDNAKA NULI
         BD=DABS(BTH)-DABS(DTH)
         IF(DABS(BD).LT.1.D-15) RETURN
C
         IF(BTH.LT.DTH) THEN
            IF(VREM.LT.BT.OR.VREM.GT.DT) IBD=1
         ELSE
            IF(VREM.GT.DT.AND.VREM.LT.BT) IBD=1
         ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
C
C======================================================================
      FUNCTION DELTA(I,J)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CRONECKER DELTA
CS.   P R O G R A M
CS.      ZA KRONEKER DELTA
C .
C ......................................................................
C
      DELTA=0.D0
      IF(I.EQ.J) DELTA=1.D0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE GLAVN3(S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA RACUNANJE EFEKTIVNOG NAPONA ZA 3/D ELEMENT
CE.   PROGRAM
CE.      TO CALCULATE EFECTIVE STRESS FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAVN3'
C
      XMY=S(1)-S(2)
      YMZ=S(2)-S(3)
      ZMX=S(3)-S(1)
      EFEK=0.5D0*(XMY*XMY+YMZ*YMZ+ZMX*ZMX)+
     13.D0*(S(4)*S(4)+S(5)*S(5)+S(6)*S(6))
      S(7)=0.
      IF(EFEK.GT.1.D-19) S(7)=DSQRT(EFEK)
C
C   GLAVNI NAPONI
C
      CALL GLAVN(S)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE GLAVN(S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA RACUNANJE GLAVNIH NAPONA ZA 3/D ELEMENT
CE.   PROGRAM
CE.      TO CALCULATE PRINCIPAL STRESSES FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /PRINCI/ PRINC(3)
      COMMON /DEFNAP/ NAPDEF
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*)
C OVU TOLERANCIJU PROVERITI 1.D-08 -STARA TOLERANCIJA
C      DATA TOL/1.D-08/
      DATA TOL/1.D-16/
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAVN'
      IZILE=0
C
C   INARIJANTE TENZORA
C
C      ONE=0.9999999999999999D0
      ONE=1.D0
      CC=(S(1)+S(2)+S(3))/3.D0
      DO 20 I=1,3                                                              
   20 PRINC(I)=S(I)-CC                                                        
      C2=(PRINC(1)*PRINC(1)+PRINC(2)*PRINC(2)+PRINC(3)*PRINC(3))*.5D0+
     1    S(4)*S(4)+S(5)*S(5)+S(6)*S(6)             
c      C2=-PRINC(1)*PRINC(2)-PRINC(2)*PRINC(3)-PRINC(3)*PRINC(1) 
c     1   +S(4)*S(4)+S(5)*S(5)+S(6)*S(6) 
      C3=PRINC(1)*(PRINC(2)*PRINC(3)-S(5)*S(5))+
     1   S(4)*(S(5)*S(6)-S(4)*PRINC(3))+        
     1   S(6)*(S(4)*S(5)-PRINC(2)*S(6))                                          
      DUM=S(1)*(S(2)*S(3)-S(5)*S(5))+
     1    S(4)*(S(5)*S(6)-S(4)*S(3))+        
     1    S(6)*(S(4)*S(5)-S(2)*S(6))                                          
C
      SQ2 =DSQRT(2.D0)
      PI23=2.D0*DACOS(0.5D0)
C      DUM=DUM*1.D-8
C      IF(DABS(DUM).LT.1.D-9) DUM=1.D-9
      DUM=DUM*TOL
      IF(DABS(DUM).LT.0.1D0*TOL) DUM=0.1D0*TOL
      DUMM=DUM
C      WRITE(3,*) 'DUMM,C2,CC',DUMM,C2,CC
      IF(C2.LT.DABS(DUM))THEN
         DO 330 I=1,3
            PRINC(I)=CC
C ZILE 
C           ZBOG DSQRT(PRINC(I))
            IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  330    CONTINUE
C ZILE
          IF(IZILE.EQ.1) THEN
             PRINC(1)=PRINC(1)+DUMM
             PRINC(2)=PRINC(2)-DUMM
          ENDIF
C ZILE
         RETURN
      ENDIF
      T = DSQRT(C2/1.5D0)                                                           
      A = C3*SQ2/(T*T*T)
C
C     IF ( A .LT. -1.D0 ) A=-ONE
C     IF ( A .GT.  1.D0 ) A= ONE
      IF(DABS(A).GT.ONE) A=DSIGN(ONE,A)
      A=DACOS(A)/3.D0                                                             
      T=T*SQ2
C
      PRINC(1)=T*DCOS(A)                                                       
      PRINC(2)=T*DCOS(A-PI23)                                                
      PRINC(3)=T*DCOS(A+PI23)                                                
C
      DO 240 I=1,3                                                              
         PRINC(I)=PRINC(I)+CC                                                    
C ZILE 
C        ZBOG DSQRT(PRINC(I))
         IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  240 CONTINUE
C  S O R T I R A N J E
      DO 220 I=1,2                                                             
      DO 215 J=I+1,3                                                             
         IF (PRINC(I).GT.PRINC(J)) GO TO 215                                     
         CC=PRINC(I)                                                              
         PRINC(I)=PRINC(J)                                                       
         PRINC(J)=CC                                                              
  215 CONTINUE
  220 CONTINUE                                                                  
C
C ZILE
       IF(IZILE.EQ.1) THEN
          IF(DABS(PRINC(1)-PRINC(2)).LT.DUMM) THEN
             PRINC(1)=PRINC(1)+DUMM
             PRINC(2)=PRINC(2)-DUMM
          ENDIF
          IF(DABS(PRINC(2)-PRINC(3)).LT.DUMM) THEN
             PRINC(2)=PRINC(2)+DUMM
             PRINC(3)=PRINC(3)-DUMM
          ENDIF
       ENDIF
C ZILE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE GLAPR3(S,V)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.      RACUNANJE GLAVNIH PRAVACA NAPONA ZA 3/D ELEMENT
CE.      CALCULATE PRINCIPAL STRESS DIRECTIONS FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /PRINCI/ PRINC(3)
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*),V(3,*),VP(3),A(6)
C OVU TOLERANCIJU PROVERITI 1.D-06 -STARA TOLERANCIJA
C ZA SAVIJANJE
C      DATA TOL/1.D-06/
C ZA ZATEZANJE
C      DATA TOL/1.D-08/
      DATA TOL/1.D-16/
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAPR3'
CE  CHECK FOR ALL SAME ROOTS
      IF(DABS(PRINC(1)-PRINC(2)).LT.TOL.AND.
     1   DABS(PRINC(2)-PRINC(3)).LT.TOL) THEN
        CALL CLEAR(V,9)
        V(1,1)=1.D0
        V(2,2)=1.D0
        V(3,3)=1.D0
        RETURN
      ENDIF
CE  CHECK FOR TWO SAME ROOTS
      IP=2
      IDB=0
      IF(DABS(PRINC(1)-PRINC(2)).LT.TOL)THEN
        IP=3
        IS=1
        IT=2
        IDB=1
      ELSEIF(DABS(PRINC(2)-PRINC(3)).LT.TOL)THEN
        IP=1
        IS=2
        IT=3
        IDB=1
      ELSEIF(DABS(PRINC(1)-PRINC(3)).LT.TOL)THEN
        IP=2
        IS=3
        IT=1
        IDB=1
      ENDIF
C
    1 CONTINUE
C
        DO 10 I=1,3
   10   A(I) = S(I)-PRINC(IP)
        DO 20 I=4,6
   20   A(I) = S(I)
C
      VP(1)=A(4)*A(5)-A(2)*A(6)
      VP(2)=A(4)*A(6)-A(1)*A(5)
      VP(3)=A(1)*A(2)-A(4)*A(4)
      AI2 =DOT(VP,VP,3)
      CALL JEDNA1(V(1,IP),VP,3)
      VP(1)=A(5)*A(6)-A(3)*A(4)
      VP(2)=A(1)*A(3)-A(6)*A(6)
      VP(3)=A(4)*A(6)-A(1)*A(5)
      AI22=DOT(VP,VP,3)
      IF(AI22.GT.AI2) THEN
        AI2=AI22
        CALL JEDNA1(V(1,IP),VP,3)
      ENDIF
      VP(1)=A(2)*A(3)-A(5)*A(5)
      VP(2)=A(5)*A(6)-A(3)*A(4)
      VP(3)=A(4)*A(5)-A(2)*A(6)
      AI22=DOT(VP,VP,3)
      IF(AI22.GT.AI2) THEN
        CALL JEDNA1(V(1,IP),VP,3)
      ENDIF
      CALL JEDV(V(1,IP),V(2,IP),V(3,IP))
C
      IF(IDB.EQ.0) THEN
        IF(IP.EQ.2) THEN
          IP=3
          GO TO 1
        ELSEIF(IP.EQ.3) THEN
          IP=2
          IS=3
          IT=1
        ENDIF
      ELSEIF(IDB.EQ.1) THEN
        VP(1)=1.D0
        VP(2)=0.D0
        VP(3)=VP(2)
        AI2 =DOT(V(1,IP),VP,3)
        IF(DABS(1.D0-DABS(AI2)).LE.TOL) THEN
C        IF(DABS(1.D0-AI2).LE.TOL) THEN
          VP(2)=VP(1)
          VP(1)=VP(3)
        ENDIF
        CALL AXBV(V(1,IP),VP,V(1,IS))
        CALL JEDV(V(1,IS),V(2,IS),V(3,IS))
      ENDIF
C     GLAVNI PRAVCI PO KOLONAMA
      CALL AXBV(V(1,IP),V(1,IS),V(1,IT))
C     TRANSPONOVANJE MATRICE V DA GLAVNI PRAVCI BUDU PO VRSTAMA
      DO 200 I=1,2
      DO 200 J=I+1,3
      DUM=V(I,J)
      V(I,J)=V(J,I)
  200 V(J,I)=DUM
C     DA PROJEKCIJA PRVOG GLAVNOG PRAVCA NA X OSU BUDE POZITIVNA 
C      IF(V(1,1).LT.0.D0) THEN
C         DO 300 I=1,3
C            V(1,I)=-V(1,I)
C            V(3,I)=-V(3,I)
C  300    CONTINUE
C      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CISTIN(TA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA CISCENJE NUMERICKIH GRESAKA ZA NAPONE
CE.   PROGRAM
CE.      TO CLEANING NUMERICAL ERRORS FOR STRESS
C . 
C ......................................................................
C 
      COMMON /FILTER/ TOLNAP,TOLPOM
      COMMON /MATIZO/ E,V
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /CDEBUG/ IDEBUG
      DIMENSION TA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CISTIN'
C
C     ODREDJIVANJE MAKSIMALNE VREDNOSTI MODULA ELASTICNOSTI
C
      TOL=E
      IF(TOL.LT.EX) TOL=EX
      IF(TOL.LT.EY) TOL=EY
      IF(TOL.LT.EZ) TOL=EZ
C      TOL=1.
      TOLNAP=TOL/1.D12
C      RETURN
C
      DO 10 I=1,N
         IF(DABS(TA(I)).LT.TOLNAP) TA(I)=0.D0
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CISTIP(TA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.   PROGRAM
CE.      TO CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C . 
C ......................................................................
C 
      COMMON /FILTER/ TOLNAP,TOLPOM
      COMMON /MAXDUZ/ XL,YL,ZL
      COMMON /CDEBUG/ IDEBUG
      DIMENSION TA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CISTIP'
C
C     ODREDJIVANJE MAKSIMALNE DIMENZIJE MODELA
C
      TOL=XL
      IF(TOL.LT.YL) TOL=YL
      IF(TOL.LT.ZL) TOL=ZL
C      TOL=1.
      TOLPOM=TOL/1.D12
      DO 10 I=1,N
         IF(DABS(TA(I)).LT.TOLPOM) TA(I)=0.D0
   10 CONTINUE
      RETURN
      END
C======================================================================
      FUNCTION TDOTAN(TAUD,E)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TAUD(*),E(*)
      DUM=TAUD(1)-TAUD(2)
      TDOTAN=E(1)*DUM*DUM
      DUM=TAUD(1)-TAUD(3)
      TDOTAN=TDOTAN+E(2)*DUM*DUM
      DUM=TAUD(2)-TAUD(3)
      TDOTAN=TDOTAN+E(3)*DUM*DUM+2.*(E(4)*TAUD(4)*TAUD(4)+
     &       E(5)*TAUD(5)*TAUD(5)+E(6)*TAUD(6)*TAUD(6))
      RETURN
      END      
C======================================================================
      SUBROUTINE XMXMAT(XMX,DEFQP,EM1DVT,CY,AN,AN1,KDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE XMX
CE     BACK STRESS CONSTITUTIVE MATRIX - ANISOTROPIC PLASTICITY AND
C                                                    VISCOPLASTICITY
C
      DIMENSION XMX(*),CY(*),AN(*),AN1(*)
      TOLP = 1.D-7
      TOLEXP = 1.D-6
      DO 10 I=1,KDIM
      AN1A = AN1(I)
      IF (DEFQP.LT.TOLP .AND. AN1A.LT.TOLEXP) THEN
         XMX(I) = EM1DVT*AN(I)*CY(I)*TOLP**TOLEXP
      ELSE
         XMX(I)=EM1DVT*AN(I)*CY(I)*DEFQP**AN1A
      ENDIF
   10 CONTINUE
C************  KONSTANTNO XM
C      XM   =EM1*DVT*AN*CY*DEFQP**AN1
C      DO 10 I=1,KDIM
C   10 XMX(I)=XM
      RETURN
      END
C======================================================================
      SUBROUTINE ENTY(TAUY,E,Y0,CY,AN,EM,DEFQP,DEFQPT,KDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
C  EN KOEFICIJENTI ZA ANIZOTROPNU PLASTICNOST  I  TAUY
C    IND = 0   EN-OVI KONSTANTNI TOKOM ANALIZE
C    IND = 1   EN-OVI KONSTANTNI U KORAKU
C    IND = 2   EN-OVI ODGOVARAJU TEKUCIM PLASTICNIM DEFORMACIJAMA
C
      DIMENSION E(*),Y0(*),CY(*),AN(*),T(6)
      IND=2
      IF (IND.EQ.0) GO TO 19
      IF  (IND.EQ.2) THEN
          EMDEFQ = EM*DEFQP
      ELSE  
          EMDEFQ = EM*DEFQPT
      ENDIF
      DO 10 I = 1,6
      IF  (DABS(AN(I)).LT.1.D-8) THEN
          T(I) = Y0(I) + CY(I)
      ELSE
          T(I) =Y0(I)+CY(I)*(EMDEFQ)**AN(I)
      ENDIF
   10 CONTINUE
      TAUY =FUNAR(T)
      TAUYT=TAUY
      GO TO 22
C
   19 DO 11 I=1,6
   11 T(I)=Y0(I)
      TAUYT=FUNAR(T)
C
   22 SIG2=TAUYT*TAUYT*2./3.
      DO 20 I=1,KDIM
      E(I)=0.5/(T(I)*T(I))
   20 CONTINUE
      F=E(2)+E(3)-E(1)
      G=E(3)+E(1)-E(2)
      H=E(1)+E(2)-E(3)
      E(1)=H
      E(2)=G
      E(3)=F
C... FORM N1 TO N6
      DO 30 I=1,KDIM
   30 E(I)=E(I)*SIG2
      RETURN
      END
C======================================================================
      FUNCTION FUNAR(T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  EFEKTIVNI NAPON ZA ANIZOTROPNU PLASTICNOST (MODEL10)
      DIMENSION T(*)
      FUNAR=DSQRT(0.5*( (T(1)*T(1)+T(2)*T(2)+T(3)*T(3))/3.+
     +                   T(4)*T(4)+T(5)*T(5)+T(6)*T(6) ))
      RETURN
      END
C=========================================================================      
      SUBROUTINE SMOOTH(F,STRESS,NN,NR,NS,NT,RR,SS,TT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   F             ARRAY STORES SAMPLING POINT VALUES
C   NR,NS,NT      INTERPOLATION ORDERS
C   RR,SS,TT      COORDINATES WHERE F IS TO BE EVALUATED
C   FP            VALUE OF F AT (RR,SS,TT)
      COMMON /INTGRA/ INCOTX,INCOTY,INCOTZ
      DIMENSION F(NN,*),STRESS(*)
      DO 50 II=1,NN
      JG=0
      FP=0.D0
      DO 10 I=1,NR
        XLR=1.D0
        IF(NR.GT.1) XLR=POLY(I,NR,RR,INCOTX)
      DO 10 J=1,NS
        XLS=1.D0
        IF(NS.GT.1) XLS=POLY(J,NS,SS,INCOTY)
      DO 10 K=1,NT
        XLT=1.D0
        IF(NT.GT.1) XLT=POLY(K,NT,TT,INCOTZ)
      JG=JG+1
      FP=FP+XLR*XLS*XLT*F(II,JG)
10    CONTINUE
      STRESS(II)=FP
50    CONTINUE
      RETURN
      END
C=========================================================================      
      FUNCTION POLY(I,N,Z,INCOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C   EVALUTE N-TH ORDER LEGENDRE POLYNOMAL OF DEGREE I AT POINT Z
      DIMENSION XG(15),XNC(15),NREF(6)
C
      DATA NREF/0,1,3,6,10,15/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0/
      DATA XNC/            0.D0,              -1.D0,               1.D0,
     1                    -1.D0,               0.D0,               1.D0,
     2                    -1.D0,-.333333333333333D0, .333333333333333D0,
     3                     1.D0,              -1.D0,             -0.5D0,
     4                     0.D0,              0.5D0,               1.D0/
C
      JGI=NREF(N)+I
      IF(INCOT.EQ.0)THEN
        ZI=XG(JGI)
      ELSE
        ZI=XNC(JGI)
      ENDIF
      POLY=1.D0
      DO 10 J=1,N
      JGJ=NREF(N)+J
      IF(J.NE.I)THEN
        IF(INCOT.EQ.0)THEN
          ZJ=XG(JGJ)
        ELSE
          ZJ=XNC(JGJ)
        ENDIF
        POLY=POLY*(Z-ZJ)/(ZI-ZJ)
      ENDIF
10    CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE TRANS4(BLT,TRLN,ND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     MATRICE BM, BB I BS U GAUS TACKI R I S U GLOBALNOM SISTEMU
CS     JEDNACINE (119), (121) I (123)
C
      DIMENSION BLT(6,*)
C ------------------------------
      DIMENSION TRLN(6,6),POM(6)
C ------------------------------
C
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRANS4'
C     DRUGI SLUCAJ
      DO 10 I=1,ND
         POM(1)=TRLN(1,1)*BLT(1,I)+TRLN(1,2)*BLT(2,I)+TRLN(1,3)*BLT(3,I)
     1         +TRLN(1,4)*BLT(4,I)+TRLN(1,5)*BLT(5,I)+TRLN(1,6)*BLT(6,I)
         POM(2)=TRLN(2,1)*BLT(1,I)+TRLN(2,2)*BLT(2,I)+TRLN(2,3)*BLT(3,I)
     1         +TRLN(2,4)*BLT(4,I)+TRLN(2,5)*BLT(5,I)+TRLN(2,6)*BLT(6,I)
         POM(3)=TRLN(3,1)*BLT(1,I)+TRLN(3,2)*BLT(2,I)+TRLN(3,3)*BLT(3,I)
     1         +TRLN(3,4)*BLT(4,I)+TRLN(3,5)*BLT(5,I)+TRLN(3,6)*BLT(6,I)
         POM(4)=TRLN(4,1)*BLT(1,I)+TRLN(4,2)*BLT(2,I)+TRLN(4,3)*BLT(3,I)
     1         +TRLN(4,4)*BLT(4,I)+TRLN(4,5)*BLT(5,I)+TRLN(4,6)*BLT(6,I)
         POM(5)=TRLN(5,1)*BLT(1,I)+TRLN(5,2)*BLT(2,I)+TRLN(5,3)*BLT(3,I)
     1         +TRLN(5,4)*BLT(4,I)+TRLN(5,5)*BLT(5,I)+TRLN(5,6)*BLT(6,I)
         POM(6)=TRLN(6,1)*BLT(1,I)+TRLN(6,2)*BLT(2,I)+TRLN(6,3)*BLT(3,I)
     1         +TRLN(6,4)*BLT(4,I)+TRLN(6,5)*BLT(5,I)+TRLN(6,6)*BLT(6,I)
          BLT(1,I)=POM(1)
          BLT(2,I)=POM(2)
          BLT(3,I)=POM(3)
          BLT(4,I)=POM(4)
          BLT(5,I)=POM(5)
          BLT(6,I)=POM(6)
   10 CONTINUE
      RETURN
      END
C=======================================================================
      INTEGER FUNCTION MODPR1( NMODM )
C
CE    SPACE NECESSARY IN MATERIAL MODEL FOR 1D ELEMENT, 
CE    ON GAUSS INTEGRATION POINT LEVEL
CS    PROSTOR POTREBAN NA NIVOU GAUSOVE TACKE U MATERIJALNOM MODELU
CS    ZA 1D ELEMENT
C
      GO TO (999,999,999,999,  5,999,999,999,999,999,
     1        11, 11,999,999,999,999,999,999,999,999,
     2       999,999,999, 24,999,999,999,999,999, 30,
     3       999,999,999,999,999,999,999,999,999,999,
     4       999,999,999,999,999,999,999,999,999,999,
     5       999,999,999,999,999,999,999,999,999,999,
     6       999,999,999,999,999,999,999,999,999,999,
     7       999,999,999,999,999,999,999,999,999,999,
     8       999,999,999,999,999,999,999,999,999,999,
     9       999,999,999,999,999,999,999,999,999,999),NMODM
    5 MODPR1=4
      RETURN
   11 MODPR1=0
      RETURN
   24 MODPR1=5
      RETURN
   30 MODPR1=21
  999 RETURN
      END
C=======================================================================
      INTEGER FUNCTION MODPR2( NMODM )
C
CE    SPACE NECESSARY IN MATERIAL MODEL FOR 2D ELEMENT, 
CE    ON GAUSS INTEGRATION POINT LEVEL
CS    PROSTOR POTREBAN NA NIVOU GAUSOVE TACKE U MATERIJALNOM MODELU
CS    ZA 2D ELEMENT
C
      GO TO (999,999,999,999,  5,  6,  6,  6,  5,  6,
     1       999,999,999,  6, 15, 16,  6, 18, 19,  6,
     2         6, 22,999,999, 25, 26,999, 28,999, 30,
     3        31,999,999,999,999,999,999,999,999, 40,
     4       999,999,999,999,999,999,999,999,999,999,
     5        51,999,999,999,999,999,999,999,999,999,
     6        61,999,999,999,999,999,999,999,999,999,
     7       999,999,999,999,999,999,999,999,999,999,
     8       999,999,999,999,999,999,999,999,999,999,
     9       999,999,999,999,999,999,999,999,999,999),NMODM
    5 MODPR2=16
      RETURN
    6 MODPR2=21
      RETURN
   15 MODPR2=26
      RETURN
   16 MODPR2=36
      RETURN
   18 MODPR2=30
      RETURN
   19 MODPR2=41
      RETURN
   22 MODPR2=19
      RETURN
   25 MODPR2=16
      RETURN
C VISKOELASTICAN ( 5 JE BROJ TACAKA)
   26 MODPR2=8+4*5
      RETURN
   28 MODPR2=30
      RETURN
   30 MODPR2=26
      RETURN
   31 MODPR2=13
      RETURN
   40 MODPR2=19
      RETURN
C ARGENTINA
   51 MODPR2=32
      RETURN
C     Isotropic Damage Model (Oliver 1996)
   61 MODPR2=11
  999 RETURN
      END
C=======================================================================
      INTEGER FUNCTION MODPRO( NMODM )
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
C
CE    SPACE NECESSARY IN MATERIAL MODEL FOR 3D ELEMENT, SHELL AND 
CE    BEAM SUPERELEMENT, ON GAUSS INTEGRATION POINT LEVEL
CS    PROSTOR POTREBAN NA NIVOU GAUSOVE TACKE U MATERIJALNOM MODELU
CS    ZA 3D I LJUSKU ELEMENT
C
      GO TO (999,999,999,999,  5,  6,  6,  6,  5,  6,
     1       999,999,  6,  6, 15, 16,  6, 18, 19,  6,
     2         6, 22, 15, 24, 25, 26, 27, 28, 29, 30,
     3        31, 32, 33, 34,999,999,999,999,999,999,
     4        41, 42, 43, 44, 45,999,999,999,999,999,
     5       999, 52, 53, 54,999,999,999,999,999,999,
     6        61,999,999,999,999,999,999,999,999,999,
     7       999,999,999,999,999,999,999,999,999,999,
     8       999,999,999,999,999,999,999,999,999,999,
     9       999,999,999,999,999,999,999,999,999,999),NMODM
    5 MODPRO=22
      RETURN
    6 MODPRO=29
      RETURN
   15 MODPRO=36
      RETURN
   16 MODPRO=52
      RETURN
   18 MODPRO=42
      RETURN
   19 MODPRO=57
      RETURN
   22 MODPRO=27
      RETURN
   24 MODPRO=12
      RETURN
   25 MODPRO=6+20*MXS
      RETURN
   26 MODPRO=12+6*5
      RETURN
   27 MODPRO=43
      RETURN
   28 MODPRO=22
      RETURN
   29 MODPRO=22
      RETURN
   30 MODPRO=30
      RETURN
   31 MODPRO=15
      RETURN
   32 MODPRO=20
      RETURN
   33 MODPRO=20
      RETURN
   34 MODPRO=26
      RETURN
C     Drucker-Prager (rakic)
   41 MODPRO=20
      RETURN
C     Mohr-Coulomb (rakic)
   42 MODPRO=20
      RETURN
C     Hoek-Brown (rakic)
   43 MODPRO=20
      RETURN
C     Generalized Hoek-Brown (rakic)
   44 MODPRO=25
      RETURN
C     Maksimovic (rakic)
   45 MODPRO=25
      RETURN
C     Crystal      
   52 MODPRO=155
      RETURN     
C     SMA      
   53 MODPRO=119
      RETURN  
C     SMA VLADA     
   54 MODPRO=57
      RETURN 
C     Isotropic Damage Model (Oliver 1996)
   61 MODPRO=15
  999 RETURN
      END
C=======================================================================
      INTEGER FUNCTION MODTMP( NMODM )
C
CE    INDICATOR FOR MODELS WITH TEMPERATURE DEPENDENCE
CS    INDIKATOR ZA TERMO MODELE
C
      GO TO (999,999,  1,  1,999,999,999,999,999,999,
     1       999,999,999,  1,  1,  1,  1,  1,  1,999,
     2       999,999,999,999,999,999,999,999,999,999,
     3       999,999,999,999,999,999,999,999,999,999,
     4       999,999,999,999,999,999,999,999,999,999,
     5       999,999,  1,999,999,999,999,999,999,999,
     6       999,999,999,999,999,999,999,999,999,999,
     7       999,999,999,999,999,999,999,999,999,999,
     8       999,999,999,999,999,999,999,999,999,999,
     9       999,999,999,999,999,999,999,999,999,999),NMODM
    1 MODTMP=1
      RETURN
  999 MODTMP=0
      RETURN
      END
C=======================================================================
      INTEGER FUNCTION MODORT( NMODM )
C
CE    INDICATOR FOR MATERIAL ORTHOTROPY
CS    INDIKATOR ZA ORTOTROPNE MODELE
C
      GO TO (999,  1,999,  1,999,999,999,999,999,  1,
     1       999,999,  1,999,999,999,  1,  1,  1,999,
     2         1,999,999,999,999,  1,999,999,999,  1,
     3       999,999,999,999,999,999,999,999,999,999,
     4       999,999,999,999,999,999,999,999,999,999,
     5       999,999,999,999,999,999,999,999,999,999,
     6       999,999,999,999,999,999,999,999,999,999,
     7       999,999,999,999,999,999,999,999,999,999,
     8       999,999,999,999,999,999,999,999,999,999,
     9       999,999,999,999,999,999,999,999,999,999),NMODM
    1 MODORT=1
      RETURN
  999 MODORT=0
      RETURN
      END
C=======================================================================
      SUBROUTINE BTAB(FUN,NTA,NFE,MATE,TGT,NL,IND,KFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    LINEAR INTERPOLATION FOR GIVEN TEMPERATURE CURVES
C
      DIMENSION FUN(KFUN,MATE,*),NTA(*)
C
CE    LEFT BOUNDARY OF INTERVAL
C
      IND=0
      NCT=NTA(NFE)
      IF (NCT.EQ.1) THEN
        IND=1
        RETURN
      END IF
      IF (TGT.LT.FUN(1,NFE,1).OR.TGT.GT.FUN(1,NFE,NCT)) THEN
        IND=2
        RETURN
      END IF
      NR=NCT
      NL=1
      NM=(NL+NR)/2
   10 IF ((NR-NL).EQ.1) RETURN
        IF (TGT.GT.FUN(1,NFE,NM)) THEN
          NL=NM
        ELSE IF (TGT.LT.FUN(1,NFE,NM)) THEN
          NR=NM
        ELSE
          NL=NM
          RETURN
        END IF
        NM=(NL+NR)/2
      GO TO 10
      END
C=======================================================================
      SUBROUTINE ORNL3(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    "O.R.N.L." - HARDENING RULE
C
      DIMENSION TAU(*),DEFC(*),OPLUS(*),OMINS(*),EPLUS(6),EMINS(6),
     &          ODIST(6)
C
      DVT   =2.D0/3.
C
      DO 10 I=1,6
        EPLUS(I)=DEFC(I)-OPLUS(I)
        EMINS(I)=DEFC(I)-OMINS(I)
        ODIST(I)=OPLUS(I)-OMINS(I)
   10 CONTINUE
      EPLEF=DSQRT(DVT*TENDOT(EPLUS))
      EMIEF=DSQRT(DVT*TENDOT(EMINS))
      EHET =DSQRT(DVT*TENDOT(ODIST))
C
      IF (IOR1.EQ.0) THEN
        IF ((TDOTAB(TAU,EPLUS).GE.0.).OR.
     &      ((EPLEF.LE.EHET).AND.(TDOTAB(TAU,EMINS).LT.0.).AND.
     &       (EPLEF.GE.EMIEF))) THEN
          DEFQC=EPLEF
        ELSE
          IOR1=-1
          DO 20 I=1,6
            OMINS(I)=DEFC(I)
   20     CONTINUE
          DEFQC=EMIEF
        END IF
      ELSE
        IF ((TDOTAB(TAU,EMINS).GE.0.).OR.
     &      ((EMIEF.LE.EHET).AND.(TDOTAB(TAU,EPLUS).LT.0.).AND.
     &       (EMIEF.GE.EPLEF))) THEN
          DEFQC=EMIEF
        ELSE
          IOR1=0
          DO 30 I=1,6
            OPLUS(I)=DEFC(I)
   30     CONTINUE
          DEFQC=EPLEF
        END IF
      END IF
      RETURN
      END
C=======================================================================
      FUNCTION TDOTAB(TENA,TENB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    THE SCALAR PRODUCT OF TWO TENSORS
C
      DIMENSION TENA(*),TENB(*)
      TDOTAB=TENA(1)*TENB(1)+TENA(2)*TENB(2)+TENA(3)*TENB(3)+
     &       2.*(TENA(4)*TENB(4)+TENA(5)*TENB(5)+TENA(6)*TENB(6))
      RETURN
      END
C=======================================================================
      SUBROUTINE BISECB(X,XL,XD,DX,F,FL,FD,AF,IB,XMIN,TOL0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    SOLVE A NONLINEAR EQUATION BY BISECTION
C
C     X   - CURRENT VALUE OF ARGUMENT IN BISECTION        (INPUT)
C     XL  - LEFT BOUNDARY                                 (INPUT)
C     XD  - RIGHT BOUNDARY                                (OUTPUT)
C     DX  - INCREMENT OF X                                (INPUT)
C     F   - F(X)                                          (INPUT)
C     FL  - F(XL)                                         (INPUT)
C     FD  - F(XD)                                         (OUTPUT)
C     AF  - ACCELERATION FACTOR FOR DX                    (INPUT)
C     IB  - BISECTION FLAG                                (INPUT)
C           0 = SEARCH FOR XL AND XD  (INITIAL)
C           1 = SEARCH FOR SOLUTION BY BISECTION
C          -1 = SOLUTION IS OBTAINED DURING THE BOUNDARY SEARCH
C     XMIN- LEFT BOUNDARY FOR ARGUMENT                    (INPUT)
C     TOL0- MINIMUM DIFERENCE OF BOUNDARIES               (INPUT)
C
      COMMON /CDEBUG/ IDEBUG
C
      IF (IDEBUG.GT.0) PRINT*,' BISECB '
C
      IF (IB.EQ.1) GO TO 100
C
      IF (FL*F.LT.0.D0) THEN
        IF (DABS(XL-X).LE.TOL0) THEN
          X=0.5*(XL+X)
          IB=-1
        ELSE
          IF (DX.LT.0.D0) THEN
            XD=XL
            FD=FL
            XL=X
            FL=F
          ELSE
            XD=X
            FD=F
          END IF
          X=0.5*(XL+XD)
          IB=1
        END IF
      ELSE IF (FL*F.GT.0.D0) THEN
        XL=X
        FL=F
        DX=AF*DX
        X=X+DX
        IF (X.LT.XMIN) X=XMIN
      ELSE
        IB=-1
      END IF
      RETURN
C
CE    BISECTION
C
  100 IF (F*FL.GT.0.D0) THEN
        XL1=X
        IF (DABS(F-FL).GT.1.D-10) THEN
          X=X+F/(FL-F)*(X-XL)
        ELSE
          X=0.5*(XL1+XD)
        END IF
        IF (X.GE.XD.OR.X.LE.XL) X=0.5*(XL1+XD)
        XL=XL1
        FL=F
        DX=X-XL
      ELSE IF (F*FL.LT.0.D0) THEN
        XD1=X
        IF (DABS(F-FD).GT.1.D-10) THEN
          X=X-F/(FD-F)*(XD-X)
        ELSE
          X=0.5*(XL+XD1)
        END IF
        IF (X.LE.XL.OR.X.GE.XD) X=0.5*(XL+XD1)
        XD=XD1
        FD=F
        DX=XD-X
      ELSE
        DX=0.D0
      END IF
      RETURN
      END
C=======================================================================
      FUNCTION EC(TEQ,TIME,THETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    EFFECTIVE CREEP STRAIN
C
      CHARACTER*250 RP1,RP2
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      GO TO (10,20,30,40) ICLAW
C
C     BAILEY-NORTON-OV ZAKON PUZANJA
C
   10 EC   = ACR(1)*(TEQ**ACR(2))*(TIME**ACR(3))
      RETURN
C
C     EKSPONENCIJALNI ZAKON PUZANJA
C
   20 IF (ICELAW.NE.1) THEN
        F    = ACR(1)*DEXP(ACR(2)*TEQ)
        R    = ACR(3)*(TEQ/ACR(4))**ACR(5)
        G    = ACR(6)*DEXP(ACR(7)*TEQ)
C
C     EKSPONENCIJALNI ZAKON PUZANJA (Snyder & Bathe)
C
      ELSE
        F    = ACR(1)*(TEQ**ACR(2))
        R    = ACR(3)*DEXP(ACR(4)*TEQ)
        G    = ACR(5)*((DSINH(ACR(6)*TEQ))**ACR(7))
      END IF
      EC   = F*(1.D0-DEXP(-R*TIME))+G*TIME
      RETURN
C
C     ZAKON PUZANJA SA 8 PARAMETARA
C
   30 EC   = ACR(1)*TEQ**ACR(2)*(TIME**ACR(3)+ACR(4)*TIME**ACR(5)+
     &       ACR(6)*TIME**ACR(7))*DEXP(-ACR(8)/(THETA+273.16D0))
      RETURN
C
C     KORISNICKI ZAKON PUZANJA
C
   40 EC   = UEC(TEQ,THETA,TIME,ACR,RP1,ISG1,ITI1,ITH1)
      RETURN
      END
C=======================================================================
      FUNCTION ECDOT(TEQ,TIME,THETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    EFFECTIVE CREEP STRAIN RATE
C
      CHARACTER*250 RP1,RP2
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      GO TO (10,20,30,40) ICLAW
C
C     BAILEY-NORTON-OV ZAKON PUZANJA
C
   10 ECDOT= ACR(1)*(TEQ**ACR(2))*ACR(3)*(TIME**(ACR(3)-1.))
      RETURN
C
C     EKSPONENCIJALNI ZAKON PUZANJA
C
   20 IF (ICELAW.NE.1) THEN
        F = ACR(1)*DEXP(ACR(2)*TEQ)
        R = ACR(3)*(TEQ/ACR(4))**ACR(5)
        G = ACR(6)*DEXP(ACR(7)*TEQ)
C
C     EKSPONENCIJALNI ZAKON PUZANJA (Snyder & Bathe)
C
      ELSE
        F = ACR(1)*(TEQ**ACR(2))
        R = ACR(3)*DEXP(ACR(4)*TEQ)
        G = ACR(5)*((DSINH(ACR(6)*TEQ))**ACR(7))
      END IF
      ECDOT= R*F*DEXP(-R*TIME)+G
      RETURN
C
C     ZAKON PUZANJA SA 8 PARAMETARA
C
   30 ECDOT= ACR(1)*TEQ**ACR(2)*(ACR(3)*TIME**(ACR(3)-1.D0)+
     &       ACR(4)*ACR(5)*TIME**(ACR(5)-1.D0)+
     &       ACR(6)*ACR(7)*TIME**(ACR(7)-1.D0))*
     &       DEXP(-ACR(8)/(THETA+273.16D0))
      RETURN
C
C     KORISNICKI ZAKON PUZANJA
C
   40 ECDOT= UEC(TEQ,THETA,TIME,BCR,RP2,ISG2,ITI2,ITH2)
      RETURN
      END
C=======================================================================
      FUNCTION TDO2AB(TENA,TENB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    THE SCALAR PRODUCT OF TWO TENSORS
C
      DIMENSION TENA(*),TENB(*)
      TDO2AB=TENA(1)*TENB(1)+TENA(2)*TENB(2)+TENA(4)*TENB(4)+
     &       2.*TENA(3)*TENB(3)
      RETURN
      END
C=======================================================================
      FUNCTION TENDO2(TENV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TENV(*)
      TENDO2=TENV(1)*TENV(1)+TENV(2)*TENV(2)+TENV(4)*TENV(4)+
     1    2.*TENV(3)*TENV(3)
      RETURN
      END
C======================================================================
      SUBROUTINE ORNL2(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    "O.R.N.L." - HARDENING RULE
C
      DIMENSION TAU(*),DEFC(*),OPLUS(*),OMINS(*),EPLUS(4),EMINS(4),
     &          ODIST(4)
C
      DVT   =2.D0/3.
C
      DO 10 I=1,4
        EPLUS(I)=DEFC(I)-OPLUS(I)
        EMINS(I)=DEFC(I)-OMINS(I)
        ODIST(I)=OPLUS(I)-OMINS(I)
   10 CONTINUE
      EPLEF=DSQRT(DVT*TENDO2(EPLUS))
      EMIEF=DSQRT(DVT*TENDO2(EMINS))
      EHET =DSQRT(DVT*TENDO2(ODIST))
C
      IF (IOR1.EQ.0) THEN
        IF ((TDO2AB(TAU,EPLUS).GE.0.).OR.
     &      ((EPLEF.LE.EHET).AND.(TDO2AB(TAU,EMINS).LT.0.).AND.
     &       (EPLEF.GE.EMIEF))) THEN
          DEFQC=EPLEF
        ELSE
          IOR1=-1
          DO 20 I=1,4
            OMINS(I)=DEFC(I)
   20     CONTINUE
          DEFQC=EMIEF
        END IF
      ELSE
        IF ((TDO2AB(TAU,EMINS).GE.0.).OR.
     &      ((EMIEF.LE.EHET).AND.(TDO2AB(TAU,EPLUS).LT.0.).AND.
     &       (EMIEF.GE.EPLEF))) THEN
          DEFQC=EMIEF
        ELSE
          IOR1=0
          DO 30 I=1,4
            OPLUS(I)=DEFC(I)
   30     CONTINUE
          DEFQC=EPLEF
        END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE DEV3(SHAT,SE,GB,CXX,CYY,DXX,DYY,CED,FP,D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SHAT(*),SE(*),GB(*)
        D       = 1. + FP*(CXX+CYY + FP*CED)
        SHAT(1) = (SE(1) + FP*DXX)/D
        SHAT(2) = (SE(2) + FP*DYY)/D
        SHAT(3) = - SHAT(1) - SHAT(2)
        SHAT(4) = SE(4)/GB(1)
        SHAT(5) = SE(5)/GB(2)
        SHAT(6) = SE(6)/GB(3)
      RETURN
      END
C======================================================================
      SUBROUTINE INIOR(SE,E,CE,CPE,GB,CXX,CYY,CXY,CYX,DXX,DYY,CED,CHM,
     &                  CHV,FP,KDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SE(*),E(*),CE(KDIM,*),CPE(3,*),GB(*),CHM(3,*),CHV(*)
C
      DVA = 2.D0
      XX1 = E(1)+DVA*E(2)
      XX2 = E(1)-E(3)
      XX3 = DVA*E(2)+E(3)
      CXX =  XX1*(CPE(1,1)+CHM(1,1))-XX2*(CPE(1,2)+CHM(1,2))
     &      -XX3*(CPE(1,3)+CHM(1,3))
      CYX = -XX1*(CPE(2,1)+CHM(2,1))+XX2*(CPE(2,2)+CHM(2,2))
     &      +XX3*(CPE(2,3)+CHM(2,3))
      XX1 = E(1)-E(2)
      XX2 = E(1)+DVA*E(3)
      XX3 = E(2)+DVA*E(3)
      CYY = -XX1*(CPE(2,1)+CHM(2,1))+XX2*(CPE(2,2)+CHM(2,2))
     &      -XX3*(CPE(2,3)+CHM(2,3))
      CXY =  XX1*(CPE(1,1)+CHM(1,1))-XX2*(CPE(1,2)+CHM(1,2))
     &      +XX3*(CPE(1,3)+CHM(1,3))
      CED = CXX*CYY - CXY*CYX
      DXX = CYY*SE(1)+CXY*SE(2)
      DYY = CYX*SE(1)+CXX*SE(2)
      IF (KDIM.EQ.6) THEN
        DO 10 I=1,KDIM-3
          I3=I+3
   10     GB(I)=1+FP*(CE(I3,I3)*DVA+CHV(I))*E(I3)
      ELSE
        GB(1)=1+FP*(CE(3,3)*DVA+CHV(1))*E(4)
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE CHATP(CHPM,CHPV,EPL,DEFQP,EM1,DVT,CY,AN,AN1,KDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     KONSTITUTIVNA MATRICA ZA "BACK STRESS"
CE     BACK STRESS CONSTITUTIVE MATRIX
C
      DIMENSION CHPM(3,*),CHPV(*),EPL(*),XMX(6)
      DIMENSION CY(*),AN(*),AN1(*)
      DO 10 I=1,KDIM
   10 XMX(I)=EM1*DVT*AN(I)*CY(I)*DEFQP**AN1(I)
      CHPM(1,1)=XMX(1)*(EPL(1)+EPL(2))
      CHPM(1,2)=-XMX(2)*EPL(1)
      CHPM(1,3)=-XMX(3)*EPL(2)
      CHPM(2,1)=-XMX(1)*EPL(1)
      CHPM(2,2)=XMX(2)*(EPL(1)+EPL(3))
      CHPM(2,3)=-XMX(3)*EPL(3)
      CHPM(3,1)=-XMX(1)*EPL(2)
      CHPM(3,2)=-XMX(2)*EPL(3)
      CHPM(3,3)=XMX(3)*(EPL(2)+EPL(3))
      DO 20 I=1,KDIM-3
   20   CHPV(I)=1.7320508*XMX(I+3)*EPL(I+3)
      RETURN
      END
C======================================================================
      SUBROUTINE CHATC(CHCM,CHCV,ECR,EMC1,DVT,CCF,KDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    KONSTIT. MAT. ZA VEKTOR TRANSLACIJE POTENCIJALNE POVRSI PUZANJA
CE    CONSTIT. MAT. FOR TRANSLATION VECTOR OF POTENTIAL CREEP SURFACE
C
      DIMENSION CHCM(3,*),CHCV(*),ECR(*)
      FAK=DVT*EMC1*CCF
      CHCM(1,1)=FAK*(ECR(1)+ECR(2))
      CHCM(1,2)=-FAK*ECR(1)
      CHCM(1,3)=-FAK*ECR(2)
      CHCM(2,1)=CHCM(1,2)
      CHCM(2,2)=FAK*(ECR(1)+ECR(3))
      CHCM(2,3)=-FAK*ECR(3)
      CHCM(3,1)=CHCM(1,3)
      CHCM(3,2)=CHCM(2,3)
      CHCM(3,3)=FAK*(ECR(2)+ECR(3))
      DO 10 I=1,KDIM-3
   10   CHCV(I)=FAK*ECR(I+3)
      RETURN
      END
C======================================================================
      SUBROUTINE INDEF3(DDEF,E,SHAT,FP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    RACUNANJE PRIRASTAJA NEELASTICNE DEFORMACIJE
CE    INELASTIC DEFORMATION INCREMENT
C
      DIMENSION DDEF(*),E(*),SHAT(*)
      DDEF(1) = FP*( (E(1)+E(2))*SHAT(1)-E(1)*SHAT(2)-E(2)*SHAT(3))
      DDEF(2) = FP*(-E(1)*SHAT(1)+(E(1)+E(3))*SHAT(2)-E(3)*SHAT(3))
      DDEF(3) =-DDEF(1)-DDEF(2)
      DDEF(4) = FP*E(4)*SHAT(4)
      DDEF(5) = FP*E(5)*SHAT(5)
      DDEF(6) = FP*E(6)*SHAT(6)
      RETURN
      END
C======================================================================
      SUBROUTINE ABACK3(ALFA,ALFAT,DDEF,CHM,CHV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    VEKTOR ZAOSTALOG NAPONA / TRANSLACIJE POVRSI PUZANJA
CE    BACK STRESS VECTOR / VECTOR OF CREEP SURFACE TRANSLATION
C
      DIMENSION ALFA(*),ALFAT(*),DDEF(*),CHM(3,*),CHV(*)
      DO 10 L=1,3
        L3 = L+3
        DUM = 0.D0
        DO 20 M=1,3
   20     DUM = DUM+CHM(L,M)*DDEF(M)
        ALFA(L)  = ALFAT(L)+DUM
        DUM      = CHV(L)*DDEF(L3)
   10   ALFA(L3) = ALFAT(L3)+DUM
      RETURN
      END
C======================================================================
      FUNCTION TDOTA2(TAUD,E)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TAUD(*),E(*)
      DUM=TAUD(1)-TAUD(2)
      TDOTA2=E(1)*DUM*DUM
      DUM=TAUD(1)-TAUD(4)
      TDOTA2=TDOTA2+E(2)*DUM*DUM
      DUM=TAUD(2)-TAUD(4)
      TDOTA2=TDOTA2+E(3)*DUM*DUM+2.*E(4)*TAUD(3)*TAUD(3)
      RETURN
      END      
C======================================================================
      SUBROUTINE DEV2(SHAT,SE,GB,CXX,CYY,DXX,DYY,CED,FP,D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SHAT(*),SE(*),GB(*)
        D       = 1. + FP*(CXX+CYY + FP*CED)
        SHAT(1) = (SE(1) + FP*DXX)/D
        SHAT(2) = (SE(2) + FP*DYY)/D
        SHAT(4) = - SHAT(1) - SHAT(2)
        SHAT(3) = SE(3)/GB(1)
      RETURN
      END
C======================================================================
      SUBROUTINE INDEF2(DDEF,E,SHAT,FP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    RACUNANJE PRIRASTAJA NEELASTICNE DEFORMACIJE
CE    INELASTIC DEFORMATION INCREMENT
C
      DIMENSION DDEF(*),E(*),SHAT(*)
      DDEF(1) = FP*( (E(1)+E(2))*SHAT(1)-E(1)*SHAT(2)-E(2)*SHAT(4))
      DDEF(2) = FP*(-E(1)*SHAT(1)+(E(1)+E(3))*SHAT(2)-E(3)*SHAT(4))
      DDEF(4) =-DDEF(1)-DDEF(2)
      DDEF(3) = FP*E(4)*SHAT(3)
      RETURN
      END
C======================================================================
      SUBROUTINE ABACK2(ALFA,ALFAT,DDEF,CHM,CHV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    VEKTOR ZAOSTALOG NAPONA / TRANSLACIJE POVRSI PUZANJA
CE    BACK STRESS VECTOR / VECTOR OF CREEP SURFACE TRANSLATION
C
      DIMENSION ALFA(*),ALFAT(*),DDEF(*),CHM(3,*),CHV(*)
      ALFA(1) = CHM(1,1)*DDEF(1)+CHM(1,2)*DDEF(2)+CHM(1,3)*DDEF(4)
     &         +ALFAT(1)
      ALFA(2) = CHM(2,1)*DDEF(1)+CHM(2,2)*DDEF(2)+CHM(2,3)*DDEF(4)
     &         +ALFAT(2)
      ALFA(3) = CHV(1)*DDEF(3)+ALFAT(3)
      ALFA(4) =-ALFA(1)-ALFA(2)
      RETURN
      END
C=======================================================================
      SUBROUTINE UCITAM(MODEL,MOD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO FIND POSITION DATA FOR MATERIALS
CS.    P R O G R A M
CS.        ZA DESIFROVANJE POLOZAJA PODATAKA O MATERIJALIMA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUPLAP/ IDVA
      COMMON /REPERM/ MREPER(4)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION MODEL(4,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UCITAM'
      IF(NETIP.EQ.10.OR.NETIP.EQ.11.OR.NETIP.EQ.12.OR.NETIP.EQ.14)RETURN
      CALL ICLEAR(MREPER,4)
      DO 130 J=1,NMATM
      MO=MODEL(1,J)
      IF(MO.EQ.MOD) GO TO 140
  130 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MOD
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MOD
      STOP
C
  140 GO TO (    1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     1          11, 11, 13, 14, 15, 14, 17, 18, 17,  6,
     2          21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     3          31, 32, 33, 34,999,999,999,999,999, 40,
     4          41, 42, 43, 44, 45,999,999,999,999,999,
     5          51, 52, 53, 54,999,999,999,999,999,999,
     6          61,999,999,999,999,999,999,999,999,999,
     7         999,999,999,999,999,999,999,999,999,999,
     8         999,999,999,999,999,999,999,999,999,999,
     9         999,999,999,999,999,999,999,999,999,999),MOD
      GO TO 999
C
    1 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
    2 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
    3 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*6*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
    4 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*24*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
    5 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+2)*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(4)=MATE
      RETURN
C
    6 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+6)*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(4)=MATE
      RETURN
C
    7 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
    8 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
    9 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   10 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   11 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*2*MAXT*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(4)=MATE
      RETURN
C
   13 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   14 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*16*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
   15 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*6*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
   17 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*72*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
   18 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*26*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN
C
   21 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   22 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   23 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
C     VISKOZNI OTPOR STAPA
   24 LREP = MODEL(4,J)
      MREPER(1) = LREP
      MREPER(4) = MODEL(2,J)
      RETURN
C
   25 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   26 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*2*(MAXT+1)*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(3)=MAXT
      MREPER(4)=MATE
      RETURN
C
   27 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   28 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C
   29 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C
   30 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C 
   31 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C
   32 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C
   33 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C
   34 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN
C
   40 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(6*MAXT+6)*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(4)=MATE
      RETURN
C     Drucker-Prager (rakic)
   41 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     Mohr-Coulomb (rakic)
   42 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     Hoek-Brown (rakic)
   43 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     Generalized Hoek-Brown (rakic)
   44 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     Maksimovic (rakic)
   45 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     ARGENTINA
   51 LREP=MODEL(4,J)
      MREPER(1)=LREP
      MREPER(4)=MODEL(2,J)
      RETURN
C     Crystal     
   52 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN    
C     SMA     
   53 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(3)=LTEM
      MREPER(4)=MATE
      RETURN         
C     SMA VLADA    
   54 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MREPER(1)=LREP
      MREPER(4)=MATE
      RETURN 
C     Isotropic Damage Model (Oliver 1996)
   61 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+2)*IDVA
      MREPER(1)=LREP
      MREPER(2)=LNTA
      MREPER(4)=MATE
c      MREPER(1)=LREP
c      MREPER(4)=MODEL(2,J)
      RETURN
C
  999 RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(///' P O D A C I   O   M O D E L U   M A T E R I J A L A
     1 B R O J  =',I5,'  N I S U   U C I T A N I')
C-----------------------------------------------------------------------
 6000 FORMAT(///' D A T A    F O R    M O D E L    O F   M A T E R I A L
     1   N U M B E R  =',I5,'  I S   N O T   R E A D')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEGV(A,B,LM,W,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(NWE)   - MATRIX
CE.       B(JJ,II) - MATRIX
CE.       LM(II)   - EQUATION NUMBER
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(NWE)   - MATRICA
CS.       B(JJ,II) - MATRICA
CS.       LM(II)   - BROJEVI JEDNACINA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A + (BT * B) * W
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(JJ,*),LM(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGV'
      IJ=0
      DO 10 I=1,II
      DO 10 J=I,II
         IJ=IJ+1
         IF(LM(I).EQ.0.OR.LM(J).EQ.0) GO TO 10
         X=0.D0
         DO 30 K=1,JJ
   30    X=X+B(K,I)*B(K,J)
         A(IJ)=A(IJ)+X*W
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE ZBIRNN(A,B,C,N)
C....  VEKTORSKI IZRAZ   A() <--- A() * C - B()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(*)
      DO 10 I=1,N
   10 A(I)=A(I)*C-B(I)
      RETURN
      END
C=======================================================================
C
C======================================================================
      DOUBLE PRECISION FUNCTION RNORM2(B,N)
C....  KVADRAT NORME VEKTORA
      DOUBLE PRECISION B
      DIMENSION B(*)
      RNORM2=0.D0
      DO 10 I=1,N
   10 RNORM2=RNORM2+B(I)*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ENERG(RTDT,FTDT,JEDN,ENE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     RACUNANJE ENERGIJE
C
      DIMENSION RTDT(*),FTDT(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ENERG'
      ENE=0.0D0
      DO 10 I=1,JEDN
      ENE=ENE+RTDT(I)*FTDT(I)
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C======================================================================
      DOUBLE PRECISION FUNCTION DOT(A,B,N)
C....  SKALARNI PROIZVOD VEKTORA
      DOUBLE PRECISION A,B
      DIMENSION A(*),B(*)
      DOT=0.D0
      DO 10 I=1,N
   10 DOT=DOT+A(I)*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      DOUBLE PRECISION FUNCTION RNORM(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....   EUKLIDSKA NORMA VEKTORA
      DIMENSION A(*)
      RNORM=0.D0
      DO 10 I=1,N
   10 RNORM=RNORM+A(I)*A(I)
      RNORM=DSQRT(RNORM)
      RETURN
      END
C=======================================================================
      SUBROUTINE ENEDOF(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(*)
C....  ENERGIJE PO STEPENIMA SLOBODE
      DO 10 I=1,N
10    B(I)=A(I)*B(I)
      RETURN
      END 
C=======================================================================
C
C======================================================================
      SUBROUTINE DIJAD(F,Q,P,P1,P2,P3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO DIADIC
CS.   P R O G R A M
CS.      ZA DIJADSKO MNOZENJE
C .
C ......................................................................
C
      DIMENSION F(3,*),Q(3,*),P(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' DIJAD '
C
      F(1,1)=P1*Q(1,1)*P(1,1)+P2*Q(2,1)*P(2,1)+P3*Q(3,1)*P(3,1)
      F(1,2)=P1*Q(1,1)*P(1,2)+P2*Q(2,1)*P(2,2)+P3*Q(3,1)*P(3,2)
      F(1,3)=P1*Q(1,1)*P(1,3)+P2*Q(2,1)*P(2,3)+P3*Q(3,1)*P(3,3)
      F(2,1)=P1*Q(1,2)*P(1,1)+P2*Q(2,2)*P(2,1)+P3*Q(3,2)*P(3,1)
      F(2,2)=P1*Q(1,2)*P(1,2)+P2*Q(2,2)*P(2,2)+P3*Q(3,2)*P(3,2)
      F(2,3)=P1*Q(1,2)*P(1,3)+P2*Q(2,2)*P(2,3)+P3*Q(3,2)*P(3,3)
      F(3,1)=P1*Q(1,3)*P(1,1)+P2*Q(2,3)*P(2,1)+P3*Q(3,3)*P(3,1)
      F(3,2)=P1*Q(1,3)*P(1,2)+P2*Q(2,3)*P(2,2)+P3*Q(3,3)*P(3,2)
      F(3,3)=P1*Q(1,3)*P(1,3)+P2*Q(2,3)*P(2,3)+P3*Q(3,3)*P(3,3)
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE DIJADS(Q,P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO DIADIC
CS.   P R O G R A M
CS.      ZA DIJADSKO MNOZENJE SIMETRICNIH TENZORA
C .
C ......................................................................
C
      DIMENSION F(6),Q(3,*),P(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' DIJADS '
C
      F(1)=P(1)*Q(1,1)*Q(1,1)+P(2)*Q(2,1)*Q(2,1)+P(3)*Q(3,1)*Q(3,1)
      F(2)=P(1)*Q(1,2)*Q(1,2)+P(2)*Q(2,2)*Q(2,2)+P(3)*Q(3,2)*Q(3,2)
      F(3)=P(1)*Q(1,3)*Q(1,3)+P(2)*Q(2,3)*Q(2,3)+P(3)*Q(3,3)*Q(3,3)
      F(4)=P(1)*Q(1,1)*Q(1,2)+P(2)*Q(2,1)*Q(2,2)+P(3)*Q(3,1)*Q(3,2)
      F(5)=P(1)*Q(1,2)*Q(1,3)+P(2)*Q(2,2)*Q(2,3)+P(3)*Q(3,2)*Q(3,3)
      F(6)=P(1)*Q(1,1)*Q(1,3)+P(2)*Q(2,1)*Q(2,3)+P(3)*Q(3,1)*Q(3,3)
      CALL JEDNA1(P,F,6)
      RETURN
      END
C=======================================================================
C
C======================================================================
      FUNCTION UEC(SIGMA,THETA,TIME,ACR,RP,ISG,ITI,ITH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER*250 RP
      CHARACTER*8 OPER
      CHARACTER*1 B
      DIMENSION STACK(20),ACR(20)
C
            DATA OPER/'-+*/^#!?'/
C
      IF (ISG.NE.21) ACR(ISG)=SIGMA
      IF (ITH.NE.21) ACR(ITH)=THETA
      IF (ITI.NE.21) ACR(ITI)=TIME
C
      REZ=0.D0
      ICR=0
      I=0
      ISTACK=0
      KRP=KRAJ(RP)
   20 I=I+1
      IF (I.GT.KRP) GO TO 999
      B=RP(I:I)
      IF (B.EQ.'A') THEN
        ICR=ICR+1
        IF (ICR.GT.9) THEN
          I=I+2
        ELSE
          I=I+1
        END IF
        IF (ICR.NE.1) THEN
          ISTACK=ISTACK+1
          STACK(ISTACK)=REZ
        END IF
        REZ=ACR(ICR)
        GO TO 20
      ELSE
        DO 30 J=1,8
   30     IF (B.EQ.OPER(J:J)) GOTO (110,120,130,140,150,160,170,20) J
        STOP 'ERROR IN THE RPN-STRING'
  110   REZ=STACK(ISTACK)-REZ
        STACK(ISTACK)=0.D0
        ISTACK=ISTACK-1
        GO TO 20
  120   REZ=STACK(ISTACK)+REZ
        STACK(ISTACK)=0.D0
        ISTACK=ISTACK-1
        GO TO 20
  130   REZ=STACK(ISTACK)*REZ
        STACK(ISTACK)=0.D0
        ISTACK=ISTACK-1
        GO TO 20
  140   REZ=STACK(ISTACK)/REZ
        STACK(ISTACK)=0.D0
        ISTACK=ISTACK-1
        GO TO 20
  150   REZ=STACK(ISTACK)**REZ
        STACK(ISTACK)=0.D0
        ISTACK=ISTACK-1
        GO TO 20
  160   REZ=DEXP(REZ)
        GO TO 20
  170   REZ=-REZ
        GO TO 20
      END IF
  999 UEC=REZ
      RETURN
      END
C=======================================================================
      FUNCTION ECTAU(TEQ,DEFQCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    EFFECTIVE CREEP STRAIN ( ONLY FOR " POWER CREEP " )
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPC/  ACR(20),BCR(20)
C
      ECTAU=((ACR(1)**(1.D0/ACR(3)))*(TEQ**(ACR(2)/ACR(3)))*DT+
     &      (DEFQCT**(1.D0/ACR(3))))**ACR(3)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZBIRKM(A,B,C,N,MAXA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE MATRICE I VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(MAXA(I))=A(MAXA(I))+C*B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),MAXA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZBIRKM'
      DO 10 I=1,N
         M=MAXA(I)
   10 A(M)=A(M)+C*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MNOZMU(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO ADD 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA SABIRANJE MATRICE I VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=A(I)+B(I)*C(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZMU'
      DO 10 I=1,N
   10 A(I)=A(I)+B(I)*C(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STERM3(ETHERM,TGT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R  O G R A M
CE.      TO EVALUATE THERMAL STRAIN
CS.   P R O G R A M
CS.      ZA IZRACUNAVANJE TERMICKE DEFORMACIJE
C .
C ......................................................................
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION ETHERM(3)
C
      DTGT0=TGT-TEMP0
      DO 10 I=1,3
   10   ETHERM(I)=ALFA(I)*DTGT0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEGKN(A,B1,C,B,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(II,JJ)   - MATRIX
CE.       B1(JJ,II) - MATRIX
CE.       C(JJ,JJ) - MATRIX
CE.       B(JJ,II) - MATRIX
CE.       LM(II)   - EQUATION NUMBERS
CE.       W        - WEIGHTS COEFFICIENTS
CS.       A(II,JJ)   - MATRICA
CS.       B1(JJ,II) - MATRICA
CS.       C(JJ,JJ) - MATRICA
CS.       B(JJ,II) - MATRICA
CS.       LM(II)   - BROJEVI JEDNACINA
CS.       W        - TEZINSKI KOEFICIJENTI
C .       A = A - B1T * C * B
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B1(JJ,*),C(JJ,*),B(JJ,*),P(30)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGKN'
      DO 10 I=1,II
            DO 20 J=1,JJ
               P(J)=0.D0
            DO 20 K=1,JJ
   20       P(J)=P(J)+B1(K,I)*C(K,J)
      DO 10 J=1,II
         X=0.D0
         DO 30 K=1,JJ
   30    X=X+P(K)*B(K,J)
         A(I,J)=A(I,J)-X
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTEGKP(A,B1,C,B,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INTEGRATE MATRIX - A 
CS.   P R O G R A M
CS.      ZA INTEGRACIJU MATRICE - A
C .
CE.       A(II,II)   - MATRIX
CE.       B1(JJ,II) - MATRIX
CE.       C(JJ,JJ) - MATRIX
CE.       B(JJ,II) - MATRIX
CS.       A(II,II)   - MATRICA
CS.       B1(JJ,II) - MATRICA
CS.       C(JJ,JJ) - MATRICA
CS.       B(JJ,II) - MATRICA
C .       A = A + B1T * C * B
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(II,*),B1(JJ,*),C(JJ,*),B(JJ,*),P(30)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTEGKP'
      DO 10 I=1,II
            DO 20 J=1,JJ
               P(J)=0.D0
            DO 20 K=1,JJ
   20       P(J)=P(J)+B1(K,I)*C(K,J)
      DO 10 J=1,II
         X=0.D0
         DO 30 K=1,JJ
   30    X=X+P(K)*B(K,J)
         A(I,J)=A(I,J)+X
   10 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE MNOZT1(A,B,C,II,JJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO MULTIPLY VECTOR AND MATRIX 
CS.   P R O G R A M
CS.      ZA MNOZENJE VEKTORA I MATRICE
C .
CE.       A(II)    - VECTOR
CE.       B(JJ)    - VECTOR
CE.       C(JJ,II) - MATRIX
CS.       A(II)    - VEKTOR
CS.       B(JJ)    - VEKTOR
CS.       C(JJ,II) - MATRICA
C .       A = A + B * C
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*),C(JJ,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MNOZI1'
      DO 10 I=1,II
      DO 10 J=1,JJ
   10 A(I)=A(I)+B(J)*C(J,I)
      RETURN
      END
CR======================================================================
CR    Subroutine ODUZBA dodata pri pravljenju oduzimanje pocetne 
CR    deformacije za geomehaniku!
CR======================================================================
      SUBROUTINE ODUZBA(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO SUBTRACT 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA ODUZIMANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .          A(I)=B(I)-A(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ODUZBA'
      DO 10 I=1,N
   10 A(I)=B(I)-A(I)
      RETURN
      END
C=======================================================================
C====Podprogram za bisekcije u materijalnim modelima za geomehaniku=====
C=======================================================================
      SUBROUTINE BISECTG (X,XM,XP,DX,F,FM,FP,AF,IB,JP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C     
      IF(IDEBUG.GT.0) PRINT *, ' BISEC '
C
C     SOLVE A NONLINEAR EQUATION BY BISECTION
C
C     X  = X-ARGUMENT (POSITIVE)
C     XM = X-MINUS
C     XP = X-PLUS
C     F  = FUNCTION
C     FM = F-MINUS
C     FP = F-PLUS
C     DX = ICREMENT OF X
C     IB = BISECTION FLAG
C          0 = SEARCH FOR XM AND XP
C          1 = XM AND XP ARE OBTAINED (SEARCH FOR SOLUTION BY BISECTION)
C     AF = ACCELERATION FACTOR FOR DX
C
      IF (IB.EQ.1) GO TO 100
C
      IF (JP.EQ.1) THEN    ! ZNAK PRVOG IZVODA POZITIVAN
C
          DX = AF*DX
          IF (F.LT.0.D0) THEN
              FM = F
              XM = X
              X  = 0.5*(XM+XP)
              IB = 1
          ELSE
              FP = F
              XP = X
              X  = X-DX
          ENDIF
          RETURN
      ENDIF
C
      IF (JP.EQ.2) THEN    ! ZNAK PRVOG IZVODA NEGATIVAN
C
          DX = AF*DX
          IF (F.GT.0.D0) THEN
              FP = F
              XP = X
              X  = 0.5*(XM+XP)
              IB = 1
          ELSE
              FM = F
              XM = X
              X  = X-DX
          ENDIF
          RETURN
      ENDIF

C
C     BISECTION
C     
  100 IF (JP.EQ.1)THEN    ! ZNAK PRVOG IZVODA POZITIVAN
          IF (F.GT.0.D0) THEN
              XP1 = X
              X   = X-(XP-X)/(FP-F)*F
              IF (X.GE.XP) X=0.5*(XM+XP1)
              XP  = XP1
              FP  = F
              DX  = X-XP
          ELSEIF (F.LT.0.D0) THEN
              XM1 = X
              X   = X-(XM-X)/(FM-F)*F
              IF (X.LE.XM) X=0.5*(XM1+XP)
              XM  = XM1
              FM  = F
              DX  = X-XM
          ENDIF
          RETURN
      ENDIF
C
      IF (JP.EQ.2)THEN    ! ZNAK PRVOG IZVODA NEGATIVAN
          IF (F.LT.0.D0) THEN 
              XM1 = X
              X   = X-(XM-X)/(FM-F)*F
              IF (X.GE.XM) X=0.5*(XM1+XP)
              XM  = XM1
              FM  = F
              DX  = X-XM
          ELSEIF (F.GT.0.D0) THEN
              XP1 = X
              X   = X-(XP-X)/(FP-F)*F
              IF (X.LE.XP) X=0.5*(XM+XP1)
              XP  = XP1
              FP  = F
              DX  = X-XP
          ENDIF
          RETURN
      ENDIF
      END
C=======================================================================
C=======================================================================
C====Podprogram za bisekcije u materijalnim modelima za geomehaniku=====
C=======================================================================
      SUBROUTINE BISECTR (X,XM,XP,DX,F,FM,FP,AF,IB,JP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C     
      IF(IDEBUG.GT.0) PRINT *, ' BISEC '
C
C     SOLVE A NONLINEAR EQUATION BY BISECTION
C
C     X  = X-ARGUMENT (POSITIVE)
C     XM = X-MINUS
C     XP = X-PLUS
C     F  = FUNCTION
C     FM = F-MINUS
C     FP = F-PLUS
C     DX = ICREMENT OF X
C     IB = BISECTION FLAG
C          0 = SEARCH FOR XM AND XP
C          1 = XM AND XP ARE OBTAINED (SEARCH FOR SOLUTION BY BISECTION)
C     AF = ACCELERATION FACTOR FOR DX
C
      IF (IB.EQ.1) GO TO 100
C
      IF (JP.EQ.1) THEN    ! ZNAK PRVOG IZVODA POZITIVAN
C
          DX = AF*DX
          IF (F.GT.0.D0) THEN
              FP = F
              XP = X
              X  = 0.5*(XM+XP)
              IB = 1
          ELSE
              FM = F
              XM = X
              X  = X+DX
          ENDIF
          RETURN
      ENDIF
C
      IF (JP.EQ.2) THEN    ! ZNAK PRVOG IZVODA NEGATIVAN
C
          DX = AF*DX
          IF (F.LT.0.D0) THEN
              FM = F
              XM = X
              X  = 0.5*(XM+XP)
              IB = 1
          ELSE
              FP = F
              XP = X
              X  = X+DX
          ENDIF
          RETURN
      ENDIF
C
C     BISECTION
C     
  100 IF (JP.EQ.1)THEN    ! ZNAK PRVOG IZVODA POZITIVAN
          IF (F.GT.0.D0) THEN
              XP = X
              FP = F
              X  = XM-(XM-XP)/(FM-FP)*FM
              DX = X-XM
          ELSEIF (F.LT.0.D0) THEN
              XM = X
              FM = F
              X  = XM-(XM-XP)/(FM-FP)*FM
              DX = X-XM
          ENDIF
          RETURN
      ENDIF
C
      IF (JP.EQ.2)THEN    ! ZNAK PRVOG IZVODA NEGATIVAN
          IF (F.LT.0.D0) THEN 
              XM = X
              FM = F
              X  = XP-(XP-XM)/(FP-FM)*FP
              DX  = X-XP
          ELSEIF (F.GT.0.D0) THEN
              XP = X
              FP = F
              X  = XP-(XP-XM)/(FP-FM)*FP
              DX = X-XP
          ENDIF
          RETURN
      ENDIF
      END
C==========================================================================
