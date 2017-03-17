C=======================================================================
C
CE       READ DATA FOR MATERIALS
CS       UCITAVANJE PODATAKA O MATERIJALIMA
C
C   SUBROUTINE UCMATE
C              UCMODE
C              MATREP
C              MATERI
C              UMOD01
C              UMOD02
C              UMOD03
C              UMOD04
C              UMOD05
C              UMOD07
C              UMOD09
C              UMOD10
C              UMOD11
C              UMOD13
C              UMOD41
C              UMOD42
C              UMOD43
C              UMOD44
C              UMOD45
C              UMOD52
C              UMOD53
C              UMOD54
C              NAPPLD
C
C=======================================================================
      SUBROUTINE UCMATE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO FORM POINTERS AND CALL ROUTINES FOR READING DATA FOR ALL
CE.        MATERIAL MODELS
CS.    P R O G R A M
CS.        ZA FORMIRANJE REPERA I POZIVANJE PROGRAMA
CS.        ZA UCITAVANJE PODATAKA O MATERIJALIMA
C .
CE.    P O I N T E R S
CE.        LMODEL - DATA FOR MATERIAL MODELS
CE.        LGUSM  - DATA FOR MATERIAL MASS DENSITY
CS.    R E P E R I
CS.        LMODEL - PODACI O MODELIMA MATERIJALA
C .
CE.    V A R I A B L E S
CE.        NMATM  - TOTAL NUMBER OF DIFFERENT MATERIAL MODELS /3/
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /MATIZO/ E,V
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UCMATE'
      E=0.D0
      EX=0.D0
      EY=0.D0
      EZ=0.D0
      LMODEL=LMAX
      LMAX=LMODEL+NMATM*4
      IF(LMAX.GT.MTOT) CALL ERROR(1)
CE    READING OF GENERAL DATA FOR MATERIAL MODELS
      CALL UCMODE(A(LMODEL))
CE    FORMING OF POINTER FOR MATERIAL MODELS
      CALL MATREP(A(LMODEL))
CE    READING OF DATA FOR ALL MATERIALS
      CALL MATERI(A(LMODEL),A(LGUSM))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UCMODE(MODEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ GENERAL DATA FOR MATERIAL MODELS
CS.    P R O G R A M
CS.        ZA UCITAVANJE OPSTIH PODATAKA O MATERIJALNIM MODELIMA
C .
CE.    I=1,NMATM  (NMATM - TOTAL NUMBER OF MATERIAL MODELS /3/)
CS.    I=1,NMATM  (NMATM - UKUPAN BROJ MODELA MATERIJALA)
CE.        /11/ USER MANUAL
C .
CE.           MODEL(1,I) - MATERIAL MODEL TYPE WITH SEQUENCE NUMBER
CE.                        =1; ELASTIC ISOTROPIC
CE.                        =2; ELASTIC ANISOTROPIC
CE.                        =3; THERMO-ELASTIC ISOTROPIC
CE.                        =4; THERMO-ELASTIC ANISOTROPIC
CE.                        =5; VON MISES ELASTIC-PLASTIC WITH
CE.                            ISOTRPIC HARDENING
CE.                        =6; VON MISES ELASTIC-PLASTIC WITH
CE.                            MIXED HARDENING
CE.                        =7; GENERALIZED CAP MODEL
CE.                        =8;
CE.                        =9; SOIL MODEL
CE.                        =10;ELASTIC-PLASTIC ANISOTROPIC
CE.                        =11;GAP-STRESS
CE.                        =12;GAP-FORCE
CE.                        =13;GENERAL ANISOTROPIC ELASTO-PLASTIC
CE.                        =14;VON MISES THERMO-ELASTIC-PLASTIC WITH
CE.                            MIXED HARDENING
CE.                        =15;THERMO-ELASTIC WITH CREEP
CE.                        =16;VON MISES THERMO-ELASTIC-PLASTIC (WITH
CE.                            MIXED HARDENING) AND CREEP
CE.                        =17;VON MISES THERMO-ELASTIC-PLASTIC WITH
CE.                            MIXED HARDENING (ORTHOTROPIC)
CE.                        =18;THERMO-ELASTIC WITH CREEP (ORTHOTROPIC)
CE.                        =19;VON MISES THERMO-ELASTIC-PLASTIC WITH
CE.                            MIXED HARDENING AND CREEP (ORTHOTROPIC)
CE.                        =20;ISOTROPIC VISCO PLASTIC
CE.                        =21;ORTHOTROPIC VISCO PLASTIC
CE.                        =22; SOIL MODEL
CE.                        =23; CAP MODEL - VISCO PLASTIC
CE.                        =25; ROCKS-CST
CE.                        =26; VISCOELASTIC MODEL
CE.                        =30; FIBER - FIBER MODEL
CE.                        =31; HILLS MUSCLE MODEL
CE.                        =32; STRESS-STRETCH MODEL FOR PASSIVE STATE
CE.                        =33; STRESS-STRETCH CREEP MODEL
CE.                        =34; STRESS-STRETCH MODEL FOR ACTIVE STATE
CE.                        =40; GURSON
CE.                        =61; Isotropic Damage Model (Oliver 1996)
C .
CS.           MODEL(1,I) - VRSTA MATERIJALNOG MODELA
CS.                        =1; ELASTICAN IZOTROPAN
CS.                        =2; ELASTICAN ANIZOTROPAN
CS.                        =3; TERMO-ELASTICAN IZOTROPAN
CS.                        =4; TERMO-ELASTICAN ANIZOTROPAN
CS.                        =5; MIZESOV ELASTO-PLASTICAN
CS.                            SA IZOTROPNIM OJACANJEM
CS.                        =6; MIZESOV ELASTO-PLASTICAN
CS.                            SA MESOVITIM OJACANJEM
CS.                        =7; GENERALISANI MODEL TLA SA KAPOM
CS.                        =9; TLO PLASTICAN
CS.                        =10;ELASTO-PLASTICAN ANIZOTROPAN
CS.                        =11;ZAZOR SA ZADATIM NAPONOM
CS.                        =12;ZAZOR SA ZADATOM SILOM
CS.                        =13;GENERALNI ANIZOTROPNI ELASTO-PLASTICAN
CS.                        =14;MIZESOV TERMO-ELASTO-PLASTICAN
CS.                            SA IZOTROPNIM OJACANJEM
CS.                        =15;TERMO-ELASTICAN SA PUZANJEM
CS.                        =16;MIZESOV TERMO-ELASTO-PLASTICAN (SA
CS.                            MESOVITIM OJACANJEM) I PUZANJE
CS.                        =17;MIZESOV TERMO-ELASTO-PLASTICAN SA
CS.                            MESOVITIM OJACANJEM (ORTOTROPAN)
CS.                        =18;TERMO-ELASTICAN SA PUZANJEM (ORTOTROPAN)
CS.                        =19;MIZESOV TERMO-ELASTO-PLASTICAN (SA
CS.                            MESOVITIM OJACANJEM) I PUZANJE (ORTOTROPAN)
CS.                        =20;IZOTROPAN VISKO PLASTICAN
CS.                        =21;ORTOTROPAN VISKO PLASTICAN
CE.                        =22; TLO VISKO-PLASTICNOST
CS.                        =23; MODEL TLA SA KAPOM VISKO-PLASTICNOST
CS.                        =25; STENA-CST
CS.                        =26; VISKOELASTICAN               
CS.                        =30; MODEL VLAKNO - VLAKNO               
CS.                        =31; HILOV MODEL MISICA      
CS.                        =32; NAPON-STREC MODEL ZA PASIVNO STANJE
CS.                        =33; NAPON-STREC MODEL PUZANJA
CS.                        =34; NAPON-STREC MODEL ZA AKTIVNO STANJE
CS.                        =40; GURSON
CS.                        =61; Izotropni Model Ostecenja (Oliver 1996)
C .
CE.           MODEL(2,I) - NUMBER OF DIFFERENT MATERIAL SETS OF MATERIAL
CE.                        MODEL TYPE (=MODEL(1,I))
CS.           MODEL(2,I) - BROJ MATERIJALA OVE VRSTE
C .
CE.           MODEL(3,I) - MAXIMUM NUMBER OF PAIRS, IF TABLES ARE USED
CE.                        TO DEFINE MATERIAL CHARACTERISTICS OF 
CE.                        MATERIAL MODEL TYPE (=MODEL(1,I))
CS.           MODEL(3,I) - MAKSIMALAN BROJ TACAKA ZA KRIVE KARAKTERI-
CS.                        STIKA MATERIJALA OVOG MODELA MATERIJALA
C .
CE.           MODEL(4,I) - POINTER (=LREP) FOR MATERIAL DATA OF MATERIAL
CE.                        MODEL TYPE (=MODEL(1,I)),(SEE ROUTINE MATERI)
CS.           MODEL(4,I) - POCETNI REPER ZA PODATKE O MATERIJALIMA
CS.                        OVOG MODELA MATERIJALA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OBNOVA/ IRUSI
      COMMON /IZADJI/ NASILU
      COMMON /SRPSKI/ ISRPS
      DIMENSION MODEL(4,*)
      COMMON /CDEBUG/ IDEBUG
C INDIKATOR ZA BRISANJE DEFORMACIJA I POMERANJA PRE RESTARTOVANJA - NE VALJA!!!
c     klasican restart =0, sa brisanjem =1
c      IRUSI=1
      IRUSI=0
C INDIKATOR ZA NASILNI IZLAZAK IZ MODELA - RADI SAMO ELASTICNO RASTOJANJE OD LOMA      
c      NASILU=0
      NASILU=0
      if(nasilu.eq.1) write(*,*) ' NASILNO RACUNANJE SAMO ELASTICNO!!!!'
      IF(IDEBUG.GT.0) PRINT *, ' UCMODE'
      DO 10 K=1,NMATM
        CALL ISPITA(ACOZ)
        IF(K.EQ.1) KARTI=KARTIC
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (MODEL(J,K),J=1,3)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) (MODEL(J,K),J=1,3)
   10 CONTINUE
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,1)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      DO 20 I=1,NMATM
        WRITE(IZLAZ,5010) (MODEL(J,I),J=1,3)
   20 CONTINUE
      RETURN
C
 1000 FORMAT(14I5)
 5010 FORMAT(9X,I5,8X,I5,19X,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'P O D A C I    O    M O D E L I M A    M A T E R I J A L A'/
     16X,58('-')///11X,
     2'MODEL    BROJ MATERIJALA    NAJVECI BROJ TACAKA KRIVE')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'D A T A   F O R   M A T E R I A L   M O D E L S'/
     16X,47('-')///11X,
     2'MODEL   NUMBER OF MATERIALS   MAX NUMBER OF POINTS ON CURVE')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MATREP(MODEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO FORM POINTER FOR MATERIAL MODELS
CS.    P R O G R A M
CS.        ZA FORMIRANJE REPERA ZA MODELE MATERIJALA
C .
CE.    P O I N T E R S
CE.        LREP   - MATERIAL DATA
CE.        LGUSM  - DATA FOR MATERIAL MASS DENSITY
CE.    V A R I A B L E S
CE.        NMATM  - TOTAL NUMBER OF DIFFERENT MATERIAL MODELS /3/
CE.          MAT  - NUMBER OF DIFFERENT MATERIAL SETS FOR ONE MATERIAL
CE.                 MODEL /11/
CE.         MATM  - MAXIMUM NUMBER OF DIFFERENT MATERIAL SETS FOR ONE
CE.                 MATERIAL MODEL
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUPLAP/ IDVA
      COMMON /MATERM/ LMODEL,LGUSM
      DIMENSION MODEL(4,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MATREP'
      LMA=LMAX
      MATM=0
C
      DO 200 I=1,NMATM
        MOD=MODEL(1,I)
        MAT=MODEL(2,I)
        MAX=MODEL(3,I)
        IF(MATM.LT.MAT) MATM=MAT
        LREP=LMAX
        CALL DELJIV(LREP,2,INDL)
        IF(INDL.EQ.0) LREP=LREP+1
        GO TO (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     1          11, 11, 13, 14, 15, 14, 17, 18, 17, 20,
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
CE      READING DATA FOR MATERIAL MODEL 1
    1   LMAX=LREP+MAT*2*IDVA
        GO TO 100
C
CE      READING DATA FOR MATERIAL MODEL 2
    2   LMAX=LREP+MAT*9*IDVA
        GO TO 100
C
    3   LMAX=LREP+MAT*(3+(6*MAX+1)*IDVA)
        GO TO 100
C
    4   LMAX=LREP+MAT*(12+(24*MAX+1)*IDVA)
        GO TO 100
C
    5   LMAX=LREP+MAT*(1+(2*MAX+2)*IDVA)
        GO TO 100
C
    6   LMAX=LREP+MAT*(1+(2*MAX+6)*IDVA)
        GO TO 100 
C
    7   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C
    8   LMAX=LREP+MAT*12*IDVA
        GO TO 100
C
    9   LMAX=LREP+MAT*9*IDVA
        GO TO 100
C
   10   LMAX=LREP+MAT*29*IDVA
        GO TO 100
C
   11   LMAX=LREP+MAT*(1+2*MAX*IDVA)
        GO TO 100
C
   13   LMAX=LREP+MAT*77*IDVA
        GO TO 100
C
   14   LMAX=LREP+MAT*(4+(16*MAX+1)*IDVA)
        GO TO 100
C
   15   LMAX=LREP+MAT*(3+(6*MAX+1)*IDVA)
        GO TO 100
C
   17   LMAX=LREP+MAT*(18+(72*MAX+1)*IDVA)
        GO TO 100
C
   18   LMAX=LREP+MAT*(13+(26*MAX+1)*IDVA)
        GO TO 100
C
   20   LMAX=LREP+MAT*(1+(2*MAX+6)*IDVA)
        GO TO 100
C
   21   LMAX=LREP+MAT*30*IDVA
        GO TO 100
C
   22   LMAX=LREP+MAT*10*IDVA
        GO TO 100
C
   23   LMAX=LREP+MAT*34*IDVA
        GO TO 100
C
C       VISKOZNI OTPOR ZA STAP
   24   LMAX = LREP + MAT*5*IDVA
        GO TO 100
C
   25   LMAX=LREP+MAT*30*IDVA
        GO TO 100
C
   26   LMAX=LREP+MAT*(1+2*(MAX+1)*IDVA)
        GO TO 100
C
   27   LMAX=LREP+MAT*6*IDVA
        GO TO 100
C
   28   LMAX=LREP+MAT*30*IDVA
        GO TO 100
C
   29   LMAX=LREP+MAT*10*IDVA
        GO TO 100
C
   30   LMAX=LREP+MAT*14*IDVA
        GO TO 100 
C 
   31   LMAX=LREP+MAT*8*IDVA
        GO TO 100 
C
   32   LMAX=LREP+MAT*5*IDVA
        GO TO 100
C
   33   LMAX=LREP+MAT*5*IDVA
        GO TO 100
C
   34   LMAX=LREP+MAT*5*IDVA
        GO TO 100
C
   40   LMAX=LREP+MAT*(1+(2*MAX+18)*IDVA)
        GO TO 100 
C       Drucker-Prager (rakic)
   41   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C       Mohr-Coulomb (rakic)
   42   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C       Hoek-Brown (rakic)
   43   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C       Generalized Hoek-Brown (rakic)
   44   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C       Maksimovic (rakic)
   45   LMAX=LREP+MAT*11*IDVA
        GO TO 100
C       ARGENTINA
   51   LMAX=LREP+MAT*33*IDVA
        GO TO 100
C       Crystal
   52   LMAX=LREP+MAT*160*IDVA
        GO TO 100   
C       SMA
   53   LMAX=LREP+MAT*24*IDVA
        GO TO 100                
C       SMA VLADA
   54   LMAX=LREP+MAT*17*IDVA
        GO TO 100
C       Isotropic Damage Model (Oliver 1996)
   61   LMAX=LREP+MAT*5*IDVA
        GO TO 100 
C
  999   STOP 'MATREP - 999'
C
  100   MODEL(4,I)=LREP
        IF(LMAX.GT.MTOT) CALL ERROR(1)
  200 CONTINUE
CE    DEFINE SPACE FOR MATERIAL MASS DENSITY
CS    DEFINISANJE PROSTORA ZA GUSTINU MATERIJALA
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LGUSM=LMAX
      NRAZ=100*MATM
      LMAX=LGUSM+NRAZ*IDVA
      CALL CLEAR(A(LGUSM),NRAZ)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MATERI(MODEL,GUSM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR ALL MATERIALS
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALIMA
C .
CE.    V A R I A B L E S
CE.     I=1,NMAT (NMAT - TOTAL NUMBER OF DIFFERENT MATERIAL SETS)
CE.      /12/ USER MANUAL
CE.          MOD  - MATERIAL MODEL FOR WHICH PHYSICAL CHARACTERISTICS
CE.                 ARE INPUT
CE.          MAT  - MATERIAL SET FOR MATERIAL MODEL (MOD), FOR WHICH
CE.                 PHYSICAL CHARACTERISTICS FOLLOW
CE.         GUST  - MATERIAL MASS DENSITY (GUSM(MOD,MAT)=GUST)
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /NIDEAS/ IDEAS
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /MATER9/ SCMS(3,100),S1S3(2,100,10),KRIT(100)
      common /kadamb/ kadmb(100)
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      DIMENSION MODEL(4,*),GUSM(50,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MATERI'
C
      CALL CLEAR(SCMS,300)
      CALL CLEAR(S1S3,2000)
      CALL ICLEAR(KRIT,100)
      CALL ICLEAR(kadmb,100)
      CALL ICLEAR(korz,100*100*3)
      CALL CLEAR(evg,100*100*3)
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005)
      ENDIF
CE    DETERMINATING TOTAL NUMBER OF DIFFERENT MATERIALS
CS    RACUNANJE UKUPNOG BROJA MATERIJALA
      NMAT=0
      DO 110 I=1,NMATM
        NMAT=NMAT+MODEL(2,I)
  110 CONTINUE
CE    CALL ROUTINES FOR READING MATERIAL DATA
CS    POZIVANJE PROGRAMA ZA UCITAVANJE PODATAKA ZA SVE MATERIJALE
      DO 125 I=1,NMAT
        CALL ISPITA(ACOZ)
        KARTI=KARTIC
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) MOD,MAT,GUST,
     1                  (korz(mod,mat,k),k=1,3),(evg(mod,mat,k),k=1,3)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1000) MOD,MAT,GUST,
     1                    (korz(mod,mat,k),k=1,3),(evg(mod,mat,k),k=1,3)
        DO 130 J=1,NMATM
          MO=MODEL(1,J)
          IF(MO.EQ.MOD) GO TO 140
  130   CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MOD
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MOD
      STOP 'PROGRAM STOP: PAK03 - MATERI'
C
  140   GO TO (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     1          11, 11, 13, 14, 15, 16, 17, 18, 19,  6,
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
CE    READING DATA FOR MATERIAL MODEL 1
    1 LREP=MODEL(4,J)
      CALL UMOD01(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
CE    READING DATA FOR MATERIAL MODEL 2
    2 LREP=MODEL(4,J)
      CALL UMOD02(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    3 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*6*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD03(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    4 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*24*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD04(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    5 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+2)*IDVA
      CALL UMOD05(A(LREP),A(LNTA),MATE,MAXT,MOD,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    6 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+6)*IDVA
      CALL UMOD05(A(LREP),A(LNTA),MATE,MAXT,MOD,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    7 LREP=MODEL(4,J)
      CALL UMOD07(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    8 LREP=MODEL(4,J)
      CALL UMOD08(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
    9 LREP=MODEL(4,J)
      CALL UMOD09(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
   10 LREP=MODEL(4,J)
C      CALL UMOD10(A(LREP),MAT,KARTI)
      CALL UMOD1U(A(LREP),MAT,KARTI)
      GO TO 120
C
   11 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*2*MAXT*IDVA
      IF(MOD.EQ.11) CALL UMOD11(A(LREP),A(LNTA),MATE,MAXT,MAT,KARTI,1)
      IF(MOD.EQ.12) CALL UMOD11(A(LREP),A(LNTA),MATE,MAXT,MAT,KARTI,2)
      GO TO 120
C
   13 LREP=MODEL(4,J)
      CALL UMOD13(A(LREP),MAT,KARTI)
      GO TO 120
C
   14 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*16*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD14(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
   15 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*6*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD15(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   16 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*16*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD16(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   17 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*72*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD17(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   18 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*26*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD18(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   19 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LTEM=LREP+MATE*72*MAXT*IDVA
      LNTA=LTEM+MATE*IDVA
      CALL UMOD19(A(LREP),A(LTEM),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   21 LREP=MODEL(4,J)
      CALL UMOD21(A(LREP),MAT,KARTI)
      GO TO 120
C
   22 LREP=MODEL(4,J)
      CALL UMOD22(A(LREP),MAT,KARTI)
      GO TO 120
C
   23 LREP=MODEL(4,J)
      CALL UMOD23(A(LREP),MAT,KARTI)
      GO TO 120
C
C     VISKOZNI ZA STAP
   24 LREP = MODEL(4,J)
      CALL UMOD24(A(LREP),MAT,KARTI)
      GO TO 120
C
   25 LREP=MODEL(4,J)
      CALL UMOD25(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
   26 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*2*(MAXT+1)*IDVA
      CALL UMOD26(A(LREP),A(LNTA),MATE,MAXT,MAT,KARTI)
      GO TO 120
C
   27 LREP=MODEL(4,J)
      CALL UMOD27(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
   28 LREP=MODEL(4,J)
      CALL UMOD28(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
	write(3,*)'lrepul',lrep
      GO TO 120
C
   29 LREP=MODEL(4,J)
      CALL UMOD29(A(LREP),MAT,KARTI)
      GO TO 120
C
   30 LREP=MODEL(4,J)
      CALL UMOD30(A(LREP),MAT,KARTI)
      GO TO 120
C
   31 LREP=MODEL(4,J)
      CALL UMOD31(A(LREP),MAT,KARTI)
      GO TO 120
C
   32 LREP=MODEL(4,J)
      CALL UMOD32(A(LREP),MAT,KARTI)
      GO TO 120
C
   33 LREP=MODEL(4,J)
      CALL UMOD33(A(LREP),MAT,KARTI)
      GO TO 120
C
   34 LREP=MODEL(4,J)
      CALL UMOD34(A(LREP),MAT,KARTI)
      GO TO 120
C
   40 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      MAXT=MODEL(3,J)
      LNTA=LREP+MATE*(2*MAXT+18)*IDVA
      CALL UMOD40(A(LREP),A(LNTA),MATE,MAXT,MOD,MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C     Drucker-Prager (rakic)
   41 LREP=MODEL(4,J)
      CALL UMOD41(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
c     Mohr-Coulomb (rakic)
   42 LREP=MODEL(4,J)
      CALL UMOD42(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
c     Hoek-Brown (rakic)
   43 LREP=MODEL(4,J)
      CALL UMOD43(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
c     Generalized Hoek-Brown (rakic)
   44 LREP=MODEL(4,J)
      CALL UMOD44(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
c     Maksimovic (rakic)
   45 LREP=MODEL(4,J)
      CALL UMOD45(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C ARGENTINA
   51 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      CALL UMOD51(A(LREP),MAT,KARTI,MATE)
      GO TO 120
C     Crystal      
   52 LREP=MODEL(4,J)
      CALL UMOD52(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120     
C     SMA      
   53 LREP=MODEL(4,J)
      CALL UMOD53(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120      
C     SMA VLADA
   54 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      CALL UMOD54(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120      
C     Isotropic Damage Model (Oliver 1996)
   61 LREP=MODEL(4,J)
      MATE=MODEL(2,J)
      CALL UMOD61(A(LREP),MAT,KARTI,
     +            GUST,NBLGR,IDEAS,I)
      GO TO 120
C
  999 STOP 'MATERI - 999'
C
CE    MATERIAL MASS DENSITY
  120 GUSM(MOD,MAT)=GUST
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) GUST
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) GUST
      if(korz(mod,mat,1).ne.0.or.korz(mod,mat,2).ne.0.or.
     1                           korz(mod,mat,3).ne.0) then
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2200) (korz(mod,mat,k),k=1,3),(evg(mod,mat,k),k=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,2200) (korz(mod,mat,k),k=1,3),(evg(mod,mat,k),k=1,3)
      ENDIF
      ENDIF
  125 CONTINUE
      RETURN
C
 1000 FORMAT(2I5,F10.2,3i5,3f10.0)
 2200 FORMAT(' tunel'/3I5,3(1pe12.4))
C-----------------------------------------------------------------------
 2100 FORMAT(/
     116X,'G U S T I N A ....................... GUST =',1PD12.5)
 2005 FORMAT(//////6X,
     1'P  O  D  A  C  I     O     M  A  T  E  R  I  J  A  L  I  M  A'/
     16X,61('-'))
 2000 FORMAT(////' P O D A C I   O   M O D E L U   M A T E R I J A L A
     1 B R O J  =',I5/' N I S U   U C I T A N I'//' PROGRAM STOP'//)
C-----------------------------------------------------------------------
 6100 FORMAT(/
     116X,'D E N S I T Y ....................... GUST =',1PD12.5)
 6005 FORMAT(//////6X,
     1'D  A  T  A     F  O  R     M  A  T  E  R  I  A  L  S'/
     16X,52('-'))
 6000 FORMAT(////' DATA IS NOT READ FOR MATERIAL MODEL =',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD01(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL TYPE 1
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 1
C .
CE.    V A R I A B L E S
CE.      /12-1/ USER MANUAL
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS  - E
CS.   FUNMAT(1,MAT) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(2,MAT) - POISSON*S RATIO - V
CS.   FUNMAT(2,MAT) - POISSONOV BROJ - V
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(2,*),AMAT(30)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD01'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(2,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(2,MAT)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     1  (LINEARAN ELASTICAN IZOTROPAN)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I ... E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J ........... V =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     1  (LINEAR ELASTIC ISOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S    M O D U L U S ........... E =',1PD12.5//
     216X,'P O I S S O N S    R A T I O ........... V =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD02(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 2
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 2
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  EX
CE.   FUNMAT(2,MAT) - YOUNG*S MODULUS IN DIRECTION Y  -  EY
CE.   FUNMAT(3,MAT) - YOUNG*S MODULUS IN DIRECTION Z  -  EZ
CS.   FUNMAT(1,MAT) - MODUL ELASTICNOSTI U PRAVCU X OSE - EX
CS.   FUNMAT(2,MAT) - MODUL ELASTICNOSTI U PRAVCU Y OSE - EY
CS.   FUNMAT(3,MAT) - MODUL ELASTICNOSTI U PRAVCU Z OSE - EZ
C .
CE.   FUNMAT(4,MAT) - POISSON*S RATIO FOR PLANE XY    -  VXY
CE.   FUNMAT(5,MAT) - POISSON*S RATIO FOR PLANE YZ    -  VYZ
CE.   FUNMAT(6,MAT) - POISSON*S RATIO FOR PLANE ZX    -  VZX
CS.   FUNMAT(4,MAT) - POISSONOV BROJ U PRAVCU X I Y OSE - VXY
CS.   FUNMAT(5,MAT) - POISSONOV BROJ U PRAVCU Y I Z OSE - VYZ
CS.   FUNMAT(6,MAT) - POISSONOV BROJ U PRAVCU Z I X OSE - VZX
C .
CE.   FUNMAT(7,MAT) - SHEAR MODULUS FOR PLANE XY  -  GXY
CE.   FUNMAT(8,MAT) - SHEAR MODULUS FOR PLANE YZ  -  GYZ
CE.   FUNMAT(9,MAT) - SHEAR MODULUS FOR PLANE ZX  -  GZX
CS.   FUNMAT(7,MAT) - MODUL SMICANJA U PRAVCU X I Y OSE - GXY
CS.   FUNMAT(8,MAT) - MODUL SMICANJA U PRAVCU Y I Z OSE - GYZ
CS.   FUNMAT(9,MAT) - MODUL SMICANJA U PRAVCU Z I X OSE - GZX
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(9,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD02'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
C
C... CHECK MATERIAL CONSTANTS
      CALL ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(2)=FUNMAT(2,MAT)
         AMAT(3)=FUNMAT(3,MAT)
         AMAT(4)=FUNMAT(7,MAT)
         AMAT(5)=FUNMAT(8,MAT)
         AMAT(6)=FUNMAT(9,MAT)
         AMAT(7)=FUNMAT(4,MAT)
         AMAT(8)=FUNMAT(5,MAT)
         AMAT(9)=FUNMAT(6,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,9)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,9)
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     2  (LINEARAN ELASTICAN ANIZOTROPAN)'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     2  (LINEAR ELASTIC ORTHOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U S'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD03(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 3
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 3
C .
CE.   M=(MAT-1)*3    (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*3    (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J)   - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J)   - YOUNGS  MODULUS - E
CS.   FUNMAT(1,M+1,J)   - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J)   - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - POISSONS  NUMBER - V
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - POISSONOV BROJ - V
C .
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - COEFFICIENT OF THERMAL EXPANSION  - A
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - A
C .
CE.         TREF(MAT) - REFERENCE TEMPERATURE
CS.         TREF(MAT) - REFERENTNA TEMPERATURA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS FOR ELASTIC MODULUS - E
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - E
C .
CE.      NTFUN(M+2)   - NUMBER OF POINT TO DEFINE POISSON*S NUMBER  - V
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - V
C .
CE.      NTFUN(M+3)   - NUMBER OF POINT FOR CURVE - A
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - A
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(2,MATE*3,*),TREF(*),NTFUN(*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD03'
      MA=(MAT-1)*3
      DO 10 K=1,3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)                              
C        write(3,*) 'NTFUN(MA+K),MATE',NTFUN(MA+K),MATE
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD03'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,2)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,2)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(2,MA+1,1)
         AMAT(7)=FUNMAT(2,MA+2,1)
         AMAT(10)=FUNMAT(2,MA+3,1)
         AMAT(14)=TREF(MAT)
         AMAT(13)=GUST
         AMAT(19)=1.
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MATG,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MA+1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MA+1)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+1,J),I=1,2),J=1,NTFUN(MA+1))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MA+2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MA+2)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+2,J),I=1,2),J=1,NTFUN(MA+2))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NTFUN(MA+3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NTFUN(MA+3)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+3,J),I=1,2),J=1,NTFUN(MA+3))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      RETURN
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     3  (LINEARAN TERMO-ELASTICAN',
     2' IZOTROPAN)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  -  E'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               E')
 2030 FORMAT(//
     116X,'P O I S S O N O V    B R O J  -  V'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               V')
 2040 FORMAT(//
     116X,'K O E F I C I J E N T    S I R E N J A  -  A'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               A')
 2050 FORMAT(//
     116X,'R E F E R E N T N A   TEMPERATURA ... TREF =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     3  (LINEAR THERMO-ELASTIC',
     2' ISOTROPIC)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S      M O D U L U S  -  E'//
     121X,'NUMBER OF POINTS FOR CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               E')
 6030 FORMAT(//
     116X,'P O I S S O N S      N U M B E R  -  V'//
     121X,'NUMBER OF POINTS FOR CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               V')
 6040 FORMAT(//
     116X,'COEFFICIENT OF THERMAL EXPANSION  -  A'//
     121X,'NUMBER OF POINTS FOR CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               A')
 6050 FORMAT(//
     116X,'R E F E R E N C E    TEMPERATURE .... TREF =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD04(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 4
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 4
C .
CE.   M=(MAT-1)*12   (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*12   (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J)   - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J)   - YOUNG*S MODULUS - EX
CE.   FUNMAT(1,M+2,J)   - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J)   - YOUNG*S MODULUS - EY
CE.   FUNMAT(1,M+3,J)   - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J)   - YOUNG*S MODULUS - EZ
CS.   FUNMAT(1,M+1,J)   - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J)   - MODUL ELASTICNOSTI - EX
CS.   FUNMAT(1,M+2,J)   - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J)   - MODUL ELASTICNOSTI - EY
CS.   FUNMAT(1,M+3,J)   - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J)   - MODUL ELASTICNOSTI - EZ
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - POISSON*S RATIO - VXY
CE.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+5,J) - POISSON*S RATIO - VYZ
CE.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+6,J) - POISSON*S RATIO - VZX
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - POISSONOV BROJ - VXY
CS.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+5,J) - POISSONOV BROJ - VYZ
CS.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+6,J) - POISSONOV BROJ - VZX
C .
CE.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+7,J) - SHEAR MODULUS - GXY
CE.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+8,J) - SHEAR MODULUS - GYZ
CE.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+9,J) - SHEAR MODULUS - GZX
CS.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+7,J) - MODUL SMICANJA - GXY
CS.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+8,J) - MODUL SMICANJA - GYZ
CS.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+9,J) - MODUL SMICANJA - GZX
C .
CE.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+10,J) - COEFFICIENT OF THERMAL EXPANSION - AX
CE.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+11,J) - COEFFICIENT OF THERMAL EXPANSION - AY
CE.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+12,J) - COEFFICIENT OF THERMAL EXPANSION  - AZ
CS.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+10,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AX
CS.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+11,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AY
CS.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+12,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AZ
C .
CE.         TREF(MAT) - REFERENCE TEMPERATURE
CS.         TREF(MAT) - REFERENTNA TEMPERATURA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS - EX
CE.      NTFUN(M+2)   - NUMBER OF POINTS - EY
CE.      NTFUN(M+3)   - NUMBER OF POINTS - EZ
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - EX
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - EY
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - EZ
C .
CE.      NTFUN(M+4)   - NUMBER OF POINTS - VXY
CE.      NTFUN(M+5)   - NUMBER OF POINTS - VYZ
CE.      NTFUN(M+6)   - NUMBER OF POINTS - VZX
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVU - VXY
CS.      NTFUN(M+5)   - BROJ TACAKA ZA KRIVU - VYZ
CS.      NTFUN(M+6)   - BROJ TACAKA ZA KRIVU - VZX
C .
CE.      NTFUN(M+7)   - NUMBER OF POINTS - GXY
CE.      NTFUN(M+8)   - NUMBER OF POINTS - GYZ
CE.      NTFUN(M+9)   - NUMBER OF POINTS - GZX
CS.      NTFUN(M+7)   - BROJ TACAKA ZA KRIVU - GXY
CS.      NTFUN(M+8)   - BROJ TACAKA ZA KRIVU - GYZ
CS.      NTFUN(M+9)   - BROJ TACAKA ZA KRIVU - GZX
C .
CE.      NTFUN(M+10)  - NUMBER OF POINTS - AX
CE.      NTFUN(M+11)  - NUMBER OF POINTS - AY
CE.      NTFUN(M+12)  - NUMBER OF POINTS - AZ
CS.      NTFUN(M+10)  - BROJ TACAKA ZA KRIVU - AX                      .
CS.      NTFUN(M+11)  - BROJ TACAKA ZA KRIVU - AY                      .
CS.      NTFUN(M+12)  - BROJ TACAKA ZA KRIVU - AZ                      .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(2,MATE*12,*),TREF(*),NTFUN(*),AMAT(30)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD04'
      DO 10 K=1,12
        MA=(MAT-1)*12
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD04'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,2)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,2)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(2,MA+1,1)
         AMAT(2)=FUNMAT(2,MA+2,1)
         AMAT(3)=FUNMAT(2,MA+3,1)
         AMAT(4)=FUNMAT(2,MA+7,1)
         AMAT(5)=FUNMAT(2,MA+8,1)
         AMAT(6)=FUNMAT(2,MA+9,1)
         AMAT(7)=FUNMAT(2,MA+4,1)
         AMAT(8)=FUNMAT(2,MA+5,1)
         AMAT(9)=FUNMAT(2,MA+6,1)
         AMAT(10)=FUNMAT(2,MA+10,1)
         AMAT(11)=FUNMAT(2,MA+11,1)
         AMAT(12)=FUNMAT(2,MA+12,1)
         AMAT(14)=TREF(MAT)
         AMAT(13)=GUST
         AMAT(19)=1.
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      MA=(MAT-1)*12
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MA+1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MA+1)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+1,J),I=1,2),J=1,NTFUN(MA+1))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) NTFUN(MA+2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) NTFUN(MA+2)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+2,J),I=1,2),J=1,NTFUN(MA+2))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060) NTFUN(MA+3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060) NTFUN(MA+3)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+3,J),I=1,2),J=1,NTFUN(MA+3))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MA+4)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MA+4)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+4,J),I=1,2),J=1,NTFUN(MA+4))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) NTFUN(MA+5)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) NTFUN(MA+5)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+5,J),I=1,2),J=1,NTFUN(MA+5))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2080) NTFUN(MA+6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6080) NTFUN(MA+6)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+6,J),I=1,2),J=1,NTFUN(MA+6))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NTFUN(MA+7)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NTFUN(MA+7)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+7,J),I=1,2),J=1,NTFUN(MA+7))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2090) NTFUN(MA+8)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6090) NTFUN(MA+8)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+8,J),I=1,2),J=1,NTFUN(MA+8))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NTFUN(MA+9)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NTFUN(MA+9)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+9,J),I=1,2),J=1,NTFUN(MA+9))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2110) NTFUN(MA+10)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6110) NTFUN(MA+10)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+10,J),I=1,2),J=1,NTFUN(MA+10))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2120) NTFUN(MA+11)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6120) NTFUN(MA+11)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+11,J),I=1,2),J=1,NTFUN(MA+11))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2130) NTFUN(MA+12)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6130) NTFUN(MA+12)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+12,J),I=1,2),J=1,NTFUN(MA+12))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2140) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6140) TREF(MAT)
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     4  (LINEARAN TERMO-ELASTICAN',
     2' ANIZOTROPAN)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  -  EX'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               EX')
 2050 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  -  EY'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               EY')
 2060 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  -  EZ'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               EZ')
 2030 FORMAT(//
     116X,'P O I S S O N O V    B R O J  -  VXY'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              VXY')
 2070 FORMAT(//
     116X,'P O I S S O N O V    B R O J  -  VYZ'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              VYZ')
 2080 FORMAT(//
     116X,'P O I S S O N O V    B R O J  -  VZX'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              VZX')
 2040 FORMAT(//
     116X,'M O D U L    S M I C A N J A  -  GXY'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              GXY')
 2090 FORMAT(//
     116X,'M O D U L    S M I C A N J A  -  GYZ'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              GYZ')
 2100 FORMAT(//
     116X,'M O D U L    S M I C A N J A  -  GZX'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA              GZX')
 2110 FORMAT(//
     116X,'K O E F I C I J E N T    S I R E N J A  -  AX'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               AX')
 2120 FORMAT(//
     116X,'K O E F I C I J E N T    S I R E N J A  -  AY'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               AY')
 2130 FORMAT(//
     116X,'K O E F I C I J E N T    S I R E N J A  -  AZ'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT              FUNKCIJA'/
     120X,'TEMPERATURA               AZ')
 2140 FORMAT(//
     116X,'R E F E R E N T N A   TEMPERATURA ... TREF =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     4  (LINEAR THERMO-ELASTIC',
     2' ORTHOTROPIC)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S      M O D U L U S  -  EX'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               EX')
 6050 FORMAT(//
     116X,'Y O U N G S      M O D U L U S  -  EY'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               EY')
 6060 FORMAT(//
     116X,'Y O U N G S      M O D U L U S  -  EZ'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               EZ')
 6030 FORMAT(//
     116X,'P O I S S O N S      R A T I O  -  VXY'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE              VXY')
 6070 FORMAT(//
     116X,'P O I S S O N S      R A T I O  -  VYZ'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURA              VYZ')
 6080 FORMAT(//
     116X,'P O I S S O N S      R A T I O  -  VZX'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE              VZX')
 6040 FORMAT(//
     116X,'S H E A R    M O D U L U S  -  GXY'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE              GXY')
 6090 FORMAT(//
     116X,'S H E A R    M O D U L U S  -  GYZ'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE              GYZ')
 6100 FORMAT(//
     116X,'S H E A R    M O D U L U S  -  GZX'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE              GZX')
 6110 FORMAT(//
     116X,'COEFFICIENT OF THERMAL EXPANSION  -  AX'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               AX')
 6120 FORMAT(//
     116X,'COEFFICIENT OF THERMAL EXPANSION  -  AY'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               AY')
 6130 FORMAT(//
     116X,'COEFFICIENT OF THERMAL EXPANSION  -  AZ'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT              FUNCTION'/
     120X,'TEMPERATURE               AZ')
 6140 FORMAT(//
     116X,'R E F E R E N C E   TEMPERATURE ..... TREF =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD05(FUNMAT,NTFUN,MATE,MAXT,MOD,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 5,6
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 5,6
C .
CE.               MAT - MATERIAL NUMBER
CE.   J=2,(NTMAX - NUMBER OF POINTS FOR CURVE) + 1
CS.               MAT - MATERIJAL BROJ
CS.   J=2,(NTMAX - BROJ TACAKA ZA KRIVU) + 1
C .
CE.   FUNMAT(1,MAT,1) - YOUNG*S MODULUS - E
CE.   FUNMAT(2,MAT,1) - POISSON*S NUMBER - V
CS.   FUNMAT(1,MAT,1) - MODUL ELASTICNOSTI - E
CS.   FUNMAT(2,MAT,1) - POISSONOV BROJ - V
C .
CE. A) NTMAX=1 (NTMAX - NUMBER OF POINTS FOR CURVE TAU(DEF))
CS. A) NTMAX=1 (NTMAX - BROJ TACAKA ZA KRIVU TAU(DEF))
C .
CE.   FUNMAT(1,MAT,2) - YIELD STRESS - TAUY
CE.   FUNMAT(2,MAT,2) - TANGENTIAL MODULUS - ET
CS.   FUNMAT(1,MAT,2) - NAPON TECENJA - TAUY
CS.   FUNMAT(2,MAT,2) - TANGENTNI MODUL - ET
C .
CE. B) NTMAX>1 (NTMAX - NUMBER OF POINTS FOR CURVE TAU(DEF))
CS. B) NTMAX>1 (NTMAX - BROJ TACAKA ZA KRIVU TAU(DEF))
C .
CE.   FUNMAT(1,MAT,J) - ARGUMENT - STRAIN - DEF
CE.   FUNMAT(2,MAT,J) - STRESS - TAU
CS.   FUNMAT(1,MAT,J) - ARGUMENT - DEFORMACIJA - DEF
CS.   FUNMAT(2,MAT,J) - NAPON - TAU
C .
CE.        NTFUN(MAT) - NUMBER OF POINTS FOR CURVE - TAU(DEF)
CE.                     =1; BILINEAR CURVE
CE.                     >1; MULTILINEAR CURVE
CS.        NTFUN(MAT) - BROJ TACAKA ZA KRIVU - TAU(DEF)
CS.                     =1; BILINEARNA KRIVA
CS.                     >1; MULTILINEARNA KRIVA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(2,MATE,*),NTFUN(*),AMAT(30)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD05'
      ISIMO=0
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NTFUN(MAT),INDEP
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NTFUN(MAT),INDEP
        IF(NTFUN(MAT).GT.MAXT.OR.NTFUN(MAT).LE.0) STOP 'UMOD05'
      NBRT=NTFUN(MAT)
      IF (NBRT.GT.1 .AND. INDEP.EQ.0 ) STOP 'INDEP - PAK03'
      IF (NBRT.GT.MAXT) STOP 'MAXT - PAK03'
      IPROM=2
      IF (NBRT.GT.1) IPROM=2+NBRT
      FUNMAT(1,MAT,4)=0.
      DO 20 J=1,IPROM
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (FUNMAT(I,MAT,J),I=1,2),DUM
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1010) (FUNMAT(I,MAT,J),I=1,2),DUM
        IF(MOD.EQ.6.AND.DABS(DUM).GT.1.D-10.AND.IPROM.EQ.2)THEN
          FUNMAT(1,MAT,4)=DUM
          ISIMO=1
        ENDIF
   20 CONTINUE
      IF((MOD.EQ.6.OR.MOD.EQ.20).AND.NTFUN(MAT).EQ.1)THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (FUNMAT(I,MAT,3),I=1,2)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1010) (FUNMAT(I,MAT,3),I=1,2)
        IF(DABS(FUNMAT(1,MAT,3)).LT.1.D-6) FUNMAT(1,MAT,3)=1.D0
      ENDIF
      IF(MOD.EQ.20.AND.NTFUN(MAT).EQ.1)THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) FUNMAT(1,MAT,4)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1010) FUNMAT(1,MAT,4)
      ENDIF
      IF(MOD.EQ.6 .AND. NBRT.GT.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MAT),INDEP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MAT),INDEP
      WRITE(IZLAZ,5020) ((FUNMAT(I,MAT,J),I=1,2),J=3,IPROM)
      CALL FICA(FUNMAT,MAT,MATE,IPROM,INDEP)
      NTFUN(MAT)=1
      ENDIF
      IF(MOD.EQ.6.AND.FUNMAT(2,MAT,2).LT.1.D-8) FUNMAT(2,MAT,2)=1.D-8
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) GO TO 90
      CALL WBROJK(KARTI,0)
      IF(MOD.EQ.5) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      ENDIF
      IF(MOD.EQ.6) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      ENDIF
      IF(MOD.EQ.20) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2006) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6006) MAT
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(I,MAT,1),I=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(I,MAT,1),I=1,2)
C
      IF(NTFUN(MAT).EQ.1.AND.MOD.EQ.5) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) (FUNMAT(I,MAT,2),I=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) (FUNMAT(I,MAT,2),I=1,2)
      ENDIF
C
      IF(NTFUN(MAT).EQ.1.AND.(MOD.EQ.6.OR.MOD.EQ.20)) THEN
      IF(ISIMO.EQ.0) THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2045) ((FUNMAT(I,MAT,J),I=1,2),J=2,3)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6045) ((FUNMAT(I,MAT,J),I=1,2),J=2,3)
C
      ELSEIF(ISIMO.EQ.1) THEN
        IF(ISRPS.EQ.0)THEN
          WRITE(IZLAZ,2075) (FUNMAT(I,MAT,2),I=1,2),FUNMAT(1,MAT,4),
     1                      (FUNMAT(I,MAT,3),I=1,2)
        ELSEIF(ISRPS.EQ.1)THEN
          WRITE(IZLAZ,6075) (FUNMAT(I,MAT,2),I=1,2),FUNMAT(1,MAT,4),
     1                      (FUNMAT(I,MAT,3),I=1,2)
        ENDIF
        FUNMAT(2,MAT,2)=FUNMAT(2,MAT,2)-FUNMAT(1,MAT,2)
      ENDIF
C
       IF(MOD.EQ.20) THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2060) FUNMAT(1,MAT,4)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6060) FUNMAT(1,MAT,4)
       ENDIF
      ENDIF
C
   90 IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT,1)
         AMAT(7)=FUNMAT(2,MAT,1)
         AMAT(13)=GUST
         AMAT(15)=FUNMAT(1,MAT,2)
         AMAT(16)=FUNMAT(2,MAT,2)
         AMAT(17)=FUNMAT(2,MAT,3)
         AMAT(19)=3.
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NTFUN(MAT).EQ.1) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MAT),INDEP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MAT),INDEP
      WRITE(IZLAZ,5020) ((FUNMAT(I,MAT,J),I=1,2),J=4,IPROM)
C
CE    TRANSFORMATION OF CURVE ON FORM STRESS-PLASTIC STRAIN
CS    PREVODJENJE KRIVE NA OBLIK NAPON-PLASTICNA DEFORMACIJA
C
      IF(INDEP.EQ.1) CALL NAPPLD(FUNMAT,NTFUN,MATE,MAT)
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     5'/6X,'(MIZESOV ELASTO-PLASTICAN',
     2' SA IZOTROPNIM OJACANJEM)'///11X,'MATERIJAL BROJ =',I5)
 2005 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     6'/6X,'(MIZESOV ELASTO-PLASTICAN',
     2' SA MESOVITIM OJACANJEM)'///11X,'MATERIJAL BROJ =',I5)
 2006 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    20'/6X,'(MIZESOV VISKO-PLASTICAN',
     2' SA MESOVITIM OJACANJEM)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5)
 2040 FORMAT(//
     116X,'N A P O N    T E C E N J A  ........  TAUY =',1PD12.5//
     216X,'T A N G E N T N I    M O D U L  ......  ET =',1PD12.5)
 2045 FORMAT(//6X,
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'N A P O N    T E C E N J A  ........  TAUY =',1PD12.5//
     216X,'M N O Z I O C  .....................    CY =',1PD12.5//
     316X,'E K S P O N E N T  .................    AN =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5)
 2060 FORMAT(//
     116X,'K O E F   V I S K O Z N O S T I  ...   ETA =',1PD12.5)
 2030 FORMAT(//
     116X,'K R I V A    N A P O N  -  D E F O R M A C I J A'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'VRSTA DEFORMACIJE , INDEP=',I5/
     126X,'.EQ.1; UKUPNA'/
     126X,'.EQ.2; PLASTICNA'//
     121X,'ARGUMENT             FUNKCIJA'/
     120X,'DEFORMACIJA             NAPON')
 2075 FORMAT(//6X,
     1'   J.C.SIMO FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'NAPON TECENJA .....................   TAUY =',1PD12.5//
     216X,'ASIMPTOTSKI NAPON ................. TAUINF =',1PD12.5//
     216X,'KOEFICIJENT LINEARNOG OJACANJA ....      H =',1PD12.5//
     316X,'EKSPONENT .........................  DELTA =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ...     EM =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =    5'/6X,'(VON MISES ELASTIC-PLASTIC',
     2' WITH ISOTROPIC HARDENING)'///11X,'MATERIAL NUMBER =',I5)
 6005 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =    6'/6X,'(VON MISES ELASTIC-PLASTIC',
     2' WITH MIXED HARDENING)'///11X,'MATERIAL NUMBER =',I5)
 6006 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =   20'/6X,'(VON MISES VISCO-PLASTIC',
     2' WITH MIXED HARDENING)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5)
 6040 FORMAT(//
     116X,'Y I E L D    S T R E S S  ..........  TAUY =',1PD12.5//
     216X,'T A N G E N T    M O D U L U S  ......  ET =',1PD12.5)
 6045 FORMAT(//6X,
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     116X,'Y I E L D    S T R E S S  ..........  TAUY =',1PD12.5//
     216X,'M U L T I P L A Y E R  .............    CY =',1PD12.5//
     316X,'E X P O N E N T  ...................    AN =',1PD12.5///
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5)
 6060 FORMAT(//
     116X,'V I S C O S I T Y   C O E F F  .....   ETA =',1PD12.5)
 6030 FORMAT(//
     116X,'S T R E S S - S T R A I N    C U R V E'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'DEFORMATION SHAPE , INDEP=',I5/
     126X,'.EQ.1; TOTAL STRAIN'/
     126X,'.EQ.2; PLASTIC STRAIN'//
     121X,'ARGUMENT             FUNCTION'/
     120X,'STRAIN                 STRESS')
 6075 FORMAT(//6X,
     1'   J.C.SIMO FORMULA MATERIAL CONSTANTS'//
     116X,'YIELD STRESS  .....................   TAUY =',1PD12.5//
     216X,'ASYMPTOTIC STRESS ................. TAUINF =',1PD12.5//
     216X,'LINEAR HARDENING COEFFICIENT ......      H =',1PD12.5//
     316X,'EXPONENT ..........................  DELTA =',1PD12.5///
     416X,'MIXED HARDENING COEFFICIENT  ......     EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD07(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 7
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 7
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.   COEFFICIENTS FOR DEFINITION CURVE F1=0.0
CE.   FUNMAT(3,MAT) - COEFFICIENT   k 
CE.   FUNMAT(4,MAT) - COEFFICIENT   B 
CE.   FUNMAT(5,MAT) - COEFFICIENT   A 
CE.   FUNMAT(6,MAT) - COEFFICIENT   alfa
CE.   FUNMAT(7,MAT) - TENSION CUTOFF   T
CE.   FUNMAT(8,MAT) - INITIAL CAP POSITION   Lo
C
CE.   DATA FOR HARDENING FUNCTION epv 
CE.   FUNMAT(9,MAT)  - COEFFICIENT   W 
CE.   FUNMAT(10,MAT) - COEFFICIENT   D 
CE.   FUNMAT(11,MAT) - SEMIAXIS RATIO  R
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD07'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,8)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,8)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=9,11)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,11)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=3,8)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=3,8)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2065) (FUNMAT(J,MAT),J=9,11)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=9,11)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     7 (TLO PLASTICNOST-MODEL SA KAPOM)'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F1=0 '/21X,'A',13X,'B ',
     112X,'C ',12X,'TETA'/16X,4(1PD12.5,2X)//
     120X,'AT  ',11X,'AI1A0 '/16X,2(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'
     1/20X,'W ',13X,'D ',11X,'   R'/16X,3(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    7 (SOIL PLASTICITY - GENERAL CAP)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F1=0 '
     1/21X,'A',13X,'B ',12X,'C ',12X,'TETA'/16X,4(1PD12.5,2X)//
     120X,'AT  ',11X,'AI1A0 ',9X/16X,2(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/20X,'W ',13X,'D ',11X,'   R'/16X,3(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD08(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 9
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 9
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - AEM
CE.   FUNMAT(2,MAT) - CONSTANT  - ALAM
CE.   FUNMAT(3,MAT) - CONSTANT  - AKA
CE.   FUNMAT(4,MAT) - CONSTANT  - AE1
CE.   FUNMAT(5,MAT) - CONSTANT  - AE0
CS.   FUNMAT(1,MAT) - KONSTANTA - AEM
CS.   FUNMAT(2,MAT) - KONSTANTA - ALAM
CS.   FUNMAT(3,MAT) - KONSTANTA - AKA
CS.   FUNMAT(4,MAT) - KONSTANTA - AE1
CS.   FUNMAT(5,MAT) - KONSTANTA - AE0
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(12,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD08'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,12)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,12)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=7,12)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=7,12)
      RETURN
 1000 FORMAT(6F10.0)
 1010 FORMAT(6F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,  
     1'MODEL MATERIJALA BROJ =     8 (TLO ANIZOTROPNI CAM-CLAY)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//6X,
     2/6X,'AEM',8X,'ALAM',9X,'AKA',8X,'AE0 ',8X,'G',8X,'BKP0'/
     3' ',6(1PD12.4))
 2045 FORMAT(//6X,
     2/6X,'AK0',8X,'D',9X,'ANN',8X,'AMM ',8X,'GAMA',8X,'SIGX0'/
     3' ',6(1PD12.4))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     8  (SOIL ANISOTROPIC PLASTICITY',
     2' (CAM-CLAY))'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//6X,
     2/6X,'AEM',8X,'ALAM',9X,'AKA',8X,'AE0 ',8X,'G',8X,'BKP0'/
     3' ',6(1PD12.4))
 6045 FORMAT(//6X,
     2/6X,'AK0',8X,'D',9X,'ANN',8X,'AMM ',8X,'GAMA',8X,'SIGX0'/
     3' ',6(1PD12.4))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD09(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 9
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 9
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
CE.   FUNMAT(3,MAT) - TENSION CUTOFF                  - AT
CE.   FUNMAT(4,MAT) - INDICATOR FOR ELASTICITY DEF.   - IEL 
C
CE.   FUNMAT(5,MAT) - CONSTANT  - AEM
CE.   FUNMAT(6,MAT) - CONSTANT  - ALAM
CE.   FUNMAT(7,MAT) - CONSTANT  - AKA
CE.   FUNMAT(8,MAT) - CONSTANT  - AE1
CE.   FUNMAT(9,MAT) - CONSTANT  - AE0
C
CS.   FUNMAT(5,MAT) - KONSTANTA - AEM
CS.   FUNMAT(6,MAT) - KONSTANTA - ALAM
CS.   FUNMAT(7,MAT) - KONSTANTA - AKA
CS.   FUNMAT(8,MAT) - KONSTANTA - AE1
CS.   FUNMAT(9,MAT) - KONSTANTA - AE0
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OBNOVA/ IRUSI
      COMMON /SRPSKI/ ISRPS
      COMMON /MATER9/ SCMS(3,100),S1S3(2,100,10),KRIT(100)
      common /kadamb/ kadmb(100)
      DIMENSION FUNMAT(9,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD09'
      IF(MAT.GT.100) STOP 'UMOD09 MAT>100'
c
c     klasican restart =0, sa brisanjem =1
c      IRUSI=1
c      IRUSI=0
c
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3),IEL,IKRIT,(SCMS(I,MAT),I=1,3) 
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010)(FUNMAT(I,MAT),I=1,3),IEL,IKRIT,(SCMS(I,MAT),I=1,3)
      FUNMAT(4,MAT)=IEL
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=5,9),kadmb(mat)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1030) (FUNMAT(I,MAT),I=5,9),kadmb(mat)
c     racunanje rastojanja od kriticne povrsi preko s1s3 a ne scms     
      IF(IKRIT.EQ.1) THEN
         CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) KRIT(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1020) KRIT(MAT)
         KRIT(MAT)=KRIT(MAT)+1
         DO 10 K=2,KRIT(MAT)
            IF(K.GT.10.OR.MAT.GT.100) STOP 'UMOD09 NTAC>10.OR.MAT>100'
            CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (S1S3(J,MAT,K),J=1,2) 
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (S1S3(J,MAT,K),J=1,2)
   10    CONTINUE
         S1S3(1,MAT,1)=-FUNMAT(3,MAT)
         S1S3(2,MAT,1)=0.
      ENDIF
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         AMAT(18)=2.
         AMAT(19)=3.
         ISUMGR=MAT
C
         AMAT(25)=FUNMAT(5,MAT)
         AMAT(26)=FUNMAT(6,MAT)
         AMAT(27)=FUNMAT(7,MAT)
         AMAT(28)=FUNMAT(8,MAT)
         AMAT(29)=FUNMAT(9,MAT)
         AMAT(30)=FUNMAT(3,MAT)
C
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
      ENDIF
C      
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,3),IEL,IKRIT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,3),IEL,IKRIT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=5,9),kadmb(mat)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=5,9),kadmb(mat)
      IF(SCMS(3,MAT).GT.1.D-10) THEN
         WRITE(IZLAZ,2046) (SCMS(J,MAT),J=1,3)
      ENDIF
      IF(IKRIT.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) KRIT(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) KRIT(MAT)
         DO 20 K=1,KRIT(MAT)
         WRITE(IZLAZ,1000) (S1S3(J,MAT,K),J=1,2)
   20    CONTINUE
      ENDIF
      RETURN
 
 1000 FORMAT(7F10.0)
 1030 FORMAT(5F10.0,I5)
 1010 FORMAT(3F10.0,2I5,3F10.0)
 1020 FORMAT(14I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,    
     1'MODEL MATERIJALA BROJ =     9  (TLO PLASTICNOST CAM-CLAY (GLINA))
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//
     216X,'G R A N I C A   Z A T E Z A NJ A ...... AT =',1PD12.5//
     216X,'IDIKATOR ZADAVANJA MODULA ELASTIC. ....IEL =',I12//
     216X,'IDIKATOR ZADAVANJA KRIVE LOMA .......IKRIT =',I12//)
 2045 FORMAT(//6X,
     2/6X,'AEM',8X,'ALAM',9X,'AKA',8X,'P0  ',8X,'AEE0',4x,'kadmb'/
     3' ',5(1PD12.4),i9)
 2046 FORMAT(//6X,' m ',8X,' s  ',9X,'Sc '/' ',3(1PD12.4))
 2020 FORMAT(//
     1' BROJ TACAKA ZA LOM ................... KRIT =',I5//
	1' SREDNJI NAPON   SQRT(3*J2D)')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     9  (SOIL ISOTROPIC PLASTICITY (CAM-CL
     2AY))'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//
     216X,'T E N S I O N      C U T O F F  ....... AT =',1PD12.5//
     216X,'I N D I C A T O R  OF  C H A N G E  E .IEL =',I12//
     216X,'IDICATOR OF CRUSH CURVE .............IKRIT =',I12//)
 6045 FORMAT(//6X,
     26X,'AEM',8X,'ALAM',9X,'AKA',8X,'AEE1',8X,'AEE0',4x,'kadmb'/
     3' ',5(1PD12.4),i9)
 6020 FORMAT(//
     1' NUMBER OF POINTS FOR CRUSH CURVE ....... KRIT =',I5//
	1' MEAN STRESS    SQRT(3*J2D)')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD10(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 10
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 10
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  EX
CE.   FUNMAT(2,MAT) - YOUNG*S MODULUS IN DIRECTION Y  -  EY
CE.   FUNMAT(3,MAT) - YOUNG*S MODULUS IN DIRECTION Z  -  EZ
C .
CE.   FUNMAT(4,MAT) - POISSON*S RATIO FOR PLANE XY    -  VXY
CE.   FUNMAT(5,MAT) - POISSON*S RATIO FOR PLANE YZ    -  VYZ
CE.   FUNMAT(6,MAT) - POISSON*S RATIO FOR PLANE ZX    -  VZX
C .
CE.   FUNMAT(7,MAT) - SHEAR MODULUS FOR PLANE XY  -  GXY
CE.   FUNMAT(8,MAT) - SHEAR MODULUS FOR PLANE YZ  -  GYZ
CE.   FUNMAT(9,MAT) - SHEAR MODULUS FOR PLANE ZX  -  GZX
C .
CE.   FUNMAT(10-15,MAT) - YIELD STRESSES X,Y,Z, XY, YZ, ZX
CE.   FUNMAT(16,MAT)    - EFFECTIVE YIELD STRESSES TAUY
CE.   FUNMAT(17,18,MAT) - RAMBERG-OSGOOD CONSTANTS CY,AN
CE.   FUNMAT(19,MAT)    - MIXED HARDENING PARAMETER  M
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(19,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD10'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=10,15)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=10,15)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=17,19)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=17,19)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      FUNMAT(16,MAT)=DSQRT(0.5*((FUNMAT(10,MAT)*FUNMAT(10,MAT)+
     1FUNMAT(11,MAT)*FUNMAT(11,MAT)+FUNMAT(12,MAT)*FUNMAT(12,MAT))/3.+
     1FUNMAT(13,MAT)*FUNMAT(13,MAT)+FUNMAT(14,MAT)*FUNMAT(14,MAT)+
     1FUNMAT(15,MAT)*FUNMAT(15,MAT)))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,19)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,19)
      DO 10 I=10,15
        FUNMAT(I,MAT)=0.5/(FUNMAT(I,MAT)*FUNMAT(I,MAT))
   10 CONTINUE
      F=FUNMAT(11,MAT)+FUNMAT(12,MAT)-FUNMAT(10,MAT)
      G=FUNMAT(12,MAT)+FUNMAT(10,MAT)-FUNMAT(11,MAT)
      H=FUNMAT(10,MAT)+FUNMAT(11,MAT)-FUNMAT(12,MAT)
      FUNMAT(10,MAT)=H
      FUNMAT(11,MAT)=G
      FUNMAT(12,MAT)=F
C... FORM N1 TO N6
      SIG2=FUNMAT(16,MAT)
      SIG2=SIG2*SIG2*2./3.
      DO 20 I=10,15
   20 FUNMAT(I,MAT)=FUNMAT(I,MAT)*SIG2
C... CHECK MATERIAL CONSTANTS
      CALL       ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    10  (ELASTO-PLASTICAN ANIZOTROPAN)'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1//16X,'N A P O N I    T E C E N J A'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'N A P O N    T E C E N J A (EFEKTIV.) TAUY =',1PD12.5//
     216X,'M N O Z I O C  .....................    CY =',1PD12.5//
     316X,'E K S P O N E N T  .................    AN =',1PD12.5//
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    10  (ELASTIC-PLASTIC ANISOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U S'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1//16X,'Y I E L D    S T R E S S E S'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     116X,'Y I E L D    S T R E S S (EFFECTIV).  TAUY =',1PD12.5//
     216X,'M U L T I P L A Y E R  .............    CY =',1PD12.5//
     316X,'E X P O N E N T  ...................    AN =',1PD12.5//
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE ANICHK(F,IZLAZ,ISRPS)
C
CE      CHECK ANISOTROPY INPUT DATA
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(*)
      IND=2
C
      IF((F(1)/F(2)-F(4)*F(4)).GT.0.D0 .AND.
     &   (F(2)/F(3)-F(5)*F(5)).GT.0.D0 .AND.
     &   (F(3)/F(1)-F(6)*F(6)).GT.0.D0) IND=IND-1
C
      CHEK=(1.-F(4)*F(4)*F(2)/F(1)-F(5)*F(5)*F(3)/F(2)-
     &         F(6)*F(6)*F(1)/F(3))*0.5
      IF((F(4)*F(5)*F(6) .LT. CHEK) .AND.
     &   (CHEK .LT. 0.5D0)) IND=IND-1
C
      IF(IND.GT.0)THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2000) IND
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6000) IND
        STOP
      ENDIF
      RETURN
 2000 FORMAT(/6X,'MATERIJALNE KONSTANTE NISU FIZICKI KONZISTENTNE ZA ',
     1 I3,' TESTA')
 6000 FORMAT(/6X,'MATERIAL CONSTANTS ARE NOT PHISICALY CONSISTENT IN ',
     1 I3,' TESTS')
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD11(FUNMAT,NTFUN,MATE,MAXT,MAT,KARTI,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 11 AND 12
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 11 I 12
C .
CE.               MAT - MATERIAL NUMBER
CE.   J=1,(NTMAX - NUMBER OF POINTS FOR CURVE)
CS.               MAT - MATERIJAL BROJ
CS.   J=1,(NTMAX - BROJ TACAKA ZA KRIVU)
C .
CE. A) MATERIAL MODEL NUMBER 11 "GAP - STRESS" (IND=1)
CS. A) MODEL MATERIJALA BROJ 11 "ZAZOR SA ZADATIM NAPONOM" (IND=1)
C .
CE.    NTMAX>=1 (NTMAX - NUMBER OF POINTS FOR CURVE TAU(DEF))
CS.    NTMAX>=1 (NTMAX - BROJ TACAKA ZA KRIVU TAU(DEF))
C .
CE.   FUNMAT(1,MAT,J) - ARGUMENT - STRAIN - DEF
CE.   FUNMAT(2,MAT,J) - STRESS - TAU
CS.   FUNMAT(1,MAT,J) - ARGUMENT - DEFORMACIJA - DEF
CS.   FUNMAT(2,MAT,J) - NAPON - TAU
C .
CE.        NTFUN(MAT) - NUMBER OF POINTS FOR CURVE - TAU(DEF)
CE.                     =1; LINEAR CURVE
CE.                         (ASSUMPTION: FIRST POINT IS (0.,0.))
CE.                     >1; MULTILINEAR CURVE
CS.        NTFUN(MAT) - BROJ TACAKA ZA KRIVU - TAU(DEF)
CS.                     =1; LINEARNA KRIVA
CS.                         (PRETPOSTAVKA: PRVA TACKA JE (0.,0.))
CS.                     >1; MULTILINEARNA KRIVA
C .
CE. B) MATERIAL MODEL NUMBER 12 (GAP - FORCES) (IND=2)
CS. B) MODEL MATERIJALA BROJ 12 (ZAZOR SA ZADATOM SILOM) (IND=2)
C .
CE.    NTMAX>=1 (NTMAX - NUMBER OF POINTS FOR CURVE FRC(EXT))
CS.    NTMAX>=1 (NTMAX - BROJ TACAKA ZA KRIVU FRC(EXT))
C .
CE.   FUNMAT(1,MAT,J) - ARGUMENT - EXTENSION - EXT
CE.   FUNMAT(2,MAT,J) - FORCE - FRC
CS.   FUNMAT(1,MAT,J) - ARGUMENT - IZDUZENJE - EXT
CS.   FUNMAT(2,MAT,J) - SILA - FRC
C .
CE.        NTFUN(MAT) - NUMBER OF POINTS FOR CURVE - FRC(EXT)
CE.                     =1; LINEAR CURVE
CE.                         (ASSUMPTION: FIRST POINT IS (0.,0.))
CE.                     >1; MULTILINEAR CURVE
CS.        NTFUN(MAT) - BROJ TACAKA ZA KRIVU - FRC(EXT)
CS.                     =1; LINEARNA KRIVA
CS.                         (PRETPOSTAVKA: PRVA TACKA JE (0.,0.))
CS.                     >1; MULTILINEARNA KRIVA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(2,MATE,*),NTFUN(*)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD11'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NTFUN(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NTFUN(MAT)
        IF(NTFUN(MAT).GT.MAXT.OR.NTFUN(MAT).LE.0) STOP 'UMOD11'
      DO 20 J=1,NTFUN(MAT)
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (FUNMAT(I,MAT,J),I=1,2)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1010) (FUNMAT(I,MAT,J),I=1,2)
   20 CONTINUE
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(IND.EQ.1) THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2000) MAT
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6000) MAT
      ELSE
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2010) MAT
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6010) MAT
      END IF
C
      IF(IND.EQ.1) THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2030) NTFUN(MAT)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6030) NTFUN(MAT)
      ELSE
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2040) NTFUN(MAT)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6040) NTFUN(MAT)
      END IF
      WRITE(IZLAZ,5020) ((FUNMAT(I,MAT,J),I=1,2),J=1,NTFUN(MAT))
C
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    11'/6X,'(ZAZOR SA ZADATIM NAPONOM',
     2')'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    12'/6X,'(ZAZOR SA ZADATOM SILOM',
     2')'///11X,'MATERIJAL BROJ =',I5)
 2030 FORMAT(//
     116X,'K R I V A    N A P O N  -  D E F O R M A C I J A'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT             FUNKCIJA'/
     120X,'DEFORMACIJA             NAPON')
 2040 FORMAT(//
     116X,'K R I V A    S I L A  -  I Z D U Z E NJ E'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'ARGUMENT             FUNKCIJA'/
     120X,'IZDUZENJE                SILA')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =   11'/6X,'(GAP - STRESS)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =   12'/6X,'(GAP - FORCE)'
     2///11X,'MATERIAL NUMBER =',I5)
 6030 FORMAT(//
     116X,'S T R E S S - S T R A I N    C U R V E'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT             FUNCTION'/
     120X,'STRAIN                 STRESS')
 6040 FORMAT(//
     116X,'F O R C E - E X T E N S I O N    C U R V E'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'ARGUMENT             FUNCTION'/
     120X,'EXTENSION               FORCE')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE NAPPLD(FUNMAT,NTFUN,MATE,MAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO TRANSFORMATION OF CURVE ON FORM STRESS-PLASTIC STRAIN
CS.    P R O G R A M
CS.        ZA PREVODJENJE KRIVE NA OBLIK NAPON-PLASTICNA DEFORMACIJA
C .
C ......................................................................
C
      DIMENSION FUNMAT(2,MATE,*),NTFUN(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' NAPPLD'
      FUNMAT(1,MAT,2)=0.0D0
      DO 10 J=3,NTFUN(MAT)+1
   10 FUNMAT(1,MAT,J)=FUNMAT(1,MAT,J)-FUNMAT(2,MAT,J)/FUNMAT(1,MAT,1)
      RETURN
      END
C=======================================================================
      SUBROUTINE UMOD1U(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 10
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 10
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1-3,MAT)   - YOUNG*S MODULUS EX,EY,EZ
CE.   FUNMAT(4-6,MAT)   - POISSON*S RATIO  VXY, VYZ, VZX
CE.   FUNMAT(7-9,MAT)   - SHEAR MODULUS GXY, GYZ, GZX
CE.   FUNMAT(10-15,MAT) - YIELD STRESSES X,Y,Z, XY, YZ, ZX
CE.   FUNMAT(16-21,MAT) - RAMBERG-OSGOOD CONSTANTS (C)  X,Y,Z,XY,YZ,ZX
CE.   FUNMAT(22-27,MAT) - RAMBERG-OSGOOD CONSTANTS (AN) X,Y,Z,XY,YZ,ZX
CE.   FUNMAT(28,MAT)    - EFFECTIVE YIELD STRESSES TAUY
CE.   FUNMAT(29,MAT)    - MIXED HARDENING PARAMETER  M
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(29,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD10'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
C... NAPONI TECENJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=10,15)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=10,15)
C... RAMBERG OSGOOD MNOZIOCI
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=16,21)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=16,21)
      DO 10 I=16,21
         IF(FUNMAT(I,MAT).LT.1.D-8) FUNMAT(I,MAT)=1.D-8
   10 CONTINUE
C... RAMBERG OSGOOD EKSPONENTI
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=22,27)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=22,27)
C... KOEFICIJENT MESOVITOG OJACANJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=29,29)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=29,29)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      FUNMAT(28,MAT)=FUNAR(FUNMAT(10,MAT))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,29)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,29)
C... CHECK MATERIAL CONSTANTS
      CALL       ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    10  (ELASTO-PLASTICAN ANIZOTROPAN)'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     1//16X,'N A P O N I    T E C E N J A'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     1//16X,'M N O Z I O C I'             //11X,'CX ',9X,'CY ',
     19X,'CZ ',9X,'CXY',9X,'CYZ',9X,'CZX'/7X,6(1PD11.4,1X)//
     1//16X,'E K S P O N E N T I'         //11X,'AX ',9X,'AY ',
     19X,'AZ ',9X,'AXY',9X,'AYZ',9X,'AZX'/7X,6(1PD11.4,1X)//
     116X,'N A P O N    T E C E N J A (EFEKTIV.) TAUY =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5) 
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    10  (ELASTIC-PLASTIC ANISOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U S'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     1//16X,'Y I E L D    S T R E S S E S'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     116X,'M U L T I P L A Y E R S'       //11X,'CX ',9X,'CY ',
     19X,'CZ ',9X,'CXY',9X,'CYZ',9X,'CZX'/7X,6(1PD11.4,1X)//
     1//16X,'E X P O N E N T S'           //11X,'AX ',9X,'AY ',
     19X,'AZ ',9X,'AXY',9X,'AYZ',9X,'AZX'/7X,6(1PD11.4,1X)//
     116X,'Y I E L D    S T R E S S (EFFECTIVE)  TAUY =',1PD12.5//
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE UMOD13(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 13
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 13
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1-3,MAT)   - YOUNG*S MODULUS EX,EY,EZ
CE.   FUNMAT(4-6,MAT)   - POISSON*S RATIO  VXY, VYZ, VZX
CE.   FUNMAT(7-9,MAT)   - SHEAR MODULUS GXY, GYZ, GZX
CE.   FUNMAT(10,MAT)    - EFFECTIVE YIELD STRESSES TAUY
CE.   FUNMAT(11,MAT)    - MIXED HARDENING PARAMETER  M
CE.   FUNMAT(12-74,MAT) - 21*3 RAMBERG-OSGOOD CONSTANTS (Y0,C,AN) -> N MATRIX
CE.   FUNMAT(75-77,MAT) - 3 RAMBERG-OSGOOD CONSTANTS BACK STRESS RATE (C BAR)
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(77,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD13'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
C... RAMBERG OSGOOD KONSTANTE
      I1=9
      DO 10 I=1,22
      CALL ISPITA(ACOZ)
      I1=I1+3
      I2=I1+2
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(J,MAT),J=I1,I2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(J,MAT),J=I1,I2)
   10 CONTINUE
C... KOEFICIJENT MESOVITOG OJACANJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(11,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(11,MAT)
C
CE  CALCULATE INITIAL EFFECTIVE YIELD STRESS
C
      DUM=0.D0
      I1=-9
C      DO 20 I=1,6
C        I1=I1+24-3*I
C        DUM=DUM+1./FUNMAT(I1,MAT)
C   20 CONTINUE
C      FUNMAT(10,MAT)=DSQRT(0.25*DUM)
C      FUNMAT(10,MAT)=FUNAR(EN)
      DO 20 I=1,6
        I1=I1+24-3*I
        DUM=DUM+FUNMAT(I1,MAT)
   20 CONTINUE
      FUNMAT(10,MAT)=DUM/6.
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,9),
     1               (JJ,(FUNMAT(J,MAT),J=12+3*(JJ-1),11+3*JJ),JJ=1,21),
     2                  FUNMAT(10,MAT),FUNMAT(11,MAT),
     3                  (FUNMAT(J,MAT),J=75,77)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,9),
     1               (JJ,(FUNMAT(J,MAT),J=12+3*(JJ-1),11+3*JJ),JJ=1,21),
     2                  FUNMAT(10,MAT),FUNMAT(11,MAT),
     3                  (FUNMAT(J,MAT),J=75,77)
C... CHECK MATERIAL CONSTANTS
      CALL       ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    13  (GENERALNI PLASTICAN ANIZOTROPAN)'
     2///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'R.B.',8X,'Y0',10X,'CY',10X,'AN'/
     1  21(/12X,I5,2X,3(1PD11.4,1X))//
     116X,'N A P O N    T E C E N J A (EFEKTIV.) TAUY =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5//
     1'   RAMBERG-OSGOOD KONSTANTE ZA KINEMATSKO OJACANJE'
     1//28X,'A0',10X,'CY',10X,'AN'//19X,3(1PD11.4,1X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    13  (GENERAL-PLASTIC ANISOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U S'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     116X,'R.B.',8X,'Y0',10X,'CY',10X,'AN'/
     1  21(/12X,I5,2X,3(1PD11.4,1X))//
     116X,'Y I E L D    S T R E S S (EFFECTIVE)  TAUY =',1PD12.5//
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5//
     1'   RAMBERG-OSGOOD CONSTANTS FOR BACKE STRESS'
     1//28X,'A0',10X,'CY',10X,'AN'//19X,3(1PD11.4,1X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD14(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 14
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 14
C .
CE.   M=(MAT-1)*4   (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*4   (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - E
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - POISSON*S RATIO - V
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - POISSON - OV BROJ - V
C .
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - COEFFICIENT OF THERMAL EXPANSION - A
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - KOEFICIJENT TERMICKOG SIRENJA - A
C .
C .
CE.       RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS
CS.       RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - YIELD STRESS - TAUY
CE.   FUNMAT(3,M+4,J) - MULTIPLAYER - CY
CE.   FUNMAT(4,M+4,J) - EXPONENT - AN
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - NAPON TECENJA - TAUY
CS.   FUNMAT(3,M+4,J) - MNOZIOC - CY
CS.   FUNMAT(4,M+4,J) - EKSPONENT - AN
C .
C .
CE.      TREF(MAT)    - REFERENCE TEMPERATURE
CS.      TREF(MAT)    - REFERENTNA TEMPERATURA
C .
CE.   FUNMAT(3,M+1,1) - MIXED HARDENING COEFFICIENT
CS.   FUNMAT(3,M+1,1) - KOEFICIJENT MESOVITOG OJACANJA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS ON CURVE - E
CE.      NTFUN(M+2)   - NUMBER OF POINTS ON CURVE - V
CE.      NTFUN(M+3)   - NUMBER OF POINTS ON CURVE - A
CE.      NTFUN(M+4)   - NUMBER OF POINTS ON CURVES - TAUY,CY,AN
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - E
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - V
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - A
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVE - TAUY,CY,AN
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION FUNMAT(4,MATE*4,*),TREF(*),NTFUN(*),AMAT(30)
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD14'
      NFUN=1
      MA=(MAT-1)*4
      DO 10 K=1,4
        IF (K.EQ.4) NFUN=3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD14'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(K.EQ.4) THEN
             IF(FUNMAT(3,MA+K,J).LT.1.D-8) FUNMAT(3,MA+K,J)=1.D-8
          ENDIF
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(3,MA+1,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(3,MA+1,1)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(2,MA+1,1)
         AMAT(7)=FUNMAT(2,MA+2,1)
         AMAT(10)=FUNMAT(2,MA+3,1)
         AMAT(14)=TREF(MAT)
         AMAT(13)=GUST
         AMAT(15)=FUNMAT(2,MA+4,1)
         AMAT(16)=FUNMAT(3,MA+4,1)
         AMAT(17)=FUNMAT(3,MA+1,1)
         AMAT(19)=3.
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MA+1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MA+1)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+1,J),I=1,2),J=1,NTFUN(MA+1))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MA+2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MA+2)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+2,J),I=1,2),J=1,NTFUN(MA+2))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NTFUN(MA+3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NTFUN(MA+3)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+3,J),I=1,2),J=1,NTFUN(MA+3))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060) NTFUN(MA+4)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060) NTFUN(MA+4)
      WRITE(IZLAZ,5030) ((FUNMAT(I,MA+4,J),I=1,4),J=1,NTFUN(MA+4))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) FUNMAT(3,MA+1,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) FUNMAT(3,MA+1,1)
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(9X,1PD12.5,6X,1PD12.5,6X,1PD12.5,6X,1PD12.5)
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =    14  (MIZESOV TERMO-ELASTO-PLASTICAN '/
     &38X,'IZOTROPAN SA MESOVITIM OJACANJEM)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E')
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V')
 2060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A ',
     &   '(MATERIJALNE CONSTANTE)'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',3(10X,'FUNKCIJA')/
     &10X,'TEMPERATURA',10X,'TAUY',15X,'CY',16X,'AN')
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A')
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2070 FORMAT(//7X,
     &'K O E F I C I J E N T    M E S O V I T O G    O J A C A NJ A'//
     &7X,'EM =',1PD12.5)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =    14  (VON MISES THERMO-ELASTIC-PLASTIC
     &)'/38X,'ISOTROPIC WITH MIXED HARDENING)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E')
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V')
 6060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A S ',
     &   '(MATERIAL CONSTANTS)'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',3(10X,'FUNCTION')/
     &10X,'TEMPERATURE',10X,'TAUY',15X,'CY',16X,'AN')
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A')
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6070 FORMAT(//7X,
     &'M I X E D    H A R D E N I N G    C O E F F I C I E N T'//
     &7X,'EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD15(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 15
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 15
C .
CE.   M=(MAT-1)*3   (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*3   (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - E
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - POISSON*S RATIO - V
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - POISSON - OV BROJ - V
C .
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - COEFFICIENT OF THERMAL EXPANSION - A
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - KOEFICIJENT TERMICKOG SIRENJA - A
C .
CE.      TREF(MAT)    - REFERENCE TEMPERATURE
CS.      TREF(MAT)    - REFERENTNA TEMPERATURA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS ON CURVE - E
CE.      NTFUN(M+2)   - NUMBER OF POINTS ON CURVE - V
CE.      NTFUN(M+3)   - NUMBER OF POINTS ON CURVE - A
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - E
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - V
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - A
C .
C .
CE.             CREEP LAW
CS.           ZAKON PUZANJA
C.
CE.      ICLAW - CREEP LAW FORM INDICATOR
CE.            = 1 BAILEY-NORTON CREEP LAW    ( NACR=3 )
CE.            = 2 EXPONENTIAL CREEP LAW      ( NACR=7 )
CE.            = 3 EIGHT PARAMETER CREEP LAW  ( NACR=8 )
CE.            = 4 CREEP LAW - USER DEFINED
CS.      ICLAW - INDIKATOR OBLIKA ZAKONA PUZANJA
CS.            = 1 BAILEY-NORTON-OV ZAKON PUZANJA   ( NACR=3 )
CS.            = 2 EKSPONENCIJALNI ZAKON PUZANJA    ( NACR=7 )
CS.            = 3 ZAKON PUZANJA SA 8 PARAMETARA    ( NACR=8 )
CS.            = 4 KORISNIKOV ZAKON PUZANJA
C.
CE.      ICELAW - EXPONENTIAL CREEP LAW FORM INDICATOR
CE.             = 1 SNYDER & BATHE EXPONENTIAL CREEP LAW
CS.      ICELAW - EKSPONENCIJALNI INDIKATOR OBLIKA ZAKONA PUZANJA
CS.             = 1 SNYDER & BATHE EKSPONENCIJALNI ZAKON PUZANJA
C.
CE.  FOR ICLAW = 1,2,3
CS.  ZA  ICLAW = 1,2,3
C.
CE.  ACR(I),I=1,NACR - MATERIAL CONSTANTS, TEMPERATURE INDEPENDENT
CS.  ACR(I),I=1,NACR - MATERIJALNE KONSTANTE, NEZAVISNE OD TEMPERATURE
C.
CE.  FOR ICLAW = 4
CS.  ZA  ICLAW = 4
C.
CE.           USTR1  - CREEP LAW (STRAIN)
CS.           USTR1  - ZAKON PUZANJA (DEFORMACIJA)
C.
CE.           USTR2  - CREEP LAW (STRAIN RATE)
CS.           USTR2  - ZAKON PUZANJA (BRZINA DEFORMACIJA)
C ......................................................................
C
      CHARACTER*250 ACOZ,USTR1,USTR2,IZRAZ,RP1,RP2
      CHARACTER*17 FRMT
      CHARACTER*41 FMTS,FMTE
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      DIMENSION FUNMAT(2,MATE*3,*),TREF(*),NTFUN(*),NACR(3)
      DATA NACR/3,7,7/
      DATA FRMT/'(5X,7(8X,1HA,I1))'/
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD15'
      MA=(MAT-1)*3
      DO 10 K=1,3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD15'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,2)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,2)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) ICLAW,ICELAW
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) ICLAW,ICELAW
C
      IF (ICLAW.NE.4) THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,*) (ACR(I),I=1,NACR(ICLAW))
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,1010) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) THEN
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     &    READ(IULAZ,*) ACR(8)
          IF(INDFOR.EQ.2)
     &    READ(ACOZ,1010) ACR(8)
        END IF
      ELSE
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR1
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR1
        DO 29 I=1,80
          RP1(I:I)=' '
   29     IZRAZ(I:I)=' '
        CALL CITAC(USTR1,IZRAZ,ACR,ISG1,ITI1,ITH1)
        RP1(1:1)='$'
        CALL KONVEB(IZRAZ,RP1)
        DO 30 I=1,80
          RP2(I:I)=' '
   30     IZRAZ(I:I)=' '
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR2
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR2
        CALL CITAC(USTR2,IZRAZ,BCR,ISG2,ITI2,ITH2)
        RP2(1:1)='$'
        CALL KONVEB(IZRAZ,RP2)
      END IF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MA+1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MA+1)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+1,J),I=1,2),J=1,NTFUN(MA+1))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MA+2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MA+2)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+2,J),I=1,2),J=1,NTFUN(MA+2))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NTFUN(MA+3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NTFUN(MA+3)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+3,J),I=1,2),J=1,NTFUN(MA+3))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      IF(ICLAW.EQ.1) THEN
        FMTS='(20X,30HBAILEY-NORTON-OV ZAKON PUZANJA//)'
        FMTE='(20X,30HBAILEY-NORTON CREEP LAW       //)'
      ELSE IF(ICLAW.EQ.2) THEN
        FMTS='(20X,30HEKSPONENCIJALNI ZAKON PUZANJA //)'
        FMTE='(20X,30HEKSPONENTIAL CREEP LAW        //)'
      ELSE IF(ICLAW.EQ.3) THEN
        FMTS='(20X,30HZAKON PUZANJA SA 8 PARAMETARA //)'
        FMTE='(20X,30HEIGHT PARAMETER CREEP LAW     //)'
      ELSE
        FMTS='(20X,30HKORISNIKOV ZAKON PUZANJA      //)'
        FMTE='(20X,30HCREEP LAW - USER DEFINED      //)'
      END IF
C
      IF (ICLAW.NE.4) THEN
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,2080)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,6080)
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,FMTS)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,FMTE)
        IF (ICLAW.EQ.1) FRMT(5:5)='3'
        WRITE(IZLAZ,FRMT)(I,I=0,NACR(ICLAW)-1)
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(43X,''A7'')')
        WRITE(IZLAZ,5030) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(35X,1PD10.3)') ACR(8)
      ELSE
        WRITE(IZLAZ,'(10X,''S T R A I N''//A80)')USTR1
        WRITE(IZLAZ,'(/10X,''S T R A I N    R A T E''//A80)')USTR2
      END IF
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(5X,7(1PD10.3))
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =     15  (TERMO-ELASTICAN SA PUZANJEM)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E')
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V')
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A')
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2080 FORMAT(//
     &7X,'Z A K O N    P U Z A NJ A   (MATERIJALNE CONSTANTE)'//)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =    15 (THERMO-ELASTIC WITH CREEP)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E')
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V')
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A')
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6080 FORMAT(//7X,'C R E E P    L A W    (MATERIAL CONSTANTS)'//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD16(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 16
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 16
C .
CE.   M=(MAT-1)*4   (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*4   (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - E
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - POISSON*S RATIO - V
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - POISSON - OV BROJ - V
C .
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - COEFFICIENT OF THERMAL EXPANSION - A
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - KOEFICIJENT TERMICKOG SIRENJA - A
C .
C .
CE.       RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS
CS.       RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - YIELD STRESS - TAUY
CE.   FUNMAT(3,M+4,J) - MULTIPLAYER - CY
CE.   FUNMAT(4,M+4,J) - EXPONENT - AN
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - NAPON TECENJA - TAUY
CS.   FUNMAT(3,M+4,J) - MNOZIOC - CY
CS.   FUNMAT(4,M+4,J) - EKSPONENT - AN
C .
C .
CE.      TREF(MAT)    - REFERENCE TEMPERATURE
CS.      TREF(MAT)    - REFERENTNA TEMPERATURA
C .
CE.   FUNMAT(3,M+1,1) - MIXED HARDENING COEFFICIENT
CS.   FUNMAT(3,M+1,1) - KOEFICIJENT MESOVITOG OJACANJA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS ON CURVE - E
CE.      NTFUN(M+2)   - NUMBER OF POINTS ON CURVE - V
CE.      NTFUN(M+3)   - NUMBER OF POINTS ON CURVE - A
CE.      NTFUN(M+4)   - NUMBER OF POINTS ON CURVES - TAUY,CY,AN
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - E
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - V
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - A
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVE - TAUY,CY,AN
C .
C .
CE.             CREEP LAW
CS.           ZAKON PUZANJA
C.
CE.      ICLAW - CREEP LAW FORM INDICATOR
CE.            = 1 BAILEY-NORTON CREEP LAW    ( NACR=3 )
CE.            = 2 EXPONENTIAL CREEP LAW      ( NACR=7 )
CE.            = 3 EIGHT PARAMETER CREEP LAW  ( NACR=8 )
CE.            = 4 CREEP LAW - USER DEFINED
CS.      ICLAW - INDIKATOR OBLIKA ZAKONA PUZANJA
CS.            = 1 BAILEY-NORTON-OV ZAKON PUZANJA   ( NACR=3 )
CS.            = 2 EKSPONENCIJALNI ZAKON PUZANJA    ( NACR=7 )
CS.            = 3 ZAKON PUZANJA SA 8 PARAMETARA    ( NACR=8 )
CS.            = 4 KORISNIKOV ZAKON PUZANJA
C.
CE.      ICELAW - EXPONENTIAL CREEP LAW FORM INDICATOR
CE.             = 1 SNYDER & BATHE EXPONENTIAL CREEP LAW
CS.      ICELAW - EKSPONENCIJALNI INDIKATOR OBLIKA ZAKONA PUZANJA
CS.             = 1 SNYDER & BATHE EKSPONENCIJALNI ZAKON PUZANJA
C.
CE.  FOR ICLAW = 1,2,3
CS.  ZA  ICLAW = 1,2,3
C.
CE.  ACR(I),I=1,NACR - MATERIAL CONSTANTS, TEMPERATURE INDEPENDENT
CS.  ACR(I),I=1,NACR - MATERIJALNE KONSTANTE, NEZAVISNE OD TEMPERATURE
C.
CE.  FOR ICLAW = 4
CS.  ZA  ICLAW = 4
C.
CE.           USTR1  - CREEP LAW (STRAIN)
CS.           USTR1  - ZAKON PUZANJA (DEFORMACIJA)
C.
CE.           USTR2  - CREEP LAW (STRAIN RATE)
CS.           USTR2  - ZAKON PUZANJA (BRZINA DEFORMACIJA)
C ......................................................................
C
      CHARACTER*250 ACOZ,USTR1,USTR2,IZRAZ,RP1,RP2
      CHARACTER*17 FRMT
      CHARACTER*41 FMTS,FMTE
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      DIMENSION FUNMAT(4,MATE*4,*),TREF(*),NTFUN(*),NACR(3)
      DATA NACR/3,7,7/
      DATA FRMT/'(5X,7(8X,1HA,I1))'/
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD16'
      NFUN=1
      MA=(MAT-1)*4
      DO 10 K=1,4
        IF (K.EQ.4) NFUN=3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD16'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(K.EQ.4) THEN
            IF(FUNMAT(3,MA+K,J).LT.1.D-8) FUNMAT(3,MA+K,J)=1.D-8
          ENDIF
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(3,MA+1,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(3,MA+1,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) ICLAW,ICELAW
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) ICLAW,ICELAW
C
      IF (ICLAW.NE.4) THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,*) (ACR(I),I=1,NACR(ICLAW))
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,1010) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) THEN
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     &    READ(IULAZ,*) ACR(8)
          IF(INDFOR.EQ.2)
     &    READ(ACOZ,1010) ACR(8)
        END IF
      ELSE
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR1
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR1
        DO 29 I=1,80
          RP1(I:I)=' '
   29     IZRAZ(I:I)=' '
        CALL CITAC(USTR1,IZRAZ,ACR,ISG1,ITI1,ITH1)
        RP1(1:1)='$'
        CALL KONVEB(IZRAZ,RP1)
        DO 30 I=1,80
          RP2(I:I)=' '
   30     IZRAZ(I:I)=' '
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR2
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR2
        CALL CITAC(USTR2,IZRAZ,BCR,ISG2,ITI2,ITH2)
        RP2(1:1)='$'
        CALL KONVEB(IZRAZ,RP2)
      END IF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MA+1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MA+1)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+1,J),I=1,2),J=1,NTFUN(MA+1))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MA+2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MA+2)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+2,J),I=1,2),J=1,NTFUN(MA+2))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NTFUN(MA+3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NTFUN(MA+3)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MA+3,J),I=1,2),J=1,NTFUN(MA+3))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060) NTFUN(MA+4)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060) NTFUN(MA+4)
      WRITE(IZLAZ,5025) ((FUNMAT(I,MA+4,J),I=1,4),J=1,NTFUN(MA+4))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) FUNMAT(3,MA+1,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) FUNMAT(3,MA+1,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      IF(ICLAW.EQ.1) THEN
        FMTS='(20X,30HBAILEY-NORTON-OV ZAKON PUZANJA//)'
        FMTE='(20X,30HBAILEY-NORTON CREEP LAW       //)'
      ELSE IF(ICLAW.EQ.2) THEN
        FMTS='(20X,30HEKSPONENCIJALNI ZAKON PUZANJA //)'
        FMTE='(20X,30HEKSPONENTIAL CREEP LAW        //)'
      ELSE IF(ICLAW.EQ.3) THEN
        FMTS='(20X,30HZAKON PUZANJA SA 8 PARAMETARA //)'
        FMTE='(20X,30HEIGHT PARAMETER CREEP LAW     //)'
      ELSE
        FMTS='(20X,30HKORISNIKOV ZAKON PUZANJA      //)'
        FMTE='(20X,30HCREEP LAW - USER DEFINED      //)'
      END IF
C
      IF (ICLAW.NE.4) THEN
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,2080)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,6080)
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,FMTS)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,FMTE)
        IF (ICLAW.EQ.1) FRMT(5:5)='3'
        WRITE(IZLAZ,FRMT)(I,I=0,NACR(ICLAW)-1)
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(43X,''A7'')')
        WRITE(IZLAZ,5030) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(35X,1PD10.3)') ACR(8)
      ELSE
        WRITE(IZLAZ,'(10X,''S T R A I N''//A80)')USTR1
        WRITE(IZLAZ,'(/10X,''S T R A I N    R A T E''//A80)')USTR2
      END IF
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5025 FORMAT(9X,1PD12.5,6X,1PD12.5,6X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(5X,7(1PD10.3))
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =    16  (MIZESOV TERMO-ELASTO-PLASTICAN '/
     &24X,'SA PUZANJEM - IZOTROPAN SA MESOVITIM OJACANJEM)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E')
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V')
 2060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A ',
     &   '(MATERIJALNE CONSTANTE)'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',3(10X,'FUNKCIJA')/
     &10X,'TEMPERATURA',10X,'TAUY',15X,'CY',16X,'AN')
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A')
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2070 FORMAT(//7X,
     &'K O E F I C I J E N T    M E S O V I T O G    O J A C A NJ A'//
     &7X,'EM =',1PD12.5)
 2080 FORMAT(//
     &7X,'Z A K O N    P U Z A NJ A   (MATERIJALNE CONSTANTE)'//)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =   16  (VON MISES THERMO-ELASTIC-PLASTI',
     &'C'/24X,'WITH CREEP - ISOTROPIC WITH MIXED HARDENING)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E')
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V')
 6060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A S ',
     &   '(MATERIAL CONSTANTS)'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',3(10X,'FUNCTION')/
     &10X,'TEMPERATURE',10X,'TAUY',15X,'CY',16X,'AN')
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A')
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6070 FORMAT(//7X,
     &'M I X E D    H A R D E N I N G    C O E F F I C I E N T'//
     &7X,'EM =',1PD12.5)
 6080 FORMAT(//7X,'C R E E P    L A W    (MATERIAL CONSTANTS)'//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD17(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 17
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 17
C .
CE.   M=(MAT-1)*18  (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*18  (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - EX
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - YOUNG*S MODULUS - EY
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - YOUNG*S MODULUS - EZ
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - EX
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - MODUL ELASTICNOSTI - EY
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - MODUL ELASTICNOSTI - EZ
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - POISSON*S RATIO - VXY
CE.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+5,J) - POISSON*S RATIO - VYZ
CE.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+6,J) - POISSON*S RATIO - VZX
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - POISSONOV BROJ - VXY
CS.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+5,J) - POISSONOV BROJ - VYZ
CS.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+6,J) - POISSONOV BROJ - VZX
C .
CE.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+7,J) - SHEAR MODULUS - GXY
CE.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+8,J) - SHEAR MODULUS - GYZ
CE.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+9,J) - SHEAR MODULUS - GZX
CS.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+7,J) - MODUL SMICANJA - GXY
CS.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+8,J) - MODUL SMICANJA - GYZ
CS.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+9,J) - MODUL SMICANJA - GZX
C .
CE.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+10,J) - COEFFICIENT OF THERMAL EXPANSION - AX
CE.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+11,J) - COEFFICIENT OF THERMAL EXPANSION - AY
CE.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+12,J) - COEFFICIENT OF THERMAL EXPANSION  - AZ
CS.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+10,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AX
CS.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+11,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AY
CS.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+12,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AZ
C .
CE.       RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (X - DIRECTION)
CS.       RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (X - PRAVAC)
C .
CE.  FUNMAT(1,M+13,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+13,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+13,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+13,J) - EXPONENT - AN
CS.  FUNMAT(1,M+13,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+13,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+13,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+13,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (Y - DIRECTION)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (Y - PRAVAC)
C .
CE.  FUNMAT(1,M+14,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+14,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+14,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+14,J) - EXPONENT - AN
CS.  FUNMAT(1,M+14,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+14,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+14,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+14,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (Z - DIRECTION)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (Z - PRAVAC)
C .
CE.  FUNMAT(1,M+15,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+15,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+15,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+15,J) - EXPONENT - AN
CS.  FUNMAT(1,M+15,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+15,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+15,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+15,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (XY - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (XY - RAVAN)
C .
CE.  FUNMAT(1,M+16,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+16,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+16,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+16,J) - EXPONENT - AN
CS.  FUNMAT(1,M+16,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+16,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+16,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+16,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (YZ - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (YZ - RAVAN)
C .
CE.  FUNMAT(1,M+17,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+17,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+17,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+17,J) - EXPONENT - AN
CS.  FUNMAT(1,M+17,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+17,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+17,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+17,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (ZX - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (ZX - RAVAN)
C .
CE.  FUNMAT(1,M+18,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+18,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+18,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+18,J) - EXPONENT - AN
CS.  FUNMAT(1,M+18,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+18,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+18,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+18,J) - EKSPONENT - AN
C .
CE.         TREF(MAT) - REFERENCE TEMPERATURE
CS.         TREF(MAT) - REFERENTNA TEMPERATURA
C .
CE.   FUNMAT(3,M+1,1) - MIXED HARDENING COEFFICIENT
CS.   FUNMAT(3,M+1,1) - KOEFICIJENT MESOVITOG OJACANJA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS - EX
CE.      NTFUN(M+2)   - NUMBER OF POINTS - EY
CE.      NTFUN(M+3)   - NUMBER OF POINTS - EZ
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - EX
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - EY
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - EZ
C .
CE.      NTFUN(M+4)   - NUMBER OF POINTS - VXY
CE.      NTFUN(M+5)   - NUMBER OF POINTS - VYZ
CE.      NTFUN(M+6)   - NUMBER OF POINTS - VZX
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVU - VXY
CS.      NTFUN(M+5)   - BROJ TACAKA ZA KRIVU - VYZ
CS.      NTFUN(M+6)   - BROJ TACAKA ZA KRIVU - VZX
C .
CE.      NTFUN(M+7)   - NUMBER OF POINTS - GXY
CE.      NTFUN(M+8)   - NUMBER OF POINTS - GYZ
CE.      NTFUN(M+9)   - NUMBER OF POINTS - GZX
CS.      NTFUN(M+7)   - BROJ TACAKA ZA KRIVU - GXY
CS.      NTFUN(M+8)   - BROJ TACAKA ZA KRIVU - GYZ
CS.      NTFUN(M+9)   - BROJ TACAKA ZA KRIVU - GZX
C .
CE.      NTFUN(M+10)  - NUMBER OF POINTS - AX
CE.      NTFUN(M+11)  - NUMBER OF POINTS - AY
CE.      NTFUN(M+12)  - NUMBER OF POINTS - AZ
CS.      NTFUN(M+10)  - BROJ TACAKA ZA KRIVU - AX
CS.      NTFUN(M+11)  - BROJ TACAKA ZA KRIVU - AY
CS.      NTFUN(M+12)  - BROJ TACAKA ZA KRIVU - AZ
C .
CE.      NTFUN(M+13)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( X )
CE.      NTFUN(M+14)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( Y )
CE.      NTFUN(M+15)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( Z )
CE.      NTFUN(M+16)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( XY )
CE.      NTFUN(M+17)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( YZ )
CE.      NTFUN(M+18)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( ZX )
CS.      NTFUN(M+13)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( X )
CS.      NTFUN(M+14)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( Y )
CS.      NTFUN(M+15)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( Z )
CS.      NTFUN(M+16)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( XY )
CS.      NTFUN(M+17)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( YZ )
CS.      NTFUN(M+18)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( ZX )
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      CHARACTER*1  CDIR(3)
      CHARACTER*2  CPLN(3),CMARK
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION FUNMAT(4,MATE*18,*),TREF(*),NTFUN(*)
      DATA CDIR/'X','Y','Z'/,CPLN/'XY','YZ','ZX'/
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD17'
      NFUN=1
      MA=(MAT-1)*18
      DO 10 K=1,18
        IF (K.GE.13) NFUN=3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD17'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(3,MA+1,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(3,MA+1,1)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      K=0
      DO 30 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2010) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6010) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   30 CONTINUE
      DO 40 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2030) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6030) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   40 CONTINUE
      DO 50 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2080) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6080) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   50 CONTINUE
      DO 60 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2040) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6040) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   60 CONTINUE
      CMARK='  '
      DO 70 L=1,6
        IF (L.LE.3) THEN
          CMARK=CDIR(L)
        ELSE
          CMARK=CPLN(L-3)
        END IF
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2060) CMARK,NTFUN(MA+K)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6060) CMARK,NTFUN(MA+K)
        WRITE(IZLAZ,5030) ((FUNMAT(I,MA+K,J),I=1,4),J=1,NTFUN(MA+K))
   70 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) FUNMAT(3,MA+1,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) FUNMAT(3,MA+1,1)
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(9X,1PD12.5,6X,1PD12.5,6X,1PD12.5,6X,1PD12.5)
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =    17  (MIZESOV TERMO-ELASTO-PLASTICAN '/
     &38X,'ORTOTROPAN SA MESOVITIM OJACANJEM)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E',A1)
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V',A2)
 2080 FORMAT(//
     &7X,'M O D U L    S M I C A NJ A  -  G',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'G',A2)
 2060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A ',
     &   '(MATERIJALNE CONSTANTE)'//
     &11X,'( ',A2,' )'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',3(10X,'FUNKCIJA')/
     &10X,'TEMPERATURA',10X,'TAUY',15X,'CY',16X,'AN')
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A',A1)
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2070 FORMAT(//7X,
     &'K O E F I C I J E N T    M E S O V I T O G    O J A C A NJ A'//
     &7X,'EM =',1PD12.5)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =   17  (VON MISES THERMO-ELASTIC-PLASTI',
     &'C'/38X,'ORTOTROPIC WITH MIXED HARDENING)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E',A1)
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V',A2)
 6080 FORMAT(//
     &7X,'S H E A R    M O D U L U S  -  G',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'G',A2)
 6060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A S ',
     &   '(MATERIAL CONSTANTS)'//
     &11X,'( ',A2,' )'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',3(10X,'FUNCTION')/
     &10X,'TEMPERATURE',10X,'TAUY',15X,'CY',16X,'AN')
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A',A1)
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6070 FORMAT(//7X,
     &'M I X E D    H A R D E N I N G    C O E F F I C I E N T'//
     &7X,'EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD18(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 18
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 18
C .
CE.   M=(MAT-1)*13  (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*13  (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - EX
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - YOUNG*S MODULUS - EY
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - YOUNG*S MODULUS - EZ
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - EX
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - MODUL ELASTICNOSTI - EY
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - MODUL ELASTICNOSTI - EZ
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - POISSON*S RATIO - VXY
CE.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+5,J) - POISSON*S RATIO - VYZ
CE.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+6,J) - POISSON*S RATIO - VZX
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - POISSONOV BROJ - VXY
CS.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+5,J) - POISSONOV BROJ - VYZ
CS.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+6,J) - POISSONOV BROJ - VZX
C .
CE.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+7,J) - SHEAR MODULUS - GXY
CE.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+8,J) - SHEAR MODULUS - GYZ
CE.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+9,J) - SHEAR MODULUS - GZX
CS.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+7,J) - MODUL SMICANJA - GXY
CS.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+8,J) - MODUL SMICANJA - GYZ
CS.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+9,J) - MODUL SMICANJA - GZX
C .
CE.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+10,J) - COEFFICIENT OF THERMAL EXPANSION - AX
CE.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+11,J) - COEFFICIENT OF THERMAL EXPANSION - AY
CE.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+12,J) - COEFFICIENT OF THERMAL EXPANSION  - AZ
CS.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+10,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AX
CS.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+11,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AY
CS.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+12,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AZ
C .
CE.       MATERIAL COEFFICIENTS (N)
CS.       KOEFICIJENTI OBLIKA (N)
C .
CE.              E(1) - N1
CE.              E(2) - N2
CE.              E(3) - N3
CE.              E(4) - NXY
CS.              E(5) - NYZ
CS.              E(6) - NXZ
C .
CE.         TREF(MAT) - REFERENCE TEMPERATURE
CS.         TREF(MAT) - REFERENTNA TEMPERATURA
C .
CE.  FUNMAT(2,M+13,1) - TRANSLATION COEFFICIENT OF CREEP POTENTIAL
CS.  FUNMAT(2,M+13,1) - KOEFICIJENT POMERANJA POTENCIJALA PUZANJA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS - EX
CE.      NTFUN(M+2)   - NUMBER OF POINTS - EY
CE.      NTFUN(M+3)   - NUMBER OF POINTS - EZ
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - EX
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - EY
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - EZ
C .
CE.      NTFUN(M+4)   - NUMBER OF POINTS - VXY
CE.      NTFUN(M+5)   - NUMBER OF POINTS - VYZ
CE.      NTFUN(M+6)   - NUMBER OF POINTS - VZX
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVU - VXY
CS.      NTFUN(M+5)   - BROJ TACAKA ZA KRIVU - VYZ
CS.      NTFUN(M+6)   - BROJ TACAKA ZA KRIVU - VZX
C .
CE.      NTFUN(M+7)   - NUMBER OF POINTS - GXY
CE.      NTFUN(M+8)   - NUMBER OF POINTS - GYZ
CE.      NTFUN(M+9)   - NUMBER OF POINTS - GZX
CS.      NTFUN(M+7)   - BROJ TACAKA ZA KRIVU - GXY
CS.      NTFUN(M+8)   - BROJ TACAKA ZA KRIVU - GYZ
CS.      NTFUN(M+9)   - BROJ TACAKA ZA KRIVU - GZX
C .
CE.      NTFUN(M+10)  - NUMBER OF POINTS - AX
CE.      NTFUN(M+11)  - NUMBER OF POINTS - AY
CE.      NTFUN(M+12)  - NUMBER OF POINTS - AZ
CS.      NTFUN(M+10)  - BROJ TACAKA ZA KRIVU - AX
CS.      NTFUN(M+11)  - BROJ TACAKA ZA KRIVU - AY
CS.      NTFUN(M+12)  - BROJ TACAKA ZA KRIVU - AZ
C .
C .
CE.      EQUIVALENT CREEP CURVE
CS.      EKVIVALENTNA KRIVA PUZANJA
C.
CE.      ICLAW - CREEP LAW FORM INDICATOR
CE.            = 1 BAILEY-NORTON CREEP LAW    ( NACR=3 )
CE.            = 2 EXPONENTIAL CREEP LAW      ( NACR=7 )
CE.            = 3 EIGHT PARAMETER CREEP LAW  ( NACR=8 )
CE.            = 4 CREEP LAW - USER DEFINED
CS.      ICLAW - INDIKATOR OBLIKA ZAKONA PUZANJA
CS.            = 1 BAILEY-NORTON-OV ZAKON PUZANJA   ( NACR=3 )
CS.            = 2 EKSPONENCIJALNI ZAKON PUZANJA    ( NACR=7 )
CS.            = 3 ZAKON PUZANJA SA 8 PARAMETARA    ( NACR=8 )
CS.            = 4 KORISNIKOV ZAKON PUZANJA
C.
CE.      ICELAW - EXPONENTIAL CREEP LAW FORM INDICATOR
CE.             = 1 SNYDER & BATHE EXPONENTIAL CREEP LAW
CS.      ICELAW - EKSPONENCIJALNI INDIKATOR OBLIKA ZAKONA PUZANJA
CS.             = 1 SNYDER & BATHE EKSPONENCIJALNI ZAKON PUZANJA
C.
CE.  FOR ICLAW = 1,2,3
CS.  ZA  ICLAW = 1,2,3
C.
CE.  ACR(I),I=1,NACR - MATERIAL CONSTANTS, TEMPERATURE INDEPENDENT
CS.  ACR(I),I=1,NACR - MATERIJALNE KONSTANTE, NEZAVISNE OD TEMPERATURE
C.
CE.  FOR ICLAW = 4
CS.  ZA  ICLAW = 4
C.
CE.           USTR1  - CREEP LAW (STRAIN)
CS.           USTR1  - ZAKON PUZANJA (DEFORMACIJA)
C.
CE.           USTR2  - CREEP LAW (STRAIN RATE)
CS.           USTR2  - ZAKON PUZANJA (BRZINA DEFORMACIJA)
C ......................................................................
C
      CHARACTER*250 ACOZ,USTR1,USTR2,IZRAZ,RP1,RP2
      CHARACTER*17 FRMT
      CHARACTER*41 FMTS,FMTE
      CHARACTER*1  CDIR(3)
      CHARACTER*2  CPLN(3)
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ECOEF/ ECR(6)
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      DIMENSION FUNMAT(2,MATE*13,*),TREF(*),NTFUN(*),NACR(3)
      DATA CDIR/'X','Y','Z'/,CPLN/'XY','YZ','ZX'/
      DATA NACR/3,7,7/
      DATA FRMT/'(5X,7(8X,1HA,I1))'/
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD18'
      MA=(MAT-1)*13
      DO 10 K=1,12
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD18'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,2)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,2)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (ECR(J),J=1,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (ECR(J),J=1,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(2,MA+13,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(2,MA+13,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) ICLAW,ICELAW
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) ICLAW,ICELAW
C
      IF (ICLAW.NE.4) THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,*) (ACR(I),I=1,NACR(ICLAW))
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,1010) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) THEN
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     &    READ(IULAZ,*) ACR(8)
          IF(INDFOR.EQ.2)
     &    READ(ACOZ,1010) ACR(8)
        END IF
      ELSE
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR1
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR1
        DO 29 I=1,80
          RP1(I:I)=' '
   29     IZRAZ(I:I)=' '
        CALL CITAC(USTR1,IZRAZ,ACR,ISG1,ITI1,ITH1)
        RP1(1:1)='$'
        CALL KONVEB(IZRAZ,RP1)
        DO 27 I=1,80
          RP2(I:I)=' '
   27     IZRAZ(I:I)=' '
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR2
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR2
        CALL CITAC(USTR2,IZRAZ,BCR,ISG2,ITI2,ITH2)
        RP2(1:1)='$'
        CALL KONVEB(IZRAZ,RP2)
      END IF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      K=0
      DO 30 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2010) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6010) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   30 CONTINUE
      DO 40 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2030) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6030) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   40 CONTINUE
      DO 50 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2080) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6080) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   50 CONTINUE
      DO 60 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2040) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6040) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   60 CONTINUE
      IF(ISRPS.EQ.0) WRITE(IZLAZ,2060) (L,L=1,6)
      IF(ISRPS.EQ.1) WRITE(IZLAZ,6060) (L,L=1,6)
      WRITE(IZLAZ,5040) (ECR(J),J=1,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) FUNMAT(2,MA+13,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) FUNMAT(2,MA+13,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      IF(ICLAW.EQ.1) THEN
        FMTS='(20X,30HBAILEY-NORTON-OV ZAKON PUZANJA//)'
        FMTE='(20X,30HBAILEY-NORTON CREEP LAW       //)'
      ELSE IF(ICLAW.EQ.2) THEN
        FMTS='(20X,30HEKSPONENCIJALNI ZAKON PUZANJA //)'
        FMTE='(20X,30HEKSPONENTIAL CREEP LAW        //)'
      ELSE IF(ICLAW.EQ.3) THEN
        FMTS='(20X,30HZAKON PUZANJA SA 8 PARAMETARA //)'
        FMTE='(20X,30HEIGHT PARAMETER CREEP LAW     //)'
      ELSE
        FMTS='(20X,30HKORISNIKOV ZAKON PUZANJA      //)'
        FMTE='(20X,30HCREEP LAW - USER DEFINED      //)'
      END IF
C
      IF (ICLAW.NE.4) THEN
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,2090)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,6090)
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,FMTS)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,FMTE)
        IF (ICLAW.EQ.1) FRMT(5:5)='3'
        WRITE(IZLAZ,FRMT)(I,I=0,NACR(ICLAW)-1)
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(43X,''A7'')')
        WRITE(IZLAZ,5030) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(35X,1PD10.3)') ACR(8)
      ELSE
        WRITE(IZLAZ,'(10X,''S T R A I N''//A80)')USTR1
        WRITE(IZLAZ,'(/10X,''S T R A I N    R A T E''//A80)')USTR2
      END IF
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(5X,7(1PD10.3))
 5040 FORMAT(5X,6(1PD12.3))
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =    18  (TERMO-ELASTICAN SA PUZANJEM '/
     &38X,'- ORTOTROPAN)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E',A1)
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V',A2)
 2080 FORMAT(//
     &7X,'M O D U L    S M I C A NJ A  -  G',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'G',A2)
 2060 FORMAT(//
     &11X,'K O E F I C I J E N T I    O B L I K A  ',
     &   '(POTENCIJAL PUZANJA)'//
     &11X,6('N',I1,10X)//)
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A',A1)
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2070 FORMAT(//7X,
     &'K O E F I C I J E N T    P O M E R A NJ A  (POTENCIJAL PUZANJA)'
     &//7X,'EMC =',1PD12.5)
 2090 FORMAT(//
     &7X,'Z A K O N    P U Z A NJ A   (MATERIJALNE CONSTANTE)'//)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =   18  (THERMO-ELASTIC WITH CREEP '/
     &38X,'- ORTOTROPIC)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E',A1)
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V',A2)
 6080 FORMAT(//
     &7X,'S H E A R    M O D U L U S  -  G',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'G',A2)
 6060 FORMAT(//
     &11X,'M A T E R I A L    C O E F F I C I E N T S ',
     &   '(CREEP POTENTIAL)'//
     &11X,6('N',I1,10X)//)
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A',A1)
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6070 FORMAT(//7X,
     &'T R A N S L A T I O N    C O E F F I C I E N T    (CREEP ',
     &'POTENTIAL)'//7X,'EMC =',1PD12.5)
 6090 FORMAT(//7X,'C R E E P    L A W    (MATERIAL CONSTANTS)'//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD19(FUNMAT,TREF,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 19
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 19
C .
CE.   M=(MAT-1)*18  (MAT - MATERIAL NUMBER)
CE.   J=1,(NUMBER OF POINTS ON CURVE)
CS.   M=(MAT-1)*18  (MAT - MATERIJAL BROJ)
CS.   J=1,(DO BROJA TACAKA ZA KRIVU)
C .
CE.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+1,J) - YOUNG*S MODULUS - EX
CE.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+2,J) - YOUNG*S MODULUS - EY
CE.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+3,J) - YOUNG*S MODULUS - EZ
CS.   FUNMAT(1,M+1,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+1,J) - MODUL ELASTICNOSTI - EX
CS.   FUNMAT(1,M+2,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+2,J) - MODUL ELASTICNOSTI - EY
CS.   FUNMAT(1,M+3,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+3,J) - MODUL ELASTICNOSTI - EZ
C .
CE.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+4,J) - POISSON*S RATIO - VXY
CE.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+5,J) - POISSON*S RATIO - VYZ
CE.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+6,J) - POISSON*S RATIO - VZX
CS.   FUNMAT(1,M+4,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+4,J) - POISSONOV BROJ - VXY
CS.   FUNMAT(1,M+5,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+5,J) - POISSONOV BROJ - VYZ
CS.   FUNMAT(1,M+6,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+6,J) - POISSONOV BROJ - VZX
C .
CE.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+7,J) - SHEAR MODULUS - GXY
CE.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+8,J) - SHEAR MODULUS - GYZ
CE.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURE
CE.   FUNMAT(2,M+9,J) - SHEAR MODULUS - GZX
CS.   FUNMAT(1,M+7,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+7,J) - MODUL SMICANJA - GXY
CS.   FUNMAT(1,M+8,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+8,J) - MODUL SMICANJA - GYZ
CS.   FUNMAT(1,M+9,J) - ARGUMENT - TEMPERATURA
CS.   FUNMAT(2,M+9,J) - MODUL SMICANJA - GZX
C .
CE.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+10,J) - COEFFICIENT OF THERMAL EXPANSION - AX
CE.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+11,J) - COEFFICIENT OF THERMAL EXPANSION - AY
CE.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+12,J) - COEFFICIENT OF THERMAL EXPANSION  - AZ
CS.  FUNMAT(1,M+10,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+10,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AX
CS.  FUNMAT(1,M+11,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+11,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AY
CS.  FUNMAT(1,M+12,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+12,J) - KOEFICIJENT TEMPERATURSKOG SIRENJA - AZ
C .
CE.       RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (X - DIRECTION)
CS.       RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (X - PRAVAC)
C .
CE.  FUNMAT(1,M+13,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+13,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+13,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+13,J) - EXPONENT - AN
CS.  FUNMAT(1,M+13,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+13,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+13,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+13,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (Y - DIRECTION)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (Y - PRAVAC)
C .
CE.  FUNMAT(1,M+14,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+14,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+14,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+14,J) - EXPONENT - AN
CS.  FUNMAT(1,M+14,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+14,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+14,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+14,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (Z - DIRECTION)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (Z - PRAVAC)
C .
CE.  FUNMAT(1,M+15,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+15,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+15,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+15,J) - EXPONENT - AN
CS.  FUNMAT(1,M+15,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+15,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+15,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+15,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (XY - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (XY - RAVAN)
C .
CE.  FUNMAT(1,M+16,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+16,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+16,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+16,J) - EXPONENT - AN
CS.  FUNMAT(1,M+16,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+16,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+16,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+16,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (YZ - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (YZ - RAVAN)
C .
CE.  FUNMAT(1,M+17,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+17,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+17,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+17,J) - EXPONENT - AN
CS.  FUNMAT(1,M+17,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+17,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+17,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+17,J) - EKSPONENT - AN
C .
CE.      RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS (ZX - PLANE)
CS.      RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE (ZX - RAVAN)
C .
CE.  FUNMAT(1,M+18,J) - ARGUMENT - TEMPERATURE
CE.  FUNMAT(2,M+18,J) - YIELD STRESS - TAUY
CE.  FUNMAT(3,M+18,J) - MULTIPLAYER - CY
CE.  FUNMAT(4,M+18,J) - EXPONENT - AN
CS.  FUNMAT(1,M+18,J) - ARGUMENT - TEMPERATURA
CS.  FUNMAT(2,M+18,J) - NAPON TECENJA - TAUY
CS.  FUNMAT(3,M+18,J) - MNOZIOC - CY
CS.  FUNMAT(4,M+18,J) - EKSPONENT - AN
C .
CE.         TREF(MAT) - REFERENCE TEMPERATURE
CS.         TREF(MAT) - REFERENTNA TEMPERATURA
C .
CE.   FUNMAT(3,M+1,1) - MIXED HARDENING COEFFICIENT
CS.   FUNMAT(3,M+1,1) - KOEFICIJENT MESOVITOG OJACANJA
C .
CE.   FUNMAT(4,M+1,1) - TRANSLATION COEFFICIENT
CS.   FUNMAT(4,M+1,1) - KOEFICIJENT TRANSLACIJE POVRSI PUZANJA
C .
CE.      NTFUN(M+1)   - NUMBER OF POINTS - EX
CE.      NTFUN(M+2)   - NUMBER OF POINTS - EY
CE.      NTFUN(M+3)   - NUMBER OF POINTS - EZ
CS.      NTFUN(M+1)   - BROJ TACAKA ZA KRIVU - EX
CS.      NTFUN(M+2)   - BROJ TACAKA ZA KRIVU - EY
CS.      NTFUN(M+3)   - BROJ TACAKA ZA KRIVU - EZ
C .
CE.      NTFUN(M+4)   - NUMBER OF POINTS - VXY
CE.      NTFUN(M+5)   - NUMBER OF POINTS - VYZ
CE.      NTFUN(M+6)   - NUMBER OF POINTS - VZX
CS.      NTFUN(M+4)   - BROJ TACAKA ZA KRIVU - VXY
CS.      NTFUN(M+5)   - BROJ TACAKA ZA KRIVU - VYZ
CS.      NTFUN(M+6)   - BROJ TACAKA ZA KRIVU - VZX
C .
CE.      NTFUN(M+7)   - NUMBER OF POINTS - GXY
CE.      NTFUN(M+8)   - NUMBER OF POINTS - GYZ
CE.      NTFUN(M+9)   - NUMBER OF POINTS - GZX
CS.      NTFUN(M+7)   - BROJ TACAKA ZA KRIVU - GXY
CS.      NTFUN(M+8)   - BROJ TACAKA ZA KRIVU - GYZ
CS.      NTFUN(M+9)   - BROJ TACAKA ZA KRIVU - GZX
C .
CE.      NTFUN(M+10)  - NUMBER OF POINTS - AX
CE.      NTFUN(M+11)  - NUMBER OF POINTS - AY
CE.      NTFUN(M+12)  - NUMBER OF POINTS - AZ
CS.      NTFUN(M+10)  - BROJ TACAKA ZA KRIVU - AX
CS.      NTFUN(M+11)  - BROJ TACAKA ZA KRIVU - AY
CS.      NTFUN(M+12)  - BROJ TACAKA ZA KRIVU - AZ
C .
CE.      NTFUN(M+13)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( X )
CE.      NTFUN(M+14)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( Y )
CE.      NTFUN(M+15)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( Z )
CE.      NTFUN(M+16)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( XY )
CE.      NTFUN(M+17)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( YZ )
CE.      NTFUN(M+18)  - NUMBER OF POINTS ON CURVES - TAUY,CY,AN ( ZX )
CS.      NTFUN(M+13)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( X )
CS.      NTFUN(M+14)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( Y )
CS.      NTFUN(M+15)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( Z )
CS.      NTFUN(M+16)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( XY )
CS.      NTFUN(M+17)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( YZ )
CS.      NTFUN(M+18)  - BROJ TACAKA ZA KRIVE - TAUY,CY,AN ( ZX )
C .
C .
CE.      EQUVALENT CREEP CURVE
CS.      EKVIVALENTNA KRIVA PUZANJA
C.
CE.      ICLAW - CREEP LAW FORM INDICATOR
CE.            = 1 BAILEY-NORTON CREEP LAW    ( NACR=3 )
CE.            = 2 EXPONENTIAL CREEP LAW      ( NACR=7 )
CE.            = 3 EIGHT PARAMETER CREEP LAW  ( NACR=8 )
CE.            = 4 CREEP LAW - USER DEFINED
CS.      ICLAW - INDIKATOR OBLIKA ZAKONA PUZANJA
CS.            = 1 BAILEY-NORTON-OV ZAKON PUZANJA   ( NACR=3 )
CS.            = 2 EKSPONENCIJALNI ZAKON PUZANJA    ( NACR=7 )
CS.            = 3 ZAKON PUZANJA SA 8 PARAMETARA    ( NACR=8 )
CS.            = 4 KORISNIKOV ZAKON PUZANJA
C.
CE.      ICELAW - EXPONENTIAL CREEP LAW FORM INDICATOR
CE.             = 1 SNYDER & BATHE EXPONENTIAL CREEP LAW
CS.      ICELAW - EKSPONENCIJALNI INDIKATOR OBLIKA ZAKONA PUZANJA
CS.             = 1 SNYDER & BATHE EKSPONENCIJALNI ZAKON PUZANJA
C.
CE.  FOR ICLAW = 1,2,3
CS.  ZA  ICLAW = 1,2,3
C.
CE.  ACR(I),I=1,NACR - MATERIAL CONSTANTS, TEMPERATURE INDEPENDENT
CS.  ACR(I),I=1,NACR - MATERIJALNE KONSTANTE, NEZAVISNE OD TEMPERATURE
C.
CE.  FOR ICLAW = 4
CS.  ZA  ICLAW = 4
C.
CE.           USTR1  - CREEP LAW (STRAIN)
CS.           USTR1  - ZAKON PUZANJA (DEFORMACIJA)
C.
CE.           USTR2  - CREEP LAW (STRAIN RATE)
CS.           USTR2  - ZAKON PUZANJA (BRZINA DEFORMACIJA)
C ......................................................................
C
      CHARACTER*250 ACOZ,USTR1,USTR2,IZRAZ,RP1,RP2
      CHARACTER*17 FRMT
      CHARACTER*41 FMTS,FMTE
      CHARACTER*1  CDIR(3)
      CHARACTER*2  CPLN(3),CMARK
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ECOEF/ ECR(6)
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /CREEPC/ ACR(20),BCR(20)
      COMMON /CREEPL/ RP1,RP2
C
      DIMENSION FUNMAT(4,MATE*18,*),TREF(*),NTFUN(*),NACR(3)
      DATA CDIR/'X','Y','Z'/,CPLN/'XY','YZ','ZX'/
      DATA NACR/3,7,7/
      DATA FRMT/'(5X,7(8X,1HA,I1))'/
C

      IF(IDEBUG.GT.0) PRINT *, ' UMOD19'
      NFUN=1
      MA=(MAT-1)*18
      DO 10 K=1,18
        IF (K.GE.13) NFUN=3
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MA+K)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MA+K)
        IF(NTFUN(MA+K).GT.MAXT.OR.NTFUN(MA+K).LE.0) STOP 'UMOD19'
        DO 20 J=1,NTFUN(MA+K)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MA+K,J),I=1,NFUN+1)
   20   CONTINUE
   10 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (ECR(J),J=1,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (ECR(J),J=1,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) TREF(MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) TREF(MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(3,MA+1,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(3,MA+1,1)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(4,MA+1,1)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) FUNMAT(4,MA+1,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) ICLAW,ICELAW
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) ICLAW,ICELAW
C
      IF (ICLAW.NE.4) THEN
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,*) (ACR(I),I=1,NACR(ICLAW))
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,1010) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) THEN
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     &    READ(IULAZ,*) ACR(8)
          IF(INDFOR.EQ.2)
     &    READ(ACOZ,1010) ACR(8)
        END IF
      ELSE
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR1
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR1
        DO 29 I=1,80
          RP1(I:I)=' '
   29     IZRAZ(I:I)=' '
        CALL CITAC(USTR1,IZRAZ,ACR,ISG1,ITI1,ITH1)
        RP1(1:1)='$'
        CALL KONVEB(IZRAZ,RP1)
        DO 27 I=1,80
          RP2(I:I)=' '
   27     IZRAZ(I:I)=' '
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     &  READ(IULAZ,'(A80)')USTR2
        IF(INDFOR.EQ.2)
     &  READ(ACOZ,'(A80)')USTR2
        CALL CITAC(USTR2,IZRAZ,BCR,ISG2,ITI2,ITH2)
        RP2(1:1)='$'
        CALL KONVEB(IZRAZ,RP2)
      END IF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      K=0
      DO 30 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2010) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6010) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   30 CONTINUE
      DO 40 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2030) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6030) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   40 CONTINUE
      DO 50 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2080) CPLN(L),NTFUN(MA+K),CPLN(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6080) CPLN(L),NTFUN(MA+K),CPLN(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   50 CONTINUE
      DO 60 L=1,3
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2040) CDIR(L),NTFUN(MA+K),CDIR(L)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6040) CDIR(L),NTFUN(MA+K),CDIR(L)
        WRITE(IZLAZ,5020) ((FUNMAT(I,MA+K,J),I=1,2),J=1,NTFUN(MA+K))
   60 CONTINUE
      CMARK='  '
      DO 70 L=1,6
        IF (L.LE.3) THEN
          CMARK=CDIR(L)
        ELSE
          CMARK=CPLN(L-3)
        END IF
        K=K+1
        IF(ISRPS.EQ.0) WRITE(IZLAZ,2060) CMARK,NTFUN(MA+K)
        IF(ISRPS.EQ.1) WRITE(IZLAZ,6060) CMARK,NTFUN(MA+K)
        WRITE(IZLAZ,5030) ((FUNMAT(I,MA+K,J),I=1,4),J=1,NTFUN(MA+K))
   70 CONTINUE
      IF(ISRPS.EQ.0) WRITE(IZLAZ,2065) (L,L=1,6)
      IF(ISRPS.EQ.1) WRITE(IZLAZ,6065) (L,L=1,6)
      WRITE(IZLAZ,5050) (ECR(J),J=1,6)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TREF(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TREF(MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) FUNMAT(3,MA+1,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) FUNMAT(3,MA+1,1)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) FUNMAT(4,MA+1,1)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) FUNMAT(4,MA+1,1)
C
CS    ZAKON PUZANJA
CE    CREEP LAW
C
      IF(ICLAW.EQ.1) THEN
        FMTS='(20X,30HBAILEY-NORTON-OV ZAKON PUZANJA//)'
        FMTE='(20X,30HBAILEY-NORTON CREEP LAW       //)'
      ELSE IF(ICLAW.EQ.2) THEN
        FMTS='(20X,30HEKSPONENCIJALNI ZAKON PUZANJA //)'
        FMTE='(20X,30HEKSPONENTIAL CREEP LAW        //)'
      ELSE IF(ICLAW.EQ.3) THEN
        FMTS='(20X,30HZAKON PUZANJA SA 8 PARAMETARA //)'
        FMTE='(20X,30HEIGHT PARAMETER CREEP LAW     //)'
      ELSE
        FMTS='(20X,30HKORISNIKOV ZAKON PUZANJA      //)'
        FMTE='(20X,30HCREEP LAW - USER DEFINED      //)'
      END IF
C
      IF (ICLAW.NE.4) THEN
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,2090)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,6090)
        IF(ISRPS.EQ.0)
     &  WRITE(IZLAZ,FMTS)
        IF(ISRPS.EQ.1)
     &  WRITE(IZLAZ,FMTE)
        IF (ICLAW.EQ.1) FRMT(5:5)='3'
        WRITE(IZLAZ,FRMT)(I,I=0,NACR(ICLAW)-1)
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(43X,''A7'')')
        WRITE(IZLAZ,5040) (ACR(I),I=1,NACR(ICLAW))
        IF (ICLAW.EQ.3) WRITE(IZLAZ,'(35X,1PD10.3)') ACR(8)
      ELSE
        WRITE(IZLAZ,'(10X,''S T R A I N''//A80)')USTR1
        WRITE(IZLAZ,'(/10X,''S T R A I N    R A T E''//A80)')USTR2
      END IF
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(9X,1PD12.5,6X,1PD12.5)
 5030 FORMAT(9X,1PD12.5,6X,1PD12.5,6X,1PD12.5,6X,1PD12.5)
 5040 FORMAT(5X,7(1PD10.3))
 5050 FORMAT(5X,6(1PD12.3))
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     &'MODEL MATERIJALA BROJ =    19  (MIZESOV TERMO-ELASTO-PLASTICAN '/
     &24X,'SA PUZANJEM - ORTOTROPAN SA MESOVITIM OJACANJEM)'///
     &11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     &7X,'M O D U L    E L A S T I C N O S T I  -  E',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'E',A1)
 2030 FORMAT(//
     &7X,'P O I S S O N O V    B R O J  -  V',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'V',A2)
 2080 FORMAT(//
     &7X,'M O D U L    S M I C A NJ A  -  G',A2//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'G',A2)
 2060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A ',
     &   '(MATERIJALNE CONSTANTE)'//
     &11X,'( ',A2,' )'//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',3(10X,'FUNKCIJA')/
     &10X,'TEMPERATURA',10X,'TAUY',15X,'CY',16X,'AN')
 2065 FORMAT(//
     &11X,'K O E F I C I J E N T I    O B L I K A  ',
     &   '(POTENCIJAL PUZANJA)'//
     &11X,6('N',I1,10X)//)
 2040 FORMAT(//
     &7X,'K O E F I C I J E N T    T E R M I C K O G    S I R E N J A',
     &'  -  A',A1//
     &11X,'BROJ TACAKA ZA KRIVU =',I5//
     &11X,'ARGUMENT',10X,'FUNKCIJA'/
     &10X,'TEMPERATURA',11X,'A',A1)
 2050 FORMAT(//7X,
     &'R E F E R E N T N A   T E M P E R A T U R A'//
     &7X,'TREF =',1PD12.5)
 2070 FORMAT(//7X,
     &'K O E F I C I J E N T    M E S O V I T O G    O J A C A NJ A'//
     &7X,'EM =',1PD12.5)
 2020 FORMAT(//7X,
     &'K O E F I C I J E N T    P O M E R A NJ A  (POTENCIJAL PUZANJA)'
     &//7X,'EMC =',1PD12.5)
 2090 FORMAT(//
     &7X,'Z A K O N    P U Z A NJ A   (MATERIJALNE CONSTANTE)'//)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     &'MATERIAL MODEL NUMBER  =   19  (VON MISES THERMO-ELASTIC-PLASTI',
     &'C'/24X,'WITH CREEP - ORTOTROPIC WITH MIXED HARDENING)'///
     &11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     &7X,'Y O U N G S      M O D U L U S  -  E',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'E',A1)
 6030 FORMAT(//
     &7X,'P O I S S O N S      N U M B E R  -  V',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'V',A2)
 6080 FORMAT(//
     &7X,'S H E A R    M O D U L U S  -  G',A2//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'G',A2)
 6060 FORMAT(//
     &7X,'R A M B E R G  -  O S G O O D    F O R M U L A S ',
     &   '(MATERIAL CONSTANTS)'//
     &11X,'( ',A2,' )'//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',3(10X,'FUNCTION')/
     &10X,'TEMPERATURE',10X,'TAUY',15X,'CY',16X,'AN')
 6065 FORMAT(//
     &11X,'M A T E R I A L    C O E F F I C I E N T S ',
     &   '(CREEP POTENTIAL)'//
     &11X,6('N',I1,10X)//)
 6040 FORMAT(//
     &7X,'C O E F F I C I E N T   O F   T H E R M A L   ',
     &    'E X P A N S I O N  -  A',A1//
     &11X,'NUMBER OF POINTS ON CURVE =',I5//
     &11X,'ARGUMENT',10X,'FUNCTION'/
     &10X,'TEMPERATURE',11X,'A',A1)
 6050 FORMAT(//7X,
     &'R E F E R E N C E    T E M P E R A T U R E'//
     &7X,'TREF =',1PD12.5)
 6070 FORMAT(//7X,
     &'M I X E D    H A R D E N I N G    C O E F F I C I E N T'//
     &7X,'EM =',1PD12.5)
 6020 FORMAT(//7X,
     &'T R A N S L A T I O N    C O E F F I C I E N T    (CREEP ',
     &'POTENTIAL)'//7X,'EMC =',1PD12.5)
 6090 FORMAT(//7X,'C R E E P    L A W    (MATERIAL CONSTANTS)'//)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE UMOD21(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 21
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 21
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1-3,MAT)   - YOUNG*S MODULUS EX,EY,EZ
CE.   FUNMAT(4-6,MAT)   - POISSON*S RATIO  VXY, VYZ, VZX
CE.   FUNMAT(7-9,MAT)   - SHEAR MODULUS GXY, GYZ, GZX
CE.   FUNMAT(10-15,MAT) - YIELD STRESSES X,Y,Z, XY, YZ, ZX
CE.   FUNMAT(16-21,MAT) - RAMBERG-OSGOOD CONSTANTS (C)  X,Y,Z,XY,YZ,ZX
CE.   FUNMAT(22-27,MAT) - RAMBERG-OSGOOD CONSTANTS (AN) X,Y,Z,XY,YZ,ZX
CE.   FUNMAT(28,MAT)    - EFFECTIVE YIELD STRESSES TAUY
CE.   FUNMAT(29,MAT)    - MIXED HARDENING PARAMETER  M
CE.   FUNMAT(30,MAT)    - VISCOSITY COEFFICIENT ETA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(30,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD10'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
C... NAPONI TECENJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=10,15)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=10,15)
C... RAMBERG OSGOOD MNOZIOCI
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=16,21)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=16,21)
C... RAMBERG OSGOOD EKSPONENTI
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=22,27)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=22,27)
C... KOEFICIJENT MESOVITOG OJACANJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(29,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(29,MAT)
C... KOEFICIJENT VISKOZNOSTI
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(30,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(30,MAT)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      FUNMAT(28,MAT)=FUNAR(FUNMAT(10,MAT))
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,30)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,30)
C... CHECK MATERIAL CONSTANTS
      CALL       ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    21  (VISKO-PLASTICAN ANIZOTROPAN)'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     1//16X,'N A P O N I    T E C E N J A'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     1//16X,'M N O Z I O C I'             //11X,'CX ',9X,'CY ',
     19X,'CZ ',9X,'CXY',9X,'CYZ',9X,'CZX'/7X,6(1PD11.4,1X)//
     1//16X,'E K S P O N E N T I'         //11X,'AX ',9X,'AY ',
     19X,'AZ ',9X,'AXY',9X,'AYZ',9X,'AZX'/7X,6(1PD11.4,1X)//
     116X,'N A P O N    T E C E N J A (EFEKTIV.) TAUY =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5/// 
     116X,'K O E F   V I S K O Z N O S T I  ...   ETA =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    21  (VISCOPLASTICITY ORTHOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U S'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X)//
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     1//16X,'Y I E L D    S T R E S S E S'//11X,' X ',9X,' Y ',
     19X,' Z ',9X,' XY',9X,' YZ',9X,' ZX'/7X,6(1PD11.4,1X)//
     116X,'M U L T I P L A Y E R S'       //11X,'CX ',9X,'CY ',
     19X,'CZ ',9X,'CXY',9X,'CYZ',9X,'CZX'/7X,6(1PD11.4,1X)//
     1//16X,'E X P O N E N T S'           //11X,'AX ',9X,'AY ',
     19X,'AZ ',9X,'AXY',9X,'AYZ',9X,'AZX'/7X,6(1PD11.4,1X)//
     116X,'Y I E L D    S T R E S S (EFFECTIVE)  TAUY =',1PD12.5///
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5///
     116X,'V I S C O S I T Y   C O E F F  .....   ETA =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD22(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 22
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 22
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - AEM
CE.   FUNMAT(2,MAT) - CONSTANT  - ALAM
CE.   FUNMAT(3,MAT) - CONSTANT  - AKA
CE.   FUNMAT(4,MAT) - CONSTANT  - AE1
CE.   FUNMAT(5,MAT) - CONSTANT  - AE0
CS.   FUNMAT(1,MAT) - KONSTANTA - AEM
CS.   FUNMAT(2,MAT) - KONSTANTA - ALAM
CS.   FUNMAT(3,MAT) - KONSTANTA - AKA
CS.   FUNMAT(4,MAT) - KONSTANTA - AE1
CS.   FUNMAT(5,MAT) - KONSTANTA - AE0
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(10,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD22'
      IEL=0
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3),IEL
      FUNMAT(4,MAT)=IEL
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,3),IEL
      FUNMAT(4,MAT)=IEL
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=5,10)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=5,10)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,3),IEL
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,3),IEL
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=5,10)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=5,10)
      RETURN
 1000 FORMAT(7F10.0)
 1010 FORMAT(3F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    22  (TLO - VISKO PLASTICNOST)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//
     216X,'G R A N I C A   Z A T E Z A NJ A ...... AT =',1PD12.5//
     216X,'IDIKATOR ZADAVANJA MODULA ELASTIC. ....IEL =',I12//)
 2045 FORMAT(//
     16X,'M A T E R I J A L N E     K O N S T A N T E          '//
     2,6X,' EM',8X,'LAMBDA',7X,'KAPA',7X,'P0  ',8X,'  E0',8X,'ETA '/
     3' ',6(1PD12.4))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    22  (CAM-CLAY - VISCOPLASTICITY)'///
     211X,'MATERIAL NUMBER =',I5)
 6045 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     2,6X,' EM',6X,'ALAMDA',8X,'ETAS',10X,'P0',10X,'E0',10X,'ETA'/
     3' ',6(1PD12.4))
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//
     216X,'T E N S I O N      C U T O F F  ....... AT =',1PD12.5//
     216X,'I N D I C A T O R  OF  C H A N G E  E .IEL =',I12//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD23(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 23
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 23
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  EX
CE.   FUNMAT(2,MAT) - YOUNG*S MODULUS IN DIRECTION Y  -  EY
CE.   FUNMAT(3,MAT) - YOUNG*S MODULUS IN DIRECTION Z  -  EZ
C .
CE.   FUNMAT(4,MAT) - POISSON*S RATIO FOR PLANE XY    -  VXY
CE.   FUNMAT(5,MAT) - POISSON*S RATIO FOR PLANE YZ    -  VYZ
CE.   FUNMAT(6,MAT) - POISSON*S RATIO FOR PLANE ZX    -  VZX
C .
CE.   FUNMAT(7,MAT) - SHEAR MODULUS FOR PLANE XY  -  GXY
CE.   FUNMAT(8,MAT) - SHEAR MODULUS FOR PLANE YZ  -  GYZ
CE.   FUNMAT(9,MAT) - SHEAR MODULUS FOR PLANE ZX  -  GZX
C .
CE.   FUNMAT(10-15,MAT) - YIELD STRESSES X,Y,Z, XY, YZ, ZX
CE.   FUNMAT(16,MAT)    - EFFECTIVE YIELD STRESSES TAUY
CE.   FUNMAT(17,18,MAT) - RAMBERG-OSGOOD CONSTANTS CY,AN
CE.   FUNMAT(19,MAT)    - MIXED HARDENING PARAMETER  M
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(34,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD23'
      IEL=0
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,6),IEL
      FUNMAT(7,MAT)=IEL
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,6),IEL
      FUNMAT(7,MAT)=IEL
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=8,14)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=8,14)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=15,20)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=15,20)
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=21,26)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=21,26)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=27,33)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=27,33)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(34,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(34,MAT)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,34)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,34)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =  23 (MODEL TLA SA KAPOM -        
     1VISKOPLASTICAN  )'/
     2//11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F1=0 '/21X,'A',13X,'B ',
     112X,'C ',12X,'TETA'/16X,4(1PD12.5,2X)//
     120X,'AT  ',11X,'AI1A0 ',9X,'IEL '/16X,3(1PD12.5,2X)//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'
     1/20X,'W ',13X,'D ',11X,'ALFA'/16X,3(1PD12.5,2X)
     1/20X,'W1',13X,'D1',11X,'D2  ',11X,'D3 ' /16X,4(1PD12.5,2X)//
     1//16X,'K O E F I C I J E N T I  K R I V E  R=R(L)'
     1/20X,'R0',13X,'R1',11X,'R2  ',11X,'R3 ' /16X,4(1PD12.5,2X)
     1/20X,'R4',13X,'R5' /16X,2(1PD12.5,2X)//
     1//16X,'M O D U L    E L A S T I C N O S T I '
     1/20X,'KEI',13X,'KS',11X,'K1 ',11X,'K2 ' /16X,4(1PD12.5,2X)
     1/20X,'BETA',11X,'DELTA' /16X,2(1PD12.5,2X)//
     1//16X,'M O D U L    S M I C A N J A'/16X,'GE1',12X,'GS ',
     111X,'G1 ',11X,'G2 ',11X,'G3 '/10X,5(1PD12.5,2X)/
     119X,'ETA',11X,'GAMA '/10X,2(1PD12.5,2X)//
     116X,'K O E F I C I J E N T  V I S K O Z N O S T I '/21X,'ETA'
     1/16X,1PD12.5,2X)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =   23  (SOIL CAP MODEL - VISCOPLASTICITY)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'C O E F I C I J E N T S  F O R  C U R V E  F1=0 '
     1/21X,'A',13X,'B ',12X,'C ',12X,'TETA'/16X,4(1PD12.5,2X)//
     120X,'AT  ',11X,'AI1A0 ',9X,'IEL '/16X,3(1PD12.5,2X)//
     1//16X,'C O E F I C I J E N T S  F O R  H A R D E N I N G  
     1C U R V E'/20X,'W ',13X,'D ',11X,'ALFA'/16X,3(1PD12.5,2X)
     1/20X,'W1',13X,'D1',11X,'D2  ',11X,'D3 ' /16X,4(1PD12.5,2X)//
     1//16X,'C O E F I C I J E N T S  F O R  C U R V E  R=R(L)'
     1/20X,'R0',13X,'R1',11X,'R2  ',11X,'R3 ' /16X,4(1PD12.5,2X)
     1/20X,'R4',13X,'R5' /16X,2(1PD12.5,2X)//
     1//16X,'Y O U N G S  M O D U L U S '
     1/20X,'KEI',13X,'KS',11X,'K1 ',11X,'K2 ' /16X,4(1PD12.5,2X)
     1/20X,'BETA',11X,'DELTA' /16X,2(1PD12.5,2X)//
     1//10X,'S H E A R   M O D U L U S'/16X,'GE1',12X,'GS ',
     111X,'G1 ',11X,'G2 ',11X,'G3 '/10X,5(1PD12.5,2X)/
     119X,'ETA',11X,'GAMA '/10X,2(1PD12.5,2X)/
     116X,'V I S C O S I T Y  C O E F F I C I J E N T '/21X,'ETA'
     1/16X,1PD12.5,2X)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD24(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 24
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 24
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - TYPE OF MATERIAL, MTYPE: = 0 (DEFAULT 1)
CE.                                            = 1 MAXVELL TYPE
CE.                                            = 2 VOIGT MODEL
CE.                                            = 3 KELVIN
CE.                                                (STANDARD MODEL) 
CE.   FUNMAT(2,MAT) - YOUNG*S MODULUS E (=ER FOR MTYPE=3)
CE.   FUNMAT(3,MAT) - VISCOUS COEFFICIENT B
CE.   FUNMAT(4,MAT) - YOUNG*S MODULUS E1 FOR MTYPE=3
CS.   FUNMAT(1,MAT) - TIP MATERIJALA,   MTYPE: = 0 (DEFAULT 1)
CS.                                            = 1 MAXVELL-OV MODEL
CS.                                            = 2 VOIGT-OV MODEL
CS.                                            = 3 KELVIN-OV
CS.                                                (STANDARDNI MODEL) 
CS.   FUNMAT(2,MAT) - JUNGOV MODUL ELASTICNOSTI E (=ER ZA MTYPE=3)
CS.   FUNMAT(3,MAT) - VISKOZNI KOEFICIJENT OTPORA B SIGV = B*VDEF
CS.   FUNMAT(4,MAT) - JUNGOV MODUL ELASTICNOSTI E1 ZA MTYPE=3
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(5,*)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD24'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) MTYPE,(FUNMAT(I,MAT),I=2,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) MTYPE,(FUNMAT(I,MAT),I=2,5)
      CALL ISPITA(ACOZ)
      IF (MTYPE.EQ.0) MTYPE = 1
      FUNMAT(1,MAT) = MTYPE
      IF (MTYPE.NE.3) FUNMAT(4,MAT)=FUNMAT(5,MAT)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF (ISRPS.EQ.0) THEN
         IF(MTYPE.EQ.3) WRITE(IZLAZ,2010) MTYPE,(FUNMAT(I,MAT),I=2,5)
         IF(MTYPE.NE.3) WRITE(IZLAZ,2011) MTYPE,(FUNMAT(I,MAT),I=2,4)
      ELSE
         IF(MTYPE.EQ.3) WRITE(IZLAZ,6010) MTYPE,(FUNMAT(I,MAT),I=2,5)
         IF(MTYPE.NE.3) WRITE(IZLAZ,6011) MTYPE,(FUNMAT(I,MAT),I=2,4)
      ENDIF
      RETURN
C
 1000 FORMAT(I5,4F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    24  (LINEARNI VISKOELASTICNI)'
     2 ///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'TIP MODELA .............................. MTYPE =',I5/
     316X,'  = 1 MAXWELL-OV (EQ.0 - ZAMENIMO U 1) '/
     416X,'  = 2 VOIGT-OV'/
     516X,'  = 3 KELVINOV (STANDARDNI LINEARNI SOLID)'/       
     116X,'MODUL ELASTICNOSTI ......................... ER =',1PD12.5/
     116X,'KOEFICIJENT VISKOZNOG OTPORA ................ B =',1PD12.5/
     116X,'MODUL ELASTICNOSTI ......................... E1 =',1PD12.5/
     116X,'POASONOV BROJ ............................... V =',1PD12.5)
 2011 FORMAT(//
     116X,'TIP MODELA .............................. MTYPE =',I5/
     316X,'  = 1 MAXWELL-OV (EQ.0 - ZAMENIMO U 1) '/
     416X,'  = 2 VOIGT-OV'/
     516X,'  = 3 KELVINOV (STANDARDNI LINEARNI SOLID)'/       
     116X,'MODUL ELASTICNOSTI .........................  E =',1PD12.5/
     116X,'KOEFICIJENT VISKOZNOG OTPORA ................ B =',1PD12.5/
     116X,'POASONOV BROJ ............................... V =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =     1  (LINEAR VISCOELASTIC)'
     2 ///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'MODEL TYPE .............................. MTYPE =',I5/
     316X,'  = 1 MAXWELL (EQ.0 - DEFAULT TO 1) '/
     416X,'  = 2 VOIGT'/
     516X,'  = 3 KELVIN (STANDARD LINEAR SOLID)'/       
     116X,'YOUNG*S MODULUS ............................ ER =',1PD12.5/
     116X,'COEFFICIENT OF VISCOUS RESISTANCE ........... B =',1PD12.5/
     116X,'YOUNG*S MODULUS ............................ E1 =',1PD12.5/
     116X,'POISSONS RATIO .............................. V =',1PD12.5)
 6011 FORMAT(//
     116X,'MODEL TYPE .............................. MTYPE =',I5/
     316X,'  = 1 MAXWELL (EQ.0 - DEFAULT TO 1) '/
     416X,'  = 2 VOIGT'/
     516X,'  = 3 KELVIN (STANDARD LINEAR SOLID)'/       
     116X,'YOUNG*S MODULUS ............................ ER =',1PD12.5/
     116X,'COEFFICIENT OF VISCOUS RESISTANCE ........... B =',1PD12.5/
     116X,'POISSONS RATIO .............................. V =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD25(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C ......................................................................
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 25
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 25
C ......................................................................
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(30,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD25'
C ELASTICNE KONSTANTE
      IEL=0
      CALL CLEAR(FUNMAT(1,MAT),30)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3),IEL,(FUNMAT(J,MAT),J=12,14)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,3),IEL,(FUNMAT(J,MAT),J=12,14)
      FUNMAT(4,MAT)=IEL
      JOINT=0
      IF(FUNMAT(14,MAT).GT.0.D0) JOINT=1
C PLASTICNE KONSTANTE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=5,11)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=5,11)
C KONSTANTE OJACANJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=15,21)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=15,21)
      IHARD=FUNMAT(15,MAT)
      IF(IHARD.EQ.2)
     1  FUNMAT(22,MAT)=FUNMAT(16,MAT)*DEXP(-FUNMAT(17,MAT))*
     1                 DEXP(FUNMAT(19,MAT)*DEXP(-FUNMAT(20,MAT)))
C KONSTANTE OSTECENJA
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=27,30)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=27,30)
C
      IF(JOINT.EQ.1)THEN
        FUNMAT(3,MAT)=0.D0
        FUNMAT(6,MAT)=0.D0
      ENDIF
c
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
C
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(JOINT.EQ.0)THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,3),IEL
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,3),IEL
      ELSE
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2015) (FUNMAT(J,MAT),J=1,3),IEL,
     2                    (FUNMAT(J,MAT),J=12,14)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6015) (FUNMAT(J,MAT),J=1,3),IEL,
     2                    (FUNMAT(J,MAT),J=12,14)
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2043) (FUNMAT(J,MAT),J=5,7)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6043) (FUNMAT(J,MAT),J=5,7)
C
      IF(IHARD.LE.1)THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2044) IHARD,(FUNMAT(J,MAT),J=16,19)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6044) IHARD,(FUNMAT(J,MAT),J=16,19)
      ELSE
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2045) IHARD,(FUNMAT(J,MAT),J=16,22)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6045) IHARD,(FUNMAT(J,MAT),J=16,22)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2055) (FUNMAT(J,MAT),J=27,30)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6055) (FUNMAT(J,MAT),J=27,30)
      RETURN
 1000 FORMAT(7F10.0)
 1010 FORMAT(3F10.0,I5,3F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,    
     1'MODEL MATERIJALA BROJ =    25  (STENA-CST / PRSLINA)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     2  16X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//
     2  16X,'G R A N I C A   Z A T E Z A NJ A ...... AT =',1PD12.5//
     2  16X,'IDIKATOR ZADAVANJA MODULA ELASTIC .....IEL =',I12//)
 2015 FORMAT(//      
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'NORMALNI MODUL ELASTICNOSTI .......... EKN =',1PD12.5//
     2  16X,'TANGENCIJALNI MODUL ELASTICNOSTI ..... EKT =',1PD12.5//
     2  16X,'GRANICA ZATEZANJA ....................  AT =',1PD12.5//
     2  16X,'IDIKATOR ZADAVANJA MODULA ELASTIC .... IEL =',I12//
     2  16X,'VERTIKALNI AZIMUT PRSLINA ............ FIV =',1PD12.5//
     2  16X,'HORIZONTALNI AZIMUT PRSLINA .......... FIH =',1PD12.5//
     2  16X,'ZAPREMINSKI UDEO PRSLINA .............   F =',1PD12.5//)
 2043 FORMAT(//16X,'K O N S T A N T E    P L A S T I C N O S T I'//
     1  16X,'                                      GAMA =',1PD12.5//
     1  16X,'                                      BETA =',1PD12.5//
     1  16X,'                                        EN =',1PD12.5//)
 2044 FORMAT(//16X,'K O N S T A N T E    O J A C A N J A'//
     1  16X,' (0=CST, 1=DESAI, 2=DIVAC)           IHARD =',I5//
     1  16X,'                                      ALAM =',1PD12.5//
     1  16X,'                                       AKA =',1PD12.5//
     1  16X,'                                        P0 =',1PD12.5// 
     1  16X,'                                      AEE0 =',1PD12.5//)
 2045 FORMAT(//16X,'K O N S T A N T E    O J A C A N J A'//
     1  16X,' (0=CST, 1=DESAI, 2=DIVAC)           IHARD =',I5//
     1  16X,'                                         A =',1PD12.5//
     1  16X,'                                     AKPSI =',1PD12.5//
     1  16X,'                                    AKPSI1 =',1PD12.5// 
     1  16X,'                                       AKN =',1PD12.5// 
     1  16X,'                                      AKN1 =',1PD12.5// 
     1  16X,'                                     USTAR =',1PD12.5//
     1  16X,'                                     ALFAR =',1PD12.5//)
 2055 FORMAT(//16X,'K O N S T A N T E    O S T E C E N J A'/
     1         16X,'    r=R0+RU*(1-EXP(-AK*(EPSPD)**R)'/
     1  16X,'EPSPD = EKVIVALENTNA DEVIJATORSKA PLASTICNA DEFORMACIJA'//
     1  16X,'INICIJALNI NIVO OSTECENJA ............  R0 =',1PD12.5//
     1  16X,'GRANICNI PRIRAST OSTECENJA ...........  RU =',1PD12.5//
     2  16X,'KOEFICIJENT MNOZIOC  .................  AK =',1PD12.5//
     2  16X,'KOEFICIJENT EKSPONENT  ...............   R =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    25  (ROCKS-CST / JOINT)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     1  16X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     2  16X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//
     2  16X,'T E N S I O N      C U T O F F  ....... AT =',1PD12.5//
     2  16X,'I N D I C A T O R  OF  C H A N G E  E .IEL =',I12//)
 6015 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     1  16X,'NORMAL ELASTIC MODULUS ............... EKN =',1PD12.5//
     2  16X,'TANGENTIAL ELASTIC MODULUS ........... EKT =',1PD12.5//
     2  16X,'TENSION CUTOFF .......................  AT =',1PD12.5//
     2  16X,'INDICATOR OF CHANGE OF E ............. IEL =',I12//
     2  16X,'VERTICAL ASIMUT OF JOINTS ............ FIV =',1PD12.5//
     2  16X,'HORISONTAL ASIMUT OF JOINTS .......... FIH =',1PD12.5//
     2  16X,'VOLUMETRIC DENSITI OF JOINTS .........   F =',1PD12.5//)
 6043 FORMAT(//16X,'P L A S T I C I T Y    C O N S T A N T S'//
     1  16X,'                                      GAMA =',1PD12.5//
     1  16X,'                                      BETA =',1PD12.5//
     1  16X,'                                        EN =',1PD12.5//)
 6044 FORMAT(//16X,'H A R D E N I N G    C O N S T A N T S'//
     1  16X,' (0=CST, 1=DESAI, 2=DIVAC)           IHARD =',I5//
     1  16X,'                                      ALAM =',1PD12.5//
     1  16X,'                                       AKA =',1PD12.5//
     1  16X,'                                        P0 =',1PD12.5// 
     1  16X,'                                      AEE0 =',1PD12.5//)
 6045 FORMAT(//16X,'H A R D E N I N G    C O N S T A N T S'//
     1  16X,' (0=CST, 1=DESAI, 2=DIVAC)           IHARD =',I5//
     1  16X,'                                         A =',1PD12.5//
     1  16X,'                                     AKPSI =',1PD12.5//
     1  16X,'                                    AKPSI1 =',1PD12.5// 
     1  16X,'                                       AKN =',1PD12.5// 
     1  16X,'                                      AKN1 =',1PD12.5// 
     1  16X,'                                     USTAR =',1PD12.5//
     1  16X,'                                     ALFAR =',1PD12.5//)
 6055 FORMAT(//16X,'D A M A G E   C O N S T A N T S'/
     1         16X,'    r=R0+RU*(1-EXP(-AK*(EPSPD)**R)'/
     1  16X,'EPSPD = EQIVALENT DEVIATORIC PLASTIC STRAIN'//
     1  16X,'INITIAL DAMAGE  ......................  R0 =',1PD12.5//
     1  16X,'LIMITING DAMAGE INCRESE  .............  RU =',1PD12.5//
     2  16X,'COEFICIENT MULTIPLAYER  ..............  AK =',1PD12.5//
     2  16X,'COEFICIENT EXPONENT   ................   R =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD26(FUNMAT,NTFUN,MATE,MAXT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 26
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 26
C .
CE.      MAT         (MAT - MATERIAL NUMBER)
CS.      MAT         (MAT - MATERIJAL BROJ)
C .
CE.   FUNMAT(1,MAT,MAXT+1)   - YOUNGS  MODULUS - E0
CE.   FUNMAT(2,MAT,MAXT+1)   - POISSONS  NUMBER - V
CS.   FUNMAT(1,MAT,MAXT+1)   - MODUL ELASTICNOSTI - E0
CS.   FUNMAT(2,MAT,MAXT+1)   - POISSONOV BROJ - V
C .
CE.   NTFUN(MAT)   - NUMBER OF POINTS FOR ELASTIC MODULUS - Et
CS.   NTFUN(MAT)   - BROJ TACAKA ZA MODUL ELATICNOSTI - Et
C .
CE.   J=1,(NUMBER OF POINTS)
CS.   J=1,(DO BROJA TACAKA)
CE.   FUNMAT(1,MAT,J)   - YOUNGS  MODULUS - Ek
CE.   FUNMAT(2,MAT,J)   - CONTANT - Ak
CS.   FUNMAT(1,MAT,J)   - MODUL ELASTICNOSTI - Ek
CS.   FUNMAT(2,MAT,J)   - KONSTANTA - Ak
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(2,MATE,*),NTFUN(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD26'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT,MAXT+1),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT,MAXT+1),I=1,2)
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) NTFUN(MAT)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) NTFUN(MAT)                              
        IF(NTFUN(MAT).GT.MAXT.OR.NTFUN(MAT).LE.0) STOP 'UMOD31'
        DO 20 J=1,NTFUN(MAT)
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*) (FUNMAT(I,MAT,J),I=1,2)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1010) (FUNMAT(I,MAT,J),I=1,2)
   20   CONTINUE
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) (FUNMAT(I,MAT,MAXT+1),I=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) (FUNMAT(I,MAT,MAXT+1),I=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) NTFUN(MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) NTFUN(MAT)
      WRITE(IZLAZ,5020) ((FUNMAT(I,MAT,J),I=1,2),J=1,NTFUN(MAT))
      RETURN
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.0)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    31  (LINEARAN VISKO-ELASTICAN',
     2' IZOTROPAN)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  -  Et=E0*Bk*EXP(-Ak*t)'
     1//21X,'BROJ TACAKA - K =',I5//
     120X,'         Bk              Ak')
 2050 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  . E0 =',1PD12.5//
     116X,'P O I S S O N O V    B R O J ........... V =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    31  (LINEAR VISCO-ELASTIC',
     2' ISOTROPIC)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S      M O D U L U S  -  Et=E0*Bk*EXP(-Ak*t)'//
     121X,'NUMBER OF POINTS - K =',I5//
     120X,'         Ek              Ak')
 6050 FORMAT(//
     116X,'Y O U N G S      M O D U L U S ........ E0 =',1PD12.5//
     116X,'P O I S S O N S      N U M B E R ....... V =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD27(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 27
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 27
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CYLINDER COMPRESSIVE STRENGTH - FC
CS.   FUNMAT(1,MAT) - JACINA CILINDRA NA PRITISAK - FC
C .
CE.   FUNMAT(2,MAT) - SHEAR RETENTION FACTOR - SRF
CS.   FUNMAT(2,MAT) - STABILIZACIONI FAKTOR SMICANJA - SRF
C .
CE.   FUNMAT(3,MAT) - INITIAL YOUNG*S MODULUS - E
CS.   FUNMAT(3,MAT) - POCETNI MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(4,MAT) - INITIAL POISSON*S RATIO - V
CS.   FUNMAT(4,MAT) - POCETNI POISSONOV BROJ - V
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /CONFAC/ FACTOR
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(6,*),ZERO3(3),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
      DATA ZERO3/0.D0,0.D0,0.D0/
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD27'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT),FACTOR
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT),FACTOR
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(2,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(2,MAT)
C
      CALL CONCON(FUNMAT(1,MAT),ZERO3,AK0,G0,AKS0,GS0,E0,ANI0,
     1            FUNMAT(3,MAT),FUNMAT(4,MAT))
c      CALL CONCON(FUNMAT(1,MAT),ZERO3,AK0,G0,AKS0,GS0,E0,ANI0,
c     1            FUNMAT(3,MAT),FUNMAT(4,MAT),SIGID)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2),E0,ANI0,FACTOR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2),E0,ANI0,FACTOR
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    27   (BETON) '///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'JACINA CILINDRA NA PRITISAK ........... FC =',1PD12.5//
     116X,'STABILIZACIONI FAKTOR SMICANJA ....... SRF =',1PD12.5//
     116X,'M O D U L    E L A S T I C N O S T I .. E0 =',1PD12.5//
     216X,'P O I S S O N O V    B R O J .......... V0 =',1PD12.5//
     216X,'FAKTOR MERNIH JEDINICA ............ FACTOR =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    27  (CONCRETE)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'CYLINDER COMPRESSIVE STRENGTH ......... FC =',1PD12.5//
     116X,'SHEAR RETENTION FACTOR ............... SRF =',1PD12.5//
     116X,'Y O U N G S    M O D U L U S .......... E0 =',1PD12.5//
     216X,'P O I S S O N S    R A T I O .......... V0 =',1PD12.5//
     216X,'FAKCTOR OF MEASURED UNIT........... FACTOR =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD28(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 1
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 1
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS  - E
CS.   FUNMAT(1,MAT) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(2,MAT) - POISSON*S RATIO - V
CS.   FUNMAT(2,MAT) - POISSONOV BROJ - V
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(3,*),AMAT(30)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD28'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(2,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(2,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(3,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(3,MAT)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,3)
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    28  (GUMA)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M A T E R I J A L N A  K O N S T A N T A   C1 =',1PD12.5//
     216X,'M A T E R I J A L N A  K O N S T A N T A   C2 =',1PD12.5//
     316X,'P O I S S O N O V  B R O J                 V  =',1PD12.5//)
C--------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    28  (RUBBER)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'M A T E R I A L  C O N S T A N T S......C1 =',1PD12.5//
     216X,'M A T E R I A L  C O N S T A N T S......C2 =',1PD12.5//
     316X,'P O I S S O N S  R A T I O        ......V  =',1PD12.5//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD29(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 50
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 50
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS  - E
CS.   FUNMAT(1,MAT) - MODUL ELASTICNOSTI - E
C .
CE.   FUNMAT(2,MAT) - POISSON*S RATIO - V
CS.   FUNMAT(2,MAT) - POISSONOV BROJ - V
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(10,*)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ANAND / AB,AKSI,AM,SL0,H0,AL,SLHET,AN
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD29'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(2,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(2,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
C     1READ(IULAZ,*) AB,AKSI,AM,SL0
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,6)
      IF(INDFOR.EQ.2)
C     1READ(ACOZ,1000) AB,AKSI,AM,SL0
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
C     1READ(IULAZ,*) H0,AL,SLHET,AN
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,10)
      IF(INDFOR.EQ.2)
C     1READ(ACOZ,1000) H0,AL,SLHET,AN
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,10)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0) WRITE(IZLAZ,2010)
C     1 (FUNMAT(J,MAT),J=1,2),AB,AKSI,AM,SL0,H0,AL,SLHET,AN
     1 (FUNMAT(J,MAT),J=1,10)
      IF(ISRPS.EQ.1) WRITE(IZLAZ,6010)
C     1 (FUNMAT(J,MAT),J=1,2),AB,AKSI,AM,SL0,H0,AL,SLHET,AN
     1 (FUNMAT(J,MAT),J=1,10)
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    29  (ANAND-OV PLASTICAN IZOTROPAN)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I ... E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J ........... V =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ............ A =',1PD12.5//
     216X,'M A T.    K O N S T A N T A .......... KSI =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ............ m =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ........... s0 =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ........... h0 =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ............ a =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ....... sTILDA =',1PD12.5//
     216X,'M A T.    K O N S T A N T A ............ n =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    29  (ANAND`S PLASTIC ISOTROPIC)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S    M O D U L U S ........... E =',1PD12.5//
     216X,'P O I S S O N S    R A T I O ........... V =',1PD12.5//
     216X,'M A T.    C O N S T A N T .............. A =',1PD12.5//
     216X,'M A T.    C O N S T A N T ............ KSI =',1PD12.5//
     216X,'M A T.    C O N S T A N T .............. m =',1PD12.5//
     216X,'M A T.    C O N S T A N T ............. s0 =',1PD12.5//
     216X,'M A T.    C O N S T A N T ............. h0 =',1PD12.5//
     216X,'M A T.    C O N S T A N T .............. a =',1PD12.5//
     216X,'M A T.    C O N S T A N T ......... sTILDA =',1PD12.5//
     216X,'M A T.    C O N S T A N T .............. n =',1PD12.5)
C------------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD30(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 30
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 30
C          ELASTICAN SA UNUTRASNJIM TRENJEM I VISKOZNOSCU
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X (FIBER MODULUS) EX
CE.   FUNMAT(2,MAT) - YOUNG*S MODULUS IN DIRECTION Y  -  EY
CE.   FUNMAT(3,MAT) - YOUNG*S MODULUS IN DIRECTION Z  -  EZ
CS.   FUNMAT(1,MAT) - MODUL ELASTICNOSTI U PRAVCU X OSE (MOD. VLAKNA) EX
CS.   FUNMAT(2,MAT) - MODUL ELASTICNOSTI U PRAVCU Y OSE - EY
CS.   FUNMAT(3,MAT) - MODUL ELASTICNOSTI U PRAVCU Z OSE - EZ
C .
CE.   FUNMAT(4,MAT) - POISSON*S RATIO FOR PLANE XY    -  VXY
CE.   FUNMAT(5,MAT) - POISSON*S RATIO FOR PLANE YZ    -  VYZ
CE.   FUNMAT(6,MAT) - POISSON*S RATIO FOR PLANE ZX    -  VZX
CS.   FUNMAT(4,MAT) - POISSONOV BROJ U PRAVCU X I Y OSE - VXY
CS.   FUNMAT(5,MAT) - POISSONOV BROJ U PRAVCU Y I Z OSE - VYZ
CS.   FUNMAT(6,MAT) - POISSONOV BROJ U PRAVCU Z I X OSE - VZX
C .
CE.   FUNMAT(7,MAT) - SHEAR MODULUS FOR PLANE XY  -  GXY (FOR FIBERS)
CE.   FUNMAT(8,MAT) - SHEAR MODULUS FOR PLANE YZ  -  GYZ
CE.   FUNMAT(9,MAT) - SHEAR MODULUS FOR PLANE ZX  -  GZX  
CS.   FUNMAT(7,MAT) - MODUL SMICANJA U PRAVCU X I Y OSE - GXY (VLAKNA)
CS.   FUNMAT(8,MAT) - MODUL SMICANJA U PRAVCU Y I Z OSE - GYZ
CS.   FUNMAT(9,MAT) - MODUL SMICANJA U PRAVCU Z I X OSE - GZX
C .
C .
CS.   FUNMAT(10,MAT) - ADHESIONI PRITISAK 
CE.   FUNMAT(10,MAT) - ADHESION PRESSURE
C
CS.   FUNMAT(11,MAT) - KOEFICIENT KULONOVOG TRENJA
CE.   FUNMAT(11,MAT) - COEEFICIENT OF COULOMB FRICTION
C .
CS.   FUNMAT(12,MAT) - KOEFICIENT VISKOZNOG OTPORA
CE.   FUNMAT(12,MAT) - VISCOUS COEEFICIENT
C .
CS.   FUNMAT(13,MAT) - DUZINA VLAKNA
CE.   FUNMAT(13,MAT) - FIBER LENGTH
C
CS.   FUNMAT(14,MAT) - DEBLJINA VLAKNA
CE.   FUNMAT(14,MAT) - FIBER THICKNESS
C
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(14,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD30'
       CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,6)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=7,9)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=7,9)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=10,14)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=10,14)
C
C... CHECK MATERIAL CONSTANTS
      CALL       ANICHK(FUNMAT(1,MAT),IZLAZ,ISRPS)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,9)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,9)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2011) (FUNMAT(J,MAT),J=10,14)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6011) (FUNMAT(J,MAT),J=10,14)
      RETURN
C
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    30  (ELASTICAN SA UNUTRASNJIM',
     1 ' TRENJEM) '///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L I    E L A S T I C N O S T I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N O V I    B R O J E V I'/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'M O D U L I    S M I C A N J A'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X))
 2011 FORMAT(//
     116X,'ADHESIONI PRITISAK ...................... P0 =',1PD12.5//
     216X,'KOEFICIENT KULONOVOG TRENJA ............ AMI =',1PD12.5//
     216X,'KOEFICIENT VISKOZNOG OTPORA ............ CNI =',1PD12.5//
     216X,'DUZINA VLAKNA ........................... EL =',1PD12.5//
     216X,'DEBLJINA VLAKNA ......................... AT =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    30  (ELASTIC MATERIAL WITH INTERNAL',
     1 ' FRICTION) '///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U I'/21X,'EX',12X,'EY',
     112X,'EZ'/16X,3(1PD12.5,2X)//
     1//16X,'P O I S S O N S      R A T I O S '/20X,'VXY',11X,'VYZ',
     111X,'VZX'/16X,3(1PD12.5,2X)//
     1//16X,'S H E A R    M O D U L U I'/20X,'GXY',11X,'GYZ',
     111X,'GZX'/16X,3(1PD12.5,2X))
 6011 FORMAT(//
     116X,'ADHESION  PRESSURE ...................... P0 =',1PD12.5//
     216X,'COEFFICIENT OF COULOMB FRICTION ........ AMI =',1PD12.5//
     216X,'COEFFICIENT OF VISCOUS RESISTANCE ...... CNI =',1PD12.5//
     216X,'FIBERS LENGTH ........................... EL =',1PD12.5//
     216X,'FIBERS THICKNESS ........................ AT =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD31(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 31
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 31
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - E
CE.   FUNMAT(2,MAT) - CONSTANT  - ANI
CE.   FUNMAT(3,MAT) - CONSTANT  - ALFA
CE.   FUNMAT(4,MAT) - CONSTANT  - BETA
CE.   FUNMAT(5,MAT) - CONSTANT  - A
CE.   FUNMAT(6,MAT) - CONSTANT  - DEFMP
CE.   FUNMAT(7,MAT) - CONSTANT  - AKA
CS.   FUNMAT(8,MAT) - KONSTANTA - FI
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CDEBUG/ IDEBUG
C
      COMMON /MISIC/ NTS,NTALF,DALM(2,2,30)
      DIMENSION FUNMAT(8,*),NTAV(2)
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD31 '
      TMAX = VREME
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,4)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=5,8)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=5,8)
C
      DO 100 N=1,2
      CALL ISPITA(ACOZ)
      KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) IBR,IMAX
      IF(INDFOR.EQ.2)
     1READ(ACOZ,3000) IBR,IMAX
C
      NTAV(IBR)=IMAX
      DO 20 J=1,IMAX
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (DALM(I,IBR,J),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,3010) (DALM(I,IBR,J),I=1,2)
   20 CONTINUE
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      IF(N.EQ.1)
     1CALL WBROJK(KARTI,0)
      IF(IBR.EQ.1) THEN 
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,4)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,4)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=5,8)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=5,8)
      ENDIF
C
      IF(IBR.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,4001) IBR,NTAV(IBR)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,7001) IBR,NTAV(IBR)
      ENDIF
      IF(IBR.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,4011)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,7011)
      GO TO 150
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,4000) IBR,NTAV(IBR)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,7000) IBR,NTAV(IBR)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,4010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,7010)
  150 CONTINUE
      IF(ISRPS.EQ.0.OR.ISRPS.EQ.1)
     1WRITE(IZLAZ,5020) ((DALM(I,IBR,J),I=1,2),J=1,IMAX)
C
  100 CONTINUE
      NTS = NTAV(1)
      NTALF = NTAV(2)
      IF(DALM(1,2,IMAX).LT.TMAX) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,1510)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,2510)
      STOP
      ENDIF
C     NTM = MAX0(NTS,NTALF)
      RETURN
C
 1000 FORMAT(4F10.0)
 3000 FORMAT(14I5)
 3010 FORMAT(2F10.0)
 5020 FORMAT(9X,1PD12.5,2X,1PD12.5)
C-----------------------------------------------------------------------
 1510 FORMAT(2X,'FUNKCIJA AKTIVACIJE NIJE DEFINISANA U CELOM INTERVALU')
 2510 FORMAT(2X,'FUNCTION ACTIVATION IS NOT DEFINED IN ALL INTERVALS  ')
C-----------------------------------------------------------------------
 4000 FORMAT(//6X,'NAPON SIG0          =',I5,5X,'BROJ TACAKA  =',I5)
 4001 FORMAT(//6X,'FUNKCIJA AKTIVACIJE =',I5,5X,'BROJ TACAKA  =',I5)
 4010 FORMAT(/12X,'STRECH',6X,'SIG0')
 4011 FORMAT(/12X,'VREME',6X,'FUNKCIJA AKTIVACIJE')
C-----------------------------------------------------------------------
 7000 FORMAT(//6X,'STRESS SIG0       =',I5,5X,'NUMBER POINTS  =',I5)
 7001 FORMAT(//6X,'FUNCTION ACTIVATION =',I5,5X,'NUMBER POINTS=',I5)
 7010 FORMAT(/12X,'STRECH',6X,'SIG0')
 7011 FORMAT(/12X,'TIME ',6X,'FUNCTION ACTIVATION')
C-----------------------------------------------------------------------
 2000 FORMAT(6X,  
     1'MODEL MATERIJALA BROJ =    26  (  MISIC   )'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//
     216X,'ELASTICNA KONSTANTA TENDONA ...... ALFA=',1PD12.5//
     216X,'ELASTICNA KONSTANTA TENDONA ...... BETA=',1PD12.5)
 2045 FORMAT(//6X,
     2/6X,'A',10X,'DEFMP',9X,'AKA',10X,'FI  '/
     3' ',4(1PD12.4))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    26  (    MUSCLE    )'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//
     216X,'ELASTIC CONSTANT TENDON ...........ALFA=',1PD12.5//
     216X,'ELASTIC CONSTANT TENDON .......... BETA=',1PD12.5)
 6045 FORMAT(//6X,
     2/6X,'A',10X,'DEFMP',9X,'AKA',10X,'FI  '/
     3' ',4(1PD12.4))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD32(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 32
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 32
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - ANI
CE.   FUNMAT(2,MAT) - CONSTANT  - AYP
CE.   FUNMAT(3,MAT) - CONSTANT  - BYP
CE.   FUNMAT(4,MAT) - CONSTANT  - AXP 
CE.   FUNMAT(5,MAT) - CONSTANT  - BXP 
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION FUNMAT(5,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD32 '
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=2,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=2,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,5)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) (FUNMAT(J,MAT),J=4,5)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) (FUNMAT(J,MAT),J=4,5)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(7I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    32  (NAPON-STREC MODEL ZA PAS',
     2'IVNO STANJE)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//16X,'POISSONOV BROJ .................... V =',1PD12.5) 
 2020 FORMAT(//
     111X,'KOEFICIJENTI ZA UZDUZNI PRAVAC                        '//
     216X,'PRVI KOEFICIJENT .............. ........ AYP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BYP =',1PD12.5)
 2030 FORMAT(//
     111X,'KOEFICIJENTI ZA POPRECNI PRAVAC                       '//
     216X,'PRVI KOEFICIJENT .............. ........ AXP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BXP =',1PD12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    32  (STRESS-STRETCH MODEL',
     2' FOR PASSIVE STATE)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//16X,'POISSONS MODULUS .................. V =',1PD12.5) 
 6020 FORMAT(//
     111X,'COEFFICIENTS FOR LONGITUDINAL DIRECTION               '//
     216X,'FIRST COEFFICIENT.............. ........ AYP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BYP =',1PD12.5)
 6030 FORMAT(//
     111X,'COEFFICIENTS FOR CIRCUMFERENTIAL DIRECTION            '//
     216X,'FIRST COEFFICIENT.............. ........ AXP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BXP =',1PD12.5/)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD33(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 33
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 33
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - ANI
CE.   FUNMAT(2,MAT) - CONSTANT  - AYP
CE.   FUNMAT(3,MAT) - CONSTANT  - BYP
CE.   FUNMAT(4,MAT) - CONSTANT  - AXP 
CE.   FUNMAT(5,MAT) - CONSTANT  - BXP 
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CDEBUG/ IDEBUG
C
      COMMON /RELAKS/ DLAMDC(10,4),ET(2),NLAM,INDREL
C
      DIMENSION FUNMAT(5,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD33 '
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=2,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=2,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,5)
C
C ZBOG RELAKSACIJE !!!
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NLAM,INDREL
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) NLAM,INDREL
      IF(NLAM.EQ.0) INDREL=100
C        IF(NTFUN(MAT).GT.MAXT.OR.NTFUN(MAT).LE.0) STOP 'UMOD33'
      DO 20 I=1,NLAM
        CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (DLAMDC(I,J),J=1,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (DLAMDC(I,J),J=1,4)
   20 CONTINUE
C
C ZBOG RELAKSACIJE !!!
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) (FUNMAT(J,MAT),J=4,5)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) (FUNMAT(J,MAT),J=4,5)
C
C ZBOG RELAKSACIJE !!!
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NLAM,INDREL
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NLAM,INDREL
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,5020) ((DLAMDC(I,J),J=1,4),I=1,NLAM)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,5020) ((DLAMDC(I,J),J=1,4),I=1,NLAM)
C
C ZBOG RELAKSACIJE !!!
C
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(7I5)
 5020 FORMAT(5X,1PD12.5,3X,1PD12.5,3X,1PD12.5,3X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    33  (NAPON-STREC MODEL PUZANJ',
     2'A)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//16X,'POISSONOV BROJ .................... V =',1PD12.5) 
 2020 FORMAT(//
     111X,'KOEFICIJENTI ZA UZDUZNI PRAVAC                        '//
     216X,'PRVI KOEFICIJENT .............. ........ AYP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BYP =',1PD12.5)
 2030 FORMAT(//
     111X,'KOEFICIJENTI ZA POPRECNI PRAVAC                       '//
     216X,'PRVI KOEFICIJENT .............. ........ AXP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BXP =',1PD12.5)
 2040 FORMAT(/
     116X,'K R I V E    N A P O N  -  V R E M E'//
     16X,'BROJ KRIVIH RELAKSACIJE =',I5,6X,'KORAK RELAKSACIJE =',I5//
     18X,'STREC            C0             C1             C2  ')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    33  (STRESS-STRETCH CREEP',
     2' MODEL)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//16X,'POISSONS MODULUS .................. V =',1PD12.5) 
 6020 FORMAT(//
     111X,'COEFFICIENTS FOR LONGITUDINAL DIRECTION               '//
     216X,'FIRST COEFFICIENT.............. ........ AYP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BYP =',1PD12.5)
 6030 FORMAT(//
     111X,'COEFFICIENTS FOR CIRCUMFERENTIAL DIRECTION            '//
     216X,'FIRST COEFFICIENT.............. ........ AXP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BXP =',1PD12.5/)
 6040 FORMAT(/
     120X,'S T R E S S - T I M E    C U R V E S'//
     16X,'NUMBER OF RELAXATION CURVES =',I5,6X,'STEP OF RELAXATION =',I5
     1//8X,'STRETCH          C0             C1             C2  ')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD34(FUNMAT,MAT,KARTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 34
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 34
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - CONSTANT  - ANI
CE.   FUNMAT(2,MAT) - CONSTANT  - AYP
CE.   FUNMAT(3,MAT) - CONSTANT  - BYP
CE.   FUNMAT(4,MAT) - CONSTANT  - AXP 
CE.   FUNMAT(5,MAT) - CONSTANT  - BXP 
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CDEBUG/ IDEBUG
C
      COMMON /AKTIV/ COEFIC(10,3),FAKTIV(10,2),NAKTIV,IAKTIV,IMAX
C
      DIMENSION FUNMAT(5,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD34 '
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) FUNMAT(1,MAT)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) FUNMAT(1,MAT)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=2,3)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=2,3)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=4,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=4,5)
C
C ZBOG AKTIVNOG STANJA !!!
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NAKTIV,IAKTIV
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) NAKTIV,IAKTIV
      IF(NAKTIV.EQ.0) IAKTIV=100
C        IF(NTFUN(MAT).GT.MAXT.OR.NTFUN(MAT).LE.0) STOP 'UMOD34'
      DO 30 I=1,NAKTIV
       CALL ISPITA(ACOZ)
       IF(INDFOR.EQ.1)
     1 READ(IULAZ,*) (COEFIC(I,J),J=1,3)
       IF(INDFOR.EQ.2)
     1 READ(ACOZ,1000) (COEFIC(I,J),J=1,3)
   30 CONTINUE
C
C FUNKCIJA AKTIVACIJE
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) IMAX
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) IMAX
C
      DO 40 I=1,IMAX
       CALL ISPITA(ACOZ)
       IF(INDFOR.EQ.1)
     1 READ(IULAZ,*) (FAKTIV(I,J),J=1,2)
       IF(INDFOR.EQ.2)
     1 READ(ACOZ,1000) (FAKTIV(I,J),J=1,2)
   40 CONTINUE
C
C ZBOG AKTIVNOG STANJA !!!
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) FUNMAT(1,MAT)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) (FUNMAT(J,MAT),J=2,3)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) (FUNMAT(J,MAT),J=4,5)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) (FUNMAT(J,MAT),J=4,5)
C
C ZBOG AKTIVNOG STANJA !!!
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) NAKTIV,IAKTIV
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) NAKTIV,IAKTIV
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,5030) ((COEFIC(I,J),J=1,3),I=1,NAKTIV)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,5030) ((COEFIC(I,J),J=1,3),I=1,NAKTIV)
C
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060) IMAX
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060) IMAX
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,5040) ((FAKTIV(I,J),J=1,2),I=1,IMAX)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,5040) ((FAKTIV(I,J),J=1,2),I=1,IMAX)
C
C ZBOG AKTIVNOG STANJA !!!
C
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(7I5)
 5030 FORMAT(5X,1PD12.5,3X,1PD12.5,3X,1PD12.5)
 5040 FORMAT(5X,1PD12.5,3X,1PD12.5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    34  (NAPON-STREC MODEL ZA AKT',
     2'IVNO STANJE)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//16X,'POISSONOV BROJ .................... V =',1PD12.5) 
 2020 FORMAT(//
     111X,'KOEFICIJENTI ZA UZDUZNI PRAVAC                        '//
     216X,'PRVI KOEFICIJENT .............. ........ AYP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BYP =',1PD12.5)
 2030 FORMAT(//
     111X,'KOEFICIJENTI ZA POPRECNI PRAVAC                       '//
     216X,'PRVI KOEFICIJENT .............. ........ AXP =',1PD12.5//
     316X,'DRUGI KOEFICIJENT ...................... BXP =',1PD12.5)
 2050 FORMAT(/
     16X,'K R I V E    N A P O N - S T R E C - B R Z I N A   D E F. '//
     16X,'BROJ KRIVIH KONTRAKCIJE =',I5,6X,'KORAK KONTRAKCIJE =',I5//
     12X,'BRZINA DEFORMACIJE     AY2A          BY2A ')
 2060 FORMAT(/
     116X,'F U N K C I J A  A K T I V A C I J E '//
     125X,'BROJ TACAKA =',I5//
     13X,'VREME     FUN. AKTIVACIJE')
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    34  (STRESS-STRETCH MODEL',
     2' FOR ACTIVE STATE)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//16X,'POISSONS MODULUS .................. V =',1PD12.5) 
 6020 FORMAT(//
     111X,'COEFFICIENTS FOR LONGITUDINAL DIRECTION               '//
     216X,'FIRST COEFFICIENT.............. ........ AYP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BYP =',1PD12.5)
 6030 FORMAT(//
     111X,'COEFFICIENTS FOR CIRCUMFERENTIAL DIRECTION            '//
     216X,'FIRST COEFFICIENT.............. ........ AXP =',1PD12.5//
     216X,'SECOND COEFFICIENT...................... BXP =',1PD12.5/)
 6050 FORMAT(//
     18X,'S T R E S S - S T R E T C H - S T R A I N  R A T E  C U R.'//
     16X,'NUMBER OF CONTRACTION CURVES=',I5,6X,'CONTRACTION STEP =',I5//
     15X,' STRAIN RATE        AY2A           BY2A ')
 6060 FORMAT(//
     122X,'F U N C T I O N  A C T I V A T I O N '//
     129X,'NUMBER OF POINTS =',I5//
     19X,'TIME      FUN. ACTIVATION')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD40(FUNMAT,NTFUN,MATE,MAXT,MOD,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 40
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 40
C .
CE.               MAT - MATERIAL NUMBER
CE.   J=2,(NTMAX - NUMBER OF POINTS FOR CURVE) + 1
CS.               MAT - MATERIJAL BROJ
CS.   J=2,(NTMAX - BROJ TACAKA ZA KRIVU) + 1
C .
CE.   FUNMAT(1,MAT,1) - YOUNG*S MODULUS - E
CE.   FUNMAT(2,MAT,1) - POISSON*S NUMBER - V
CS.   FUNMAT(1,MAT,1) - MODUL ELASTICNOSTI - E
CS.   FUNMAT(2,MAT,1) - POISSONOV BROJ - V
C .
CE. A) NTMAX=1 (NTMAX - NUMBER OF POINTS FOR CURVE TAU(DEF))
CS. A) NTMAX=1 (NTMAX - BROJ TACAKA ZA KRIVU TAU(DEF))
C .
CE.   FUNMAT(1,MAT,2) - YIELD STRESS - TAUY
CE.   FUNMAT(2,MAT,2) - TANGENTIAL MODULUS - ET
CS.   FUNMAT(1,MAT,2) - NAPON TECENJA - TAUY
CS.   FUNMAT(2,MAT,2) - TANGENTNI MODUL - ET
C .
CE. B) NTMAX>1 (NTMAX - NUMBER OF POINTS FOR CURVE TAU(DEF))
CS. B) NTMAX>1 (NTMAX - BROJ TACAKA ZA KRIVU TAU(DEF))
C .
CE.   FUNMAT(1,MAT,J) - ARGUMENT - STRAIN - DEF
CE.   FUNMAT(2,MAT,J) - STRESS - TAU
CS.   FUNMAT(1,MAT,J) - ARGUMENT - DEFORMACIJA - DEF
CS.   FUNMAT(2,MAT,J) - NAPON - TAU
C .
CE.        NTFUN(MAT) - NUMBER OF POINTS FOR CURVE - TAU(DEF)
CE.                     =1; BILINEAR CURVE
CE.                     >1; MULTILINEAR CURVE
CS.        NTFUN(MAT) - BROJ TACAKA ZA KRIVU - TAU(DEF)
CS.                     =1; BILINEARNA KRIVA
CS.                     >1; MULTILINEARNA KRIVA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION FUNMAT(2,MATE,*),NTFUN(*),AMAT(30)
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      
      IF(IDEBUG.GT.0) PRINT *, ' UMOD40'
      ISIMO=0
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NTFUN(MAT),INDEP
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NTFUN(MAT),INDEP
      NBRT=NTFUN(MAT)
      IF (NBRT.GT.1 .AND. INDEP.EQ.0 ) STOP 'INDEP - PAK03'
      IF (NBRT.GT.MAXT) STOP 'MAXT - PAK03'
      IPROM=10
      IF (NBRT.GT.1) IPROM=2+NBRT
      DO 20 J=1,2
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (FUNMAT(I,MAT,J),I=1,2),DUM
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1010) (FUNMAT(I,MAT,J),I=1,2),DUM
   20 CONTINUE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT,3),I=1,2),FUNMAT(2,MAT,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT,3),I=1,2),FUNMAT(2,MAT,4)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(2,MAT,J),J=5,10)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(2,MAT,J),J=5,10)

      IF(DABS(FUNMAT(1,MAT,3)).LT.1.D-6) FUNMAT(1,MAT,3)=1.D0
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) GO TO 90
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(I,MAT,1),I=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(I,MAT,1),I=1,2)
C
      IF(NTFUN(MAT).EQ.1) THEN
      IF(ISIMO.EQ.0) THEN
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2045) ((FUNMAT(I,MAT,J),I=1,2),J=2,3),
     1   FUNMAT(2,MAT,4)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6045) ((FUNMAT(I,MAT,J),I=1,2),J=2,3)
        IF(ISRPS.EQ.0)
     1  WRITE(IZLAZ,2145) (FUNMAT(2,MAT,J),J=5,10)
        IF(ISRPS.EQ.1)
     1  WRITE(IZLAZ,6145) (FUNMAT(2,MAT,J),J=5,10)
C
      ELSEIF(ISIMO.EQ.1) THEN
        IF(ISRPS.EQ.0)THEN
          WRITE(IZLAZ,2075) (FUNMAT(I,MAT,2),I=1,2),FUNMAT(1,MAT,4),
     1                      (FUNMAT(I,MAT,3),I=1,2)
        ELSEIF(ISRPS.EQ.1)THEN
          WRITE(IZLAZ,6075) (FUNMAT(I,MAT,2),I=1,2),FUNMAT(1,MAT,4),
     1                      (FUNMAT(I,MAT,3),I=1,2)
        ENDIF
        FUNMAT(1,MAT,2)=FUNMAT(2,MAT,2)-FUNMAT(1,MAT,2)
      ENDIF
      ENDIF
C      
   90 IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT,1)
         AMAT(7)=FUNMAT(2,MAT,1)
         AMAT(13)=GUST
         AMAT(15)=FUNMAT(1,MAT,2)
         AMAT(16)=FUNMAT(2,MAT,2)
         AMAT(17)=FUNMAT(2,MAT,3)
         AMAT(19)=3.
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT(AMAT,MATG,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NTFUN(MAT).EQ.1) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030) NTFUN(MAT),INDEP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) NTFUN(MAT),INDEP
      WRITE(IZLAZ,5020) ((FUNMAT(I,MAT,J),I=1,2),J=4,IPROM)
C
CE    TRANSFORMATION OF CURVE ON FORM STRESS-PLASTIC STRAIN
CS    PREVODJENJE KRIVE NA OBLIK NAPON-PLASTICNA DEFORMACIJA
C
      CALL NAPPLD(FUNMAT,NTFUN,MATE,MAT)
      RETURN
C
 1000 FORMAT(14I5)
 1010 FORMAT(7F10.3)
 5020 FORMAT(19X,1PD12.5,10X,1PD12.5)
C-----------------------------------------------------------------------
 2005 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    40'/6X,'(GURSON MODEL',
     2' SA MESOVITIM OJACANJEM)'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     116X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5)
 2040 FORMAT(//
     116X,'N A P O N    T E C E N J A  ........  TAUY =',1PD12.5//
     216X,'T A N G E N T N I    M O D U L  ......  ET =',1PD12.5)
 2045 FORMAT(//6X,
     1'   RAMBERG-OSGOOD FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'N A P O N    T E C E N J A  ........  TAUY =',1PD12.5//
     216X,'M N O Z I O C  .....................    CY =',1PD12.5//
     316X,'E K S P O N E N T  .................    AN =',1PD12.5//
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ....    EM =',1PD12.5)
 2145 FORMAT(//6X,
     1'   MATERIALNE KONSTANTE - GURSON '//
     116X,' F0 =',1PD12.5/
c     216X,' K  =',1PD12.5/
     316X,' Ff =',1PD12.5/
     416X,' q1 =',1PD12.5/
     516X,' q2 =',1PD12.5/
     616X,' q3 =',1PD12.5/
     716X,' Fc =',1PD12.5)
 2030 FORMAT(//
     116X,'K R I V A    N A P O N  -  D E F O R M A C I J A'//
     121X,'BROJ TACAKA ZA KRIVU =',I5//
     121X,'VRSTA DEFORMACIJE , INDEP=',I5/
     126X,'.EQ.1; UKUPNA'/
     126X,'.EQ.2; PLASTICNA'//
     121X,'ARGUMENT             FUNKCIJA'/
     120X,'DEFORMACIJA             NAPON')
 2075 FORMAT(//6X,
     1'   J.C.SIMO FORMULA -  MATERIJALNE KONSTANTE'//
     116X,'NAPON TECENJA .....................   TAUY =',1PD12.5//
     216X,'ASIMPTOTSKI NAPON ................. TAUINF =',1PD12.5//
     216X,'KOEFICIJENT LINEARNOG OJACANJA ....      H =',1PD12.5//
     316X,'EKSPONENT .........................  DELTA =',1PD12.5///
     416X,'KOEFICIJENT MESOVITOG OJACANJA  ...     EM =',1PD12.5)
C-----------------------------------------------------------------------
 6005 FORMAT(6X,
     1'MATERIAL MODEL NUMBER  =    40'/6X,'(GURSON MODEL',
     2' WITH MIXED HARDENING)'///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5)
 6040 FORMAT(//
     116X,'Y I E L D    S T R E S S  ..........  TAUY =',1PD12.5//
     216X,'T A N G E N T I A L    M O D U L U S    ET =',1PD12.5)
 6045 FORMAT(//6X,
     1'   RAMBERG-OSGOOD FORMULAS MATERIAL CONSTANTS'//
     116X,'Y I E L D    S T R E S S  ..........  TAUY =',1PD12.5//
     216X,'M U L T I P L A Y E R  .............    CY =',1PD12.5//
     316X,'E X P O N E N T  ...................    AN =',1PD12.5///
     416X,'MIXED HARDENING COEFFICIENT  .......    EM =',1PD12.5)
 6145 FORMAT(//6X,
     1'   MATERIAL CONSTANTS - GURSON '//
     116X,' F0 =',1PD12.5/
c     216X,' K  =',1PD12.5/
     316X,' Ff =',1PD12.5/
     416X,' q1 =',1PD12.5/
     516X,' q2 =',1PD12.5/
     616X,' q3 =',1PD12.5/
     716X,' Fc =',1PD12.5)
 6030 FORMAT(//
     116X,'S T R E S S - S T R A I N    C U R V E'//
     121X,'NUMBER OF POINTS ON CURVE =',I5//
     121X,'DEFORMATION SHAPE , INDEP=',I5/
     126X,'.EQ.1; TOTAL STRAIN'/
     126X,'.EQ.2; PLASTIC STRAIN'//
     121X,'ARGUMENT             FUNCTION'/
     120X,'STRAIN                 STRESS')
 6075 FORMAT(//6X,
     1'   J.C.SIMO FORMULA FORMULAS MATERIAL CONSTANTS'//
     116X,'YIELD STRESS  .....................   TAUY =',1PD12.5//
     216X,'ASYMPTOTIC STRESS ................. TAUINF =',1PD12.5//
     216X,'LINEAR HARDENING COEFFICIENT ......      H =',1PD12.5//
     316X,'EXPONENT ..........................  DELTA =',1PD12.5///
     416X,'MIXED HARDENING COEFFICIENT  ......     EM =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE FICA(FUNMAT,MAT,MATE,IPROM,INDEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAUY0/ Y0
      DIMENSION X(50),Y(50),SIG(50),A(2),LISTA(2),COVAR(2,2)
      DIMENSION ALPHA(2,2)
      DIMENSION FUNMAT(2,MATE,*)
      IF ((IPROM-2).GT.50) STOP 'NBRT > 50'
      E=FUNMAT(1,MAT,1)
      Y0=FUNMAT(1,MAT,2)
      NBRTAC=IPROM-2
      DO 4 J=3,IPROM
      I=J-2
      X(I)=FUNMAT(1,MAT,J)
      Y(I)=FUNMAT(2,MAT,J)
   4  CONTINUE
      DO 5 I=1,NBRTAC
      SIG(I)=1.
      IF (INDEP.EQ.1)
     1X(I)=X(I)-Y(I)/E
 5    CONTINUE
      IF (DABS(X(NBRTAC)-X(1)).LT.1.D-10) STOP
      CY=(Y(NBRTAC)-Y(1))/(X(NBRTAC)-X(1))
      A(1)=CY
C     A(1)=1.
      A(2)=0.29447
        LISTA(1)=1
        LISTA(2)=2
      ALAMDA=-1.
14    CALL MRQMIN(X,Y,SIG,NBRTAC,A,2,LISTA,2,COVAR,ALPHA,2,CHISQ,
     1ALAMDA)
      IF (ALAMDA.GE.1.D13) STOP 'DIVERGIRA'
      IF (ABS(CHISQ).GE.1.D-7) GOTO 14
      ANULL=0.0
      CALL MRQMIN(X,Y,SIG,NBRTAC,A,2,LISTA,2,COVAR,ALPHA,2,CHISQ,ANULL)
      AM=FUNMAT(2,MAT,2)
      FUNMAT(2,MAT,2)=A(1)
      FUNMAT(1,MAT,3)=A(2)
      FUNMAT(2,MAT,3)=AM
      END
C============================================================================
      SUBROUTINE FUNCS(X,A,Y,DYDA,NA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAUY0/ Y0
      DIMENSION A(NA),DYDA(NA)
      IF(X.LT.1.D-20) THEN
      XA=0.
      ELSE
      XA=X**A(2)
      ENDIF
C     XA=DEXP(A(2)*DLOG(X))
      Y=Y0+A(1)*XA
      DYDA(1)=XA
C     IF(X.EQ.0.) THEN
      IF(X.LT.1.D-20) THEN
C     DYDA(2)=1.D20 
C     DYDA(2)=A(1)
      DYDA(2)=0.
      ELSE
      DYDA(2)=A(1)*XA*DLOG(X)
      ENDIF
      WRITE(2,102) A(1),A(2),X,Y,DYDA(1),DYDA(2)
  102 FORMAT ('A(1),A(2),X,Y,DYDA(1),DYDA(2)',6D10.3)
      RETURN
      END
C============================================================================
C============================================================================
C     SUBROUTINE FUNCS(X,A,Y,DYDA,NA)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DIMENSION A(NA),DYDA(NA)
C     Y=0.
C     DO 11 I=1,NA-1,3
C      ARG=(X-A(I+1))/A(I+2)
C      EX=EXP(-ARG**2)
C      FAC=A(I)*EX*2.*ARG
C      Y=Y+A(I)*EX
C      DYDA(I)=EX
C      DYDA(I+1)=FAC/A(I+2)
C      DYDA(I+2)=FAC*ARG/A(I+2)
C 11  CONTINUE
C     RETURN
C     END
C============================================================================
      SUBROUTINE MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *  COVAR,ALPHA,NCA,CHISQ,ALAMDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      PARAM=2.
      IF(ALAMDA.LT.0.)THEN
        KK=MFIT+1
        DO 12 J=1,MA
          IHIT=0
          DO 11 K=1,MFIT
            IF(LISTA(K).EQ.J)IHIT=IHIT+1
11        CONTINUE
          IF (IHIT.EQ.0) THEN
            LISTA(KK)=J
            KK=KK+1
          ELSE IF (IHIT.GT.1) THEN
            PAUSE 'Improper permutation in LISTA'
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) PAUSE 'Improper permutation in LISTA'
        ALAMDA=1.D-3
        CALL MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,CHISQ)
C       WRITE(*,*)'IZASAO IZ POTPROGRAMA MRQCOF PRVI PUT'
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
C     WRITE(*,*)'PRE POZIVA GAUSSJ'
      CALL GAUSSJ(COVAR,MFIT,NCA,DA,1,1)
C     WRITE(*,*)'POSLE POZIVA GAUSSJ'
      IF(DABS(ALAMDA).LT.1.D-10)THEN
        CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=ATRY(LISTA(J))+DA(J)
16    CONTINUE
      CALL MRQCOF(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ)
C       WRITE(*,*)'IZASAO IZ POTPROGRAMA MRQCOF DRUGI PUT'
      IF (OCHISQ-CHISQ.LT.1.D-3) THEN
      CHISQ=1.D-8
      RETURN
      ENDIF
      IF(CHISQ.LT.OCHISQ)THEN
        ALAMDA=1./PARAM*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MFIT
          DO 17 K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          A(LISTA(J))=ATRY(LISTA(J))
18      CONTINUE
      ELSE
        ALAMDA=PARAM*ALAMDA
        CHISQ=OCHISQ
      ENDIF
      RETURN
      END
C============================================================================
      SUBROUTINE MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,CH
     *ISQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT)
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        CALL FUNCS(X(I),A,YMOD,DYDA,MA)
        SIG2I=1./(SIG(I)*SIG(I))
        DY=Y(I)-YMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))*SIG2I
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY*SIG2I
15    CONTINUE
      DO 17 J=2,MFIT
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END
C============================================================================
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (DABS(A(ICOL,ICOL)).LT.1.D-10) PAUSE 'Singular matrix-PAK03'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
C============================================================================
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END
C============================================================================
      SUBROUTINE CITAC(USTR,IZRAZ,ACR,ISG,ITI,ITH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER*250 USTR,IZRAZ
      CHARACTER*20 BROJ
      CHARACTER*11 CIFRE
      CHARACTER*6  KCIF,OPER
      DIMENSION ACR(20)
C
      DATA CIFRE/'0123456789.'/
      DATA KCIF/'-+*/) '/
      DATA OPER/'-+*/()'/
C
      BROJ='                    '
      ISG=21
      ITH=21
      ITI=21
      ICR=0
      I=0
      KEND=0
   20 I=I+1
      IF (I.GT.80) GO TO 900
      IF (USTR(I:I).EQ.' ') GO TO 20
      IF (USTR(I:I).EQ.'E') THEN
        I=I+2
        KEND=KEND+1
        IZRAZ(KEND:KEND)='#'
        GO TO 20
      END IF
      IF (USTR(I:I).EQ.'S') THEN
        ICR=ICR+1
        ISG=ICR
        I=I+4
        GO TO 70
      ELSE IF (USTR(I:I).EQ.'T') THEN
        ICR=ICR+1
        IF (USTR(I+1:I+1).EQ.'H') THEN
          ITH=ICR
          I=I+4
        ELSE
          ITI=ICR
          I=I+3
        END IF
        GO TO 70
      END IF
      DO 80 J=1,6
        IF (USTR(I:I).EQ.OPER(J:J)) THEN
          KEND=KEND+1
          IF (J.EQ.3.AND.USTR(I+1:I+1).EQ.'*') THEN
            I=I+1
            IZRAZ(KEND:KEND)='^'
            GO TO 20
          END IF
          IZRAZ(KEND:KEND)=USTR(I:I)
          GO TO 20
        END IF
   80 CONTINUE
      DO 30 J=1,11
        IF (USTR(I:I).EQ.CIFRE(J:J)) THEN
          ISTART=I
          GO TO 40
        END IF
   30 CONTINUE
   40 I=I+1
      DO 50 J=1,6
        IF (USTR(I:I).EQ.KCIF(J:J)) THEN
          IF ((USTR(I:I).EQ.'+'.OR.USTR(I:I).EQ.'-').AND.
     &        (USTR(I-1:I-1).EQ.'E'.OR.USTR(I-1:I-1).EQ.'e')) GO TO 40
          I=I-1
          IEND=I
          GO TO 60
        END IF
   50 CONTINUE
      GO TO 40
   60 BROJ='                    '
      BROJ(20-IEND+ISTART:20)=USTR(ISTART:IEND)
      ICR=ICR+1
      READ(BROJ,'(F20.0)')ACR(ICR)
   70 KEND=KEND+1
      IZRAZ(KEND:KEND)='A'
      IF (ICR.GT.9) THEN
        WRITE(IZRAZ(KEND+1:KEND+2),'(I2)')ICR
        KEND=KEND+2
      ELSE
        WRITE(IZRAZ(KEND+1:KEND+1),'(I1)')ICR
        KEND=KEND+1
      END IF
      GO TO 20
  900 KEND=KEND+1
      IZRAZ(KEND:KEND)='$'
      RETURN
      END
C============================================================================
      SUBROUTINE KONVEB(IZRAZ,RP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER*250 IZRAZ,RP,PRP
      CHARACTER*20 STACK,POM
      CHARACTER*8 OPER
      CHARACTER*1 B
      DIMENSION IPR(20),IOPER(7)
      DATA OPER/'-+*/^#()'/
      DATA IOPER/1,1,2,2,3,4,0/
C
      STACK='$'
      DO 10 I=1,20
   10   IPR(I)=0
      KIZRAZ=KRAJ(IZRAZ)
      IO=1
      IX=1
      DO 20 I=1,KIZRAZ
        B=IZRAZ(I:I)
        IF (B.EQ.' ') GO TO 20
        K=0
        DO 30 J=1,8
          IF (B.EQ.OPER(J:J)) THEN
            K=J
            GO TO 40
          END IF
   30   CONTINUE
   40   IF (K.EQ.8) THEN
   50     KSTACK=KRAJ(STACK)
          IF (KSTACK.EQ.0) GO TO 20
          B=STACK(1:1)
          STACK=STACK(2:KSTACK+1)
          STACK(KSTACK+1:KSTACK+1)=' '
          IX=IX-1
          IF (B.EQ.'(') GO TO 20
          KRP=KRAJ(RP)
          IF (KRP.NE.0) THEN
            PRP=RP
            RP=PRP(1:KRP)//B//'$'
          ELSE
            RP=B//'$'
          END IF
          GO TO 50
        ELSE IF (K.EQ.0) THEN
          IO=0
          KRP=KRAJ(RP)
          IF (KRP.NE.0) THEN
            PRP=RP
            RP=PRP(1:KRP)//B//'$'
          ELSE
            RP=B//'$'
          END IF
        ELSE
          IS=IOPER(K)
          IF (IS.EQ.1.AND.IO.NE.0) THEN
            IS=6
            IF (B.EQ.'-') B='!'
            IF (B.EQ.'+') B='?'
          END IF
          IF (IS.NE.0) THEN
            IO=1
   70       IF (IX.EQ.1) GO TO 60
            IF (IPR(IX-1).LT.IS) GO TO 60
            KRP=KRAJ(RP)
            IF (KRP.NE.0) THEN
              PRP=RP
              RP=PRP(1:KRP)//STACK(1:1)//'$'
            ELSE
              RP=STACK(1:1)//'$'
            END IF
            IX=IX-1
            KSTACK=KRAJ(STACK)+1
            STACK=STACK(2:KSTACK)
            STACK(KSTACK:KSTACK)=' '
            GO TO 70
          END IF
   60     IPR(IX)=IS
          POM=STACK
          STACK=B//POM
          IX=IX+1
        END IF
   20 CONTINUE
      KSTACK=KRAJ(STACK)
      IF (KSTACK.NE.0) THEN
        KRP=KRAJ(RP)
        PRP=RP
        RP=PRP(1:KRP)//STACK(1:KSTACK)//'$'
      END IF
      RETURN
      END
C============================================================================
      FUNCTION KRAJ(A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 A
      DO 10 I=1,250
   10   IF (A(I:I).EQ.'$') GO TO 20
   20 KRAJ=I-1
      RETURN
      END
C=======================================================================
C     (rakic)
C=======================================================================
      SUBROUTINE UMOD41(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 41 Drucker-Prager
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 41 
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.   COEFFICIENTS FOR DEFINITION CURVE F1=0.0
CE.   FUNMAT(3,MAT) - COEFFICIENT   k 
CE.   FUNMAT(6,MAT) - COEFFICIENT   alfa
CE.   FUNMAT(7,MAT) - TENSION CUTOFF   T
CE.   FUNMAT(8,MAT) - INITIAL CAP POSITION   Lo=Xo
C
CE.   DATA FOR HARDENING FUNCTION epv 
CE.   FUNMAT(9,MAT)  - COEFFICIENT   W 
CE.   FUNMAT(10,MAT) - COEFFICIENT   D 
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD41'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,8)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,8)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=9,11)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,11)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) FUNMAT(3,MAT),(FUNMAT(J,MAT),J=6,8)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) FUNMAT(3,MAT),(FUNMAT(J,MAT),J=6,8)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2065) (FUNMAT(J,MAT),J=9,10)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=9,10)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     41 (TLO PLASTICNOST - DRUCKER-PRAGER)
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F1=0 '/
     121X,'k',11X,'Alfa',12X,'T ',12X,'Lo'/
     116X,4(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    41 (SOIL PLASTICITY - DRUCKER-PRAGER)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F1=0 '/
     121X,'k',11X,'Alfa',12X,'T ',12X,'Lo'/
     116X,4(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C     (rakic)
C=======================================================================
      SUBROUTINE UMOD42(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 42 Mohr-Coulomb
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 42 
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.     FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.     FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.     COEFFICIENTS FOR DEFINITION CURVE F=0.0
CE.     FUNMAT(3,MAT) - Cohesion - ce 
CE.     FUNMAT(4,MAT) - Friction angle - phi
CE.     FUNMAT(5,MAT) - Dilatation angle - psi
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD42'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,5)
cr      CALL ISPITA(ACOZ)
cr      IF(INDFOR.EQ.1)
cr     1READ(IULAZ,*) (FUNMAT(I,MAT),I=9,11)
cr      IF(INDFOR.EQ.2)
cr     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,11)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) FUNMAT(3,MAT),(FUNMAT(J,MAT),J=6,8)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) FUNMAT(3,MAT),(FUNMAT(J,MAT),J=6,8)
c      IF(ISRPS.EQ.0)
c     1WRITE(IZLAZ,2065) (FUNMAT(J,MAT),J=9,10)
c      IF(ISRPS.EQ.1)
c     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=9,10)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     42 (TLO PLASTICNOST - MOHR-COULOMB)
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F=0 '/
     121X,'ce',11X,'phi',12X,'psi'/
     116X,3(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    42 (SOIL PLASTICITY - MOHR-COULOMB)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F1=0 '/
     121X,'ce',11X,'phi',12X,'psi'/
     116X,3(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C     (rakic)
C=======================================================================
      SUBROUTINE UMOD43(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 43 Hoek-Brown
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 43 
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.   COEFFICIENTS FOR DEFINITION CURVE F=0.0
CE.   FUNMAT(3,MAT) - COEFFICIENT   em
CE.   FUNMAT(4,MAT) - COEFFICIENT   emd
CE.   FUNMAT(5,MAT) - COEFFICIENT   smc
CE.   FUNMAT(6,MAT) - COEFFICIENT   es
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD43'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,6)
CR    Sledecu linija pod komentarom da ne bi citao prazan red iz DAT
CR      CALL ISPITA(ACOZ)
c      IF(INDFOR.EQ.1)
c     1READ(IULAZ,*) (FUNMAT(I,MAT),I=9,11)
c      IF(INDFOR.EQ.2)
c     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,11)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=3,6)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=3,6)
c      IF(ISRPS.EQ.0)
c     1WRITE(IZLAZ,2065) (FUNMAT(J,MAT),J=9,10)
c      IF(ISRPS.EQ.1)
c     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=9,10)
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     43 (TLO PLASTICNOST - Hoek-Brown)
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F=0 '/
     121X,'smc',11X,'em',12X,'emd',12X,'es'/
     116X,4(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    43 (SOIL PLASTICITY - Hoek-Brown)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F=0 '/
     121X,'smc',11X,'em',12X,'emd',12X,'es'/
     116X,4(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C     (rakic)
C=======================================================================
      SUBROUTINE UMOD44(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 44 Generalized Hoek-Brown
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 44 
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.   COEFFICIENTS FOR DEFINITION CURVE F=0.0
CE.   FUNMAT(3,MAT) - COEFFICIENT   smc
CE.   FUNMAT(4,MAT) - COEFFICIENT   em
CE.   FUNMAT(6,MAT) - COEFFICIENT   emd
CE.   FUNMAT(6,MAT) - COEFFICIENT   es
CE.   FUNMAT(7,MAT) - COEFFICIENT   ahb
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD44'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,7)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,7)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=3,7)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=3,7)
c
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     44 (TLO PLASTICNOST - Gen.Hoek-Brown)
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F=0 '/
     121X,'smc',11X,'em',12X,'emd',12X,'es',12X,'ahb'/
     116X,5(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    44 (SOIL PLASTICITY - Gen.Hoek-Brown)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F=0 '/
     121X,'smc',11X,'em',12X,'emd',12X,'es',12X,'ahb'/
     116X,5(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
C     (rakic)
C=======================================================================
      SUBROUTINE UMOD45(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 45 Maksimovic
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 45 
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
CE.   FUNMAT(1,MAT) - YOUNG*S MODULUS IN DIRECTION X  -  E
CE.   FUNMAT(2,MAT) - POISSON*S RATIO                 -  V
C
CE.   COEFFICIENTS FOR DEFINITION CURVE F=0.0
CE.   FUNMAT(3,MAT) - COEFFICIENT   fib
CE.   FUNMAT(4,MAT) - COEFFICIENT   dfi
CE.   FUNMAT(5,MAT) - COEFFICIENT   pn
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(11,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD45'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,5)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,5)
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
C            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
C            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2045) (FUNMAT(J,MAT),J=3,5)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=3,5)
c
      RETURN
C
 1000 FORMAT(7F10.0)
 1010 FORMAT(6F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =     44 (TLO PLASTICNOST - Maksimovic)
     2'///11X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//
     16X,'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'
     1//16X,'M O D U L    E L A S T I C N O S T I  .  E =',1PD12.5//
     216X,'P O I S S O N O V    B R O J  .........  V =',1PD12.5//)
 2045 FORMAT(//
     116X,'K O E F I C I J E N T I  K R I V E F=0 '/
     121X,'fib',11X,'dfi',12X,'pn'/
     116X,3(1PD12.5,2X)//)
 2065 FORMAT(//
     1//16X,'K O E F I C I J E N T I  K R I V E  O J A C A N J A'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    44 (SOIL PLASTICITY - Maksimovic)'
     2///11X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S       M O D U L U S  ......  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O  .......  V =',1PD12.5//)
 6045 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  C U R V E  F=0 '/
     121X,'fib',11X,'dfi',12X,'pn'/
     116X,3(1PD12.5,2X)//)
 6065 FORMAT(//
     112X,'C O E F F I C I J E N T S  F O R  H A R D E N I N G  ',  
     1 'C U R V E'/
     120X,'W ',13X,'D '/16X,2(1PD12.5,2X))
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE UMOD52(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 52
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 52
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(160,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD52'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,24)
      IF(INDFOR.EQ.2) THEN
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,8)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,16)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=17,24)
      END IF  
C      
C---- Provera vrste kristala
      CHECK=0.0D0
      DO J=10,21
         CHECK=CHECK+DABS(FUNMAT(J,MAT))
      END DO

      IF (CHECK.EQ.0.0D0) THEN
         DO J=4,9
            CHECK=CHECK+DABS(FUNMAT(J,MAT))
         END DO

         IF (CHECK.EQ.0.0D0) THEN
            IF ((FUNMAT(3,MAT)).EQ.0.0D0) THEN
C-----  Isotropic material
                ivrsmat=1
            ELSE

C-----  Cubic material
                ivrsmat=2
            END IF

         ELSE

C-----  Orthotropic metarial                
                ivrsmat=3
         END IF

      ELSE

C-----  General anisotropic material
                ivrsmat=4
      END IF
         
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=25,32)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=25,32)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=33,40)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=33,40)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=41,48)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=41,48)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=49,56)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=49,56)          
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=57,64)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=57,64)          
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=65,72)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (FUNMAT(I,MAT),I=65,72)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=73,80)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=73,80)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=81,88)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=81,88)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=89,96)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=89,96)          
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=97,104)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=97,104)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=105,112)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=105,112)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=113,120)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=113,120)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=121,128)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=121,128)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=129,136)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=129,136)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=137,144)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=137,144)          
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=145,152)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=145,152)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=153,160)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=153,160)
    
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
C         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF

      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0.and.ivrsmat.eq.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.0.and.ivrsmat.eq.2)
     1WRITE(IZLAZ,6020) (FUNMAT(J,MAT),J=1,3)
      IF(ISRPS.EQ.0.and.ivrsmat.eq.3)
     1WRITE(IZLAZ,6030) (FUNMAT(J,MAT),J=1,9)
      IF(ISRPS.EQ.0.and.ivrsmat.eq.4)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,21)
      IF(ISRPS.EQ.1.and.ivrsmat.eq.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,2)
      IF(ISRPS.EQ.1.and.ivrsmat.eq.2)
     1WRITE(IZLAZ,6020) (FUNMAT(J,MAT),J=1,3)
      IF(ISRPS.EQ.1.and.ivrsmat.eq.3)
     1WRITE(IZLAZ,6030) (FUNMAT(J,MAT),J=1,9)
      IF(ISRPS.EQ.1.and.ivrsmat.eq.4)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,21)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=25,25)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6045) (FUNMAT(J,MAT),J=25,25)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6051) (FUNMAT(J,MAT),J=33,38)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6051) (FUNMAT(J,MAT),J=33,38)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6052) (FUNMAT(J,MAT),J=41,46)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6052) (FUNMAT(J,MAT),J=41,46)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6053) (FUNMAT(J,MAT),J=49,54)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6053) (FUNMAT(J,MAT),J=49,54)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6060) (FUNMAT(J,MAT),J=57,62)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060) (FUNMAT(J,MAT),J=57,62)               
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6061) (FUNMAT(J,MAT),J=65,70)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6061) (FUNMAT(J,MAT),J=65,70)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=73,74)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6065) (FUNMAT(J,MAT),J=73,74)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6066) (FUNMAT(J,MAT),J=81,82)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6066) (FUNMAT(J,MAT),J=81,82)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6067) (FUNMAT(J,MAT),J=89,90)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6067) (FUNMAT(J,MAT),J=89,90)               
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6070) (FUNMAT(J,MAT),J=97,99)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) (FUNMAT(J,MAT),J=97,99)     
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6071) (FUNMAT(J,MAT),J=105,106)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6071) (FUNMAT(J,MAT),J=105,106)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6072) (FUNMAT(J,MAT),J=113,115)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6072) (FUNMAT(J,MAT),J=113,115)    
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6073) (FUNMAT(J,MAT),J=121,122)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6073) (FUNMAT(J,MAT),J=121,122)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6074) (FUNMAT(J,MAT),J=129,131)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6074) (FUNMAT(J,MAT),J=129,131)     
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6075) (FUNMAT(J,MAT),J=137,138)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6075) (FUNMAT(J,MAT),J=137,138)     
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6085) (FUNMAT(J,MAT),J=145,146)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6085) (FUNMAT(J,MAT),J=145,146)     
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6095) (FUNMAT(J,MAT),J=153,155)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6095) (FUNMAT(J,MAT),J=153,155)
  
      RETURN
C
 1000 FORMAT(8F10.0)
 1001 FORMAT(F10.0)
 1010 FORMAT(8F5.0) 
C-----------------------------------------------------------------------   
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    52 (CRISTAL MODEL)'///
     111X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'YOUNG S  MODULUS......  E =',1PD12.5//
     116X,  'POISSON S  RATIO......  V =',1PD12.5//)
 6020 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'ELASTIC  CONSTANT  C11 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  C12 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  C44 = ',1PD12.5//)
 6030 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'ELASTIC  CONSTANT  D1111 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1122 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2222 = ',1PD12.5// 
     116X,'ELASTIC  CONSTANT  D1133 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2233 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D3333 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1212 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1313 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2233 = ',1PD12.5//)
 6040 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'ELASTIC  CONSTANT  D1111 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1122 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2222 = ',1PD12.5//  
     116X,'ELASTIC  CONSTANT  D1133 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2233 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D3333 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1112 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2212 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D3312 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1212 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1113 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2213 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D3313 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1213 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1313 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1123 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2223 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D3323 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1223 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D1323 = ',1PD12.5//
     116X,'ELASTIC  CONSTANT  D2323 = ',1PD12.5//)      
 6045 FORMAT(//
     112X,'NUMBER OF SETS OF SLIP SYSTEMS = '
     121X,1PD12.5//)
 6051 FORMAT(//
     112X,'NORMAL TO SLIP PLANE IN THE 1st SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//     
     112X,'SLIP DIRECTION IN THE 1st SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//)
 6052 FORMAT(//
     112X,'NORMAL TO SLIP PLANE IN THE 2nd SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//     
     112X,'SLIP DIRECTION IN THE 2nd SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//)
 6053 FORMAT(//
     112X,'NORMAL TO SLIP PLANE IN THE 3rd SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//
     112X,'SLIP DIRECTION IN THE 2nd SET : ',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//)        
 6060 FORMAT(//     
     112X,'DIRECTION IN LOCAL SYSTEM OF THE 1st VECTOR :',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//
     112X,'DIRECTION IN GLOBAL SYSTEM OF THE 1st VECTOR:',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//)
 6061 FORMAT(//     
     112X,'DIRECTION IN LOCAL SYSTEM OF THE 2nd VECTOR :',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//
     112X,'DIRECTION IN GLOBAL SYSTEM OF THE 2nd VECTOR:',2X,1PD12.5,2X,
     11PD12.5,2X,1PD12.5//)
 6065 FORMAT(//
     112X,'POWER HARDENING EXP. AND HARD. COEFF. 1st SET(n,adot)'//
     112X,'n=',2X,1PD12.5,2X,'adot=',2X,1PD12.5,2X)
 6066 FORMAT(//
     112X,'POWER HARDENING EXP. AND HARD. COEFF. 2nd SET(n,adot)'//
     112X,'n=',2X,1PD12.5,2X,'adot=',2X,1PD12.5,2X)
 6067 FORMAT(//
     112X,'POWER HARDENING EXP. AND HARD. COEFF. 3rd SET(n,adot)'//
     112X,'n=',2X,1PD12.5,2X,'adot=',2X,1PD12.5,2X)  
 6070 FORMAT(//
     112X,'Initial hardening modulus, saturation stress and initial crit
     1ical resolved shear stress 1st set'//
     112X,'h0=',2X,1PD12.5,2X,'taus=',2X,1PD12.5,2X,'tau0=',2X,1PD12.5) 
 6071 FORMAT(//
     112X,'Ratios of latent to self-hardening in the same and different
     1sets of slip systems 1st set'//
     112X,'q=',2X,1PD12.5,2X,'q1=',2X,1PD12.5,2X)
 6072 FORMAT(//
     112X,'Initial hardening modulus, saturation stress and initial crit
     1ical resolved shear stress 2nd set'//
     112X,'h0=',2X,1PD12.5,2X,'taus=',2X,1PD12.5,2X,'tau0=',2X,1PD12.5) 
 6073 FORMAT(//
     112X,'Ratios of latent to self-hardening in the same and different
     1sets of slip systems 2nd set'//
     112X,'q=',2X,1PD12.5,2X,'q1=',2X,1PD12.5,2X)
 6074 FORMAT(//
     112X,'Initial hardening modulus, saturation stress and initial crit
     1ical resolved shear stress 3rd set'//
     112X,'h0=',2X,1PD12.5,2X,'taus=',2X,1PD12.5,2X,'tau0=',2X,1PD12.5) 
 6075 FORMAT(//
     112X,'Ratios of latent to self-hardening in the same and different
     1sets of slip systems 3rd set'//
     112X,'q=',2X,1PD12.5,2X,'q1=',2X,1PD12.5,2X)          
 6085 FORMAT(//
     112X,'theta=',2X,1PD12.5,2X,'nlgeom=',2X,1PD12.5,2X)
 6095 FORMAT(//
     112X,'ITRATN=',2X,1PD12.5,2X,'ITRMAX=',2X,1PD12.5,2X,'GAMERR=',2X,
     11PD12.5)         
C     
C-----------------------------------------------------------------------
      END
C=======================================================================  
C=======================================================================
      SUBROUTINE UMOD53(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 53
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 53
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(24,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD53'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,24)
      IF(INDFOR.EQ.2) THEN
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,8)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,16)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=17,24)
      END IF
    
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF

      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,24)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,24)
  
      RETURN
C
 1000 FORMAT(8F10.0)
 1001 FORMAT(F10.0)
 1010 FORMAT(8F5.0) 
C-----------------------------------------------------------------------   
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    53 (SMA MODEL)'///
     111X,'MATERIAL NUMBER =',I5)
 6040 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'MATERIAL  CONSTANT  1 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  2 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  3 = ',1PD12.5//  
     116X,'MATERIAL  CONSTANT  4 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  5 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  6 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  7 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  8 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  9 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  10 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  11 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  12 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  13 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  14 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  15 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  16 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  17 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  18 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  19 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  20 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  21 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  22 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  23 = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  24 = ',1PD12.5//)
C     
C-----------------------------------------------------------------------
      END
C=======================================================================          
C=======================================================================
      SUBROUTINE UMOD54(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 54
CS.    P R O G R A M
CS.        ZA UCITAVANJE PODATAKA O MATERIJALNOM MODELU BROJ 54
C .
CE.             MAT - MATERIAL NUMBER
CS.             MAT - MATERIJAL BROJ
C .
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(17,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD54'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,17)
      IF(INDFOR.EQ.2) THEN
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,4)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=5,8)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=9,10)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=11,14)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=15,16)
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) (FUNMAT(I,MAT),I=17,17)                 
      END IF
    
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
            CALL TGRMAT9(AMAT,MATG,MAT,IEL,49)
            CALL TGRMAT9(AMAT,MAT,MAT,IEL,49)
            CALL TGRMAT(AMAT,MAT,49)
      ENDIF

      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,17)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) (FUNMAT(J,MAT),J=1,17)
  
      RETURN
C
 1000 FORMAT(8F10.0)
 1001 FORMAT(F10.0)
 1010 FORMAT(8F5.0) 
C-----------------------------------------------------------------------   
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    54 (SMA VLADA MODEL)'///
     111X,'MATERIAL NUMBER =',I5)
 6040 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'
     1//16X,'MATERIAL  CONSTANT  dE = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dNi = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dEpsilonL = ',1PD12.5//  
     116X,'MATERIAL  CONSTANT  dT = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dTsas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dTfas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dTssa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dTfsa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dCas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dCsa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dSsas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dSfas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dSssa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dSfsa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dBas = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dBsa = ',1PD12.5//
     116X,'MATERIAL  CONSTANT  dAlpha = ',1PD12.5//)
C     
C-----------------------------------------------------------------------
      END
C=======================================================================     
