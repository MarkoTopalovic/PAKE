C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C                                                     .             .
C                AUTOMATSKO OPTERECIVANJE             .  VER.  1.01 .
C                                                     .  05-02-1990 .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE AUTMNW (D,U,DU1,INDC,LL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C       AUTOMATSKO OPTERECIVANJE ; ALGORITAM KONSTANTAN LUK
C..................................................................
C       U    INKREMENT VEKTORA POMERANJA = U(T+DT) - U(T)    REPER N9
C       D    POMERNAJA KOJA ODGOVARAJU REZIDUALNIM SILAMA
C            ( DU SA CRTOM )                                 REPER LRTDT
C       DU1  POMERNAJA KOJA ODGOVARAJU SPOLJASNJIM SILAMA
C            ( DU SA DVE CRTE )                              REPER N10
C       RESID  REZIDUALNE SILE       REPER LFTDT
C
C       KUSLOV  TIP USLOVA ZA IZBOR ZNAKA PROMENE LAMDA (DLAM)
C               = 1    DETERMINANTA SISTEMA
C               = 2    POSTOJANJE NEGATIVNOG PIVOTA
C
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /AUTST1/ TLAM,DLAM,TLAM0,TLAM00,TW,W00,DTUNR2,UNORM2,
     1                IWORK,ITE0
      COMMON /AUTST2/ PARAM(5)
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ENERGI/ ENE1,ENE2,ENE
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /PROBAS/ IILS
      COMMON /GRUPER/ LIGRUP
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK0/ IRESCT
      COMMON /RESTAR/ TSTART,IREST
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      EQUIVALENCE (PARAM(1),ARCG20), (PARAM(2),ARCMAX), (PARAM(3),DET0),
     1            (PARAM(4),SG0)   , (PARAM(5),ZLAM)
C
      DIMENSION D(*),U(*),DU1(*)
C      DFLOAT(M)=DBLE(FLOAT(M))
C
      IF(IDEBUG.GT.0) PRINT *, ' AUTMNW'
      KUSLOV=2
      IUPMAX=10
      IUPDT=0
      ITRR=-1
    1 ITER=-1
      INDC=0
C....  INICIJALIZOVANJE JEDINICNIH VEKTORA  V = V(T)
      CALL REWDRA(1)     
C
      CALL RSTAZ(A(LIPODS),N10,38)
      IILS=-1
      CALL BACKSB ( DU1 )
      IILS=0
C.................................................
C                PRVI KORAK                      .
C.................................................
      IF(KOR.GT.1) GO TO 100
      IRESCT=0
      IF(IUPDT.GT.0) IRESCT=1
      IF(IUPDT.GT.0) WRITE(IZLAZ,1021) IUPDT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,1031) DTUK
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6031) DTUK
C
      IWORK=0
      CALL CLEAR(U,JEDN)
   30 ITER=ITER+1
      ITRR=ITRR+1
C
C.....  NULTA ITERACIJA, PRVI KORAK
C
      IF(ITER.EQ.0) THEN
        IF(DABS(DU1(KEQ)).LT.1.D-10)THEN
          WRITE(IZLAZ,1040)
          STOP'  AUTMN3'
        ENDIF
                    TLAM=DTUK/DU1(KEQ)
                    SG0=TLAM/DABS(TLAM)
                    DO 10 I=1,JEDN
   10               D(I)=TLAM*DU1(I)
                    CALL RSTAZ(A(LIPODS),LFTDT,38)
C
             ELSE
C
C.....  ITERACIJA 1,2,3..., PRVI KORAK
C
           IF(LL.EQ.3)THEN
C 	write(3,*) ' iter,tlam',iter,tlam
             CALL TANGK8 (TLAM)
             CALL RSTAZ(A(LIPODS),N10,38)
C      call wrr6(a(n10),jedn,'n10u')
C      call wrr6(du1,jedn,'du1u')
             IILS=-1
             CALL BACKSB ( DU1 )
C      call wrr6(du1,jedn,'du1i')
             IILS= 0
           ELSE
             CALL FRESAU(TLAM)
           ENDIF
             CALL BACKSB (D)
             DLAM=-D(KEQ)/DU1(KEQ)
             TLAM=TLAM+DLAM
             DO 20 I=1,JEDN
   20        D(I)=D(I)+DLAM*DU1(I)
                    ENDIF
CE    AZURIRANJE KONTAKTNIH SILA
CS    ARANGE CONTACT FORCES
      IF(ICONT.NE.0) CALL CONCOR(A(LIGRUP),9)
C
C	write(3,*) ' iter,tlam,dlam',iter,tlam,dlam
C      call wrr6(a(lftdt),jedn,'ftdt')
C      call wrr6(a(n10),jedn,' n10')
C      call wrr6(u,jedn,'   u')
C      call wrr6(d,jedn,'   d')
C      call wrr6(du1,jedn,' du1')
      CALL KONVER(INDC)
      DO 62 I=1,JEDN
   62 U(I)=U(I)+D(I)
C
      IF(INDC.GT.0.AND.ITER.GT.0) THEN
C
C.....  MEMORISANJE PARAMETARA  -  KRAJ PRVOG KORAKA
C
                                 DTUNR2=ALFG*ALFG*RNORM2(U,JEDN)
                                 CALL POMREC
C
                                 ZLAM=DSQRT(DTUNR2)/TLAM
                                 TLAM=TLAM*ZLAM
C.....                           SKALIRANJE SPOLJASNJIH SILA
                                 CALL RSTAZ(A(LIPODS),N10,38)
                                 DO 70 I=1,JEDN
   70                            DU1(I)=DU1(I)/ZLAM
                                 CALL WSTAZ(A(LIPODS),N10,38)
C
                                 ARCG20=TLAM*TLAM+DTUNR2
                                 ARCMAX=ARCG20*DELS*DELS
C
C                                TW=ENE1/2.
C                                TDTW=TW
                                 TLAM0=TLAM
                                 ITE0=ITER
                                 ITER=ITRR
                                 DET0=DETER
          VREME=TLAM/ZLAM
CC          TTTTT=TLAM/ZLAM
CC          WRITE(IZLAZ,1010) TTTTT
                                 RETURN
                                 ENDIF
      IF(ITER.EQ.MAXIT) THEN
      IF(IUPDT.EQ.IUPMAX)THEN
                         WRITE(IZLAZ,1001)
                         STOP
                         ENDIF
                       DTUK=DTUK/2.D0
                       IUPDT=IUPDT+1
                       GO TO 1
                       ENDIF
C
      GO TO 30
C.....................................................................
C                KORACI      2 , 3 ....                              .
C.....................................................................
  100 IF(IUPDT.GT.0) WRITE(IZLAZ,1020) IUPDT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,1030) DSQRT(ARCG20)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) DSQRT(ARCG20)
C
      TLAM=TLAM0
      SG=SG0
      ARCLG2=ARCG20
      IF(KUSLOV.EQ.1) THEN
        IF(DETER*DET0.LT.0.D0) SG=-SG0
      ENDIF
      IF(KUSLOV.EQ.2) THEN
        SG=1.0D0
        IF(IPIVOT.EQ.-1) SG=-SG
      ENDIF
      CALL CLEAR(U,JEDN)
C
      IF(ICONT.NE.0.AND.IUPDT.GT.0)THEN
        DUM=ARCLG2
        IREST=1
        IRESCT=1
        CALL BACKUP(2)
        IRESCT=0
        IREST=0
        ARCLG2=DUM
        ARCG20=DUM
      ENDIF
C--------------------------------------------
  120 ITER=ITER+1
      ITRR=ITRR+1
      IF(ITER.EQ.0) THEN
C
C.....  NULTA ITERACIJA ---------------------
C
            CALL RSTAZ(A(LIPODS),LFTDT,38)
            IF(IWORK.EQ.0) THEN
                           DU12=DOT(DU1,DU1,NEQ)
                           DLAM=DSQRT(ARCLG2/(DU12+1.D0))
                             ELSE
C                            A1=DOT(RESID,DU1,NEQ)
C                            DLAM=-TLAM0+DSQRT(TLAM0*TLAM0+DABS(TW/A1))
C                            A1=A1*2.
                           ENDIF
                           DLAM=SG*DLAM
                    CALL JEDNAK(D,DU1,DLAM,JEDN)
      ELSE
C
C.....  ITERACIJE  1 ,2 , 3 ... ---------------
C
           IF(LL.EQ.3)THEN
             CALL TANGK8 (TLAM)
             CALL RSTAZ(A(LIPODS),N10,38)
             CALL BACKSB ( DU1 )
             DU12=DOT(DU1,DU1,NEQ)
           ELSE
             CALL FRESAU(TLAM)
           ENDIF
      CALL BACKSB (D)
C
       C11=DOT(U,DU1,NEQ)
       C12=DOT(U,U  ,NEQ)
       C13=DOT(U,D  ,NEQ)
       GAMA0=C12+C13
C---------------------------------------------------------
C....  KOEFICIJENTI KVADRATNE JEDNACINE
      IF(IWORK.EQ.0) THEN
         CALL KOEFAL (A1,A2,A3,U,D,DU1,DU12,C11,TLAM,TLAM0,ARCLG2,NEQ)
             ELSE
C            CALL RSTAZ(A(LIPODS),LFTDT,38)
C            CALL KOEFW (A1,A2,A3,RESID,D,TLAM,TDTW,JEDN,ITER)
                     ENDIF
C....  DISKRIMINANTA KVADRATNE JEDNACINE
      DSKR=A2*A2-2.*A1*A3
C....  IMAGINARNI KORENOVI
      IF(DSKR.LT.0.D0)THEN
      WRITE(IZLAZ,*)'     ***  IMAGINARNI KORENOVI ***'
C     IF(IWORK.EQ.1) STOP '***  DISKRIMINANTA .LT. 0.  ; IWORK .EQ. 1'
      IF(IUPDT.GT.IUPMAX)THEN
                         WRITE(IZLAZ,1000)
                         STOP
                         ENDIF
                      ARCG20=ARCG20/4.D0
                      IUPDT=IUPDT+1
C                     IWORK=1
                      GO TO 1
                      ENDIF
C....     RESENJE ZA DSKR=0.
      DLAM=-A2/A1
C
      IF(DSKR.GT.0.D0) THEN
                       DSKR=DSQRT(DSKR)
                       DLAM2=DLAM-DSKR/A1
                       DLAM=DLAM+DSKR/A1
                       ENDIF
C....    PROVERA PREKO GAMA PREMA  CRISFIELDU
C....    DLAML  -  LINEARNO RESENJE
      GAMA1=GAMA0+DLAM *C11
      GAMA2=GAMA0+DLAM2*C11
      IF(GAMA2.GT.0.D0.AND.GAMA1.LT.0.D0) DLAM=DLAM2
      IF(GAMA1.GT.0.D0.AND.GAMA2.GT.0.D0) THEN
             DLAML=-A3/A2
             IF(DABS(DLAM2-DLAML).LT.DABS(DLAM-DLAML)) DLAM=DLAM2
                                          ENDIF
C-------------------------------------------------
      DO 90 I=1,JEDN
   90 D(I)=D(I)+DLAM*DU1(I)
CE    AZURIRANJE KONTAKTNIH SILA
CS    ARANGE CONTACT FORCES
      IF(ICONT.NE.0) CALL CONCOR(A(LIGRUP),9)
C
      ENDIF
C-------------------------------------------------
      CALL KONVER(INDC)
      TLAM=TLAM+DLAM
      DO 92 I=1,JEDN
   92 U(I)=U(I)+D(I)
C
      IF(INDC.GT.0.AND.ITER.GT.0) THEN
                                  UNORM2=RNORM2(D,JEDN)
C-------------------------------- ALFA KRITERIJUM
      IF(IWORK.EQ.0.AND.
     1   UNORM2.GT.DTUNR2) THEN
      IF(IUPDT.GT.IUPMAX)THEN
                         WRITE(IZLAZ,1000)
                         STOP
                         ENDIF
                                     ARCG20=ARCLG2*DTUNR2/UNORM2
                                     IUPDT=IUPDT+1
                                     GO TO 1
                                     ENDIF
C
C.....  MEMORISANJE PARAMETARA  -  KRAJ 2 ,3 ,4 ... KORAKA
C
                    CALL POMREC
                    TLAM00=TLAM0
                    TLAM0=TLAM
C                   W00=TW
C                   TW=TDTW
                    ITE0=ITER
                    DET0=DETER
                    SG0=SG
C....        TEKUCA VREDNOST LUKA  (KVADRAT)
                    FLT=FLOAT(ITEOPT)
                    FL0=FLOAT(ITE0)
                    ARCLG2=ARCLG2*FLT/FL0
                    ARCG20=DMIN1(ARCLG2,ARCMAX)
                    ITER=ITRR
C....        TEKUCA VREDNOST RADA
C                   TDTW=TW*TLAM0/TLAM00
C....        PROVERA ZA PRELAZ NA EXT. WORK
C                   WW=DABS(DABS(TW/W00)-1.)
C                   IWORK=0
C                   IF(WW.LT.DELS) IWORK=1
          VREME=TLAM/ZLAM
CC          TTTTT=TLAM/ZLAM
CC          WRITE(IZLAZ,1010) TTTTT
                    RETURN
                    ENDIF
      IF(ITER.EQ.MAXIT) THEN
      IF(IUPDT.EQ.IUPMAX)THEN
                         WRITE(IZLAZ,1000)
                         STOP
                         ENDIF
                      ARCG20=ARCG20/4.D0
                      IUPDT=IUPDT+1
C                       IF(IWORK.EQ.1) RETURN
C                       IWORK=1
                        GO TO 1
                        ENDIF
C
      GO TO 120
C-----------------------------------------------------------------
 1000 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA LUCNE DUZINE U TEK
     1UCEM KORAKU')
 1001 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA ZADATOG POMERANJA
     1U TEKUCEM KORAKU ')
CC 1010 FORMAT(//' SKALARNI MNOZIOC VEKTORA OPTERECENJA (LAMDA) = ',G10.4
CC     1///)
 1020 FORMAT(///' PROCEDURA DIVERGIRA - KOREKCIJA ',I3,'  LUCNE DUZINE')
 1021 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ ITERACIJA - KOREKCIJA ',I3,
     1'  ZADATOG POMERANJA')
 1030 FORMAT(/////' LUCNA DUZINA (ARC LENGTH) = ',G10.4)
 1031 FORMAT(/////' TEKUCE ZADATO POMERANJE (DTUK) = ',G10.4)
 1040 FORMAT(//' METODE AUTOMATSKOG OPTERECIVANJA'/' DEFINISATI ZADATO',
     2' POMERANJE (DTUK) ZA DRUGI STEPEN SLOBODE')
C-----------------------------------------------------------------
 6000 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA LUCNE DUZINE U TEK
     1UCEM KORAKU')
 6001 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA ZADATOG POMERANJA
     1U TEKUCEM KORAKU ')
CC 6010 FORMAT(//' SKALARNI MNOZIOC VEKTORA OPTERECENJA (LAMDA) = ',G10.4
CC     1///)
 6020 FORMAT(///' PROCEDURA DIVERGIRA - KOREKCIJA ',I3,'  LUCNE DUZINE')
 6021 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ ITERACIJA - KOREKCIJA ',I3,
     1'  ZADATOG POMERANJA')
 6030 FORMAT(/////' ARC LENGTH = ',G10.4)
 6031 FORMAT(/////' CURRENT PRESCRIBED DISPLACEMENT (DTUK) = ',G10.4)
 6040 FORMAT(//' METODE AUTOMATSKOG OPTERECIVANJA'/' DEFINISATI ZADATO',
     2' POMERANJE (DTUK) ZA DRUGI STEPEN SLOBODE')
      END
C======================================================================
C
C======================================================================
      SUBROUTINE AUTPRD ( D,U,DU1,INDC,LL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C       KONTROLISANA POMERANJA
C
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /AUTST1/ TLAM,DLAM,TLAM0,TLAM00,TW,W00,DTUNR2,UNORM2,
     1                IWORK,ITE0
      COMMON /AUTST2/ PARAM(5)
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /GRUPER/ LIGRUP
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK0/ IRESCT
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      EQUIVALENCE (PARAM(1),DTUK0),  (PARAM(2),DTUKMX)
C
      DIMENSION D(*),U(*),DU1(*)
C      DFLOAT(M)=DBLE(FLOAT(M))
C
      IUPMAX=10
      IUPDT=0
      IF(KOR.GT.1)THEN
        CALL FRESA2(TLAM)
             IILS=-1
        CALL BACKSB(D)
             IILS= 0
        DRUK=D(KEQ)
      ELSE
        DRUK=0.D0
      ENDIF
      ITRR=-1
    1 ITER=-1
      INDC=0
C
      IF(KOR.EQ.1.AND.IUPDT.EQ.0) THEN
        DTUKMX=DTUK*ALFG
        DTUK0=DTUK
        TLAM0=0.D0
      ENDIF
C
      DTUK=DTUK0
      IRESCT=0
      IF(IUPDT.GT.0) IRESCT=1
      IF(IUPDT.GT.0) WRITE(IZLAZ,1020) IUPDT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,1030) DTUK
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030) DTUK
C
C....  INICIJALIZOVANJE JEDINICNIH VEKTORA  V = V(T)
      CALL REWDRA(1)     
C
      CALL RSTAZ(A(LIPODS),N10,38)
             IILS=-1
      CALL BACKSB ( DU1 )
             IILS= 0
      CALL CLEAR(U,JEDN)
      TLAM=TLAM0
   30 ITER=ITER+1
      ITRR=ITRR+1
      IF(ITER.EQ.0) THEN
        IF(DABS(DU1(KEQ)).LT.1.D-10)THEN
          WRITE(IZLAZ,1040)
          STOP'  AUTPRD'
        ENDIF
      DLAM=(DTUK-DRUK)/DU1(KEQ)
      IF(DABS(DLAM).GT.DELS)THEN
        IF(IUPDT.EQ.IUPMAX)THEN
           WRITE(IZLAZ,1000)
           STOP
        ENDIF
        DTUK0=DTUK0/2.D0
        IUPDT=IUPDT+1
        GO TO 1
      ENDIF
      TLAM=TLAM+DLAM
      DO 10 I=1,JEDN
   10 D(I)=DLAM*DU1(I)
      CALL RSTAZ(A(LIPODS),LFTDT,38)
             ELSE
C....              ITER.GT.0
C
           IF(LL.EQ.3)THEN
             CALL TANGK8 (TLAM)
             CALL RSTAZ(A(LIPODS),N10,38)
             IILS=-1
             CALL BACKSB ( DU1 )
             IILS= 0
           ELSE
             CALL FRESAU(TLAM)
           ENDIF
C
             CALL BACKSB (D)
C
             DLAM=-D(KEQ)/DU1(KEQ)
             TLAM=TLAM+DLAM
             DO 20 I=1,JEDN
   20        D(I)=D(I)+DLAM*DU1(I)
                    ENDIF
C
CE    AZURIRANJE KONTAKTNIH SILA
CS    ARANGE CONTACT FORCES
      IF(ICONT.NE.0) CALL CONCOR(A(LIGRUP),9)
C
      CALL KONVER(INDC)
      DO 62 I=1,JEDN
   62 U(I)=U(I)+D(I)
C
      IF(INDC.GT.0.AND.ITER.GT.0) THEN
                                  CALL POMREC
C
                                  TLAM0=TLAM
C....  DEFINISANJE OPTIMALNOG ZADATOG POMERANJA
                    FLT=FLOAT(ITEOPT)
                    FL0=FLOAT(ITER)
                    DTUK0=DTUK*FLT/FL0
                    IF(DABS(DTUKMX).LT.DABS(DTUK0))DTUK0=DTUKMX
                    ITER=ITRR
C
                                  VREME=TLAM
CC                                  WRITE(IZLAZ,1010) TLAM
                                  RETURN
                                  ENDIF
      IF(ITER.EQ.MAXIT) THEN
      IF(IUPDT.EQ.IUPMAX)THEN
                         WRITE(IZLAZ,1000)
                         STOP
                         ENDIF
                        DTUK0=DTUK0/2.D0
                        IUPDT=IUPDT+1
                        GO TO 1
                        ENDIF
C
      GO TO 30
C-----------------------------------------------------------------
 1000 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA ZADATOG POMERANJA
     1U TEKUCEM KORAKU ')
CC 1010 FORMAT(//' SKALARNI MNOZIOC VEKTORA OPTERECENJA (LAMDA) = ',G10.4
CC     1///)
 1020 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ ITERACIJA - KOREKCIJA ',I3,
     1'  ZADATOG POMERANJA')
 1030 FORMAT(/////' TEKUCE ZADATO POMERANJE (DTUK) = ',G10.4)
 1040 FORMAT(//' METODE AUTOMATSKOG OPTERECIVANJA'/' DEFINISATI ZADATO',
     2' POMERANJE (DTUK) ZA DRUGI STEPEN SLOBODE')
C-----------------------------------------------------------------
 6000 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ KOREKCIJA ZADATOG POMERANJA
     1U TEKUCEM KORAKU ')
CC 6010 FORMAT(//' SKALARNI MNOZIOC VEKTORA OPTERECENJA (LAMDA) = ',G10.4
CC     1///)
 6020 FORMAT(///' POSTIGNUT MAKSIMALAN BROJ ITERACIJA - KOREKCIJA ',I3,
     1'  ZADATOG POMERANJA')
 6030 FORMAT(/////' CURRENT PRESCRIBED DISPLACEMENT (DTUK) = ',G10.4)
 6040 FORMAT(//' METODE AUTOMATSKOG OPTERECIVANJA'/' DEFINISATI ZADATO',
     2' POMERANJE (DTUK) ZA DRUGI STEPEN SLOBODE')
C-----------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE RBJED(ID,NELCV,NPI,ICVEL)
C....   REDNI BROJ JEDNACINE KOJA ODGOVARA KONTROLISANOM POMERANJU
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
C
      DIMENSION ID(NP,*),NELCV(*)
      KPNN=KPNOD
      IF(ICVEL.NE.0)KPNN=NELCV(KPNOD-NPI+1)
      KEQ=ID(KPNN,KPDIR)
      IF(KEQ.EQ.0)THEN
        WRITE(IZLAZ,2000)KPDIR,KPNOD
        STOP'  RBJED'
      ENDIF
 2000 FORMAT(///' ZA STEPEN SLOBODE',I5,'  CVORA',I5,'  POSTOJI OGRANICE
     1NJE'//' METODE AUTOMATSKOG OPTERECIVANJA'/' DEFINISATI ZADATO',
     2' POMERANJE (DTUK) ZA DRUGI STEPEN SLOBODE')
      RETURN
      END
C======================================================================
      SUBROUTINE POMREC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....  ZAPISIVANJE POMERANJA NA KRAJU KORAKA
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
C
      CALL UKUP2
      CALL WSTAZ(A(LIPODS),LRTDT,52)
C....  PREPISIVANJE JEDINICNIH VEKTORA
      CALL REWDRA(2)
      RETURN
      END
C======================================================================
      SUBROUTINE UKUP2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....  UKUPNA POMERANJA ( ITERATIVNI METODI  6,7,8,9 )
      include 'paka.inc'
      
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /UPRIRI/ LUPRI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
C
C....   PRIRASTAJ POMERANJA
                    CALL JEDNA1(A(LUPRI),A(LRTDT),JEDN)
C....   ROTIRANJE JEDINICNIH VEKTORA CVOROVA
                    CALL ROTDRV(1)
                    CALL REWDRV
C
      CALL RSTAZ(A(LIPODS),LRTDT,52)
C	call wrr6(a(lupri),jedn,'upr2')
C	call wrr6(a(lrtdt),jedn,'rt52')
      CALL ZBIRM1(A(LRTDT),A(N9),JEDN)
C	call wrr6(a(n9),jedn,'  n9')
C	call wrr6(a(lrtdt),jedn,'rtd2')
      RETURN
      END
C======================================================================
      SUBROUTINE FRESAU(TLAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....   REZIDUALNE SILE ( ITERATIVNI METODI  6  I  8 )
      include 'paka.inc'
      
      COMMON /GRUPER/ LIGRUP
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' FRESAU'
      CALL UKUP2
      CALL INTMNJ(A(LIGRUP),A(LIPODS))
      CALL RSTAZ(A(LIPODS),LRTDT,38)
      CALL ZBIRNN(A(LRTDT),A(LFTDT),TLAM,JEDN)
      IF(ICONT.NE.0)THEN
         CALL ZBIRM1(A(LRTDT),A(LRCTDT),JEDN)
         CALL CLEAR(A(LRCTDT),JEDN)
      ENDIF
      CALL JEDNA1(A(LFTDT),A(LRTDT),JEDN)
      RETURN
      END
C======================================================================
      SUBROUTINE TANGK8(TLAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....   MATRICA KRUTOSTI I REZIDUALNE SILE ( ITERATIVNI METOD  7 I 9)
      include 'paka.inc'
      
      COMMON /GRUPER/ LIGRUP
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TANGK8'
      CALL UKUP2
C	call wrr6(a(lsk),60,'skpr')
      CALL INTFNJ(A(LIGRUP),A(LIPODS))
C	call wrr6(a(lsk),60,'skpo')
      CALL RSTAZ(A(LIPODS),LRTDT,38)
C	call wrr6(a(lrtdt),jedn,'rtk8')
C	call wrr6(a(lftdt),jedn,'ftk8')
      CALL ZBIRNN(A(LRTDT),A(LFTDT),TLAM,JEDN)
C	call wrr6(a(lrtdt),jedn,'rtu8')
      IF(ICONT.NE.0)THEN
         CALL ZBIRM1(A(LRTDT),A(LRCTDT),JEDN)
         CALL CLEAR(A(LRCTDT),JEDN)
      ENDIF
      CALL JEDNA1(A(LFTDT),A(LRTDT),JEDN)
C	call wrr6(a(lftdt),jedn,'ftd8')
      RETURN
      END
C======================================================================
      SUBROUTINE FRESA2(TLAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....   REZIDUALNE SILE ( ITERATIVNI METODI  6,7,8  PRE ITERACIJA )
      include 'paka.inc'
      
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' FRESA2'
      CALL RSTAZ(A(LIPODS),LRTDT,38)
      CALL ZBIRNN(A(LRTDT),A(LFTDT),TLAM,JEDN)
      IF(ICONT.NE.0)THEN
         CALL ZBIRM1(A(LRTDT),A(LRCTDT),JEDN)
         CALL CLEAR(A(LRCTDT),JEDN)
      ENDIF
      CALL JEDNA1(A(LFTDT),A(LRTDT),JEDN)
      RETURN
      END
C======================================================================
      SUBROUTINE KOEFAL (A1,A2,A3,U,D,DU1,DU12,C11,TLAM,TLAM0,ARCLG2,
     1 JEDN)
C....   KOEF. KVADRATNE J-NE ZA "ARC-LENGTH"
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      DIMENSION D(*),U(*),DU1(*)
C
      C1=DOT(D,DU1,JEDN)
      C2=DOT(U,U  ,JEDN)
      C3=0.D0
      DO 50 I=1,JEDN
   50 C3=C3+(D(I)+2*U(I))*D(I)
C... U SLUCAJU KONTAKTA PRELAZI NA CILINDRICNI POSTUPAK
      DD=0.D0
      D1=0.D0
CC      IF(ICONT.EQ.0)THEN
        DD=TLAM-TLAM0
        D1=1.D0
CC      ENDIF
C
C... A1 IMA DVOSTRUKU VREDNOST , KASNIJE JE UZETO U OBZIR
      A1=2.D0*(DU12+D1)
      A2=(C1+C11+DD)*2.D0
      A3=C2-ARCLG2 +C3+DD*DD
      RETURN
C======================================================================
C     SUBROUTINE KOEFW (A1,A2,A3,R,D,TLAM,TDTW,JEDN,ITER)
C....   KOEF. KVADRATNE J-NE ZA "EXT-WORK"
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     DIMENSION D(*),R(*)
C
C     C2=DOT(R,D,JEDN)
C
C     A2=TLAM*A1+C2
C     A3=2.*TLAM*C2
C     IF(ITER.EQ.0) A3=A3-2.*TDTW
C     RETURN
C     END
      END
