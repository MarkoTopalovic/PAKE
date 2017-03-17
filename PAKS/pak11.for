C=======================================================================
C
C   SUBROUTINE UL1EGL
C              UL1EK
C              UL1EK2
C              LIMMH1
C              TGRAF1
C
C=======================================================================
      SUBROUTINE UL1EGL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM ZA ULAZNIE PODATAKE STAPOVA
CE     MAIN PROGRAM FOR INPUT DATA ABOUT 1/D ELEMENTS
C
      include 'paka.inc'
      
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
C
CS     REPERI U VEKTORU ZA ULAZNE PODATKE
CE     POINTERS IN INPUT VECTOR
C
      NCVE=2
      LAPRS=1
      LTBTH=LAPRS+NE*IDVA
      LTDTH =LTBTH
      IF(INDBTH.EQ.1) LTDTH=LTBTH+NE*IDVA
      LNMAT =LTDTH
      IF(INDDTH.EQ.1) LNMAT=LTDTH+NE*IDVA
      LISNA=LNMAT+NE
      LNEL=LISNA+NE
      LCEL=LNEL
      IF(ICVEL.EQ.1) LNEL=LCEL+NE
      LAU=LMAX
C
CS     POZIVANJE PROGRAMA ZA ULAZNE PODATKE U VEKTOR AU
CE     CALL ROUTINES FOR INPUT DATA IN VECTOR  AU
C
      CALL UL1EK(A(LAU))
      NCVE3=NCVE*3
      MXAE = NCVE3+(NCVE3+NCVE3*(NCVE3+1)/2)*IDVA
      LMAX = LMAX + MXAE
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UL1EK(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
CE     MENAGEMENT PROGRAM FOR INPUT DATA IN  AU  VECTOR
C
      include 'paka.inc'
      
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /INCONF/ LINDBEL,LBIRTHC
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /RESTAR/ TSTART,IREST
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION AU(*)
      REAL AU
C
CS     POZIVANJE PROGRAMA ZA ULAZNE PODATKE O STAPOVIMA
CE     CALL ROUTINE FOR INPUT DATA ABOUT TRUSSES
C
      LL=NE*NCVE
      CALL ICLEAR(AU(LNEL),LL)
      NMA=0
      NMI=0
      CALL UL1EK2(AU(LNEL),AU(LNMAT),AU(LISNA),AU(LAPRS),
     1            AU(LTBTH),AU(LTDTH),AU(LCEL),NMA,NMI,ICVEL)
      IF(NBLGR.GE.0) 
     +          CALL TGRAF1(AU(LNEL),AU(LCEL),ICVEL,AU(LNMAT),AU(LAPRS))
      IF(NBLGR.GE.0) 
     +       CALL TGRAU1(AU(LNEL),AU(LCEL),ICVEL,AU(LNMAT),AU(LAPRS),49)
      IF(ICVEL.EQ.1) CALL VEZACE(AU(LNEL),A(LELCV),NE,NCVE)
C
      NCVE3=NCVE*3
      LLMEL=LNEL+NE*NCVE
      LELC=LLMEL
      IF(ICVEL.EQ.1) THEN
        LLMEL=LELC+NMA-NMI+1
        CALL ICLEAR(AU(LELC),NMA-NMI+1)
        CALL VEZAEL(AU(LCEL),AU(LELC),NE,NMI)
      ENDIF
      LMXAU=LLMEL+NE*NCVE3
      CALL DELJIV(LMXAU,2,INDL)
      IF(INDL.EQ.0) LMXAU=LMXAU+1
CE    pointer for initial configuration for birth elements
      LBIRTHC=LMXAU
      LINDBEL=LMXAU
      IF(INDBTH.EQ.1) THEN
         LINDBEL=LBIRTHC+NE*NCVE*3*IDVA
         LMXAU=LINDBEL+NE
      ENDIF
      MXAU = LMXAU - 1
      CALL DELJIV(MXAU,2,INDL)
      IF(INDL.EQ.1) MXAU=MXAU+1
      LMAX=LAU+MXAU
      IF (LMAX.LE.MTOT) GO TO 5
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) LMAX,MTOT
      STOP
C
CS     FORMIRANJE VEKTORA LM
CE     FORM  LM  VECTOR
C
    5 LL=NE*NCVE3
      CALL ICLEAR(AU(LLMEL),LL)
      CALL LMIMH1(A(LID),AU(LNEL),AU(LLMEL),NCVE3)
C
      LMAX8=LMAX8+1
      WRITE(IELEM,REC=LMAX8)
     1 NCVE,ITERME,MXAU,LNEL,LNMAT,LAPRS,LISNA,LLMEL,LCEL,LELC,NMA,NMI,
     1 INDBTH,INDDTH,LTBTH,LTDTH,LINDBEL,LBIRTHC
      CALL WRITDD(AU(LAPRS),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      IF(IATYP.EQ.0) GO TO 10
      IF(NMODM.LE.4.OR.NMODM.EQ.11.OR.NMODM.EQ.12) GO TO 10
CS....  PLASTICNOST
CS....  PROSTOR ZA VELICINE PRETHODNOG KORAKA
CE     PLASTICITY
      LPLAST=LMAX
      NPROS=NE*MODPR1( NMODM )*IDVA
      LPLAS1=LPLAST+NPROS
      CALL DELJIV(LPLAS1,2,INDL)
      IF(INDL.EQ.0) LPLAS1=LPLAS1+1
CS..... PROSTOR ZA VELICINE TEKUCEG KORAKA
      NPROS1 = NPROS
      LMAX=LPLAS1+NPROS1
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      IF(LMAX.GT.MTOT) CALL ERROR(1)
CS....  INICIJALIZACIJA
      IF(IREST.NE.1) THEN
         CALL INCPL1(AU(LNMAT),NE,NMODM,LMODEL)
         CALL WRITDD(A(LPLAST),NPROS/IDVA,IELEM,LMAX8,LDUZI)
         CALL WRITDD(A(LPLAS1),NPROS1/IDVA,IELEM,LMAX8,LDUZI)
      ELSE
         CALL READDD(A(LPLAST),NPROS/IDVA,IELEM,LMAX8,LDUZI)
         CALL READDD(A(LPLAS1),NPROS1/IDVA,IELEM,LMAX8,LDUZI)
      ENDIF
      RETURN
C
   10 LSIGMA=LMAX
      NPROS=2*NE*IDVA
      LMAX=LSIGMA+NPROS
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      IF(IREST.NE.1) THEN
         CALL CLEAR(A(LSIGMA),NPROS/IDVA)
         CALL WRITDD(A(LSIGMA),NPROS/IDVA,IELEM,LMAX8,LDUZI)
      ELSE
         CALL READDD(A(LSIGMA),NPROS/IDVA,IELEM,LMAX8,LDUZI)
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2005 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA PODATKE O STAPOVIM
     1A'/' POTREBNA DIMENZIJA , LMAX=',I10/' RASPOLOZIVA DIMENZIJA,
     2NTOT =',I10)
 6005 FORMAT(///' NOT ENOUGH SPACE IN WORKING VECTOR  A'
     1/' REQUESTED DIMENSION , LMAX=',I10
     2/' AVAILABLE DIMENSION , MTOT=',I10)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE INCPL1(NMAT,NE,NMODM,LMODEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C....  INICIJALIZACIJA PLASTICNOSTI
      COMMON /ITERBR/ ITER
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      DIMENSION NMAT(*)
C
      ITER=0
      IRAC=0
      TGT =0.D0
      STRAIN =0.D0
      STRESS =0.D0
      CALL UCITAM(A(LMODEL),NMODM)
      DO 50 NLM=1,NE
      MAT=NMAT(NLM)
   50 CALL MODMA1(STRAIN,STRESS,TGT,NMODM,NLM,IRAC)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UL1EK2(NEL,MATV,ISNAP,APRS,
     1                  TBTH,TDTH,MCVEL,NMA,NMI,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PODPROGRAM ZA UCITAVANJE PODATAKA O STAPOVIMA
CE     SUBROUTINE FOR READING DATA ABOUT TRUSSES
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /STAPRO/ THIDP(100),NPROPR(100)
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION NEL(NE,*),MATV(*),ISNAP(*),APRS(*),
     1          TBTH(*),TDTH(*),MCVEL(*)
C
CS     P  O  D  A  C  I      O      E  L  E  M  E  N  T  I  M  A
CE     D  A  T  A      A  B  O  U  T      E  L  E  M  E  N  T  S
C
      CALL ICLEAR(NPROPR,100)
      IF(ISRPS.EQ.0.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,6000)
      NMATS=  1
      ISNAA=  0
      APRSS=  0.D0
      I = 0
      NAUT=0
    5 I = I + 1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NN,NMAT,ISNA,APR,KORC,BTH,DTH
      IF(INDFOR.EQ.2.and.ind56.eq.0)
     1READ(ACOZ,1010) NN,NMAT,ISNA,APR,KORC,BTH,DTH
      IF(INDFOR.EQ.2.and.ind56.eq.1)
     1READ(ACOZ,4010) NN,NMAT,ISNA,APR,KORC,BTH,DTH
      IF(NPROPR(NMAT).EQ.0) THEN
         IF(NMAT.GT.100) STOP 'NMAT.GT.100'
         NPROPR(NMAT)=NMAT
      ENDIF
      IF(ICVEL.EQ.1) THEN
         MCVEL(I)=NN
         NI=NN
         IF(I.EQ.1) THEN
            NMA=NN
            NMI=NN
         ELSE
            IF(NMA.LT.NN) NMA=NN
            IF(NMI.GT.NN) NMI=NN
         ENDIF
         NN=I
      ENDIF
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (NEL(NN,J),J=1,2)
      IF(INDFOR.EQ.2.and.ind56.eq.0)
     1READ(ACOZ,1000) (NEL(NN,J),J=1,2)
      IF(INDFOR.EQ.2.and.ind56.eq.1)
     1READ(ACOZ,4000) (NEL(NN,J),J=1,2)
      IF(NMAT.EQ.0) NMAT=NMATS
      NMATS=NMAT
      MATV(NN)= NMAT
      IF(ISNA.EQ.0) ISNA=ISNAA
      ISNAA=ISNA
      IF(ISNA.LT.0) ISNA = 0
      ISNAP(NN)=ISNA
      IF(DABS(APR).LT.1.0D-10) APR=APRSS
      APRSS=APR
      APRS(NN)=APR
      IF(INDBTH.EQ.1) TBTH(NN)=BTH
      IF(INDDTH.EQ.1) TDTH(NN)=DTH
      IF(ICVEL.EQ.1) IN=MCVEL(NN)
      IF(ICVEL.EQ.0) IN=NN
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)THEN
         WRITE(IZLAZ,5001) IN,(NEL(NN,J),J=1,2),NMAT,ISNA,APR,KORC
         IF(INDBTH.EQ.1) THEN
         IF(DABS(TBTH(NN)).GT.1.0D-10) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) TBTH(NN)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) TBTH(NN)
         ENDIF
         ENDIF
         IF(INDDTH.EQ.1) THEN 
         IF(DABS(TDTH(NN)).GT.1.0D-10) THEN 
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TDTH(NN)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TDTH(NN)
         ENDIF
         ENDIF
      ENDIF
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(I.EQ.NE) GO TO 50
      GO TO 5
C
   20 NAUT=1
      IF(ICVEL.EQ.1) N1=NI
      IF(ICVEL.EQ.0) NN1=NN
      KORA=KORC
      GO TO 5
C
CS     AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
CE     AUTOMATIC  GENERATION BETWEN NODES  N1  AND  N2
C
   30 NN2=NN
      IF(ICVEL.EQ.1) THEN
         NN1=NN-1
         N2=NI
      ENDIF
      N11=NEL(NN1,1)
      N22=NEL(NN2,1)
      N12=N22-N11
      CALL DELJIV(N12,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      N11=NEL(NN1,2)
      N22=NEL(NN2,2)
      N21=N22-N11
      IF(N12.NE.N21) GO TO 100
      IF(ICVEL.EQ.1) NNN=N2-N1-1
      IF(ICVEL.EQ.0) NNN=NN2-NN1-1
      NGG=N12/KORA-1
      IF(NNN.NE.NGG) GO TO 150
      IAUT=N12/KORA-1
      IF(ICVEL.EQ.1) THEN
         NM=NN+IAUT
         MCVEL(NM)=MCVEL(NN)
         MATV(NM) = MATV(NN)
         ISNAP(NM)=ISNAP(NN)
         APRS(NM)=APRS(NN)
         IF(INDBTH.EQ.1) TBTH(NM) = TBTH(NN)
         IF(INDDTH.EQ.1) TDTH(NM) = TDTH(NN)
         DO 32 II=1,2
            NEL(NM,II)=NEL(NN,II)
   32    CONTINUE
      ENDIF
      DO 34 J=1,IAUT
         JJ=NN1+J
         IF(ICVEL.EQ.1) THEN
            N1=N1+1
            MCVEL(JJ)=N1
         ENDIF
         MATV(JJ) = MATV(NN1)
         ISNAP(JJ)=ISNAP(NN1)
         APRS(JJ)=APRS(NN1)
         IF(INDBTH.EQ.1) TBTH(JJ) = TBTH(NN1)
         IF(INDDTH.EQ.1) TDTH(JJ) = TDTH(NN1)
         DO 35 II=1,2
            NODP=NEL(JJ-1,II)
            NEL(JJ,II)=NODP+KORA
            IF(NODP.EQ.0) NEL(JJ,II) = 0
   35    CONTINUE
   34 CONTINUE
      IF(ICVEL.EQ.1) I=NM
      IF(ICVEL.EQ.0) I=I+IAUT
      IF(I.EQ.NE) GO TO 50
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      IF(ICVEL.EQ.1) N1=N2
      IF(ICVEL.EQ.0) NN1=NN2
      GO TO 5
C
CS     STAMPANJE UCITANIH I GENERISANIH PODATAKA
CE     PRINT INPUT AND GENERATED DATA
C
   50 IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2140)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6140)
      DO 70 I=1,NE
         IF(ICVEL.EQ.1) IN=MCVEL(I)
         IF(ICVEL.EQ.0) IN=I
         WRITE(IZLAZ,5001) IN,(NEL(I,J),J=1,2),MATV(I),ISNAP(I),APRS(I)
         IF(INDBTH.EQ.1) THEN
         IF(DABS(TBTH(I)).GT.1.0D-10) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) TBTH(I)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) TBTH(I)
         ENDIF
         ENDIF
         IF(INDDTH.EQ.1) THEN 
         IF(DABS(TDTH(I)).GT.1.0D-10) THEN 
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2050) TDTH(I)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6050) TDTH(I)
         ENDIF
         ENDIF
   70 CONTINUE
      RETURN
C
  150 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2110) NN1,NN2,NGG
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6110) NN1,NN2,NGG
      STOP
  100 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) N22,N11,KORA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) N22,N11,KORA
      STOP
C
 1000 FORMAT(14I5)
 4000 FORMAT(14I10)
 1010 FORMAT(3I5,F10.3,I5,2F10.0)
 4010 FORMAT(I10,2I5,F10.3,I5,2F10.0)
 5001 FORMAT(1X,I5,4X,2I6,3X,I6,2X,I6,14X,1PD10.3,I6)
C-----------------------------------------------------------------------
 2000 FORMAT(///6X,'U C I T A N I    P O D A C I    O    S T A P O V I M
     1 A'/6X,54('-')///
     111X,'(MOGUCE GENERISANJE PODATAKA O ELEMENTIMA,'/
     212X,'UCITATI POTREBAN BROJ KARTICA ZA DEF. SVIH ELEMENATA)'///1X,
     3'ELEMENT','   CVOR1 CVOR2      ','MAT.',4X,
     4'NAPON',5X,'     ',5X,'POVRSINA',4X,'KORAK'/2X,'BROJ',3X,
     4'             ',4X,'BROJ',4X,'STAMPA')
 2100 FORMAT(///' BROJ CVORA N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJ
     1A CVORA N1=',I5,' I KONACNOG BROJA KORAKA KORA=',I5)
 2110 FORMAT(///' IZMEDJU ELEMENATA N1=',I5,' I N2=',I5,' NEMA NG=',I5,
     1' ELEMENATA KOJE TRBA GENERISATI')
 2140 FORMAT(6X,'G E N E R I S A N I    P O D A C I    O    S T A P O V
     1I M A'/6X,60('-')///1X,
     3'ELEMENT','   CVOR1 CVOR2      ','MAT.',4X,
     4'NAPON',5X,'     ',5X,'POVRSINA'           /2X,'BROJ',3X,
     4'             ',4X,'BROJ',4X,'STAMPA')
 2040 FORMAT(/' VREME NASTAJANJA ELEMENTA - TBTH =',1PD10.3)
 2050 FORMAT(/' VREME NESTAJANJA ELEMENTA - TDTH =',1PD10.3)
C-----------------------------------------------------------------------
 6000 FORMAT(///6X,'I N P U T  D A T A   F O R  T R U S S  E L E M E N T
     1 S'/6X,54('-')///1X,
     3'ELEMENT','   NODE1 NODE2      ','MAT.',4X,
     4'STRESS',5X,'    ',5X,'AREA',4X,' STEP'/2X,'NUMBER',3X,
     4'            ',4X,'NUMBER',4X,'PRINT ')
 6100 FORMAT(///' NODE NUMBER N2=',I5,' CAN NOT BE DEFINED BY SUPERPOSIT
     1ION OF NODE N1=',I5,' WITH A NUMBER OF STEPS  KORA=',I5)
 6110 FORMAT(///' AMONG ELEMENTS N1=',I5,' AND N2=',I5,' THERE IS NO  ',
     1'NG=',I5,' ELEMENTS TO GENERATE')
 6140 FORMAT(6X,'G E N E R A T E D  D A T A  F O R  T R U S S  E L E M E
     1 N T S'/6X,61('-')///1X,
     3'ELEMENT','   NODE1 NODE2      ','MAT.',4X,
     4'STRESS',5X,'    ',5X,'AREA    '/2X,'NUMBER',3X,
     4'            ',4X,'NUMBER',4X,'PRINT ')
 6040 FORMAT(/' ELEMENT BIRTH TIME - TBTH =',1PD10.3)
 6050 FORMAT(/' ELEMENT DEATH TIME - TDTH =',1PD10.3)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE LMIMH1(ID,NEL,LMEL,NCVE3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE VEKTORA LM I VISINA STUBOVA (MHT)
CE     OBTAIN VECTOR  LM   AND COLUMN HEIGHTS
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
C
      DIMENSION ID(NP,*),NEL(NE,*),LMEL(NCVE3,*)
C     PETLJA PO ELEMENTIMA
      DO 100 NLM=1,NE
      DO 10 NC=1,NCVE
      IF(NEL(NLM,NC).EQ.0) GO TO 10
      NNC=3*(NC-1)
      JJ=NEL(NLM,NC)
      DO 20 I=1,3
   20 LMEL(NNC+I,NLM)=ID(JJ,I)
   10 CONTINUE
C     FORMIRANJE VISINA STUBOVA
      CALL VISINE(A(LMHT),NCVE3,LMEL(1,NLM))
      IF (IABS(ICCGG).EQ.1) THEN
         ND=NCVE3
         WRITE(IDRAKCE) ND,(LMEL(I,NLM),I=1,ND)
         NELUK=NELUK+1
      ENDIF
  100 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE TGRAF1(NEL,MCVEL,ICVEL,NMAT,THID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PROGRAM ZA STAMPANJE STAPA  U UNIVERZALNI FILE
CE     PROGRAM FOR PRINTOUT DATA IN UNIVERSAL GRAPHICS FILE
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /NIDEAS/ IDEAS
C
      DIMENSION NEL(NE,*),MCVEL(*),NMAT(*),THID(*)
      DIMENSION FIZ(40)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TGRAF1'
      IF(ideas.eq.-1) return
      ISUMGR=ISUMGR+1
C
C     FIZICKE OSOBINE
C
      NULA=0
      ZERO=0.
      ONE=1.
      JEDAN=1
      INDPR=25
      INDPD=2
      I11=11
      I2=2
      I8=8
      I14=14
C
      CALL CLEAR(FIZ,40)
      ICL=1
      IPH=10
      NPR=1
      FIZ(1)=THID(1)
      IND=-1
      IF(IDEAS.GT.6) THEN
        ITYP=2437
      ELSE
        ITYP=772
      ENDIF
      IF(IDEAS.GT.6) THEN
         WRITE(IGRAF,1100) IND
         WRITE(IGRAF,1100) ITYP
         WRITE(IGRAF,1000) ISUMGR,IPH,NPR
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2000) ISUMGR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6000) ISUMGR
         WRITE(IGRAF,1000) INDPR,INDPD,JEDAN
         WRITE(IGRAF,1300) ZERO
C         WRITE(IGRAF,1200) (FIZ(I),I=ICL,18)
         WRITE(IGRAF,1100) IND
C
         ITYP=776
         WRITE(IGRAF,1100) IND
         WRITE(IGRAF,1100) ITYP
         WRITE(IGRAF,1000) JEDAN,NULA,NULA
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2000) ISUMGR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6000) ISUMGR
         WRITE(IGRAF,1200) (FIZ(I),I=31,40)
         WRITE(IGRAF,1200) (FIZ(I),I=ICL,40)
         WRITE(IGRAF,1000) I11,I2,I8,I14,JEDAN,JEDAN
         WRITE(IGRAF,1001) NULA,JEDAN,JEDAN,I11,ONE
         WRITE(IGRAF,1100) IND
      ENDIF
C
C     M A T E R I J A L I
C
      MAT=NMAT(1)
      IF(IDEAS.LT.7) CALL MELGR(MAT)
C
C     E L E M E N T I
C
      IF(ISRPS.EQ.0.AND.NCVE.GT.2) WRITE(IZLAZ,2200) NGE
      IF(ISRPS.EQ.1.AND.NCVE.GT.2) WRITE(IZLAZ,6200) NGE
C     GRAFICKI OPIS STAPA: SA 2 CVORA = 1
      IFGD=1
C     VRSTA  ELEMENTA: STAP SA 2 CVORA  = 11
      IFDI=11
C     TABELA FIZICKIH OSOBINA
      IPTN=ISUMGR
C     TABELA MATERIJALA
      IF(IDEAS.LT.7) THEN
         MPTN=ISUMGR
      ELSE
         MPTN=MAT
      ENDIF         
C     BOJA  
      ICOL=7
C     BROJ CVOROVA NA ELEMENTU
      NNODS=2
      IND=-1
      IF(IDEAS.GT.6) THEN
         ITYP=2412
      ELSEIF(IDEAS.EQ.6) THEN
         ITYP=780
         ICOL=11
      ELSE
         ITYP=71
      ENDIF
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      JEDAN=1
      NULA=0
      DO 10 I=1,NE
C        REDNI BROJ ELEMENTA
         IEL=I+ISUMEL
         IF(ICVEL.EQ.1) IEL=MCVEL(I)
         IF(IDEAS.GT.6) THEN
            WRITE(IGRAF,1000) IEL,IFDI,IPTN,MPTN,ICOL,NNODS
            WRITE(IGRAF,1000) NULA,JEDAN,JEDAN
         ELSEIF(IDEAS.EQ.6) THEN
            WRITE(IGRAF,1000) IEL,IFDI,JEDAN,IPTN,MPTN,JEDAN,ICOL,NNODS
         ELSE
            WRITE(IGRAF,1000) IEL,IFGD,IFDI,IPTN,MPTN,ICOL,NNODS
         ENDIF
         WRITE(IGRAF,1000) (NEL(I,J),J=1,2)
   10 CONTINUE
      WRITE(IGRAF,1100) IND
      ISUMEL=ISUMEL+NE
      RETURN
C
 1000 FORMAT(8I10)
 1001 FORMAT(4I10,1PE13.5)
 1100 FORMAT(I6)
 1200 FORMAT(6(1PE13.5))
 1300 FORMAT(1PE15.7)
C-----------------------------------------------------------------------
 2000 FORMAT('STAP (',I3,')')
 2200 FORMAT(//' PROGRAM ZA GRAFICKO PRIKAZIVANJE REZULTATA "IDEAS"'/
     1' ZAHTEVA STAP SA 2 CVORA PO ELEMENTU U GRUPI ELEMENATA NGE ='
     1,I5)
C-----------------------------------------------------------------------
 6000 FORMAT('ROD (',I3,')')
 6200 FORMAT(//' GRAPHIC PACKAGE   "IDEAS"'/
     1' PERMITS ONLY ROD WITH 2 NODES PER ELEMENT IN GROUP',
     1'  NGE =',I5)
C-----------------------------------------------------------------------
      END