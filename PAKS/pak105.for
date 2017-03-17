C=======================================================================
C
C=======================================================================
      SUBROUTINE MIKA(NPODS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER TIME PERIODS AND STEPS
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /GRUPER/ LIGRUP
      COMMON /DSTAZE/ NSTAZ
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /LANCZO/ LA,LB,LTT,LW
      COMMON /RADIZA/ INDBG
      COMMON /CDEBUG/ IDEBUG
      INCLUDE 'mpif.h'
      INTEGER ierr,myid
      DIMENSION NPODS(JPS1,*)
      DIMENSION CJ(1000)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      write(*,*) 'iz mike iz procesa',myid
C
      IF (myid.eq.0) THEN
         IF(IDEBUG.GT.0) PRINT *, ' MIKA'
C
         IST=0
         IF(IST.EQ.1) PRINT *,'MMAXA,LRAD,LMAX,LID',JMAXA,LRAD,LMAX,LID
         LRAD=JMAXA
      ENDIF
C
C        INVERZIJA MATRICE KRUTOSTI PODSTRUKTURE
C
         DO 250 JPBB=1,JPS
            IF (myid.eq.0) THEN
               JPBR=JPBB 
               CALL PODDAT(NPODS,1)
               CALL PODDAT(NPODS,3)
               LSK=LRAD
               LRTDT=LSK+NWK*IDVA
               LMAX=LRTDT+JEDN*IDVA
               CALL RSTAZK(NPODS,LSK,35)
C              CALL IWRR(A(LMAXA),JEDN,'MAX1')
C              CALL WRR(A(LSK),NWG,'SK1 ')
C              LDLT - Kj
            ENDIF
            CALL RESEN(A(LSK),A(LRTDT),A(LMAXA),JEDN,1)
            IF (myid.eq.0) THEN
               CALL WSTAZK(NPODS,LSK,60)
               LRAD=LID
            ENDIF
  250    CONTINUE
C
C     VADJENJE Kbb IZ Kab+Kbb
C
      IF (myid.eq.0) THEN
         DO 300 JPBR=1,JPS
            CALL PODDAT(NPODS,1)
            LSK=LRAD
            IF(IST.EQ.1) WRITE(3,*) 'JPBR,LRAD',JPBR,LRAD
            LSKG=LSK+(NWK-NWP)*IDVA
            JED=JEDN-JEDNP
            NWKP=JED*(JED+1)/2
            LIGRUP=LSKG+NWKP*IDVA
            LMAXA=LIGRUP+NGELEM*5
            LMAX=LMAXA+JEDN+1
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            NP6=NGELEM*5+JEDN+1 
            LMAX13=NPODS(JPBR,12)-1
            CALL IREADD(A(LIGRUP),NP6,IPODS,LMAX13,LDUZI)
            LMAX13=NPODS(JPBR,61)-1
            CALL READDD(A(LSK),NWK-NWP,IPODS,LMAX13,LDUZI)
            CALL TROUGO(A(LSK),A(LSKG),A(LMAXA+JEDNP),JED,NWP)
            LMAX13=NPODS(JPBR,37)-1
C           Kbb
            CALL WRITDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
            IF(IST.EQ.1) CALL WRR(A(LSK),NWK-NWP,'SK  ')
            IF(IST.EQ.1) CALL WRR(A(LSKG),NWKP,'SKG ')
  300    CONTINUE
      ENDIF
C
C     FORMIRANJE MATRICE  Kg
C
      IF (myid.eq.0) THEN
         LRAD=JMAXA
         JPBR=JPS1
         CALL PODDAT(NPODS,1)
         CALL PODDAT(NPODS,3)
C
         LSK=LRAD
         LRAD=LSK+NWG*IDVA
         CALL CLEAR(A(LSK),NWG)
C
         DO 400 JPBR=1,JPS
            JEDN=NPODS(JPBR,6)
            JEDNP=NPODS(JPBR,23)
            JED=JEDN-JEDNP
            NWKP=JED*(JED+1)/2
            LSKG=LRAD
            LLMG=LSKG+NWKP*IDVA
            LMAX=LLMG+JED
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            LMAX13=NPODS(JPBR,37)-1
            CALL READDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
            LMAX13=NPODS(JPBR,27)-1
            CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
            CALL SPAKUJ(A(LSK),A(JMAXA),A(LSKG),A(LLMG),JED)
  400    CONTINUE
         JPBR=JPS1
         JEDN=JEDNG
C        LDLT - Kg
         IF(IST.EQ.1) CALL WRR(A(LSK),NWG,'KG1 ')
C         CALL IWRR(A(JMAXA),JEDN,'MAXA')
C         CALL WRR(A(LSK),NWG,'A   ')
      ENDIF
      CALL RESEN(A(LSK),A(LRTDT),A(JMAXA),JEDN,1)
      IF (myid.eq.0) THEN
         IF(IST.EQ.1)  CALL WRR(A(LSK),NWG,'KGIN')
         CALL WSTAZK(NPODS,LSK,60)
         LRAD=JMAXA
         IF(IST.EQ.1) WRITE(3,*) 'LRAD100',LRAD
C
         IF(NSOPV.GT.1000) THEN
            STOP '(NSOPV.GT.1000), PROGRAM STOP'
         ELSE
            DO 10 I=1,NSOPV
   10       CJ(I)=0.
         ENDIF
      ENDIF
C
C        FORMIRANJE RITZOVIH VEKTORA
C  
      DO 100 IX=1,NSOPV
C
         IF (myid.eq.0) NST=NSTAZ*(IX-1)
C
C        5.   REDUKCIJA DESNE STRANE PODSTRUKTURE - Fj(k)
C
         DO 200 JPBB=1,JPS
            IF (myid.eq.0) THEN
               JPBR=JPBB 
               CALL PODDAT(NPODS,1)
               CALL PODDAT(NPODS,3)
               LSK=LRAD
               LRTDT=LSK+NWK*IDVA
               LMAX=LRTDT+JEDN*IDVA
C              LDLT - Kg
               CALL RSTAZK(NPODS,LSK,60)
               LMAX13=NPODS(JPBR,39)-1+NST
               CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               IF(IX.EQ.1) THEN
                  CALL POCETI(A(LRTDT),1)
               ENDIF
               IF(IST.EQ.1) CALL WRR(A(LRTDT),JEDN,'POCE')
            ENDIF
      IF(IST.EQ.1) CALL WRR(A(LRTDT),JEDN,'POCE')
            CALL RESEN(A(LSK),A(LRTDT),A(LMAXA),JEDN,2)
            IF (myid.eq.0) THEN
               IF(IST.EQ.1) CALL WRR(A(LRTDT),JEDN,'RDES')
               LMAX13=NPODS(JPBR,39)-1+NST
               CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               LRAD=LID
            ENDIF
  200    CONTINUE
C
         IF (myid.eq.0) THEN
            JPBR=JPS1
            CALL PODDAT(NPODS,1)
            CALL PODDAT(NPODS,3)
            LSK=LRAD
            LRTDT=LSK+NWG*IDVA
            LRAD=LRTDT+JEDNG*IDVA
            CALL CLEAR(A(LRTDT),JEDNG)
C
            DO 410 JPBR=1,JPS
               JEDN=NPODS(JPBR,6)
               JEDNP=NPODS(JPBR,23)
               JED=JEDN-JEDNP
               LRTG=LRAD
               LLMG=LRTG+JEDN*IDVA
               LMAX=LLMG+JED
               IF(LMAX.GT.MTOT) CALL ERROR(1)
               LMAX13=NPODS(JPBR,39)-1+NST
               CALL READDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,27)-1
               CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
               LRTG=LRTG+JEDNP*IDVA
C              Fg(k)
               CALL SPAKUD(A(LRTDT),A(LRTG),A(LLMG),JED)
  410       CONTINUE
C
            JPBR=JPS1
            JEDN=JEDNG
C           LDLT - Kg
            CALL RSTAZK(NPODS,LSK,60)
C
CE          BACKSUBSTITUTION - SOLUTION OF SYSTEM EQUATIONS
CS          ZAMENA UNAZAD - RESAVANJE :SISTEMA JEDNACINA
C
C           6.           Kg . Xg(k) = Fg(k)
            IF(IST.EQ.1) CALL WRR(A(LSK),NWG,'KG  ')
            IF(IST.EQ.1) CALL WRR(A(LRTDT),JEDN,'RTDT')
            IF(IST.EQ.1) CALL IWRR(A(JMAXA),JEDN+1,'MAXA')
         ENDIF
         CALL RESEN(A(LSK),A(LRTDT),A(JMAXA),JEDN,2)
         IF (myid.eq.0) THEN
            IF(IST.EQ.1) CALL WRR(A(LRTDT),JEDN,'UTDT')
C
C
CE          L I N E A R     A N A L Y S I S
CS          L I N E A R N A      A N A L I Z A
C
            LID=LRAD
            LCVEL=LID+NP*6
            NPM=NPA-NPI+1
            LELCV=LCVEL+NP
            LMAX=LELCV+NPM
            NP6=NP*7+NPM
            LMAX13=NPODS(JPBR,2)-1
            CALL IREADD(A(LID),NP6,IPODS,LMAX13,LDUZI)
C
CE          PRINT CONTOUR DISPLACEMENTS 
CS          STAMPANJE KONTURNIH POMERANJA
C
C           CALL STAMPA
C
            DO 600 JPBR=1,JPS
               JEDN=NPODS(JPBR,6)
               JEDNP=NPODS(JPBR,23)
               JED=JEDN-JEDNP
               LRTG=LRAD
               LLMG=LRTG+JEDN*IDVA
               LMAX=LLMG+JED
               IF(LMAX.GT.MTOT) CALL ERROR(1)
               LMAX13=NPODS(JPBR,27)-1
               CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,39)-1+NST
               CALL READDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
               LRTGG=LRTG+JEDNP*IDVA
               CALL SPAKUU(A(LRTDT),A(LRTGG),A(LLMG),JED)
               LMAX13=NPODS(JPBR,39)-1+NST
C              MESOVITI VEKTOR F I X
               CALL WRITDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
  600       CONTINUE
            LRAD=JMAXA
         ENDIF
C
C        7.   NALAZANJE Xj(k)**
C
         DO 800 JPBB=1,JPS
C
            IF (myid.eq.0) THEN
               JPBR=JPBB
               CALL PODDAT(NPODS,1)
               CALL PODDAT(NPODS,3)
               LSK=LRAD
               LRTDT=LSK+NWK*IDVA
               LMAX=LRTDT+JEDN*IDVA
C              LDLT - Kg
               CALL RSTAZK(NPODS,LSK,60)
               LMAX13=NPODS(JPBR,39)-1+NST
               CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
            ENDIF
            CALL RESEN(A(LSK),A(LRTDT),A(LMAXA),JEDN,3)
            IF (myid.eq.0) THEN
               LMAX13=NPODS(JPBR,39)-1+NST
C              7.    Xj(k)**
               CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               IF(IX.EQ.1) THEN
                  LMAX13=NPODS(JPBR,52)-1
C                 7.    Xj(k)*
                  CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               ENDIF
C
CE             R E S U L T S
CS             R E Z U L T A T I
C
CE             PRINT DISPLACEMENTS AND STRESSES
CS             STAMPANJE POMERANJA I NAPONA
C
C               CALL STAMPA
C
               LRAD=LID
            ENDIF
  800    CONTINUE
C
C           ORTOGONALIZACIJA
C        9.    RACUNANJE KOEFICIJENATA C(1) DO C(k-1)
C
         IF (myid.eq.0) THEN
            IF(IX.GT.1) THEN
               DO 700 JPBB=1,JPS
C
                  JPBR=JPBB
                  CALL PODDAT(NPODS,1)
                  CALL PODDAT(NPODS,3)
                  LRTDT=LRAD
                  LRM=LRTDT+JEDN*IDVA
                  LMAX=LRM+JEDN*IDVA
                  LMAX13=NPODS(JPBR,39)-1+NST
                  CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               DO 700 JK=1,IX-1
                  LMAX13=NPODS(JPBR,39)-1+NSTAZ*(JK-1)
                  CALL READDD(A(LRM),JEDN,IPODS,LMAX13,LDUZI)
                  CALL SCALAR(A(LRM),A(LRTDT),CE,JEDN)
                  CJ(JK)=CJ(JK)-CE
C
                  LRAD=LID
  700          CONTINUE
               DO 710 JPBB=1,JPS
C
                  JPBR=JPBB
                  CALL PODDAT(NPODS,1)
                  CALL PODDAT(NPODS,3)
                  LSK=LRAD
                  LRTDT=LSK+NWK*IDVA
                  LRM=LRTDT+JEDN*IDVA
                  LMAX=LRM+JEDN*IDVA
                  LMAX13=NPODS(JPBR,39)-1+NST
                  CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
                  DO 720 JK=1,IX-1
                     LMAX13=NPODS(JPBR,52)-1+NSTAZ*(JK-1)
                     CALL READDD(A(LRM),JEDN,IPODS,LMAX13,LDUZI)
                     CALL ZBIRM(A(LRTDT),A(LRM),CJ(JK),JEDN)
  720             CONTINUE
                  LMAX13=NPODS(JPBR,52)-1+NST
                  CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
                  CALL CLEAR(A(LRM),JEDN)
C                 M - 
                  CALL RSTAZK(NPODS,LSK,54)
C                 CALL IWRR(A(LMAXA),JEDN,'MAX ')
C                 CALL WRR(A(LSK),NWK,'M   ')
                  LMAX13=NPODS(JPBR,39)-1+NST
                  CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
                  CALL MAXAPR(A(LSK),A(LRTDT),A(LRM),A(LMAXA),JEDN)
                  LMAX13=NPODS(JPBR,39)-1+NST
                  CALL WRITDD(A(LRM),JEDN,IPODS,LMAX13,LDUZI)
C
                  LRAD=LID
  710          CONTINUE
            ENDIF
C
C           NORMALIZACIJA
C           11.    RACUNANJE ALFA
C
            ALFA=0.D0
            DO 810 JPBB=1,JPS
C
               JPBR=JPBB
               CALL PODDAT(NPODS,1)
               CALL PODDAT(NPODS,3)
               LSK=LRAD
               LRTDT=LSK+NWK*IDVA
               LRM=LRTDT+JEDN*IDVA
               LMAX=LRM+JEDN*IDVA
               CALL CLEAR(A(LRM),JEDN)
C              M - 
               CALL RSTAZK(NPODS,LSK,54)
               LMAX13=NPODS(JPBR,52)-1+NST
               CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               CALL MAXAPR(A(LSK),A(LRTDT),A(LRM),A(LMAXA),JEDN)
               CALL SCALAR(A(LRM),A(LRTDT),ALF,JEDN)
               ALFA=ALFA+ALF
C
               LRAD=LID
  810       CONTINUE
            ALFA=1.D0/DSQRT(ALFA)
C
C           13.  NALAZENJE  Xj(k)
C
            DO 820 JPBB=1,JPS
C
               JPBR=JPBB
               CALL PODDAT(NPODS,1)
               CALL PODDAT(NPODS,3)
               LSK=LRAD
               LRTDT=LSK+NWK*IDVA
               LRM=LRTDT+JEDN*IDVA
               LMAX=LRM+JEDN*IDVA
C              M - 
               CALL RSTAZK(NPODS,LSK,54)
               LMAX13=NPODS(JPBR,52)-1+NST
               CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               CALL JEDNAK(A(LRM),A(LRTDT),ALFA,JEDN)
               LMAX13=NPODS(JPBR,52)-1+NST
               CALL WRITDD(A(LRM),JEDN,IPODS,LMAX13,LDUZI)
C
C              15.NALAZENJE  Fj(k+1)
C
               CALL CLEAR(A(LRTDT),JEDN)
               CALL MAXAPR(A(LSK),A(LRM),A(LRTDT),A(LMAXA),JEDN)
               IF(IX.LT.NSOPV) THEN
                  LMAX13=NPODS(JPBR,39)-1+NST+NSTAZ
                  CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
               ENDIF
               LRAD=LID
  820       CONTINUE
C
            LRAD=JMAXA
         ENDIF
  100  CONTINUE
C
C        NALAZENJE  MATRICE K
C
       IF (myid.eq.0) THEN
         LSKG=LRAD
         LRAD=LSKG+NSOPV*NSOPV*IDVA
         CALL CLEAR(A(LSKG),NSOPV*NSOPV)
         DO 900 JPBB=1,JPS
C
            JPBR=JPBB
            CALL PODDAT(NPODS,1)
            CALL PODDAT(NPODS,3)
            LSK=LRAD
            LRTDT=LSK+NWK*IDVA
            LRM=LRTDT+JEDN*NSOPV*IDVA
            LMAX=LRM+JEDN*IDVA
C           M - 
            CALL RSTAZK(NPODS,LSK,35)
            LMAX13=NPODS(JPBR,52)-1
            DO 910 I=1,NSOPV
               LRTD=LRTDT+JEDN*(I-1)*IDVA
               CALL READDD(A(LRTD),JEDN,IPODS,LMAX13,LDUZI)
  910       CONTINUE
C            WRITE(3,*) 'JPBR,JEDN,NSOPV',JPBR,JEDN,NSOPV
C            CALL WRR(A(LSK),NWK,'KPOD')
C            CALL WRR(A(LRTDT),JEDN*NSOPV,'RTDT')
            IJ=-1
            DO 920 I=1,NSOPV
               LRTD=LRTDT+JEDN*(I-1)*IDVA
               CALL CLEAR(A(LRM),JEDN)
               CALL MAXAPR(A(LSK),A(LRTD),A(LRM),A(LMAXA),JEDN)
            DO 920 J=1,NSOPV
               IJ=IJ+1
               LS=LSKG+IJ*IDVA
               LRT=LRTDT+JEDN*(J-1)*IDVA
               CALL SCALAR(A(LRM),A(LRT),CE,JEDN)
               CALL ZBIRM1(A(LS),CE,1)
  920       CONTINUE
C            CALL WRR(A(LSKG),NSOPV*NSOPV,'LSKG')
  900    CONTINUE
         LRAD=LID
         IF(INDBG.EQ.0) THEN
            IF(NDIN.GT.0) THEN
               IF(ISRPS.EQ.0)
     1WRITE(*,2020)
               IF(ISRPS.EQ.1)
     1WRITE(*,6020)
            ELSE
               IF(ISRPS.EQ.0)
     1WRITE(*,2030)
               IF(ISRPS.EQ.1)
     1WRITE(*,6030)
            ENDIF
         ENDIF
         NROOT=(NSOPV+1)/2
         NN=NSOPV
         NNM=NN+1
C        NWM=NWK
C        IF(NDIN.GT.0.AND.IMASS.EQ.2) NWM=NN
         RTOL=1.0D-06
         NC=MIN0(2*NROOT,NROOT+8)
         IF(NC.GT.NN) THEN
            NROT=NN/2
            NRO8=NROT
            IF(NN.GT.8) NRO8=NN-8
            NRO=MAX0(NROT,NRO8)
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NROOT,NRO
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NROOT,NRO
            STOP 'PROGRAM STOP - MIKA '
         ENDIF
         NITEM=30
         IFSS=IPROV
         IFPR=ISTSV
         IOUT = IZLAZ
         NNC = NC*(NC+1)/2
         NWK=NSOPV*(NSOPV+1)/2
         NWM=NSOPV
         LA=LRAD
         CALL DELJIV(LA,2,INDL)
         IF(INDL.EQ.0) LA=LA+1
         LB = LA + NWK*IDVA
         LMAXA= LB + NWM*IDVA
         LMAX=LMAXA+NSOPV+1
C      CALL WRR(A(LSKG),NSOPV*NSOPV,'SKG ')
         CALL PREBAC(A(LSKG),A(LA),A(LB),A(LMAXA),NSOPV)
         LR = LMAX
         LEIGV = LR + NN *NC*IDVA
         LTT = LEIGV + NC * IDVA
         LW = LTT + NN * IDVA
         LAR = LW + NN * IDVA
         LBR = LAR + NNC * IDVA
         LVEC = LBR + NNC * IDVA
         LD = LVEC + NC * IDVA
         LRTOLV = LD + NC * IDVA
         LBUP = LRTOLV + NC * IDVA
         LBLO = LBUP + NC * IDVA
         LBUPC = LBLO + NC * IDVA
         LMAX = LBUPC + NC * IDVA
         IF(LMAX.GT.MTOT) CALL ERROR(1)
C
         CALL SOPSTP(A(LA),A(LB),A(LMAXA),A(LR),A(LEIGV),A(LTT),A(LW)
     1               ,A(LAR),A(LBR),A(LVEC),A(LD),A(LRTOLV),A(LBUP)
     2               ,A(LBLO),A(LBUPC),NN,NWK,NWM,NROOT,RTOL,NC
     3               ,NITEM,IFSS,IFPR)
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' GRESKA U ULAZNIM PODACIMA'/
     1'  TRAZENI BROJ SOPSTVENIH VREDNOSTI JE NSOPV =',I5/
     2'  A MOZE SE MAKSIMALNO NACI',I5//'  PROGRAM STOP'//)
 2020 FORMAT(' *** RACUNANJE SOPSTVENIH VREDNOSTI ***')
 2030 FORMAT(' *** RACUNANJE KRITICNIH SILA ***')
C-----------------------------------------------------------------------
 6000 FORMAT(//' ERROR IN INPUT DATA'/
     1'  TRAZENI BROJ SOPSTVENIH VREDNOSTI JE NSOPV =',I5/
     2'  A MOZE SE MAKSIMALNO NACI',I5//'  PROGRAM STOP'//)
 6020 FORMAT(' *** CALCULATE EIGEN VALUE ***')
 6030 FORMAT(' *** CALCULATE CRITICAL LOADS ***')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE POCETI(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      ZA POCETNO POBUDJIVANJE
CS.   P R O G R A M
CS.      ZA POCETNO POBUDJIVANJE
C .
CE.       A(I)    - VECTOR
CS.       A(I)    - VECTOR
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' POCETI'
      A(N)=1.D0
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SCALAR(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO SCALAR MULTIPLY TWO VECTORS 
CS.   P R O G R A M
CS.      ZA SKALARNO MNOZENJE DVA VEKTORA 
C .
CE.       A(I)    - VECTOR
CE.       B(I)    - VECTOR
CE.       C       - CONSTANT
CS.       A(I)    - VECTOR
CS.       B(I)    - VECTOR
CS.       C       - KONSTANTA
C .       C = A * B
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SCALAR'
      C=0.D0
      DO 10 I=1,N
   10 C=C+A(I)*B(I)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE PREBAC(S,A,B,MAXA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      ZA PREPISIVANJE KVADRATNE MATRICE U MATRICU PO STUPCIMA 
CS.   P R O G R A M
CS.      ZA PREPISIVANJE KVADRATNE MATRICE U MATRICU PO STUPCIMA 
C .
CE.       S(I,J)  - MATRIX
CE.       A(I)    - VECTOR
CE.       B(I)    - VECTOR
CS.       S(I,J)  - MATRICA
CS.       A(I)    - VECTOR
CS.       B(I)    - VECTOR
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(N,*),A(*),B(*),MAXA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' PREBAC'
      MAXA(1)=1
      DO 10 I=1,N
         MAXA(I+1)=MAXA(I)+I
         B(I)=1.D0
         N1=MAXA(I)
         N2=MAXA(I+1)-1
         K=-1
      DO 10 J=N1,N2
         K=K+1
         IJ=I-K
         A(J)=S(IJ,I)
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SOPSTP(A9,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,
     1BUPC,NN,NWK,NWM,NROOT,RTOL,NC,NITEM,IFSS,IFPR)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      ZA POZIVANJE PROGRAMA ZA: UCITAVANJE MATRICA,
CE.      RESAVANJE SOPSTVENIH VREDNOSTI I STAMPANJE REZULTATA
CS.   PROGRAM
CS.      ZA POZIVANJE PROGRAMA ZA: UCITAVANJE MATRICA,
CS.      RESAVANJE SOPSTVENIH VREDNOSTI I STAMPANJE REZULTATA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /IMPERF/ NMODS,LIDIM,LSCIM,MODES
      COMMON /LANCZO/ LA,LB,LTT,LW
      DIMENSION A9(*),B(*),R(NN,*),TT(*),W(*),EIGV(*),
     1D(*),VEC(NC,*),AR(*),BR(*),RTOLV(*),BUP(*),BLO(*),BUPC(*),MAXA(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' SOPSTP'
C      CALL RSTAZK(A(LIPODS),LA,35)
C      CALL RSTAZK(A(LIPODS),LB,54)
C
CE    SLUCAJ ISTYP=2 , SAMO ZA NDIN=0
CS    SLUCAJ ISTYP=2 , SAMO ZA NDIN=0
C           KNL = (T-DT)K - (T)K
C
C      IF(ISTYP.EQ.2.AND.NDIN.EQ.0) CALL TDTKTK(A9,B,NWK)
C
CE    KORIGOVANJE MATRICA MAKSIMALNIM BROJEM NA DIAGONALI
CS    KORIGOVANJE MATRICA MAKSIMALNIM BROJEM NA DIAGONALI
C      CALL MAKSIM(A9,MAXA,NWK,NN,AMAXK)
C      CALL MAKSIM(B,MAXA,NWM,NN,AMAXF)
C      IF(AMAXK.GT.AMAXF) THEN
C         CALL MAXDEL(A9,AMAXK,NWK)
C         CALL MAXDEL(B,AMAXK,NWM)
C      ELSE
C         CALL MAXDEL(A9,AMAXF,NWK)
C         CALL MAXDEL(B,AMAXF,NWK)
C      ENDIF
C
C      CALL IWRR(MAXA,NC+1,'MA  ')
C      CALL WRR(A9,NWK,'A9  ')
C      CALL WRR(B,NWM,'B   ')
      CALL SSPACE (A9,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,
     1BUPC,NN,NWK,NWM,NROOT,RTOL,NC,NITEM,IFSS,IFPR)
C
CE    PRINT RESULTS
CS    STAMPANJE REZULTATA
C
C      IF(NDIN.GT.0) THEN
         PI=4.0D0*ATAN(1.0D0)
         DO 10 I=1,NC
   10    EIGV(I)=DSQRT(EIGV(I))/PI/2.0D0
         CALL WRR(EIGV,NC,'EIGV')
C         DO 580 J=1,NROOT
C            CALL MAKSIS(R(1,J),A(LID),SIL)
C            CALL STAPSV(R(1,J),A(LID),A(LCVEL),ICVEL,NP,J,EIGV(J),SIL,0)
C            IF(NBLGR.GE.0)
C     1    CALL STAGSV(R(1,J),A(LID),A(LCVEL),ICVEL,NP,J,EIGV(J),IGRAF,0)
C  580    CONTINUE
C      ELSE
C
CE       SLUCAJ ISTYP=1 , RADI PREKO GAMA
CS       SLUCAJ ISTYP=1 , RADI PREKO GAMA
C
C         IF(ISTYP.EQ.1)THEN
C            DO 20 I=1,NC
C   20       EIGV(I)=1.0D0/(1.0D0-EIGV(I))
C         ENDIF
CC
C         CALL RSTAZ(A(LIPODS),LTT,53)
C         CALL RSTAZ(A(LIPODS),LW,39)
C
C        OPEN FILE  'ZMODES'
C         CALL MODEOP(MODES)
C         REWIND MODES
CC
C         DO 30 J=1,NROOT
C            B1=EIGV(J)
C            A1=1.0D0-B1
C            CALL ZBIR2(B,TT,W,A1,B1,NN)
C            CALL STAPSV(B,A(LID),A(LCVEL),ICVEL,NP,J,EIGV(J),SIL,1)
C            CALL SUMSIL(B,A(LID),IOPGL,NP,IZLAZ,1)
C            CALL MAKSIS(R(1,J),A(LID),SIL)
C            CALL STAPSV(R(1,J),A(LID),A(LCVEL),ICVEL,NP,J,EIGV(J),SIL,2)
C            IF(NBLGR.GE.0)
C     1    CALL STAGSV(R(1,J),A(LID),A(LCVEL),ICVEL,NP,J,EIGV(J),IGRAF,1)
C
CE          SAVE FOR INITIAL IMPERFECTION
CS          PODACI ZA INICIJALNU IMPERFEKCIJU
C
C            WRITE(MODES) EIGV(J),(R(I,J),I=1,NN)
C
C   30    CONTINUE
C         CLOSE (MODES,STATUS='KEEP')
C      ENDIF
      RETURN
      END
