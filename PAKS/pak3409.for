C=======================================================================
C
CS   PLASTICNOST 3/D ELEMENT
CE  PLASTICITY CAM-CLAY MATERIAL  MODEL 3D ELEMENT
C
C    SUBROUTINE D3M9
C               TAUI35
C               TEQBI3
C               PRILA3
C               DEVEQ3
C
      SUBROUTINE D3M9(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC,lpoc0)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
CS     NA NIVOU INTEGRACIONE TACKE
CE     PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
      include 'paka.inc'
      
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /OBNOVA/ IRUSI
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M9'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LPOT=LDEFPP + 6
      LEMP=LPOT + 1
      LOCR=LEMP + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LPOT1=LDEFP1 + 6
      LEMP1=LPOT1 + 1
      LOCR1=LEMP1 + 1
C
      IF(IRUSI.EQ.1) THEN
      LTAU0=LPOC0
      LDEFT0=LTAU0 + 6
      LDEFP0=LDEFT0 + 6
      LPOT0=LDEFP0 + 6
      LEMP0=LPOT0 + 1
      LOCR0=LEMP0 + 1
      ELSE
      LTAU0=1
      LDEFT0=LTAU0 + 6
      ALLOCATE (PLAS0(12))
      ENDIF

C
      CALL TAUI39(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),PLAST(LPOT),
     1            PLAST(LEMP),PLAST(LOCR),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),PLAS1(LPOT1),
     1            PLAS1(LEMP1),PLAS1(LOCR1),
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC,
     1            PLAS0(LTAU0),PLAS0(LDEFT0))
!      CALL TAUI39(A(LTAU),A(LDEFT),A(LDEFPP),A(LPOT),A(LEMP),A(LOCR),
!     1       A(LTAU1),A(LDEFT1),A(LDEFP1),A(LPOT1),A(LEMP1),A(LOCR1),
!     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC,
!     1       A(LTAU0),A(LDEFT0))
C
      IF(IRUSI.NE.1) DEALLOCATE(PLAS0)
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TAUI39(TAUT,DEFT,DEFPP,P0T,EMP,OCR,
     1                  TAU1,DEF1,DEFP1,P0TDT,EMP1,OCR1,
     1                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC,
     1                  tau0,def0)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS     CAM-CLAY MATERIJALNI MMODEL
CE     PROGRAM FOR STRESS INTEGRATION FOR CAM-CLAY CAP MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1               DETAU(6),DDEF(6)
      COMMON /MAT2D/ EE,ANI,ET,TEQY0
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
C
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
C
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /RESTAR/ TSTART,IREST
      COMMON /restap/ irestp,lplas0,lstaz0
      COMMON /OBNOVA/ IRUSI
      COMMON /IZADJI/ NASILU
      common /kadamb/ kadmb(100)
      common /crklie/ icrkli(100000)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(*),DEFT(*),DEFPP(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP1(*),tau0(6),def0(6),TAUD0(6),TAUDP(6)
      DIMENSION FUN(9,MATE),NTA(*)
      IF(IDEBUG.EQ.1) PRINT *, 'TAUI39'
C
      IF(DABS(DEF(1)).GT.10..OR.DABS(DEF(2)).GT.10.
     1                      .OR.DABS(DEF(3)).GT.10.) THEN
         WRITE(*,*) 'NLM,MAT,DEF',NLM,MAT,(DEF(I),I=1,3)
         WRITE(3,*) 'NLM,MAT,DEF',NLM,MAT,(DEF(I),I=1,3)
         STOP 'TAUI39'
      ENDIF
C
CS    OSNOVNE KONSTANTE
CE    BASIC KONSTANTS
C
      IBISE=0
      ISOFT=0
      IPL1=0
      ICSS=0
C
      DVT=2.0D0/3.0D0
      SQ2=DSQRT(2.0D0)
      DJP=DSQRT(1.5D0)
C
CS    MATERIJALNE KONSTANTE
CE    MATERIAL CONSTANTS
C
      EE    = FUN(1,MAT)
      ANI  = FUN(2,MAT)
      AT   = FUN(3,MAT)
      IEL  = INT(FUN(4,MAT))
C
      AEM  = FUN(5,MAT)
      ALAM = FUN(6,MAT)
      AKA  = FUN(7,MAT)
      P0   = -FUN(8,MAT)
      AE0  = FUN(9,MAT)
C
      AMM=AEM

C     TOLERANCIJE
c      TOLFY=1.0D-8
C      TOLFY=EE*TOLFY
      TOLFY=1.0D-5
      TOLBIS=1.0D-5
c      TOLP0=1.0D-8
C      TOLP0=EE*TOLP0
      TOLP0=1.0D-3
C
      moje=0
c za moje=1 radi dobro pre restarta ali ne valja posle restarta f=con
c      IPOC=1
      ipoc=irusi
      SM0=0.0D0
      IF(IPOC.EQ.1) SM0=(TAU0(1)+TAU0(2)+TAU0(3))/3.0D0
      IF(KOR.EQ.1) THEN
         P0T=P0
c proveri kako radi sa nulom
         SMT=P0
      ELSE
         if(moje.eq.0) SMT=(TAUT(1)+TAUT(2)+TAUT(3))/3.0D0
         if(moje.eq.1) SMT=(TAU1(1)+TAU1(2)+TAU1(3))/3.0D0
      ENDIF
C
      EM0=0.D0
c posle restarta = pocetno EE         
      IF(IPOC.EQ.1) EM0=(DEF0(1)+DEF0(2)+DEF0(3))/3.0D0
      if(moje.eq.0) EMT=(DEFT(1)+DEFT(2)+DEFT(3))/3.0D0
      if(moje.eq.1) EMT=(DEF(1)+DEF(2)+DEF(3))/3.0D0
c posle restarta = pocetno Et         
      ET=(1.0D0+AE0)*DEXP(3.*(EMT+EM0))-1.0D0
      ALK = ALAM-AKA 
      EAL=-(1.+ET)/ALK
C
      IF(IEL.EQ.1) THEN
c         SMTA=DABS(SMT)
         SMTA=-SMT+AT
c         SMTA=-P0T+AT
         AMS=(1.0D0+ET)*SMTA/AKA
c posle restarta = pocetno EE         
         EE=3.0D0*AMS*(1.0D0-2.0D0*ANI)
         if(EE.lt.1.e-8) then
            write(*,*) ' NLM,MAT,EE',NLM,MAT,EE
            write(3,*) ' NLM,MAT,EE',NLM,MAT,EE
            stop ' EE<0'
         endif
      ENDIF
C
      IF(KOR.EQ.1) THEN   
         SMT=SM0
      ENDIF
      P0TDT=P0T
C
CZILE
C OVO SE KORISTI SAMO ZA BETON KOD TUNELA ! NECE=0
      nece=1
c menja beton u prvom koraku posle restarta (za bogovinu)
C      matfix=1
C      if(mat.eq.matfix.and.irestp.eq.1.and.nece.eq.0) then
c menja beton u drugom koraku posle restarta (za tunel)
      matfix=2
      if((mat.eq.matfix.and.irestp.eq.1.and.irest.eq.0.and.nece.eq.0)
     1   .or.(kadmb(mat).gt.0.and.kadmb(mat).le.kor))then
         icss=1
         do i=1,6
c	      def(i)=def(i)-deft(i)
            defpp(i)=0.
            defp1(i)=0.
            if(dabs(defpp(i)).gt.1.e-9) then
               write(3,*) ' pp',i,defpp(i),nlm,mat,irest
               stop
            endif
            if(dabs(defp1(i)).gt.1.e-9) then
               write(3,*) ' p1',i,defp1(i),nlm,mat,irest
               stop
            endif
            if(dabs(tau0(i)).gt.1.e-9) then
               write(3,*) ' t0',i,tau0(i),nlm,mat,irest
               stop
            endif
         enddo   
         emp=0.
         emp1=0.
         ocr=0.
         ocr1=0.
         if(dabs(emp).gt.1.e-9.or.dabs(emp1).gt.1.e-9) then
            write(3,*) ' em',emp,emp1,nlm,mat,irest
            stop
         endif
         if(dabs(ocr).gt.1.e-9.or.dabs(ocr1).gt.1.e-9) then
            write(3,*) ' em',ocr,ocr1,nlm,mat,irest
            stop
         endif
c        bogovina	      
c         EE=15.0e6
c         ANI=0.18
c        tunel
!         EE=30.0e6
!         ANI=0.2
      endif
      CALL MEL3EL(ELAST,EE,ANI)
      AE=(1.0D0+ANI)/EE
      AM=(1.0D0-2.0D0*ANI)/EE
      CM=1.0D0/AM
C
C
      IF(IRAC.EQ.2) return
      DO 10 I=1,6
         TAU(I)=0.0D0
         DDEF(I)=DEF(I)-DEFPP(I)
   10 CONTINUE
C
CE    SOLUTION IS ELASTOPLASTIC , DETERMINE PLASTIC STRAIN AND STRESS
CE    INCREMENT BY USING EFFECTIVE STRESS FUNCTION
CS    RESENJE JE ELASTOPLASTICNO , TREBA ODREDITI PRIRASTAJ PLASTICNE
CS    DEFORMACIJE IPRIRASTAJ NAPONA PRIMENOM FUNKCIJE EFEKTIVNOG NAPONA
C
CS    DEVIJATOR UKUPNE DEFORMACIJE
CE    TOTAL STRAIN DEVIATOR 
C
      EMT=(DEF(1)+DEF(2)+DEF(3))/3.0D0
      DO 180 I=1,6
         IF(I.LE.3) THEN
            DEFDPR(I)=DEF(I)-EMT
         ELSE
            DEFDPR(I)=0.5D0*DEF(I)
         ENDIF 
  180 CONTINUE
C
CS    DEVIJATOR PROBNE ELASTICNE DEFORMACIJE
CE    ELASTIC STRAIN DEVIATOR 
C
      DO 30 I=1,6
         IF(I.LE.3) THEN
            DEFDS(I)=DEFDPR(I)-DEFPP(I)+EMP
         ELSE
            DEFDS(I)=DEFDPR(I)-DEFPP(I)
         ENDIF
   30 CONTINUE
C
CS    ODREDIVANJE POJEDINIH CLANOVA FUNKCIJE EFEKTIVNOG NAPONA
CE    EFFECTIVE STRESS FUNCTION DEFINITION
C
C
CS    SKALARNI PROIZVOD DEVIJATORA ELASTICNE DEFORMACIJE I NORMA
C
      DKV=(DEFDS(1)*DEFDS(1)+DEFDS(2)*DEFDS(2)+DEFDS(3)*DEFDS(3)
     1+2.0D0*(DEFDS(4)*DEFDS(4)+DEFDS(5)*DEFDS(5)+DEFDS(6)*DEFDS(6)))
      DD=DSQRT(DKV)
C
CS    PROVERA ELASTICNOG RESENJA
CS    ODREDIVANJE NAPONA KOJI ODGOVARA ELASTICNOM RESENJU
CE    PLASTIC YIELDING CHECK
CE    ELASTIC SOLUTION
C
      DO 32 I=1,6
         if(i.le.3) then
            taud0(i)=tau0(i)-sm0
         else
            taud0(i)=tau0(i)
         endif
         IF(IPOC.EQ.1) THEN
            TAUD(I)=DEFDS(I)/AE+taud0(i)
         ELSE
            TAUD(I)=DEFDS(I)/AE
         ENDIF
   32 continue
C
C     KOREN IZ J2D
C
      AJ2DE=0.5*(TAUD(1)*TAUD(1)+TAUD(2)*TAUD(2)+TAUD(3)*TAUD(3))+
     1           TAUD(4)*TAUD(4)+TAUD(5)*TAUD(5)+TAUD(6)*TAUD(6)
      AJ2DQ=DSQRT(3.0D0*AJ2DE)
C
C     SREDNJA ELASTICNA DEFORMACIJA=SREDNJA UKUPNA-SREDNJA PLASTICNA
C
      EMS=EMT-EMP
C
C     SREDNJI NAPON
C
      SMTE=EMS/AM+sm0
C     SMTE=EMS/AM
      SMTDT=SMTE
C
CZILE
C     PRESKACE PLASTICNOST UVEK
c     if(mat.eq.matfix.and.nece.eq.0) then
C     BOGOVINA PRESKACE U PRVOM KORAKU POSLE RESTARTA
c      if(mat.eq.matfix.and.irestp.eq.1.and.nece.eq.0) then
C     TUNEL PRESKACE U DRUGOM KORAKU POSLE RESTARTA
c      if(mat.eq.matfix.and.irestp.eq.1.and.irest.eq.0.and.nece.eq.0)then
      if((mat.eq.matfix.and.irestp.eq.1.and.irest.eq.0.and.nece.eq.0)
     1   .or.(kadmb(mat).gt.0.and.kadmb(mat).le.kor))then
         go to 100
      endif
C
C     NASILNI IZLAZAK IZ MODELA, RACUNA SAMO ELASTICNO ODSTOJANJE OD LOMA
      IF(NASILU.EQ.1) GO TO 100
CZILE
C      
C OVO PROVERITI NA KRAJU
C     RAZE=2.0D0*SMTE-P0T
      RAZE=2.0D0*SMTE-P0T-AT
      RAZ=RAZE 
      DELEMP=0.0D0
      DLAM=0.0D0
      AER=AE
C
CS    KOREKCIJA FY ZBOG POMERENOSTI POVRSI TECENJA ZA AT
CE    CORRECTION OF FY ACCORDING TO YIELD SURFACE TRANSLATION AT
      RAZE=2.0D0*SMTE-P0T-AT
      PPT=SMTE-(P0T+AT)/2.0D0
      PMT=(P0T-AT)/2.0D0
      FYE=PPT*PPT- PMT*PMT+ 3.* AJ2DE/(AMM*AMM)
C
      IF(DABS(SMTE).GT.TOLFY) THEN
         PC=SMTE+3.0D0*AJ2DE/(AMM*AMM*SMTE)
         IF(DABS(PC).GT.TOLFY) OCR1=P0T/PC
      ELSE
         OCR1=1.0D3
      ENDIF
C
      TANAE=0.
      IF(DABS(SMTE).GT.1.E-3) TANAE=-AJ2DQ/(SMTE-AT)
C
C     ELASTICNO RESENJE
C
      ist=0
      if(jg.eq.1.and.ist.eq.1) then
      if(jg.eq.1) write(3,*) ' NOVI ULAZ U MODEL'
      write(3,*) ' JG,iter,nlm',JG,iter,nlm
      write(3,*) ' p0t,emp,ocr',p0t,emp,ocr
      write(3,*) ' p0tdt,emp1,ocr1',p0tdt,emp1,ocr1
      call wrr6(tau0,6,'tau0')
      call wrr6(taut,6,'taut')
      call wrr6(tau1,6,'tau1')
      call wrr6(def0,6,'def0')
      call wrr6(deft,6,'deft')
      call wrr6(def1,6,'def1')
      call wrr6(def,6,'def ')
      call wrr6(defpp,6,'dept')
      call wrr6(defp1,6,'dep1')
	write(3,*) ' fye,tolfy,raze,tanae',fye,tolfy,raze,tanae
	WRITE(3,*) 'ULAZ-SMTDT,P0TDT',SMTDT,P0TDT
	CALL WRR6(tauD,6,'TAD0')
      IF(FYE.LE.TOLFY) then
	   WRITE(3,*)' ELASTICNO'
	ELSE
	   WRITE(3,*)' PLASTICNO'
	ENDIF
      endif
C
C
      IF(FYE.LE.TOLFY) GO TO 100
C
C     ELASTOPLASTICNO RESENJE
C
      IPL1=1
      ICSS=1
C
CS    KRITICNO STANJE IZNAD TEMENA A ELIPSE 
CE    CRITICAL STATE
C
      TOLPT=1.0D0-2.0D0*SMTE/(P0T+AT)      
C KOJIC
C      TOLPT=DABS(2.0D0*SMTE-P0T-AT)/(P0T+AT)      
      IF(DABS(TOLPT).LT.TOLP0) THEN
         DAST=DABS(SMTE-AT)*AMM
c naredni uslov treba najverovatnije izbaciti!!!!!!
cn         IF(AJ2DQ.GT.DAST) THEN
cn            write(3,*) ' A iznad'
            OCR1=1.0D0
            SMTDT=SMTE
            P0TDT=2.0D0*SMTDT-AT
            DLAM= (DJP*DD/AMM/dabs(SMTDT-AT)-AE)*AMM*AMM/3.0D0
            if(dlam.lt.-1.d-10) then
               write(3,*) ' lamda1<0',dlam
               write(3,*) ' JG,iter,nlm',JG,iter,nlm
c               DLAM= (DJP*DD/AMM/(SMTDT-AT)-AE)*AMM*AMM/3.0D0
c               if(dlam.lt.-1.d-10) then
c                  write(3,*) ' lamda2<0',dlam
c                  stop
c               endif
            endif
            AER=AE+3.0D0*DLAM/AMM/AMM
            DELEMP=0.0D0     
cn	   ELSE
cn            write(3,*) ' A ispod'
cn            stop ' A ispod'
cn         ENDIF 
         GO TO 330
      ENDIF
C 
CS    PODRUCJE OJACANJA
CE    HARDENING REGION
C
C
CS       ODREDIVANJE NULE FUNKCUJE EFEKTIVNOG NAPONA
CE       SOLUTION OF NONLINEAR EQUATIN FN=0.0 - NEWTON METHOD
C
         IBISE=0
         IBISM=100
         ICSS=0
         TOLBIS=1.0D-5
         DELEMP=0.0D0
         DELE1=0.0D0
C        POROZNOST
         ETDT=(1.0D0+AE0)*DEXP(3.*(EMT+EM0))-1.0D0
         ALK = ALAM-AKA 
         EAL=-(1.0D0+ETDT)/ALK
         BVTDT=ALK/(1.0D0+ETDT)/3.0D0
         OCR1=1.0D0
C
   50    IBISE=IBISE+1
         DELE1=DELEMP
C
C        KOREKCIJA SREDNJEG NAPONA
         SMTDT=SMTE-DELEMP/AM
C        KOREKCIJA SREDNJE ELASTICNE DEFORMACIJE
         EME=EMS-DELEMP
C        NOVA VELICINA ELIPSE
cb         P0TDT=P0T*DEXP(EAL*3.0D0*DELEMP)
         P0TDT=P0T*DEXP(-DELEMP/BVTDT)
         RAZ=2.0D0*SMTDT-P0TDT-AT
         TOLPT=1.0D0-2.0D0*SMTDT/(P0TDT+AT)      
         IF(DABS(TOLPT).LT.TOLP0) THEN
cs           write(*,*)NLM,MAT,delemp,rAZ,TOLPT,'NLM,MAT,delemp,rAZ,TOLPT'
cs           write(3,*)NLM,MAT,delemp,RAZ,TOLPT,'NLM,MAT,delemp,rAZ,TOLPT'
C            STOP
C            SMTDT=SMTE
C            P0TDT=2.0D0*SMTDT-AT
c            DLAM=0.
            IF(dabs(raz).gt.tolp0) then
               DLAM=3.0D0*DELEMP/RAZ
               if(dlam.lt.-1.d-10) then
                  write(3,*) ' lamda3<0',dlam
                  write(*,*) ' lamda3<0',dlam
cccccc	         dlam=0
c                  stop
               endif
            else
               DLAM=3.0D0*DELEMP/RAZ
               write(*,*) NLM,MAT,Dlam,'NLM,MAT,Dlam,OPASNO!!!'
               write(3,*) NLM,MAT,Dlam,'NLM,MAT,Dlam,OPASNO!!!'
c               stop 'raz=0'
            endif
C            DELEMP=0.0D0     
            AER=AE+3.0D0*DLAM/AMM/AMM
c            GO TO 330
c                  OCR1=1.0D0
c                  SMTDT=SMTE
c                  P0TDT=2.0D0*SMTDT-AT
c                  DLAM= (DJP*DD/AMM/(SMTDT-AT)-AE)*AMM*AMM/3.0D0
c                  AER=AE+3.0D0*DLAM/AMM/AMM
c                  DELEMP=0.0D0     
c            GO TO 330
         ELSE
            DLAM=3.0D0*DELEMP/RAZ
            AER=AE+3.0D0*DLAM/AMM/AMM
         ENDIF
C
         AER=AE+3.0D0*DLAM/AMM/AMM
C	
CS       KOREKCIJA FY ZBOG POMERENOSTI POVRSI TECENJA ZA AT
CE       CORRECTION OF FY ACCORDING TO YIELD SURFACE TRANSLATION AT
         PPT=SMTDT-(P0TDT+AT)/2.0D0
         PMT=(P0TDT-AT)/2.0D0
         TFQ=PPT*PPT- PMT*PMT+ 3.0*DKV/(AER*AER)/(2.*AMM*AMM)
C
cb         TFQP=-RAZ/AM-0.5D0*(RAZ+P0TDT-AT)*3.0D0*EAL*P0TDT-
         TFQP=-RAZ/AM+0.5D0*(RAZ+P0TDT-AT)*P0TDT/BVTDT-
     1         3.0D0*DKV/AMM**4/AER**3*9.*
cb     2        (1.0D0/RAZ+DELEMP*(2.0/AM+3.0D0*EAL*P0TDT)/(RAZ*RAZ))
     2        (1.0D0/RAZ+DELEMP*(2.0/AM-P0TDT/BVTDT)/(RAZ*RAZ))
C 
         DELEMP=DELEMP-TFQ/TFQP
      if(delemp.gt.1.) then
         write(3,*) delemp,TDEL,TFQ,NLM,MAT,'delemp,TDEL,TFQ,NLM,MAT'
         write(*,*) delemp,TDEL,TFQ,NLM,MAT,'delemp,TDEL,TFQ,NLM,MAT'
      ENDIF
         TDEL=DABS((DELEMP-DELE1)/DELEMP)
C
C
C        PROVERA KONVERGENCIJE RESENJA ZA FUNKCIJU TECENJA
C
c         IF(TDEL.LT.TOLBIS.OR.DABS(TFQ).LT.TOLBIS) THEN
C         IF(TDEL.LT.TOLBIS.OR.DABS(TFQ).LT.TOLFY) THEN
         IF(TDEL.LT.TOLBIS.and.DABS(TFQ).LT.TOLFY) THEN
C*          IZLAZI POSTO JE KONVERGIRAO
            GO TO 330
C*
            DO 61 I=1,6
               TAUDP(I)=DEFDS(I)/AER
   61       continue
            AJ2DP=0.5*
     1           (TAUDP(1)*TAUDP(1)+TAUDP(2)*TAUDP(2)+TAUDP(3)*TAUDP(3))
     1           +TAUDP(4)*TAUDP(4)+TAUDP(5)*TAUDP(5)+TAUDP(6)*TAUDP(6)
            AJ2DP=DSQRT(3.0D0*AJ2DP)
            TANAT=-AJ2DP/SMTDT
C*          IZLAZI POSTO JE KONVERGIRAO
            GO TO 330
c*
	      IF(TANAT.GT.TANAE) THEN
	         write(3,*) ' desno iznad'
               DAST=DABS(SMTDT-AT)*AMM
               IF(AJ2DP.GT.DAST) THEN
                  OCR1=1.0D0
C                 SMTDT=SMTE
                  P0TDT=2.0D0*SMTDT-AT
                  DLAM= (DJP*DD/AMM/(SMTDT-AT)-AE)*AMM*AMM/3.0D0
                  AER=AE+3.0D0*DLAM/AMM/AMM
                  DELEMP=0.0D0     
               ENDIF 
            else
	         write(3,*) ' desno ispod'
	      ENDIF
            GO TO 330
         ENDIF
C
         IF(IBISE.GT.IBISM) THEN
            IF(TDEL.LT.TOLBIS) GO TO 330
            IF(DABS(TFQ).LT.TOLFY) GO TO 330
            WRITE(*,1000) NLM,MAT,tdel,tolbis,tfq,tolfy
            WRITE(3,1000) NLM,MAT,tdel,tolbis,tfq,tolfy
 1000       FORMAT(' ','DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUINT5'/
     1                 'NLM,MAT,tdel,tolbis,tfq,tolfy'/
     1             2I5,4(1PD12.5))
            icrkli(nlm)=1
            if(icrkli(nlm).eq.1) then
               call clear(tau,6)
               call clear(tau1,6)
               return
            endif
c            STOP
         ENDIF
C
         GO TO 50
C
C      ELSE
C
CS       TACKA IZNAD KRITICNE LINIJE - OMEKSANJE
CE       POINT ABOVE CRITICAL LINE - SOFTENING
C         DAST=DABS(SMTE-AT)*AMM
C         IF(AJ2DE.GT.DAST) THEN
c      	  write(3,*) ' levo iznad'
C            OCR1=1.0D0
C            SMTDT=SMTE
C            P0TDT=2.0D0*SMTDT-AT
C            DLAM= (DJP*DD/AMM/(SMTDT-AT)-AE)*AMM*AMM/3.0D0
C           DLAM= (DD/(sq2*AMM*SMTDT)-AE)*AMM*AMM/3.0D0
C            AER=AE+3.0D0*DLAM/AMM/AMM
C            DELEMP=0.0D0     
C	   ELSE
C            write(3,*) ' levo ispod'
C        ENDIF 
C      ENDIF
  
C
C
CS    DEVIJATOR NAPONA
CE    STRESS DEVIATOR COMPONENTS
C
  330 DO 65 I=1,6
         TAUD(I)=DEFDS(I)/AER
   65 continue
C
CS    SREDNJA PLASTICNA DEFORMACIJA
C
      EMP1=EMP+DELEMP
C
CS    UKUPNA PLASTICNA DEFORMACIJA = STARA PLASTICNA DEFORMACIJA +
CS    PRIRASTAJ DEVIJATORA + PRIRASTAJ SREDNJE PLASTICNE DEFORMACIJE
CE    TOTAL PLASTIC STRAIN = OLD PLASTIC STRAIN +
CE    DEVIATORIC INCREMENT + MEAN INCREMENT PLASTIC STRAIN 
C
c      IF(IRESTP.EQ.1.AND.IPOC.EQ.1) THEN
c         DO 75 I=1,6
c            TAUD(I)=TAUD(I)+taud0(i)
c   75    continue
c      ENDIF
       DO 70 I=1,6
         IF(I.LE.3) THEN
            DEFP1(I)=DEFPP(I)+3.0D0*DLAM/(AMM*AMM)*TAUD(I)+DELEMP
         ELSE
            DEFP1(I)=DEFPP(I)+3.0D0*DLAM/(AMM*AMM)*TAUD(I)
         ENDIF
   70 CONTINUE
C
C     DODAVANJE POCENTOG DEVIJATORA NAPONA
C
      IF(IRESTP.EQ.1.AND.IPOC.EQ.1) THEN
         DO 75 I=1,6
            TAUD(I)=TAUD(I)+taud0(i)
  75     continue
      ENDIF
      if(DABS(delemp).gt.1..OR.DABS(DLAM).GT.10000.) then
         write(*,*) delemp,Dlam,NLM,MAT, 'delemp,Dlam,NLM,MAT'
         write(3,*) delemp,Dlam, NLM,MAT, 'delemp,Dlam,NLM,MAT'
         stop 'lamda>>'
      endif
C
C
CS     NAPON
CE     STRESS
C
  100 DO 260 I=1,6
         IF(I.LE.3) THEN
            TAU(I)=TAUD(I)+SMTDT
         ELSE
            TAU(I)=TAUD(I)
         ENDIF
  260 CONTINUE
C
CS    ELASTOPLASTICNA MATRICA
CE    ELASTIC-PLASTIC CONSTITUTIVE MATRIX
C
C     PROVERI USLOV DA NE ULAZI ZA ELASTICNO RESENJE
      IF(ITER.EQ.0.OR.ICSS.EQ.1) GO TO 110
      IF((METOD.EQ.4).OR.(METOD.EQ.3).OR.(METOD.EQ.63).OR.
     *(METOD.EQ.73)) THEN
C
       CALL MEL39(TAU,DEFDS,SMTDT,P0TDT,DLAM,EAL,AMM,AER,AM,
     1           RAZ,DKV,DELEMP,ICSS)
C
      ENDIF
C
CS    KORIGOVANJE VELICINA IZ PRETHODNOG KORAKA KAD JE POSTIGNUTA
CS    KONVERGENCIJA
CE    CORECTION OF VALUES FROM PREVIOUS STEP WHEN CONVERGENCE IS
CE    REATCHED
C
  110 CONTINUE
      DO 290 I=1,6
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
      if(jg.eq.1.and.ist.eq.1) then
	write(3,*) ' IZLAZ'
	write(3,*) ' izlaz-tfq,tolfy,raz',tfq,tolfy,raz
	write(3,*) ' izlaz-Ip0t,emp,ocr',p0t,emp,ocr
	write(3,*) ' izlaz-Ip0tdt,emp1,ocr1',p0tdt,emp1,ocr1
      call wrr6(tau0,6,'taI0')
      call wrr6(taut,6,'taIt')
      call wrr6(tau1,6,'taI1')
      call wrr6(def0,6,'def0')
      call wrr6(deft,6,'deft')
      call wrr6(def1,6,'def1')
      call wrr6(def,6,'def ')
      call wrr6(defpp,6,'dept')
      call wrr6(defp1,6,'dep1')
C     provera funkcije tecenja pre izlaska
C
C     KOREN IZ J2D
C
      AJ2DE=0.5*(TAUD(1)*TAUD(1)+TAUD(2)*TAUD(2)+TAUD(3)*TAUD(3))+
     1           TAUD(4)*TAUD(4)+TAUD(5)*TAUD(5)+TAUD(6)*TAUD(6)
      AJ2DQ=DSQRT(3.0D0*AJ2DE)
      PPT=SMTDT-(P0TDT+AT)/2.0D0
      PMT=(P0TDT-AT)/2.0D0
      FYE=PPT*PPT- PMT*PMT+ 3.* AJ2DE/(AMM*AMM)
	WRITE(3,*) ' izlaz-SMTDT,P0TDT',SMTDT,P0TDT
	CALL WRR6(tauD,6,'TAUD')
	write(3,*) ' izlaz-fye,tolfy',fye,tolfy
	IF(FYE.GT.TOLFY) WRITE(3,*)' NE VALJA',JG,NLM
 	IF(FYE.GT.TOLFY) WRITE(*,*)' NE VALJA',JG,NLM
      endif
      RETURN
      END
C======================================================================
C
      SUBROUTINE MEL39(TAU,DEFDS,SMTDT,P0TDT,DLAM,EAL,AMM,AER,AM,
     1           RAZ,DKV,DELEMP,ICSS)
C
C----------------------------------------------------------------------
CS    PROGRAM ZA FORMIRANJE MATRICE C<E> ILI C<EP> ZA CAM-CLAY MODEL
CE    PROGRAM FOR CONSTITUTIVE MATRIX C<E> OR C<EP>
C----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION C(6,6),EPSD(6),TAU(6),DEFDS(6)
      IF(IDEBUG.EQ.1) PRINT *, 'MEL09 '
C
      DVT=2.0/3.0
      TR=1.0/3.0
      CM=1.0/AM
C
      DO 1 I=1,3
      EPSD(I)=DEFDS(I)
    1 EPSD(I+3)=2.0*DEFDS(I+3)
C
      D2=EPSD(1)**2+EPSD(2)**2+EPSD(3)**2+
     *(EPSD(4)**2+EPSD(5)**2+EPSD(6)**2)/2.
C
       AA=1.0/AER
       AAM=AA/AMM**2
       A2M2=1.0/(AER*AER*AMM*AMM)
       AL2=3.0*P0TDT*EAL  
       SA2=SMTDT*AL2
C
        IF(ICSS.EQ.0) THEN
C
       AL3=(3.+DLAM*(2.0*CM+AL2))/RAZ
       A43=AMM**4*AER**3
       D9=9.0*DKV/A43
C  AL4
       AL4=CM*RAZ+SA2+D9*AL3
C  AL15
       AL5=3.0*AL3*A2M2/AL4
       AL5=3.0*AL5*A2M2
        ELSE
         AL5=DSQRT(DVT*DKV)*AMM/SMTDT
        ENDIF 
C
       RAZ2=RAZ*RAZ
         IF(ICSS.EQ.0) THEN
       AL14=-2.0*DLAM*CM/RAZ
       AL15= 3.0/RAZ+3.0*DELEMP/RAZ2*(2.*CM+AL2)
       AL16=RAZ*CM-D9*AL14
       AL17=RAZ*CM+SA2+D9*AL15
       A67 = AL16/AL17
       AL18=AL14+AL15*A67
       AL19=CM*(1.0-A67)
C
        ELSE
        AL18=AL5
        AL19=CM 
        ENDIF
C
      DO 4 I=1,6
      DO 4 J=1,6
      ELAST(I,J)=0.0D0
    4 C(I,J)=0.0D0
C
CS   MATRICA __Cik' NADVUCENO
C
      DO 5 I=1,6
      DO 5 J=1,6
      C(I,J)=-AL5*DEFDS(I)*EPSD(J)
      IF(I.EQ.J)C(I,J)=C(I,J)+AA
    5 CONTINUE
C
CS    MATRICA C'ij
      ACMP=-A2M2*AL18
      DO 6 I=1,3
      DO 6 J=1,3
      DO 16 K=1,3
       IF(J.EQ.K) THEN
       TT=DVT
       ELSE
       TT=-TR
       ENDIF
      ELAST(I,J)=ELAST(I,J)+C(I,K)*TT
   16 ELAST(I+3,J)=ELAST(I+3,J)+C(I+3,K)*TT
    6 CONTINUE
C
      DO 40 I=1,3
      DO 40 J=1,3
      ELAST(I,J)=ELAST(I,J)+ACMP*DEFDS(I)/3.
   40 ELAST(I+3,J)=ELAST(I+3,J)+ACMP*DEFDS(I+3)/3.
      DO 7 I=1,6
      DO 7 J=4,6
      ELAST(I,J)=C(I,J)/2.0
    7 CONTINUE
C
      AM78=(1.0-AL17/AL18)*CM
       IF(ICSS.EQ.0) THEN
        AML4=-3.0/AL4*CM*A2M2
       ELSE
        AML4=0.0D0
       ENDIF
C  
      AL19=AL19/3.0
      DO 8 I=1,3
      DO 8 J=1,3
      DO 26 K=1,3
       IF(J.EQ.K) THEN
       TT=DVT
       ELSE
       TT=-TR
       ENDIF
   26  ELAST(I,J)=ELAST(I,J)+AML4*EPSD(K)*TT
       ELAST(I,J+3)=ELAST(I,J+3)+AML4*EPSD(J+3)*0.5
       ELAST(I,J)=ELAST(I,J)+AL19
    8 CONTINUE
C
      DO 30 I=1,6
      DO 30 J=I,6
30    ELAST(I,J)=(ELAST(I,J)+ELAST(J,I))/2.
C
      DO 31 I=1,6
      DO 31 J=1,I
31    ELAST(I,J)=ELAST(J,I)
      RETURN
      END

