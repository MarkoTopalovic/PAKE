C=======================================================================
C
C  MODEL 24 : VISKOZNI MODEL ZA 3D
C             SILA PROPORCIONALNA RELATIVNOJ BRZINI U GL. PRAVCIMA
C
C               D3M24
C               TAU324
C
C======================================================================
      SUBROUTINE D3M24(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     REPERI 
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' D3M24'
C
      LFUN=MREPER(1)
C
      LTAUT = LPOCG
      LEPST = LTAUT + 6
C
      LTAU1 = LPOC1
      LEPS1 = LTAU1 + 6
C
      CALL TAU324(DEF,TAU,PLAST(LTAUT),PLAST(LEPST),
     1                    PLAS1(LTAU1),PLAS1(LEPS1),
     2                    A(LFUN),IRAC)
      RETURN
      END
C======================================================================
      SUBROUTINE TAU324(DEF,TAU,TAUT,EPST,
     1                          TAU1,EPS1,FUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     RACUNANJE SILE PO JEDINICI POVRSINE : ELASTICNA + VISKOZNI OTPOR
C     MODEL = 24 ZA 3D
C
C        DEF    = TEKUCA DEFORMACIJA STAPA = EPST1
C        TAU    = TEKUCI NAPON = TAUT1
C        TAUT   = NAPON NA KRAJU PRETHODNOG KORAKA
C        EPST   = DEFORMACIJA NA KRAJU PRETHODNOG KORAKA
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION TAUT(*),EPST(*),TAU1(*),EPS1(*),TAU(*),DEF(*)
C
      DIMENSION FUN(5,*)
C
      TOL = 1.D-10
CC      ITYPE = FUN(1,MAT)
C... USVOJEN KELVINOV MODEL
      ITYPE = 3
C
      EE = FUN(2,MAT)
      BB = FUN(3,MAT)
      E1 = FUN(4,MAT)
      ANI= FUN(5,MAT)
      IF (IRAC.EQ.0) RETURN
      BBT = BB/DT
      IF  (ITYPE.EQ.1) THEN
          IF (BBT.LT.TOL .OR. EE.LT.TOL) THEN
             WRITE (IZLAZ,2000) EE,BB,DT
             STOP
          ENDIF
          ER = 1./(1./BBT+1./EE)
      ELSE
          IF(ITYPE.EQ.2) ER= EE + BBT
      ENDIF
      IF(ITYPE.NE.3) GO TO 10
      IF  (EE.LT.TOL.OR.E1.LT.TOL) THEN
          WRITE(IZLAZ,2001) EE,E1
          STOP
      ENDIF
      TAUEPS = BB/E1
      TAUSIG = BB*(1./EE+1./E1)
CC      DEN = DT+TAUEPS
CC      ER = EE*(DT+TAUSIG)/DEN
   10 IF (IRAC.EQ.2) RETURN
C
      GO TO (100,200,300),ITYPE
C     MAKSVELOV MODEL
  100 CONTINUE
      GO TO 400
C     VOIGTOV MODEL
  200 CONTINUE
      GO TO 400
C
C     KELVINOV MODEL ---- RADI SAMO OVO
C
300   CALL MEL3EL(ELAST,EE,ANI)
      CC1=-TAUSIG/(DT+TAUEPS)
      CALL JEDNAK(EPS1,EPST,CC1,6)
      CALL CLEAR(TAU1,6)
      CALL MNOZI1(TAU1,ELAST,EPS1,6,6)
      CC2=(DT+TAUSIG)/(DT+TAUEPS)
      CALL JEDNAK(ELAST,ELAST,CC2,36)
      CALL MNOZI1(TAU1,ELAST,DEF,6,6)
      CC3=TAUEPS/(DT+TAUEPS)
      CALL ZBIRM(TAU1,TAUT,CC3,6)
C
  400 CALL JEDNA1(TAU,TAU1,6)
      CALL JEDNA1(EPS1,DEF,6)
      RETURN
C
 2000 FORMAT(//' STOP U PAK3424 - MODUL ELASTICNOSTI ILI KOEFICIJENT'/
     1         '                  VISKOZNOG OTPORA SU ISUVISE MALI'/
     2         '                  ZA MAXWELL-OV MODEL'/
     2         '                  MORAJU BITI:  E.GT.1.D-10'/
     3         '                                (B/DT).GT.1.D-10'//
     4         '  E  =',D12.5/
     5         '  B  =',D12.5/
     6         ' DT  =',D12.5)                                
 2001 FORMAT(//' STOP U PAK3424 - JEADN OD MODULA ELASTICNOSTI JE',
     1                            ' MALI'/
     2         '                  ZA KELVINOV-OV MODEL'/
     2         '                  MORAJU BITI:  ER.GT.1.D-10'/
     3         '                                E1.GT.1.D-10'//
     4         '  ER  =',D12.5/
     5         '  E1  =',D12.5)
C
      END
