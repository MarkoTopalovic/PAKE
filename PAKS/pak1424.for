C=======================================================================
C
C  MODEL 24 : VISKOZNI MODEL ZA STAP
C             SILA PROPORCIONALNA RELATIVNOJ BRZINI IZMEDJU CVOROVA 
C             STAPA
C
C               D1M24
C               TAU124
C
C======================================================================
      SUBROUTINE D1M24(DEF,TAU,NLM,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA POZIVANJE PODPROGRAMA ODREDJIVANJE VISKOZNE SILE
C     (NAPONA) I KRUTOSTI
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
C
      LFUN=MREPER(1)
C      write(3,4020) (mreper(i),i=1,4)
C 4020 format(/' d1m24: mreper(i)=',4i5)
C      write(3,4021) lfun
C 4021 format(' lfun=',5i5) 
C
      NLL=NE*IDVA
      NLM1=(NLM-1)*IDVA
C
      LTAUT =LPLAST+NLM1
C      LEPST=LPLAST+NLL+NLM1
      LTAUVT = LTAUT + NLL
      LEPST = LTAUVT + NLL
      LEPSVT = LEPST + NLL
      LFT = LEPSVT + NLL  
C
      LTAUT1  =LPLAS1+NLM1
C      LEPST1 =LPLAS1+NLL+NLM1
      LTAUV1 = LTAUT1 + NLL
      LEPST1 = LTAUV1 + NLL
      LEPSV1 = LEPST1 + NLL
      LFT1 = LEPSV1 + NLL  
C      write(3,4010)ltaut,ltauvt,lepst,lepsvt,lft,ltaut1,ltauv1,
C     1  lepst1,lepsv1,lft1
C 4010 format(' reperi u d1m24: ltaut,ltauvt,lepst,lepsvt,lft=',5i5/
C     1 ' ltaut1,ltauv1,lepst1,lepsv1,lft1=',5i5)
C
      CALL TAU124(DEF,TAU,A(LTAUT),A(LTAUVT),A(LEPST),A(LEPSVT),A(LFT),
     1            A(LTAUT1),A(LTAUV1),A(LEPST1),A(LEPSV1),A(LFT1),
     2            A(LFUN),IRAC)
      RETURN
      END
C======================================================================
      SUBROUTINE TAU124(DEF,TAU,TAUT,TAUVT,EPST,EPSVT,FT,
     1                  TAUT1,TAUV1,EPST1,EPSV1,FT1,FUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     RACUNANJE SILE PO JEDINICI POVRSINE : ELASTICNA + VISKOZNI OTPOR
C     MODEL = 24 ZA STAP
C
C        DEF    = TEKUCA DEFORMACIJA STAPA = EPST1
C        TAU    = TEKUCI NAPON = TAUT1
C        TAUT   = NAPON NA KRAJU PRETHODNOG KORAKA
C        EPST   = DEFORMACIJA NA KRAJU PRETHODNOG KORAKA
C        TAUVT  = VISKOZNI NAPON NA KRAJU PRETHODNOG KORAKA
C        EPSVT  = BRZINA DEFORMACIJE NA KRAJU PRETHODNOG KORAKA
C        FT     = REAKCIJA STAPA NA KRAJU PRETHODNOG KORAKA
C        TAUV1,EPSV1,FT1 = TAUV,EPSV,F - NA KRAJU KORAKA
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /VISK1/ APR
C
      DIMENSION FUN(5,*)
C
C      write(3,4000) irac,mat,dt,(fun(i,mat),i=1,4)
C 4000 format(/' tau124: irac,mat,dt=',2i5,d15.7/' fun(1,2,3;mat)=',
C     1 4d15.7)
      TOL = 1.D-10
      ITYPE = FUN(1,MAT)
      EE = FUN(2,MAT)
      BB = FUN(3,MAT)
      E1 = FUN(4,MAT)
C      write(3,4053) itype
C 4053 format(' itype =',i5)
      IF (IRAC.EQ.0) RETURN
      BBT = BB/DT
      IF  (ITYPE.EQ.1) THEN
          IF (BBT.LT.TOL .OR. EE.LT.TOL) THEN
             WRITE (IZLAZ,2000) EE,BB,DT
             STOP
          ENDIF
          ELAST(1,1) = 1./(1./BBT+1./EE)
      ELSE
          IF(ITYPE.EQ.2) ELAST(1,1)= EE + BBT
      ENDIF
      IF(ITYPE.NE.3) GO TO 10
      IF  (EE.LT.TOL.OR.E1.LT.TOL) THEN
          WRITE(IZLAZ,2001) EE,E1
          STOP
      ENDIF
      TAUEPS = BB/E1
      TAUSIG = BB*(1./EE+1./E1)
      DEN = DT+TAUEPS
      ELAST(1,1) = EE*(DT+TAUSIG)/DEN
   10 IF (IRAC.EQ.2) RETURN
C      write(3,4070)bb,e1,ee,taueps,tausig 
C 4070 format(' bb,e1,ee=',3d15.7/' taueps,tausig=',2d15.7)      
C      write(3,4001) elast(1,1)
C 4001 format(' elast(1,1)=',d15.7)
C
      IF (KOR.EQ.1) THEN
         EPST = 0.D0
         TAUT = 0.D0
      ENDIF
      EPST1 = DEF
      EPSV1 = (EPST1-EPST)/DT
C
      GO TO (100,200,300),ITYPE
C     MAKSVELOV MODEL
  100 TAUT1 = ELAST(1,1)*(EPSV1*DT+TAUT/EE)
          TAUV1 = TAUT1
      GO TO 400
C     VOIGTOV MODEL
  200 TAUT1 = EE*EPST1
      TAUV1 = BB*EPSV1 
      TAUT1 = TAUT1 + TAUV1
      GO TO 400
C     KELVINOV MODEL 
  300 TAUV1 = BB*EPSV1
      TAUT1 = ELAST(1,1)*EPST1 + (TAUEPS*TAUT-TAUSIG*EE*EPST)/DEN
C
  400 TAU = TAUT1
      FT1 = TAU*APR
C      write(3,4003)epst1,epst,epsv1,taut1,tauv1,tau,apr,ft1
C 4003 format('epst1,epst,epsv1=',3d15.7/
C     1 ' taut1,tauv1,tau=',3d15.7/' apr,ft1=',2d15.7 )
C
C      write(3,4004)epst,epsvt,taut,tauvt,ft
C 4004 format(' vraca se iz tau1424:'/'epst,epsvt=',2d15.7/
C     1 ' taut,tauvt,ft=',3d15.7)
      RETURN
C
 2000 FORMAT(//' STOP U PAK1424 - MODUL ELASTICNOSTI ILI KOEFICIJENT'/
     1         '                  VISKOZNOG OTPORA SU ISUVISE MALI'/
     2         '                  ZA MAXWELL-OV MODEL'/
     2         '                  MORAJU BITI:  E.GT.1.D-10'/
     3         '                                (B/DT).GT.1.D-10'//
     4         '  E  =',D12.5/
     5         '  B  =',D12.5/
     6         ' DT  =',D12.5)                                
 2001 FORMAT(//' STOP U PAK1424 - JEADN OD MODULA ELASTICNOSTI JE',
     1                            ' MALI'/
     2         '                  ZA KELVINOV-OV MODEL'/
     2         '                  MORAJU BITI:  ER.GT.1.D-10'/
     3         '                                E1.GT.1.D-10'//
     4         '  ER  =',D12.5/
     5         '  E1  =',D12.5)
C
      END
