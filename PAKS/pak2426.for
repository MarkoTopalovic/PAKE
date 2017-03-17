C=======================================================================
C
C   VISKOELASTICNOST 2/D ELEMENT    (15.02.1998)
C
C=======================================================================
      SUBROUTINE D2M26(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MAXT=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 4*IDVA
      LTAUKT=LDEFT  + 4*IDVA
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 4*IDVA
      LTAUK =LDEFT1 + 4*IDVA
C
      CALL TAU226(A(LTAU),A(LDEFT),A(LTAU1),A(LDEFT1),A(LTAUKT),
     1            A(LTAUK),A(LFUN),A(LNTA),MAXT,MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAU226(TAUT,DEFT,TAU1,DEF1,TAUKT,TAUK,
     1                  FUN,NTA,MAXT,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     INTEGRACIJA KONSTITUTIVNIH RELACIJA ZA VISKOELASTICAN MATERIJAL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
C
      DIMENSION TAU(*),DEF(*),TAUT(*),DEFT(*),TAU1(*),DEF1(*)
      DIMENSION FUN(2,MATE,*),NTA(*),BDEF(4),TAUKT(*),TAUK(*)
C
      IF(IRAC.EQ.2) RETURN
C
CE  E,V
C
      WRITE(3,*) 'NLM,NGR,NGS,ITER,KOR',NLM,NGR,NGS,ITER,KOR
      WRITE(3,*) 'VREME,DT',VREME,DT
      E=FUN(1,MAT,MAXT+1)
      V=FUN(2,MAT,MAXT+1)
      NTAC=NTA(MAT)
C
      DO 20 I=1,KK
         BDEF(I)=DEF(I)-DEFT(I)
   20 CONTINUE
C
      EK=0.D0
      DO 10 K=1,NTAC
         BK=FUN(1,MAT,K)
         AK=FUN(2,MAT,K)
         DUM=-AK*DT
         EK=EK+BK*(1.-EXP(DUM))/AK
      WRITE(3,*) 'BK,AK,EK',BK,AK,EK
         LI=(K-1)*KK
      DO 10 I=1,KK
         TAUK(LI+I)=TAUKT(LI+I)*EXP(DUM)
   10 CONTINUE
      E=E*EK/DT
      CALL MEL226(E,V)
C
      WRITE(3,*) 'E,V',E,V
      CALL WRR(DEF,KK,'DEF ')
      CALL WRR(DEFT,KK,'DEFT')
      CALL WRR(BDEF,KK,'BDEF')
      CALL WRR(TAUT,KK,'TAUT')
C
      CALL CLEAR(TAU,KK)
      DO 40 K=1,NTAC
         LI=(K-1)*KK
         DO 30 I=1,KK
C            WRITE(3,*) 'ELAST',(ELAST(I,N),N=1,KK)
         DO 30 J=1,KK
            TAUK(LI+I)=TAUK(LI+I)+ELAST(I,J)*BDEF(J)
   30    CONTINUE
      DO 40 M=1,KK
          TAU(M)=TAU(M)+TAUK(LI+M)
   40 CONTINUE
         CALL WRR(TAU,KK,'TAU ')
C
      DO 290 I=1,KK
         DEF1(I)=DEF(I)
         TAU1(I)=TAU(I)
  290 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE MEL226(E,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      DO 15 I=1,4
      DO 15 J=I,4
   15 ELAST(I,J)=0.0D0
CS     RAVANSKO STANJE NAPONA
CE     PLANE STRESS
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
        ELAST(1,1)=E/(1.-V*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V
        ELAST(3,3)=ELAST(1,1)*(1.-V)/2.
      ELSE
CS     RAVANSKA DEFORMACIJA
CE     PLANE STRAIN
        ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V/(1.-V)
        ELAST(3,3)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
        IF(IETYP.NE.1) GO TO 40
        ELAST(4,4)=ELAST(1,1)
        ELAST(1,4)=ELAST(1,2)
        ELAST(2,4)=ELAST(1,2)
      ENDIF
C
   40 DO 50 I=1,4
      DO 50 J=I,4
   50 ELAST(J,I)=ELAST(I,J)
      RETURN
      END

