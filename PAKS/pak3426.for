C=======================================================================
C
C   VISKOELASTICNOST 3/D ELEMENT    (15.02.1998)
C
C=======================================================================
      SUBROUTINE D3M26(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
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
      LDEFT =LTAU   + 6
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6
C
      CALL TAUT26(PLAST(LTAU),PLAST(LDEFT),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),
     1            A(LFUN),A(LNTA),MAXT,MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUT26(TAUT,DEFT,TAU1,DEF1,
     1                  FUN,NTA,MAXT,MATE,TAU,DEF,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     INTEGRACIJA KONSTITUTIVNIH RELACIJA ZA VISKOELASTICAN MATERIJAL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      DIMENSION TAU(*),DEF(*),TAUT(*),DEFT(*),TAU1(*),DEF1(*)
      DIMENSION FUN(2,MATE,*),NTA(*),BDEF(6)
C
      IF(IRAC.EQ.2) RETURN
C
CE  E,V
C
      E=FUN(1,MAT,MAXT+1)
      V=FUN(2,MAT,MAXT+1)
      NTAC=NTA(MAT)
      EK=0.D0
      DO 10 K=1,NTAC
         BK=FUN(1,MAT,K)
         AK=FUN(2,MAT,K)
         DUM=-AK*VREME
         EK=EK+BK*(1.-EXP(DUM))
   10 CONTINUE
      E=E*EK
      CALL MELK26(E,V)
      DO 20 I=1,6
         BDEF(I)=(DEF(I)-DEFT(I))/DT
         TAU(I)=TAUT(I)
   20 CONTINUE
      CALL MNOZI1(TAU,ELAST,BDEF,6,6)
      DO 290 I=1,6
         DEF1(I)=DEF(I)
         TAU1(I)=TAU(I)
  290 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE MELK26(E,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
CE     FORM ( ELAST ) MATRIX
C
CE  NULL ( ELAST )
      DO 15 I=1,6
      DO 15 J=I,6
   15 ELAST(I,J)=0.0D0
C
      ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
      ELAST(2,2)=ELAST(1,1)
      ELAST(3,3)=ELAST(1,1)
      ELAST(1,2)=ELAST(1,1)*V/(1.-V)
      ELAST(1,3)=ELAST(1,2)
      ELAST(2,3)=ELAST(1,2)
      ELAST(4,4)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
      ELAST(5,5)=ELAST(4,4)
      ELAST(6,6)=ELAST(4,4)
      DO 50 I=1,6
      DO 50 J=I,6
   50 ELAST(J,I)=ELAST(I,J)
      RETURN
      END
