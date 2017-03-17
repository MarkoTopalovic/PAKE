C=======================================================================
      SUBROUTINE UMOD61(FUNMAT,MAT,KARTI,
     +                  GUST,NBLGR,IDEAS,MATG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 61; Isotropic Damage Model (Oliver 1996)
C .
CE.             MAT - MATERIAL NUMBER
C .
CE.   FUNMAT(1,MAT) - YOUNG'S MODULUS  -  E
CE.   FUNMAT(2,MAT) - POISSON'S RATIO  -  v
CE.   FUNMAT(3,MAT) - ULTIMATE STRESS   - Su 
CE.   FUNMAT(4,MAT) - FRACTURE ENERGY   - Gf 
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(5,*),AMAT(30)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD61'
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=1,2)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(I,MAT),I=3,4),idam
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(I,MAT),I=3,4),idam
      funmat(5,mat)=idam
C
      IF(NBLGR.GE.0) THEN
         CALL CLEAR(AMAT,30)
         AMAT(1)=FUNMAT(1,MAT)
         AMAT(7)=FUNMAT(2,MAT)
         AMAT(13)=GUST
         AMAT(15)=FUNMAT(3,MAT)
         AMAT(16)=FUNMAT(4,MAT)
         ISUMGR=MAT
         IF(IDEAS.EQ.8) THEN
            CALL MIDEAS(AMAT,ISUMGR,MAT,IGRAF) 
         ELSEIF(IDEAS.EQ.7) THEN
            CALL MIDEA7(AMAT,ISUMGR,MAT,IGRAF) 
         ENDIF
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
     1WRITE(IZLAZ,2010) (FUNMAT(J,MAT),J=1,4),idam
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) (FUNMAT(J,MAT),J=1,4),idam
      RETURN
C
 1000 FORMAT(2F10.0,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    61 (IZOTROPAN MODEL OSTECENJA)'///
     211X,'MATERIJAL BROJ =',I5)
 2010 FORMAT(//6X,
     1'M  A  T  E  R  I  J  A  L  N  E     K  O  N  S  T  A  N  T  E'//
     116X,'M O D U L     E L A S T I C N O S T I .  E =',1PD12.5//
     216X,'P O A S O N O V      O D N O S ........  V =',1PD12.5// 
     316X,'G R A N I C N I      N A P O N ........  Su=',1PD12.5// 
     416X,'E N E R G I J A      L O M A ..........  Gf=',1PD12.5//
     416X,'OMEKSANJE (=0-LINEARNO; =1-EKSPONENCIJALNO) idam=',i5//)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =   61 (Isotropic Damage Model)'///
     211X,'MATERIAL NUMBER =',I5)
 6010 FORMAT(//
     16X,'M  A  T  E  R  I  A  L      C  O  N  S  T  A  N  T  S'//
     116X,'Y O U N G S      M O D U L U S ........  E =',1PD12.5//
     216X,'P O I S S O N S      R A T I O ........  V =',1PD12.5// 
     316X,'U L T I M A T E      S T R E S S ......  Su=',1PD12.5// 
     416X,'F R A C T U R E      E N E R G Y ......  Gf=',1PD12.5//
     416X,'SOFTENING (=0-LINEAR; =1-EXPONENTIAL ) idam=',i5//)
C-----------------------------------------------------------------------
      END
