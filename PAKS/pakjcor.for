	
      SUBROUTINE  KOREKCIJA(IPRSLINA,PP,PC,A0,YA)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ******************************************************************
      IF(IPRSLINA.EQ.1) THEN
C        CENTRALNA PRSLINA
         YA=1.0
      ELSEIF(IPRSLINA .EQ.2) THEN
C        IVICNA PRSLINA
         YA=1.12
      ELSEIF(IPRSLINA.EQ.3) THEN
C        PRSLINA NA KRUZNOM OTVORU
         STEP=-1./6.
         YA=1.42*(A0/PP)**STEP
      ELSEIF(IPRSLINA.EQ.4) THEN
C        ELIPTICNO IVICNA PRSLINA
         YA=A0/PC
         YA=0.896
      ELSEIF(IPRSLINA.EQ.5) THEN
C        UNUTRASNJA ELIPTICNA PRSLINA
         YA=A0/PC
         YA=0.823   
      ELSEIF(IPRSLINA.EQ.6) THEN
C        UNUTRASNJA ELIPTICNA PRSLINA
         YA=1.0

      ENDIF
	
      RETURN
      END	   	   	   	       	  
