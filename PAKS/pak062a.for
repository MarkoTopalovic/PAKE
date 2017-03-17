C===============================================
CE    PRORAMED BY DRAKCE
C
CE    ITERSP - SOLVING SIMETRIC MATRIX A
CE    ITSPNE - SOLVING NONSIMETRICAL MATRIX A
C
C      (n+1)     1              (n)
C     Xj     = ----- (Rj - Aji*Xi  ), (i=1,j-1 i=j+1,NN) j=1,NN
C               Ajj
C
CE    This iterative solver can be applied if the matirix of coeficients
ce    is diagonal orieantaded
C===============================================
      SUBROUTINE ITERSP(A,MAXA,IROW,NN,NWK,X,R,NITER,TOL)
C=======================================================================
C
C     A() - MATRICA KOEFICIJENATA
C     MAXA() - VEKTOR INDEKSA GLAVNE DIJAGONALE
C     IROW() - VEKTOR INDEKSA VRSTE
C     NN - BROJ JEDNACINA
C     NWK - BROJ CLANOWA MATRICE A
C     X() - VEKTOR RESENJA
C     R() - VEKTOR DESNE STRANE
C     NITER - MAKSIMALAN BROJ DOZVOLJENIH ITERACIJA
C     TOL - TOLERANCIJA
C
C     V() - POMOCNI VEKTOR
C
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       SOLVING 
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(NWK),MAXA(NN+1),V(NN),IROW(NWK),X(NN),R(NN)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
	ITER=0
100	CONTINUE
        ITER=ITER+1
	  PRINT *, ' ITERATION NO. ',ITER
      DO 10 I=1,NN
        V(I)=X(I)
	  X(I)=R(I)
10	CONTINUE
      DO 20 J=1,NN
	  DO 30 I=MAXA(I)+1,MAXA(I+1)-1
	    X(J)=X(J)-A(I)*V(IROW(I))
	    X(IROW(I))=X(IROW(I))-A(I)*V(J)
30	  CONTINUE
20	CONTINUE
      XX=0.0D0
      DO 40 I=1,NN
	  X(I)=X(I)/A(MAXA(I))
	  XX=XX+(X(I)-V(I))*(X(I)-V(I))
40	CONTINUE
      XX=SQRT(XX)
	IF ((XX.GT.TOL).OR.(ITER.LT.NITER)) GOTO 100
      RETURN
      END

C===============================================
CE    PRORAMED BY DRAKCE
C===============================================
      SUBROUTINE ITSPNE(A,MAXA,IROW,NN,NWK,X,R,NITER,TOL)
C=======================================================================
C
C     A() - MATRICA KOEFICIJENATA
C     MAXA() - VEKTOR INDEKSA GLAVNE DIJAGONALE
C     IROW() - VEKTOR INDEKSA VRSTE
C     NN - BROJ JEDNACINA
C     NWK - BROJ CLANOWA MATRICE A
C     X() - VEKTOR RESENJA
C     R() - VEKTOR DESNE STRANE
C     NITER - MAKSIMALAN BROJ DOZVOLJENIH ITERACIJA
C     TOL - TOLERANCIJA
C
C     V() - POMOCNI VEKTOR
C
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       SOLVING 
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(2*NWK),MAXA(NN+1),V(NN),IROW(NWK),X(NN),R(NN)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
	ITER=0
100	CONTINUE
        ITER=ITER+1
	  PRINT *, ' ITERATION NO. ',ITER
      DO 10 I=1,NN
        V(I)=X(I)
	  X(I)=R(I)
10	CONTINUE
      DO 20 J=1,NN
	  DO 30 I=MAXA(I)+1,MAXA(I+1)-1
	    X(J)=X(J)-A(I+NWK)*V(IROW(I))
	    X(IROW(I))=X(IROW(I))-A(I)*V(I)
30	  CONTINUE
20	CONTINUE
      XX=0.0D0
      DO 40 I=1,NN
	  X(I)=X(I)/A(MAXA(I))
	  XX=XX+(X(I)-V(I))*(X(I)-V(J))
40	CONTINUE
      XX=SQRT(XX)
	IF ((XX.GT.TOL).OR.(ITER.LT.NITER)) GOTO 100
      RETURN
      END
