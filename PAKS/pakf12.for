C=======================================================================
C=======================================================================
C    SUBROUTINE  STRM2D
C                STREM2
C                FILLID
C                STRUNV
C=======================================================================
C=======================================================================
      SUBROUTINE STRM2D(A,GNODE,CCORD,NEL,NDIM,NET,LMAX,NPT,NTOT,NASLOV,
     &KOR,VREME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /BROJFK/ INDFOR,NULAZ
      COMMON /CDEBUG/ IDEBUG

      CHARACTER*250 NASLOV
      DIMENSION GNODE(2,5,*),NEL(NDIM,*)
      DIMENSION CCORD(3,*)
      dimension a(*)
      real a
      
      IIZLAZ=56
      IUNV=59
      MAXVEC=NTOT
      NETIP=2
      ISRPS=1
      INDFOR=2
      IDEBUG=0


      LMAXOL=LMAX


C      NEQ=3*NPT
C==========================================================================
      CALL MEMORY(LID,LMAX,3*NPT,1,MAXVEC,IIZLAZ)
      CALL FILLID(A(LID),NPT,NEQ)
      CALL MEMORY(LMHT,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LMAXA,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)

      CALL MAXATE(A(LMAXA),A(LMHT),A(LID),NEL,NET,NDIM,NEQ,NWK,3)
C      CALL MAXATF(A(LMAXA),A(LMHT),A(LID),NEL,NET,NDIM,NEQ,NWK,3,NDIM)

      LMAX0=LMAX
      NDES=3*NDIM
      MEM=0.5*NDES*(NDES+1)
      CALL MEMORY(LSKE,LMAX,MEM,2,NTOT,IIZLAZ)
      CALL MEMORY(LSKEF,LMAX,NDES*NDES,2,NTOT,IIZLAZ)
      CALL MEMORY(LALEVO,LMAX,NWK,2,NTOT,IIZLAZ)
C      CALL MEMORY(LDESNO,LMAX,NWK,2,NTOT,IIZLAZ)
      CALL MEMORY(LSILE,LMAX,NEQ,2,NTOT,IIZLAZ)
      CALL MEMORY(LTT1,LMAX,NEQ,2,NTOT,IIZLAZ)

      CALL CLEARR(A(LMAX0),LMAX-LMAX0+1)
      
     
      CALL STREM2(GNODE,A(LID),CCORD,NEL,A(LMAXA),A(LTT1),A(LALEVO),
     &A(LSKE),A(LSILE),A(LSKEF),NPT,NET,NDIM,NDES,NEQ,IIZLAZ,NWK)

      CALL STRUNV(A(LID),A(LTT1),59,NASLOV,KOR,VREME,NPT)

      LMAX=LMAXOL
C End of subroutine STRM2D
      END
C=======================================================================
C=======================================================================
      SUBROUTINE STREM2(GNODE,ID,CORD,NEL,MAXA,TT1,ALEVO,SKE,SILE,
     &SKEF,NPT,NET,NDIM,NDES,NEQ,IIZLAZ,NWK)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION TT1(*),ALEVO(*),SKE(*),SILE(*),CORD(3,*),SKEF(NDES,*)
      DIMENSION GNODE(2,5,*)
      DIMENSION NEL(NDIM,*),ID(3,*),MAXA(*)

      DIMENSION X(9),Y(9),TT21(3*9),B(3,3*9),BT(3*9,3),BS(2,3*9)
C      DIMENSION A11(9,9),ADX(9,9),ADY(9,9)
C      DIMENSION ADXX(9,9),ADYY(9,9)
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION LM2(3*9)
      DIMENSION H(9),ZVHX(9),ZVHY(9),F36(3*9),DT(3,3),DTS(2,2)
      



      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556


      R(2,1)=0.57735026918963
      R(2,2)=-0.5773502691896
      S(2,1)=0.57735026918963
      S(2,2)=-0.5773502691896
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      W(1,1)=2.0




      ITER=0
      EPSTR=1.D-2
      ISRPS=1


C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIM.EQ.4) IBRGT=2

      IAXIS=3


C      EEE=10.92D0
C      ANI=0.3D0
C      THICK=0.01D0
      EEE=1.D0
      ANI=0.D0
      THICK=1.D0

      EMNOZI=(EEE*THICK**3)/(12.D0*(1.D0-ANI*ANI))
      DT(1,1)=EMNOZI
      DT(1,2)=ANI*EMNOZI
      DT(1,3)=0.D0
      DT(2,1)=DT(1,2)
      DT(2,2)=1.D0*EMNOZI
      DT(2,3)=0.D0
      DT(3,1)=DT(1,3)
      DT(3,2)=DT(2,3)
      DT(3,3)=EMNOZI*(1.D0-ANI)/2.D0

      DTS(1,1)=(EEE*THICK)/(2.D0*(1.D0+ANI))
      DTS(2,2)=DTS(1,1)
      DTS(1,2)=0.D0
      DTS(2,1)=0.D0

 100  CALL CLEAR(SILE,NEQ) 
      CALL CLEAR(ALEVO,NWK) 
C      CALL CLEAR(DESNO,NWK) 

CE MAIN LOOP OVER ELEMENTS
      DO 400 NBREL=1,NET

C      DO I=1,NDIM
C       DO J=1,NDIM
C        A11(I,J)=0.D0
C        ADX(I,J)=0.D0
C        ADY(I,J)=0.D0
C        ADXX(I,J)=0.D0
C        ADYY(I,J)=0.D0
C       ENDDO
C      ENDDO

	DO I=1,3*NDIM
         LM2(I)=0
        ENDDO


      DO 125 I=1,NDES
C      TT210(I)=0.D0
 125  TT21(I)=0.D0
C=========================================================================
      JJ=-2
       DO 130 KLM=1,NDIM
        NODE=NEL(KLM,NBREL)
        X(KLM)=CORD(1,NODE)
        Y(KLM)=CORD(2,NODE)
       JJ=JJ+3
       JX=ID(1,NODE)
       JY=ID(2,NODE)
       JZ=ID(3,NODE)
       LM2(JJ)=JX
       LM2(JJ+1)=JY
       LM2(JJ+2)=JZ
C        JU=ID(1,NODE)
C        JV=ID(2,NODE)
C        JPSI=ID(3,NODE)
C        LM2(KLM)=JU
C        LM2(KLM+NDIM)=JV
C        LM2(KLM+2*NDIM)=JPSI
C        IF (JU.NE.0) TT21(KLM)=TT1(JU)
C        IF (JV.NE.0) TT21(KLM+NDIM)=TT1(JV)
C        IF (JPSI.NE.0) TT21(KLM+2*NDIM)=TT1(JPSI)
 130  CONTINUE

C INCIJALIZACIJA MATRICE SKEF I F36
      DO 260 I=1,NDES
      DO 258 J=1,NDES
       SKEF(I,J)=0.D0
 258  CONTINUE
       F36(I)=0.D0
 260  CONTINUE

       DO I=1,NDES*(NDES+1)*0.5
        SKE(I)=0.D0
       ENDDO

C INTEGRACIJA U GAUSOVIM TACKAMA
      NGAUS=0
      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT
       NGAUS=NGAUS+1
      CALL INTERS(R(IBRGT,I),S(IBRGT,J),TT21,DETJ,X,Y,B,BT,NDIM,NBREL,
     &ZVHX,ZVHY,H,BS)
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

C      DO 165 K=1,NDIM
C      DO 164 N=1,NDIM
C Calculation of stifness matrix Kuu
C       A11(K,N)=A11(K,N)+H(K)*H(N)*WDT
C       ADXX(K,N)=ADXX(K,N)+ZVHX(K)*ZVHX(N)*WDT
C       ADYY(K,N)=ADYY(K,N)+ZVHY(K)*ZVHY(N)*WDT
C       ADX(K,N)=ADX(K,N)+H(K)*ZVHX(N)*WDT
C       ADY(K,N)=ADY(K,N)+H(K)*ZVHY(N)*WDT
C  164 CONTINUE
C  165 CONTINUE


      DO 165 K=1,NDIM*3
      DO 164 N=1,NDIM*3

C Calculation of stifness matrix Kuu
      DO II=1,IAXIS
      DO JJ=1,IAXIS
       SKEF(K,N)=SKEF(K,N)+B(II,K)*DT(II,JJ)*B(JJ,N)*WDT
      ENDDO
      ENDDO

      DO II=1,2
      DO JJ=1,2
       SKEF(K,N)=SKEF(K,N)+BS(II,K)*DTS(II,JJ)*BS(JJ,N)*WDT
      ENDDO
      ENDDO


  164 CONTINUE
  165 CONTINUE

  170 CONTINUE
  180 CONTINUE


C      DO 305 I=1,NDIM
C      DO 300 J=1,NDIM
C       SKEF(I,J)=-A11(I,J)
C       SKEF(I,J+2*NDIM)=ADY(I,J)
C       SKEF(I+NDIM,J+NDIM)=A11(I,J)
C       SKEF(I+NDIM,J+2*NDIM)=ADX(I,J)
C       SKEF(I+2*NDIM,J)=ADX(I,J)
C       SKEF(I+2*NDIM,J+NDIM)=ADY(I,J)
C 300  CONTINUE
C 305  CONTINUE


C       DO 320 I=1,NDES
C        DO 320 J=1,NDES
C        F36(I)=F36(I)-SKEF(I,J)*TT21(J)
C 320   CONTINUE

       CALL PSKEFN(SKEF,SKE,NDES)	  
       CALL MATSTE (ALEVO,MAXA,SILE,SKE,F36,LM2,NDES,1)


C=========================================================================

C      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F36,MAXA,LM2,NDES,1)

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE



C PLATE BENDING EXAMPLE 1 - HINTON,OWEN
C       SILE(ID(1,5))=0.25D0
C
      RADIUS=1.D0

      DO 405 NODE=1,NPT
C==========================================
C FOR AXI-SYMMETRIC FLOW:
      IF (INDAX.EQ.1) RADIUS=CORD(1,NODE)
C==========================================
       JX=ID(2,NODE)       
       IF (JX.NE.0) SILE(JX)=GNODE(2,2,NODE)*1.0D35
         ALEVO(MAXA(JX))=1.0D35
       JY=ID(3,NODE)       
       IF (JY.NE.0) SILE(JY)=-GNODE(2,1,NODE)*1.0D35
         ALEVO(MAXA(JY))=1.0D35
405   CONTINUE

C       DO 2423 I=1,NWK
C 2423    WRITE(IIZLAZ,*)I,'ALEVO',ALEVO(I)

C      DO 424 I=1,NEQ
C 424    WRITE(IIZLAZ,*)I,'SILE',SILE(I)


      CALL RESENF(ALEVO,SILE,MAXA,NEQ,NWK,1)
      CALL RESENF(ALEVO,SILE,MAXA,NEQ,NWK,2)

C        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,NEQ,1)
C        CALL UACTCF(ALEVO,DESNO,SILE,MAXA,NEQ,2)

       DO I=1,NEQ
        TT1(I)=TT1(I)+SILE(I)
       ENDDO 
      
C      WRITE(IIZLAZ,*)'RESENJA BRE:'
C      DO I=1,NPT
C       WRITE(IIZLAZ,3000) I,TT1(ID(1,I)),TT1(ID(2,I)),TT1(ID(3,I))
C      ENDDO
C 3000   FORMAT(I5,3D13.5)

      WRITE(*,*)'ITER_S= ',ITER
       CALL KONVTF(TT1,SILE,KONVV1,1,ID,ITER,NPT,EPSTR,ISRPS,3)
       CALL KONVTF(TT1,SILE,KONVV2,2,ID,ITER,NPT,EPSTR,ISRPS,3)
       CALL KONVTF(TT1,SILE,KONVV3,3,ID,ITER,NPT,EPSTR,ISRPS,3)
       IF (KONVV1*KONVV2*KONVV3.EQ.0) THEN
        ITER=ITER+1
C        GO TO 100
       ENDIF

      END

C=======================================================================
C=======================================================================
      SUBROUTINE FILLID(ID,NPT,NEQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ID(3,*)


      DO I=1,NPT
       ID(1,I)=0
       ID(2,I)=0
       ID(3,I)=0
      ENDDO

C       ID(1,17)=1

C PLATE BENDING EXAMPLE 1 - HINTON,OWEN
C      ID(1,1)=1
C      ID(2,1)=1
C      ID(3,1)=1
C      ID(1,2)=1
C      ID(2,2)=1
C      ID(3,2)=0
C      ID(1,3)=1
C      ID(2,3)=1
C      ID(3,3)=0
C      ID(1,4)=1
C      ID(2,4)=0
C      ID(3,4)=1
C      ID(1,5)=0
C      ID(2,5)=0
C      ID(3,5)=0
C      ID(1,6)=0
C      ID(2,6)=1
C      ID(3,6)=0
C      ID(1,7)=1
C      ID(2,7)=0
C      ID(3,7)=1
C      ID(1,8)=0
C      ID(2,8)=0
C      ID(3,8)=1
C      ID(1,9)=0
C      ID(2,9)=1
C      ID(3,9)=1
C

      K=0
      DO I=1,NPT
       DO J=1,3
        IF (ID(J,I).EQ.0) THEN
         K=K+1
         ID(J,I)=K
        ELSE
         ID(J,I)=0
        ENDIF
       ENDDO
      ENDDO

      NEQ=K

      END
C=======================================================================
C=======================================================================
      SUBROUTINE STRUNV(ID,TT1,II,NASLOV,KOR,VREME,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ID(3,*)
      DIMENSION TT1(*)

      CHARACTER*250 NASLOV

      IND1=-1


C
C  FOR PRINTING STREAMLINES
C
      WRITE(II,5100) IND1
      WRITE(II,5100) 55
      WRITE(II,5003) NASLOV
      WRITE(II,7006)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,35,2,1
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) VREME
      DO 110 I=1,NPT
       WRITE(II,3006) I
       JJ=ID(1,I)
       PSI=0.D0
       IF (JJ.NE.0) PSI=TT1(ID(1,I))
       WRITE(II,5202) PSI
  110 CONTINUE
      WRITE(II,5100) IND1


 3006 FORMAT(6I10)
 5100 FORMAT(I6)
 5202 FORMAT(1PE13.2)
 5003 FORMAT(A80)
 7006 FORMAT(' NODAL STREAMLINES')
 2001 FORMAT(' DATE'/
     1       ' EMPTY'/
     1       ' LOAD CASE',I10)


      END
