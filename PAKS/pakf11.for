C==========================================================================
C==========================================================================
C    SUBROUTINE  MOVMSH
C                MOVE2D
C                SOL2D
C                INTERS
C                MATSTE
C                RESENF
C                ULAZS1
C                MAXATE
C                PSKEFN
C                PRNTSS
C                ZADNOS
C==========================================================================
C==========================================================================
      SUBROUTINE MOVMSH(A,CCORD,NZAD,ZADVRE,NEL,IDSS,IDENT,IDS,IDF,
     &BRZ,TT1S,CORDF,NDIM,NET,LMAX,NP,NTOTF,NPTF,NETIP,NPTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IDENT(2,*),IDS(NP,*),IDF(4,*),NEL(NDIM,*)
      DIMENSION BRZ(*),TT1S(*)
      DIMENSION CCORD(3,*),CORDF(3,*)
      DIMENSION IDSS(3,*)
	DIMENSION V(3),D(3),NQS(3),NQF(3),DD(3),VOLD(3)
      DIMENSION NZAD(3,*),ZADVRE(*)
 
      DIMENSION A(*)    
      REAL A

      
      K=0
      DO 10 NODEI=1,NPTI 
	    NODES=IDENT(1,NODEI)
        NODEF=IDENT(2,NODEI)

	  DO I=1,NETIP
C WE PUT NODEI INSTEAD NODES BECAUSE FREE NUMERATION NODES
        NQS(I)=IDS(NODES,I)
C        NQS(I)=IDS(NODEI+50,I)

        NQF(I)=IDF(I,NODEF)
        V(I)=0.D0
        D(I)=0.D0
        
        
        IF(NQS(I).NE.0) THEN
          IF (NSTAC.EQ.0) V(I)=BRZ(NQS(I))
          D(I)=TT1S(NQS(I))
        ENDIF
          DD(I)=D(I)-(CCORD(I,NODEF)-CORDF(I,NODEF))
          K=K+1
          ZADVRE(K)=DD(I)
          NZAD(1,K)=NODES
          NZAD(2,K)=I
          NZAD(3,K)=1
        ENDDO
 10    CONTINUE
   

      NUMZAD=K
      EEE=2.0D6
      ANI=0.3D0

      CALL MOVE2D(A,CCORD,NZAD,ZADVRE,NEL,IDSS,NDIM,NUMZAD,NET,LMAX,
     &NPTF,EEE,ANI,NTOTF)
      
      END
C==========================================================================
C=======================================================================
      SUBROUTINE MOVE2D(A,CCORD,NZAD,ZADVRE,NEL,ID,NDIM,NUMZAD,NET,LMAX,
     &NPT,EEE,ANI,NTOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      PARAMETER (NTOT = 100000000)
      COMMON /BROJFK/ INDFOR,NULAZ
      COMMON /CDEBUG/ IDEBUG

      DIMENSION NZAD(3,*),NEL(NDIM,*),ID(3,*)
      DIMENSION CCORD(3,*),ZADVRE(*)
      DIMENSION A(*)    
      REAL A
      
C      ISOLID=37
      IIZLAZ=32
      IUNV=33
      MAXVEC=NTOT
      
C           WRITE(2,*)'NTOTF',NTOT
C           WRITE(2,*)'MAXVEC',MAXVEC
      
      NETIP=2
      ISRPS=1
      INDFOR=2
      IDEBUG=0
C      EEE=2.0D6
C      ANI=0.3D0

C==========================================================================
C Opening the files
C      OPEN(ISOLID,FILE='SOLID.DAT')
C      OPEN(IIZLAZ,FILE='SOLID.OUT')
C      OPEN(IUNV,FILE='SOLID.UNV')
C==========================================================================
C Reading numberation of NODES,ID-matrix, COORDinates of nodes
C      REWIND (ISOLID)
C      READ(ISOLID,1000) NPT
C      KK=0
C      DO I=1,NPT
C       READ(ISOLID,1001)N,(ID(J,N),J=1,3)
C      DO JJ=1,3
C       IF (ID(JJ,N).EQ.0) THEN
C        KK=KK+1
C        ID(JJ,N)=KK
C       ELSE
C        ID(JJ,N)=0
C       ENDIF
C      ENDDO
C      ENDDO
C      NEQ=KK
C============================
C  1001   FORMAT(I5,1X,3I2)
      CALL MEMORY(LID,LMAX,3*NPT,1,MAXVEC,IIZLAZ)
C      CALL MEMORY(LCORD,LMAX,3*NPT,2,MAXVEC,IIZLAZ)
C      CALL ULAZS1(A(LID),A(LCORD),NPT,ISOLID,IIZLAZ,3,NEQ)
C==========================================================================
C Reading finite elements
C      READ(ISOLID,1000) NET,NDIM
C      CALL MEMORY(LNEL,LMAX,NDIM*NET,1,MAXVEC,IIZLAZ)
C      CALL ULAZF2(A(LNEL),NET,NETIP,NDIM,ISOLID,IIZLAZ,ISRPS)
C==========================================================================
C Reading prescribed values
C      READ(ISOLID,1000) NUMZAD
C      CALL MEMORY(LNZAD,LMAX,3*NUMZAD,1,MAXVEC,IIZLAZ)
C      CALL MEMORY(LZADVR,LMAX,NUMZAD,2,MAXVEC,IIZLAZ)
C      CALL ULAZF3(A(LNZAD),A(LZADVR),NUMZAD,ISOLID,IIZLAZ)
C==========================================================================
      CALL MEMORY(LMHT,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)
      CALL MEMORY(LMAXA,LMAX,NEQ+1,1,MAXVEC,IIZLAZ)
CZ ODBRAVIO
c      CALL MAXATE(A(LMAXA),A(LMHT),A(LID),NEL,NET,NDIM,NEQ,NWK,3)
C      CALL MAXATF(A(LMAXA),A(LMHT),A(LID),A(LNEL),NET,NDIM,NEQ,NWK,3,
C     &NDIM)

      LMAX0=LMAX
      NDES=2*NDIM
      MEM=0.5*NDES*(NDES+1)
      CALL MEMORY(LSKE,LMAX,MEM,2,NTOT,IIZLAZ)
      CALL MEMORY(LSKEF,LMAX,NDES*NDES,2,NTOT,IIZLAZ)
      CALL MEMORY(LALEVO,LMAX,NWK,2,NTOT,IIZLAZ)
CZ ODBRAVIO
      LDESNO=LMAX
C      CALL MEMORY(LDESNO,LMAX,NWK,2,NTOT,IIZLAZ)
C      CALL MEMORY(LGNODE,LMAX,NPT*3,2,NTOT,IIZLAZ)
      CALL MEMORY(LSILE,LMAX,NEQ,2,NTOT,IIZLAZ)
      CALL MEMORY(LTT1,LMAX,NEQ,2,NTOT,IIZLAZ)

      CALL CLEARR(A(LMAX0),LMAX-LMAX0+1)
      
      CALL SOL2D(CCORD,ID,CCORD,NEL,A(LMAXA),A(LTT1),
     &A(LALEVO),A(LSILE),NZAD,ZADVRE,A(LSKEF),A(LSKE),NPT,NET,
     &NDIM,NDES,NEQ,NUMZAD,EEE,ANI,A(LDESNO),IIZLAZ,NWK)


C       CALL RESEN(A(LALEVO),A(LSILE),A(LMAXA),NEQ,NWK,1)
C       CALL RESEN(A(LALEVO),A(LSILE),A(LMAXA),NEQ,NWK,2)
C  
C       CALL UACTCF(A(LALEVO),A(LDESNO),A(LSILE),A(LMAXA),NEQ,1)
C       CALL UACTCF(A(LALEVO),A(LDESNO),A(LSILE),A(LMAXA),NEQ,2)

       CALL PRNTSS(CCORD,A(LTT1),ID,NPT,IIZLAZ)


1000  FORMAT (2I5)

C End of program SOLID
      END
C=======================================================================
      SUBROUTINE SOL2D(GNODE,ID,CORD,NEL,MAXA,TT1,ALEVO,SILE,NZAD,
     &ZADVRE,SKEF,SKE,NPT,NET,NDIM,NDES,NEQ,NUMZAD,EEE,ANI,DESNO,
     &IIZLAZ,NWK)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION TT1(*),ALEVO(*),SILE(*),CORD(3,*),SKEF(NDES,*),SKE(*)
      DIMENSION GNODE(2,5,*),ZADVRE(*),DESNO(*)
      DIMENSION NEL(NDIM,*),ID(3,*),MAXA(*),NZAD(3,*)

      DIMENSION X(9),Y(9),TT21(18),B(4,18),BT(18,4)
      DIMENSION AKUU(18,18),F36(18),STRAIN(4),TAU(4)
      DIMENSION R(3,3),S(3,3),W(3,3),DT(4,4),H(9)
      DIMENSION LM2(18)
      



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

C      CALL ZADNOS(GNODE,ZADVRE,NZAD,NUMZAD)



C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIM.EQ.4) IBRGT=2


CE MAIN LOOP OVER ELEMENTS
      DO 400 NBREL=1,NET

	DO I=1,2*NDIM
         LM2(I)=0
        ENDDO

       IAXIS=4

      EMNOZI=EEE*(1.D0-ANI)/((1.D0+ANI)*(1.D0-2.D0*ANI))
      DT(1,1)=1.D0*EMNOZI
      DT(1,2)=ANI*EMNOZI/(1.D0-ANI)
      DT(1,3)=0.D0
      DT(2,1)=DT(1,2)
      DT(2,2)=1.D0*EMNOZI
      DT(2,3)=0.D0
      DT(3,1)=DT(1,3)
      DT(3,2)=DT(2,3)
      DT(3,3)=EMNOZI*(1.D0-2.D0*ANI)/(2.D0*(1.D0-ANI))
      DT(1,4)=DT(1,2)
      DT(2,4)=DT(1,2)
      DT(3,4)=0.D0
      DT(4,4)=EMNOZI*1.D0
      DT(4,1)=DT(1,4)
      DT(4,2)=DT(2,4)
      DT(4,3)=DT(3,4)

      DO 125 I=1,NDES
C      TT210(I)=0.D0
 125  TT21(I)=0.D0
C=========================================================================
      JJ=-1
       DO 130 KLM=1,NDIM
       NODE=NEL(KLM,NBREL)
       JJ=JJ+2
       DO NR=1,2
       IF (ID(NR,NEL(KLM,NBREL)) .NE. 0) THEN
        TT21(JJ+NR-1)=GNODE(2,NR,NODE)
C        TT21(JJ+NR-1)=TT1(ID(NR,NEL(KLM,NBREL)))
       ENDIF
      ENDDO
       X(KLM)=CORD(1,NODE)
       Y(KLM)=CORD(2,NODE)
       JX=ID(1,NODE)
       JY=ID(2,NODE)
       LM2(JJ)=JX
       LM2(JJ+1)=JY
       TT21(KLM)=GNODE(2,1,NODE)        
       TT21(KLM+NDIM)=GNODE(2,2,NODE)        
       IF (JX.NE.0) TT21(KLM)=TT1(JX)
       IF (JY.NE.0) TT21(KLM+NDIM)=TT1(JY)
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
      CALL INTERS(R(IBRGT,I),S(IBRGT,J),TT21,DETJ,X,Y,B,BT,NDIM,NBREL)

C      CALL INTERS(R(IBRGT,I),S(IBRGT,J),TT21,DETJ,X,Y,B,BT,NDIM,NBREL,
C     &ZVHX,ZVHY,H,BS)
     
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

C        CALL CLEAR(STRAIN,4) 
C        CALL CLEAR(TAU,4) 

    
C       DO II=1,IAXIS
C        DO NN=1,NDIM*2
C         STRAIN(II)=STRAIN(II)+B(II,NN)*TT21(NN)
C        ENDDO
C       ENDDO
      

      

      DO 165 K=1,NDIM*2
      DO 164 N=1,NDIM*2

C Calculation of stifness matrix Kuu
      DO II=1,IAXIS
      DO JJ=1,IAXIS
       SKEF(K,N)=SKEF(K,N)+BT(K,II)*DT(II,JJ)*B(JJ,N)*WDT
      ENDDO
      ENDDO
  164 CONTINUE
  165 CONTINUE
  170 CONTINUE
  180 CONTINUE








C=========================================================================

C       DO I=1,NDES
C        DO J=1,NDES
C         F36(I)=F36(I)-SKEF(I,J)*TT21(J)
C        ENDDO
C       ENDDO       


       CALL PSKEFN(SKEF,SKE,NDES)	  
       CALL MATSTE (ALEVO,MAXA,SILE,SKE,F36,LM2,NDES,1)

C      CALL ADDSTF(ALEVO,SILE,DESNO,SKEF,F36,MAXA,LM2,NDES,1)

C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE


      DO 405 I=1,NUMZAD
       JJ=ID(NZAD(2,I),NZAD(1,I))       
       IF (JJ.NE.0) THEN
        SILE(JJ)=ZADVRE(I)*1.0D35
        ALEVO(MAXA(JJ))=1.0D35
       ENDIF
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
      

      END

C=======================================================================
C======================================================================
C      SUBROUTINE INTERS(R,S,TT21,DETJ,X,Y,B,BT,NDIM,NBREL,ZVHX,ZVHY,H,
C     &BS)
      SUBROUTINE INTERS(R,S,TT21,DETJ,X,Y,B,BT,NDIM,NBREL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION X(9),Y(9),TT21(3*9),B(3,3*9),BT(3*9,3),BS(2,3*9)
      DIMENSION H(9),ZVHX(9),ZVHY(9),ZVHR(9),ZVHS(9),AJ(2,2)


      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
      IF (NDIM.EQ.8) THEN
        H(9)=0.0
      ELSE
        H(9)=RR*SS
      ENDIF
      H(8)=0.5*RP*SS-0.5*H(9)
      H(7)=0.5*RR*SM-0.5*H(9)
      H(6)=0.5*RM*SS-0.5*H(9)
      H(5)=0.5*RR*SP-0.5*H(9)
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP

      ENDIF



      IF (NDIM.GT.4) THEN
      IF (NDIM.EQ.8) THEN
        ZVHR(9)=0.D0
        ZVHS(9)=0.D0
      ELSE
        ZVHR(9)=-2.D0*R*SS
        ZVHS(9)=-2.D0*S*RR
      ENDIF
      ZVHR(8)=0.5*SS-0.5*ZVHR(9)
      ZVHS(8)=-RP*S-0.5*ZVHS(9)
      ZVHR(7)=-R*SM-0.5*ZVHR(9)
      ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
      ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
      ZVHS(6)=-RM*S-0.5*ZVHS(9)
      ZVHR(5)=-R*SP-0.5*ZVHR(9)
      ZVHS(5)=0.5*RR-0.5*ZVHS(9)
      ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
      ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
      ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
      ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
      ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
      ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)

      ELSE

      ZVHR(4)=0.25*SM
      ZVHS(4)=-0.25*RP
      ZVHR(3)=-0.25*SM
      ZVHS(3)=-0.25*RM
      ZVHR(2)=-0.25*SP
      ZVHS(2)=0.25*RM
      ZVHR(1)=0.25*SP
      ZVHS(1)=0.25*RP

      ENDIF

      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)


      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
       WRITE(*,*)'DETERMINANT LESS THEN ZERO!!!'
       WRITE(*,*)'FOR ELEMENT NUMBER ',NBREL
      STOP
      ENDIF

      DO 20 I=1,NDIM
       ZVHX(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
       ZVHY(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
  20  CONTINUE  

      

      DO J=1,3*NDIM
        BS(1,J)=0.D0
        BS(2,J)=0.D0
       DO I=1,3
        B(I,J)=0.D0
       ENDDO
      ENDDO


C BATHE:
C      JG=0
C      DO I=1,NDIM
C       IG=JG+1
C       BS(1,IG)=ZVHY(I)
C       BS(2,IG)=ZVHX(I)
C       IG=IG+1
C       JG=IG+1
C       BS(2,IG)=H(I)
C       BS(1,JG)=-H(I)
C        B(1,IG)=ZVHX(I)
C        B(3,IG)=ZVHY(I)
C        B(2,JG)=-ZVHY(I)
C        B(3,JG)=-ZVHX(I)
C      ENDDO


C HINTON:
      JG=0
      DO I=1,NDIM
       IG=JG+1
       BS(1,IG)=ZVHX(I)
       BS(2,IG)=ZVHY(I)
       IG=IG+1
       JG=IG+1
       BS(1,IG)=-H(I)
       BS(2,JG)=-H(I)
        B(1,IG)=-ZVHX(I)
        B(3,IG)=-ZVHY(I)
        B(2,JG)=-ZVHY(I)
        B(3,JG)=-ZVHX(I)
      ENDDO



      DO I=1,3
       DO J=1,3*NDIM
        BT(J,I)=B(I,J)
       ENDDO
      ENDDO

C End of subroutine INTERS
      END
C======================================================================


C=======================================================================
      SUBROUTINE MATSTE(SK,MAXA,F,SKE,FE,LM,NCV,INDSK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS  PAKOVANJE MATRICA I VEKTORA ELEMENATA U MATRICE I VEKTORE SISTEMA
CE  INSERTING ELEMENT MATRIXES AND VECTORS
CE  INTO SYSTEM MATRIXES AND VECTORS
      DIMENSION SK(*),MAXA(*),F(*),SKE(*),FE(*),LM(*)
C
      K=0
      DO 200 I=1,NCV
      IVR=LM(I)
      IF(INDSK.EQ.0) GO TO 110
      DO 100 J=I,NCV
      K=K+1
      KOL=LM(J)
C     IF(IVR.EQ.0.OR.KOL.EQ.0) GO TO 100
      IF(IVR.LE.0.OR.KOL.LE.0) GO TO 100
      IF(IVR-KOL) 10,10,20
   10 KS=MAXA(KOL) + KOL - IVR
      GO TO 50
   20 KS=MAXA(IVR)+IVR-KOL
   50 SK(KS)=SK(KS)+SKE(K)
C
  100 CONTINUE
C
C 110 IF(IVR.EQ.0) GO TO 200
  110 IF(IVR.LE.0) GO TO 200
      F(IVR)=F(IVR)+FE(I)
C
  200 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE RESENF(A,V,MAXA,NN,NWK,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C        WRITE(3,*) 'NN,NWK',NN,NWK
C        CALL IWRR(MAXA,NN,'JDIA')
C        CALL WRR(V,NN,'V   ')
C        CALL WRR(A,NWK,'AA  ')
C
CS     L*D*L(T) FAKTORIZACIJA
CE     L*D*L(T) FACTORIZATION
C
      IF(KKK-2)40,150,150
   40 DO 140 N=1,NN
      KN=MAXA(N)
      KL=KN+1
      KU=MAXA(N+1)-1
      KH=KU-KL
      IF(KH)110,90,50
   50 K=N-KH
      IC=0
      KLT=KU
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF(ND)80,80,60
   60 KK=MIN0(IC,ND)
      C=0.
      DO 70 L=1,KK
   70 C=C+A(KI+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
   90 K=N
      B=0.
      DO 100 KK=KL,KU
      K=K-1
      KI=MAXA(K)
C OVDE UBACENO:
      IF (A(KI).EQ.0.) A(KI)=1.D-15
C
      C=A(KK)/A(KI)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
  110 IF(A(KN).LE. 0.) THEN
        IF(ISRPS.EQ.0)
     *  WRITE(*,2000)N,A(KN)
        IF(ISRPS.EQ.1)
     *  WRITE(*,6000)N,A(KN)
C        STOP
      ENDIF
  140 CONTINUE
      RETURN
C
CS     REDUKOVANJE SLOBODNOG VEKTORA
CE     FORWARD REDUCTION
C
  150 DO 180 N=1,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL)180,160,160
  160 K=N
      C=0.
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
C
CS     ZAMENA UNAZAD
CE     BACK SUBSTITUTION
C
      DO 200 N=1,NN
      K=MAXA(N)
  200 V(N)=V(N)/A(K)
      IF(NN.EQ.1)RETURN
      N=NN
      DO 230 L=2,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF (KU-KL)230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
      RETURN
C
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2D20.12)
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITION'
     1//' ','PIVOT IS NOT POSITIVE FOR EQUATION NUM.',I4,//' ','PIVOT=',
     2D20.12)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZS1(ID,CORD,NPT,IULAZ,IIZLAZ,IDIM,NEQ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CHARACTER*250 ACOZ
      DIMENSION ID(IDIM,*),CORD(3,*)

      KK=0
      DO 10 I=1,NPT
      CALL ISPITF(ACOZ,IULAZ)
      READ(ACOZ,1006) N,(ID(II,N),II=1,IDIM),(CORD(J,N),J=1,3)
      DO JJ=1,IDIM
       IF (ID(JJ,N).EQ.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
        ID(JJ,N)=0
       ENDIF
      ENDDO

   10 CONTINUE
      NEQ=KK

C 1005 FORMAT(I5,4(3X,I2),3F10.6,2I2)
 1006 FORMAT(I5,1X,3(I2),2X,3F10.6,I5)
      WRITE(IIZLAZ,2000)
 2000 FORMAT(//
     *11X,'REDNI BROJEVI,OGRANICENJA I KOORDINATE CVOROVA'/)
      DO 12 I=1,NPT
      WRITE(IIZLAZ,1006) I,(ID(II,I),II=1,IDIM),(CORD(J,I),J=1,3)
   12 CONTINUE
      END
C==========================================================================
C==========================================================================
      SUBROUTINE MAXATE(MAXA,MHT,ID,NEL,NE,NTE,JEDN,NWK,IDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    PODPROGRAM ZA FORMIRANJE VEKTORA VISINA STUBOVA I MAXA
CS    KONACNO SE SMESTAJU U ISTI PROSTOR
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
C
C
      DIMENSION MAXA(*),MHT(*),NEL(NTE,*),LM(63),ID(IDIM,*)
C
CS    PETLJA PO ELEMENTIMA
CE    ELEMENT LOOP
C
      CALL IWRR(ID,66*IDIM,'ID  ')
      CALL IWRR(NEL,NE*NTE,'NEL    ')
      
      DO I=1,63
       LM(I)=0
      ENDDO

      DO I=1,JEDN+1
       MHT(I)=0
       MAXA(I)=0
      ENDDO

      DO 100 NLM=1,NE
         KK=0
         DO 2 I=1,NTE
            IF(NEL(I,NLM).EQ.0) GO TO 2
            N=NEL(I,NLM)
               DO 1 J=1,IDIM
                  IF(ID(J,N).LE.0) GO TO 1
                  KK=KK+1
                  LM(KK)=ID(J,N)
    1          CONTINUE
C            ENDIF
    2    CONTINUE
C
         LS=JEDN+1
         DO 10 I=1,KK
            IF (LM(I).LT.LS) LS=LM(I)
   10    CONTINUE

C
         DO 20 I=1,KK
            II=LM(I)
            ME=II-LS
            IF(ME.GT.MHT(II)) MHT(II)=ME
   20    CONTINUE
C
  100 CONTINUE
C
CS    VEKTOR MAXA
CE    VECTOR MAXA
C
      MAXA(1)=1
      MAXA(2)=2
      DO 200 I=2,JEDN
       MAXA(I+1)=MAXA(I)+MHT(I)+1
 200  CONTINUE
      NWK=MAXA(JEDN+1)-1
      LS = JEDN+1
      DO 210 I=1,LS
 210  MHT(I)=MAXA(I)

      RETURN
      END
C==========================================================================
C=========================================================================
      SUBROUTINE PSKEFN(SKEF,SKEFN,NDES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION SKEF(NDES,*),SKEFN(*)

	 K=0
       DO I=1,NDES
        DO J=I,NDES
          K=K+1
         	SKEFN(K)=SKEF(I,J)
        ENDDO
       ENDDO

      END
C==========================================================================
C=========================================================================
      SUBROUTINE PRNTSS(GNODE,TT1,ID,NPT,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(2,5,*),TT1(*),ID(3,*)


       WRITE(IIZLAZ,*)'RESULTS '
       DO NODE=1,NPT
         JX=ID(1,NODE)
         JY=ID(2,NODE)
         JZ=ID(3,NODE)
        IF (JX.NE.0) GNODE(2,1,NODE)=GNODE(2,1,NODE)+TT1(JX)
        IF (JY.NE.0) GNODE(2,2,NODE)=GNODE(2,2,NODE)+TT1(JY)
        IF (JZ.NE.0) GNODE(2,3,NODE)=GNODE(2,3,NODE)+TT1(JZ)
        WRITE(IIZLAZ,1000)NODE,(GNODE(2,J,NODE),J=1,3)
       ENDDO

1000   FORMAT(I5,3(D13.5))

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZADNOS(GNODE,ZADVRE,NZAD,NUMZAD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(3,*),ZADVRE(*)
       DIMENSION NZAD(3,*)

C
CE Subroutine ZADNOS is used for inclusion prescribed values
C

      
      DO 425 I=1,NUMZAD
        GNODE(NZAD(2,I),NZAD(1,I))=ZADVRE(I)
  425  CONTINUE 
      
      END