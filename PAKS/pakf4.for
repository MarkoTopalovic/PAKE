C==========================================================================
C==========================================================================
C    SUBROUTINE  JACT
C                STRANA
C                KOORD1
C                INTEF2
C                INTER1
C                KOORD
C                RAZL1
C                CROSS
C                TANNOR
C                CROSSC
C                CROSXY
C                VVNEW
C                CORNEW
C                FGRAF1
C                TGRAXY
C                STAGPP
C                STAGP2
C                TGRA22
C                ZPLANE
C                TGRA23
C                EXPAND
C                COLLEA
C                INITIA
C                RETPRI
C                FILLN
C                ZADNOD
C                SHEARS
C                VECMUL
C                SSTRES
C                PRTIME
C                MAXMEM
C==========================================================================
C==========================================================================
       SUBROUTINE JACT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION H(21),P(3,21),NEL(NDIM,*)
       DIMENSION XJ(3,3),LE(3),ME(3),CK(21,3),PJ(3,21),HP(8)
       DIMENSION TT21(*),V1(21),V2(21),V3(21),PP(8),XJJ(3,3),TEMP(21)
       DIMENSION VTANG(3,21),AN(3),VSHX(3),VSHY(3),VSHZ(3),SHEAR(*)
       DIMENSION IPERM(8),NOD9(13)
C
CE  Subroutine JACT is used for integration 3D finite elements
C
C       COMMON /VISKOZ/ AMI,INDAMI
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C       COMMON /SRPSKI/ ISRPS
C       COMMON /ULAZNI/ IULAZ,IIZLAZ
C       COMMON /AJVVV/ HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW
C       COMMON /UPWIND/ IUPWIN

C      H(21)           - Interolacione funkcije
C      P(3,21)         - Izvodi interpolacionih funkcija
C      NEL(NDIM,*)     - Koordinate cvorova elementa
C      XJ(3,3)         - Jakobijan
C      XJJ(3,3)        - Inverzni Jakobijan
C      LE(3)           - 
C      ME(3)           - 
C      CK(21,3)        - 
C      PJ(3,21)        - Matrica izvoda interpolacionih funkcija 
C                        po globalnim koordinatama
C      HP(8)           - 
C      TT21(*)         - 
C      V1(21)          - Brzina u X pravcu
C      V2(21)          - Brzina u Y pravcu
C      V3(21)          - Brzina u Z pravcu
C      PP(8)           - Pritisak
C      TEMP(21)        - Temperatura
C      VTANG(3,21)     - 
C      AN(3)           - 
C      VSHX(3)         - 
C      VSHY(3)         -  
C      VSHZ(3)         - 
C      SHEAR(*)        - 
C      IPERM(8)        - Oznacavanje cvorova od 1 do 8
C      NOD9(13)        - Oznacavanje cvorova od 9 do 21

       DATA IPERM/2,3,4,1,6,7,8,5/
       DATA NOD9/9,10,11,12,13,14,15,16,17,18,19,20,21/

       IELX=NDIM
       INDUP=0
       NND9=13
       CALL CLEAR(SHEAR,3)
C (zagrade kod interpolacionih funckija: P=+, M=-)
      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
      DO 82 I=1,21
      H(I)=0.
      DO 82 J=1,3
   82 P(J,I)=0.
C
C     INTERPOLACIJSKE FUNKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
C
CS    PRVIH 8 CVOROVA
CE    FIRST 4 NODES
C
      H(1)=0.125*RP*SP*TP
      H(2)=0.125*RM*SP*TP
      H(3)=0.125*RM*SM*TP
      H(4)=0.125*RP*SM*TP
      H(5)=0.125*RP*SP*TM
      H(6)=0.125*RM*SP*TM
      H(7)=0.125*RM*SM*TM
      H(8)=0.125*RP*SM*TM
      
 100  DO JJ=1,8
       HP(JJ)=H(JJ)
      ENDDO
C hi,1 za prvih 8 cvorova
      P(1,1)=0.125*SP*TP
      P(1,2)=-P(1,1)
      P(1,3)=-0.125*SM*TP
      P(1,4)=-P(1,3)
      P(1,5)=0.125*SP*TM
      P(1,6)=-P(1,5)
      P(1,7)=-0.125*SM*TM
      P(1,8)=-P(1,7)
C hi,2 za prvih 8 cvorova
      P(2,1)=0.125*RP*TP
      P(2,2)=0.125*RM*TP
      P(2,3)=-P(2,2)
      P(2,4)=-P(2,1)
      P(2,5)=0.125*RP*TM
      P(2,6)=0.125*RM*TM
      P(2,7)=-P(2,6)
      P(2,8)=-P(2,5)
C hi,3 za prvih 8 cvorova 
      P(3,1)=0.125*RP*SP
      P(3,2)=0.125*RM*SP
      P(3,3)=0.125*RM*SM
      P(3,4)=0.125*RP*SM
      P(3,5)=-P(3,1)
      P(3,6)=-P(3,2)
      P(3,7)=-P(3,3)
      P(3,8)=-P(3,4)
C
      IF (IELX.EQ.8) GO TO 50
C
CS    STEPENI SLOBODE ZA CVOROVE PREKO 8
CE    DEGREES OF FREADOM FOR NODES OVER 8
C
      I=0
    2 I=I+1
      IF (I.GT.NND9) GO TO 30
      NN=NOD9(I)-8
      GO TO (9,10,11,12,13,14,15,16,17,18,19,20,21),NN
C
    9 H(9)=0.25*RR*SP*TP
      P(1,9)=-0.50*R*SP*TP
      P(2,9)=0.25*RR*TP
      P(3,9)=0.25*RR*SP
      GO TO 2
   10 H(10)=0.25*RM*SS*TP
      P(1,10)=-0.25*SS*TP
      P(2,10)=-0.50*RM*S*TP
      P(3,10)=0.25*RM*SS
      GO TO 2
   11 H(11)=0.25*RR*SM*TP
      P(1,11)=-0.50*R*SM*TP
      P(2,11)=-0.25*RR*TP
      P(3,11)=0.25*RR*SM
      GO TO 2
   12 H(12)=0.25*RP*SS*TP
      P(1,12)=0.25*SS*TP
      P(2,12)=-0.50*RP*S*TP
      P(3,12)=0.25*RP*SS
      GO TO 2
   13 H(13)=0.25*RR*SP*TM
      P(1,13)=-0.50*R*SP*TM
      P(2,13)=0.25*RR*TM
      P(3,13)=-0.25*RR*SP
      GO TO 2
   14 H(14)=0.25*RM*SS*TM
      P(1,14)=-0.25*SS*TM
      P(2,14)=-0.50*RM*S*TM
      P(3,14)=-0.25*RM*SS
      GO TO 2
   15 H(15)=0.25*RR*SM*TM
      P(1,15)=-0.50*R*SM*TM
      P(2,15)=-0.25*RR*TM
      P(3,15)=-0.25*RR*SM
      GO TO 2
   16 H(16)=0.25*RP*SS*TM
      P(1,16)=0.25*SS*TM
      P(2,16)=-0.50*RP*S*TM
      P(3,16)=-0.25*RP*SS
      GO TO 2
   17 H(17)=0.25*RP*SP*TT
      P(1,17)=0.25*SP*TT
      P(2,17)=0.25*RP*TT
      P(3,17)=-0.50*RP*SP*T
      GO TO 2
   18 H(18)=0.25*RM*SP*TT
      P(1,18)=-0.25*SP*TT
      P(2,18)=0.25*RM*TT
      P(3,18)=-0.50*RM*SP*T
      GO TO 2
   19 H(19)=0.25*RM*SM*TT
      P(1,19)=-0.25*SM*TT
      P(2,19)=-0.25*RM*TT
      P(3,19)=-0.50*RM*SM*T
      GO TO 2
   20 H(20)=0.25*RP*SM*TT
      P(1,20)=0.25*SM*TT
      P(2,20)=-0.25*RP*TT
      P(3,20)=-0.50*RP*SM*T
      GO TO 2
   21 H(21)=RR*SS*TT
      P(1,21)=-2.0*R*SS*TT
      P(2,21)=-2.0*S*RR*TT
      P(3,21)=-2.0*T*RR*SS
      GO TO 2
C
CS    KOREKCIJE PRVIH 20 FUNKCIJA AKO JE UPOTREBLJEN CVOR 21
CE    CORECTION OF FIRST 20 FUNCTIONS IF NODE 21 EXISTS
C
   30 IN=NOD9(NND9)
      IF(IN.NE.21) GO TO 40
      DO 36 I=1,8
      H(I)=H(I)-0.125*H(21)
      DO 36 J=1,3
   36 P(J,I)=P(J,I)-0.125*P(J,21)
      IF(NND9.EQ.1) GO TO 51
      DO 37 I=1,NND9-1
      IN=NOD9(I)
      H(IN)=H(IN)-0.25*H(21)
      DO 37 J=1,3
   37 P(J,IN)=P(J,IN)-0.25*P(J,21)
C
CS    KOREKCIJE PRVIH 8 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 8
CE    CORECTION OF FIRST 8 FUNCTIONS IF NODES OVER 8 EXISTS
C
   40 IH=0
   41 IH=IH+1
      IF(IH.GT.NND9) GO TO 50
      IN=NOD9(IH)
      IF(IN.GT.16) GO TO 46
      I1=IN-8
      I2=IPERM(I1)
C
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 45 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   45 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   46 IF(IN.EQ.21) GO TO 51
      I1=IN-16
      I2=I1+4
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 47 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   47 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   51 H(NND9+8)=H(21)
      DO 39 J=1,3
   39 P(J,NND9+8)=P(J,21)

      
   50 HH=0.D0
      DO I=1,NDIM
       HH=HH+H(I)
      ENDDO
C      WRITE(IIZLAZ,*)'ZBIR H=',HH

      DO I=1,NDIM
C Brzina I-tog cvora u X-pravcu
      V1(I)=TT21(I)
C Brzina I-tog cvora u Y-pravcu
      V2(I)=TT21(I+NDIM)
C Brzina I-tog cvora u Z-pravcu
      V3(I)=TT21(I+2*NDIM)
C Temperatura I-tog cvora
      TEMP(I)=TT21(I+3*NDIM+8)
      IF (I.LE.8) THEN
C Pritisak I-tog cvora za  IPERM cvorove
        PP(I)=TT21(I+3*NDIM)
      ENDIF
      ENDDO
      
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C

      DO I=1,3
      DO J=1,3
       XJ(I,J)=0.D0
       XJJ(I,J)=0.D0
        DO  KK=1,NDIM
C Matrica Jakobijana
         XJ(I,J)=XJ(I,J)+P(I,KK)*CK(KK,J)
C Inicijalizacija Inverzne Matrice Jakobijana
         XJJ(I,J)=XJJ(I,J)+P(I,KK)*CK(KK,J)
        ENDDO    
C      WRITE(IIZLAZ,*)'XJ=',XJ(I,J)
      ENDDO    
      ENDDO    

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)

C Odredjivanje Inverzne Matrice Jakobijana
      CALL MINV(XJJ,3,DET1,LE,ME)

      IF (DET1.LT.1.D-15) THEN
       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
       WRITE(IIZLAZ,*)'NODES COORDINATES'
       DO I=1,NDIM
        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
       ENDDO
       STOP
      ELSE
C       WRITE(IIZLAZ,*)'DET1=',DET1
      ENDIF
        
 1000 FORMAT(I5,3(D13.5))

C     Matrica izvoda interpolacionih funkcija 
C     po globalnim koordinatama
      DO 85 I=1,3
      DO 85 JJ=1,NDIM
      PJ(I,JJ)=0.D0
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJJ(I,K)*P(K,JJ)
   85 CONTINUE


      IF (IUPWIN.EQ.1.AND.INDUP.EQ.0) THEN
       CALL INTER1(CK,V1,V2,V3,H,PJ)
       INDUP=1
       GOTO 100
      ENDIF
      

      HXU=0.D0
      HYU=0.D0
      HZU=0.D0
      HXV=0.D0
      HYV=0.D0
      HZV=0.D0
      HXW=0.D0
      HYW=0.D0
      HZW=0.D0
      ZVXT=0.D0
      ZVYT=0.D0
      ZVZT=0.D0
      DO L =1,NDIM
        HXU=HXU+PJ(1,L)*V1(L)
        HYU=HYU+PJ(2,L)*V1(L)
        HZU=HZU+PJ(3,L)*V1(L)
        HXV=HXV+PJ(1,L)*V2(L)
        HYV=HYV+PJ(2,L)*V2(L)
        HZV=HZV+PJ(3,L)*V2(L)
        HXW=HXW+PJ(1,L)*V3(L)
        HYW=HYW+PJ(2,L)*V3(L)
        HZW=HZW+PJ(3,L)*V3(L)
        ZVXT=ZVXT+PJ(1,L)*TEMP(L)
        ZVYT=ZVYT+PJ(2,L)*TEMP(L)
        ZVZT=ZVZT+PJ(3,L)*TEMP(L)
      ENDDO
      
      IF(KFIX.GT.0) GO TO 70
      RETURN

C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
  70   GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      ANX=R*(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))
      ANY=R*(XJ(3,1)*XJ(2,3)-XJ(3,3)*XJ(2,1))
      ANZ=R*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      ANX=-S*(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))
      ANY=S*(XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1))
      ANZ=-S*(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
      ANX=T*(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))
      ANY=-T*(XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1))
      ANZ=T*(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))
   74 DET=DSQRT(DET)
      ANX=ANX/DET
      ANY=ANY/DET
      ANZ=ANZ/DET
      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
C      WRITE(IIZLAZ,*)'ANX=',ANX
C      WRITE(IIZLAZ,*)'ANY=',ANY
C      WRITE(IIZLAZ,*)'ANZ=',ANZ
      CALL SHEARS(AN,TT21,VTANG,NDIM)
      
      V1X=0.D0
      V1Y=0.D0
      V1Z=0.D0
      V2X=0.D0
      V2Y=0.D0
      V2Z=0.D0
      V3X=0.D0
      V3Y=0.D0
      V3Z=0.D0
      VSHX(1)=0.D0
      VSHX(2)=0.D0
      VSHX(3)=0.D0
      VSHY(1)=0.D0
      VSHY(2)=0.D0
      VSHY(3)=0.D0
      VSHZ(1)=0.D0
      VSHZ(2)=0.D0
      VSHZ(3)=0.D0
      DO I=1,NDIM
       V1X=V1X+PJ(1,I)*V1(I)
       V1Y=V1Y+PJ(2,I)*V1(I)
       V1Z=V1Z+PJ(3,I)*V1(I)
       V2X=V2X+PJ(1,I)*V2(I)
       V2Y=V2Y+PJ(2,I)*V2(I)
       V2Z=V2Z+PJ(3,I)*V2(I)
       V3X=V3X+PJ(1,I)*V3(I)
       V3Y=V3Y+PJ(2,I)*V3(I)
       V3Z=V3Z+PJ(3,I)*V3(I)
       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,3)
      SHEAR(2)=DOT(AN,VSHY,3)
      SHEAR(3)=DOT(AN,VSHZ,3)
      PRIT=DOT(HP,PP,8)
C      SF1=-PRIT*ANX
C      SF2=-PRIT*ANY
C      SF3=-PRIT*ANZ
      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C=======================================================================
       SUBROUTINE STRANA(NEL,NDIM,NODES,NPOV)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
       DIMENSION NODES(5),NEL(NDIM,*),IPM(8,3),IBRP(3),IBRM(3)
C
CE Subroutine STRANA is used for finding which face(side) of 3D element 
CE is used for surface boundary conditions
C
C       DATA ITP/1,2,3,4/,ITM/5,6,7,8/
C       DATA ISP/1,5,6,2/,ISM/3,7,8,4/
C       DATA IRP/1,4,8,5/,IRM/3,2,6,7/

       DATA IPM/1,4,8,5,3,2,6,7,1,5,6,2,3,7,8,4,1,2,3,4,5,6,7,8/
C      DATA IPM/1,4,8,5,2,6,7,3,1,5,6,2,4,3,7,8,1,2,3,4,5,8,7,6/


        NBREL=NODES(1)
        DO K=1,3
        IBRP(K)=0
        IBRM(K)=0
       DO I=2,5
         DO J=1,4
          IF(NODES(I).EQ.NEL(IPM(J,K),NBREL)) IBRP(K)=IBRP(K)+1
          IF(NODES(I).EQ.NEL(IPM(J+4,K),NBREL)) IBRM(K)=IBRM(K)+1
         ENDDO      
       ENDDO      
        
       IF (IBRP(K).EQ.4) NPOV=2*K-1
       IF (IBRM(K).EQ.4) NPOV=2*K
 
       ENDDO
        
       END

C=======================================================================
C==========================================================================
       SUBROUTINE KOORD1(X,Y,H,CORD,NDIM)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION H(*),CORD(*),X(*),Y(*)
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

       DO I=1,2
        CORD(I)=0.D0
       ENDDO

       DO I=1,NDIM
          CORD(1)=CORD(1)+H(I)*X(I)
          CORD(2)=CORD(2)+H(I)*Y(I)
       ENDDO


      END
C=======================================================================
C==========================================================================
       SUBROUTINE INTEF2(X,Y,V1,V2,HH,ZVHX,ZVHY,NDIM,AKT,GUSM,AMI)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C       COMMON /VISKOZ/ AMI,INDAMI
C       COMMON /KONST/ GUSM(1),CC,AKT
C       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

        

       DIMENSION V1(*),V2(*),HH(*),H(9)
       DIMENSION RP1(3),RM1(3),SP1(3),SM1(3)
       DIMENSION RST0(3),EX(3),EY(3),X(*),Y(*),ZVHX(*),ZVHY(*)



      DO I=1,3
       RP1(I)=0.D0
       RM1(I)=0.D0
       SP1(I)=0.D0
       SM1(I)=0.D0
       RST0(I)=0.D0
       EX(I)=0.D0
       EY(I)=0.D0
      ENDDO
       
	DO 20 IBR=0,4
      R=0.D0
      S=0.D0

      IF (IBR.EQ.1) R=1.D0
      IF (IBR.EQ.2) S=1.D0
      IF (IBR.EQ.3) R=-1.D0
      IF (IBR.EQ.4) S=-1.D0

      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
      H(9)=RR*SS
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



C
C     INTERPOLACIJSKE FUNKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
       HV1=DOT(H,V1,NDIM)
       HV2=DOT(H,V2,NDIM)
  

      IF (IBR.EQ.0) CALL KOORD1(X,Y,H,RST0,NDIM)

      IF (IBR.EQ.1) CALL KOORD1(X,Y,H,RP1,NDIM)
      IF (IBR.EQ.2) CALL KOORD1(X,Y,H,SP1,NDIM)
      IF (IBR.EQ.3) CALL KOORD1(X,Y,H,RM1,NDIM)
      IF (IBR.EQ.4) CALL KOORD1(X,Y,H,SM1,NDIM)
     
20    CONTINUE       
       
      HX=RAZL(RP1,RM1)
      HY=RAZL(SP1,SM1)

      CALL RAZL1(RP1,RST0,EX)
      CALL RAZL1(SP1,RST0,EY)


       UX=EX(1)*HV1+EX(2)*HV2
       UY=EY(1)*HV1+EY(2)*HV2

       VV=HV1**2+HV2**2
C       CALL KOORD(CK,H,RST0)       
       IF (DABS(VV).LT.1.D-10) RETURN


	AK=2.D0*AMI/GUSM

      ALFAX=UX*HX/(AK)
      ALFAY=UY*HY/(AK)
      
      AX=0.D0
      AY=0.D0
      IF (ALFAX.NE.0.D0) AX=1.D0/DTANH(ALFAX)-1.D0/ALFAX
      IF (ALFAY.NE.0.D0) AY=1.D0/DTANH(ALFAY)-1.D0/ALFAY
     
C       AX=HUGHES(ALFAX)
C       AY=HUGHES(ALFAY)
C       AZ=HUGHES(ALFAZ)

           
      AKK=((AX*UX*HX+AY*UY*HY)/2.D0)



      DO I=1,NDIM
        PP=AKK*(HV1*ZVHX(I)+HV2*ZVHY(I))/VV
        HH(I)=HH(I)+PP
      ENDDO
      
            
      END
C=======================================================================
C==========================================================================
       SUBROUTINE INTER1(CK,V1,V2,V3,HH,PJ)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON /VISKOZ/ AMI,INDAMI
       COMMON /KONST/ GUSM,CC,AKT
       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

        

       DIMENSION CK(21,*),XJ(3,3),V1(*),V2(*),V3(*),HH(*),H(21),P(3,21)
       DIMENSION RP1(3),RM1(3),SP1(3),SM1(3),TP1(3),TM1(3),PJ(3,*)
       DIMENSION RST0(3),EX(3),EY(3),EZ(3)

C
CE Subroutine INTER1 is used for integration 3D finite element with UPWIND addition
C

      TOL=1.D-8
       
	DO 20 IBR=0,6
      R=0.D0
      S=0.D0
      T=0.D0

      IF (IBR.EQ.1) R=1.D0
      IF (IBR.EQ.2) S=1.D0
      IF (IBR.EQ.3) T=1.D0

      IF (IBR.EQ.4) R=-1.D0
      IF (IBR.EQ.5) S=-1.D0
      IF (IBR.EQ.6) T=-1.D0

      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
      DO 82 I=1,21
      H(I)=0.D0
   82 CONTINUE
C
C     INTERPOLACIJSKE FUNKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
C
CS    PRVIH 8 CVOROVA
CE    FIRST 4 NODES
C
      H(1)=0.125*RP*SP*TP
      H(2)=0.125*RM*SP*TP
      H(3)=0.125*RM*SM*TP
      H(4)=0.125*RP*SM*TP
      H(5)=0.125*RP*SP*TM
      H(6)=0.125*RM*SP*TM
      H(7)=0.125*RM*SM*TM
      H(8)=0.125*RP*SM*TM

      IF (IBR.EQ.0) THEN
      P(1,1)=0.125*SP*TP
      P(1,2)=-P(1,1)
      P(1,3)=-0.125*SM*TP
      P(1,4)=-P(1,3)
      P(1,5)=0.125*SP*TM
      P(1,6)=-P(1,5)
      P(1,7)=-0.125*SM*TM
      P(1,8)=-P(1,7)
C
      P(2,1)=0.125*RP*TP
      P(2,2)=0.125*RM*TP
      P(2,3)=-P(2,2)
      P(2,4)=-P(2,1)
      P(2,5)=0.125*RP*TM
      P(2,6)=0.125*RM*TM
      P(2,7)=-P(2,6)
      P(2,8)=-P(2,5)
C
      P(3,1)=0.125*RP*SP
      P(3,2)=0.125*RM*SP
      P(3,3)=0.125*RM*SM
      P(3,4)=0.125*RP*SM
      P(3,5)=-P(3,1)
      P(3,6)=-P(3,2)
      P(3,7)=-P(3,3)
      P(3,8)=-P(3,4)

       HV1=DOT(H,V1,NDIM)
       HV2=DOT(H,V2,NDIM)
       HV3=DOT(H,V3,NDIM)
  
      DO I=1,3
      DO J=1,3
       XJ(I,J)=0.D0
        DO  KK=1,NDIM
         XJ(I,J)=XJ(I,J)+P(I,KK)*CK(KK,J)
        ENDDO    
      ENDDO    
      ENDDO    
C       UX=XJ(1,1)*HV1+XJ(1,2)*HV2+XJ(1,3)*HV3
C       UY=XJ(2,1)*HV1+XJ(2,2)*HV2+XJ(2,3)*HV3
C       UZ=XJ(3,1)*HV1+XJ(3,2)*HV2+XJ(3,3)*HV3

C       VV=HV1**2+HV2**2+HV3**2
       CALL KOORD(CK,H,RST0)       
C       IF (DABS(VV).LT.1.D-10) RETURN
      ENDIF

      IF (IBR.EQ.1) CALL KOORD(CK,H,RP1)
      IF (IBR.EQ.2) CALL KOORD(CK,H,SP1)
      IF (IBR.EQ.3) CALL KOORD(CK,H,TP1)

      IF (IBR.EQ.4) CALL KOORD(CK,H,RM1)
      IF (IBR.EQ.5) CALL KOORD(CK,H,SM1)
      IF (IBR.EQ.6) CALL KOORD(CK,H,TM1)
     
20    CONTINUE       
       
      HX=RAZL(RP1,RM1)
      HY=RAZL(SP1,SM1)
      HZ=RAZL(TP1,TM1)

      CALL RAZL1(RP1,RST0,EX)
      CALL RAZL1(SP1,RST0,EY)
      CALL RAZL1(TP1,RST0,EZ)


       UX=EX(1)*HV1+EX(2)*HV2+EX(3)*HV3
       UY=EY(1)*HV1+EY(2)*HV2+EY(3)*HV3
       UZ=EZ(1)*HV1+EZ(2)*HV2+EZ(3)*HV3

       VV=HV1**2+HV2**2+HV3**2
       CALL KOORD(CK,H,RST0)       
       IF (DABS(VV).LT.TOL) RETURN


	AK=2.D0*AMI/GUSM
C	AK=2.D0/AMI

      ALFAX=UX*HX/AK
      ALFAY=UY*HY/AK
      ALFAZ=UZ*HZ/AK
      
      AX=0.D0
      AY=0.D0
      AZ=0.D0
      IF (DABS(ALFAX).GT.0.D0) AX=1.D0/DTANH(ALFAX)-1.D0/ALFAX
      IF (DABS(ALFAY).GT.0.D0) AY=1.D0/DTANH(ALFAY)-1.D0/ALFAY
      IF (DABS(ALFAZ).GT.0.D0) AZ=1.D0/DTANH(ALFAZ)-1.D0/ALFAZ
     
C       AX=HUGHES(ALFAX)
C       AY=HUGHES(ALFAY)
C       AZ=HUGHES(ALFAZ)

           
      AKK=((AX*UX*HX+AY*UY*HY+AZ*UZ*HZ)/2.D0)*1.D2



      DO I=1,NDIM
        PP=AKK*(HV1*PJ(1,I)+HV2*PJ(2,I)+HV3*PJ(3,I))/VV
        HH(I)=HH(I)+PP
      ENDDO
      
            
      END
C=======================================================================

C==========================================================================
       SUBROUTINE KOORD(CK,H,CORD)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CK(21,*),H(*),CORD(*)
       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C
CE Subroutine KOORD is used for finding coordinate of nodes
C

       DO I=1,3 
        CORD(I)=0.D0
       ENDDO

      DO J=1,3
       DO I=1,NDIM
         CORD(J)=CORD(J)+H(I)*CK(I,J)
        ENDDO
       ENDDO


      END
C=======================================================================
C==========================================================================
       FUNCTION RAZL(A,B)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION A(*),B(*)

       RAZL=DSQRT((A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2)

      END
C=======================================================================

C==========================================================================
       SUBROUTINE RAZL1(A,B,C)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION A(*),B(*),C(*)
C
CE Subroutine RAZL1 is used for finding difference between two vectors
C

       AINTEZ=DSQRT((A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2)
       C(1)=(A(1)-B(1))/AINTEZ
       C(2)=(A(2)-B(2))/AINTEZ
       C(3)=(A(3)-B(3))/AINTEZ

      END
C=======================================================================
C==========================================================================
       FUNCTION HUGHES(ALFA)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)


       IF (ALFA.LT.-1.D0.AND.DABS(ALFA).GT.0.D0) HUGHES=-1.D0-1.D0/ALFA
       IF (ALFA.GE.-1.D0.AND.ALFA.LE.1.D0) HUGHES=0.D0
       IF (ALFA.GT.1.D0.AND.DABS(ALFA).GT.0.D0) HUGHES=1.D0-1.D0/ALFA
      END
C=======================================================================
C==========================================================================
      SUBROUTINE CROSS(TT1,ID,NEL,CORD,VVREME,KKORAK,PRIT,NASLOV,NDIM,
     &IDPRIT,NETIP,PENALT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /IMEULZ/ PAKLST,PAKUNV,IDUZIN
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM

      DIMENSION TT1(*),ID(5,*),NEL(NDIM,*),CORD(3,*),PRIT(IDPRIT,*)
      DIMENSION TANG(3),ANORM(3),A(3),B(3),AN(3)
C
CE Subroutine CROSS is used for special postprocessing calculation for example:
CE Blood flow through the carotid bifurcation artery
CE This subroutine is used for calculation axial and secondary velocities
C
      
      CHARACTER*3 DOD
      CHARACTER*250 NASLOV
      CHARACTER *24 PAKLST,STRING


      PAKLST='BIFURCAT................'

      KNODE=66
      KELEM=50
      NUMSEC=3


C       NSTART=-(KNODE-1)
C       ISEL=-(KELEM-1)
       NSTART=1
       ISEL=1
	 IDEV=0

       DO 10 I=1,NUMSEC
        IDEV=IDEV+1
        II=30+IDEV
        DOD=CHAR(I+ICHAR('A')-1)
	  STRING=PAKLST(1:7)//'.'//DOD
        OPEN(II,FILE=STRING)
         
        NSTART=NSTART+KNODE*10
        NEND=NSTART+KNODE-1
        ISEL=ISEL+KELEM*10
        IEEL=ISEL+KELEM-1



      NODE1=NEL(1,ISEL)
      NODE2=NEL(2,ISEL)
      NODE3=NEL(3,ISEL)

      B(1)=CORD(1,NODE2)-CORD(1,NODE1)
      B(2)=CORD(2,NODE2)-CORD(2,NODE1)
      B(3)=CORD(3,NODE2)-CORD(3,NODE1)
      
      A(1)=CORD(1,NODE3)-CORD(1,NODE1)
      A(2)=CORD(2,NODE3)-CORD(2,NODE1)
      A(3)=CORD(3,NODE3)-CORD(3,NODE1)
      AN(1)=A(2)*B(3)-A(3)*B(2)
      AN(2)=-(A(1)*B(3)-A(3)*B(1))
      AN(3)=A(1)*B(2)-A(2)*B(1)
      BB=DSQRT(B(1)**2+B(2)**2+B(3)**2)
      ANN=DSQRT(AN(1)**2+AN(2)**2+AN(3)**2)

      TANG(1)=B(1)/BB
      TANG(2)=B(2)/BB
      TANG(3)=B(3)/BB

      ANORM(1)=AN(1)/ANN
      ANORM(2)=AN(2)/ANN
      ANORM(3)=AN(3)/ANN

        CALL FGRAF1(CORD,NSTART,NEND,II,ANORM)
        CALL TGRA22(NEL,4,ISEL,IEEL,1,II,NETIP,NDIM)
C        CALL STAGPP(TT1,ID,NASLOV,VVREME,KKORAK,1,NPT,II,1,NET,NEL,PRIT,
C     &CORD,NSTART,NEND,ISEL,IEEL,ANORM,TANG,NDIM,NETIP,IDPRIT,PENALT)
 10    CONTINUE        


      END
C=======================================================================
C=======================================================================
        CHARACTER *3 FUNCTION  INTSTR (INT)
         CHARACTER *1 J,D,S
         J='0'
         D='0'
         S='0'

         NULA=ICHAR('0')


         IS=INT/100
         NN=INT-IS*100
         ID=NN/10
         IJ=MOD(NN,10)
         J=CHAR(IJ+NULA)
         D=CHAR(ID+NULA)
         S=CHAR(IS+NULA)
         INTSTR=S//D//J
        END

C==========================================================================
C==========================================================================
      SUBROUTINE TANNOR(NEL,NDIM,ISEL,CORD,TANG,ANORM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION NEL(NDIM,*),TANG(*),ANORM(*),CORD(3,*)
      DIMENSION A(3),B(3),AN(3)

C
CE Subroutine TANNOR is used for finding tangent and normal vector 
CE for presribed cross-section
C

      NODE1=NEL(1,ISEL)
      NODE2=NEL(2,ISEL)
      NODE3=NEL(3,ISEL)

      B(1)=CORD(1,NODE2)-CORD(1,NODE1)
      B(2)=CORD(2,NODE2)-CORD(2,NODE1)
      B(3)=CORD(3,NODE2)-CORD(3,NODE1)
      
      A(1)=CORD(1,NODE3)-CORD(1,NODE1)
      A(2)=CORD(2,NODE3)-CORD(2,NODE1)
      A(3)=CORD(3,NODE3)-CORD(3,NODE1)
      AN(1)=A(2)*B(3)-A(3)*B(2)
      AN(2)=-(A(1)*B(3)-A(3)*B(1))
      AN(3)=A(1)*B(2)-A(2)*B(1)
      BB=DSQRT(B(1)**2+B(2)**2+B(3)**2)
      ANN=DSQRT(AN(1)**2+AN(2)**2+AN(3)**2)

      TANG(1)=B(1)/BB
      TANG(2)=B(2)/BB
      TANG(3)=B(3)/BB

      ANORM(1)=AN(1)/ANN
      ANORM(2)=AN(2)/ANN
      ANORM(3)=AN(3)/ANN
      END
C==========================================================================

C==========================================================================
      SUBROUTINE CROSSC(GNODE,ID,NEL,CORD,VVREME,KKORAK,PRIT,NASLOV,
     &NDIM,IDPRIT,NETIP,PENALT,NUMSEC,KNODE,KELEM,IBLOOD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /IMEULZ/ PAKLST,PAKUNV,IDUZIN
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM

      DIMENSION GNODE(2,5,*),ID(5,*),NEL(NDIM,*),CORD(3,*)
      DIMENSION PRIT(IDPRIT,*)
      DIMENSION TANG(3),ANORM(3)
      
C
CE Subroutine CROSSC is used for special postprocessing calculation for
CE example: Blood flow through the carotid bifurcation artery
C
      CHARACTER*4 DOD,STEP
      CHARACTER*250 NASLOV
      CHARACTER *12 STRING
      CHARACTER *3 INTSTR

C      KNODE=66
C      KELEM=50
C      NUMSEC=21


C       NSTART=-(KNODE-1)
C       ISEL=-(KELEM-1)
       NSTART=1-KNODE
       ISEL=1-KELEM
	 IDEV=0

C        STEP=INTSTR(KKORAK)
 


       DO 10 I=1,NUMSEC
        IDEV=IDEV+1
        II=65+IDEV
	  STEP='S'//INTSTR(KKORAK)
        DOD='C'//INTSTR(I)
        STRING=STEP//DOD//'.'//'UNV'
        OPEN(II,FILE=STRING)
         

        NSTART=NSTART+KNODE

         IF(I.EQ.11.AND.IBLOOD.EQ.1) THEN
            KNODE=71
         ELSEIF (I.GT.11.AND.IBLOOD.EQ.1) THEN
            KNODE=72
         ELSE
            KNODE=66
         ENDIF

        NEND=NSTART+KNODE-1

        ISEL=ISEL+KELEM
        IEEL=ISEL+KELEM-1
			   
C        CALL TANNOR(NEL,NDIM,ISEL,CORD,TANG,ANORM)
        CALL TANNOR(NEL,NDIM,IEEL,CORD,TANG,ANORM)

        CALL FGRAF1(CORD,NSTART,NEND,II,ANORM)
        CALL TGRA22(NEL,4,ISEL,IEEL,1,II,NETIP,NDIM)
      CALL STAGPP(GNODE,ID,NASLOV,VVREME,KKORAK,1,NPT,II,1,NET,NEL,PRIT,
     &CORD,NSTART,NEND,ISEL,IEEL,ANORM,TANG,NDIM,NETIP,IDPRIT,PENALT)


        CLOSE(II)
 10    CONTINUE        


      END
C=======================================================================
C==========================================================================
      SUBROUTINE CROSXY(GNODE,ID,NEL,CORD,VVREME,KKORAK,PRIT,NASLOV,
     &NDIM,IDPRIT,NETIP,PENALT,NPT,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /IMEULZ/ PAKLST,PAKUNV,IDUZIN
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM

      DIMENSION GNODE(2,5,*),ID(5,*),NEL(NDIM,*),CORD(3,*)
      DIMENSION PRIT(IDPRIT,*)
C
CE Subroutine CROSXY is used for special postprocessing calculation for
CE example: Blood flow through the carotid bifurcation artery
C
      
      CHARACTER*250 NASLOV
      CHARACTER *12 STRING
      CHARACTER *3 INTSTR


        STRING='PLANE'//INTSTR(KKORAK)//'.UNV'
        II=31
        OPEN(II,FILE=STRING)
        NSTART=1
        NEND=NPT
        ISEL=1
        IEEL=NET

        ZCORD=0.D0
        CALL TGRAXY(CORD,NSTART,NEND,II,ZCORD)
        CALL TGRA23(NEL,4,ISEL,IEEL,1,II,NETIP,NDIM,ZCORD,CORD)
        CALL STAGP2(GNODE,ID,NASLOV,VVREME,KKORAK,1,NPT,II,1,NET,NEL,
     &PRIT,CORD,NSTART,NEND,ISEL,IEEL,NDIM,NETIP,IDPRIT,PENALT,ZCORD)
        CLOSE(II)


      END
C=======================================================================
C==========================================================================
      SUBROUTINE VVNEW(NODE,ID,ANORM,TANG,VNEW,GNODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION GNODE(2,5,*),ID(5,*),ANORM(*),TANG(*),V(3),VNEW(*)

C
CE Subroutine VVNEW is used for finding axial anr radial veclocity vectors
CE when normal and tangent vector is prescribed
C
        
        DO J=1,3
           V(J)=GNODE(2,J,NODE)
        ENDDO
        VAXIAL=0.D0
        VRADI=0.D0
        DO J=1,3
          VAXIAL=VAXIAL+V(J)*ANORM(J)
          VRADI=VRADI+V(J)*TANG(J)
         ENDDO
	 VNEW(1)=VRADI
	 VNEW(2)=V(3)
	 VNEW(3)=VAXIAL

      END
C=======================================================================

C==========================================================================
      SUBROUTINE CORNEW(NODE,CORD,ANORM,CNEW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION CORD(3,*),P(2,2),E(2),CNEW(*),ANORM(*)
C
CE Subroutine CORNEW is used for calculation new coordinate of nodes
CE For example: Blood flow through the carotid bifurcation artery
C

        E(1)=CORD(1,NODE)
        E(2)=CORD(2,NODE)
         
        DC=ANORM(2)
        DS=ANORM(1)

    	  P(1,1)=DC
   	  P(1,2)=DS
   	  P(2,1)=-DS
C   	  P(1,2)=-DS
C   	  P(2,1)=DS
   	  P(2,2)=DC
             
         DO I=1,2
          CNEW(I)=0.D0
           DO J=1,2
            CNEW(I)=CNEW(I)+P(I,J)*E(J)
          ENDDO
         ENDDO
         
      END
C=======================================================================

C======================================================================
      SUBROUTINE FGRAF1(CORD,NSTART,NEND,II,ANORM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine FGRAF1 is used for printing data in output *.UNV file
CE For example: Blood flow through the carotid bifurcation artery
C

      DIMENSION CORD(3,*),ANORM(*),CNEW(2)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' TGRAFC'
      REWIND II
      IND=-1
      ITYP=15
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
 1100 FORMAT(I6)
      IT1=0
      IT2=0
      ICOL=8
      DO 10 NODE=NSTART,NEND
      CALL CORNEW(NODE,CORD,ANORM,CNEW)
C      WRITE(II,1000) NODE,IT1,IT2,ICOL,CNEW(2),CORD(3,NODE),0.D0
      WRITE(II,1000) NODE,IT1,IT2,ICOL,CNEW(1),CORD(3,NODE),0.D0
 1000 FORMAT(4I10,3E13.5)
   10 CONTINUE
      WRITE(II,1100) IND
      RETURN
      END
C=======================================================================
C======================================================================
      SUBROUTINE TGRAXY(CORD,NSTART,NEND,II,ZCORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
      COMMON /CDEBUG/ IDEBUG

C
CE Subroutine TGRAXY is used for printing data in output *.UNV file
CE For example: Blood flow through the carotid bifurcation artery
C

      REWIND II
      IND=-1
      ITYP=15
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
 1100 FORMAT(I6)
      IT1=0
      IT2=0
      ICOL=8
      DO 10 NODE=NSTART,NEND
      IF (DABS(CORD(3,NODE)-ZCORD).GT.1.D-8) GOTO 10
       WRITE(II,1000) NODE,IT1,IT2,ICOL,(CORD(J,NODE),J=1,3)
 1000 FORMAT(4I10,3E13.5)
   10 CONTINUE
      WRITE(II,1100) IND
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE STAGPP(GNODE,ID,NASLOV,VREME,KOR,NDT,NP,II,IND,NET,NEL,
     1PRITP,CORD,NSTART,NEND,ISEL,IEEL,ANORM,TANG,NDIM,NETIP,IDPRIT,
     &PENALT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS

C
CE Subroutine STAGPP is used for printout velocities,pressure,temperature,
CE Printing moving of fluid mesh
C
C      COMMON /TIPEL/ NP2DMX
C      COMMON /TIPELM/ NETIP
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C      COMMON /REPEA2/ NEL2D
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET1
C
C ......................................................................
C
      DIMENSION NEL(NDIM,*)
      DIMENSION PRITP(IDPRIT,*),CORD(3,*)
      CHARACTER*250 NASLOV
      DIMENSION GNODE(2,5,*),ID(5,*)
      DIMENSION ANORM(*),TANG(*),FSP(6),VNEW(3)

      IF(IDEBUG.GT.0) PRINT *, ' STAGP1'
      IMOTY=1
      NNDIM=4
      IF (NDIM.EQ.9) NNDIM=8
      IANTY=1
      IFAT1=1
      IFAT2=1					    
      FATY8=0.0D0
      IF(NDT.GT.1) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
      IDACH=3
      IF(IND.EQ.0) ISDTY=8
      IF(IND.EQ.1) ISDTY=11
      IF(IND.EQ.2) ISDTY=12
      IDATY=3
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=2
      IND1=-1
      ITYP=55
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,2
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) VREME
      DO 10 I=NSTART,NEND
        WRITE(II,5000) I
        CALL VVNEW(I,ID,ANORM,TANG,VNEW,GNODE)         
        WRITE(II,5200) (VNEW(JJ),JJ=1,3)
   10 CONTINUE
      WRITE(II,5100) IND1
      NNII=4
      IF (NDIM.EQ.9) NNII=8
C
C     FOR PRINTING PRESSURE
C
      WRITE(II,5100) IND1
      WRITE(II,5100) 57
      WRITE(II,5003) NASLOV
      WRITE(II,3005)
      WRITE(II,2000) KOR
      WRITE(II,3006) 1,1,4,20,2,6
      WRITE(II,3006) 1,1,KOR
      WRITE(II,5202) VREME
      IF (PENALT.GT.1.D0) THEN
      DO 30 I=ISEL,IEEL
      WRITE(II,3006) I,1,4,6
      DO 29 J=1,4
   29 WRITE(II,5200) PRITP(1,I),0.,0.,0.,0.,0.
   30 CONTINUE     
C   38 CONTINUE 
      CONTINUE 
C 40  CONTINUE     
      CONTINUE     

      ENDIF

      WRITE(II,5100) IND1


C FOR DISPLACEMENT OF SOLID
C
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
C     WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      WRITE(II,5000) IMOTY,IANTY,IDACH,8,IDATY,6
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 130 I=NSTART,NEND
        CALL VVNEW(I,ID,ANORM,TANG,VNEW,GNODE)         
            WRITE(II,5000) I
      DO IBROJ=1,6
       FSP(IBROJ)=0.D0
      ENDDO 
         FSP(3)=VNEW(3)
         WRITE(II,5200) (FSP(J),J=1,6)
  130 CONTINUE
      WRITE(II,5100) IND1



      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
 5202 FORMAT(1PE13.2)
C-----------------------------------------------------------------------
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 3005 FORMAT('PRITISCI U FLUIDU')
 3006 FORMAT(6I10)
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C 2001 FORMAT(' DATE'/
C     1       ' EMPTY'/
C     1       ' LOAD CASE',I10)
C-----------------------------------------------------------------------
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOCITIES AND ANGLE VELOCITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
C 7006 FORMAT(' NODAL TEMPERATURE')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C=======================================================================
      SUBROUTINE STAGP2(GNODE,ID,NASLOV,VREME,KOR,NDT,NP,II,IND,NET,NEL,
     1PRITP,CORD,NSTART,NEND,ISEL,IEEL,NDIM,NETIP,IDPRIT,PENALT,ZCORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS

C
CE Subroutine STAGP2 is used for special postprocessing printing for
CE example: Blood flow through the carotid bifurcation artery
C

C      COMMON /TIPEL/ NP2DMX
C      COMMON /TIPELM/ NETIP
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C      COMMON /REPEA2/ NEL2D
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET1
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS, VELOSITIES AND ACCELERATIONS
CE.      IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA U UNIVERZALNI FILE
C .
C ......................................................................
C
      DIMENSION NEL(NDIM,*),NN(4)
      DIMENSION PRITP(IDPRIT,*),CORD(3,*)
      CHARACTER*250 NASLOV
      DIMENSION GNODE(2,5,*),ID(5,*),FSP(6)

      IF(IDEBUG.GT.0) PRINT *, ' STAGP1'
      IMOTY=1
      NNDIM=4
      IF (NDIM.EQ.9) NNDIM=8
      IANTY=1
      IFAT1=1
      IFAT2=1					    
      FATY8=0.0D0
      IF(NDT.GT.1) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
      IDACH=3
      IF(IND.EQ.0) ISDTY=8
      IF(IND.EQ.1) ISDTY=11
      IF(IND.EQ.2) ISDTY=12
      IDATY=3
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=2
      IND1=-1
      ITYP=55
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,2
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) VREME
      DO 10 I=NSTART,NEND
      IF (DABS(CORD(3,I)-ZCORD).GT.1.D-8) GOTO 10
            WRITE(II,5000) I
         DO 20 J=1,NETIP
            FSP(J) = 0.0D0
            FSP(J)=GNODE(2,J,I)
   20    CONTINUE
         WRITE(II,5200) (FSP(J),J=1,2)
   10 CONTINUE
      WRITE(II,5100) IND1
      NNII=4
      IF (NDIM.EQ.9) NNII=8
C
C     ZA ISPISIVANJE PRITISAKA
C
      WRITE(II,5100) IND1
      WRITE(II,5100) 57
      WRITE(II,5003) NASLOV
      WRITE(II,3005)
      WRITE(II,2000) KOR
      WRITE(II,3006) 1,1,4,20,2,6
      WRITE(II,3006) 1,1,KOR
      WRITE(II,5202) VREME
      IF (PENALT.GT.1.D0) THEN
      DO 30 I=ISEL,IEEL
      
      CALL ZPLANE(NEL,NDIM,CORD,ZCORD,NN,K,I)
	IF(K.LT.4) GOTO 30
     

      WRITE(II,3006) I,1,4,6
      DO 29 J=1,4
   29 WRITE(II,5200) PRITP(1,I),0.,0.,0.,0.,0.
   30 CONTINUE     
C   38 CONTINUE 
      CONTINUE 
C 40  CONTINUE     
      CONTINUE     

      ENDIF

      WRITE(II,5100) IND1

C ZA ISPISIVANJE POMERANJA SOLIDA
C
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
C     WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      WRITE(II,5000) IMOTY,IANTY,IDACH,8,IDATY,6
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 130 I=NSTART,NEND
      IF (DABS(CORD(3,I)-ZCORD).GT.1.D-8) GOTO 130
            WRITE(II,5000) I
      DO IBROJ=1,6
       FSP(IBROJ)=0.D0
      ENDDO 
         DO 120 J=1,NETIP
            FSP(J)=GNODE(2,J,I)
 120    CONTINUE
         WRITE(II,5200) (FSP(J),J=1,6)
  130 CONTINUE
      WRITE(II,5100) IND1

      RETURN


 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
 5202 FORMAT(1PE13.2)
C-----------------------------------------------------------------------
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 3005 FORMAT('PRITISCI U FLUIDU')
 3006 FORMAT(6I10)
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C 2001 FORMAT(' DATE'/
C     1       ' EMPTY'/
C     1       ' LOAD CASE',I10)
C-----------------------------------------------------------------------
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOCITIES AND ANGLE VELOCITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
C 7006 FORMAT(' NODAL TEMPERATURE')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C=======================================================================
      SUBROUTINE TGRA22(NEL,NBR2,ISTART,IEND,NGE,II,NETIP,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
C      COMMON /TIPELM/ NETIP

C
CE Subroutine TGRA22 is used for special postprocessing printing for
CE example: Blood flow through the carotid bifurcation artery
C

C      COMMON /TIPEL/ NP2DMX

C      CHARACTER*250 ACOZ

      COMMON /SRPSKI/ ISRPS
      DIMENSION NEL(NDIM,*)
C
C     E L E M E N T I   2/D
C
C      IF(ISRPS.EQ.0.AND.(NBR2.NE.4.AND.NBR2.NE.8))
C     1WRITE(IIZLAZ,2200) NGE
C      IF(ISRPS.EQ.1.AND.(NBR2.NE.4.AND.NBR2.NE.8))
C     1WRITE(IIZLAZ,6200) NGE


      NBR=NBR2
      IF(NBR2.LT.8) NBR2=4
      IF(NBR2.EQ.9) NBR2=8
C     GRAFICKI OPIS RAVANSKOG ELEMENTA: SA 4 CVORA = 27, SA 8 CVOROVA = 28
      ITYPE=27 
      IF(NBR2.EQ.8) ITYPE=28
C     VRSTA 2/D ELEMENTA: 
      IE1=44
      IF(NBR2.EQ.8) IE1=45
C     TABELA FIZICKIH OSOBINA
      IE2=1
C     TABELA MATERIJALA
      IE3=1
C     BOJA  
      ICOL=8
      IND=-1
      ITYP=71
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 10 I=ISTART,IEND
         WRITE(II,1000) I,ITYPE,IE1,IE2,IE3,ICOL,NBR2
         IF(NBR2.EQ.4) THEN
            WRITE(II,1000) (NEL(J,I),J=1,4)
         ELSE
            WRITE(II,1000) (NEL(J,I),NEL(J+4,I),J=1,4)
         ENDIF
   10 CONTINUE
      WRITE(II,1100) IND
      NBR2=NBR
      RETURN
C
 1100 FORMAT(I6)
 1000 FORMAT(8I10)
C-----------------------------------------------------------------------
C 2200 FORMAT(//' PROGRAM ZA GRAFICKO PRIKAZIVANJE REZULTATA "IDEAS"'/
C     1' ZAHTEVA 2/D ELEMENT SA 4 ILI 8 CVOROVA U GRUPI ELEMENATA NGE ='
C     1,I5)
C-----------------------------------------------------------------------
C 6200 FORMAT(//' GRAPHIC PACKAGE   "IDEAS"'/
C     1' PERMITS ONLY 2/D ELEMENTS WITH 4 OR 8 NODES PER ELEMENT IN',
C     1' GROUP   NGE =',I5)
C-----------------------------------------------------------------------
C=======================================================================
C=======================================================================
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ZPLANE(NEL,NDIM,CORD,ZCORD,NN,K,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NEL(NDIM,*),CORD(3,*),NN(4)
C
CE Subroutine ZPLANE is used for special postprocessing calculation for
CE example: Blood flow through the carotid bifurcation artery
C
        
         K=0
         DO JJ=1,NDIM
		NODE=NEL(JJ,I)
          IF (DABS(CORD(3,NODE)-ZCORD).LT.1.D-5) THEN
           K=K+1
           NN(K)=NODE
          ENDIF
         ENDDO   

		IF (K.EQ.4) THEN
    	     II=NN(3)
           NN(3)=NN(4)
           NN(4)=II
          ENDIF

      END
C=======================================================================
C=======================================================================
      SUBROUTINE TGRA23(NEL,NBR2,ISTART,IEND,NGE,II,NETIP,NDIM,ZCORD,
     &CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /SRPSKI/ ISRPS
      DIMENSION NEL(NDIM,*),CORD(3,*),NN(4)

C
CE Subroutine TGRA23 is used for special postprocessing printing for
CE example: Blood flow through the carotid bifurcation artery
C

      NBR=NBR2
      IF(NBR2.LT.8) NBR2=4
      IF(NBR2.EQ.9) NBR2=8
C     GRAFICKI OPIS RAVANSKOG ELEMENTA: SA 4 CVORA = 27, SA 8 CVOROVA = 28
      ITYPE=27 
      IF(NBR2.EQ.8) ITYPE=28
C     VRSTA 2/D ELEMENTA: 
      IE1=44
      IF(NBR2.EQ.8) IE1=45
C     TABELA FIZICKIH OSOBINA
      IE2=1
C     TABELA MATERIJALA
      IE3=1
C     BOJA  
      ICOL=8
      IND=-1
      ITYP=71
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 10 I=ISTART,IEND

      CALL ZPLANE(NEL,NDIM,CORD,ZCORD,NN,K,I)
	 IF(K.LT.4) GOTO 10

         WRITE(II,1000) I,ITYPE,IE1,IE2,IE3,ICOL,NBR2
         IF(NBR2.EQ.4) THEN
            WRITE(II,1000) (NN(J),J=1,4)
         ELSE
            WRITE(II,1000) (NEL(J,I),NEL(J+4,I),J=1,4)
         ENDIF
   10 CONTINUE
      WRITE(II,1100) IND
      NBR2=NBR
      RETURN
C
 1100 FORMAT(I6)
 1000 FORMAT(8I10)
C-----------------------------------------------------------------------
C 2200 FORMAT(//' PROGRAM ZA GRAFICKO PRIKAZIVANJE REZULTATA "IDEAS"'/
C     1' ZAHTEVA 2/D ELEMENT SA 4 ILI 8 CVOROVA U GRUPI ELEMENATA NGE ='
C     1,I5)
C-----------------------------------------------------------------------
C 6200 FORMAT(//' GRAPHIC PACKAGE   "IDEAS"'/
C     1' PERMITS ONLY 2/D ELEMENTS WITH 4 OR 8 NODES PER ELEMENT IN',
C     1' GROUP   NGE =',I5)
C-----------------------------------------------------------------------
C=======================================================================
C=======================================================================
      END
C=======================================================================
C==========================================================================
      SUBROUTINE EXPAND(CCORD,CORD,GTIME,NPT,ZADVRE,ZADVR1,NUMZAD,
     &NZAD,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       DIMENSION CCORD(3,*),CORD(3,*),ZADVRE(*),ZADVR1(*)
       DIMENSION NZAD(3,*)

C
CE Subroutine EXAPND is used for special calculation for
CE example: Flow through expandable tube
C

       PI=4.D0*DATAN(1.D0)
       C=0.25D0
       FI=(1.D0+C)**(1.D0/3.D0)
       P=3.D0
       AN=2*PI/P
       AK=(FI-1.D0)/(FI+1.D0)
       FT=1.D0+AK*DSIN(AN*GTIME)
       BB=3.7125D-3
       AA=250.D-6
	 XX=0.D0

       DO I=1,NUMZAD
C        IF (I.LE.10) THEN
C          ZADVR1(I)=ZADVRE(I)*DCOS(AN*GTIME)
C        ELSE IF(I.GT.10) THEN
          IP=NZAD(2,I)
          NODE=NZAD(1,I)
          IF (IP.EQ.2) XX=CORD(2,NODE)
          IF (IP.EQ.1) XX=CORD(1,NODE)
C          IF (IP.EQ.2) XX=BB
C          IF (IP.EQ.1) XX=AA
C          ZADVR1(I)=AN*AK*XX*DCOS(AN*GTIME)
           ZADVR1(I)=AN*AK*XX*DSIN(AN*GTIME)
          IF (IP.EQ.3) ZADVR1(I)=ZADVRE(I)
C          WRITE(IIZLAZ,100)NODE,IP,XX,ZADVR1(I),GTIME
C        ENDIF
       ENDDO
      
        DO I=1,NPT
          CCORD(1,I)=CORD(1,I)*FT
          CCORD(2,I)=CORD(2,I)*FT
        ENDDO

C 100   FORMAT(2I5,2(1PE13.5),F10.3)
      END
C==========================================================================
C==========================================================================
      SUBROUTINE COLLEA(CORD,CCORD,TT1,NPT,ID,IIZLAZ,VVREME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION CORD(3,*),CCORD(3,*),TT1(*),ID(5,*)
      DIMENSION XXLM(3)
C
CE Subroutine COLLEA is used for special calculation for some example
C

      XXLM(1)=0.D0
      XXLM(2)=1.9305E-3
      XXLM(3)=3.7125E-3

      TOL=1.D-7

      WRITE(IIZLAZ,*)'TIME(S)= ',VVREME
      DO II=1,3
      WRITE(IIZLAZ,*)'X/XLM= ',II
      DO NODE=1,NPT
       XX=CCORD(2,NODE) 
       X=CORD(2,NODE) 
       IF (DABS(X-XXLM(II)).LT.TOL) THEN
        RR=CORD(1,NODE)
        JAXI=ID(2,NODE)
        JRAD=ID(1,NODE)
        AXI=0.D0
        RAD=0.D0
        IF (JAXI.NE.0) AXI=TT1(JAXI)
        IF (JRAD.NE.0) RAD=TT1(JRAD)
        WRITE(IIZLAZ,100) NODE,AXI,RR ,RAD,RR
       ENDIF
      ENDDO
      ENDDO

 100  FORMAT(I5,4(1PE13.5))
       

      END
C==========================================================================
C==========================================================================
      SUBROUTINE INITIA(CORD,CCORD,NEL,TT1,INDEL,SPSIL,NET,NPT,NEQF,
     &NDIM,KKORAK,VVREME,NETIP,GNODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CCORD(3,*),CORD(3,*),TT1(*),SPSIL(NETIP,*),GNODE(2,5,*)
      DIMENSION NEL(NDIM,*),INDEL(*)
     
C
CE Subroutine INITIA is used for initialisation all global variable
C

      VVREME=0.D0
      KKORAK=0
       
C INITIALISATION:
       DO I=1,NPT
        DO J=1,NETIP
         SPSIL(J,I)=0.D0
        ENDDO
        INDEL(I)=0
        DO J=1,5
         GNODE(1,J,I)=0.D0
         GNODE(2,J,I)=0.D0
        ENDDO
        DO J=1,3
         CCORD(J,I)=CORD(J,I)
        ENDDO
       ENDDO

       DO NBREL=1,NET
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         INDEL(NODE)=INDEL(NODE)+1
        ENDDO        
       ENDDO

c TT1-vector of unknowns values
      CALL CLEAR(TT1,NEQF)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE RETPRI(GNODE,INDEL,NPT,TT1,ID,NDIM,PENALT,PRIT,
     &IDPRIT,NEL,NET,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION INDEL(*),ID(5,*),NEL(NDIM,*)
      DIMENSION GNODE(2,5,*),TT1(*),N(21),PRIT(IDPRIT,*),P(21)
C      
CE Subroutine RETPRI is used for calculation pressure at all nodes per element
CE with prescribed values
CS RACUNANJE PRITISAKA U SVIM CVOROVIMA AKO SE KORISTE ZADATE VREDNOSTI
C 23.02.2011.
      DO I=1,NPT
       INDEL(I)=0
      ENDDO
       DO NBREL=1,NET
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         INDEL(NODE)=INDEL(NODE)+1
        ENDDO        
       ENDDO
      DO NBREL=1,NET
       DO I=1,NDIM
        N(I)=NEL(I,NBREL)       
        P(I)=GNODE(2,4,N(I))
       ENDDO
      IF(NDIM.GE.8.AND.NETIP.EQ.2) THEN
        P(5)=0.5D0*(P(1)+P(2))
        P(6)=0.5D0*(P(2)+P(3))
        P(7)=0.5D0*(P(3)+P(4))
        P(8)=0.5D0*(P(4)+P(1))
        IF(NDIM.EQ.9) P(9)=0.25D0*(P(1)+P(2)+P(3)+P(4))
            ELSEIF(NETIP.EQ.3.AND.NDIM.GE.20)THEN
            P(9)=0.5D0*(P(1)+P(2))
            P(10)=0.5D0*(P(2)+P(3))
            P(11)=0.5D0*(P(3)+P(4))
            P(12)=0.5D0*(P(4)+P(1))
            P(13)=0.5D0*(P(5)+P(6))
            P(14)=0.5D0*(P(6)+P(7))
            P(15)=0.5D0*(P(7)+P(8))
            P(16)=0.5D0*(P(8)+P(5))
            P(17)=0.5D0*(P(1)+P(5))
            P(18)=0.5D0*(P(2)+P(6))
            P(19)=0.5D0*(P(3)+P(7))
            P(20)=0.5D0*(P(4)+P(8))
       IF(NDIM.EQ.21)   
     1  P(21)=0.125D0*(P(1)+P(2)+P(3)+P(4)+P(5)+P(6)+P(7)+P(8))
       ENDIF
      IF(NETIP.EQ.2) IBRCVR=4
      IF(NETIP.EQ.3) IBRCVR=8
       DO I=IBRCVR+1,NDIM
          N(I)=NEL(I,NBREL)
          GNODE(2,4,N(I))=P(I)
       ENDDO
      ENDDO
      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE RETPRIP(GNODE,INDEL,NPT,TT1,ID,NDIM,PENALT,PRIT,
     &IDPRIT,NEL,NET,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION INDEL(*),ID(5,*),NEL(NDIM,*)
      DIMENSION GNODE(2,5,*),TT1(*),N(21),PRIT(IDPRIT,*),P(21)
C
CE Subroutine RETPRIP is used for calculation pressure at element level
C
CS RACUNANJE PRITISAKA AKO JE UKLJUCENA PENALTY METODA
CE CALCULATION OF PRESSURE IF PENALTY METHOD IS ENABLED
C 23.02.2011.

      DO I=1,NPT
       GNODE(2,4,I)=0.D0
       INDEL(I)=0
      ENDDO
       DO NBREL=1,NET
        DO I=1,NDIM
         NODE=NEL(I,NBREL)
         INDEL(NODE)=INDEL(NODE)+1
        ENDDO        
       ENDDO
      DO NBREL=1,NET
       DO I=1,NDIM
        P(I)=0.D0
        N(I)=NEL(I,NBREL)       
        JJ=ID(4,N(I))
        IF (JJ.NE.0) P(I)=TT1(JJ)
       ENDDO
       IF (PENALT.GT.1.D0) THEN
        DO I=1,NDIM
         P(I)=PRIT(1,NBREL)
        ENDDO
       ENDIF
        DO I=1,NDIM
         NODE=N(I)
          GNODE(2,4,NODE)=GNODE(2,4,NODE)+P(I)/(INDEL(NODE)*1.D0)
C                    WRITE(3,*)'INDEL,NODE ',NODE,INDEL(NODE)
       ENDDO
      ENDDO
      END
C==========================================================================
C==========================================================================
      SUBROUTINE FILLN(GNODE,TT1,ID,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(2,5,*),TT1(*),ID(5,*)
C
CE Subroutine FILLN is used for re-writing values at the end of time step
C
       DO NODE=1,NPT 
        DO I=1,5
         JJ=ID(I,NODE)
C         GNODE(I,NODE)=0.D0
         IF (JJ.NE.0) GNODE(2,I,NODE)=TT1(JJ)
        ENDDO
       ENDDO
      END
C==========================================================================
C==========================================================================
      SUBROUTINE FILZAD(GNODE,ID,NPT,NZAD,ZADVRE,NUMZAD,ITFMAX,NTABFT,
     &IIZLAZ,TABF,VVREME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,5,*),ID(5,*),NZAD(3,*),ZADVRE(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
       
CS    PROGRAM
CS           ZA UBACIVANJE ZADATIH VREDNOSTI U GNODE(2,I,NODE)
CE    PROGRAM
CE           FOR INSERT OF PRESCRIBED VALIEAS IN GNODE(2,I,NODE)
CE    LAST UPDATE: 03. 02. 2011
      IF(NUMZAD.GT.0)THEN
      
      DO 425 I=1,NUMZAD
       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I),
     &NTABFT,IIZLAZ)
C====================================
C FOR AKIRA'S EXAMPLE
C        FK1=VSR
C        IF(NZAD(1,I).GT.(0.75*NPT)) FK1=VSR1
C====================================
        GNODE(2,NZAD(2,I),NZAD(1,I))=ZADVRE(I)*FK1
  425  CONTINUE
  
      ENDIF
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ZADNOD(GNODE,ZADVRE,NZAD,TABF,VVREME,ITFMAX,NTABFT,
     &IIZLAZ,NUMZAD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION GNODE(2,5,*),ZADVRE(*),TABF(2,NTABFT,*)
       DIMENSION ITFMAX(*),NZAD(3,*)

      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE
      COMMON /AKIRA1/ VOLT1,PBALL1,VSR1,PTUBE1
C
CE Subroutine ZADNOD is used for inclusion prescribed values
C
      DO 425 I=1,NUMZAD
       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I),
     &NTABFT,IIZLAZ)
C====================================
C FOR AKIRA'S EXAMPLE
C        FK1=VSR
C        IF(NZAD(1,I).GT.(0.75*NPT)) FK1=VSR1
C====================================
        GNODE(2,NZAD(2,I),NZAD(1,I))=ZADVRE(I)*FK1
  425  CONTINUE 
      END
C==========================================================================
C==========================================================================
      SUBROUTINE SHEARS(AN,TT21,VTANG,NDIM) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION TT21(*),VTANG(3,*),AN(*),V(3),P(3),T(3)
 
C
CE Subroutine SHEARS is used for calculation shear stresses
C
       DO I=1,NDIM
         V(1)=TT21(I)          
         V(2)=TT21(I+NDIM)          
         V(3)=TT21(I+2*NDIM)
         CALL VECMUL(AN,V,P)
         CALL VECMUL(P,AN,T)
         TT=SQRT(T(1)**2+T(2)**2+T(3)**2) 
         IF (DABS(TT).LT.1.D-9) GOTO 30
         T(1)=T(1)/TT
         T(2)=T(2)/TT
         T(3)=T(3)/TT
  30     VV=DOT(T,V,3)
         VTANG(1,I)=VV*T(1)
         VTANG(2,I)=VV*T(2)
         VTANG(3,I)=VV*T(3)
       ENDDO

      END
C==========================================================================
C==========================================================================
      SUBROUTINE VECMUL(A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(*),C(*)

C
CE Subroutine VECMUL is used for calculation multiple of vectors
C

        C(1)=A(2)*B(3)-A(3)*B(2)
        C(2)=A(3)*B(1)-A(1)*B(3)
        C(3)=A(1)*B(2)-A(2)*B(1)
       

      END
C==========================================================================
C==========================================================================
      SUBROUTINE SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,NUMZAD,
     &IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,ZADVRE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NEL(NDIM,*),ID(5,*),INDEL(*),NZAD(3,*)
      DIMENSION CK(21,*),PJ(3,*),H(*),HP(*),TT21(*),PRES(3,*),ZADVRE(*)
	DIMENSION NID(21),N(21),NWALL(6),ITR(6,4)
      DIMENSION XG(*),WGT(*),NREF(*),SHEAR(3)

C
CE Subroutine SSTRES is used for calculation shear stresses
C
      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT

C      WRITE(IIZLAZ,*)'ELEMENT= ',NBREL
      
	DO I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.0.AND.ID(2,N(I)).EQ.0.AND.ID(3,N(I)).EQ.0)
     &NID(I)=1
      ENDDO

      ITR(1,1)=N(8)
      ITR(1,2)=N(4)
      ITR(1,3)=N(5)
      ITR(1,4)=N(1)

      ITR(2,1)=N(7)
      ITR(2,2)=N(3)
      ITR(2,3)=N(6)
      ITR(2,4)=N(2)

      ITR(3,1)=N(6)
      ITR(3,2)=N(2)
      ITR(3,3)=N(5)
      ITR(3,4)=N(1)

      ITR(4,1)=N(7)
      ITR(4,2)=N(3)
      ITR(4,3)=N(8)
      ITR(4,4)=N(4)

      ITR(5,1)=N(3)
      ITR(5,2)=N(2)
      ITR(5,3)=N(4)
      ITR(5,4)=N(1)

      ITR(6,1)=N(7)
      ITR(6,2)=N(6)
      ITR(6,3)=N(8)
      ITR(6,4)=N(5)

      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300
C
CS  PETLJA PO GAUSOVIM TACKAMA
CE  GAUSS POINTS LOOP
C
      IF(NPOV.GT.2) GO TO 45
      KFIX=1
      NGXP=NGAUSY
      NGYP=NGAUSZ
      R=1.
      IF(NPOV.EQ.2) R=-1.
      GO TO 60
C
   45 IF(NPOV.GT.4) GO TO 50
      KFIX=2
      NGXP=NGAUSX
      NGYP=NGAUSZ
      S=1.
      IF(NPOV.EQ.4) S=-1.
      GO TO 60
C
   50 NGXP=NGAUSX
      NGYP=NGAUSY
      KFIX=3
      T=1.
      IF(NPOV.EQ.6) T=-1.
C
   60 KK=0
      DO 210 NGX=1,NGXP
      JR=NREF(NGXP) + NGX
      XX = XG(JR)
      WX=WGT(JR)
      DO 210 NGY=1,NGYP
      JS=NREF(NGYP) + NGY
      YY = XG(JS)
      WY=WGT(JS)
      GO TO (71,72,73),KFIX
   71 S=XX
      T=YY
      GO TO 75
   72 R=XX
      T=YY
      GO TO 75
   73 R=XX
      S=YY
   75 WT=WX*WY

      KK=KK+1
      NODE=ITR(NPOV,KK)
      DO JJ=1,NUMZAD
       IF(NODE.EQ.NZAD(1,JJ).AND.DABS(ZADVRE(JJ)).LT.1.D-12) GOTO 150
      ENDDO
C     GOTO 300
 150   CALL JACT(R,S,T,DET1,CK,KFIX,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT)
C      WRITE(3,*)'SHEAR STRESS= ',PRES(1,NODE),PRES(2,NODE),PRES(3,NODE)

      PRES(1,NODE)=PRES(1,NODE)-AMI*SHEAR(1)/INDEL(NODE)
      PRES(2,NODE)=PRES(2,NODE)-AMI*SHEAR(2)/INDEL(NODE)
      PRES(3,NODE)=PRES(3,NODE)-AMI*SHEAR(3)/INDEL(NODE)

C      INDEL(NODE)=INDEL(NODE)+1

  210 CONTINUE    
  300 CONTINUE    
      END
C==========================================================================
C==========================================================================
      SUBROUTINE PRTIME(NNPER,KKORAK,VVREME,ISRPS,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine PRTIME is used for printing output data 
CE for number of periods, time steps
C
      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,100) NNPER         
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,101) NNPER         

 101   FORMAT (11X,'PERIOD NUMBER ......................... NNPER =',I5)
 100   FORMAT (11X,'PERIOD BROJ ........................... NNPER =',I5)


      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,200) KKORAK
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,201) KKORAK        

 201   FORMAT (11X,'STEP NUMBER ........................... NSTEP =',I5)
 200   FORMAT (11X,'KORAK BROJ ............................ NSTEP =',I5)



      IF (ISRPS.EQ.0)
     *WRITE(IIZLAZ,300) VVREME
      IF (ISRPS.EQ.1)
     *WRITE(IIZLAZ,301) VVREME        

 301   FORMAT (11X,'TIME ................................ TIME =',E10.3)
 300   FORMAT (11X,'VREME ............................... TIME =',E10.3)


      END
C=======================================================================
C=======================================================================
      SUBROUTINE MAXMEM(LMAX,NTOT,IZLAZ,KSTEP,NWK)
C
      IF (KSTEP.NE.1) RETURN

      WRITE(IZLAZ,100) LMAX,NTOT,NWK

 100  FORMAT(//6X,'D A T A   F O R   M E M O R Y   O C C U P A T I O N'
     1/6X,51('-')///
     18X,'MAXIMUM OCCUPATION MEMORY FOR SOLVING PROBLEM.. LMAX= ',I10//
     18X,'MAXIMUM DEFINITION MEMORY ..................... NTOT= ',I10//
     18X,'MAX. NUM. OF TERMS IN STIFFNESS MATRIX ........  NWK= ',I10//)

       WRITE(IZLAZ,200) LMAX*4./1.D6
       WRITE(IZLAZ,300) NTOT*4./1.D6
       WRITE(*,200) LMAX*4./1.D6
       WRITE(*,300) NTOT*4./1.D6
 200  FORMAT(' YOU NEED FOR SOLVING THIS PROBLEM  ',
     *F5.2,'  MB RAM MEMORY')
 300  FORMAT(' THIS VERSION OF THE PROGRAM PAK-F HAS MAXIMUM  ',
     *F5.2,'  MB RAM MEMORY')

      END
