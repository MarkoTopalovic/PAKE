C======================================================================
CE        ALL SUBROUTINES IN THIS FILE ARE PROGRAMD BY DRAKCE
C======================================================================
      SUBROUTINE ISPAKUJ(IROW,ISK,MAXA,NWK,JEDN)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      include 'paka.inc'
      
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /CDEBUG/ IDEBUG
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      INTEGER*8 MAXA,NWK,I8,NZERO8
      INTEGER*1 ISK(1+NWK/7)
      DIMENSION MAXA(*),LM(100),IROW(*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUJ'
      NBLOCK=1
      IF(NBLOCK.EQ.1) THEN
CE       WITHOUT BLOCKS
CS       BEZ BLOKOVA
CS     UNULJIVANJE MATRICE
CE     CLEARING MATRIX ISK 
      DO 10 I8=1,1+NWK/7
         ISK(I8)=0
  10  CONTINUE
      DO 11 I=1, JEDN
         n256=MAXA(I)-((MAXA(I)-1)/7)*7
         if(n256.LE.0) THEN
          WRITE(3,*) 'N256,I,MAXA(I),((MAXA(I)-1)/7)*7',
     1                N256,I,MAXA(I),((MAXA(I)-1)/7)*7
          stop 'ISPAKUJ - n256=MAXA(I)-((MAXA(I)-1)/7)*7'
          ENDIF
         CALL POST1(ISK(1+(MAXA(I)-1)/7),MAXA(I)-((MAXA(I)-1)/7)*7)
  11  CONTINUE
      REWIND(IDRAKCE)
      DO 20 J=1,NELUK
c         READ(IDRAKCE,1000) ND,(LM(I),I=1,ND)
         READ(IDRAKCE) ND,(LM(I),I=1,ND)
 1000    FORMAT(101I5)
         LMNQ=1
         LLREC=1
         LR=1
C
C         CALL IIISWR(LM,ND,' LM ')  
C
         CALL ISPAKUA(IROW,ISK,MAXA,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC),
     &               NWK,JEDN)
  20  CONTINUE
      CALL ICOUNT(ISK,NWK,JEDN,NZERO8,NEED1,NEED2,NEED3,NNZERO)
C PROVERI TREBA DA SE IZBACI      NNZERO=NWK-NZERO8
      ELSE
CE       WITH BLOCKS
CS       SA BLOKOVIMA
         STOP 'NERADI SA BLOKOVIMA'
C         WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
      ENDIF
      RETURN
      END
C====================
      SUBROUTINE ISPAKUA(IROW,ISK,MAXA,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC,NWK,JEDN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
      include 'paka.inc'
      
!      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      INTEGER*8 MAXA,NWK,MI,KK,KKNWK
      DIMENSION LM(*),MAXA(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(*)
      INTEGER*1 ISK(1+NWK/7)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
C      IF(IILS.NE.-1) CALL KORIGC(A(LCORUL),A(LCMPC),MMP,NP,NEZAV)
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        IF(INDD.EQ.1)GO TO 9
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(KB0)
      MNQ1=MNQ(KB0+1)-1
      MXMN=MAXA(MNQ(KB0))-1
      LDB=MAXA(MNQ(KB0+1))-MAXA(MNQ(KB0))
      CALL RBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   IF(ICCGG.EQ.2) THEN
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ELSE
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ENDIF
C-----------------------------------------------
   11 NDI=0
99999 FORMAT(5I5)
      DO 200 I=1,ND
      II=LM(I)
C
C
      IF(II.LT.0)THEN
        IIP=-II
        ICM=MPC(1,IIP)
        DO 320 L=1,NEZAV
          II=MPC(L+1,IIP)
          IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 320
          CMI=CMPC(ICM,L)
          MI=MAXA(II)-MXMN
          KS=I
          DO 310 J=1,ND
            JJ=LM(J)
            IF(JJ)303,310,307
  303       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           CMJ=CMPC(JCM,K)
                IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  KKNWK=KK+NWK
                  CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                  ISK(KK)=1
C                  ISK(KKNWK)=1
                ELSE          
                  KK=MI+IJ
                  CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                  ISK(KK)=1
                ENDIF
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         IF(ICCGG.EQ.2) THEN
                KK=MI-IJ
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
                CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                ISK(KK)=1
C                ISK(KKNWK)=1
              ELSE
                KK=MI+IJ
                KSS=KS
                IF(J.GE.I)KSS=J+NDI
                CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                ISK(KK)=1
              ENDIF
  310         KS=KS+ND-J
  320   CONTINUE
        GO TO 200
      ENDIF
      IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 200
      MI=MAXA(II)-MXMN
      KS=I
      DO 220 J=1,ND
      JJ=LM(J)
      IF(JJ)420,220,110
C
C
  420       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  KKNWK=KK+NWK
                  CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                  ISK(KK)=1
C                  ISK(KKNWK)=1
                ELSE
                  KK=MI+IJ
                  CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C                  ISK(KK)=1
                ENDIF
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 IF(ICCGG.EQ.2) THEN
         KK=MI-IJ
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
         CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C         ISK(KK)=1
C         ISK(KKNWK)=1
      ELSE
         KK=MI+IJ
         KSS=KS
         IF(J.GE.I)KSS=J+NDI
         CALL POST1(ISK(1+(KK-1)/7),KK-((KK-1)/7)*7)
C         ISK(KK)=1
      ENDIF
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
   20 CONTINUE
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C===========================================================
      SUBROUTINE ICOUNT(ISK,NWK,NN,NZERO,NEED1,NEED2,NEED3,NNZERO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      INTEGER*8 NWK,I8,NZERO
      INTEGER*1 ISK(1+NWK/7)
      NZERO=0
      DO 10 I8=1,NWK
	  CALL UZMI1(ISK(1+(I8-1)/7),I8-((I8-1)/7)*7,JEDAN)
        IF(JEDAN.EQ.0) NZERO=NZERO+1
10    CONTINUE
      NNZERO=NWK-NZERO
      A1=FLOAT(NZERO)/FLOAT(NWK)*100.
      A2=FLOAT(NNZERO)/FLOAT(NWK)*100.
      NEED1=NNZERO*2+NN*5
      NEED2=NNZERO
      NEED3=NEED1+NEED2
      IF(ISRPS.EQ.0)
     1WRITE(*,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
      IF(ISRPS.EQ.1)
     1WRITE(*,30)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,30)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
20    FORMAT(//' DIMENZIJA MATRICE    : ',I12/
     *         ' UKUPAN BROJ CLANOVA  : ',I12/
     *         ' BROJ CLANOVA = 0     : ',I12,'  ILI ',F6.0,' %'/
     *         ' BROJ NE NULA CLANOVA : ',I12,'  ILI ',F6.0,' %'/
     *         ' POTREBNA VELICINA GLAVNOG VEKTORA  : ',I12/
     *         ' POTREBNA VELICINA POMOCNOG VEKTORA : ',I12/
     *' UKUPNA POTREBNA VELICINA GLAVNOG I POMOCNOG VEKTORA : ',I9//)
30    FORMAT(//' DIMENSION OF SYSTEM MATRIX: ',I12/
     *         ' TOTAL NUMBER OF MEMBERS   : ',I12/
     *         ' NUMBER OF MEMBERS = 0     : ',I12,'  OR ',F6.0,' %'/
     *         ' NUMBER OF NONZERO MEMBERS : ',I12,'  OR ',F6.0,' %'/
     *         ' REQUIRED DIMENSION IN MAIN ARRAY     : ',I12/
     *         ' REQUIRED DIMENSION IN AUXILIARY ARRAY: ',I12/
     *' TOTAL REQUIRED DIMENSION IN MAIN AND AUXILIARY ARRAY: ',I12//)
      RETURN
      END
C=======================================================================
      SUBROUTINE ISPAK0(SK,IROW,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
C      SUBROUTINE ISPAK(SK,IROW,MAXA,LM,ND,INDD,
C     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      DIMENSION LM(*),MAXA(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(NWK),SK(NWK),SKE(*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
C      IF(IILS.NE.-1) CALL KORIGC(A(LCORUL),A(LCMPC),MMP,NP,NEZAV)
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        IF(INDD.EQ.1)GO TO 9
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(KB0)
      MNQ1=MNQ(KB0+1)-1
      MXMN=MAXA(MNQ(KB0))-1
      LDB=MAXA(MNQ(KB0+1))-MAXA(MNQ(KB0))
      CALL RBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   IF(ICCGG.EQ.2) THEN
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ELSE
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ENDIF
C-----------------------------------------------
   11 NDI=0
99999 FORMAT(5I5)
      DO 200 I=1,ND
      II=LM(I)
C
C
      IF(II.LT.0)THEN
        IIP=-II
        ICM=MPC(1,IIP)
        DO 320 L=1,NEZAV
          II=MPC(L+1,IIP)
          IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 320
          CMI=CMPC(ICM,L)
          MI=MAXA(II)-MXMN
          KS=I
          DO 310 J=1,ND
            JJ=LM(J)
            IF(JJ)303,310,307
  303       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           CMJ=CMPC(JCM,K)
                IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
  666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 666
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSSU)
                   SK(KKNWK)=SK(KKNWK)+CMI*CMJ*SKE(KSSL)
                ELSE          
                  KK=MI+IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
 1666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 1666
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSS)
                ENDIF
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         IF(ICCGG.EQ.2) THEN
                KK=MI-IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
 6166               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6166
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
                SK(KK)=SK(KK)+CMI*SKE(KSSU)
                SK(KKNWK)=SK(KKNWK)+CMI*SKE(KSSL)
              ELSE
                KK=MI+IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
 6616               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6616
                KSS=KS
                IF(J.GE.I)KSS=J+NDI
                SK(KK)=SK(KK)+CMI*SKE(KSS)
              ENDIF
  310         KS=KS+ND-J
  320   CONTINUE
        GO TO 200
      ENDIF
      IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 200
      MI=MAXA(II)-MXMN
      KS=I
      DO 220 J=1,ND
      JJ=LM(J)
      IF(JJ)420,220,110
C
C
  420       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
 6661               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6661
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMJ*SKE(KSSU)
                  SK(KKNWK)=SK(KKNWK)+CMJ*SKE(KSSL)
                ELSE
                  KK=MI+IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 2666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 2666
                  SK(KK)=SK(KK)+CMJ*SKE(KSS)
                ENDIF
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 IF(ICCGG.EQ.2) THEN
         KK=MI-IJ
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
                  KK=KK+1
 6266               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6266
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
         SK(KK)=SK(KK)+SKE(KSSU)
         SK(KKNWK)=SK(KKNWK)+SKE(KSSL)
      ELSE
         KK=MI+IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6626               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6626
         KSS=KS
         IF(J.GE.I)KSS=J+NDI
         SK(KK)=SK(KK)+SKE(KSS)
      ENDIF
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
   20 CONTINUE
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C=======================================================================
      SUBROUTINE ISWR(IA,MAXA,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       PRINT PROFILED MATRIX
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION MAXA(*),IV(30)
      INTEGER*1 IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      WRITE(IZLAZ,5010) CHAR
      DO 20 I=1,N
      DO 10 J=I,N
        IV(J)=0.D0
        IJ=J-I
        KK=MAXA(J)+IJ
        IF(KK.LT.MAXA(J+1)) IV(J)=IA(KK)
10    CONTINUE
        WRITE(IZLAZ,5000) (IV(J),J=1,N)
        IV(I)=0
20    CONTINUE
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(30I5)
      END
C=======================================================================
      SUBROUTINE IISWR(IA,MAXA,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       PRINT PROFILED MATRIX
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION MAXA(*)
      INTEGER*1 IV(77500),IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      OPEN(99,FILE='FOR')
      WRITE(99,5010) CHAR
      DO 20 I=1,N
      DO 10 J=I,N
        IV(J)=0
        IJ=J-I
        KK=MAXA(J)+IJ
        IF(KK.LT.MAXA(J+1)) IV(J)=IA(KK)
10    CONTINUE
        WRITE(99,5000) (IV(J),J=1,N)
        WRITE(99,*) ' '
        IV(I)=0
20    CONTINUE
      CLOSE(99)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(80000I1)
      END
C=======================================================================
      SUBROUTINE IIISWR(LM,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       PRINT PROFILED MATRIX
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION LM(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      WRITE(IZLAZ,5010) CHAR
        WRITE(IZLAZ,5000) (LM(J),J=1,N)
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(30I5)
      END
C=======================================================================
      SUBROUTINE FORM(IROW,ISK,MAXA,JEDN,NWK,NWK1)
      COMMON /smumps/ imumps,ipar
      DIMENSION IROW(*),MAXA(JEDN+1),IVRS(JEDN)
      INTEGER*1 ISK(1+NWK1/7)
      DO 10 I=1,JEDN
         IVRS(I)=MAXA(I)+I+1-MAXA(I+1)
         MAXA(I)=MAXA(I+1)-1
   10 CONTINUE
      IPOZ=1
      DO 20 I=1,JEDN
         IMAX=IPOZ
         DO 30 JJ=I,JEDN
            IF(IVRS(JJ).EQ.I) THEN
               J=MAXA(JJ)
               CALL UZMI1(ISK(1+(J-1)/7),J-((J-1)/7)*7,IJ)
               IF(IJ.EQ.1)THEN
                  IROW(IPOZ)=JJ
                  IPOZ=IPOZ+1
               ENDIF
               MAXA(JJ)=MAXA(JJ)-1
               IVRS(JJ)=IVRS(JJ)+1
            ENDIF
   30    CONTINUE
         MAXA(I)=IMAX
   20 CONTINUE
      MAXA(JEDN+1)=IPOZ
      RETURN
      END
C=======================================================================
      SUBROUTINE ISP0(SK,MAXA,NWK,JEDN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SK(NWK),MAXA(JEDN+1)
      OPEN(7,FILE='0.SK')
      JE=JEDN
      IF (JEDN.GT.100) JE=100
      DO 10 I=1,JE
         DO 20 J=MAXA(I),MAXA(I+1)-1
            WRITE(7,1000)SK(J),I,I+MAXA(I)-J
   20    CONTINUE
   10 CONTINUE
      CLOSE(7)
      RETURN
 1000 FORMAT(E13.5,2I6)
      END

      SUBROUTINE ISP1(SK,IROW,MAXA,NWK,JEDN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SK(NWK),MAXA(JEDN+1),IROW(NWK)
      OPEN(9,FILE='1.SK')
      JE=JEDN
      IF (JEDN.GT.100) JE=100
      DO 10 I=1,JE
         DO 20 J=MAXA(I),MAXA(I+1)-1
            WRITE(9,1000) SK(J), I, IROW(J)
   20    CONTINUE
   10 CONTINUE
      CLOSE(9)
      RETURN
 1000 FORMAT(E13.5,2I6)
      END
C=======================================================================
      SUBROUTINE POST1(ISK,NP)
      INTEGER*1 N,ISK,IV(7)
      IV(7)= 1
      IV(6)= 2
      IV(5)= 4
      IV(4)= 8
      IV(3)=16
      IV(2)=32
      IV(1)=64
      N=ISK
      DO 10 I=1,NP-1
         IF(N.GE.IV(I)) N=N-IV(I)
   10 CONTINUE
      IF (N.LT.IV(NP)) ISK=ISK+IV(NP)
      RETURN
      END
C=======================================================================
      SUBROUTINE UZMI1(ISK,NP,IJ)
      INTEGER*1 N,ISK,IV(7)
      IV(7)= 1
      IV(6)= 2
      IV(5)= 4
      IV(4)= 8
      IV(3)=16
      IV(2)=32
      IV(1)=64
      N=ISK
      DO 10 I=1,NP-1
         IF(N.GE.IV(I)) N=N-IV(I)
   10 CONTINUE
      IJ=0
      IF (N.GE.IV(NP)) IJ=1
      RETURN
      END
C================================================================	
      SUBROUTINE FORM0(IROW,ISK,MAXA,MAXA8,JEDN,NWK,NWK1)
      COMMON /smumps/ imumps,ipar
      INTEGER*8 MAXA8(*),NWK1,J8
      DIMENSION IROW(*),MAXA(*),IVRS(JEDN)
      INTEGER*1 ISK(1+NWK1/7)
      IPOZ=1
      DO 20 I=1,JEDN
         IMAX=IPOZ
         DO 30 J8=MAXA8(I),MAXA8(I+1)-1
            CALL UZMI1(ISK(1+(J8-1)/7),J8-((J8-1)/7)*7,IJ)
C            CALL UZMI1(ISK(1+(JJ-1)/7),J-((JJ-1)/7)*7,IJ)
            IF(IJ.EQ.1)THEN
               IROW(IPOZ)=I+MAXA8(I)-J8
               if(imumps.eq.1) IROW(NWK+IPOZ)=I
               IPOZ=IPOZ+1
            ENDIF
C         DO 30 JJ=I,JEDN
C            IF(IVRS(JJ).EQ.I) THEN
C               J=MAXA(JJ)
C               CALL UZMI1(ISK(1+(J-1)/7),J-((J-1)/7)*7,IJ)
C               IF(IJ.EQ.1)THEN
C                  IROW(IPOZ)=JJ
C                  IPOZ=IPOZ+1
C               ENDIF
C               MAXA(JJ)=MAXA(JJ)-1
C               IVRS(JJ)=IVRS(JJ)+1
C            ENDIF
   30    CONTINUE
         MAXA(I)=IMAX
   20 CONTINUE
      MAXA(JEDN+1)=IPOZ
C      CALL IWRR(IROW,NN,'IROW')
      RETURN
      END
C=======================================================================
      SUBROUTINE ISPAKG(SK,IROW,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
C      SUBROUTINE ISPAK(SK,IROW,MAXA,LM,ND,INDD,
C     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      DIMENSION LM(*),MAXA(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(*),SK(*),SKE(*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
C      IF(IILS.NE.-1) CALL KORIGC(A(LCORUL),A(LCMPC),MMP,NP,NEZAV)
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        IF(INDD.EQ.1)GO TO 9
C        OPEN(99,FILE='PRISTUPI',STATUS='NEW',
C     &          FORM='FORMATTED',ACCESS='SEQENTIAL')
C     &          FORM='FORMATTED',ACCESS='DIRECT')
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(KB0)
      MNQ1=MNQ(KB0+1)-1
      MXMN=MAXA(MNQ(KB0))-1
      LDB=MAXA(MNQ(KB0+1))-MAXA(MNQ(KB0))
      CALL RBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   IF(ICCGG.EQ.2) THEN
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*ND)
        ELSE
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
        ENDIF
C-----------------------------------------------
   11 NDI=0
99999 FORMAT(5I5)
      DO 200 I=1,ND
      II=LM(I)
C
C
      IF(II.LT.0)THEN
        IIP=-II
        ICM=MPC(1,IIP)
        DO 320 L=1,NEZAV
          II=MPC(L+1,IIP)
          IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 320
          CMI=CMPC(ICM,L)
          MI=MAXA(II)-MXMN
          KS=I
          DO 310 J=1,ND
            JJ=LM(J)
            IF(JJ)303,310,307
  303       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           CMJ=CMPC(JCM,K)
                IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
  666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 666
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSSU)
C                  WRITE(99,99999)KSSU,KK
C                  WRITE(*,99999)KSSU,KK 
C                  PRINT *,KSSU,KK
                   SK(KKNWK)=SK(KKNWK)+CMI*CMJ*SKE(KSSL)
C                  WRITE(99,99999)KSSL,KKNWK
C                  WRITE(*,99999)KSSL,KKNWK
C                  PRINT *,KSSL,KKNWK
                ELSE          
                  KK=MI+IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 1666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 1666
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSS)
C                  WRITE(99,99999)KSSL,KK
C                  WRITE(*,99999)KSSL,KK 
C                  PRINT *,KSSL,KK
                ENDIF
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         IF(ICCGG.EQ.2) THEN
                KK=MI-IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6166               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6166
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
                SK(KK)=SK(KK)+CMI*SKE(KSSU)
C                WRITE(99,99999)KSSU,KK
C                WRITE(*,99999)KSSU,KK 
C                PRINT *,KSSU,KK
                SK(KKNWK)=SK(KKNWK)+CMI*SKE(KSSL)
C                WRITE(99,99999)KSSL,KKNWK
C                WRITE(*,99999)KSSL,KKNWK
C                PRINT *,KSSL,KKNWK
              ELSE
                KK=MI+IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6616               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6616
                KSS=KS
                IF(J.GE.I)KSS=J+NDI
                SK(KK)=SK(KK)+CMI*SKE(KSS)
C                WRITE(99,99999)KSS,KK
C                WRITE(*,99999)KSS,KK 
C                PRINT *,KSS,KK
              ENDIF
  310         KS=KS+ND-J
  320   CONTINUE
        GO TO 200
      ENDIF
      IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 200
      MI=MAXA(II)-MXMN
      KS=I
      IVRS=0
      DO 220 J=1,ND
      JJ=LM(J)
      IF(JJ.GT.0)IVRS=IVRS+1
      IF(JJ)420,220,110
C
C
  420       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6661               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6661
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMJ*SKE(KSSU)
                  SK(KKNWK)=SK(KKNWK)+CMJ*SKE(KSSL)
C                  WRITE(99,99999)KSSL,KK
C                  WRITE(*,99999)KSSL,KK 
C                  PRINT *,KSSL,KK
C                  WRITE(99,99999)KSSL,KKNWK
C                  WRITE(*,99999)KSSL,KKNWK
C                  PRINT *,KSSL,KKNWK
                ELSE
                  KK=MI+IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 2666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 2666
                  SK(KK)=SK(KK)+CMJ*SKE(KSS)
C                  WRITE(99,99999)KSSL,KK
C                  WRITE(*,99999)KSSL,KK 
C                  PRINT *,KSSL,KK
                ENDIF
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 IF(ICCGG.EQ.2) THEN
         KK=MI-IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6266               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6266
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
         SK(KK)=SK(KK)+SKE(KSSU)
C         WRITE(99,99999)KSSL,KK
C         WRITE(*,99999)KSSL,KK 
C         PRINT *,KSSL,KK
         SK(KKNWK)=SK(KKNWK)+SKE(KSSL)
C         WRITE(99,99999)KSSL,KK
C         WRITE(*,99999)KSSL,KK
C         PRINT *,KSSL,KK
      ELSE
         KK=MI+IJ+1
                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
 6626               CONTINUE
                       KK=KK-1
                    IF(IROW(KK).NE.II-IJ)GOTO 6626
         KSS=KS
         IF(J.GE.I)KSS=J+NDI
         SK(KK)=SK(KK)+SKE(KSS)
C         WRITE(99,99999)KSSL,KK
C         WRITE(*,99999)KSSL,KK
C         PRINT *,KSSL,KK
      ENDIF
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
C
C      CLOSE(99,STATUS='KEEP')
C
   20 CONTINUE
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C===================================================================
      SUBROUTINE ISPAK(SK,IROW,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
C      SUBROUTINE ISPAK(SK,IROW,MAXA,LM,ND,INDD,
C     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      DIMENSION LM(*),MAXA(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(NWK),SK(NWK),SKE(*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
C      IF(IILS.NE.-1) CALL KORIGC(A(LCORUL),A(LCMPC),MMP,NP,NEZAV)
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        IF(INDD.EQ.1)GO TO 9
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(KB0)
      MNQ1=MNQ(KB0+1)-1
      MXMN=MAXA(MNQ(KB0))-1
      LDB=MAXA(MNQ(KB0+1))-MAXA(MNQ(KB0))
      CALL RBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   IF(ICCGG.EQ.2) THEN
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*ND)
        ELSE
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
           call wrr6(ske,ND*(ND+1)/2,'skeR')
           call iwrr(lm,nd,'lmR ')
        ENDIF
c   10   IF(ICCGG.EQ.2) THEN
c           READ(ISCRC,END=15,ERR=999)
c     &     ND,(LM(I),I=1,ND)
c        ELSE
c           READ(ISCRC,END=15,ERR=999)
c     &     ND,(LM(I),I=1,ND)
c        ENDIF
C-----------------------------------------------
   11 NDI=0
99999 FORMAT(5I5)
      DO 200 I=1,ND
      II=LM(I)
C
C
      IF(II.LT.0)THEN
        IIP=-II
        ICM=MPC(1,IIP)
        DO 320 L=1,NEZAV
          II=MPC(L+1,IIP)
          IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 320
          CMI=CMPC(ICM,L)
          MI=MAXA(II)-MXMN
          KS=I
          DO 310 J=1,ND
            JJ=LM(J)
            IF(JJ)303,310,307
  303       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           CMJ=CMPC(JCM,K)
                IF(ICCGG.EQ.2) THEN
                  KK=MAXA(JJ+1)
  666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 666
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSSU)
                   SK(KKNWK)=SK(KKNWK)+CMI*CMJ*SKE(KSSL)
                ELSE
                  KK=MAXA(JJ+1)
 1666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 1666
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSS)
                ENDIF
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         IF(ICCGG.EQ.2) THEN
                KK=MAXA(JJ+1)
 6166               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 6166
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
                SK(KK)=SK(KK)+CMI*SKE(KSSU)
                SK(KKNWK)=SK(KKNWK)+CMI*SKE(KSSL)
              ELSE
                KK=MAXA(JJ+1)
 6616               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 6616
                KSS=KS
                IF(J.GE.I)KSS=J+NDI
                SK(KK)=SK(KK)+CMI*SKE(KSS)
              ENDIF
  310         KS=KS+ND-J
  320   CONTINUE
        GO TO 200
      ENDIF
      IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 200
      MI=MAXA(II)-MXMN
      KS=I
      DO 220 J=1,ND
      JJ=LM(J)
      IF(JJ)420,220,110
C
C
  420       JJP=-JJ
            JCM=MPC(1,JJP)
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           IF(ICCGG.EQ.2) THEN
                  KK=MAXA(JJ+1)
 6661               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 6661
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMJ*SKE(KSSU)
                  SK(KKNWK)=SK(KKNWK)+CMJ*SKE(KSSL)
                ELSE
                  KK=MAXA(JJ+1)
 2666               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 2666
                  SK(KK)=SK(KK)+CMJ*SKE(KSS)
                ENDIF
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 IF(ICCGG.EQ.2) THEN
              KK=MAXA(JJ+1)
6266               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 6266
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
         SK(KK)=SK(KK)+SKE(KSSU)
         SK(KKNWK)=SK(KKNWK)+SKE(KSSL)
      ELSE
                    KK=MAXA(JJ+1)
 6626               CONTINUE
                    KK=KK-1
                    IF(IROW(KK).NE.II)GOTO 6626
         KSS=KS
         IF(J.GE.I)KSS=J+NDI
         SK(KK)=SK(KK)+SKE(KSS)
      ENDIF
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
   20 CONTINUE
C      WRITE(3,*)'LSK',LSK
C	CALL DRSWRR(SK,MAXA,IROW,JEDN,'SK 1')
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C=======================================================================
C=======================================================================
      SUBROUTINE DRSWRR(A,MAXA,IROW,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       PRINT PROFILED MATRIX
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),MAXA(*),V(N),IROW(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      WRITE(IZLAZ,5010) CHAR
      DO 20 I=1,N
         DO 10 J=I,N
            V(J)=0.D0
10       CONTINUE
         DO 30 J=MAXA(I),MAXA(I+1)-1
            V(IROW(J))=A(J)
30       CONTINUE
        WRITE(IZLAZ,5000) (V(J),J=1,N)
        V(I)=0.D0
20    CONTINUE
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(30G10.3)
      END
	