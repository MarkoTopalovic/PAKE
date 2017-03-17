C=======================================================================
C
CE    BLOCKS
CS    BLOKOVI
C
C   SUBROUTINE   SPAKUJ
C                SPAKUA
C                RESEN
C                RESENA
C                ZADLEV
C                BLKDEF
C                RESENB
C                CLEARB
C                RBLOCK
C                WBLOCK
C                READDB
C                WRITEB
C
C=======================================================================
      SUBROUTINE SPAKUJ(SK,MAXA,SKE,LM,ND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - GLAVNI
C .
C ......................................................................
C
      include 'paka.inc'
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      common /elejob/ indjob
      DIMENSION SKE(*),LM(*),MAXA(*),SK(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUJ'
      IF(MDVI.EQ.3.AND.indjob.EQ.2) return
      IF(NBLOCK.EQ.1) THEN
CE       WITHOUT BLOCKS
CS       BEZ BLOKOVA
         IF(IABS(ICCGG).EQ.1)THEN
          IF(ICCGG.EQ.1)THEN
           CALL ISPAK(SK,A(IROWS),MAXA,SKE,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
          ELSE
           CALL ISPAKG(SK,A(IROWS),MAXA,SKE,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
          ENDIF
         ELSE
           CALL SPAKUA(SK,MAXA,SKE,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
         ENDIF
         ND2=ND*(ND+1)/2
         IF(ICCGG.EQ.2) ND2=ND*ND
         IF(ISKDSK.NE.0)
     &   WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND2)
      ELSE
CE       WITH BLOCKS
CS       SA BLOKOVIMA
         WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SPAKUA(SK,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
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
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      DIMENSION SKE(*),LM(*),MAXA(*),SK(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*)
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
        ENDIF
C-----------------------------------------------
   11 NDI=0
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
                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSSU)
                  SK(KKNWK)=SK(KKNWK)+CMI*CMJ*SKE(KSSL)
                ELSE
                  KK=MI+IJ
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
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
                SK(KK)=SK(KK)+CMI*SKE(KSSU)
                SK(KKNWK)=SK(KKNWK)+CMI*SKE(KSSL)
              ELSE
                KK=MI+IJ
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
                  KKNWK=KK+NWK
                  SK(KK)=SK(KK)+CMJ*SKE(KSSU)
                  SK(KKNWK)=SK(KKNWK)+CMJ*SKE(KSSL)
                ELSE
                  KK=MI+IJ
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
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
         SK(KK)=SK(KK)+SKE(KSSU)
         SK(KKNWK)=SK(KKNWK)+SKE(KSSL)
      ELSE
         KK=MI+IJ
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
C
C=======================================================================
      SUBROUTINE RESEN(B,V,MAXA,NN,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO SOLVE SYMMETRIC SYSTEM OF EQUATIONS
CS.   P R O G R A M
CS.      ZA RESAVANJE SISTEMA JEDNACINA
C .
C ......................................................................
C
      include 'paka.inc'
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      COMMON /RADIZA/ INDBG
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DUPLAP/ IDVA
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /FORSPA/ Lr1,Lz1,Lp1,LM0,LMaxM1,LColM
      COMMON /PRIVRM/ NM
      INCLUDE 'mpif.h'
      integer ierr, myid, k
      COMMON /smumps/ imumps,ipar
      DIMENSION B(*),V(*),MAXA(*)
      C(I1)=B(NWK+I1)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      if ( myid .ne. 0 ) goto 10
      k=kkk
C
      write(3,*)'kkk',kkk
      IF(IDEBUG.GT.0) PRINT *, ' RESEN'
      IF(ISRPS.EQ.0.AND.INDBG.EQ.0)
     1WRITE(*,2000) KKK
      IF(ISRPS.EQ.1.AND.INDBG.EQ.0)
     1WRITE(*,6000) KKK
C      DO 1 I=1,JEDN
C         WRITE(IZLAZ,12)I,B(MAXA(I)),V(I)
C         WRITE(IZLAZ,11)(B(J),J=MAXA(I),MAXA(I+1)-1)
C    1 CONTINUE
C      CALL SWRR(B,MAXA,JEDN,'SK  ')

   11 FORMAT(8E9.2)
   12 FORMAT(I5,2E9.2)
C     CALL SWRR(A(LSK),MAXA,JEDN,'SK U')
10    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iccgg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(imumps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(k,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(IABS(ICCGG).EQ.1) THEN
         IF(ICCGG.EQ.1) THEN
            if (myid.ne.0) goto 20
            IF(kkk.EQ.1)THEN
c              PSI=1.D-6
               psi=tolg
               MAXA(JEDN+1)=MAXA(JEDN)+1
C	         CALL DRSWRR(A(LSK),MAXA,A(IROWS),JEDN,'SK U')
               CALL ICM(JEDN,NWK,NNZERO,NM,NMREJ,PSI,
     +                  A(LSK),MAXA,A(IROWS),A(LM0),A(LMaxM1),A(LColM),
     +                  A(Lr1),A(Lz1),A(Lp1))
               PRINT*, NWK, NNZERO, NM
C	         CALL DRSWRR(A(LSK),MAXA,A(IROWS),JEDN,'SK I')
            ELSE
               CALL JEDNA1(A(Lr1),V,JEDN)
C              EPSILON=1.D-10
               EPSILON=1.D-8
               CALL ICM_CG(V,A(LSK),MAXA,A(IROWS),A(LM0),A(LMaxM1),
     +                  A(LColM),A(Lr1),A(Lz1),A(Lp1),
     +                  JEDN,NWK,NM,EPSILON,TOL,NITER)
20          ENDIF
         ELSE
            if(imumps.eq.1) then
c	          if(k.eq.2) then
                    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                    CALL dmumps1(A(IROWS),A(IROWS+nwk),B,V,nwk,nn,k)
                    IF (myid.ne.0) return
c	          end if
            else
                IF (myid.eq.0) CALL ICCGMA(B,V,MAXA,NWK,NN,
     1                                      K,IZLAZ,TOLG,ALFAG)
            endif
         ENDIF
      ELSE
         IF (myid.ne.0) goto 30
         IF(NBLOCK.EQ.1) THEN
CE          WITHOUT BLOCKS
CS          BEZ BLOKOVA
            IF(JPS.GT.1.AND.JPBR.LT.JPS1) THEN
               CALL RESENP(B,V,MAXA,NN,JEDNP,IZLAZ,K)
            ELSE
C              if(nn.le.30) call swrr(B,MAXA,NN,'B-  ')
C              call wrr6(B,NWK,'B-  ')
C              call wrr6(B,36,'B-  ')
C              if(nn.le.30) call wrr6(V,NN,'V-  ')
               IF(ICCGG.EQ.2) THEN
                  CALL UACTCF(B,C(1),V,MAXA,NN,K)
               ELSE
                  CALL RESENA(B,V,MAXA,NN,IZLAZ,K)
               ENDIF
c              if(nn.le.30) call swrr(B,MAXA,NN,'B+  ')
c              if(nn.le.30) call wrr6(V,NN,'V+  ')
            ENDIF
         ELSE
CE          WITH BLOCKS
CS          SA BLOKOVIMA
            CALL RESENB(B,V,MAXA,A(LMNQ),A(LICPL),
     1                  A(LLREC),NN,KC,NBLOCK,LR,IBLK,K)
30       ENDIF
C
      ENDIF
C
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(' *** RESAVANJE SISTEMA JEDNACINA (',I1,') ***')
C----------------------------------------------------------------------
 6000 FORMAT(' *** SOLUTION EQUATIONS SYSTEM (',I1,') ***')
C----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RESENP(A,V,MAXA,NN,JEDNP,IZLAZ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE------- RESEN FOR ONE PART SYSTEM -------------
CS------- RESEN ZA JEDNODELNI SITEM -------------
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RESENP'
      IF(JEDNP.EQ.0) RETURN
C          WRITE(3,*) 'KKK,ITER',KKK,ITER
C          call wrr6(A,36,'A-  ')
C          call wrr6(V,JEDNP+NN,'W-  ')
      IF(KKK-2)40,150,950
  40  DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      AMAXK=A(1)
      AMINK=A(1)
      AMAXF=A(1)
      AMINF=A(1)
C
C     FAKTORIZACIJA BLOKA 1 PREKO BLOKA 1
C
      DO 140 N=1,JEDNP
        KN=MAXA(N)
        KL=KN+1
        KU=MAXA(N+1)-1
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0'
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KH)110,90,50
  50    K=N-KH
        IC=0
        KLT=KU
        DO 80 J=1,KH
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)80,80,60
  60    KK=MIN0(IC,ND)
        C=0.D0
        DO 70 L=1,KK
  70    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
  80    K=K+1
  90    K=N
        B=0.D0
        DO 100 KK=KL,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 100    A(KK)=C
        A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
 115  IF(IPIVOT.EQ.1)IPIVOT=-1
      IF(IPIVOT)125,120,125
 120    CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)N,A(KN)
      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
      STOP
 125  IF(IDETER.EQ.0) GO TO 140
C     DETER=DETER*A(KN)
 140   CONTINUE
C
C     FAKTORIZACIJA BLOKA 2A PREKO BLOKA 1  I BLOKA 2B
C
      DO 540 N=JEDNP+1,NN
        KN=MAXA(N)
        KLA=KN+N-JEDNP
        KLB=KN+1
        KU=MAXA(N+1)-1
        KUB=KLB+N-JEDNP-1
        KHA=KU -KLA
        KHB=KUB-KLB
        KH =KU -KLB
C....  BLOK 2A
      IF(KHA)540,545,450
 450    K=JEDNP+1-KHA
        IC=0
        KLT=KU
        DO 480 J=1,KHA
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)480,480,460
 460    KK=MIN0(IC,ND)
        C=0.D0
        DO 470 L=1,KK
 470    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
 480    K=K+1
C....  BLOK 2B
 545  IF(KHB)490,490,550
 550    K=N-KHB
        KLT=KUB
        KJJ=K-JEDNP-1
        DO 580 J=1,KHB
        KJJ=KJJ+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-KJJ
        IF(ND)580,580,560
 560    KK=MIN0(ND,KHA+1)
        C=0.D0
        DO 570 L=KJJ,KJJ+KK-1
 570    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
 580    K=K+1
C
 490    K=JEDNP+1
        B=0.D0
        DO 500 KK=KLA,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 500    A(KK)=C
        A(KN)=A(KN)-B
 540   CONTINUE
       RETURN
CE     REDUCE  FREE  VECTOR
CS     REDUKOVANJE SLOBODNOG VEKTORA
C
C      STANDARDNA REDUKCIJA
C
 150    DO 180 N=1,JEDNP
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)180,160,160
 160    K=N
        C=0.
        DO 170 KK=KL,KU
        K=K-1
 170    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 180    CONTINUE
C
C      MODIFIKOVANA REDUKCIJA UPOTREBOM BLOKA 2A
C
        DO 780 N=JEDNP+1,NN
        KL=MAXA(N)+N-JEDNP
        KU=MAXA(N+1)-1
        IF(KU-KL)780,760,760
 760    K=JEDNP+1
        C=0.
        DO 770 KK=KL,KU
        K=K-1
 770    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 780    CONTINUE
CP
        IF(KKK.EQ.2)RETURN
CP
C       ZAMENA UNAZAD
 950    DO 200 N=1,JEDNP
        K=MAXA(N)
 200    V(N)=V(N)/A(K)
CCC***        IF(JEDNP.EQ.1)RETURN
C.....   KOREKCIJA OD POMERANJA GLOBALNIH CVOROVA 
        N=NN
        DO 230 L=JEDNP+1,NN
        KL=MAXA(N)+N-JEDNP
        KU=MAXA(N+1)-1
        IF(KU-KL)230,210,210
 210    K=JEDNP+1
        DO 220 KK=KL,KU
        K=K-1
 220    V(K)=V(K)-A(KK)*V(N)
 230    N=N-1
C
        N=JEDNP
        DO 930 L=2,JEDNP
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)930,910,910
 910    K=N
        DO 920 KK=KL,KU
        K=K-1
 920    V(K)=V(K)-A(KK)*V(N)
 930    N=N-1
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL CISTIP(V,JEDNP)
CZILE
        RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'//
     1  ' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I8,//' ','PIVOT=',
     2  D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I8)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1  ' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I8,//' ','PIVOT=',
     2  D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I8)
C-----------------------------------------------------------------------
        END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RESENA(A,V,MAXA,NN,IZLAZ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE       TO L*D*LT FAKTORIZE WITHOUT BLOCKS
CS.   P R O G R A M
CS.      ZA L*D*LT FAKTORIZACIJU BEZ BLOKOVA
C .
C ......................................................................
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RESENA'
      IF(NN.EQ.0) RETURN
      IF(KKK-2)40,150,150
  40  DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      AMAXK=A(1)
      AMINK=A(1)
      AMAXF=A(1)
      AMINF=A(1)
      DO 140 N=1,NN
        KN=MAXA(N)
        KL=KN+1
        KU=MAXA(N+1)-1
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15.AND.N.LE.NEQ) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0'
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KH)110,90,50
  50    K=N-KH
        IC=0
        KLT=KU
        DO 80 J=1,KH
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)80,80,60
  60    KK=MIN0(IC,ND)
        C=0.D0
        DO 70 L=1,KK
  70    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
  80    K=K+1
  90    K=N
        B=0.D0
        DO 100 KK=KL,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 100    A(KK)=C
        A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
CC 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.KN.LE.NEQ)))
 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.N.LE.NEQ)))
     $   IPIVOT=-1
      IF(IPIVOT)125,116,125
CC 116  IF(ICONT.NE.0.AND.KN.GT.NEQ)GO TO 125
 116  IF(ICONT.NE.0.AND.N.GT.NEQ)GO TO 125
 120    CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)N,A(KN)
      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
      STOP
 125  IF(IDETER.EQ.0) GO TO 140
C     DETER=DETER*A(KN)
 140   CONTINUE
C
C       WRITE(3,*) 'AMAXK,AMINK',AMAXK,AMINK
C       WRITE(3,*) 'AMAXF,AMINF',AMAXF,AMINF
       RETURN
CE     REDUCE  FREE  VECTOR
CS     REDUKOVANJE SLOBODNOG VEKTORA
 150    CONTINUE
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL CISTIN(V,NN)
CZILE
        DO 180 N=1,NN
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)180,160,160
 160    K=N
        C=0.D0
        DO 170 KK=KL,KU
        K=K-1
 170    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 180    CONTINUE
C       ZAMENA UNAZAD
        DO 200 N=1,NN
        K=MAXA(N)
CZILE
C        IF(DABS(A(K)-1.D0).LT.1.D-15.AND.N.LE.NEQ) V(N)=0.D0
CZILE
 200    V(N)=V(N)/A(K)
        IF(NN.EQ.1)RETURN
        N=NN
        DO 230 L=2,NN
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)230,210,210
 210    K=N
        DO 220 KK=KL,KU
        K=K-1
 220    V(K)=V(K)-A(KK)*V(N)
 230    N=N-1
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL CISTIP(V,NN)
CZILE
        RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'//
     1  ' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I8,//' ','PIVOT=',
     2  D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I8)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1  ' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I8,//' ','PIVOT=',
     2  D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I8)
C-----------------------------------------------------------------------
        END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RESENB(A,V,MAXA,MNQ,ICPL,LREC,NN,KC,NBLOCK,LR,IBLK,
     1  KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE       TO L*D*LT FAKTORIZE WITH BLOCKS
CS.   P R O G R A M
CS.      ZA L*D*LT FAKTORIZACIJU SA BLOKOVIMA
C .
C ......................................................................
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),V(*),MAXA(*),MNQ(*),ICPL(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RESENB'
      KC1=KC+1
      IF(KKK-2)40,150,150
   40 DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      DO 145 NB=1,NBLOCK
      IF(ISRPS.EQ.0)
     1WRITE(*,2002) NB,NBLOCK
      IF(ISRPS.EQ.1)
     1WRITE(*,6002) NB,NBLOCK
      CALL RBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      IF(NB.EQ.1)THEN
        AMAXK=A(1)
        AMINK=A(1)
        AMAXF=A(1)
        AMINF=A(1)
      ENDIF
      KB=ICPL(NB)
      KHB=NB-KB+1
      KI0=1
      IF(KHB.GT.1) KI0=KC1
CE    LOOP OVER COUPLED BLOCKS
CS    PETLJA PO SPREGNUTIM BLOKOVIMA
        DO 142 JB=1,KHB
        IF(KB.NE.NB)
     1CALL RBLOCK(A(KC1),LREC,KB,LR,MAXA(MNQ(KB+1))-MAXA(MNQ(KB)),IBLK)
CE    LOOP OVER EQUATIONS OF BLOCK NB
CS    PETLJA PO JEDNACINAMA BLOKA NB
      DO 140 N=MNQ(NB),MNQ(NB+1)-1
        KN=MAXA(N)-MAXA(MNQ(NB))+1
        KL=KN+1
        KU=MAXA(N+1)-MAXA(MNQ(NB))
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15.AND.N.LE.NEQ) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0'
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KB.LT.NB) THEN
                   IF(KH) 140,140,50
                   ELSE
                   IF(KH) 110,90,50
                   ENDIF
C
   50 K1=N-KH
CE    AKO JE K1.GE.MNQ(KB+1)  JEDNACINA N NIJE POVEZANA SA BLOKOM KB
CS    AKO JE K1.GE.MNQ(KB+1)  JEDNACINA N NIJE POVEZANA SA BLOKOM KB
         IF(K1-MNQ(KB+1)) 52,140,140
CE    K  JE PRVA JEDNACINA IZ BLOKA KB SA KOJOM JE SPREGNUTA J-NA N
CS    K  JE PRVA JEDNACINA IZ BLOKA KB SA KOJOM JE SPREGNUTA J-NA N
   52    K=MAX0(K1,MNQ(KB))
         IC=K-K1
CE     KLT  ADRESA PRVOG CLANA G (KASNIJE L)  KOJI JE POVEZAN SA  KB
CS     KLT  ADRESA PRVOG CLANA G (KASNIJE L)  KOJI JE POVEZAN SA  KB
       KLT=KU-IC
CE     KH   BROJ CLANOVA POVEZANIH SA  KB
CS     KH   BROJ CLANOVA POVEZANIH SA  KB
          IF(KB.LT.NB) THEN
          KH=MNQ(KB+1)-K
                ELSE
                KI0=1
                KH=N-K
                ENDIF
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF(ND)80,80,60
   60 KK=MIN0(IC,ND)
CE    KIB  ADRESSES OF DIAGONAL ELEMENT OF BLOCK KB
CS    KIB  ADRESA DIJAGONALNOG CLANA BLOKA KB
       KIB=KI-MAXA(MNQ(KB))+KI0
      C=0.0D0
      DO 70 L=1,KK
   70 C=C+A(KIB+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
C
   90 IF(KB.LT.NB)GO TO 140
      K=N
      B=0.0D0
      DO 100 KK=KL,KU
      K=K-1
      C=A(KK)/V(K)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.N.LE.NEQ)))
     $   IPIVOT=-1
      IF(IPIVOT)125,116,125
 116  IF(ICONT.NE.0.AND.N.GT.NEQ)GO TO 125
 120  CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)N,A(KN)
      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
      STOP
  125 IF(IDETER.EQ.0) GO TO 139
C     DETER=DETER*A(KN)
  139 V(N)=A(KN)
  140 CONTINUE
  142 KB=KB+1
      CALL WBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
CE    PAMTI DIJAGONALNE CLANOVE U DRUGOJ POLOVINI VEKTORA  A
CS    PAMTI DIJAGONALNE CLANOVE U DRUGOJ POLOVINI VEKTORA  A
      DO 147 N=1,NN
  147 A(KC1+N)=V(N)
  145 CONTINUE
      LMAXRC=LREC(NBLOCK+1)-1
      CALL WRITDD(A(KC1+1),NN,IBLK,LMAXRC,LDUZI)
      RETURN
C
CE    REDUKOVANJE SLOBODNOG VEKTORA
CS    REDUKOVANJE SLOBODNOG VEKTORA
  150 DO 185 NB=1,NBLOCK
      CALL RBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      DO 180 N=MNQ(NB),MNQ(NB+1)-1
      KL=MAXA(N)-MAXA(MNQ(NB))+2
      KU=MAXA(N+1)-MAXA(MNQ(NB))
      IF(KU-KL)180,160,160
  160 K=N
      C=0.0D0
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
  185 CONTINUE
C
CE    BACKSUBSTITUTION
CS    ZAMENA UNAZAD
C
      LMAXRC=LREC(NBLOCK+1)-1
      CALL READDD(A(KC1+1),NN,IBLK,LMAXRC,LDUZI)
      DO 200 N=1,NN
  200 V(N)=V(N)/A(KC1+N)
      IF(NN.EQ.1)RETURN
C
      N=NN
      NB=NBLOCK
      DO 235 LB=1,NBLOCK
      CALL RBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      L0=MNQ(NB)
      IF(NB.EQ.1) L0=2
      DO 230 L=L0,MNQ(NB+1)-1
      KL=MAXA(N)-MAXA(MNQ(NB))+2
      KU=MAXA(N+1)-MAXA(MNQ(NB))
      IF (KU-KL)230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
  235 NB=NB-1
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I8,//' ','PIVOT=',
     2D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I8)
 2002 FORMAT(' *** *** BLOK MATRICE SISTEMA',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
 6000 FORMAT(//' ',' MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I8,//' ','PIVOT=',
     2D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I8)
 6002 FORMAT(' *** *** BLOCK OF SYSTEM MATRIX',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZADLEV(B,MAXA,NZADJ,NZADP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRESCRIBED DISPLACEM. - UPDATING OF LEFT SIDE WITH BLOCKS
CS.   P R O G R A M
CS.      ZA ZADATA POMERANJA - KOREKCIJA LEVE STRANE SA BLOKOVIMA
C .
C ......................................................................
C
      include 'paka.inc'
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),NZADJ(*),IA(1)
      EQUIVALENCE (A(1),IA(1))
      MNQ(NN)=IA(LMNQ+NN-1)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZADLEB'
      VBR=1.0D35
      DO 20 NB=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(NB)
      MNQ1=MNQ(NB+1)-1
      MXMN=MAXA(MNQ(NB))-1
      LDB=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL RBLOCK(B,A(LLREC),NB,LR,LDB,IBLK)
   11   DO 10 I=1,NZADP
          NJ=NZADJ(I)
          IF(NJ.LT.MNQ0.OR.(NJ.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 10
          NDIJ=MAXA(NJ)-MXMN
          B(NDIJ) = VBR
   10   CONTINUE
      IF(NBLOCK.EQ.1)RETURN
      CALL WBLOCK(B,A(LLREC),NB,LR,LDB,IBLK)
   20 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CLEARB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,NUL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY WITH BLOCKS
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA SA BLOKOVIMA
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEARB'
      CALL CLEAR(B,NUL)
      IF(NBLOCK.EQ.1)RETURN
      DO 10 NB=1,NBLOCK
  10  CALL WBLOCK(B,LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RBLOCK(B,LREC,NB,LR,KC,IBLK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO READING OF BLOCK  NB
CS.   P R O G R A M
CS.      ZA CITANJE  BLOKA  NB
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RBLOCK'
      I=0
      DO 10 N=LREC(NB),LREC(NB+1)-1
      LL=LR
      IF(N.EQ.LREC(NB+1)-1) LL=KC-I
      READ(IBLK,REC=N) (B(I+J),J=1,LL)
   10 I=I+LR
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WBLOCK(B,LREC,NB,LR,KC,IBLK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO WRITE OF BLOCK  NB
CS.   P R O G R A M
CS.      ZA ZAPISIVANJE BLOKA  NB
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WBLOCK'
      I=0
      DO 10 N=LREC(NB),LREC(NB+1)-1
      LL=LR
      IF(N.EQ.LREC(NB+1)-1) LL=KC-I
      WRITE(IBLK,REC=N) (B(I+J),J=1,LL)
   10 I=I+LR
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READDB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,IU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSFER SEQUENTIAL FILE INTO BLOCKS
CS.   P R O G R A M
CS.      ZA PREBACIVANJE SEKVENCIJALNE DATOTEKE U BLOKOVE
C .
C ......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READDB'
      DO 10 NB=1,NBLOCK
      LL=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL READDD(B,LL,IPODS,IU,LDUZI)
   10 CALL WBLOCK(B,LREC,NB,LR,LL,IBLK)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE WRITEB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,IU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSFER BLOCKS INTO SEQUENTIAL FILE
CS.   P R O G R A M
CS.      ZA PREBACIVANJE  BLOKOVA U SEKVENCIJALNU DATOTEKU
C .
C ......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITEB'
      DO 10 NB=1,NBLOCK
      LL=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL RBLOCK(B,LREC,NB,LR,LL,IBLK)
   10 CALL WRITDD(B,LL,IPODS,IU,LDUZI)
      RETURN
      END
C=======================================================================
      SUBROUTINE SILPAK(F,FE,LM,ND,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE      ASSEMBLE INTERNAL FORCES
CS      RAZMESTANJE UNUTRASNJIH SILA
C
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION F(*),FE(*),LM(*),CMPC(MMP,*),MPC(NEZA1,*)
C
      DO 50 L=1,ND
      LL=LM(L)
      IF(LL.GT.0)THEN
        F(LL)=F(LL)+FE(L)
C
      ELSEIF(LL.LT.0)THEN
        LLP=-LL
        LCM=MPC(1,LLP)
        DO 20 K=1,NEZAV
          LL=MPC(K+1,LLP)
          IF(LL.EQ.0)GO TO 20
          CML=CMPC(LCM,K)
          F(LL)=F(LL)+CML*FE(L)
   20   CONTINUE
      ENDIF
C
   50 CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION CONDOF(U,CMPC,MPC,III)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE      TO CALCULATE CONSTRAIND DOF
CS      IZRACUNAVANJE ZAVISNOG STEPENA SLOBODE (MULTIPOINT CONSTRAINT)
C
      COMMON /CDEBUG/ IDEBUG
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION U(*),CMPC(MMP,*),MPC(NEZA1,*)
      IF(IDEBUG.GT.0) PRINT *, ' CONDOF'
C
      II=III
      IIP=-II
      ICM=MPC(1,IIP)
      CONDOF=0.D0
      DO 20 L=1,NEZAV
      II=MPC(L+1,IIP)
   20 IF(II.GT.0) CONDOF=CONDOF+U(II)*CMPC(ICM,L)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZADATL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO FORM POINTERS OF PRESCRIBED DISPLACEMENTS
CS.   P R O G R A M
CS.      ZA FORMIRANJE REPERA ZA ZADATA POMERANJA
C .
C ......................................................................
C
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DUPLAP/ IDVA
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' FORMCT'
      LZADFM = LRAD
      LNZADF = LZADFM + NZADP*IDVA
      LNZADJ = LNZADF + NZADP
      LMAX   = LNZADJ + NZADP
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      NPRO=LMAX-LRAD
      CALL DELJIV(NPRO,2,INDL)
      IF(INDL.EQ.1) NPRO=NPRO+1
C
      CALL READDD(A(LZADFM),NPRO/IDVA,IPODS,LMAX13,LDUZI)
      CALL ZADLEV(A(LSK),A(LMAXA),A(LNZADJ),NZADP)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SWRR(A,MAXA,N,CHAR)
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
      DIMENSION A(*),MAXA(*),V(30)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      WRITE(IZLAZ,5010) CHAR
      DO 20 I=1,N
      DO 10 J=I,N
        V(J)=0.D0
        IJ=J-I
        KK=MAXA(J)+IJ
        IF(KK.LT.MAXA(J+1)) V(J)=A(KK)
10    CONTINUE
        WRITE(IZLAZ,5000) (V(J),J=1,N)
        V(I)=0.D0
20    CONTINUE
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(30G10.3)
      END

