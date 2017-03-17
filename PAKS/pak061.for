      SUBROUTINE ICCGMA(AA,B,MAXA,NWK,NN,KKK,IOUT,TOL,ALFA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'paka.inc'
      
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      DIMENSION AA(NWK),B(NN),MAXA(NN+1)
C********************************************************************
C.......A() su clanovi glavne matrice
C.......B() su clanovi desne strane jednacine
C.......MAXA() su dijagonalni clanovi matrice A()
C.......NWK je ukupan broj clanova matrice A()
C.......NN je dimenzija matrice A() i vektora B()
C.......NEED2 je potrebna dimenzija vektora IROW() i ne 0 clanova
C.......KKK=1 za nepotpunu faktorizaciju matrice A();
C       KKK=2 za resavanje desne strane jednacine
C.......IOUT je izlazna jedinica za poruke
C.......TOL je zadata tolerancija;Preporucena vrednost TOL=0.001
C.......ALFA je koficijent za PP. SHIFT; normalno je ALFA=0.0
C********************************************************************
      IF(KKK.EQ.1) THEN
CC
CC        WRITE(*,*)(A(IFF),IFF=1+NEED2*2+NN*5,1+NEED2*2+NN*5+2),IFF
        IF(DABS(ALFA).GT.1.D-20) CALL SHIFT(AA,MAXA,NN,ALFA)
C
C        CALL SWRR(AA,MAXA,NN,' SKA')
C        CALL IISWR(A(IROWS),MAXA,NN,'IROW')
        CALL ILU(AA,MAXA,NN,IOUT,A(LAILU),A(IROWS))
C        CALL SWRR(AA,MAXA,NN,' SKA')
C        CALL IISWR(A(IROWS),MAXA,NN,'IROW')
C        CALL SWRR(A(LAILU),MAXA,NN,'AILU')
C
CC
CC        WRITE(*,*)(A(IFF),IFF=1+NEED2*2+NN*5,1+NEED2*2+NN*5+2),IFF
        RETURN
      ELSEIF(KKK.EQ.2) THEN
CC
CC        WRITE(*,*)(A(IFF),IFF=1+NEED2*2+NN*5,1+NEED2*2+NN*5+2),IFF
        CALL ICCG(NN,TOL,B,A(LUCG),A(LVCG),A(LWCG),A(LPCG),A(LRCG),
     *            AA,A(LAILU),MAXA,A(IROWS),IOUT)
C
      ELSE
        WRITE(IOUT,*)' ***STOP (ICCGMAIN): KKK.NE.1.OR.KKK.NE.2 '
        STOP' ***ERROR 2 : (ICCGMAIN)'
      ENDIF
      RETURN
      END
C=====================================================================
      SUBROUTINE BWD(NN,BBAR,V,A,MAXA,IROW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),MAXA(*),V(*),BBAR(*),IROW(*)
C
      DO 40 I=1,NN
40    V(I)=BBAR(I)
      DO 10 N=NN,1,-1
      V(N)=V(N)/A(MAXA(N))
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL) 10,20,20
20    CONTINUE
      DO 30 L=KL,KU
      K=IROW(L)
30    V(K)=V(K)-A(L)*V(N)
10    CONTINUE
      RETURN
      END
C=====================================================================
      SUBROUTINE COUNT(A,NWK,NN,NZERO,NEED1,NEED2,NEED3,IOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
      NZERO=0
      DO 10 I=1,NWK
      IF(DABS(A(I)).LT.1.D-20) NZERO=NZERO+1
10    CONTINUE
      NNZERO=NWK-NZERO
      A1=FLOAT(NZERO)/FLOAT(NWK)*100.
      A2=FLOAT(NNZERO)/FLOAT(NWK)*100.
      NEED1=NNZERO*2+NN*5
      NEED2=NNZERO
      NEED3=NEED1+NEED2
C      WRITE(IOUT,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
      WRITE(*,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3
20    FORMAT(//' DIMENZIJA MATRICE    : ',I9/
     *         ' UKUPAN BROJ CLANOVA  : ',I9/
     *         ' BROJ CLANOVA = 0     : ',I9,'  ILI ',F6.0,' %'/
     *         ' BROJ NE NULA CLANOVA : ',I9,'  ILI ',F6.0,' %'/
     *         ' POTREBNA VELICINA GLAVNOG VEKTORA  : ',I9/
     *         ' POTREBNA VELICINA POMOCNOG VEKTORA : ',I9/
     *' UKUPNA POTREBNA VELICINA GLAVNOG I POMOCNOG VEKTORA : ',I9//)
      RETURN
      END
C=====================================================================
C=====================================================================
      SUBROUTINE FWD(NN,W,V,A,MAXA,IROW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),MAXA(*),V(*),W(*),IROW(*)
C
      DO 10 N=1,NN
      V(N)=W(N)
      C=0.
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL) 30,20,20
20    CONTINUE
      DO 40 KK=KL,KU
      K=IROW(KK)
40    C=C+A(KK)*V(K)
30    V(N)=(V(N)-C)/A(MAXA(N))
10    CONTINUE
      RETURN
      END
C=====================================================================
      SUBROUTINE ICCG(N,TOL,B,U,V,W,P,R,A,ALT,IADD,IROW,IOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION B(*),A(*),ALT(*),IADD(*),U(*),V(*),
     *W(*),P(*),R(*),IROW(*)
C
C
      CALL FWD(N,B,B,ALT,IADD,IROW)
      ANORMB=0.0
C
      DO 10 I=1,N
      R(I)=B(I)
      P(I)=R(I)
c      ANORMB=ANORMB+DABS(B(I))
      ANORMB=ANORMB+B(I)*B(I)
      B(I)=0.0
10    CONTINUE
      if(anormb.lt.1.e-10) stop 'anormb.lt.1.e-10'
      ANORMB=DSQRT(ANORMB)
C
      K=0
100   K=K+1
      WRITE(*,*) 'ITERACIJA : ',K
      ANORMR=0.0
      RR=0.0
      DO 20 I=1,N
C      ANORMR=ANORMR+DABS(R(I))
      ANORMR=ANORMR+R(I)*R(I)
      RR=RR+R(I)*R(I)
20    CONTINUE
      if(anormr.lt.1.e-10) stop 'anormr.lt.1.e-10'
      ANORMR=DSQRT(ANORMR)
      ANORC=ANORMR/ANORMB
C
      IF(ANORC.LT.TOL) GOTO 120
      IF(K-1) 70,70,50
50    BETA=RR/RR1
      DO 60 I=1,N
      P(I)=BETA*P(I)+R(I)
60    CONTINUE
70    CALL BWD(N,P,V,ALT,IADD,IROW)
      CALL MULTPL(N,W,V,A,IADD,IROW)
      CALL FWD(N,W,U,ALT,IADD,IROW)
      H=0.0
      DO 80 I=1,N
      H=H+P(I)*U(I)
80    CONTINUE
      ALPHA=RR/H
      DO 90 I=1,N
      B(I)=B(I)+ALPHA*P(I)
      R(I)=R(I)-ALPHA*U(I)
90    CONTINUE
      RR1=RR
      GOTO 100
C
120   WRITE(IOUT,160)K
160   FORMAT(//' (ICCG): TOTAL ITERATIONS = ',I8///)
      CALL BWD(N,B,B,ALT,IADD,IROW)
      RETURN
      END
C=====================================================================
      SUBROUTINE ILU(A,MAXA,NN,IOUT,AILU,IROW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),MAXA(*),AILU(*),IROW(*)
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
     1WRITE(IOUT,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IOUT,6001)N
          A(KN)=1.D0
c          DO 10 IO=KL,KU
c   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
C******************************************************************
         IF(KH.LT.0) THEN
            AILU(KN)=A(KN)
            IF(AILU(KN).LT.1.D-20) THEN
               WRITE(IOUT,2000) N,AILU(KN)
               STOP ' ***ERROR 1 : (ILU)'
            ENDIF
            AILU(KN)=DSQRT(AILU(KN))
C******************************************************************
         ELSEIF(KH.GE.0) THEN
            K=IROW(KU)
            IC=0
            KLT=KU
C
            AILU(KU)=A(KU)/AILU(MAXA(K))
            B=AILU(KU)*AILU(KU)
C
            DO 80 J=1,KH
               IC=IC+1
               KLT=KLT-1
               K=IROW(KLT)
               KI=MAXA(K)
               ND=MAXA(K+1)-KI-1
               IF(ND.GT.0) THEN
                  KK=MIN0(IC,ND)
                  C=0.
                  MML=KI+1
                  MMR=KLT+1
                  DO 70 L=1,KK
45                   IF(MML.GT.(KI+ND)) GOTO 75
                     IF(IROW(MML).EQ.IROW(MMR)) THEN
                        C=C+AILU(MML)*AILU(MMR)
                        MMR=MMR+1
                        IF(MMR.GT.KU) GOTO 75
                        GOTO 70
                     ENDIF
                     IF(IROW(MML).LT.IROW(MMR)) THEN
                        MMR=MMR+1
                        IF(MMR.GT.KU) GOTO 75
                        GOTO 45
                     ENDIF
                     MML=MML+1
                     GOTO 45
70                CONTINUE
75                AILU(KLT)=A(KLT)-C
               ELSE
                  AILU(KLT)=A(KLT)
               ENDIF
C
               AILU(KLT)=AILU(KLT)/AILU(MAXA(K))
               B=B+AILU(KLT)*AILU(KLT)
C
80          CONTINUE
C
            AILU(KN)=A(KN)-B
            aminus=1.
            IF(AILU(KN).LT.1.D-20) THEN
               WRITE(IOUT,2000) N,AILU(KN)
               WRITE(*,*) ' PIVOT < 0 '
CZILE
               IF(dabs(AILU(KN)).LT.1.D-20) STOP ' ***ERROR 2 : (ILU)'
               AILU(KN)=DABS(AILU(KN))
               aminus=-1.
CZILE
            ENDIF
            IF(AILU(KN).GT.1.D-10) AILU(KN)=aminus*DSQRT(AILU(KN))
         ENDIF
C
         IF(A(KN).GT.AMAXF) AMAXF=A(KN)
         IF(A(KN).LT.AMINF) AMINF=A(KN)
C
140   CONTINUE
      RETURN
2000  FORMAT(//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//
     *' NONPOSITIVE PIVOT FOR EQUATION ',I6,//
     *' PIVOT = ',E20.12)
C-----------------------------------------------------------------------
C 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'//
C     1  ' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I8,//' ','PIVOT=',
C     2  D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I8)
C-----------------------------------------------------------------------
C 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
C     1  ' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I8,//' ','PIVOT=',
C     2  D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I8)
C-----------------------------------------------------------------------
      END
C=====================================================================
      SUBROUTINE MULTPL(N,R,V,A,MAXA,IROW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION A(*),MAXA(*),V(*),R(*),IROW(*)
C
      DO 10 IQ=1,N
10    R(IQ)=0.0
C
      DO 20 IV =1,N
C
CC      R(IROW(MAXA(IV)))=R(IROW(MAXA(IV)))+A(MAXA(IV))*V(IV)
      R(IV)=R(IV)+A(MAXA(IV))*V(IV)
C
      DO 30 I=MAXA(IV)+1,MAXA(IV+1)-1
C
CE    RIGHT FROM DIAGONAL
CS    DESNO OD DIJAGONALE
C
      R(IROW(I))=R(IROW(I))+A(I)*V(IV)
C
CE    LEFT FROM DIAGONAL
CS    LEVO OD DIJAGONALE
C
      R(IV)=R(IV)+A(I)*V(IROW(I))
30    CONTINUE
20    CONTINUE
      RETURN
      END
C=====================================================================
      SUBROUTINE PACKFI(A,MAXA,NN,IROW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),MAXA(*),IROW(*)
      OPEN(50,FILE='ICCG.TMP',ACCESS='SEQUENTIAL',
     *     STATUS='UNKNOWN',FORM='UNFORMATTED')
      J=1
      DO 10 I=1,NN
      KL=MAXA(I)+1
      KU=MAXA(I+1)-1
      II=MAXA(I)
      A(J)=A(MAXA(I))
      WRITE(50) I
      MAXA(I)=J
      J=J+1
      DO 30 K=KL,KU
      IF(DABS(A(K)).LT.1.D-20) GOTO 30
        A(J)=A(K)
        WRITE(50) I-(K-II)
        J=J+1
30    CONTINUE
10    CONTINUE
      MAXA(NN+1)=J
      REWIND(50)
      DO 20 I=1,J-1
20    READ(50) IROW(I)
      CLOSE(50,STATUS='DELETE')
      RETURN
      END
C=====================================================================
      SUBROUTINE SHIFT(A,MAXA,NN,ALFA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),MAXA(*)
      COEF=1./(1.+ALFA)
      DO 10 N=1,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL) 10,20,20
20    DO 30 KK=KL,KU
      A(KK)=A(KK)*COEF
30    CONTINUE
10    CONTINUE
      RETURN
      END
