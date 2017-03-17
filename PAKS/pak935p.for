C=======================================================================
C
C   DINAMIKA 3/D KONTAKTA
C
C=======================================================================
      SUBROUTINE MASE93P(AM,MAXA,ID,EMAS,NEL,NGLOB,NP,NE,NTTOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PRIPADAJUCE MASE CVOROVA
CE     NODES MASSES
C
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION AM(*),MAXA(*),EMAS(*),ID(NP,*),NGLOB(*),NEL(NE,*)
C       WRITE(IZLAZ,*) 'MASS, NE'
      DO 10 NLM=1,NE
        NCK=NEL(NLM,1)
        EMAS(NLM)=0.D0
        DO 5 JJ=1,3
         IEQ=ID(NCK,JJ)
         IF(IEQ.GT.0)THEN
           IF(IMASS.EQ.1)THEN
             EMAS(NLM)=AM(MAXA(IEQ))
           ELSE
             EMAS(NLM)=AM(IEQ)
           ENDIF
         ENDIF
C      WRITE(IZLAZ,*) 'MASS, NLM', EMAS(NLM)
    5   CONTINUE
   10 CONTINUE
C
C      WRITE(IZLAZ,*) 'MASS, NTTOT'
      DO 30 NCKT=1,NTTOT
      IC=NGLOB(NCKT)
      NN=NE+NCKT
      EMAS(NN)=0.D0
        DO 25 JJ=1,3
         IEQ=ID(IC,JJ)
         IF(IEQ.GT.0)THEN
           IF(IMASS.EQ.1)THEN
             EMAS(NN)=AM(MAXA(IEQ))
           ELSE
             EMAS(NN)=AM(IEQ)
           ENDIF
         ENDIF
C      WRITE(IZLAZ,*) 'MASS, NN', EMAS(NN)
   25   CONTINUE
   30 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE VAUPD3P(VA,VA0,ID,RC,EMAS,AA,CC,BB,ISNA,IK,IK1,ALFK,
     &                  NELSF,ITSRF,IDC,NEL,NCVSF,SIL,LM,U,
     &                  NP,NE,NTTOT,NTSF,DT,IVA,IKV,JDIAG,NGLOB,
     &                  NELAB,XYZ,HE,TRA,S,LD,NCVE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     KOREKCIJA UBRZANJA I BRZINA PRI UDARU
CE     UPDATE ACCELERATION AND VELOCITY AFFTER IMPACT
C
C   IVA = 7  UBRZANJA
C   IVA = 8  BRZINE

      include 'paka.inc'
       
      LOGICAL AFAC,BACK
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      DIMENSION VA(*),ID(NP,*),EMAS(*),AA(*),CC(*),BB(NTTOT,*),ISNA(*),
     &          IK(*),IK1(*),ALFK(5,*),ITSRF(*),NEL(NE,*),NCVSF(*),
     &          SIL(3,*),IDC(NE,*),VA0(*),RC(*),LM(*),U(*),DA(3),
     &          IKV(*),JDIAG(*),NGLOB(*),NELSF(NCVE,*),NELAB(NCVE,*),
     &          XYZ(3,*),HE(9,*),TRA(3,*),S(NCVE,*),LD(*),ILM(3)
      DIMENSION VAL(3),VALT(3),VALC(3),VAGT(3),VAGC(3)
      DATA IT1/1000000/
CE....  IF SLIDING, USE ONLY NORMAL COMPONENTS
C
      NCVE2=NCVE*NCVE
      DT2=-1.D0
      IF(IVA.EQ.7)DT2=DT*.5
C
C.. MATRICE   AA, CC   I   BB
C
C.. FORMIRANJE MATRICE AA,CC, BB   BEZ KONTAKTORA
      I=1
      DO 10 N=1,NTTOT
      EM=EMAS(NE+N)
      JD=JDIAG(N)
      AA(JD)=EM
      CC(JD)=1.D0
      IC=NGLOB(N)
C      WRITE(IZLAZ,*) 'AA, JD', AA(JD),JD
C      WRITE(IZLAZ,*) 'CC, JD', CC(JD),JD
       DO 5 II=1,3
        IEQ=ID(IC,II)
        IF(IEQ.GT.0) BB(I,II)=EM*VA0(IEQ)
C      WRITE(IZLAZ,*) 'BB, I, II', BB(I,II), I,II
    5  CONTINUE
   10 I=I+1
C.. DODATAK MATRICAMA  AA, CC, BB   OD KONTAKTORA
      DO 30 NLM=1,NE
       IC  =NEL(NLM,1)
       ISRF=NEL(NLM,2)
       MAXCE=ISNA(ISRF)
       NODMAX=MAXCE+1
       LL  =IABS(ITSRF(ISRF))-1
C
       IND0=IK(NLM) /IT1
       IND1=IK1(NLM)/IT1
C
C
      IF(IND1.GT.0)THEN
       NCAA=IK1(NLM)-IND1*IT1
       NCA =NCAA+LL
C
       EM  =EMAS(NLM)
       R0  =ALFK(1,NLM)
       S0  =ALFK(2,NLM)
C
C  KOREKCIJA KOORDINATA
C
       ICAL=2
       CALL UPXYZ(XYZ,NLM,NELSF,A(LCORD),U,ID,NP,MAXCE,ICAL)
       ICAL=1
       CALL UPXYZ(XYZ,NLM,NELSF(1,NCA),A(LCORD),U,ID,NP,MAXCE,ICAL)
C
C   INTERPOLACIJE HE()  I   MATRICA TRANSFORMACIJE  TRA()
C
       CALL TRAHE(XYZ,HE,TRA,R0,S0,MAXCE)
       ICT = NODTRG(NELAB,HE,LD,KLR,NCA,NCVE)
       CALL CLEAR(S,NCVE2)
       DO 17 KM=1,NCVE
        IF(LD(KM).EQ.0)GO TO 17
        DUM=EM*HE(KM,1)
        DO 16 KK=1,NCVE
        IF(LD(KK).EQ.0)GO TO 16
        S(KM,KK)=DUM*HE(KK,1)
      WRITE(IZLAZ,*) 'SS, KM, KK', S(KM,KK), KM,KK
   16  CONTINUE
   17  CONTINUE
	   write (IZLAZ,*) "935p"
         do ii=1,18
	   write (IZLAZ,*) "i",ii,ld(ii)
	   enddo
        call WRR(AA,NTTOT,'AA   ')
C        call WRR(CC,NTTOT,'CC   ')
       CALL ADDS93(AA,CC,S,JDIAG,LD,NCVE,MAXCE)
C        call WRR(AA,NTTOT,'AA   ')
C        call WRR(CC,NTTOT,'CC   ')
      ENDIF
C.. MATRICA   BB
      DO 21 KM=1,NCVE
        ICT=LD(KM)
        IF(ICT.EQ.0)GO TO 21
        DUM=EM*HE(KM,1)
       DO 20 II=1,3
        IEQ=ID(IC,II)
        IF(IEQ.GT.0.AND.IND1.GT.0)THEN
C**        IF(IEQ.GT.0.AND.NCAA.GT.0)THEN
C*         IF(IND1.GT.0) BB(ICT,II)=BB(ICT,II)+EM*VA0(IEQ)
         BB(ICT,II)=BB(ICT,II)+DUM*VA0(IEQ)
C      WRITE(IZLAZ,*) 'BB, ICT, II', BB(ICT,II), ICT,II
         IF(IVA.EQ.7.OR.(IVA.EQ.8.AND.IND1.EQ.0.AND.IND0.GT.0))
     1   BB(ICT,II)=BB(ICT,II)+DT2*HE(KM,1)*SIL(II,NLM)
C      WRITE(IZLAZ,*) 'SILE, II,NLM', SIL(II,NLM), II,NLM
C      WRITE(IZLAZ,*) 'BB, ICT, II', BB(ICT,II), ICT,II
        ENDIF
   20  CONTINUE
   21  CONTINUE
   30 CONTINUE
C----------------------------------------
C.. NALAZENJE KORIGOVANIH VREDNOSTI
C
      AFAC = .TRUE. 
      BACK = .TRUE.
      CALL UACTCL(AA,CC,BB(1,1),JDIAG,NTTOT,AFAC,BACK)
      AFAC = .FALSE. 
      CALL UACTCL(AA,CC,BB(1,2),JDIAG,NTTOT,AFAC,BACK)
      CALL UACTCL(AA,CC,BB(1,3),JDIAG,NTTOT,AFAC,BACK)
C----------------------------------------
       call WRR(AA,NTTOT,'AA   ')
       call WRR(BB,NTTOT,'BB   ')
       call WRR(CC,NTTOT,'CC   ')
C.. RASPOREDJIVANJE KORIGOVANIH VREDNOSTI
C
C
C.. CILJ
C
      DO 50 ICT=1,NTTOT
      IC=NGLOB(ICT)
       DO 35 II=1,3
        IEQ=ID(IC,II)
        IF(IEQ.GT.0)THEN
C010795
C... KORIGOVANA VREDNOST ILI IZRACUNATA PREKO NJUMARKA
C    (AKO POSTOJI RAZLIKA U ODNOSU NA NEKORIGOVANU)
        IF(DABS(VA0(IEQ)-BB(ICT,II)).GT.1.D-15)THEN
          VA(IEQ)=BB(ICT,II)
        ELSE
          BB(ICT,II)=VA(IEQ)
        ENDIF
C**          IF(IVA.EQ.7) VA(IEQ)=BB(ICT,II)
        ELSE
          BB(ICT,II)=0.D0
        ENDIF
   35  CONTINUE
   50 CONTINUE
C
C.. KONTAKTOR
C
      DO 80 NLM=1,NE
       NCK=NEL(NLM,1)
       ITS=NEL(NLM,2)
       IND0=IK(NLM) /IT1
       IND1=IK1(NLM)/IT1
       INDV=IKV(NLM)/IT1
C.. NO CHANGE OF STATUSS
       IMPC=0
C.. IMPACT
       IF(IND0.EQ.0.AND.IND1.GT.0)IMPC=1
C.. RELEASE
       IF(IND0.GT.0.AND.IND1.EQ.0)IMPC=-1
C
       NCAA=IK1(NLM)-IND1*IT1
       LL  =IABS(ITSRF(ITS))-1
       NCA =NCAA+LL
C
       EM  =EMAS(NLM)
       R0  =ALFK(1,NLM)
       S0  =ALFK(2,NLM)
C
       ICAL=2
       CALL UPXYZ(XYZ,NLM,NELSF,A(LCORD),U,ID,NP,MAXCE,ICAL)
       ICAL=1
C*270894
      IF(NCAA.GT.0)THEN
       CALL UPXYZ(XYZ,NLM,NELSF(1,NCA),A(LCORD),U,ID,NP,MAXCE,ICAL)
       CALL TRAHE(XYZ,HE,TRA,R0,S0,MAXCE)
       ICT = NODTRG(NELAB,HE,LD,KLR,NCA,NCVE)
      ENDIF
C.. PODACI O PRIPADAJUCIM JEDNACINAMA
C Lagrange
C      LMILM=0
C      DO 135 II=1,3
C        KK=II*(MAXCE+2)-MAXCE-1
C        LM(KK)=ID(NCK,II)
C        ILM(II)=IDC(NLM,II)
C        LM(KK+1)=ILM(II)
C        IF(LMILM.EQ.0) LMILM=ILM(II)
C        K1=KK+1
C        IF(NCAA.GT.0)THEN
C          DO 132 L=1,MAXCE
C          NCV=NELSF(L,NCA)
C  132     LM(K1+L)=ID(NCV,II)
C        ENDIF
C  135 CONTINUE
C KRAJ Lagrange
C PENALTY
         IF(MAXCE.EQ.2) THEN
           K2=2
         ELSE
           K2=3
         ENDIF
      DO 135 K=1,K2        
         LM(K) =ID(NCK,K)  
C       WRITE(3,*)'LM(KK),KK',LM(K),K
  135 CONTINUE
       IF(NCAA.GT.0)THEN
        DO 132 L=1,MAXCE
          NCV=NELSF(L,NCA)
          K1=K2*L
          DO 133 K=1,K2
           KK=K1+K
  133        LM(KK)=ID(NCV,K)
C   33   WRITE(3,*)'LM(KK),KK',LM(KK),KK
  132   CONTINUE
       ENDIF
C KRAJ Penalty
C
      IF(IND1.NE.0)THEN
        DO 55 II=1,3
         IEQ=ID(NCK,II)
         DA(II)=0.D0
         IF(IEQ.GT.0) DA(II) =VA(IEQ)
   55   CONTINUE
C        IF(IND1.EQ.2) CALL RSN123(VAL,DA,TRA,2)
        CALL RSN123(VAL,DA,TRA,2)
C
        IF(NCAA.GT.0)THEN
          DO 61 II=1,3
   61     VALC(II)=0.D0
          DO 145 KK=1,MAXCE
            DO 62 II=1,3
               NN=NELAB(KK,NCA)
               NG=NGLOB(NN)
               IEQKK=ID(NG,II)
               VAGT(II)=BB(NN,II)
               IF(IEQKK.GT.0) VA(IEQKK)=VAGT(II)
   62       CONTINUE
C            IF(IND1.EQ.2) CALL RSN123(VALT,VAGT,TRA,2)
            CALL RSN123(VALT,VAGT,TRA,2)
            DO 63 II=1,3
   63       VALC(II)=VALC(II)+HE(KK,1)*VALT(II)
  145     CONTINUE
          IF(IND1.EQ.2)THEN
            VALC(1)=VAL(1)
            VALC(2)=VAL(2)
C            CALL RSN123(VALC,VAGC,TRA,1)
          ENDIF
            CALL RSN123(VALC,VAGC,TRA,1)
        ENDIF
C
        DO 64 II=1,3
         IEQ=ID(NCK,II)
         IF(IEQ.GT.0) VA(IEQ)=VAGC(II)
         DA(II) =DA(II)-VAGC(II)
   64   CONTINUE
      ENDIF 
C
C
C
      IF(IVA.EQ.7)THEN
        IF(IMPC.EQ.-1)THEN
          DO 65 II=1,3
C           KK=II*(MAXCE+2)-MAXCE
           IEQ=ID(NCK,II)
           LMKK=LM(II)
C           LMKK=LM(KK)
           IF(IEQ.GT.0.AND.LMKK.GT.0)VA(IEQ)=VA(IEQ)-U(LMKK)/EM
C           IF(IEQ.GT.0) VA(IEQ)=VA(IEQ)+SIL(II,NLM)/EM
           SIL(II,NLM)=0.D0
   65     CONTINUE
        ELSEIF(IMPC.EQ.1)THEN
          IF(IND1.NE.0)THEN
            DO 70 II=1,3
   70       SIL(II,NLM)=SIL(II,NLM)+EM*DA(II)          
          ENDIF
      IF(NCAA.GT.0)THEN
         DO 72 JJ=1,MAXCE
          NCT=NE+NELAB(JJ,NCA)
          DO 71 II=1,3
   71     SIL(II,NCT)=SIL(II,NCT)-HE(JJ,1)*SIL(II,NLM)
   72   CONTINUE
      ENDIF
C
C  3)  RASPOREDJIVANJE KONTAKTNIH SILA
C
      WRITE(IZLAZ,*) 'KONTAKTNE SILE, DINAMKA'
      DO 93 II=1,3
C       KK=II*(MAXCE+2)-MAXCE
C       KKLL=LM(KK)
       KKLL=LM(II)
       IF(KKLL.NE.0) THEN
        KKL=LM(II)
C        KKL=LM(KK-1)
        IF(KKL.NE.0) RC(KKL)=SIL(II,NLM)
        WRITE(IZLAZ,*) 'KKL, II' ,RC(KKL), KKL, II
        DO 90 L=1,MAXCE
          KKL=LM(KK+L)
          NCT=NE+NELAB(L,NCA)
          IF(KKL.NE.0) RC(KKL)=SIL(II,NCT)
          WRITE(IZLAZ,*) 'RC, KKL' ,RC(KKL), KKL
   90   CONTINUE
       ENDIF
   93 CONTINUE
C.. ENDIF   IVA.EQ.7
        ENDIF
      ENDIF
   80 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE UACTCLP(A,C,B,JDIAG,NEQ,AFAC,BACK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C.... UNSYMMETRIC,ACTIVE COLUMN PROFILE EQUATION SOLVER
C
      LOGICAL AFAC,BACK
      DIMENSION A(*),B(*),JDIAG(*),C(*)
C
C.... FACTOR A TO UT*D*U REDUCE B TO Y
C
      IF(.NOT.AFAC) GO TO 1
      DO 100 J = 1,NEQ
        JD = JDIAG(J)
        C(JD) = 1.D0
  100 CONTINUE
    1 JR = 0
      DO 300 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD - JR
      IF(JH.LE.1) GO TO 300
      IS = J + 1 - JH
      IE = J - 1
      IF(.NOT.AFAC) GO TO 250
      K = JR + 1
      ID = 0
C
C.... REDUCE ALL EQUATIONS EXCEPT DIAGONAL
C
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID - IR - 1,I - IS)
      IF(IH.EQ.0) GO TO 150
      A(K) = A(K) - DOT(A(K - IH),C(ID - IH),IH)
      C(K) = C(K) - DOT(C(K - IH),A(ID - IH),IH)
  150 IF(DABS(A(ID)).GT.1.D-10) C(K) = C(K)/A(ID)
  200 K = K + 1
C
C.... REDUCE DIAGONAL TERM
C
      A(JD) = A(JD) - DOT(A(JR + 1),C(JR + 1),JH - 1)
C
C.... FORWARD REDUCE THE R.H.S.
C
  250 IF(BACK) B(J) = B(J) - DOT(C(JR + 1),B(IS),JH - 1)
  300 JR = JD
      IF(.NOT.BACK) RETURN
C
C.... BACKSUBSTITUTION
C
      J = NEQ
      JD = JDIAG(J)
  500 IF(DABS(A(JD)).GT.1.D-10) B(J) = B(J)/A(JD)
      D = B(J)
      J = J - 1
      IF(J.LE.0) RETURN
      JR = JDIAG(J)
      IF(JD - JR.LE.1) GO TO 700
      IS = J - JD + JR + 2
      K = JR - IS + 1
      DO 600 I = IS,J
  600 B(I) = B(I) - A(I + K)*D
  700 JD = JR
      GO TO 500
      END                 
C=======================================================================
      SUBROUTINE ADDS93P(A,C,S,JDIAG,LD,NST,NEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C.... ASSEMBLE GLOBAL ARRAYS
C
      DIMENSION A(*),C(*),JDIAG(*),S(NST,*),LD(*)
      DO 200 J = 1,NEL
      K = LD(J)
      IF(K.EQ.0) GO TO 200
      L = JDIAG(K) - K
      DO 100 I = 1,NEL
      M = LD(I)
      IF(M.GT.K.OR.M.EQ.0) GO TO 100
      M0 = L + M
      A(M0) = A(M0) + S(I,J)
      IF(M.NE.K) C(M0) = C(M0) + S(J,I)
  100 CONTINUE
  200 CONTINUE
      RETURN
      END    
C=======================================================================
C=======================================================================
      SUBROUTINE PROF93P(JDIAG,NELAB,NELSF,NGLOB,NTTOT,NUMEL,NCVE,
     1NWKCDY)
C
C  PROFIL SISTEMA ZA DINAMIKU (IMPULSE-MOMENTUM BALANCE)
C
      DIMENSION JDIAG(*),NELAB(NCVE,*),NELSF(NCVE,*),NGLOB(*)
      CALL ICLEAR(JDIAG,NTTOT)
      NTT=1
C... COLUMN HEIGHTS
      DO 30 N=1,NUMEL
        DO 20 I=1,NCVE
          II=NELAB(I,N)
C
          IF(II.EQ.0)THEN
             GO TO 20
             ELSEIF(II.EQ.NTT)THEN
               NGLOB(NTT)=NELSF(I,N)
               NTT=NTT+1
          ENDIF
C
          DO 10 J=1,NCVE
            JJ=NELAB(J,N)
            IF(JJ.EQ.0) GO TO 10
            M=MAX0(II,JJ)
            JDIAG(M)=MAX0(JDIAG(M),IABS(II-JJ))
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C... DIAGONAL POINTERS FOR PROFILE
      NWKCDY=1
      JDIAG(1)=1
      IF(NTTOT.EQ.1) RETURN
      DO 40 N=2,NTTOT
   40 JDIAG(N)=JDIAG(N)+JDIAG(N-1)+1
      NWKCDY=JDIAG(NTTOT)
      RETURN
      END
