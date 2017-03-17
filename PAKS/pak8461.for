C=======================================================================
C
C     Isotropic Damage Model (Oliver 1996)  / 3D element
C
C     SUBROUTINE D8M61
C                TAUIL61
C                DYADSS
C
C=======================================================================
      SUBROUTINE D8M61(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    SPACE IN WORKING VECTOR A(*)
CE    POINTERS FOR VARIABLES IN AN INTEGRATION POINT FOR THE START AND END OF TIME STEP
C
      include 'paka.inc'
      
C      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
C 
      LTAU=LPOCG
      LDEFT=LTAU + 6*IDVA
      LDEFPT=LDEFT + 6*IDVA
      LTEQT=LDEFPT + 1*IDVA
      LTEQYT=LTEQT + 1*IDVA
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6*IDVA
      LDEFP1=LDEFT1 + 6*IDVA
      LTEQT1=LDEFP1 + 1*IDVA
      LTEQY1=LTEQT1 + 1*IDVA
C
      CALL TAUIL61(A(LTAU),A(LDEFT),A(LDEFPT),A(LTEQT),A(LTEQYT),
     1             A(LTAU1),A(LDEFT1),A(LDEFP1),A(LTEQT1),A(LTEQY1),
     1             A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE TAUIL61(S_N,EPSN,DN,RN,QN,
     1                   S_N1,EPSN1,DN1,RN1,QN1,
     1                   FUN,MATE,S_I1,EPS_I1,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    Isotropic Damage Model (Oliver 1996)
C
      COMMON/IZLE4B/Hja(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERBR/ ITER
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /PRINTS/ IPRS
C
      DIMENSION S_N(6),EPSN(6),S_N1(6),EPSN1(6),S_I1(6),EPS_I1(6)
      DIMENSION FUN(5,*),S_EFF(6),SS(6,6),DUM(6)
      DIMENSION CA(3,3,3,3),SA(3,3)
C
CE    MATERIAL DATA
C
      E=FUN(1,MAT)
      V=FUN(2,MAT)
      S_U=FUN(3,MAT)
      GF=FUN(4,MAT)
      idam=INT(FUN(5,MAT))
C
CE	Parameter controls type of damage model
CE	idam: 0 (linear), 1 (exponential)
C
c	idam=1
C  
CE    CALCULATION ELASTIC CONSTITUTIVE MATRIX
C
      CALL MEL861(E,V,ETP)
      IF(IRAC.EQ.2) RETURN
C
CE    INITIALIZATION
C
      R0=S_U/DSQRT(E)
      ZERO_Q=1.E-6*R0
      IF(RN.LT.1.E-10) THEN
         RN=R0
         QN=R0
      ENDIF
C
CE    EVALUATION OF EFFECTIVE STRESS
C          
      CALL CLEAR(S_EFF,6)
      CALL MNOZI1(S_EFF,ETP,EPS_I1,6,6)
C
CE    PREDICT ELASTIC TRIAL STRAIN
C
      R_TRIAL=dsqrt(DOT(S_EFF,EPS_I1,6))
C
      FLOAD=0.
C
CE    UPDATE STRESSES
C
      IF(R_TRIAL.GT.RN) THEN
         FLOAD=1.
         DELTA_R=R_TRIAL-RN
         RN1=R_TRIAL
      ELSE
         FLOAD=0.
         RN1=RN
      ENDIF
C
CE    SOFTENING REGULARIZATION
C
      H=-0.1
	IF(H.LT.1.E-10) THEN
	   HINTR=-S_U*S_U/(2.*E*GF)
c         A=-HINTR
c     approximatively
         A=R0*R0/GF
	   IF(idam.EQ.0) THEN
		 ALENGTH=1.0
           H=HINTR*ALENGTH
         ELSE
           H=A*((ZERO_Q-R0)/R0)*EXP(A*(1.-RN1/R0))
	   ENDIF
	ENDIF

      IF(R_TRIAL.GT.RN) THEN
         QN1=QN+H*DELTA_R
         IF(QN1.LT.ZERO_Q) QN1=ZERO_Q
      ELSE
         QN1=QN
      ENDIF
C
      DN1=1.-QN1/RN1

      CALL JEDNAK(S_I1,S_EFF,(1.-DN1),6)
      CALL JEDNA1(S_N1,S_I1,6)
      CALL JEDNA1(EPSN1,EPS_I1,6)
C
CE    TRANSFORM STRESSES TO GLOBAL COORDINATE SYSTEM
C
      IF(IPRS.EQ.2) THEN
         CALL CLEAR(DUM,6) 
         CALL MNOZI2(DUM,TSGD,S_N1,6,6)
         CALL JEDNA1(S_N1,DUM,6)
      ENDIF
c      if(nlm.eq.2) then
c        write(3,*) 'nlm,kor,iter,r_trial,rn',nlm,kor,iter,r_trial,rn
c        write(3,*) 'a,r0,gf,hintr,fload',a,r0,gf,hintr,fload
c        write(3,*) 'qn,h,delta_r,zero_q',qn,h,delta_r,zero_q
c        write(3,*) 'dn1,qn1,rn1,idam',dn1,qn1,rn1,idam
c        call wrr6(eps_i1,6,'eps1')
c        call wrr6(s_eff,6,'seff')
c        call wrr6(s_n1,6,'s_n1')
c      endif
C
CE    EVALUATION OF TANGENTIAL CONSTITUTIVE MATRIX
C
      CALL JEDNAK(ETP,ETP,(1.-DN1),36)
      IF(FLOAD.GT.1.E-10) THEN
         COEF=(QN1-H*RN1)/(RN1*RN1*RN1)
         CALL TENVEK(SA,S_EFF,0)
         CALL CEPMT(ETP,CA,0)
         CALL DYADSS(CA,SA,COEF)
         CALL CEPMT(ETP,CA,1)
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE DYADSS(C,S,COEF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO DYADIC PRODUCT C=SxS 
C .
C .       Cijkl = Sij x Skl
C .
C ......................................................................
C
      DIMENSION S(3,3),C(3,3,3,3)
C
      
      DO 10 I=1,3
      DO 10 J=1,3
      DO 10 K=1,3
      DO 10 L=1,3
         C(I,J,K,L)=C(I,J,K,L)-COEF*S(I,J)*S(K,L)
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE MEL861(E,V,ETP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    ELASTICITY MATRIX FOR SHELL ELEMENT
C
      DIMENSION ETP(6,*)
C
      EA=E/(1.D0-V*V)
      E1=0.5D0*(1.D0-V)
C
      CALL CLEAR(ETP,36)
      ETP(1,1)=EA
      ETP(2,2)=EA
      ETP(2,1)=EA*V
      ETP(1,2)=ETP(2,1)
      ETP(4,4)=EA*E1
      ETP(5,5)=ETP(4,4)
      ETP(6,6)=ETP(4,4)
      RETURN
      END
