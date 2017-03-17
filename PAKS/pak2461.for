C=======================================================================
C
C     Isotropic Damage Model (Oliver 1996) / 3D element
C
C     SUBROUTINE D3M61
C                TAUI361
C                
C
C=======================================================================
      SUBROUTINE D2M61(TAU,DEF,IRAC,LPOCG,LPOC1)
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
      LDEFT=LTAU + 4*IDVA
      LDEFPT=LDEFT + 4*IDVA
      LTEQT=LDEFPT + 1*IDVA
      LTEQYT=LTEQT + 1*IDVA
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LTEQT1=LDEFP1 + 1*IDVA
      LTEQY1=LTEQT1 + 1*IDVA
C
      CALL TAUI261(A(LTAU),A(LDEFT),A(LDEFPT),A(LTEQT),A(LTEQYT),
     1             A(LTAU1),A(LDEFT1),A(LDEFP1),A(LTEQT1),A(LTEQY1),
     1             A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE TAUI261(S_N,EPSN,DN,RN,QN,
     1                   S_N1,EPSN1,DN1,RN1,QN1,
     1                   FUN,MATE,S_I1,EPS_I1,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    Isotropic Damage Model (Oliver 1996)
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERBR/ ITER
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /PRINTS/ IPRS
C
      DIMENSION S_N(6),EPSN(6),S_N1(6),EPSN1(6),S_I1(6),EPS_I1(6)
      DIMENSION FUN(5,*),S_EFF(6),SS(6,6),DUM(6)
      DIMENSION CA(3,3,3,3),SA(3,3),ELAS(3,3),ELA(6,6)
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
      CALL MEL261(E,V,ELAST,IETYP)
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
      DO I=1,3
      DO J=1,3
      ELAS(I,J)=ELAST(I,J)
      ENDDO
      ENDDO          
      CALL CLEAR(S_EFF,6)
      CALL MNOZI1(S_EFF,ELAS,EPS_I1,KK,KK)
C
CE    PREDICT ELASTIC TRIAL STRAIN
C
      R_TRIAL=dsqrt(DOT(S_EFF,EPS_I1,KK))
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
C
      IF(R_TRIAL.GT.RN) THEN
         QN1=QN+H*DELTA_R
         IF(QN1.LT.ZERO_Q) QN1=ZERO_Q
      ELSE
         QN1=QN
      ENDIF
C
      DN1=1.-QN1/RN1

      CALL JEDNAK(S_I1,S_EFF,(1.-DN1),KK)
      CALL JEDNA1(S_N1,S_I1,KK)
      CALL JEDNA1(EPSN1,EPS_I1,KK)
c      if(nlm.eq.2) then
c        write(3,*) 'nlm,kor,iter,r_trial,rn',nlm,kor,iter,r_trial,rn
c        write(3,*) 'a,r0,gf,hintr,fload',a,r0,gf,hintr,fload
c        write(3,*) 'qn,h,delta_r,zero_q',qn,h,delta_r,zero_q
c        write(3,*) 'dn1,qn1,rn1,idam',dn1,qn1,rn1,idam
c        call wrr6(eps_i1,KK,'eps1')
c        call wrr6(s_eff,KK,'seff')
c        call wrr6(s_n1,KK,'s_n1')
c      endif
C
CE    EVALUATION OF TANGENTIAL CONSTITUTIVE MATRIX
C
      CALL JEDNAK(ELA,ELAST,(1.-DN1),36)
      ELA(4,4)=ELA(3,3)
      ELA(3,3)=0.
      IF(FLOAD.GT.1.E-10) THEN
         COEF=(QN1-H*RN1)/(RN1*RN1*RN1)
         CALL TENVEK(SA,S_EFF,0)
         CALL CEPMT(ELA,CA,0)
         CALL DYADSS(CA,SA,COEF)
         CALL CEPMT(ELAST,CA,1)
         ELAST(3,3)=ELAST(4,4)
         ELAST(4,4)=0.
      ENDIF
      RETURN
      END
C=======================================================================
C
C======================================================================
      SUBROUTINE MEL261(E,V,ELAST,IETYP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    ELASTICITY MATRIX FOR 3D ELEMENT
C
      DIMENSION ELAST(6,*)
C
      CALL CLEAR(ELAST,36)
C
      IF(IETYP.EQ.0.OR.IETYP.EQ.3) THEN
CS      RAVANSKO STANJE NAPONA
CE      PLANE STRESS
        ELAST(1,1)=E/(1.D0-V*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V
        ELAST(3,3)=ELAST(1,1)*(1.D0-V)/2.D0
      ELSE
CS      RAVANSKA DEFORMACIJA
CE      PLANE STRAIN
        ELAST(1,1)=E*(1.D0-V)/(1.D0+V)/(1.D0-2.D0*V)
        ELAST(2,2)=ELAST(1,1)
        ELAST(1,2)=ELAST(1,1)*V/(1.D0-V)
        ELAST(3,3)=ELAST(1,1)*(1.D0-2.D0*V)/(1.D0-V)/2.D0
        IF(IETYP.NE.1) GO TO 40
        ELAST(4,4)=ELAST(1,1)
        ELAST(1,4)=ELAST(1,2)
        ELAST(2,4)=ELAST(1,2)
      ENDIF
   40 DO 50 I=1,6
      DO 50 J=I,6
   50 ELAST(J,I)=ELAST(I,J)
      RETURN
      END
