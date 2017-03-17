C=====================================================================
C    01.12. (GPM sa CEP1)
C        PLASTICNOST 2D ELEMENT
C
C       SUBROUTINE D2M08
C
C
C
      SUBROUTINE D2M8(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       PROGRAM ZA ODREDJIVANJE LOKACIJA VELICINA KOJE SE 
C       CUVAJU NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(4),DEF(4)
C
      IF(IDEBUG.GT.0) PRINT *, ' D2M08'
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU+4*IDVA
      LDEFPP=LDEFT+4*IDVA
      LALFMT=LDEFPP+4*IDVA
      LEMP=LALFMT+1*IDVA
      LDMI=LEMP+1*IDVA
      LALFDT=LDMI+1*IDVA      
C
      LTAU1=LPOC1
      LDEFT1=LTAU1+4*IDVA
      LDEFP1=LDEFT1+4*IDVA
      LALFM1=LDEFP1+4*IDVA
      LEMP1=LALFM1+1*IDVA
      LDMI1=LEMP1+1*IDVA
      LALFD1=LDMI1+1*IDVA
C
C
      CALL TI2408(A(LTAU),A(LDEFT),A(LDEFPP),A(LALFMT),A(LEMP),
     1             A(LDMI),A(LALFDT),A(LTAU1),A(LDEFT1),A(LDEFP1),     
     1             A(LALFM1),A(LEMP1),A(LDMI1),A(LALFD1),
     1             A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C
C======================================================================
C     
      SUBROUTINE TI2408(TAUT,DEFT,EPSP,KK1T,EPSPM,DMIT,ALF,
     1                    TAU1,DEF1,EPSP1,KK11,EPSPM1,DMI,ALF1,
     1                    FUN,MATE,TAU,DEF,IRAC)            
C**********************************************************************
C NUMERICKA INTEGRACIJA ANIZOTROPNOG MODELA ZA CIKLICNA OPTERECENJA
C**********************************************************************
C  DATUM POSLEDNJE VERZIJE 11.1998.
C**********************************************************************
C    ULAZNI PODACI:
C
C	AEM - NAGIB LCS
C	ALAM - NAGIB LNK
C	AKA - NAGIB LB
C	AE0 - POCETNI KOEFIJENT POROZNOSTI
C	G - MODUL SMICANJA
C	BKP0 - KP0 MODUL OJACANJA
C	OCR - STEPEN PREKONSOLIDACIJE
C	PIN - INICIJALNI NAPON SM
C	D - ODNOS POLUPRECNIKA A1/A0
C	ANN,AMM -NAGIBI LINIJE KRITICNOG STANJA AKO JE ELIPSA POMERENA PO P OSI
C	GAMA
C	SIGX0 - POCETNI VERTIKALNI NAPON
C	EPSD1(I) - DEVIJATORSKE KOMPONENTE UKUPNE DEFORMACIJE
C	EPSDM - ZAPREMINSKA KOMPONENTA UKUPNE DEFORMACIJE
C	ALFD(I) - DEVIJATORSKE KOORDINATE CENTRA ELIPSE TECENJA
C	ALFM - ZAPREMINSKA KOORDINATA CENTRA ELIPSE TECENJA
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD2/ TAUD(4),DEFDPR(4),DEFDS(4),DDEFP(4),
     1               DETAU(4),DDEF(4)
      COMMON /MAT2D/ EE,ANI,ET,TEQY0
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON/PLASTI/LPLAST,LPLAS1,LSIGMA
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERBR/ ITER
      COMMON /CDEBUG/ IDEBUG
      COMMON/DEVI/ S1(4),EPSD(4),EPSD1(4),EPSQJ(4),S(4),EPSQ1,EPSQS
      COMMON/KONSTANT/ ALAM,AKA,AEM,AE0,BE,G,RK,KK1      
      COMMON/PROM1/ DLAM,DLAM1,A1,AT,D,AM,AQ,ALFQ1,OCR,K1
      COMMON/MIRA/CPE(4,4),K6
C
      DIMENSION EPSP(4),SIG1(4),DELPD(4),
     *EPSP1(4),DEPSD(4),TAUT(4),DEFT(4),ALF(4),ALF1(4),
     *TAU(4),DEF(4),TAU1(4),DEF1(4),BPE(4,4),ELM(4,4),
     *EN(4),DEN(4),DEPSP(4)
      DIMENSION FUN(12,MATE),FUN2(6)
      DATA TOL1/1.D-8/
C     
      WRITE(3,1002)
1002  FORMAT(//'KOMPONENTE DEFORMACIJE I NAPONA'///)
      WRITE(3,1001)AEM,ALAM,AKA,AE0,AK0
1001  FORMAT('ULAZNI PODACI'//'M=',F5.3//'LAMDA=',F5.3//
     *'K=',F5.3,//'E0=',F5.3//'AK0=',F8.2//)   
C
      TOLER=1.0D-06
      KOR1=KOR	
C
C       INICIJALIZACIJA OSNOVNIH VELICINA
C
      AEM=FUN(1,MAT)
      ALAM=FUN(2,MAT)
      AKA=FUN(3,MAT)
      AE0=FUN(4,MAT)
      G=FUN(5,MAT)
      BKP0=FUN(6,MAT)
      AK0=FUN(7,MAT)
      D=FUN(8,MAT)
      ANN=FUN(9,MAT)
      AMM=FUN(10,MAT)
      GAMA=FUN(11,MAT)
      SIGY0=FUN(12,MAT)
C
      Z=0.0D0
      RK=ANN/AMM
      AEE=1/2./G
      AEM=AEM/(3.**0.5)
      BK0=(1+AE0)/(ALAM-AKA)
C      
C     RAVANSKO STANJE DEFORMACIJE 
C
      IF(IETYP.EQ.2)THEN
      SM0=SIGY0*(1.D0+AK0+0.3*(1.D0+AK0))/3.
      S01=SIGY0*AK0-SM0
      S02=SIGY0-SM0
      S03=Z
      S04=0.3*SIGY0*(AK0+1.0D0)-SM0                  
      ELSE
C
C     ROTACIONA SIMETRIJA
C      
      SM0=SIGY0*(1.D0+2*AK0)/3.
      S01=SIGY0*AK0-SM0
      S02=SIGY0-SM0
      S03=Z
      S04=S01
      ENDIF                        
C
      FM0=(S01*S01+S02*S02+S04*S04)/2.+S03*S03
      Q0=(FM0)**0.5      
C
      DIS=(4*SM0**2*RK**2-4*(RK**2-1)*(SM0**2+Q0**2/AEM**2))**0.5
      IF((RK-1.).LT.TOLER)THEN
      AI=(SM0**2+Q0**2/AEM**2)/(2*SM0)
      ELSE
      AI=(2*SM0-DIS)/(2*(RK**2-1.))
      ENDIF      
C-------------------------------------------------
      SMS=-(TAUT(1)+TAUT(2)+TAUT(4))/3.+SM0
      EVS=DEFT(1)+DEFT(2)+DEFT(4)
      ETT=(1.+AE0)*DEXP(EVS)-1.
      AMS=(1+ETT)*SMS/AKA
      EJ=3.0D0*AMS*(1.0D0-2.0D0*0.3D0)
C      
C      
      ANI=0.3D0
      FUN2(1)=EJ
      FUN2(2)=ANI
C
       MAT=1
      CALL MEL01(FUN2)
      IF(IRAC.EQ.2) RETURN
C
      IF(ITER.EQ.0)THEN
      DO 166 I=1,4
      TAU(I)=TAUT(I)
166   CONTINUE
      GO TO 20
      ENDIF
C      
C------------------------------------------------------------------------      
C	POTREBNE VELICINE IZ KONFIG. T  KOJE SE NE CUVAJU 
C	NA NIVOU INTEGRACIONE TACKE
C-------------------------------------------------------------------------
      DO 223 I=1,4
      DEFT(I)=-DEFT(I)
      TAUT(I)=-TAUT(I)
223   DEF(I)=-DEF(I)
      DEFT(3)=-DEFT(3)
      DEF(3)=-DEF(3)
      TAUT(3)=-TAUT(3)
C------------------------------------------------------------------------
C    RAVANSKO STANJE DEFORMACIJE I ROTACIONA SIMETRIJA (NAPONI U T) 
C------------------------------------------------------------------------
      IF(IETYP.EQ.2)THEN
      SM=(TAUT(1)+TAUT(2)+TAUT(4))/3.0D0+SM0
      S(1)=TAUT(1)-SM+SIGY0*AK0
      S(2)=TAUT(2)-SM+SIGY0
      S(3)=TAUT(3)
      S(4)=TAUT(4)-SM+0.3*SIGY0*(AK0+1.0D0)
      ELSE
      SM=(TAUT(1)+TAUT(2)+TAUT(4))/3.0D0+SM0
      S(1)=TAUT(1)-SM+SIGY0*AK0
      S(2)=TAUT(2)-SM+SIGY0
      S(3)=TAUT(3)
      S(4)=TAUT(4)-SM+SIGY0*AK0
      ENDIF      
C-----------------------------------------
      AK1=3**0.5*(1-AK0)/(1+2*AK0)
      AK0T=AK1
      Q4=AK1*SM      
C      
C	VELICINE ZA PRVI KORAK
C
      IF(KOR.EQ.1)THEN
      DMIT=Z
      DELMT=Z
      EPSPM=Z
      EPSQPT=Z
      KK1T=2
      KK1=KK1T
      DO 81 I=1,4
      DEFT(I)=Z
81    EPSP(I)=Z            
C      
      ALFM=(SM0*(D-1)+AI*RK)/D
      ALFQ=Q0*(D-1.)/D
      Q=Q0
C      
C      
      DP=(ALFQ**2/AEM**2+(ALFM-AI*RK)**2)**0.5
      XI1=AI*(ALFM-AI*RK)/(D*DP)
      YI11=DABS((AI/D)**2-XI1**2)
      YI1=AEM*YI11**0.5
      YI2=-AEM*YI11**0.5
      XI2=-AI*(ALFM-AI*RK)/(D*DP)
      YI33=DABS((AI/D)**2-XI2**2)
      YI3=AEM*YI33**0.5
      YI4=-AEM*YI33**0.5
      DELT1=(((D-1)*YI1-ALFQ)**2/AEM**2+((D-1)*XI1+AI*RK-ALFM)**2)**0.5
      DELT2=(((D-1)*YI2-ALFQ)**2/AEM**2+((D-1)*XI1+AI*RK-ALFM)**2)**0.5
      DELT3=(((D-1)*YI3-ALFQ)**2/AEM**2+((D-1)*XI2+AI*RK-ALFM)**2)**0.5
      DELT4=(((D-1)*YI4-ALFQ)**2/AEM**2+((D-1)*XI2+AI*RK-ALFM)**2)**0.5
      DELT0=MAX(DELT1,DELT2,DELT3,DELT4)
      KK1=2
      ELSE
C  
      KK1=KK1T
      ALFQ=ALF(1)
      Q=ALF(2)
      ALFM=ALF(3)
      EPSQPT=ALF(4)
      ENDIF
C------------------------------------------------------------------------      
C------------------------------------------------------------------------
      EPSM=(DEFT(1)+DEFT(2)+DEFT(4))/3.0D0
      EPSDM=(DEF(1)+DEF(2)+DEF(4))/3.0D0
      DO 21 I=1,4
      EPSD1(I)=DEF(I)-EPSDM
21    EPSD(I)=DEFT(I)-EPSM
      EPSD1(3)=DEF(3)
      EPSD(3)=DEFT(3)
C
C
      EPSVT=-3.*EPSM
      ETT=(1.+AE0)*DEXP(EPSVT)-1.
      BKT=3.0D0*(1.0D0+AE0)*EPSPM/(ALAM-AKA)
      AT=AI*DEXP(BKT)
      A0T=AT/D
      CT=RK*AT
C
C   DELT0
C   
      DETSM=(CT**2-(1+(AK1/AEM)**2)*(CT**2-AT**2))**0.5
      SMD0=(CT+DETSM)/(1+(AK1/AEM)**2)
      QD0=AK1*SMD0
      ALFM0=(SMD0*(D-1)+CT)/D
      ALFQ0=QD0*(D-1)/D
      DP=(ALFQ0**2/AEM**2+(ALFM0-CT)**2)**0.5
      XI1=AT*(ALFM0-CT)/(D*DP)
      YI11=DABS((AT/D)**2-XI1**2)
      YI1=AEM*YI11**0.5
      YI2=-AEM*YI11**0.5
      XI2=-AT*(ALFM0-CT)/(D*DP)
      YI33=DABS((AT/D)**2-XI2**2)
      YI3=AEM*YI33**0.5
      YI4=-AEM*YI33**0.5
      DELT1=(((D-1)*YI1-ALFQ0)**2/AEM**2+((D-1)*XI1+CT-ALFM0)**2)**0.5
      DELT2=(((D-1)*YI2-ALFQ0)**2/AEM**2+((D-1)*XI1+CT-ALFM0)**2)**0.5
      DELT3=(((D-1)*YI3-ALFQ0)**2/AEM**2+((D-1)*XI2+CT-ALFM0)**2)**0.5
      DELT4=(((D-1)*YI4-ALFQ0)**2/AEM**2+((D-1)*XI2+CT-ALFM0)**2)**0.5
      DELT0=MAX(DELT1,DELT2,DELT3,DELT4)
C            
C
      WRITE(3,7999)(TAUT(I),I=1,4)
7999  FORMAT('TAUT(I)',4E15.7/)      
      
C--------------------------------------------------------------------------

C	
C      
      NBIS=200
      IBISEM=50
      TOLBIS=1.0D-04
      TOLB1=1.0D-04
C      
      WRITE(3,1003) KOR           
1003  FORMAT( 'K O R A K  B R O J =',I5//)      
C
C---------------------------------------------------------      
      EPSQ1=((EPSD1(1)*EPSD1(1)+EPSD1(2)*EPSD1(2)+
     *EPSD1(4)*EPSD1(4))/2.+EPSD1(3)*EPSD1(3))**0.5
      EPSQ=((EPSD(1)*EPSD(1)+EPSD(2)*EPSD(2)+
     *EPSD(4)*EPSD(4))/2.+EPSD(3)*EPSD(3))**0.5      
      DEPSQ=EPSQ1-EPSQ
C
      EPSV=-3*EPSDM
      E1=(1.+AE0)*DEXP(EPSV)-1.
      K1=1
C
117   DEPSM=EPSDM-EPSM

      WRITE(3,4566)AT,SM,Q
4566  FORMAT(/'AT,SM,Q',3E15.7)      
C
      BKS=(1.+E1)/AKA
      BK1=AKA/(3.*(1.+E1))
      BK2=3*BK1
C                 
C     RAZLIKA DEVIJATORA DEFORMACIJE U (T+DT) I DEVIJATORA PLASTICNE
C     DEFORMACIJE U T
C
      DO 3 I=1,4
3     DEPSD(I)=EPSD1(I)-EPSD(I)
      DEPSD(3)=DEPSD(3)/2.
      DEPSQ=EPSQ1-EPSQ
      EPSQS=DEPSQ+Q*AEE
C      
C###################################################################
C     ODREDJIVANJE NAPONA U T+DT
C###################################################################
C
      IM=0
      IP=0
      TM=0.0D0
      TP=0.0D0
      TFM=0.0D0
      TFP=0.0D0
      IBISE=0
      BE=SM*BKS
      OCR=2.
      K4=0
      K6=0
      K5=0
      K3=0
C
C    ELASTICNO PREDVIDJANJE
C      
      SM1=SM+3*BE*DEPSM
      Q1=Q+DEPSQ/AEE
      IF(K1.EQ.1)GO TO 88
C--------------------------------------------------------------------    
4     IBISE=IBISE+1
C
C--------------------------------------------------------
C	OCR=1
C--------------------------------------------------------      
      IF((DABS(DELTT0).LE.1.0D-05).AND.(KK1.NE.1))THEN
C      
      OCR=1.0D0
      IF(IBISE.EQ.1)THEN
C
      GI1=BKPT*(XT**2+YT**2/AEM**4)**0.5
      G1=1./GI1
      SK0=(XT**2+YT**2/4./AEM**4)**0.5
      SK1=(XT**2+YT**2/AEM**4)**0.5
      G2=G1/(2*SK0)
      DSME=BE*3*DEPSM
      DQE=2*G*DEPSQ
      DELP=G1*(XT*DSME+YT*DQE/AEM**2)/
     *(1+2*G2*(BE*XT**2+G*YT**2/AEM**4))
      DELP=DABS(DELP)
C
      DELVP=DELP*XT/SK0
      DELM=DELVP/3.
      DELQP=DELP*YT/AEM**2/(2*SK0)
C
      SM10=DSME+SM-BE*DELVP
      Q10=DQE+Q-2*G*DELQP
C
      C10=1+(D-1)*BK0*DELVP
      DALFM0=((D-1)*BK0*DELVP*(SM10-ALFM)+AT*BK0*DELVP)/C10
      DALFQ0=(D-1)*BK0*DELVP*(Q10-ALFQ)/C10
      SMR=D*(SM10-ALFM-DALFM0)+AT*RK
      QR=D*(Q10-ALFQ-DALFQ0)
      BETAM=SMR-SM10
      BETAQ=QR-Q10
      DELTT=(BETAQ**2/AEM**2+BETAM**2)**0.5
      DELTT0=DELTT/DELT0
      K6=1
      IF(DABS(DELTT0).GT.1.0D-08)GO TO 120
C      
      SKT=(4*XT**2+YT**2/AEM**2)**0.5
      DELPB=XT*BE*3*DEPSM+YT*2*G*DEPSQ/AEM**2
      DELPI=2*XT**2*(BE+(D-1)*BK0*XT+AT*BK0*RK)/SKT+
     *2*YT*(G*YT/AEM**2+(D-1)*BK0*YT*XT)/AEM**2/SKT+
     *2*A0T**2*BK0*XT/SKT 
      DELP=DELPB/DELPI
      DELP=DABS(DELP)
      DX=DELP/10.
C      
      DELP0=Z
      EPSVP0=3*EPSPM
      EPSQP0=EPSQPT
C
      ALFM0=ALFM
      ALFQ0=ALFQ      
C
      SM10=SME1
      Q10=QE1
C
      IB=0
      XP=Z
      FP=FE
C      
      WRITE(3,347)XP,FP
347   FORMAT( 'XP0,FP0',2E15.7/)      
C--------------------------------------------------------------------
C       ZA CEP1
C-------------------------------------------------      
      S10=2*G*DEPSD(1)+S(1)      
      S20=2*G*DEPSD(2)+S(2)      
      S30=2*G*DEPSD(4)+S(4)      
      S40=2*G*DEPSD(3)+S(3)      
C
      YI0=Q10-ALFQ0
      XI0=SM10-ALFM0
C
      IF(DABS(Q10-Z).LT.TOL1)THEN
      EN01=2*XI0/3.
      EN02=2*XI0/3.
      EN03=2*XI0/3.
      EN04=Z
      SK00=(XI0**2+YI0**2/4./AEM**4)**0.5
      FS0=2*SK00
C      FS0=(EN01**2+EN02**2+EN03**2+2*EN04**2)**0.5
      EN01=EN01/FS0
      EN02=EN02/FS0
      EN03=EN03/FS0
      GO TO 77
      ENDIF
C      
      EN01=(YI0*(2*S10-S20-S30)/Q10/AEM**2+2*XI0)/3.
      EN02=(YI0*(-S10+2*S20-S30)/Q10/AEM**2+2*XI0)/3.
      EN03=(YI0*(-S10-S20+2*S30)/Q10/AEM**2+2*XI0)/3.
      EN04=YI0*S40/Q10/AEM**2
      SK00=(XI0**2+YI0**2/4./AEM**4)**0.5
      FS0=2*SK00
C      FS0=(EN01**2+EN02**2+EN03**2+2*EN04**2)**0.5
      EN01=EN01/FS0
      EN02=EN02/FS0
      EN03=EN03/FS0
      EN04=EN04/FS0
C
77    SK0X0=XI0/SK00
      SK0Y0=YI0/AEM**2./SK00      
      SK10=(XI0**2+YI0**2/AEM**4)**0.5
      SK1X0=XI0/SK10
      SK1Y0=YI0/AEM**2./SK10
C
      ALFC0=SK10/SK00            
C            
      ENDIF
C----------------------------------------------------------
      WRITE(3,102)IBISE
102   FORMAT( 'IBISE=',I5)      
C
      X0=SM10-ALFM0
      Y0=Q10-ALFQ0
      SK0=(X0**2+Y0**2/(4.*AEM**4))**0.5
C---------------------------------------------------------
C      
C    PLASTICNE DEFORMACIJE U T+DT
C---------------------------------------------------------      
      EPSVP=EPSVP0+(DELP-DELP0)*X0/SK0
      EPSQP=EPSQP0+(DELP-DELP0)*Y0/(2.*AEM**2*SK0)      
C---------------------------------------------------
C   NAPONI U T+DT
C---------------------------------------------------
C
      SM1=SM10-(DELP-DELP0)*BE*X0/SK0
      Q1=Q10-(DELP-DELP0)*G*Y0/AEM**2/SK0
      DELP2=(XT*(SM1-SM)+YT*(Q1-Q)/AEM**2)/
     *(BKPT*(XT**2+YT**2/AEM**4)**0.5)
C      
C
      WRITE(3,4276)SM1,Q1,DELP2,KK1,DELTT0
4276  FORMAT(/'SM1,Q1,DELP2',3E15.7/'KK1,DELTT0',I5,E15.7/)   
C      
C-----------------------------------------------------------
C	UNUTRASNJE PROMENLJIVE ALFM1,ALFQ1
C-----------------------------------------------------------
C
      SK1=(X0**2+Y0**2/AEM**4)**0.5
      ALFC=SK1/SK0
      DALFM=-(DELP-DELP0)*ALFC*BK0*(AT*RK+(D-1)*XT)*X0/SK1
      DALFQ=-(DELP-DELP0)*ALFC*AEM**2*(D-1)*BK0*XT*Y0/(AEM**2*SK1)
      ALFM1=ALFM0-DALFM
      ALFQ1=ALFQ0-DALFQ
C
      X=SM1-ALFM1
      Y=Q1-ALFQ1
C
      DLAM=DELP/((4*X**2+Y**2/AEM**4)**0.5)
      DLI=BKPT*((X**2+Y**2/AEM**4)*(4*X**2+Y**2/AEM**4))**0.5
      DLAMP=(X*(SM1-SM)+Y*(Q1-Q)/AEM**2)/DLI
      DELM=DLAM*2*X/3.
      DELQP=DLAM*Y/AEM**2      
C
      EPSPV1=3.*(EPSPM+DELM)
      BB1=BK0*EPSPV1
      A1=AI*DEXP(BB1)      
      C=A1*RK
      A0=A1/D
      C1=1+(D-1)*BK0*3*DELM      
C
      F=X**2.+Y**2./AEM**2.-A0**2.   
C
      WRITE(3,163)F,DELP,DLAM,DLAMP,OCR,SM1,Q1,X,Y
163   FORMAT( 'F,DELP',2E15.7/'DLAM,DLAMP,OCR',3E15.7/
     *'SM1,Q1,X,Y',4E15.6//)      
C
      IF(DABS(F).LE.TOLBIS)GO TO 90
C
      IF(IBISE.EQ.IBISEM)THEN
      WRITE(3,159)
159   FORMAT( 'MAX. BROJ BISEKCIJA')      
      STOP
      ENDIF
C
      DELP0=DELP
      EPSVP0=EPSVP
      EPSQP0=EPSQP
      SM10=SM1
      Q10=Q1
      ALFM0=ALFM1
      ALFQ0=ALFQ1
C      
       AF=1.0
      CALL BISEC1(DELP,XM,XP,DX,F,FM,FP,AF,IB)
      GO TO 4
C           
      ENDIF
C      
C
C========================================================
C      OCR>1
C---------------------------------------------------------
C      
120   IF(IBISE.EQ.1)THEN
      Q10=QE1
      SM10=SME1
      ALFM0=ALFM
      ALFQ0=ALFQ
      X0=SM10-ALFM0
      Y0=Q10-ALFQ0
      A00=A0T
      F0=FE
      DELPD0=Z

C
      DSME=BE*3*DEPSM
      DQE=2*G*DEPSQ
C--------------------------------------------------------------------
C       ZA CEP1
C-------------------------------------------------      
      S10=2*G*DEPSD(1)+S(1)      
      S20=2*G*DEPSD(2)+S(2)      
      S30=2*G*DEPSD(4)+S(4)      
      S40=2*G*DEPSD(3)+S(3)      
C
      YI0=Q-ALFQ0
      XI0=SM-ALFM0
C
      IF(DABS(Q-Z).LT.TOL1)THEN
      EN01=2*XI0/3.
      EN02=2*XI0/3.
      EN03=2*XI0/3.
      EN04=Z
      SK00=(XI0**2+YI0**2/4./AEM**4)**0.5
      FS0=2*SK00
C      FS0=(EN01**2+EN02**2+EN03**2+2*EN04**2)**0.5
      EN01=EN01/FS0
      EN02=EN02/FS0
      EN03=EN03/FS0
      GO TO 11
      ENDIF
C      
      EN01=(YI0*(2*S10-S20-S30)/Q/AEM**2+2*XI0)/3.
      EN02=(YI0*(-S10+2*S20-S30)/Q/AEM**2+2*XI0)/3.
      EN03=(YI0*(-S10-S20+2*S30)/Q/AEM**2+2*XI0)/3.
      EN04=YI0*S40/Q/AEM**2
      SK00=(XI0**2+YI0**2/4./AEM**4)**0.5
      FS0=2*SK00
C      FS0=(EN01**2+EN02**2+EN03**2+2*EN04**2)**0.5
      EN01=EN01/FS0
      EN02=EN02/FS0
      EN03=EN03/FS0
      EN04=EN04/FS0
C
11    SK0X0=XI0/SK00
      SK0Y0=YI0/AEM**2./SK00      
      SK10=(XI0**2+YI0**2/AEM**4)**0.5
      SK1X0=XI0/SK10
      SK1Y0=YI0/AEM**2./SK10
C
      ALFC0=SK10/SK00            
C            
      ENDIF
C--------------------------------------------------------      
C    PLASTICNE DEFORMACIJE U T+DT
C---------------------------------------------------------      
C      
      GI1=BKPT*(X0**2+Y0**2/AEM**4)**0.5
      G1=1./GI1
      SK0=(X0**2+Y0**2/4./AEM**4)**0.5
      G2=G1/(2*SK0)
C
      DELP=G1*(X0*DSME+Y0*DQE/AEM**2)/
     *(1+2*G2*(BE*X0**2+G*Y0**2/AEM**4))
      DELP=DABS(DELP)
C      
      SK1=(X0**2+Y0**2/AEM**4)**0.5
      ALFC=SK1/SK0
C
      DELVP=DELP*X0/SK0
      DELM=DELVP/3.
      DELQP=DELP*Y0/AEM**2/(2*SK0)
      EPSQP=EPSQPT+DELQP     
C---------------------------------------------------------
C	NAPONI U T+DT
C---------------------------------------------------------
C
      SM1=DSME+SM-BE*DELVP
      Q1=DQE+Q-2*G*DELQP
C---------------------------------------------------------
C	UNUTRASNJE PROMENLJIVE
C---------------------------------------------------------
C
      C1=1.+(D-1)*BK0*3*DELM
      EPSPV1=3.*(EPSPM+DELM)
      BB1=BK0*EPSPV1
      A1=AI*DEXP(BB1)      
      C=A1*RK
      A0=A1/D
      DA0=A0-A00
      DA1=A1-D*A00
      BETAM=D*A00*RK-SM10+D*X0
      BETAQ=-Q10+D*Y0      
C      
      AM=A1*RK+SM10*(D-1)-D*ALFM0
      AQ=Q10*(D-1)-D*ALFQ0
C
      FX0=2*X0*(SM1-SM10-(D-1)*DA0*X0/A00-DA1*RK)
      FY0=2*Y0*(Q1-Q10-(D-1)*DA0*Y0/A00)/AEM**2
      DMIB=F0+FX0+FY0-2*A00*DA0
      DMII=2*(X0*BETAM+Y0*BETAQ/AEM**2)
      DMI=DMIB/DMII
C
      WRITE(3,913)DMI,DMIB,DMII,BETAM,BETAQ
913   FORMAT( 'DMI,DMIB,DMII',3E15.7/'BETAM,BETAQ',2E15.7/)
C
      DALFM=BETAM*DMI+(D-1)*DA0*X0/A00+DA1*RK
      DALFQ=BETAQ*DMI+(D-1)*DA0*Y0/A00
      ALFM1=ALFM0+DALFM
      ALFQ1=ALFQ0+DALFQ
      X=SM1-ALFM1
      Y=Q1-ALFQ1
      DLAM=DELP/((4*X**2+Y**2/AEM**4)**0.5)
      GI1=BKPT*(X**2+Y**2/AEM**4)**0.5
      G1=1./GI1
      SKK=(4*X**2+Y**2/AEM**4)**0.5
      G21=G1/SKK
      DLAMP=G21*(X*(SM1-SM)+Y*(Q1-Q)/AEM**2)
      FY=X**2+Y**2/AEM**2-A0**2
      F=FY
C
      WRITE(3,911)FY,DELP,DMI,SM1,Q1,X,Y
911   FORMAT( 'FY,DELP,DMI',3E15.7/'SM1,Q1,X,Y',4E15.7//)
C-----------------------------------------------------      
C
      WRITE(3,102)IBISE
C
      DDELP=DABS(DELP-DELPD0)/DELP
      F4=DELP*BKPT-
     *(X*(SM1-SM)-Y*(Q1-Q)/AEM**2)/(X**2+Y**2/AEM**4)**0.5
      WRITE(3,912)F4
912   FORMAT( 'F4',E15.7/)
C
      IF(DABS(FY).LE.TOLBIS)GO TO 90
C
      IF(IBISE.EQ.IBISEM)THEN
      WRITE(3,159)
      STOP
      ENDIF
C
      SM10=SM1
      Q10=Q1
      ALFM0=ALFM1
      ALFQ0=ALFQ1
      X0=X
      Y0=Y
      A00=A0
      F0=FY
      DELPD0=DELP
      GO TO 4
C
C-----------------------------------------------------
C	ELASTICNO RESENJE
C=====================================================
C
88    SME1=SM1
      QE1=Q1      
      FE=(SME1-ALFM)**2+(QE1-ALFQ)**2/AEM**2-(AT/D)**2
C
      XT=SM-ALFM
      YT=Q-ALFQ
C
      IF((FE.LT.Z).AND.(K1.EQ.1)) GO TO 303
      IF((KK1.EQ.1).AND.(K1.EQ.1)) GO TO 113      
C
C    POCETNA VREDNOST ZA DELM
C
      SMRT=(SM-ALFM)*D+CT
      QRT=(Q-ALFQ)*D
      GFFT=((SMRT-CT)**2+QRT**2/(AEM**4))**0.5
      BKPRT=AT*BK0*((SMRT-CT)*RK+AT)*(SMRT-CT)/GFFT**2
      DELTT=(((QRT-Q)/AEM)**2+(SMRT-SM)**2)**0.5
      BKPT=BKPRT+BKP0*(DELTT/DELT0)**(GAMA+1)      
      DELTT0=DELTT/DELT0
      WRITE(3,1696)DELTT0,BKPT
1696  FORMAT( 'DELTT0,BKPT',2E15.7/)      
C      
      G1T=BKPT*((XT**2+YT**2/AEM**4)*(4*XT**2+YT**2/AEM**4))**0.5
      G1=1/G1T
      DQA=G*YT/(XT*AEM**2*BE)      
      DQB=2*G*DEPSQ-3*DEPSM*G*YT/(XT*AEM**2)
      WRITE(3,1697)XT,YT
1697  FORMAT( 'XT,YT',2E15.7/)            
C      
      G3=1./BE+2*XT**2*G1+2*YT*XT*G1*DQA/AEM**2
      WRITE(3,1698)G3
1698  FORMAT( 'G3',E15.7/)                  
      DSM0B=3*DEPSM-2*XT*YT*G1*DQB/AEM**2
      DSM0=DSM0B/G3
      DQ0=DQA*DSM0+DQB

      DELM=2*G1*XT*(XT*DSM0+YT*DQ0/(AEM**2))/3.  
      EPSPV1=3.*(EPSPM+DELM)
      BB1=BK0*EPSPV1
      A1=AI*DEXP(BB1)      
      A0=A1/D      
      FE1=(SM+DSM0-ALFM)**2+(Q+DQ0-ALFQ)**2/AEM**2-A0**2
C
      DELM00=DELM
C      

C      
      DLAMF=G1*(XT*DSM0+YT*DQ0/(AEM**2))
      DLAM0=DLAMF
C      
      DELQP=G1*YT*(XT*DSM0+YT*DQ0/(AEM**2))/AEM**2  

      DELMP1=DEPSM-DSM0/BE/3.
C      
      WRITE(3,1691)FE,FE1,DELM,DELMP1,SME1,QE1,DEPSM,DEPSQ
1691  FORMAT( 'FE,FE1,DELM,DELMP1',4E15.7/'SME1,QE1',2E15.7/
     *'DEPSM,DEPSQ',2E15.7//)

C
      K1=2
      GO TO 4
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
90    DMIE=((ALFM1-ALFM)-(D-1)*BK0*DELVP*X-(A1-AT)*RK)/
     *(C-SM1+D*X)
      WRITE(3,161)DMIE,DMI,DELP,DLAM,DLAMP,OCR
161   FORMAT( 'DMIE,DMI,DELP',3E15.7/'DLAM,DLAMP,OCR',3E15.7/)
C
      DO 327 I=1,4
      S1(I)=Q1*(DEPSD(I)+S(I)*AEE)/(Q1*AEE+DLAM*Y/(AEM**2*C1))
327   CONTINUE
C
C------------------------------------------------------------------
C	ZA CEP1
C------------------------------------------------------------------
      IF(DABS(Q1-Z).LT.TOL1)THEN
      EN(1)=2*X/3.
      EN(2)=2*X/3.
      EN(3)=2*X/3.
      EN(4)=Z
      FS=2*SK0
C      FS=(EN(1)**2+EN(2)**2+EN(3)**2+2*EN(4)**2)**0.5
      EN(1)=EN(1)/FS
      EN(2)=EN(2)/FS
      EN(3)=EN(3)/FS
      GO TO 71
      ENDIF
C      
      EN(1)=(Y*(2*S1(1)-S1(2)-S1(4))/(C1*Q1)/AEM**2+2*X)/3.
      EN(2)=(Y*(-S1(1)+2*S1(2)-S1(4))/(C1*Q1)/AEM**2+2*X)/3.
      EN(3)=(Y*(-S1(1)-S1(2)+2*S1(4))/(C1*Q1)/AEM**2+2*X)/3.
      EN(4)=Y*S1(3)/(C1*Q1)/AEM**2
      FS=2*SK0
      FS1=(EN(1)**2+EN(2)**2+EN(3)**2+2*EN(4)**2)**0.5
      EN(1)=EN(1)/FS
      EN(2)=EN(2)/FS
      EN(3)=EN(3)/FS      
      EN(4)=EN(4)/FS
C 
71    DELPR=DELP
C
C      
      DEN(1)=(EN(1)-EN01)/DELPR
      DEN(2)=(EN(2)-EN02)/DELPR
      DEN(3)=(EN(3)-EN03)/DELPR
      DEN(4)=(EN(4)-EN04)/DELPR
C           
      SK0X1=X/SK0
      DSK0X=(SK0X1-SK0X0)/DELPR
      SK0Y1=Y/AEM**2/SK0
      DSK0Y=(SK0Y1-SK0Y0)/DELPR
      SK1X1=X/SK1
      DSK1X=(SK1X1-SK1X0)/DELPR
      SK1Y1=Y/AEM**2/SK1
      DSK1Y=(SK1Y1-SK1Y0)/DELPR
      DALFC=(ALFC-ALFC0)/DELPR
C------------------------------------------------------------------
C      
      GO TO 301
C      
C*******************************************************************
C    RESENJE ZA SLUCAJ DA JE STANJE U T BILO ELASTICNO, A U (T+DT)
C       ELASTOPLASTICNO
C*******************************************************************
113   IF(DABS(DEPSQ-Z).LT.TOL1)THEN
      SMP1=(A0T**2-YT**2/AEM**2)**0.5+ALFM
      SMP2=-(A0T**2-YT**2/AEM**2)**0.5+ALFM
      D1SM=SM-SMP1
      D2SM=SM-SMP2
      DD1SM=DABS(D1SM)
      DD2SM=DABS(D2SM)
      DSM=MIN(DD1SM,DD2SM)
      IF(DABS(DD1SM-DSM).LT.TOL1)THEN
      SMPRIM=SMP1
      ELSE
      SMPRIM=SMP2
      ENDIF
      QPRIM=Q
      DO 903 I=1,4
903   S1(I)=S(I)
      DEPSME=(SMPRIM-SM)/(SM*BKS)
      GO TO 172
      ENDIF
C
      IF(DABS(DEPSM-Z).LT.TOL1)THEN
      QP1=AEM*(A0T**2-XT**2)**0.5+ALFQ
      QP2=-AEM*(A0T**2-XT**2)**0.5+ALFQ
      D1Q=Q-QP1
      D2Q=Q-QP2
      DD1Q=DABS(D1Q)
      DD2Q=DABS(D2Q)
      DQP=MIN(DD1Q,DD2Q)
      IF(DABS(DD1Q-DQP).LT.TOL1)THEN
      QPRIM=QP1
      ELSE
      QPRIM=QP2
      ENDIF      
      SMPRIM=SM
      DEPSME=Z
      GO TO 219
      ENDIF
C      
      D1=3*DEPSM/DEPSQ
      D2=2*BK2*G/SM/D1
      AA=1+(D2/AEM)**2
      BB=2*(SM-ALFM)+2*D2*(Q-ALFQ)/AEM**2
      CC=(SM-ALFM)**2+(Q-ALFQ)**2/AEM**2-(AT/D)**2
      DSM1=(-BB+(BB**2-4*AA*CC)**0.5)/(2*AA)
      DSM1A=DABS(DSM1)
      DSM2=(-BB-(BB**2-4*AA*CC)**0.5)/(2*AA)
      DSM2A=DABS(DSM2)
      DSM=MIN(DSM1A,DSM2A)
C      
      IF(DABS(DSM1A-DSM).LT.TOL1)THEN
      SMPRIM=SM+DSM1
      QPRIM=Q+DSM1*D2
      DEPSME=DSM1*BK2/SM/3.
      ELSE
      SMPRIM=SM+DSM2
      QPRIM=Q+DSM2*D2
      DEPSME=DSM2*BK2/SM/3.
      ENDIF 
C
219   CONTINUE
      DO 211 I=1,4
      S1(I)=S(I)+DEPSD(I)*(QPRIM-Q)/DEPSQ
211   EPSD(I)=EPSD(I)+(S1(I)-S(I))*AEE
172   FEP=(SMPRIM-ALFM)**2+(QPRIM-ALFQ)**2/AEM**2-(AT/D)**2
      
C
C      
      EPSM=EPSM+DEPSME
C
      DO 213 I=1,4
213   S(I)=S1(I)
      EPSQ=EPSQ+(QPRIM-Q)*AEE
      DEPSQ=EPSQ1-EPSQ
C        
C
      Q=QPRIM
      SM=SMPRIM      
      KK1=2
C
C
      GO TO 117 
C*********************************************************************      
C
C
C     KONACNI NAPONI I DEFORMACIJE
C
C     NAPONI
C-----------------------------------------------------------
C     ELASTICNO STANJE
C-----------------------------------------------------------
C	
303   KK1=1
      DELM=Z
      ALFM1=ALFM
      ALFQ1=ALFQ
      EPSPV1=3.*(EPSPM+DELM)
      BB1=BK0*EPSPV1
      A1=AI*DEXP(BB1)      
      A0=A1/D      
C      
      DO 310 I=1,4
      DELPD(I)=Z
310   S1(I)=2*G*DEPSD(I)+S(I)      
C
      GO TO 305
C-------------------------------------------------------------------
C     ELASTO-PLASTICNO STANJE
C-------------------------------------------------------------------
C      
301   KK1=2
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(3,1004)IBISE    
1004  FORMAT(/'BROJ BISEKCIJA =',I5//)    
C
C
305   SIG1(1)=S1(1)+SM1
      SIG1(2)=S1(2)+SM1
      SIG1(3)=S1(3)
      SIG1(4)=S1(4)+SM1 
C
C
C     DEFORMACIJE
C
      EPSPM1=EPSPM+DELM
C
      IF(KK1.EQ.1)GO TO 309      
C      
      DO 31 I=1,4
      DELPD(I)=DEPSD(I)-(S1(I)-S(I))*AEE
31    CONTINUE      
C
309   DO 32 I=1,4
32    EPSP1(I)=EPSP(I)+DELPD(I)
C      
      DELQP=DEPSQ-(Q1-Q)*AEE
C      
      IF(KK1.EQ.1)THEN
      EPSQP=EPSQPT+DELQP
      ENDIF
C      
      EPSQP1=EPSQP
C      
      EPSQ2=EPSQ1*2./3.**0.5
      EPSV1=3*EPSDM
      WRITE(3,1006)
1006  FORMAT( 'DEFORMACIJE:'//'UKUPNE'//)
      WRITE(3,1007)(EPSD1(I),I=1,4)
1007  FORMAT( 'EPSD1',4E15.7//)                  
      WRITE(3,1307)EPSV1,EPSQ2,EPSDM
1307  FORMAT( 'EPSV1,EPSQ2,EPSDM',3E15.7//)                  
      WRITE(3,1008)
1008  FORMAT( 'PLASTICNE'//)
      WRITE(3,1009)EPSPM1,(EPSP1(I),I=1,4)
1009  FORMAT( 'ZAPREMINSKA',E15.7,/'DEVIJATORSKE',4E15.7//)      
      WRITE(3,1010)
1010  FORMAT( 'NAPONI:'//)
C

      SIGX1=S1(1)+SM1
      DSIGX=S1(1)+SM1-S(1)-SM
      SIGY1=S1(2)+SM1
      DSIGY=S1(2)+SM1-S(2)-SM
      SIGZ1=S1(4)+SM1
      DSIGZ=S1(4)+SM1-S(4)-SM
C
      Q2=Q1*3.**0.5
      D2I=(S1(1)**2+S1(2)**2+S1(4)**2+2*S1(3)**2)/2.
      Q5=D2I**0.5
      DELE1=DEPSD(1)+DEPSM
      DELE2=DEPSD(2)+DEPSM
      DELE3=DEPSD(4)+DEPSM
      DELE4=DEPSD(3)
C      
      WRITE(3,1011)SM1,(S1(I),I=1,4),Q1,Q2,Q5,SIGX1,SIGY1,SIGZ1,EPSQP1
1011  FORMAT( 'SFERNA KOMP.',E15.7,/'DEVIJATORSKE',4E15.7/'Q1,Q2,Q3 = ',
     *3E15.7/'SIGX1,SIGY1,SIGZ1',3E15.7//'EPSQP1',E15.7//)
C     
      WRITE(3,1012)A1,A0
1012  FORMAT('POLUPRECNICI A1 I A0'//'POLUPRECNIK A1 = ',E15.7/
     *'POLUPRECNIK A0 = ',E15.7//)      
      WRITE(3,4356) ALFM1,ALFQ1
4356  FORMAT( 'ALFM1,ALFQ1',2E15.7)      
      WRITE(3,4388) E1,AMS,BE
4388  FORMAT( 'E,AMS,BE',3E15.7//)
C
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C	VELICINE KOJE SE VRACAJU U GLAVNI PROGRAM
C----------------------------------------------------------------
C
      KK11=KK1
      ALF1(1)=ALFQ1
      ALF1(2)=Q1
      ALF1(3)=ALFM1
      ALF1(4)=EPSQP1
C      
      IF(IETYP.EQ.2)THEN
      SIGX0 = SIGY0*AK0
      TAU(1)=-(SIG1(1)-SIGX0)
      TAU(2)=-(SIG1(2)-SIGY0)
      TAU(3)=SIG1(3)
      TAU(4)=-(SIG1(4)-0.3*SIGY0*(AK0+1.0D0))
      ELSE
      SIGX0 = SIGY0*AK0
      TAU(1)=-(SIG1(1)-SIGX0)
      TAU(2)=-(SIG1(2)-SIGY0)
      TAU(3)=SIG1(3)
      TAU(4)=-(SIG1(4)-SIGX0)
      ENDIF      
C      
      DO 75 I=1,4
      DEF(I)=-DEF(I)
      DEF1(I)=DEF(I)
75    TAU1(I)=TAU(I)
      DEF(3)=-DEF(3)
      DEF1(3)=DEF(3)
C
      IF(KK1.EQ.1)GO TO 73
C------------------------------------------------------------
C	ZA CEP1
C------------------------------------------------------------      
C
      BF1=-2*X*BE*(X/SK0+DELP*DSK0X)-2*Y*G*(Y/AEM**2/SK0+DELP*DSK0Y)/
     *AEM**2
      BF2=2*X*BK0*((D-1)*XT+AT*RK)*(ALFC*(X/SK1+DELP*DSK1X)+
     *DELP*DALFC*X/SK1)+
     *2*Y*BK0*(D-1)*XT*(ALFC*(Y/AEM**2/SK1+DELP*DSK1Y)+
     *DELP*DALFC*Y/AEM**2/SK1)
      BF2=2*BK0*(X/SK0+DELP*DSK0X)*(X*(X*(D-1)+A1*RK)+
     *Y**2*(D-1)/AEM**2)
      BF=BF1-BF2
      BF=-1/BF
C
      WRITE(3,852)BF1,BF2,BF
852   FORMAT( 'AF1,AF2,AF',3E15.7/)      
      IF(DABS(Q1-Z).LT.TOL1)THEN
      DEPSP(1)=BF*2*BE*X      
      DEPSP(2)=BF*2*BE*X      
      DEPSP(3)=BF*2*BE*X      
      DEPSP(4)=Z
      GO TO 72
      ENDIF
C            
      DEPSP(1)=BF*2*(BE*X+G*Y*(2*S1(1)-S1(2)-S1(4))/(3*Q1*AEM**2))      
      DEPSP(2)=BF*2*(BE*X+G*Y*(-S1(1)+2*S1(2)-S1(4))/(3*Q1*AEM**2))      
      DEPSP(3)=BF*2*(BE*X+G*Y*(-S1(1)-S1(2)+2*S1(4))/(3*Q1*AEM**2))      
      DEPSP(4)=BF*2*G*Y*S1(3)/AEM**2/Q1
72    CONTINUE
C
      DO 78 I=1,4
      DO 78 J=1,4
      ELM(I,J)=Z
      CPE(I,J)=Z
78    BPE(I,J)=Z
C
      DO 80 I=1,3
      DO 80 J=1,3
      ELM(I,J)=BE-2*G/3.
      IF(I.EQ.J)THEN
      ELM(I,J)=BE+4*G/3.
      ENDIF
80    CONTINUE
      ELM(4,4)=2*G
C
      DO 70 I=1,4
      DO 70 J=1,4
      BPE(I,J)=(EN(I)+DELP*DEN(I))*DEPSP(J)
      IF(I.EQ.J)THEN
      BPE(I,J)=1.-BPE(I,J)
      ELSE
      BPE(I,J)=-BPE(I,J)
      ENDIF
70    CONTINUE      
C
      DO 79 I=1,4
      DO 79 J=1,4
      CPE(I,J)=ELM(I,1)*BPE(1,J)+ELM(I,2)*BPE(2,J)+
     *ELM(I,3)*BPE(3,J)+ELM(I,4)*BPE(4,J)    
79    CONTINUE     
      DO 93 I=1,4
93    WRITE(3,555)(CPE(I,J),J=1,4)
555   FORMAT( 'CEP1'/4E15.6//)
C
      C11=(1.-BF*2*X*BE*(X/SK0+DELP*DSK0X))*BE
      C12=-BF*4*G*Y*(X/SK0+DELP*DSK0X)*BE/AEM**2
      C21=-BF*2*X*BE*(Y/SK0/AEM**2+DELP*DSK0Y)*G
      C22=(1.-BF*2*G*Y*(Y/SK0/AEM**2+DELP*DSK0Y)/AEM**2)*2*G
      DSM3=SM1-SM
      DQ3=Q1-Q
      DSM4=C11*3*DEPSM+C12*DEPSQ
      DQ4=C21*3*DEPSM+C22*DEPSQ
      WRITE(3,556)SK0,FS1,DSM3,DSM4,DQ3,DQ4
556   FORMAT( 'SK0,FS',2E15.6/'DSM1,DSM2',2E15.6/'DQ1,DQ2',2E15.6//)      
C
      DSIG1=CPE(1,1)*DELE1+CPE(1,2)*DELE2+CPE(1,3)*DELE3+
     *CPE(1,4)*DELE4
      DSIG2=CPE(2,1)*DELE1+CPE(2,2)*DELE2+CPE(2,3)*DELE3+
     *CPE(2,4)*DELE4
      DSIG4=CPE(3,1)*DELE1+CPE(3,2)*DELE2+CPE(3,3)*DELE3+
     *CPE(3,4)*DELE4
      DSM5=(DSIG1+DSIG2+DSIG4)/3.
      WRITE(3,459)DSIG1,DSIG2,DSIG4,DSIGX,DSIGY,DSIGZ,DSM5
C--------------------------------------------------------------------        
73    CONTINUE      
C	MATRICA CEP  	L2=2 ZA CE, L2=1 ZA CEP
C      IF(ISKNP.NE.2)THEN
C      CALL CEP(DELM,SM1,Q1,L2)
C      ENDIF     
C      CALL CEP(DELM,SM1,Q1,L2)
C
C
C      WRITE(3,459)DSIG1,DSIG2,DSIG4,DSIGX,DSIGY,DSIGZ,DSM5
459   FORMAT( 'DSIG1,2,3',3E15.7/'DSIGX,Y,Z',3E15.7/'DSM5',E15.7//)      
      WRITE(3,460)DLAM,DLAM1
460   FORMAT( 'DLAM,DLAM1',2E15.7/)      
C
C
20    CONTINUE
      RETURN
      END
C
C==============================================================
C	ELASTO-PLASTICNA  MATRICA
C==============================================================
      SUBROUTINE CEP(DELM,SM1,Q1,L2)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/DEVI/S1(4),EPSD(4),EPSD1(4),EPSQJ(4),S(4),EPSQ1,EPSQS
      COMMON/KONSTANT/ALAM,AKA,AEM,AE0,BE,G,RK,KK1      
      COMMON/PROM1/ DLAM,DLAM1,A1,AT,D,AM,AQ,ALFQ1,OCR,K1
      COMMON/MIRA/CPE(4,4),K6
      DATA TOL1/1.D-8/
C
C      
      S3=S1(4)
      S4=S1(3)
      ST3=S(4)
      ST4=S(3)
      EPSD3=EPSD1(4)
      EPS3=EPSD(4)
      EPSD4=EPSD1(3)
      EPS4=EPSD(3)
      S1(3)=S3
      S1(4)=S4
      S(3)=ST3
      S(4)=ST4
      EPSD1(3)=EPSD3
      EPSD1(4)=EPSD4
      EPSD(3)=EPS3
      EPSD(4)=EPS4
      Z=0.0D0
      DLAM1=DLAM
C
C	INICIJALIZACIJA CEP
C      
      DO 1 I=1,4
      DO 1 J=1,4
1     CPE(I,J)=Z      
C
      
C----------------------------------------------------------------
C              ELASTICNA MATRICA
C----------------------------------------------------------------      
C
      IF((KK1.EQ.1).OR.(DABS(Q1-Z).LT.TOL1).OR.(DABS(EPSQ1-Z)).LT.TOL1)
     +          THEN
      DO 4 I=1,3
      DO 4 J=1,3
      IF(I.EQ.J)THEN
      CPE(I,J)=4*G/3.+BE
      ELSE
      CPE(I,J)=-2*G/3.+BE
      ENDIF
4     CONTINUE
      CPE(4,4)=2*G
      GO TO 20
      ENDIF
C
C---------------------------------------------------------------------
C            ELASTOPLASTICNA MATRICA
C----------------------------------------------------------------------
C
      BK0=(1+AE0)/(ALAM-AKA)
      C1=1+3*(D-1)*BK0*DELM
      C1M=3*BK0*(D-1)
      BKM=-3*BE
      G2=Q1/2/G+DLAM*(Q1-ALFQ1)/(AEM**2)
      G3=AEM**2*C1/(2*G*DLAM)+1.
C      
      IF(DABS(OCR-1.D0).LT.TOL1)THEN
      ADELM=3*((A1/D)**2*BK0-3*DELM/(4*DLAM1**2))
      GO TO 11
      ENDIF
C      
      AMI=(BKM-3*(C1M*DELM+C1)/(2*DLAM)-3*AT*BK0*RK)/AM
      AQJ=(AEM**2*C1M*(EPSQS-Q1/2./G)/DLAM+AQ*AMI)/G3
      ADELM=3*BK0*(A1/D)**2+AEM**2*AQJ*(EPSQS-Q1/2/G)/(2*G*DLAM**2)-
     *9*DELM/(4*DLAM**2) 
C
11    EPSQJ(1)=(2*EPSD1(1)/3.-EPSD1(2)/3.-EPSD1(3)/3.)/(2*EPSQ1)
      EPSQJ(2)=(-EPSD1(1)/3.+2*EPSD1(2)/3.-EPSD1(3)/3.)/(2*EPSQ1)
      EPSQJ(3)=(-EPSD1(1)/3.-EPSD1(2)/3.+2*EPSD1(3)/3.)/(2*EPSQ1)
      EPSQJ(4)=2*EPSD1(4)/(2*EPSQ1)
C
      DO 6 I=1,3
      DO 6 J=1,3
C      
      IF(DABS(OCR-1.D0).LT.TOL1)THEN
      QJ=EPSQJ(J)/(1./2/G+DLAM1/D/AEM**2)
      ADELM=3*(BK0*(A1/D)**2-3*DELM/4/DLAM**2)
      DELMJ=(Q1-ALFQ1)*EPSQJ(J)/(ADELM*DLAM1*(1.+D*AEM**2/2/G/DLAM1))
      GO TO 12
      ENDIF
C      
      BQJ=(AEM**2*C1*EPSQJ(J)/DLAM+AQ*BE/AM)/G3
      DELMJ=AEM**2*(EPSQS-Q1/2/G)*(EPSQJ(J)-BQJ/2/G)/(DLAM**2*ADELM)
      QJ=BQJ+AQJ*DELMJ
12    SMJ=BE+BKM*DELMJ
C
C      
      IF(I.EQ.J)THEN
      SIJ=(QJ*(EPSD1(I)-EPSD(I)+S(I)/2/G)+Q1*2/3.-S1(I)*(EPSQJ(J)-
     *3*DLAM*(Q1-ALFQ1)*(D-1)*BK0*DELMJ/(AEM**2*C1**2)))/G2
      ELSE      
      SIJ=(QJ*(EPSD1(I)-EPSD(I)+S(I)/2/G)-Q1/3.-S1(I)*(EPSQJ(J)-
     *3*DLAM*(Q1-ALFQ1)*(D-1)*BK0*DELMJ/(AEM**2*C1**2)))/G2
      ENDIF
C
6     CPE(I,J)=CPE(I,J)+SIJ+SMJ
C
      DO 7 I=1,3
C      
      IF(DABS(OCR-1.D0).LT.TOL1)THEN
      QJ=EPSQJ(4)/(1./2/G+DLAM1/D/AEM**2)
      ADELM=3*(BK0*(A1/D)**2-3*DELM/4/DLAM**2)
      DELMJ=(Q1-ALFQ1)*EPSQJ(4)/(ADELM*DLAM1*(1.+D*AEM**2/2/G/DLAM1))
      GO TO 13
      ENDIF   
C         
      BQJ=(AEM**2*C1*EPSQJ(4)/DLAM)/G3
      DELMJ=AEM**2*(EPSQS-Q1/2/G)*(EPSQJ(4)-BQJ/2/G)/(DLAM**2*ADELM)
      QJ=BQJ+AQJ*DELMJ
13    SMJ=BKM*DELMJ
C
C      
      SIJ=(QJ*(EPSD1(I)-EPSD(I)+S(I)/2/G)-S1(I)*(EPSQJ(4)-
     *3*DLAM*(Q1-ALFQ1)*(D-1)*BK0*DELMJ/(AEM**2*C1**2)))/G2
C      
7     CPE(I,4)=CPE(I,4)+SIJ+SMJ
C
      DO 8 J=1,3
C
      IF(DABS(OCR-1.D0).LT.TOL1)THEN
      ADELM=3*(BK0*(A1/D)**2-3*DELM/4/DLAM**2)
      DELMJ=(Q1-ALFQ1)*EPSQJ(J)/(ADELM*DLAM1*(1.+D*AEM**2/2/G/DLAM1))      
      QJ=EPSQJ(J)/(1./2/G+DLAM1/D/AEM**2)
      GO TO 14
      ENDIF   
C            
      BQJ=(AEM**2*C1*EPSQJ(J)/DLAM+AQ*BE/AM)/G3
      DELMJ=AEM**2*(EPSQS-Q1/2/G)*(EPSQJ(J)-BQJ/2/G)/(DLAM**2*ADELM)
      QJ=BQJ+AQJ*DELMJ
C
C      
14    SIJ=(-S1(4)*(EPSQJ(J)-3*DLAM*(Q1-ALFQ1)*
     *(D-1)*BK0*DELMJ/(AEM**2*C1**2))+
     *QJ*(EPSD1(4)-EPSD(4)+S(4)/2/G))/G2
C      
8     CPE(4,J)=CPE(4,J)+SIJ
C
      IF(DABS(OCR-1.D0).LT.TOL1)THEN
      ADELM=3*(BK0*(A1/D)**2-3*DELM/4/DLAM1**2)
      DELMJ=(Q1-ALFQ1)*EPSQJ(4)/(ADELM*DLAM1*(1.+D*AEM**2/2/G/DLAM1))      
      QJ=EPSQJ(4)/(1./2/G+DLAM1/D/AEM**2)
      GO TO 15
      ENDIF   
C      
      BQJ=(AEM**2*C1*EPSQJ(4)/DLAM)/G3
      DELMJ=AEM**2*(EPSQS-Q1/2/G)*(EPSQJ(4)-BQJ/2/G)/(DLAM**2*ADELM)
      QJ=BQJ+AQJ*DELMJ
      
C
15    SIJ=(QJ*(EPSD1(4)-EPSD(4)+S(4)/2/G)+Q1-S1(4)*(EPSQJ(4)-
     *3*DLAM*(Q1-ALFQ1)*(D-1)*BK0*DELMJ/(AEM**2*C1**2)))/G2
C      
      CPE(4,4)=CPE(4,4)+SIJ
C
20    DO 9 I=1,4
9     WRITE(3,555)(CPE(I,J),J=1,4)
555   FORMAT( 'CEP'//4E15.6)
C
      EPSD1(3)=EPSD4
      EPSD1(4)=EPSD3
      EPSD(3)=EPS4
      EPSD(4)=EPS3
      S1(3)=S4
      S1(4)=S3
      S(3)=ST4
      S(4)=ST3
CS    TRANSFORMACIJA ELAST
C      CALL CLEAR(ELAST,36)
C      DO 81 I=1,4
C      DO 81 J=1,4
C   81 ELAST(I,J)=CPE(I,J)
C
      RETURN
      END
C===========================================================
C     PROGRAM  ZA RESAVANJE NELINEARNE JEDNACINE
C============================================================
      SUBROUTINE TEQBIS1(X,DX,F,XM2,FM2,XP2,FP2,IM,IP,IBISE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA RESAVANJE NELINEARNE JEDNACINE F(X)=0.
C
C
      FUBRZ=1.0D0
      IF(IBISE.GT.2) GO TO 10
      IF(IBISE.GT.1) GO TO 5
      X1=X
      F1=F
    5 IF(IBISE.GT.2) GO TO 10
      IR=1
      IF(X.GT.X1.AND.F.GT.F1) IR=2
C
10    IF(F)20,100,30
20    IM=IM+1
      IF(IM.GT.1) GO TO 22
      XM2=X
      FM2=F
      GO TO 50
   22 IF(F.LT.FM2) GO TO 50
      XM1=XM2
      FM1=FM2
      XM2=X
      FM2=F
C  25 X=XM1+(XM2-XM1)/(FM2-FM1)*FM1
C     GO TO 100
      GO TO 50
   30 IP=IP+1
      IF(IP.GT.1) GO TO 32
      XP2=X
      FP2=F
      GO TO 50
   32 IF(F.GT.FP2) GO TO 50
      XP1=XP2
      FP1=FP2
      XP2=X
      FP2=F
C  35 X=XP1+(XP2-XP1)/(FP2-FP1)*FP1
C     GO TO 100
   50 CONTINUE
      IF(IBISE.EQ.1) GO TO 60
      IF(IM.GT.0.AND.IP.GT.0) GO TO 80
      IF(IP.GT.0.) GO TO 60
      DX=DABS(0.01D0*X)*FUBRZ
      IF(IR.EQ.1) GO TO 55
      X=X+DX
      GO TO 100
   55 X=X-DX
      GO TO 100
   60 DX=DABS(0.01D0*X)*FUBRZ
      IF(IR.EQ.1) GO TO 65
      X=X-DX
      GO TO 100
   65 X=X+DX
      GO TO 100
C
C     KORENI SU RAZDVOJENI - METOD BISEKCIJE
C
   80 IF(IM.GT.5.AND.IP.GT.5) GO TO 90
C  80 CONTINUE
      X=(XP2+XM2)/2.0D0
      DX=XP2-XM2
      GO TO 100
C
C     METOD SECICE
C
C  90 IF(IBISE.GT.10) GO TO 0
   90 CONTINUE
      X=XP2-(XM2-XP2)/(FM2-FP2)*FP2
      DX=XP2-XM2
  100 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------C
C===================================================================
      SUBROUTINE BISEC1(X,XM,XP,DX,F,FM,FP,AF,IB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BISEC '
C            
C	SOLVE A NONLINEAR EQUATION BY BISECTION
C
C	X = ARGUMENT (POSITIVE)
C	XM = X-MINUS
C	XP = X-PLUS
C	F = FUNCTION
C	FM = F-MINUS
C	FP = F-PLUS
C	DX = INCREMENT OF X
C	IB = BISECTION FLAG OF X
C	     0 = SEARCH FOR XM AND XP
C	     1 = XM AND XP OBTAIND (SEARC FOR SOLUTION BY BISECTION)
C
C	AF = ACCELERATION FACTOR FOR DX
C
      IF(IB.EQ.1) GO TO 100
C
      DX=AF*DX
      IF(F.LE.0D0)THEN
      FM=F
      XM=X
      X=0.5*(XM+XP)
      IB=1
      ELSE
      FP=F
      XP=X
      X=X+DX
      ENDIF
      RETURN
C
C	BISECTION
C
100   IF(F.LT.0D0)THEN
      XM1=X
      X=X-F/(FM-F)*(XM-X)
      IF(X.LE.XP) X=0.5*(XM1+XP)
      XM=XM1
      FM=F
      DX=XM-X
      ELSEIF (F.GT.0.D0) THEN
      XP1=X
      X=X+F/(FP-F)*(X-XP)
      IF(X.GE.XM) X=0.5*(XM+XP1)
      XP=XP1
      FP=F
      DX=X-XP
      ELSE
      DX=0.D0
      ENDIF
      RETURN
      END


