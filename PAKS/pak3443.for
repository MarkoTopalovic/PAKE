C==========================================================================
C==========================================================================
CS                      HOEK-BROWN MATERIJALNI MODEL
CE                      HOEK-BROWN MATERIAL MODEL
C==========================================================================
C==========================================================================
CE    SUBROUTINE D3M43
CE               TI3443
C
      SUBROUTINE D3M43(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
C
      include 'paka.inc'
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M43'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LEMP=LDEFPP + 6
      LXT=LEMP + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LEMP1=LDEFP1 + 6
      LXTDT=LEMP1 + 1
C
      CALL TI3443(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     &            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     &            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  ========================================================================
C
      SUBROUTINE TI3443(TAUT,DEFT,DEFPP,EMP,XT,
     &                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     &                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    HOEK-BROWN MODEL 
CE    PROGRAM FOR STRESS INTEGRATION FOR HOEK-BROWN MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     &               DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ITERBR/ ITER
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(6),DEFT(6),DEFPP(6),TAU(6),DEF(6),TAU1(6),DEF1(6), 
     &          DEFP1(6)
      DIMENSION FUN(11,*),NTA(*),DSIG(6),DEPS(6),DFDS(6),DGDS(6),ALAM(6)
     &         ,DDEFE(6)
      dimension di1ds(6),dj2dds(6),dthds(6),dj3dds(6),dgds1(6)
     &         ,dgds2(6),dtds1(6),dtds2(6),dfds1(6),dfds2(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3443'
C
      IF(IRAC.EQ.2) RETURN
c
CE    BASIC KONSTANTS
      TOLL =  1.D-6               ! tolerance
      MAXT =  500                 ! max. no. of iterations
      atriq=dsqrt(3.d0)
c==========================================================================
c     Material constants
      E    = FUN(1,MAT)           ! young's modulus
      ANI  = FUN(2,MAT)           ! poisson's ratio 
c
      smc  = FUN(3,MAT)           ! material constant
      em   = FUN(4,MAT)           ! material constant
      emd  = FUN(5,MAT)           ! material constant
      es   = FUN(6,MAT)           ! material constant
C==========================================================================
c
c==========================================================================
c       PLASTIC DOMAIN
c==========================================================================
c     elasticity matrix for 3d
      call mel3el(elast,e,ani)    ! formiranje matrice elasticnosti
c
c     {de}
      call jedna1(deps,def,6)     ! deps(i)=def(i)
      call zbirmm(deps,deft,6)    ! deps(i)=deps(i)-deft(i)
c
c     deltaSigma_E
      call clear(dsig,6)
      call mnozi1(dsig,elast,deps,6,6) !dsig(i)=elast(i,k)*deps(k)
c
c     Sigma_E
      call zbir2b(tau,taut,dsig,6)       ! tau(i)=taut(i)+dsig(i)
c
c     Prva invarijanta napona I1
      ai1=tau(1)+tau(2)+tau(3)           ! i1=sigma1+sigma2+sigma3
c
c     Druga invarijanta devijatora napona J2D
      aj2d=1.d0/6.d0*((tau(1)-tau(2))**2  +
     &                (tau(2)-tau(3))**2  +
     &                (tau(3)-tau(1))**2) +
     &                 tau(4)**2+tau(5)**2+tau(6)**2
c
c     Sqrt(J2D)
      AJ2D=DSQRT(AJ2D*AJ2D)
      aj2dq=dsqrt(aj2d)
c
C     Tension cutoff limit
      sigt=es*smc/em
      if(ai1.gt.sigt) then
        do i=1,6
          IF(I.LE.3)THEN
            tau(i)=sigt
          ELSE
            tau(i)=0.d0
          ENDIF
        enddo
       go to 400
      endif 
c
c     Increment of plastic strain is zerro in elastic domain
      call clear(DDEFP,6)
c
      demp=0.d0
c
c     Pan-Hudson yield curve
      Fph =ai1*em/3.d0+atriq/2.d0*em*aj2dq-es*smc+3.d0/smc*aj2d
c
      Fphe=Fph
c
c     Yielding check 
c     Fph>0
      if(Fphe.gt.toll) goto 100
c     Fph<0
      if(Fphe.le.toll) goto 400
c
c     U slucaju prolaska svih uslova
      stop 'ERROR! Hoke-Brown (Pan-Hudson) passed all conditions!!!'
c==========================================================================
c     PLASTIC DOMAIN
c==========================================================================
  100   continue
c-------------------------------------------------         
c ***** {dG/dSigma}T *****************************
c-------------------------------------------------
c ***** dG/dI1 ***********************************
        dgdi1=emd/3.d0
c
c ***** dG/dJ2D **********************************
        if(aj2dq.le.0.d0) stop 'aj2dq<0!, 1'
        dgdj2d =atriq/4.d0*emd/aj2dq+3.d0/smc
c
c ***** {dI1/dSigma}T ****************************
        di1ds(1)=1.d0
        di1ds(2)=1.d0
        di1ds(3)=1.d0
        di1ds(4)=0.d0
        di1ds(5)=0.d0
        di1ds(6)=0.d0
c
c ***** {dJ2D/dSigma}T ***************************
        dj2dds(1)=(2.d0*tau(1)      -tau(2)      -tau(3))/3.d0
        dj2dds(2)=(    -tau(1) +2.d0*tau(2)      -tau(3))/3.d0
        dj2dds(3)=(    -tau(1)      -tau(2) +2.d0*tau(3))/3.d0
        dj2dds(4)=2.d0*tau(4)
        dj2dds(5)=2.d0*tau(5)
        dj2dds(6)=2.d0*tau(6)
c
c ***** {dG/dSigma}T (suma)***********************
        call jednak(dgds1,di1ds,dgdi1,6)
        call jednak(dgds2,dj2dds,dgdj2d,6)
        call zbir2b(dgds,dgds1,dgds2,6)
c
c-------------------------------------------------         
c ***** {dF/dSigma}T *****************************
c-------------------------------------------------
c
c ***** dF/dI1 ***********************************
        dfdi1=em/3.d0
c
c ***** dF/dJ2D **********************************
        if(aj2dq.le.toll) stop 'aj2dq<0!, 2'
        dfdj2d =atriq/4.d0*em/aj2dq+3.d0/smc
c
c ***** {dF/dSigma}T (sum)************************
        call jednak(dfds1,di1ds,dfdi1,6)
        call jednak(dfds2,dj2dds,dfdj2d,6)
c
        call zbir2b(dfds,dfds1,dfds2,6)
c
c--------------------------------------------------------------------------  
c ***** dLambda **********************************
        call clear(alam,6)
        call mnozt1(alam,dfds,elast,6,6)      ! {dF/dSigma}T*[Ce]
        call clear(fcg,1)
        fcg=dot(alam,dgds,6)    ! {dF/dSigma}T*[Ce]*{dG/dSigma}
        call clear(fce,1)
        fce=dot(alam,deps,6)    ! {dF/dSigma}T*[Ce]*{de}
        dlam=fce/fcg
c
c--------------------------------------------------------------------------  
cr      Inicijalizacija za bisekcije
c--------------------------------------------------------------------------  
        call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
        call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
        call clear(dsig,6)
        call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
c
c       {Sigma_t+dt}={Sigma_t}+{dSigma}
        call zbir2b(tau,taut,dsig,6)
c
c       Prva invarijanta napona
        ai1=tau(1)+tau(2)+tau(3)       
c
c       J2D
        aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                  (tau(2)-tau(3))**2 +
     &                  (tau(3)-tau(1))**2)+
     &                   tau(4)**2+tau(5)**2+tau(6)**2
c
c       sqrt(J2D)
        AJ2D=DSQRT(AJ2D*AJ2D)
        aj2dq=dsqrt(aj2d)
c
        if(aj2dq.le.toll) stop 'aj2dq<0!, 3'
c
        Fphm =ai1*em/3.d0+atriq/2.d0*em*aj2dq-es*smc+3.d0/smc*aj2d
        Fphp=Fphe
        Fph=Fphe
        dlamm=dlam
        dlamp=0.d0
        dlam=0.d0
        dx=0.1d0*dlamm
        ib=0
        jp=2
        af=2.d0
        I=0
        write(3,*)'dlam, Fphe, Fphm',dlam,Fphe,Fphm
C=================================================
CD      BISECTION LOOP                            
        write(3,*)' I,      dlam,      dlamp,      dlamm,  
     &      Fph,       Fphp,       Fphm,kor'
  110   I = I + 1
c-------------------------------------------------          
c-------------------------------------------------          
            call BISECTG (dlam,dlamm,dlamp,dx,Fph,Fphm,Fphp,af,ib,jp)
c
            call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
            call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
            call clear(dsig,6)
            call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
c
            demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.d0 
            emp1=emp+demp
c
c           {Sigma_t+dt}={Sigma_t}+{dSigma}
            call zbir2b(tau,taut,dsig,6)
c
c           Prva invarijanta napona
            ai1=tau(1)+tau(2)+tau(3)       
c
c           J2D
            aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                      (tau(2)-tau(3))**2 +
     &                      (tau(3)-tau(1))**2)+
     &                       tau(4)**2+tau(5)**2+tau(6)**2
c
c           sqrt(J2D)
            AJ2D=DSQRT(AJ2D*AJ2D)
            aj2dq=dsqrt(aj2d)
c
c           Pan-Hudson yield curve
            if(aj2dq.le.0.d0) stop 'aj2dq<0!, 5'
            Fph =ai1*em/3.d0+atriq/2.d0*em*aj2dq-es*smc+3.d0/smc*aj2d
c
c           write(3,*)'dlam, Fphe, Fph',dlam,Fphe, Fph
            write(3,1001) I,dlam,dlamp,dlamm,Fph,Fphe,Fphm,kor
c
            if(I.gt.maxt) then
                stop 'Max. num. of bisection in Hoek-Brown (PH) model!'
            endif
c
 1001 FORMAT(I3,6E12.4,I3)
c
      if(dabs(Fph).gt.toll) goto 110
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.d0
      emp1=emp+demp
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      goto 400 
C==========================================================================
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
C
c========================================================================
c     Corection of values from previous step when convergence is reatched
      CALL ZBIR2B(DEFP1,DEFPP,DDEFP,6)
      call jedna1(def1,def,6)
      call jedna1(tau1,tau,6)
      return   
      end
C==========================================================================
