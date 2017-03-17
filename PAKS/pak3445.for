C==========================================================================
C==========================================================================
CS                      MOHR-COULOMB MATERIJALNI MODEL
CE                      MOHR-COULOMB MATERIAL MODEL
C==========================================================================
C==========================================================================
CE    SUBROUTINE D3M45
CE               TI3445
C
      SUBROUTINE D3M45(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
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
      IF(IDEBUG.GT.0) PRINT *, ' D3M45'
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
      CALL TI3445(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3445(TAUT,DEFT,DEFPP,EMP,XT,
     1                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     1                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    HIPERBOLICKI MOHR-COULOMB MODEL 
CE    PROGRAM FOR STRESS INTEGRATION FOR HUYPERBOLIC MOHR-COULOMB MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1               DETAU(6),DDEF(6)
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
      dimension di1ds(6),dj2dds(6),dthds(6),dj3dds(6),dgdi1T(6)
     &         ,dgdj2dT(6),dgdthT(6),dtdj3dT(6),dtdj2dT(6),dgdi11(6)
     &         ,dfdi1T(6),dfdj2dT(6),dfdthT(6),dfdi11(6),dphids(6)
     &         ,dfdphiT(6),dfds1(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3445'
C
      IF(IRAC.EQ.2) RETURN
c==========================================================================
CE    BASIC KONSTANTS
      TOLL = 1.D-6                ! tolerance
      MAXT = 500                  ! max. no. of iterations
      pi   = 4.d0*atan(1.d0)
c==========================================================================
c     Material constants
      E    = FUN(1,MAT)           ! young's modulus
      ANI  = FUN(2,MAT)           ! poisson's ratio 
c
      fib     = FUN(3,MAT)        ! basic angle
      dfi     = FUN(4,MAT)        ! angle diference
      pn      =-FUN(5,MAT)        !
C==========================================================================
      fib  = fib*pi/180.d0        ! deg. to rad.
      dfi  = dfi*pi/180.d0        ! deg. to rad.
C==========================================================================
c     elasticity matrix for 3d
      call mel3el(elast,e,ani)    ! formiranje matrice elasticnosti
c
c     {de}
      call jedna1(deps,def,6)     ! deps(i)=def(i)
      call zbirmm(deps,deft,6)    ! deps(i)=deps(i)-deft(i)
c      write(3,*)'def=',def
c      write(3,*)'deft=',deft
c      write(3,*)'deps=',deps
c
c     deltaSigma_E 
      call clear(dsig,6)
      call mnozm1(dsig,elast,deps,6,1,6) !dsig(i,j)=elast(i,k)*deps(k,j)
c      write(3,*)'dsig=',dsig
c
c     Sigma_E
      call zbir2b(tau,taut,dsig,6)       ! tau(i)=taut(i)+dsig(i)
c      write(3,*)'tau=',tau
c
      call clear(ddefp,6)
      if(kor.eq.1.and.iter.eq.0) goto 400
c
c     Prva invarijanta napona I1
      ai1=tau(1)+tau(2)+tau(3)           ! i1=sigma1+sigma2+sigma3
c      write(3,*)'ai1=',ai1
c
c     Tension cutoff
      if(ai1.gt.0.d0) then
        do i=1,6
            if(i.le.6)then
                tau(i)=0.d0
            endif
        enddo
        go to 400
      endif 
c
c     Druga invarijanta napona I2
      ai2= tau(1)*tau(2)+tau(2)*tau(3)+tau(3)*tau(1)
     &    -tau(4)**2-tau(5)**2-tau(6)**2
c      write(3,*)'ai2=',ai2
c
c     Treca invarijanta napona I3
      ai3= tau(1)*tau(2)*tau(3)
     &    -tau(1)*tau(5)**2-tau(2)*tau(6)**2-tau(3)*tau(4)**2
     &    +2.d0*tau(4)*tau(5)*tau(6)
c      write(3,*)'ai3=',ai3
c
c     Druga invarijanta devijatora napona J2D
      aj2d=1.d0/6.d0*((tau(1)-tau(2))**2  +
     &                (tau(2)-tau(3))**2  +
     &                (tau(3)-tau(1))**2) +
     &                 tau(4)**2+tau(5)**2+tau(6)**2
c      write(3,*)'aj2d=',aj2d
c
c     Treca invarijanta devijatora napona J3D
      aj3d=ai3-1.d0/3.d0*ai1*ai2+2.d0/27.d0*ai1**3
c      write(3,*)'aj3d=',aj3d
c
c     Sqrt(J2D)
      aj2dq=dsqrt(aj2d)
c
c     Sqrt(3)
      atriq=dsqrt(3.d0)
c
c     J3D/J2D^3/2  
      if(dabs(aj2dq).lt.toll) then
        aj2d3d=1.d0
      else
        aj2d3d=aj3d/(aj2dq**3)
      endif
c      write(3,*)'aj2d3d=',aj2d3d
c
c     Lode's angle argument  
      alode=-3.d0*atriq/2.d0*aj2d3d
      if(alode.gt. 4.d0) then
           alode= 1.d0
      endif
      if(alode.lt.-4.d0) then
           alode=-1.d0
      endif
      if(alode.gt. 3.d0) alode= 4.00d0-alode
      if(alode.lt.-3.d0) alode=-4.00d0-alode
      if(alode.gt. 2.d0) alode=-2.00d0+alode
      if(alode.lt.-2.d0) alode= 2.00d0+alode
      if(alode.gt. 1.d0) alode= 2.00d0-alode
      if(alode.lt.-1.d0) alode=-2.00d0-alode
c      alode= 1.d0
c
c     Lode's angle (Theta)
      theta=1.d0/3.d0*dasin(alode)
c
c     SigmaM
      sigmam=ai1/3.d0
c
c     fiM
      phim=fib+dfi/2.d0
c
c     pAV
      pav=pn/3.d0*(3.d0-dsin(phim))/(1.d0-(dsin(phim))**2)
c
c     fi
      phi=fib+dfi/(1.d0+sigmam/pav)
c
c     Increment of plastic strain is zerro in elastic domain
      call clear(DDEFP,6)
c
c      write(3,*)'sigmam,phim,pav,phi',sigmam,phim,pav,phi
c     Mohr-Coulomb yield curve for Maksimovic model
      Fmc =ai1*dsin(phi)/3.d0+aj2dq*(dcos(theta)-dsin(theta)*dsin(phi)
     &      /atriq)
      Fmce=Fmc
c      write(3,*)'Fmc,theta',Fmc,theta
c
c     Yielding check 
c     Fmc>0
      if(Fmc.gt.toll) goto 100
c     Fmc<0
      if(Fmc.le.toll) goto 400
c
c     U slucaju prolaska svih uslova
      stop 'Hiperbolicni Mohr-Coulomb prosao sve uslove!!!'
c==========================================================================
c         PLASTIC DOMAIN
c==========================================================================
  100     continue
c
c ******* dF/dI1 *********************************
          dfdi1=dsin(phi)/3.d0
c
c ******* {dI1/dSigma}T **************************
          di1ds(1)=1.d0
          di1ds(2)=1.d0
          di1ds(3)=1.d0
          di1ds(4)=0.d0
          di1ds(5)=0.d0
          di1ds(6)=0.d0
c
c ******* dF/dJ2D ********************************
          dfdj2d =1.d0/(2.d0*aj2dq)*(dcos(theta)
     &           -1.d0/atriq*dsin(theta)*dsin(phi))
c
c ******* {dJ2D/dSigma}T *************************
          dj2dds(1)=(2.d0*tau(1)      -tau(2)      -tau(3))/3.d0
          dj2dds(2)=(    -tau(1) +2.d0*tau(2)      -tau(3))/3.d0
          dj2dds(3)=(    -tau(1)      -tau(2) +2.d0*tau(3))/3.d0
          dj2dds(4)=2.d0*tau(4)
          dj2dds(5)=2.d0*tau(5)
          dj2dds(6)=2.d0*tau(6)
c
c ******* dF/dphi ********************************
          dfdphi=dcos(phi)*(ai1/3.d0-aj2dq/atriq*dsin(theta))
c
c ******* {dphi/dSigma}T *************************
          dphids(1)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(2)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(3)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(4)=0.d0
          dphids(5)=0.d0
          dphids(6)=0.d0
c
c--------------------------------------------------------------------------          
c         {dF/dSigma}T
          call jednak(dfdi1T,di1ds,dfdi1,6)
          call jednak(dfdj2dT,dj2dds,dfdj2d,6)
          call jednak(dfdphiT,dphids,dfdphi,6)
c
          call zbir2b(dfds1,dfdi1T,dfdj2dT,6)
          call zbir2b(dfds,dfds1,dfdphiT,6)
c
          call jedna1(dgds,dfds,6)
c--------------------------------------------------------------------------  
c         lambda
          call clear(alam,6)
          call mnozt1(alam,dfds,elast,6,6)      ! {dFmc/dSigma}T*[Ce]
          fcg=dot(alam,dgds,6)    ! {dFmc/dSigma}T*[Ce]*{dGmc/dSigma}
          fce=dot(alam,deps,6)    ! {dFmc/dSigma}T*[Ce]*{de}
          dlam=fce/fcg
c          write(3,*)'dlam',dlam
c
C==========================================================================
c         Inicijalizacija za bisekcije
          call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
          call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
          call clear(dsig,6)
          call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
c
c         {Sigma_t+dt}={Sigma_t}+{dSigma}
          call zbir2b(tau,taut,dsig,6)
c
c         Prva invarijanta napona
          ai1=tau(1)+tau(2)+tau(3)       
c
c         J2D
          aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                    (tau(2)-tau(3))**2 +
     &                    (tau(3)-tau(1))**2)+
     &                     tau(4)**2+tau(5)**2+tau(6)**2
c
c         sqrt(J2D)
          aj2dq=dsqrt(aj2d)
c
c         Sqrt(3)
          atriq=dsqrt(3.d0)
c
c         SigmaM
          sigmam=ai1/3.d0
c
c         fiM
          phim=fib+dfi/2.d0
c
c         pAV
          pav=pn/3.d0*(3.d0-dsin(phim))/(1.d0-(dsin(phim))**2)
c
c         fi
          phi=fib+dfi/(1.d0+sigmam/pav)
c
          Fmcm=ai1*dsin(phi)/3.d0+aj2dq*(dcos(theta)
     &           -dsin(theta)*dsin(phi)/atriq)
          Fmcp=Fmce
c
          I=0
          dlamp=0.d0
          dlamm=dlam
          dx=0.01*dlamm
          Fmc=Fmcp
          dlam=0.d0
          af=2
          ib=0
          jp=2
c            
C==========================================================================
          call clear(ddefp,6)
CD        BISECTION LOOP
c          write(3,*)'alode, theta,Fmce,Fmcp'
c     &              ,alode, theta*180.d0/pi,Fmce,Fmcm
c          write(3,*)' I,      dlam,      dlamp,      dlamm,  
c     &      Fmc,       Fmcp,       Fmcm,'
  110     I = I + 1
            call BISECTG (dlam,dlamm,dlamp,dx,fmc,fmcm,fmcp,af,ib,jp)
c--------------------------------------------------------------------------          
            call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
            call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
            call clear(dsig,6)
            call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
c
c           SIGMA
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
            aj2dq=dsqrt(aj2d)
c
c           Sqrt(3)
            atriq=dsqrt(3.d0)
c
c           SigmaM
            sigmam=ai1/3.d0
c
c           fiM
            phim=fib+dfi/2.d0
c
c           pAV
            pav=pn/3.d0*(3.d0-dsin(phim))/(1.d0-(dsin(phim))**2)
c
c           fi
            phi=fib+dfi/(1.d0+sigmam/pav)
c
c           Mohr-Coulomb yield curve
            Fmc=ai1*dsin(phi)/3.d0+aj2dq*(dcos(theta)
     &             -dsin(theta)*dsin(phi)/atriq)
c
c            write(3,1001) I,dlam,dlamp,dlamm,Fmc,Fmcp,Fmcm
c
            if(I.gt.maxt) then
                stop 'Max. num. of bisection in Mohr-Coulomb model!'
            endif
c
            if(dabs(Fmc).gt.toll) goto 110
c
c 1001       FORMAT(I3,6E12.4)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            demp=(ddefp(1)+ddefp(2)+ddefp(3))/3
            emp1=emp+demp
            evp1=3*emp1
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      goto 400 
C==========================================================================

CE    UPDATES FOR NEXT STEP
  400 CONTINUE
C
      CALL ZBIR2B(DEFP1,DEFPP,DDEFP,6)
c      write(3,*)'DEFP1,DEFPP,DDEFP',DEFP1,DEFPP,DDEFP
C
c========================================================================
c     Corection of values from previous step when convergence is reatched
      call jedna1(def1,def,6)
      call jedna1(tau1,tau,6)
c      write(3,*)'def1,def',def1,def
c      write(3,*)'tau1,tau',tau1,tau
      return   
      end
C==========================================================================

