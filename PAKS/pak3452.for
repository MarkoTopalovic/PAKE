C=======================================================================
C
CS    MODEL KRISTALA (IZ ABAQUS-a)
CE
      SUBROUTINE D3M52(TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
CS    NA NIVOU INTEGRACIONE TACKE
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
      include 'paka.inc'
      
C     DODAT COMMON ZA ELAST
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      
C     DODAT COMMON ZA DT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR   
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      COMMON /ITERBR/ ITER
         
C     DODAT DIMENSION ZA DROT,NKONS
      DIMENSION TAU(6),DEF(6),DROT(3,3),nkons(4)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M52'
      STOP 'PAK 3452'
      
C      
C     DODELA VREDNOSTI MATRICI DROT=I
C
C      do 10 i=1,3
C      do 20 j=1,3
C      if (i.eq.j) then
C      drot(i,j)=1.0D0
C      else
C      drot(i,j)=0.0D0
C      endif
C   20 continue
C   10 continue
C     
C     NKONS - VEKTOR KONSTANTI POTREBNIH LOKALNO (NDI,NSHR,NTENS,ND)   
C  
      nkons(1)=3
      nkons(2)=3
      nkons(3)=6
      nkons(4)=12
C  
C     REPERI ZA MATERIJALE U PAK-U
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C       
C     REPERI ZA VELICINE VEKTORA STATEV U TRENUTKU T       
C
      LGSLIPT=LPOCG
      LGAMMAT=LGSLIPT+12*IDVA
      LTAUSLPT=LGAMMAT+12*IDVA
      LSLPNORT=LTAUSLPT+12*IDVA
      LSLPDIRT=LSLPNORT+3*12*IDVA
      LGMSLTLT=LSLPDIRT+3*12*IDVA
      LGAMTOLT=LGMSLTLT+12*IDVA
      LNSLIPT=LGAMTOLT+1*IDVA
      LNSLPTLT=LNSLIPT+3*IDVA
C
C     REPERI ZA DEFORMACIJU, PRIRASTAJ DEFORMACIJE I NAPON U TRENUTKU T
C
      LDEFT=LNSLPTLT+1*IDVA
      LDSTRANT=LDEFT+6*IDVA      
      LSTREST=LDSTRANT+6*IDVA
      LFSLIPT=LSTREST+6*IDVA
C     Ukupno 155*IDVA      
C    
C     REPERI ZA VELICINE VEKTORA STATEV U TRENUTKU t+dt
C
      LGSLIPTT=LPOC1
      LGAMMATT=LGSLIPTT+12*IDVA
      LTAUSLPTT=LGAMMATT+12*IDVA
      LSLPNORTT=LTAUSLPTT+12*IDVA
      LSLPDIRTT=LSLPNORTT+3*12*IDVA
      LGMSLTLTT=LSLPDIRTT+3*12*IDVA
      LGAMTOLTT=LGMSLTLTT+12*IDVA
      LNSLIPTT=LGAMTOLTT+1*IDVA
      LNSLPTLTT=LNSLIPTT+3*IDVA
C      
C     REPERI ZA DEFORMACIJU, PRIRASTAJ DEFORMACIJE I NAPON I TRENUTKU T+dT      
C     
      LDEFTT=LNSLPTLTT+1*IDVA
      LDSTRANTT=LDEFTT+6*IDVA
      LSTRESTT=LDSTRANTT+6*IDVA
      LFSLIPTT=LSTRESTT+6*IDVA
C
C     POZIV SUBROUTINE TI3452
C
      CALL TI3452(
     1A(LGSLIPT),A(LGSLIPTT),
     1A(LGAMMAT),A(LGAMMATT),
     1A(LTAUSLPT),A(LTAUSLPTT),
     1A(LSLPNORT),A(LSLPNORTT),
     1A(LSLPDIRT),A(LSLPDIRTT),
     1A(LGMSLTLT),A(LGMSLTLTT),
     1A(LGAMTOLT),A(LGAMTOLTT),
     1A(LNSLIPT),A(LNSLIPTT),
     1A(LNSLPTLT),A(LNSLPTLTT),
     1A(LDSTRANT),A(LDSTRANTT),
     1A(LFSLIPT),A(LFSLIPTT),
     1nkons(1),nkons(2),nkons(3),nkons(4),
     1tau,A(LSTRESTT),A(LSTREST),
     1def,A(LDEFTT),A(LDEFT),
     1elast,dt,A(LFUN),drot,IRAC)
      
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3452(
     1gslipt,gslip,
     1gammat,gamma,
     1tauslpt,tauslp,
     1slpnort,slpnor,
     1slpdirt,slpdir,
     1gmsltlt,gmsltl,
     1gamtolt,gamtol,
     1nslipt,nslip,
     1nslptlt,nslptl,
     1dstrant,dstran,
     1fslipt,fslip,
     1ndi,nshr,ntens,nd,
     1stress,tautt,taut,
     1def,deftt,deft,
     1ddsdde,dtime,props,drot,IRAC)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
C     DODAT COMMON ZA BROJANJE ITERACIJA ZBOG USLOVA ZA NULTU ITERACIJU
      COMMON /ITERBR/ ITER 
      COMMON /PRINCI/ PRINC(3)  
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR   
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK   

      EXTERNAL F
C    
      dimension stress(ntens),taut(ntens),tautt(ntens),
     1          ddsdde(ntens,ntens),props(160),
     1          drot(3,3),def(ntens),deft(ntens),deftt(ntens),
     1          dstran(ntens),dstrant(ntens),
     1          gmsltl(nslptl),gmsltlt(nslptl)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), nslipt(3),
     2          SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND),
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND),
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3),
     5          FSLIP(ND), DFDXSP(ND), DDEMSD(6,ND),
     6          H(ND,ND), DDGDDE(ND,6),
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND),
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3)

      DIMENSION FSLIPT(ND),GAMMA(ND),gammat(ND),TAUSLP(ND),tauslpt(ND),
     1          GSLIP(ND),gslipt(ND),slpnort(3,ND),slpdirt(3,ND),
     1          DGAMOD(ND),DHDGDG(ND,ND),POM(3,3)

      IF(IRAC.EQ.2) RETURN
      
C     GLAVNE VREDNOSTI
C     INVERZNO LAMBDA
      write (3,*) 'princ',(princ(i),i=1,3)
      write (3,*) 'QP',((qp(i,j),j=1,3),i=1,3)
      write (3,*) 'XJ',((XJ(i,j),j=1,3),i=1,3)
      P1=1.D0/DSQRT(PRINC(1))
      P2=1.D0/DSQRT(PRINC(2))
      P3=1.D0/DSQRT(PRINC(3))
C     INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
      CALL DIJAD(POM,QP,QP,P1,P2,P3)
C     TENZOR ROTACIJE R
C     R = Fe * Ue**-1 
      CALL MNOZM1(DROT,XJ,POM,3,3,3)
      write (3,*) 'drot',((drot(i,j),j=1,3),i=1,3)
      
CVLA  Prirastaj deformacije u koraku
CVLA  Nulta iteracija iskljucena zbog istih vrednosti u DEF i DEFT 

      if (iter.gt.0) then
      do 100 I=1,6
         dstran(I)=DEF(I)-DEFT(I)
  100 continue
      end if
CVLA  Odredjivanje vrste materijala na osnovu broja yadatih materijalnih 
CVLA  konstanti u ulaznom fajlu
C
C-----  Elastic matrix in local cubic crystal system: DLOCAL 

      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.0D0
         END DO
      END DO
      
      CHECK=0.0D0
      DO J=10,21
         CHECK=CHECK+DABS(PROPS(J))
      END DO

      IF (CHECK.EQ.0.0D0) THEN
         DO J=4,9
            CHECK=CHECK+DABS(PROPS(J))
         END DO

         IF (CHECK.EQ.0.0D0) THEN
            IF (PROPS(3).EQ.0.0D0) THEN

C-----  Isotropic material

               GSHEAR=PROPS(1)/2.0D0/(1.0D0+PROPS(2))
               E11=2.0D0*GSHEAR*(1.0D0-PROPS(2))/(1.0D0-2.0D0*PROPS(2))
               E12=2.0D0*GSHEAR*PROPS(2)/(1.0D0-2.0D0*PROPS(2))      
               DO J=1,3
                  DLOCAL(J,J)=E11

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=E12
                  END DO

                  DLOCAL(J+3,J+3)=GSHEAR
               END DO
            ELSE

C-----  Cubic material

               DO J=1,3
                  DLOCAL(J,J)=PROPS(1)

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=PROPS(2)
                  END DO

                  DLOCAL(J+3,J+3)=PROPS(3)
               END DO
            END IF

         ELSE

C-----  Orthotropic metarial

            DLOCAL(1,1)=PROPS(1)
            DLOCAL(1,2)=PROPS(2)
            DLOCAL(2,1)=PROPS(2)
            DLOCAL(2,2)=PROPS(3)

            DLOCAL(1,3)=PROPS(4)
            DLOCAL(3,1)=PROPS(4)
            DLOCAL(2,3)=PROPS(5)
            DLOCAL(3,2)=PROPS(5)
            DLOCAL(3,3)=PROPS(6)

            DLOCAL(4,4)=PROPS(7)
            DLOCAL(5,5)=PROPS(8)
            DLOCAL(6,6)=PROPS(9)

         END IF

      ELSE

C-----  General anisotropic material

         ID=0
         DO J=1,6
            DO I=1,J
               ID=ID+1
               DLOCAL(I,J)=PROPS(ID)
               DLOCAL(J,I)=DLOCAL(I,J)
            END DO
         END DO       

      END IF

C-----  Rotation matrix: ROTATE, i.e. direction cosines of [100], [010]
C     and [001] of a cubic crystal in global system
C
      CALL ROTATION (PROPS(57), ROTATE)
C      
C-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL
C     to global elastic matrix D
C
      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.0D0*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

C-----  Elastic matrix in global system: D
C     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose
C
      DO J=1,6
         DO I=1,6
            D(I,J)=0.0D0
         END DO
      END DO

      DO J=1,6
         DO I=1,J
            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO
            D(J,I)=D(I,J)
         END DO
      END DO
      
C-----  Total number of sets of slip systems: NSET

      NSET=NINT(PROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         STOP
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*)
     2     '***ERROR - more than three sets of slip systems'
         STOP
      END IF

C-----  Implicit integration parameter: THETA

      THETA=PROPS(145)

C-----  Finite deformation ?
C-----  NLGEOM = 0,   small deformation theory
C       otherwise, theory of finite rotation and finite strain, Users
C     must declare "NLGEOM" in the input file, at the *STEP card
C
      IF (PROPS(146).EQ.0.0D0) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF

C-----  Iteration?
C-----  ITRATN = 0, no iteration
C       otherwise, iteration (solving increments of stresses and
C     solution dependent state variables)
C
      IF (PROPS(153).EQ.0.0D0) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF
      ITRMAX=NINT(PROPS(154))
      GAMERR=PROPS(155)
      
      NITRTN=-1

      DO J=1,ND
         DGAMOD(J)=0.0D0
      END DO

C-----  Increment of spin associated with the material element: DSPIN
C     (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO

            TERM(J,J)=TERM(J,J)+1.D0
            TRM0(J,J)=TRM0(J,J)-1.D0
         END DO

         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP)

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)

      END IF
      
C-----  Increment of dilatational strain: DEV

      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO
C
C-----  Iteration starts (only when iteration method is used)
1000  CONTINUE

C-----  Parameter NITRTN: number of iterations
C       NITRTN = 0 --- no-iteration solution
C
      NITRTN=NITRTN+1
C-----  Check whether the current stress state is the initial state
C
      IF (gslipt(1).EQ.0.0D0) THEN
C
C-----  Initial state
C
C-----  Generating the following parameters and variables at initial
C     state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Unit vectors in initial slip directions SLPDIR
C          Unit normals to initial slip planes SLPNOR
C
         NSLPTL=0
         DO I=1,NSET
            ISPNOR(1)=NINT(PROPS(25+8*I))
            ISPNOR(2)=NINT(PROPS(26+8*I))
            ISPNOR(3)=NINT(PROPS(27+8*I))

            ISPDIR(1)=NINT(PROPS(28+8*I))
            ISPDIR(2)=NINT(PROPS(29+8*I))
            ISPDIR(3)=NINT(PROPS(30+8*I))
            
            CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), slpdir,
     2                    slpnor, ROTATE)
            call JEDNA1 (slpdirt, slpdir, 36)
            call JEDNA1 (slpnort, slpnor, 36)  
     
            NSLPTL=NSLPTL+NSLIP(I)
             
         END DO
         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*)
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
         END IF

C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
         SLPDEF(1,J)=slpdir(1,J)*slpnor(1,J)
         SLPDEF(2,J)=slpdir(2,J)*slpnor(2,J)
         SLPDEF(3,J)=slpdir(3,J)*slpnor(3,J)
         SLPDEF(4,J)=slpdir(1,J)*slpnor(2,J)+slpdir(2,J)*slpnor(1,J)
         SLPDEF(5,J)=slpdir(1,J)*slpnor(3,J)+slpdir(3,J)*slpnor(1,J)
         SLPDEF(6,J)=slpdir(2,J)*slpnor(3,J)+slpdir(3,J)*slpnor(2,J)
         END DO

C-----  Initial value of the current strength for all slip systems

        CALL GSLPINIT (gslip, NSLIP, NSLPTL, NSET, PROPS(97))
        CALL JEDNA1 (gslipt, gslip, 12)

C-----  Initial value of shear strain in slip systems
C-----  Initial value of cumulative shear strain in each slip systems

         DO I=1,NSLPTL
            gammat(I)=0.0D0
            gamma(I)=0.0D0
            gmsltlt(I)=0.0D0
            gmsltl(I)=0.0D0
         END DO
         gamtolt=0.0D0
         gamtol=0.0D0
         
C-----  Initial value of the resolved shear stress in slip systems
         DO I=1,NSLPTL
            TERM1=0.0D0
            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO
            tauslpt(I)=TERM1
            tauslp(I)=TERM1
         END DO

      ELSE

C-----  Slip deformation tensor: SLPDEF (Schmid factors)

         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

      END IF
      
C-----  Slip spin tensor: SLPSPN (only needed for finite rotation)

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=5.0D-1*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=5.0D-1*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=5.0D-1*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF

C-----  Double dot product of elastic moduli tensor with the slip
C     deformation tensor (Schmid factors) plus, only for finite
C     rotation, the dot product of slip spin tensor with the stress:
C     DDEMSD
C
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.0D0
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF

C-----  Shear strain-rate in a slip system at the start of increment:
C     FSLIP, and its derivative: DFDXSP
C         
      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL STRAINRATE (gamma, tauslp,
     2                    gslip, NSLIP(I), FSLIP, DFDXSP,
     3                    PROPS(65+8*I))
      END DO
      
C-----  Self- and latent-hardening laws

       CALL LATENTHARDEN (gamma, tauslp,
     2                   gslip, gmsltl,
     3                   gamtol, NSLIP, NSLPTL,
     4                   NSET, H(1,1), PROPS(97), ND)

C-----  LU decomposition to solve the increment of shear strain in a
C     slip system
C     
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
         X=TAUSLP(I)/GSLIP(I)
         TERM2=TERM1*DFDXSP(I)/GSLIP(I)
         TERM3=TERM1*X*DFDXSP(I)/GSLIP(I)

         DO J=1,NSLPTL
            TERM4=0.0D0
            DO K=1,6
               TERM4=TERM4+DDEMSD(K,I)*SLPDEF(K,J)
            END DO

            WORKST(I,J)=TERM2*TERM4+H(I,J)*TERM3*DSIGN(1.D0,FSLIP(J))
            
            IF (NITRTN.GT.0) WORKST(I,J)=WORKST(I,J)+TERM3*DHDGDG(I,J)

         END DO

         WORKST(I,I)=WORKST(I,I)+1.0D0
      END DO
      
      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP)     

C-----  Increment of shear strain in a slip system: DGAMMA

      TERM1=THETA*DTIME
      DO I=1,NSLPTL
         IF (NITRTN.EQ.0) THEN
            X=TAUSLP(I)/GSLIP(I)
            TERM2=TERM1*DFDXSP(I)/GSLIP(I)
            FSLIPT(I)=0.0D0
            DGAMMA(I)=0.0D0
            DO J=1,NDI
               DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF

            DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME

         ELSE

            DGAMMA(I)=TERM1*(FSLIP(I)-FSLIPT(I))+FSLIPT(I)*DTIME
     2                -DGAMOD(I)

         END IF
      END DO
      
      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)
       
      DO I=1,NSLPTL
         DGAMMA(I)=DGAMMA(I)+DGAMOD(I)
      END DO

C-----  Update the shear strain in a slip system: GAMMA
C
      DO I=1,NSLPTL
         gamma(I)=gammat(I)+DGAMMA(I)
      END DO

C-----  Increment of current strength in a slip system: DGSLIP
      DO I=1,NSLPTL
         DGSLIP(I)=0.0D0
         DO J=1,NSLPTL
            DGSLIP(I)=DGSLIP(I)+H(I,J)*DABS(DGAMMA(J))
         END DO
      END DO

C-----  Update the current strength in a slip system: GSLIP
C
      DO I=1,NSLPTL
         gslip(I)=gslipt(I)+DGSLIP(I)
      END DO

C-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.0D0
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

C-----  Increment of deformation gradient associated with lattice
C     stretching in the current state, i.e. the velocity gradient
C     (associated with lattice stretching) times the increment of time:
C     DVGRAD (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.0D0
                  ELSE
                     TERM1=-1.0D0
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                  END DO
               END IF

            END DO
         END DO

      END IF

C-----  Increment of resolved shear stress in a slip system: DTAUSP
      DO I=1,NSLPTL
         DTAUSP(I)=0.0D0
         DO J=1,6
            DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
         END DO
      END DO

C-----  Update the resolved shear stress in a slip system: TAUSLP

      DO I=1,NSLPTL
         tauslp(I)=tauslpt(I)+DTAUSP(I)
      END DO

C-----  Increment of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.0D0
         END DO
      ELSE
         DO I=1,NTENS
            DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF

      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR
            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
            END DO

         END DO
      END IF
      
C-----  Update the stress: STRESS

      DO I=1,NTENS
         STRESS(I)=TAUT(I)+DSTRES(I)
      END DO
 
C-----  Increment of normal to a slip plane and a slip direction (only
C     needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.0D0
               DSPDIR(I,J)=0.0D0

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO

C-----  Update the normal to a slip plane and a slip direction (only
C     needed for finite rotation)

         DO J=1,NSLPTL
            DO I=1,3
               slpnor(I,J)=slpnort(I,J)+DSPNOR(I,J)
               slpdir(I,J)=slpdirt(I,J)+DSPDIR(I,J)
            END DO
         END DO
      END IF

C-----  Derivative of shear strain increment in a slip system w.r.t.
C     strain increment: DDGDDE
C
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
            X=TAUSLP(J)/GSLIP(J)
            TERM2=TERM1*DFDXSP(J)/GSLIP(J)
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO
C-----  Derivative of stress increment w.r.t. strain increment, i.e.
C     Jacobian matrix
C
C-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.0D0
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF

C-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2                                DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2                            DDEMSD(I,K)*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2                            DDEMSD(J+3,K)*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF

      IF (ITRATN.NE.0) THEN
         DO J=1,NTENS
            DO I=1,NTENS
               DDSDDE(I,J)=DDSDDE(I,J)/(1.0D0+DEV)
            END DO
         END DO
      END IF
      
C-----  Iteration ?

      IF (ITRATN.NE.0) THEN

C-----  Solution dependent state
C     variables DGAMOD (for the next
C     iteration)

         DO J=1,NSLPTL
            DGAMOD(J)=DGAMMA(J)
         END DO

C-----  Check if the iteration solution converges

         IDBACK=0
         ID=0
         DO I=1,NSET        
            DO J=1,NSLIP(I)
               ID=ID+1
               X=tauslp(ID)/gslip(ID)
               RESIDU=THETA*DTIME*F(X,PROPS(65+8*I))+DTIME*(1.0-THETA)*
     2                FSLIPT(ID)-DGAMMA(ID)
               IF (DABS(RESIDU).GT.GAMERR) IDBACK=1
            END DO
         END DO

         IF (IDBACK.NE.0.AND.NITRTN.LT.ITRMAX) THEN

C-----  Iteration: arrays for iteration

            CALL ITERATION (gamma, tauslp,
     2                      gslip, gmsltl,
     3                      gamtol, NSLPTL,
     4                      NSET, NSLIP, ND, PROPS(97), DGAMOD,
     5                      DHDGDG)

CFIXB     
            GO TO 1000

         ELSE IF (NITRTN.GE.ITRMAX) THEN

         write (3,*) 'Solution not converge within maximum number of 
     1 iteration (the solution without iteration will be used)'
            
         END IF
      END IF

C-----  Total cumulative shear strains on all slip systems (sum of the
C       absolute values of shear strains in all slip systems)
CFIX--  Total cumulative shear strains on each slip system (sum of the
CFIX    absolute values of shear strains in each individual slip system)
C
      DO I=1,NSLPTL
         gamtol=gamtolt+DABS(DGAMMA(I))
         gmsltl(I)=gmsltlt(I)+DABS(DGAMMA(I))
      END DO     
C     
      do 200 I=1,6
         deftt(I)=DEF(I)
         tautt(I)=STRESS(I)
  200 continue
  
      write (3,*) 'tauslp na kraju',(tauslp(i),i=1,12)
      write (3,*) 'stress na kraju',(tautt(i),i=1,6)
      write (3,*) 'stran',(deftt(i),i=1,6)
     
      CALL TRANSS(TSG,QP)
C     Pg=Qs*Pd
      CALL CLEAR(DEF,6)
      CALL MNOZI1(DEF,TSG,DDEFPS,6,6)
        
      RETURN
      END

C----------------------------------------------------------------------


      SUBROUTINE ROTATION (PROP, ROTATE)

C-----  This subroutine calculates the rotation matrix, i.e. the
C     direction cosines of cubic crystal [100], [010] and [001]
C     directions in global system

C-----  The rotation matrix is stored in the array ROTATE.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3)

C-----  Subroutines:
C
C       CROSS  -- cross product of two vectors
C
C       LUDCMP -- LU decomposition
C
C       LUBKSB -- linear equation solver based on LU decomposition
C                 method (must call LUDCMP first)


C-----  PROP -- constants characterizing the crystal orientation
C               (INPUT)
C
C            PROP(1) - PROP(3) -- direction of the first vector in
C                                 local cubic crystal system
C            PROP(4) - PROP(6) -- direction of the first vector in
C                                 global system
C
C            PROP(9) - PROP(11)-- direction of the second vector in
C                                 local cubic crystal system
C            PROP(12)- PROP(14)-- direction of the second vector in
C                                 global system
C
C-----  ROTATE -- rotation matrix (OUTPUT):
C
C            ROTATE(i,1) -- direction cosines of direction [1 0 0] in
C                           local cubic crystal system
C            ROTATE(i,2) -- direction cosines of direction [0 1 0] in
C                           local cubic crystal system
C            ROTATE(i,3) -- direction cosines of direction [0 0 1] in
C                           local cubic crystal system

C-----  local matrix: TERM1
      CALL CROS (PROP(1), PROP(9), TERM1, ANGLE1)

C-----  LU decomposition of TERM1

      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP)

C-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.0D0
            ELSE
               TERM2(I,J)=0.0D0
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

C-----  global matrix: TERM1
      CALL CROS (PROP(4), PROP(12), TERM1, ANGLE2)

C-----  Check: the angle between first and second vector in local and
C     global systems must be the same.  The relative difference must be
C     less than 0.1%.
C
      IF (DABS(ANGLE1/ANGLE2-1.0D0).GT.1.0D-3) THEN
         WRITE (6,*)
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF

C-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.0D0
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END


C-----------------------------------


           SUBROUTINE CROS (A, B, C, ANGLE)

C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.0D0) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.0D0) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.0D0
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*)
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE SLIPSYS (ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR,
     2                    ROTATE)

C-----  This subroutine generates all slip systems in the same set for
C     a CUBIC crystal.  For other crystals (e.g., HCP, Tetragonal,
C     Orthotropic, ...), it has to be modified to include the effect of
C     crystal aspect ratio.

C-----  Denote s as a slip direction and m as normal to a slip plane.
C     In a cubic crystal, (s,-m), (-s,m) and (-s,-m) are NOT considered
C     independent of (s,m).

C-----  Subroutines:  LINE1 and LINE

C-----  Variables:
C
C     ISPDIR -- a typical slip direction in this set of slip systems
C               (integer)  (INPUT)
C     ISPNOR -- a typical normal to slip plane in this set of slip
C               systems (integer)  (INPUT)
C     NSLIP  -- number of independent slip systems in this set
C               (OUTPUT)
C     SLPDIR -- unit vectors of all slip directions  (OUTPUT)
C     SLPNOR -- unit normals to all slip planes  (OUTPUT)
C     ROTATE -- rotation matrix (INPUT)
C          ROTATE(i,1) -- direction cosines of [100] in global system
C          ROTATE(i,2) -- direction cosines of [010] in global system
C          ROTATE(i,3) -- direction cosines of [001] in global system
C
C     NSPDIR -- number of all possible slip directions in this set
C     NSPNOR -- number of all possible slip planes in this set
C     IWKDIR -- all possible slip directions (integer)
C     IWKNOR -- all possible slip planes (integer)

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISPDIR(3), ISPNOR(3), SLPDIR(3,50), SLPNOR(3,50),
     *          ROTATE(3,3), IWKDIR(3,24), IWKNOR(3,24), TERM(3) 
     
      NSLIP=0
      NSPDIR=0
      NSPNOR=0

C-----  Generating all possible slip directions in this set
C
C       Denote the slip direction by [lmn].  I1 is the minimum of the
C     absolute value of l, m and n, I3 is the maximum and I2 is the
C     mode, e.g. (1 -3 2), I1=1, I2=2 and I3=3.  I1<=I2<=I3.

      I1=MIN(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I3=MAX(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I2=IABS(ISPDIR(1))+IABS(ISPDIR(2))+IABS(ISPDIR(3))-I1-I3

      RMODIR=SQRT(FLOAT(I1*I1+I2*I2+I3*I3))

C     I1=I2=I3=0
      IF (I3.EQ.0) THEN
         WRITE (3,*) '***ERROR - slip direction is [000]'
         STOP

C     I1=I2=0, I3>0   ---   [001] type
      ELSE IF (I2.EQ.0) THEN
         NSPDIR=3
         DO J=1,3
            DO I=1,3
               IWKDIR(I,J)=0
               IF (I.EQ.J) IWKDIR(I,J)=I3
            END DO
         END DO

C     I1=0, I3>=I2>0
      ELSE IF (I1.EQ.0) THEN

C        I1=0, I3=I2>0   ---   [011] type
         IF (I2.EQ.I3) THEN
            NSPDIR=6
            DO J=1,6
               DO I=1,3
                  IWKDIR(I,J)=I2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKDIR(I,J)=0
                  IWKDIR(1,6)=-I2
                  IWKDIR(2,4)=-I2
                  IWKDIR(3,5)=-I2
               END DO
            END DO

C        I1=0, I3>I2>0   ---   [012] type
         ELSE
            NSPDIR=12
            CALL LINE1 (I2, I3, IWKDIR(1,1), 1)
            CALL LINE1 (I3, I2, IWKDIR(1,3), 1)
            CALL LINE1 (I2, I3, IWKDIR(1,5), 2)
            CALL LINE1 (I3, I2, IWKDIR(1,7), 2)
            CALL LINE1 (I2, I3, IWKDIR(1,9), 3)
            CALL LINE1 (I3, I2, IWKDIR(1,11), 3)

         END IF

C     I1=I2=I3>0   ---   [111] type
      ELSE IF (I1.EQ.I3) THEN
         NSPDIR=4
         CALL LINE (I1, I1, I1, IWKDIR)

C     I3>I2=I1>0   ---   [112] type
      ELSE IF (I1.EQ.I2) THEN
         NSPDIR=12
         CALL LINE (I1, I1, I3, IWKDIR(1,1))
         CALL LINE (I1, I3, I1, IWKDIR(1,5))
         CALL LINE (I3, I1, I1, IWKDIR(1,9))

C     I3=I2>I1>0   ---   [122] type
      ELSE IF (I2.EQ.I3) THEN
         NSPDIR=12
         CALL LINE (I1, I2, I2, IWKDIR(1,1))
         CALL LINE (I2, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I2, I1, IWKDIR(1,9))

C     I3>I2>I1>0   ---   [123] type
      ELSE
         NSPDIR=24
         CALL LINE (I1, I2, I3, IWKDIR(1,1))
         CALL LINE (I3, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I3, I1, IWKDIR(1,9))
         CALL LINE (I1, I3, I2, IWKDIR(1,13))
         CALL LINE (I2, I1, I3, IWKDIR(1,17))
         CALL LINE (I3, I2, I1, IWKDIR(1,21))

      END IF

C-----  Generating all possible slip planes in this set
C
C       Denote the normal to slip plane by (pqr).  J1 is the minimum of
C     the absolute value of p, q and r, J3 is the maximum and J2 is the
C     mode, e.g. (1 -2 1), J1=1, J2=1 and J3=2.  J1<=J2<=J3.

      J1=MIN(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J3=MAX(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J2=IABS(ISPNOR(1))+IABS(ISPNOR(2))+IABS(ISPNOR(3))-J1-J3

      RMONOR=SQRT(FLOAT(J1*J1+J2*J2+J3*J3))

      IF (J3.EQ.0) THEN
         WRITE (3,*) '***ERROR - slip plane is [000]'
         STOP

C     (001) type
      ELSE IF (J2.EQ.0) THEN
         NSPNOR=3
         DO J=1,3
            DO I=1,3
               IWKNOR(I,J)=0
               IF (I.EQ.J) IWKNOR(I,J)=J3
            END DO
         END DO

      ELSE IF (J1.EQ.0) THEN

C     (011) type
         IF (J2.EQ.J3) THEN
            NSPNOR=6
            DO J=1,6
               DO I=1,3
                  IWKNOR(I,J)=J2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKNOR(I,J)=0
                  IWKNOR(1,6)=-J2
                  IWKNOR(2,4)=-J2
                  IWKNOR(3,5)=-J2
               END DO
            END DO

C     (012) type
         ELSE
            NSPNOR=12
            CALL LINE1 (J2, J3, IWKNOR(1,1), 1)
            CALL LINE1 (J3, J2, IWKNOR(1,3), 1)
            CALL LINE1 (J2, J3, IWKNOR(1,5), 2)
            CALL LINE1 (J3, J2, IWKNOR(1,7), 2)
            CALL LINE1 (J2, J3, IWKNOR(1,9), 3)
            CALL LINE1 (J3, J2, IWKNOR(1,11), 3)

         END IF

C     (111) type
      ELSE IF (J1.EQ.J3) THEN
         NSPNOR=4
         CALL LINE (J1, J1, J1, IWKNOR)

C     (112) type
      ELSE IF (J1.EQ.J2) THEN
         NSPNOR=12
         CALL LINE (J1, J1, J3, IWKNOR(1,1))
         CALL LINE (J1, J3, J1, IWKNOR(1,5))
         CALL LINE (J3, J1, J1, IWKNOR(1,9))

C     (122) type
      ELSE IF (J2.EQ.J3) THEN
         NSPNOR=12
         CALL LINE (J1, J2, J2, IWKNOR(1,1))
         CALL LINE (J2, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J2, J1, IWKNOR(1,9))

C     (123) type
      ELSE
         NSPNOR=24
         CALL LINE (J1, J2, J3, IWKNOR(1,1))
         CALL LINE (J3, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J3, J1, IWKNOR(1,9))
         CALL LINE (J1, J3, J2, IWKNOR(1,13))
         CALL LINE (J2, J1, J3, IWKNOR(1,17))
         CALL LINE (J3, J2, J1, IWKNOR(1,21))

      END IF

C-----  Generating all slip systems in this set
C
C-----  Unit vectors in slip directions: SLPDIR, and unit normals to
C     slip planes: SLPNOR in local cubic crystal system
C
      WRITE (3,*) '          '
      WRITE (3,*) ' #          Slip plane          Slip direction'

      DO J=1,NSPNOR
         DO I=1,NSPDIR

            IDOT=0
            DO K=1,3
               IDOT=IDOT+IWKDIR(K,I)*IWKNOR(K,J)
            END DO

            IF (IDOT.EQ.0) THEN
               NSLIP=NSLIP+1
               DO K=1,3
                  SLPDIR(K,NSLIP)=IWKDIR(K,I)/RMODIR
                  SLPNOR(K,NSLIP)=IWKNOR(K,J)/RMONOR
               END DO

               WRITE (3,10) NSLIP,
     2                      (IWKNOR(K,J),K=1,3), (IWKDIR(K,I),K=1,3)

            END IF

         END DO
      END DO
10    FORMAT(1X,I2,9X,'(',3(1X,I2),1X,')',10X,'[',3(1X,I2),1X,']')

      WRITE (3,*) 'Number of slip systems in this set = ',NSLIP
      WRITE (3,*) '          '

      IF (NSLIP.EQ.0) THEN
         WRITE (3,*)
     *      'There is no slip direction normal to the slip planes!'
         STOP

      ELSE

C-----  Unit vectors in slip directions: SLPDIR, and unit normals to
C     slip planes: SLPNOR in global system
C
         DO J=1,NSLIP
            DO I=1,3
               TERM(I)=0.0D0
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPDIR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPDIR(I,J)=TERM(I)
            END DO

            DO I=1,3
               TERM(I)=0.0D0
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPNOR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPNOR(I,J)=TERM(I)
            END DO
         END DO

      END IF

      RETURN
      END


C----------------------------------


           SUBROUTINE LINE (I1, I2, I3, IARRAY)

C-----  Generating all possible slip directions <lmn> (or slip planes
C     {lmn}) for a cubic crystal, where l,m,n are not zeros.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,4)

           DO J=1,4
              IARRAY(1,J)=I1
              IARRAY(2,J)=I2
              IARRAY(3,J)=I3
           END DO

           DO I=1,3
              DO J=1,4
                 IF (J.EQ.I+1) IARRAY(I,J)=-IARRAY(I,J)
              END DO
           END DO

           RETURN
           END


C-----------------------------------


           SUBROUTINE LINE1 (J1, J2, IARRAY, ID)

C-----  Generating all possible slip directions <0mn> (or slip planes
C     {0mn}) for a cubic crystal, where m,n are not zeros and m does
C     not equal n.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,2)

           IARRAY(ID,1)=0
           IARRAY(ID,2)=0

           ID1=ID+1
           IF (ID1.GT.3) ID1=ID1-3
           IARRAY(ID1,1)=J1
           IARRAY(ID1,2)=J1

           ID2=ID+2
           IF (ID2.GT.3) ID2=ID2-3
           IARRAY(ID2,1)=J2
           IARRAY(ID2,2)=-J2

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE GSLPINIT (GSLIP0, NSLIP, NSLPTL, NSET, PROP)

C-----  This subroutine calculates the initial value of current
C     strength for each slip system in a rate-dependent single crystal.
C     Two sets of initial values, proposed by Asaro, Pierce et al, and
C     by Bassani, respectively, are used here.  Both sets assume that
C     the initial values for all slip systems are the same (initially
C     isotropic).

C-----  These initial values are assumed the same for all slip systems
C     in each set, though they could be different from set to set, e.g.
C     <110>{111} and <110>{100}.

C-----  Users who want to use their own initial values may change the
C     function subprogram GSLP0.  The parameters characterizing these
C     initial values are passed into GSLP0 through array PROP.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GSLP0
      DIMENSION GSLIP0(NSLPTL), NSLIP(NSET), PROP(16,NSET)

C-----  Function subprograms:
C
C       GSLP0 -- User-supplied function subprogram given the initial
C                value of current strength at initial state

C-----  Variables:
C
C     GSLIP0 -- initial value of current strength (OUTPUT)
C
C     NSLIP  -- number of slip systems in each set (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C
C     PROP   -- material constants characterizing the initial value of
C               current strength (INPUT)
C
C               For Asaro, Pierce et al's law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C
C               For Bassani's law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of
C                            slip systems (or the breakthrough stress
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C

      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            GSLIP0(ID)=GSLP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET)
         END DO
      END DO

      RETURN
      END


C----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION GSLP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET)

C-----     User-supplied function subprogram given the initial value of
C        current strength at initial state

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION NSLIP(NSET), PROP(16)

           GSLP0=PROP(3)

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE STRAINRATE (GAMMA, TAUSLP, GSLIP, NSLIP, FSLIP,
     2                       DFDXSP, PROP)

C-----  This subroutine calculates the shear strain-rate in each slip
C     system for a rate-dependent single crystal.  The POWER LAW
C     relation between shear strain-rate and resolved shear stress
C     proposed by Hutchinson, Pan and Rice, is used here.

C-----  The power law exponents are assumed the same for all slip
C     systems in each set, though they could be different from set to
C     set, e.g. <110>{111} and <110>{100}.  The strain-rate coefficient
C     in front of the power law form are also assumed the same for all
C     slip systems in each set.

C-----  Users who want to use their own constitutive relation may
C     change the function subprograms F and its derivative DFDX,
C     where F is the strain hardening law, dGAMMA/dt = F(X),
C     X=TAUSLP/GSLIP.  The parameters characterizing F are passed into
C     F and DFDX through array PROP.

C-----  Function subprograms:
C
C       F    -- User-supplied function subprogram which gives shear
C               strain-rate for each slip system based on current
C               values of resolved shear stress and current strength
C
C       DFDX -- User-supplied function subprogram dF/dX, where x is the
C               ratio of resolved shear stress over current strength

C-----  Variables:
C
C     GAMMA  -- shear strain in each slip system at the start of time
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in each slip system (INPUT)
C     GSLIP  -- current strength (INPUT)
C     NSLIP  -- number of slip systems in this set (INPUT)
C
C     FSLIP  -- current value of F for each slip system (OUTPUT)
C     DFDXSP -- current value of DFDX for each slip system (OUTPUT)
C
C     PROP   -- material constants characterizing the strain hardening
C               law (INPUT)
C
C               For the current power law strain hardening law
C               PROP(1) -- power law hardening exponent
C               PROP(1) = infinity corresponds to a rate-independent
C               material
C               PROP(2) -- coefficient in front of power law hardening


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F, DFDX
      DIMENSION GAMMA(NSLIP), TAUSLP(NSLIP), GSLIP(NSLIP),
     2          FSLIP(NSLIP), DFDXSP(NSLIP), PROP(8)

      DO I=1,NSLIP
         X=TAUSLP(I)/GSLIP(I)
         FSLIP(I)=F(X,PROP)
         DFDXSP(I)=DFDX(X,PROP)
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION F(X,PROP)

C-----     User-supplied function subprogram which gives shear
C        strain-rate for each slip system based on current values of
C        resolved shear stress and current strength
C
C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           F=PROP(2)*(DABS(X))**PROP(1)*DSIGN(1.D0,X)

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
C
           REAL*8 FUNCTION DFDX(X,PROP)

C-----     User-supplied function subprogram dF/dX, where x is the
C        ratio of resolved shear stress over current strength

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           DFDX=PROP(1)*PROP(2)*(DABS(X))**(PROP(1)-1.0D0)

           RETURN
           END


C----------------------------------------------------------------------

CFIXA
      SUBROUTINE LATENTHARDEN (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL,
     2                         NSLIP, NSLPTL, NSET, H, PROP, ND)
CFIXB

C-----  This subroutine calculates the current self- and latent-
C     hardening moduli for all slip systems in a rate-dependent single
C     crystal.  Two kinds of hardening law are used here.  The first
C     law, proposed by Asaro, and Pierce et al, assumes a HYPER SECANT
C     relation between self- and latent-hardening moduli and overall
C     shear strain.  The Bauschinger effect has been neglected.  The
C     second is Bassani's hardening law, which gives an explicit
C     expression of slip interactions between slip systems.  The
C     classical three stage hardening for FCC single crystal could be
C     simulated.

C-----  The hardening coefficients are assumed the same for all slip
C     systems in each set, though they could be different from set to
C     set, e.g. <110>{111} and <110>{100}.

C-----  Users who want to use their own self- and latent-hardening law
C     may change the function subprograms HSELF (self hardening) and
C     HLATNT (latent hardening).  The parameters characterizing these
C     hardening laws are passed into HSELF and HLATNT through array
C     PROP.


C-----  Function subprograms:
C
C       HSELF  -- User-supplied self-hardening function in a slip
C                 system
C
C       HLATNT -- User-supplied latent-hardening function

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip 
C      system
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems
C               (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C
C     H      -- current value of self- and latent-hardening moduli
C               (OUTPUT)
C               H(i,i) -- self-hardening modulus of the ith slip system
C                         (no sum over i)
C               H(i,j) -- latent-hardening molulus of the ith slip
C                         system due to a slip in the jth slip system
C                         (i not equal j)
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets
C                            of slip systems to self-hardening in the
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of
C                            slip systems (or the breakthrough stress
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given
C                            interaction between slip systems in the
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given
C                            interaction between slip systems in the
C                            ith set and jth set (i not equal j)
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets
C                            of slip systems to self-hardening in the
C                            ith set of slip systems Q1
C
C     ND     -- leading dimension of arrays defined in subroutine UMAT
C               (INPUT)


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL HSELF, HLATNT
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET),
     3          H(ND,NSLPTL)
CFIXB

      CHECK=0.0D0
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+DABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO LATENT=1,NSLPTL
               IF (LATENT.EQ.ISELF) THEN
CFIXA
                  H(LATENT,ISELF)=HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                  NSET,NSLIP,PROP(1,I),CHECK,
     3                                  ISELF,ISET)
CFIXB
               ELSE
CFIXA
                  H(LATENT,ISELF)=HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                   NSET,NSLIP,PROP(1,I),CHECK,
     3                                   ISELF,ISET,LATENT)
CFIXB

               END IF
            END DO

         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP,CHECK,ISELF,ISET)
CFIXB

C-----     User-supplied self-hardening function in a slip system

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           IF (CHECK.EQ.0.0D0) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              HSELF=PROP(1)*TERM2**2

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.0D0
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT)
CFIXB

C-----     User-supplied latent-hardening function

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.0D0) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              HLATNT=PROP(1)*TERM2**2*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

              ID=0
              G=1.0D0
              DO I=1,NSET
                 IF (I.EQ.ISET) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 DO J=1,NSLIP(I)
                    ID=ID+1
                    IF (ID.NE.ISELF) THEN
CFIXA
		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
		    END IF

                 END DO
              END DO

              HLATNT=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------

CFIXA
      SUBROUTINE ITERATION (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL,
     2                      NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD,
     3                      DHDGDG)
CFIXB

C-----  This subroutine generates arrays for the Newton-Rhapson
C     iteration method.

C-----  Users who want to use their own self- and latent-hardening law
C     may change the function subprograms DHSELF (self hardening) and
C     DHLATN (latent hardening).  The parameters characterizing these
C     hardening laws are passed into DHSELF and DHLATN through array
C     PROP.


C-----  Function subprograms:
C
C       DHSELF -- User-supplied function of the derivative of self-
C                 hardening moduli
C
C       DHLATN -- User-supplied function of the derivative of latent-
C                 hardening moduli

C-----  Variables:
C
C     GAMMA  -- shear strain in all slip systems at the start of time
C               step  (INPUT)
C     TAUSLP -- resolved shear stress in all slip systems (INPUT)
C     GSLIP  -- current strength (INPUT)
CFIX  GMSLTL -- total cumulative shear strains on each individual slip 
C     system
CFIX            (INPUT)
C     GAMTOL -- total cumulative shear strains over all slip systems
C               (INPUT)
C     NSLPTL -- total number of slip systems in all the sets (INPUT)
C     NSET   -- number of sets of slip systems (INPUT)
C     NSLIP  -- number of slip systems in each set (INPUT)
C     ND     -- leading dimension of arrays defined in subroutine UMAT
C               (INPUT)
C
C     PROP   -- material constants characterizing the self- and latent-
C               hardening law (INPUT)
C
C               For the HYPER SECANT hardening law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- saturation stress TAUs in the ith set of
C                            slip systems
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets
C                            of slip systems to self-hardening in the
C                            ith set of slip systems Q1
C
C               For Bassani's hardening law
C               PROP(1,i) -- initial hardening modulus H0 in the ith
C                            set of slip systems
C               PROP(2,i) -- stage I stress TAUI in the ith set of
C                            slip systems (or the breakthrough stress
C                            where large plastic flow initiates)
C               PROP(3,i) -- initial critical resolved shear stress
C                            TAU0 in the ith set of slip systems
C               PROP(4,i) -- hardening modulus during easy glide Hs in
C                            the ith set of slip systems
C               PROP(5,i) -- amount of slip Gamma0 after which a given
C                            interaction between slip systems in the
C                            ith set reaches peak strength
C               PROP(6,i) -- amount of slip Gamma0 after which a given
C                            interaction between slip systems in the
C                            ith set and jth set (i not equal j)
C                            reaches peak strength
C               PROP(7,i) -- representing the magnitude of the strength
C                            of interaction in the ith set of slip
C                            system
C               PROP(8,i) -- representing the magnitude of the strength
C                            of interaction between the ith set and jth
C                            set of system
C               PROP(9,i) -- ratio of latent to self-hardening Q in the
C                            ith set of slip systems
C               PROP(10,i)-- ratio of latent-hardening from other sets
C                            of slip systems to self-hardening in the
C                            ith set of slip systems Q1
C
C-----  Arrays for iteration:
C
C       DGAMOD (INPUT)
C
C       DHDGDG (OUTPUT)
C

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL DHSELF, DHLATN
CFIXA
      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET),
     3          DGAMOD(NSLPTL), DHDGDG(ND,NSLPTL)
CFIXB

      CHECK=0.0D0
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+DABS(PROP(J,I))
         END DO
      END DO

C-----  CHECK=0   --  HYPER SECANT hardening law
C       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO KDERIV=1,NSLPTL
               DHDGDG(ISELF,KDERIV)=0.0D0

               DO LATENT=1,NSLPTL
                  IF (LATENT.EQ.ISELF) THEN
CFIXA
                     DHDG=DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           KDERIV)
CFIXB
                  ELSE
CFIXA
                     DHDG=DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           LATENT,KDERIV)
CFIXB
                  END IF

                  DHDGDG(ISELF,KDERIV)=DHDGDG(ISELF,KDERIV)+
     2                                 DHDG*DABS(DGAMOD(LATENT))
               END DO

            END DO
         END DO
      END DO

      RETURN
      END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of self-hardening
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL),
     2               NSLIP(NSET), PROP(16)
CFIXB

           IF (CHECK.EQ.0.0D0) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHSELF=-2.0D0*PROP(1)*TERM2**2*DTANH(TERM1)*TERM3

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.0D0*(PROP(1)-PROP(4))*TERM2**2*DTANH(TERM1)*TERM3
                 ID=0
                 G=1.0D0
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.0D0*DEXP(-TERM4)/(1.0D0+DEXP(-2.0D0*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHSELF=F*G

           END IF

           RETURN
           END


C-----------------------------------


C-----  Use single precision on cray
CFIXA
           REAL*8 FUNCTION DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT,
     3                            KDERIV)
CFIXB

C-----  User-supplied function of the derivative of latent-hardening
C     moduli

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
CFIXA
           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), NSLIP(NSET),
     2               PROP(16)
CFIXB

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.0D0) THEN

C-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHLATN=-2.0D0*PROP(1)*TERM2**2*DTANH(TERM1)*TERM3*Q

           ELSE

C-----  Bassani's hardening law
CFIXA
              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
CFIXB
              TERM2=2.0D0*DEXP(-TERM1)/(1.0D0+DEXP(-2.0D0*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.0D0*(PROP(1)-PROP(4))*TERM2**2*DTANH(TERM1)*TERM3
                 ID=0
                 G=1.0D0
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1
CFIXA
                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)
CFIXB
                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF
CFIXA
                 TERM4=GMSLTL(KDERIV)/GAMMA0
CFIXB
                 TERM5=2.0D0*DEXP(-TERM4)/(1.0D0+DEXP(-2.0D0*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHLATN=F*G*Q

           END IF

           RETURN
           END


C----------------------------------------------------------------------


      SUBROUTINE LUDCMP (A, N, NP, INDX, D)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0D-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.0D0
      DO I=1,N
         AAMAX=0.0D0

         DO J=1,N
            IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.0D0) PAUSE 'Singular matrix.'
         VV(I)=1.0D0/AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.0D0

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*DABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.0D0) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1.0D0/A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

      RETURN
      END


C----------------------------------------------------------------------


      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.0D0) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END
**
**
**
******************** Load step follows *****************************
