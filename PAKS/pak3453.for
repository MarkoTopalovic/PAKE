C=======================================================================
C
CS    MODEL SHAPE MEMORY ALLOY (LAGOUDAS)  UGRADJENO 07.11.2011
CE
      SUBROUTINE D3M53(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
C
CS    PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
CS    NA NIVOU INTEGRACIONE TACKE
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION POINT LEVEL

C     DODATI COMMON zbog ELAST
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6),STATEV(100)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M53'
C  
C     REPERI ZA MATERIJALE U PAK-U
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)                  
C           
C     REPERI ZA VELICINE VEKTORA STATEV U TRENUTKU T       
C
      LTAU = LPOCG
      LDEF = LTAU + 6
      LSTATEV = LDEF + 6
      LTEMP = LSTATEV + 100
C    
C     REPERI ZA VELICINE VEKTORA STATEV U TRENUTKU T+dT
C
      LTAU1 = LPOC1
      LDEF1 = LTAU1 + 6
      LSTATEV1 = LDEF1 + 6
      LTEMP1 = LSTATEV1 + 100
C
C     POZIV SUBROUTINE TI3453
C      
      CALL TI3453(STATEV,PLAS1(LSTATEV1),PLAST(LSTATEV),
     1  TGT,PLAS1(LTEMP1),PLAST(LTEMP),
     1  TAU,PLAS1(LTAU1),PLAST(LTAU),
     1  DEF,PLAS1(LDEF1),PLAST(LDEF),
     1  ELAST,A(LFUN),IRAC,DT,VREME,ITER)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3453(
     1STATEV,STATEVTT,STATEVT,
     1TGT,TEMPTT,TEMPT,
     1TAU,TAUTT,TAUT,
     1DEF,DEFTT,DEFT,
     1DT,PROPS,IRAC,DTIME,VREME,ITERPAK)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXMTP=20000)
C
      DIMENSION STRESS(6),STATEV(100),STATEVT(100),STATEVTT(100),
     1TAUT(6),TAUTT(6),DDSDDE(6,6),DDSDDT(6),DRPLDE(6),
     2STRAN(6),DEFT(6),DEFTT(6),DSTRAN(6),TIME(2),
     3PROPS(24),DROT(3,3),POM(3,3),DEF(6),TAU(6),
     4DT(6,6),SD(6,2),ST(6),SDDT(6),RPLE(6),yumma(6),  
     5PropsUM(40), StatevUM(90)
C
C     Sve velicine u DATA01,DATA02,DATA03 i PUNGA se izracunavaju
C     u okviru UMAT-a, tj nisu ulazne velicine
C
      COMMON /DATA01/DSA(6,6),DSM(6,6),DA1(6,6),DA2(3),TEMPAL(MAXMTP)
      COMMON /DATA02/MAXMTP1,NELMTP,IMTP,ITERAB,KEYN,ISOLVE,KFLAG
      COMMON /DATA03/TIME1,DTIME1,TIME0,NUMMIN
      COMMON /PUNGA/ imli,MODEL
C
C     DODATi COMMONi
      
      COMMON /PRINCI/ PRINC(3)
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP

C
C     MOZE DA TREBA ZA VELIKE DEFORMACIJE
C     COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
c     Here we initialize PropsUM using Props
C     NPropsUM = 40
      Call Clear(PropsUM, 40)
      PropsUM(1) = Props(1) !IPHASE (1-austenit, 2-martensite)
      PropsUM(3) = Props(2) !MODEL (1-tanaka; 2-boyd, lagoudas; 3-liang, rogers)
      PropsUM(5) = Props(3) !TOL (convergence criterion tolerance)
      PropsUM(6) = Props(4) !xi0 (initial value for martensitic volume fraction)
      PropsUM(7) = Props(5) !NELMTP (number of integration points)
      Do i=1,12
         PropsUM(10+i) = Props(5+i) ! EA,EM,nu,alphaA,alphaM,H,rDsOA,rDsOM
      End Do
C
      Do i=1,4
         PropsUM(15+i) = Props(10+i) - 273.E+0 !Mos,Mof,Aos,Aof
      end do
C
      Do i=1, 6
         PropsUM(30+i) = Props(17+i) !epstr11,epstr22,epstr33,2epstr23,2epstr13,2epstr12
      End Do
C
      PropsUM(38) = Props(24) !FRULE - flag for the form of lambda 1-eq5, 2-eq7
C
      model = Int(Props(2)+0.1E0)
C
      if(model.eq.1) then
         PropsUM(23) = (0.01E0)/(Props(13)-Props(14))
         PropsUM(24) = (0.01E0)/(Props(11)-Props(12))
      else if(model.eq.2) then
         PropsUM(23) = -Props(16)*(Props(14)-Props(13))
         PropsUM(24) = -Props(17)*(Props(11)-Props(12))
      else
         PropsUM(23) = 3.14159265358979/(Props(14)-Props(13))
         PropsUM(24) = 3.14159265358979/(Props(11)-Props(12))
      end if
C
      PSI=STATEVT(2)
C
      CALL FORMDD(PROPSUM,DSM,DSA,DA1,DA2)
C
      CALL ELASTF(PSI,DSM,DSA,DT)
C
      IF(IRAC.EQ.2) RETURN
C      
C     DODELA VREDNOSTI ZA UKUPNO VREME
C
      TIME(2)=VREME
C
C     GLAVNE VREDNOSTI
C     INVERZNO LAMBDA
C
C      IVELIKE = 1
C      if (IVELIKE.eq.0) then
C      write (3,*) 'princ',(princ(i),i=1,3)
C      write (3,*) 'QP',((qp(i,j),j=1,3),i=1,3)
c      write (3,*) 'XJ',((XJ(i,j),j=1,3),i=1,3)
c      P1=1.D0/DSQRT(PRINC(1))
C      P2=1.D0/DSQRT(PRINC(2))
C      P3=1.D0/DSQRT(PRINC(3))
C     INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
C      CALL DIJAD(POM,QP,QP,P1,P2,P3)
C     TENZOR ROTACIJE R, a potrebno je izracunati prirastaj tenzora
C     rotacije u koraku
C     R = Fe * Ue**-1
c      CALL MNOZM1(DROT,XJ,POM,3,3,3)
C      else
C      do 30 i=1,3
C      do 20 j=1,3
C      if (i.eq.j) then
C      drot(i,j)=1.0D0
C      else
C      drot(i,j)=0.0D0
c      endif
c   20 continue
C   30 continue
c      endif
C
C---- Prirastaj deformacije u koraku
C
      CALL CLEAR(DSTRAN,6)
      CALL CLEAR(STRAN,6)
      CALL CLEAR(STRESS,6)
      CALL CLEAR(STATEV,100)
C
      do 100 I=1,6
  100    DSTRAN(I)=DEF(I)-DEFT(I)
C
      CALL MTASSIGN(DEFT,STRAN,6,1)
      CALL MTASSIGN(TAUT,STRESS,6,1)
C
      IF (TIME(2).EQ.1.0D0) THEN
         TEMPT=TGT
      ENDIF
C
      TEMP=TEMPT
C
c     initiating data only for the 1st material point of the 1st iter. 
c     of the first increment     
      IF(KFLAG.NE.920716)THEN
          KFLAG=920716
          imli=0
          IMTP=0    ! initializing the counter for material points
          MAXMTP1=MAXMTP
          NELMTP=INT(PROPSUM(7)+0.1)  ! given # of material points
C
          IF(NELMTP.GT.MAXMTP)THEN  ! incase max. # of mat. pts. is low
              WRITE(6,4000)MAXMTP,NELMTP
              WRITE(*,4000)MAXMTP,NELMTP
              STOP
          ENDIF       
C
          MODEL=INT(PROPSUM(3)+0.1D0)
C
          IF (MODEL.NE.1.AND.MODEL.NE.2.AND.MODEL.NE.3) THEN
              WRITE(6,*) 'THE MODEL NUMBER IS WRONG IN THE INPUT'
              WRITE(*,*) 'THE MODEL NUMBER IS WRONG IN THE INPUT'
              stop 'Wrong model number'
          ENDIF    
C
      CALL CLEAR(TEMPAL, MAXMTP)
C          
      CALL FORMDD(PROPSUM,DSM,DSA,DA1,DA2) ! calc. mat.ten.
C
      TIME1=-999999999.D0 ! initializing time and inc. of time
      DTIME1=TIME1
C
      ENDIF !KRAJ INICIJALIZACIJE ZA PRVU TACKU, ITERACIJU I INKREMENT
C
C     Ne postoji propsum(2) nigde pre ovoga
C
      IOPT=INT(PROPSUM(2)+0.1) ! info. about plane stress/strain 
C      
c     count the material point number
C
      IMTP=IMTP+1
      IF(IMTP.GT.NELMTP) IMTP=1 ! beginning of a new inc.
      REALTEMP=TGT-273.0D0
      TIME0=TIME(2)
C
      if(imtp.eq.1) imli=imli+1
C
C     Inicijalizacija STATEV promenljivih
      IF (TIME(2).EQ.1.0D0) THEN
C
         CALL CLEAR(STATEVT,100)
         STATEVT(1)=PROPSUM(1) !direction of transformation 1-forward 2-reverse
         STATEVT(2)=PROPSUM(6) !martensitic volume fraction
         STATEVT(17)=0.0D0
         STATEVT(26)=0.0D0
C
         TEMPAL(IMTP)=REALTEMP !temperatura u gausovim tackama za TIME(2)=0
C
         CALL MTASSIGN(PROPSUM(31),STATEVT(5),6,1) !Prepisuje iz PROPSUM u STATEV transformation strain
C
         DO 5 I=1,6
    5       STATEVT(19+I)=PROPSUM(10)*PROPSUM(20) !Transformation direction tensor
C
      ENDIF
C
C     IZRACUNAVANJE PRIRASTAJA TEMPERATURE
C
      REALTEMPT=TEMPT-273.0D0
C
      DTEMP=REALTEMP-REALTEMPT
C
C     PREPISIVANJE STATEVT IZ TRENUTKA T U TEKUCI VEKTOR STATEV
C
      CALL MTASSIGN(STATEVT,STATEV,100,1)                   
C
c     check if this is the first iteration cycle of an increment
C
      ISOLVE=3 ! if the same increment
      KEYN=0
      ITERAB=2
C
      IF(TIME1.NE.TIME(2))THEN  ! if new increment
          ITERAB=1
          KEYN=1
          TIME1=TIME(2)
          DTIME1=DTIME
      ELSEIF(DTIME.NE.DTIME1)THEN  ! if time increment changes
          ITERAB=1
          KEYN=2
          TIME1=TIME(2)
          DTIME1=DTIME
      ENDIF
C
c     converting the input quantities into second order tensorial form
C     transform 3/2-D stress and strain tensors into 1-D ones
c     Vraca ST(i)-stress, SD(i,1)-strain i SD(i,2)-incremet of strain
C
      CALL STD3_2(ST,SD,STRESS,STRAN,DSTRAN)
C
c     calculate stress increment for given strain and temperature
c     increments, and update the tangent stiffness matrix
C
      CALL SMA(ST,SD,STATEV,SDDT,RPLE,RPL,DRPLDT,PROPSUM,DT,
     *         REALTEMPT,DTEMP,DSA,DSM,DA1,DA2,TEMPAL(IMTP),DTIME,
     *         yumma,tumma,NOEL)
C
c     convert 3-D quantities into proper dimensions
C
      CALL STF3_2(DT,DDSDDE,STRESS,ST,STRAN,SD,DDSDDT,DRPLDE,SDDT,
     *            RPLE)

      STATEV(48)=DSQRT(2.0D0)/2.0D0*DSQRT((STRESS(1)-STRESS(2))**2
     *           +(STRESS(1)-STRESS(3))**2+(STRESS(2)-STRESS(3))**2)
      STATEV(49)=DSQRT(2.0D0)/3.0D0*DSQRT((STRAN(1)+DSTRAN(1)
     *              -STRAN(2)-DSTRAN(2))**2+(STRAN(1)+DSTRAN(1)
     *              -STRAN(3)-DSTRAN(3))**2+(STRAN(2)+DSTRAN(2)
     *              -STRAN(3)-DSTRAN(3))**2)
      STATEV(50)=(STRESS(1)+STRESS(2)+STRESS(3))/3.0D0
C
      CALL MTASSIGN(stran,def,6,1)
      CALL MTASSIGN(stress,tau,6,1)
      CALL MTASSIGN(stran,deftt,6,1)
      CALL MTASSIGN(stress,tautt,6,1)
C
      TEMPTT=REALTEMP+273.0D0
C
      CALL MTASSIGN(STATEV,STATEVTT,100,1)
C
 4000 FORMAT(//2X
     *' THE REAL MATERIAL POINTS ARE MORE THAN THE MAXIUM MATERIAL ',
     *' POINTS ALLOWED!'//2X,
     *' PLEASE ENLARGE THE PARAMETER, MAXMTP, IN SUBROUTINE, UMAT TO',
     *' ENSURE THAT IT IS LARGER THAN THE REAL MATERIAL POINTS.'//2X,
     *' THE MAXIUM MATERIAL POINTS ALLOWED NOW .... (MAXMTP)=',I8/2X,
     *' THE REAL MATERIAL POINTS INPUTTED. . . . . . . (NELMTP)=',I8)
  10  FORMAT(1X, I2, 7E12.5)
  11  FORMAT(1X, /)
C      
      RETURN
      END 
C      
c     ******************************************************************
c     *************************** END OF UMAT **************************
c     ******************************************************************
C
c     ****************************************************************** 
c     *********************** 3-D CASE STARTS HERE *********************
c     ******************************************************************
C
c     this subroutine is to perform an integration scheme to find the
c     real stresses under the given strain input, and then to update the
c     state variables, for a 2-d/3-d case.
C
      SUBROUTINE SMA(ST0,SD,STATEV,SDDT,RPLE,RPL,DRPLDT,PROPS,
     *               DT,TEMP,DTEMP,DSA,DSM,DA1,DA2,TEMPINI,DTIME,
     *               yumma,tumma,NOEL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      COMMON /DATA02/MAXMTP1,NELMTP,IMTP,ITERAB,KEYN,ISOLVE,KFLAG
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION ST0(6),SD(6,2),STATEV(100),PROPS(40),DT(6,6),
     * DSA(6,6),DSM(6,6),ST(6),SDT(6),RL(6),DA1(6,6),DA2(3),
     * SDDT(6),RPLE(6),yumma(6),FOLD(6),FNEW(6),DEV_ST(6),
     * STRAN(6),RLNEW(6)
C
c ********************* A R R A Y S ********************
C
c     ST0(1) : stress vector from MAIN 
c     ST(6)  : local stress vector 
c     DD(6)  : local total strain tensor 
c     RL(6)  : lamda vector 
c     SP(6)  : deviatoric stress vector 
c     SDT(6) : transformation strain vector 
c     SDDT(6): thermal stifness vector 
c     RPLE(6) :
C
      STR_FLG=0.00001D0
      IOPT=INT(PROPS(2)+0.1) ! info about plane stress/strain
      PSI=STATEV(2) ! martensitic volume fraction
      IPHASE=INT(STATEV(1)+0.1) ! info. about no yield/for./rev. trans.
      SEF0 =STATEV(13) ! von-mises stress
      ALPH=PROPS(20)
      BETA=PROPS(27)
      GAMA=PROPS(28)
C
      IRULE = INT(PROPS(38)+0.1)
C
      CALL MTASSIGN(STATEV(5),SDT,1,6) ! assigning state variables to SDT transformation strain vector
C
      NRL=INT(STATEV(17)+0.1D0)
C
      CALL MTASSIGN(STATEV(20),RL,1,6) ! local variables LAMBDE
C
C     ULAZIMO SA TEMPRATUROM IZ PRETHODNOG KORAKA
C
      TEMP0=TEMPINI
      DTEMP0=TEMP+DTEMP-TEMP0 ! start temp.-To
      TEMP1=TEMP+DTEMP ! present temperature
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of thermal (formula 2 u radu Qidwai)
                                          ! expansion
c     calculate the elastic predictor if no time inc. change
C     this subroutine finds the elastic stiffness matrix DT
C
      CALL ELASTF(PSI,DSM,DSA,DT) !(formula 2 u radu Qidwai racuna S(DT) od Sa i Sm)
C
      DO 15 I=1,3
         ST(I)=ST0(I) ! use local stress vector
         I3=I+3
         ST(I3)=ST0(I3)+DT(I3,I3)*SD(I3,2) ! update stress
         DO 15 J=1,3
 15         ST(I)=ST(I)+DT(I,J)*(SD(J,2)-ALPHA*DTEMP)
C
      IF(STRS_DOT(ST,ST).LT.STR_FLG) THEN
         DO 17 I=1,6
 17         RLNEW(I)=PROPS(10)*PROPS(20)         
      ELSE
         CALL GETRL(RLNEW,ST,SDT,ALPH,BETA,GAMA,IPHASE, IRULE)
      ENDIF
C
c     check yield of the material point
C
      CALL CLEAR(FOLD,6)
      CALL CLEAR(FNEW,6)
C
      DO 20 I=1,3
         I3=I+3
         FOLD(I3)=DA1(I3,I3)*ST0(I3)
         FNEW(I3)=DA1(I3,I3)*ST(I3)
         DO 20 J=1,3
            FOLD(I)=FOLD(I)+DA1(I,J)*ST0(J)
 20         FNEW(I)=FNEW(I)+DA1(I,J)*ST(J)
C      
      DOLD=DOT(ST0,RL,6)+0.5D0*DOT(FOLD,ST0,6)+DOT(ST0,DA2,3)
     * *(TEMP-TEMP0)+PROPS(22)*(TEMP-PROPS(16))
C
C
      DNEW=DOT(ST,RLNEW,6)+0.5D0*DOT(FNEW,ST,6)+DOT(ST,DA2,3)
     * *(TEMP1-TEMP0)+PROPS(22)*(TEMP1-PROPS(16))
C
      DECIDE=DNEW-DOLD
C      
C
      IF(DECIDE.GT.0.0D0)THEN
         IPHASE=1
         IFYD=1
      ELSEIF(DECIDE.LT.0.0D0)THEN
         IPHASE=2
         IFYD=1
      ELSE
         IFYD=0
         GOTO 70
      ENDIF
C
c     if yielding happened, return mapping integration scheme will 
c     be performed to find the real stress increment under the given 
c     strain and temperature increments and current state input.
C
      CALL BEINT1(ST0,SDT,PSI,SD,TEMP,DTEMP,TEMP0,PROPS,DSA,DSM,
     *            IPHASE,DA1,DA2,DT,SDDT,RPLE,RPL,DRPLDT,DTIME,
     *            yumma,tumma,RL,NRL,IFYD,NOEL)
C
c     update stresses, and all state vaiables
C
 70   STATEV(1)=DFLOAT(IPHASE)  ! whether forward/reverse trans.
      STATEV(2)=PSI
      STATEV(3)=DFLOAT(IFYD)  
C
      CALL MTASSIGN(SDT,STATEV(5),6,1) ! assign transform. strain 
C
      STATEV(11)=SD(2,1)+SD(2,2) ! 33-component of strain
      STATEV(12)=ST0(2)  ! 33-component of stress
C
      CALL DEV_STRS(ST0,DEV_ST)    
C
      SEF=STRS_DOT(DEV_ST,DEV_ST) ! von-mises stress
      STATEV(13)=DSQRT(1.5D0*SEF)
      STATEV(14)=DSQRT(2.D0/3.D0*STRN_DOT(SDT,SDT))! inelastic fiber strain
      STATEV(15)=TEMP1  ! assigning temp. 
C
c     get lamda(RL), effective stress(U)      
C                                   
      STATEV(16)=DOT(RL,ST0,6)/PROPS(20) ! modified effective stress
      STATEV(17)=DFLOAT(NRL)
      STATEV(26)=DTEMP0
C
      CALL MTASSIGN(RL,STATEV(20),6,1)
C
      DO 80 I=1,6
 80      STRAN(I)=SD(I,1)+SD(I,2)
C
      RETURN
      END   
C
c     ******************************************************************
C
c     this subroutine is to integrate the stress, inelastic strain, and
c     internal state variable, increment for gaven strain and temperature
c     increments by using Return Mapping Method (elastic predictor-plastic
c     corrector scheme), for a 2-d/3-d case.
C
      SUBROUTINE BEINT1(ST0,SDT0,PSI0,SD,TEMP,DTEMP,TEMP0,PROPS,
     *        DSA,DSM,IPHASE,DA1,DA2,DT,SDDT,RPLE,RPL,DRPLDT,DTIME,
     *           yumma,tumma,RL,NRL,IFYD,NOEL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /DATA02/MAXMTP1,NELMTP,IMTP,ITERAB,KEYN,ISOLVE,KFLAG
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION SD(6,2),ST0(6),SDT0(6),PROPS(40),DSA(6,6),DSM(6,6),
     *          DT(6,6),ST(6),SDT(6),DD1(6),DD(6),Q(6),RL(6),R(6),
     *          DA1(6,6),DA2(3),SDDT(6),RPLE(6),yumma(6),F(6) 
C
c     ********************* A  R   R   A   Y   S *********************** 
C
c     SDT0(6)     : transformation strain, actually SDT(6) from sma
c     DD1(6)      : local total strain vector 
c     DD(6)       : local strain inc. vector
c     Q(6)        : see papers
c     R(6)        : see papers
c     D(6,6)      : elastic stiffness matrix
c     ST(6)       : local stress vector
c     SDT(6)      : local transformation strain vector
c     RL(6)       : lamda vector, see papers
c     ES1(6)      : elastic strain at previous increment
c     Q1(6)       : see papers
C
      TOL=PROPS(5)              ! criteria of convergence
      DTMP=DTEMP                ! temp. inc. for each inc. within iter.
      TMP=TEMP                  ! starting temp.
      PSI=PSI0                  ! martensitic volume frac.
      ALPH=PROPS(20)            
      BETA=PROPS(27)                
      GAMA=PROPS(28)
C
      IRULE = INT(PROPS(38)+0.1)
C
      DO 10 I=1,6           ! assigining stress, strain and strain inc. 
          ST(I)=ST0(I)      ! values from sma to local vectors
          SDT(I)=SDT0(I)
          DD1(I)=SD(I,1)
 10       DD(I)=SD(I,2)
C
      FLAG1F=1.0D0
      FLAG1R=0.0D0
      FLAG2R=0.001D0
      FLAG2F=0.99D0
C
      if(psi.le.0.0001D0) then
        iphase=1
      end if
C
c     to implement the return mapping integration scheme, calculate 
c     elastic predictor
C
      TMP=TMP+DTMP          ! present temp.
      DTMP0=TMP-TEMP0           ! present temp-To
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C
      CALL ELASTF(PSI,DSM,DSA,DT) !get elas. stiff. mat.     
      CALL ELATTF(DT,ALPHA,SDDT) 
C
      DO 20 I=1,6               ! updating strain and stress
         DD1(I)=DD1(I)+DD(I)
         DO 20 J=1,6
 20         ST(I)=ST(I)+DT(I,J)*DD(J)
C
      DO 30 I=1,3     ! adding the thermal part
 30      ST(I)=ST(I)+SDDT(I)*DTMP   
C
c     iteration for plastic corrector
C
      DO 50 LOOP=1,100
C
c     calculate rl and the value of the yield function
C
      CALL GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE,IRULE)
C
         CALL YDFUN(FYD,ST,SDT,RL,TMP,DTMP0,PSI,PROPS,IPHASE,
     *              DA1,DA2,F) 
C
         IF(LOOP.EQ.1)THEN
            IF(FYD.LE.0.0D0)THEN
               IFYD=0
               GOTO 9040  
            ELSEIF(IPHASE.EQ.1.AND.PSI.GE.(FLAG1F-0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO 9040
            ELSEIF(IPHASE.EQ.2.AND.PSI.LE.(FLAG1R+0.0001).AND.
     *              DABS(FYD).GT.0.0D0)THEN
               IFYD=0
               GOTO 9040
            ENDIF
         ENDIF
C         
c     get rl, and b, also see papers/scheme              
C
         CALL GETQ(Q,DA1,DA2,DTMP0,RL,ST)
C
         CALL GETB(B,PSI,Q,PROPS,DT,IPHASE)
C
c     update PSI, SDT, and ST     
C
         DPSI=-FYD/B            ! inc. of mart. vol. frac.
         psii=psi
         PSI=PSI+DPSI           ! updating mart. vol. frac.
C
         if(psi.gt.FLAG1F)then
            if(iphase.eq.1)then
               psi=FLAG1F
               dpsi=psi-psii
            elseif(iphase.eq.2.and.psii.lt.FLAG1F)then
               psi=FLAG2R
               dpsi=psi-psii
            endif
         elseif(psi.lt.FLAG1R)then
            if(iphase.eq.1.and.psii.gt.FLAG1R)then
               psi=FLAG2F
               dpsi=psi-psii
            elseif(iphase.eq.2)then
               psi=FLAG1R
               dpsi=psi-psii
            endif
         endif
C
         DO 40 I=1,6            ! see scheme
 40         SDT(I)=SDT(I)+RL(I)*DPSI 
C
         CALL ELASTF(PSI,DSM,DSA,DT)
         ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C
         CALL CLEAR(ST,6)
C
c         pause"is the problem here?"
C
         DO 45 I=1,3
            I3=I+3
            ST(I3)=DT(I3,I3)*(DD1(I3)-(SDT(I3)-PROPS(30+I3)))
            DO 45 J=1,3
 45            ST(I)=ST(I)+DT(I,J)*(DD1(J)-ALPHA*DTMP0
     *                -(SDT(J)-PROPS(30+J)))
C
c     check convergence
C
         IF(DABS(DPSI).LT.TOL) GOTO 9020
C
 50   CONTINUE
C
      WRITE(*,*)'ITERATION FAILS TO CONVERGE!Del_xi=',DPSI,'xi=',PSI
C
 9020 CONTINUE
  100 CONTINUE  
C
c     update tangent stifness matrix
C
      if(iphase.eq.1.and.psi.lt.FLAG1F
     *     .or.iphase.eq.2.and.psi.gt.FLAG1R)then
      CALL TANSTF(DT,SDDT,PSI,ST,SDT,DTMP0,DSA,DSM,PROPS,IPHASE,DA1,
     *            DA2,RL)
      endif
C
      if(psi.le.0.0001D0) then
        iphase=1
      end if
C
      GOTO 9040   ! return
C
 9030 CONTINUE  
C
c     if PSI is out of bound, forward Euler integration scheme is applied. 
C
      IF(PSI.LT.1.0E-3.AND.IPHASE.EQ.2.OR.PSI.GT.99.99E-2.AND.
     *     IPHASE.EQ.1)THEN
         NSS=100
         DTMP=DTEMP/DFLOAT(NSS)
         DTMP0=TEMP-TEMP0
         PSI=PSI0
         DO 110 I=1,6           ! assigning local variables with global values
            ST(I)=ST0(I)
            SDT(I)=SDT0(I)
            DD1(I)=SD(I,1)
 110        DD(I)=SD(I,2)/DFLOAT(NSS) 
C
         CALL ELASTF(PSI,DSM,DSA,DT)    
C
         DO 200 N=1,NSS 
            DTMP1=DTMP0+DTMP 
C
            CALL GETRL(RL,ST,SDT,H,B2,IPHASE,IRULE)
            CALL GETQ(Q,DA1,DA2,DTMP1,RL,ST)
            CALL GETB(B,PSI,Q,PROPS,DT,IPHASE)
            CALL GETRS(R,S,Q,DT,DA2,ST,PROPS,IPHASE)
C
            DPSI=-(DOT(R,DD,6)+S*DTMP)/B
            PSI=PSI+DPSI      
C
            IF(PSI.GT.1.0D-5.AND.PSI.LT.0.999999D0)THEN ! if not near the 
               DTMP0=DTMP0+DTMP ! finish
               DO 120 I=1,6
                  DD1(I)=DD1(I)+DD(I) ! strain
 120              SDT(I)=SDT(I)+RL(I)*DPSI ! transformation strain
C
               ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C
               CALL ELASTF(PSI,DSM,DSA,DT) ! elas. stiff. mat. 
               CALL ELATTF(DT,ALPHA,SDDT)
               CALL CLEAR(ST,6)
C
               DO 130 I=1,3     ! stress
                  I3=I+3
                  ST(I3)=DT(I3,I3)*(DD1(I3)-(SDT(I3)-PROPS(30+I3)))
                  DO 130 J=1,3
 130                 ST(I)=ST(I)+DT(I,J)*(DD1(J)-ALPHA*DTMP0
     *                     -(SDT(J)-PROPS(30+J)))
            ELSE       ! near the finish and is done elastic
               DO 140 I=1,6
 140              DD(I)=SD(I,1)+SD(I,2)-DD1(I) ! strain inc.
C
               PSI=PSI-DPSI ! why
C
               ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C
               CALL ELASTF(PSI,DSM,DSA,DT) ! elastic stiffness
               CALL ELATTF(DT,ALPHA,SDDT)
C
               DO 150 I=1,3     ! stress
                  I3=I+3
                  ST(I3)=ST(I3)+DT(I3,I3)*DD(I3)
                  ST(I )=ST(I )-(DT(I,1)+DT(I,2)+DT(I,3))*ALPHA
     *                   *(TEMP+DTEMP-DTMP0-TEMP0)
                  DO 150 J=1,3
 150                 ST(I)=ST(I)+DT(I,J)*DD(J)
C
               GOTO 9040
            ENDIF      
 200     CONTINUE
      ELSE
         WRITE(*  ,*)'IPHASE-PSI WRONG!  IPHASE,PSI=',IPHASE,PSI
         WRITE(MSG,*)'IPHASE-PSI WRONG!  IPHASE,PSI=',IPHASE,PSI
      ENDIF
C
 9040 CONTINUE 
C
      IF(INT(PROPS(8)+0.1D0).EQ.1.AND.IFYD.EQ.1)THEN
         if(iphase.eq.1.and.psi.lt.(FLAG1F-0.0001d0)
     *        .or.iphase.eq.2.and.psi.gt.(FLAG1R+0.0001d0))then         
            T=TEMP+DTEMP        ! current temperature
            DTEMP0=TEMP+DTEMP-TEMP0 ! current temperature-reference temp.
            DPSI=PSI-PSI0
C
c     get lamda(RL), and then Q      
C
            CALL GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE,IRULE) ! model 1,3
            CALL GETQ(Q,DA1,DA2,DTEMP0,RL,ST)
C
c     getting the variables necessary for coupled heat transfer analysis
C
            CALL CALCRPL(RPL,RPLE,DRPLDT,PROPS,IPHASE,T,SD,ST,st0,DT,
     *                   DA2,SDDT,Q,DTEMP,DTIME,DPSI,PSI,yumma,
     *                   tumma)
         else
            CALL CLEAR(RPLE,6)
            
            RPL=0.0D0
            DRPLDT=0.0D0
         endif
      ELSE
         CALL CLEAR(RPLE,6)
C        
         RPL=0.0D0
         DRPLDT=0.0D0          
      ENDIF
C
c     assign local state variables to global ones
C
      PSI0=PSI
      DO 90 I=1,6
         SD(I,1)=DD1(I)
         SD(I,2)=DD(I)
         SDT0(I)=SDT(I)         ! transformation strain
 90      ST0(I)=ST(I)           ! stress
C      
      RETURN
      END   
C
c     ******************************************************************
C
c     this subroutine is to get the tangent STIFFNESS matrix DT, and
c     thermal STIFFNESS coefficients, ALS.
C
      SUBROUTINE TANSTF(DT,SDDT,PSI,ST,SDT,DTEMP0,DSA,DSM,PROPS,
     *                  IPHASE,DA1,DA2,RL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION DT(6,6),SDDT(6),ST(6),SDT(6),DSM(6,6),DSA(6,6),
     *          TEMP(6),PROPS(40),R(6),Q(6),RL(6),DA1(6,6),DA2(3)   
     
C
c     *********************  A    R   R   A   Y   S ********************
C
c     R(6)        : see scheme
c     Q(6)        : see scheme
c     RL(6)       : lamda vector          

      CALL CLEAR(SDDT,6)
      CALL CLEAR(TEMP,6)
C
      ALPH=PROPS(20)
      BETA=PROPS(27)    
      GAMA=PROPS(28)
C
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI) ! coeff. of ther. exp.
C
c     get q, r, s, and b, also see scheme                        
C
      CALL GETQ(Q,DA1,DA2,DTEMP0,RL,ST)
      CALL ELASTF(PSI,DSM,DSA,DT)
      CALL GETRS(R,DS,Q,DT,DA2,ST,PROPS,IPHASE)
      CALL GETB(B,PSI,Q,PROPS,DT,IPHASE)
C
      DO 10 I=1,3               ! thermal stiffness vector
         I3=I+3
         SDDT(I+3)=DT(I3,I3)*Q(I3)*DS/B
         DO 10 J=1,3
 10         SDDT(I)=SDDT(I)+DT(I,J)*(Q(J)*DS/B-ALPHA)
C
      DO 15 I=1,3
         I3=I+3
         TEMP(I3)=DT(I3,I3)*Q(I3)
         DO 15 J=1,3
 15         TEMP(I)=TEMP(I)+DT(I,J)*Q(J)                
C
      DO 20 I=1,6
         DO 20 J=1,6
 20         DT(I,J)=DT(I,J)+TEMP(I)*R(J)/B ! tangent stiffness matrix
C
      RETURN
      END   
C
c     ******************************************************************
C
c     this subroutine calculates q for a 2-d/3-d case
c     also see scheme.           
C
      SUBROUTINE GETQ(Q,DA1,DA2,DTEMP0,RL,ST) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION Q(6),ST(6),DA1(6,6),DA2(3),RL(6)
C
      CALL CLEAR(Q,6)   
C
      DO 10 I=1,3
      I3=I+3
      Q(I3)=DA1(I3,I3)*ST(I3)+RL(I3)
          Q(I)=DA2(I)*DTEMP0+RL(I)
      DO 10 J=1,3
 10       Q(I)=Q(I)+DA1(I,J)*ST(J)
C
      RETURN   
      END   
C
c     ************************************************************
C
c     this subroutine calculates b for a 2-d/3-d case
c     also see scheme.           
C
      SUBROUTINE GETB(B,PSI0,Q,PROPS,D,IPHASE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION Q(6),PROPS(40),D(6,6),TEMP(6)
C
      PSI=PSI0  
C
      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007
C
      CALL CLEAR(TEMP,6)    
C
      DO 10 I=1,3
      I3=I+3
      TEMP(I3)=Q(I3)*D(I3,I3)
      DO 10 J=1,3   
 10       TEMP(I)=TEMP(I)+Q(J)*D(J,I)       
C  
c     see scheme again                       
C                 
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
          IF (MODEL.EQ.1) THEN  ! model 1
              B=-DOT(Q,TEMP,6)-PROPS(22)/PROPS(24)/(1.0D0-PSI)
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B=-DOT(Q,TEMP,6)-(2.0D0*PROPS(26)*PSI+PROPS(24))
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B=-DOT(Q,TEMP,6)+2.0*PROPS(22)/PROPS(24)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)               
          ENDIF        
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
          IF (MODEL.EQ.1) THEN  ! model 1
              B=DOT(Q,TEMP,6)-PROPS(21)/PROPS(23)/PSI
          ELSEIF (MODEL.EQ.2) THEN  ! model 2
              B=DOT(Q,TEMP,6)+(2.0D0*PROPS(25)*PSI+PROPS(23))    
          ELSEIF (MODEL.EQ.3) THEN  ! model 3
              B=DOT(Q,TEMP,6)-2.0*PROPS(21)/PROPS(23)
     *          /DSQRT(1.0D0-(2.0*PSI-1.0D0)**2)    
          ENDIF  
      ENDIF
C
      RETURN
      END
C
c***********************************************************************
C
c     this subroutine calculates r and s for a 2-d/3-d case
c     also see scheme. 
C
C
      SUBROUTINE GETRS(R,DS,Q,D,DA2,ST,PROPS,IPHASE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION R(6),Q(6),D(6,6),DA2(3),ST(6),PROPS(40)
C
      CALL CLEAR(R,6)   
C
      DO 20 I=1,3
      I3=I+3
      R(I3)=Q(I3)*D(I3,I3)
      IF(IPHASE.EQ.2) R(I3)=-R(I3)
      DO 10 J=1,3
 10       R(I)=R(I)+Q(J)*D(J,I)
 20   IF(IPHASE.EQ.2) R(I)=-R(I)        
C
        ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C       
c     see scheme again                       
C                       
      IF(IPHASE.EQ.1)THEN  ! forward transforamtion 
      DS=DOT(DA2,ST,3)+PROPS(22)-(R(1)+R(2)+R(3))*ALPHA
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation 
      DS=-DOT(DA2,ST,3)-PROPS(21)-(R(1)+R(2)+R(3))*ALPHA
      ENDIF
C
      RETURN
      END
C
c***********************************************************************
C
c     this subroutine calculates lamda(RL) and effective stress(U) for a
c     2-d/3-d case
C
      SUBROUTINE GETRL(RL,ST,SDT,ALPH,BETA,GAMA,IPHASE, IRULE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)              
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION RL(6),ST(6),SDT(6),DST(6),DSTINV(3,3),DSTTEMP(6) 
C      
c     *********************** A   R   R   A   Y   s ********************
C
c     SP(6)       : deviatoric stress vector
C
      IF((IPHASE.EQ.2).AND.(IRULE.NE.2))THEN
         DEN=DSQRT(2.0D0/3.0D0*STRN_DOT(SDT,SDT))
C
         IF(DEN.EQ.0.0D0)THEN
            DO 10 I=1,6
 10            RL(I)=ALPH
         ELSE
            DEN=ALPH/DEN
C
            DO 20 I=1,6
 20            RL(I)=DEN*SDT(I)
         ENDIF            
      ELSE
         CALL DEV_STRS(ST,DST)
C
         DJ2=0.5D0*STRS_DOT(DST,DST)
C
         IF(DJ2.EQ.0.D0)THEN    ! anomaly
            DO 30 I=1,6
 30            RL(I)=1.0D0
         ELSE
            DEN=ALPH*(3.0d0/2.0D0)/DSQRT(3.0D0*DJ2)
C
            DO 40 I=1,3
               RL(I)=DEN*DST(I)+GAMA
 40            RL(I+3)=2.0D0*DEN*DST(I+3)
         ENDIF
      ENDIF
C
      RETURN
      END
C
c     ******************************************************************  
C
c     calculates the value of the yield function for a material point of 
c     a 2-d/3-d case.
C
      SUBROUTINE YDFUN(FYD,ST,SDT,RL,TEMP1,DTEMP1,PSI0,PROPS,IPHASE,
     *                 DA1,DA2,F)   
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION ST(6),SDT(6),PROPS(40),RL(6),DA1(6,6),
     *      DA2(3),F(6) 
C
c     ********************* A   R   R   A   Y   S **********************
C
c     Rl(6)       : local lamda vector      
C
      FYD=0.D0  ! yield value of the function
      PSI=PSI0  ! martensitic volume fraction
C
      IF(PSI.GT.0.99993) PSI=0.99993
      IF(PSI.LT.0.00007) PSI=0.00007
C
      ALPH=PROPS(20)  ! max. transformation strain
      BETA=PROPS(27)
      GAMA=PROPS(28)
C      
c     calculate lamda vector and effective stress      
C                                                      
      CALL CLEAR(F,6)
C
      DO 10 I=1,3
      I3=I+3
          F(I3)=DA1(I3,I3)*ST(I3)
          DO 10 J=1,3
  10          F(I)=F(I)+DA1(I,J)*ST(J) 
C
      IF(IPHASE.EQ.1)THEN   ! forward transforamtion
          IF (MODEL.EQ.1) THEN  ! model 1
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*TEMP1+PROPS(22)/PROPS(24)
     *        *DLOG(1.D0-PSI)-PROPS(22)*PROPS(16)
          ELSEIF (MODEL.EQ.2) THEN   ! model 2 
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*(TEMP1-PROPS(16))-(PROPS(26)*PSI
     *            +PROPS(24))*PSI
          ELSEIF (MODEL.EQ.3) THEN   ! model 3 
              FYD=DOT(ST,RL,6)+(1.0D0/2.0D0)*DOT(F,ST,6)+DOT(DA2,ST,3)
     *            *DTEMP1+PROPS(22)*TEMP1-PROPS(22)/PROPS(24)
     *            *(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)-PROPS(22)
     *        *PROPS(16)
          ENDIF               
      ELSEIF(IPHASE.EQ.2)THEN  ! reverse transformation
          IF (MODEL.EQ.1) THEN  ! model 1
             FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*TEMP1-PROPS(21)/PROPS(23)
     *        *DLOG(PSI)+PROPS(21)*PROPS(18)
        ELSEIF (MODEL.EQ.2) THEN   ! model 2 
             FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*(TEMP1-PROPS(19))+(PROPS(25)*PSI
     *            +PROPS(23))*PSI
          ELSEIF (MODEL.EQ.3) THEN   ! model 3
             FYD=-DOT(ST,RL,6)-(1.0D0/2.0D0)*DOT(F,ST,6)-DOT(DA2,ST,3)
     *            *DTEMP1-PROPS(21)*TEMP1+PROPS(21)/PROPS(23)
     *        *(DACOS(2.0D0*PSI-1.0D0)-22.0D0/7.0D0)+PROPS(21)
     *        *PROPS(19) 
          ENDIF
      ENDIF
C                                                    
      RETURN
      END
C
c     ******************************************************************
C
c     this subroutine finds the elastic stiffness matrix
C
      SUBROUTINE ELASTF(PSI,DSM,DSA,D) 
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
C      
      DIMENSION DSM(6,6),DSA(6,6),D(6,6),TEMP(3,3)
C
      DO 10 I=1,6
          DO 10 J=1,6
 10           D(I,J)=DSA(I,J)+PSI*(DSM(I,J)-DSA(I,J))
C
      DO 20 I=1,3
      DO 20 J=1,3
 20       TEMP(I,J)=D(I,J)
C
      CALL GET_INV_DET(TEMP,GARB,0)
C
      DO 30 I=1,3
      DO 30 J=1,3
 30       D(I,J)=TEMP(I,J)
C    
      D(4,4)=1.0D0/D(4,4)           
      D(5,5)=1.0D0/D(5,5)           
      D(6,6)=1.0D0/D(6,6)           
C           
      RETURN
      END
C
c     ******************************************************************
c     this subroutine finds the elastic thermal stiffness matrix
C
      SUBROUTINE ELATTF(D,ALPHA,TST)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION D(6,6),TST(6)
C
      CALL CLEAR(TST,6)
C
      DO 10 I=1,3
 10   TST(I)=-ALPHA*(D(I,1)+D(I,2)+D(I,3))
C
      RETURN
      END
C        
c***********************************************************************        
c     this subroutine is used to calculate the variables associated with
c     the coupled mechanical-thermal analysis
C
      SUBROUTINE CALCRPL(RPL,RPLE,DRPLDT,PROPS,IPHASE,T,SD,ST,st0,DT,
     *                   DA2,SDDT,Q,DTEMP,DTIME,DPSI,PSI,yumma,
     *                   tumma)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION RPLE(6),PROPS(40),ST(6),SD(6,2),DT(6,6),SDDT(6),
     *          G(6),DA2(3),Q(6),DRPE(6),yumma(6),st0(6)
C
      CALL CLEAR(DRPE,6)    
      RPL=0.0D0
      TEMP=T
C
c     calling the subroutine to calculate the a and g vectors beside h
C
      CALL RPLEXT(G,HT,DA2,TEMP,ST,IPHASE,PROPS,Q,DPSI,RPLE,
     *            DRPLDT,DT,PSI,SDDT,st0,rpl)
C
      DO 10 I=1,6
         RPLE(I)=RPLE(I)/DTIME  
         DO 10 J=1,6
 10         DRPE(I)=DRPE(I)+G(J)*DT(J,I)
C
      DRPT=DOT(G,SDDT,6)+HT 
C
      DRPLDT=DRPLDT/DTIME   
      RPL=RPL/DTIME 
C
      do 40 i=1,6
 40      yumma(i)=drpe(i)
C
      tumma=drpt
C
      RETURN
      END
C
c     ******************************************************************
      SUBROUTINE RPLEXT(G,HT,DA2,T,ST,IPHASE,PROPS,Q,DPSI,RPLE,
     *              DRPLDT,DT,PSI,SDDT,st0,rpl)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION A(6),G(6),DA2(3),ST(6),PROPS(40),Q(6),RPLE(6),
     *        DT(6,6),SDDT(6),st0(6),dst(6)         
C
c     only 2nd model now
C
      CALL CLEAR(RPLE,6)    
      CALL CLEAR(A,6)
      CALL CLEAR(G,6)
C
      T=T+273.0D0
      DALPH=DA2(1)
      ALPHA=ALFA(PROPS(14),PROPS(15),PSI)
C
      IF (IPHASE.EQ.1) THEN   ! forward transformation
         CALL MTASSIGN(Q,A,6,1) 
         B=DOT(DA2,ST,3)+PROPS(22)
         C=-(2.0D0*PROPS(26)*PSI+PROPS(24))
         DMU2=1.0D0/6.0D0*(PROPS(25)-PROPS(26))+1.0D0/4.0D0
     *        *(PROPS(23)-PROPS(24))
         FYD=-1.0D0/2.0D0*PROPS(22)*(PROPS(19)-PROPS(16))
     *        -DMU2
         HELLO=FYD-T*DALPH*(ST(1)+ST(2)+ST(3))-PROPS(22)*T
         YELLO=PROPS(22)*T-FYD
         DELLO=PROPS(22)
      ELSEIF(IPHASE.EQ.2) THEN  ! reverse transformation
         DO 15 I=1,6
 15         A(I)=-Q(I)      
         B=-DOT(DA2,ST,3)-PROPS(21)
         C=2.0D0*PROPS(25)*PSI+PROPS(23)
         DMU2=1.0D0/6.0D0*(PROPS(25)-PROPS(26))+1.0D0/4.0D0
     *        *(PROPS(23)-PROPS(24))
         FYD=1.0D0/2.0D0*PROPS(21)*(PROPS(19)-PROPS(16))
     *        +DMU2
         HELLO=FYD-T*DALPH*(ST(1)+ST(2)+ST(3))-PROPS(21)*T
         YELLO=PROPS(21)*T-FYD
         DELLO=PROPS(21)
      ENDIF     
C
      DO 20 I=1,6
         DO 20 J=1,6
 20         RPLE(I)=RPLE(I)+A(J)*DT(J,I)
C               
      DRPLDT=YELLO/C*(DOT(A,SDDT,6)+B)-DELLO*DPSI               
C
      DO 25 I=1,3
         RPLE(I)=RPLE(I)*YELLO/C    
         G(I)=-(T*ALPHA+A(I)*HELLO/C)
         RPLE(I+3)=RPLE(I+3)*YELLO/C    
 25      G(I+3)=-A(I+3)*HELLO/C
C
      HT=-B*HELLO/C
C
      do 26 i=1,6
 26      dst(i)=st(i)-st0(i)
C
      rpl=hello*dpsi-t*alpha*(dst(1)+dst(2)+dst(3))
C
      RETURN    
      END 
C
c     ******************************************************************
c     this subroutine is to calculate the deviatoric stress tensor, SP
c     and the effective stress, SEF.
C
      SUBROUTINE DEV_STRS(ST,SP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION ST(6),SP(6)
C
      S1=(ST(1)+ST(2)+ST(3))/3.D0  ! calculating diagonal elements
      SP(1)=ST(1)-S1
      SP(2)=ST(2)-S1
      SP(3)=ST(3)-S1
C      
      DO 10 I=1,3
 10       SP(I+3)=ST(I+3) ! calculating off-diagonal elements
C
      RETURN
      END
C
c     ******************************************************************
c     this function calculates the product of strain components
C
      FUNCTION STRN_DOT(Q,S) 
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      DIMENSION Q(6),S(6)
C
      D=0.D0
C      
      DO 10 I=1,3
 10       D=D+Q(I)*S(I)+0.5D0*Q(I+3)*S(I+3) 
C 
      STRN_DOT=D 
C      
      RETURN
      END   
C      
c     ******************************************************************
c     this function calculates the product of stress components      
C
      FUNCTION STRS_DOT(Q,S)    
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      DIMENSION Q(6),S(6)
C
      D=0.D0 
C      
      DO 10 I=1,3
 10       D=D+Q(I)*S(I)+2.0D0*Q(I+3)*S(I+3)   
C 
      STRS_DOT=D    
C
      RETURN
      END   
C      
c     ******************************************************************      
c   this subroutine inverts a given matrix and return it  
C        
      SUBROUTINE GET_INV_DET(XS,DET,NFLAG)
C                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                         
C                 
      DIMENSION XS(3,3),A(3,3) 
C    
      DO 10 I=1,3                                                       
      I1=I+1                                                            
      I2=I+2                                                            
      IF(I1.GT.3) I1=I1-3                                               
      IF(I2.GT.3) I2=I2-3                                               
      DO 10 J=1,3                                                       
          J1=J+1                                                            
          J2=J+2                                                            
          IF(J1.GT.3) J1=J1-3                                               
          IF(J2.GT.3) J2=J2-3                                               
   10         A(I,J)=XS(I1,J1)*XS(I2,J2)-XS(I1,J2)*XS(I2,J1)
C                      
      DET=0.D0 
C                                                             
      DO 20 I=1,3                                                       
   20     DET=DET+A(1,I)*XS(1,I) 
C   
      NE=1  
C
      if (NFLAG.EQ.0.AND.det.le.0.0d0) then
      write(*,*) 'error in the jacobain, DETERMINANT', det
      pause'error in the jacobian' 
      endif         
C                                           
      DO 30 I=1,3                                                      
      DO 30 J=1,3                                                       
   30         XS(I,J)=A(J,I)/DET 
      RETURN                                                            
      END           
C
c***********************************************************************
C     
c     this subroutine is to initiate the material data for a 2-3/3-d case.      
C
      SUBROUTINE FORMDD(PROPS,DSM,DSA,DA1,DA2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
C      
      COMMON /PUNGA/ imli,MODEL
C
      DIMENSION PROPS(40),DSM(6,6),DSA(6,6),DA1(6,6),DA2(3)
C
      CALL CLEAR(DSA,36) ! initializing 
      CALL CLEAR(DSM,36)
      CALL CLEAR(DA1,36) ! initializing
      CALL CLEAR(DA2,3)  ! initializing
C      
      V=PROPS(13)  ! poisson's ratio
      EA=PROPS(11) ! young's modulii
      EM=PROPS(12) 
C      
c     working on martenstic material matrix, see any mechanics book for
c     isotropic material matrix      
C
      DO 30 I=1,3  ! diagonal terms
      I3=I+3
      DSA(I,I)=1.0D0/EA
      DSM(I,I)=1.0D0/EM
      DA1(I,I)=DSM(I,I)-DSA(I,I)
C
      DSA(I3,I3)=2.0D0*(1+V)/EA
      DSM(I3,I3)=2.0D0*(1+V)/EM
 30   DA1(I3,I3)=DSM(I3,I3)-DSA(I3,I3)
C                             
C     off-diagonal terms                             
C                             
      DSA(1,2)=-V/EA
      DSM(1,2)=-V/EM
      DA1(1,2)=DSM(1,2)-DSA(1,2)
C
      DSA(1,3)=DSA(1,2)
      DSA(2,1)=DSA(1,2)
      DSA(2,3)=DSA(1,2)
      DSA(3,1)=DSA(1,2)
      DSA(3,2)=DSA(1,2)
C
      DSM(1,3)=DSM(1,2)
      DSM(2,1)=DSM(1,2)
      DSM(2,3)=DSM(1,2)
      DSM(3,1)=DSM(1,2)
      DSM(3,2)=DSM(1,2)
C
      DA1(1,3)=DA1(1,2)
      DA1(2,1)=DA1(1,2)
      DA1(2,3)=DA1(1,2)
      DA1(3,1)=DA1(1,2)
      DA1(3,2)=DA1(1,2)
C                              
      DO 40 I=1,3
 40   DA2(I)=PROPS(15)-PROPS(14)
      RETURN
      END   
C      
c     ******************************************************************
C
c     this subroutine is to transform 3/2-D stress and strain tensors 
c     into 1-D ones.      
C
      SUBROUTINE STD3_2(ST,SD,STRESS,STRAIN,DSTRAN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C      
      DIMENSION ST(6),SD(6,2),STRESS(6),STRAIN(6),DSTRAN(6),
     *      STATEV(100)
C                                               
      CALL CLEAR(ST,6) ! initializing
      CALL CLEAR(SD,12)  
C                                                
      DO 10 I=1,6
         ST(I)=STRESS(I)     ! stress
         SD(I,1)=STRAIN(I)   ! strain
 10      SD(I,2)=DSTRAN(I)   ! inc. of strain
      RETURN
      END   
C      
c     ******************************************************************
C
c     this subroutine is to convert 1-D stiffness matrix and stress
c     tensor to 3/2-D ones.

      SUBROUTINE STF3_2(D3,D2,STRESS,ST,STRAIN,SD,DDSDDT,DRPLDE,SDDT,
     *                  RPLE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C      
      DIMENSION D3(6,6),D2(6,6),STRESS(6),ST(6),
     *          DDSDDT(6),DRPLDE(6),SDDT(6),RPLE(6),SD(6,2),
     *          STRAIN(6),DSTRAN(6)
C
c     assigning local strain and strain inc. components to global 
c     variables      
C      
      CALL MTASSIGN(D3,D2,6,6)  
C                   
c     doing the same thing for stress                   
C 
      DO 30 I=1,6
         DDSDDT(I)=SDDT(I)
         DRPLDE(I)=RPLE(I)
         STRAIN(I)=SD(I,1)   ! strain
         DSTRAN(I)=SD(I,2)   ! inc. of strain
 30      STRESS(I)=ST(I)
      RETURN
      END  
C      
c     ******************************************************************
C
c     this subroutine assigns matrix a's values to matrix b      
C
      SUBROUTINE MTASSIGN(A,B,M,N)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
C      
      DIMENSION A(M,N),B(M,N)
C
      DO 10 I=1,M
          DO 10 J=1,N
 10           B(I,J)=A(I,J) 
C     
      RETURN
      END    
C
c     ******************************************************************      
C
c     this function calculates the present thermal expansion coefficient
C
      FUNCTION ALFA(ALPHAA,ALPHAM,PSI)
C   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    
      ALFA=ALPHAA+PSI*(ALPHAM-ALPHAA)
C    
      RETURN
      END 
C
c     ******************************************************************
c     ******************************************************************