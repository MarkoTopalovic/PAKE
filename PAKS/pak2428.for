C     3D OR PLANE STRAIN NEARLY INCOMPRESIBLE NEO-HOOKEAN
C
C
      SUBROUTINE D2M28 (TAU,DEF,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE

      include 'paka.inc'
      
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D2M28'
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6*IDVA
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6*IDVA
C
      CALL TI228(A(LTAU),A(LDEFT),
     1           A(LTAU1),A(LDEFT1),
     1           A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI228(STREST,STRAIT,
     1                 STRES1,STRAI1,
     1                 FUN,MATE,STRESS,STRAIN,IRAC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     MATERIJALNI MMODEL GUME (MOONEY-RIVLIN)
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION ctens(3,3,3,3),sigma(3,3),STRESS(*),STRAIN(*),btens(3,3)
      DIMENSION FUN(3,MATE)
C
      DIMENSION STREST(*),STRAIT(*),
     1          STRES1(*),STRAI1(*)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI228'
C
      IF(IRAC.EQ.2) RETURN
      call clear(ctens,81)
      call clear(sigma,9)
C     MATERIJALNE KONSTANTE
C     PRVA JE MI_10,DRUGA JE MI_01 A TRECA JE K
      S1=FUN(1,MAT)
      S2=FUN(2,MAT)
      S3=FUN(3,MAT)
      call jedna1(btens,grae,9)
!     Material 5 nearly incompressible neo-Hookean
!
      ndime=2
      detf=detg
	ak=2*(S1+S2)/(1-2*S3)
      call wrr3(btens,9,'b   ')
      call mooney(ndime,S1,S2,ak, detf, btens,press,sigma) ! napon 
      call wrr3(sigma,9,'sig ')
      call cdevia (ndime, S1,S2, detf, btens,press,ctens,1)!tenzor elasticnosti
      call cpress (ndime, press, ctens,0)
      call ckapa (ndime,ak, detf, ctens)

      CALL TENVEK(sigma,STRESS,1)
      CALL CEPMT(ELAST,ctens,1)
      call jedna1(strai1,strain,4)
      call jedna1(stres1,stress,4)
      call wrr(strai1,4,'stra')
      call wrr(stres1,4,'stre')
      RETURN
      END
       
	!-----------------------------------------------------------------------
      SUBROUTINE mooney(ndime, ami10,ami01,ak, detf, btens,press, sigma)
!-----------------------------------------------------------------------
!
!    Racuna  Kosijev  tenzor napona za Mooney-Rivlin materijalni model
!
!     ndime  -->  dimenzija 
!     ami10,ami01    -->  materijalne konstante
!     ak      -->  bulk moduo 
!     detf   -->  determinanta F
!     btens  -->  levi Kosi-Grinov tenzor
!     sigma  -->  Kosijev tenyor napona
!
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION sigma (3, 3), btens (3, 3),bb(3,3)
!	DATA ( (delta (i, j), i = 1, 3), j = 1, 3) 
!     &    / 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
!
      call invar(btens,bI,bII)

      ami1 = 2*ami10 / (detf** (5.d0 / 3.d0) )
      ami0 = 2*ami01 / (detf** (7.d0 / 3.d0) )

!     prva dva clana u izrazu za devijatorski deo

      DO id = 1, ndime
         DO jd = 1, ndime
            sigma(id,jd)=ami1*(btens(id,jd)-bI*delta(id,jd)/3.d0)
     &	                +ami0*bI*(btens(id,jd)-bI*delta(id,jd)/3)
         END DO
      END DO

      call bbproiz(btens,bb)

!      treci clan koji se oduzima

      DO id=1,ndime
       DO jd=1,ndime
        sigma(id,jd)=sigma(id,jd)-ami0*(bb(id,jd)-bII*delta(id,jd)/3.d0)
       END DO
      END DO

!     Sferni deo 

      press=ak*(detf-1)

      DO id=1,3
         sigma(id,id)=sigma(id,id)+press
      END DO
      


      END SUBROUTINE mooney
 

!-----------------------------------------------------------------------
      SUBROUTINE invar(btens,bI,bII)
!-----------------------------------------------------------------------
!     Racuna invarijante tenzora 3x3 
!     
!     btens --> tenzor 3x3
!     bI --> prva invarijanta bI=tr b 
!     bII --> druga invarijanta bII
!     bIII --> treca invarijanta bIII=det(b)
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION btens (3, 3)
      
      bI = 0.d0
      DO id = 1, 3
        bI = bI + btens(id, id)
      END DO
!     druga invarijanta u obliku bII=b:b=BklBkl
      DO id=1,3
        DO jd=1,3

          bII=bII+btens(id,jd)*btens(id,jd)

        END DO
      END DO


      END SUBROUTINE invar

      SUBROUTINE bbproiz(btens,bb)
!-----------------------------------------------------------------------
!     Racuna proizvod oblika bb=btens*btens
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION bb (3, 3), btens (3, 3),bo(3,3)
!	DATA ( (bo (i, j), i = 1, 3), j = 1, 3) 
!     &    / 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
!      racuna bb=F*C*FT -->moze se koristiti i subroutine PIOKOS 
!      ali prethodno treba transformisati matricu C u vektor 
      call clear(bb,9)
      DO id=1,3
        DO jd=1,3
          DO kd=1,3
            bb(id,jd)=bb(id,jd)+btens(id,kd)*btens(kd,jd)
          END DO
        END DO
      END DO
      
      END SUBROUTINE bbproiz
!-----------------------------------------------------------------------
       SUBROUTINE izvinvar(aI,IN)
!-----------------------------------------------------------------------
!    Racuna izvod inverznog tenzora tj. - d ct ^(-1)/d ct pri cemu je primenjena 
!    operacija unapred tj. tenzor i 
!    IN je indikator 
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  aI (3, 3 , 3 ,3)
       IF (IN.EQ.0) THEN
      DO id=1,3
        DO jd=1,3
          DO kd=1,3
	      DO ld=1,3
             aI(id,jd,kd,ld)=0.5*(delta(id,kd)*delta(jd,ld)
     &		   +delta(id,ld)*delta(jd,kd))
          END DO
        END DO
      END DO
	END DO
	END IF 
       IF (IN.EQ.1) THEN
      DO id=1,3
        DO jd=1,3
          DO kd=1,3
	      DO ld=1,3
             aI(id,jd,kd,ld)=delta(id,kd)*delta(jd,ld)
          END DO
        END DO
      END DO
	END DO
	END IF 
      
      END SUBROUTINE izvinvar
!-----------------------------------------------------------------------

      SUBROUTINE cdevia (ndime, ami10,ami01,detf, btens,press,ctens,IND)
!-----------------------------------------------------------------------
!
!     Racuna komponentu tenzora elasticnosti 
!
!
!     ndime  --> dimenzija
!     ami10,ami01    -->  materijalne konstante 
!     detf   -->  determinant  F
!     btens  -->  levi Kosi-Grinov tenzor 
!     ctens  -->  tenzor elasticnosti  c(i,j,k,l)
!
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION ctens (3, 3, 3, 3), btens (3, 3),bb(3,3),aI(3,3,3,3)
!	DATA ( (delta(i, j), i = 1, 3), j = 1, 3) 

!    &    / 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
 
!
!     
!
      ami1=4*ami10/(3*detf** (5.d0 / 3.d0) )
      ami0=4*ami01/(detf** (7.d0 / 3.d0) )
      trace = 0.d0
      

      call invar(btens,bI,bII)
      call bbproiz(btens,bb)
	call izvinvar(aI,IND)

!
!     
!    slucaj kada je tenzor identiteta jednak delta(i,k)*delta(j,l)
!      if (IND.EQ.0) THEN

      DO i = 1, ndime
        DO j = 1, ndime
          DO k = 1, ndime
            DO l = 1, ndime
              ct1 = ami1*(bI*aI(i,j,k,l) 
     &            -btens(i,j)*delta(k,l)-delta(i,j)*btens(k,l)
     &            +bI*delta(i,j)*delta(k,l)/3.d0)
               ct2 = ami0*(btens(i,j)*btens(k,l)
     &               -bI*btens(i,j)*delta(k,l)/3.d0         
     &               -2*bI*delta(i,j)*btens(k,l)/3.d0
     &              +2*delta(i,j)*delta(k,l)*(bI*bI)/9.d0        
     &               -bI*btens(i,j)*delta(k,l)/3 
     &                +aI(i,j,k,l)*(bI*bI)/3.d0)
                ct3 = ami0*(-2*delta(i,j)*bb(k,l)/3.d0+
     &                  2*bII*delta(i,j)*delta(k,l)/9.d0
     &                   +btens(i,k)*btens(j,l)          
     &                   -2*bb(i,j)*delta(k,l)/3
     &                    +bII*aI(i,j,k,l)/3.d0)
              ctens(i,j,k,l)=ct1+ct2-ct3                  
            END DO
          END DO
        END DO
      END DO
!     END IF
!    slucaj kada je tenzor identiteta jednak (0.5*delta(i,k)*delta(j,l)+0.5delta(i,l)*delta(j,k))
!      IF (IND.EQ.1) THEN
!      DO i = 1, ndime
!        DO j = 1, ndime
!          DO k = 1, ndime
!           DO l = 1, ndime
!
!              ct1 = ami1*(bI*(0.5*delta(i,k)*delta(j,l)
!     &			  +0.5*delta(i,l)*delta(j,k)) 
!     &              -delta(i,j)*btens(k,l)-btens(i,j)*delta(k,l)
!     &            +bI*delta(i,j)*delta(k,l)/3.d0)
!               ct2 = ami0*(btens(i,j)*btens(k,l)
!     &               -bI*delta(i,j)*btens(k,l)/3.d0         
!     &               -2*bI*btens(i,j)*delta(k,l)/3.d0
!     &              +2*delta(i,j)*delta(k,l)*(bI**2)/9.d0        
!    &               -bI*delta(i,j)*btens(k,l)/3 
!     &                +(0.5*delta(i,k)*delta(j,l)
!     &                +0.5*delta(i,l)*delta(j,k))*(bI**2.d0)/3.d0)
!                ct3 = ami0*(-2*delta(i,j)*bb(k,l)/3.d0+
!     &                  2*bII*delta(i,j)*delta(k,l)/9.d0
!     &                   +btens(i,j)*btens(k,l)          
!     &                   -2*bb(i,j)*delta(k,l)/3
!     &                    +bII*(0.5*delta(i,k)*delta(j,l)
!     &                    +0.5*delta(i,l)*delta(j,k))/3.d0)
!             ctens(i,j,k,l)=ct1+ct2-ct3                  
!            END DO
!          END DO
!        END DO
!      END DO
!     END IF




      END SUBROUTINE cdevia
 
!-----------------------------------------------------------------------
      SUBROUTINE cpress (ndime, press, ctens,IND)
!-----------------------------------------------------------------------
!
!     Racuna komponentu tenzora elasticnosti oblika Cp=p(I*I-2i)
!
!
!     ndime  -->  dimenzija
!     press  -->  pritisak
!     ctens  -->  tenzor elasticnosti c(i,j,k,l)
!
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION ctens (3, 3, 3, 3),aI(3,3,3,3)
!      DATA ( (delta (i, j), i = 1, 3), j = 1, 3) 
!     &    / 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
!
      
	call izvinvar(aI,IND)
!    slucaj kada je tenzor identiteta jednak delta(i,k)*delta(j,l)
!      if (IND.EQ.0) THEN
      DO i = 1, ndime
        DO j = 1, ndime
          DO k = 1, ndime
            DO l = 1, ndime

               ctens (i, j, k, l) = ctens (i, j, k, l) 
     &                            + press * delta (i, j) * delta (k, l) 
     &                            -2* press * aI(i,j,k,l) 
            END DO
          END DO
        END DO
      END DO
!     END IF
!    slucaj kada je tenzor identiteta jednak (0.5*delta(i,k)*delta(j,l)+0.5delta(i,l)*delta(j,k))
!     IF (IND.EQ.1) THEN
!      DO i = 1, ndime
!        DO j = 1, ndime
!          DO k = 1, ndime
!            DO l = 1, ndime
!               ctens (i, j, k, l) = ctens (i, j, k, l) 
!     &                            + press * delta (i, j) * delta (k, l) 
!     &                            - press * delta (i, k) * delta (j, l) 
!     &                            - press * delta (i, l) * delta (j, k)
!     
!            END DO
!          END DO
!        END DO
!
!      END DO
!      END IF

      END SUBROUTINE cpress
 
C  =====================================================================


!-----------------------------------------------------------------------
      SUBROUTINE ckapa (ndime,ak, detf, ctens)
!-----------------------------------------------------------------------
!
!     Racuna komponentu tenzora elasticnosti oblika Ck=kJ I*I
!
!
!     ndime  -->  dimenzija
!     k  -->  bulk moduo
!     detf--> determinanta tenzora gradijenta deformacije
!     ctens  -->  tenzor elasticnosti c(i,j,k,l)
!
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION ctens (3, 3, 3, 3)
!      DATA ( (delta (i, j), i = 1, 3), j = 1, 3) 
!     &    / 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /
!

!
      DO i = 1, ndime
        DO j = 1, ndime
          DO k = 1, ndime
            DO l = 1, ndime
              ctens (i, j, k, l) = ctens (i, j, k, l) 
     &                           + ak * detf*delta (i, j) * delta (k, l) 

            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE ckapa


C======================================================================
      SUBROUTINE KOSI3(PLAST,CORGT,AU,ISNA,IPRC,MCVEL,ICVEL,LL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      
C
C     PODPROGRAM ZA TRANSFORMACIJU PIOLA-KIRCHHOFF NAPONA
C     U CAUCHY NAPON (3D)
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /SRPSKI/ ISRPS
      COMMON /GRADPOLOS/ E11,E22,E33,E12,E13,E23,E21,E31,E32
      DIMENSION PLAST(*),ISNA(*),IPRC(*),MCVEL(*)
      DIMENSION CORGT(3,NGS12,*),SOK(9),BOK(9),R1(9),R2(9),R3(9)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAN39'
         SOK(1)=PLAST(LL)
         SOK(2)=PLAST(LL+3)
         SOK(3)=PLAST(LL+5)
         SOK(4)=PLAST(LL+3)
         SOK(5)=PLAST(LL+1)
         SOK(6)=PLAST(LL+4)
         SOK(7)=PLAST(LL+5)
         SOK(8)=PLAST(LL+4)
         SOK(9)=PLAST(LL+2)
! VELICINE E NISU DEFINISANE
!STA RADI OVA SUBROUTINA
         BOK(1)=E11+1
         BOK(2)=E21
         BOK(3)=E31
         BOK(4)=E12
         BOK(5)=E22+1
         BOK(6)=E32
         BOK(7)=E13
         BOK(8)=E23
         BOK(9)=E33+1
         CALL MTRA(BOK,R1,3,3)
         CALL GMPRD (SOK,R1,R2,3,3,3)
         CALL GMPRD (BOK,R2,R3,3,3,3)
         PLAST(LL)=R3(1)
         PLAST(LL+1)=R3(5)
         PLAST(LL+2)=R3(9)
         PLAST(LL+3)=R3(2)
         PLAST(LL+4)=R3(3)
         PLAST(LL+5)=R3(7)
                 RETURN
                 END

