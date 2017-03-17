C     3D OR PLANE STRAIN NEARLY INCOMPRESIBLE NEO-HOOKEAN
C
C
      SUBROUTINE D3M28 (TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
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
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D2M28'
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
C
      CALL TI328(PLAST(LTAU),PLAST(LDEFT),
     1           PLAS1(LTAU1),PLAS1(LDEFT1),
     1           A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI328(STREST,STRAIT,
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
      COMMON /BOJANA/ CTENS(3,3,3,3),P2KStress(3,3)
      COMMON /CDEBUG/ IDEBUG
      COMMON /PRINCI/ PRINC(3)
C
      DIMENSION sigma(3,3),STRESS(*),STRAIN(*),btens(3,3)!definisi sve matrice i vektore
      DIMENSION tenzC(3,3),Cinv(3,3),vekC(3),alamdw(3),alamd(3)
      DIMENSION disowOverdlambda(3),d2isowOdlambda1dlambda2(3,3),
     1          d2isowOverdlamd2(3)
      DIMENSION wa(3),Yval(3,3)
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
c  	s1=184.84375;
c 	s2=26.40625;
c 	s3=0.0;
	ndime=3;
	write(3,*)'s1,s2,s3',s1,s2,s3
	call wrr6(strain,6,'stri')
c     Racunamo desni Kosi-Grinov tenzor
	call desniKG(strain,tenzC)
	call wrr3(tenzc,9,'tenc')
c     Racunamo njegov inverzni tenzor i odgovarajucu determinantu
      CALL JEDNA1(Cinv,tenzC,9)! prepisujemo u Cinv tenzC da ga ne bi izgubili
      call MINV3(Cinv,detC)
	call wrr3(cinv,9,'cinv')
c      Racunamo vrednost J jna (4.61) B.Jeremic
      aJ=dsqrt(detC);
	write(3,*)'aJ',aJ
	CALL TENVEK(tenzC,vekC,1)
	call wrr6(vekc,6,'vekc')
c     Trazimo sopstvene vrednosti tenzora C
c	call GLAVNBOJ(vekC)
      call GLAVN(vekC)
	call wrr3(princ,3,'c123')
	
      call lamda_wave(princ,aJ,alamdw,alamd,index)
c	call wrr3(alamd,3,'lamd')
c	call wrr3(alamdw,3,'lamw')
c	write(3,*) 'index',index

c     KORELC=0
ckorelc
c      IF(KORELC.EQ.1) THEN
c         call Testmr(S1,S2,alamd,tenzC,Ctens,Sigma)
c        prebacivanje tenzora u vektor
c         CALL TENVEK(Sigma,Stress,1)
c         CALL CEPMT(ELAST,ctens,1)
c         call jedna1(strai1,strain,6)
c         call jedna1(stres1,stress,6)
c         return
c      ENDIF
ckorelc

c     odredjivanje energije W
      call MooneyRivlinWenergy(s1,s2,s3,aJ,alamdw,w_iso,w_vol)
c	write(3,*) 'w_iso,w_vol',w_iso,w_vol

c     pozivi drugih podprograma za racunanje odgovarajucih izvoda energije
      call disoWpodLamda(s1,s2,alamdw,disowOverdlambda)
c	call wrr3(disowOverdlambda,3,'dis1')
      call dvolWpodJ(s1,s2,s3,aJ,dvolwOdJ)
c	write(3,*) 'dvolwOdJ',dvolwOdJ
      call d2volWpod2J(s1,s2,s3,aJ,d2volwOdJ2)
c	write(3,*) 'd2volwOdJ2',d2volwOdJ2
      call diso2Wpod2Lamda(s1,s2,alamdw,d2isowOverdlamd2)
c	call wrr3(d2isowOverdlamd2,3,'dis2')
      call diso2WpodLamda1dLamda2(s1,s2,d2isowOdlambda1dlambda2)
c	call wrr3(d2isowOdlambda1dlambda2,9,'dis0')
	call waDEF(alamdw,disowOverdlambda,wa)
c	call wrr3(wa,3,'wa  ')
      call yab(alamdw,d2isowOdlambda1dlambda2,
     &                disowOverdlambda,d2isowOverdlamd2,Yval)
c	call wrr3(Yval,9,'Yval')
c     racunanje Piola-Kirhofovog tenzora napona II reda 
c	call wrr3(P2KStress,9,'P2K0')
      call cP2KStress(tenzC,aJ,dvolwOdJ,alamd,Cinv,index,P2KStress,wa)
	call wrr3(P2KStress,9,'P2K3')
c     prebacivanje tenzora u vektor
      CALL TENVEK(P2KStress,Stress,1)
c	call wrr6(Stress,6,'STRE')
c      pretvaranje P2K napona u Kosijev 
c      call PIOKOS(GRAD,Stress)
c	CALL TENVEK(Stress,sigma,1)
c      call wrr3(sigma,9,'sig ') 
c     racunanje tenzora elasticnosti 
      call Stiffness(tenzC,Cinv,ndime,alamd,alamdw,index,aJ,ctens,wa
     &,d2isowOverdlamd2,disowOverdlambda,dvolwOdJ,d2volwOdJ2,Yval)
c      call Stiffness1(ctens,grad,ctens)      
c      CALL TENVEK(sigma,STRESS,1)
C OVO IZBACITI?????
      CALL CEPMT(ELAST,ctens,1)
c	call wrr6(ELAST,36,'elas')
      call jedna1(strai1,strain,6)
      call jedna1(stres1,stress,6)
c      call wrr6(strai1,6,'stra')
c      call wrr6(stres1,6,'stre')
      RETURN
      END
      
!     call jedna1(btens,grap,9)
!     Material 5 nearly incompressible neo-Hookean
!
!      ndime=2
!      detf=detg
	
!     call wrr3(btens,9,'b   ')
!      call mooney(ndime,S1,S2,S3, detf, btens,press,sigma) ! napon 
!      call wrr3(sigma,9,'sig ')
!     call cdevia (ndime, S1,S2, detf, btens,press,ctens,0)!tenzor elasticnosti
!      call cpress (ndime, press, ctens,0)
!      call ckapa (ndime,S3, detf, ctens)
!-----------------------------------------------------------------------
      SUBROUTINE TENVEK5(T,V)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     TRANSFORMISANJE : VEKTOR -> TENZOR  (IND=0)
CS                       TENZOR  -> VEKTOR (IND=1)

      DIMENSION T(3,*),V(*)
C
C     T(I,J) = V(K) 
C
 
         T(1,1)=V(1)
         T(2,2)=V(2)
         T(3,3)=V(3)
         T(1,2)=0.5*V(4)
         T(2,3)=0.5*V(5)
         T(1,3)=0.5*V(6)
         T(2,1)=T(1,2)
         T(3,2)=T(2,3)
         T(3,1)=T(1,3)
 
C

      END SUBROUTINE TENVEK5

!-----------------------------------------------------------------------
      SUBROUTINE desniKG(strain,tenzC)
!-----------------------------------------------------------------------
!
!    Racuna desni Kosi-Grinov tenzor C
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  tenzC(3,3),tenStr(3,3)
         
      call TENVEK5(tenStr,strain)
	DO id=1,3
	Do jd=1,3
	tenzC(id,jd)=2.*tenStr(id,jd)+delta(id,jd);
	END DO
	END DO
     
            
      END SUBROUTINE desniKG
 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE desniKG1(grad,tenzC)
!-----------------------------------------------------------------------
!
!    Racuna desni Kosi-Grinov tenzor C=FT*F
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  tenzC(3,3),grad(3,3)
         
      
	DO id=1,3
	Do jd=1,3
	   tenzC(id,jd)=0
	DO k=1,3
	   tenzC(id,jd)=tenzC(id,jd)+grad(k,id)*grad(k,jd)
	end do
	END DO
	END DO
     
            
      END SUBROUTINE desniKG1
 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      SUBROUTINE GLAVNBOJ(S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
!-----------------------------------------------------------------------

C .
CS.   PROGRAM
CS.      ZA RACUNANJE GLAVNIH NAPONA ZA 3/D ELEMENT
CE.   PROGRAM
CE.      TO CALCULATE PRINCIPAL STRESSES FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /PRINCI/ PRINC(3)
      COMMON /DEFNAP/ NAPDEF
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*)
C OVU TOLERANCIJU PROVERITI 1.D-08 -STARA TOLERANCIJA
C      DATA TOL/1.D-08/
      DATA TOL/1.D-16/
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAVN'
      IZILE=0
C
C   INARIJANTE TENZORA
C
C      ONE=0.9999999999999999D0
      ONE=1.D0
      CC=(S(1)+S(2)+S(3))/3.D0
      DO 20 I=1,3                                                              
   20 PRINC(I)=S(I)-CC                                                        
      C2=(PRINC(1)*PRINC(1)+PRINC(2)*PRINC(2)+PRINC(3)*PRINC(3))*.5D0+
     1    S(4)*S(4)+S(5)*S(5)+S(6)*S(6)             
c      C2=-PRINC(1)*PRINC(2)-PRINC(2)*PRINC(3)-PRINC(3)*PRINC(1) 
c     1   +S(4)*S(4)+S(5)*S(5)+S(6)*S(6) 
      C3=PRINC(1)*(PRINC(2)*PRINC(3)-S(5)*S(5))+
     1   S(4)*(S(5)*S(6)-S(4)*PRINC(3))+        
     1   S(6)*(S(4)*S(5)-PRINC(2)*S(6))                                          
      DUM=S(1)*(S(2)*S(3)-S(5)*S(5))+
     1    S(4)*(S(5)*S(6)-S(4)*S(3))+        
     1    S(6)*(S(4)*S(5)-S(2)*S(6))                                          
C
      SQ2 =DSQRT(2.D0)
      PI23=2.D0*DACOS(0.5D0)
C      DUM=DUM*1.D-8
C      IF(DABS(DUM).LT.1.D-9) DUM=1.D-9
      DUM=DUM*TOL
      IF(DABS(DUM).LT.0.1D0*TOL) DUM=0.1D0*TOL
      DUMM=DUM
C      WRITE(3,*) 'DUMM,C2,CC',DUMM,C2,CC
      IF(C2.LT.DABS(DUM))THEN
         DO 330 I=1,3
            PRINC(I)=CC
C ZILE 
C           ZBOG DSQRT(PRINC(I))
            IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  330    CONTINUE
C ZILE
          IF(IZILE.EQ.1) THEN
             PRINC(1)=PRINC(1)+DUMM
             PRINC(2)=PRINC(2)-DUMM
          ENDIF
C ZILE
         RETURN
      ENDIF
      T = DSQRT(C2/1.5D0)                                                           
      A = C3*SQ2/(T*T*T)
C
C     IF ( A .LT. -1.D0 ) A=-ONE
C     IF ( A .GT.  1.D0 ) A= ONE
      IF(DABS(A).GT.ONE) A=DSIGN(ONE,A)
      A=DACOS(A)/3.D0                                                             
      T=T*SQ2
C
      PRINC(1)=T*DCOS(A)                                                       
      PRINC(2)=T*DCOS(A-PI23)                                                
      PRINC(3)=T*DCOS(A+PI23)                                                
C
      DO 240 I=1,3                                                              
         PRINC(I)=PRINC(I)+CC                                                    
C ZILE 
C        ZBOG DSQRT(PRINC(I))
         IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  240 CONTINUE
      END SUBROUTINE GLAVNBOJ
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      SUBROUTINE d_macheps(eps)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
            
         eps = 1.0;
         icounter = 1;
	DO 10 WHILE ( (1.0 + eps) .GT. (1.0) )
        eps = eps/2.0;
        icounter=icounter+1;
   10  END DO
        eps = eps*2.0;
       icounter=icounter-1;


       END SUBROUTINE d_macheps
!-----------------------------------------------------------------------


      SUBROUTINE lamda_wave(eigen,aJ,alamdw,alamd,index)
!-----------------------------------------------------------------------
!      Jednacina (4.83) iz B.Jeremic
!      Racuna  lamde i odredjuje index u zavisnosti od odnosa lamdi
!      eigen  -> sopstvene vrednosti tenzora C
!      alamd  - > predstavlja vektor (lamda1,lamda2,lamda3)
!      alamdw-> predstavlja vektor(lamda_wave1,lamda_wave2,lamda_wave3)
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  alamd(3),alamdw(3),eigen(3)
         
      DO id=1,3
          alamd(id)=dsqrt(eigen(id));
      END DO
       
      aJJJ = 1./(aJ** (1./3));
      
      DO id=1,3
       alamdw(id)=alamd(id)*aJJJ;
      END DO
      
      diff12 = abs(alamd(1)-alamd(2));
      diff23 = abs(alamd(2)-alamd(3));
      call d_macheps(eps)
      perturbation = eps**0.4;
c	perturbation = 0.01;
c	write(3,*)'tol,l1,l2,l3',perturbation ,alamd(1),alamd(2),alamd(3)
       if (( diff12.GE. perturbation) .AND. (diff23 .GE. perturbation))
     &  index = 0;
	
       if ((diff12 .GE. perturbation ).AND.(diff23 .LT. perturbation))
     &  index = 11;
	
       if ((diff12 .LT. perturbation).AND.(diff23 .GE. perturbation)) 
     &  index = 13;
	
       if ((diff12 .LT. perturbation).AND.(diff23 .LT. perturbation)) 
     &  index = 2;
c	write(3,*)'index',index,diff12,diff23
         
      END SUBROUTINE lamda_wave
 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE MooneyRivlinWenergy(c1,c2,c3,aJ,alamdw,w_iso,w_vol)
!-----------------------------------------------------------------------
!
!     Energijska funkcija data jednacinom (4.124) sa str. 113  B.Jeremica 
!     Dati model predstavlja model koji se sastoji iz izohorskog dela datog jednacinom (4.124)
!     iz B.Jeremica i zapreminskog dela koji je dat Simo-Pister-ovim modelom (4.137)
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)

      DIMENSION alamdw(3)
         
      temp1 = alamdw(1) * alamdw(1)
     &                      + alamdw(2) * alamdw(2)
     &                      + alamdw(3) * alamdw(3) - 3.0;
      temp2 = 1.0 / (alamdw(1)*alamdw(1))
     &                      + 1.0 / (alamdw(2)*alamdw(2))
     &                      + 1.0 / (alamdw(3)*alamdw(3)) - 3.0;
      w_iso = c1 * temp1 + c2 * temp2;
      
      ! zapreminski deo energije -Simo-Pister model(4.118)
c      w_vol=c3*(aJ*aJ-1.-2*log(aJ))/4.
       w_vol=0;
      
      
      END SUBROUTINE MooneyRivlinWenergy
 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE disoWpodLamda(c1,c2,alamdw,disowOverdlambda)
!-----------------------------------------------------------------------
!
!    Jednacina (4.125) B.Jeremic
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  disowOverdlambda(3),alamdw(3)
         
 
        disowOverdlambda(1) = 2.0 * c1 *alamdw(1) 
     &                        - 2.0 * c2 /(alamdw(1)**3.0);
        disowOverdlambda(2) = 2.0 * c1 *alamdw(2)
     &                        - 2.0 * c2 /(alamdw(2)**3.0);
        disowOverdlambda(3) = 2.0 * c1 *alamdw(3) 
     &                        - 2.0 * c2 /(alamdw(3)**3.0);
            
      END SUBROUTINE disoWpodLamda
 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE dvolWpodJ(c1,c2,c3,aJ,dvolwOdJ)
!-----------------------------------------------------------------------
!
!   
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      
         
        dvolwOdJ =c3*(aJ-1./aJ)/2. ;
            
      END SUBROUTINE dvolWpodJ
 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE d2volWpod2J(c1,c2,c3,aJ,d2volwOdJ2)
!-----------------------------------------------------------------------
!
!   
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      
         
        d2volwOdJ2 =c3*(1+1/(aJ*aJ))/2. ;
            
      END SUBROUTINE d2volWpod2J
 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE diso2Wpod2Lamda(c1,c2,alamdw,d2isowOverdlambda2)
!-----------------------------------------------------------------------
!
!     Jednacina (4.126) B.Jeremic 
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  d2isowOverdlambda2(3),alamdw(3)
         
 
       d2isowOverdlambda2(1) = 2.0 * c1  + 6.0 * c2/(alamdw(1)**4.0);
       d2isowOverdlambda2(2) = 2.0 * c1  + 6.0 * c2/(alamdw(2)**4.0);
       d2isowOverdlambda2(3) = 2.0 * c1  + 6.0 * c2/(alamdw(3)**4.0);
            
      END SUBROUTINE diso2Wpod2Lamda
 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE diso2WpodLamda1dLamda2(c1,c2
     &                                  ,d2isowOdlambda1dlambda2)
!-----------------------------------------------------------------------
!
!     Jednacina (4.127) B.Jeremic 
!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  d2isowOdlambda1dlambda2(3,3)
         
      
      DO i = 1,3
        DO j = 1,3
          
       d2isowOdlambda1dlambda2(i,j) = 0.0;!proveriti

       END DO 
	END DO
            
      END SUBROUTINE diso2WpodLamda1dLamda2
 

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      SUBROUTINE waDEF(alamdw,disowOverlambda,wa)
!-----------------------------------------------------------------------
!     Jednacina (D.25) iz B.Jeremic
!-----------------------------------------------------------------------
!  
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  wa(3),alamdw(3),disowOverlambda(3)

       temp = disowOverlambda(1) * alamdw(1) +
     &        disowOverlambda(2) * alamdw(2) +
     &        disowOverlambda(3) * alamdw(3) ;
       temp = temp * (-0.3333333333333333333333333333);
       wa(1) = temp + disowOverlambda(1) * alamdw(1);
       wa(2) = temp + disowOverlambda(2) * alamdw(2);
       wa(3) = temp + disowOverlambda(3) * alamdw(3);


       end SUBROUTINE waDEF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE yab(alamdw,d2isowOdlambda1dlambda2,
     &                disowOdlambda,d2isowOdlambda2,Yval)
!-----------------------------------------------------------------------

!  
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  d2(3,3),d1(3),d11(3),tempi(3),Yval(3,3)
	DIMENSION  alamdw(3)
	DIMENSION d2isowOdlambda1dlambda2(3,3),disowOdlambda(3)
	DIMENSION d2isowOdlambda2(3)
      ! popraviti!
       
         call jedna1(d2,d2isowOdlambda1dlambda2,9)
         call jedna1(d1,disowOdlambda,3)
         call jedna1(d11,d2isowOdlambda2,3)

	  d2(1,1) = d11(1);
        d2(2,2) = d11(2);
        d2(3,3) = d11(3);
c	call wrr3(d2,9,'d2  ')
        tempd = d1(1)*alamdw(1) + d1(2)*alamdw(2) + d1(3)*alamdw(3) ;
        tempcd = 0.0;
         DO id=1,3
          tempi(id) = 0.0;
           DO jd=1,3

           tempi(id)=tempi(id)+ d2(id,jd) * alamdw(id) * alamdw(jd);
           tempcd= tempcd+ d2(id,jd) * alamdw(id) * alamdw(jd);
        
          END DO
          END DO

!        Jednacina (D.33) B.Jeremic   
          
        DO id=1,3
        DO jd=1,3
                
              Yval(id,jd) = d1(id)*delta(id,jd)*alamdw(jd) 
     &         + d2(id,jd)*alamdw(id)*alamdw(jd) -
     &         (  tempi(id) + tempi(jd) + d1(id)*alamdw(id) 
     &         + d1(jd)*alamdw(jd) ) / 3.0 +
     &                      ( tempcd + tempd ) / 9.0;

	END DO
	END DO

      


      end SUBROUTINE yab
!-----------------------------------------------------------------------
      SUBROUTINE cP2KStress(C,aJ,dvolwOdJ,alamd,Cinv,index,PK2Stress,wa)
!-----------------------------------------------------------------------
!
!      Racuna  drugi Piola-Kirhfofov tenzor napona 
!      alamdw-> predstavlja vektor(lamda_wave1,lamda_wave2,lamda_wave3)

!-----------------------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (a - h, o - z)
      DIMENSION  aM1(3,3),aM2(3,3),aM3(3,3),PK2Stress(3,3)
      DIMENSION  aMr(3,3),volPK2Stress(3,3)
	DIMENSION  alamd(3),wa(3),Cinv(3,3),C(3,3)

      aI1 = alamd(1)*alamd(1)+alamd(2)*alamd(2)+alamd(3)*alamd(3);   ! j-na (5.71) FlagShyp 
       
      if (index == 0)then
      
      !     jednacina (5.77) B.Jeremic
       d1 = (alamd(1)+alamd(2))*(alamd(1)+alamd(3))*(alamd(1)
     &                                 -alamd(2))*(alamd(1)-alamd(3));
       d2 = (alamd(2)+alamd(3))*(alamd(2)+alamd(1))*(alamd(2)
     &                                 -alamd(3))*(alamd(2)-alamd(1));
       d3 = (alamd(3)+alamd(1))*(alamd(3)+alamd(2))*(alamd(3)
     &                                 -alamd(1))*(alamd(3)-alamd(2));

      alamd1=alamd(1);
      alamd2=alamd(2);
      alamd3=alamd(3);
      
      ! jednacina (5.86) iz B. Jeremic
      
      DO id=1,3
        DO jd=1,3

        aM1(id,jd) = ( C(id,jd) - delta(id,jd)*(aI1-alamd1*alamd1) 
     &                + Cinv(id,jd)*(aJ*aJ/(alamd1*alamd1)) )/d1;
        aM2(id,jd) = ( C(id,jd) - delta(id,jd)*(aI1-alamd2*alamd2)
     &                + Cinv(id,jd)*(aJ*aJ/(alamd2*alamd2)) )/d2;
        aM3(id,jd) = ( C(id,jd) - delta(id,jd)*(aI1-alamd3*alamd3) 
     &                + Cinv(id,jd)*(aJ*aJ/(alamd3*alamd3)) )/d3;

	  PK2Stress (id,jd)= aM1(id,jd)*Wa(1) + aM2(id,jd)*Wa(2) !drugi clan od (4.99)
     &                     + aM3(id,jd)*Wa(3);
        END DO
      END DO
     
c      call wrr3(aM1,9,'aM1')
c      call wrr3(aM2,9,'aM2')
c	call wrr3(aM3,9,'aM3')
c	call wrr3(PK2Stress,9,'PK2Stress')
      
      end if

       if (index == 11) then 
! isto kao i slucaj kada su lamda1=lamda2 a razlicito od lamda3
       DO id=1,3
        DO jd=1,3

          aM1(id,jd) = (delta(id,jd) - Cinv(id,jd) 
     &       * (alamd2*alamd2)) * (1.0/(alamd1+alamd2)/(alamd1-alamd2));
          aMr(id,jd) = Cinv(id,jd) - aM1(id,jd);

          PK2Stress(id,jd) = aMr(id,jd)*Wa(3) + aM1(id,jd)*Wa(1);

         END DO
        END DO
     
      end if

       if (index == 13) then
	!j-na (12) rad B.Jeremica
	!kada je lamda1=lamda2razlicito od lamda3

        DO id=1,3
         DO jd=1,3
          aM3(id,jd) = (delta(id,jd) - Cinv(id,jd) * (alamd2*alamd2)) 
     &                * (1.0/(alamd3+alamd2)/(alamd3-alamd2));

          aMr(id,jd) = Cinv(id,jd) - aM3(id,jd);

           PK2Stress(id,jd) = aMr(id,jd)*Wa(1) + aM3(id,jd)*Wa(3);

          END DO
         END DO

       
       end if
 
       if (index == 2)then
	! slucaj kada je lamda1=lamda2=lamda3,j-na (12) rad B.Jeremica
c      DO id=1,3
c        DO jd=1,3
c          PK2Stress(id,jd) = 0;
c         end do
c        end do
        end if   
   
c  zapreminski deo 
      
      
             dWdJ = dvolwOdJ;!izvod zapreminskog dela po J
        DO id=1,3
          DO jd=1,3
             volPK2Stress(id,jd) = Cinv(id,jd) * aJ * dWdJ;
c ukupna vrednost napona 
           end do
        end do
       DO id=1,3
        DO jd=1,3
          PK2Stress(id,jd) = volPK2Stress(id,jd) + PK2Stress(id,jd);
         end do
        end do

c	call wrr3(PK2Stress,9,'PK2Stress')
  
      END SUBROUTINE cP2KStress
 



!-----------------------------------------------------------------------
      SUBROUTINE Stiffness(C,Cinv,ndime,alamd,alamdw,index,aJ,ctens,wa
     &	     ,d2isowOverdlamd2,disowOdlambda,dvolwOdJ,d2volwOdJ2,Yval)
!-----------------------------------------------------------------------
!     Racuna tenzor elasticnosti 
!-----------------------------------------------------------------------

       IMPLICIT DOUBLE PRECISION (a - h, o - z)
       DIMENSION  aI(3,3,3,3),CinvCinv(3,3,3,3),aI4s(3,3,3,3)
	 DIMENSION  aICinv(3,3,3,3),tempI(3,3,3,3),CinvCinv_ICinv(3,3,3,3)
	 DIMENSION aM1(3,3),aM2(3,3),aM3(3,3),aMr(3,3),wa(3)
       DIMENSION Cm1M1M1Cm1(3,3,3,3),Cm1M2M2Cm1(3,3,3,3)
       DIMENSION Cm1M3M3Cm1(3,3,3,3),dM1M1d(3,3,3,3),dM2M2d(3,3,3,3)
	DIMENSION  aM1M1(3,3,3,3),aM2M2(3,3,3,3),aM3M3(3,3,3,3)
	DIMENSION calM1(3,3,3,3),calM2(3,3,3,3),calM3(3,3,3,3)
	DIMENSION aL_iso(3,3,3,3),dM3M3d(3,3,3,3),calMr(3,3,3,3)
	DIMENSION aL_iso_1(3,3,3,3),aL_iso_2(3,3,3,3),ctens(3,3,3,3)
	DIMENSION aL_vol(3,3,3,3),d11(3),C(3,3),dd1(3),Cinv(3,3)
	DIMENSION alamd(3),d2isowOverdlamd2(3),disowOdlambda(3),Yval(3,3)
	DIMENSION alamdw(3),cll(3,3,3,3)
! Racuna deo aL_iso

   

      DO i = 1, ndime
        DO j = 1, ndime
          DO k = 1, ndime
            DO l = 1, ndime
              aI(i,j,k,l)=delta(i,j)*delta(k,l);
	         CinvCinv(i,j,k,l)=Cinv(i,j)*Cinv(k,l);
            END DO
          END DO
         END DO
        END DO
        
      DO i = 1, ndime
        DO j = 1, ndime
          DO k = 1, ndime
            DO l = 1, ndime
              aI4s(i,j,k,l)=(aI(i,k,j,l)+aI(i,l,j,k))/2.;! tenzor identiteta
c             cll(i,j,k,l)=(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))/2.;
c	        clk=cll(i,j,k,l)-aI4s(i,j,k,l);
	        aICinv(i,j,k,l)=(CinvCinv(i,k,j,l)+CinvCinv(i,l,j,k))/2;
	        tempI(i,j,k,l) = aI4s(i,j,k,l) - aI(i,j,k,l);
	        CinvCinv_ICinv(i,j,k,l) = CinvCinv(i,j,k,l)
     &                                   - aICinv(i,j,k,l); 
            END DO
          END DO
         END DO
        END DO
    
       

       aI1 = alamd(1)*alamd(1) + alamd(2)*alamd(2) + alamd(3)*alamd(3);


        if(index == 0)then
      !sve lamde su razlicite
    !     jednacina (4.77) B.Jeremic
        d1 = (alamd(1)+alamd(2))*(alamd(1)+alamd(3))*(alamd(1)
     &        -alamd(2))*(alamd(1)-alamd(3));
        d2 = (alamd(2)+alamd(3))*(alamd(2)
     &        +alamd(1))*(alamd(2)-alamd(3))*(alamd(2)-alamd(1));
        d3 = (alamd(3)+alamd(1))*(alamd(3)
     &        +alamd(2))*(alamd(3)-alamd(1))*(alamd(3)-alamd(2));
	  
     
      alamd1=alamd(1);
      alamd2=alamd(2);
      alamd3=alamd(3);
      
      ! jednacina (5.86) iz B. Jeremic
      
      DO id=1,3
        DO jd=1,3
   
        aM1(id,jd)=(C(id,jd)- delta(id,jd)*(aI1-alamd1*alamd1)
     &               + Cinv(id,jd)*aJ*aJ/(alamd1*alamd1))/d1;
        aM2(id,jd)=(C(id,jd)-delta(id,jd)*(aI1-alamd2*alamd2) 
     &               + Cinv(id,jd)*(aJ*aJ/(alamd2*alamd2)))/d2;
        aM3(id,jd)=(C(id,jd)-delta(id,jd)*(aI1-alamd3*alamd3) 
     &               +Cinv(id,jd)*(aJ*aJ/(alamd3*alamd3)))/d3;
        END DO
      END DO
    
    ! jednacina (D.15) B. Jeremic 
    
       d1p = 4.0 *alamd1*alamd1*alamd1*alamd1 
     &       - aI1*alamd1*alamd1 - aJ*aJ /(alamd1*alamd1);
       d2p = 4.0 *alamd2*alamd2*alamd2*alamd2 
     &       - aI1*alamd2*alamd2 - aJ*aJ /(alamd2*alamd2);
       d3p = 4.0 *alamd3*alamd3*alamd3*alamd3 
     &       - aI1*alamd3*alamd3 - aJ*aJ /(alamd3*alamd3);
   
   ! jednacina (D.17) iz B.Jeremic

      DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
         Cm1M1M1Cm1(i,j,k,l)=Cinv(i,j)*aM1(k,l) + aM1(i,j)*Cinv(k,l);
         Cm1M2M2Cm1(i,j,k,l)= Cinv(i,j)*aM2(k,l) + aM2(i,j)*Cinv(k,l);
         Cm1M3M3Cm1 (i,j,k,l)= Cinv(i,j)*aM3(k,l) + aM3(i,j)*Cinv(k,l);
         dM1M1d(i,j,k,l) = delta(i,j)*aM1(k,l) + aM1(i,j)*delta(k,l);
         dM2M2d(i,j,k,l) = delta(i,j)*aM2(k,l) + aM2(i,j)*delta(k,l);
         dM3M3d(i,j,k,l) = delta(i,j)*aM3(k,l) + aM3(i,j)*delta(k,l);
         aM1M1(i,j,k,l) = aM1(i,j) * aM1(k,l);
         aM2M2(i,j,k,l) = aM2(i,j) * aM2(k,l);
         aM3M3(i,j,k,l) = aM3(i,j) * aM3(k,l);
             END DO
            END DO
           END DO
       end do
    ! JEDNACINA (D.17) B.Jeremic
      DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
      calM1(i,j,k,l) = (tempI(i,j,k,l) + (CinvCinv_ICinv(i,j,k,l)
     & -Cm1M1M1Cm1(i,j,k,l))*(aJ*aJ/(alamd1*alamd1)) + 
     & dM1M1d(i,j,k,l)*(alamd1*alamd1) - aM1M1(i,j,k,l)*d1p)/d1;
      calM2(i,j,k,l) = ( tempI(i,j,k,l) + (CinvCinv_ICinv(i,j,k,l)
     &  -Cm1M2M2Cm1(i,j,k,l))*(aJ*aJ/(alamd2*alamd2)) +
     & dM2M2d(i,j,k,l)*(alamd2*alamd2) - aM2M2(i,j,k,l)*d2p )/d2;
      calM3(i,j,k,l) = ( tempI(i,j,k,l) + (CinvCinv_ICinv(i,j,k,l)
     &  -Cm1M3M3Cm1(i,j,k,l))*(aJ*aJ/(alamd3*alamd3)) + 
     & dM3M3d(i,j,k,l)*(alamd3*alamd3) - aM3M3(i,j,k,l)*d3p )/d3;
           END DO
          END DO
        END DO
	END DO

! JEDNACINA (4.105) iz B.Jeremic za aL_iso
      DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
        aL_iso_1(i,j,k,l) = (calM1(i,j,k,l)*Wa(1) + calM2(i,j,k,l)*Wa(2)
     &                       + calM3(i,j,k,l)*Wa(3)) * 2.0;
         aL_iso_2 (i,j,k,l)=  aM1(i,j) * aM1(k,l) * Yval(1,1)  
     &                     + aM1(i,j) * aM2(k,l) * Yval(1,2)  + aM1(i,j) 
     &                    * aM3(k,l) * Yval(1,3)  + aM2(i,j) * aM1(k,l)
     &                    * Yval(2,1)  + aM2(i,j) * aM2(k,l) * Yval(2,2)  
     &                      + aM2(i,j) * aM3(k,l) * Yval(2,3)
     &                 +aM3(i,j) * aM1(k,l) * Yval(3,1)  + aM3(i,j) 
     &      * aM2(k,l) * Yval(3,2)  + aM3(i,j) *aM3(k,l) * Yval(3,3);
         aL_iso(i,j,k,l)= aL_iso_1(i,j,k,l) + aL_iso_2(i,j,k,l) ;
           END DO
         END DO
       END DO
      END DO
  
      END IF
  
      if(index == 11) then
 
       d1 = (alamd1+alamd2)*(alamd1+alamd3)
     &      *(alamd1-alamd2)*(alamd1-alamd3);
      DO i=1,3
        DO j=1,3
       
            aM1(i,j) = (delta(i,j) - Cinv(i,j) * (alamd2*alamd2))
     &                 * (1.0/(alamd1+alamd2)/(alamd1-alamd2));! PROVERITI OVE IZRAZE!
            aMr(i,j) = Cinv(i,j) - aM1(i,j);
         end do
      end do
      
           
           d1p = 4.0 *alamd1*alamd1*alamd1*alamd1 
     &            - aI1*alamd1*alamd1 - aJ*aJ /(alamd1*alamd1);
       DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3    
          Cm1M1M1Cm1(i,j,k,l) = Cinv(i,j)*aM1(k,l) + aM1(i,j)*Cinv(k,l);
            dM1M1d(i,j,k,l) = delta(i,j)*aM1(k,l) + aM1(i,j)*delta(k,l);
            aM1M1(i,j,k,l) = aM1(i,j) * aM1(k,l);
            end do
           end do
         end do
       end do
       
      DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3    
       
           calM1(i,j,k,l) = (tempI(i,j,k,l) + (CinvCinv_ICinv(i,j,k,l) 
     &     -Cm1M1M1Cm1(i,j,k,l))*(aJ*aJ/(alamd1*alamd1)) 
     &     + dM1M1d(i,j,k,l)*(alamd1*alamd1) 
     &     - aM1M1(i,j,k,l)*d1p ) *(1.0/d1);
           calMr(i,j,k,l) = (aICinv (i,j,k,l)+ calM1(i,j,k,l)) * (-1.0);
           aL_iso_1(i,j,k,l) = ( calM1(i,j,k,l)*Wa(1) 
     &     + calMr(i,j,k,l)*Wa(3) ) * 2.0;
           aL_iso_2(i,j,k,l) =  aM1(i,j) * aM1(k,l) * Yval(1,1) 
     &      + aM1(i,j) * aMr(k,l) * Yval(1,3)+ aMr(i,j) 
     &     * aM1(k,l) * Yval(3,1)  + aMr(i,j) * aMr(k,l) * Yval(3,3);
           aL_iso (i,j,k,l)= aL_iso_1(i,j,k,l) + aL_iso_2(i,j,k,l) ;
            end do
           end do
          end do
       end do
      end if 



       if(index == 13)then 
  
       d3 = (alamd3+alamd1)*(alamd3+alamd2)
     &       *(alamd3-alamd1)*(alamd3-alamd2);
       DO i=1,3
       DO j=1,3
        aM3(i,j) = (delta(i,j) - Cinv(i,j) * (alamd2*alamd2)) 
     &                /(alamd3+alamd2)/(alamd3-alamd2);
        aMr(i,j) = Cinv (i,j)- aM3(i,j);
        end do
      end do 
      d3p = 4.0 *alamd3*alamd3*alamd3*alamd3                  
     &         - aI1*alamd3*alamd3 - aJ*aJ /(alamd3*alamd3);
      DO i=1,3
        DO j=1,3 
         DO k=1,3
           DO l=1,3
           Cm1M3M3Cm1(i,j,k,l) = Cinv(i,j)*aM3(k,l) + aM3(i,j)*Cinv(k,l);
           dM3M3d(i,j,k,l) = delta(i,j)*aM3(k,l) + aM3(i,j)*delta(k,l);
           aM3M3(i,j,k,l) = aM3(i,j) *aM3(k,l);
            end do
          end do
        end do
      end do
        
       DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
        
         calM3(i,j,k,l) = ( tempI(i,j,k,l) + (CinvCinv_ICinv(i,j,k,l)
     &                   -Cm1M3M3Cm1(i,j,k,l))*(aJ*aJ/(alamd3*alamd3)) 
     &                    + dM3M3d(i,j,k,l)*(alamd3*alamd3) 
     &                    - aM3M3(i,j,k,l)*d3p ) *(1.0/d3);
         calMr(i,j,k,l) = (aICinv(i,j,k,l) + calM3(i,j,k,l))*(-1.0);
        aL_iso_1(i,j,k,l)= ( calM3(i,j,k,l)*Wa(3) 
     &   +calMr(i,j,k,l)*Wa(1) ) * 2.0;
         aL_iso_2(i,j,k,l) = aM3(i,j) * aM3(k,l) * Yval(3,3)  
     & + aM3(i,j) * aMr(k,l) * Yval(3,1)+aMr(i,j) * aM3(k,l) * Yval(1,3) 
     &     + aMr(i,j) * aMr(k,l) * Yval(1,1);
       aL_iso (i,j,k,l)= aL_iso_1(i,j,k,l) + aL_iso_2(i,j,k,l);
            end do
         end do
       end do
      end do
       end if
       if(index == 2) then 
      
      shearmoduo=0.4225;
	do id=1,3
      d11(id)= d2isowOverdlamd2(id);
      dd1(id)= disowOdlambda(id);
	end do
c        G2 = d11(1)*alamdw(2)*alamdw(2)+ dd1(1)*alamdw(2);//greskom prevedeno
        G2 = d11(2)*alamdw(2)*alamdw(2)+ dd1(2)*alamdw(2);
       DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
        
c      Izvadjeno iz OS programa 
       aL_iso(i,j,k,l) = (aICinv(i,j,k,l) 
     &       - CinvCinv(i,j,k,l) /3.0) * G2;
cn      clan=(aICinv(i,j,k,l)
cn     &       - CinvCinv(i,j,k,l)/3.0);
cn	clan2=aI4s(i,j,k,l);
cn	clan1=clan*2*shearmoduo*alamd(2)**(-4.);
cn      aL_iso(i,j,k,l) = (aI4s(i,j,k,l) 
cn     &       - delta(i,j)*delta(k,l)/3.0)*2*shearmoduo*alamd(2)**(-4.);!ISPRAVLJENO PREMA RADU!
	
c      if (aL_iso(i,j,k,l).lt.0.0 ) then
cc	write(3,*)'stao je kod L_iso',i,j,k,l,aL_iso(i,j,k,l)
c	endif
	
            end do
         end do
        end do
       end do
       end if
  

! racuna zapremisnski deo tenzora elasticnosti jednacina (4.102) iz B.jeremica
 
       dWdJ = dvolwOdJ;  
       d2WdJ2 = d2volwOdJ2;
        wj = d2WdJ2*aJ*aJ + aJ*dWdJ;
	
       DO i=1,3
        DO j=1,3
         DO k=1,3
           DO l=1,3
            aL_vol (i,j,k,l)= CinvCinv(i,j,k,l)*wj
     &                 - aICinv (i,j,k,l)*2.0*aJ*dWdJ;
            ctens(i,j,k,l)=aL_iso(i,j,k,l)+aL_vol(i,j,k,l);

c	     if (ctens(i,j,k,l).lt.0.0 ) then
c	write(3,*)'stao je kod ctensa',i,j,k,l,ctens(i,j,k,l)
c	endif
	   
            end do
            end do
            end do
            end do
c	write(3,*)'tenzor aL_iso'
c      call wcten(aL_iso)
c	write(3,*)'tenzor ctens'
c	call wcten(ctens)

		  
            
       end SUBROUTINE Stiffness
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE Stiffness1(Stiff,gr,Stiff1)
!-----------------------------------------------------------------------
!     Ova subroutina je ubacena samo iz razloga da se prilagodi OpenSeesu
!-----------------------------------------------------------------------

       IMPLICIT DOUBLE PRECISION (a - h, o - z)
       DIMENSION  Stiff(3,3,3,3),Stiff1(3,3,3,3)
	 DIMENSION  gr(3,3)
! Racuna deo aL_iso

       ndime=3;

      DO  m= 1, ndime
        DO b = 1, ndime
          DO n = 1, ndime
            DO d = 1, ndime
	      suma=0;
	       DO a = 1, ndime
               DO c = 1, ndime
              suma=suma+gr(m,a)
     &                *gr(n,c)*Stiff(a,b,c,d);
	    
                 END DO
              END DO
	     Stiff1(m,b,n,d)=suma;
             END DO
           END DO
	  END DO
	END DO
        
      
		  
            
       end SUBROUTINE Stiffness1
!-----------------------------------------------------------------------
      subroutine wcten(c)
	implicit double precision(a-h,o-z)
	dimension c(3,3,3,3)
	write(3,2000) c(1,1,1,1),c(1,1,2,2),c(1,1,3,3),
     1              c(1,1,1,2),c(1,1,2,3),c(1,1,3,1)
      write(3,2000) c(2,2,1,1),c(2,2,2,2),c(2,2,3,3),
     1              c(2,2,1,2),c(2,2,2,3),c(2,2,3,1)
	write(3,2000) c(3,3,1,1),c(3,3,2,2),c(3,3,3,3),
     1              c(3,3,1,2),c(3,3,2,3),c(3,3,3,1)
	write(3,2000) c(1,2,1,1),c(1,2,2,2),c(1,2,3,3),
     1              c(1,2,1,2),c(1,2,2,3),c(1,2,3,1)
      write(3,2000) c(2,3,1,1),c(2,3,2,2),c(2,3,3,3),
     1              c(2,3,1,2),c(2,3,2,3),c(2,3,3,1)
	write(3,2000) c(3,1,1,1),c(3,1,2,2),c(3,1,3,3),
     1              c(3,1,1,2),c(3,1,2,3),c(3,1,3,1)
 2000 format(6(1pe12.5))
      return
      end subroutine wcten

