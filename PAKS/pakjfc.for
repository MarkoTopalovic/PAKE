
C=======================================================================
      SUBROUTINE FATCRACK

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ******************************************************************
C     ULAZNE VELICINE						/ JEDINICE /

C     AINT-POCETNA DUZINA PRSLINE			/ m /
C     AFIN-KONACNA DUZINA PRSLINE			/ m /
C     STRMAX-MAXIMALNA VREDNOST NAPONA	/ MPa /
C     STRMIN-MINIMALNA VREDNOST NAPONA	/ MPa /
C     DN-dN U PARISOVOM ZAKONU            / - /
C     PEC-KONSTANTA C U PARISOVOM ZAKONU	/ m**(11/2)/(MN)**3 /
C     PEM-KONSTANTA m U PARISOVOM ZAKONU	/ - /
C     UKC-LOMNA ZILAVOST	                / MN/m**(3/2) / 
C     SS-STEPEN SIGURNOSTI                / - /
C     YA-KOREKCIONI FAKTOR                / - /   
C     ******************************************************************

      DIMENSION A(1000000),DA(1000000),
     1          UN(1000000),DK(1000000),
     1          AZ(100),KZ(100)


C     ======================================================
C     ULAZNI PODACI
C     ======================================================
      !!!!!KARAKTERISTIKE PRSLINE
      AINT=0.01	    !POCETNA DUZINA PRSLINE	   (m)
      AFIN=0.05	    !MAKSIMALNA DUZINA PRSLINE (m)

      !!!!!ISTORIJA OPTERECENJA
      STRMAX=10.0    !MAKSIMALNI NAPON	(MPa)
      STRMIN=0.0     !MINIMALNI NAPON	(MPa)

      !!!!!! PARIS-ov ZAKON
      PEC=0.005E-11   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
      PEM=4.          !KOEFICIJENT m U PARISOVOM ZAKONU (/)

      !!!!!! KRITERIJUM OTKAZA 
      UKC=90.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
      SS=1.5          !STEPEN SIGURNOSTI	 (/)

      !!!!!!!!NACIN ODREDJIVANJA FAKTORA K    
      INTER=1	! INTER=0 ANALITICKIM PUTEM 
                ! INTER=1 INTERPOLACIJA NUMERICKIH VREDNOSTI K
      IF(INTER .EQ.1) THEN
         NNIZ=11  !UNETI VREDNOST
C         DO L=1,NNIZ
         AZ(1)=0.02	!UNETI VREDNOST
         KZ(1)=7.6 !UNETI VREDNOST
         AZ(2)=0.025
         KZ(2)=8.48 
         AZ(3)=0.030
         KZ(3)=9.32 
         AZ(4)=0.035
         KZ(4)=10.09 
         AZ(5)=0.040
         KZ(5)=11.06 
         AZ(6)=0.045
         KZ(6)=11.92 
         AZ(7)=0.050
         KZ(7)=12.96 
         AZ(8)=0.055
         KZ(8)=14.00 
         AZ(9)=0.060
         KZ(9)=15.67 
         AZ(10)=0.065
         KZ(10)=17.90 
         AZ(11)=0.070
         KZ(11)=20.8 

C         ENDDO
         AINT=AZ(1)
         AFIN=AZ(NNIZ)
      ENDIF

      ! IF(INTER.EQ.0) THEN
      !!!!!! IZBOR OBLIKA PRSLINE
      !!!!!! IPRSLINA=1 CENTRALNA 
      !!!!!!         =2 IVICNA
      !!!!!!         =3 NA KRUZNOM OTVORU
      !!!!!!         =4 IVICNA ELIPTICNA PRSLINA
      !!!!!!         =5 UNUTRASNJA ELIPTICNA PRSLINA 
      ! IPRSLINA=3
      ! IF(IPRSLINA.EQ.3) THEN
!    PP=0.08       !UNETI VREDNOST POLUPRECNIKA KRUZNOG OTVORA/M/
!	    PC=0.00        !VECA POLUOSA ELIPSE 
!	  ELSEIF(IPRSLINA .EQ.4 .OR. IPRSLINA.EQ.5) THEN
!	    PP=0.080
!	    PC=0.004      !UNETI VREDNOST VECE POLUOSE ELIPSE
!        ELSE 
!	    PP=0.000		!ZA SLUCAJEVE IPRSLINA=1,2
!	    PC=0.000
!        ENDIF	   
!	ENDIF
	

	!!!!!! VRSTA ZAKONA KOJA SE KORISTI
	!!!!!!    IN=0	PARISOV ZAKON
	!!!!!!    IN=1  KORIGOVANI PARISOV ZAKON
	!!!!!!    IN=2  WALKER-ov ZAKON
	
      IN=0     
      IF(IN.EQ.1) THEN
         DKTH=6.0  !UNETI VREDNOST PRAGA SIF-a (MN/m**(3/2))
      ENDIF
      IF(IN.EQ.2) THEN
         GAMA=0.7
         R=STRMIN/STRMAX
      ENDIF 

      !!!!!! PROVERA VALIDNOSTI LEML
      !!!!!! DA LEML=1
      !!!!!! NE LEML=0
      write(6,100)
  100 FORMAT(/' UNESITE BROJ PRIMERA (4,9,7,10)')
      READ(5,200) IPRIMER
  200 FORMAT(I2)    
	
C     IPRIMER=10
      IF(IPRIMER.EQ.4.OR. IPRIMER.EQ.9) THEN
         INTER=0 !ODREDJIVANJE FAK. K ANALITICKIM PUTEM
         IPRSLINA=5
         DKTH=6.
         LEML=1
         PC=0.002   !VECA POLUOSA ELIPSE

         AINT=0.001	    !POCETNA DUZINA PRSLINE	   (m)
         AFIN=0.055	    !MAKSIMALNA DUZINA PRSLINE (m)

         STRMAX=400.0	    !MAKSIMALNI NAPON	(MPa)
         STRMIN=200.0      !MINIMALNI NAPON	(MPa)

         PEC=0.05E-11   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
         PEM=4.        !KOEFICIJENT m U PARISOVOM ZAKONU (/)
         SIGMAY=1200
         !!!!!! KRITERIJUM OTKAZA 
         UKC=60.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
         SS=1.4          !STEPEN SIGURNOSTI	 (/)

         IF(IPRIMER.EQ.4) THEN
            IN=0  !PARISOV YAKON 
         ELSEIF(IPRIMER.EQ.9) THEN
            IN=2
            GAMA=0.7
            R=STRMIN/STRMAX	       
         ENDIF

      ENDIF

      IF(IPRIMER.EQ.7) THEN
         INTER=0 !ODREDJIVANJE FAK. K ANALITICKIM PUTEM
         IPRSLINA=3
         DKTH=6.
         LEML=1
         PP=0.08	!POLUPRECNIK OTVORA
         PC=0.0   !VECA POLUOSA ELIPSE

         AINT=0.08	    !POCETNA DUZINA PRSLINE	   (m)
         AFIN=0.1 	    !MAKSIMALNA DUZINA PRSLINE (m)

         SIGMAY=600

         STRMAX=60.0	    !MAKSIMALNI NAPON	(MPa)
         STRMIN=12.0      !MINIMALNI NAPON	(MPa)

         PEC=1.0E-11   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
         PEM=3.        !KOEFICIJENT m U PARISOVOM ZAKONU (/)

	 !!!!!! KRITERIJUM OTKAZA 
         UKC=90.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
         SS=1.5          !STEPEN SIGURNOSTI	 (/)
         IN=0		  !PARISOV YAKON 
      ENDIF

      IF(IPRIMER.EQ.10) THEN
         INTER=0 !ODREDJIVANJE FAK. K ANALITICKIM PUTEM
         IPRSLINA=2
         DKTH=6.
         LEML=0
         PP=0.00	!POLUPRECNIK OTVORA
         PC=0.0   !VECA POLUOSA ELIPSE

         AINT=0.001	    !POCETNA DUZINA PRSLINE	   (m)
         AFIN=0.081 	    !MAKSIMALNA DUZINA PRSLINE (m)

         SIGMAY=600

         STRMAX=100.0	    !MAKSIMALNI NAPON	(MPa)
         STRMIN=0.0      !MINIMALNI NAPON	(MPa)

         PEC=0.250E-11   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
         PEM=3.2        !KOEFICIJENT m U PARISOVOM ZAKONU (/)

         !!!!!! KRITERIJUM OTKAZA 
         UKC=90.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
         SS=1.4          !STEPEN SIGURNOSTI	 (/)
         IN=1		  !MODIFIKOVANI PARISOV YAKON 
      ENDIF
	     

      IF(LEML .EQ.1.AND. INTER.EQ.0) THEN
         Rpz=2.5*(UKC/(SS*SIGMAY))**2.
         A0=AINT 
         CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)
         Acrit=((UKC/(SS*YA*STRMAX))**2.)/3.14
         IF(IPRSLINA.EQ.3) THEN
             Acrit=(UKC/(SS*YA*STRMAX*(PP**(1./6.))*SQRT(3.14)))**3.
         ENDIF
         WRITE(6,*) 'POLUPRECNIK PLASTICNE ZONE Rpz',Rpz
         WRITE(6,*) 'KRITICNA DUZINA PRSLINE Acrit',Acrit 
         IF(Acrit.GT.Rpz) THEN 
            WRITE(6,*) 'VAZI LEML JER Acrit>Rpz'
            PAUSE
         ELSE
            WRITE(6,*) 'NE VAZI LEML JER Acrit<Rpz'
            PAUSE
         ENDIF
   
         IF(Acrit .GT. AFIN .AND. INTER.EQ.0) THEN
            WRITE(6,*)'KRITICNA DUZINA JE VECA OD KONACNO ZADATE'
            WRITE(6,*)'DA LI HOCETE KOREKCIJU KONACNE DUZINE PRSLINE'
            !!!!!!! KOR=1 DA / KOR=0 NE
            KOR=1
            IF(KOR.EQ.1) THEN
               AFIN=Acrit
            ENDIF
         ELSEIF(Acrit .LT. AFIN .AND. INTER.EQ.0) THEN
            WRITE(6,*)'KRITICNA DUZINA JE MANJA OD KONACNO ZADATE'
            WRITE(6,*)'NEOPHODNA JE KOREKCIJA KONACNE DUZINE PRSLINE'
            AFIN=Acrit
         ENDIF
      ENDIF
      !!!!!! KORAK BROJA CIKLUSA
      !ZADAJE SE	    (/)
   
C     ======================================================
C     ======================================================
      DSTR=STRMAX-STRMIN
   
      IF(DSTR.LT.99) THEN
         DN=100.
      ELSE
         DN=10.
      ENDIF	 	   
C     INICIJALIZACIJA POCETNIH VREDNOSTI
      NCO=1
      A(NCO)=AINT
      UN(NCO)=0.0
      A0=A(NCO)
      IF(INTER.EQ.0) THEN
         CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)
         DK(NCO)=DSTR*YA*SQRT(3.14*A(NCO))
      ELSEIF(INTER.EQ.1) THEN
         CALL INTERPOLACIJA(A0,DKK,AZ,KZ,NNIZ)
         DK(NCO)=DKK 
      ENDIF
      IF(IN.EQ.0)THEN
         DA(NCO)=PEC*DN*(DK(NCO)**PEM)
      ELSEIF(IN.EQ.1) THEN 
         DA(NCO)=PEC*DN*(DK(NCO)**PEM-DKTH**PEM)
      ELSEIF(IN.EQ.2) THEN 
         DA(NCO)=PEC*DN*DK(NCO)**PEM/((1.-R)**(PEM-PEM*GAMA))
      ENDIF
   
      IF(IN.EQ.1 .AND.DK(1).GT.DKTH) THEN
         WRITE (6,*) 'DOLAZI DO RASTA PRSLINE K(a1)>Kth',DK(1),DKTH
         PAUSE
      ENDIF    
C     ======================================================
C     ======================================================
      
   10 IF(DK(NCO).LT.UKC/SS .AND. NCO.LT.1000000) THEN
         GOTO 20
   
      ELSE
         GOTO 30
   
      ENDIF
	  
   20 IF(A(NCO) .LT. AFIN)  THEN
         NCO=NCO+1
         A(NCO)=A(NCO-1)+DA(NCO-1)
         A0=A(NCO)
         IF(INTER.EQ.0) THEN
            CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)
            DK(NCO)=DSTR*YA*SQRT(3.14*A(NCO))
         ELSEIF(INTER.EQ.1) THEN
            CALL INTERPOLACIJA(A0,DKK,AZ,KZ,NNIZ)
            DK(NCO)=DKK 
         ENDIF
   
         UN(NCO)=UN(NCO-1)+DN       
         IF(IN.EQ.0) THEN
            DA(NCO)=PEC*DN*(DK(NCO)**PEM)
         ELSEIF(IN.EQ.1) THEN 
            DA(NCO)=PEC*DN*(DK(NCO)**PEM-DKTH**PEM)
         ELSEIF(IN.EQ.2) THEN 
            DA(NCO)=PEC*DN*DK(NCO)**PEM/((1.-R)**(PEM-PEM*GAMA))
         ENDIF
         GOTO 10
      
      ELSE 
         GOTO 30  
   
      ENDIF
	  
   30 WRITE(6,*) 'KONACNA DUZINA PRSLINE I BROJ CIKLUSA' 
      WRITE(6,*) 'NA OSNOVU KRITERIJUMA OTKAZA'
      write(6,*) 'NCO A(NCO) UN(NCO) DK(NCO) UKC'	
   
      
      write(6,*)  NCO,A(NCO),UN(NCO),DK(NCO),UKC
	  
      NNCO=NCO/20
      WRITE(6,*) 'NNCO',NNCO
      NDIJ=0
      WRITE(3,*)'RAST PRSLINE        N                            a'
      DO N=1,NCO,NNCO
         WRITE(3,2100) UN(N), A(N)
         NDIJ=NDIJ+1
      ENDDO
C     UPISIVANJE POSLEDNJEG CLANA NIZA
      WRITE(3,2100) UN(NCO), A(NCO)
      NDIJ=NDIJ+1
      write(6,*)'BROJ UPISANIH CLANOVA NIZA',NDIJ
	  
      NDIJ1=0
      WRITE(3,*)'BRZINA RASTA PRSLINE   dK                     da/dN'
      DO N=1,NCO,NNCO
         WRITE(3,2100) DK(N),DA(N)/DN 
         NDIJ1=NDIJ1+1
      ENDDO
C     UPISIVANJE POSLEDNJEG CLANA NIZA
      NDIJ1=NDIJ1+1
      WRITE(3,2100) DK(NCO), DA(NCO)/DN
      write(6,*)'BROJ UPISANIH CLANOVA NIZA',NDIJ1
      pause
C     ===================================================================
C     ===================================================================
c     PRAVLJENJE DATOTEKE ZA TEC-PLOT
      OPEN(UNIT=10,FILE='RAST PRSLINE.DAT') 
      WRITE(10,*)'TITLE = "RAST PRSLINE USLED MONOTONOG OPTERECENJA"' 
      WRITE(10,*)'VARIABLES = "N", "a"'
      WRITE(10,1030) NDIJ
      DO N=1,NCO,NNCO
         WRITE(10,2100) UN(N),A(N)
      ENDDO 	
      WRITE(10,2100) UN(NCO),A(NCO)
      CLOSE(10)
      OPEN(UNIT=10,FILE='BRZINA RASTA PRSLINE.DAT') 
      WRITE(10,*)'TITLE ="BRZINA RASTA PRI MONOTONOM OPTERECENJU"'
      WRITE(10,*)'VARIABLES = "dK", " da/dN"'
      WRITE(10,1040) NDIJ1 	
      DO N=1,NCO,NNCO
         WRITE(10,2100) DK(N),DA(N)/DN
      ENDDO 	
      WRITE(10,2100) DK(NCO),DA(NCO)/DN
1030  FORMAT('ZONE T="ZONE N-a", I=',I2,' F=POINT')
1040  FORMAT('ZONE T="ZONE dK-da/dN", I=',I2,' F=POINT')
2100  FORMAT(20X,1PE10.3,20X,1PE10.3)
      CLOSE(10)      
C     ===================================================================
C     ===================================================================
      return

      end
C============================================================
C============================================================
      SUBROUTINE INTERPOLACIJA(A0,DKK,AZ,KZ,NNIZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION AZ(100),KZ(100)

      DO L=1,NNIZ-1
         IF(A0.GT.AZ(L) .AND. A0.LT. AZ(L+1)) THEN
            DKK=KZ(L)+(KZ(L+1)-KZ(L))*(A0-AZ(L))/(AZ(L+1)-AZ(L))
         ELSEIF(A0.EQ.AZ(L)) THEN
            DKK=KZ(L)
         ELSEIF(A0.EQ.AZ(L+1)) THEN
            DKK=KZ(L+1)
         ENDIF
      ENDDO

      RETURN
      END	  		 		  		  
C============================================================
C============================================================
