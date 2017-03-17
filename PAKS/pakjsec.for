
C=======================================================================
      SUBROUTINE FATSEC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ******************************************************************
C     ******************************************************************
C     ZAMORNI RAST RAST PRSLINE USLED KOMBINOVANOG OPTERECENJA
C     ******************************************************************
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
     1          SIG0(4),SIG(4),Ath0(4),Ath(4),
     1          APOC(5),ACON(5),USEC(5),IBB(5),NCRT(5),
     1          NMIN(4),NS0(4),NS(4),RR(4)
      
      DIMENSION UNN(4),UN1(100),A1(100),DK1(100),DA1(100)

C     ======================================================
C     ULAZNI PODACI
C     ======================================================
      write(6,100)
  100 FORMAT(/' UNESITE BROJ PRIMERA (8,11,12)')
      READ(5,200) IPRIMER
  200 FORMAT(I2)    

      !!!!!KARAKTERISTIKE PRSLINE
      AINT=0.001  	    !POCETNA DUZINA PRSLINE	   (m)
      AFIN=0.005  	    !MAKSIMALNA DUZINA PRSLINE (m)

      !!!!!ISTORIJA OPTERECENJA
      NSEC=4
      SIG0(1)=285.0
      SIG0(2)=270.0
      SIG0(3)=200.0
      SIG0(4)=92.0
      NS0(1)=1
      NS0(2)=3
      NS0(3)=12
      NS0(4)=44

      !!!!!! PARIS-ov ZAKON
      PEC=0.25E-11   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
      PEM=3.2        !KOEFICIJENT m U PARISOVOM ZAKONU (/)

      !!!!!! KRITERIJUM OTKAZA 
      UKC=90.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
      SS=1.5          !STEPEN SIGURNOSTI	 (/)
      DKTH=8.0        !UNETI VREDNOST PRAGA SIF-a (MN/m**(3/2))
      IN=0
      IPRSLINA=2
      !!!!!! IZBOR OBLIKA PRSLINE
      !!!!!! IPRSLINA=1 CENTRALNA 
      !!!!!!	   =2 IVICNA
      !!!!!!       =3 NA KRUZNOM OTVORU
      !!!!!!       =4 IVICNA ELIPTICNA PRSLINA
      !!!!!!       =5 UNUTRASNJA ELIPTICNA PRSLINA 

      IF(IPRIMER.EQ.11.OR.IPRIMER.EQ.8) THEN
         !!!!!KARAKTERISTIKE PRSLINE
         AINT=0.002  	    !POCETNA DUZINA PRSLINE	   (m)
         AFIN=0.006  	    !MAKSIMALNA DUZINA PRSLINE (m)
         
         !!!!!ISTORIJA OPTERECENJA
         NSEC=2
         SIG0(1)=200.0
         SIG0(2)=400.0
         NS0(1)=1
         NS0(2)=1
         
         !!!!!! PARIS-ov ZAKON
         PEC=1.0 E-13   !KOEFICIJENT C U PARISOVOM ZAKONU (e-11)
         PEM=4.        !KOEFICIJENT m U PARISOVOM ZAKONU (/)
         
         !!!!!! KRITERIJUM OTKAZA 
         UKC=70.0        !LOMNA ZILAVOST		 (MN/m**(3/2))
         SS=1.5          !STEPEN SIGURNOSTI	 (/)
         DKTH=6.0        !UNETI VREDNOST PRAGA SIF-a (MN/m**(3/2))
         IN=0
         IPRSLINA=4
         PC=0.004
      ENDIF

      IF(IPRSLINA.EQ.3) THEN
         PP=0.08       !UNETI VREDNOST POLUPRECNIKA KRUZNOG OTVORA/M/
         PC=0.00        !VECA POLUOSA ELIPSE 
      ELSEIF(IPRSLINA .EQ.4 .OR. IPRSLINA.EQ.5) THEN
         PP=0.000
         PC=0.004      !UNETI VREDNOST VECE POLUOSE ELIPSE
      ELSE 
         PP=0.000		!ZA SLUCAJEVE IPRSLINA=1,2
         PC=0.000
      ENDIF	   

      !!!!!! VRSTA ZAKONA KOJA SE KORISTI
      !!!!!!    IN=0	PARISOV ZAKON
      !!!!!!    IN=1  KORIGOVANI PARISOV ZAKON
      !!!!!!    IN=2  WALKER-ov ZAKON
	
	     
      IF(IN.EQ.1) THEN
         DKTH=8.0	  !UNETI VREDNOST PRAGA SIF-a (MN/m**(3/2))
      ENDIF
      IF(IN.EQ.2) THEN
         GAMA=0.7
         DO I=1,NSEC
            RR(I)=0.0     !!!!!KOEFICIJENT SREDNJEG NAPONA
         ENDDO
      ENDIF     

C     DEFINISANJE POCETNE VREDNOSTI KOREKCIONE FUNKCIJE  
      A0=AINT 
      CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)

	
      !!!DEFINISANJE MINIMALNE DUZINE PRSLINE PRI KOJOJ DOLAZI DO 
      !!!RASTA 
      PES=2.
      DO I=1,NSEC
         Ath0(I)=((DKTH/(SIG0(I)*YA))**PES)/3.14
      ENDDO
      WRITE(6,*) 'I,Ath0(I)' 
      DO I=1,NSEC
         WRITE(6,*)  I,Ath0(I)
      ENDDO
C      PAUSE   

      DO J=1,NSEC
         AMIN=Ath0(J)
         Ath(J)=AMIN 
         NMIN(J)=J
         DO I=1,NSEC
            IF(AMIN.GT.Ath0(I)) THEN
               AMIN=Ath0(I)
               Ath(J)=AMIN
               NMIN(J)=I
            ENDIF
         ENDDO
         Ath0(NMIN(J))=1000000.
      ENDDO

      WRITE(6,*) 'I,Ath(I),NMIN(I)' 
      DO I=1,NSEC
         WRITE(6,*)  I,Ath(I),NMIN(I)
      ENDDO
C      PAUSE
	   
      DO I=1,NSEC
         SIG(I)=SIG0(NMIN(I))
         NS(I)=NS0(NMIN(I))
         WRITE(6,*) I,SIG(I),NS(I)
      ENDDO
C      PAUSE       
      DO I=1,NSEC
         IF(Ath(I).GT.AFIN) THEN
            AFIN=1.05*Ath(I)
            WRITE(6,*)'KOREKCIJA KRAJNJE DUZINE PRSLINE NA',AFIN
         ENDIF
      ENDDO  
C     ======================================================
C     ======================================================
C     DEFINISANJE POJEDINACNIH SEKVENCI
      IB=0
      APOC(1)=AINT
      DO I=1,NSEC
         IF(Ath(I).GT.AINT .AND. Ath(I).LT.AFIN) THEN
            IB=IB+1
            IF(IB.EQ.1 .AND.I.LT.NSEC) THEN
               APOC(IB)=AINT
               ACON(IB)=Ath(I)
               IBB(IB)=I
            ELSEIF(IB.GT.1 .AND.I.LT.NSEC) THEN
               APOC(IB)=ACON(IB-1)
               ACON(IB)=Ath(I)
               IBB(IB)=I
            ELSEIF(IB.GT.1 .AND.I.EQ.NSEC) THEN
               APOC(IB)=ACON(IB-1)
               ACON(IB)=Ath(I)
               IBB(IB)=I
               IB=IB+1
               APOC(IB)=ACON(IB-1)
               ACON(IB)=AFIN
               IBB(IB)=NSEC+1
            ELSEIF(IB.EQ.1 .AND.I.EQ.NSEC) THEN
               APOC(IB)=AINT
               ACON(IB)=Ath(I)
               IBB(IB)=I
               IB=IB+1
               APOC(IB)=ACON(IB-1)
               ACON(IB)=AFIN
               IBB(IB)=NSEC+1
            ENDIF
         ENDIF
      ENDDO
      write(6,*)'IB',IB
      IF(IB.EQ.0) THEN
         IB=IB+1
         APOC(IB)=AINT
         ACON(IB)=AFIN
         IBB(IB)=NSEC+1
      ENDIF
	
      WRITE(6,*) 'BROJ PODSEKVENCI',IB		      	          
      WRITE(6,*) 'I,APOC(I),ACON(I),IBB(I)'
      DO I=1,IB
         WRITE(6,*) I,APOC(I),ACON(I),IBB(I)
      ENDDO
      PAUSE

C     ======================================================
C     ======================================================
      DO KK=1,IB 
C        PODELA NA PODSEKVENCE
         AINT=APOC(KK)
         AFIN=ACON(KK)
C        EKVIVALENTNI NAPON UKUPNI I BROJ CIKLUSA PO SEKVENCI
         DSTR=0.
         NUK=0
         R=0.0
         DO I=1,IBB(KK)-1                          !!!!!!!NSEC
            DSTR=DSTR+NS(I)*(SIG(I)**PEM)
            NUK=NUK+NS(I)
            IF(IN.EQ.2) THEN
               R=R+NS(I)*(RR(I)**PEM)
            ENDIF
         ENDDO
         DSTR=(DSTR/NUK)**(1./PEM)
         R=(R/NUK)**(1./PEM)
         WRITE(6,*) '          '
         WRITE(6,*) 'PODSEKVENCA',KK
         WRITE(6,*) 'EKV NAPON=',DSTR
C         PAUSE

         IF(DSTR.LT.99) THEN
            DN=100.
         ELSE
            DN=10.
         ENDIF	 	   
C        INICIJALIZACIJA POCETNIH VREDNOSTI
C         IF(KK.EQ.1) THEN
         NCO=1         !!!!!ZA PRVU PODSEKVENCU
         UN(NCO)=0.0
C         ELSE
C            NCO=NCO+1
C            UN(NCO)=UN(NCO)
C         ENDIF
         A(NCO)=AINT
         A0=A(NCO) 
         CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)
C         UN(NCO)=0.0

         DK(NCO)=DSTR*YA*SQRT(3.14*A(NCO))

         IF(IN.EQ.0) THEN
            DA(NCO)=PEC*DN*(DK(NCO)**PEM)
         ELSEIF(IN.EQ.1) THEN 
            DA(NCO)=PEC*DN*(DK(NCO)**PEM-DKTH**PEM)
         ELSEIF(IN.EQ.2) THEN 
            DA(NCO)=PEC*DN*DK(NCO)**PEM/((1.-R)**(PEM-PEM*GAMA))
         ENDIF

C        ======================================================
C        ======================================================
   
   10    IF(DK(NCO).LT.UKC/SS) THEN
            GOTO 20

         ELSE
            GOTO 30

         ENDIF
	
   20    IF(A(NCO) .LT. AFIN)  THEN
            NCO=NCO+1
            A(NCO)=A(NCO-1)+DA(NCO-1)
            A0=A(NCO) 
            CALL KOREKCIJA(IPRSLINA,PP,PC,A0,YA)
            DK(NCO)=DSTR*YA*SQRT(3.14*A(NCO))
            UN(NCO)=UN(NCO-1)+DN       
C            DA(NCO)=PEC*DN*(DK(NCO)**PEM)
            IF(IN.EQ.0)	THEN
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
C        BROJ PODSEKVENCI
   30    USEC(KK)=UN(NCO)/NUK 
         UNN(KK)=UN(NCO)
       WRITE(6,*)'KONACNA DUZINA PRSLINE,BROJ CIKLUSA PODSEKVENCE I SIF' 
       write(6,*)' A(NCO)                 UN(NCO)               DK(NCO)'	
         write(6,*)   A(NCO),UN(NCO),DK(NCO)
         WRITE(6,*) 'BROJ CIKLUSA U PODSEKVENCI',USEC(KK)

C        PRIPREMA NIZA ZA UPISIVANJE VREDNOSTI U DATOTEKU	
	 NNCO=NCO/20
         WRITE(6,*) 'NNCO',NNCO
         NDIJ1=0
         IF(KK.EQ.1) THEN
            NDIJ=0
         ENDIF

         WRITE(3,*)'PODSEKVENCA',KK
       WRITE(3,*)'BRZINA RASTA PRSLINE     dK                     da/dN'
         DO N=1,NCO,NNCO
            NDIJ=NDIJ+1
            NDIJ1=NDIJ1+1
            IF(KK.EQ.1) THEN
               UN1(NDIJ)=UN(N)
            ELSE
               UN1(NDIJ)=UN(N)+UNN(KK-1)
            ENDIF
            A1(NDIJ)=A(N)
            DK1(NDIJ)=DK(N)
            DA1(NDIJ)=DA(N)/DN
            WRITE(3,2100) DK(N),DA(N)/DN
         ENDDO
C        FORMIRANJE POSLEDNJEG CLANA NIZA
         NDIJ=NDIJ+1
         NDIJ1=NDIJ1+1
         WRITE(3,2100) DK(NCO),DA(NCO)/DN	
         IF(KK.EQ.1) THEN
            UN1(NDIJ)=UN(NCO)
         ELSE
            UN1(NDIJ)=UN(NCO)+UNN(KK-1)
         ENDIF
         A1(NDIJ)=A(NCO)
         DK1(NDIJ)=DK(NCO)
         DA1(NDIJ)=DA(NCO)/DN
	
C        BROJ CLANOVA NIZA ZA CRTANJE PO PODSEKVENCAMA
         NCRT(KK)=NDIJ1 
         write(6,*)'BROJ UPISANIH CLANOVA NIZA',NDIJ
         pause

	
      ENDDO   !!!ZAVRSENA PETLJA PO PODSEKVENCAMA


C     UKUPAN BROJ SEKVENCI
      UBS=0.0
      DO KK=1,IB
         UBS=UBS+USEC(KK)
      ENDDO
      WRITE(6,*) '                    '
      WRITE(6,*) 'UKUPAN BROJ SEKVENCI',UBS

      write(6,*)'BROJ CLANOVA NIZA',NDIJ
      pause

C     UPISIVANJE VREDNOSTI U DATOTEKU
      WRITE(3,*)'RAST PRSLINE        N                             a'
      DO N=1,NDIJ
         WRITE(3,2100) UN1(N), A1(N)
      ENDDO

C     ===================================================================
C     ===================================================================
c     PRAVLJENJE DATOTEKE ZA TEC-PLOT
      OPEN(UNIT=10,FILE='RAST PRS SPEC.DAT') 
      WRITE(10,*)'TITLE = "RAST PRSLINE USLED SPEKTRA OPTERECENJA"' 	
      WRITE(10,*)'VARIABLES = "N", "a"'
      WRITE(10,1030) NDIJ
      DO N=1,NDIJ
         WRITE(10,2100) UN1(N),A1(N)
      ENDDO 	
      CLOSE(10)

      OPEN(UNIT=10,FILE='BRZ RAS PRS SPEC.DAT') 
      WRITE(10,*)'TITLE ="BRZINA RASTA PRSLINE PRI SPEKTRU OPTERECENJA"' 	
      WRITE(10,*)'VARIABLES = "dK", "da/dN"'
      DO KK=1,IB
         IP=NCRT(KK)
         WRITE(10,1040) KK,IP 
         IF(KK.EQ.1) THEN
            IPOC=1
            IKRAJ=NCRT(1)
         ELSE
            IPOC=NCRT(KK-1)+1
            IKRAJ=IPOC+NCRT(KK)-1
         ENDIF
         DO N=IPOC,IKRAJ
            WRITE(10,2100) DK1(N),DA1(N)
         ENDDO
      ENDDO
      CLOSE(10)

2100  FORMAT(20X,1PE10.3,20X,1PE10.3)
1030  FORMAT('ZONE T="ZONE N-a", I=',I2,' F=POINT')
1040  FORMAT('ZONE T="PODSEKVENCA',I2,'" I=',I2,' F=POINT')
      return
      end
C============================================================
