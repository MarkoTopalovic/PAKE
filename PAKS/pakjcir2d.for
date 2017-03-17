
C=======================================================================
      SUBROUTINE S_FUN(CORD,NUM,PSI,NERING,MIE,QST,Q,
     1           NP,NE,NCVE,XL,YL,XGG,YGG,NNOD,MODEL)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ******************************************************************
C     DEFINISANJE ELEMENATA PO KOJIMA SE VRSI INTEGRACIJA J-INTEGRALA I
C     RACUNANJE TEZINSKE FUNK.      
C     ******************************************************************
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /SRPSKI/ ISRPS
      DIMENSION VX1(10,2),VX2(10,2),PP(4)
      DIMENSION CORD(NP,3),NUM(NE,*),PSI(NCRACK,MAXRIN,NE,MAXNOD),
     1          NERING(NCRACK,MAXSEG,MAXRIN),
     1          MIE(NCRACK,MAXSEG,MAXRIN,N100),
     1          QST(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          Q(NP,MAXRIN),
     1          XL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          XGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),NNOD(NE)

C     ======================================================
C     Ulazne velicine u potprogram su :
C        NP-broj cvorova u domenu 
C        CORD[NP,1],CORD[NP,2]-koordinate cvorova;
C        NUM[NE,NCVE]-numeracija cvorova po elementima;
C        NE-broj elemenata;
C        NCVE-broj cvorova po elementu
C        Xc i Yc-koordinate centra kruznice;
C        R-je poluprecnik kruznice-integracione putanje;          
C        NCR - TEKUCI BROJ PRSLINE
C        NCS - TEKUCI BROJ SEGMENTA PO DUBINI PRSLINE (U 2D NCS=1)
C        NSEG-MAX. BROJ SEGMENTATA (u 2D NSEG=1)
C        N100-MAX. BROJ ELEMENATA U SEGMENTU (N100=100)
C        IRNIG-BROJ PRSTENOVA PO KOJIMA SE VRSI INTEGRACIJA
C     =======================================================
C     =======================================================
C     Izlazne velicine potprograma su: 
C     NERING(NCR,NCS,KR)-broj elemenata u datom prstenu KR
C     MIE(NCR,NCS,KR,NERING(NCR,NCS,KR))-globalna numeracija elemenata
c         koji pripadaju prstenu KR
C     QST(NCR,NCS,KR,i,k)-vrednost funkcije Q po nodovima k,
c         elemenata i, u KR-tom prstenu
C     ======================================================
C     INICIJALIZACIJA VREDNOSTI NIZA 
      CALL ICLEAR(NERING,NCRACK*MAXSEG*MAXRIN)

      DO NCR=1,NCRACK	      !PETLJA PO PRSLINAMA
         MODEL=NODCR(NCR,10)
         NSEG=NODCR(NCR,11)
         IRING=NODCR(NCR,12)

C        DUZINA PRSLINE (VRH PRSLINE NODCR(*,3) KRAJ LICA NODCR(*,4))
         DUZ=DSQRT((CORD(NODCR(NCR,3),1)-CORD(NODCR(NCR,4),1))**2+
     1             (CORD(NODCR(NCR,3),2)-CORD(NODCR(NCR,4),2))**2)
         WRITE(3,*) ' '
         if (isrps.eq.0) then
            WRITE(6,*) ' PRSLINA',NCR,' DUZINA',DUZ
            WRITE(3,*) ' PRSLINA',NCR,' DUZINA',DUZ
         else
            WRITE(6,*) ' CRACK',NCR,' LENGTH',DUZ
            WRITE(3,*) ' CRACK',NCR,' LENGTH',DUZ
         endif
C	  PAUSE          	  
C        POLUPRECNICI PRSTENOVA PO KOJIMA SE VRSI INTEGRACIJA	  
c              pd1-pd4
C         PP(1)=0.20*DUZ   !NA 20% DUZINE PRSLINE
C         PP(2)=0.30*DUZ   !NA 30% DUZINE PRSLINE
C         PP(3)=0.40*DUZ   !NA 40% DUZINE PRSLINE 
C         PP(4)=0.50*DUZ   !NA 50% DUZINE PRSLINE
         PP(1)=DUZ*NODCR(NCR,6)/100   !NA 20% DUZINE PRSLINE
         PP(2)=DUZ*NODCR(NCR,7)/100   !NA 30% DUZINE PRSLINE
         PP(3)=DUZ*NODCR(NCR,8)/100   !NA 40% DUZINE PRSLINE 
         PP(4)=DUZ*NODCR(NCR,9)/100   !NA 50% DUZINE PRSLINE

         !JEDINICNI VEKTORI PRAVCA-X1 PRSLINE NCR
         VX1(NCR,1)=CORD(NODCR(NCR,3),1)-CORD(NODCR(NCR,4),1)
         VX1(NCR,2)=CORD(NODCR(NCR,3),2)-CORD(NODCR(NCR,4),2)
         VV=dsqrt(VX1(NCR,1)*VX1(NCR,1)+VX1(NCR,2)*VX1(NCR,2))
         VX1(NCR,1)=VX1(NCR,1)/VV
         VX1(NCR,2)=VX1(NCR,2)/VV
         CONTE=VX1(NCR,1)
         SINTE=VX1(NCR,2)
         !JEDINICNI VEKTORI PRAVCA-X2 PRSLINE NCR
         VX2(NCR,1)=-1.D0*VX1(NCR,2)
         VX2(NCR,2)=VX1(NCR,1)

         WRITE(6,*)'NCR,VX1(1),VX1(2)',NCR,VX1(NCR,1),VX1(NCR,2)
         WRITE(6,*)'NCR,VX2(1),VX2(2)',NCR,VX2(NCR,1),VX2(NCR,2)
C        KOORDINATE VRHA PRSLINE 
         Xc=CORD(NODCR(NCR,3),1)
         Yc=CORD(NODCR(NCR,3),2)
C         dex=dabs(VX1(NCR,1)*(CORD(NODCR(NCR,4),1)-Xc)+
C     1            VX2(NCR,1)*(CORD(NODCR(NCR,4),2)-Yc))
C         dey=dabs(VX1(NCR,2)*(CORD(NODCR(NCR,6),1)-Xc)+
C     1            VX2(NCR,2)*(CORD(NODCR(NCR,6),2)-Yc))
         dex=DSQRT((CORD(NODCR(NCR,4),1)-Xc)**2.+
     1             (CORD(NODCR(NCR,4),2)-Yc)**2.) 
c         dey=DSQRT((CORD(NODCR(NCR,6),1)-Xc)**2.+
c     1             (CORD(NODCR(NCR,6),2)-Yc)**2.) 

C         R1=dex

         WRITE(3,*)
         DO NCS=1,NSEG	      !PETLJA PO SEGMENTIMA
            DO KR=1,IRING     !PETLJA PO PRSTENOVIMA
	     
C               CALL IWRR(NNOD,NE,'NCV2') 
C               R=PP(KR)*R1-0.001D0*R1    !Tekuci poluprecnik int. putanje
               R=PP(KR)
               if(isrps.eq.0.AND.NCS.EQ.1) then
                  WRITE(6,*) ' PRSTEN',KR,' POLUPRECNIK',R
                  WRITE(3,*) ' PRSTEN',KR,' POLUPRECNIK',R
               else
                  WRITE(6,*) ' RING',KR,' RADIUS',R
                  WRITE(3,*) ' RING',KR,' RADIUS',R
               endif

               do i=1,NE
                  do k=1,NNOD(i)
                     PSI(NCR,KR,i,k)=
     1            dsqrt(DABS((CORD(NUM(i,k),1)-Xc)*(CORD(NUM(i,k),1)-Xc)
     1                  +(CORD(NUM(i,k),2)-Yc)*(CORD(NUM(i,k),2)-Yc)))-R
                  enddo
               enddo

	
C              Identifikacija elemenata koji su preseceni sa int.putanjom
               do i=1,NE
                  PSIMIN=PSI(NCR,KR,i,1)
                  PSIMAX=PSI(NCR,KR,i,1)
                  do k=1,NNOD(i)
                     IF(PSIMAX.LT.PSI(NCR,KR,i,k))PSIMAX=PSI(NCR,KR,i,k)
                     IF(PSIMIN.GT.PSI(NCR,KR,i,k))PSIMIN=PSI(NCR,KR,i,k)
                  enddo 
                  !)THEN  .OR.
                  IF(PSIMAX*PSIMIN.LT.0.0 .OR.
     1               (PSIMAX.LT.0.D0 .AND. PSIMIN.LT.0.D0)) THEN
                     NERING(NCR,NCS,KR)=NERING(NCR,NCS,KR)+1
                     MIE(NCR,NCS,KR,NERING(NCR,NCS,KR))=i
                  ELSE
C                     write(6,*) KR,NERING(NCR,NCS,KR),
C     1	              MIE(NCR,NCS,KR,NERING(NCR,NCS,KR))
C                     pause
                  ENDIF
               enddo              


C              DODELJIVANJE VREDNOSTI tezinske funkcije QST ELEMENTIMA PRSTENA
C              DEF.LOKALNIH KOORDINATA NODOVA ELEMENATA PRESECENIH SA PRSTENOM 
	       DO KE=1,NERING(NCR,NCS,KR)
	          JG1=MIE(NCR,NCS,KR,KE)

	          do k=1,NNOD(JG1)
	             IF(PSI(NCR,KR,JG1,k).LT.0.000D0) THEN
	                QST(NCR,NCS,KR,KE,k)=1.D0
C                   write(6,*)'KR,JG1,k,PSI1',KR,JG1,k,PSI1(NCR,KR,JG1,k)
C                        write(6,*)'KR,i,MIE',KR,i,MIE(NCR,NCS,KR,i)
C                       write(6,*)'KR,i,k,QST',KR,i,k,QST(NCR,NCS,KR,i,k)
C                        PAUSE
                     ELSE
                        QST(NCR,NCS,KR,KE,k)=0.D0
C                   write(6,*)'KR,JG1,k,PSI1',KR,JG1,k,PSI1(NCR,KR,JG1,k)
C                        write(6,*)'KR,i,MIE',KR,i,MIE(NCR,NCS,KR,i)
C                       write(6,*)'KR,i,k,QST',KR,i,k,QST(NCR,NCS,KR,i,k)
C                        PAUSE
                     ENDIF
                  enddo

                  do K=1,NNOD(JG1)
c                    GLOBALNE KOORDINATE
                     XGG(NCR,NCS,KR,KE,K)=CORD(NUM(JG1,K),1)
                     YGG(NCR,NCS,KR,KE,K)=CORD(NUM(JG1,K),2)

C                    LOKALNE KOORDINATE
                     XL(NCR,NCS,KR,KE,K)=
     1                                VX1(NCR,1)*(CORD(NUM(JG1,K),1)-Xc)
     1                               +VX1(NCR,2)*(CORD(NUM(JG1,K),2)-Yc)
                     YL(NCR,NCS,KR,KE,K)=
     1                          -1.D0*VX1(NCR,2)*(CORD(NUM(JG1,K),1)-Xc)
     1                               +VX1(NCR,1)*(CORD(NUM(JG1,K),2)-Yc)
C      write(6,*)'JG1,NCR,KR,KE,K,XL',JG1,NCR,KR,KE,K,XL(NCR,NCS,KR,KE,K)
C      write(6,*)'JG1,NCR,KR,KE,K,YL',JG1,NCR,KR,KE,K,YL(NCR,NCS,KR,KE,K)             
C                     PAUSE
                  enddo

               ENDDO


            ENDDO !ZAVRSENA PETLJA PO PRSTENOVIMA (KR=1-IRING)
         ENDDO !ZAVRSENA PETLJA PO SEGMENTIMA (NCS=1-NSEG)
      ENDDO !ZAVRSENA PETLJA PO PRSLINAMA (NCR=1-NCRACK)


C     PRAVLJENJE DATOTEKE ZA CRTANJE
C     INICIJALIZACIJA POCETNIH VREDNOSTI ZA Q
      DO KR=1,IRING
         DO J=1,NE
            DO M=1,NNOD(J)
               Q(NUM(J,M),KR)=101.D0
            ENDDO
         ENDDO
      ENDDO
      !NCS=1 U SLEDECIM REDOVIMA ZA 2D SLUCAJ 
      !U 3D SLUCAJU Q MORA DA IMA DIMENZIJU I SEGMENTA!!!!!!
      DO KR=1,IRING
         DO NCR=1,NCRACK
            DO J=1,NERING(NCR,1,KR)
               JG1=MIE(NCR,1,KR,J)
               DO M=1,NNOD(JG1)
                  Q(NUM(JG1,M),KR)=QST(NCR,1,KR,J,M)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO KR=1,IRING
         DO NCR=1,NCRACK
            DO J=1,NE
               DO M=1,NNOD(J)
                  IF(Q(NUM(J,M),KR).GT.100.D0) Q(NUM(J,M),KR)=0.D0
               ENDDO		 	   
            ENDDO
         ENDDO
      ENDDO

      return
      end
C============================================================
