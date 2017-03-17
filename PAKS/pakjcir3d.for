C=======================================================================
C
C=======================================================================
      SUBROUTINE S_FUN_3D(CORD,NUM,PSI,NERING,MIE,QST,Q,
     1                    NP,NE,XL,YL,ZL,XGG,YGG,ZGG,NNOD,MODEL,DEB)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ******************************************************************
C     DEFINISANJE ELEMENATA PO KOJIMA SE VRSI INTEGRACIJA J-INTEGRALA I
C     RACUNANJE TEZINSKE FUNK.      
C     ******************************************************************
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /SRPSKI/ ISRPS

      DIMENSION Xpf(3),Xvp(3),VX1(10,3),VX2(10,3),VX3(10,3),PP(4)
      DIMENSION Tl(3),Tf(3),Bf(3)
      DIMENSION CORD(NP,3),NUM(NE,*),
     1          PSI(NCRACK,MAXRIN,NE,MAXNOD),
     1          NERING(NCRACK,MAXSEG,MAXRIN),
     1          MIE(NCRACK,MAXSEG,MAXRIN,N100),
     1          QST(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          Q(NP,MAXRIN),
     1          XL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          ZL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          XGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          ZGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          NNOD(NE)

      DIMENSION MFRONT(1000),		  !DIMENZIJE BROJA ELEMENATA DUZ FRONTA MAXSEG
     *          DEB(1000),		  !DIMENZIJE BROJA ELEMENATA DUZ FRONTA MAXSEG
     *          RAST(1000),
     *          R(1000),
     *          PSI1(10000,21),	  !DIMENZIJE UKUPNO BROJA ELEMENATA
     *          FI1(10000,21),	  !DIMENZIJE UKUPNO BROJA ELEMENATA
     *          PSI2(10000,21),	  !DIMENZIJE UKUPNO BROJA ELEMENATA
     *          FI2(10000,21)       !DIMENZIJE BROJA ELEMENAATA U PRSTENU N100
C	data PP /0.8D0,1.D0,1.1D0,1.2D0/
c	data PP /1.3D0,1.5D0,1.6D0,1.8D0/
C       data PP /2.D0,2.2D0,2.5D0,2.8D0/      
c      DATA PP /3.D0,4.D0,5.D0,6.D0/	     
c      DATA PP /7.D0,8.D0,9.D0,10.D0/
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
      CALL CLEAR(QST,NCRACK*MAXSEG*MAXRIN*N100*MAXNOD)



      DO NCR=1,NCRACK	      !PETLJA PO PRSLINAMA
         MODEL=NODCR(NCR,10)
         NSEG=NODCR(NCR,11)
         IRING=NODCR(NCR,12)

         !JEDINICNI VEKTORI PRAVCA-X1 (LICE PRSLINE) NCR
         VX1(NCR,1)=CORD(NODCR(NCR,3),1)-CORD(NODCR(NCR,4),1)
         VX1(NCR,2)=CORD(NODCR(NCR,3),2)-CORD(NODCR(NCR,4),2)
         VX1(NCR,3)=CORD(NODCR(NCR,3),3)-CORD(NODCR(NCR,4),3)
         VV=dsqrt(VX1(NCR,1)*VX1(NCR,1)+VX1(NCR,2)*VX1(NCR,2)+
     1            VX1(NCR,3)*VX1(NCR,3) )
         VX1(NCR,1)=VX1(NCR,1)/VV
         VX1(NCR,2)=VX1(NCR,2)/VV
         VX1(NCR,3)=VX1(NCR,3)/VV

         CONTE=VX1(NCR,1)
         SINTE=VX1(NCR,2)

         !DUZINA LICA PRSLINE
         DUZ=VV
         WRITE(3,*) ' '
         if (isrps.eq.0) then
            WRITE(6,*) ' PRSLINA',NCR,' DUZINA',DUZ
            WRITE(3,*) ' PRSLINA',NCR,' DUZINA',DUZ
         else
            WRITE(6,*) ' CRACK',NCR,' LENGTH',DUZ
            WRITE(3,*) ' CRACK',NCR,' LENGTH',DUZ
         endif
C         PAUSE
         !POLUPRECNICI PRSTENOVA PO KOJIMA SE VRSI INTEGRACIJA
c              pd1-pd4
C         PP(1)=0.1*DUZ   !NA 20% DUZINE PRSLINE
C         PP(2)=0.15*DUZ   !NA 30% DUZINE PRSLINE
C         PP(3)=0.20*DUZ   !NA 40% DUZINE PRSLINE 
C         PP(4)=0.25*DUZ   !NA 50% DUZINE PRSLINE
         PP(1)=DUZ*NODCR(NCR,6)/100   !NA 20% DUZINE PRSLINE
         PP(2)=DUZ*NODCR(NCR,7)/100   !NA 30% DUZINE PRSLINE
         PP(3)=DUZ*NODCR(NCR,8)/100   !NA 40% DUZINE PRSLINE 
         PP(4)=DUZ*NODCR(NCR,9)/100   !NA 50% DUZINE PRSLINE

         !JEDINICNI VEKTORI PRAVCA-X3 (FRONT PRSLINE) NCR
         VX3(NCR,1)=CORD(NODCR(NCR,3),1)-CORD(NODCR(NCR,5),1)
         VX3(NCR,2)=CORD(NODCR(NCR,3),2)-CORD(NODCR(NCR,5),2)
         VX3(NCR,3)=CORD(NODCR(NCR,3),3)-CORD(NODCR(NCR,5),3)
         VV=dsqrt(VX3(NCR,1)*VX3(NCR,1)+VX3(NCR,2)*VX3(NCR,2)+
     1            VX3(NCR,3)*VX3(NCR,3) )
         VX3(NCR,1)=VX3(NCR,1)/VV
         VX3(NCR,2)=VX3(NCR,2)/VV
         VX3(NCR,3)=VX3(NCR,3)/VV

         !JEDINICNI VEKTORI PRAVCA-X2 (BINORMALA PRSLINE) NCR
         VX2(NCR,1)=VX3(NCR,2)*VX1(NCR,3)-VX3(NCR,3)*VX1(NCR,2)
         VX2(NCR,2)=VX3(NCR,3)*VX1(NCR,1)-VX3(NCR,1)*VX1(NCR,3)
         VX2(NCR,3)=VX3(NCR,1)*VX1(NCR,2)-VX3(NCR,2)*VX1(NCR,1)

         WRITE(6,*)'VX1(1),VX1(2)',VX1(NCR,1),VX1(NCR,2),VX1(NCR,3)
         WRITE(6,*)'VX2(1),VX2(2)',VX2(NCR,1),VX2(NCR,2),VX2(NCR,3)
         WRITE(6,*)'VX3(1),VX3(2)',VX3(NCR,1),VX3(NCR,2),VX3(NCR,3)

C        KOORDINATE VRHA PRSLINE 
         Xvp(1)=CORD(NODCR(NCR,3),1)
         Xvp(2)=CORD(NODCR(NCR,3),2)
         Xvp(3)=CORD(NODCR(NCR,3),3) 


C        KOORDINATE PRAVCA FRONTA PRSLINE
         Xpf(1)=CORD(NODCR(NCR,5),1)
         Xpf(2)=CORD(NODCR(NCR,5),2)
         Xpf(3)=CORD(NODCR(NCR,5),3) 

         Tl(1)=VX1(NCR,1)
         Tl(2)=VX1(NCR,2)
         Tl(3)=VX1(NCR,3)

         Tf(1)=VX3(NCR,1)
         Tf(2)=VX3(NCR,2)
         Tf(3)=VX3(NCR,3)

         Bf(1)=VX2(NCR,1)
         Bf(2)=VX2(NCR,2)
         Bf(3)=VX2(NCR,3)

         CALL FRONT(MFRONT,DEB,RAST,R,NFRONT,PSI1,PSI2,FI1,
     *              CORD,NUM,Tl,Tf,Bf,Xpf,Xvp,NP,NE,NNOD)


         DO NCS=1,NFRONT	    !PETLJA PO SEGMENTIMA NSEG=NFRONT

            !CENTAR TEKUCEG SEGMENTA
            Xc=CORD(MFRONT(NCS),1)
            Yc=CORD(MFRONT(NCS),2)
            Zc=CORD(MFRONT(NCS),3)


            WRITE (6,*) 'NFRONT,MFRONT(NCS)',NFRONT,MFRONT(NCS)
            WRITE (6,*) 'DEB(NCS),R(NCS)',   DEB(NCS),R(NCS)
            WRITE (6,*) 'CENTAR TEKUCEG SEGMENTA'
            WRITE (6,*) Xc,Yc,Zc
            !debljina elementa
            THI=DEB(NCS)	 ! dez=1.
	   
            DO KR=1,IRING     !PETLJA PO PRSTENOVIMA
	     
C               CALL IWRR(NNOD,NE,'NCV2') 

C               R1=PP(KR)*R(NCS)    !Tekuci poluprecnik int. putanje
               R1=PP(KR)
               if(isrps.eq.0.AND.NCS.EQ.1) then
                  WRITE(6,*) ' PRSTEN',KR,' POLUPRECNIK',R1
                  WRITE(3,*) ' PRSTEN',KR,' POLUPRECNIK',R1
               else
                  WRITE(6,*) ' RING',KR,' RADIUS',R1
                  WRITE(3,*) ' RING',KR,' RADIUS',R1
               endif
	     
               do i=1,NE
                  do k=1,20     !NNOD(i)
C                     PSI(NCR,KR,i,k)=
C     1            dsqrt(DABS((CORD(NUM(i,k),1)-Xc)*(CORD(NUM(i,k),1)-Xc)
C     1                 +(CORD(NUM(i,k),2)-Yc)*(CORD(NUM(i,k),2)-Yc)
C     1                 +(CORD(NUM(i,k),3)-Zc)*(CORD(NUM(i,k),3)-Zc)))-R1
                   PSI(NCR,KR,i,k)=DSQRT(PSI1(i,k)**2.+PSI2(i,k)**2.)-R1
                  enddo
               enddo

	
C              Identifikacija elemenata koji su preseceni sa int.putanjom
               do i=1,NE
                  PSIMIN=PSI(NCR,KR,i,1)
                  PSIMAX=PSI(NCR,KR,i,1)

                  do k=1,NNOD(i)     !NNOD(i) 20
                     IF(PSIMAX.LT.PSI(NCR,KR,i,k))PSIMAX=PSI(NCR,KR,i,k)
                     IF(PSIMIN.GT.PSI(NCR,KR,i,k))PSIMIN=PSI(NCR,KR,i,k)
                  enddo 
                  !USLOV DA JE ELEMENT PRESECEN SA KRUZNICOM ILI SE NALAZI UNUTAR NJE
                  IF(PSIMAX*PSIMIN.LT.0.D0 .OR.
     1               (PSIMAX.LT.0.D0 .AND. PSIMIN.LT.0.D0)) THEN
                     NERING(NCR,NCS,KR)=NERING(NCR,NCS,KR)+1
                     MIE(NCR,NCS,KR,NERING(NCR,NCS,KR))=i
                  ENDIF

               enddo              


C              DODELJIVANJE VREDNOSTI tezinske funkcije QST ELEMENTIMA PRSTENA
               DO KE=1,NERING(NCR,NCS,KR)
                  JG1=MIE(NCR,NCS,KR,KE)
                  do m=1,NNOD(JG1)
                     FI2(KE,m)=(CORD(NUM(JG1,m),1)-Xc)*Tf(1)+
     *                         (CORD(NUM(JG1,m),2)-Yc)*Tf(2)+
     *                         (CORD(NUM(JG1,m),3)-Zc)*Tf(3)
                  enddo
	 
                  do k=1,NNOD(JG1)    
C                     write(6,*)'FI2(JG1,k)',FI2(JG1,k)
                     !ULSOV DA SE CVOR NALAZI UNUTAR CILINDRA
                     IF(PSI(NCR,KR,JG1,k).LT.0.0D0 .AND.
     1                  DABS(FI2(KE,k)).LT.0.5001D0*THI) THEN
C     1	         ((FI2(JG1,k).LT.0.5001D0*THI.AND.FI2(JG1,k).GT.0.D0).OR.
C     1        (FI2(JG1,k).GT.-0.5001D0*THI.AND.FI2(JG1,k).LT.0.D0)))THEN 
C     1             ABS(ABS(CORD(NUM(JG1,k),3))-ABS(Zc)).LT.0.5*dez) THEN

                        QST(NCR,NCS,KR,KE,k)=1.D0

C           write(6,*)'KR,JG1,k,PSI', KR,JG1,NUM(JG1,k),PSI(NCR,KR,JG1,k)
C                        write(6,*)'KR,KE,MIE',KR,KE,MIE(NCR,NCS,KR,KE)
C         write(6,*)'KR,KE,k,QST',  KR,KE,NUM(JG1,k),QST(NCR,NCS,KR,KE,k)
C                        PAUSE
                     ELSE
                        QST(NCR,NCS,KR,KE,k)=0.D0
C                    write(6,*)'KR,JG1,k,PSI', KR,JG1,k,PSI(NCR,KR,JG1,k)
C                        write(6,*)'KR,KE,MIE', KR,KE,MIE(NCR,NCS,KR,KE)
C                  write(6,*)'KR,KE,k,QST',  KR,KE,k,QST(NCR,NCS,KR,KE,k)
C                        PAUSE
                     ENDIF
                  enddo
C                  PAUSE

                  do K=1,NNOD(JG1)     !NNOD(JG1)  20

c                    GLOBALNE KOORDINATE
                     XGG(NCR,NCS,KR,KE,K)=CORD(NUM(JG1,K),1)
                     YGG(NCR,NCS,KR,KE,K)=CORD(NUM(JG1,K),2)
                     ZGG(NCR,NCS,KR,KE,K)=CORD(NUM(JG1,K),3)

C                    LOKALNE KOORDINATE
                     XL(NCR,NCS,KR,KE,K)=
     1                                VX1(NCR,1)*(CORD(NUM(JG1,K),1)-Xc)
     1                               +VX1(NCR,2)*(CORD(NUM(JG1,K),2)-Yc)
     1                               +VX1(NCR,3)*(CORD(NUM(JG1,K),3)-Zc)

                     YL(NCR,NCS,KR,KE,K)=
     1                                VX2(NCR,1)*(CORD(NUM(JG1,K),1)-Xc)
     1                               +VX2(NCR,2)*(CORD(NUM(JG1,K),2)-Yc)
     1                               +VX2(NCR,3)*(CORD(NUM(JG1,K),3)-Zc)

                     ZL(NCR,NCS,KR,KE,K)=
     1                                VX3(NCR,1)*(CORD(NUM(JG1,K),1)-Xc)
     1                               +VX3(NCR,2)*(CORD(NUM(JG1,K),2)-Yc)
     1                               +VX3(NCR,3)*(CORD(NUM(JG1,K),3)-Zc)
	
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
         DO M=1,NP
            Q(M,KR)=0.D0
         ENDDO
      ENDDO

      !NCS=1 U SLEDECIM REDOVIMA ZA 2D SLUCAJ 
      !U 3D SLUCAJU Q MORA DA IMA DIMENZIJU I SEGMENTA!!!!!!
      DO KR=1,IRING	            !PETLJA PO PRSTENOVIMA
         DO NCR=1,NCRACK		    !PETLJA PO PRSLINAMA
            DO J=1,NERING(NCR,1,KR)	!PETLJA PO ELEMENTIMA
               JG1=MIE(NCR,1,KR,J)
               DO M=1,NNOD(JG1)      !NNOD(JG1)  20
C                  IF(QST(NCR,1,KR,J,M).LT.0.D0)THEN
C                     QST(NCR,1,KR,J,M)=0.D0
C                  ENDIF
                  Q(NUM(JG1,M),KR)=QST(1,1,KR,J,M)
	
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      WRITE(6,*) 'STAMPANJE UPOZORENJA'
      DO KR=1,IRING
         DO I=I,NP
            IF(Q(I,KR).LT.0.D0) THEN
               WRITE(6,*) I,KR,Q(I,KR)
               PAUSE
            ENDIF
         ENDDO
      ENDDO

      return
      end
C============================================================
