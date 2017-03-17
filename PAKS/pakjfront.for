
      SUBROUTINE FRONT(MFRONT,DEB,RAST,R,NFRONT,PSI1,PSI2,FI1,
     *           CORD,NUM,TL,TF,BF,Xpf,Xvp,NP,NE,NNOD)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION TF(3),TL(3),BF(3)
      DIMENSION Xpf(3),Xvp(3)
      DIMENSION CORD(NP,3),
     *          PSI1(22000,21),
     *          FI1(22000,21),
     *          PSI2(22000,21),
     *          NUM(NE,*),NNOD(NE),
     *          MFRONT(1000),
     *          DEB(1000),
     *          RAST(1000),
     *          R(1000)

C     MFRONT-SREDISNJI NODOVI DUZ FRONTA ZA INTEGRACIJU
C     DEB-DEBLJINA ELEMENTA DUZ FRONTA
C     RAST-RASTOJANJE SREDISNJEG NODA OD VRHA PRSLINE  
C     T1-DUZINA FRONTA
c     R-POLUPRECNIK KRUZNICE ZA INTEGRACIJU

C     DUZINA FRONTA
      T1=dsqrt((Xpf(1)-Xvp(1))*(Xpf(1)-Xvp(1))+
     1         (Xpf(2)-Xvp(2))*(Xpf(2)-Xvp(2))+
     1         (Xpf(3)-Xvp(3))*(Xpf(3)-Xvp(3)))

c     TANGENTA FRONTA (TF)-X3
      tx1=TF(1)
      ty1=TF(2)
      tz1=TF(3)

C     NORMALA NA FRONT (TANGENTA-LICA-TL)-X1
      bx1=TL(1)
      by1=TL(2)
      bz1=TL(3)

C     BINORMALA FRONTA (BF)-X2
      bbx1=BF(1)
      bby1=BF(2)
      bbz1=BF(3)

C      WRITE(6,*) 'duzina fronta T1',T1
C      WRITE(6,*) 'fornt tx1,ty1,tz1',tx1,ty1,tz1
C      WRITE(6,*) 'normala bx1,by1,bz1',bx1,by1,bz1
C      pause
C      DEFINISANJE FI i PSI FUNKCIJA
      do i=1,NE
         do m=1,NNOD(i)
            FI1(i,m)=(CORD(NUM(i,m),1)-Xpf(1))*TF(1)+
     *               (CORD(NUM(i,m),2)-Xpf(2))*TF(2)+
     *               (CORD(NUM(i,m),3)-Xpf(3))*TF(3)

            PSI1(i,m)=(CORD(NUM(i,m),1)-Xpf(1))*TL(1)+
     *                (CORD(NUM(i,m),2)-Xpf(2))*TL(2)+
     *                (CORD(NUM(i,m),3)-Xpf(3))*TL(3)

            PSI2(i,m)=(CORD(NUM(i,m),1)-Xpf(1))*BF(1)+
     *                (CORD(NUM(i,m),2)-Xpf(2))*BF(2)+
     *                (CORD(NUM(i,m),3)-Xpf(3))*BF(3)

         enddo
      enddo
      WRITE(6,*) 'FRONT'
C     DEFINISANJE MEDJUCVOROVA ZA INTEGRACIJU PO SEGMENTIMA
      NBRO=0
      DO I=1,NE
         DO M=9,NNOD(I)

C	   IF( FI1(I,M).LT.0.D0        .AND.
C     *	  -0.000001D0.GT.PSI1(I,M) .AND.
C     *       PSI1(I,M).LT.0.000001D0) THEN

            IF( PSI1(I,M).GT.-0.0000001D0 .AND.
     *          PSI1(I,M).LT.0.0000001D0  .AND.
     *          PSI2(I,M).GT.-0.0000001D0 .AND.
     *          PSI2(I,M).LT.0.0000001D0  .AND.
     *          FI1(I,M).GT.-0.0000001D0 ) THEN

               WRITE(6,*) 'I,M,PSI1(I,M),NUM(I,M)'
               WRITE(6,*) I,M,PSI1(I,M),NUM(I,M)
               PAUSE
               NF1=0
               IF(NBRO.GT.0) THEN
                  DO K=1,NBRO
                     IF(NUM(I,M).EQ.MFRONT(K)) THEN
                        NF1=NF1+1
                     ENDIF
                  ENDDO
               ENDIF

               IF(NF1.EQ.0) THEN
                  NBRO=NBRO+1
                  MFRONT(NBRO)=NUM(I,M)
                  II=I
                  MM=M
                  CALL DEBLJINA(DEB,R,FI1,PSI1,PSI2,NBRO,II,MM)
                  RAST(NBRO)=T1-DABS(FI1(II,MM))
                  WRITE(6,*)'NBRO,DEB(NBRO),R(NBRO),RAST(NBRO)'
                  WRITE(6,*) NBRO,DEB(NBRO),R(NBRO),RAST(NBRO)
                  PAUSE
               ENDIF

            ENDIF

         ENDDO
      ENDDO

      NFRONT=NBRO

      RETURN
      END
C     ******************************************************

C     ******************************************************

      SUBROUTINE DEBLJINA(DEB,R,FI1,PSI1,PSI2,NBRO,II,MM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION FI1(22000,21),PSI1(22000,21),PSI2(22000,21),
     *          DEB(1000),R(1000)

      IF(MM.EQ.9) THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,1))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,11)**2.+PSI2(II,11)**2.),
     1               DSQRT(PSI1(II,13)**2.+PSI2(II,13)**2.))

      ELSE IF(MM.EQ.10)THEN
          DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,2))) 
          R(NBRO)=MAX(DSQRT(PSI1(II,12)**2.+PSI2(II,12)**2.),
     1                DSQRT(PSI1(II,14)**2.+PSI2(II,14)**2.))
	
      ELSE IF(MM.EQ.11)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,3))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,9)**2.+PSI2(II,9)**2.),
     1               DSQRT(PSI1(II,15)**2.+PSI2(II,15)**2.))
	
      ELSE IF(MM.EQ.12)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,4))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,10)**2.+PSI2(II,10)**2.),
     1               DSQRT(PSI1(II,16)**2.+PSI2(II,16)**2.))
	
      ELSE IF(MM.EQ.13)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,5))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,9)**2.+PSI2(II,9)**2.),
     1               DSQRT(PSI1(II,15)**2.+PSI2(II,15)**2.))
	
      ELSE IF(MM.EQ.14)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,6))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,10)**2.+PSI2(II,10)**2.),
     1               DSQRT(PSI1(II,16)**2.+PSI2(II,16)**2.))
	
      ELSE IF(MM.EQ.15)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,7))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,11)**2.+PSI2(II,11)**2.),
     1               DSQRT(PSI1(II,13)**2.+PSI2(II,13)**2.))
	
      ELSE IF(MM.EQ.16)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,8))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,12)**2.+PSI2(II,12)**2.),
     1               DSQRT(PSI1(II,14)**2.+PSI2(II,14)**2.))
	
      ELSE IF(MM.EQ.17)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,1))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,18)**2.+PSI2(II,18)**2.),
     1	             DSQRT(PSI1(II,20)**2.+PSI2(II,20)**2.))
C         write(6,*) 'PSI1(II,18),PSI2(II,18)'
C         write(6,*)  PSI1(II,18),PSI2(II,18)
C         write(6,*) 'PSI1(II,20),PSI2(II,20)'
C         write(6,*)  PSI1(II,20),PSI2(II,20)
C         pause
      ELSE IF(MM.EQ.18)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,2))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,17)**2.+PSI2(II,17)**2.),
     1               DSQRT(PSI1(II,19)**2.+PSI2(II,19)**2.))
C         write(6,*) 'PSI1(II,17),PSI2(II,17)'
C         write(6,*)  PSI1(II,17),PSI2(II,17)
C         write(6,*) 'PSI1(II,19),PSI2(II,19)'
C         write(6,*)  PSI1(II,19),PSI2(II,19)
C         pause
	
      ELSE IF(MM.EQ.19)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,3))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,18)**2.+PSI2(II,18)**2.),
     1	             DSQRT(PSI1(II,20)**2.+PSI2(II,20)**2.))
	
      ELSE IF(MM.EQ.20)THEN
         DEB(NBRO)=2.D0*DABS(DABS(FI1(II,MM))-DABS(FI1(II,4))) 
         R(NBRO)=MAX(DSQRT(PSI1(II,17)**2.+PSI2(II,17)**2.),
     1               DSQRT(PSI1(II,19)**2.+PSI2(II,19)**2.))
      
      ENDIF

	
      RETURN
      END
