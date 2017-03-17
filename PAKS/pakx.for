C
C    *****************************************************************
C    *                                                               *
C    *   IDENTIFIKACIJA DISKONTINUITETA PRIMENOM LS-FUNKCIJA         *
C    *                                                               *
C    *   Developed by G.Jovicic & M.Zivkovic                         *
C    *                                                               *
C    *                                                               *
C    *   University of Kragujevac                                    *
C    *   Faculty of Mechanical Engineering                           *
C    *   Department of Applaed Mechanic                              *
C    *   Kragujevac, YUGOSLAVIA, September 2002 A.D.                 *
C    *                           November  2003 A.D.                 *
C    *                                                               *
C    *****************************************************************
C
      SUBROUTINE LS_EDGE(NUM,FI,PSI,KELEM,X,Y,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
      
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /VRHP/ VRHX,VRHY
      DIMENSION NUM(NE,*),X(*),Y(*),FI(NE,*),PSI(NE,*),KELEM(*)
      DIMENSION DUZ(100)
C      write(6,*)'START LS_EDGE.FOR'
C     ******************************************************************
C     Ovaj potprogram sluzi za definisanje LS-funkcija ivicne prsline (Edge Crack)
C     ******************************************************************
C     ==================================================================
C     Ulazne velicine:

C             NSEGX      - broj segmenata prsline
C             Xp(NSEGX+1)- x-koordinate tacaka segmenata
C             Yp(NSEGX+1)- y-koordinate tacaka segmenata
C    ===================================================================
C     Potrebne velicine  iz potprograma SETUP.FOR su:
C             NB        - broj cvorova mreze  
C             X(NB)     - x-koordinate cvorova mreze
C             Y(NB)     - y-koordinate cvorova mreze                        
C             NE        - broj elemenata Narrow-Band-a
C             NUM(NE,4) - globalna numeracija cvorova

C    ===================================================================   
C     Pomocne velicine LS_EDGE.FOR potprograma su:
C             FIend,FIstart- FI funkcija pocetka i kraja segmenta
C             tx1,ty1,t1   - komponente tangente i njen intenzitet za prvi segment
C             bx1,by1      - komponente normale prvog ili narednog segmenta
C             tx,ty,t      - komponente tangente i njen intenzitet tekuci segment
C             bx,by        - komponente normale za tekuci segment
C             FIMAX,FIMIN  - maks.i min. vrednost funkcije FI u tekucem elementu  
C             PSIMAX,PSIMIN- maks.i min. vrednost funkcije PSI u tekucem elementu 
C    ===================================================================
C     Izlazne velicine potprograma LS_EDGE.FOR su:
C             FI(NE,4)  - LS-funkcija kojom se identifikuje vrh prsline
C             PSI(NE,4) - LS-funkcija kojom se identifikuje polozaj prsline (PSI=0)
C             NEC       - Broj elemenata presecenih sa prslinom
C             KELEM(NEC)- Niz koji pamti elemente koji su preseceni sa prslinom   
C     ==================================================================
C     Upozorenje!
C             1. Nizovi Xp i Yp su postavljeni u pakxfem.inc i dimenzionisani 
C                na maksimalnu vrednost NSEGX+1=101, DUZ(NSEGX)      
C             2. U potprogram se unose vrednosti velicina NSEGX,Xp(NSEGX+1),Yp(NSEGX+1)
C     ==================================================================
C     NAPOMENA!!!!!!!!
c             LS-FUNKCIJA FI SE RACUNA OD VRHA PRIHVACENOG SEGMENTA
C     ==================================================================
C     ==================================================================
C     BLOK ZA UNOSENJE ULAZNIH PODATAKA  GEOMETRIJE PRSLINE 
c      write(6,10)
c   10 FORMAT(/' UNESITE BROJ SEGMENATA PRSLINE')
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) NSEGX,NNB,MNB
      IF(INDFOR.EQ.2)
     1READ(ACOZ,20) NSEGX,NNB,MNB
   20 FORMAT(3i5)    
      write(3,90) NSEGX,NNB,MNB
   90 format(/' NSEGX=',I5,' NNB=',I5,' MNB=',I5)      
      if(nsegX.gt.100) stop 'nsegX > 100'
C      write(6,30)
      do i=1,NSEGX+1
c       write(6,40)i
c       read(5,50)Xp(i)
c       write(6,60)i
c       read(5,70)Yp(i) 
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1read(1,*)Xp(i),Yp(i) 
      IF(INDFOR.EQ.2)
     1read(ACOZ,70)Xp(i),Yp(i) 
      enddo
CV    PRENOS VRHA PRSLINE	
	VRHX=Xp(NSEGX+1)
	VRHY=Yp(NSEGX+1)

c     crtanje linije
C
      II=49
      MJ=-1
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
C
      WRITE(II,1000) MJ
      M470=470
      WRITE(II,1001) M470
      do i=1,nsegx+1
      write(ii,2000) i,0,0,1,24,0,0,0,0,0,0,xp(i),yp(i),0.
 2000 format(i5,',',10(i3,','),3(1pe13.5,','))     
      enddo       
 1000 FORMAT(10I5)
 1001 FORMAT(I6,F12.4)
      WRITE(II,1000) MJ
      WRITE(II,1000) MJ
      M471=471
      WRITE(II,1001) M471
      do i=1,nsegx
      write(ii,2002) i,4,0,1,0,0.
      write(ii,2001) i,i+1,0,0,0,0,0,0,0,0
      write(ii,2001) 0,0,0,0,0,0,0,0,0,0
      write(ii,2001) 0,0,0,0,0,0,0,0,0,0
 2001 format(10(i5,','))     
 2002 format(5(i5,','),1pe13.5,',')     
      enddo       
      WRITE(II,1000) MJ
C
      
   30 format(/' UNESITE KOORDINATE i-te TACKE SEGMENTA')   
   40 format(' Xp(',I3,')=')   
   50 format(F15.10)          
   60 format(' Yp(',I3,')=')   
   70 format(2F10.4)        
      do i=1,NSEGX+1
         write(3,80)i,Xp(i),i,Yp(i)
      enddo
   80 format(/' Xp(',I3,')=',F15.10,' ,Yp(',I3,')=',F15.10)      
C     ==================================================================   
C     ==================================================================
c     10. 10. 2003 DUZINA SVAKOG SEGMENTA
      DO I=1,NSEGX
       DUZ(I)=DSQRT((Xp(I+1)-Xp(I))**2+(Yp(I+1)-Yp(I))**2) 
      ENDDO

C     Tangenta i normala prvog segmenta (komponente normale su oznacene sa bx1,by1) 
      t1=Dsqrt((Xp(2)-Xp(1))**2+(Yp(2)-Yp(1))**2)
      tx1=(Xp(2)-Xp(1))/t1
      ty1=(Yp(2)-Yp(1))/t1
      bx1=-ty1
      by1=tx1
C    ===================================================================
C     DEFINISANJE LS-funkcija: FI i PSI po cvorovima elemenata
C     BLOK KOJI SE IZVRSAVA UKOLIKO IMAMO VISE SEGMENATA      
      IF(NSEGX.GT.1) THEN
       DO k=1,NSEGX-1
C       RASTOJANJE DO VRHA POSLEDNJEG SEGMENTA
        RAST=0.0
        DO ICE=NSEGX,k+1,-1
         RAST=RAST-DUZ(ICE)
        ENDDO
C       Tangenta i normala tekuceg segmenta
        t=Dsqrt((Xp(k+1)-Xp(k))**2+(Yp(k+1)-Yp(k))**2)
        tx=(Xp(k+1)-Xp(k))/t
        ty=(Yp(k+1)-Yp(k))/t
        bx=-ty
        by=tx
C       Tangenta i normala narednog segmenta         
        t1=sqrt((Xp(k+2)-Xp(k+1))**2+(Yp(k+2)-Yp(k+1))**2)
        tx1=(Xp(k+2)-Xp(k+1))/t1
        ty1=(Yp(k+2)-Yp(k+1))/t1
        bx1=-ty1
        by1=tx1

c        WRITE(6,*)'SEGMENT'
c        WRITE(6,*)k
c        pause     
C       Petlja po elementima         
        do i=1,NE
         do m=1,4
        
          FIstart=(X(NUM(i,m))-Xp(k))*tx+(Y(NUM(i,m))-Yp(k))*ty 
          FIend=(X(NUM(i,m))-Xp(k+1))*tx+(Y(NUM(i,m))-Yp(k+1))*ty
          PSIend=(X(NUM(i,m))-Xp(k+1))*bx+(Y(NUM(i,m))-Yp(k+1))*by
                                        
          FIstart1=(X(NUM(i,m))-Xp(k+1))*tx1+(Y(NUM(i,m))-Yp(k+1))*ty1
          FIend1=(X(NUM(i,m))-Xp(k+2))*tx1+(Y(NUM(i,m))-Yp(k+2))*ty1                               
C         ==============================================================
C         PRVI SEGMENT            
          if(k.EQ.1)then
C          Uslov ako se cvor nalazi u zoni prvog segmenta           
c           if(FIstart.GE.0.0.AND.FIend.LE.0.0) then
           if(FIend.LE.0.0) then
            FI(i,m)=FIend+RAST
            PSI(i,m)=PSIend
           endif
C          Uslov ako je tup ugao izmedju prvog i drugog segmenta
           if(FIstart.GT.0.0.AND.FIend.GE.0.0.AND.FIstart1.LE.0.0) then
            FI(i,m)=RAST !!!!!!FIend1+
            PSI(i,m)=Dsign(Dsqrt((X(NUM(i,m))-Xp(k+1))**2+
     *                           (Y(NUM(i,m))-Yp(k+1))**2),PSIend)
           endif                              
c          =============================================================
          endif     
C         ==============================================================
C         UNUTRASNJI SEGMENTI       
          if(k.GT.1) then                   
           
C          Uslov ako se cvor nalazi na startnoj liniji tekuceg segmenta           
C           if(FIstart.EQ.0.0 .AND.Dabs(PSI(i,m)).GE.Dabs(PSIend))then
           if(Dabs(FIstart).lt.1.e-9
     *                          .AND.Dabs(PSI(i,m)).GE.Dabs(PSIend))then
            FI(i,m)=FIend+RAST
            PSI(i,m)=PSIend
           endif   
                                                                   
C          Uslov ako se cvor poklapa sa pocetnom tackom tekuceg i krajnjom tackom predhodnog segmenta                                                                    
           if(Dabs(FIstart).lt.1.e-9.AND.Dabs(FI(i,m)).lt.1.e-9
     *                              .AND.Dabs(PSI(i,m)).lt.1.e-9) then
             FI(i,m)=FIend+RAST
           endif

C          Uslov ako se cvor nalazi unutar zone tekuceg segmenta                                
           if(FIstart.GT.0.0.AND.FIend.LE.0.0) then
            FI(i,m)=FIend+RAST
            PSI(i,m)=PSIend
           endif 
                     
C          Uslov ako je tup ugao izmedju segmenata (ako se cvor nalazi u zoni tupog ugla)                    
           if(FIstart.GT.0.0.AND.FIend.GE.0.0.AND.FIstart1.LE.0.0) then
            FI(i,m)=RAST !!!!+RAST !!!!!FIend1 !!!!!!FIend+(druga var)
            PSI(i,m)=sign(sqrt((X(NUM(i,m))-Xp(k+1))**2+
     *                       (Y(NUM(i,m))-Yp(k+1))**2),PSIend)
           endif                              

c          =============================================================

          endif   
C         ==============================================================            
         enddo
        enddo
         
       ENDDO !PETLJA OD PRVOG DO PREDPOSLEDNJEG SEGMENTA
C      =================================================================
C      POSLEDNJI SEGMENT
C      Tangenta i normala poslednjeg segmenta
       t=dsqrt((Xp(NSEGX+1)-Xp(NSEGX))**2+(Yp(NSEGX+1)-Yp(NSEGX))**2)
       tx=(Xp(NSEGX+1)-Xp(NSEGX))/t
       ty=(Yp(NSEGX+1)-Yp(NSEGX))/t
       bx=-ty
       by=tx
C      Tangenta i normala predposlednjeg segmenta
       t1=dsqrt((Xp(NSEGX)-Xp(NSEGX-1))**2+(Yp(NSEGX)-Yp(NSEGX-1))**2)
       tx1=(Xp(NSEGX)-Xp(NSEGX-1))/t1
       ty1=(Yp(NSEGX)-Yp(NSEGX-1))/t1
       bx1=-ty1
       by1=tx1
       
C      Petlja po elementima         
       do i=1,NE
        do m=1,4
        FIstart=(X(NUM(i,m))-Xp(NSEGX))*tx+(Y(NUM(i,m))-Yp(NSEGX))*ty 
        FIend=(X(NUM(i,m))-Xp(NSEGX+1))*tx+(Y(NUM(i,m))-Yp(NSEGX+1))*ty
        PSIend=(X(NUM(i,m))-Xp(NSEGX+1))*bx+(Y(NUM(i,m))-Yp(NSEGX+1))*by
        FIend1=(X(NUM(i,m))-Xp(NSEGX))*tx1+(Y(NUM(i,m))-Yp(NSEGX))*ty1
C        Uslov ako se cvor nalazi na startnoj liniji poslednjeg segmenta
         if(Dabs(FIstart).lt.1.e-9
     *                        .AND.dabs(PSI(i,m)).GE.dabs(PSIend))then
          FI(i,m)=FIend
          PSI(i,m)=PSIend
         endif   
          
C        Uslov ako se cvor nalazi unutar zone poslednjeg segmenta         
         if(FIstart.GT.0.0) then
          FI(i,m)=FIend
          PSI(i,m)=PSIend
         endif                                                       
          
c        ===============================================================
c         Uslov ako dolazi do preklapanja
c          if(FIstart.GT.0.0.AND.FIend1.EQ.0.0 ) then
c            FI(i,m)=FIend                                           
c            PSI(i,m)=sign(sqrt((X(NUM(i,m))-Xp(NSEGX))**2+
c      *                       (Y(NUM(i,m))-Yp(NSEGX))**2),PSIend)     
c          endif
c        ===============================================================
        enddo
       enddo
C      =================================================================           
      ENDIF   !PETLJA AKO IMA VISE SEGMENATA OD JEDNOG
C     Zavrsen uslov ako imamo vise segmenata
C     ==================================================================
C     ==================================================================
C     BLOK KOJI SE IZVRSAVA UKOLIKO POSTOJI SAMO JEDAN SEGMENT
      IF(NSEGX.EQ.1) THEN
      
       do i=1,NE
        do m=1,4
         FI(i,m)=(X(NUM(i,m))-Xp(2))*tx1+(Y(NUM(i,m))-Yp(2))*ty1
         PSI(i,m)=(X(NUM(i,m))-Xp(2))*bx1+(Y(NUM(i,m))-Yp(2))*by1
        enddo
       enddo
        
      ENDIF
c      do i=1,NE
c       do m=1,4
c        write(3,*) 'i,m,num(i,m)',i,m,num(i,m)
c        write(3,*) 'fi,psi',fi(i,m),psi(i,m)
c       enddo 
c      enddo 
C     ============================================================================
C     Identifikacija elemenata koji su preseceni sa prslinom

      NEC=0
      do i=1,NE
        PSIMIN=PSI(i,1)
        PSIMAX=PSI(i,1)
        FIMIN=FI(i,1)
        FIMAX=FI(i,1)	 !!!!!FIMMAX
        do k=1,4
          IF(PSIMAX.LE.PSI(i,k)) PSIMAX=PSI(i,k)
          IF(PSIMIN.GT.PSI(i,k)) PSIMIN=PSI(i,k)
          IF(FIMAX.LE.FI(i,k))   FIMAX=FI(i,k)
          IF(FIMIN.GT.FI(i,k))   FIMIN=FI(i,k)
        enddo 

        IF(PSIMAX*PSIMIN.LE.0.0 .AND.
     *    (FIMAX*FIMIN.LE.0.0 .OR.(FIMAX.LT.0.0.AND.FIMIN.LT.0.0))) THEN
           NEC=NEC+1					   !LE 3.11.03  	!LE	3.11.03
           KELEM(NEC)=i
cc           write(6,*) NEC,KELEM(NEC)
C           pause
        ENDIF
      enddo              
C      EKSTENZIJA PRSLINE I REDEFINISANJE LS-funkcija NA OSNOVU NOVOG 
C      POLSEDNJEG SEGMENTA 
C      write(6,*)'UNESITE VREDNOST 1 UKOLIKO POSTOJI EKSTENZIJA PRSLINE I 
C     *            0 UKOLIKO NE POSTOJI EKSTENZIJA PRSLINE'
     
C      WRITE(6,*) ' INDIKATOR EKSTENZIJE='
C      READ(1,110)IND_EKST
C      WRITE(3,*) ' INDIKATOR EKSTENZIJE=',IND_EKST
C  110 FORMAT(I5)

C      IF(IND_EKST.EQ.1) THEN
C       POZIV POTPROGRAMA ZA EKSTENZIJU PRSLINE EKST_EDG      
C       call EKST_EDG(NUM,FI,PSI,KELEM,X,Y,NE)
C      ENDIF                 
C      =================================================================
C      =================================================================
C      write(6,*)'END LS_EDGE.FOR'                                                                
      return                       
      end
C=======================================================================
C
C==========================================================================
C
C==========================================================================
      SUBROUTINE XEC_YEC(NUM,FI,PSI,KELEM,X,Y,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
            
      include 'pakxfem.inc'
      
      DIMENSION NUM(NE,*),X(*),Y(*),FI(NE,*),PSI(NE,*),KELEM(*)

c      write(6,*)'START XEC_YEC.FOR'
C    ***************************************************************************      
C     PODPROGRAM XEC_YEC DEFINISE:
c        1.) BROJ PRESECNIH TACAKA PRSLINE I ELEMENATA
C        2.) KOORDINATE PRESECNIH TACAKA
c        3.) NIZOVE KOJI PAMTE CVOROVE (globalno numerisane) IZMEDJU KOJIH 
C                     SE PRESECNA TACKA NALAZI U x- i y-PRAVCU RESPEKTIVNO       
C        4.) NIZ KOJI PAMTI ELEMENTE U KOJIMA SE NALAZI PRESECNA TACKA
C    ***************************************************************************
c    ***************************************************************************
C     ULAZNE VELICINE U POTPROGRAM XEC_YEC SU:
C        NEC          -Broj elemenata presecenih sa diskontinuitetom
c        KELEM(NEC)   -Niz elemenata presecenih sa diskontinuitetom
C        NE           -Broj elemnata Narrow-Band-a
c        NUM(NE,4)    -Niz globalne numeracije cvorova po elementima 
c        X(NUM(NE,4)) -X-koordinate cvorova Narrow-Band-a  
c        Y(NUM(NE,4)) -Y-koordinate cvorova Narrow-Band-a        
c        PSI(NE,4)    -LS-funkcija PSI po cvorovima elemenata
c    ***************************************************************************
c    ***************************************************************************
c     IZLAZNE VELICINE POTPROGRAMA XEC_YEC SU:
c        NES                  -BROJ PRESECNIH TACAKA PRSLINE I ELEMENATA
C        XEC(NES),YEC(NES)    -KOORDINATE PRESECNIH TACAKA
c        NFX(NES,2),NFY(NES,2)-NIZOVE KOJI PAMTE CVOROVE (globalno numerisane) IZMEDJU KOJIH 
C                              SE PRESECNA TACKA NALAZI U x- i y-PRAVCU RESPEKTIVNO       
C        KELC(NES,4)          -NIZ KOJI PAMTI ELEMENTE U KOJIMA SE NALAZI PRESECNA TACKA
C    ***************************************************************************
C    ****************************************************************************         
C    *Inicijalizacija nizova
      do i=1,NEC+1
       NFX(i,1)=0
       NFX(i,2)=0
       NFY(i,1)=0
       NFY(i,2)=0
       KELC(i,1)=0
       KELC(i,2)=0
       KELC(i,3)=0
       KELC(i,4)=0
       XEC(i)=0.0
       YEC(i)=0.0
      enddo
C    *Nulovanje broja presecnih tacaka NES=0       
      NES=0                                 
      DO MM=1,NEC
       KEL=KELEM(MM)
C      PROVERI DA LI OVO MOZE DA RADI ZA ELEMNTE KOJI NISU CETVOROCVORNI?
       do m=1,4
         NX=0
         NY=0
         if(m.EQ.1) then
           k=1
           l=2
           NX=1
         endif
         if(m.EQ.2) then
           k=2
           l=3
           NY=1
         endif
         if(m.EQ.3) then
           k=3
           l=4
           NX=1
         endif
         if(m.EQ.4) then
           k=4
           l=1
           NY=1
         endif                             
         
c         IF(PSI(KEL,k)*PSI(KEL,l).LE.0.0 .AND. 
c     *     (FI(KEL,k)*FI(KEL,l).LE.0.0 .OR. 
c     *     (FI(KEL,k).LE.0.0 .AND.FI(KEL,l).LE.0.0))) THEN
         IF((PSI(KEL,k)*PSI(KEL,l).LT.0.0.AND.((FI(KEL,k).LE.0.0.AND.
     *      FI(KEL,l).LE.0.0).OR.FI(KEL,k)*FI(KEL,l).LE.0.0)).OR. 
     *     (FI(KEL,k).LE.0.0 .AND.DABS(PSI(KEL,k)).LT.1.E-9).OR.
     *     (FI(KEL,l).LE.0.0 .AND.DABS(PSI(KEL,l)).LT.1.E-9)) THEN
                  
           DEL=PSI(KEL,l)-PSI(KEL,k)
           
           IBRO=0
           KBRO=k
           LBRO=l
           call INVEST(NUM,NE)
           
           IF(IBRO.EQ.0) THEN
C            ================================================             
C            Definisanje XEC,YEC,NFX,NFY kada je PSI(KEL,k)=0.0 ili PSI(KEL,l)=0.0             
             if(DABS(PSI(KEL,k)).LT.1.E-9.AND.FI(KEL,k).LE.0.0) then
     
               NES=NES+1
               XEC(NES)=X(NUM(KEL,k))
               NFX(NES,1)=NUM(KEL,k)
               NFX(NES,2)=NUM(KEL,k)
               YEC(NES)=Y(NUM(KEL,k))
               NFY(NES,1)=NUM(KEL,k)
               NFY(NES,2)=NUM(KEL,k)
             endif

             if(DABS(PSI(KEL,l)).LT.1.E-9.AND.FI(KEL,l).LE.0.0) then
     
               NES=NES+1
               XEC(NES)=X(NUM(KEL,l))
               NFX(NES,1)=NUM(KEL,l)
               NFX(NES,2)=NUM(KEL,l)
               YEC(NES)=Y(NUM(KEL,l))
               NFY(NES,1)=NUM(KEL,l)
               NFY(NES,2)=NUM(KEL,l)    
             endif
C            ===============================================              
           ENDIF
C          =================================================           
C          =================================================              
C          Definisanje XEC i YEC kada DEL nije 0.0             
           IF(DABS(DEL).GT.1.E-9) THEN
              IBRO=0
              KBRO=k
              LBRO=l
              call INVEST(NUM,NE)
              if(IBRO.EQ.0)THEN
              
                 NES=NES+1
                 RAZX=X(NUM(KEL,l))-X(NUM(KEL,k))
                 RAZY=Y(NUM(KEL,l))-Y(NUM(KEL,k))         
               
                 XEC(NES)=X(NUM(KEL,k))-PSI(KEL,k)*RAZX/DEL   
                 YEC(NES)=Y(NUM(KEL,k))-PSI(KEL,k)*RAZY/DEL
                                                   
                 if(NX.EQ.1) then
                   NFX(NES,1)=NUM(KEL,k)
                   NFX(NES,2)=NUM(KEL,l)
                 endif

                if(NY.EQ.1) then
                   NFY(NES,1)=NUM(KEL,k)
                   NFY(NES,2)=NUM(KEL,l)
                endif
                
              endif

           ENDIF       
C          =========================================================================
C        Zatvaranje IF po PSI*PSI           
         ENDIF
C      Zatvaranje petlje po m           
       enddo    
C     Zatvaranje petlje po MM     
      ENDDO 
      IF(NES.GT.101) STOP 'NES.GT.101'   
C     ==========================================================================
C     DEFINISANJE ELEMENATA KOJIMA PRIPADA DATA PRESECNA TACKA
      DO N=1,NES
        IBR=0
        DO MM=1,NEC
          KEL=KELEM(MM)
          do m=1,4
            if(m.EQ.1) then
              k=1
              l=2
            else if(m.EQ.2) then
              k=2
              l=3
            else if(m.EQ.3) then
              k=3
              l=4
            else if(m.EQ.4) then
              k=4
              l=1
            endif

            IF(KELC(N,1).NE.KEL.AND.KELC(N,2).NE.KEL.AND.
     *         KELC(N,3).NE.KEL.AND.KELC(N,4).NE.KEL)THEN  
     
              IF(NFX(N,1).EQ.NUM(KEL,k).AND.NFX(N,2).EQ.NUM(KEL,l)) THEN
                 IBR=IBR+1
                 KELC(N,IBR)=KEL

              ELSE IF(NFX(N,1).EQ.NUM(KEL,l).AND.NFX(N,2).EQ.NUM(KEL,k))
     *        THEN           
                 IBR=IBR+1
                 KELC(N,IBR)=KEL

              ELSE IF(NFY(N,1).EQ.NUM(KEL,k).AND.NFY(N,2).EQ.NUM(KEL,l))
     *        THEN           
                 IBR=IBR+1
                 KELC(N,IBR)=KEL      
          
              ELSE IF(NFY(N,1).EQ.NUM(KEL,l).AND.NFY(N,2).EQ.NUM(KEL,k))
     *        THEN           
                 IBR=IBR+1
                 KELC(N,IBR)=KEL                                 

              ELSE IF(NFX(N,1).EQ.NUM(KEL,k).AND.NFX(N,2).EQ.NUM(KEL,k)
     *           .AND.NFY(N,1).EQ.NUM(KEL,k).AND.NFY(N,2).EQ.NUM(KEL,k))
     *        THEN         
                 IBR=IBR+1
                 KELC(N,IBR)=KEL                                 
                
              ELSE IF(NFX(N,1).EQ.NUM(KEL,l).AND.NFX(N,2).EQ.NUM(KEL,l)
     *           .AND.NFY(N,1).EQ.NUM(KEL,l).AND.NFY(N,2).EQ.NUM(KEL,l))
     *        THEN         
                 IBR=IBR+1
                 KELC(N,IBR)=KEL                                 
                 
              ENDIF
            ENDIF

C         zavrsetak petlje po m      
          enddo
c       zavrsetak petlje po MM          
        ENDDO
        IF(IBR.GT.4) STOP 'IBR.GT.4'
C     zavrsetak petlje po N
      ENDDO  

c      WRITE(6,*) NES,  'PRESECNIH TACAKA'       
C     Stampanje rezultata

c      do i=1,NES           
c       write(6,*) '  KOORDINATE PRESECNIH TACAKA'
c       write(6,10) i,XEC(i),i,YEC(i)             
C       write(6,*) ' GLOBALNO NUMERISANI CVOROVI IZMEDJU KOJIH SE NALAZI'
C       write(6,20) i,NFX(i,1),i,NFX(i,2) 
C       write(6,30) i,NFY(i,1),i,NFY(i,2)                               
c       write(6,*) ' ELEMENTI U KOJIMA SE NALAZI PRESECNA TACKA'
c       write(6,40) i,KELC(i,1),i,KELC(i,2),i,KELC(i,3),i,KELC(i,4)
c       pause
c      enddo
c   10 format(' XEC(',I3,')=',F15.10,' YEC(',I3,')=',F15.10)
c   20 format(' NFX(',I3,',1)=',I3,' NFX(',I3,',2)=',I3)     
c   30 format(' NFY(',I3,',1)=',I3,' NFY(',I3,',2)=',I3)    
c   40 format(' KELC(',I3,',1)=',I3,' KELC(',I3,',2)=',I3,
c     *       ' KELC(',I3,',3)=',I3,' KELC(',I3,',4)=',I3)     
     
      return
      end       
C==========================================================================
C
C==========================================================================
      SUBROUTINE INVEST(NUM,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      include 'pakxfem.inc'

      DIMENSION NUM(NE,*)
            
C      write(6,*)' KBRO=  LBRO='
C      WRITE(6,*) KBRO,LBRO
C      PAUSE  
      
      do N=1,NES 
      
        IF( NFX(N,1).EQ.NUM(KEL,KBRO).AND.NFX(N,2).EQ.NUM(KEL,LBRO).OR.
     *      NFX(N,1).EQ.NUM(KEL,LBRO).AND.NFX(N,2).EQ.NUM(KEL,KBRO).OR.
     *      NFY(N,1).EQ.NUM(KEL,KBRO).AND.NFY(N,2).EQ.NUM(KEL,LBRO).OR.
     *      NFY(N,1).EQ.NUM(KEL,LBRO).AND.NFY(N,2).EQ.NUM(KEL,KBRO))
     *      THEN
            IBRO=IBRO+1
            
        ELSE IF(NFX(N,1).EQ.NUM(KEL,KBRO).AND.NFX(N,2).EQ.NUM(KEL,KBRO)
     *     .AND.NFY(N,1).EQ.NUM(KEL,KBRO).AND.NFY(N,2).EQ.NUM(KEL,KBRO))
     *      THEN
            IBRO=IBRO+1
            
        ELSE IF(NFX(N,1).EQ.NUM(KEL,LBRO).AND.NFX(N,2).EQ.NUM(KEL,LBRO)
     *     .AND.NFY(N,1).EQ.NUM(KEL,LBRO).AND.NFY(N,2).EQ.NUM(KEL,LBRO)) 
     *      THEN
            IBRO=IBRO+1   
            
        ENDIF
      enddo 
      
      return
      end
C==========================================================================
C
C==========================================================================
      SUBROUTINE CRT_LS(FIE,PSIE,NUM,FI,PSI,X,Y,NE,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
	
      DIMENSION FIE(*),PSIE(*)
      DIMENSION NUM(NE,*),X(*),Y(*),FI(NE,*),PSI(NE,*)
      
C     ======================================================
C     PP SLUZI ZA PRAVLJENJE DATOTEKE ZA CRTANJE U TEC_PLOT_u
C     CRTA SE:
C     LS-FUNKCIJA FI
C     LS-FUNKCIJA PSI
C     NES BROJ PRESECNIH TACAKA PRSLINE I ELEMENATA        
C     ======================================================
C     CRTANJE TEC-PLOT-u
      IGR=99
      OPEN(UNIT=IGR,FILE='OUTPUT.DAT') 
      WRITE(IGR,*) 'TITLE = "Example: LEVEL SET 2D Plot"'
C      WRITE(IGR,*) 'VARIABLES = "X", "Y", "Psi"'
C      WRITE(IGR,1030) NNB+1,MNB+1
      DO J=1,NE
       DO M=1,4
         PSIE(NUM(J,M))=PSI(J,M) 
         FIE(NUM(J,M))=FI(J,M) 
       ENDDO
      ENDDO
C      DO I=1,NP
C         WRITE(IGR,1000) X(I),Y(I),PSIE(I)
C      ENDDO	   	   
C      FORMAT(I6,10I3,3(1PE13.5),I2)
      WRITE(IGR,*) 'VARIABLES = "X", "Y", "Z"'
      WRITE(IGR,1020) NNB+1,MNB+1 
      DO I=1,NP
	 IF(X(I).GT.129.999 .AND. X(I).LT.245.001) THEN
         WRITE(IGR,1000) X(I),Y(I),Z
	 ENDIF !1. uslov za xfem3, xfem31, xfem4, xfem41
      ENDDO  
      WRITE(IGR,*) 'VARIABLES = "X", "Y", "CRACK"'
      WRITE(IGR,1010) NSEGX+1

      Z=0.
c      DO I=1,NES
c         WRITE(IGR,1000) Xec(I),Yec(I),Z
c      ENDDO	   	   
      DO I=1,NSEGX+1
         WRITE(IGR,1000) Xp(I),Yp(I),Z
      ENDDO	   	   
C      WRITE(IGR,1000) Xp(NSEGX+1),Yp(NSEGX+1),Z
 1000 FORMAT(F15.10,F15.10,F15.10)
 1010 FORMAT('ZONE T="ZONE XecYec", I=',I2,', J=1, F=POINT')
 1020 FORMAT('ZONE T="ZONE FI", I=',I2,',J=',I2,' F=POINT')
 1030 FORMAT('ZONE T="ZONE PSI", I=',I2,',J=',I2,' F=POINT')
      CLOSE(IGR)      
C     ======================================================
      RETURN
      END       
C=======================================================================
C
C=======================================================================
      SUBROUTINE ENRICH1(CORD,NODTIP,NSSN,HZNAK,HNOD,
     1                  NSSE,FI,PSI,KELEM,NUM,NE,NP,K1,K2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
      
      DIMENSION CORD(NP,*),NUM(NE,*),NODTIP(*),NSSN(*),HNOD(*)
      DIMENSION HZNAK(NE,4),NSSE(*),FI(NE,*),PSI(NE,*),KELEM(*)
	DIMENSION V21(3),V31(3), V41(3)
	DIMENSION T1(3),T2(3),T3(3), T4(3)
	DIMENSION V1(3),V2(3),V3(3),V4(3),D01(3), D02(3) 
CCCCCCCCCCC     BOGACENJE SAMO SA HEVISAJDOVOM FUNKCIOM    CCCCCCCCCCCC      
C     ==================================================================
C     ULAZNE VELICINE:
C     NB          BROJ CVOROVA U DOMENU
C     NEC         BROJ ELEMENATA PRESECENIH SA PRSLINOM
C     KELEM(NEC) 	NIZ KOJI PAMTI ELEMENTE PRESECENE SA PRSLINOM
C     FIE         FUNKCIJA FI PO CVOROVIMA DOMENA
C     PSIE        FUNKCIJA PSI PO CVOROVIMA DOMENA
C     ======================================================================
C     ======================================================================
C     IZLAZNE VELICINE:
C     DEFINISANJE VRSTE NODA SVI NODOVI U DOMENU IMAJU SVOJU VREDNOST NODTIP(NP)
C     NODTIP(I)=0 OBICAN NOD
C     NODTIP(I)<0 OBOGACEN SA H
C     HZNAK(NE,4) HEVISAJDOVA FUNKCIJA PO CVOROVIMA ELEMENTA
C     NSSN(NB)    BROJ STEPENI SLOBODE PO CVOROVIMA DOMENA
C     NSSE(NE)    BROJ STEPENI SLOBODE U ELEMENTU
C     ======================================================================= 
      NB=NP
c     DUZINA PRSLINE
      DPRS=0.0 
      DO I=1,NSEGX
	DPRS=DPRS+DSQRT((Xp(I+1)-Xp(I))*(Xp(I+1)-Xp(I))+
     1                (Yp(I+1)-Yp(I))*(Yp(I+1)-Yp(I)) )
	ENDDO

	DPRS=1.05*DPRS
	DPRS=-DPRS
c     DEFINISANJE ELEMENTA U KOJEM JE VRH IVRHELEM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO I=1,NEC

	J0=KELEM(I)

 
      V21(1)=CORD(NUM(J0,2),1)-CORD(NUM(J0,1),1)
      V21(2)=CORD(NUM(J0,2),2)-CORD(NUM(J0,1),2)
	V21(3)=0.0
	V31(1)=CORD(NUM(J0,3),1)-CORD(NUM(J0,1),1)
      V31(2)=CORD(NUM(J0,3),2)-CORD(NUM(J0,1),2)
	V31(3)=0.0
      V41(1)=CORD(NUM(J0,4),1)-CORD(NUM(J0,1),1)
      V41(2)=CORD(NUM(J0,4),2)-CORD(NUM(J0,1),2)
	V41(3)=0.0
 

      T1(1)=CORD(NUM(J0,1),1)-Xp(NSEGX+1)
      T1(2)=CORD(NUM(J0,1),2)-Yp(NSEGX+1)

      T2(1)=CORD(NUM(J0,2),1)-Xp(NSEGX+1)
      T2(2)=CORD(NUM(J0,2),2)-Yp(NSEGX+1)

      T3(1)=CORD(NUM(J0,3),1)-Xp(NSEGX+1)
      T3(2)=CORD(NUM(J0,3),2)-Yp(NSEGX+1)

      T4(1)=CORD(NUM(J0,4),1)-Xp(NSEGX+1)
      T4(2)=CORD(NUM(J0,4),2)-Yp(NSEGX+1)

C     UNULJAVANJE TRECEG CLANA
	T1(3)=0.0
	T2(3)=0.0
	T3(3)=0.0
	T4(3)=0.0

C     UNULJAVANJE SVIH CLANOVA NIZA
	DO M=1,3
	   D01(M)=0.0
	   D02(M)=0.0
	   V1(M)=0.0
	   V2(M)=0.0
	   V3(M)=0.0
	   V4(M)=0.0
      ENDDO
C     VEKTORSKI PROIZVODI
	CALL AXBV(V21,V31,D01)
	CALL AXBV(V31,V41,D02)
      CALL AXBV(T2,T1,V1)
	CALL AXBV(T3,T2,V2)
	CALL AXBV(T4,T3,V3)
	CALL AXBV(T1,T4,V4)

C     SKALARNI PROIZVOD
      P01=DABS(D01(3))
	P02=DABS(D02(3))
	Pelem=0.5*(P01+P02)
	P1=0.5*DABS(V1(3))
	P2=0.5*DABS(V2(3))
	P3=0.5*DABS(V3(3))
	P4=0.5*DABS(V4(3))

	SUMA=P1+P2+P3+P4
C      WRITE(6,*) 'J0',J0
C      WRITE(6,*) 'Pelem,Suma',Pelem,SUMA
C      PAUSE

	IF( DABS(Pelem-SUMA).LT.10.E-9 
     1    ) THEN
	   IVRHELEM=J0 
	ENDIF
 
      ENDDO
c	WRITE(6,*) 'IVRHELEM', IVRHELEM
c	PAUSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO I=1,NEC


	!DEFINISANJE STRANICA ELEMENATA
	J0=KELEM(I)
	DE1=DSQRT((CORD(NUM(J0,1),1)-CORD(NUM(J0,2),1))**2+
     1          (CORD(NUM(J0,1),2)-CORD(NUM(J0,2),2))**2 )

	DE2=DSQRT((CORD(NUM(J0,1),1)-CORD(NUM(J0,4),1))**2+
     1          (CORD(NUM(J0,1),2)-CORD(NUM(J0,4),2))**2 )
	!DEFINISANJE DIJAGONALE ELEMENTA KOJI JE PRESECEN SA PRSLINOM
c	J0=KELEM(I)
c	DE1=DSQRT((CORD(NUM(J0,1),1)-CORD(NUM(J0,3),1))**2+
c     1          (CORD(NUM(J0,1),2)-CORD(NUM(J0,3),2))**2 )

c	DE2=DSQRT((CORD(NUM(J0,2),1)-CORD(NUM(J0,4),1))**2+
c     1          (CORD(NUM(J0,2),2)-CORD(NUM(J0,4),2))**2 )

	IF(DE1.gt.DE2) THEN
	   D1=0.800*DE1	 !0.80D0  0.750D0
      ELSE
	   D1=0.800*DE2
      ENDIF

        PSIMIN=PSI(KELEM(I),1)
        PSIMAX=PSI(KELEM(I),1)
        FIMIN=FI(KELEM(I),1)
        FIMAX=FI(KELEM(I),1)

        DO J=1,4
          IF(PSIMAX.LE.PSI(KELEM(I),J))  PSIMAX=PSI(KELEM(I),J)
          IF(PSIMIN.GT.PSI(KELEM(I),J))  PSIMIN=PSI(KELEM(I),J)
          IF(FIMAX.LE.FI(KELEM(I),J))    FIMAX=FI(KELEM(I),J)
          IF(FIMIN.GT.FI(KELEM(I),J))    FIMIN=FI(KELEM(I),J)
        ENDDO	    							   

        DO J=1,4          !LE LT                         !LT LE
	    KALAM=0
          IF(PSIMAX*PSIMIN.LT.0.0 .AND. FI(KELEM(I),J).LE.0.0
     1       .AND. DABS(PSI(KELEM(I),J)).LT.D1  
     1       .AND. FI(KELEM(I),J).GT.DPRS 
     1       .AND.  KELEM(I).NE. IVRHELEM
     1       ) THEN	
             NODTIP(NUM(KELEM(I),J))=-1      
          ENDIF

	    !PROVERA OBOGACENOSTI  7.10.2004. 
          IF(NODTIP(NUM(KELEM(I),J)).NE.0) THEN
            DO MM=1,NES
	       IF( NUM(KELEM(I),J).EQ.NFX(MM,1).OR.
     1	       NUM(KELEM(I),J).EQ.NFX(MM,2).OR.
     1           NUM(KELEM(I),J).EQ.NFY(MM,1).OR.
     1           NUM(KELEM(I),J).EQ.NFY(MM,2)) THEN
                KALAM=1
             ENDIF
            ENDDO
          ENDIF
          IF(NODTIP(NUM(KELEM(I),J)).NE.0 .AND. KALAM.EQ.0    
     1	   .AND. DABS(PSI(KELEM(I),J)).GT.D1) THEN
	      NODTIP(NUM(KELEM(I),J))=0	  !AKO JE CVOR PO DIJAGONALI DALEKO VRACA
          ENDIF                         !SE NA NULU OBOGACENOST  
				        	  		  
        ENDDO

	!SLUCAJ KADA SE PRSLINA POKLAPA SA ELEMENTOM
        DO J=1,4                
          KALAM=0
          IF((PSIMAX*PSIMIN).EQ.0.0.AND.    !LT.1.E-9  
     1       (PSI(KELEM(I),J)).EQ.0.0 .AND. !LT.1.E-9 
     1       FI(KELEM(I),J).LT.0.0 .AND.
     1       FI(KELEM(I),J).GT.DPRS) THEN
             NODTIP(NUM(KELEM(I),J))=-1
          ENDIF
	    !PROVERA OBOGACENOSTI  7.10.2004. 
          IF(NODTIP(NUM(KELEM(I),J)).NE.0) THEN
            DO MM=1,NES
	       IF( NUM(KELEM(I),J).EQ.NFX(MM,1).OR.
     1	       NUM(KELEM(I),J).EQ.NFX(MM,2).OR.
     1           NUM(KELEM(I),J).EQ.NFY(MM,1).OR.
     1           NUM(KELEM(I),J).EQ.NFY(MM,2)) THEN
                KALAM=1
             ENDIF
            ENDDO
          ENDIF
          IF(NODTIP(NUM(KELEM(I),J)).NE.0 .AND. KALAM.EQ.0    
     1	   .AND. DABS(PSI(KELEM(I),J)).GT.D1) THEN
	      NODTIP(NUM(KELEM(I),J))=0	  !AKO JE CVOR PO DIJAGONALI DALEKO VRACA
          ENDIF                         !SE NA NULU OBOGACENOST  
		
				        	  		  
        ENDDO

      ENDDO ! ZAVRSENA PETLJA PO ELEMENTIMA KOJI SU PRESECENI SA PRSLINOM

	!7.10.2004.	KOREKCIJA OBOGACENIH CVOROVA
C	DO J=1,NB
C	 KALAM=0 
C	 IF(NODTIP(J).NE.0) THEN
C          DO I=1,NES
C	       IF( J.EQ.NFX(I,1).OR.J.EQ.NFX(I,2).OR.
C     1           J.EQ.NFY(I,1).OR.J.EQ.NFY(I,2)) THEN
C                KALAM=1
C             ENDIF
C          ENDDO
C        ENDIF
C        IF(NODTIP(J).NE.0 .AND. KALAM.EQ.0    !!AND.DABS(PSI(J)).GT.D1
C     1	  ) THEN
C	     NODTIP(J)=0
C        ENDIF
C      ENDDO	  		 	  				   			  		      


C     BROJ STEPENI SLOBODE PO CVORU I UKUPAN BROJ STEPENI SLOBODE
      KTIP1=0
      KTIP2=0
      NSS=0
      DO I=1,NB
        IF(NODTIP(I).EQ.0)  THEN
          NSSN(I)=2
          NSS=NSS+2
        ELSEIF(NODTIP(I).LT.0)  THEN
          KTIP1=KTIP1+1
          NODTIP(I)=-KTIP1
          NSSN(I)=4
          NSS=NSS+4
        ELSEIF(NODTIP(I).GT.0)  THEN
          KTIP2=KTIP2+1
          NIZ2(KTIP2)=I
          NODTIP(I)=KTIP2
          NSSN(I)=10
          NSS=NSS+10
        ENDIF
      ENDDO


c      DO I=1,NB
c        IF(NODTIP(I).NE.0)  THEN
c          WRITE(6,*) I,NODTIP(I)
c        ENDIF 		   
c      ENDDO
c	PAUSE
c     ==================================================================
c     ==================================================================
C     HEVISAJDOVA FUNKCIJA BOGACENJA
C     PETLJA PO SVIM ELEMENTIMA U DOMENU

      DO I=1,NE
C       BROJ STEPENI SLOBODE PO ELEMENTU
        NSSE(I)=0
C       NALAZENJE MAKSIMALNE I MINIMALNE VREDNOSTI PO ELEMENTU
        PSIMIN=PSI(I,1)
        PSIMAX=PSI(I,1)
        FIMIN=FI(I,1)
        FIMAX=FI(I,1)
        DO J=1,4
          NSSE(I)=NSSE(I)+NSSN(NUM(I,J))
C
          IF(PSIMAX.LE.PSI(I,J))  PSIMAX=PSI(I,J)
          IF(PSIMIN.GT.PSI(I,J))  PSIMIN=PSI(I,J)
          IF(FIMAX.LE.FI(I,J))    FIMAX=FI(I,J)
          IF(FIMIN.GT.FI(I,J))    FIMIN=FI(I,J)
        ENDDO	    							   

        DO J=1,4

         IF( NODTIP(NUM(I,J)).LT.0) THEN
C OVI USLOVI DA SE PROVERE, SIGURNO NE VALJA ZA =0? ZILE         
          IF(DABS(PSI(I,J)).EQ.0.0 .AND.	!LT.1.E-9
     1       PSIMAX.GT.0.0 .AND. DABS(PSIMAX*PSIMIN).EQ.0.0) THEN
            ZNAK=1.D0
            ZNAK1=PSIMAX
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
	                 !LT.1.E-9
          ELSEIF(DABS(PSI(I,J)).EQ.0.0 .AND.			   !LT.1.E-9
     1           PSIMIN.LT.0.0 .AND. DABS(PSIMAX*PSIMIN).EQ.0.0) THEN
            ZNAK=1.D0
            ZNAK1=PSIMIN
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
            	                  !LT.1.E-9
          ELSEIF(DABS(PSI(I,J)).EQ.0.0 .AND. PSIMAX*PSIMIN.LT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSIMIN    !!!DVOZNACNO JE MOZE DA SE STAVI I PSIMAX
            VRED=DSIGN(ZNAK,ZNAK1)     !!!CVOR 1,15,29 PRIMER 4 i PRIMER 8
            HZNAK(I,J)=-1.D0    !!! PRSLINA PROLAZI PO CVORU 4.4.2004 (VRED)
	                          !GT.1.E-9
          ELSEIF(DABS(PSI(I,J)).GT.0.0 .AND. PSIMAX*PSIMIN.LT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
	                           !GT.1.E-9
          ELSEIF(DABS(PSI(I,J)).GT.0.0 .AND. PSIMAX*PSIMIN.GT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
	                           !GT.1.E-9
          ELSEIF(DABS(PSI(I,J)).GT.0.0 .AND.
     1           DABS(PSIMAX*PSIMIN).EQ.0.0) THEN   !LT.1.E-9
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
	    
          ENDIF
         ENDIF
        ENDDO
      ENDDO 

C     DEFINISANJE HEVISAJDOVE FUNKCIJE PO CVOROVIMA
      DO I=1,NE
        DO J=1,4
CV         KONTROLNA STAMPA - VLADE
c          WRITE(3,*)'HZNAK(I,J),I',HZNAK(I,J),I
          HNOD(NUM(I,J))=HZNAK(I,J)
        ENDDO 
      ENDDO		  

C     ==========================================================
C     ==========================================================
C     ! KOSINUSI PRAVCA PRSLINE U ODNOSU NA GLOBALNI KOORD. SIS.
      VX1(1)=Xp(NSEGX+1)-Xp(NSEGX)
      VX1(2)=Yp(NSEGX+1)-Yp(NSEGX)
      VV=dsqrt(VX1(1)*VX1(1)+VX1(2)*VX1(2))
      VX1(1)=VX1(1)/VV
      VX1(2)=VX1(2)/VV
C     ! JEDINICNI VEKTORI PRAVCA-X2 PRSLINE NCR
      VX2(1)=-1.D0*VX1(2)	   !OVO NIJE DIRANO
      VX2(2)=VX1(1)
C     ======================================================
	
      K1=KTIP1
      K2=KTIP2
c         DO N=1,NE
c            WRITE(3,*) N,(NUM(N,JJ),JJ=1,4)
c            WRITE(3,*) N,(NODTIP(NUM(N,JJ)),JJ=1,4)
c         ENDDO
      RETURN
      END
C=======================================================================
      SUBROUTINE ENRICH2(CORD,FIE,PSIE,NODTIP,NSSN,HZNAK,HNOD,
     1                  NSSE,FI,PSI,KELEM,NUM,NE,NP,K1,K2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
      
      DIMENSION CORD(NP,*),NUM(NE,*),FIE(*),PSIE(*),NODTIP(*)
      DIMENSION HZNAK(NE,*),NSSE(*),FI(NE,*),PSI(NE,*),KELEM(*)
      DIMENSION NSSN(*),HNOD(*)      
C     ==================================================================
C     ULAZNE VELICINE:
C     NB          BROJ CVOROVA U DOMENU
C     NEC         BROJ ELEMENATA PRESECENIH SA PRSLINOM
C     KELEM(NEC) 	NIZ KOJI PAMTI ELEMENTE PRESECENE SA PRSLINOM
C     FIE         FUNKCIJA FI PO CVOROVIMA DOMENA
C     PSIE        FUNKCIJA PSI PO CVOROVIMA DOMENA
C     ======================================================================
C     ======================================================================
C     IZLAZNE VELICINE:
C     DEFINISANJE VRSTE NODA SVI NODOVI U DOMENU IMAJU SVOJU VREDNOST NODTIP(NP)
C     NODTIP(I)=0 OBICAN NOD
C     NODTIP(I)<0 OBOGACEN SA H
C     NODTIP(I)>0 OBOGACEN SA NEAR TIP FUNKCIJAMA
C     NIZ2(10)    PAMTI GLOBALNU NUMERACIJU CVOROVA OBOGACENIH SA TIPOM 2
C     KTIP2       BROJ CVOROVA U DOMENU OBOGACENIH SA TIPOM 2
C     HZNAK(NE,4) HEVISAJDOVA FUNKCIJA PO CVOROVIMA ELEMENTA
C     NSSN(NB)    BROJ STEPENI SLOBODE PO CVOROVIMA DOMENA
C     NSSE(NE)    BROJ STEPENI SLOBODE U ELEMENTU
C     ======================================================================= 
      NB=NP
      DPRS=0.0 
      DO I=1,NSEGX
	DPRS=DPRS+DSQRT((Xp(I+1)-Xp(I))*(Xp(I+1)-Xp(I))+
	1                (Yp(I+1)-Yp(I))*(Yp(I+1)-Yp(I)) )
	ENDDO
	DPRS=1.005*DPRS
	DPRS=-DPRS

	
      DO I=1,NEC

	!DEFINISANJE DIJAGONALE ELEMENTA KOJI JE PRESECEN SA PRSLINOM
c	J0=KELEM(I)
c	DE1=DSQRT((CORD(NUM(J0,1),1)-CORD(NUM(J0,3),1))**2+
c     1          (CORD(NUM(J0,1),2)-CORD(NUM(J0,3),2))**2 )
c
c	DE2=DSQRT((CORD(NUM(J0,2),1)-CORD(NUM(J0,4),1))**2+
c     1          (CORD(NUM(J0,2),2)-CORD(NUM(J0,4),2))**2 )
c
	J0=KELEM(I)
	DE1=DSQRT((CORD(NUM(J0,1),1)-CORD(NUM(J0,2),1))**2+
     1          (CORD(NUM(J0,1),2)-CORD(NUM(J0,2),2))**2 )

	DE2=DSQRT((CORD(NUM(J0,3),1)-CORD(NUM(J0,4),1))**2+
     1          (CORD(NUM(J0,3),2)-CORD(NUM(J0,4),2))**2 )

	IF(DE1.GT.DE2) THEN
	   D1=0.800D0*DE2
      ELSE
	   D1=0.80D0*DE1
      ENDIF


        PSIMIN=PSI(KELEM(I),1)
        PSIMAX=PSI(KELEM(I),1)
        FIMIN=FI(KELEM(I),1)
        FIMAX=FI(KELEM(I),1)
C       PROVERI DA LI OVO RADI ZA BROJ CVOROVA RAZLICIT OD 4?        
        DO J=1,4
          IF(PSIMAX.LE.PSI(KELEM(I),J))  PSIMAX=PSI(KELEM(I),J)
          IF(PSIMIN.GT.PSI(KELEM(I),J))  PSIMIN=PSI(KELEM(I),J)
          IF(FIMAX.LT.FI(KELEM(I),J))    FIMAX=FI(KELEM(I),J)
          IF(FIMIN.GT.FI(KELEM(I),J))    FIMIN=FI(KELEM(I),J)
        ENDDO	    							   
	!BOGACENJE NT FUNKCIJOM                    !LE                       
	IF(PSIMAX*PSIMIN.LT.0.D0 
	1   .AND. FIMAX*FIMIN.LT.0.D0 
     1   )THEN
          DO J=1,4
c           bogacenje samo ispred vrha
	      IF( !FI(KELEM(I),J).LT.DE1.AND.
     1       FI(KELEM(I),J).GT.0.0 	
     1		  ) THEN  
               NODTIP(NUM(KELEM(I),J))=2
            ENDIF
          ENDDO
      ENDIF
C	BOGACENJE HEVISAJDOVOM FUNKCIJOM
        DO J=1,4          !LE                          !LE
          IF(PSIMAX*PSIMIN.LT.0.D0 .AND. FI(KELEM(I),J).LT.0.0
     1       .AND. NODTIP(NUM(KELEM(I),J)).LT.2 
     1       .AND. DABS(PSI(KELEM(I),J)).LT.D1  
     1       .AND. FI(KELEM(I),J).GT.DPRS ) THEN	
             NODTIP(NUM(KELEM(I),J))=-1      
          ENDIF		        	  		  
        ENDDO

C        !!!!KOREKCIJA DA LI JE VRSITI PROVERITI U TOKU TESTIRANJA
C        !!!!OVO MOZE DA SE VRSI KADA JE VRH NA IVICI ELEMENTA a i ne mora
C        DO J=1,4                                        !LE
C          IF(PSIMAX*PSIMIN.LT.0.0 .AND.DABS(FIMAX*FIMIN).LT.1.E-9 .AND.
C     1       FI(KELEM(I),J).LT.0.0 .AND.
C     1                               NODTIP(NUM(KELEM(I),J)).EQ.2) THEN
C            NODTIP(NUM(KELEM(I),J))=-1
C          ENDIF
C        ENDDO

ccccccccccccccccccccPRSLINA DUZ IVICE ELEMENTAcccccccccccccccccccccccc

C       !  BOGACENJE SA NT U DOMENU (-90, 90 STEPENI) ISPRED VRHA
        !  ZA PRSLINU KOJA SE NALAZI DUZ IVICE ELEMENTA 
        DO J=1,4
     
          IF(DABS(PSIMAX*PSIMIN).LT.1.E-9 .AND.
     1       ((DABS(PSI(KELEM(I),J)).GT.1.E-9 .AND.
     1                              DABS(FI(KELEM(I),J)).LT.1.E-9) .OR. 
     1        (DABS(PSI(KELEM(I),J)).GT.1.E-9 .AND.
     1                                   FI(KELEM(I),J).GT.0.0) .OR. 
     1        (DABS(PSI(KELEM(I),J)).LT.1.E-9 .AND.
     1                                   FI(KELEM(I),J).GT.0.0) )) THEN
            NODTIP(NUM(KELEM(I),J))=2    !TREBA 2 
          ENDIF		        	  		  
        ENDDO
	!     BOGACENJE SA H FUNKCIJOM ZA PRSLINU DUZ IVICA ELEMENATA
        DO J=1,4                  !LE
          IF(DABS(PSIMAX*PSIMIN).LT.1.E-9 .AND. 
     1       DABS(PSI(KELEM(I),J)).LT.1.E-9 .AND.
     1       NODTIP(NUM(KELEM(I),J)).LT.2
     1            ) THEN
            NODTIP(NUM(KELEM(I),J))=-1
          ENDIF		        	  		  
        ENDDO

	!SPECIJALAN SLUCAJ KADA SE VRH NALAZI U TEMENU
	!BOGACENJE SA NT SAMO ISPRED VRHA 7.4.2004.
C        DO J=1,4          
C          IF(DABS(FI(KELEM(I),J)).LT.1.E-9 .AND.
C     1      DABS(PSI(KELEM(I),J)).GT.1.E-9 .AND.
C     1      NODTIP(NUM(KELEM(I),J)).EQ.2)THEN
C            NODTIP(NUM(KELEM(I),J))=0 
C          ENDIF	                    
C        ENDDO

C  !    POSEBAN SLUCAJ KOREKCIJE 29.10.2003. PRIMER 8
C  !    NE DOZVOLJAVA DA SE NT-FUNKCIJE NADJU U CVORU NA PRSLINI
        DO J=1,4
          IF(NODTIP(NUM(KELEM(I),J)).EQ.2 .AND.
     1       DABS(PSI(KELEM(I),J)).LT.1.E-9 .AND.
     1       FI(KELEM(I),J).LE.0.0 ) THEN
            NODTIP(NUM(KELEM(I),J))=0    !0 NIJE OBOGACEN VRH
          ENDIF						   
        ENDDO
	!BOGACENJE SA NT FUNKCIJOM IZA VRHA 
	!KADA SE PRSLINA NALAZI DUZ IVICE ELEMENTA 
CC        IF(DABS(PSIMAX*PSIMIN).LT.1.E-9) THEN 
CC	     NT=0
CC	     DO J=1,4
CC	       IF(NODTIP(NUM(KELEM(I),J)).EQ.2) THEN
CC	          NT=NT+1
CC		   ENDIF	     
CC	     ENDDO
CC	     IF(NT.GT.0) THEN
CC            DO J=1,4
CC	         IF(NODTIP(NUM(KELEM(I),J)).EQ.0) THEN
CC	            NODTIP(NUM(KELEM(I),J))=2
CC		     ENDIF	     
CC	 	   ENDDO
CC		  ENDIF
CC	   ENDIF


ccccccccccccccccccccPRSLINA DUZ IVICE ELEMENTAcccccccccccccccccccccccc
cc	!NEMA BOGACENJA  VRHA  KADA
cc	! SE PRSLINA NALAZI DUZ IVICA ELEMENATA    
c        IF(PSIMAX*PSIMIN.LT.0.0 .OR. DABS(PSIMAX*PSIMIN).LT.1.E-9) THEN 
c 	     DO J=1,4
c	   	  IF( DABS(FI(KELEM(I),J)).LT.1.E-9 .AND.
c     1          DABS(PSI(KELEM(I),J)).LT.1.E-9 ) THEN	     
c	      	    NODTIP(NUM(KELEM(I),J))=0
c             ENDIF
c		  ENDDO
c         ENDIF
	   
	   		     				   	 
      ENDDO ! ZAVRSENA PETLJA PO ELEMENTIMA


C     BROJ STEPENI SLOBODE PO CVORU I UKUPAN BROJ STEPENI SLOBODE
      KTIP1=0
      KTIP2=0
      NSS=0
      DO I=1,NB
        IF(NODTIP(I).EQ.0)  THEN
          NSSN(I)=2
          NSS=NSS+2
        ELSEIF(NODTIP(I).LT.0)  THEN
          KTIP1=KTIP1+1
          NODTIP(I)=-KTIP1
          NSSN(I)=4
          NSS=NSS+4
        ELSEIF(NODTIP(I).GT.0)  THEN
          KTIP2=KTIP2+1
          NIZ2(KTIP2)=I
          NODTIP(I)=KTIP2
          NSSN(I)=10
          NSS=NSS+10
        ENDIF
      ENDDO
      IF(KTIP2.GT.10) STOP 'KTIP2.GT.10'

c     DO I=1,NB
c       IF(NODTIP(I).NE.0)  THEN
c         WRITE(6,*) I,NODTIP(I)
c       ENDIF 		   
c     ENDDO

c     ==================================================================
c     ==================================================================
C     HEVISAJDOVA FUNKCIJA BOGACENJA
C     PETLJA PO SVIM ELEMENTIMA U DOMENU

      DO I=1,NE
C       BROJ STEPENI SLOBODE PO ELEMENTU
        NSSE(I)=0
C       NALAZENJE MAKSIMALNE I MINIMALNE VREDNOSTI PO ELEMENTU
        PSIMIN=PSI(I,1)
        PSIMAX=PSI(I,1)
        FIMIN=FI(I,1)
        FIMAX=FI(I,1)
        DO J=1,4
          NSSE(I)=NSSE(I)+NSSN(NUM(I,J))
C
          IF(PSIMAX.LE.PSI(I,J))  PSIMAX=PSI(I,J)
          IF(PSIMIN.GT.PSI(I,J))  PSIMIN=PSI(I,J)
          IF(FIMAX.LE.FI(I,J))    FIMAX=FI(I,J)
          IF(FIMIN.GT.FI(I,J))    FIMIN=FI(I,J)
        ENDDO	    							   

        DO J=1,4

         IF( NODTIP(NUM(I,J)).LT.0) THEN
         
          IF(DABS(PSI(I,J)).LT.1.E-9 .AND.
     1       PSIMAX.GT.0.0 .AND. DABS(PSIMAX*PSIMIN).LT.1.E-9) THEN
            ZNAK=1.D0
            ZNAK1=PSIMAX
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED

          ELSEIF(DABS(PSI(I,J)).LT.1.E-9 .AND.
     1           PSIMIN.LT.0.0 .AND. DABS(PSIMAX*PSIMIN).LT.1.E-9) THEN
            ZNAK=1.D0
            ZNAK1=PSIMIN
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED

          ELSEIF(DABS(PSI(I,J)).LT.1.E-9 .AND. PSIMAX*PSIMIN.LT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSIMIN    !!!DVOZNACNO JE MOZE DA SE STAVI I PSIMAX
            VRED=DSIGN(ZNAK,ZNAK1)     !!!CVOR 1,15,29 PRIMER 4 i PRIMER 8
            HZNAK(I,J)=VRED  !!!0.D0, VRED, VRED PRSLINA PRELAZI PO CVORU

          ELSEIF(DABS(PSI(I,J)).GT.1.E-9 .AND. PSIMAX*PSIMIN.LT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED

          ELSEIF(DABS(PSI(I,J)).GT.1.E-9 .AND. PSIMAX*PSIMIN.GT.0.0)THEN
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED

          ELSEIF(DABS(PSI(I,J)).GT.1.E-9 .AND.
     1           DABS(PSIMAX*PSIMIN).LT.1.E-9) THEN
            ZNAK=1.D0
            ZNAK1=PSI(I,J)
            VRED=DSIGN(ZNAK,ZNAK1)
            HZNAK(I,J)=VRED
	    
          ENDIF
         ENDIF

        ENDDO
      ENDDO 

C     DEFINISANJE HEVISAJDOVE FUNKCIJE PO CVOROVIMA
      DO I=1,NE
        DO J=1,4
	     HNOD(NUM(I,J))=HZNAK(I,J)
        ENDDO 
      ENDDO		  

C     ==========================================================
C     ==========================================================
C     ! KOSINUSI PRAVCA PRSLINE U ODNOSU NA GLOBALNI KOORD. SIS.
      VX1(1)=Xp(NSEGX+1)-Xp(NSEGX)
      VX1(2)=Yp(NSEGX+1)-Yp(NSEGX)
      VV=dsqrt(VX1(1)*VX1(1)+VX1(2)*VX1(2))
      VX1(1)=VX1(1)/VV
      VX1(2)=VX1(2)/VV
C     ! JEDINICNI VEKTORI PRAVCA-X2 PRSLINE NCR
      VX2(1)=-1.D0*VX1(2)	   !OVO NIJE DIRANO
      VX2(2)=VX1(1)
c     ==================================================================
C     BOGACENJE SA PSI FUNKCIJAMA
      DO K=1,KTIP2
        I=NIZ2(K)
c	  WRITE(6,*)'K,NIZ2(K),I',K,NIZ2(K),I
        ROO(K)=DSQRT((PSIE(I))**2+(FIE(I))**2)

        IF(FIE(I).GT.0.0 .AND. PSIE(I).GE.0.0) THEN 	    
          TETA(K)=180.D0*DATAN(PSIE(I)/FIE(I))/PI             ! 1. KV.

        ELSEIF(FIE(I).GT.0.0 .AND. PSIE(I).LT.0.0) THEN 	    
          TETA(K)=180.D0*DATAN(PSIE(I)/FIE(I))/PI      ! 4. KV.

        ELSEIF(FIE(I).LT.0.0 .AND. PSIE(I).GT.0.0) THEN
          TETA(K)=180.D0+180.D0*DATAN(PSIE(I)/FIE(I))/PI      ! 2. KV.

        ELSEIF(FIE(I).LT.0.0 .AND. PSIE(I).LT.0.0) THEN
          TETA(K)=-180.D0+180.D0*DATAN(PSIE(I)/FIE(I))/PI    ! 3. KV.	 

        ELSEIF(DABS(FIE(I)).LT.1.E-9 .AND. PSIE(I).GT.0.0) THEN
          TETA(K)=90.D0                                       ! +y OSA

        ELSEIF(DABS(FIE(I)).LT.1.E-9 .AND. PSIE(I).LT.0.0) THEN
          TETA(K)=-90.D0    ! 270.D0                          ! -y OSA

        ENDIF
	  TETA(K)=TETA(K)*PI/180.D0
        PSI1(1,K)=DSQRT(ROO(K))*DSIN(0.5*TETA(K))
        PSI1(2,K)=DSQRT(ROO(K))*DCOS(0.5*TETA(K))
        PSI1(3,K)=DSQRT(ROO(K))*DSIN(0.5*TETA(K))*DSIN(TETA(K))
        PSI1(4,K)=DSQRT(ROO(K))*DCOS(0.5*TETA(K))*DSIN(TETA(K))

        PSI1R(1,K)=0.5*DSIN(0.5*TETA(K))/DSQRT(ROO(K))
        PSI1R(2,K)=0.5*DCOS(0.5*TETA(K))/DSQRT(ROO(K))
        PSI1R(3,K)=0.5*DSIN(0.5*TETA(K))*DSIN(TETA(K))/DSQRT(ROO(K))
        PSI1R(4,K)=0.5*DCOS(0.5*TETA(K))*DSIN(TETA(K))/DSQRT(ROO(K))       
	
        PSI1T(1,K)=0.5*DSQRT(ROO(K))*DCOS(0.5*TETA(K))
        PSI1T(2,K)=-0.5*DSQRT(ROO(K))*DSIN(0.5*TETA(K))
        PSI1T(3,K)=DSQRT(ROO(K))*(0.5*DCOS(0.5*TETA(K))*DSIN(TETA(K))+
     1	         DSIN(0.5*TETA(K))*DCOS(TETA(K)))
        PSI1T(4,K)=DSQRT(ROO(K))*(-0.5*DSIN(0.5*TETA(K))*DSIN(TETA(K))+
     1           DCOS(0.5*TETA(K))*DCOS(TETA(K)))

        PSI1X1(1,K)=DCOS(TETA(K))*PSI1R(1,K)-
     1              DSIN(TETA(K))*PSI1T(1,K)/ROO(K)	  
        PSI1X1(2,K)=DCOS(TETA(K))*PSI1R(2,K)-
     1              DSIN(TETA(K))*PSI1T(2,K)/ROO(K)
        PSI1X1(3,K)=DCOS(TETA(K))*PSI1R(3,K)-
     1              DSIN(TETA(K))*PSI1T(3,K)/ROO(K)
        PSI1X1(4,K)=DCOS(TETA(K))*PSI1R(4,K)-
     1              DSIN(TETA(K))*PSI1T(4,K)/ROO(K)
	
        PSI1X2(1,K)=DSIN(TETA(K))*PSI1R(1,K)+  
     1              DCOS(TETA(K))*PSI1T(1,K)/ROO(K)	  
        PSI1X2(2,K)=DSIN(TETA(K))*PSI1R(2,K)+
     1              DCOS(TETA(K))*PSI1T(2,K)/ROO(K)
        PSI1X2(3,K)=DSIN(TETA(K))*PSI1R(3,K)+
     1              DCOS(TETA(K))*PSI1T(3,K)/ROO(K)
        PSI1X2(4,K)=DSIN(TETA(K))*PSI1R(4,K)+
     1              DCOS(TETA(K))*PSI1T(4,K)/ROO(K)
	
	  !AKO JE ZAOKRENUTA PRSLINA


        PSI1X(1,K)=PSI1X1(1,K)*VX1(1)+PSI1X2(1,K)*VX1(2)
        PSI1X(2,K)=PSI1X1(2,K)*VX1(1)+PSI1X2(2,K)*VX1(2)
        PSI1X(3,K)=PSI1X1(3,K)*VX1(1)+PSI1X2(3,K)*VX1(2)
        PSI1X(4,K)=PSI1X1(4,K)*VX1(1)+PSI1X2(4,K)*VX1(2)
	
        PSI1Y(1,K)=PSI1X2(1,K)*VX2(1)+PSI1X2(1,K)*VX2(2)
        PSI1Y(2,K)=PSI1X2(2,K)*VX2(1)+PSI1X2(2,K)*VX2(2)
        PSI1Y(3,K)=PSI1X2(3,K)*VX2(1)+PSI1X2(3,K)*VX2(2)
        PSI1Y(4,K)=PSI1X2(4,K)*VX2(1)+PSI1X2(4,K)*VX2(2)

      ENDDO

c      DO I=1,KTIP2
c        WRITE(6,*) I,NIZ2(I),ROO(I),TETA(I)*180.0/PI
c      ENDDO
c	PAUSE


      K1=KTIP1
      K2=KTIP2
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE JXFEM(ID1,ID2,K1,K2,JEDN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      
C     PROGRAM ZA DODAVANJE JEDNACINA ZA XFEM      
C
      COMMON /GREECE/ JEQ
      DIMENSION ID1(2,*),ID2(8,*)
C      
      JEQ=JEDN
      JED=JEDN
      DO I=1,K1
        DO J=1,2
          JED=JED+1
          ID1(J,I)=JED
        ENDDO
      ENDDO    
C
      DO I=1,K2
        DO J=1,8
          JED=JED+1
          ID2(J,I)=JED
        ENDDO
      ENDDO
      JEDN=JED    
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE LMXFEM(LM,NODTIP,ID1,ID2,NEL,
     1                  NLM,NE,NCVE,ND,ND1,ND2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      
C     PROGRAM ZA DODAVANJE JEDNACINA ZA XFEM      
C
      include 'pakxfem.inc'
      
      DIMENSION NEL(NE,*),LM(*),NODTIP(*),ID1(2,*),ID2(8,*)
C      
C
      DO 40 NC=1,NCVE
         IF(NEL(NLM,NC).EQ.0) GO TO 40              
         JJ=NODTIP(NEL(NLM,NC))
         IF(JJ.LT.0) THEN
            JJ=-JJ
            DO J=1,2
               LM(ND+ND1+J)=ID1(J,JJ)
            ENDDO
            ND1=ND1+2
         ENDIF
   40 CONTINUE                                                        
C
      DO 50 NC=1,NCVE
         IF(NEL(NLM,NC).EQ.0) GO TO 50              
         JJ=NODTIP(NEL(NLM,NC))
         IF(JJ.GT.0) THEN
            DO J=1,8
               LM(ND+ND1+ND2+J)=ID2(J,JJ)
            ENDDO
            ND2=ND2+8
         ENDIF
   50 CONTINUE       

c      CALL IWRR(LM,ND+ND1+ND2,'LM  ')
c      CALL IWRR(ID1,KTIP1*2,'ID1 ')
c      CALL IWRR(ID2,KTIP2*8,'ID2 ')
      RETURN
      END
C======================================================================
C
C=======================================================================
      SUBROUTINE BETBEX(CORD,H,BLT,NOP,X1,NODTIP,HZNAK,KELEM,JG,NP,r,s)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     MATRICE BLT (BE TRANSPONOVANO)
CE     MATRIX  BLT ( B TRANSPOSED )
C
      include 'pakxfem.inc'
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
	COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2

      DIMENSION H(NCVE,*),BLT(KK,*),CORD(NP,*),NOP(NE,*)
      DIMENSION NODTIP(*),HZNAK(NE,*),KELEM(*),XC(4),YC(4)
	DIMENSION PSG(4),PSGX(4),PSGY(4)
      COMMON /CDEBUG/ IDEBUG


      IF(IDEBUG.GT.0) PRINT *, ' BETBED'
      JJ=-1
      DO 40 I=1,NCVE
         JJ=JJ+2
         IF(NOP(NLM,I).EQ.0) GO TO 40
         BLT(1,JJ)=XJ(1,1)*H(I,2)+XJ(1,2)*H(I,3)
         BLT(3,JJ)=XJ(2,1)*H(I,2)+XJ(2,2)*H(I,3)
         BLT(2,JJ+1)=BLT(3,JJ)
         BLT(3,JJ+1)=BLT(1,JJ)
         IF(IETYP.NE.1) GO TO 40
         IF(X1.GT.1.0D-8) GO TO 30
         BLT(4,JJ)=BLT(1,JJ)
         GO TO 40
   30    BLT(4,JJ)=H(I,1)/X1
   40 CONTINUE

C     XFEM
      ND=JJ+1

      ND1=0 !LOKALNE PROMENLJIVE
      ND2=0

	DO LL=1,NCVE
	  XC(LL)=CORD(NOP(NLM,LL),1)
        YC(LL)=CORD(NOP(NLM,LL),2)
      ENDDO


c      if(nlm.eq.9.and.jg.eq.1) then
c      write(3,*) 'xi,eta',r,s
c      write(3,*) 'rn',(h(ij,1),ij=1,4)
c      write(3,*) 'dnr',(h(ij,2),ij=1,4)
c      write(3,*) 'dns',(h(ij,3),ij=1,4)
c      write(3,*) 'rj',((rj(ij,ji),ij=1,2),ji=1,2)
c      write(3,*) 'xj',((xj(ij,ji),ij=1,2),ji=1,2)
c      write(3,*) 'x',(xc(ij),ij=1,4)
c      write(3,*) 'y',(yc(ij),ij=1,4)
c      endif

      DO 45 I=1,NCVE
         IF(NOP(NLM,I).EQ.0) GO TO 45              
         II=NODTIP(NOP(NLM,I))
         IF(II.LT.0) THEN
	      CALL GAUSHH(H,XC,YC,HZNAK,KELEM,NLM,JG,HH)
            II=-II
            JJ=ND+ND1+1         
            BLT(1,JJ)=(XJ(1,1)*H(I,2)+XJ(1,2)*H(I,3))*HH 
            BLT(3,JJ)=(XJ(2,1)*H(I,2)+XJ(2,2)*H(I,3))*HH 
            BLT(2,JJ+1)=BLT(3,JJ)
            BLT(3,JJ+1)=BLT(1,JJ)
            ND1=ND1+2
c            WRITE(3,*) ' NOP(NLM,I),HH,II,JJ',NOP(NLM,I),HH,II,JJ
         ENDIF
   45 CONTINUE                                                        
C
      DO 50 I=1,NCVE
         IF(NOP(NLM,I).EQ.0) GO TO 50              
         II=NODTIP(NOP(NLM,I))
         IF(II.GT.0) THEN
            DO J=1,4
	         CALL GAUSNT(H,XC,YC,PSG,PSGX,PSGY)
               JJ=ND+ND1+ND2+1         
               BLT(1,JJ)=(XJ(1,1)*H(I,2)+XJ(1,2)*H(I,3))*PSG(J)+ 
     *                   H(I,1)*PSGX(J)          
               BLT(3,JJ)=(XJ(2,1)*H(I,2)+XJ(2,2)*H(I,3))*PSG(J)+
     *                   H(I,1)*PSGY(J)              
               BLT(2,JJ+1)=BLT(3,JJ)
               BLT(3,JJ+1)=BLT(1,JJ)
               ND2=ND2+2
c      if(nlm.eq.9.and.jg.eq.1) then
c               WRITE(3,*) ' NOP(NLM,I),II,JJ',NOP(NLM,I),II,JJ
c               WRITE(3,*) ' blt1,3',BLT(1,JJ),BLT(3,JJ)
c      endif
            ENDDO   
         ENDIF
   50 CONTINUE                                


      RETURN
      END


C=======================================================================
C=======================================================================
	SUBROUTINE GAUSHH(H,XC,YC,HZNAK,KELEM,NLM,JG,HH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
	COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
	DIMENSION H(NCVE,*)
	DIMENSION XC(4),YC(4),VL1(2),VL2(2),VEKTOR(2)
	DIMENSION HZNAK(NE,*),KELEM(*)


	INDIK=0	 !INDIKATOR DA JE ELEMENT PRESECEN SA PRSLINOM
	IF(JG.EQ.1) THEN
	   IKP=0
	   IKN=0
      ENDIF
	DO LL=1,NEC
	  IF(KELEM(LL).EQ.NLM) THEN
	     INDIK=1
        ENDIF
	ENDDO

!     DEFINISANJE LOKALNE ORIJENTACIJE PRSLINE U ELEMENTU
      IPRES=0
      IF(INDIK.GT.0) THEN
	  DO KK=1,NES
	   DO LL=1,4
	     
		 IF(KELC(KK,LL) .EQ.NLM) THEN
	       IPRES=IPRES+1
	       IF(IPRES.EQ.1) THEN
	         XPRES1=XEC(KK)
	         YPRES1=YEC(KK)
              ELSEIF(IPRES.EQ.2) THEN
	         XPRES2=XEC(KK)
	         YPRES2=YEC(KK)
              ENDIF
           ENDIF

         ENDDO
        ENDDO	   		  
	ENDIF

	IF(INDIK.GT.0 .AND.     !IPRES.EQ.1 .AND. 
	1  DSQRT((XPRES1-XPRES2)**2+
     1        (YPRES1-YPRES2)**2).LT.1.E-9) THEN
	   INDIK=0


!     ILI STAVITI   OVO SU KOSINUSI PRAVCA OD POCETKA DO VRHA
C       VL1(1)=VX1(1)
C       VL1(2)=VX1(2)
C	  INDIK=1
C       VL2(1)=-1.D0*VL1(2)	   
C       VL2(2)=VL1(1)
C	 GOTO 100
      ENDIF

	IF(INDIK.GT.0) THEN
C    ! KOSINUSI PRAVCA PRSLINE U ELEMENTU NLM

       VL1(1)=XPRES2-XPRES1
       VL1(2)=YPRES2-YPRES1
       VV=dsqrt(VL1(1)*VL1(1)+VL1(2)*VL1(2))
       VL1(1)=VL1(1)/VV
       VL1(2)=VL1(2)/VV


	 ORIJ=VL1(1)*VX1(1)+VL1(2)*VX1(2) !SKALARNI PROIZVOD
	 IF(ORIJ .LT.0.D0) THEN
	   POMX=XPRES2
	   POMY=YPRES2
	   XPRES2=XPRES1
	   YPRES2=YPRES1
	   XPRES1=POMX
	   YPRES1=POMY
         VL1(1)=XPRES2-XPRES1
         VL1(2)=YPRES2-YPRES1
         VV=dsqrt(VL1(1)*VL1(1)+VL1(2)*VL1(2))
         VL1(1)=VL1(1)/VV
         VL1(2)=VL1(2)/VV
	 ENDIF

       VL2(1)=-1.D0*VL1(2)	   !OVO NIJE DIRANO
       VL2(2)=VL1(1)

	ENDIF !ZAVRSENA PETLJA PO (INDIK.GT.0)

  100   IF(INDIK.GT.0) THEN
	  GX=0.D0
	  GY=0.D0
	  DO LL=1,NCVE
	!   GLOBALNE KOORDINATE GAUSOVE TACKE GX,GY
	    GX=GX+H(LL,1)*XC(LL)
	    GY=GY+H(LL,1)*YC(LL)
	  ENDDO
	    VEKTOR(1)=GX-XPRES2       !Xp(NSEGX+1)
		VEKTOR(2)=GY-YPRES2             !Yp(NSEGX+1)
		ORIJ=VEKTOR(1)*VL2(1)+VEKTOR(2)*VL2(2) !SKALARNI PROIZVOD


		IF(ORIJ.LT.0.D0) THEN
		   HH=-1.D0
	       IKN=IKN+1

		ELSEIF(ORIJ.GT.0.D0)THEN
		   HH=1.D0
	       IKP=IKP+1
		ELSE
	       IF(IKN.GT.IKP) THEN
		       HH=1.D0  		         
	           IKP=IKP+1
             ELSE  
		       HH=-1.D0 
			   IKN=IKN+1             
             ENDIF
             
	    ENDIF
	    

      ENDIF
	
CC	U ELEMENTIMA KOJI NISU PRESECENI SA PRSLINOM H ZNAK U GAUSOVIM   
CC	TACKAMA SE ODREDJUJE INTERPOLACIJOM H ZNAKOVA U CVOROVIMA
	IF(INDIK.EQ.0) THEN
        HH1=0.D0
        DO I=1,NCVE
	   	HH1=HH1+H(I,1)*HZNAK(NLM,I)
	  ENDDO
	  ZN=1.D0
	  VRED=DSIGN(ZN,HH1)
	  HH=VRED
C	  IF(HH.GT.0.D0) THEN
C	     IKP=IKP+1
C        ELSEIF(HH.LT.0.D0) THEN
C	     IKN=IKN+1
C        ENDIF 
	ENDIF


c     KONTROLNA STAMPA   
c      IF(JG.EQ.36) THEN
c      IF(NLM.EQ.703.OR.NLM.EQ.677.OR.NLM.EQ.729.OR.NLM.EQ.755.OR.
c     1   NLM.EQ.781.OR.NLM.EQ.807 ) THEN
c	WRITE(6,*) 'NLM',NLM
c	WRITE(6,*) 'IKP,IKN ',IKP,IKN
c	PAUSE
c	ENDIF         
c      ENDIF

	RETURN
	END
C=======================================================================
C=======================================================================
	SUBROUTINE GAUSNT(H,XC,YC,PSG,PSGX,PSGY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
	COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
	DIMENSION H(NCVE,*)
	DIMENSION XC(4),YC(4)
	DIMENSION PSG(4),PSGX(4),PSGY(4)

	GX=0.D0
	GY=0.D0
	DO LL=1,NCVE
	! GLOBALNE KOORDINATE GAUSOVE TACKE GX,GY
	  GX=GX+H(LL,1)*XC(LL)
	  GY=GY+H(LL,1)*YC(LL)
	ENDDO

C      Tangenta i normala poslednjeg segmenta
      t=dsqrt((Xp(NSEGX+1)-Xp(NSEGX))**2+(Yp(NSEGX+1)-Yp(NSEGX))**2)
      tx=(Xp(NSEGX+1)-Xp(NSEGX))/t
      ty=(Yp(NSEGX+1)-Yp(NSEGX))/t
      bx=-ty
      by=tx
      FIG=(GX-Xp(NSEGX+1))*tx+(GY-Yp(NSEGX+1))*ty
      PSIG=(GX-Xp(NSEGX+1))*bx+(GY-Yp(NSEGX+1))*by

      RR=DSQRT((PSIG)**2+(FIG)**2)

        IF(FIG.GT.0.0 .AND. PSIG.GE.0.0) THEN 	    
          TEG=180.D0*DATAN(PSIG/FIG)/PI             ! 1. KV.

        ELSEIF(FIG.GT.0.0 .AND. PSIG.LT.0.0) THEN 	    
          TEG=180.D0*DATAN(PSIG/FIG)/PI            ! 4. KV.

        ELSEIF(FIG.LT.0.0 .AND. PSIG.GT.0.0) THEN
          TEG=180.D0+180.D0*DATAN(PSIG/FIG)/PI     ! 2. KV.

        ELSEIF(FIG.LT.0.0 .AND. PSIG.LT.0.0) THEN
          TEG=-180.D0+180.D0*DATAN(PSIG/FIG)/PI    ! 3. KV.	 

        ELSEIF(DABS(FIG).LT.1.E-9 .AND. PSIG.GT.0.0) THEN
          TEG=90.D0                                ! +y OSA

        ELSEIF(DABS(FIG).LT.1.E-9 .AND. PSIG.LT.0.0) THEN
          TEG=-90.D0                               ! -y OSA

        ENDIF

	TEG=TEG*PI/180.D0
      PSG(1)=DSQRT(RR)*DSIN(0.5*TEG)
      PSG(2)=DSQRT(RR)*DCOS(0.5*TEG)
      PSG(3)=DSQRT(RR)*DSIN(0.5*TEG)*DSIN(TEG)
      PSG(4)=DSQRT(RR)*DCOS(0.5*TEG)*DSIN(TEG)

      PSG1R=0.5*DSIN(0.5*TEG)/DSQRT(RR)
      PSG2R=0.5*DCOS(0.5*TEG)/DSQRT(RR)
      PSG3R=0.5*DSIN(0.5*TEG)*DSIN(TEG)/DSQRT(RR)
      PSG4R=0.5*DCOS(0.5*TEG)*DSIN(TEG)/DSQRT(RR)       
	

      PSG1T=0.5*DSQRT(RR)*DCOS(0.5*TEG)
      PSG2T=-0.5*DSQRT(RR)*DSIN(0.5*TEG)
      PSG3T= DSQRT(RR)*(0.5*DCOS(0.5*TEG)*DSIN(TEG)+
     1                      DSIN(0.5*TEG)*DCOS(TEG))

      PSG4T=DSQRT(RR)*(-0.5*DSIN(0.5*TEG)*DSIN(TEG)+
     1                      DCOS(0.5*TEG)*DCOS(TEG))



      PSG1X1=DCOS(TEG)*PSG1R-DSIN(TEG)*PSG1T/RR	  
      PSG2X1=DCOS(TEG)*PSG2R-DSIN(TEG)*PSG2T/RR
      PSG3X1=DCOS(TEG)*PSG3R-DSIN(TEG)*PSG3T/RR
      PSG4X1=DCOS(TEG)*PSG4R-DSIN(TEG)*PSG4T/RR
	
      PSG1X2=DSIN(TEG)*PSG1R+DCOS(TEG)*PSG1T/RR	  
      PSG2X2=DSIN(TEG)*PSG2R+DCOS(TEG)*PSG2T/RR
      PSG3X2=DSIN(TEG)*PSG3R+DCOS(TEG)*PSG3T/RR
      PSG4X2=DSIN(TEG)*PSG4R+DCOS(TEG)*PSG4T/RR
	
	  !AKO JE ZAOKRENUTA PRSLINA
	  PSGX(1)=PSG1X1*VX1(1)+PSG1X2*VX1(2)
        PSGX(2)=PSG2X1*VX1(1)+PSG2X2*VX1(2)
        PSGX(3)=PSG3X1*VX1(1)+PSG3X2*VX1(2)
        PSGX(4)=PSG4X1*VX1(1)+PSG4X2*VX1(2)
	
        PSGY(1)=PSG1X1*VX2(1)+PSG1X2*VX2(2)
        PSGY(2)=PSG2X1*VX2(1)+PSG2X2*VX2(2)
        PSGY(3)=PSG3X1*VX2(1)+PSG3X2*VX2(2)
        PSGY(4)=PSG4X1*VX2(1)+PSG4X2*VX2(2)

	RETURN 
	END



C=======================================================================
C=======================================================================
      SUBROUTINE EKST_EDG(NUM,FI,PSI,KELEM,X,Y,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'pakxfem.inc'
      
      DIMENSION NUM(NE,*),X(*),Y(*),FI(NE,*),PSI(NE,*),KELEM(*)

cc      write(6,*)'START EKST_EDG.FOR'


C     EKSTENZIJA PRSLINE I REDEFINISANJE LS-funkcija NA OSNOVU NOVOG POLSEDNJEG SEGMENTA 
   10 FORMAT(F15.10)
C     PRIRASTAJ PRSLINE SE DAJE U POLARNOM SISTEMU VEZANOM ZA VRH PRSLINE      
      
c      write(6,*) 'PRIRASTAJ PRSLINE Dea'
      read(5,10)Dea
C     UGAO PRIRASTAJA PRSLINE SE KRECE U GRANICAMA OD -70 DO +70      
c      write(6,*) 'UGAO PRIRASTAJA PRSLINE TE U STEPENIMA'
      read(5,10)TET
c     Ogranicenja za ugao prirastaja         
      if(TET.GT.70.) then
        TET=70.
      endif   
      if(TET.LT.-70.) then
        TET=-70.
      endif  
c     Prevodjenje ugla prirastaja u radijane                
      TET=TET*PI/180.
C     Tangenta i normala postojeceg poslednjeg segmenta
      t=Dsqrt( (Xp(NSEGX+1)-Xp(NSEGX))**2 +(Yp(NSEGX+1)-Yp(NSEGX))**2 )
      tx=(Xp(NSEGX+1)-Xp(NSEGX))/t
      ty=(Yp(NSEGX+1)-Yp(NSEGX))/t
      bx=-ty
      by=tx        
c      WRITE(6,10)TET
c      pause
      Fx1=Dea*Dcos(TET)
      Fy1=Dea*Dsin(TET)
      Fx=tx*Fx1+ty*Fy1
      Fy=bx*Fx1+by*Fy1                                        
      NSEGX=NSEGX+1
c      WRITE(6,*)'NOVI BROJ SEGMENATA JE'
c      WRITE(6,*)NSEGX
      Xp(NSEGX+1)=Xp(NSEGX)+Fx
      Yp(NSEGX+1)=Yp(NSEGX)+Fy
cc     do i=NSEGX,NSEGX+1
cc        write(6,20)i,Xp(i),i,Yp(i)
c      enddo
c   20 format(/' Xp(',I3,')=',F15.10,' ,Yp(',I3,')=',F15.10)         
C     Tangenta i normala NOVOG poslednjeg segmenta
      t=Dsqrt( (Xp(NSEGX+1)-Xp(NSEGX))**2 +(Yp(NSEGX+1)-Yp(NSEGX))**2 )
      tx=(Xp(NSEGX+1)-Xp(NSEGX))/t
      ty=(Yp(NSEGX+1)-Yp(NSEGX))/t
      bx=-ty
      by=tx  
c     Tangenta i normala predposlednjeg segmenta      
      t1=Dsqrt( (Xp(NSEGX)-Xp(NSEGX-1))**2 +(Yp(NSEGX)-Yp(NSEGX-1))**2 )
      tx1=(Xp(NSEGX)-Xp(NSEGX-1))/t1
      ty1=(Yp(NSEGX)-Yp(NSEGX-1))/t1
      bx1=-ty1
      by1=tx1

C     Petlja po elementima         
      DO i=1,NE
       do m=1,4
        FIstart=(X(NUM(i,m))-Xp(NSEGX))*tx+(Y(NUM(i,m))-Yp(NSEGX))*ty 
        FIend=(X(NUM(i,m))-Xp(NSEGX+1))*tx+(Y(NUM(i,m))-Yp(NSEGX+1))*ty
        FIend1=(X(NUM(i,m))-Xp(NSEGX))*tx1+(Y(NUM(i,m))-Yp(NSEGX))*ty1
        PSIend=(X(NUM(i,m))-Xp(NSEGX+1))*bx+(Y(NUM(i,m))-Yp(NSEGX+1))*by
C       ****************************************************************
c       OSTAJANJE CVORA U MEDJU-ZONI
c       ================================================================
C       a) Uslov ako se cvor nalazi na kraju predposlednjeg segmenta
        IF(DABS(FIend1).LT.1.E-9) THEN
           FI(i,m)=FIend
        ENDIF
C       b) Uslov ako je tup ugao izmedju predposlednjeg i poslednjeg segmenta
        IF(FIend1.GT.0.0 .AND. FIstart.LT.0.0) THEN
           FI(i,m)=FIend
           PSI(i,m)=Dsign(Dsqrt((X(NUM(i,m))-Xp(NSEGX))**2+
     *                          (Y(NUM(i,m))-Yp(NSEGX))**2),PSIend)
        ENDIF
c       ****************************************************************
c       USLOV ZA NOVI SEGMENT   
        IF(FIstart.GE.0.0) THEN                               
C         ==============================================================
C         Uslov ako se cvor nalazi na startnoj liniji NOVOG poslednjeg segmenta
C         a) ostar ugao izmedju predposlednjeg i polsednjeg segmenta
          if(DABS(FIstart).LT.1.E-9 .AND. FIend1.GT.0.0) then
             FI(i,m)=FIend
             PSI(i,m)=PSIend
          endif   
C         b) tup ugao izmedju predposlednjeg i polsednjeg segmenta
          if(DABS(FIstart).LT.1.E-9 .AND. FIend1.LT.0.0 .AND.
     *       Dabs(PSI(i,m)).GE.Dabs(PSIend)) then   
             FI(i,m)=FIend
             PSI(i,m)=PSIend
          endif 
C         ==============================================================
C         Uslov ako se cvor nalaci unutar zone NOVOG poslednjeg segmenta         
          if(FIstart.GT.0.0) then
             FI(i,m)=FIend
             PSI(i,m)=PSIend
          endif
        ENDIF  
C     *******************************************************************************************        
       enddo
      ENDDO                                                                                    
C     ===========================================================================================
C     ===========================================================================================
C     REDEFINISANJE ELEMENATA PRESECENIH PRSLINOM                       
C     ===========================================================================================      
C     Blok za stampanje LS-funkcija: FI i PSI po cvorovima elemenata      
c      write(6,*)' FI i PSI po cvorovima elemenata'
c      do i=1,NE
c       do m=1,4
c         write(6,100) i,m,FI(i,m),i,m,PSI(i,m) 
c         pause
c       enddo
c      enddo
c  100 format(' FI(',I3,','I1,')=',F15.10,', PSI(',I3,',',I1,')=',F15.10)     
C     ===========================================================================================  
C     ===========================================================================================
C     Identifikacija elemenata koji su preseceni sa prslinom
c      write(6,*)' ELEMENTI PRESECENI SA PRSLINOM'
      NEC=0
      do i=1,NE
        PSIMIN=PSI(i,1)
        PSIMAX=PSI(i,1)
        FIMIN=FI(i,1)
        FIMMAX=FI(i,1)
        do k=1,4
          IF(PSIMAX.LE.PSI(i,k)) PSIMAX=PSI(i,k)
          IF(PSIMIN.GT.PSI(i,k)) PSIMIN=PSI(i,k)
          IF(FIMAX.LE.FI(i,k))   FIMAX=FI(i,k)
          IF(FIMIN.GT.FI(i,k))   FIMIN=FI(i,k)
        enddo 
        IF(PSIMAX*PSIMIN.LE.0.0.AND.FIMAX*FIMIN.LE.0.0) THEN
          NEC=NEC+1
          KELEM(NEC)=i
c          write(6,*) NEC,KELEM(NEC)
c          pause
        ENDIF
      enddo              
C     ==================================================================           
      return                       
      end
C=======================================================================
C
C=======================================================================
c      SUBROUTINE KOR_POM_X(NODTIP,NOP,NLM,NE,NCVE,RTDTL,RGOCA,HNOD)
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      include 'pakxfem.inc'
c      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
c     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
c
c      DIMENSION NOP(NE,*),NODTIP(*),HNOD(*)
c      DIMENSION RTDTL(40)
c      DIMENSION RGOCA(40)
      
!     KOREKCIJA POMERANJA U OBOGACENIM CVOROVIMA H FUNKCIJOM

c      K11=0
	     
c      DO I=1,NCVE
c         RGOCA(I+I-1)=RTDTL(I+I-1)
c         RGOCA(I+I)=RTDTL(I+I)			 
c         IF(NODTIP(NOP(NLM,I)).LT.0 ) THEN
	              
C            UH=RTDTL(9+K11)
C            RTDTL(I+I-1)=RTDTL(I+K)+UH*HNOD(NEL(NLM,I)) !X-pravac
C            UH=RTDTL(10+K11)
C            RTDTL(I+I)=RTDTL(I+1+K)+UH*HNOD(NEL(NLM,I)) !Y-pravac
c            RGOCA(9+K11)=RTDTL(9+K11)*HNOD(NOP(NLM,I))
c            RGOCA(10+K11)=RTDTL(10+K11)*HNOD(NOP(NLM,I))
c            K11=K11+2
c
c         ENDIF
c      ENDDO	!ZAVRSENA PETLJA PO KLASICNIM STEPENIMA SLOBODE U CVORU


c      RETURN
c      END  
