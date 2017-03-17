      SUBROUTINE J2DX(NUM,NERING,MIE,QST,XL,YL,XGG,YGG,
     1               SIF1AVG,TAU,N45,NDI,LM,LMEL,RTDT,DEF,NNOD,
     1               ID1,ID2,NODTIP,HNOD,MODEL,NGG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'paka.inc'

      include 'pakxfem.inc'
      COMMON /MATIZO/ E,V
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /JINTEG/ SQ(9),COORDL(2,9),COORDG(2,9),NCR,NCS,KR,KE,JGE
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
      COMMON /SRPSKI/ ISRPS
  
C    **************************************************************
C    *   IZRACUNAVANJE FAKTORA INTENZIVNOSTI NAPONA J-INTEGRALOM  *
C    **************************************************************
      DIMENSION LM(ND),LMEL(ND,NE),RTDT(*),
     1          NERING(NCRACK,MAXSEG,MAXRIN),
     1          MIE(NCRACK,MAXSEG,MAXRIN,N100),
     1          QST(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          XL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YL(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),      
     1          XGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),
     1          YGG(NCRACK,MAXSEG,MAXRIN,N100,MAXNOD),NNOD(NE),
     1          SIF1AVG(NCRACK,NSEG,3),
     1          TAU(N45,NGS12,NE,*),DEF(N45,NGS12,NE,*),NUM(NE,*)
      DIMENSION ID1(2,*),ID2(8,*)
      DIMENSION NODTIP(*),HNOD(*),NGG(*)
      DIMENSION VEDI1(10,1,4),SIF1(10,1,4),EDI1(10,1,4,1000)
      DIMENSION VEDI2(10,1,4),SIF2(10,1,4),EDI2(10,1,4,1000)
C     ================================================================
C     ULAZNE VELIZINE
C     NERING(10,1,4)-BROJ ELEMENATA U KR-tom PRSTENU I-tog seg.
C     MIE(10,1,4,100)-GLOBALNA NUM.J-tog EL. U KR-tom PRST. I-tog SEG. 
C     QST(10,1,4,100,8)-FUN.S PO NOD. J-tog EL.KR-tog PRST.,I-tog SEG.  
C     ===============================================================
C     EDI(10,1,4,100)-J INTEGRAL U SVAKOM ELEMENTU PRSTENA-KR
C     VEDI(10,1,4)-UKUPNA VREDNOST J INTEGRALA U PRSTENU KR (KR=1-4)
C     SIF1(10,1,4)-FAKTOR K1 U PRSTENU-KR (KR=1-4) 
C     SIF1AVG(NCRACK,NSEG,1)-OSREDNJENO K1 PO SVIM PRSTENOVIMA SEGMENTA
C     SIF1AVG(NCRACK,NSEG,2)-OSREDNJENA K2 PO SVIM PRSTENOVIMA SEGMENTA
C     ==============================================================
	
      DO 400 NCR=1,NCRACK        !PETLJA PO PRSLINAMA
         WRITE(3,*) ' '
         if(isrps.eq.0) then
            WRITE(6,*) ' PRSLINA',NCR
            WRITE(3,*) ' PRSLINA',NCR
         else
            WRITE(6,*) ' CRACK',NCR
            WRITE(3,*) ' CRACK',NCR
         endif
         FAK1=0.
         FAK2=0.
         DO 300 NCS=1,NSEG	         !PETLJA PO SEGMENTIMA
            if(isrps.eq.0) then
               WRITE(6,*) ' SLOJ',NCS
               WRITE(3,*) ' SLOJ',NCS
            else
               WRITE(6,*) ' LAYER',NCS
               WRITE(3,*) ' LAYER',NCS
            endif
            DO 200 KR=1,IRING	         !PETLJA PO PRSTENOVIMA
C     INICIJ.UKUPNE VREDNOSTI U SVIM PRSTENOVIMA FAKTORA 
C     K1 (SIF-a) PRSLINE-NCR I SEGMENTA-NCS 
C      IF(KR.EQ.1) THEN
C         SIF1AVG(NCR,NCS,1)=0.D0 
C         SIF2AVG(NCR,NCS,1)=0.D0 
C      ENDIF

               VEDI1(NCR,NCS,KR)=0.0   !J-INTEGRAL U KR-tom PRSTENU
               VEDI2(NCR,NCS,KR)=0.0   !J-INTEGRAL U KR-tom PRSTENU
               NEKR=NERING(NCR,NCS,KR) !BROJ ELEM. U KR-tom PRSTENU
               if(isrps.eq.0) then
                  WRITE(6,*) ' PRSTEN',KR,' BROJ CVOROVA',NEKR
                  WRITE(3,*) ' PRSTEN',KR,' BROJ CVOROVA',NEKR
               else
                  WRITE(6,*) ' RING',KR,' NUMBER OF NODES',NEKR
                  WRITE(3,*) ' RING',KR,' NUMBER OF NODES',NEKR
               endif
C      WRITE(6,*) KR,NEKR
C      WRITE(3,*) 'J-INTEGRAL i BROJ ELEM. U KR-tom PRSTENU', KR,NEKR
C      PAUSE

               DO 100 KE=1,NEKR		!PETLJA PO ELEMENTIMA KR-tog PRSTENA
                  JGE=MIE(NCR,NCS,KR,KE) !GLOBALNO NUMERISANI ELEMENT
C         WRITE(6,*) 'KR,KE,JGE',KR,KE,JGE 
                  DO I=1,NNOD(JGE)
                     COORDL(1,I)=XL(NCR,NCS,KR,KE,I)  !LOKALNA KOORD.x1
                     COORDL(2,I)=YL(NCR,NCS,KR,KE,I)  !LOKALNA KOORD.x2
                     COORDG(1,I)=XGG(NCR,NCS,KR,KE,I) !GLOBALNA KOORD.X1
                     COORDG(2,I)=YGG(NCR,NCS,KR,KE,I) !GLOBALNA KOORD.X2
                     SQ(I)=QST(NCR,NCS,KR,KE,I)       !TEZ. FUN. PO NODOVIMA
C            WRITE(6,*)'SQ I',SQ(I),I
                  ENDDO
   
         !DEFINISANJE EDI-INTEGRALA U JGE ELEMENTU
            CALL INTEGX(NUM,EDI1,EDI2,TAU,N45,LM,LMEL,RTDT,NDI,DEF,NNOD,
     1                  ID1,ID2,NODTIP,HNOD,NGG,A(LPLAS1))


         !SABIRANJE J-INTEGRALA PO ELEMENTIMA TEKUCEG KR PRSTENA
                  VEDI1(NCR,NCS,KR)=VEDI1(NCR,NCS,KR)+
     1                              EDI1(NCR,NCS,KR,KE)       

                  VEDI2(NCR,NCS,KR)=VEDI2(NCR,NCS,KR)+
     1                              EDI2(NCR,NCS,KR,KE)


  100          CONTINUE

               F=1.D0 !ZA 2D-SLUCAJ	
               VEDI1(NCR,NCS,KR)=VEDI1(NCR,NCS,KR)/F
               VEDI2(NCR,NCS,KR)=VEDI2(NCR,NCS,KR)/F


      !KADA JE MODEL=2 RADI SE POLA MODELA I J-INTEGRAL TREBA DA SE MNOZI SA 2 
      !KADA JE MODEL=1 RADI SE CEO  MODEL  I J-INTEGRAL TREBA DA SE MNOZI SA 1 
               IF(MODEL.EQ.2) THEN
                  CIN=2.D0
               ELSE
                  CIN=1.D0 
               ENDIF

               IF(IETYP.EQ.0) THEN
                  PMOD=E           !PLANE STRESS
               ELSEIF(IETYP.GT.0) THEN
                  PMOD=E/(1-V*V)   !PLANE STRAIN
               ENDIF
C OVO PROVERI ZASTO JE PONOVO STAVLJENO BEZ USLOVA, ZASTO NIJE OSTAO GORNJI USLOV      
c      CIN=1.D0    !ZA CEO MODEL svi primeri radjeni kao ceo model
               DUM=CIN*VEDI1(NCR,NCS,KR)*PMOD
               SIF1(NCR,NCS,KR)=SQRT(DABS(DUM))

C     SIF1(NCR,NCS,KR)=0.5*SQRT(PMOD)*
C    1                 (SQRT(DABS(VEDI1(NCR,NCS,KR)-VEDI2(NCR,NCS,KR)))+
C    1                  SQRT(DABS(VEDI1(NCR,NCS,KR)+VEDI2(NCR,NCS,KR))))

C      WRITE(6,*) 'KR,SIF1',KR,SIF1(NCR,NCS,KR)
C      WRITE(3,*) 'Broj prstena, KI  ',KR,SIF1(NCR,NCS,KR)


               DUM=CIN*VEDI2(NCR,NCS,KR)*PMOD
               SIF2(NCR,NCS,KR)=SQRT(DABS(DUM))

               SIF2(NCR,NCS,KR)=0.5D0*SQRT(PMOD)*
     1                 (SQRT(DABS(VEDI1(NCR,NCS,KR)-VEDI2(NCR,NCS,KR)))-
     1                  SQRT(DABS(VEDI1(NCR,NCS,KR)+VEDI2(NCR,NCS,KR))))

C      IF(VX1(2).LT.0.000 .AND. SIF2(NCR,NCS,KR).LT.0.000 ) THEN 
C         SIF2(NCR,NCS,KR)=-SIF2(NCR,NCS,KR)
C      ENDIF

C      WRITE(6,*) 'KR,SIF2',KR,SIF2(NCR,NCS,KR)
C      WRITE(3,*) 'Broj prstena, KII ',KR,SIF2(NCR,NCS,KR)
C      WRITE(6,*) 'KR,SIF1,SIF2',KR,SIF1(NCR,NCS,KR),SIF2(NCR,NCS,KR)
C      WRITE(3,*) 'KR,SIF1,SIF2',KR,SIF1(NCR,NCS,KR),SIF2(NCR,NCS,KR)
C      WRITE(3,*) ' '

C     UKUPNA VREDNOST SIF-a PO SVIM PRSTENOVIMA NCS SEGMENTA
               SIF1AVG(NCR,NCS,1)=SIF1AVG(NCR,NCS,1)+SIF1(NCR,NCS,KR)
c      IF(KR.LT.4) THEN
               SIF1AVG(NCR,NCS,2)=SIF1AVG(NCR,NCS,2)+SIF2(NCR,NCS,KR)
c      ENDIF

            if(isrps.eq.0) then
             WRITE(6,*) ' FAKTOR INTENZITETA NAPONA K1',SIF1(NCR,NCS,KR)
             WRITE(6,*) ' FAKTOR INTENZITETA NAPONA K2',SIF2(NCR,NCS,KR)
             WRITE(3,*) ' FAKTOR INTENZITETA NAPONA K1',SIF1(NCR,NCS,KR)
             WRITE(3,*) ' FAKTOR INTENZITETA NAPONA K2',SIF2(NCR,NCS,KR)
            else
               WRITE(6,*) ' STRESS INTENSITY FACTOR K1',SIF1(NCR,NCS,KR)
               WRITE(6,*) ' STRESS INTENSITY FACTOR K2',SIF2(NCR,NCS,KR)
               WRITE(3,*) ' STRESS INTENSITY FACTOR K1',SIF1(NCR,NCS,KR)
               WRITE(3,*) ' STRESS INTENSITY FACTOR K2',SIF2(NCR,NCS,KR)
            endif
  200       CONTINUE
	
C     OSREDNJENA VREDNOST SIF-a U NCS SEGMENTU
            SIF1AVG(NCR,NCS,1)=SIF1AVG(NCR,NCS,1)/IRING
            FAK1=FAK1+SIF1AVG(NCR,NCS,1)
c      WRITE(6,*)NCR,NCS,SIF1AVG(NCR,NCS,1)
C      WRITE(3,*)'Prslina br., NCS,  KI',NCR,NCS,SIF1AVG(NCR,NCS,1)

            SIF1AVG(NCR,NCS,2)=SIF1AVG(NCR,NCS,2)/IRING
            FAK2=FAK2+SIF1AVG(NCR,NCS,2)
c      WRITE(6,*)NCR,NCS,SIF1AVG(NCR,NCS,2)
C      WRITE(3,*)'Prslina br., NCS, KII',NCR,NCS,SIF1AVG(NCR,NCS,2)


  300    CONTINUE
         WRITE(3,*) ' '
         if(isrps.eq.0) then
            WRITE(6,*) ' PROSECNI FAKTOR INTENZITETA NAPONA K1',FAK1
            WRITE(6,*) ' PROSECNI FAKTOR INTENZITETA NAPONA K2',FAK2
            WRITE(3,*) ' PROSECNI FAKTOR INTENZITETA NAPONA K1',FAK1
            WRITE(3,*) ' PROSECNI FAKTOR INTENZITETA NAPONA K2',FAK2
         else
            WRITE(6,*) ' AVERAGE STRESS INTENSITY FACTOR K1',FAK1
            WRITE(6,*) ' AVERAGE STRESS INTENSITY FACTOR K2',FAK2
            WRITE(3,*) ' AVERAGE STRESS INTENSITY FACTOR K1',FAK1
            WRITE(3,*) ' AVERAGE STRESS INTENSITY FACTOR K2',FAK2
         endif
         WRITE(3,*) ' '
  400 CONTINUE

c      CALL FATCRACK 
C      CALL FATSEC  
c     DIREKTNO ODREDJIVANJE FAKTORA K
CCC        CALL DIREKTK(NUM,LM,LMEL,RTDT,NDI,ID2,NODTIP)
      RETURN
      end

C============================================================  

C============================================================
      SUBROUTINE INTEGX(NUM,EDI1,EDI2,TAU,N45,LM,LMEL,RTDT,NDI,DEF,NNOD,
     1                  ID1,ID2,NODTIP,HNOD,NGG,PLAST)
	        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'pakxfem.inc'

      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /MATIZO/ E,V
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /JINTEG/ SQ(9),COORDL(2,9),COORDG(2,9),NCR,NCS,KR,KE,JGE
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2

C     ********************************************************
C    *       ODREDJIVANJE J-INTEGRALA KE-tog ELEMENTA		  *
C     ********************************************************
      DIMENSION EDI1(10,1,4,100),EDI2(10,1,4,100),
     1          LM(ND),LMEL(ND,NE),RTDT(*)
      DIMENSION DER(NCVE,3),DSHAP(NCVE,2),DQRST(100,2),QF(100),NNOD(NE)
C     ,THID(*)

      DIMENSION TERMEDI1(100,2),TERMEDI2(100,2),WEDI(4),
     1          SIGMA(2,2),DQX(100,2),DUX1(100,2),DUX2(100,2),
     1          DISP(NCVE,2),STRS(3),STRN(3),
     1          TAU(N45,NGS12,NE,*),DEF(N45,NGS12,NE,*),PLAST(*),
     1          NUM(NE,*)
      DIMENSION ID1(2,*),ID2(8,*)
      DIMENSION NODTIP(*),HNOD(*),NGG(*)
      DIMENSION XG(55),WGT(55),NREF(11)
C
      DATA NREF/0,1,3,6,10,15,21,28,36,45,55/
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0,
     5       .171324492379170D0, .360761573048139D0, .467913934572691D0,
     6       .467913934572691D0, .360761573048139D0, .171324492379170D0,
     7       .129484966168870D0, .279705391489277D0, .381830050505119D0,
     8       .417959183673469D0, .381830050505119D0, .279705391489277D0,
     9       .129484966168870D0, .101228536290376D0, .222381034453374D0,
     9       .313706645877887D0, .362683783378362D0, .362683783378362D0,
     1       .313706645877887D0, .222381034453374D0, .101228536290376D0,
     2       .081274388361574D0, .180648160694857D0, .260610696402935D0,
     3       .312347077040003D0, .330239355001260D0, .312347077040003D0,
     4       .260610696402935D0, .180648160694857D0, .081274388361574D0,
     5       .066671344308688D0, .149451349150581D0, .219086362515982D0,
     6       .269266719309996D0, .295524224714753D0, .295524224714753D0,
     7       .269266719309996D0, .219086362515982D0, .149451349150581D0,
     8       .066671344308688D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0,
     5      -.932469514203152D0,-.661209386466265D0,-.238619186083197D0,
     6       .238619186083197D0, .661209386466265D0, .932469514203152D0,
     7      -.949107912342759D0,-.741531185599394D0,-.405845151377397D0,
     8                     0.D0, .405845151377397D0, .741531185599394D0,
     9       .949107912342759D0,-.960289856497536D0,-.796666477413627D0,
     9      -.525532409916329D0,-.183434642495650D0, .183434642495650D0,
     1       .525532409916329D0, .796666477413627D0, .960289856497536D0,
     2      -.968160239507626D0,-.836031107326636D0,-.613371432700590D0,
     3      -.324253423403809D0,               0.D0, .324253423403809D0,
     4       .613371432700590D0, .836031107326636D0, .968160239507626D0,
     5      -.973906528517172D0,-.865063366688985D0,-.679409568299024D0,
     6      -.433395394129247D0,-.148874338981631D0, .148874338981631D0,
     7       .433395394129247D0, .679409568299024D0, .865063366688985D0,
     8       .973906528517172D0/

C
      NGAUS=NGG(JGE)       !ZA 2D NGAUS=4 ili 36,NDI=2 
      IF(NGAUS.EQ.36) THEN
         NGAUSX=6
         NGAUSY=6
      ELSEIF(NGAUS.EQ.100) THEN
	   NGAUSX=10
         NGAUSY=10 
      ELSE
         NGAUSX=2
         NGAUSY=2
      ENDIF
		 
C      write(6,*) 'NGG(JGE),N45,NGS12 NGAUS'
C      write(6,*)  JGE,NGG(JGE),N45,NGS12,NGAUS
C      PAUSE
      !INICIJALIZACIJA J-INTEGRALA ZA TEKUCI ELEMENT
      EDI1(NCR,NCS,KR,KE)=0.D0  
      EDI2(NCR,NCS,KR,KE)=0.D0  

      DO K=1,NGAUS            !PETLJA PO GAUSOVIM TACKAMA
         DQX(K,1)=0.D0        !INIC.dS/dx1,dS/dx2 PO GAUS.TAC.
         DQX(K,2)=0.D0 
         DUX1(K,1)=0.D0       !INIC.dU1/dx1,dU2/dx1 PO GAUS.TAC.
         DUX1(K,2)=0.D0  
         DUX2(K,1)=0.D0       !INIC.dU1/dx2,dU2/dx2 PO GAUS.TAC.
         DUX2(K,2)=0.D0
         QF(K)=0.D0           !INIC.VREDNOSTI TEZINSKE FUNKCIJE
         DQRST(K,1)=0.D0      !INIC.POCETNE VREDNOSTI dQ/dr
         DQRST(K,2)=0.D0      !INIC.POCETNE VREDNOSTI dQ/ds
      ENDDO

      DO I=1,NNOD(JGE)
         DO J=1,2         !NDI
            DISP(I,J)=0.D0
         ENDDO
      ENDDO

      !UCITAVANJE POMERANJA (DISP) ZA DATI ELEMENT JGE
      CALL IJEDN1(LM,LMEL(1,JGE),ND)
      K=0
      DO I=1,4           !4 BROJ CVOROVA PO ELEMENTU
         DO J=1,NDI    !NDI=2
            K=K+1
            IF(LM(K).GT.0) THEN
               L=LM(K)
               DISP(I,J)=RTDT(L)
            ENDIF

            IF(NODTIP(NUM(JGE,I)).LT.0 ) THEN
               JJ=-NODTIP(NUM(JGE,I))
               L1=ID1(J,JJ)
               DISP(I,J)=DISP(I,J)+RTDT(L1)*HNOD(NUM(JGE,I))
            ENDIF   

            IF(NODTIP(NUM(JGE,I)).GT.0 ) THEN
               JJ=NODTIP(NUM(JGE,I))
               IF(J.EQ.1) THEN
                  KPOC=1
               ELSE
                  KPOC=2
               ENDIF	
               M=0 					   			 
               DO KK1=KPOC,8,2
                  M=M+1
                  L1=ID2(KK1,JJ) 
                  DISP(I,J)=DISP(I,J)+RTDT(L1)*PSI1(M,JJ)
               ENDDO
            ENDIF			 

         ENDDO
      ENDDO

C      DO I=1,NNOD(JGE)
C         DO J=1,NDI
C            K=K+1
C            IF(LM(K).GT.0) THEN
C               L=LM(K)
C               DISP(I,J)=RTDT(L)
C            ENDIF
C         ENDDO
C      ENDDO



      NLM=JGE
C      THI=THID(NLM)
C     GLAVNA PETLJA PO GAUSOVIM TACKAMA
      LI=0
      NPR56=MODPR2( NMODM )
      NGXYZ=NGS12*MXS
      NPROS=(NLM-1)*NGXYZ-1
      DO 1199 NGR=1,NGAUSX
         JGR=NREF(NGAUSX)+NGR
         R=XG(JGR)
         WR=WGT(JGR)
C
      DO 1199 NGS=1,NGAUSY
         JGR=NREF(NGAUSY)+NGS
         S=XG(JGR)
         WS=WGT(JGR)
C
         LI=LI+1

	 !UCITATI STRS(3),STRN(3) ZA DATU GAUSOVU TACKU-LI
         IF(NMODM.GT.4) THEN
            IBTC=LI
            LL=1+(NPROS+IBTC)*NPR56
            CALL JEDNA1(STRS,PLAST(LL),3)
            CALL JEDNA1(STRN,PLAST(LL+4),3)
         ELSE
            DO I=1,3
               STRS(I)=TAU(I,LI,JGE,1) !(JGE-GLOBALNI ELEM,LI-GAUS.TAC U TOM ELEM.)
               STRN(I)=DEF(I,LI,JGE,1)
            ENDDO
         ENDIF

C        KONSTRUKCIJA TENZORA NAPONA
         SIGMA(1,1)=STRS(1)
         SIGMA(2,2)=STRS(2)
         SIGMA(1,2)=STRS(3)
         SIGMA(2,1)=STRS(3)

C        INTERPOLACIJSKE FUNKCIJE I JAKOBIJAN U TACKI R, S, T
         CALL JACTEL(NUM,COORDG,DER,R,S,0)
  
         DO I=1,NNOD(JGE)
            DSHAP(I,1)=XJ(1,1)*DER(I,2)+XJ(1,2)*DER(I,3)
            DSHAP(I,2)=XJ(2,1)*DER(I,2)+XJ(2,2)*DER(I,3)
         ENDDO
  
         DO I=1,NNOD(JGE)  !PETLJA PO NODOVIMA
            DUX1(LI,1)=DUX1(LI,1)+DSHAP(I,1)*DISP(I,1)	 !DEF.DUX1<dU1/dx1>  
            DUX1(LI,2)=DUX1(LI,2)+DSHAP(I,1)*DISP(I,2)	 !DEF.DUX1<dU2/dx1>
            DUX2(LI,1)=DUX2(LI,1)+DSHAP(I,2)*DISP(I,1)	 !DEF.DUX2<dU1/dx2>
            DUX2(LI,2)=DUX2(LI,2)+DSHAP(I,2)*DISP(I,2)	 !DEF.DUX2<dU2/dx2>
         ENDDO

         DO I=1,NNOD(JGE)
            QF(LI)=QF(LI)+DER(I,1)*SQ(I)
            DQRST(LI,1)=DQRST(LI,1)+DER(I,2)*SQ(I)   
            DQRST(LI,2)=DQRST(LI,2)+DER(I,3)*SQ(I)   
         ENDDO
       
C        DEFINISANJE VREDNOSTI DQX (dQ/dx1,dQ/dx2)
         DQX(LI,1)=XJ(1,1)*DQRST(LI,1)+XJ(1,2)*DQRST(LI,2)
         DQX(LI,2)=XJ(2,1)*DQRST(LI,1)+XJ(2,2)*DQRST(LI,2)

 
C        PRVI SABIRAK EDI-INTEGRALA ZA DATU GAUSOVU TACKU
         TERMEDI1(LI,1)=DUX1(LI,1)*SIGMA(1,1)*DQX(LI,1)+
     1                  DUX1(LI,1)*SIGMA(1,2)*DQX(LI,2)+
     1                  DUX1(LI,2)*SIGMA(2,1)*DQX(LI,1)+
     1	              DUX1(LI,2)*SIGMA(2,2)*DQX(LI,2)

         TERMEDI2(LI,1)=DUX2(LI,1)*SIGMA(1,1)*DQX(LI,1)+
     1                  DUX2(LI,1)*SIGMA(1,2)*DQX(LI,2)+
     1                  DUX2(LI,2)*SIGMA(2,1)*DQX(LI,1)+
     1	              DUX2(LI,2)*SIGMA(2,2)*DQX(LI,2)


C        DEFINISANJE SPECIFICNOG RADA U DATOJ GAUSOVOJ TACKI
        WEDI(LI)=0.5*STRS(1)*STRN(1)+0.5*STRS(2)*STRN(2)+STRS(3)*STRN(3)
 
C        PRORACUN DRUGOG SABIRKA EDI-INTEGRALA
         TERMEDI1(LI,2)=WEDI(LI)*DQX(LI,1)
         TERMEDI2(LI,2)=WEDI(LI)*DQX(LI,2)

C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         Xsr=0.0
c         Ysr=0.0
c         DO I=1,NNOD(JGE)
c            Xsr=Xsr+COORDL(1,I)
c            Ysr=Ysr+COORDL(2,I)
c         ENDDO
c         Xsr=Xsr/NNOD(JGE)
c         Ysr=Ysr/NNOD(JGE)
c         SINTE=Ysr/SQRT(Xsr*Xsr+Ysr*Ysr)
c         CONTE=Xsr/SQRT(Xsr*Xsr+Ysr*Ysr)	  
C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         EDI1GL=(TERMEDI1(LI,1)-TERMEDI1(LI,2))*WR*WS*DET
         EDI2GL=(TERMEDI2(LI,1)-TERMEDI2(LI,2))*WR*WS*DET
C     ZILE PROMENIO
            EDI1LOC=CONTE*EDI1GL+SINTE*EDI2GL
            EDI2LOC=-SINTE*EDI1GL+CONTE*EDI2GL  
C         IF(VX1(2).LT.0.00D0) THEN
C            EDI1LOC=CONTE*EDI1GL-DABS(SINTE)*EDI2GL
C            EDI2LOC=DABS(SINTE)*EDI1GL+CONTE*EDI2GL  
C         ELSE
C            EDI1LOC=CONTE*EDI1GL+DABS(SINTE)*EDI2GL
C            EDI2LOC=-DABS(SINTE)*EDI1GL+CONTE*EDI2GL  
C         ENDIF

C        PRORACUN EDI INTEGRALA LOC ,GL
         EDI1(NCR,NCS,KR,KE)=EDI1(NCR,NCS,KR,KE)+EDI1LOC
         EDI2(NCR,NCS,KR,KE)=EDI2(NCR,NCS,KR,KE)+EDI2LOC

C        *THI

 1199 CONTINUE	 !ZAVRSENA PETLJA PO GAUSOVIM TACKAMA

      RETURN
      END
C============================================================  
C============================================================  
      SUBROUTINE DIREKTK(NUM,LM,LMEL,RTDT,NDI,ID2,NODTIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'pakxfem.inc'
      COMMON /MATIZO/ E,V
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
      DIMENSION LM(ND),LMEL(ND,NE),RTDT(*) 
      DIMENSION ID2(8,*)
      DIMENSION NODTIP(*)
      DIMENSION NUM(NE,*),DIS(4,2),U(2)

C     DEFINISANJE MODULA KLIZANJA 
      Gmod=0.5D0*E/(1+V)
C     DEFINISANJE KOLOSOVE KONSTANTE
      IF(IETYP.EQ.0) THEN
         AKAPA=3.D0-4.D0*V       !PLANE STRESS
      ELSEIF(IETYP.GT.0) THEN
         AKAPA=(3.D0-V)/(1.D0+V) !PLANE STRAIN
      ENDIF
C     DELILAC
      DEL=2.D0*Gmod*DSQRT(2.D0*PI)	

      DO NN=1,NE
         JGE=NN
         !UCITAVANJE POMERANJA (DISP) ZA DATI ELEMENT JGE
         CALL IJEDN1(LM,LMEL(1,JGE),ND)
         K=0
         DO N=1,4  !4 BROJ CVOROVA PO ELEMENTU
            U(1)=0.D0
            U(2)=0.D0
            DO J=1,NDI    !NDI=2
               K=K+1
               IF(LM(K).GT.0) THEN
                  L=LM(K)
                  DIS(N,J)=RTDT(L)
                  U(J)=RTDT(L)
               ENDIF

CC              IF(NODTIP(NUM(JGE,N)).LT.0 ) THEN
CC                 JJ=-NODTIP(NUM(JGE,N))
CC                 L1=ID1(J,JJ)
CC                 DIS(N,J)=DIS(N,J)+RTDT(L1)*HNOD(NUM(JGE,N))
CC              ENDIF   
               IF(NODTIP(NUM(NN,N)).GT.0) THEN
                  JJ=NODTIP(NUM(JGE,N))
                  IF(J.EQ.1) THEN
                     KPOC=1
                  ELSE
                     KPOC=2
                  ENDIF	
                  M=0 					   			 
                  DO KK=KPOC,8,2
                     M=M+1
                     L1=ID2(KK,JJ) 
                     DIS(N,J)=DIS(N,J)+RTDT(L1)*PSI1(M,JJ)
c                     U(J)=U(J)+RTDT(L1)*PSI1(M,JJ)
                  ENDDO
               ENDIF

            ENDDO !ZAVRSENA PETLJA PO J (STEPENIMA SLOBODE)
            IF(NODTIP(NUM(NN,N)).GT.0) THEN
               D11=((AKAPA-1.D0)*PSI1(1,JJ)+PSI1(3,JJ))/DEL
               D12=((AKAPA+1.D0)*PSI1(2,JJ)+PSI1(4,JJ))/DEL
               D21=((AKAPA+1.D0)*PSI1(2,JJ)-PSI1(4,JJ))/DEL
               D22=((-AKAPA+1.D0)*PSI1(1,JJ)+PSI1(3,JJ))/DEL
c               DK1=(dabs(U(1)*D22)-dabs(U(2)*D12))/(D11*D22-D21*D12)
cc              DK2=(dabs(U(2)*D11)-dabs(U(1)*D21))/(D11*D22-D21*D12)
               DK1=dabs(U(1)*D22-U(2)*D12)/dabs(D11*D22-D21*D12)
               DK2=(dabs(U(2)*D11)-dabs(U(1)*D21))/dabs(D11*D22-D21*D12)
c               WRITE(6,*) 'NOD',NUM(JGE,N)
c               WRITE(6,*) 'DK1,DK2'
c               WRITE(6,*)  DK1,DK2
c               PAUSE
            ENDIF	!ZAVRSEN USLOV PO OBOGACENOM CVORU

         ENDDO !ZAVRSENA PETLJA PO -N- CVOROVIMA ELEMENTA
      ENDDO !ZAVRSENA PETLJA PO -NN- ELEMENTIMA MREZE      
      RETURN
      END
C==========================================================
C==========================================================



