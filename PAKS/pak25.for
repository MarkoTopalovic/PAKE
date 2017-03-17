C=======================================================================
C
C        DINAMIKA 2/D ELEMENT
C
C   SUBROUTINE K09MAS
C              READM
C              SISTTM
C              ELTM
C              JACGMU
C              DIMAS
C              PODMAS
C
C=======================================================================
      SUBROUTINE K09MAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C    GLAVNI PROGRAM ZA POZIVANJE PROGRAMA ZA RACUNANJE MATRICA ELEMENATA
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' K09MAS'
C
      LAU=LMAX
      CALL READM(A(LAU))
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LAE=LMAX
      IF(IETYP.LE.2) THEN
        NCVE2=NCVE*2
        MXAE = NCVE2+(7*NCVE+NCVE2*(NCVE2+1)/2)*IDVA
      ELSE
        NCVE3=NCVE*3
        MXAE = NCVE3+(7*NCVE+15*NCVE*NCVE+NCVE3*(NCVE3+1)/2)*IDVA
      ENDIF
      LMAX = LAE + MXAE
      IF(LMAX.LT.MTOT) GO TO 70
      WRITE(IZLAZ,2009) LMAX,MTOT
      STOP
C
C     FORMIRANJE MATRICE KRUTOSTI ELEMENATA I PAKOVANJE U SISTEM
C
   70 CALL SISTTM(A(LAE),A(LAU))
C
      RETURN
 2009 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA MATRICE ELEMENATA'
     1/' POTREBNA DIMENZIJA, LMAX=',I10/
     2' RASPOLOZIVA DIMENZIJA, MTOT=',I10)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READM(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
C
      include 'paka.inc'
      
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /DUPLAP/ IDVA
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /GAUSVR/ LTEMGT,LCORGT,ICORGT
      COMMON /CVSILE/ NSILA,LESILA
      COMMON /MIXEDM/ MIXED,IOPGS(6),NDS
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /CRACKS/ CONTE,SINTE,FK123(10,3),NODCR(10,14),NCRACK,LQST,
     1                LNERING,LMIE,LPSI,LQ,N100,IRING,NSEG,MAXRIN,MAXSEG
     1                ,MAXNOD,LXL,LYL,LZL,LSIF1,LXGG,LYGG,LZGG,LNNOD
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION AU(*)
      REAL AU
C
C     POZIVANJE PROGRAMA ZA ULAZNE PODATKE .
      IF(IDEBUG.GT.0) PRINT *, ' READM'
C
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     1IETYP,NGAUSX,NGAUSY,BETA,NCVE,ITERME,ICORGT,LCEL,LELC,NMA,NMI,
     1MXAU,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,LLMEL,IPODT,ISHEAR,
     1ND,NGS12,MSLOJ,MXS,MSET,NSILA,LESILA,LNNOD,
     1LNSLOJ,LMATSL,LDSLOJ,LBBET,INCOTX,INCOTY,INCOTZ,
     1LBET0,(CPP(J,1),J=1,3),IBB0,LALFE,LHAEM,LHINV,LGEEK,IALFA,
     1INDBTH,INDDTH,LTBTH,LTDTH,NDS,ILEDE,NLD,ICPM1
      LSTAZA(2)=LMAX8+1
      CALL READDD(AU,MXAU/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAU+MXAU
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SISTTM(AE,AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA(KCAL=1)
C     RACUNANJE NAPONA (KCAL=2)
C
      include 'paka.inc'
      
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /DUPLAP/ IDVA
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /CDEBUG/ IDEBUG
C
      common /ielmkt/ inmkt
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1       IOPGL(6),KOSI,NDIN,ITEST
C
      DIMENSION AE(*),AU(*)
      REAL AE,AU
C
C     REPERI U VEKTORU ELEMENATA AE
C
      IF(IDEBUG.GT.0) PRINT *, ' SISTTM'
C
      if(nkrt.gt.0) then
          if(inmkt.gt.0) go to 20
          ncv1=3
          llm=1
          lske=llm+ncv1
          mxae1=lske+(ncv1*(ncv1+1)/2)*idva
            do 10 i=1,nkrt
              call eltmkt(ae(lske),ae(llm),a(lncgl),ncv1,a(lid),
     1                    a(lxmass),np,nkrt,izlaz,i)

              call spakuj(a(lsk),a(lmaxa),ae(lske),ae(llm),ncv1)
c             call stakr(a(lsk),a(lmaxa))
   10       continue
   20    continue
      endif 
      NCVE2=NCVE*2
      LHE=1
      LBET=LHE+3*NCVE*IDVA
      LSKE=LBET+NCVE2*2*IDVA
      IF(IETYP.LE.2) THEN
        NCVEN=NCVE2
        NWE=NCVE2*(NCVE2+1)/2
        LLM=LSKE+NWE*IDVA
        MXAE1=LLM+NCVE2-1
        LSKEP=LLM
        LSKEP1=LSKEP
      ELSE
        NCVEN=3*NCVE
        NWE=NCVEN*(NCVEN+1)/2
        LSKEP=LSKE+NWE*IDVA
        LSKEP1=LSKEP+NCVE2*NCVEN*IDVA
        LLM=LSKEP1+NCVEN*NCVEN*IDVA
        MXAE1=LLM+NCVEN-1
      ENDIF
      IZBR=7*NCVE+NWE
      IF(IETYP.EQ.3) IZBR= 7*NCVE+15*NCVE*NCVE+NWE
C
C     OSNOVNA PETLJA PO ELEMENTIMA
C
      DO 100 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 100
C
        CALL CLEAR(AE,IZBR)
        KORD=LCORD
        IF(IATYP.EQ.3) KORD=LCORUL
      IF(IMASS.EQ.1) GO TO 150
      CALL DIMAS(A(LSK),A(LMAXA),AE(LLM),AU(LNEL),AU(LNMAT),AU(LTHID),
     1           AU(LLMEL),AE(LHE),A(KORD),A(LGUSM),NCVE2)
      GO TO 100
C
C     RACUNANJE KONZISTENTNE MATRICE MASA
C
  150 CALL ELTM(AE(LSKE),AE(LLM),AU(LNEL),AU(LNMAT),
     1AU(LTHID),AE(LHE),AE(LBET),A(KORD),AU(LIPGC),AU(LLMEL),A(LGUSM),
     2AE(LSKEP),AE(LSKEP1),NCVEN,NCVE2)
C
C     PAKOVANJE MATRICA ELEMENATA U MATRICU SISTEMA
C
      CALL SPAKUJ(A(LSK),A(LMAXA),AE(LSKE),AE(LLM),NCVEN)
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ELTM(SKE,LM,NEL,NMAT,THID,HE,BET,CORD,
     1IPGC,LMEL,GUSM,SKP,SKP1,NCVEN,NCVE2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C                                                                       
C     INTEGRACIJA PO POVRSINI MATRICA ELEMENATA                         
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION SKE(*),LM(*),NEL(NE,*),NMAT(*),                         
     1THID(*),CORD(NP,*),HE(NCVE,*),BET(NCVE2,*),IPGC(*),
     1LMEL(NCVEN,*),GUSM(50,*),SKP(NCVEN,*),SKP1(NCVEN,*)
      DIMENSION XG(55),WGT(55),NREF(11),XGG(15)                          
      DIMENSION TTE(2,3),CORDL(2,9),A12(3),A13(3),EN(3),Y(3)
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
      DATA XGG/            0.D0,              -1.D0,               1.D0,
     1                    -1.D0,               0.D0,               1.D0,
     2                    -1.D0,-.333333333333333D0, .333333333333333D0,
     3                     1.D0,              -1.D0,              -.5D0,
     4                     0.D0,               .5D0,               1.D0/
C
C     FORMIRANJE VEKTORA LM                                             
C                                                                       
      IF(IDEBUG.GT.0) PRINT *, ' ELTM'
C
      DO 10 NC=1,NCVEN                                                  
      LM(NC)=LMEL(NC,NLM)                                               
   10 CONTINUE                                                          
C
C    FORMIRANJE MATRICE TRANSFORMACIJE
C
      IF(IETYP.EQ.3) THEN
        I1=NEL(NLM,1)
        I2=NEL(NLM,2)
        I3=NEL(NLM,3)
        DO 620 I=1,3
          A12(I) = CORD(I2,I) - CORD(I1,I)
  620     A13(I) = CORD(I3,I) - CORD(I1,I)
C     VEKTOR NORMALE
        EN(1) = A12(2)*A13(3) - A12(3)*A13(2)
        EN(2) = A12(3)*A13(1) - A12(1)*A13(3)
        EN(3) = A12(1)*A13(2) - A12(2)*A13(1)
C     VEKTOR Y
        Y(1) = EN(2)*A12(3) - EN(3)*A12(2)
        Y(2) = EN(3)*A12(1) - EN(1)*A12(3)
        Y(3) = EN(1)*A12(2) - EN(2)*A12(1)
        YI = DSQRT(Y(1)*Y(1)+Y(2)*Y(2)+Y(3)*Y(3))
        A12I = DSQRT(A12(1)*A12(1)+A12(2)*A12(2)+A12(3)*A12(3))
        DO 640 I=1,3
          TTE(1,I) = A12(I)/A12I
  640     TTE(2,I) = Y(I)/YI
C
C     FORMIRANJE LOKALNIH KOORDINATA
C
        I1=NEL(NLM,1)
        CORDL(1,1)=0.0D0
        CORDL(2,1)=0.0D0
        DO 650 I=2,NCVE
          IN=NEL(NLM,I)
          DO 660 J=1,2
            CORDL(J,I)=0.0D0
            DO 670 K =1,3
  670       CORDL(J,I) = CORDL(J,I) + TTE(J,K)*(CORD(IN,K) - CORD(I1,K))
  660     CONTINUE
  650   CONTINUE
      ELSE
        DO 600 I=1,NCVE
          II = NEL(NLM,I)
          DO 610 J=1,2
  610     CORDL(J,I)=CORD(II,J)
  600   CONTINUE
      ENDIF
C                                                                       
C     PETLJA PO GAUSOVIM TACKAMA                                        
C                                                                       
      IPGCT=IPGC(NLM)
      MAT=NMAT(NLM)
      GUST=GUSM(NMODM,MAT)
      THICK=THID(NLM)                                                   
C                                                                       
      DO 500 NGR=1,NGAUSX                                               
      JGR=NREF(NGAUSX)+NGR                                              
      R=XG(JGR)                                                         
      IF(IPGCT.EQ.1) R=XGG(JGR)                                         
      WR=WGT(JGR)                                                       
C                                                                       
      DO 50 NGS=1,NGAUSY                                                
      JGS=NREF(NGAUSY)+NGS                                              
      S=XG(JGS)                                                         
      IF(IPGCT.EQ.1) S=XGG(JGS)                                         
      WS=WGT(JGS)                                                       
      WRWS=WR*WS                                                        
C                                                                       
C     JAKOBIJAN U TACKI                                                 
C                                                                       
      CALL JACTEL(NEL,CORDL,HE,R,S,0)                                    
      IF(IETYP.EQ.1)                                                    
     1CALL JACGMU(NEL,CORDL,HE,X1)                                       
      WD=WRWS*DET                                                       
      CONAX=WD*THICK*GUST                                               
      IF(IETYP.EQ.1) CONAX=X1*WD*GUST                                   
      NCV2=NCVE*2-1
      DO 150 I=1,NCV2,2
      NCI=(I+1)/2
      BET(I,1)=HE(NCI,1)
      BET(I,2)=0.0D0
      BET(I+1,1)=0.0D0
      BET(I+1,2)=HE(NCI,1)
  150 CONTINUE
C                                                                       
C     INTEGRACIJA SKE                                                   
C
      IF(IETYP.EQ.3) THEN
        DO 450 I=1,NCVE*2
          DO 450 J=I,NCVE*2
            XX=0.0
            DO 460 K=1,2
  460       XX=XX+BET(I,K)*BET(J,K)
          SKP1(I,J)=SKP1(I,J)+XX*CONAX
          SKP1(J,I)=SKP1(I,J)
  450   CONTINUE
      ELSE
      IJ=0                                                              
      DO 430 I=1,NCVE*2                                                 
      DO 430 J=I,NCVE*2                                                 
      IJ=IJ+1                                                           
      IF(LM(I).EQ.0.OR.LM(J).EQ.0) GO TO 430                            
      XX=0.0D00
      DO 440 K=1,2
      XX=XX+BET(I,K)*BET(J,K)                                           
  440 CONTINUE
      SKE(IJ)=SKE(IJ)+XX*CONAX                                          
  430 CONTINUE
      ENDIF
   50 CONTINUE                                                          
  500 CONTINUE
      IF(IETYP.EQ.3) THEN
        DO 800 N=1,NCVE
          DO 800 I=1,3
            II=3*(N-1)+I
            DO 800 J=1,NCVE*2
              SKP(II,J)=0.0D0
              DO 800 K=1,2
                MM=2*(N-1)+K
                SKP(II,J)=SKP(II,J)+TTE(K,I)*SKP1(J,MM)
  800   CONTINUE 
        DO 850 N=1,NCVE
          DO 850 I=1,3
            II=3*(N-1)+I
            DO 850 J1=1,NCVE
              DO 850 J2=1,3
                J=3*(J1-1)+J2
                SKP1(II,J)=0.0D0
                DO 850 K=1,2
                  MM=2*(J1-1)+K
                  SKP1(II,J)=SKP1(II,J)+TTE(K,J2)*SKP(II,MM)
  850   CONTINUE
        IJ=0
        DO 900 I=1,NCVEN
          DO 900 J=I,NCVEN
            IJ=IJ+1
            SKE(IJ)=SKP1(I,J)
  900   CONTINUE
      ENDIF
      RETURN                                                            
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE JACGMU(NEL,CORDL,H,X1)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C                                                                       
C     JAKOBIJAN I INTERPOLACIJA                                         
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NEL(NE,*),CORDL(2,9),H(NCVE,*)                          
C                                                                       
C     KOORDINATE TACKE                                                  
C                                                                       
      IF(IDEBUG.GT.0) PRINT *, ' JACGMU'
C
      X1=0.0                                                            
      DO 80 I=1,NCVE                                                    
      IF(NEL(NLM,I).EQ.0) GO TO 80                                      
      X1=X1+H(I,1)*CORDL(1,I)                                           
   80 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C======================================================================
      SUBROUTINE DIMAS(SK,MAXA,LM,NEL,NMAT,THID,LMEL,HE,CORD,GUSM,
     1NCVE2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     INTEGRACIJA PO POVRSINI MATRICA ELEMENATA
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION CMC(9),SK(*),MAXA(*),LM(*),NEL(NE,*),NMAT(*),
     1THID(*),CORD(NP,*),HE(NCVE,*),GUSM(50,*),LMEL(NCVE2,*),CORDL(2,9)
      DIMENSION XG(15),WGT(15),NREF(6)
      DATA NREF/0,1,3,6,10,15/                                          
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0/
C
C     FORMIRANJE VEKTORA LM                                             
C                                                                       
      IF(IDEBUG.GT.0) PRINT *, ' DIMAS'
C
      DO 10 NC=1,NCVE2                                                 
      LM(NC)=LMEL(NC,NLM)                                               
   10 CONTINUE                                                          
        DO 600 I=1,NCVE
          II = NEL(NLM,I)
          DO 610 J=1,2
  610     CORDL(J,I)=CORD(II,J)
  600   CONTINUE
C                                                                       
C     PETLJA PO GAUSOVIM TACKAMA                                        
C                                                                       
      AMASA=0.0D0
      MAT=NMAT(NLM)
      GUST=GUSM(NMODM,MAT)
      THICK=THID(NLM)                                                   
C                                                                       
      DO 500 NGR=1,NGAUSX                                               
      JGR=NREF(NGAUSX)+NGR                                              
      R=XG(JGR)                                                         
      WR=WGT(JGR)                                                       
C                                                                       
      DO 50 NGS=1,NGAUSY                                                
      JGS=NREF(NGAUSY)+NGS                                              
      S=XG(JGS)                                                         
      WS=WGT(JGS)                                                       
      WRWS=WR*WS                                                        
C                                                                       
C     JAKOBIJAN U TACKI                                                 
C                                                                       
      CALL JACTEL(NEL,CORDL,HE,R,S,0)                                    
      IF(IETYP.EQ.1)                                                    
     1CALL JACGMU(NEL,CORD,HE,X1)                                       
      WD=WRWS*DET                                                       
      CONAX=WD*THICK*GUST           
      IF(IETYP.EQ.1) CONAX=X1*WD*GUST
C                                    
C     INTEGRACIJA SKE                
C                                    
      AMASA= AMASA+CONAX
   50 CONTINUE                       
  500 CONTINUE                       
      CALL PODMAS(NEL,CMC)
      DO 80 NC=1,NCVE
         IF(NEL(NLM,NC).EQ.0) GO TO 80
         NNC=2*(NC-1)
         DO 75 I=1,2
            NJ=LM(NNC+I)
            IF(NJ.EQ.0) GO TO 75
C            NJM=MAXA(NJ)
C            SK(NJM)=SK(NJM)+AMASA*CMC(NC)
            SK(NJ)=SK(NJ)+AMASA*CMC(NC)
   75    CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE PODMAS(NEL,CMC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     PODELA MASA ELEMENTA PO CVOROVIMA
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NEL(NE,*),CMC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' PODMAS'
C
      IF(NCVE.GT.4) GO TO 100
      CMC(1)=0.25D0
      CMC(2)=CMC(1)
      CMC(3)=CMC(1)
      CMC(4)=CMC(1)
      CMC(5)=0.D0
      CMC(6)=CMC(5)
      CMC(7)=CMC(5)
      CMC(8)=CMC(5)
      CMC(9)=CMC(5)
      GO TO 200
C
C     KOREKCIJA FUNKCIJA KADA JE BROJ CVOROVA VECI OD 4
C
  100 CMC(1)=0.0625D0
      CMC(2)=CMC(1)
      CMC(3)=CMC(1)
      CMC(4)=CMC(1)
      CMC(5)=0.125D0
      CMC(6)=CMC(5)
      CMC(7)=CMC(5)
      CMC(8)=CMC(5)
      CMC(9)=0.25D0
      C112=1.D0/12.D0
      C124=1.D0/24.D0
      C148=1.D0/48.D0
      IF(NEL(NLM,9).NE.0) GO TO 110
      CMC(1)=CMC(1)+C148
      CMC(2)=CMC(2)+C148
      CMC(3)=CMC(3)+C148
      CMC(4)=CMC(4)+C148
      CMC(5)=CMC(5)+C124
      CMC(6)=CMC(6)+C124
      CMC(7)=CMC(7)+C124
      CMC(8)=CMC(8)+C124
      CMC(9)=0.D0
C     PETI CVOR
      IF(NEL(NLM,5).NE.0) GO TO 10
      CMC(1)=CMC(1)+C112
      CMC(2)=CMC(2)+C112
C     SESTI CVOR
   10 IF(NEL(NLM,6).NE.0) GO TO 20
      CMC(2)=CMC(2)+C112
      CMC(3)=CMC(3)+C112
C     SEDMI CVOR
   20 IF(NEL(NLM,7).NE.0) GO TO 30
      CMC(3)=CMC(3)+C112
      CMC(4)=CMC(4)+C112
C     OSMI CVOR
   30 IF(NEL(NLM,8).NE.0) GO TO 40
      CMC(1)=CMC(1)+C112
      CMC(4)=CMC(4)+C112
   40 GO TO 200
C     PETI CVOR
  110 IF(NEL(NLM,5).NE.0) GO TO 11
      CMC(1)=CMC(1)+C124
      CMC(2)=CMC(2)+C124
      CMC(9)=CMC(9)+C124
C     SESTI CVOR
   11 IF(NEL(NLM,6).NE.0) GO TO 21
      CMC(2)=CMC(2)+C124
      CMC(3)=CMC(3)+C124
      CMC(9)=CMC(9)+C124
C     SEDMI CVOR
   21 IF(NEL(NLM,7).NE.0) GO TO 31
      CMC(3)=CMC(3)+C124
      CMC(4)=CMC(4)+C124
      CMC(9)=CMC(9)+C124
C     OSMI CVOR
   31 IF(NEL(NLM,8).NE.0) GO TO 200
      CMC(1)=CMC(1)+C124
      CMC(4)=CMC(4)+C124
      CMC(9)=CMC(9)+C124
  200 CONTINUE
      RETURN
      END

