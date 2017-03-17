C=======================================================================
C
CE       READ DATA: COORDINATES
CS       UCITAVANJE PODATAKA: KOORDINATE
C
C   SUBROUTINE UCKORD
C              UCKOR
C              ULAZE1
C              ULAZE2
C              ULAZE3
C              ULAZE4
C              ULAZE5
C              VEZACV
C              POLAR
C              TGRAFK
C
C=======================================================================
      SUBROUTINE UCKORD
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .                                                                     
CE.    P R O G R A M                                                    
CE.        TO CALL ROUTINE FOR READING NODAL POINT DATA
CS.    P R O G R A M                                                    
CS.        ZA FORMIRANJE REPERA BROJEVA STAZA I OPSTIH PODATAKA O       
CS.        PODSTRUKTURAMA, I ZA POZIVANJE PROGRAMA KOJI FORMIRA REPERE  
CS.        ZA KOORDINATE I OGRANICENJA CVORNIH TACAKA                   
C .                                                                     
CE.    P O I N T E R S                                                  
CE.        LIPODS- IS POINTER FOR ARRAY NPODS(JPS+1,100). 
CE.                NPODS()- IS BASIC ARRAY WITH GENERAL PROGRAM CONTROL
CE.                         PARAMETERS, POINTERS AND RECORD NUMBERS. 
CE.                JPS    - IS NUMBER OF SUBSTRUCTURES.
CE.                         (JPS=1 - SUBSTRUCTURES ARE NOT USED)
CS.    R E P E R I                                                      
CS.        LIPODS - OPSTI PODACI O PODSTRUKTURAMA I BROJEVI STAZA       
C .
CE.    V A R I A B L E S
CE.       LMAX13- RECORD NUMBER OF DIRECT ACCESS FILE (ZIPODS)
C .                                                                     
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /MAXDUZ/ XL,YL,ZL
C
      XL=0.D0
      YL=0.D0
      ZL=0.D0
CZxxx
      CALL CZINIT
CZxxx
      LMAX13=0
      LIPODS=LMAX
      JP=JPS1*100
      LMAX=LIPODS+JP
      IF(LMAX.GT.MTOT) CALL ERROR(1)
CZxxx
      CALL CZPROV(LIPODS,1,34)
CZxxx
      CALL ICLEAR(A(LIPODS),JP)
CE    CALLING OF ROUTINE FOR READING NODAL POINT DATA
      CALL UCKOR(A(LIPODS))
CZxxx
      CALL ZASTIT
CZxxxx
C 
C     U L A Z N I    P O D A C I    Z A    M E H A N I K U    L O M A
C
C     IXFEM > 0, UCITAVANJE PODATAKA ZA XFEM
c      CALL UCXFEM
C     IGBM > 0, UCITAVANJE PODATAKA ZA MESLESS
      CALL UCGBM
C     NCRACK > 0, UCITAVANJE PODATAKA O PRSLINAMA ZA J INTEGRAL
      CALL UCRACK
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UCKOR(NPODS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO FORM POINTERS AND CALL ROUTINES FOR READING NODAL
CE.        COORDINATES AND LOCAL/GLOBAL CONSTRAINS
CS.    P R O G R A M
CS.        ZA FORMIRANJE REPERA I POZIVANJE PROGRAMA
CS.        ZA UCITAVANJE KOORDINATA I OGRANICENJA CVORNIH TACAKA
C .
CE.    P O I N T E R S
CE.        LCORD - NODAL POINT COORDINATES
CE.        LID   - NODAL POINT CONSTRAINS
CE.        LLJUS - INDICATORS FOR LOCAL/GLOBAL ROTATION  
CS.    R E P E R I
CS.        LCORD - KOORDINATE CVORNIH TACAKA
CS.        LID   - INDIKATORI OGRANICENJA PA BROJEVI JEDNACINA
CE.    V A R I A B L E S
CE.        NP    - TOTAL NUMBER OF NODAL POINTS
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /MASINA/ INDPC,ICRTA
      COMMON /SRPSKI/ ISRPS
      COMMON /KAKO6O/ LLJUS,LKAKO6
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /MIXEDM/ MIXED,IOPGS(6),NDS
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /krut3d/ kt3d,mmpogr
      common /ispakk/ liracv
      COMMON /CRXFEM/ NCXFEM,LNODTIP,LNSSN,LPSIE,LFIE,LHNOD,
     1                LPSIC,LFI,LHZNAK,LNSSE,LKELEM,LID1,LID2
      COMMON /MESLESS/ IGBM,ndif,idif(50),NKI,IKI(10)
      COMMON /NIDEAS/ IDEAS
      COMMON /CDEBUG/ IDEBUG
      COMMON/VERSION/ IVER
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UCKOR'
      if(IVER.EQ.1.AND.np.gt.300) STOP 
     1 'PROGRAM STOP - STUDENTSKA VERZIJA PROGRAMA, BROJ CVOROVA VECI JE
     1 OD 300'
CZxxx
      CALL CZPROV(NP,1,34)
CZxxx
      NP6=NP*6
      IF(MIXED.EQ.1) NP6=NP6*2
      NP3=NP*3
      LNCVP=LMAX
      IF(NCVPR.GT.0) LMAX=LNCVP+NP
      LLJUS=LMAX
      LMAX=LLJUS+NP
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LCORD=LMAX
      LID=LCORD+NP3*IDVA
      LMAX=LID+NP6
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      if(nkrt.gt.0) then
          lnkt=lmax
          llncvz=lnkt+nlnc
          lnbv=llncvz+nlnc
          lncgl=lnbv+nlnc
          lnvez=lncgl+nkrt
          lncvez=lnvez+nkrt
          livez=lncvez+nkrt*mxvez
          lmax=livez+nkrt-1
          CALL DELJIV(LMAX,2,INDL)
          IF(INDL.EQ.0) LMAX=LMAX+1
          lxmass=lmax
          lmax=lxmass+nkrt*3*idva
          if(kt3d.eq.1)  lmax=lxmass+nkrt*6*idva
      endif  
C
CE    INPUT DATA FOR NODAL POINT COORDINATES AND CONSTRAINTS
CS    ULAZNI PODACI ZA KOORDINATE I OGRANICENJA
C
      IF(ITEST.GT.0) THEN
         CALL ISPITA(ACOZ)
         IF(INDFOR.EQ.1)
     1   READ(IULAZ,*) JEDN,NWK
         IF(INDFOR.EQ.2)
     1   READ(ACOZ,1000) JEDN,NWK
         IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) JEDN,NWK
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) JEDN,NWK
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
         ENDIF
         CALL ISPITA(ACOZ)
         CALL READTT(A(LID),NP6)
         RETURN
      ENDIF
C
      IF(ICVEL.EQ.0) THEN
         if(nkrt.gt.0) then
            call ulazkt(A(LCORD),A(LID),A(LLJUS),a(lnkt),a(llncvz),
     +           a(lnbv),a(lncgl),a(lnvez),a(lncvez),a(livez),a(lxmass))
         else
CE          READING OF NODAL POINT COORDINATES AND CONSTRAINTS
            CALL ULAZE1(A(LCORD),A(LID),A(LLJUS))
         endif
         CALL DIMENZ(A(LCORD),NP)
         LCVEL=LMAX
         LELCV=LMAX
         NPA=0
         NPI=0
      ELSE
         LCVEL=LMAX
         LMAX=LCVEL+NP
         IF(LMAX.GT.MTOT) CALL ERROR(1)
         CALL ULAZE2(A(LCORD),A(LID),A(LCVEL),NPA,NPI,NP,A(LLJUS))
         CALL DIMENZ(A(LCORD),NP)
         NPM=NPA-NPI+1
         LELCV=LMAX
         LMAX=LELCV+NPM
         IF(LMAX.GT.MTOT) CALL ERROR(1)
         CALL ICLEAR(A(LELCV),NPM)
         CALL VEZACV(A(LCVEL),A(LELCV),NP,NPI)
         LCVELN=LMAX
         LELCVN=LCVELN+NP
         LMAX=LELCVN+NPM
         CALL DELJIV(LMAX,2,INDL)
         IF(INDL.EQ.0) LMAX=LMAX+1
         LCORDN=LMAX
         LIDN=LCORDN+NP*3*IDVA
         LLJUSN=LIDN+NP*6
         LMAX=LLJUSN+NP
         IF(LMAX.GT.MTOT) CALL ERROR(1)
         CALL ULAZE3(A(LCVEL),A(LCVELN),A(LELCV),A(LELCVN),A(LID),
     1   A(LIDN),A(LCORD),A(LCORDN),NP,NPM,A(LLJUS),A(LLJUSN))
         LMAX=LCVELN
         NP6=NP*7+NPM
      ENDIF
      NPODS(JPS1,1)=LMAX13+1
      CALL WRITDD(A(LCORD),NP3,IPODS,LMAX13,LDUZI)
      NPODS(JPS1,2)=LMAX13+1
      IF(NCVPR.GT.0) THEN
         CALL ICLEAR(A(LNCVP),NP)
         CALL INDST(A(LNCVP),A(LMAX),A(LELCV),A(LCVEL),
     +              ICVEL,NP,NPI,NCVPR)
      ENDIF
C...  INICIJALIZACIJA REPERA I DIMENZIJA
      if(nkrt.gt.0) go to 521
C OVU GLUPOST DOBRO PROVERITI
         MMP=0
C         MMP=1
         NEZAV=1
         NEZA1=NEZAV+1
         LCMPC=LMAX
         LMPC =LMAX
  521 IF(NMPC.GT.0) THEN
         if(nkrt.eq.0) then 
            IND=1
            CALL MPCONS(A(LCMPC),A(LMPC),A(LID),A(LELCV),
     &                  ICVEL,NPI,NP,LMAX,IULAZ,IZLAZ,IND)
            IND=2
            CALL MPCONS(A(LCMPC),A(LMPC),A(LID),A(LELCV),
     &                  ICVEL,NPI,NP,LMAX,IULAZ,IZLAZ,IND)
         else
          lcmpc=lmax
          call deljiv(lcmpc,2,indl)
          if(indl.eq.0) lcmpc=lcmpc+1
          neza1=nezav+1
          lmpc=lcmpc+mmp*nezav*idva
          lmax=lmpc+neza1*mmp
          call deljiv(lmax,2,indl)
          if(indl.eq.0) lmax=lmax+1
          if(kt3d.eq.0) then
           lac=lmax
           nkrt2=nkrt*2
           limpc=lmax+nkrt2*nezav*idva
           ljmpc=limpc+nkrt*nezav
           lmax=ljmpc+nkrt
           call deljiv(lmax,2,indl)
           if(indl.eq.0) lmax=lmax+1
c
c
c
c         formiranje matrice koefic. glavnih cvorova krutih tela
           call formac(a(lcord),a(lnkt),a(llncvz),a(lnbv),a(lncgl),
     1                a(lnvez),a(lncvez),izlaz,np,a(lac),nkrt2,nezav)
c
c
c
c
c        formiranje matrice koeficijenata za zavisna pomeranja
                call accmpc(a(lcord),a(lnkt),a(llncvz),a(lnbv),a(lncgl),
     1                a(lnvez),a(lncvez),izlaz,np,a(lac),nkrt2,
     2                nezav,a(lcmpc),mmp)
c
c
c
c
c         vektor jednacina zavisnih pomeranja u x i y pravcu
               call formim(a(lid),np,a(lnkt),a(llncvz),a(lnbv),a(lncgl),
     1                a(lnvez),a(lncvez),izlaz,a(limpc),a(ljmpc),nezav)
c
c

c
c
c     izracunavanje broja korekcija ukoliko ima ogranicenja zadnjeg krutog tela
               call brkorekc(a(lid),np,nkrt,a(lncgl),ikor)
c
c
          	if(ikor.eq.0) go to 612
             	lijogr=lmax
             	lmax=lijogr+ikor
             	call ogrcvor(a(lid),np,a(lijogr),ikor)
c
c
c        formiranje matrice jednacina zavisnih pomeranja
  612          call formmpc(a(lid),np,a(lnkt),a(lnbv),a(lncgl),a(lnvez),
     1                 a(lncvez),izlaz,a(limpc),a(ljmpc),nezav,a(lmpc),
     2                 neza1,mmp)
c
c
         	if(ikor.eq.0) go to 613
           	do 614 ik=1,ikor
c
c
c             korigovanje matrica id,mpc,cmpc
             	 call korig(a(lid),np,a(lijogr),ik,ikor,a(lmpc),
     1                     a(lcmpc),izlaz)
c
c
  614      continue
  613    continue
c
          	if(mxrac.gt.0) then
             	liracv=lmax+nkrt
             	lmax=liracv+mxvc*mxrac
             	call deljiv(lmax,2,indl)
             	if(indl.eq.0) lmax=lmax+1
c
c
c            matrica racvanja u tackama veze krutih tela za zadate krutosti
             	call formirac(a(lid),a(lncgl),a(lnvez),a(lncvez),izlaz
     1                         ,np,a(livez),a(liracv))
       		endif
           else
             lro=lmax
             lnid=lro+3*idva
             lmax=lnid+3
             call deljiv(lmax,2,indl)
             if(indl.eq.0) lmax=lmax+1
             call cmpc3d(a(lcord),a(lncgl),a(lnvez),a(lncvez),izlaz,
     1                   np,nezav,a(lcmpc),mmp,a(lro))
             call mpc3d(a(lid),np,a(lncgl),a(lnvez),a(lncvez),izlaz,
     1                  nezav,a(lmpc),neza1,mmp,a(lnid))
           endif
        endif
      ENDIF
      IF(IDEAS.LT.7.and.ideas.gt.-1) THEN
         IF(NBLGR.GE.0.AND.NP.GT.0) THEN
            CALL TGRAFK(A(LCORD),A(LCVEL),ICVEL,NP,IGRAF)
            CALL TGRAFB(A(LID),A(LCVEL),ICVEL,NP,IGRAF)
            CALL TGRAFR(A(LID),A(LCVEL),ICVEL,NP,IGRAF)
         ENDIF
      ENDIF
      CALL IWRITD(A(LID),NP6,IPODS,LMAX13,LDUZI)
C
      JIDG=LID
      JCORG=LCORD
      JCVEL=LCVEL
      JELCV=LELCV
      NGA=NPA
      NGI=NPI
      NPG=NP
      JEDNG=JEDN
      IF(JPS.GT.1) THEN
         JMAXA=LMAX
         LMAX=JMAXA+JEDNG+1
         CALL ICLEAR(A(JMAXA),JEDNG+1)
         NPODS(JPS1,3)=NPG
         NPODS(JPS1,4)=NGA
         NPODS(JPS1,5)=NGI
         NPODS(JPS1,6)=JEDNG
         NPODS(JPS1,69)=1
         IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2060)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6060)
         ENDIF
         DO 100 I=1,JPS
            CALL ISPITA(ACOZ)
            IF(INDFOR.EQ.1)
     1      READ(IULAZ,*) JPBR,NGELEM,NP,NPK,NPUP
            IF(INDFOR.EQ.2)
     1      READ(ACOZ,1000) JPBR,NGELEM,NP,NPK,NPUP
            IF(JPBR.LT.1.OR.JPBR.GT.JPS) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) JPBR,JPS
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) JPBR,JPS
               STOP 'PROGRAM STOP: PAK02 - UCKOR'
            ENDIF
            IF(NP.LT.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020) NP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020) NP
               STOP 'PROGRAM STOP: PAK02 - UCKOR'
            ENDIF 
            IF(NPK.LE.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040) NPK
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040) NPK
               STOP 'PROGRAM STOP: PAK02 - UCKOR'
            ENDIF
            IF(NGELEM.EQ.0) NGELEM = 1
CZxxx
            CALL CZPROV(NGELEM,1,34)
CZxxx
            IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
               CALL WBROJK(KARTIC,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2070) JPBR,NGELEM,NP,NPK,NPUP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6070) JPBR,NGELEM,NP,NPK,NPUP
            ENDIF
            NPP=NP+NPK
            CALL DELJIV(LMAX,2,INDL)
            IF(INDL.EQ.0) LMAX=LMAX+1
            LCORD=LMAX
            LID=LCORD+NPP*3*IDVA
            LCVEL=LID+NPP*6
            LCVE=LCVEL+NP
            LMAX=LCVE+NPK
            IF(NP.EQ.0) THEN
               NPA=-1
               NPI=0
               JEDNP=0
               JEDN=0
               NPM=NPA-NPI+1
               LELCV=LMAX
               IF(LMAX.GT.MTOT) CALL ERROR(1)
               GO TO 50
            ENDIF
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            CALL ULAZE2(A(LCORD),A(LID),A(LCVEL),NPA,NPI,NPP,A(LLJUS))
            CALL DIMENZ(A(LCORD),NPP)
            NPM=NPA-NPI+1
            LELCV=LMAX
            LMAX=LELCV+NPM
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            CALL ICLEAR(A(LELCV),NPM)
            CALL VEZACV(A(LCVEL),A(LELCV),NP,NPI)
            LCVELN=LMAX
            LELCVN=LCVELN+NP
            LMAX=LELCVN+NPM
            CALL DELJIV(LMAX,2,INDL)
            IF(INDL.EQ.0) LMAX=LMAX+1
            LCORDN=LMAX
            LIDN=LCORDN+NP*3*IDVA
            LMAX=LIDN+NP*6
            IF(LMAX.GT.MTOT) CALL ERROR(1)
CZxxx
            CALL CZPROV(NP,1,34)
CZxxx
            CALL ULAZE3(A(LCVEL),A(LCVELN),A(LELCV),A(LELCVN),A(LID),
     1            A(LIDN),A(LCORD),A(LCORDN),NPP,NPM,A(LLJUS),A(LLJUSN))
            LMAX=LCVELN
            JEDNP=JEDN
C
   50       CALL ULAZE4(A(LCVE),NAA,NII)
            NMM=NAA-NII+1
            LECV=LMAX
            LMAX=LECV+NMM
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            CALL ICLEAR(A(LECV),NMM)
            CALL VEZACG(A(LCVE),A(LECV),NPK,NII,NP)
            LCVELN=LMAX
            LELCVN=LCVELN+NPK
            LMAX=LELCVN+NMM
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            CALL ULAZE5(A(LCVEL),A(LCVELN),A(LECV),A(LELCVN),A(LID),
     1      A(JIDG),A(LCORD),A(JCORG),A(JELCV),NPP,NMM)
            LMAX=LCVELN
            NPODS(JPBR,1)=LMAX13+1
            CALL WRITDD(A(LCORD),NPP*3,IPODS,LMAX13,LDUZI)
            NPODS(JPBR,2)=LMAX13+1
            NP6=NPP*7+NPM+NMM
            CALL IWRITD(A(LID),NP6,IPODS,LMAX13,LDUZI)
            NPODS(JPBR,3)=NP
            NPODS(JPBR,4)=NPA
            NPODS(JPBR,5)=NPI
            NPODS(JPBR,6)=JEDN
            NPODS(JPBR,7)=NPK
            NPODS(JPBR,8)=NAA
            NPODS(JPBR,9)=NII
            NPODS(JPBR,10)=NPUP
            NPODS(JPBR,11)=NGELEM
            NPODS(JPBR,23)=JEDNP
            IF(NBLGR.GE.0.AND.NP.GT.0) THEN
              CALL TGRAFK(A(LCORD),A(LCVEL),ICVEL,NPP,IGRAF)
              CALL TGRAFB(A(LID),A(LCVEL),ICVEL,NPP,IGRAF)
            ENDIF
            LMG=LMAX
            NP6=JEDN-JEDNP
            LMAX=LMG+NP6
            IF(LMAX.GT.MTOT) CALL ERROR(1)
            CALL LMGMHT(A(JIDG),A(LCVE),A(JELCV),A(LMG))
            NPODS(JPBR,27)=LMAX13+1
            CALL IWRITD(A(LMG),NP6,IPODS,LMAX13,LDUZI)
            JEDN=JEDNG
            CALL VISINE(A(JMAXA),NP6,A(LMG))
            LMAX=LCORD
  100    CONTINUE
         LMAX=LCORD+JEDNG+1
         IF(LMAX.GT.MTOT) CALL ERROR(1)
         LMAXM=0
         JPBR=JPS1
CZxxx
         CALL CZPROV(LCORD,1,34)
CZxxx
         CALL ADRESEP(A(LCORD),A(JMAXA))
         NWG=NWK
         NPODS(JPS1,24)=NWG
         NPODS(JPS1,26)=LMAXM
         NPODS(JPS1,12)=LMAX13+1
         CALL IWRITD(A(JMAXA),JEDNG+1,IPODS,LMAX13,LDUZI)
      ENDIF
      LMAX=JCORG
C..   KOREKCIJA ZA MULTIPOINT CONSTRAINTS
      IF(NMPC.GT.0)THEN
         if(nkrt.gt.0)then
           CALL IJEDN1(A(LMAX),A(lncgl),nkrt) 
           lncgl=lmax
           lmax=lncgl+nkrt
           CALL DELJIV(LMAX,2,INDL)
           IF(INDL.EQ.0) LMAX=LMAX+1
cc
           CALL IJEDN1(A(LMAX),A(lnvez),nkrt) 
           CALL IJEDN1(A(LMAX+nkrt),A(lncvez),mxvez*nkrt) 
           lnvez=lmax
           lncvez=lnvez+nkrt
           lmax=lncvez+mxvez*nkrt
           CALL DELJIV(LMAX,2,INDL)
           IF(INDL.EQ.0) LMAX=LMAX+1
cc
           if(ndin.gt.0) then
               ll3=nkrt*3
               if(kt3d.eq.1) ll3=ll3*2
               CALL JEDNA1(A(LMAX),A(lxmass),ll3)
               lxmass=lmax
               lmax=lxmass+ll3*idva
               CALL DELJIV(LMAX,2,INDL)
               IF(INDL.EQ.0) LMAX=LMAX+1
           endif
         endif
         LL=MMP*NEZAV
         LL1=NMPC*(NEZAV+1)
         CALL JEDNA1(A(LMAX),A(LCMPC),LL)
         CALL IJEDN1(A(LMAX+LL*IDVA),A(LMPC),LL1) 
         LCMPC=LMAX
         LMPC =LCMPC+MMP*NEZAV*IDVA
         LMAX=LMPC+NMPC*(NEZAV+1)
         CALL DELJIV(LMAX,2,INDL)
         IF(INDL.EQ.0) LMAX=LMAX+1
      ENDIF
      MAX13=LMAX13
      RETURN
C
 1000 FORMAT(14I5)
C-----------------------------------------------------------------------
 2000 FORMAT(///
     1'  BROJ JEDNACINA ................................... JEDN =',I5//
     1'  BROJ CLANOVA U MATRICI ............................ NWK =',I5)
 2010 FORMAT(///' V E K T O R   I D')
 2060 FORMAT(///'1'///6X,
     1  'P  O  D  A  C  I     O     P  O  D  S  T  R  U  K  T  U  R  A  
     1M  A'/6X,67('-'))
 2100 FORMAT(//' UCITANI BROJ PODSTRUKTURE JPBR =',I5/' MORA SE NALAZITI
     1 U GRANICAMA OD    1 DO JPS =',I5//' PROGRAM STOP'//)
 2020 FORMAT(/' GRESKA U ULAZNIM PODACIMA :'/' BROJ CVOROVA - NPU =',
     1I5/' A NE SME BITI MANJI OD NULE !'//' PROGRAM STOP'/)
 2040 FORMAT(/' GRESKA U ULAZNIM PODACIMA :'/' BROJ CVOROVA - NPK =',
     1I5/' A MORA BITI VECI OD NULE !'//' PROGRAM STOP'/)
 2070 FORMAT(6X,
     1'P  O  D  S  T  R  U  K  T  U  R  A     B  R  O  J .... JPBR =',I5
     1/6X,66('-')///
     111X,'BROJ GRUPA ELEMENATA .......................... NGELEP =',I5/
     116X,'EQ.0; JEDNA GRUPA ELEMENATA'///
     111X,'BROJ UNUTRASNJIH CVOROVA.......................... NPU =',I5/
     116X,'LT.0; PREKIDA SE IZVRSAVANJE PROGRAMA'///
     111X,'BROJ KONTURNIH CVOROVA ........................... NPK =',I5/
     116X,'LE.0; PREKIDA SE IZVRSAVANJE PROGRAMA'///
     111X,'INDIKATOR RACUNANJA UNUTRASNJIH POMERANJA ....... IPUP =',I5/
     116X,'EQ.0; RACUNAJU SE'/
     116X,'EQ.1; NE RACUNAJU SE'/)
C-----------------------------------------------------------------------
 6000 FORMAT(///
     1'  NUMBER OF EQUATIONS .............................. JEDN =',I5//
     1'  NUMBER OF MATRIX ELEMENTS ......................... NWK  ',I5)
 6010 FORMAT(///' A R R A Y    I D')
 6060 FORMAT(///'1'///6X,
     1  'D  A  T  A    A  B  O  U  T    S  U  B  S  T  R  U  C  T  U  R
     1 E  S'/6X,68('-'))
 6100 FORMAT(//' SUBSTRUCTURE NUMBER, JPBR =',I5/' MUST BE IN BOUNDARY 
     1FROM   1 TO JPS =',I5//' PROGRAM STOP'//)
 6020 FORMAT(/' ERROR IN INPUT DATA:'/
     1' NUMBER OF NODAL POINT - NP =',I5/
     1' MUST NOT BE LEATER THAN ZERO !'//' PROGRAM STOP'/)
 6040 FORMAT(/' ERROR IN INPUT DATA:'/
     1' NUMBER OF NODAL POINT - NPK =',I5/
     1' MUST BE GREATER THAN ZERO !'//' PROGRAM STOP'/)
 6070 FORMAT(6X,
     1'S U B S T R U C T U R E    N U M B E R ............... JPBR =',I5
     1/6X,66('-')///
     111X,'NUMBER OF ELEMENT GROUPS ...................... NGELEP =',I5/
     116X,'EQ.0; DEFAULT SET "1"'///
     111X,'NUMBER OF INTERNAL NODES ......................... NPU =',I5/
     116X,'LT.0; PROGRAM STOP'///
     111X,'NUMBER OF CONTOUR NODES .......................... NPK =',I5/
     116X,'LE.0; PROGRAM STOP'///
     111X,'INDICATOR FOR INTERNAL DISPLACEMENTS ............ IPUP =',I5/
     116X,'EQ.0; CALCULATED'/
     116X,'EQ.1; NOT CALCULATED'/)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE LMGMHT(ID,NCVEL,NELCV,LM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C .......................................................................
C .
CE.   P R O G R A M
CE.       TO FORM VECTOR LM AND HEIGHT COLUMNS FOR SUBSTRUCTURE
CS.   P R O G R A M
CS.       ZA FORMIRANJE VEKTORA LM I VISINA STUBOVA ZA PODSTRUKTURE 
C .
C .......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      DIMENSION ID(NPG,*),NCVEL(*),NELCV(*),LM(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' LMGMHT'
CE    LOOP OVER CONTOUR NODES
CS    PETLJA PO KONTURNIM CVOROVIMA
      IJ=0
      DO 10 NC=1,NPK                             
         NN=NCVEL(NC)
         NN=NN-NGI+1                               
         JJ=NELCV(NN)                             
         DO 20 I=1,IDF                                
            N=ID(JJ,I)
            IF(N.GT.0) THEN
               IJ=IJ+1
               LM(IJ)=N
            ENDIF
   20    CONTINUE
   10 CONTINUE                                   
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZE1(CORD,ID,LJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 CH
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ NODAL POINT COORDINATES AND CONSTRAINTS
CS.    P R O G R A M
CS.        ZA UCITAVANJE KOORDINATA I OGRANICENJA CVORNIH TACAKA
C .
CE.           I=1,NP (NP - TOTAL NUMBER OF NODAL POINTS /3/)
CS.           I=1,NP (NP - UKUPAN BROJ CVORNIH TACAKA)
CE.            /10/ USER MANUAL
C .
CE.            CH - INDICATOR OF POLAR COORDINATES SYSTEM
CS.            CH - KARAKTER INDIKATOR POLARNOG KOORD. SISTEMA
C .
CE.            CORD(I,1) - X COORDINATE
CE.            CORD(I,2) - Y COORDINATE
CE.            CORD(I,3) - Z COORDINATE
CS.            CORD(I,1) - KOORDINATA  X
CS.            CORD(I,2) - KOORDINATA  Y
CS.            CORD(I,3) - KOORDINATA  Z
C .
CE.              ID(I,1) - CONSTRAINT OF TRANSLATION IN DIRECTION X AXIS
CE.              ID(I,2) - CONSTRAINT OF TRANSLATION IN DIRECTION Y AXIS
CE.              ID(I,3) - CONSTRAINT OF TRANSLATION IN DIRECTION Z AXIS
CE.              ID(I,4) - CONSTRAINT OF ROTATION AROUND X AXIS
CE.              ID(I,5) - CONSTRAINT OF ROTATION AROUND Y AXIS
CE.              ID(I,6) - CONSTRAINT OF ROTATION AROUND Z AXIS
CS.              ID(I,1) - OGRANICENJE POMERANJA U PRAVCU - X OSE
CS.              ID(I,2) - OGRANICENJE POMERANJA U PRAVCU - Y OSE
CS.              ID(I,3) - OGRANICENJE POMERANJA U PRAVCU - Z OSE
CS.              ID(I,4) - OGRANICENJE ROTACIJA OKO - X OSE
CS.              ID(I,5) - OGRANICENJE ROTACIJA OKO - Y OSE
CS.              ID(I,6) - OGRANICENJE ROTACIJA OKO - Z OSE
C .
CE.              LJ(I)   - INDICATOR FOR LOCAL/GLOBAL ROTATION
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SRPSKI/ ISRPS
      COMMON /MIXEDM/ MIXED,IOPGS(6),NDS
      DIMENSION CORD(NP,*),ID(NP,*),LJ(*)
      DIMENSION DXYZ(3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZE1'
CZxxx
      CALL CZPROV(NMATM,1,34)
CZxxx
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      ENDIF
C
      N6=6
      IF(MIXED.EQ.1) N6=12
      I = 0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) N,(ID(N,J),J=1,6),CH,(CORD(N,J),J=1,3),KORC,LJ(N)
     1                  ,(ID(N,J),J=7,N6)
      IF(INDFOR.EQ.2.and.ind56.eq.0)
!     1READ(ACOZ,1000) N,(ID(N,J),J=1,6),CH,(CORD(N,J),J=1,3),KORC,LJ(N)
!     1                    ,(ID(N,J),J=7,N6)
     1READ(ACOZ,1000) N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
     1                    ,(ID(N,J),J=7,N6)
      IF(INDFOR.EQ.2.and.ind56.eq.1)
     1READ(ACOZ,4000) N,(ID(N,J),J=1,6),CH,(CORD(N,J),J=1,3),KORC,LJ(N)
     1                    ,(ID(N,J),J=7,N6)
      IF(N.GT.NP) GO TO 110
      IF((NULAZ.EQ.1.OR.NULAZ.EQ.3).and.ind56.eq.0)
     1WRITE(IZLAZ,5010)N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
     1                     ,(ID(N,J),J=7,N6)
      IF((NULAZ.EQ.1.OR.NULAZ.EQ.3).and.ind56.eq.1)
     1WRITE(IZLAZ,8010)N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
     1                     ,(ID(N,J),J=7,N6)
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(N,1),CORD(N,2))
      IF(I.EQ.NP) GO TO 50
      GO TO 5
C
   20 NAUT=1
      N1=N
      KORA=KORC
      GO TO 5
C
CE    AUTOMATIC GENERATE DATA BETWEEN NODAL POINT N1 AND N2
CS    AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
C
   30 N2=N
      CALL DELJIV(N2-N1,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      RKORA = KORA
      RN1N2 = N2-N1
      DD = RKORA/RN1N2
      IF(DD.LT.0.0D00) GO TO 100
      DO 31 J=1,3
   31 DXYZ(J)=DD*(CORD(N2,J)-CORD(N1,J))
      IAUT=(N2-N1)/KORA-1
      N = N1
      DO 34 J=1,IAUT
      DO 35 K=1,3
   35 CORD(N+KORA,K)=CORD(N,K)+DXYZ(K)
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(N,1),CORD(N,2))
      DO 36 K=1,N6
   36 ID(N+KORA,K)=ID(N1,K)
      LJ(N+KORA)=LJ(N1)
      N=N+KORA
   34 CONTINUE
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(N,1),CORD(N,2))
      I=I+IAUT
      IF(I.EQ.NP)THEN
        IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(N2,1),CORD(N2,2))
        GO TO 50
      ENDIF
      NAUT=0
      IF(KORC.EQ.0)THEN
        IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(N2,1),CORD(N2,2))
        GO TO 5
      ENDIF
      KORA=KORC
      NAUT=1
      N1=N2
      GO TO 5
C
CE    UPDATE MATRIX ID
CS    KORIGOVANJE MATRICE ID
C
   50 JEDN=0
      NMPC=0
      DO 60 N=1,NP
      DO 60 J=1,6
C      DO 60 J=1,N6
         IF(J.GT.6) THEN
            IF(IOPGS(J-6).EQ.1) THEN
            ID(N,J)=0
            GO TO 60
            ENDIF
         ELSE
            IF(IOPGL(J).EQ.1) THEN
            ID(N,J)=0
            GO TO 60
            ENDIF
         ENDIF
      IF(ID(N,J)) 64,66,65
   64 NMPC=NMPC+1
      GO TO 60
   65 ID(N,J)=0
      GO TO 60
   66 JEDN=JEDN+1
      ID(N,J)=JEDN
   60 CONTINUE
C
      DO 160 N=1,NP
      DO 160 J=7,N6
         IF(J.GT.6) THEN
            IF(IOPGS(J-6).EQ.1) THEN
            ID(N,J)=0
            GO TO 160
            ENDIF
         ELSE
            IF(IOPGL(J).EQ.1) THEN
            ID(N,J)=0
            GO TO 160
            ENDIF
         ENDIF
      IF(ID(N,J)) 164,166,165
  164 NMPC=NMPC+1
      GO TO 160
  165 ID(N,J)=0
      GO TO 160
  166 JEDN=JEDN+1
      ID(N,J)=JEDN
  160 CONTINUE
C
CE    PRINT GENERATED DATA FOR NODAL POINTS
CS    STAMPANJE GENERISANIH PODATAKA O CVOROVIMA
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,1)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      DO 70 I=1,NP
      WRITE(IZLAZ,5030) I,(ID(I,J),J=1,6),(CORD(I,J),J=1,3),LJ(I)
     1                   ,(ID(I,J),J=7,N6)
   70 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) N2,N1,KORA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) N2,N1,KORA
      STOP 'PROGRAM STOP: PAK021 - ULAZE1'
C
  110 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2200) N,NP
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6200) N,NP
      STOP 'PROGRAM STOP: PAK021 - ULAZE1'
C
! 1000 FORMAT(I5,1X,6I2,1X,A1,3F10.0,I5,I5,8X,6I2)
 1000 FORMAT(I5,A1,6I2,2X,3F10.0,I5,I5,8X,6I2)
 4000 FORMAT(I10,1X,6I2,1X,A1,3F15.0,I5,I5,8X,6I2)
 5010 FORMAT(1X,I5,4X,A1,6I3,2X,3F12.6,I5,I5,3X,6I5)
 8010 FORMAT(1X,I10,4X,A1,6I3,2X,3F15.6,I5,I5,3X,6I5)
 5030 FORMAT(1X,I10,6I10,3F15.9,I5,3X,6I5)
C-----------------------------------------------------------------------
 2000 FORMAT(//////1X,
     1'U C I T A N I   P O D A C I   O   C V O R O V I M A'/
     11X,51('-')//
     21X,'(MOGUCE GENERISANJE PODATAKA O CVOROVIMA,'/
     32X,'UCITATI POTREBAN BROJ KARTICA ZA DEFINISANJE SVIH CVOROVA)'//
     11X,' CVOR KORD    OGRANICENJA          K  O  O  R  D  I  N  A  T 
     1E     KOR KOSI'/
     11X,' BROJ SIST NX NY NZ WX WY WZ        X           Y           Z
     1      GEN SIST')
 2020 FORMAT(1X,
     1'G E N E R I S A N I   P O D A C I   O   C V O R O V I M A'/
     11X,57('-')///
     11X,' CVOR  J E D N A C I N A    B R O J     K  O  O  R  D  I  N  A
     1  T  E    KOSI'/
     11X,' BROJ                                    X           Y        
     1   Z      SIST')
 2100 FORMAT(///' BROJ N2=',I9,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
     1,I9,' I KONACNOG BROJA KORAKA KORA=',I5)
 2200 FORMAT(///' UCITANI CVOR BROJ N =',I10,' VECI JE OD UKUPNOG BROJA
     1CVORNIH TACAKA NP =',I10)
C-----------------------------------------------------------------------
 6000 FORMAT(//////1X,
     1'I N P U T   N O D A L   D A T A'/
     11X,31('-')///
     11X,' NODE COOR    CONSTRAINTS        C  O  O  R  D  I  N  A  T  E 
     1 S   STEP SKEW'/
     11X,' NUMB SYST NX NY NZ WX WY WZ        X           Y           Z 
     1     GENE SYST')
 6020 FORMAT(1X,
     1'G E N E R A T E D   N O D A L   D A T A'/
     11X,39('-')///
     11X,' NODE  E Q U A T I O N  N U M B E R   C  O  O  R  D  I  N  A  
     1T  E  S   SKEW'/
     11X,' NUMB                                    X           Y        
     1   Z      SYST')
 6100 FORMAT(///' NUMBER N2=',I10,' HAS NOT OBTAIN BY THE SUM OF NUBMBER
     1N1 =',I10,'AND FINITE NUMBER OF STEPS KORA=',I5)
 6200 FORMAT(///' NODE NUMBER N =',I10,' IS GREATER THAN TOTAL NUMBER OF
     1NODAL POINT NP =',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZE2(CORD,ID,NCVEL,NPA,NPI,NPP,LJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ NODAL POINT COORDINATES AND CONSTRAINTS
CS.    P R O G R A M
CS.        ZA UCITAVANJE KOORDINATA I OGRANICENJA CVORNIH TACAKA
CS.        SLOBODNO NUMERISANIH
C .
CE.           I=1,NPP (NPP - TOTAL NUMBEL OF NODAL POINT
CS.           I=1,NPP (NPP - UKUPAN BROJ CVORNIH TACAKA)
C .
CE.            CH - INDICATOR OF POLAR COORDINATES SYSTEM
CS.            CH - KARAKTER INDIKATOR POLARNOG KOORD. SISTEMA
C .
CE.            CORD(I,1) - X COORDINATE
CE.            CORD(I,2) - Y COORDINATE
CE.            CORD(I,3) - Z COORDINATE
CS.            CORD(I,1) - KOORDINATA  X
CS.            CORD(I,2) - KOORDINATA  Y
CS.            CORD(I,3) - KOORDINATA  Z
C .
CE.              ID(I,1) - CONSTRAINT OF TRANSLATION IN DIRECTION X AXIS
CE.              ID(I,2) - CONSTRAINT OF TRANSLATION IN DIRECTION Y AXIS
CE.              ID(I,3) - CONSTRAINT OF TRANSLATION IN DIRECTION Z AXIS
CE.              ID(I,4) - CONSTRAINT OF ROTATION AROUND X AXIS
CE.              ID(I,5) - CONSTRAINT OF ROTATION AROUND Y AXIS
CE.              ID(I,6) - CONSTRAINT OF ROTATION AROUND Z AXIS
CS.              ID(I,1) - OGRANICENJE POMERANJA U PRAVCU - X OSE
CS.              ID(I,2) - OGRANICENJE POMERANJA U PRAVCU - Y OSE
CS.              ID(I,3) - OGRANICENJE POMERANJA U PRAVCU - Z OSE
CS.              ID(I,4) - OGRANICENJE ROTACIJA OKO - X OSE
CS.              ID(I,5) - OGRANICENJE ROTACIJA OKO - Y OSE
CS.              ID(I,6) - OGRANICENJE ROTACIJA OKO - Z OSE
C .
C ......................................................................
C
      include 'paka.inc'
      CHARACTER*1 CH
      CHARACTER*250 ACOZ
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      COMMON /MIXEDM/ MIXED,IOPGS(6),NDS
      DIMENSION CORD(NPP,*),ID(NPP,*),NCVEL(*),LJ(*)
      DIMENSION DXYZ(3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZE2'
CZxxx
      CALL CZPROV(NP,1,34)
CZxxx
      N6=6
      IF(MIXED.EQ.1) N6=12
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         IF(JPS.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
         ELSE
            IF(JPBR.EQ.JPS1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)
            ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002) JPBR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002) JPBR
            ENDIF
         ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005)
      ENDIF
C
      I = 0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) N,CH,(ID(I,J),J=1,6),(CORD(I,J),J=1,3),KORC,LJ(I)
C     1                    ,(ID(I,J),J=7,N6)
      IF(INDFOR.EQ.2.and.ind56.eq.0)
     1READ(ACOZ,1000) N,(ID(I,J),J=1,6),CH,(CORD(I,J),J=1,3),KORC,LJ(I)
C     1                    ,(ID(I,J),J=7,N6)
      IF(INDFOR.EQ.2.and.ind56.eq.1)
     1READ(ACOZ,4000) N,(ID(I,J),J=1,6),CH,(CORD(I,J),J=1,3),KORC,LJ(I)
C     1                    ,(ID(I,J),J=7,N6)
      IF(I.EQ.1) THEN
        NPA=N
        NPI=N
      ELSE
        IF(NPA.LT.N) NPA=N
        IF(NPI.GT.N) NPI=N
      ENDIF
      NCVEL(I)=N
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,5010)N,CH,(ID(I,J),J=1,6),(CORD(I,J),J=1,3),KORC,LJ(I)
C     1                    ,(ID(I,J),J=7,N6)
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(I,1),CORD(I,2))
      IF(I.EQ.NP) GO TO 50
      GO TO 5
C
   20 NAUT=1
      N1=N
      I1=I
      KORA=KORC
      GO TO 5
C
CE    AUTOMATIC GENERATE DATA BETWEEN NODAL POINT N1 AND N2
CS    AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
C
   30 N2=N
      CALL DELJIV(N2-N1,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      RKORA = KORA
      RN1N2 = N2-N1
      DD = RKORA/RN1N2
      IF(DD.LT.0.0D00) GO TO 100
      DO 31 J=1,3
   31 DXYZ(J)=DD*(CORD(I,J)-CORD(I-1,J))
      IAUT=(N2-N1)/KORA-1
      I2=I+IAUT
      DO 32 K=1,3
   32 CORD(I2,K)=CORD(I,K)
      DO 33 K=1,6
   33 ID(I2,K)=ID(I,K)
      IF(JPS.EQ.1) LJ(I2)=LJ(I)
      NCVEL(I2)=NCVEL(I)
      I = I1
      N = N1
      DO 34 J=1,IAUT
      I=I+1
      N=N+KORA
      NCVEL(I)=N
      DO 35 K=1,3
   35 CORD(I,K)=CORD(I-1,K)+DXYZ(K)
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(I-1,1),CORD(I-1,2))
      DO 36 K=1,6
   36 ID(I,K)=ID(I1,K)
      IF(JPS.EQ.1) LJ(I)=LJ(I1)
   34 CONTINUE
      IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(I,1),CORD(I,2))
      I=I2
      IF(I.EQ.NP)THEN
        IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(I,1),CORD(I,2))
        GO TO 50
      ENDIF
      NAUT=0
      IF(KORC.EQ.0)THEN
        IF(CH.EQ.'p'.OR.CH.EQ.'P') CALL POLAR(CORD(I2,1),CORD(I2,2))
        GO TO 5
      ENDIF
      KORA=KORC
      NAUT=1
      N1=N2
      GO TO 5
C
   50 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) CALL WBROJK(KARTI,1)
      RETURN
C
  100 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) N2,N1,KORA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) N2,N1,KORA
      STOP
C
 1000 FORMAT(I5,1X,6I2,1X,A1,3F10.0,I5,I5,8X,6I2)
 4000 FORMAT(I10,1X,6I2,1X,A1,3F15.0,I5,I5,8X,6I2)
 5010 FORMAT(1X,I10,4X,A1,6I3,2X,3F15.9,I5,I5)
C-----------------------------------------------------------------------
 2000 FORMAT(//////1X,
     1'U C I T A N I   P O D A C I   O   C V O R O V I M A'/
     11X,51('-')//)
 2001 FORMAT(//////1X,
     1'U C I T A N I   P O D A C I   O   K O N T U R N I M   C V O R O V
     1 I M A'/1X,71('-')//)
 2002 FORMAT(//////1X,
     1'UCITANI PODACI O UNUTRASNJIM CVOROVIMA PODSTRUKTURE BROJ',I5/
     11X,61('-')//)
 2005 FORMAT(1X,'(MOGUCE GENERISANJE PODATAKA O CVOROVIMA,'/
     32X,'UCITATI POTREBAN BROJ KARTICA ZA DEFINISANJE SVIH CVOROVA)'//
     11X,' CVOR KORD    OGRANICENJA          K  O  O  R  D  I  N  A  T  
     1E     KOR KOSI'/
     11X,' BROJ SIST NX NY NZ WX WY WZ        X           Y           Z 
     1      GEN SIST')
 2100 FORMAT(///' BROJ N2=',I10,' NE MOZE SE DOBITI SABIRANJEM BROJA N1=
     1',I10,' I KONACNOG BROJA KORAKA KORA=',I5)
C-----------------------------------------------------------------------
 6000 FORMAT(//////1X,
     1'I N P U T   N O D A L   D A T A'/
     11X,31('-')//)
 6001 FORMAT(//////1X,
     1'I N P U T   D A T A   A B O U T   C O N T O U R   N O D E S'/
     11X,59('-')//)
 6002 FORMAT(//////1X,
     1'INPUT DATA ABOUT INTERNAL NODES FOR SUBSTRUCTURE NUMBER',I5/
     11X,60('-')//)
 6005 FORMAT(
     11X,' NODE COOR    CONSTRAINTS        C  O  O  R  D  I  N  A  T  E 
     1 S   STEP SKEW'/
     11X,' NUMB SYST NX NY NZ WX WY WZ        X           Y           Z 
     1     GENE SYST')
 6100 FORMAT(///' NUMBER N2=',I10,' HAS NOT OBTAIN BY THE SUM OF NUBMBER 
     1N1 =',I10,'AND FINITE NUMBER OF STEPS KORA=',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZE3(NCVEL,NCVELN,NELCV,NELCVN,ID,IDN,CORD,CORDN,NPP,
     1NPM,LJ,LJN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO UPDATE MATRIX ID
CS.    P R O G R A M
CS.       ZA KORIGOVANJE MATRICE ID
CS.       ZA UREDJIVANJE BROJEVA CVOROVA PO RASTUCEM REDOSLEDU,      
CS.       ZA STAMPANJE GENERISANIH PODATAKA O KOORDINATAMA I        
CS.       BROJEVA JEDNACINA CVORNIH TACAKA       
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      COMMON /REDJA / IREDJA
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /RESTAR/ TSTART,IREST
      DIMENSION CORD(NPP,*),ID(NPP,*),NCVEL(*),NELCV(*),LJ(*)
      DIMENSION CORDN(NP,*),IDN(NP,*),NCVELN(*),NELCVN(*),LJN(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZE3'
CZxxx
      CALL CZPROV(NP,1,34)
CZxxx
      JEDN=0
      NMPC=0
C
CE    UPDATE MATRIX ID
CS    KORIGOVANJE MATRICE ID, BROJEVI JEDNACINA SU U UCITANOM REDOSLEDU
C
      IREDJ=0
      IF(IREDJA.EQ.-1) IREDJ=-1
      IF(IREDJA.EQ.-1) IREDJA=0
      IF(IREDJA.EQ.0) THEN
         DO 50 N=1,NP
         DO 50 J=1,6
            IF(IOPGL(J).EQ.1) THEN
               ID(N,J)=0
               GO TO 50
            ENDIF
            IF(ID(N,J)) 54,56,55
   54       NMPC=NMPC+1
            GO TO 50
   55       ID(N,J)=0
            GO TO 50
   56       JEDN=JEDN+1
            ID(N,J)=JEDN
   50    CONTINUE
      ENDIF
C
CE    REDJANJE CVOROVA PO RASTUCEM REDU
CS    REDJANJE CVOROVA PO RASTUCEM REDU
C
      IF(IREST.EQ.2.AND.IOPGL(6).EQ.1.AND.IREDJ.EQ.0) GO TO 99
      J=0
      DO 10 I=1,NPM
         N=NELCV(I)
         IF(N.EQ.0) THEN
            NELCVN(I)=0
         ELSE
            J=J+1
            NELCVN(I)=J
            NCVELN(J)=NCVEL(N)
            DO 20 K=1,3
   20       CORDN(J,K)=CORD(N,K)
            DO 30 K=1,6
   30       IDN(J,K)=ID(N,K)
            IF(JPS.EQ.1) LJN(J)=LJ(N)
        ENDIF
   10 CONTINUE
C
CE    EQUALIZING VECTORS I MATRICES
CS    IZJEDNACAVANJE MATRICA I VEKTORA
C
      DO 80 I=1,NP
         NCVEL(I)=NCVELN(I)
         DO 85 K=1,3
   85    CORD(I,K)=CORDN(I,K)
         DO 90 K=1,6
   90    ID(I,K)=IDN(I,K)
         IF(JPS.EQ.1) LJ(I)=LJN(I)
   80 CONTINUE
      DO 95 I=1,NPM
   95 NELCV(I)=NELCVN(I)
C
CE    UPDATE MATRIX ID
CS    KORIGOVANJE MATRICE ID, BROJEVI JEDNACINA SU U RASTUCEM REDOSLEDU
C
   99 IF(IREDJA.EQ.1) THEN
         DO 60 N=1,NP
         DO 60 J=1,6
            IF(IOPGL(J).EQ.1) THEN
               ID(N,J)=0
               GO TO 60
            ENDIF
            IF(ID(N,J)) 64,66,65
   64       NMPC=NMPC+1
            GO TO 60
   65       ID(N,J)=0
            GO TO 60
   66       JEDN=JEDN+1
            ID(N,J)=JEDN
   60    CONTINUE
      ENDIF
C
CE    PRINT GENERATED DATA FOR NODAL POINT
CS    STAMPANJE GENERISANIH PODATAKA O CVOROVIMA
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      IF(JPS.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      ELSE
         IF(JPBR.EQ.JPS1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2021)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6021)
         ELSE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2022) JPBR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6022) JPBR
         ENDIF 
      ENDIF 
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2025)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6025)
      DO 70 I=1,NP
         IF(JPS.EQ.1) 
     1WRITE(IZLAZ,5030) NCVEL(I),(ID(I,J),J=1,6),(CORD(I,J),J=1,3),LJ(I)
         IF(JPS.GT.1) 
     1WRITE(IZLAZ,5030) NCVEL(I),(ID(I,J),J=1,6),(CORD(I,J),J=1,3)
   70 CONTINUE
      RETURN
C
 5030 FORMAT(1X,I10,6I10,3F12.6,I5)
C-----------------------------------------------------------------------
 2020 FORMAT(1X,'G E N E R I S A N I   P O D A C I   O   C V O R O V I M
     1 A'/1X,57('-')//)
 2021 FORMAT(1X,'G E N E R I S A N I   P O D A C I   O   K O N T U R N I
     1 M   C V O R O V I M A'/1X,77('-')//)
 2022 FORMAT(1X,'GENERISANI PODACI O UNUTRASNJIM CVOROVIMA PODSTRUKTURE 
     1BROJ',I5/1X,64('-')//)
 2025 FORMAT(
     11X,' CVOR  J E D N A C I N A    B R O J     K  O  O  R  D  I  N  A
     1  T  E    KOSI'/
     11X,' BROJ                                    X           Y        
     1   Z      SIST')
C-----------------------------------------------------------------------
 6020 FORMAT(1X,'G E N E R A T E D   N O D A L   D A T A'/1X,39('-')///)
 6021 FORMAT(1X,'G E N E R A T E D   D A T A   A B O U T   C O N T O U R
     1   N O D E S'/1X,67('-')//)
 6022 FORMAT(1X,'GENERATED DATA ABOUT INTERNAL NODES SUBSTRUCTURE NUMBER
     1',I5/1X,60('-')//)
 6025 FORMAT(
     11X,' NODE  E Q U A T I O N  N U M B E R   C  O  O  R  D  I  N  A  
     1T  E  S   SKEW'/
     11X,' NUMB                                    X           Y        
     1   Z      SYST')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZE4(NCVEL,NPA,NPI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .                                                                     
CS.    P R O G R A M                                                    
CS.        ZA UCITAVANJE KONTURNIH CVOROVA PODSTRUKTURE                 
C .                                                                     
CS.    I=1,NPK     (NPK - BROJ KONTURNIH CVOROVA PODSTRUKTURE)          
CS.            NCVEL(I) - VEKTOR BROJEVA KONTURNIH CVOROVA              
CS.                 NPA - NAJVECI KONTURNI CVOR PODSTRUKTURE        
CS.                 NPI - NAJMANJI KONTURNI CVOR PODSTRUKTURE           
C .                                                                     
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      DIMENSION NCVEL(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZE4'
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) JPBR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) JPBR
      ENDIF
C
      NK=1+(NPK-1)/10
      DO 10 K=1,NK
        IK=(K-1)*10+1
        NN=K*10
        IF(K.EQ.NK) NN=NPK
        CALL ISPITA(ACOZ)
        IF(K.EQ.1) KARTI=KARTIC
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*) (NCVEL(J),J=IK,NN)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) (NCVEL(J),J=IK,NN)
   10 CONTINUE 
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,5010) (NCVEL(J),J=1,NPK)
      DO 20 K=1,NPK
      N=NCVEL(K)
      IF(N.LT.NGI.OR.N.GT.NGA) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) N,NGI,NGA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) N,NGI,NGA
      STOP 'PROGRAM STOP: PAK021 - ULAZE4'
      ENDIF
        IF(K.EQ.1) THEN
          NPA=N
          NPI=N
        ELSE
          IF(NPA.LT.N) NPA=N
          IF(NPI.GT.N) NPI=N
        ENDIF
   20 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) CALL WBROJK(KARTI,0)
      RETURN
C
 1000   FORMAT(10I5)
 5010 FORMAT(6X,10I5)
C-----------------------------------------------------------------------
 2000 FORMAT(//////6X,
     1'UCITANI KONTURNI CVOROVI PODSTRUKTURE BROJ',I9/6X,47('-')//)
 2100 FORMAT(//' UCITANI KONTURNI CVOR N=',I9/' MORA SE NALAZITI U ',
     1'GRANICAMA IZMEDJU NAJMANJEG NGI =',I9/' I NAJVECEG KONTURNOG',
     1' CVORA NGA =',I9//' PROGRAM STOP'//)
C-----------------------------------------------------------------------
 6000 FORMAT(//////6X,
     1'READ CONTOUR NODES SUBSTRUCTURE NUMBER',I9/6X,43('-')//)
 6100 FORMAT(//' READ CONTOUR NODE,  N=',I9/' MUST BE IN BAUNDARY', 
     1' BETWEEN MINIMUM, NGI =',I9/' AND MAXIMUM CONTOUR',
     1' NODE, NGA =',I9//' PROGRAM STOP'//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZE5(NCVEL,NCVELN,NELCV,NELCVN,ID,IDN,CORD,CORDN,
     1NELCVG,NPP,NPM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .                                                                     
CE.    P R O G R A M
CE.       TO UPDATE MATRIX ID FOR CONTOUR NODES 
CS.    P R O G R A M                                                    
CS.        ZA UREDJIVANJE KONTURNIH CVOROVA PO RASTUCEM REDOSLEDU,      
CS.        ZA IZDVAJANJE PODATAKA O KONTURNIM CVOROVIMA PODSTRUKTURE,   
CS.        ZA STAMPANJE GENERISANIH PODATAKA O KOORDINATAMA I        
CS.        BROJEVA JEDNACINA KONTURNIH CVORNIH TACAKA            
C .                     
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      DIMENSION CORD(NPP,*),ID(NPP,*),NCVEL(*),NELCV(*),NELCVG(*)
      DIMENSION CORDN(NPG,*),IDN(NPG,*),NCVELN(*),NELCVN(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZE5'
C
CE    COORDINATE CONTOURAL NODES FOR SUBSTRUCTURE
CS    KOORDINATE KONTURNIH CVOROVA ZA PODSTRUKTURU
CS    REDJANJE CVOROVA U RASTUCEM REDU
C
      J=NP
      DO 10 I=1,NPM
      N=NELCV(I)
      IF(N.EQ.0) THEN
        NELCVN(I)=0
      ELSE
        J=J+1
        NELCVN(I)=J
        NN=NCVEL(N)      
        NCVELN(J-NP)=NN
        NM=NN-NGI+1
        MN=NELCVG(NM)
          DO 20 K=1,3
   20     CORD(J,K)=CORDN(MN,K)
          DO 30 K=1,6
   30     ID(J,K)=IDN(MN,K)
      ENDIF
   10 CONTINUE
C
CE    EQUALIZING VECTORS I MATRICES
CS    IZJEDNACAVANJE MATRICA I VEKTORA
C
      DO 80 I=1,NPK
   80 NCVEL(NP+I)=NCVELN(I)
      DO 95 I=1,NPM
   95 NELCV(I)=NELCVN(I)
C
CE    UPDATE MATRIX ID
CS    KORIGOVANJE MATRICE ID, BROJEVI JEDNACINA SU U RASTUCEM REDOSLEDU
C
      DO 60 N=NP+1,NPP
      DO 60 J=1,6
      IF(IOPGL(J).EQ.1) THEN
        ID(N,J)=0
        GO TO 60
      ENDIF
      IF(ID(N,J)) 65,65,66
   65 ID(N,J)=0
      GO TO 60
   66 JEDN=JEDN+1
      ID(N,J)=JEDN
   60 CONTINUE
C
CE    PRINT GENERATED DATA FOR CONTOUR NODAL POINT
CS    STAMPANJE GENERISANIH PODATAKA O KONTURNIM CVOROVIMA
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)JPBR
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)JPBR
      DO 70 I=NP+1,NPP
      WRITE(IZLAZ,5030) NCVEL(I),(ID(I,J),J=1,6),(CORD(I,J),J=1,3)
   70 CONTINUE
      RETURN
C
 5030 FORMAT(1X,I10,6I10,3F12.6,I5)
C-----------------------------------------------------------------------
 2020 FORMAT(1X,'GENERISANI PODACI O KONTURNIM CVOROVIMA PODSTRUKTURE BR     
     1OJ',I5/1X,62('-')///
     11X,' CVOR  J E D N A C I N A    B R O J     K  O  O  R  D  I  N  A
     1  T  E    KOSI'/
     11X,' BROJ                                    X           Y        
     1   Z      SIST')
C-----------------------------------------------------------------------
 6020 FORMAT(1X,'GENERATED DATA ABOUT CONTOUR NODES SUBSTRUCTURE NUMBER'
     1,I5/1X,59('-')///
     11X,' NODE  E Q U A T I O N  N U M B E R   C  O  O  R  D  I  N  A  
     1T  E  S   SKEW'/
     11X,' NUMB                                    X           Y        
     1   Z      SYST')
C-----------------------------------------------------------------------
      END
C======================================================================
C
C======================================================================
      SUBROUTINE VEZACV(NCVEL,NELCV,NP,NPI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS NODES IN FREE NUMERATION
CS.   P R O G R A M
CS.        ZA VEZU BROJEVA CVOROVA I NJIHOVOG REDNOG BROJA U VEKTORU
CS.        BROJEVA CVOROVA U SLOBODNOJ NUMERACIJI
C .
CE.    V E C T O R S
CE.               NCVEL(NP) - VECTOR NUMBER OF NODES
CE.        NELCV(NPA-NPI+1) - VECTOR CONECT NUMBER OF NODES
CS.    V E K T O R I
CS.               NCVEL(NP) - VEKTOR BROJEVA CVOROVA
CS.        NELCV(NPA-NPI+1) - VEKTOR VEZE BROJEVA CVOROVA I NJIHOVOG
CS.                           REDNOG BROJA U VEKTORU BROJEVA CVOROVA
C .
C ......................................................................
C
      DIMENSION NCVEL(*),NELCV(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' VEZACV'
      DO 10 I=1,NP
         J=NCVEL(I)
         J=J-NPI+1
         NELCV(J)=I
   10 CONTINUE
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE VEZACG(NCVEL,NELCV,NPK,NII,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS CONTOUR NODES IN FREE NUMERATION 
CS.   P R O G R A M
CS.        ZA VEZU BROJEVA KONTURNIH CVOROVA I NJIHOVOG REDNOG BROJA U 
CS.        VEKTORU BROJEVA CVOROVA U SLOBODNOJ NUMERACIJI
C .
CE.    V E C T O R S
CE.               NCVEL(NPK) - VECTOR NUMBER OF CONTOUR NODES
CE.        NELCV(NAA-NII+1) - VECTOR CONECT NUMBER OF CONTOUR NODES
CS.    V E K T O R I
CS.               NCVEL(NP) - VEKTOR BROJEVA KONTURNIH CVOROVA
CS.        NELCV(NAA-NII+1) - VEKTOR VEZE BROJEVA CVOROVA I NJIHOVOG
CS.                           REDNOG BROJA U VEKTORU BROJEVA CVOROVA
C .
C ......................................................................
C
      DIMENSION NCVEL(*),NELCV(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' VEZACG'
      DO 10 I=1,NPK
         J=NCVEL(I)
         J=J-NII+1
         NELCV(J)=I+NP
   10 CONTINUE
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE POLAR(CORD1,CORD2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSFORMATION OF POLAR COORDINATES
CS.   P R O G R A M
CS.      ZA TRANSFORMACIJU POLARNIH KOORDINATA NA DEKARTOVE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' POLAR '
      CORD2=CORD2*1.745329252D-2
      X=CORD1*DCOS(CORD2)
      CORD2=CORD1*DSIN(CORD2)
      CORD1=X
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE TGRAFK(CORD,NCVEL,ICVEL,NPP,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT COORDINATES IN UNIVERSAL FILE
CS.    P R O G R A M
CS.       ZA STAMPANJE KOORDINATA CVOROVA U UNIVERZALNI FILE
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /NIDEAS/ IDEAS
      DIMENSION CORD(NPP,*),NCVEL(*)
C
      IND=-1
      IF(IDEAS.GT.6) THEN
         ITYP=2411
         IT1=0
         IT2=0
         ICOL=11
      ELSEIF(IDEAS.EQ.6) THEN
         ITYP=781
         IT1=0
         IT2=0
         ICOL=11
      ELSE
         ITYP=15
         IT1=0
         IT2=0
         ICOL=8
      ENDIF
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      DO 10 I=1,NP
      IF(ICVEL.EQ.0) THEN
        NI=I
      ELSE
        NI=NCVEL(I)
      ENDIF
      IF(IDEAS.GT.0) THEN
         WRITE(II,5411) NI,IT1,IT2,ICOL,(CORD(I,J),J=1,3)
      ELSE
         WRITE(II,5000) NI,IT1,IT2,ICOL,(CORD(I,J),J=1,3)
      ENDIF
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I6)
 5000 FORMAT(4I10,3(1PE13.5))
 5411 FORMAT(4I10/3(1PD25.16))
      END
C======================================================================
C
C======================================================================
      SUBROUTINE TGRAFB(ID,NCVEL,ICVEL,NPP,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT D.O.F. IN UNIVERSAL FILE
CS.    P R O G R A M
CS.       ZA STAMPANJE OGRANICENJA CVOROVA U UNIVERZALNI FILE
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /NIDEAS/ IDEAS
      DIMENSION ID(NPP,*),NCVEL(*)
      DIMENSION IDD(6)
C
      IND=-1
      ITYP=757
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      IND1=JPBR-1
      WRITE(II,5000) IND1
      WRITE(II,5200) IND1
C     COLOR
      ICOL=4
      DO 10 I=1,NP
         IF(ICVEL.EQ.0) THEN
            NI=I
         ELSE
            NI=NCVEL(I)
         ENDIF
         IST=0
         DO 20 IJ=1,6
            IF(ID(I,IJ).GT.0) THEN
               IDD(IJ)=0
            ELSE
               IDD(IJ)=1
               IST=1
            ENDIF
   20    CONTINUE
         IF(IST.EQ.1.OR.IDEAS.LT.6)
     1   WRITE(II,5000) NI,ICOL,(IDD(J),J=1,6)
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I6)
 5000 FORMAT(2I10,6I2)
 5200 FORMAT(' DOF SET',I10)
      END
C======================================================================
C
C======================================================================
      SUBROUTINE TGRAFR(ID,NCVEL,ICVEL,NPP,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT D.O.F. IN UNIVERSAL FILE
CS.    P R O G R A M
CS.       ZA STAMPANJE OGRANICENJA CVOROVA U UNIVERZALNI FILE
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /NIDEAS/ IDEAS
      DIMENSION ID(NPP,*),NCVEL(*)
      DIMENSION IDD(6)
C
      IF(IDEAS.LT.6) RETURN
      IND=-1
      ITYP=791
      NTYP=7
      JEDAN=1
      NULA=0
      ZERO=0.
C     COLOR
      ICOL=4
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      IND1=JPBR-1
      WRITE(II,5000) IND1,JEDAN
      WRITE(II,5200) IND1 
      DO 10 I=1,NP
         IF(ICVEL.EQ.0) THEN
            NI=I
         ELSE
            NI=NCVEL(I)
         ENDIF
         IST=0
         DO 20 IJ=1,6
            IF(ID(I,IJ).GT.0) THEN
               IDD(IJ)=0
            ELSE
               IDD(IJ)=1
               IST=1
            ENDIF
   20    CONTINUE
         IF(IST.EQ.1) THEN
            WRITE(II,5000) NI,ICOL,(IDD(J),J=1,6),NULA,NTYP,NULA
            WRITE(II,5300) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     1                     NULA,NULA,NULA,NULA,NULA,NULA
         ENDIF
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I6)
 5000 FORMAT(2I10,9I2)
 5200 FORMAT(' RESTRAIN SET',I10)
 5300 FORMAT(3(1PD25.16)/3(1PD25.16)/6I10)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UCPUVA(NPODS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO FORM POINTERS AND CALL SUBROUTINE FOR READ DATA
CE.        ABOUT INITIAL CONDITIONS
CS.    P R O G R A M
CS.        ZA FORMIRANJE REPERA I POZIVANJE PROGRAMA
CS.        ZA UCITAVANJE PODATAKA O POCETNIM USLOVIMA
C .
CE.    POINTERS
CE.        LPUU  -  INITIAL DISPLACEMENTS
CE.        LPUV  -  INITIAL VELOCITIES
CE.        LPUA  -  INITIAL ACCELERATIONS
CS.    R E P E R I
CS.        LPUU  -  POCETNA POMERANJA
CS.        LPUV  -  POCETNE BRZINE
CS.        LPUA  -  POCETNA UBRZANJA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DUPLAP/ IDVA
      COMMON /REZREP/ LPUUU,LPUUV,LPUUA
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' UCPUVA'
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         CALL WBROJK(KARTIC,1)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      ENDIF
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
C*      LPUU=LRAD
      LPUU=LID+NP*6
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
      JEDN3=JEDN*3
      CALL CLEAR(A(LPUU),JEDN3)
      IF(IPUU.EQ.0) GO TO 20
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      ENDIF
      CALL UPUVA(A(LPUU),A(LID),IPUU,JEDN,NP)
   20 IF(IPUV.EQ.0) GO TO 30
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2030)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6030)
      ENDIF
      CALL UPUVA(A(LPUV),A(LID),IPUV,JEDN,NP)
   30 IF(IPUA.EQ.0) GO TO 40
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2040)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6040)
      ENDIF
      CALL UPUVA(A(LPUA),A(LID),IPUA,JEDN,NP)
   40 CONTINUE
      LMAX13=MAX13
      NPODS(JPBR,59)=LMAX13+1
      CALL WRITDD(A(LPUU),JEDN,IPODS,LMAX13,LDUZI)
      NPODS(JPBR,55)=LMAX13+1
      CALL WRITDD(A(LPUV),JEDN,IPODS,LMAX13,LDUZI)
      NPODS(JPBR,57)=LMAX13+1
      CALL WRITDD(A(LPUA),JEDN,IPODS,LMAX13,LDUZI)
      MAX13=LMAX13
      LMAX=LPUU
      RETURN
C
 1000 FORMAT(8I5)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,'P O D A C I    O    P O C E T N I M     U S L O V I M A
     1'/6X,55('-')///
     111X,'BROJ CVOROVA SA POMERANJIMA RAZLICITIM OD NULE .. IPUU =',I5/
     116X,'EQ.0; SVA POCETNA POMERANJA JEDNAKA NULI'///
     111X,'BROJ CVOROVA SA BRZINAMA RAZLICITIM OD NULE ..... IPUV =',I5/
     116X,'EQ.0; SVE POCETNE BRZINE JEDNAKE NULI'///
     111X,'BROJ CVOROVA SA UBRZANJIMA RAZLICITIM OD NULE ... IPUA =',I5/
     116X,'EQ.0; SVA POCETNA UBRZANJA JEDNAKA NULI'///
     111X,'INDIKATOR STAMPANJA POMERANJA ................... ISUU =',I5/
     116X,'LT.0; ZELJENO STAMPANJE(PODPROGRAM KRSTAMP)'/
     116X,'EQ.0; STAMPAJU SE'/
     116X,'GT.0; NE STAMPAJU SE'///
     111X,'INDIKATOR STAMPANJA BRZINA ...................... ISUV =',I5/
     116X,'LT.0; ZELJENO STAMPANJE(PODPROGRAM KRSTAMP)'/
     116X,'EQ.0; STAMPAJU SE'/
     116X,'GT.0; NE STAMPAJU SE'///
     111X,'INDIKATOR STAMPANJA UBRZANJA .................... ISUA =',I5/
     116X,'LT.0; ZELJENO STAMPANJE(PODPROGRAM KRSTAMP)'/
     116X,'EQ.0; STAMPAJU SE'/
     116X,'GT.0; NE STAMPAJU SE')
 2020 FORMAT(//////6X,'P O D A C I    O    P O C E T N I M',
     1'    P O M E R A NJ I M A'/6X,60('-')///)
 2030 FORMAT(//////6X,'P O D A C I    O    P O C E T N I M',
     1'    B R Z I N A M A'/6X,55('-')///)
 2040 FORMAT(//////6X,'P O D A C I    O    P O C E T N I M',
     1'    U B R Z A NJ I M A'/6X,58('-')///)
C-----------------------------------------------------------------------
 6000 FORMAT(6X,'D A T A    A B O U T    I N I T I A L    C O N D I T I
     *O N'/6X,58('-')///
     111X,'NUMBER OF NODES WITH NON-ZERO DISPLACEMENTS ..... IPUU =',I5/
     116X,'EQ.0; ALL INITIAL DISPLACEMENTS EQUAL TO ZERO'///
     111X,'NUMBER OF NODES WITH NON-ZERO VELOCITIES ........ IPUV =',I5/
     116X,'EQ.0; ALL INITIAL VELOCITIES ARE EQUAL TO ZERO'///
     111X,'NUMBER OF NODES WITH NON-ZERO ACCELERATIONS ..... IPUA =',I5/
     116X,'EQ.0; ALL INITIAL ACCELERATIONS ARE EQUAL TO ZERO'///
     111X,'PRINT DISPLACEMENT CODE ......................... ISUU =',I5/
     116X,'LT.0; DISPLACEMENTS ARE PRINTED (PROGRAM KRSTAMP)'/
     116X,'EQ.0; DISPLACEMENTS ARE PRINTED'/
     116X,'GT.0; DISPLACEMENT ARE NOT PRINTED'///
     111X,'PRINT VELOCITIES CODE ........................... ISUV =',I5/
     116X,'LT.0; VELOCITIES ARE PRINTED (PROGRAM KRSTAMP)'/
     116X,'EQ.0; VELOCITIES ARE PRINTED'/
     116X,'GT.0; VELOCITIES ARE NOT PRINTED'///
     111X,'PRINT ACCELERATIONS CODE ........................ ISUA =',I5/
     116X,'LT.0; ACCELERATIONS ARE PRINTED (PROGRAM KRSTAMP)'/
     116X,'EQ.0; ACCELERATIONS ARE PRINTED'/
     116X,'GT.0; ACCELERATIONS ARE NOT PRINTED')
 6020 FORMAT(//////6X,'D A T A   A B O U T   I N I T I A L',
     1'   D I S P L A C E M E N T S'/6X,64('-')///)
 6030 FORMAT(//////6X,'D A T A   A B O U T   I N I T I A L',
     1'   V E L O C I T I E S'/6X,57('-')///)
 6040 FORMAT(//////6X,'D A T A   A B O U T   I N I T I A L',
     1'   A C C E L E R A T I O N S'/6X,63('-')///)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UPUVA(PU,ID,IPU,JEDN,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO READ INITIAL DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
CS.   PROGRAM
CS.      ZA UCITAVANJE POCENTIH POMERANJA, BRZINE I UBRZANJA
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION PU(*),ID(NP,*)
      DIMENSION POM(6),DXYZ(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' UPUVA'
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      ENDIF
      I=0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) N,KN,(POM(JI),JI=1,6)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) N,KN,(POM(JI),JI=1,6)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,5000) N,KN,(POM(JI),JI=1,6)
      DO 6 J=1,6
      IF(ID(N,J).EQ.0) GO TO 6
      NJ=ID(N,J)
      PU(NJ)=POM(J)
    6 CONTINUE
      IF(NAUT.GT.0) GO TO 30
      IF(KN.NE.0) GO TO 20
      IF(I.EQ.IPU) GO TO 50
      GO TO 5
   20 NAUT=1
      N1=N
      KORA = KN
      GO TO 5
   30 N2=N
      CALL DELJIV(N2-N1,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      RKORA = KORA
      RN1N2 = N2 - N1
      DD = RKORA/RN1N2
      DO 31 J=1,6
      NJ12=0
      NJ1=ID(N1,J)
      IF(NJ1.GT.0) NJ12=NJ12+1
      NJ2=ID(N2,J)
      IF(NJ2.GT.0) NJ12=NJ12+1
      IF(NJ12-1) 31,110,32
   32 DXYZ(J)=DD*(PU(NJ2)-PU(NJ1))
   31 CONTINUE
      IAUT=(N2-N1)/KORA-1
      N=N1
      DO 35 L=1,IAUT
      N=N+KORA
      DO 35 J=1,6
      NJ12=0
      NJ1=ID(N,J)
      IF(NJ1.GT.0) NJ12=NJ12+1
      NJ2=ID(N-KORA,J)
      IF(NJ2.GT.0) NJ12=NJ12+1
      IF(NJ12-1) 35,110,34
   34 PU(NJ1)=PU(NJ2)+DXYZ(J)
   35 CONTINUE
      I=I+IAUT
      IF(I.EQ.IPU) GO TO 50
      NAUT=0
      IF(KN.EQ.0) GO TO 5
      KORA = KN
      NAUT=1
      N1=N2
      GO TO 5
   50 IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,1)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      DO 70 IN=1,JEDN
      PP=PU(IN)
      IPP=1
      IF(DABS(PP).LT.0.1D-10) IPP=0
      IF(IPP.EQ.0) GO TO 70
      WRITE(IZLAZ,5100) IN,PP
   70 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2105) N2,N1,KORA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6105) N2,N1,KORA
      STOP
C
  110 NN=N-KORA
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2110) NN,N
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6110) NN,N
      STOP
C
 1000 FORMAT(2I5,6F10.2)
 5000 FORMAT(6X,2I10,4X,6D12.4)
 5100 FORMAT(10X,I10,15X,D12.4)
C-----------------------------------------------------------------------
 2010 FORMAT(6X,'U C I T A N I    P O D A C I'/6X,28('-')///11X,'CVOR BR
     1OJ',3X,'KORAK',9X,'X',11X,'Y',11X,'Z',10X,'WX',10X,'WY',10X,'WZ')
 2020 FORMAT(6X,'G E N E R I S A N I    P O D A C I'/6X,34('-')///11X,
     1'   JEDNACINA BROJ',9X,'V R E D N O S T')
 2105 FORMAT(///' BROJ N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
     1,I5,' I KONACNOG BROJA KORAKA KORA=',I5)
 2110 FORMAT(///' CVOROVI   N1 =',I5,'  I   N2 =',I5,' NEMAJU ISTA OGRAN
     1ICENJA PA GENERISANJE PODATAKA NIJE MOGUCE')
C-----------------------------------------------------------------------
 6010 FORMAT(6X,'I N P U T     D A T A'/6X,21('-')///11X,'NODE NO ',
     *3X,'STEP',10X,'X',11X,'Y',11X,'Z',10X,'WX',10X,'WY',10X,'WZ')
 6020 FORMAT(6X,'G E N E R A T E D    D A T A'/6X,28('-')///11X,
     1'   EQUATION NO',13X,'V A L U E')
 6105 FORMAT(///' BROJ N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
     1,I5,' I KONACNOG BROJA KORAKA KORA=',I5)
 6110 FORMAT(///' CVOROVI   N1 =',I5,'  I   N2 =',I5,' NEMAJU ISTA OGRAN
     1ICENJA PA GENERISANJE PODATAKA NIJE MOGUCE')
C-----------------------------------------------------------------------
      END
C=========================================================
      SUBROUTINE MPCONS(CMPC,MPC,ID,NELCV,  ICVEL,NPI,NP,
     &                  LMAX,IULAZ,IZLAZ,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   M U L T Y P O I N T    C O N S T R A I N T
CS    MMP    = UKUPAN BROJ RAZLICITIH LINEARNIH VEZA
CS    NMPC   = UKUPAN BROJ ZAVISNIH STEPENI SLOBODE (-1 U MATRICI ID)
CS    NEZAV  = MAKSIMALAN BROJ NEZAVISNIH U LINEARNIM VEZAMA
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /SRBAMP/ IPO(15,342)
      DIMENSION CMPC(MMP,*),MPC(NEZA1,*),ID(NP,*),NELCV(*),IPOM(15)
C
      GO TO (1,2,3),IND
C
    1 CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) MMP,NEZAV
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) MMP,NEZAV
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      CALL WBROJK(KARTIC,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2080) MMP,NEZAV
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6080) MMP,NEZAV
      ENDIF
      IF(NEZAV.EQ.0) NEZAV=1
      NEZA1=NEZAV+1
      IF(NEZAV.GT.6) STOP'***   MPCONS   ***'
      LCMPC=LMAX
      LMPC =LCMPC+MMP*NEZAV*IDVA
      RETURN
C
C
    2 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      ENDIF
      DO 10 K=1,MMP
      CALL ISPITA(ACOZ)
      IF(K.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) IRB,(CMPC(IRB,J),J=1,NEZAV)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) IRB,(CMPC(IRB,J),J=1,NEZAV)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,1005) IRB,(CMPC(IRB,J),J=1,NEZAV)
      IF(IRB.GT.MMP) STOP'***   MPCONS   ***'
   10 CONTINUE
C
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      ENDIF
      KI=0
      NMPC0=NMPC
      CALL ICLEAR(IPOM,15)
!      IF((3+2*NEZAV).GT.15) STOP 'SRBA PAK021-3+2*NEZAV.GT.15'
!      IF(NMPC0.GT.342) STOP 'SRBA PAK021-NMPC0.GT.342'
      DO 50 I=1,NMPC0
       CALL ISPITA(ACOZ)
CS   IPOM(1) = REDNI BROJ LINEARNE VEZE
CS   IPOM(2)   ZAVISNA JEDNACINA ISKAZANA PREKO CVORA I PRAVCA  DOF
CS   IPOM(3)
CS   IPOM(4,5) ... = PODACI O NEZAVISNIM STEPENIMA SLOBODE
       IF(INDFOR.EQ.1)
     1 READ(IULAZ,*) (IPOM(J),J=1,3+2*NEZAV)
       IF(INDFOR.EQ.2.and.ind56.eq.0)
     1 READ(ACOZ,1010) (IPOM(J),J=1,3+2*NEZAV)
       IF(INDFOR.EQ.2.and.ind56.eq.1)
     1 READ(ACOZ,4010) (IPOM(J),J=1,3+2*NEZAV)
c ovo je samo za srbu za geometrijsku nelinearnost kod sfere?????     
c       CALL IJEDN1(IPO(1,I),IPOM,3+2*NEZAV)
       IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1 WRITE(IZLAZ,1015) (IPOM(J),J=1,3+2*NEZAV)
CS.  DIREKTNO SPREZANJE DVA STEPENA SLOBODE   ( DOF1 = DOF2 ) 
CS   TADA SE KORIGUJE SAMO MATRICA ID
       IF((NEZAV.EQ.1.OR.IPOM(6).EQ.0).AND.
     1     DABS(CMPC(IPOM(1),1)-1.D0).LT.1.D-10)THEN
         NN1=IPOM(2)
         NN2=IPOM(4)
         N1=NN1
         N2=NN2
         IF(ICVEL.NE.0)N1=NELCV(NN1-NPI+1)
         IF(ICVEL.NE.0)N2=NELCV(NN2-NPI+1)
         IF(ID(N1,IPOM(3)).GE.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN1
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN1
      STOP'***   MPCONS   ***'
         ENDIF 
         IF(ID(N2,IPOM(5)).LE.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN2
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN2
      STOP'***   MPCONS   ***'
         ENDIF 
         ID(N1,IPOM(3))=ID(N2,IPOM(5))
         NMPC=NMPC-1
       ELSE
CS   STANDARDNA LINEARNA VEZA, VEZE ODGOVARAJUCIH JEDNACINA
       KI=KI+1
       MPC(1,KI)=IPOM(1)
       NN1=IPOM(2)
       N1=NN1
       IF(ICVEL.NE.0)N1=NELCV(NN1-NPI+1)
       IF(ID(N1,IPOM(3)).GE.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN1
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN1
      STOP'***   MPCONS   ***'
       ENDIF 
       ID(N1,IPOM(3))=-KI
       KK=2
       DO 40 K=1,NEZAV
         KK=KK+2
         NN2=IPOM(KK)
         N2=NN2
         IF(ICVEL.NE.0)N2=NELCV(NN2-NPI+1)
       IF(ID(N2,IPOM(KK+1)).LE.0)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN2
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN2
      STOP'***   MPCONS   ***'
       ENDIF 
   40    MPC(K+1,KI)=ID(N2,IPOM(KK+1))
      ENDIF
   50 CONTINUE         
C
CS   KOREKCIJA REPERA
C
      IF(NMPC.LT.NMPC0)THEN
       IF(NMPC.GT.0)THEN
        LMAX=LMPC+NMPC*(NEZAV+1)
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1
       ELSE
        LMAX=LCMPC
       ENDIF
      ENDIF
      RETURN
C
C
    3 CONTINUE
      RETURN
C
 1000 FORMAT(I5,6F10.0)
 1005 FORMAT(I5,1X,6(1PD12.5))
 1010 FORMAT(15I5)
 4010 FORMAT(I5,7(I10,I5))
 1015 FORMAT(I5,21X,I5,7X,I5/(14I5))
C-----------------------------------------------------------------------
 2000 FORMAT(///6X,
     1'P O D A C I   O   K O F I C I J E N T I M A   V E Z A'/
     16X,53('-')///3X,'VEZA',28X,'KOEFICIJENTI'/) 
 2010 FORMAT(//3X,'VEZA',10X,'ZAVISAN:   CVOR - STEPEN SLOBODE'/
     13X,'NEZAVISNI:    CVOR - STEPEN SLOBODE'/) 
 2080 FORMAT(///6X,
     1'P O D A C I   O   S P R E G N U T I M   P O M E R A N J I M A'/
     16X,61('-')///
     111X,'UKUPAN BROJ RAZLICITIH LINEARNIH VEZA ..........   MMP =',I5/
     111X,'BROJ NEZAVISNIH (MAKSIMUM 6) ................... NEZAV =',I5)
 2100 FORMAT(//6X,'NEPRAVILNI PODACI O OGRANICENJIMA CVORA  N =',I10,
     1'  U MATRICI  ID  '/)
C-----------------------------------------------------------------------
 6000 FORMAT(///6X,
     1'D A T A    A B O U T   K O E F I C I J E N T S'/
     16X,46('-')///3X,'SET ',28X,'KOEFICIENTS'/) 
 6010 FORMAT(//3X,'SET ',10X,'CONSTRAIND: NOD - D.O.F.'/
     13X,'INDEPENDENT:   NOD - D.O.F.'/) 
 6080 FORMAT(///6X,
     1'D A T A    A B O U T   C O N S T R A I N D   D. O. F'/
     16X,52('-')///
     111X,'NUMBER OF DIFFERENT MULTIPOINT CONSTRAINTS .....   MMP =',I5/
     111X,'NUMBER OF INDEPENDENT D.O.F. (MAXIMUM 6) ,...... NEZAV =',I5)
 6100 FORMAT(//6X,'IMPROPER DATA ABOUT CONSTRAINTS FOR NOD  N =',I10,
     1'  IN MATRIX  ID  '/)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE DIMENZ(CORD,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA TOLERANCIJU ZA CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.   PROGRAM
CE.      TO CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C . 
C ......................................................................
C 
      COMMON /MAXDUZ/ XL,YL,ZL
      COMMON /CDEBUG/ IDEBUG
      DIMENSION CORD(NP,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' DIMENZ'
C
      XMAX=-1.D10
      YMAX=-1.D10
      ZMAX=-1.D10
      XMIN= 1.D10
      YMIN= 1.D10
      ZMIN= 1.D10
      DO 10 I=1,NP
         IF(XMAX.LT.CORD(I,1)) XMAX=CORD(I,1)
         IF(YMAX.LT.CORD(I,2)) YMAX=CORD(I,2)
         IF(ZMAX.LT.CORD(I,3)) ZMAX=CORD(I,3)
         IF(XMIN.GT.CORD(I,1)) XMIN=CORD(I,1)
         IF(YMIN.GT.CORD(I,2)) YMIN=CORD(I,2)
         IF(ZMIN.GT.CORD(I,3)) ZMIN=CORD(I,3)
   10 CONTINUE
      XLL=DABS(XMAX-XMIN)
      YLL=DABS(YMAX-YMIN)
      ZLL=DABS(ZMAX-ZMIN)
      IF(XL.LT.XLL) XL=XLL
      IF(YL.LT.YLL) YL=YLL
      IF(ZL.LT.ZLL) ZL=ZLL
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INDST(ISKC,KSK,NELCV,NCVEL,ICVEL,NP,NPI,NCVPR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
C
CS  PODPROGRAM UCITAVA I GENERISE INDIKATORE ZA STAMPU ZELJENIH CVOROVA
CE  READ AND GENERATE INDICATORS FOR PRINTOUT NODES
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION ISKC(*),KSK(*),NELCV(*),NCVEL(*)
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(IZLAZ,2007)
      IF(ISRPS.EQ.1)
     *WRITE(IZLAZ,6007)
      ENDIF
C
      I = 0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) N,KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) N,KORC
      IF(ICVEL.EQ.0) THEN
         ISKC(N)=1
      ELSE
         J=N-NPI+1
         NI=NELCV(J)
         ISKC(NI)=1
      ENDIF
      IF(I.GT.NCVPR) STOP ' STOP - INDST - (I.GT.NCVPR) - PAK021' 
      IF(N.EQ.0) GO TO 50
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,1008) N,KORC
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(I.EQ.NCVPR) GO TO 50
      GO TO 5
C
   20 NAUT=1
      N1=N
      I1=I
      KORA=KORC
      GO TO 5
C
CE    AUTOMATIC GENERATE DATA BETWEEN NODAL POINT N1 AND N2
CS    AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
C
   30 N2=N
      CALL DELJIV(N2-N1,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      IAUT=(N2-N1)/KORA-1
      I2=I+IAUT
      I = I1
      N = N1
      DO 34 J=1,IAUT
         I=I+1
         N=N+KORA
         IF(ICVEL.EQ.0) THEN
            ISKC(N)=1
         ELSE
            JI=N-NPI+1
            NI=NELCV(JI)
            ISKC(NI)=1
         ENDIF
   34 CONTINUE
      I=I2
      IF(I.GT.NCVPR) STOP ' STOP - INDST - (I.GT.NCVPR) - PAK021' 
      IF(I.EQ.NCVPR) GO TO 50
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      N1=N2
      GO TO 5
C
   50 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      ENDIF
      J=0
      DO 23 I=1,NP
         IF(ISKC(I).EQ.0) GO TO 23
         J=J+1
         IF(ICVEL.EQ.0) THEN
            NI=I
         ELSE
            NI=NCVEL(I)
         ENDIF
            KSK(J)=NI
   23 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         CALL WBROJK(KARTI,1)
      IF(ISRPS.EQ.0)
     *WRITE(IZLAZ,2020) J
      IF(ISRPS.EQ.1)
     *WRITE(IZLAZ,6020) J
         WRITE(IZLAZ,1000) (KSK(I),I=1,J)
      ENDIF
      RETURN
C
  100 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) N2,N1,KORA
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) N2,N1,KORA
      STOP
C
 1000 FORMAT(14I5)
 1008 FORMAT(11X,I5,12X,I5)
C-----------------------------------------------------------------------
 2007 FORMAT(///6X,
     1'UCITANI ZELJENI CVOROVI ZA KOJE SE STAMPAJU REZULTATI'/6X,53('-')
     1//11X,'CVOR BROJ',10X,'KORAK')
 2020 FORMAT(///6X,
     1'GENERISANI ZELJENI CVOROVI ZA KOJE SE STAMPAJU REZULTATI',I9/6X,
     161('-')//)
 2100 FORMAT(///' BROJ N2=',I9,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
     1,I9,' I KONACNOG BROJA KORAKA KORA=',I5)
C-----------------------------------------------------------------------
 6007 FORMAT(///6X,
     1'NODES FOR RESULTS PRINTOUT'/6X,26('-')
     1//11X,'NODE NUMB.',9X,'STEP')
 6020 FORMAT(///6X,
     1'GENERATED NODES FOR RESULTS PRINTOUT',I5/6X,41('-')//)
 6100 FORMAT(///' NUMBER N2=',I9,' HAS NOT OBTAIN BY THE SUM OF NUBMBER
     1N1 =',I9,'AND FINITE NUMBER OF STEPS KORA=',I5)
C-----------------------------------------------------------------------
      END
