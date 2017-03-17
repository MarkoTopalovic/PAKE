C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZkt(CORD,ID,LJ,nktlnc,lncvz,nbvez,ncgl,nvez,ncvez,
     +                  ivez,xmass)
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
CE.           I=1,NP (NP - TOTAL NUMBEL OF NODAL POINT
CS.           I=1,NP (NP - UKUPAN BROJ CVORNIH TACAKA)
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
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /krut3d/ kt3d,mmpogr
      DIMENSION DXYZ(3)
      dimension nktlnc(*),lncvz(*),nbvez(*),ncgl(*),nvez(*),
     1          ncvez(nkrt,*),id(np,*),cord(np,*),lj(*),ivez(*),
     2          xmass(nkrt,*)
      common /ielmkt/ inmkt
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ULAZEkt'
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
      I = 0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
      IF(N.GT.NP) GO TO 110
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,5010)N,CH,(ID(N,J),J=1,6),(CORD(N,J),J=1,3),KORC,LJ(N)
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
      DO 36 K=1,6
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
      if(nkrt.gt.0) then
          do 412 ij=1,nkrt
              ncgl(ij)=0
              nvez(ij)=0
              do 412 ji=1,mxvez
                   ncvez(ij,ji)=0
  412 continue
cc
cc
         n3d=2
         if(kt3d.eq.1) n3d=3
cc
cc
         do 531 iv11=1,nkrt-1
  531    ivez(iv11)=0
        mmp0=0
        mmpogr=0
        nezav0=1
        ii=0
        icv=1
        nkrt0=0
        do 123 l=1,nlnc
         kk=0
         call ispita(acoz)
         IF(INDFOR.EQ.1)
     1   read(iulaz,*) nktlnc(l),lncvz(l),nbvez(l)
         IF(INDFOR.EQ.2)
     1   read(acoz,1100) nktlnc(l),lncvz(l),nbvez(l)
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2102) L,nktlnc(l),lncvz(l),nbvez(l)
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6102) L,nktlnc(l),lncvz(l),nbvez(l)
         ln=lncvz(l)
         nktl=nktlnc(l)
         do 120 i=1,nktl
         if(i.le.nbvez(l)) go to 120
           ii=ii+1
           kk=kk+1
           call ispita(acoz)
           IF(INDFOR.EQ.1)
     1     read(iulaz,*) ncgl(ii),nvez(ii)
           IF(INDFOR.EQ.2)
     1     read(acoz,1101) ncgl(ii),nvez(ii)
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2103) II,ncgl(ii),nvez(ii)
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6103) II,ncgl(ii),nvez(ii)
c
c
c  formiranje podataka za matricu masa za glavne cvorove
c
c
           if(ndin.eq.1) then
                inmkt=0
                if(kt3d.eq.0) then
                    call ispita(acoz)
                    IF(INDFOR.EQ.1)
     1              read(iulaz,*) xmasa,xmi3
                    IF(INDFOR.EQ.2)
     1              read(acoz,1102) xmasa,xmi3
                    xmass(ii,1)=xmasa
                    xmass(ii,2)=xmasa
                    xmass(ii,3)=xmi3
                    IF(ISRPS.EQ.0)
     1              WRITE(IZLAZ,2104) xmasa,xmi3
                    IF(ISRPS.EQ.1)
     1              WRITE(IZLAZ,6104) xmasa,xmi3
                else
                    call ispita(acoz)
                    IF(INDFOR.EQ.1)
     1              read(iulaz,*) xmasa,xmi1,xmi2,xmi3
                    IF(INDFOR.EQ.2)
     1              read(acoz,1102) xmasa,xmi1,xmi2,xmi3
                    xmass(ii,1)=xmasa
                    xmass(ii,2)=xmasa
                    xmass(ii,3)=xmasa
                    xmass(ii,4)=xmi1
                    xmass(ii,5)=xmi2
                    xmass(ii,5)=xmi3
                    IF(ISRPS.EQ.0)
     1              WRITE(IZLAZ,2105) xmasa,xmi1,xmi2,xmi3
                    IF(ISRPS.EQ.1)
     1              WRITE(IZLAZ,6105) xmasa,xmi1,xmi2,xmi3

c                   write(izlaz,*) (xmass(ii,ik),ik=1,6)
                endif 
           endif 
           nv=nvez(ii)
           ncgl1=ncgl(ii)
             
c           
c           if(kk.ne.i) then
c              do 124 k=1,2
c  124         id(ncgl(ii),k)=-1
c           endif
c         mmp0=mmp0+nvez(ii)
         nkrt0=nkrt0+1
         IDOKLE=1+(NV-1)/14
         DO 10 MENJAJ=1,IDOKLE
            IKRENI=1+(MENJAJ-1)*14
            ISTANI=MENJAJ*14
            IF(ISTANI.GT.NV) ISTANI=NV
         call ispita(acoz)
         IF(INDFOR.EQ.1)
     1   read(iulaz,*) (ncvez(ii,j),j=IKRENI,ISTANI)
         IF(INDFOR.EQ.2)
     1   read(acoz,1100) (ncvez(ii,j),j=IKRENI,ISTANI)
   10    CONTINUE
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2106) (ncvez(ii,j),j=1,nv)
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6106) (ncvez(ii,j),j=1,nv)
c       
c        odre|ivanje cvorova veza izme|u krutih tela vektor ivez
c
         if(mxrac.eq.0.or.kt3d.eq.1) go to 513
         if(i.gt.1) then
             iv=ncvez(ii,1)
             l11=icv
             do 512 l1=1,l11
             iv1=ivez(l1)
             if(iv.eq.iv1) go to 513
             if(l1.eq.l11) then
                ivez(icv)=iv
                mxvc=icv
                icv=icv+1
             endif
  512        continue
         endif
  513    continue
c
c
c
         do 121 j=1,nvez(ii)
            do 121 k=1,n3d
               korg=id(ncvez(ii,j),k)
               if(korg.eq.0.or.korg.eq.-1) then
                    id(ncvez(ii,j),k)=-1
                    mmp0=mmp0+1
c                    write(izlaz,*) ncvez(ii,j),k,mmp0
                    if(mmpogr.eq.0) then
                         if(k.gt.3.and.korg.eq.0) then
                             mmpogr=1
                         endif
				   endif
                endif
  121    continue
         if(kt3d.eq.0) then
             iogx=id(ncgl(ii),1)
             iogy=id(ncgl(ii),2)
             if(i.eq.nktl.and.iogx.eq.1.and.iogy.eq.1) then
                 id(ncgl(ii),1)=-2
                 id(ncgl(ii),2)=-2
             elseif(i.gt.1.and.iogx.eq.0.and.iogy.eq.0) then
                 id(ncgl(ii),1)=-1
                 id(ncgl(ii),2)=-1
             elseif(i.eq.1) then
                 go to 120
             else
             stop '** ULAZEkt ** greska u ogranicenjima krutih tela ***'
             endif
         endif
  120    continue
         nezav0=max(nezav0,nktlnc(l))
  123    continue
         nezav=nezav0+1
cc
cc
         if(kt3d.eq.1) nezav=nezav+1
         mmp=mmp0
c         if(mmpogr.eq.1) mmp=mmp0*2
cc
cc
         if((nkrt0-nkrt).ne.0) stop '*** nkrt.ne.nkrt0 ***'
C      if(mxrac.gt.0) write(izlaz,612) (ivez(ij),ij=1,mxvc)
      endif
 1101 format(2i5,2d15.4)
 1102 format(4d15.5)
 1100 format(14i5)
c  612 format(///' vektor veznih cvorova ivez'/10i5)
      NMPC=0
      DO 60 N=1,NP
      DO 60 J=1,6
      IF(IOPGL(J).EQ.1) THEN
        ID(N,J)=0
        GO TO 60
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
 1000 FORMAT(I5,A1,6I2,2X,3F10.2,I5,I5)
 5010 FORMAT(1X,I5,4X,A1,6I3,2X,3F12.6,I5,I5)
 5030 FORMAT(1X,I5,6I5,3F12.6,I5)
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
 2102 FORMAT(///7X,
     1'L A N A C................................................L =',I5/
     111X,'BROJ KRUTIH TELA U LANCU........................NKTLNC =',I5/
     111X,'BROJ KRUTOG TELA ZA KOJE JE VEZAN LANAC..........LNCVZ =',I5/
     111X,'BROJ ZAJEDNICKIH KRUTIH TELA ZA POVEZANE LANCE.  NBVEZ =',I5)
 2103 FORMAT(//
     111X,'KRUTO TELO..........................................II =',I5/
     111X,'GLAVNI CVOR KRUTOG TELA /TEZISNI/.................NCGL =',I5/
     111X,'BROJ KONTURNIH CVOROVA NA KRUTOM TELU.............NVEZ =',I5)
 2104 FORMAT(/
     111X,'MASA KRUTOG TELA..............................XMASA =',D14.4/
     111X,'MOMENT INERCIJE ZA Z-OSU.......................XMI3 =',D14.4)
 2105 FORMAT(/
     111X,'MASA KRUTOG TELA..............................XMASA =',D14.4/
     111X,'MOMENT INERCIJE ZA X-OSU.......................XMI1 =',D14.4/
     111X,'MOMENT INERCIJE ZA Y-OSU.......................XMI2 =',D14.4/
     111X,'MOMENT INERCIJE ZA Z-OSU.......................XMI3 =',D14.4)
 2106 FORMAT(/'K O N T U R N I    C V O R O V I'/11X,14I5)
 2020 FORMAT(1X,
     1'G E N E R I S A N I   P O D A C I   O   C V O R O V I M A'/
     11X,57('-')///
     11X,' CVOR  J E D N A C I N A    B R O J     K  O  O  R  D  I  N  A
     1  T  E    KOSI'/
     11X,' BROJ                                    X           Y        
     1   Z      SIST')
 2100 FORMAT(///' BROJ N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
     1,I5,' I KONACNOG BROJA KORAKA KORA=',I5)
 2200 FORMAT(///' UCITANI CVOR BROJ N =',I5,' VECI JE OD UKUPNOG BROJA
     1CVORNIH TACAKA NP =',I5)
C-----------------------------------------------------------------------
 6000 FORMAT(//////1X,
     1'I N P U T   N O D A L   D A T A'/
     11X,31('-')///
     11X,' NODE COOR    CONSTRAINTS        C  O  O  R  D  I  N  A  T  E
     1 S   STEP SKEW'/
     11X,' NUMB SYST NX NY NZ WX WY WZ        X           Y           Z
     1     GENE SYST')
 6102 FORMAT(///7X,
     1'C H A I N................................................L =',I5/
     111X,'NUMBER OF RIGID BODIS IN CHAIN..................NKTLNC =',I5/
     111X,'NR. OF THE RIGID BODY THAT CHAIN IS CONNECT. TO.NLNCVZ =',I5/
     111X,'NR. OF COMMON RIGID BODIES FOR CONNECTED CHAINS..NBVEZ =',I5)
 6103 FORMAT(/
     111X,'RIGID BODY..........................................II =',I5/
     111X,'MAIN NODE OF THE RIGID BODY.......................NCGL =',I5/
     111X,'NUMBER OF BOUNDARY NODES OF THE RIGID BODY........NVEZ =',I5)
 6104 FORMAT(/
     111X,'MASS OF THE RIGID BODY........................XMASA =',D14.4/
     111X,'MOMENTUM OF INERTIA FOR Z-AXIS.................XMI3 =',D14.4)
 6105 FORMAT(/
     111X,'MASS OF THE RIGID BODY........................XMASA =',D14.4/
     111X,'MOMENTUM OF INERTIA FOR X-AXIS.................XMI1 =',D14.4/
     111X,'MOMENTUM OF INERTIA FOR Y-AXIS.................XMI2 =',D14.4/
     111X,'MOMENTUM OF INERTIA FOR Z-AXIS.................XMI3 =',D14.4)
 6106 FORMAT(11X,'B O U N D A R Y    N O D E S'/11X,14I5)
 6020 FORMAT(1X,
     1'G E N E R A T E D   N O D A L   D A T A'/
     11X,39('-')///
     11X,' NODE  E Q U A T I O N  N U M B E R   C  O  O  R  D  I  N  A  
     1T  E  S   SKEW'/
     11X,' NUMB                                    X           Y        
     1   Z      SYST')
 6100 FORMAT(///' NUMBER N2=',I5,' HAS NOT OBTAIN BY THE SUM OF NUBMBER
     1N1 =',I5,'AND FINITE NUMBER OF STEPS KORA=',I5)
 6200 FORMAT(///' NODE NUMBER N =',I5,' IS GREATER THAN TOTAL NUMBER OF
     1NODAL POINT NP =',I5)
C-----------------------------------------------------------------------
      END
c
c
c     program za izracunavanje broja korekcija
      subroutine brkorekc(id,np,nkrt,ncgl,ikor)
      dimension id(np,*),ncgl(*)
c
c
      ikor=0
      do 10 i=1,nkrt
          kt=ncgl(i)
          ix=id(kt,1)
          iy=id(kt,2)
          if(ix.eq.-2.and.iy.eq.-2) then
             ikor=ikor+1
          endif
   10 continue 
      return
      end
      subroutine cmpc3d(cord,ncgl,nvez,ncvez,izlaz,np,nezav,
     1                  cmpc,mmp,ro)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /krut3d/ kt3d,mmpogr
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension ncgl(*),nvez(*),ncvez(nkrt,*),cord(np,*),cmpc(mmp,*),
     1          ro(*)
c     nezav --- maksimalan broj nezavisnih u linearnim vezama 
c     mmp   --- ukupan broj razli~itih linearnih veza
c	 ro    --- vektor relativnih koordinata konturni cvorova
c
c   formiranje matrice cmpc
c
      irb=0
c
      do 10 ii=1,nkrt
c
         ncg=ncgl(ii)
         nbkcv=nvez(ii)
c   koordinate tezista krutog tela  
         xc=cord(ncg,1) 
         yc=cord(ncg,2) 
         zc=cord(ncg,3) 
c
         do 20 jj=1,nbkcv
c
            nkcv=ncvez(ii,jj)
c   koordinate konturnih cvorova krutog tela 
            xk=cord(nkcv,1) 
            yk=cord(nkcv,2) 
            zk=cord(nkcv,3) 
            ro(1)=xk-xc
            ro(2)=yk-yc
            ro(3)=zk-zc
c
            do 30 i=1,6
c      
               if(i.gt.3.and.mmpogr.eq.0) go to 30
c
               irb=irb+1
               cmpc(irb,1)=1.d0
c
               if(i.gt.3) go to 50

               do 40 j=2,3
                  i1=i+j-1
                  if(i1.gt.3) i1=i1-3
                  cmpc(irb,j)=-ro(i1)*(-1)**j
   40          continue
               go to 30
c
c
   50          do 60 j=2,3
                     cmpc(irb,j)=0.d0
   60          continue
c
c
   30       continue
   20    continue
   10 continue
c
c
c
c      write(izlaz,121)
c  121 format(///' matrica cmpc'/)
c      do 700 i=1,mmp
c  700 write(izlaz,125) i,(cmpc(i,j),j=1,nezav)
c  125 format(i5,8f10.1)
      return
      end
C
C
C
C
      subroutine eltmkt(ske,lm,ncgl,ncv1,id,xmass,np,nkrt,izlaz,i2)
      implicit double precision (a-h,o-z)
      include 'paka.inc'
      
      common /ielmkt/ inmkt
      dimension ske(*),lm(*),id(np,6),xmass(nkrt,*),ncgl(*)
      inmkt=1
      nkt=ncgl(i2)
      do 5 i=1,ncv1
         i1=i
         if(i.eq.ncv1) i1=i+3
    5 lm(i)=id(nkt,i1)
C      write(izlaz,112) (lm(i),i=1,ncv1)
C  112 format(' vektor lm u eltmkt'/3i5)
      ij=0
      do 10 i=1,ncv1
      do 20 j=i,ncv1
      ij=ij+1
      ske(ij)=0.d00
      if(i.eq.j) then
           amass=xmass(i2,i)
           ske(ij)=amass
c           write(izlaz,*) i,j,ij,amass
      endif
   20 continue
   10 continue
c      write(izlaz,1201) i2
c 1201 format(//' matrica masa elementa ',i3)
      ij1=1
      ij2=0
      do 1100 i1=1,ncv1
      ij=0
      do 1200 j1=i1,ncv1
      ij=ij+1
 1200 continue
      ij2=ij2+ij
C      write(izlaz,1010) i1,(ske(ij3),ij3=ij1,ij2)
      ij1=ij1+ij
 1100 continue
C 1010 format(i2,3(d10.3,2x))      
      return
      end
C
C
      subroutine eltmkt3d(ske,lm,ncgl,ncv1,id,xmass,np,nkrt,izlaz,i2)
      implicit double precision (a-h,o-z)
      include 'paka.inc'
      
      common /ielmkt/ inmkt
      dimension ske(*),lm(*),id(np,6),xmass(nkrt,*),ncgl(*)
      inmkt=1
C
C     formiranje vektora LM
C 
            do 250 j=1,6
               ii=ncgl(i2)
               k=id((ii),j)
  250          lm(j)=k
C
C     formiranje matrice masa
C 
c           write(izlaz,*) (xmass(i2,ik),ik=1,6)
            ij=0
            do 300 ii=1,6                                                          
               do 300 jj=ii,6
                  ij=ij+1
                  ske(ij)=0.d00
                  if(ii.eq.jj)  then
                       amass=xmass(i2,ii)
                       ske(ij)=amass
                  endif
  300       continue
C
C     pakovanje u matricu sistema
C
c      write(izlaz,112) (lm(i1),i1=1,6)
c  112 format(' vektor lm u eltmkt3d'/6i5)
c      write(izlaz,1201) i2
c 1201 format(//' matrica masa elementa ',i3)
      ij1=1
      ij2=0
      do 1100 i1=1,6
      ij=0
      do 1200 j1=i1,6
      ij=ij+1
 1200 continue
      ij2=ij2+ij
c      write(izlaz,1010) i1,(ske(ij3),ij3=ij1,ij2)
      ij1=ij1+ij
 1100 continue
c 1010 format(i2,6(d10.3,2x))      
      return
      end
C
      subroutine accmpc(cord,nktlnc,lncvz,nbvez,ncgl,nvez,
     1                 ncvez,izlaz,np,ac,nkrt2,nezav,cmpc,mmp)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension nktlnc(*),lncvz(*),nbvez(*),ncgl(*),nvez(*),
     1          ncvez(nkrt,*),cord(np,*),ac(nkrt2,*),cmpc(mmp,*)
c     nezav --- maksimalan broj nezavisnih u linearnim vezama 
c     mmp --- ukupan broj razli~itih linearnih veza
c
c
c   formiranje matrice cmpc
c
      ii=0
      irb=1
      do 100 l=1,nlnc
      lnv=lncvz(l)
      nktl=nktlnc(l)
      do 80 i=1,nktl
      ikol=i+1
      if(i.le.nbvez(l)) go to 80
      ii=ii+1
      iibr=nvez(ii)
      do 90 m=1,iibr
      do 95 j=1,nezav
      cmpc(irb,j)=ac(ii,j)
      cmpc(irb+1,j)=ac(ii+nkrt,j)
      if(m.eq.1.and.i.gt.1) go to 95
      if(j.eq.ikol) then
         acoefx=-(cord(ncvez(ii,m),2)-cord(ncgl(ii),2))
         acoefy=(cord(ncvez(ii,m),1)-cord(ncgl(ii),1))
         cmpc(irb,j)=ac(ii,j)+acoefx
         cmpc(irb+1,j)=ac(ii+nkrt,j)+acoefy
      endif
   95 continue
      irb=irb+2
   90 continue
   80 continue
  100 continue
c      write(izlaz,121)
C  121 format(///' matrica cmpc'/)
C      do 700 i=1,mmp
C  700 write(izlaz,125) i,(cmpc(i,j),j=1,nezav)
C  125 format(i3,7f10.5)
      return
      end
C
      subroutine formac(cord,nktlnc,lncvz,nbvez,ncgl,nvez,
     1                 ncvez,izlaz,np,ac,nkrt2,nezav)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      common /ielmkt/ inmkt
      dimension nktlnc(*),lncvz(*),nbvez(*),ncgl(*),nvez(*),
     1          ncvez(nkrt,*),cord(np,*),ac(nkrt2,*)
c     nezav --- maksimalan broj nezavisnih u linearnim vezama 
c     mmp --- ukupan broj razli~itih linearnih veza
c
c
c     formiranje pomocne matrice ac sa koeficijentima uz pomeranja tezi
c     sta krutih tela i to prvo x pa y
c
c      print*,'usao u krttel '
      inmkt=0
      irb=0
      do 5 k=2,1,-1
      ii=0
      do 10 l=1,nlnc
      lnv=lncvz(l)
      nktl=nktlnc(l)
      do 20 i=1,nktl
      ikol=i+1
      if(i.le.nbvez(l)) go to 20
      ii=ii+1
      irb=irb+1
      ac(irb,1)=1.d0
      do 40 j=2,nezav
      ac(irb,j)=0.d0
   40 continue
c      print*,irb,j,ac(irb,j),ii
      if(i.eq.1) go to 1000
      if(i.eq.(nbvez(l)+1).and.lncvz(l).ne.0) then
c     sracunavanje repera za prve clanove lanca koji je vezan za drugi
cc         irep=0
cc         do 60 i1=1,lnv
cc         irep1=nktlnc(i1)-nbvez(i1)
cc   60    irep=irep+irep1
cc         irep=irep-(nktlnc(lncvz(l))-nbvez(l))
c     sracunavanje koeficijenata za prvo telo vezano za drugi lanac
         do 70 i2=2,ikol
         irep=lncvz(l)
         irep1=irep
         if(k.eq.1) irep1=irep+nkrt
         ac(irb,i2)=ac(irep1,i2)
         if(i2.eq.i) then
            acoef1=(-cord(ncvez(ii,1),k)+cord(ncgl(irep),k))
            if(k.eq.1) acoef1=-acoef1
            ac(irb,i2)=ac(irb,i2)+acoef1
         elseif(i2.eq.i+1) then
            ac(irb,i2)=cord(ncvez(ii,1),k)-cord(ncgl(ii),k)
            if(k.eq.1) ac(irb,i2)=-ac(irb,i2)
         endif
   70   continue
        go to 1000
      endif
      do 30 j=2,ikol
      if(j.lt.i) then
          ac(irb,j)=ac(irb-1,j)
       elseif(j.eq.i) then
          ac(irb,j)=-(cord(ncvez(ii,1),k)-cord(ncgl(ii-1),k))
          if(k.eq.1) ac(irb,j)=-ac(irb,j)
          ac(irb,j)=ac(irb,j)+ac(irb-1,j)
       else
          ac(irb,j)=cord(ncvez(ii,1),k)-cord(ncgl(ii),k)
          if(k.eq.1) ac(irb,j)=-ac(irb,j)
      endif
c      if(k.eq.1) ac(irb,j)=-ac(irb,j)
c      write(izlaz,*) irb,k,i,j,ac(irb,j)
   30 continue
 1000 continue
c      stop '** iza petlje 30'
   20 continue
   10 continue
    5 continue
c      write(izlaz,1112)
c 1112 format(///' matrica ac'//)
      do 1115 i=1,nkrt2
c      write(izlaz,125) (ac(i,j),j=1,nezav)
 1115 continue
c  125 format(6f10.2)
c      stop '*** zaustavljeno iza ac ***'
      return
      end
      subroutine formim(id,np,nktlnc,lncvz,nbvez,ncgl,nvez,
     1                 ncvez,izlaz,impc,jmpc,nezav)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension nktlnc(*),lncvz(*),nbvez(*),ncgl(*),nvez(*),
     1          ncvez(nkrt,*),id(np,*),impc(nkrt,*),jmpc(*)
c
c    formiranje pomocne matrice impc i vektora jmpc
c
c
      ki=0
      ii=0
      do 120 l=1,nlnc
      lnv=lncvz(l)
      nktl=nktlnc(l)
      kk=0
      do 130 i=1,nktl
      kt=i+1
      if(i.le.nbvez(l)) go to 130
      kk=kk+1
      ii=ii+1
c        odredjivanje repera za impc
      if(i.eq.(nbvez(l)+1).and.lncvz(l).ne.0) irep=lncvz(l)
cc           irep=0
cc           do 200 i1=1,lnv
cc           irep1=nktlnc(i1)-nbvez(i1)
cc  200      irep=irep+irep1
cc           irep=irep-(nktlnc(lncvz(l))-nbvez(l))
cc           irep2=irep
cc           if(k.eq.1) irep2=irep+nkrt
cc       endif
       do 160 ij=1,nezav

       impc(ii,ij)=0
       if(ij.le.(kt)) impc(ii,ij)=id(ncgl(ii),6)
       if(kk.eq.1.and.lncvz(l).eq.0) then
            if(ij.eq.1) then
                impc(ii,ij)=id(ncgl(ii),1)
                jmpc(ii)=id(ncgl(ii),2)
            endif
        endif
        if(kk.gt.1) then
            if(ij.le.i) impc(ii,ij)=impc(ii-1,ij)
            jmpc(ii)=jmpc(ii-1)
        endif
        if(kk.eq.1.and.lncvz(l).ne.0) then
            if(ij.le.i) then
                 impc(ii,ij)=impc(irep,ij)
                 jmpc(ii)=jmpc(irep)
            endif
        endif
  160 continue
  130 continue
  120 continue
c      write(izlaz,124)
c      do 500 i=1,nkrt
c      write(izlaz,*) (impc(i,j),j=1,nezav)      
c  500 write(izlaz,*) jmpc(i)
c  124 format(///' matrica impc'/)
      return
      end

      subroutine formirac(id,ncgl,nvez,ncvez,izlaz,np,ivez,iracv)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /ispakk/ liracv
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension ncgl(*),nvez(*),ncvez(nkrt,*),id(np,*),iracv(mxvc,*),
     1          ivez(*)
c
c
c
c     formiranje matrice iracv sa glavnim cvorovima oko cvora racvanja
c
c      write(izlaz,5501)
c 5501 format(' matrica IRACV')
      do 510 ir=1,mxvc
      do 520 jjr=1,mxrac
  520 iracv(ir,jjr)=0
      irv=ivez(ir)
      jr=0
      do 510 iik=1,nkrt
      iik1=iik
      iikv=ncgl(iik)
      lcvor=id(iikv,6)
      nv=nvez(iik)
      do 540 kr=1,nv
      krv=ncvez(iik,kr)
      if(irv.eq.krv) then
          jr=jr+1
          iracv(ir,jr)=lcvor
c          write(izlaz,5100) iik,irv,krv,ir,jr
c 5100 format(' iik=',i3,'   irv=',i3,'  krv=',i3,'  ir=',i3,'  jr=',i3)
      endif
  540 continue
c      if(jr.gt.mxrac) stop ' **** mxrac ****',mxrac
  510 continue
c      do 6000 ir1=1,mxvc
c 6000 write(izlaz,5500) (iracv(ir1,jrr),jrr=1,mxrac)
c 5500 format(10i5)
c
c
c

      return
      end
      subroutine formmpc(id,np,nktlnc,nbvez,ncgl,nvez,
     1                 ncvez,izlaz,impc,jmpc,nezav,mpc,neza1,mmp)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension nktlnc(*),nbvez(*),ncgl(*),nvez(*),ncvez(nkrt,*),
     1          id(np,*),impc(nkrt,*),jmpc(*),mpc(neza1,*)
c
      ki=0
      ii=0
      do 120 l=1,nlnc
      nktl=nktlnc(l)
      kk=0
      do 130 i=1,nktl
      if(i.le.nbvez(l)) go to 130
      kk=kk+1
      ii=ii+1
      iibr=nvez(ii)
c
c
c     formiranje matrice mpc
c
      do 210 j=1,iibr
      do 220 k=1,2
      ki=ki+1
      mpc(1,ki)=ki
      do 230 ij=1,nezav
      mpc(ij+1,ki)=impc(ii,ij)
  230 continue
      if(k.eq.2) mpc(2,ki)=jmpc(ii)
c
c     formiranje negativnih jednacina za zavisna pomeranja
c
      if(i.ne.1.and.j.eq.1) then
          id(ncgl(ii),k)=-ki
      else
          id(ncvez(ii,j),k)=-ki
      endif
  220 continue
  210 continue
  130 continue
  120 continue
c
c

C      write(izlaz,122)
C      do 400 j=1,mmp
C  400 write(izlaz,124) (mpc(i,j),i=1,neza1)      
C      write(izlaz,123)
c      do 600 i=1,np
C  600 write(izlaz,124) i,(id(i,j),j=1,6)
C  122 format(///' matrica mpc'/)
C  123 format(///' korigovana matrica id'/)
C  124 format(20i3)
      return
      end
c
c
c     program za korigovanje id mpc cmpc
      subroutine korig(id,np,ijogr,ik,ikor,mpc,cmpc,izlaz)
      implicit double precision(a-h,o-z)
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      dimension id(np,*),ijogr(*)
      dimension mpc(neza1,*),cmpc(mmp,*)
c
c
      nc=ijogr(ik)
      negx=-id(nc,1)
      negy=-id(nc,2)
      id(nc,1)=0
      id(nc,2)=0
      ifn=0
      kk=0
      do 10 i=2,neza1
         j=mpc(i,negx)
         if(j.eq.0) go to 10
         ifn1=ifn
         ifn=mpc(i,negx)
         ir=i
         kk=kk+1
   10 continue
         if(kk.lt.2) stop ' *** korig *** =1'
         k=ir
      do 20 i=1,np
         if11=id(i,6)
         if(ifn1.eq.if11) then
             nc1=i
             go to 30
          endif
   20 continue 
c
c     korekcija matrice id
   30 id(nc1,6)=-negx
      id(nc,6)=-negy
c
c     korekcija matrice mpc izbacivanje jednacina za fn i fn1
      do 40 i=k-1,k
         mpc(i,negx)=0
   40    mpc(i,negy)=0
c
c     odredjivanje glavnih koeficijenata za korekciju
      ckn1=cmpc(negx,k-2)
      ckn=cmpc(negx,k-1)
      cn1=cmpc(negy,k-2)
      cn=cmpc(negy,k-1)
      isim=0
      aaa=(cn1-ckn1/ckn*cn)
      if(abs(aaa).lt.1.d-10) isim=1
      ckn=ckn+1.d-05
      aaa=(cn1-ckn1/ckn*cn)
      irep9=k-3
c      if(isim.eq.0) go to 100
c      do 105 j=1,nezav
c        ai=0.d00
c        bi=0.d00
c        mpx=mpc(j+1,negx)
c        mpy=mpc(j+1,negy)
c        mpc(j+1,negx)=0
c        mpc(j+1,negy)=0
c        if(j.eq.3) then
c           ai=1.
c           mpc(j+1,negx)=mpx
c        endif
c        if(j.eq.2) then
c           bi=1.
c           mpc(j+1,negy)=mpy
c        endif
c        cmpc(negx,j)=ai
c        cmpc(negy,j)=bi
c  105 continue
      aa=-1./aaa
      do 50 j=1,nezav
         if(j.eq.1.or.j.gt.irep9) then
             cmpc(negx,j)=0.d00
             cmpc(negy,j)=0.d00
         endif
         cki=cmpc(negx,j)
         ci=cmpc(negy,j)
         ai=(ci-cki/ckn*cn)*aa
         bi=-(cki+ai*ckn1)/ckn
         cmpc(negx,j)=ai
         cmpc(negy,j)=bi
   50 continue
      do 55 kor=1,2
         if0=ifn1
         if(kor.eq.2) then
             if0=ifn
             negx=negy
          endif
          do 60 j=1,mmp
             koef=mpc(1,j)
             do 70 i=2,neza1
                k1=mpc(i,j)
                if(k1.eq.if0) then
                     mpc(i,j)=0
                     jlok=i-1
                     cnn=cmpc(koef,jlok)
                     cmpc(koef,jlok)=0.d00
                     cmpc(koef,1)=0.d00
                     do 80 jj=2,jlok-1
                        ai=cmpc(negx,jj)
                        cmpc(koef,jj)=cmpc(koef,jj)+cnn*ai
   80                continue
                     go to 60
                endif
   70        continue
   60     continue
   55 continue
C      write(izlaz,122)
C      do 400 j=1,mmp
C  400 write(izlaz,124) (mpc(i,j),i=1,neza1)      
C      write(izlaz,123)
C      do 600 i=1,np
C  600 write(izlaz,124) i,(id(i,j),j=1,6)
C  122 format(///' matrica ** korig ** mpc'/)
C  123 format(///' korigovana matrica posle ** korig ** id'/)
C      write(izlaz,121)
C  121 format(///' matrica cmpc ** korig **'/)
C      do 700 i=1,mmp
C  700 write(3,125) i,(cmpc(i,j),j=1,nezav)
C  125 format(i5,8f10.1)
C  124 format(20i3)
      return
      end
      subroutine mpc3d(id,np,ncgl,nvez,
     1                 ncvez,izlaz,nezav,mpc,neza1,mmp,nid)
      implicit double precision(a-h,o-z)
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /krut3d/ kt3d,mmpogr
      COMMON /DUPLAP/ IDVA
      COMMON /RADIZA/ INDBG
      dimension ncgl(*),nvez(*),ncvez(nkrt,*),id(np,*),mpc(neza1,*),
     1          nid(*)
c
c   formiranje matrice mpc
c
      irb=0
c
      do 10 ii=1,nkrt
c
         k=3
         ncg=ncgl(ii)
         nbkcv=nvez(ii)
         nid(1)=id(ncg,4)
         nid(2)=id(ncg,5)
         nid(3)=id(ncg,6)
c
         do 20 jj=1,nbkcv
c
c   koordinate konturnih cvorova krutog tela 
c
            do 30 i=1,6
c      
               if(i.gt.3.and.mmpogr.eq.0) go to 30
c
               k2=ncvez(ii,jj)
               k3=id(k2,i)
               irb=irb+1
               mpc(1,irb)=irb
               mpc(2,irb)=id(ncg,i)
c
               if(i.gt.3.and.k3.eq.-1) go to 50
c
               if(i.gt.3.and.k3.eq.0) then
                     mpc(2,irb)=0
                     go to 50
               endif
c
c
               do 40 j=3,neza1
c
                  if(k.lt.1) k=3
                  k1=nid(k)
                  mpc(j,irb)=k1
                  k=k-1
   40          continue
               go to 60
c
   50          do 70 j=3,neza1
                     mpc(j,irb)=0
   70          continue
c

c           formiranje negativnih jednacina u matrici id
c
   60       id(k2,i)=-irb
c
   30       continue
   20    continue
   10 continue
c
c
C      write(izlaz,122)
C      do 400 j=1,mmp
C  400 write(izlaz,*) (mpc(i,j),i=1,neza1)      
C      write(izlaz,123)
C      do 600 i=1,np
C  600 write(izlaz,*) i,(id(i,j),j=1,6)
C  122 format(///' matrica mpc'/)
C  123 format(///' korigovana matrica id'/)
      return
      end

c
c
c     program za odredjivane cvorova koji su ograniceni sa -2 
c     /zadnji cvorovi u lancu, ako su ogranicena pomeranja u x i y/
      subroutine ogrcvor(id,np,ijogr,ikor)
      dimension id(np,*),ijogr(*)
c
c
      ik=0
      do 10 i=1,np
          ix=id(i,1)
          iy=id(i,2)
          if(ix.eq.-2.and.iy.eq.-2) then
             ik=ik+1
             ijogr(ik)=i
          endif
   10 continue 
      if(ik.ne.ikor) stop '*** ogrcvor **** greska u ogran. krutih tela'
      return
      end
C=======================================================================
C
C=======================================================================
      SUBROUTINE SPAKKC(SK,MAXA,IRACV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   P R O G R A M
CS.      ZA RAZMESTANJE KONSTANTNIH KRUTOSTI U CVOROVIMA U SISTEM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      common /ispakk/ liracv
      COMMON /CDEBUG/ IDEBUG
      DIMENSION SK(*),MAXA(*),IRACV(MXVC,*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKKC'
C
C
      skkc=0.d00
      skkc1=0.d00
C      do 6000 ir1=1,mxvc
C 6000 write(izlaz,5500) (iracv(ir1,jrr),jrr=1,mxrac)
C 5500 format(10i5)
C      write(izlaz,1001)
C 1001 format(' kv  i  j ii jj kk   sk')
      DO 10 KV=1,MXVC
          MNQ0=1
          MNQ1=JEDN
          MXMN=0
        DO 20 I=1,MXRAC
            II=IRACV(KV,I)
            IF(II.LT.MNQ0) GO TO 20
            MI=MAXA(II)-MXMN
          DO 30 J=1,MXRAC
             JJ=IRACV(KV,J)
             IF(JJ) 35,30,40
   35        STOP 'JEDNACINA NE MOZE BITI MANJA OD NULE'
   40        IJ=II-JJ
             IF(IJ) 30,50,60
   50        KK=MI
c             SK(KK)=SK(KK)+SKKC
             SK(KK)=SKKC
C             write(izlaz,1000) kv,i,j,ii,jj,kk,sk(kk)
             GO TO 30
   60        KK=MI+IJ
c             SK(KK)=SK(KK)-SKKC
             SK(KK)=-SKKC1
C             write(izlaz,1000) kv,i,j,ii,jj,kk,sk(kk)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C 1000 format(6i3,f10.3)
C      write(izlaz,125)
C      do 524 isk=1,jedn
C      isk1=maxa(isk)
C      isk2=maxa(isk+1)-1
C  524 write(izlaz,126) isk,(sk(iskk),iskk=isk1,isk2)
C  125 format(' matrica  S K     od dijagonale navise u SPAKKC') 
C  126 format(i2,3(10d12.4/2x)/)
      RETURN
      END
             
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAKR(SK,MAXA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   P R O G R A M
CS.      ZA RAZMESTANJE KONSTANTNIH KRUTOSTI U CVOROVIMA U SISTEM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION SK(*),MAXA(*)
      IF(IDEBUG.GT.0) PRINT *, ' STAKR'
C
C
C      write(izlaz,125)
C  125 format(' matrica  S K     od dijagonale navise u STAKR') 
C      do 523 isk=1,jedn
C      isk1=maxa(isk)
C      isk2=maxa(isk+1)-1
C  523 write(izlaz,126) isk,(sk(iskk),iskk=isk1,isk2)
C  126 format(i2,3(10d12.4/2x)/)
      return
      END
             
C=======================================================================
C=======================================================================
      SUBROUTINE TGRAFKT(NCGL,NVEZ,NCVEZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PROGRAM ZA STAMPANJE KRUTIH TELA  U UNIVERZALNI FILE
C
      common /krutot/ nlnc,nkrt,mxvez,lnkt,llncvz,lnbv,lncgl,lnvez,
     1                lncvez,lac,limpc,ljmpc,livez,mxvc,mxrac,lxmass
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION NCGL(*),NVEZ(*),NCVEZ(NKRT,*)
      ISUMGR=ISUMGR+1
C     GRAFICKI OPIS KRUTOG TELA = 34
      IFGD=34
C     VRSTA  ELEMENTA: KRUTO TELO  = 1
      IFDI=1
C     TABELA FIZICKIH OSOBINA
      IPTN=ISUMGR
C     TABELA MATERIJALA
      MPTN=ISUMGR
C     BOJA  
      ICOL=7
      IND=-1
      ITYP=71
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 I=1,NKRT
C        REDNI BROJ ELEMENTA
         IEL=I+ISUMEL
C        BROJ CVOROVA NA ELEMENTU
         NNODS=NVEZ(I)+1
         nn1=nnods-1
         WRITE(IGRAF,1000) IEL,IFGD,IFDI,IPTN,MPTN,ICOL,NNODS
         WRITE(IGRAF,1000) NCGL(I),(NCVEZ(I,J),J=1,NN1)
   10 CONTINUE
      WRITE(IGRAF,1100) IND
      ISUMEL=ISUMEL+NKRT
      RETURN
C
 1000 FORMAT(8I10)
 1100 FORMAT(I6)
      END



