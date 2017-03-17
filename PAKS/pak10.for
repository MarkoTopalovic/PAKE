C=======================================================================
C
CE       DINAMICS ANALYSIS - GENERAL SUBROUTINE
CS       DINAMIKA - OPSTI POTPROGRAMI
C
C   SUBROUTINE UCPUVA
C              UPUVA
C              INTKMM
C              INTKM
C              INTKMC
C              FORMEF
C              RPRUVA
C              VILNEW
C              CENRAZ
C              VILSON
C              NEWMAR
C              CENDIF
C              CENPOC
C              CONSTD
C              STAMP
C              KSTAMP
C              SAVREZ
C              SILAAN
C
C=======================================================================
      SUBROUTINE INTKMM
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   LOOP OVER GROUP OF ELEMNTS TO INTEGRATION MATRICES M AND C
CS.   PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE M I C
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GRUPER/ LIGRUP
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /EXPLICITNA/ INDEXPL
      COMMON /EXPLALOK/ IUBRZ,IBRZINA,IBRZINAIPO,IPOMAK,IMASA,IPRIGUSEN
C
      IF(IDEBUG.GT.0) PRINT *, ' INTKMM'
      NWM=NWK
      IF(IMASS.EQ.2) NWM=JEDN
CSKDISK....
      LSKP=LSK
      LSK =LSKP+NWK*IDVA
CSKDISK
      IF(ITEST.EQ.0) GO TO 10
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010)
      CALL ISPITA(ACOZ)
      CALL READTE(A(LSK),NWM)
      GO TO 20
   10 CALL CLEAR(A(LSK),NWM)
      CALL INTKM(A(LIGRUP))

      IF(INDEXPL.EQ.1) CALL JEDNA1(A(IMASA),A(LSK),JEDN)
C
CE    STORE LINEAR MATRIX M ON UNIT 11
CS    ZAPISIVANJE LINEARNE MATRICE M NA DISK 11
C
   20 IF(IMASS.NE.2) CALL WSTAZK(A(LIPODS),LSK,54)
      IF(IMASS.EQ.2) CALL WSTAZ(A(LIPODS),LSK,54)
C      WRITE(3,*) 'LSK',LSK
c      CALL WRR(A(LSK),NWM,'MASA')
c      CALL WRR(A(IMASA),JEDN,'IMASA')
CSKDISK....
      LSK =LSKP
CSKDISK
C
CE    LOOP FOR GROUPS OF ELEMENTS TO INTEGRATION MATRIX C
CS    PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE  C
C
      IF(IDAMP.EQ.0) RETURN
      IF(ITEST.EQ.0) GO TO 30
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2020)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6020)
      CALL ISPITA(ACOZ)
      CALL READTE(A(LSK),NWM)
      GO TO 40
   30 CALL CLEAR(A(LSK),NWM)
      CALL INTKMC(A(LIGRUP))
C      IF(INDEXPL.EQ.1) CALL JEDNA1(A(IPRIGUSEN),A(LSK),JEDN)
C
CE    TO STORE LINERAR MATRIX C ON DISK 11
CS    ZAPISIVANJE LINEARNE MATRICE C NA DISK 11
C
   40 IF(IDAMP.NE.2) CALL WSTAZK(A(LIPODS),LSK,56)
      IF(IDAMP.EQ.2) CALL WSTAZ(A(LIPODS),LSK,56)
      RETURN
C-----------------------------------------------------------------------
 2010 FORMAT(///' M A T R I C A   M')
 2020 FORMAT(///' M A T R I C A   C')
C-----------------------------------------------------------------------
 6010 FORMAT(///' M A T R I X   M')
 6020 FORMAT(///' M A T R I X   C')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTKM(IGRUP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   LOOP FOR GROUP OF ELEMENTS TO INTEGRATION MATRIX M
CS.   PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE  M
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IGRUP(NGELEM,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTKM'
      DO 100 NGE = 1,NGELEM
      IATYP = IGRUP(NGE,3)
      LMAX=LRAD
      NETIP = IGRUP(NGE,1)
      NE = IGRUP(NGE,2)
      NMODM = IGRUP(NGE,4)
      LMAX8 = IGRUP(NGE,5)
C
      CALL ELEME(NETIP,5)
C
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTKMC(IGRUP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   LOOP FOR GROUPS OF ELEMENTS TO INTEGRATION MATRIX C
CS.   PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE  C
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IGRUP(NGELEM,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTKMC'
      DO 100 NGE = 1,NGELEM
      IATYP = IGRUP(NGE,3)
      LMAX=LRAD
      NETIP = IGRUP(NGE,1)
      NE = IGRUP(NGE,2)
      NMODM = IGRUP(NGE,4)
      LMAX8 = IGRUP(NGE,5)
C
      CALL ELEME(NETIP,6)
C
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE FORMEF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO FORM EFECTIVE OF RIGHT SIDE VECTOR
CS.   PROGRAM
CS.      ZA FORMIRANJE EFEKTIVNOG VEKTORA DESNE STRANE
C .         (T+DT)REF = (T+DT)R + (T+DT)RM + (T+DT)RC
C .
C ......................................................................
C
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' FORMEF'
      IF(MDVI-2) 10,10,20
   10 CALL VILNEW
      RETURN
   20 CALL CENRAZ
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE RPRUVA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO CALCULATE INCREMENT OF DISPLACEMENTS, NEW VELOCITIES AND
CE.      ACCELERATIONS IN T+DT
CS.   PROGRAM
CS.      ZA RACUNANJE PRIRASTAJA POMERANJA, NOVIH BRZINA I UBRZANJA
CS.      U TRENUTKU T+DT
C .
C ......................................................................
C
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' RPRUVA'
      IF(MDVI-2)  10,20,30
   10 CALL VILSON
      RETURN
   20 CALL NEWMAR
      RETURN
   30 CALL CENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VILNEW
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
C .   PROGRAM
CE.      TO FORM ON RIGHT SIDE VECTOR FROM MASSES AND DAMPING
CS.   PROGRAM
CS.      ZA FORMIRANJE VEKTORA DESNE STRANE USLED MASA I PRIGUSENJA
C .
C .   - WILSON
C .     (T+�T)R=�*(T+T)R+(1-�)*(T)R+M*(A0*(T)UU+A2*(T)UV+A3*(T)UA)+
C .                                  +C*(A1*(T)UU+A4*(T)UV+A5*(T)UA)
C .   - NEWMARK
C .     (T+T)R=(T+T)R+M*(A0*(T)UU+A2*(T)UV+A3*(T)UA)+
C .                    +C*(A1*(T)UU+A4*(T)UV+A5*(T)UA)
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /OPTERE/ NCF,NPP2,NPP3,NPGR,NPGRI,NPLJ,NTEMP
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DUPLAP/ IDVA
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /GRUPER/ LIGRUP
      COMMON /ITERBR/ ITER
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' VILNEW'
      NWM=NWK
      IF(IMASS.EQ.2) NWM=JEDN
      IF(ICONT.GT.0) JEDN=NEQ
CSKDISK
      LSKP=LSK+NWK*IDVA
CSKDISK
      LMAX=LRAD
      LUMC=LMAX
      LMAX=LUMC+JEDN*IDVA
C
CE    FORMIRANJE VEKTORA DESNE STRANE USLED ZADATIH SPOLJ. OPTERECENJA
CS    FORMIRANJE VEKTORA DESNE STRANE USLED ZADATIH SPOLJ. OPTERECENJA	
C
      IF(MDVI.EQ.1) THEN
C       WILSON
        LRT=LMAX
        LMAX=LRT+JEDN*IDVA
C       (T+T)R
        IF(ITERGL.EQ.0) CALL RSTAZ(A(LIPODS),LRTDT,38)
C???    NA 53 NISU ZAPISANE SILE OD NELINEARNIH PRITISAKA
        IF(ITERGL.GT.0.AND.(NPP2.GT.0.OR.NPP3.GT.0.OR.NPLJ.GT.0)) STOP
     1     'STOP - PAK10 - VILNEW - VILSON'
C       (T)R
        CALL RSTAZ(A(LIPODS),LRT,53)
        TETA=PIP
        TETAM=1.0-PIP
C       (T+�T)R=�*(T+T)R+(1-�)*(T)R
        CALL ZBIR2(A(LUMC),A(LRTDT),A(LRT),TETA,TETAM,JEDN)
        CALL JEDNA1(A(LRTDT),A(LUMC),JEDN)
        LMAX=LRT
      ELSE
C       NEWMARK
C       (T+T)R
        IF(ITERGL.EQ.0) CALL RSTAZ(A(LIPODS),LRTDT,38)
C         CALL WRR(A(LRTDT),JEDN,'RUF0')
      ENDIF
C
CE    FORMIRANJE VEKTORA DESNE STRANE USLED MASA
CS    FORMIRANJE VEKTORA DESNE STRANE USLED MASA
C
      LPUU=LMAX
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
      LTDT=LMAX
      IF(ITERGL.GT.0) LMAX=LTDT+JEDN*IDVA
C     (T)UU, (T)UV, (T)UA, (T+DT)UU
      CALL RSTAZ(A(LIPODS),LPUU,59)
      CALL RSTAZ(A(LIPODS),LPUV,55)
      CALL RSTAZ(A(LIPODS),LPUA,57)
C         CALL WRR(A(LPUU),JEDN,'PUU ')
C         CALL WRR(A(LPUV),JEDN,'PUV ')
C         CALL WRR(A(LPUA),JEDN,'PUA ')
      IF(ITERGL.GT.0)THEN
        IF(ITER.EQ.0)CALL JEDNA1(A(LTDT),A(LPUU),JEDN)
        IF(ITER.GT.0)CALL RSTAZ(A(LIPODS),LTDT,52)
      ENDIF
C     UM=A0*(T)UU+A2*(T)UV+A3*(T)UA
      IF(ITERGL.EQ.0)
     1CALL ZBIR3(A(LUMC),A(LPUU),A(LPUV),A(LPUA),A0,A2,A3,JEDN)
      IF(ITERGL.GT.0) THEN
C        UM=A2*(T)UV+A3*(T)UA
         IF(ITER.EQ.0) CALL ZBIR2(A(LUMC),A(LPUV),A(LPUA),A2,A3,JEDN)
C        UM=-A0*(T+DT)UU+A0*(T)UU+A2*(T)UV+A3*(T)UA
         A0M=-A0
         IF(ITER.GT.0)
     1   CALL ZBIR4(A(LUMC),A(LTDT),A(LPUU),A(LPUV),A(LPUA),
     1              A0M,A0,A2,A3,JEDN)
C.. KOREKCIJA UBRZANJA ZA KONTAKTNE PROBLEME
         IF(ICONT.GT.0)THEN
           IDUM=LPUU
           LPUU=LUMC
           CALL CONCOR(A(LIGRUP),7)
           LPUU=IDUM
         ENDIF
      ENDIF
      LMAX=LPUU
C
CSKDISK      LSKP=LMAX
CSKDISK      LMAX=LSKP+NWK*IDVA
CE    READ LINEAR MATRIX M FROM UNIT 11
CS    UCITAVANJE LINEARNE MATRICE M SA DISKA 11
      IF(IMASS.NE.2) CALL RSTAZK(A(LIPODS),LSKP,54)
      IF(IMASS.EQ.2) CALL RSTAZ(A(LIPODS),LSKP,54)
C     (T+T)R=(T+T)R+M*UM
      if (IABS(ICCGG).EQ.1) then 
            IF(IMASS.NE.2) CALL MAXAPRI(A(LSKP),A(LUMC),A(LRTDT),
     1                     JEDN,A(IROWS),A(IROWS+nwk),nwk)
      else
            IF(IMASS.NE.2) CALL MAXAPR(A(LSKP),A(LUMC),A(LRTDT),
     1                     A(LMAXA),JEDN)
      endif
C         CALL WRR(A(LRTDT),JEDN,'RPR ')
C         CALL WRR(A(LUMC),JEDN,'UMC ')
C         WRITE(3,*) 'LSKP,JEDN',LSKP,JEDN
C         CALL WRR(A(LSKP),JEDN,'SKP ')
      IF(IMASS.EQ.2) CALL MNOZMU(A(LRTDT),A(LSKP),A(LUMC),JEDN)
C         CALL WRR(A(LRTDT),JEDN,'RPO ')
CSKDISK      LMAX=LSKP
C
CE    FORM ON RIGHT SIDE VECTOR BECAUSE DAMPERS
CS    FORMIRANJE VEKTORA DESNE STRANE USLED PRIGUSENJA
C
      IF(IDAMP.EQ.0) GO TO 10
      LPUU=LMAX
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
      LTDT=LMAX
      IF(ITERGL.GT.0) LMAX=LTDT+JEDN*IDVA
C     (T)UU, (T)UV, (T)UA, (T+DT)UU
      CALL RSTAZ(A(LIPODS),LPUU,59)
      CALL RSTAZ(A(LIPODS),LPUV,55)
      CALL RSTAZ(A(LIPODS),LPUA,57)
      IF(ITERGL.GT.0.AND.ITER.GT.0) CALL RSTAZ(A(LIPODS),LTDT,52)
C     UC=A1*(T)UU+A4*(T)UV+A5*(T)UA
      IF(ITERGL.EQ.0)
     1CALL ZBIR3(A(LUMC),A(LPUU),A(LPUV),A(LPUA),A1,A4,A5,JEDN)
      IF(ITERGL.GT.0) THEN
C        UC=A4*(T)UV+A5*(T)UA
         IF(ITER.EQ.0)
     1   CALL ZBIR2(A(LUMC),A(LPUV),A(LPUA),A4,A5,JEDN)
C        UC=-A1*(T+DT)UU+A1*(T)UU+A4*(T)UV+A5*(T)UA
         A1M=-A1
         IF(ITER.GT.0)
     1   CALL ZBIR4(A(LUMC),A(LTDT),A(LPUU),A(LPUV),A(LPUA),
     1              A1M,A1,A4,A5,JEDN)
      ENDIF
      LMAX=LPUU
C
CSKDISK      LSKP=LMAX
CSKDISK      LMAX=LSKP+NWK*IDVA
C     UCITAVANJE LINEARNE MATRICE C SA DISKA 11
      IF(IDAMP.NE.2) CALL RSTAZK(A(LIPODS),LSKP,56)
      IF(IDAMP.EQ.2) CALL RSTAZ(A(LIPODS),LSKP,56)
C     (T+T)R=(T+T)R+C*UC
      if (IABS(ICCGG).EQ.1) then 
         IF(IDAMP.NE.2) CALL MAXAPRI(A(LSKP),A(LUMC),A(LRTDT),
     1                       JEDN,A(IROWS),A(IROWS+nwk),nwk)
      else
         IF(IDAMP.NE.2) CALL MAXAPR(A(LSKP),A(LUMC),A(LRTDT),
     &                       A(LMAXA),JEDN)
      endif
      IF(IDAMP.EQ.2) CALL MNOZMU(A(LRTDT),A(LSKP),A(LUMC),JEDN)
C
C     (T+T)R
   10 IF(ITERGL.EQ.0) CALL WSTAZ(A(LIPODS),LRTDT,38)
      LMAX=LUMC
      IF(ICONT.GT.0) JEDN=NEQ+NEQC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CENRAZ
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO FORM ON RIGHT SIDE VECTORS - CENTRAL DIFFERENS
CS.   PROGRAM
CS.      ZA FORMIRANJE VEKTORA DESNE STRANE - CENTRALNE RAZLIKE
C .
C .      (T)R=(T)R-K*(T)UU+(A0*M+A1*C)*(T-T)UU+A2*M*((T)UU-(T-T)UU)
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /DUPLAP/ IDVA
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SCRATC/ ISCRC
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' CENRAZ'
C???  NIJE UGRADJEN NELINEARNI ALGORITAM ZA CENTRALNE RAZLIKE
      IF(ITERGL.GT.0) STOP 'STOP - PAK10 - CENRAZ'
      LSKP=LSK+NWK*IDVA
      LMAX=LRAD
      LUMC=LMAX
      LPUU=LUMC+JEDN*IDVA
      LMAX=LPUU+JEDN*IDVA
C     (T)R
      CALL RSTAZ(A(LIPODS),LRTDT,38)
C     (T)UU
      CALL RSTAZ(A(LIPODS),LPUU,59)
      LSKP=LMAX
      LMAX=LSKP+NWK*IDVA
CE    READ LINEAR MATRIX K FROM DISK 4
CS    UCITAVANJE LINEARNE MATRICE K SA DISKA 4
C!!!! PROVERI DA LI JE U (LSK) MATRICA KRUTOSTI !!!!!!!???????
C     CALL RSTAZK(A(LIPODS),LSKP,35)
C     (T)R=(T)R-K*(T)UU
cc      IF(KOR.GT.1) THEN
        NUL=NWK
        IF(NBLOCK.GT.1) NUL=KC
       CALL CLEARB(A(LSKP),A(LMAXA),A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,NUL)
        OPEN (ISCRC,FILE='ZSKLIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND ISCRC
        LLM =LMAX
        LSKE=LLM+100
        CALL SPAKUA(A(LSKP),A(LMAXA),A(LSKE),A(LLM),ND,1,
     &              A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
        CLOSE (ISCRC,STATUS='KEEP')
cc      ENDIF
      CALL MAXAPM(A(LSKP),A(LPUU),A(LRTDT),A(LMAXA),JEDN)
c      call wrr6(a(LskP),NWK,'sk  ')
c      call wrr6(a(LRTDT),jedn,'R1  ')

C
C     (T-T)UU
      CALL RSTAZ(A(LIPODS),LUMC,53)
c      call wrr6(a(LUMC),jedn,'-DTU')

CE    READ LINEAR EFECTIVE MATRIX (A0*M+A1*C) FROM DISK 7
CS    UCITAVANJE LINEARNE EFEKTIVNE MATRICE (A0*M+A1*C) SA DISKA 7
      IF(IMASS.NE.2) CALL RSTAZK(A(LIPODS),LSKP,58)
c      IF(IMASS.EQ.2) CALL RSTAZ(A(LIPODS),LSKP,58)
C sneza proba
      IF(IMASS.EQ.2) CALL RSTAZK(A(LIPODS),LSKP,58)
c      call wrr6(a(LSKP),NWK,'MEFE')
c      call iwrr(a(LMAXA),jedn,'MAXA')

C     (T)R=(T)R+(A0*M+A1*C)*(T-T)UU
      if (IABS(ICCGG).EQ.1) then 
         IF(IMASS.NE.2) CALL MAXAPRI(A(LSKP),A(LUMC),A(LRTDT),
     1                       JEDN,A(IROWS),A(IROWS+nwk),nwk)
      else
         IF(IMASS.NE.2) CALL MAXAPR(A(LSKP),A(LUMC),A(LRTDT),
     &                       A(LMAXA),JEDN)
      endif
C      IF(IMASS.EQ.2) CALL MNOZMU(A(LRTDT),A(LSKP),A(LUMC),JEDN)
C sneza proba
      IF(IMASS.EQ.2) CALL MAXAPR(A(LSKP),A(LUMC),A(LRTDT),A(LMAXA),JEDN)
c      call wrr6(a(LUMC),jedn,'LUMC  ')
c      call wrr6(a(LRTDT),jedn,'R2  ')
C
      CALL JEDNA1(A(LSKP),A(LUMC),JEDN)
      AM2=-A2
C     A2*(T)UU-A2*(T-T)UU
      CALL ZBIR2(A(LUMC),A(LPUU),A(LSKP),A2,AM2,JEDN)
CE    READ LINEAR MATRIX M FROM DISK 11
CS    UCITAVANJE LINEARNE MATRICE M SA DISKA 11
      IF(IMASS.NE.2) CALL RSTAZK(A(LIPODS),LSKP,54)
      IF(IMASS.EQ.2) CALL RSTAZ(A(LIPODS),LSKP,54)
c      call wrr6(a(LSKP),NWK,'MMAS')
c      call wrr6(a(LUMC),jedn,'A2DU')
C     (T)R=(T)R+A2*M*((T)UU-(T-T)UU)
      if (IABS(ICCGG).EQ.1) then 
            IF(IMASS.NE.2) CALL MAXAPRI(A(LSKP),A(LUMC),A(LRTDT),
     1                          JEDN,A(IROWS),A(IROWS+nwk),nwk)
      else
            IF(IMASS.NE.2) CALL MAXAPR(A(LSKP),A(LUMC),A(LRTDT),
     &                          A(LMAXA),JEDN)
      endif
      IF(IMASS.EQ.2) CALL MNOZMU(A(LRTDT),A(LSKP),A(LUMC),JEDN)
C
C     (T)R
      CALL WSTAZ(A(LIPODS),LRTDT,38)
c      call wrr6(a(LRTDT),jedn,'RTDT')
      LMAX=LUMC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VILSON
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO CALCULATE INCREMENT OF DISPLACEMENTS, NEW VELOCITIES AND
CE.      ACCELERATIONS IN (T+DT) - WILSON
CS.   PROGRAM
CS.      ZA RACUNANJE POMERANJA, NOVIH BRZINA I UBRZANJA
CS.      U TRENUTKU (T+T) - WILSON
C .
C .         (T+T)UA=A6*((T+�T)UU-(T)UU)+A7*(T)UV+A8*(T)UA
C .         (T+T)UV=(T)UV+A9*((T+T)UA+(T)UA)
C .         (T+T)UU=(T)UU+T*(T)UV+A10*(2*(T)UA+(T+T)UA)
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' VILSON'
      ONE=1.D0
      LMAX=LRAD
      LUMC=LMAX
      LPUU=LUMC+JEDN*IDVA
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
C     (T)UU, (T)UV, (T)UA
      CALL RSTAZ(A(LIPODS),LPUU,59)
      CALL RSTAZ(A(LIPODS),LPUV,55)
      CALL RSTAZ(A(LIPODS),LPUA,57)
C
CE    ACCELERATIONS IN (T+T)
CS    UBRZANJA U TRENUTKU (T+T)
C     (T+T)UA=A6*((T+�T)UU-(T)UU)+A7*(T)UV+A8*(T)UA
      AM6=-A6
      CALL ZBIR4(A(LUMC),A(LRTDT),A(LPUU),A(LPUV),A(LPUA),
     1A6,AM6,A7,A8,JEDN)
C
CE    VELOSITIES IN (T+T)
CS    BRZINE U TRENUTKU (T+T)
C     (T+T)UV=(T)UV+A9*((T+T)UA+(T)UA)
      CALL ZBIR3(A(LRTDT),A(LPUV),A(LUMC),A(LPUA),ONE,A9,A9,JEDN)
C
CE    DISPLACEMENTS IN (T+T)
CS    POMERANJA U TRENUTKU (T+T)
C     (T+T)UU=(T)UU+T*(T)UV+A10*(2*(T)UA+(T+T)UA)
      DA10=2.0*A10
      CALL ZBIRM3(A(LPUU),A(LPUV),A(LUMC),A(LPUA),DT,A10,DA10,JEDN)
C
C     (T+T)UV
      CALL JEDNA1(A(LPUV),A(LRTDT),JEDN)
C     (T+T)UA
      CALL JEDNA1(A(LPUA),A(LUMC),JEDN)
C
C     (T+T)UU, (T+T)UV, (T+T)UA
      CALL WSTAZ(A(LIPODS),LPUU,59)
      CALL WSTAZ(A(LIPODS),LPUV,55)
      CALL WSTAZ(A(LIPODS),LPUA,57)
      LMAX=LUMC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE NEWMAR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO CALCULATE VELOCITIES AND ACCELERATIONS
CE.      IN (T+DT) - NEWMARK
CS.   PROGRAM
CS.      ZA RACUNANJE BRZINA I UBRZANJA
CS.      U TRENUTKU (T+T) - NEWMARK
C .
C .         (T+T)UA=A6*((T+T)UU-(T)UU)+A7*(T)UV+A8*(T)UA
C .         (T+T)UV=(T)UV+A9*(T)UA+A10*(T+T)UA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GRUPER/ LIGRUP
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' NEWMAR'
      IF(ICONT.GT.0) JEDN=NEQ
      LMAX=LRAD
      LUMC=LMAX
      LPUU=LUMC+JEDN*IDVA
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
C     (T)UU, (T)UV, (T)UA
      CALL RSTAZ(A(LIPODS),LPUU,59)
      CALL RSTAZ(A(LIPODS),LPUV,55)
      CALL RSTAZ(A(LIPODS),LPUA,57)
C
CE    ACCELERATIONS IN (T+T)
CS    UBRZANJA U TRENUTKU (T+T)
C     (T+T)UA=A6*((T+T)UU-(T)UU)+A7*(T)UV+A8*(T)UA
      AM6=-A6
      CALL ZBIR4(A(LUMC),A(LRTDT),A(LPUU),A(LPUV),A(LPUA),
     1A6,AM6,A7,A8,JEDN)
C
CE    VELOSITIES IN (T+T)
CS    BRZINE U TRENUTKU (T+T)
C     (T+T)UV=(T)UV+A9*(T)UA+A10*(T+T)UA
      CALL ZBIRM2(A(LPUV),A(LPUA),A(LUMC),A9,A10,JEDN)
C
C     (T+T)UA
C.. KOREKCIJA UBRZANJA I BRZINA ZA KONTAKTNE PROBLEME
      IF(ICONT.GT.0)THEN
        IDUM=LPUU
        LPUU=LUMC
        CALL CONCOR(A(LIGRUP),7)
        CALL JEDNA1(A(LPUA),A(LUMC),JEDN)
        CALL RSTAZ(A(LIPODS),LUMC,55)
        CALL CONCOR(A(LIGRUP),8)
        LPUU=IDUM
      ELSE
        CALL JEDNA1(A(LPUA),A(LUMC),JEDN)
      ENDIF
C     (T+T)UU
      CALL JEDNA1(A(LPUU),A(LRTDT),JEDN)
C
C     (T+T)UU, (T+T)UV, (T+T)UA
      CALL WSTAZ(A(LIPODS),LPUU,59)
      CALL WSTAZ(A(LIPODS),LPUV,55)
      CALL WSTAZ(A(LIPODS),LPUA,57)
      LMAX=LUMC
      IF(ICONT.GT.0) JEDN=NEQ+NEQC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CENDIF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO CALCULATE VELOCITIES AND ACCELERATIONS
CE.      IN (T+DT) - CENTRAL DIFFERENS
CS.   PROGRAM
CS.      ZA RACUNANJE BRZINA I UBRZANJA
CS.      U TRENUTKU (T) - CENTRALNE RAZLIKE
C .
C .         (T)UA=A0*((T-T)UU-2.*(T)UU+(T+T)UU)
C .         (T)UV=A1*(-(T-T)UU+(T+T)UU)
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' CENDIF'
      LMAX=LRAD
      LUMC=LMAX
      LPUU=LUMC+JEDN*IDVA
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
C     (T)UU
      CALL RSTAZ(A(LIPODS),LPUU,59)
c      call wrr6(a(lpuu),jedn,'TU  ')
C     (T-T)UU
      CALL RSTAZ(A(LIPODS),LUMC,53)
c      call wrr6(a(lumc),jedn,'-TU ')
C
CE    ACCELERATIONS IN (T)
CS    UBRZANJA U TRENUTKU (T)
      AM2=-2.0*A0
C     (T)UA=A0*((T-T)UU-2.*(T)UU+(T+T)UU)
      CALL ZBIR3(A(LPUA),A(LUMC),A(LPUU),A(LRTDT),A0,AM2,A0,JEDN)
C
CE    VELOCITIES IN (T)
CS    BRZINE U TRENUTKU (T)
      AM1=-A1
C     (T)UV=A1*(-(T-T)UU+(T+T)UU)
      CALL ZBIR2(A(LPUV),A(LUMC),A(LRTDT),AM1,A1,JEDN)
C
C     (T-T)UU
      CALL WSTAZ(A(LIPODS),LPUU,53)
C     (T)UU
      CALL JEDNA1(A(LPUU),A(LRTDT),JEDN)
C
C     (T)UU, (T)UV, (T)UA
      CALL WSTAZ(A(LIPODS),LPUU,59)
      CALL WSTAZ(A(LIPODS),LPUV,55)
      CALL WSTAZ(A(LIPODS),LPUA,57)
      LMAX=LUMC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CENPOC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO CALCULATE OF NEGATIVE INCREMENT DISPLACEMENTS,
CE.      STORES ON UNIT 12
CS.   PROGRAM
CS.      ZA RACUNANJE POMERANJA U TRENUTKU (-T),
CS.      ZAPISANO NA DISK IDINA=12
C .
C .         (-T)U=1.*(0)UU-T*(0)UV+A3*(0)UA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /DUPLAP/ IDVA
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' CENPOC'
      ONE=1.D0
      LMAX=LRAD
      LUMC=LMAX
      LPUU=LUMC+JEDN*IDVA
      LPUV=LPUU+JEDN*IDVA
      LPUA=LPUV+JEDN*IDVA
      LMAX=LPUA+JEDN*IDVA
C     (0)UU, (0)UV, (0)UA
      CALL RSTAZ(A(LIPODS),LPUU,59)
      CALL RSTAZ(A(LIPODS),LPUV,55)
      CALL RSTAZ(A(LIPODS),LPUA,57)
      DM=-DT
C     (-T)U=1.*(0)UU-T*(0)UV+A3*(0)UA
      CALL ZBIR3(A(LUMC),A(LPUU),A(LPUV),A(LPUA),ONE,DM,A3,JEDN)
c      call wrr6(a(lumc),jedn,'lumc')
C
      CALL WSTAZ(A(LIPODS),LUMC,53)
      LMAX=LUMC
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CONSTD
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO FORM DYNAMICS CONSTANTS
CS.   PROGRMA
CS.      ZA FORMIRANJE DINAMICKIH KONSTANTI
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /EDPCON/ A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' CONSTD'
      IF(MDVI.NE.1) GO TO 20
      TAU=PIP*DT
      A0=6.0D0/(TAU*TAU)
      A1=3.0D0/TAU
      A2=2.0D0*A1
      A3=2.0D0
      A4=2.0D0
      A5=TAU/2.0D0
      A6=A0/PIP
      A7=-A2/PIP
      A8=1.0D0-3.0D0/PIP
      A9=DT/2.0D0
      A10=(DT*DT)/6.0D0
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      RETURN
C
   20 IF(MDVI.NE.2) GO TO 30
      A0=1.0D0/(DIP*DT*DT)
      A1=PIP/(DIP*DT)
      A2=1.0D0/(DIP*DT)
      A3=1.0D0/(2.0D0*DIP)-1.0D0
      A4=PIP/DIP-1.0D0
      A5=DT*(PIP/DIP-2.0D0)/2.0D0
      A6=A0
      A7=-A2
      A8=-A3
      A9=DT*(1.0D0-PIP)
      A10=PIP*DT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      RETURN
C
   30 IF(MDVI.NE.3) RETURN
      A0=1.0D0/(DT*DT)
      A1=1.0D0/(2.0D0*DT)
      A2=2.0D0*A0
      A3=1.0D0/A2
      A4=0.0D0
      A5=0.0D0
      A6=0.0D0
      A7=0.0D0
      A8=0.0D0
      A9=0.0D0
      A10=0.0D0
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(///' D I N A M I C K I     K O E F I C I J E N T I'//
     16X,'A0',11X,'A1',11X,'A2',11X,'A3',11X,'A4',11X,'A5'/6(1PE13.6)/
     1/6X,'A6',11X,'A7',11X,'A8',11X,'A9',11X,'A10'/5(1PE13.6))
C-----------------------------------------------------------------------
 6000 FORMAT(///' D Y N A M I C S    C O E F F I C I E N T S'//
     16X,'A0',11X,'A1',11X,'A2',11X,'A3',11X,'A4',11X,'A5'/6(1PE13.6)/
     1/6X,'A6',11X,'A7',11X,'A8',11X,'A9',11X,'A10'/5(1PE13.6))
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAMP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.     TO PRINT DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
CS.   PROGRAM
CS.     ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /REZREP/ LPUUU,LPUUV,LPUUA
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAMP'
      IF(ISUU) 30,20,30
C      IF(ISUU) 10,20,30
C   10 CALL SAVREZ(A(LPUU),A(LPUUU),A(LID))
C      GO TO 30
   20 CALL STAPO1(A(LPUU),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,0,
     +            A(LNCVP),NCVPR,A(LCORD),A(LELCV))
   30 IF(ISUV) 60,50,60
C   30 IF(ISUV) 40,50,60
C   40 CALL SAVREZ(A(LPUV),A(LPUUV),A(LID))
C      GO TO 60
   50 CALL STAPO1(A(LPUV),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,1,
     +            A(LNCVP),NCVPR,A(LCORD),A(LELCV))
   60 IF(ISUA) 90,80,90
C   60 IF(ISUA) 70,80,90
C   70 CALL SAVREZ(A(LPUA),A(LPUUA),A(LID))
C      GO TO 90
   80 CALL STAPO1(A(LPUA),A(LID),A(LCVEL),ICVEL,NP,KOR,VREME,2,
     +            A(LNCVP),NCVPR,A(LCORD),A(LELCV))
   90 RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE STAGP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
CE.      IN UNIVERZAL FILE - FOR GRAPHIC
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA
CS.      U UNIVERZALNI FILE - ZA GRAFIKU
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAGP'
      IF(ISUU) 30,20,30
   20 CALL STAGP1(A(LPUU),A(LID),A(LCVEL),ICVEL,NP,IGRAF,0,
     +            A(LNCVP),NCVPR)
      CALL STAU09(A(LPUU),A(LID),A(LCVEL),ICVEL,NP,49,1,
     +            A(LNCVP),NCVPR)
   30 IF(ISUV) 60,50,60
   50 CALL STAGP1(A(LPUV),A(LID),A(LCVEL),ICVEL,NP,IGRAF,1,
     +            A(LNCVP),NCVPR)
      CALL STAU09(A(LPUV),A(LID),A(LCVEL),ICVEL,NP,49,11,
     +            A(LNCVP),NCVPR)
   60 IF(ISUA) 90,80,90
   80 CALL STAGP1(A(LPUA),A(LID),A(LCVEL),ICVEL,NP,IGRAF,2,
     +            A(LNCVP),NCVPR)
      CALL STAU09(A(LPUA),A(LID),A(LCVEL),ICVEL,NP,49,21,
     +            A(LNCVP),NCVPR)
   90 RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KSTAMP(A,VDT,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT RESULTS PER STEPS
CS.   PROGRAM
CS.      ZA STAMPANJE REZULTATA PO KORACIMA - TABELARNO
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(NDT,*),VDT(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' KSTAMP'
      IF(II.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      ENDIF
      IF(II.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6001)
      ENDIF
      IF(II.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2002)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6002)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2003)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6003)
      VREME=0.0
      NPI=NP*1
      DO 50 I=1,NDT
      DT= VDT(I)
      VREME=VREME+DT
      WRITE(IZLAZ,5000) VREME,(A(I,J),J=1,NPI)
   50 CONTINUE
      RETURN
C
 5000 FORMAT(11D12.4)
C-----------------------------------------------------------------------
 2000 FORMAT(///' K O M P O N E N T E   P O M E R A NJ A'/)
 2001 FORMAT(///' K O M P O N E N T E   B R Z I N E'/)
 2002 FORMAT(///' K O M P O N E N T E   U B R Z A NJ A'/)
 2003 FORMAT('    VREME        N1          N2          N3          N4
     1       N5          N6          N7          N8          N9
     1N10')
C-----------------------------------------------------------------------
 6000 FORMAT(///' D I S P L A C E M E N T S'/)
 6001 FORMAT(///' V E L O C I T I E S'/)
 6002 FORMAT(///' A C C E L E R A T I O N S'/)
 6003 FORMAT('    TIME         N1          N2          N3          N4
     1       N5          N6          N7          N8          N9
     1N10')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SAVREZ(A,B,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO SAVE RESULTS PER STEPS
CS.   PROGRAM
CS.      ZA CUVANJE REZULTATA PO KORACIMA
C .
C ......................................................................
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(NDT,*),ID(NP,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SAVREZ'
      K=0
      DO 10 I=1,NP
      DO 20 J=1,1
      K=K+1
      NJ=ID(I,J)
      IF(NJ.EQ.0) GO TO 15
      B(KOR,K)=A(NJ)
      GO TO 20
   15 B(KOR,K)=0.0
   20 CONTINUE
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SILAAN(RTDT,T,T1,F,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      ZA FORMIRANJE HARMONIJSKIH SPOLJASNJIH SILA
CE.      CISTO ZATEZANJE
CS.   PROGRAM
CS.      ZA FORMIRANJE HARMONIJSKIH SPOLJASNJIH SILA
CS.      CISTO ZATEZANJE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION RTDT(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SILAAN'
      IF(T1.GT.0.5D0.AND.T1.LT.1.5D0) THEN
      RTDT(1)=DSIN(T)
      ENDIF
      IF(T1.GT.1.5D0.AND.T1.LT.2.5D0) THEN
      RTDT(1)=DSIN(T)/2.0D0
      RTDT(2)=DSIN(T)/2.0D0
      ENDIF
      IF(T1.GT.2.5D0.AND.T1.LT.3.5D0) THEN
      RTDT(1)=DSIN(T)/4.0D0
      RTDT(2)=DSIN(T)/4.0D0
      RTDT(3)=DSIN(T)/4.0D0
      RTDT(4)=DSIN(T)/4.0D0
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE SILAA(RTDT,T,T1,F,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      ZA FORMIRANJE HARMONIJSKIH SPOLJASNJIH SILA
CE.      POBUDA NA TOCKOVIMA
CS.   PROGRAM
CS.      ZA FORMIRANJE HARMONIJSKIH SPOLJASNJIH SILA
CS.      POBUDA NA TOCKOVIMA
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION RTDT(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SILAAN'
      W=F*2.*3.14159265358979
      X1=A*DSIN(W*T)
      X2=0.0
      TT1=T-T1
      IF(TT1.LT.1.0D-10) GO TO 10
      X2=A*DSIN(W*TT1)
   10 RTDT(1)= 32000.0*(X1+X2)
      RTDT(2)= 32000.0*(0.88*X1-1.556*X2)
      RTDT(3)= -32000.0*0.311*X1
      RTDT(4)= -32000.0*X2
      RETURN
      END
C=======================================================================
      SUBROUTINE CONCOR(IGRUP,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   LOOP FOR GROUP OF ELEMENTS TO INTEGRATION MATRIX M
CS.   PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE  M
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IGRUP(NGELEM,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CONCOR'
      DO 100 NGE = 1,NGELEM
       NETIP = IGRUP(NGE,1)
C.. SAMO ZA KONTAKTNE ELEMENTE ( 10 = 2D , 11 = 3D )
       IF(NETIP.EQ.10.OR.NETIP.EQ.11)THEN
         NE    = IGRUP(NGE,2)
         IATYP = IGRUP(NGE,3)
         NMODM = IGRUP(NGE,4)
         LMAX8 = IGRUP(NGE,5)
         IF(IND.EQ.9) LMAX  = LRAD
         CALL ELEME(NETIP,IND)
       ENDIF
  100 CONTINUE
      RETURN
      END
