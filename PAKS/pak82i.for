C=======================================================================
C
C=======================================================================
      SUBROUTINE MATEIL(TAU,SKE,FTDT,T,WTU,GKS,LM,MSL,NOP,DRGTM0,DEB,
     1                  TGT,LPLAS,LPLA1,ISNA,MTR,TTR,IND6,R,S,COR,GMT,
     1                  CORT,GNT,GERS,ALFE,HAEM,HINV,GEEK,
     1                  N45,LA,CEGE,ESILA,IPG,DEF,INTGL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO FORM ELEMENT MATRIX FOR ISOPARAMETRIC SHELL ELEMENT 
CS.   P R O G R A M
CS.      ZA FORMIRANJE MATRICA ELEMENATA IZOPARAMETARSKE LJUSKE
C .
C .      IPODT.EQ.0; (4-9 CVOROVA)
C .      IPODT.EQ.1; (3-6 CVOROVA)
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /PLASTV/ LPLAVT,LPLAV1,LSIGMV
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /SIMO2D/ LALFE,LHAEM,LHINV,LGEEK,IALFA
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /DUPLAP/ IDVA
      COMMON /ITERBR/ ITER
      COMMON /INDNAP/ NAPON
      COMMON /FILTER/ TOLNAP,TOLPOM
      COMMON /IKOVAR/ INDKOV
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /CUVAJN/ TAUCEN(6,5),INDSEL,NGAUSU
      COMMON /MATIZO/ E,V
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      COMMON /DEBLJG/ THICK,THICT,THIC0,NNSL
      COMMON /TRANDN/ TSGD(6,6),TSGN(6,6)
      COMMON /DEFNAP/ NAPDEF
      COMMON /DRILLI/ EL14(4),EN14(3,4),IDRIL
C
      DIMENSION NOP(NE,*),SKE(*),LM(*),FTDT(*),TAU(N45,NGS12,NE,*),
     1          TA(6),STRAIN(6),STRESS(6),DEF(N45,NGS12,NE,*),
     1          ALFE(LA,*),HAEM(LA,*),HINV(LA,LA,*),GEEK(LA,ND,*)
      DIMENSION GKS(3,2,*),ESILA(ND,*)
      DIMENSION DRGTM0(9,*),MTR(*),TTR(3,3,*),DEB(NE,*),BNL(9,54)
      DIMENSION TSS(6,6),
     1          TRLN(6,6),COR(9,*),CON(9,3),GMT(3,*),CORT(9,*),GNT(3,*),
     1          GRAN(3,3),XJT(3,3),XJ0(3,3),CT(3,3,3,3)
      DIMENSION GERS(6,*),CEGE(LA,*),SKOV(6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MATEIL'
C     INDIKATOR KONTROLNE STAMPE
      IST=0
      IF(IST.EQ.1) WRITE(3,*)'NLM,NGR,NGS,NGT',NLM,NGR,NGS,NGT
C
CS    JAKOBIJAN - XJJ, INVERTOVAN JAKOBIJAN - XJ
CE    JACOBIAN - XJJ, INVERSE JACOBIAN - XJ
C
      CALL JEDNA1(XJJ,XJ,9)  
      CALL MINV3(XJ,DET)
C
CS    MATRICA TRANSFORMACIJE - TRLN (KOVARIJANTNI - GLOBALNI DEKARTOV) 
CE    TRANSFORMATION MATRIX - TRLN (COVARIANT - GLOBAL CARTESIAN)
C
      IF(INDKOV.LT.0) CALL TRANSE(TRLN,XJ)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR(TRLN,36,'TRLN') 
C
CS    FORMIRANJE MATRICE BL LINEARNO  (BLT)
CE    FORM LINEAR BL MATRIX (BLT) 
C
      IF(INDSEL.EQ.1) CALL CLEAR(BLT,6*54)
C
      CALL BETBE8(H,BLT,BNL,NOP,GKS,T,MTR)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR(BLT,6*ND,'BLT ') 
C
CS    TRANSFORMACIJA MATRICE - BLT (KOVARIJANTNI - GLOBALNI DEKARTOV) 
CE    TRANSFORM MATRIX - BLT (COVARIANT - GLOBAL CARTESIAN)
C
      IF(INDKOV.LT.0) CALL TRANS4(BLT,TRLN,ND)
C
      IF(IDRIL.EQ.1) CALL BDRIL8(H,BLT)
C
CS    FORMIRANJE MATRICE BL1 GEOMETRIJSKA NELINEARNOST - T.L.
CE    FORM MATRIX BL1 GEOMETRICAL NONLINEARITY -  T.L.
C
      IF(IATYP.EQ.2.AND.INDKOV.GE.0)
     1CALL BETL18(BLT,NOP,LM,UEL,STRAIN,H,DRGTM0,GKS,T,THICK,0)
C
CS    TRANSFORMACIJA MATRICE - BLT ZA LOKALNI SISTEM CVORA
CE    TRANSFORM MATRIX - BLT FOR LOCAL NODAL SISTEM
C
      IF(IND6.GT.0) THEN
         ICD=6
         CALL TRAB(BLT,ICD,MTR,TTR,NCVE)
      ENDIF
C
CS    MATRICA TRANSFORM. DEFORMACIJA - TSG (GLOBALNI - LOKALNI DEKARTOV)
CE    STRAIN TRANSFORMATION MATRIX - TSG (GLOBAL - LOCAL CARTESIAN) 
C
C          IF(IST.EQ.1) CALL WRR3(XJJ,9,'XJJ ')
      CALL JEDNA1(XJ0,XJJ,9)
      CALL TRANAL(XJ0,TSG,0)
      CALL TRANSS(TSGN,XJ0)
      CALL JEDNA1(TSGD,TSG,36)
C  MATRICA KOSINUSA BETA IZ MAGISTARSKOG ZA TRASFORMACIJU TENZ. KOMP.
C      CALL JEDNA1(TBET,XJ0,9)
C      IF(IST.EQ.1) CALL WRR(TBET,9,'TBET')
C
      IF(ISKNP.NE.1) THEN
C
C
CS       RACUNANJE DEFORMACIJA U GLOBALNOM DEKARTOVOM
CE       CALCULATE STRAINS IN GLOBAL CARTESIAN 
C
C
         CALL CLEAR(STRAIN,6) 
C
CS       LINEARNOST I M.N.O.
CE       LINEAR PART AND M.N.O.
C
         IF(IATYP.EQ.0.OR.IATYP.EQ.1)
     1   CALL MNOZI1(STRAIN,BLT,UEL,6,ND)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(uel,nd,'uel0') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(blt,6*nd,'blt0') 
C
CS       GEOMETRIJSKA NELINEARNOST - UKUPNE DEFORMACIJE ZA T.L. I U.L.
CE       GEOMETRICAL NONLINEARITY - TOTAL STRAIN FOR T.L. AND U.L.
C
          IF(INDKOV.GE.0) THEN
            IF(IATYP.EQ.2)
     1      CALL BETL18(BLT,NOP,LM,UEL,STRAIN,H,DRGTM0,GKS,T,THICK,1)
            IF(IATYP.EQ.3)
     1      CALL STRUL8(H,NOP,LM,DRGTM0,UEL,STRAIN,T,THICK)
C      IF(NLM.EQ.1.AND.JG.EQ.1) write(3,*) 't,thick',t,thick 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL iWRR(lm,nd,'  lm') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(h,27,'   h') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(drgtm0,27,'drgt') 
          ELSE
            IF(IATYP.EQ.2.OR.IATYP.EQ.3)
     1      CALL UKDEFI(COR,CON,GM,GMT,UEL,H,R,S,T,STRAIN,TRLN,TSS,
     1                  IATYP,NCVE)
          ENDIF
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(strain,6,'str0') 
C
CS       VELIKE DEFORMACIJE
CE       LARGE STRAIN
C
         IF(IATYP.GE.4) THEN
C
CS          STARO - Be, Cp**-1, Fp**-1
CE          OLD - B
C
            IF(ILEDE.EQ.0.AND.ICPM1.LE.1) THEN
               CALL JEDNA1(STRESS,TAU(1,JG,NLM,MSL),NLD)
               IF(IST.EQ.1) CALL WRR(STRESS,NLD,'BOLD')
            ENDIF
            IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN 
               CALL JEDNA1(GRAP,TAU(1,JG,NLM,MSL),NLD)
               CALL MINV3(GRAP,DUM)
               IF(IST.EQ.1) CALL WRR(GRAP,NLD,'BOLD')
            ENDIF
C
CS          JAKOBIJEVA MATRICA U TACKI (R,S,T) U TRENUTKU - T,0
CE          JACOBIAN MATRIX IN POINT (R,S,T) IN TIME - T,0
C
CS          KOORDINATE U TRENUTKU - 0
CE          COORDINATE IN TIME - 0
            II=0
            DO 7 I=1,NCVE
            DO 7 K=1,3
               II=II+1
               CON(I,K)=COR(I,K)-UEL(II)
    7       CONTINUE
CS          KOVARIJANTNI BAZNI VEKTORI U KONFIGURACIJI - XJT,XJ0
CE          COVARIANT BASIS VECTORS IN CONFIGURATION - XJT,XJ0
            CALL GRGSGT(CORT,GNT,XJT,H,T,NCVE)
            CALL GRGSGT(CON,GMT,XJ0,H,T,NCVE)
C
CS          INVERTOVAN JAKOBIJAN - XJT,XJ0
CE          INVERSE JACOBIAN - XJT,XJ0
C
            CALL MINV3(XJT,DUM)
            CALL MINV3(XJ0,DUM)
C
CS          RACUNANJE GRADIJENTA DEFORMACIJE OD T DO T+DT
CE          CALCULATE DEFORMATION GRADIENT FROM T TO T+DT
C
            IF(ILEDE.EQ.0) THEN
               IF(ICPM1.EQ.0) CALL MNOZM4(GRAD,XJJ,XJT,3,3,3)
               IF(ICPM1.GE.1) CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
            ENDIF
            IF(ILEDE.EQ.1) CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
C            IF(IST.EQ.1) CALL WRR3(XJJ,9,'XJJ ')
C            IF(IST.EQ.1) CALL WRR3(XJT,9,'XJT ')
C            IF(IST.EQ.1) CALL WRR3(GRAD,9,'GRAD')
CALFA            IF(IALFA.GE.0) CALL ZBIRM1(GRAD,GRAN,9)
CS          RACUNANJE NORMIRANOG GRADIJENTA DEFORMACIJE OD T DO T+DT
            CALL JEDNA1(GRAN,GRAD,9)
C            CALL DETER3(GRAD,DETG)
C            IF(DETG.GT.DETGM) DETGM=DETG
C            DETG=DEXP(-1.D0/3.D0*DLOG(DETG))
C            CALL JEDNAK(GRAN,GRAD,DETG,9)
C            IF(IST.EQ.1) CALL WRR3(GRAN,9,'GRAN')
CS          RACUNANJE GRADIJENTA DEFORMACIJE OD 0 DO T+DT
CE          CALCULATE DEFORMATION GRADIENT FROM 0 TO T+DT
            CALL MNOZM4(GRAD,XJJ,XJ0,3,3,3)
            CALL DETER3(GRAD,DETG)
            IF(DETG.LT.1.D-15) STOP 'DETG.LE.0, PAK82I.FOR'
C            IF(IST.EQ.1) CALL WRR3(XJJ,9,'XJJ ')
C            IF(IST.EQ.1) CALL WRR3(XJ0,9,'XJ0 ')
C            IF(IST.EQ.1) CALL WRR3(GRAD,9,'GRAD')
C            IF(IST.EQ.1) WRITE(3,*) 'DETG',DETG
C
CS          TRANSF. ELASTICNOG LEVOG KOSI - GRINOVOG TENZORA DEFORMACIJE
CE          TRANSFORM. ELASTIC LEFT CAUCHY - GREEN DEFORMATION TENSOR
C           Be* = Ft * Be * FtT  (ICPM1=0)
C           Be* = F * (Cp**-1) * FT  (ICPM1=1)
            IF(ILEDE.EQ.0.AND.ICPM1.LE.1) CALL PIOKOS(GRAN,STRESS)
CS          RACUNANJE ELASTICNOG GRADIJENTA DEFORMACIJE
CE          ELASTIC DEFORMATION GRADIENT
C           Fe = F * Fp**-1         Be = Fe * FeT
            IF(ILEDE.EQ.0.AND.ICPM1.EQ.2) THEN
C  GRAD ILI GRAN
               CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
               CALL MNOZM3(GRAP,GRAE,GRAE,3,3,3)
               CALL TENVEK(GRAP,STRESS,1)
               CALL JEDNA1(XJ,GRAE,9)
            ENDIF
C           Fe = F * Fp**-1         Ce = FeT * Fe
            IF(ILEDE.EQ.1) THEN
C  GRAD ILI GRAN
               CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
               CALL MNOZM2(GRAP,GRAE,GRAE,3,3,3)
               CALL TENVEK(GRAP,STRESS,1)
               CALL JEDNA1(XJ,GRAE,9)
            ENDIF
C
C           GLAVNE VREDNOSTI I PRAVCI
C
C            IF(IST.EQ.1) CALL WRR6(STRESS,6,'BECE')
            IF(IGLPR.EQ.1) THEN
               NAPDEF=1
               CALL GLAVN(STRESS)
               NAPDEF=0
               CALL GLAPR3(STRESS,QP)
               CALL GLASOR(PRINC,QP,XJJ(3,1),XJJ(3,2),XJJ(3,3),IMAX)
C               IF(IST.EQ.1) CALL WRR3(PRINC,3,'PP  ') 
C               IF(IST.EQ.1) CALL WRR3(QP,9,'QP  ') 
            ENDIF
C
CS          PROBNE ELASTICNE DEFORMACIJE
CE          TRIAL ELASTIC STRAIN
C
            IF(IATYP.EQ.4) THEN
               IF(IGLPR.EQ.1) THEN
                  CALL CLEAR(STRESS,6)
                  STRESS(1)=PRINC(1)
                  STRESS(2)=PRINC(2)
                  STRESS(3)=PRINC(3)
               ENDIF
               CALL UKDEFV(STRAIN,STRESS)
            ENDIF
            IF(IATYP.EQ.5) THEN
C              GLAVNE VREDNOSTI
C              LAMBDA
               P1=DSQRT(PRINC(1))
               P2=DSQRT(PRINC(2))
               P3=DSQRT(PRINC(3))
               IF(IST.EQ.1) WRITE (3,*) 'L1,L2,L3',P1,P2,P3
               STRAIN(1)=DLOG(P1)
               STRAIN(2)=DLOG(P2)
               STRAIN(3)=DLOG(P3)
C               IF(IST.EQ.1) CALL WRR6(STRAIN,6,'STLO')
            ENDIF
            IF(IGLPR.EQ.1) THEN
C              TRANSFORMACIJA DEFORMACIJA: GLAVNI PRAVCI - DEKARTOV SISTEM
               CALL DIJADS(QP,STRAIN)
               STRAIN(4)=2.D0*STRAIN(4)
               STRAIN(5)=2.D0*STRAIN(5)
               STRAIN(6)=2.D0*STRAIN(6)
C               IF(IST.EQ.1) CALL WRR6(STRAIN,6,'STDE')
            ENDIF
         ENDIF
C
CS       INKOMPATIBILNE DEFORMACIJE
CE       INCOMPATIBILE STRAIN
C        ovaj koncept ne valja za velike deformacije
         IF(IALFA.GE.0.AND.IATYP.LT.4)
     +   CALL MNOZI1(STRAIN,GERS,ALFE(1,NLM),6,LA)
C
CS       TRANSFORMACIJA DEFORMACIJA (GLOBALNI - LOKALNI DEKARTOV)
CE       TRANSFORM STRAIN (GLOBAL - LOCAL CARTESIAN) 
C
C        El=Qd*Ed
         CALL CLEAR(TA,6) 
         CALL MNOZI1(TA,TSG,STRAIN,6,6)
         TA(3)=0.D0
         CALL JEDNA1(STRAIN,TA,6)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(tsg,36,' tsg') 
C
CS       TERMICKE DEFORMACIJE   ETH=ALFA*(T-T0)
CE       THERMAL STRAINS
C
         IF(NMODM.EQ.3.OR.NMODM.EQ.4.AND.NGAUSU.EQ.NGS12) THEN
            DTGT0=TGT-TEMP0
            STRAIN(1)=STRAIN(1)-ALFA(1)*DTGT0
            STRAIN(2)=STRAIN(2)-ALFA(2)*DTGT0
            IF(NAPON.EQ.1) STRAIN(3)=ALFA(3)*DTGT0
            STRAIN(4)=STRAIN(4)-ALFA(4)*DTGT0
         ENDIF
C         IF(IST.EQ.1) CALL WRR6(STRAIN,6,'STRA')
C
CS       RACUNANJE NAPONA
CE       CALCULATE STRESS 
C
         CALL CLEAR(STRESS,6)
C
         IF(NMODM.LT.5) THEN
C
CS          MATERIJALNA LINEARNOST
CE          MATERIAL LINEAR
C
CS          KOSIJEVI NAPONI U GLOBAL. ILI LOKAL. DEKART. SISTEMU
CE          CAUCHY STRESS IN GLOBAL OR LOCAL CARTESIAN SYSTEM 
C
            CALL MNOZI1(STRESS,ETP,STRAIN,6,6)
C           CALL WRR6(STRESS,6,'STRE')
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(strain,6,'stra') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(etp,36,' etp') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(stress,6,'str1') 
C
CS          TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM
CE          TRANSFORM STRESS TO GLOBAL COORDS.
C
            CALL CLEAR(TA,6) 
            CALL MNOZI2(TA,TSG,STRESS,6,6)
            IF(ISNA.EQ.2) CALL JEDNA1(STRESS,TA,6)
C
CS          RACUNANJE UNUTRASNJIH SILA
CE          CALCULATE INTERNAL FORCES
C           r = BT * S 
C
            IF(NZADP.GT.0.AND.ISKNP.EQ.2.AND.NGENL.EQ.0)
     1      CALL INTEGF(FTDT,BLT,TA,LM,WTU,ND,6)
C
CS          CISCENJE NUMERICKIH GRESAKA ZA NAPONE
CE          CLEANING NUMERICAL ERRORS FOR STRESS
C
C           IF(IATYP.EQ.0) CALL CISTIN(STRESS,6)
C
            IF(INDSEL.EQ.1) THEN
               IF(NGAUSU.LT.NGS12) THEN
                  CALL JEDNA1(TAUCEN(1,NGT),STRESS,6)
               ELSE
                  CALL ZBIR2B(TAU(1,JG,NLM,MSL),TAUCEN(1,NGT),STRESS,6)
               ENDIF
            ELSE
               CALL JEDNA1(TAU(1,JG,NLM,MSL),STRESS,6)
            ENDIF
C            IF(NGAUSU.LT.NGS12) CALL WRR(TAUCEN(1,NGT),6,'TAUC')
C            IF(NGAUSU.LT.NGS12) CALL WRR(TAUCEN(1,1),6,'TAUC')
C
            IF(NAPON.EQ.1.AND.NGAUSU.EQ.NGS12) THEN
               IF(NMODM.EQ.1.OR.NMODM.EQ.3) THEN
                  STRAIN(3)=STRAIN(3)-V*(STRAIN(1)+STRAIN(2))/(1.-V)
               ENDIF
               IF(NMODM.EQ.2.OR.NMODM.EQ.4) THEN
                  DUM=(1.-VXY*VXY*EY/EX)
                  DUM1=(VZX*EX/EZ+VXY*VYZ)/DUM
                  DUM2=(VYZ+VXY*VZX*EY/EZ)/DUM
                  STRAIN(3)=STRAIN(3)-DUM1*STRAIN(1)-DUM2*STRAIN(2)
               ENDIF
            ENDIF
C
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(STRESS,6,'STRE')
         ELSE
C
CS          NAPONI ZA PLASTICAN MODEL
CE          STRESS FOR MATERIAL NONLINEARITY 
C
            IRAC=1
            CALL MODMA8(STRAIN,STRESS,NMODM,IRAC,LPLAS,LPLA1,TGT)
C            IF(IST.EQ.1) CALL WRR6(STRAIN,6,'DEFL')
C            IF(IST.EQ.1) CALL WRR6(STRESS,6,'STRL')
C
CS          TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM
CE          TRANSFORM STRESS TO GLOBAL COORDS.
C
            CALL CLEAR(TA,6) 
            CALL MNOZI2(TA,TSG,STRESS,6,6)
C
CS          PROMENA ELASTICNOG LEVOG KOSI - GRINOVOG TENZORA DEFORMACIJE
CE          UPDATE ELASTIC LEFT CAUCHY - GREEN DEFORMATION TENSOR
C
C  OVO PROVERITI ????????????????????????????
            IF(NAPON.EQ.1.AND.IATYP.GE.4) THEN
               IF(ILEDE.EQ.0) THEN
                  IF(ICPM1.EQ.0) THEN
                     CALL JEDNA1(TAU(1,JG,NLM,MSL),STRAIN,NLD)
                  ENDIF
                  IF(ICPM1.EQ.1) THEN
CS                   TRANSF. ELAS. LEVOG KOSI - GRINOVOG TENZORA DEFOR.
CS                   U DESNI PLASTICNI INVERZNI
CE                   TRANSFORM. ELASTIC LEFT CAUCHY - GREEN DEFOR. TENS.
C                    CP**-1 = F**-1 * Be * F**-T
                     CALL MINV3(GRAN,DUM)
C                     CALL WRR(GRAN,9,'GR-1')
                     CALL PIOKOS(GRAN,STRAIN)
C                     CALL WRR6(STRAIN,6,'CP-1')
C                     CALL WRR6(TAU(1,JG,NLM,MSL),6,'BST ')
                     CALL JEDNA1(TAU(1,JG,NLM,MSL),STRAIN,6)
C                     CALL WRR6(TAU(1,JG,NLM,MSL),6,'BNOV')
C                    PROVERA PLASTICNE ROTACIJE
C                     IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                        WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C                        CALL WRR3(PRINC,3,'PRI0')
C                        CALL WRR3(QP,9,'QP0 ')
C                        CALL GLAVN(STRAIN)
C                        CALL GLAPR3(STRAIN,QP)
C                        CALL WRR3(PRINC,3,'PRIN')
C                        CALL WRR3(QP,9,'QP  ')
C                     ENDIF
                  ENDIF
               ENDIF
               IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN 
                  CALL JEDNA1(GRAD,TAU(1,JG,NLM,MSL),NLD)
                  P1=DEXP(STRAIN(1))
                  P2=DEXP(STRAIN(2))
                  P3=DEXP(STRAIN(3))
                  CALL DIJAD(GRAN,QP,QP,P1,P2,P3)
                  CALL MNOZM1(GRAP,GRAN,GRAD,3,3,3)
                  CALL JEDNA1(TAU(1,JG,NLM,MSL),GRAP,NLD)
C                 PROVERA PLASTICNE ROTACIJE
C                  IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                     WRITE(3,*) 'NLM,JG,ITER,KOR',NLM,JG,ITER,KOR
C                     CALL WRR3(GRAD,9,'GRAD')
C                     CALL WRR3(PRINC,3,'PRI0')
C                     CALL WRR3(QP,9,'QP0 ')
C                     CALL WRR3(GRAP,9,'GRAP')
C                     CALL MNOZM2(GRAN,GRAP,GRAP,3,3,3)
C                     CALL WRR3(GRAN,9,'GRAN')
C                     CALL TENVEK(GRAN,STRESS,1)
C                     CALL GLAVN(STRESS)
C                     CALL GLAPR3(STRESS,QP)
C                     CALL WRR3(PRINC,3,'PRIN')
C                     CALL WRR3(QP,9,'QP  ')
C                     CALL MNOZM3(GRAN,GRAP,GRAP,3,3,3)
C                     CALL WRR3(GRAN,9,'GRA1')
C                     CALL TENVEK(GRAN,STRESS,1)
C                     CALL GLAVN(STRESS)
C                     CALL GLAPR3(STRESS,QP)
C                     CALL WRR3(PRINC,3,'PRI1')
C                     CALL WRR3(QP,9,'QP1 ')
C                  ENDIF
                  NAPKO=0
                  IF(NAPKO.EQ.1) THEN
                     IF(IATYP.EQ.4.AND.ILEDE.EQ.1) THEN
                        CALL MINV3(GRAP,DETP)
C                        IF(DETP.GT.DETGP) DETGP=DETP
                        CALL MNOZM1(GRAE,GRAD,GRAP,3,3,3)
CS                      TRANSF. PIOLA KIRKOFOV - KOSIJEV NAPON 
CE                      TRANSFORM. PIOLA KIRCKOF - CAUCHY STRESS
C                       s = F * S * FT
                        CALL PIOKOS(GRAE,TA)
CS                      TRANSFORM. ELAST. - ETP (LOKALNI - GLOBALNI DEKART.)
CE                      TRANSFORM ELAST. - ETP (LOCAL - GLOBAL CARTESIAN)
                        CALL TRAETP(ETP,ELAST,TSGD)
                        CALL CEPMT(ELAST,CT,0)
C                       Cmnop = Fmi Fnj Fok Fpl Cijkl
                        CALL RRRRC(ELAST,CT,GRAE,1)
                     ENDIF
                  ENDIF
               ENDIF
C               IF(IST.EQ.1) CALL WRR6(TAU(1,JG,NLM,MSL),NLD,'BNOV')
            ENDIF
         ENDIF
C
C        DEFORMACIJE ZBOG KOREKCIJE DEBLJINE
C          CALL WRR(ETP,36,'ELAS')
         IF(NAPON.EQ.1.AND.NMODM.LT.5) THEN
            CALL JEDNA1(DEF(1,JG,NLM,MSL),STRAIN,6)
C            CALL WRR(STRAIN,6,'STRA')
         ENDIF
C
C
CS       RACUNANJE UNUTRASNJIH SILA
CE       CALCULATE INTERNAL FORCES
C
C
         IF(NGENL.GT.0) THEN
C
CS          INTEGRACIJA VEKTORA - h
CE          INTEGRATE VEKTOR - h
C           h = GT * S
C
            IF(IALFA.GE.0) 
     +      CALL INTEV2(HAEM(1,NLM),GERS,TA,WTU,LA,6)
C
CS          RACUNANJE UNUTRASNJIH SILA
CE          CALCULATE INTERNAL FORCES
C           r = BT * S 
C
            CALL INTEGF(FTDT,BLT,TA,LM,WTU,ND,6)
C            IF(IST.EQ.1) CALL WRR3(TA,6,'STRF')
C            CALL WRR(FTDT,JEDN,'F82I')
         ENDIF
C
         IF(IPG.GT.0) THEN
            CALL JEDNAK(STRAIN,TA,WTU,6)
            CALL MNOZI2(ESILA(1,IPG),BLT,STRAIN,ND,6)
         ENDIF
C
      ENDIF
C
C
CS    INTEGRACIJA MATRICE KRUTOSTI ELEMENTA - SKE
CE    INTEGRATE ELEMENT STIFFNESS MATRIX - SKE
C
C
      IF(ISKNP.EQ.2) RETURN
C
CS    TRANSFORM. MATRICE ELAST. - ETP (LOKALNI - GLOBALNI DEKART.)
CE    TRANSFORM ELASTICITY MATRIX - ETP (LOCAL - GLOBAL CARTESIAN)
C
      CALL TRAETP(ETP,ELAST,TSG)
C
CS    INTEGRACIJA MATRICE KRUTOSTI ELEMENTA - SKE
CE    INTEGRATE ELEMENT STIFFNESS MATRIX - SKE
C     SKE = BT * C * B
C
c                  IF(IST.EQ.1) CALL WRR6(TSG,36,'TSG ')
c                  IF(IST.EQ.1) CALL WRR6(ELAST,36,'ELAS')
      CALL INTEGK(SKE,BLT,ELAST,LM,WTU,ND,6)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(elast,36,'elas') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(blt,6*nd,' blt') 
C
CS    GEOMETRIJSKI NELINEARAN DEO MATRICE KRUTOSTI
CE    GEOMETRIC NONLINEAR PART OF STIFFNESS MATRIX   
C
      IF(IATYP.GT.1) THEN
      IF(INDKOV.EQ.0.OR.INDKOV.EQ.1)
     1CALL KNL8(SKE,H,NOP,LM,TA,T,WTU,BNL,GKS,ND,MTR,TTR,IND6)
C      IF(indkov.ge.0.and.NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(ta,6,'  ta') 
C      IF(indkov.ge.0.and.NLM.EQ.1.AND.JG.EQ.1)CALL WRR6(bnl,9*nd,'bnl1') 
      IF(INDKOV.EQ.-1) THEN
         CALL CLEAR(BE,81)
         CALL CLEAR(SKOV,6) 
C
CS       TRANSFORMACIJA NAPONA SA GLOBALNOG NA KONTRAVARIJANTNI SISTEM
CE       TRANSFORM STRESS FROM GLOBAL TO CONTRAVARIANT COORD. SYSTEM
C
         CALL MNOZI2(SKOV,TRLN,TA,6,6)
C         CALL WRR(SKOV,6,'SKOV')
         SKOV(3)=0.D0
         II=-3
         DO 39 I=1,3
            II=II+3
         DO 39 J=1,3
   39    BE(II+J,II+J)=SKOV(I)
         BE(1,4)=SKOV(4)
         BE(2,5)=SKOV(4)
         BE(3,6)=SKOV(4)
         BE(1,7)=SKOV(6)
         BE(2,8)=SKOV(6)
         BE(3,9)=SKOV(6)
         BE(4,7)=SKOV(5)
         BE(5,8)=SKOV(5)
         BE(6,9)=SKOV(5)
         DO 41 I=1,9
         DO 41 J=I,9
   41    BE(J,I)=BE(I,J)
C
         CALL INTEGK(SKE,BNL,BE,LM,WTU,ND,9)
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(be,81,' be') 
C      IF(NLM.EQ.1.AND.JG.EQ.1) CALL WRR6(bnl,9*nd,' bnl') 
      ENDIF
      ENDIF
C
      IF(IALFA.GE.0) THEN
C
CS       PROIZVOD MATRICA - GT * C
CE       MULTIPLY MATRIX - GT * C
C
         CALL MNOZM2(CEGE,GERS,ELAST,LA,6,6)
C
CS       INTEGRACIJA MATRICE - H
CE       INTEGRATE MATRIX - H
C        H = GT * C * G
C
         CALL INTEM1(HINV(1,1,NLM),CEGE,GERS,WTU,LA,LA,6)
C
CS       INTEGRACIJA MATRICE - GAMA
CE       INTEGRATE MATRIX - GAMA
C        GAMA = GT * C * B
C
         CALL INTEM1(GEEK(1,1,NLM),CEGE,BLT,WTU,LA,ND,6)
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE BETBE8(H,BET,BNL,NEL,GKS,T,MTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    MATRICE BET (BE TRANSPONOVANO)
CE    MATRIX BET (BE TRANSPOSED)
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /IKOVAR/ INDKOV
      COMMON /CUVAJN/ TAUCEN(6,5),INDSEL,NGAUSU
      COMMON /CDEBUG/ IDEBUG
      DIMENSION H(9,*),BET(6,*),NEL(NE,*),GKS(3,2,*),MTR(*),RJ(2,3),
     1          BNL(9,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' BETBE8'
      JJ=-2
      JF=NCVE*3-2
C
      DO 40 I=1,NCVE
         JJ=JJ+3
         JF=JF+3
         IF(NEL(NLM,I).EQ.0) GO TO 40
         JJ1=JJ+1
         JJ2=JJ+2
         JF1=JF+1
         JF2=JF+2
C
         IF(INDKOV.GE.0) THEN
C
CS          DEO VEZAN ZA TRANSLACIJE
CE          TRANSLATION PART
C
C           D/DX
            BET(1,JJ )=XJ(1,1)*H(I,2)+XJ(1,2)*H(I,3)
            BET(1,JJ1)=0.D0
            BET(1,JJ2)=0.D0
C           D/DY
            BET(2,JJ )=0.D0
            BET(2,JJ1)=XJ(2,1)*H(I,2)+XJ(2,2)*H(I,3)
            BET(2,JJ2)=0.D0
C           D/DZ
            BET(3,JJ )=0.D0
            BET(3,JJ1)=0.D0
            BET(3,JJ2)=XJ(3,1)*H(I,2)+XJ(3,2)*H(I,3)
C           D/DY , D/DX
            BET(4,JJ )=BET(2,JJ1)
            BET(4,JJ1)=BET(1,JJ)
            BET(4,JJ2)=0.D0
C           D/DZ , D/DY
            BET(5,JJ )=0.D0
            BET(5,JJ1)=BET(3,JJ2)
            BET(5,JJ2)=BET(2,JJ1)
C           D/DZ , D/DX
            BET(6,JJ )=BET(3,JJ2)
            BET(6,JJ1)=0.D0
            BET(6,JJ2)=BET(1,JJ)
C
CS          ROTACIJE
CE          ROTATION PART
C
            G1=T*BET(1,JJ) +XJ(1,3)*H(I,1)
            G2=T*BET(2,JJ1)+XJ(2,3)*H(I,1) 
            G3=T*BET(3,JJ2)+XJ(3,3)*H(I,1)
            BET(1,JF )=GKS(1,1,I)*G1
            BET(1,JF1)=GKS(1,2,I)*G1
            BET(2,JF )=GKS(2,1,I)*G2
            BET(2,JF1)=GKS(2,2,I)*G2
            BET(3,JF )=GKS(3,1,I)*G3
            BET(3,JF1)=GKS(3,2,I)*G3
            BET(4,JF )=GKS(1,1,I)*G2+GKS(2,1,I)*G1
            BET(4,JF1)=GKS(1,2,I)*G2+GKS(2,2,I)*G1
            BET(5,JF )=GKS(2,1,I)*G3+GKS(3,1,I)*G2
            BET(5,JF1)=GKS(2,2,I)*G3+GKS(3,2,I)*G2
            BET(6,JF )=GKS(1,1,I)*G3+GKS(3,1,I)*G1
            BET(6,JF1)=GKS(1,2,I)*G3+GKS(3,2,I)*G1
         ELSE
C
CS          DEO VEZAN ZA TRANSLACIJE
CE          TRANSLATION PART
C
C           PT=1.D0-T*T
            PT=1.D0
C
            IF(INDSEL.EQ.1.AND.NGAUSU.LT.NGS12) GO TO 991
C           D/DR
            BET(1,JJ )=XJJ(1,1)*H(I,2)
            BET(1,JJ1)=XJJ(1,2)*H(I,2)
            BET(1,JJ2)=XJJ(1,3)*H(I,2)
C           D/DS
            BET(2,JJ )=XJJ(2,1)*H(I,3)
            BET(2,JJ1)=XJJ(2,2)*H(I,3)
            BET(2,JJ2)=XJJ(2,3)*H(I,3)
C           D/DT
            BET(3,JJ )=0.D0
            BET(3,JJ1)=0.D0
            BET(3,JJ2)=0.D0
C           D/DR , D/DS
            BET(4,JJ )=XJJ(1,1)*H(I,3)+XJJ(2,1)*H(I,2)
            BET(4,JJ1)=XJJ(1,2)*H(I,3)+XJJ(2,2)*H(I,2)
            BET(4,JJ2)=XJJ(1,3)*H(I,3)+XJJ(2,3)*H(I,2)
  991       IF(INDSEL.EQ.1.AND.NGAUSU.EQ.NGS12) GO TO 992
C           D/DS , D/DT
            BET(5,JJ )=XJJ(3,1)*H(I,3)*PT
            BET(5,JJ1)=XJJ(3,2)*H(I,3)*PT
            BET(5,JJ2)=XJJ(3,3)*H(I,3)*PT
C           D/DT , D/DR
            BET(6,JJ )=XJJ(3,1)*H(I,2)*PT
            BET(6,JJ1)=XJJ(3,2)*H(I,2)*PT
            BET(6,JJ2)=XJJ(3,3)*H(I,2)*PT
C
CS          ROTACIJE
CE          ROTATION PART
C
  992       CALL CLEAR(RJ,6)
            DO 20 JC=1,3
            DO 20 JZ=1,3
               RJ(1,JC)=RJ(1,JC)+GKS(JZ,1,I)*XJJ(JC,JZ)
               RJ(2,JC)=RJ(2,JC)+GKS(JZ,2,I)*XJJ(JC,JZ)
   20       CONTINUE
            IF(INDSEL.EQ.1.AND.NGAUSU.LT.NGS12) GO TO 993
C           D/DR
            BET(1,JF )=T*H(I,2)*RJ(1,1)
            BET(1,JF1)=T*H(I,2)*RJ(2,1)
C           D/DS
            BET(2,JF )=T*H(I,3)*RJ(1,2)
            BET(2,JF1)=T*H(I,3)*RJ(2,2)
C           D/DT
            BET(3,JF )=0.D0
C           BET(3,JF )=  H(I,1)*RJ(1,3)
            BET(3,JF1)=0.D0
C           BET(3,JF1)=  H(I,1)*RJ(2,3)
C           D/DR , D/DS
            BET(4,JF )=T*H(I,3)*RJ(1,1)+T*H(I,2)*RJ(1,2)
            BET(4,JF1)=T*H(I,3)*RJ(2,1)+T*H(I,2)*RJ(2,2)
  993       IF(INDSEL.EQ.1.AND.NGAUSU.EQ.NGS12) GO TO 994
C           D/DS , D/DT
            BET(5,JF )=(T*H(I,3)*RJ(1,3)+  H(I,1)*RJ(1,2))*PT
            BET(5,JF1)=(T*H(I,3)*RJ(2,3)+  H(I,1)*RJ(2,2))*PT
C           D/DT , D/DR
            BET(6,JF )=(T*H(I,2)*RJ(1,3)+  H(I,1)*RJ(1,1))*PT
            BET(6,JF1)=(T*H(I,2)*RJ(2,3)+  H(I,1)*RJ(2,1))*PT
C
  994       IF(IATYP.GT.1.AND.INDKOV.EQ.-1) THEN
            DO 25 KZ=1,9
            BNL(KZ,JJ )=0.D0
            BNL(KZ,JJ1)=0.D0
            BNL(KZ,JJ2)=0.D0
   25       BNL(KZ,JF2)=0.D0
C
CS          DEO VEZAN ZA TRANSLACIJE
CE          TRANSLATION PART
C
C           D/DR
            BNL(1,JJ )=H(I,2)
            BNL(2,JJ1)=H(I,2)
            BNL(3,JJ2)=H(I,2)
C           D/DS
            BNL(4,JJ )=H(I,3)
            BNL(5,JJ1)=H(I,3)
            BNL(6,JJ2)=H(I,3)
C
CS          ROTACIJE
CE          ROTATION PART
C
C           D/DR
            BNL(1,JF )=T*H(I,2)*GKS(1,1,I)
            BNL(1,JF1)=T*H(I,2)*GKS(1,2,I)
            BNL(2,JF )=T*H(I,2)*GKS(2,1,I)
            BNL(2,JF1)=T*H(I,2)*GKS(2,2,I)
            BNL(3,JF )=T*H(I,2)*GKS(3,1,I)
            BNL(3,JF1)=T*H(I,2)*GKS(3,2,I)
C           D/DS
            BNL(4,JF )=T*H(I,3)*GKS(1,1,I)
            BNL(4,JF1)=T*H(I,3)*GKS(1,2,I)
            BNL(5,JF )=T*H(I,3)*GKS(2,1,I)
            BNL(5,JF1)=T*H(I,3)*GKS(2,2,I)
            BNL(6,JF )=T*H(I,3)*GKS(3,1,I)
            BNL(6,JF1)=T*H(I,3)*GKS(3,2,I)
C           D/DT
            BNL(7,JF )=  H(I,1)*GKS(1,1,I)
            BNL(7,JF1)=  H(I,1)*GKS(1,2,I)
            BNL(8,JF )=  H(I,1)*GKS(2,1,I)
            BNL(8,JF1)=  H(I,1)*GKS(2,2,I)
            BNL(9,JF )=  H(I,1)*GKS(3,1,I)
            BNL(9,JF1)=  H(I,1)*GKS(3,2,I)
            ENDIF
         ENDIF
C         IF(MTR(I).EQ.1)THEN
            DO 30 KZ=1,6
   30       BET(KZ,JF2)=0.D0
C         ENDIF
   40 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UKDEFI(COR,CON,GM,GMT,UEL,H,R,S,T,
     1                  STRAIN,TRLN,TSS,IATYP,NCVE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE TOTAL STRAIN FOR T.L. AND U.L.
CS.   P R O G R A M
CS.      ZA RACUNANJE UKUPNIH DEFORMACIJA ZA T.L. I U.L.
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /CUVAJN/ TAUCEN(6,5),INDSEL,NGAUSU
C
      DIMENSION COR(9,*),CON(9,3),GM(3,*),GMT(3,*),H(9,*),UEL(*),
     1          STRAIN(*),TRLN(6,*),TSS(6,*),STRAIK(6)
      DIMENSION XT(3,3),XJJ(3,3)
      DIMENSION X0(3,3)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UKDEFI'
C
      IF(IATYP.EQ.2) THEN
CS       ZA T.L. KOORDINATE U TRENUTKU - T+DT
CE       FOR T.L. COORDINATE IN TIME - T+DT
         II=0
         DO 6 I=1,NCVE
         DO 6 K=1,3
            II=II+1
            CON(I,K)=COR(I,K)+UEL(II)
    6    CONTINUE
CS       KOVARIJANTNI BAZNI VEKTORI U KONFIGURACIJI - 0
CE       COVARIANT BASIS VECTORS IN CONFIGURATION - 0
         CALL GRGSGT(COR,GM,X0,H,T,NCVE)
CS       KOVARIJANTNI BAZNI VEKTORI U KONFIGURACIJI - T+DT
CE       COVARIANT BASIS VECTORS IN CONFIGURATION - T+DT
         CALL GRGSGT(CON,GMT,XT,H,T,NCVE)
         CALL JEDNA1(XJJ,XT,9)
CS       MATRICA TRANSF. DEFORMACIJA - TSS (GLOBALNI - LOKALNI DEKARTOV)
CE       STRAIN TRANSFORMATION MATRIX - TSS (GLOBAL - LOCAL CARTESIAN) 
         CALL TRANAL(XJJ,TSS,0)
      ELSE
CS       KOORDINATE U TRENUTKU - 0
CE       COORDINATE IN TIME - 0
         II=0
         DO 7 I=1,NCVE
         DO 7 K=1,3
            II=II+1
            CON(I,K)=COR(I,K)-UEL(II)
    7    CONTINUE
CS       KOVARIJANTNI BAZNI VEKTORI U KONFIGURACIJI - 0
CE       COVARIANT BASIS VECTORS IN CONFIGURATION - 0
         CALL GRGSGT(CON,GMT,X0,H,T,NCVE)
CS       KOVARIJANTNI BAZNI VEKTORI U KONFIGURACIJI - T+DT
CE       COVARIANT BASIS VECTORS IN CONFIGURATION - T+DT
         CALL GRGSGT(COR,GM,XT,H,T,NCVE)
      ENDIF
C
CS    IZOPARAMETARSKE UKUPNE KOVARIJANTNE DEFORMACIJE U GAUSOVOJ TACKI
CE    ISOPARAMETRIC TOTAL COVARIANT STRAIN IN GAUS POINT 
      CALL CLEAR(STRAIK,6)
      IF(INDSEL.EQ.1.AND.NGAUSU.LT.NGS12) GO TO 991
      STRAIK(1)=.5D0*(XT(1,1)*XT(1,1)+XT(1,2)*XT(1,2)+XT(1,3)*XT(1,3)
     1               -X0(1,1)*X0(1,1)-X0(1,2)*X0(1,2)-X0(1,3)*X0(1,3))
      STRAIK(2)=.5D0*(XT(2,1)*XT(2,1)+XT(2,2)*XT(2,2)+XT(2,3)*XT(2,3)
     1               -X0(2,1)*X0(2,1)-X0(2,2)*X0(2,2)-X0(2,3)*X0(2,3))
C      STRAIK(3)=.5D0*(XT(3,1)*XT(3,1)+XT(3,2)*XT(3,2)+XT(3,3)*XT(3,3)
C     1               -X0(3,1)*X0(3,1)-X0(3,2)*X0(3,2)-X0(3,3)*X0(3,3))
      STRAIK(3)=0.D0
      STRAIK(4)=XT(1,1)*XT(2,1)+XT(1,2)*XT(2,2)+XT(1,3)*XT(2,3)
     1         -X0(1,1)*X0(2,1)-X0(1,2)*X0(2,2)-X0(1,3)*X0(2,3) 
  991 IF(INDSEL.EQ.1.AND.NGAUSU.EQ.NGS12) GO TO 992
      STRAIK(5)=XT(3,1)*XT(2,1)+XT(3,2)*XT(2,2)+XT(3,3)*XT(2,3)
     1         -X0(3,1)*X0(2,1)-X0(3,2)*X0(2,2)-X0(3,3)*X0(2,3) 
      STRAIK(6)=XT(1,1)*XT(3,1)+XT(1,2)*XT(3,2)+XT(1,3)*XT(3,3)
     1         -X0(1,1)*X0(3,1)-X0(1,2)*X0(3,2)-X0(1,3)*X0(3,3) 
C     PT=1.D0-T*T
C     STRAIK(5)=STRAIK(5)*PT
C     STRAIK(6)=STRAIK(6)*PT
C
CS    TRANSFORMACIJA DEFORMACIJA (KOVARIJANTNE - GLOBALNI DEKARTOV)
CE    TRANSFORM STRAIN (COVARIANT - GLOBAL CARTESIAN) 
  992 CALL MNOZI1(STRAIN,TRLN,STRAIK,6,6)
      RETURN
      END
C=======================================================================
      SUBROUTINE GLASOR(PRINC,QP,P1,P2,P3,IMAX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION PRINC(*),QP(3,*),PRAV(3),POM(3)
C
      PRAV(1)=P1
      PRAV(2)=P2
      PRAV(3)=P3
      IMAX=0
      SMAX=0.D0
      DO 20 J=1,3
         DO 10 I=1,3
            POM(I)=QP(J,I)
  10     CONTINUE
         SKAL=DOT(POM,PRAV,3)
         SKAL=DABS(SKAL)
         IF(SKAL.GT.SMAX) THEN
            IMAX=J
            SMAX=SKAL
         ENDIF
  20  CONTINUE
C      IF(IMAX.EQ.0.OR.IMAX.EQ.1) WRITE(*,*) 'PROVERI IMAX=',IMAX
      IF(IMAX.EQ.3) RETURN
C
      IVAR=3
C     TRECA  VARIJANTA
      IF(IVAR.EQ.3) THEN
         DO 33 I=1,3
            DUM=QP(IMAX,I)
            QP(IMAX,I)=QP(3,I)
            QP(3,I)=DUM
   33    CONTINUE
         DUM=PRINC(IMAX)
         PRINC(IMAX)=PRINC(3)
         PRINC(3)=DUM
      ENDIF
C     DRUGA  VARIJANTA
      IF(IVAR.EQ.2) THEN
         IF(IMAX.EQ.2) THEN
            DO 32 I=1,3
               DUM=QP(2,I)
               QP(2,I)=QP(3,I)
               QP(3,I)=DUM
   32       CONTINUE
            DUM=PRINC(2)
            PRINC(2)=PRINC(3)
            PRINC(3)=DUM
         ENDIF
         IF(IMAX.EQ.1) THEN
            DO 42 I=1,3
               DUM=QP(1,I)
               QP(1,I)=QP(2,I)
               QP(2,I)=QP(3,I)
               QP(3,I)=DUM
   42       CONTINUE
            DUM=PRINC(1)
            PRINC(1)=PRINC(2)
            PRINC(2)=PRINC(3)
            PRINC(3)=DUM
         ENDIF
      ENDIF
C
C     PRVA VARIJANTA
      IF(IVAR.EQ.1) THEN
         IF(IMAX.EQ.2) THEN
            DO 31 I=1,3
               DUM=QP(2,I)
               QP(2,I)=QP(1,I)
               QP(1,I)=QP(3,I)
               QP(3,I)=DUM
   31       CONTINUE
            DUM=PRINC(2)
            PRINC(2)=PRINC(1)
            PRINC(1)=PRINC(3)
            PRINC(3)=DUM
         ENDIF
         IF(IMAX.EQ.1) THEN
            DO 41 I=1,3
               DUM=QP(1,I)
               QP(1,I)=QP(2,I)
               QP(2,I)=QP(3,I)
               QP(3,I)=DUM
   41       CONTINUE
            DUM=PRINC(1)
            PRINC(1)=PRINC(2)
            PRINC(2)=PRINC(3)
            PRINC(3)=DUM
         ENDIF
      ENDIF
      RETURN
      END
