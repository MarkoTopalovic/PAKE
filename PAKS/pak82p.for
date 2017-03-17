C=======================================================================
C
C=======================================================================
      SUBROUTINE PINSKY(TAU,SKE,FTDT,T,WTU,GKS,LM,MSL,NOP,DRGTM0,DEB,
     1                  TGT,LPLAS,LPLA1,ISNA,MTR,TTR,IND6,R,S,COR,
     1                  N45,DEF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO FORM ELEMENT MATRIX FOR SHELL ELEMENT - PINSKY & JUNG
CS.   P R O G R A M
CS.      ZA FORMIRANJE MATRICA ELEMENATA LJUSKE - PINSKY & JUNG
C .
C .      IPODT.EQ.2; (9 CVOROVA)
C .
C ......................................................................
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /IZLE4B/ H(9,3),GM(3,9),BLT(6,54),BE(9,54),ETP(6,6),UEL(54)
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /INDNAP/ NAPON
      COMMON /FILTER/ TOLNAP,TOLPOM
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /PINSKI/ B0(10,54),B1(10,54)
      COMMON /MATIZO/ E,V
      COMMON /MATANI/ EX,EY,EZ,VXY,VYZ,VZX,GXY,GYZ,GZX
C
      DIMENSION NOP(NE,*),SKE(*),LM(*),FTDT(*),DEF(N45,NGS12,NE,*),
     2TAU(N45,NGS12,NE,*),TA(6),STRAIN(6),STRESS(6)
      DIMENSION GKS(3,2,*)
      DIMENSION DRGTM0(9,*),DEB(*),BNL(9,54)
      DIMENSION BM(3,54),BB(3,54),BS(2,54),HN(3,4),TRLN(6,6)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' PINSKY'
C
CS    JAKOBIJAN - XJJ, INVERTOVAN JAKOBIJAN - XJ
CE    JACOBIAN - XJJ, INVERSE JACOBIAN - XJ
C
      CALL JEDNA1(XJJ,XJ,9)  
      CALL MINV3(XJ,DET)
CZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
CS    RACUNANJE KONSTANTNOG DELA MATRICE B
      IF(JG.EQ.1) CALL B0B1(B0,B1,COR,GKS,GM)
CS     MATRICA TRANSFORMACIJE TRLN I INTERPOLACIJSKE FUNKCIJE HN
CS     JEDNACINE (51), (67) I (19), (25)
       CALL TRLNHN(HN,R,S,TRLN,XJ)
CS     MATRICE BM, BB I BS U GAUS TACKI R I S U LOKALNOM SISTEMU
CS     JEDNACINE (59),(60),(81); (90),(91),(92); (96),(97),(98),(99)
       CALL BMBBBS(BLT,BM,BB,BS,R,S,T,B0,B1,HN)
CS     MATRICE BM, BB I BS U GAUS TACKI R I S
CS     JEDNACINE (119), (89) I (95)
        CALL TRANSB(BLT,TRLN)
CZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ       
C
CS    FORMIRANJE MATRICE BL1 GEOMETRIJSKA NELINEARNOST - T.L.
CE    FORM MATRIX BL1 GEOMETRICAL NONLINEARITY -  T.L.
C
      IF(IATYP.EQ.2)
     1CALL BETL18(BLT,NOP,LM,UEL,STRAIN,H,DRGTM0,GKS,T,DEB(NLM),0)
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
      CALL TRANAL(XJJ,TSG,0)
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
C
CS       GEOMETRIJSKA NELINEARNOST - UKUPNE DEFORMACIJE ZA T.L. I U.L.
CE       GEOMETRICAL NONLINEARITY - TOTAL STRAIN FOR T.L. AND U.L.
C
         IF(IATYP.EQ.2)
     1   CALL BETL18(BLT,NOP,LM,UEL,STRAIN,H,DRGTM0,GKS,T,DEB(NLM),1)
         IF(IATYP.EQ.3)
     1   CALL STRUL8(H,NOP,LM,DRGTM0,UEL,STRAIN,T,DEB(NLM))
C
CS       TRANSFORMACIJA DEFORMACIJA (GLOBALNI - LOKALNI DEKARTOV)
CE       TRANSFORM STRAIN (GLOBAL - LOCAL CARTESIAN) 
C
         CALL CLEAR(TA,6) 
         CALL MNOZI1(TA,TSG,STRAIN,6,6)
         TA(3)=0.D0
         CALL JEDNA1(STRAIN,TA,6)
C
CS       TERMICKE DEFORMACIJE   ETH=ALFA*(T-T0)
CE       THERMAL STRAINS
C
         IF(NMODM.EQ.3.OR.NMODM.EQ.4) THEN
            DTGT0=TGT-TEMP0
            STRAIN(1)=STRAIN(1)-ALFA(1)*DTGT0
            STRAIN(2)=STRAIN(2)-ALFA(2)*DTGT0
            IF(NAPON.EQ.1) STRAIN(3)=ALFA(3)*DTGT0
            STRAIN(4)=STRAIN(4)-ALFA(4)*DTGT0
         ENDIF
C          CALL WRR(STRAIN,6,'STRA')
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
C           CALL WRR(STRESS,6,'STRE')
C
CS          TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM
CE          TRANSFORM STRESS TO GLOBAL COORDS.
C
            CALL CLEAR(TA,6) 
            CALL MNOZI2(TA,TSG,STRESS,6,6)
            IF(ISNA.EQ.2) CALL JEDNA1(STRESS,TA,6)
C
CS          CISCENJE NUMERICKIH GRESAKA ZA NAPONE
CE          CLEANING NUMERICAL ERRORS FOR STRESS
C
C           IF(IATYP.EQ.0) CALL CISTIN(STRESS,6)
C
            CALL JEDNA1(TAU(1,JG,NLM,MSL),STRESS,6)
C
            IF(NAPON.EQ.1) THEN
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
         ELSE
C
CS          NAPONI ZA PLASTICAN MODEL
CE          STRESS FOR MATERIAL NONLINEARITY 
C
            IRAC=1
            CALL MODMA8(STRAIN,STRESS,NMODM,IRAC,LPLAS,LPLA1,TGT)
C
CS          TRANSFORMACIJA NAPONA NA GLOBALNI SISTEM
CE          TRANSFORM STRESS TO GLOBAL COORDS.
C
            CALL CLEAR(TA,6) 
            CALL MNOZI2(TA,TSG,STRESS,6,6)
         ENDIF
C
C        DEFORMACIJE ZBOG KOREKCIJE DEBLJINE
C          CALL WRR(ETP,36,'ELAS')
         IF(NAPON.EQ.1.AND.NMODM.LT.5) THEN
            CALL JEDNA1(DEF(1,JG,NLM,MSL),STRAIN,6)
C            CALL WRR(STRAIN,6,'STRA')
         ENDIF
C
        IF(NGENL.GT.0.OR.(NZADP.GT.0.AND.ISKNP.EQ.2.AND.NGENL.EQ.0))THEN
C
CS          RACUNANJE UNUTRASNJIH SILA
CE          CALCULATE INTERNAL FORCES
C           r = BT * S 
C
            CALL INTEGF(FTDT,BLT,TA,LM,WTU,ND,6)
C
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
      CALL INTEGK(SKE,BLT,ELAST,LM,WTU,ND,6)
C
CS    GEOMETRIJSKI NELINEARAN DEO MATRICE KRUTOSTI
CE    GEOMETRIC NONLINEAR PART OF STIFFNESS MATRIX   
C
      IF(IATYP.GT.1)
     1CALL KNL8(SKE,H,NOP,LM,TA,T,WTU,BNL,GKS,ND,MTR,TTR,IND6)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE INTERP(HRR,HRL,HSR,HSL,HR,HL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION HRR(3,*),HRL(3,*),HSR(3,*),HSL(3,*),HR(3,*),HL(3,*)
      Q=1.D0/DSQRT(3.D0)
      QP=Q*(1+Q)
      QM=Q*(1-Q)
      Q2=1-Q*Q
      HR(1,1)=-.5D0*QM
      HR(2,1)=Q2
      HR(3,1)= .5D0*QP
      HR(1,2)=-.5D0+Q
      HR(2,2)=-2.D0*Q
      HR(3,2)= .5D0+Q
      HRR(1,1)= .75D0*QM
      HRR(1,2)=      -QM
      HRR(1,3)= .25D0*QM
      HRR(1,4)=-1.5D0*Q2
      HRR(1,5)= 2.0D0*Q2
      HRR(1,6)= -.5D0*Q2
      HRR(1,7)=-.75D0*QP
      HRR(1,8)=       QP
      HRR(1,9)=-.25D0*QP
      HRR(2,1)= HRR(1,3) 
      HRR(2,2)= 0.0D0
      HRR(2,3)=-HRR(1,3)
      HRR(2,4)= HRR(1,6)
      HRR(2,5)= 0.0D0
      HRR(2,6)=-HRR(1,6)
      HRR(2,7)= HRR(1,9)
      HRR(2,8)= 0.0D0
      HRR(2,9)=-HRR(1,9)
      HRR(3,1)=-HRR(1,3) 
      HRR(3,2)=-HRR(1,2)
      HRR(3,3)=-HRR(1,1)
      HRR(3,4)=-HRR(1,6)
      HRR(3,5)=-HRR(1,5)
      HRR(3,6)=-HRR(1,4)
      HRR(3,7)=-HRR(1,9)
      HRR(3,8)=-HRR(1,8)
      HRR(3,9)=-HRR(1,7)
      HSR(1,1)= HRR(1,1) 
      HSR(1,2)= HRR(1,4)
      HSR(1,3)= HRR(1,7)
      HSR(1,4)= HRR(1,2)
      HSR(1,5)= HRR(1,5)
      HSR(1,6)= HRR(1,8)
      HSR(1,7)= HRR(1,3)
      HSR(1,8)= HRR(1,6)
      HSR(1,9)= HRR(1,9)
      HSR(2,1)= HRR(2,1) 
      HSR(2,2)= HRR(2,4)
      HSR(2,3)= HRR(2,7)
      HSR(2,4)= HRR(2,2)
      HSR(2,5)= HRR(2,5)
      HSR(2,6)= HRR(2,8)
      HSR(2,7)= HRR(2,3)
      HSR(2,8)= HRR(2,6)
      HSR(2,9)= HRR(2,9)
      HSR(3,1)= HRR(3,1) 
      HSR(3,2)= HRR(3,4)
      HSR(3,3)= HRR(3,7)
      HSR(3,4)= HRR(3,2)
      HSR(3,5)= HRR(3,5)
      HSR(3,6)= HRR(3,8)
      HSR(3,7)= HRR(3,3)
      HSR(3,8)= HRR(3,6)
      HSR(3,9)= HRR(3,9)
      Q=-1.D0/DSQRT(3.D0)
      QP=Q*(1+Q)
      QM=Q*(1-Q)
      Q2=1-Q*Q
      HL(1,1)=-.5D0*QM
      HL(2,1)=Q2
      HL(3,1)= .5D0*QP
      HL(1,2)=-.5D0+Q
      HL(2,2)=-2.D0*Q
      HL(3,2)= .5D0+Q
      HRL(1,1)= .75D0*QM
      HRL(1,2)=      -QM
      HRL(1,3)= .25D0*QM
      HRL(1,4)=-1.5D0*Q2
      HRL(1,5)= 2.0D0*Q2
      HRL(1,6)= -.5D0*Q2
      HRL(1,7)=-.75D0*QP
      HRL(1,8)=       QP
      HRL(1,9)=-.25D0*QP
      HRL(2,1)= HRL(1,3) 
      HRL(2,2)= 0.0D0
      HRL(2,3)=-HRL(1,3)
      HRL(2,4)= HRL(1,6)
      HRL(2,5)= 0.0D0
      HRL(2,6)=-HRL(1,6)
      HRL(2,7)= HRL(1,9)
      HRL(2,8)= 0.0D0
      HRL(2,9)=-HRL(1,9)
      HRL(3,1)=-HRL(1,3) 
      HRL(3,2)=-HRL(1,2)
      HRL(3,3)=-HRL(1,1)
      HRL(3,4)=-HRL(1,6)
      HRL(3,5)=-HRL(1,5)
      HRL(3,6)=-HRL(1,4)
      HRL(3,7)=-HRL(1,9)
      HRL(3,8)=-HRL(1,8)
      HRL(3,9)=-HRL(1,7)
      HSL(1,1)= HRL(1,1) 
      HSL(1,2)= HRL(1,4)
      HSL(1,3)= HRL(1,7)
      HSL(1,4)= HRL(1,2)
      HSL(1,5)= HRL(1,5)
      HSL(1,6)= HRL(1,8)
      HSL(1,7)= HRL(1,3)
      HSL(1,8)= HRL(1,6)
      HSL(1,9)= HRL(1,9)
      HSL(2,1)= HRL(2,1) 
      HSL(2,2)= HRL(2,4)
      HSL(2,3)= HRL(2,7)
      HSL(2,4)= HRL(2,2)
      HSL(2,5)= HRL(2,5)
      HSL(2,6)= HRL(2,8)
      HSL(2,7)= HRL(2,3)
      HSL(2,8)= HRL(2,6)
      HSL(2,9)= HRL(2,9)
      HSL(3,1)= HRL(3,1) 
      HSL(3,2)= HRL(3,4)
      HSL(3,3)= HRL(3,7)
      HSL(3,4)= HRL(3,2)
      HSL(3,5)= HRL(3,5)
      HSL(3,6)= HRL(3,8)
      HSL(3,7)= HRL(3,3)
      HSL(3,8)= HRL(3,6)
      HSL(3,9)= HRL(3,9)
      RETURN
      END
C=======================================================================
      SUBROUTINE B0B1(B0,B1,COR,GKS,GM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     KOEFICIJENTI T I PRIRODNE KOORDINATE U CVOROIMA I NA LINIJAMA
CS     JEDNACINE (45)
C
      DIMENSION B0(10,*),B1(10,*),COR(9,*),GKS(3,2,*),GM(3,*) 
      DIMENSION RCZR(3,3),RCZL(3,3),SCZR(3,3),SCZL(3,3)
      DIMENSION ZRCR(3,3),ZRCL(3,3),ZSCR(3,3),ZSCL(3,3)
      DIMENSION ZCRR(3,3),ZCRL(3,3),ZCSR(3,3),ZCSL(3,3)
      DIMENSION RCR(3,3),RCL(3,3),SCR(3,3),SCL(3,3)
      DIMENSION CRR(3,3),CRL(3,3),CSR(3,3),CSL(3,3)
      DIMENSION HRR(3,9),HRL(3,9),HSR(3,9),HSL(3,9),HR(3,2),HL(3,2)
C
      CALL INTERP(HRR,HRL,HSR,HSL,HR,HL)
C
      DO 1 J=1,3
         J1=(J-1)*3+1
         J2=J1+1
         J3=J1+2
         J4=J +3
         J7=J +6
      DO 1 K=1,3
         RCR(K,J)=HR(1,2)*COR(J1,K)+HR(2,2)*COR(J2,K)+HR(3,2)*COR(J3,K)
         RCL(K,J)=HL(1,2)*COR(J1,K)+HL(2,2)*COR(J2,K)+HL(3,2)*COR(J3,K)
         SCR(K,J)=HR(1,2)*COR(J ,K)+HR(2,2)*COR(J4,K)+HR(3,2)*COR(J7,K)
         SCL(K,J)=HL(1,2)*COR(J ,K)+HL(2,2)*COR(J4,K)+HL(3,2)*COR(J7,K)
         ZRCR(K,J)=HR(1,2)*GM(K,J1)+HR(2,2)*GM(K,J2)+HR(3,2)*GM(K,J3)
         ZRCL(K,J)=HL(1,2)*GM(K,J1)+HL(2,2)*GM(K,J2)+HL(3,2)*GM(K,J3)
         ZSCR(K,J)=HR(1,2)*GM(K ,J)+HR(2,2)*GM(K,J4)+HR(3,2)*GM(K,J7)
         ZSCL(K,J)=HL(1,2)*GM(K ,J)+HL(2,2)*GM(K,J4)+HL(3,2)*GM(K,J7)
         RCZR(K,J)=HR(1,1)*GM(K,J1)+HR(2,1)*GM(K,J2)+HR(3,1)*GM(K,J3)
         RCZL(K,J)=HL(1,1)*GM(K,J1)+HL(2,1)*GM(K,J2)+HL(3,1)*GM(K,J3)
         SCZR(K,J)=HR(1,1)*GM(K ,J)+HR(2,1)*GM(K,J4)+HR(3,1)*GM(K,J7)
         SCZL(K,J)=HL(1,1)*GM(K ,J)+HL(2,1)*GM(K,J4)+HL(3,1)*GM(K,J7)
         CSR(K,J)=0.D0
         CSL(K,J)=0.D0
         CRR(K,J)=0.D0
         CRL(K,J)=0.D0
         ZCSR(K,J)=0.D0
         ZCSL(K,J)=0.D0
         ZCRR(K,J)=0.D0
         ZCRL(K,J)=0.D0
      DO 1 I=1,9
         CSR(K,J)=CSR(K,J)+HSR(J,I)*COR(I,K)
         CSL(K,J)=CSL(K,J)+HSL(J,I)*COR(I,K)
         CRR(K,J)=CRR(K,J)+HRR(J,I)*COR(I,K)
         CRL(K,J)=CRL(K,J)+HRL(J,I)*COR(I,K)
         ZCSR(K,J)=ZCSR(K,J)+HSR(J,I)*GM(K,I)
         ZCSL(K,J)=ZCSL(K,J)+HSL(J,I)*GM(K,I)
         ZCRR(K,J)=ZCRR(K,J)+HRR(J,I)*GM(K,I)
         ZCRL(K,J)=ZCRL(K,J)+HRL(J,I)*GM(K,I)
    1 CONTINUE
      Q32=DSQRT(3.D0)/2.D0
      II=0
      IJ=0
      DO 10 I=1,3
      DO 10 J=1,3
         IJ=IJ+1
      DO 20 K=1,3
         II=II+1
         I27=II+27
         B0(1,I27)=0.D0
         B1(1,I27)=0.D0
         B0(2,I27)=0.D0
         B1(2,I27)=0.D0
         B0(3,I27)=0.D0
         B1(3,I27)=0.D0
         B0(4,I27)=0.D0
         B1(4,I27)=0.D0
         B0(1,II)=.5D0 *(HR(J,2)*RCR(K,I)+HL(J,2)*RCL(K,I))
         B1(1,II)= Q32 *(HR(J,2)*RCR(K,I)-HL(J,2)*RCL(K,I))
         B0(2,II)=.5D0 *(HR(I,2)*SCR(K,J)+HL(I,2)*SCL(K,J))
         B1(2,II)= Q32 *(HR(I,2)*SCR(K,J)-HL(I,2)*SCL(K,J))
         B0(3,II)=.5D0 *(HR(J,2)*CSR(K,I)+HL(J,2)*CSL(K,I))
         B1(3,II)= Q32 *(HR(J,2)*CSR(K,I)-HL(J,2)*CSL(K,I))
         B0(4,II)=.5D0 *(HR(I,2)*CRR(K,J)+HL(I,2)*CRL(K,J))
         B1(4,II)= Q32 *(HR(I,2)*CRR(K,J)-HL(I,2)*CRL(K,J))
         B0(5,II)=.5D0 *(HR(J,2)*ZRCR(K,I)+HL(J,2)*ZRCL(K,I))
         B1(5,II)= Q32 *(HR(J,2)*ZRCR(K,I)-HL(J,2)*ZRCL(K,I))
         B0(6,II)=.5D0 *(HR(I,2)*ZSCR(K,J)+HL(I,2)*ZSCL(K,J))
         B1(6,II)= Q32 *(HR(I,2)*ZSCR(K,J)-HL(I,2)*ZSCL(K,J))
         B0(7,II)=.5D0 *(HR(J,2)*ZCSR(K,I)+HL(J,2)*ZCSL(K,I))
         B1(7,II)= Q32 *(HR(J,2)*ZCSR(K,I)-HL(J,2)*ZCSL(K,I))
         B0(8,II)=.5D0 *(HR(I,2)*ZCRR(K,J)+HL(I,2)*ZCRL(K,J))
         B1(8,II)= Q32 *(HR(I,2)*ZCRR(K,J)-HL(I,2)*ZCRL(K,J))
         B0(9,II)=.5D0 *(HR(J,2)*RCZR(K,I)+HL(J,2)*RCZL(K,I))
         B1(9,II)= Q32 *(HR(J,2)*RCZR(K,I)-HL(J,2)*RCZL(K,I))
         B0(10,II)=.5D0*(HR(I,2)*SCZR(K,J)+HL(I,2)*SCZL(K,J))
         B1(10,II)= Q32*(HR(I,2)*SCZR(K,J)-HL(I,2)*SCZL(K,J))
   20 CONTINUE
      I27=II+25
      I28=I27+1
      I29=I27+2
      RRR=GKS(1,1,IJ)*RCR(1,I)+GKS(2,1,IJ)*RCR(2,I)+GKS(3,1,IJ)*RCR(3,I)
      RRL=GKS(1,1,IJ)*RCL(1,I)+GKS(2,1,IJ)*RCL(2,I)+GKS(3,1,IJ)*RCL(3,I)
         B0(5,I27)=.5D0 *(HR(J,2)*RRR+HL(J,2)*RRL)
         B1(5,I27)= Q32 *(HR(J,2)*RRR-HL(J,2)*RRL)
         B0(9,I27)=.5D0 *(HR(J,1)*RRR+HL(J,1)*RRL)
         B1(9,I27)= Q32 *(HR(J,1)*RRR-HL(J,1)*RRL)
      RRR=GKS(1,2,IJ)*RCR(1,I)+GKS(2,2,IJ)*RCR(2,I)+GKS(3,2,IJ)*RCR(3,I)
      RRL=GKS(1,2,IJ)*RCL(1,I)+GKS(2,2,IJ)*RCL(2,I)+GKS(3,2,IJ)*RCL(3,I)
         B0(5,I28)=.5D0 *(HR(J,2)*RRR+HL(J,2)*RRL)
         B1(5,I28)= Q32 *(HR(J,2)*RRR-HL(J,2)*RRL)
         B0(5,I29)=0.D0
         B1(5,I29)=0.D0
         B0(9,I28)=.5D0 *(HR(J,1)*RRR+HL(J,1)*RRL)
         B1(9,I28)= Q32 *(HR(J,1)*RRR-HL(J,1)*RRL)
         B0(9,I29)=0.D0
         B1(9,I29)=0.D0
      RRR=GKS(1,1,IJ)*SCR(1,J)+GKS(2,1,IJ)*SCR(2,J)+GKS(3,1,IJ)*SCR(3,J)
      RRL=GKS(1,1,IJ)*SCL(1,J)+GKS(2,1,IJ)*SCL(2,J)+GKS(3,1,IJ)*SCL(3,J)
         B0(6,I27)=.5D0 *(HR(I,2)*RRR+HL(I,2)*RRL)
         B1(6,I27)= Q32 *(HR(I,2)*RRR-HL(I,2)*RRL)
         B0(10,I27)=.5D0*(HR(I,1)*RRR+HL(I,1)*RRL)
         B1(10,I27)= Q32*(HR(I,1)*RRR-HL(I,1)*RRL)
      RRR=GKS(1,2,IJ)*SCR(1,J)+GKS(2,2,IJ)*SCR(2,J)+GKS(3,2,IJ)*SCR(3,J)
      RRL=GKS(1,2,IJ)*SCL(1,J)+GKS(2,2,IJ)*SCL(2,J)+GKS(3,2,IJ)*SCL(3,J)
         B0(6,I28)=.5D0 *(HR(I,2)*RRR+HL(I,2)*RRL)
         B1(6,I28)= Q32 *(HR(I,2)*RRR-HL(I,2)*RRL)
         B0(6,I29)=0.D0
         B1(6,I29)=0.D0
         B0(10,I28)=.5D0*(HR(I,1)*RRR+HL(I,1)*RRL)
         B1(10,I28)= Q32*(HR(I,1)*RRR-HL(I,1)*RRL)
         B0(10,I29)=0.D0
         B1(10,I29)=0.D0
      RRR=GKS(1,1,IJ)*CSR(1,I)+GKS(2,1,IJ)*CSR(2,I)+GKS(3,1,IJ)*CSR(3,I)
      RRL=GKS(1,1,IJ)*CSL(1,I)+GKS(2,1,IJ)*CSL(2,I)+GKS(3,1,IJ)*CSL(3,I)
         B0(7,I27)=.5D0 *(HR(J,2)*RRR+HL(J,2)*RRL)
         B1(7,I27)= Q32 *(HR(J,2)*RRR-HL(J,2)*RRL)
      RRR=GKS(1,2,IJ)*CSR(1,I)+GKS(2,2,IJ)*CSR(2,I)+GKS(3,2,IJ)*CSR(3,I)
      RRL=GKS(1,2,IJ)*CSL(1,I)+GKS(2,2,IJ)*CSL(2,I)+GKS(3,2,IJ)*CSL(3,I)
         B0(7,I28)=.5D0 *(HR(J,2)*RRR+HL(J,2)*RRL)
         B1(7,I28)= Q32 *(HR(J,2)*RRR-HL(J,2)*RRL)
         B0(7,I29)=0.D0
         B1(7,I29)=0.D0
      RRR=GKS(1,1,IJ)*CRR(1,J)+GKS(2,1,IJ)*CRR(2,J)+GKS(3,1,IJ)*CRR(3,J)
      RRL=GKS(1,1,IJ)*CRL(1,J)+GKS(2,1,IJ)*CRL(2,J)+GKS(3,1,IJ)*CRL(3,J)
         B0(8,I27)=.5D0 *(HR(I,2)*RRR+HL(I,2)*RRL)
         B1(8,I27)= Q32 *(HR(I,2)*RRR-HL(I,2)*RRL)
      RRR=GKS(1,2,IJ)*CRR(1,J)+GKS(2,2,IJ)*CRR(2,J)+GKS(3,2,IJ)*CRR(3,J)
      RRL=GKS(1,2,IJ)*CRL(1,J)+GKS(2,2,IJ)*CRL(2,J)+GKS(3,2,IJ)*CRL(3,J)
         B0(8,I28)=.5D0 *(HR(I,2)*RRR+HL(I,2)*RRL)
         B1(8,I28)= Q32 *(HR(I,2)*RRR-HL(I,2)*RRL)
         B0(8,I29)=0.D0
         B1(8,I29)=0.D0
   10 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE TRLNHN(HN,R,S,TRLN,XJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     MATRICA TRANSFORMACIJE TRLN I INTERPOLACIJSKE FUNKCIJE HN
CS     JEDNACINE (51), (67) I (19), (25)
C
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
C
      DIMENSION HN(3,*),TRLN(6,*),XJ(3,*)
C
      DO 10 II=1,2
         RS=R
         IF(II.EQ.2) RS=S
         I=(II-1)*2+1
         I1=I+1
         HN(1,I)=-.5D0*RS*(1.D0-RS)
         HN(2,I)=1.-RS*RS
         HN(3,I)=.5D0*RS*(1.D0+RS)
         HN(1,I1)=1.D0/6.D0-.5D0*RS
         HN(2,I1)=2.D0/3.D0
         HN(3,I1)=1.D0/6.D0+.5D0*RS
   10 CONTINUE
C
      TRLN(1,1)=XJ(1,1)*XJ(1,1)
      TRLN(2,1)=XJ(2,1)*XJ(2,1)
      TRLN(3,1)=XJ(3,1)*XJ(3,1)
      TRLN(4,1)=2.D0*XJ(1,1)*XJ(2,1)
      TRLN(5,1)=2.D0*XJ(2,1)*XJ(3,1)
      TRLN(6,1)=2.D0*XJ(3,1)*XJ(1,1)
C
      TRLN(1,2)=XJ(1,2)*XJ(1,2)
      TRLN(2,2)=XJ(2,2)*XJ(2,2)
      TRLN(3,2)=XJ(3,2)*XJ(3,2)
      TRLN(4,2)=2.D0*XJ(1,2)*XJ(2,2)
      TRLN(5,2)=2.D0*XJ(2,2)*XJ(3,2)
      TRLN(6,2)=2.D0*XJ(3,2)*XJ(1,2)
C
      TRLN(1,3)=XJ(1,3)*XJ(1,3)
      TRLN(2,3)=XJ(2,3)*XJ(2,3)
      TRLN(3,3)=XJ(3,3)*XJ(3,3)
      TRLN(4,3)=2.D0*XJ(1,3)*XJ(2,3)
      TRLN(5,3)=2.D0*XJ(2,3)*XJ(3,3)
      TRLN(6,3)=2.D0*XJ(3,3)*XJ(1,3)
C
      TRLN(1,4)=XJ(1,1)*XJ(1,2)
      TRLN(2,4)=XJ(2,1)*XJ(2,2)
      TRLN(3,4)=XJ(3,1)*XJ(3,2)
      TRLN(4,4)=XJ(1,1)*XJ(2,2)+XJ(1,2)*XJ(2,1)
      TRLN(5,4)=XJ(2,1)*XJ(3,2)+XJ(2,2)*XJ(3,1)
      TRLN(6,4)=XJ(3,1)*XJ(1,2)+XJ(3,2)*XJ(1,1)
C
      TRLN(1,5)=XJ(1,2)*XJ(1,3)
      TRLN(2,5)=XJ(2,2)*XJ(2,3)
      TRLN(3,5)=XJ(3,2)*XJ(3,3)
      TRLN(4,5)=XJ(1,2)*XJ(2,3)+XJ(1,3)*XJ(2,2)
      TRLN(5,5)=XJ(2,2)*XJ(3,3)+XJ(2,3)*XJ(3,2)
      TRLN(6,5)=XJ(3,2)*XJ(1,3)+XJ(3,3)*XJ(1,2)
C
      TRLN(1,6)=XJ(1,1)*XJ(1,3)
      TRLN(2,6)=XJ(2,1)*XJ(2,3)
      TRLN(3,6)=XJ(3,1)*XJ(3,3)
      TRLN(4,6)=XJ(1,1)*XJ(2,3)+XJ(1,3)*XJ(2,1)
      TRLN(5,6)=XJ(2,1)*XJ(3,3)+XJ(2,3)*XJ(3,1)
      TRLN(6,6)=XJ(3,1)*XJ(1,3)+XJ(3,3)*XJ(1,1)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE BMBBBS(BLT,BM,BB,BS,R,S,T,B0,B1,HN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    MATRICE BM, BB I BS U GAUS TACKI R I S
CS    JEDNACINE (59),(60),(81); (90),(91),(92); (96),(97),(98),(99)
C
      COMMON /INDCEL/ IND3D
C
      DIMENSION BM(3,*),BB(3,*),BS(2,*),B0(10,*),B1(10,*),HN(3,*),
     1          BLT(6,*)
C
CS    JEDNACINE (59); (90); (96)
      II=0
      DO 10 I=1,3
      DO 10 J=1,3
      DO 10 K=1,3
         II=II+1
         I27=II+27
         BM(1,II)=HN(I,3)*(B0(1,II)+R*B1(1,II))
         BM(2,II)=HN(J,1)*(B0(2,II)+S*B1(2,II))
         BM(3,II)=HN(I,4)*(B0(3,II)+R*B1(3,II))
     $           +HN(J,2)*(B0(4,II)+S*B1(4,II))
         BM(1,I27)=0.D0
         BM(2,I27)=0.D0
         BM(3,I27)=0.D0
         BB(1,II)=HN(I,3)*(B0(5,II)+R*B1(5,II))*T
         BB(2,II)=HN(J,1)*(B0(6,II)+S*B1(6,II))*T
         BB(3,II)=HN(I,4)*(B0(7,II)+R*B1(7,II))*T
     $           +HN(J,2)*(B0(8,II)+S*B1(8,II))*T
         BB(1,I27)=HN(I,3)*(B0(5,I27)+R*B1(5,I27))*T
         BB(2,I27)=HN(J,1)*(B0(6,I27)+S*B1(6,I27))*T
         BB(3,I27)=HN(I,4)*(B0(7,I27)+R*B1(7,I27))*T
     $            +HN(J,2)*(B0(8,I27)+S*B1(8,I27))*T
         BS(1,II)=HN(I,3)*(B0(9,II)+R*B1(9,II))
         BS(2,II)=HN(J,1)*(B0(10,II)+S*B1(10,II))
         BS(1,I27)=HN(I,3)*(B0(9,I27)+R*B1(9,I27))
         BS(2,I27)=HN(J,1)*(B0(10,I27)+S*B1(10,I27))
   10 CONTINUE
C
      DO 20 I=1,54
         BLT(1,I)=BM(1,I)+BB(1,I)
         BLT(2,I)=BM(2,I)+BB(2,I)
         BLT(3,I)=0.D0
         BLT(4,I)=BM(3,I)+BB(3,I)
         BLT(5,I)=BS(2,I)
         BLT(6,I)=BS(1,I)
         IF(IND3D.NE.0) THEN
            BLT(5,I)=0.D0
            BLT(6,I)=0.D0
         ENDIF
         IF(IND3D.EQ.2) BLT(1,I)=0.D0
         IF(IND3D.EQ.2) BLT(2,I)=0.D0
   20 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE TRANSB(BLT,TRLN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     MATRICE BM, BB I BS U GAUS TACKI R I S U GLOBALNOM SISTEMU
CS     JEDNACINE (119), (121) I (123)
C
      DIMENSION BLT(6,*)
C ------------------------------
      DIMENSION TRLN(6,6),POM(6)
C ------------------------------
C
C     DRUGI SLUCAJ
      DO 10 I=1,54
         POM(1)=TRLN(1,1)*BLT(1,I)+TRLN(1,2)*BLT(2,I)+TRLN(1,3)*BLT(3,I)
     1         +TRLN(1,4)*BLT(4,I)+TRLN(1,5)*BLT(5,I)+TRLN(1,6)*BLT(6,I)
         POM(2)=TRLN(2,1)*BLT(1,I)+TRLN(2,2)*BLT(2,I)+TRLN(2,3)*BLT(3,I)
     1         +TRLN(2,4)*BLT(4,I)+TRLN(2,5)*BLT(5,I)+TRLN(2,6)*BLT(6,I)
         POM(3)=TRLN(3,1)*BLT(1,I)+TRLN(3,2)*BLT(2,I)+TRLN(3,3)*BLT(3,I)
     1         +TRLN(3,4)*BLT(4,I)+TRLN(3,5)*BLT(5,I)+TRLN(3,6)*BLT(6,I)
         POM(4)=TRLN(4,1)*BLT(1,I)+TRLN(4,2)*BLT(2,I)+TRLN(4,3)*BLT(3,I)
     1         +TRLN(4,4)*BLT(4,I)+TRLN(4,5)*BLT(5,I)+TRLN(4,6)*BLT(6,I)
         POM(5)=TRLN(5,1)*BLT(1,I)+TRLN(5,2)*BLT(2,I)+TRLN(5,3)*BLT(3,I)
     1         +TRLN(5,4)*BLT(4,I)+TRLN(5,5)*BLT(5,I)+TRLN(5,6)*BLT(6,I)
         POM(6)=TRLN(6,1)*BLT(1,I)+TRLN(6,2)*BLT(2,I)+TRLN(6,3)*BLT(3,I)
     1         +TRLN(6,4)*BLT(4,I)+TRLN(6,5)*BLT(5,I)+TRLN(6,6)*BLT(6,I)
          BLT(1,I)=POM(1)
          BLT(2,I)=POM(2)
          BLT(3,I)=POM(3)
          BLT(4,I)=POM(4)
          BLT(5,I)=POM(5)
          BLT(6,I)=POM(6)
   10 CONTINUE
      RETURN
      END
