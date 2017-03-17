C==========================================================================
C==========================================================================
C  SUBROUTINE NEWID1
C             NEWIDF
C             ZAGLAF
C             STAGPF
C             WRITEF
C             READF
C             CLEARF
C             CLEARR
C             DETVSR
C             DVSR1
C             BALLON
C             ZBALL
C             FNTERP
C             NEWID
C             NEWIDN
C              TUBE
C              RTUBE
C              RTUBEC
C              RTUBE1
C              RCAROT
CC             IDENTI
C              OUTFLU
C              OUTPAK
C              OUTPAF
C==========================================================================
C==========================================================================
      SUBROUTINE NEWID1(ID,JEDN,NPT,PENALT,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /PENALL/ PENALT,PRESS			 
C      COMMON /TIPELM/ NETIP

C
CE Subroutine NEWID1 is used for definition numeration of equations 
CE at ID matrix
C
      DIMENSION ID(5,*),IS(3),IE(3),IS3(2),IE3(2)

      DATA IS/1,4,3/,IE/4,4,3/
      DATA IS3/1,4/,IE3/3,4/
 
       JEDN=0
C      IPOCU=0
C      IPOCV=0
C      IPOCP=0
C      IPOCT=0
 

      IF (NETIP.EQ.2) THEN

      IF (PENALT.GT.0.D0) THEN
       DO 10 I=1,NPT
   10  ID(3,I)=1
       ENDIF 
      NBROJ=0
C      DO 90 II=1,3
      DO 90 II=1,1
      DO 81 I=1,NPT
      DO 80 J=IS(II),IE(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
C      IF (J.EQ.1) IPOCU=NBROJ
C      IF (J.EQ.2) IPOCV=NBROJ
C      IF (J.EQ.3) IPOCP=NBROJ
C      IF (J.EQ.4) IPOCT=NBROJ
  80  CONTINUE
  81  CONTINUE
  90  CONTINUE
      ELSEIF (NETIP.EQ.3) THEN
      IF (PENALT.GT.0.D0) THEN
       DO 110 I=1,NPT
  110  ID(4,I)=1
       ENDIF 
      NBROJ=0
      DO 190 II=1,2
      DO 181 I=1,NPT
      DO 180 J=IS3(II),IE3(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
 180  CONTINUE
 181  CONTINUE
 190  CONTINUE
      ENDIF

      
C       CALL IWRR(ID,4*NPT,' ID3')
      JEDN=NBROJ
      END
C==========================================================================
C==========================================================================
      SUBROUTINE NEWIDF(ID,JEDN,NPT,PENALT,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /PENALL/ PENALT,PRESS			 
C      COMMON /TIPELM/ NETIP

C
CE Subroutine NEWIDF is used for definition numeration of equations 
CE at ID matrix
C

      DIMENSION ID(5,*),IS(3),IE(3),IS3(2),IE3(2)

      DATA IS/1,5,4/,IE/3,5,4/
      DATA IS3/1,4/,IE3/3,4/
 
       JEDN=0
C      IPOCU=0
C      IPOCV=0
C      IPOCP=0
C      IPOCT=0
 

      IF (NETIP.EQ.2) THEN

      IF (PENALT.GT.0.D0) THEN
       DO 10 I=1,NPT
   10  ID(4,I)=1
       ENDIF 
      NBROJ=0
      DO 90 II=1,3
C      DO 90 II=1,2
      DO 81 I=1,NPT
      DO 80 J=IS(II),IE(II)
C      IF (ID(J,I).LE.0) THEN
      IF (ID(J,I).GT.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
C      IF (J.EQ.1) IPOCU=NBROJ
C      IF (J.EQ.2) IPOCV=NBROJ
C      IF (J.EQ.3) IPOCP=NBROJ
C      IF (J.EQ.4) IPOCT=NBROJ
  80  CONTINUE
  81  CONTINUE
  90  CONTINUE
      ELSEIF (NETIP.EQ.3) THEN
      IF (PENALT.GT.0.D0) THEN
       DO 110 I=1,NPT
  110  ID(4,I)=1
       ENDIF 
      NBROJ=0
      DO 190 II=1,2
      DO 181 I=1,NPT
      DO 180 J=IS3(II),IE3(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
 180  CONTINUE
 181  CONTINUE
 190  CONTINUE
      ENDIF

      
C       CALL IWRR(ID,4*NPT,' ID3')
      JEDN=NBROJ
      END
C==========================================================================
C=======================================================================
      SUBROUTINE ZAGLAF(ISRPS,IIZLAZ)

C      COMMON /SRPSKI/ ISRPS
C      COMMON /CDEBUG/ IDEBUG
C      COMMON /ULAZNI/ IULAZ,IIZLAZ

C
CE Subroutine ZAGLAF is used for printing head at output file *.LST
C

C      IF(IDEBUG.GT.0) PRINT *, ' ZAGLAV'
C
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,2050)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2051)
      WRITE(IIZLAZ,6100)
      WRITE(IIZLAZ,2103)
      WRITE(IIZLAZ,2115)
      RETURN
C----------------------------------------------------------------------
 2050 FORMAT('1'/////////33X,'P R O G R A M'/
     1 10X,'       FOR FINITE ELEMENT ANALYSIS OF FLUID FLOW'/ 
     1 26X,'     WITH HEAT TRANSFER')
 2051 FORMAT('1'/////////33X,'P R O G R A M'/
     1 10X,'ZA STACIONARNU I NELINEARNU ANALIZU STRUJANJA NESTISLJIVOG'/
     1 19X,'VISKOZNOG FLUIDA SA PROVODJENJEM TOPLOTE'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 6100 FORMAT(/////////////
     1 11X,' PPPPPPPPP         AAAA        KK      KK     FFFFFFFFFT'/
     2 11X,' PPPPPPPPPP       AAAAAA       KK     KK      FFFFFFFFFT'/
     3 11X,' PP      PP      AA    AA      KK    KK       FF'/
     4 11X,' PP      PP     AA      AA     KK   KK        FF'/
     5 11X,' PPPPPPPPPP     AA      AA     KK  KK         FFFFFFFF'/
     6 11X,' PPPPPPPPP      AA      AA     KKKKK          FFFFFFFF'/
     7 11X,' PP             AAAAAAAAAA     KKKKKKK        FF'/
     8 11X,' PP             AAAAAAAAAA     KK    KK       FF'/
     1 11X,' PP             AA      AA     KK     KK      FF'/
     2 11X,' PP             AA      AA     KK      KK     FF')
 2103 FORMAT(//////3X,'---------------------------  VERSION  6.00  ',
     1                '---------------------------')
 2115 FORMAT(//
     1 30X,'LABORATORY FOR ENGINEERING SOFTWARE'/
     1 30X,'FACULTY OF ENGINEERING'/
     1 30X,'UL. S. JANJICA 6'/
     1 30X,'34000 KRAGUJEVAC'/
     1 30X,'SERBIA'/
     1 30X,'E-mail: zile@kg.ac.rs'/)
       END
C==========================================================================
C=======================================================================
      SUBROUTINE STAGPF(GNODE,NASLOV,VREME,KOR,NDT,NP,II,IND,NET,NEL,
     1PRITP,CCORD,CORD,NDIM,NETIP,PENALT,IDPRIT,ISRPS,SPSIL,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine STAGPF is used for printing velocities, pressure, temperature
CE in output file *.UNV for graphic postprocessing
C

C      COMMON /TIPEL/ NP2DMX
C      COMMON /TIPELM/ NETIP
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM
C      COMMON /SRPSKI/ ISRPS
C      COMMON /REPEA2/ NEL2D
      COMMON /CDEBUG/ IDEBUG
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET1
C      COMMON /EL2K09/NMAT2D,MAT2D,NE2D,NGAU2X,NGAU2Y,NP2DMX,
C     1IPR2DC,IPTG2,N2DGU,MEL2D,IQE2,NBR2,NCV2,LM2(36),IDEL2(5),NBR2F,
C     2NOD5(5),NND5,IELX2,NLM2,NTHIC,THIC
C
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT DISPLACEMENTS, VELOSITIES AND ACCELERATIONS
CE.      IN UNIVERSAL FILE
CS.   PROGRAM
CS.      ZA STAMPANJE POMERANJA, BRZINA I UBRZANJA U UNIVERZALNI FILE
C .
C ......................................................................
C
      DIMENSION NEL(NDIM,*)
      DIMENSION PRITP(IDPRIT,*),CCORD(3,*),CORD(3,*)
      CHARACTER*250 NASLOV
C     COMMON A(17000)
C      DIMENSION RTH(*),ID(5,*)
      DIMENSION GNODE(2,5,*)
      DIMENSION FSP(6)
      DIMENSION PRIT(20),PRIT1(8),SPSIL(NETIP,*),PRES(3,*)
      DIMENSION IORD(20)
      DATA IORD/1,9,2,10,3,11,4,12,17,18,19,20,5,13,6,14,7,15,8,16/


C     WRITE(*,*)'IDEBUG= ',IDEBUG                     
      IF(IDEBUG.GT.0) PRINT *, ' STAGP1'
CE    STRUCTURAL ANALYSIS = 1
CS    STRUKTURNA ANALIZA = 1
      IMOTY=1
CE    STEADY STATE = 1; TRANSIENT = 4
CS    STACIONARAN = 1; NESTACIONARAN = 4
      NNDIM=4
      IF (NDIM.EQ.9) NNDIM=8
      IANTY=1
      IFAT1=1
      IFAT2=1					    
      FATY8=0.0D0
      IF(NDT.GT.1) THEN
         IANTY=4
         IFAT1=2
         FATY8=VREME
      ENDIF
C     6 DOF = 3
      IDACH=3
CE    DISPLACEMENTS = 8
CS    POMERANJA = 8
      IF(IND.EQ.0) ISDTY=8
CE    VELOCITIES = 11
CS    BRZINE = 11
      IF(IND.EQ.1) ISDTY=11
CE    ACCELERATIONS = 12
CS    UBRZANJA = 12
      IF(IND.EQ.2) ISDTY=12
CE    SINGLE PRECISION = 2; DOUBLE PRECISION = 4
CS    PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4
      IDATY=3
CE    NUMBER DATA = 6
CS    BROJ PODATAKA = 6
      NDVPN=2
      IND1=-1
      ITYP=55
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2006)
      IF(ISRPS.EQ.1)
     1WRITE(II,6006)
      ENDIF
      IF(IND.EQ.2) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2007)
      IF(ISRPS.EQ.1)
     1WRITE(II,6007)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
      WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NETIP
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) VREME
      DO 10 I=1,NP
            WRITE(II,5000) I
         DO 20 J=1,NETIP
            FSP(J) = 0.0D0
C            IF(ID(J,I).EQ.0) GO TO 20
C            K = ID(J,I)
C               FSP(J)=RTH(K)
               FSP(J)=GNODE(2,J,I)
   20    CONTINUE
         WRITE(II,5200) (FSP(J),J=1,NETIP)
   10 CONTINUE
      WRITE(II,5100) IND1
      NNII=4
      IF (NDIM.EQ.9) NNII=8
C
C     ZA ISPISIVANJE PRITISAKA
C
      WRITE(II,5100) IND1
      WRITE(II,5100) 57
      WRITE(II,5003) NASLOV
      WRITE(II,3005)
      WRITE(II,2000) KOR
      WRITE(II,3006) 1,1,4,20,2,6
      WRITE(II,3006) 1,1,KOR
      WRITE(II,5202) VREME

      IF (PENALT.GT.1.D0) THEN
       DO 30 I=1,NET
        WRITE(II,3006) I,1,NDIM,6
        IF(NDIM.EQ.9) THEN 
         PRIT(1)=PRITP(1,I)
         PRIT(2)=PRITP(3,I)
         PRIT(3)=PRITP(4,I)
         PRIT(4)=PRITP(2,I)
         PRIT1(1)=PRIT(1)
         PRIT1(2)=(PRIT(1)+PRIT(2))/2.
         PRIT1(3)=PRIT(2)
         PRIT1(4)=(PRIT(2)+PRIT(3))/2.
         PRIT1(5)=PRIT(3)
         PRIT1(6)=(PRIT(3)+PRIT(4))/2.
         PRIT1(7)=PRIT(4)
         PRIT1(8)=(PRIT(1)+PRIT(4))/2.
          DO 28 K=1,8
  28       WRITE(II,5200) PRIT1(K),0.,0.,0.,0.,0.
        ELSE
          DO 29 J=1,NDIM
   29       WRITE(II,5200) PRITP(1,I),0.,0.,0.,0.,0.
        ENDIF
   30  CONTINUE     
      ELSE IF (NETIP.EQ.2) THEN
       DO 40 I=1,NET
        WRITE(II,3006) I,1,NNDIM,6
         DO 31 J=1,NNDIM
          PRIT(J)=0.0D0
   31    CONTINUE
         DO 32 J=1,NNDIM
C      IF(ID(3,NEL(J,I)).EQ.0) GO TO 32
C      PRIT(J)=RTH(ID(3,NEL(J,I)))    
          PRIT(J)=GNODE(2,4,NEL(J,I))
   32    CONTINUE
        PRIT1(1)=PRIT(1)
        PRIT1(2)=(PRIT(1)+PRIT(2))/2.
        PRIT1(3)=PRIT(2)
        PRIT1(4)=(PRIT(2)+PRIT(3))/2.
        PRIT1(5)=PRIT(3)
        PRIT1(6)=(PRIT(3)+PRIT(4))/2.
        PRIT1(7)=PRIT(4)
        PRIT1(8)=(PRIT(1)+PRIT(4))/2.
         NBROJ=4
         IF (NDIM.EQ.9) NBROJ=8
        DO 38 K=1,NBROJ
         IF (NDIM.EQ.9) THEN 
            WRITE(II,5200) PRIT1(K),0.,0.,0.,0.,0.
         ELSE 
           WRITE(II,5200) PRIT(K),0.,0.,0.,0.,0.
         END IF
   38    CONTINUE 
   40  CONTINUE     
C-----------------------------------------------------------------------
C 3D ELEMENT 21-NODES
	ELSE IF (NETIP.EQ.3.AND.NDIM.EQ.21) THEN 
         DO I=1,NET
          WRITE(II,3006) I,1,20,6
           DO K=1,8
            PRIT(K)=GNODE(2,4,NEL(K,I))
           ENDDO
            PRIT(9)=0.5D0*(PRIT(1)+PRIT(2))
            PRIT(10)=0.5D0*(PRIT(2)+PRIT(3))
            PRIT(11)=0.5D0*(PRIT(3)+PRIT(4))
            PRIT(12)=0.5D0*(PRIT(4)+PRIT(1))
            PRIT(13)=0.5D0*(PRIT(5)+PRIT(6))
            PRIT(14)=0.5D0*(PRIT(6)+PRIT(7))
            PRIT(15)=0.5D0*(PRIT(7)+PRIT(8))
            PRIT(16)=0.5D0*(PRIT(8)+PRIT(5))
            PRIT(17)=0.5D0*(PRIT(1)+PRIT(5))
            PRIT(18)=0.5D0*(PRIT(2)+PRIT(6))
            PRIT(19)=0.5D0*(PRIT(3)+PRIT(7))
            PRIT(20)=0.5D0*(PRIT(4)+PRIT(8))
           DO K=1,20
            WRITE(II,5200) PRIT(IORD(K)),0.,0.,0.,0.,0.
           ENDDO
         ENDDO
C-----------------------------------------------------------------------
      ENDIF
      
      
      WRITE(II,5100) IND1
C
C  ZA ISPISIVANJE TEMPERATURA
C
      WRITE(II,5100) IND1
      WRITE(II,5100) 55
      WRITE(II,5003) NASLOV
      WRITE(II,7006)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,5,2,1
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) VREME
      DO 110 I=1,NP
            WRITE(II,5000) I
            FSP(4) = 0.0D0
C            IF(ID(4,I).EQ.0) GO TO 120
C            K = ID(4,I)
C               FSP(4)=RTH(K)
               FSP(4)=GNODE(2,5,I)
  120      WRITE(II,5202) FSP(4)
  110 CONTINUE
      WRITE(II,5100) IND1


C
C  FOR PRINTING SHEAR STRESS
C
C NUMBER '12' IS ANALOGY NODAL ACCELERATIONS

      WRITE(II,5100) IND1
      WRITE(II,5100) 55
      WRITE(II,5003) NASLOV
      WRITE(II,7007)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,12,2,3
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) VREME
      DO  I=1,NP
            WRITE(II,5000) I
            WRITE(II,5200) (PRES(J,I),J=1,3)
	ENDDO
      WRITE(II,5100) IND1


C
C  FOR PRINTING PRESSURE AT NODES
C
C NUMBER '13' IS FOR PRESSURE AT NODES

      WRITE(II,5100) IND1
      WRITE(II,5100) 55
      WRITE(II,5003) NASLOV
      WRITE(II,7008)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,13,2,1
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) VREME
      DO  I=1,NP
            WRITE(II,5000) I
            WRITE(II,5200) GNODE(2,4,I)
	ENDDO
      WRITE(II,5100) IND1




C ZA ISPISIVANJE POMERANJA SOLIDA
C
      WRITE(II,5100) IND1
      WRITE(II,5100) ITYP
      WRITE(II,5003) NASLOV
      IF(IND.EQ.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(II,2005)
      IF(ISRPS.EQ.1)
     1WRITE(II,6005)
      ENDIF
      IF(ISRPS.EQ.0)
     1WRITE(II,2000) KOR
      IF(ISRPS.EQ.1)
     1WRITE(II,6000) KOR
C     WRITE(II,5000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      WRITE(II,5000) IMOTY,IANTY,IDACH,8,IDATY,6
      IF(NDT.EQ.1) WRITE(II,5000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(II,5000) IFAT1,IFAT2,KOR,KOR
      WRITE(II,5200) FATY8
      DO 130 I=1,NP
            WRITE(II,5000) I
      DO IBROJ=1,6
       FSP(IBROJ)=0.D0
      ENDDO 
C SAMO PRIVREMENO ZA PRIKAZIVANJE PRITISAKA NA MESTU Z-POMERANJA    
C        FSP(3)=PRES(I)
C       IWALL=0      
C       DO JJ=1,NETIP
C        IF(ID(JJ,I).EQ.0) IWALL=IWALL+1
C       ENDDO
C       IF (IWALL.NE.NETIP) GOTO 150
         DO J=1,NETIP
              FSP(J)=CCORD(J,I)-CORD(J,I)
C           FSP(J)=SPSIL(J,I)
        ENDDO
C150    WRITE(II,5200) (FSP(J),J=1,6)
          WRITE(II,5200) (FSP(J),J=1,6)
  130 CONTINUE
      WRITE(II,5100) IND1

C



      RETURN
C
 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
 5202 FORMAT(1PE13.2)
C-----------------------------------------------------------------------
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 3005 FORMAT('PRITISCI U FLUIDU')
 3006 FORMAT(6I10)
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
 2001 FORMAT(' DATE'/
     1       ' EMPTY'/
     1       ' LOAD CASE',I10)
C-----------------------------------------------------------------------
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOCITIES AND ANGLE VELOCITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
 7006 FORMAT(' NODAL TEMPERATURE')
 7007 FORMAT(' SHEAR STRESSES')
 7008 FORMAT(' PRESSURES AT NODES')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C===========================================================================
      SUBROUTINE WRITEF(NASLOV,A,IDATF,DATF,IDATA,IVECT,II,IR,LMAX)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG

      DIMENSION A(*),IDATF(*),DATF(*)
      REAL A
      CHARACTER*250 NASLOV
      CHARACTER *6 STR(42),STRR(11)

C
CE Subroutine WRITEF is used for printing all necessary data to files
CE ASCII and binary files are used
C
			   
      DATA STR/'LVREME','LNTABF','LID   ','LCORD ','LNEL  ','LNZAD ',
     &'LZADVR','LTABF ','LITFMA','LNGPSI','LMAXA ','LMHT  ','LTT1  ',
     &'LTT10 ','LCCORD','LSPAR1','LAKc  ','LIKc  ','LROW  ','LINDEL',
     &'LMAX  ','NTOT  ','MAXIT ','NETIP ','IULAZ ','IIZLAZ','INDAMI',
     &'NUMZAD','NPT   ','NDIM  ','MAXSIL','NEQF  ','NWKF  ','NET   ',
     &'NPER  ','NTABFT','NDES  ','IDPRIT','IFORM ','NSTAC ','INDAX ',
     &'KKORAK'/

      DATA STRR/'GUSM  ','CC    ','AKT   ','EPSTR ','AMI   ','BETA  ',
     &'TETAO ','FB2   ','FB3   ','PENALT','VVREME'/
      

      OPEN (IDATA,FILE='DATA.FLU',STATUS='UNKNOWN',FORM='FORMATTED',
     & ACCESS='SEQUENTIAL')
      OPEN (IVECT,FILE='VECTOR.FLU',STATUS='UNKNOWN',FORM='UNFORMATTED',
     & ACCESS='SEQUENTIAL')

C==========================================================================
C WRITING HEADING OF PROBLEM
      WRITE(IDATA,100) NASLOV
C==========================================================================
C WRITING INTEGER DATA
      DO I=1,II
       WRITE(IDATA,101) STR(I),'=',IDATF(I)
      ENDDO
C==========================================================================
C WRITING REAL DATA
      DO I=1,IR
       WRITE(IDATA,200) STRR(I),'=',DATF(I)
      ENDDO
C==========================================================================
C WRITING VECTOR 
	CALL WRITED(A(1),LMAX,IVECT)
C==========================================================================

	 REWIND(IDATA)
	 REWIND(IVECT)

 100   FORMAT(A80)
 101   FORMAT(A6,A1,I10)
 200   FORMAT(A6,A1,D13.5)
  


      
            
C      CLOSE (IDATA,STATUS='KEEP')
C      CLOSE (IVECT,STATUS='KEEP')



      END
C=======================================================================
C===========================================================================
      SUBROUTINE READF(NASLOV,A,IDATF,DATF,IDATA,IVECT,II,IR,LMAX)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG

      DIMENSION A(*),IDATF(*),DATF(*)
      REAL A
      CHARACTER*250 NASLOV
      CHARACTER *7 STR
      
C
CE Subroutine READF is used for reading all necessary data from files
CE ASCII and binary files are used
C


C==========================================================================
C READING HEADING OF PROBLEM
      READ (IDATA,100) NASLOV
C==========================================================================
C READING INTEGER DATA
      DO I=1,II
       READ(IDATA,101) STR,IDATF(I)
      ENDDO
C==========================================================================
C READING REAL DATA
      DO I=1,IR
       READ(IDATA,200) STR,DATF(I)
      ENDDO
C==========================================================================
C READING VECTOR 
	CALL READD(A(1),LMAX,IVECT)
C==========================================================================



 100   FORMAT(A80)
 101   FORMAT(A7,I10)
 200   FORMAT(A7,D13.5)
  


      
            
      CLOSE (IDATA,STATUS='KEEP')
      CLOSE (IVECT,STATUS='KEEP')



      END
C=======================================================================
C===========================================================================
      SUBROUTINE CLEARF(NASLOV,A,IDATF,DATF,II,IR,LMAX)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C
CE Subroutine CLEARF is used for clearing all data in the working vectors
C
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG

      DIMENSION A(*),IDATF(*),DATF(*)
      REAL A
      CHARACTER*250 NASLOV


C==========================================================================
C CLEARING HEADING OF PROBLEM
      NASLOV='                                                         '
C==========================================================================
C CLEARING INTEGER DATA
      DO I=1,II
       IDATF(I)=0
      ENDDO
C==========================================================================
C CLEARING REAL DATA
      DO I=1,IR
       DATF(I)=0.D0
      ENDDO
C==========================================================================
C CLEARING VECTOR 
	CALL CLEAR(A(1),LMAX)
C==========================================================================

      END
C=======================================================================
C=======================================================================
      SUBROUTINE CLEARR(A,N)
C
CE Subroutine CLEARR is used for clearing a floating-point array
C


C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR -  A)
CE.         A(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA -  A)
CS.         A(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
      REAL A
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEAR'
      DO 10 I=1,N
   10 A(I)=0.0
      RETURN
      END
C=======================================================================
      SUBROUTINE DETVSR(GNODE,NS,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,5,*)
C
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE
C
CE Subroutine DETVSR is used for some example
C

      K=0
      VSR=0.D0

      DO I=NS,NE
        V=GNODE(2,2,I)
        V0=GNODE(1,2,I)
        IF (DABS(V).GT.1.D-10) K=K+1
        VSR=VSR+0.5D0*(V+V0)
      ENDDO
       IF (K.EQ.0) RETURN
       VSR=VSR/K
       

      END
C=======================================================================
C=======================================================================
      SUBROUTINE DVSR1(GNODE,NODE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,5,*)
C
      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE

C
CE Subroutine DVSR1 is used for some example
C

        VX=GNODE(2,1,NODE)
        VY=GNODE(2,2,NODE)
        VZ=GNODE(2,3,NODE)
        VX0=GNODE(1,1,NODE)
        VY0=GNODE(1,2,NODE)
        VZ0=GNODE(1,3,NODE)
        V=DSQRT(VX**2+VY**2+VZ**2)
        V0=DSQRT(VX0**2+VY0**2+VZ0**2)
C        VSR=0.25D0*(V+V0)
        VSR=0.25D0*(VY+VY0)
C        VSR=0.25D0*(VZ+VZ0)
       

      END
C=======================================================================
      SUBROUTINE BALLON(ID,TT1,TT10,CCORD,FK1,TIME,VVREME,KKORAK,
     &ITER,NPT,IIZLAZ,PRES,ZADVRE,NUMZAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TT1(*),TT10(*),CCORD(3,*),PRES(*),ZADVRE(*)
      DIMENSION ID(5,*)

      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE

C
CE Subroutine BALLON is used for some example
C

      IF(ITER.GT.0) GOTO 100

      VOL0=8.18D-8
      P0=122.5D0

       
      IF (KKORAK.EQ.1.AND.ITER.EQ.0) THEN
       VOLT=VOL0
C       PBALL=P0
       VSR=0.D0
       PTUBE=0.D0
C       GOTO 90
      ENDIF
  

      PI=4.D0*DATAN(1.D0)
      R=250.D-6
      T=5.D0
      W=2.D0*PI/T


C      NQ1=ID(2,1351)
C      V=TT1(NQ1)/2.D0
C      V0=TT10(NQ1)/2.D0
C      V1=TT1(NQ1)
C      VSR=(V+V0)/2.D0
       
       R=CCORD(1,1359)
	 A=R**2*PI
C      VOL=VOLT+VSR*A*TIME
C      PTUBE=(VOL-VOL0)*1.5D2*98*1.D6
C      PBOX=P0*(DCOS(W*VVREME)-1.D0)
C      PBALL=PBOX-PTUBE

       VOL=VOLT+VSR*A*TIME
       PBALL=-(VOL-VOL0)*1.5D2*98.D0*1.D6
       PBOX=P0*(DCOS(W*VVREME)-1.D0)
       PTUBE=PBOX-PBALL
      
      VOLT=VOL

 90   WRITE(IIZLAZ,*)' TIME      PBALL     VOLT          QSR
     &PBOX        PTUBE,  VSR,  A'
      WRITE(IIZLAZ,200) VVREME,PBALL/98.D0,VOLT*1.D6,(VSR*A)*1.D6,
     &PBOX/98.D0,PTUBE/98.D0,VSR*1.D2,A*1.D4
C     WRITE(IIZLAZ,*)'Y coord (cm)   Pressure (cmH2O)'
C     DO I=1,NPT,18
C      Y=CORD(2,I)
C      JJ=ID(3,I)
C      PP=0.D0
C      IF(JJ.NE.0) PP=TT1(JJ)
C      WRITE(IIZLAZ,300) Y*1.D2,PP/98.D0
C     ENDDO
C 100  FK1=PBALL
 100  FK1=PTUBE
C      DO I=1,NUMZAD
C       ZADVRE(I)=FK1   
C      ENDDO           
      RETURN

200   FORMAT (F10.3,7D13.5)
300   FORMAT (2D13.5)
      END
C=======================================================================
C=======================================================================
      FUNCTION DDIST(N1,N2,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)
C
      X=CORD(1,N1)-CORD(1,N2)
      Y=CORD(2,N1)-CORD(2,N2)
      Z=CORD(3,N1)-CORD(3,N2)
      DDIST=DSQRT(X**2+Y**2+Z**2)
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ZBALL(GNODE,CCORD,TIME,VVREME,KKORAK,ITER,NPT,IIZLAZ,
     &CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,5,*),CCORD(3,*),CORD(3,*)
      DIMENSION NODE(6),RMIN(6),DIST(6)

      COMMON /AKIRA/ VOLT,PBALL,VSR,PBOX,PTUBE
      COMMON /AKIRA1/ VOLT1,PBALL1,VSR1,PTUBE1
      COMMON /POSIT/ PART(6,3)

C
CE Subroutine ZBALL is used for some example
C

      IF(ITER.GT.0) RETURN

      VOL0=8.18D-2
      P0=1225.D0
      PIN=-3.D0*980.D0
C
      DO I=1,NPT
       DO J=1,3
        CCORD(J,I)=CCORD(J,I)+GNODE(2,J,I)*TIME 
       ENDDO
      ENDDO
       
      IF (KKORAK.EQ.1.AND.ITER.EQ.0) THEN
       VOLT=VOL0
       VOLT1=VOL0
       VSR=0.D0
       VSR1=0.D0
       PTUBE=0.D0
       PTUBE1=0.D0
      ENDIF
  

      PI=4.D0*DATAN(1.D0)
      T=5.D0
      W=2.D0*PI/T

      R=DDIST(588,586,CORD)
      WRITE(IIZLAZ,*)'R1= ',R
      R1=DDIST(1125,1137,CORD)
      WRITE(IIZLAZ,*)'R2= ',R1

       A=R**2*PI
       A1=R1**2*PI

C
C      VOL=(PTUBE/PENALT)*TIME+VOLT
C      PBALL=-(VOL-VOL0)*1.5D2*980.D0

C       PBOX=P0*(DCOS(W*VVREME)-1.D0)
       PBOX=PIN+P0*(DCOS(W*VVREME)-1.D0)
       
C       PBALL=PBOX-PTUBE
C       PBALL1=PBOX-PTUBE1

       PBALL=0.5D0*(PBALL+PBOX-PTUBE)
       PBALL1=0.5D0*(PBALL1+PBOX-PTUBE1)

       VOL=VOL0-PBALL/(1.5D2*980.D0)
       VOL1=VOL0-PBALL1/(1.5D2*980.D0)

      VSR=0.5D0*(VSR+(VOL-VOLT)/(A*TIME))
      VSR1=0.5D0*(VSR1+(VOL1-VOLT1)/(A1*TIME))

C      VSR=2.D0*((VOL-VOLT)/(A*TIME))-VSR
C      VSR1=2.D0*((VOL1-VOLT1)/(A*TIME))-VSR1

      
      VOLT=VOL
      VOLT1=VOL1

 90   WRITE(IIZLAZ,*)'    TIME      PBALL1     VOLT1        QSR1
     &PBOX      PTUBE1      VSR1     A1'
      WRITE(IIZLAZ,200) VVREME,PBALL/980.D0,VOLT,(VSR*A),
     &PBOX/980.D0,PTUBE/980.D0,VSR,A
      WRITE(IIZLAZ,*)'    TIME      PBALL2     VOLT2        QSR2
     &PBOX      PTUBE2      VSR2     A2'
      WRITE(IIZLAZ,200) VVREME,PBALL1/980.D0,VOLT1,(VSR1*A1),
     &PBOX/980.D0,PTUBE1/980.D0,VSR1,A1
C
      IF (KKORAK.EQ.1) THEN
           DO I=1,3
            PART(1,I)=CORD(I,5)
            PART(2,I)=CORD(I,6)
            PART(3,I)=CORD(I,11)
            PART(4,I)=CORD(I,305)
            PART(5,I)=CORD(I,306)
            PART(6,I)=CORD(I,311)
           ENDDO
      ELSE
        DO J=1,6
        RMIN(J)=1.D5
        DO I=1,NPT
         DIST(J)=DSQRT((CORD(1,I)-PART(J,1))**2+
     &(CORD(2,I)-PART(J,2))**2+(CORD(3,I)-PART(J,3))**2)
         IF (DIST(J).LT.RMIN(J)) THEN
                 RMIN(J)=DIST(J)
                 NODE(J)=I
         ENDIF
        ENDDO
        ENDDO
        DO I=1,3
         DO J=1,3
          IF (PART(J,2).LT.0.2) 
     &PART(J,I)=PART(J,I)+GNODE(2,I,NODE(J))*TIME
         ENDDO
         DO J=4,6
          IF (PART(J,2).LT.0.38) 
     &PART(J,I)=PART(J,I)+GNODE(2,I,NODE(J))*TIME
         ENDDO
        ENDDO
      ENDIF
      WRITE(IIZLAZ,*)'POSITION OF POINT '
      DO J=1,6
       WRITE(IIZLAZ,400) (PART(J,I),I=1,3),NODE(J)
      ENDDO

200   FORMAT (F10.3,7D13.5)
300   FORMAT (2D13.5)
400   FORMAT (3e13.5,2X,I5)
      END
C=======================================================================
C======================================================================
C======================================================================
      SUBROUTINE FNTERP(R,S,NPARAM,TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,
     &ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    Subroutine FNTERP is used for integration 2D finite element
C          
C      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
C     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /VISKOZ/ AMI,INDAMI
C      COMMON /PENALL/ PENALT,PRESS

      DIMENSION V2(9),V3(9),T(9),ZVHR(9),ZVHS(9)
      DIMENSION H(9),ZVHX(9),ZVHY(9)
      DIMENSION HP(4),X(9),Y(9)
      DIMENSION TT21(36)

      DIMENSION AJ(2,2)

      DO 10 I=1,NDIM
      V2(I)=TT21(I)
      V3(I)=TT21(I+NDIM)
      T(I)=TT21(I+2*NDIM+4)
  10  CONTINUE
      

      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
      H(9)=RR*SS
      H(8)=0.5*RP*SS-0.5*H(9)
      H(7)=0.5*RR*SM-0.5*H(9)
      H(6)=0.5*RM*SS-0.5*H(9)
      H(5)=0.5*RR*SP-0.5*H(9)
      H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
      H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
      H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
      H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE
      H(4)=0.25*RP*SM
      H(3)=0.25*RM*SM
      H(2)=0.25*RM*SP
      H(1)=0.25*RP*SP

      ENDIF


      HP(4)=0.25*RP*SM
      HP(3)=0.25*RM*SM
      HP(2)=0.25*RM*SP
      HP(1)=0.25*RP*SP


      IF (NDIM.GT.4) THEN
      ZVHR(9)=-2.D0*R*SS
      ZVHS(9)=-2.D0*S*RR
      ZVHR(8)=0.5*SS-0.5*ZVHR(9)
      ZVHS(8)=-RP*S-0.5*ZVHS(9)
      ZVHR(7)=-R*SM-0.5*ZVHR(9)
      ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
      ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
      ZVHS(6)=-RM*S-0.5*ZVHS(9)
      ZVHR(5)=-R*SP-0.5*ZVHR(9)
      ZVHS(5)=0.5*RR-0.5*ZVHS(9)
      ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
      ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
      ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
      ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
      ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
      ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
      ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)
      ELSE

      ZVHR(4)=0.25*SM
      ZVHS(4)=-0.25*RP
      ZVHR(3)=-0.25*SM
      ZVHS(3)=-0.25*RM
      ZVHR(2)=-0.25*SP
      ZVHS(2)=0.25*RM
      ZVHR(1)=0.25*SP
      ZVHS(1)=0.25*RP

      ENDIF
      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)
      
      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
C      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      WRITE(*,*)'DETERMINANTE LESS THEN ZERO!!!'
      STOP
      ENDIF

      

      DO 20 I=1,NDIM
      ZVHX(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
      ZVHY(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
  20  CONTINUE  

      IF (IUPWIN.EQ.1)
     &CALL INTEF2(X,Y,V2,V3,H,ZVHX,ZVHY,NDIM,AKT,GUSM,AMI)

      IF (NPARAM.NE.0) THEN
      VXX=0.D0
      VXY=0.D0
      VYX=0.D0
      VYY=0.D0
      PRESS=0.D0
      TX=0.D0
      TY=0.D0
          IF (DABS(R-1.).LT.1.E-5 .OR. DABS(R-(-1.)).LT.1.E-5) THEN
C          IF (R.EQ.1. .OR. R.EQ.-1.) THEN
           DETJS=DSQRT(AJ(2,1)**2+AJ(2,2)**2)
           ELX=R*AJ(2,2)/DETJS
           ELY=-R*AJ(2,1)/DETJS
           ETX=R*AJ(2,1)/DETJS
           ETY=R*AJ(2,2)/DETJS
          ELSE
           DETJS=DSQRT(AJ(1,2)**2+AJ(1,1)**2)
           ELX=-S*AJ(1,2)/DETJS
           ELY=S*AJ(1,1)/DETJS
           ETX=S*AJ(1,1)/DETJS
           ETY=S*AJ(1,2)/DETJS
          ENDIF

      DO 30 I=1,NDIM
      I1=I
      I2=I+NDIM
      I3=2*NDIM+I
      I4=2*NDIM+4+I
      VXX=VXX+ZVHX(I)*TT21(I1)
      VXY=VXY+ZVHY(I)*TT21(I1)
      VYX=VYX+ZVHX(I)*TT21(I2)
      VYY=VYY+ZVHY(I)*TT21(I2)
      TX=TX+ZVHX(I)*TT21(I4)
      TY=TY+ZVHY(I)*TT21(I4)
      IF (I.LE.4.AND. PENALT.LT.1.D0) PRESS=PRESS+HP(I)*TT21(I3)
  30  CONTINUE

      IF (PENALT.GT.1.D0) THEN
       IF(INDAX.EQ.1) THEN
        RR=0.D0
        DO I=1,4
         RR=RR+HP(I)*X(I)
        ENDDO
         HV2=DOT(H,V2,NDIM)
          PRESS=-PENALT*(VXX+VYY+HV2/RR)
       ELSE
          PRESS=-PENALT*(VXX+VYY)
       ENDIF
      ENDIF
     

       FS2=0.D0
       FS3=0.D0
       TAU=0.D0

      ENDIF

      IF (NPARAM.EQ.1) THEN
       FS2=AMI*(ELX*VXX+ELY*VXY)-ELX*PRESS
       FS3=AMI*(ELX*VYX+ELY*VYY)-ELY*PRESS
      ENDIF

      IF (NPARAM.EQ.2) THEN
C       TAU=-AMI*((VXX*ETX+VYX*ETY)*ELX+(VXY*ETX+VYY*ETY)*ELY)
C       FLUX=-AKT*(ELX*TX+ELY*TY)
       TAU=-AKT*(ELX*TX+ELY*TY)
      ENDIF
      



      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)
C      HT=DOT(H,T,NDIM)
      ZVXT=DOT(ZVHX,T,NDIM)
      ZVYT=DOT(ZVHY,T,NDIM)
      ZVXV2=DOT(ZVHX,V2,NDIM)
      ZVYV3=DOT(ZVHY,V3,NDIM)
      ZVYV2=DOT(ZVHY,V2,NDIM)
      ZVXV3=DOT(ZVHX,V3,NDIM)

      END
C======================================================================













C======================================================================
C=======================================================================
C==========================================================================

C==========================================================================
      SUBROUTINE NEWID(ID,JEDN,NPT,PENALT,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /PENALL/ PENALT,PRESS			 
C      COMMON /TIPELM/ NETIP
C
CE    Subroutine NEWID is used for determination matrix ID-numbers of equation
C
      DIMENSION ID(5,*),IS(3),IE(3),IS3(2),IE3(2)

      DATA IS/1,4,3/,IE/2,4,3/
      DATA IS3/1,4/,IE3/3,4/
 
C      IPOCU=0
C      IPOCV=0
C      IPOCP=0
C      IPOCT=0
 

      IF (NETIP.EQ.2) THEN

      IF (PENALT.GT.0.D0) THEN
       DO 10 I=1,NPT
   10  ID(3,I)=1
       ENDIF 
      NBROJ=0
      DO 90 II=1,3
C      DO 90 II=1,2
      DO 81 I=1,NPT
      DO 80 J=IS(II),IE(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
C      IF (J.EQ.1) IPOCU=NBROJ
C      IF (J.EQ.2) IPOCV=NBROJ
C      IF (J.EQ.3) IPOCP=NBROJ
C      IF (J.EQ.4) IPOCT=NBROJ
  80  CONTINUE
  81  CONTINUE
  90  CONTINUE
      ELSEIF (NETIP.EQ.3) THEN
      IF (PENALT.GT.0.D0) THEN
       DO 110 I=1,NPT
  110  ID(4,I)=1
       ENDIF 
      NBROJ=0
      DO 190 II=1,2
      DO 181 I=1,NPT
      DO 180 J=IS3(II),IE3(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
 180  CONTINUE
 181  CONTINUE
 190  CONTINUE
      ENDIF

      
C       CALL IWRR(ID,4*NPT,' ID3')
      JEDN=NBROJ
      END
C==========================================================================
C==========================================================================
      SUBROUTINE NEWIDN(ID,JEDN,NPT,PENALT,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT
C      COMMON /PENALL/ PENALT,PRESS			 
C      COMMON /TIPELM/ NETIP
C
CE    Subroutine NEWIDN is used for determination matrix ID - numbers of equation
C
      DIMENSION ID(5,*),IS(3),IE(3),IS3(2),IE3(2)

C      DATA IS/1,4,3/,IE/4,4,3/
      DATA IS/1,4,3/,IE/2,4,3/
      DATA IS3/1,4/,IE3/3,4/
 
C      IPOCU=0
C      IPOCV=0
C      IPOCP=0
C      IPOCT=0
 

      IF (NETIP.EQ.2) THEN

      IF (PENALT.GT.0.D0) THEN
       DO 10 I=1,NPT
   10  ID(4,I)=1
       ENDIF 
      NBROJ=0
C      DO 90 II=1,1
      DO 90 II=1,3
      DO 81 I=1,NPT
      DO 80 J=IS(II),IE(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
C      IF (J.EQ.1) IPOCU=NBROJ
C      IF (J.EQ.2) IPOCV=NBROJ
C      IF (J.EQ.3) IPOCP=NBROJ
C      IF (J.EQ.4) IPOCT=NBROJ
  80  CONTINUE
  81  CONTINUE
  90  CONTINUE
      ELSEIF (NETIP.EQ.3) THEN
      IF (PENALT.GT.0.D0) THEN
       DO 110 I=1,NPT
  110  ID(4,I)=1
       ENDIF 
      NBROJ=0
      DO 190 II=1,2
      DO 181 I=1,NPT
      DO 180 J=IS3(II),IE3(II)
      IF (ID(J,I).LE.0) THEN
      NBROJ=NBROJ+1
      ID(J,I)=NBROJ
      ELSE
      ID(J,I)=0
      ENDIF
 180  CONTINUE
 181  CONTINUE
 190  CONTINUE
      ENDIF

      
C       CALL IWRR(ID,4*NPT,' ID3')
      JEDN=NBROJ
      END
C==========================================================================
C==========================================================================
C     SUBROUTINE NEWID(ID,JEDN,NPT)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     DIMENSION ID(4,*)
C     NBROJ=0
C     DO 81 I=1,NPT
C     DO 80 J=1,4
C     IF (ID(J,I).LE.0) THEN
C     NBROJ=NBROJ+1
C     ID(J,I)=NBROJ
C     ELSE
C     ID(J,I)=0
C     ENDIF
C 80  CONTINUE
C 81  CONTINUE
C      CALL IWRR(ID,4*NPT,' ID3')
C     JEDN=NBROJ
C     END
C==========================================================================
C=======================================================================
C==========================================================================
      SUBROUTINE TUBE(CCORD,IZLAZF,NU,NV,LAYERS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CCORD(3,*)

       DO I=0,LAYERS
         NSTART=I*(NU+1)*(NV+1)
         CALL RTUBE(CCORD,IZLAZF,NU,NV,NSTART)
       ENDDO

       END
C==========================================================================
      SUBROUTINE RTUBE(CCORD,IZLAZF,NU,NV,NSTART)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CCORD(3,*)


         DU=1.D0/NU
         DV=1.D0/NV
         KK=NSTART
         
         DY=CCORD(3,KK+36)/5.D0
         DX=CCORD(1,KK+1)/5.D0
         DO I=2,5
          CCORD(3,I*6+KK)=DY*(I-1)
          CCORD(1,7-I+KK)=DX*(I-1)
         ENDDO


C         DO K=NSTART,NEND
C         KK=K*(NU+1)*(NV+1)
         DO I=1,3
         F00=CCORD(I,1+KK)
         F01=CCORD(I,NV+1+KK)
         F10=CCORD(I,NU*(NV+1)+1+KK)
         F11=CCORD(I,(NU+1)*(NV+1)+KK)
        DO IU=0,NU-1
         U=DU*IU
         DO IV=2,NV+1
         V=DV*(IV-1)
         NODE=IU*(NV+1)+IV+KK
         NPSI1=IV+KK
         NPSI2=NPSI1+NU*(NV+1)
         NSI1=IU*(NV+1)+1+KK
         NSI2=NSI1+NV

         PSI1=CCORD(I,NSI1)
         PSI2=CCORD(I,NSI2)
         SI1=CCORD(I,NPSI1)
         SI2=CCORD(I,NPSI2)
         P1=(1.D0-V)*PSI1+V*PSI2
         P2=(1.D0-U)*SI1+U*SI2
         PP=-(1.D0-U)*(1.D0-V)*F00-(1.D0-U)*V*F01-U*V*F11-U*(1.D0-V)*F10
         CORD=P1+P2+PP
         CCORD(I,NODE)=CORD
          ENDDO
         ENDDO
        ENDDO
C       ENDDO


      END
C==========================================================================
C==========================================================================
      SUBROUTINE RTUBEC(CCORD,IZLAZF,NU,NV,NSTART)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CCORD(3,*)


         DU=1.D0/NU
         DV=1.D0/NV


         

         KK=NSTART
         DO I=1,3
         F00=CCORD(I,1+KK)
C         F01=CCORD(I,NV+1+KK)
         F01=CCORD(I,NSTART-30)
         F10=CCORD(I,NU*(NV+1)+1+KK-1)
         F11=CCORD(I,(NU+1)*(NV+1)+KK-1)
        DO IU=0,NU-1
         U=DU*IU
         DO IV=2,NV+1
         V=DV*(IV-1)
         NODE=IU*(NV+1)+IV+KK-1
         IF (IU.EQ.0) NODE=NODE+1
         NPSI1=IV+KK
         NPSI2=NPSI1+NU*(NV+1)-1
         NSI1=IU*(NV+1)+1+KK-1
         NSI2=NSI1+NV
           IF (IV.EQ.NV+1) NPSI1=NSTART-30
           IF (IU.EQ.0) NSI2=NSTART-30
           IF (IU.EQ.0) NSI1=NSI1+1
           IF (IU.EQ.0.AND.IV.EQ.NV+1) NODE=NSTART-30

         PSI1=CCORD(I,NSI1)
         PSI2=CCORD(I,NSI2)
         SI1=CCORD(I,NPSI1)
         SI2=CCORD(I,NPSI2)
         P1=(1.D0-V)*PSI1+V*PSI2
         P2=(1.D0-U)*SI1+U*SI2
         PP=-(1.D0-U)*(1.D0-V)*F00-(1.D0-U)*V*F01-U*V*F11-U*(1.D0-V)*F10
         CORD=P1+P2+PP
         CCORD(I,NODE)=CORD
          ENDDO
         ENDDO
        ENDDO
C       ENDDO


      END
C==========================================================================
C==========================================================================
      SUBROUTINE RTUBE1(CCORD,IZLAZF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CCORD(3,*)




         NU=5
         NV=5
         DU=1.D0/NU
         DV=1.D0/NV


         DO K=0,50
         KK=K*36
         DO I=1,3
        F00=CCORD(I,1+KK)
        F01=CCORD(I,6+KK)
        F10=CCORD(I,31+KK)
        F11=CCORD(I,36+KK)
        DO IU=0,4
         U=DU*IU
         DO IV=2,6
         V=DV*(IV-1)
         NODE=IU*6+IV+KK
         NPSI1=IV+KK
         NPSI2=NPSI1+30
         NSI1=IU*6+1+KK
         NSI2=NSI1+5

         PSI1=CCORD(I,NSI1)
         PSI2=CCORD(I,NSI2)
         SI1=CCORD(I,NPSI1)
         SI2=CCORD(I,NPSI2)
         P1=(1.D0-V)*PSI1+V*PSI2
         P2=(1.D0-U)*SI1+U*SI2
         PP=-(1.D0-U)*(1.D0-V)*F00-(1.D0-U)*V*F01-U*V*F11-U*(1.D0-V)*F10
         CORD=P1+P2+PP
         CCORD(I,NODE)=CORD
          ENDDO
         ENDDO
        ENDDO
       ENDDO


      END
C==========================================================================
C==========================================================================
      SUBROUTINE RCAROT(CCORD,IZLAZF,NU,NV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION CCORD(3,*)


       NSTART=0
       DO I=0,9
         CALL RTUBE(CCORD,IZLAZF,NU,NV,NSTART)
          NSTART=NSTART+36
         CALL RTUBE(CCORD,IZLAZF,NU,NV-1,NSTART)
          NSTART=NSTART+30
       ENDDO


         CALL RTUBE(CCORD,IZLAZF,NU,NV,NSTART)
          NSTART=NSTART+36

         CALL RTUBEC(CCORD,IZLAZF,NU,NV,NSTART)
          NSTART=NSTART+35

       DO I=0,21
          CALL RTUBE(CCORD,IZLAZF,NU,NV,NSTART)
          NSTART=NSTART+36
       ENDDO


      END
C==========================================================================
C==========================================================================
C     SUBROUTINE IDENTI(NPODS,LSTART,IPAKF,KOLKF)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     include 'paka.inc'
C     
C     COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
C    1                IOPGL(6),KOSI,NDIN,ITEST
C     COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C    1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C     COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
C    1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
C    1                NWP,NWG,IDF,JPS1
C
C
C     COMMON /INTERA/ IINTER,NPTI
C     DIMENSION NPODS(JPS1,*)
C     CHARACTER*6    FIPAKI
C
C
C       IOUTS=64
C       IOUTF=63
C       IINTER=62
C
C     
C  IDVA=2  
C
C
C
C       LMAX13=NPODS(JPS1,1)-1
C       NPP=NP
C       LIDS=LSTART
C       LMAX=LIDS+NPP*6
C       CALL DELJIV(LMAX,2,INDL)
C       IF(INDL.EQ.0) LMAX=LMAX+1
C       LCORDS=LMAX
C       LMAX=LCORDS+NPP*3*IDVA
C        IF(LMAX.GT.MTOT) CALL ERROR(1)
C       CALL READDD(A(LCORDS),NPP*3,IPODS,LMAX13,LDUZI)
C       CALL IREADD(A(LIDS),NP*6,IPODS,LMAX13,LDUZI)
C  
C        NPS=NP
C        CALL OUTPAK(IOUTS,A(LCORDS),A(LIDS),NPS)
C
C
C        
C     
C
C        REWIND IPAKF
C        LPAKF=LMAX
C        LADD=LMAX-LSTART
C        CALL READD(A(LPAKF),KOLKF,IPAKF)
C        CALL INTREA(LIDATF,A(LPAKF+2))
C        CALL OUTFLU(A(LIDATF+LADD),LIDF,LCORDF,NPF)
C
C        CALL OUTPAF(IOUTF,A(LCORDF+LADD),A(LIDF+LADD),NPF)
C
C        
C       CALL DELJIV(LMAX,2,INDL)
C       IF(INDL.EQ.0) LMAX=LMAX+1
C
C        LIDENT=LMAX						 
C        CALL IDENSF(A(LIDENT),NPTI,A(LCORDS),A(LCORDF+LADD),NPS,NPF)
C       
C           FIPAKI='ZIPAKI'
C           OPEN (IINTER,FILE=FIPAKI,STATUS='UNKNOWN',
C    1                  FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
C        KOLKI=(NPTI*2)/IDVA
C        IF(KOLKI.GT.0) CALL WRITED(A(LIDENT),KOLKI,IINTER)
C
C
C     END
C==========================================================================
        SUBROUTINE OUTFLU(IFFLU,LIDF,LCORDF,NPF)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION IFFLU(*)

        LIDF=IFFLU(3)
        LCORDF=IFFLU(4)
        NPF=IFFLU(30)
        

      END
C==========================================================================
        SUBROUTINE OUTPAK(IOUT,CORD,ID,NP)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION CORD(NP,*),ID(NP,*)

        OPEN (IOUT,FILE='CHECK')
         DO I=1,NP
           WRITE(IOUT,100) I,(ID(I,J),J=1,6),(CORD(I,J),J=1,3)
         ENDDO
      
         CLOSE (IOUT)
 100   FORMAT(I5,1X,6I5,3F10.3)
      END
C==========================================================================
C==========================================================================
        SUBROUTINE OUTPAF(IOUT,CORD,ID,NP)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION CORD(3,*),ID(5,*)

        OPEN (IOUT,FILE='CHECKF')
         DO I=1,NP
           WRITE(IOUT,100) I,(ID(J,I),J=1,5),(CORD(J,I),J=1,3)
         ENDDO
      
         CLOSE (IOUT)
 100   FORMAT(I5,1X,5I5,3F10.3)
      END