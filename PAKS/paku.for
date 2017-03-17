C============================================================================
      SUBROUTINE MIDEA7(AMAT,ISUMGR,MAT,IGRAF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AMAT(*)
C
      E=AMAT(1)
      V=AMAT(2)
      DEN=AMAT(3)
      G=AMAT(4)
      IND=-1
      ITYP=1710
      N=0
      N14=14
      J=1
      Z=0.
      WRITE(IGRAF,5000) IND
      WRITE(IGRAF,5000) ITYP
      WRITE(IGRAF,6100) ISUMGR,MAT,N,N,N,N,N,N,N14
 6100 FORMAT(
     1'=================================================================
     1==============='/
     2'MATERIAL'/
     1'=================================================================
     1==============='/
     4I10,'  MATERIAL (',I2,')'/
     5I10,'  LINE(S) OF TEXT'/
     6I10,'  MATERIAL CLASS(ES)'/
     7I10,'  MATERIAL ATTRIBUTE(S)'/
     8I10,'  MATERIAL COMPONENT(S)'/
     9I10,'  MATERIAL SPECIFICATION(S)'/
     9'-----------------------------------------------------------------
     9---------------'/
     1I10,'  MATERIAL VARIABLE(S)'/
     9'-----------------------------------------------------------------
     9---------------'/
     3I10,'  MATERIAL PROPERT(IES)'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6200) J,N,E,J,N,Z,J,N,V
 6200 FORMAT(
     1'MODULUS OF ELASTICITY'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'SHEAR MODULUS'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'POISSONS RATIO'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6300) J,N,Z,J,N,Z,J,N,DEN,J,N,Z
 6300 FORMAT(
     1'COEFFICIENT OF THERMAL EXPANSION'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'COEFFICIENT OF THERMAL EXPANSION        1/KELVIN                 
     5               '/                                
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'THERMAL CONDUCTIVITY'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     2'CONDUCTIVITY                            MN MM/MM/C/SEC           
     2               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'MASS DENSITY'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'MASS DENSITY                            KILOGRAM/MILLIMETER^3    
     8               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     1'SPECIFIC HEAT'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'SPECIFIC HEAT                           MN MM/KG/C               
     5               '/                                
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6400) J,N,Z,J,N,Z
 6400 FORMAT(
     8'STRUCTURAL ELEMENT DAMPING COEFFICIENT'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'THERMAL EXPANSION REFERENCE TEMPERATURE'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'TEMPERATURE                             CELSIUS                  
     9               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6500) J,N,Z,J,N,Z,J,N,Z
 6500 FORMAT(
     1'ALLOWABLE STRESS IN TENSION X-DIR'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'ALLOWABLE STRESS IN COMPRESSION X-DIR'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'ALLOWABLE STRESS IN TENSION Y-DIR'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6600) J,N,Z,J,N,Z
 6600 FORMAT(
     1'ALLOWABLE STRESS IN COMPRESSION Y-DIR'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'ALLOWABLE IN-PLANE SHEAR STRESS'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6700) J,J,J,J,J,J,J,J,J,J,J,J,J,J,J,J
 6700 FORMAT(
     1'DEFAULT MATERIAL PROPERT(IES):'/
     2'MODULUS OF ELASTICITY                   VERSION : ',I10/
     3'SHEAR MODULUS                           VERSION : ',I10/
     4'POISSONS RATIO                          VERSION : ',I10/
     5'COEFFICIENT OF THERMAL EXPANSION        VERSION : ',I10/
     6'THERMAL CONDUCTIVITY                    VERSION : ',I10/
     7'MASS DENSITY                            VERSION : ',I10/
     8'SPECIFIC HEAT                           VERSION : ',I10/
     9'STRUCTURAL ELEMENT DAMPING COEFFICIENT  VERSION : ',I10/
     9'THERMAL EXPANSION REFERENCE TEMPERATURE VERSION : ',I10/
     1'ALLOWABLE STRESS IN TENSION X-DIR       VERSION : ',I10/
     2'ALLOWABLE STRESS IN COMPRESSION X-DIR   VERSION : ',I10/
     3'ALLOWABLE STRESS IN TENSION Y-DIR       VERSION : ',I10/
     4'ALLOWABLE STRESS IN COMPRESSION Y-DIR   VERSION : ',I10/
     5'ALLOWABLE IN-PLANE SHEAR STRESS         VERSION : ',I10/
     9'-----------------------------------------------------------------
     9---------------'/
     7I10,'  REFERENCE ENTITIES'/
     8I10,'  MATERIAL TYPES'/
     9'FEM                                     ISOTROPIC MATERIALS'/
     1'=================================================================
     1===============')
      WRITE(IGRAF,5000) IND
 5000 FORMAT(I6)     
      RETURN
      END
C============================================================================
      SUBROUTINE MIDEAS(AMAT,ISUMGR,MAT,IGRAF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AMAT(*)
C
      E=AMAT(1)
      V=AMAT(2)
      DEN=AMAT(3)
      G=AMAT(4)
      IND=-1
      ITYP=1710
      N=0
      N14=27
      J=1
      Z=0.
      WRITE(IGRAF,5000) IND
      WRITE(IGRAF,5000) ITYP
      WRITE(IGRAF,6100) ISUMGR,MAT,N,N,N,N,N,N,N14
 6100 FORMAT(
     1'=================================================================
     1==============='/
     2'MATERIAL'/
     1'=================================================================
     1==============='/
     4I10,'  MATERIAL (',I2,')'/
     5I10,'  LINE(S) OF TEXT'/
     6I10,'  MATERIAL CLASS(ES)'/
     7I10,'  MATERIAL ATTRIBUTE(S)'/
     8I10,'  MATERIAL COMPONENT(S)'/
     9I10,'  MATERIAL SPECIFICATION(S)'/
     9'-----------------------------------------------------------------
     9---------------'/
     1I10,'  MATERIAL VARIABLE(S)'/
     9'-----------------------------------------------------------------
     9---------------'/
     3I10,'  MATERIAL PROPERT(IES)'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6200) J,N,E,J,N,Z,J,N,V
 6200 FORMAT(
     1'MODULUS OF ELASTICITY'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'SHEAR MODULUS'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'POISSONS RATIO'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6300) J,N,Z,J,N,Z,J,N,Z
 6300 FORMAT(
     1'COEFFICIENT OF THERMAL EXPANSION'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'COEFFICIENT OF THERMAL EXPANSION        1/KELVIN                 
     5               '/                                
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'THERMAL CONDUCTIVITY'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     2'CONDUCTIVITY                            MN MM/MM/C/SEC           
     2               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     1'SPECIFIC HEAT'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'SPECIFIC HEAT                           MN MM/KG/C               
     5               '/                                
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6350) J,N,Z,J,N,DEN
 6350 FORMAT(
     5'HEAT GENERATION RATE'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'HEAT FLUX PER UNIT VOLUME               MN MM/MM^3/SEC           
     9               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'MASS DENSITY'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'MASS DENSITY                            KILOGRAM/MILLIMETER^3    
     8               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6400) J,N,Z,J,N,Z
 6400 FORMAT(
     8'STRUCTURAL ELEMENT DAMPING COEFFICIENT'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'THERMAL EXPANSION REFERENCE TEMPERATURE'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'TEMPERATURE                             CELSIUS                  
     9               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6500) J,N,Z,J,N,Z,J,N,Z
 6500 FORMAT(
     1'ALLOWABLE STRESS IN TENSION X-DIR'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'ALLOWABLE STRESS IN COMPRESSION X-DIR'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     5'ALLOWABLE STRESS IN TENSION Y-DIR'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     9'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6600) J,N,Z,J,N,Z
 6600 FORMAT(
     1'ALLOWABLE STRESS IN COMPRESSION Y-DIR'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     6'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------'/
     8'ALLOWABLE IN-PLANE SHEAR STRESS'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     3'CONSTANT'/1PE25.16/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6610) J,N,J,N,J,N
 6610 FORMAT(
     1'YIELD STRESS'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'PRESSURE                                MILLINEWTON/MILLIMETER^2 
     5               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     8'CONVECTIVE FILM COEFFICIENT'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     2'CONVECTION COEFFICIENT                  MN MM/MM^2/C/SEC         
     2               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     5'THERMAL CAPACITY PER UNIT AREA'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'THERMAL CAPACITY PER UNIT AREA          MN MM/MM^2/C             
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6620) J,N,J,N,J,N
 6620 FORMAT(
     1'SURFACE HEAT FLUX RATE'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'HEAT FLUX PER UNIT AREA                 MN MM/MM^2/SEC           
     5               '/   
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     8'VISCOSITY'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     2'VISCOSITY                               KG/MM/SEC                
     2               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     5'COEFFICIENT OF FRICTION'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6630) J,N,J,N,J,N
 6630 FORMAT(
     1'AREA FACTOR'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     8'EMISSIVITY'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     5'ABSORPTIVITY'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6640) J,N,J,N,J,N
 6640 FORMAT(
     1'HEAT FLUX RATE'/
     2I10,'  VERSION NUMBER'/
     3I10,'  LINE(S) OF TEXT'/
     4'DIMENSIONS AND UNITS:'/
     5'HEAT FLUX PER UNIT LENGTH               MN MM/MM/SEC             
     5               '/   
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     8'INTERACTION TERM FOR TSAI-WU'/
     9I10,'  VERSION NUMBER'/
     9I10,'  LINE(S) OF TEXT'/
     1'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------'/
     5'SWELLING COEFFICIENT'/
     6I10,'  VERSION NUMBER'/
     7I10,'  LINE(S) OF TEXT'/
     8'DIMENSIONS AND UNITS:'/
     9'DIMENSIONLESS                           UNITLESS                 
     9               '/
     6'NULL_PROPERTY'/
     9'-----------------------------------------------------------------
     9---------------')
      WRITE(IGRAF,6700) J,J,J,J,J,J,J,J,J,J,J,J,J,J,J
 6700 FORMAT(
     1'DEFAULT MATERIAL PROPERT(IES):'/
     2'MODULUS OF ELASTICITY                   VERSION : ',I10/
     3'SHEAR MODULUS                           VERSION : ',I10/
     4'POISSONS RATIO                          VERSION : ',I10/
     5'COEFFICIENT OF THERMAL EXPANSION        VERSION : ',I10/
     6'THERMAL CONDUCTIVITY                    VERSION : ',I10/
     8'SPECIFIC HEAT                           VERSION : ',I10/
     8'HEAT GENERATION RATE                    VERSION : ',I10/
     7'MASS DENSITY                            VERSION : ',I10/
     9'STRUCTURAL ELEMENT DAMPING COEFFICIENT  VERSION : ',I10/
     9'THERMAL EXPANSION REFERENCE TEMPERATURE VERSION : ',I10/
     1'ALLOWABLE STRESS IN TENSION X-DIR       VERSION : ',I10/
     2'ALLOWABLE STRESS IN COMPRESSION X-DIR   VERSION : ',I10/
     3'ALLOWABLE STRESS IN TENSION Y-DIR       VERSION : ',I10/
     4'ALLOWABLE STRESS IN COMPRESSION Y-DIR   VERSION : ',I10/
     5'ALLOWABLE IN-PLANE SHEAR STRESS         VERSION : ',I10)
      WRITE(IGRAF,6750) J,J,J,J,J,J,J,J,J,J,J,J,J,J
 6750 FORMAT(
     1'YIELD STRESS                            VERSION : ',I10/
     1'CONVECTIVE FILM COEFFICIENT             VERSION : ',I10/
     1'THERMAL CAPACITY PER UNIT AREA          VERSION : ',I10/
     1'SURFACE HEAT FLUX RATE                  VERSION : ',I10/
     1'VISCOSITY                               VERSION : ',I10/
     1'COEFFICIENT OF FRICTION                 VERSION : ',I10/
     1'AREA FACTOR                             VERSION : ',I10/
     1'EMISSIVITY                              VERSION : ',I10/
     1'ABSORPTIVITY                            VERSION : ',I10/
     1'HEAT FLUX RATE                          VERSION : ',I10/
     1'INTERACTION TERM FOR TSAI-WU            VERSION : ',I10/
     1'SWELLING COEFFICIENT                    VERSION : ',I10/
     9'-----------------------------------------------------------------
     9---------------'/
     7I10,'  REFERENCE ENTITIES'/
     8I10,'  MATERIAL TYPES'/
     9'FEM                                     ISOTROPIC MATERIALS'/
     1'=================================================================
     1===============')
      WRITE(IGRAF,5000) IND
 5000 FORMAT(I6)     
      RETURN
      END
C======================================================================
      SUBROUTINE STAM31(NMAT,AU,ISNA,MCVEL,ICVEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRINT NUMBER OF MATERIALS FOR 3/D ELEMENTS IN UNIVERSAL FILE
CS.   P R O G R A M
CS.      ZA STAMPANJE BROJA MATERIJALA 3/D ELEMENATA U UNIVERZALNI FILE
C .
C ......................................................................
C
      CHARACTER*250 NASLOV
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /NASLOV/ NASLOV
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
C
      DIMENSION MCVEL(*),ISNA(*),SR(6),NMAT(*)
      DIMENSION AU(*)
      REAL AU
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STAG31'
      NNCVE=NCVE
      IF(NCVE.LT.20) NCVE=8
      IF(NCVE.EQ.21) NCVE=20
C     STRUKTURNA ANALIZA = 1
      IMOTY=1
C     STACIONARAN = 1; NESTACIONARAN = 4   
      IANTY=1
      IFAT1=1
      IFAT2=1
      FATY8=0.0D0
      IF(NDT.GT.1) THEN 
        IANTY=4
        IFAT1=2
        FATY8=VREME
      ENDIF
C     SIMETRICAN TENZOR = 4
      IDACH=4
C     NAPONI = 2
      ISDTY=2
C     PRECIZNOST JEDNOSTRUKA = 2; DVOSTRUKA = 4 
      IDATY=2
C     BROJ PODATAKA = 6
      NDVPN=6
      IEXP=1
      NNODS=NCVE
      NNODS=-NNODS
      NVPN=6
      IND1=-1
      INA1=57
      INDSN=0
      KOR=1
      VREME=0.
      DO 20 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 20
C
         ISN=ISNA(NLM)
         IF(ISN.EQ.1) GO TO 20
         IF(INDSN.EQ.0) THEN
      WRITE(IGRAF,5000) IND1
      WRITE(IGRAF,5000) INA1
      WRITE(IGRAF,5003) NASLOV
      KOR=1
      IF(ISRPS.EQ.0)
     1WRITE(IGRAF,2004) KOR
      IF(ISRPS.EQ.1)
     1WRITE(IGRAF,6004) KOR
      WRITE(IGRAF,1000) IMOTY,IANTY,IDACH,ISDTY,IDATY,NDVPN
      IF(NDT.EQ.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR
      IF(NDT.GT.1) WRITE(IGRAF,1000) IFAT1,IFAT2,KOR,KOR
      WRITE(IGRAF,5001) FATY8
         ENDIF
         INDSN=1
         NUM=NLM+ISUMEL
         IF(ICVEL.EQ.1) NUM=MCVEL(NLM)
         WRITE(IGRAF,1000) NUM,IEXP,NNODS,NVPN
         CALL CLEAR(SR,6)
           SR(1)=NMAT(NLM)
           WRITE(IGRAF,5001) (SR(I),I=1,6) 
   20 CONTINUE
      IF(INDSN.EQ.1) WRITE(IGRAF,5000) IND1
      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
 1000 FORMAT(8I10)
 5000 FORMAT(I6)
 5001 FORMAT(6(1PE13.5))
 5003 FORMAT(A80)
C-----------------------------------------------------------------------
 2004 FORMAT('NAPONI 3/D ELEMENATA'/
     1       'DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
C-----------------------------------------------------------------------
 6004 FORMAT('STRESS OF 3/D ELEMENTS'/
     1       'DATE'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
