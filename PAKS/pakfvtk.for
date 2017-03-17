C     ******************************************************************
C     *                                                                *
C     *      SUBROUTINES TO PRINT:                                     *
C     *                               * MESH IN VTK FILE               *
C     *                               * RESULTS IN VTK FILE            *
C     *                                                                *
C     *      DATE CREATED:     22. 02. 2011 BY BLAGOJEVIC MILAN        *
C     *      LAST UPDATE:      17. 03. 2011 BY BLAGOJEVIC MILAN        *
C     *                                                                *
C     *      LIST OF SUBROUTINES:                                      *
C     *           /not used/ OTVVTK  - Open *.vtk file                 *
C     *                      OTVTK   - Open *.vtk file for every step  *
C     *                      VTKHD   - Write header of *.vtk file      *
C     *                      VTKPHD  - Write header for POINT DATA     *
C     *                      VTKEHD  - Write header for CELL DATA      *
C     *                      VTKMSH  - Write MESH in *.vtk file        *
C     *                      VTKRES  - Write RESULTS in  *.vtk file    *
C     *                      VTKSTR  - Write STREAMLINES in *.vtk file *
C     *                                                                *
C     ******************************************************************

C=======================================================================
      SUBROUTINE OTVVTK(IME,PAKVTK,IDUZIN,ISRPS)
CS OTVARANJE DATOTEKE ZA VTK GRAFIKU
CE OPEN FILE FOR VTK GRAPHIC
      LOGICAL OLDNEW
      CHARACTER*24 PAKVTK
      CHARACTER*3  STAT
      CHARACTER*20 IME

CS IZLAZ ZA GRAFIKU
CE OUTPUT FOR GRAPHIC
C
   15 CONTINUE
      IF(ISRPS.EQ.0)
     *WRITE(*,*)' UNETI IME IZLAZNE DATOTEKE ZA VTK GRAFIKU /"'
     1,PAKVTK(1:IDUZIN),'"'
      IF(ISRPS.EQ.1)
     *WRITE(*,*)' ENTER NAME OF OUTPUT FILE FOR VTK GRAPHICS /"' 
     1,PAKVTK(1:IDUZIN),'"'
C   15 WRITE(*,820)
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = PAKVTK
  910 FORMAT (A)
C
   20 STAT='NEW'
      INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(OLDNEW) STAT='OLD'
      IF(STAT.EQ.'NEW') THEN
      OPEN (89,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ELSE
      OPEN (89,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ENDIF
C
      IND=0
      IF(STAT.EQ.'OLD') CALL BRISF (IME,89,IND)
      IF(IND.EQ.1)GO TO 20
      IF(IND.EQ.2)GO TO 15
      CALL VTKHD(89)
      RETURN
      END
C======================================================================

C=======================================================================
      SUBROUTINE OTVTK(IMEF,PAKVTK,IDUZIN,ISRPS,NASLOV,KORAK,IDISK)
      CHARACTER*250 NASLOV
      CHARACTER*24 PAKVTK
      CHARACTER*3  STAT
      CHARACTER*20 IMEF,IME     
      CHARACTER*4 SKORAK
C
C Subroutine to write header for the *.vtk file
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IMEF     - FILE NAME
C          PAKVTK   - USER DEFINED STRING (CURENTLY NOT USED)
C          IDUZIF   - LENGTH OF FILE NAME
C          ISRPS    - LANGUAGE VRSION
C          NASLOV   - TITLE OF PROBLEM
C          KORAK    - TIME SPEP
C          IDISK    - NUMBER OF FILE
C     
      WRITE (SKORAK,'(I0)') KORAK
      IDISK = 89
      IME = IMEF(1:IDUZIN-4)
C      OPEN (IDISK,FILE=TRIM(IME)//'_'//TRIM(SKORAK)//'.vtk',
C     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      OPEN (IDISK,FILE='VTK_'//TRIM(SKORAK)//'.vtk',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      CALL VTKHD(IDISK,NASLOV,SKORAK)
      RETURN
      END
C======================================================================

C**********************************************************************
      SUBROUTINE VTKHD(IDISK,NASLOV,SKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*4 SKORAK
      CHARACTER*250 NASLOV
C
C Subroutine to write header for the *.vtk file
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NASLOV   - TITLE OF PROBLEM
C          SKORAK   - CURRENT STEP
C
      WRITE(IDISK,2500)
      WRITE(IDISK,2501) TRIM(NASLOV)//' - Step: '//TRIM(SKORAK)
      WRITE(IDISK,2502)

 2500	FORMAT('# vtk DataFile Version 2.0')
 2501	FORMAT(A)
 2502	FORMAT('ASCII')
      RETURN
	END
C**********************************************************************

C**********************************************************************
      SUBROUTINE VTKPHD(IDISK,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 SNPT
C
C Subroutine to write header for the RESULTS per NODES
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NPT      - TOTAL NUMBER OF NODES
C
      WRITE (SNPT,'(I0)') NPT
      WRITE(IDISK,2503) 'POINT_DATA '//TRIM(SNPT)

 2503	FORMAT(A)
      RETURN
	END
C**********************************************************************

C**********************************************************************
      SUBROUTINE VTKEHD(IDISK,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 SNET
C
C Subroutine to write header for the RESULTS per ELEMENTS
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NET      - TOTAL NUMBER OF ELEMENTS
C
      WRITE (SNET,'(I0)') NET
      WRITE(IDISK,2504) 'CELL_DATA '//TRIM(SNET)

 2504	FORMAT(A)
      RETURN
	END
C**********************************************************************

C**********************************************************************
      SUBROUTINE VTKMSH(IDISK,NNODES,CORD,NETIP,NEL,NCVE,NET)
C
C Subroutine to write MESH to *.vtk file
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NNODES   - TOTAL NUMBER OF NODES
C          CORD     - NODES COORDINATES
C          NETIP    - TYPE OF FINITE ELEMENT (2-2D or 3-3D)
C          NEL      - NODES PER ELEMENTS
C          NCVE     - TOTAL NUMBER OF NODES PER ELEMENT
C          NET      - TOTAL NUMBER OF ELEMENTS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*),NEL(NCVE,*)      
      CHARACTER*10 SNNODES,SNET,SNET1,SCELLT*2
      WRITE (SNNODES,'(I0)') NNODES
      WRITE (SNET,'(I0)') NET
      WRITE (SNET1,'(I0)') NET*(NCVE+1)
      
      NNCVE = NCVE
      WRITE(IDISK,2600)
      WRITE(IDISK,2601) 'POINTS '//TRIM(SNNODES)//' double'
	DO NODE=1,NNODES
		   WRITE(IDISK,2602) (CORD(I,NODE),I=1,3)
      ENDDO

      IF(NETIP.EQ.2.AND.NCVE.GT.8) NCVE = 8 
      IF(NETIP.EQ.3.AND.NCVE.GT.20) NCVE = 20  
       
      WRITE(IDISK,2603) 'CELLS '//TRIM(SNET)//' '//TRIM(SNET1)
      DO NE=1,NET    

C     2D ELEMENT
	IF(NETIP.EQ.2.AND.NCVE.EQ.4) 
     1WRITE(IDISK,2604) NCVE,((NEL(NODE,NE)-1),NODE=1,4)
      IF(NETIP.EQ.2.AND.NCVE.EQ.8) 
     1WRITE(IDISK,2604) NCVE,((NEL(NODE,NE)-1),NODE=1,8)

C     3D ELEMENT
	IF(NETIP.EQ.3.AND.NCVE.EQ.8) 
     1WRITE(IDISK,2604) NCVE,((NEL(NODE,NE)-1),NODE=1,8)
C	2                       ((NEL(NODE,NE)-1),NODE=1,4)
      IF(NETIP.EQ.3.AND.NCVE.EQ.20)
     1WRITE(IDISK,2604) NCVE,((NEL(NODE,NE)-1),NODE=1,20)
C (** IZ NEKOG RAZLOGA NIJE POTREBNO ISPRATITI NUMERACIJU **)
C    	2                      ((NEL(NODE,NE)-1),NODE=1,4),
C     3                      ((NEL(NODE,NE)-1),NODE=13,20),
C     4                      ((NEL(NODE,NE)-1),NODE=9,12)

      ENDDO
      
      WRITE(IDISK,2605)'CELL_TYPES '//SNET
      IF(NETIP.EQ.2.AND.NCVE.EQ.4) SCELLT = '9'
      IF(NETIP.EQ.2.AND.NCVE.EQ.8) SCELLT = '23' 
      IF(NETIP.EQ.3.AND.NCVE.EQ.8) SCELLT = '12'
      IF(NETIP.EQ.3.AND.NCVE.EQ.20) SCELLT = '25'
           
      DO NE=1,NET
           WRITE(IDISK,2605) TRIM(SCELLT)
      ENDDO

      NCVE = NNCVE
 2600 FORMAT('DATASET UNSTRUCTURED_GRID')
 2601 FORMAT(A)
 2602 FORMAT(E13.6,2X,E13.6,2X,E13.6)
 2603 FORMAT(A)
 2604 FORMAT(I2,20(2X,I8))
 2605 FORMAT(A)
      RETURN
	END
C**********************************************************************

C**********************************************************************
      SUBROUTINE VTKRES(IDISK,NNODES,GNODE,IVS,IRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,5,*)
      CHARACTER VECTOR
C
C Subroutine to write SCALAR or VECTOR NODAL RESULTS IN *.VTK FILE
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 17. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NNODES   - TOTAL NUMBER OF NODES
C          GNODE    - MATRIX VITH RESULTS IN CURRENT STEP
C          IVS      - VECTOR or SKALAR
C                     IVS.EQ.1 - SKALAR
C                     IVS.EQ.3 - VEKTOR
C          IRES     - TYPE OF RESULT
C                     IRES.EQ.1 - VELOCITY IN X,Y,Z DIRECTION
C                     IRES.EQ.4 - PRESSURES AT NODES
C                     IRES.EQ.5 - TEMPERATURE
C
	IF (IVS.EQ.1) THEN 
	   IF (IRES.EQ.4) WRITE(IDISK,2714)
	   IF (IRES.EQ.5) WRITE(IDISK,2715)
	   WRITE(IDISK,2702)
	ENDIF
	IF (IVS.EQ.3) THEN
	   IF (IRES.EQ.1) WRITE(IDISK,2731)
	ENDIF
	
	DO NODID=1,NNODES 
      IF(IVS.EQ.1)
     1WRITE(IDISK,2703) GNODE(2,IRES,NODID)
      IF(IVS.EQ.3) THEN
         X = GNODE(2,IRES,NODID)
         Y = GNODE(2,IRES+1,NODID)
         Z = GNODE(2,IRES+2,NODID)
         WRITE(IDISK,2705) X,Y,Z                        
      ENDIF
	ENDDO


 2700 FORMAT(A)
 2714 FORMAT('SCALARS pressures_at_nodes double')
 2715 FORMAT('SCALARS temperature double')
 2731 FORMAT('VECTORS velocity double')
 2702 FORMAT('LOOKUP_TABLE default')
 2703 FORMAT(E13.6)
 2704 FORMAT(A)
 2705 FORMAT(E13.6,2X,E13.6,2X,E13.6)
      RETURN
	END
C**********************************************************************

C**********************************************************************
      SUBROUTINE VTKSHS(IDISK,NNODES,PRES,IVS,IRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION PRES(3,*)
      CHARACTER VECTOR
C
C Subroutine to write SHEAR STRESS RESULTS IN *.VTK FILE
C
C     Created: 31. 03. 2011 by Blagojevic Milan
C     Last update: 31. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IDISK    - NUMBER OF FILE
C          NNODES   - TOTAL NUMBER OF NODES
C          PRES     - MATRIX VITH RESULTS IN CURRENT STEP
C          IVS      - VECTOR or SKALAR
C                     IVS.EQ.1 - SKALAR
C                     IVS.EQ.3 - VEKTOR
C          IRES     - TYPE OF RESULT
C                     IRES.EQ.1 - SHEAR STRESS IN X,Y,Z DIRECTION
C
	IF (IVS.EQ.1) THEN 
	   IF (IRES.EQ.4) WRITE(IDISK,2814)
	   IF (IRES.EQ.5) WRITE(IDISK,2815)
	   WRITE(IDISK,2802)
	ENDIF
	IF (IVS.EQ.3) THEN
	   IF (IRES.EQ.1) WRITE(IDISK,2831)
	ENDIF
	
	DO NODID=1,NNODES 
      IF(IVS.EQ.1)
     1WRITE(IDISK,2803) PRES(IRES,NODID)
      IF(IVS.EQ.3) THEN
         X = PRES(IRES,NODID)
         Y = PRES(IRES+1,NODID)
         Z = PRES(IRES+2,NODID)
         WRITE(IDISK,2805) X,Y,Z                        
      ENDIF
	ENDDO


 2800 FORMAT(A)
 2814 FORMAT('SCALARS pressures_at_nodes double')
 2815 FORMAT('SCALARS temperature double')
 2831 FORMAT('VECTORS shear_stress double')
 2802 FORMAT('LOOKUP_TABLE default')
 2803 FORMAT(E13.6)
 2804 FORMAT(A)
 2805 FORMAT(E13.6,2X,E13.6,2X,E13.6)
      RETURN
	END
C**********************************************************************

C======================================================================
      SUBROUTINE VTKSTM(IDISK,ID,TT1,VREME,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ID(3,*)
      DIMENSION TT1(*)
C
C  FOR PRINTING PAKF STREAMLINES
C
      DO I=1,NPT
       JJ=ID(1,I)
       PSI=0.D0
       IF (JJ.NE.0) PSI=TT1(ID(1,I))
       WRITE(IDISK,5202) PSI
      ENDDO

 3006 FORMAT(6I10)
 5100 FORMAT(I6)
 5202 FORMAT(1PE13.2)
 5003 FORMAT(A80)
 7006 FORMAT(' NODAL STREAMLINES')
 2001 FORMAT(' DATE'/
     1       ' EMPTY'/
     1       ' LOAD CASE',I10)


      END
C======================================================================

C**********************************************************************
      SUBROUTINE VTKSTP(IMEF,PAKUSR,IDUZIF,ISRPS,NASLOV,GNODE,PRES,
     1KKORAK,NPT,CORD,NETIP,NEL,NDIM,NET)
C
C Subroutine that calls all the subroutines to write *.vtk file
C
C     Created: 14. 03. 2011 by Blagojevic Milan
C     Last update: 31. 03. 2011 by Blagojevic Milan
C
C     Variables:
C          IMEF     - FILE NAME
C          PAKUSR   - USER DEFINED STRING (CURENTLY NOT USED)
C          IDUZIF   - LENGTH OF FILE NAME
C          ISRPS    - LANGUAGE VRSION
C          NASLOV   - TITLE OF PROBLEM
C          GNODE    - MATRIX VIT RESULTS IN CURRENT STEP
C          KKORAK   - TIME SPEP
C          NPT      - TOTAL NUMBER OF NODES
C          CORD     - NODES COORDINATES
C          NETIP    - TYPE OF FINITE ELEMENT (2-2D or 3-3D)
C          NEL      - NODES PER ELEMENTS
C          NDIM     - TOTAL NUMBER OF NODES PER ELEMENT
C          NET      - TOTAL NUMBER OF ELEMENTS
C          IDISK    - NUMBER OF FILE
C     
      CHARACTER*250 NASLOV
      CHARACTER*1 IMEF*20
      CHARACTER *24 PAKUSR
       
      CALL OTVTK(IMEF,PAKUSR,IDUZIF,ISRPS,NASLOV,KKORAK,IDISK)
      CALL VTKMSH(IDISK,NPT,CORD,NETIP,NEL,NDIM,NET)
      CALL VTKPHD(IDISK,NPT)
      CALL VTKRES(IDISK,NPT,GNODE,1,4)
      CALL VTKRES(IDISK,NPT,GNODE,1,5)
      CALL VTKRES(IDISK,NPT,GNODE,3,1)
      CALL VTKSHS(IDISK,NPT,PRES,3,1)
      CLOSE (IDISK, STATUS="KEEP")
           
      RETURN
      END
C**********************************************************************
