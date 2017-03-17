C=======================================================================
      SUBROUTINE EXPLINTGR(LIPODS,RTDT,UBRZANJE,BRZINA,BRZINAIPO,
     1          POMERANJE,AMASA,PRIGUSENJE,DT,DTP,DTCR,VREM0,KOR,KOJPAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        FOR EXPLICIT INTEGRATION
CS.    P R O G R A M
CS.        ZA EKSPLICITNU INTEGRACIJU
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      include 'paka.inc'
      DIMENSION UBRZANJE(*),BRZINA(*),BRZINAIPO(*),
     1          POMERANJE(*),RTDT(*),AMASA(*),PRIGUSENJE(*)
C     

      WRITE(*,*) 'SA POCETKA EXPLINTGR'
      WRITE(*,*) 'kor=',KOR
      IF(KOR.EQ.1) then
        CALL UZMISILE(A(LIPODS),POMERANJE,KOJPAK)
        IF (DT>DTCR) DT=0.9*DTCR

        DO I = 1,JEDN
c           A(n)=A(0)+M^(-1)*(F(n)-C*V(n-1/2)
            UBRZANJE(I)=UBRZANJE(I)+RTDT(I)/AMASA(I)
c           V(n+1/2)=V(n-1/2)+0.5*DT*A(n)
            BRZINAIPO(I)=BRZINAIPO(I)+0.5*(DT)*UBRZANJE(I)
c           D(n+1) = D(0) + 0.5*(DTP+DT)*V(n+1/2)
            POMERANJE(I) = POMERANJE(I) + 0.5*(DT+DTP)*BRZINAIPO(I)
        ENDDO
        KOR = KOR + 1
        VREM0 = VREM0 + DT
        DTP = DT
      endif
c    
      CALL UZMISILE(A(LIPODS),POMERANJE,KOJPAK)
    
      DO I = 1,JEDN
c       A(n+1)=M^(-1)*(F(n)-C*V(n-1/2)
        UBRZANJE(I)=RTDT(I)/AMASA(I)
c       V(n+1)=V(n+1/2)+1/2*DT*A(N)
        BRZINAIPO(I) = BRZINAIPO(I) + 0.5*(DTP+DT)*UBRZANJE(I)
c       D(n+1) = D(n) + 1/2*DT*V(n+1/2)
        POMERANJE(I) = POMERANJE(I) + 0.5*(DT+DTP)*BRZINAIPO(I)
      ENDDO
c     
      WRITE(3,*) 'UBRZANJE ZA KORAK BROJ ',KOR
      CALL WRR(UBRZANJE,JEDN,'UBRZ')
      WRITE(3,*) 'BRZINA ZA KORAK BROJ ',KOR
      CALL WRR(BRZINAIPO,JEDN,'BRZI')
      WRITE(3,*) 'POMERANJA ZA KORAK BROJ ',KOR
      CALL WRR(POMERANJE,JEDN,'POME')
c
      RETURN
      END      
C=======================================================================
      SUBROUTINE UZMISILE(NPODS,POMERANJE,KOJPAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        FOR GETING SUM OF FORCES
CS.    P R O G R A M
CS.        ZA DOBIJANJE SUME SILA
C .
C ......................................................................
      include 'paka.inc'
      COMMON /GRUPER/ LIGRUP
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      DIMENSION NPODS(*),POMERANJE(*)

      CALL JEDNA1(A(LRTDT),POMERANJE,JEDN)
      CALL INTNAP(A(LIGRUP))
      CALL DESNAL(NPODS,KOJPAK)
      CALL DESSTR(NPODS)
      
      RETURN
      END