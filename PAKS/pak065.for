      subroutine dmumps1(irn,jcn,a,rhs,nz,n,kkk)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) mumps_par
c      common /mumps/ mumps_par
      INTEGER IERR,I,j,n,nz,irn(nz),jcn(nz),kkk
      double precision a,rhs,amax
      dimension a(nz),rhs(n)
c      amax=a(1)
c      do i = 1,nz
c          if (irn(i).eq.jcn(i).and.amax.le.a(i)) amax=a(i)
c      enddo
c      write(3,*) 'najveci clan dijagonale matrice krutosti amax=',amax
      do i = 1,nz
        if (irn(i).eq.jcn(i).and.a(i).lt.1.d-15) then
            a(i)=1.
        endif
      enddo
      if (kkk.eq.1) then
         if (mumps_par%JOB.eq.3) then
            mumps_par%JOB = -2
            CALL DMUMPS(mumps_par)
            IF ( mumps_par%MYID .eq. 0 ) then
                  DEALLOCATE( mumps_par%IRN )
                  DEALLOCATE( mumps_par%JCN )
                  DEALLOCATE( mumps_par%A   )
                  DEALLOCATE( mumps_par%RHS )
            endif
         endif
         mumps_par%COMM = MPI_COMM_WORLD
         mumps_par%SYM = 1
         mumps_par%PAR = 1
         mumps_par%JOB = -1
         CALL DMUMPS(mumps_par)
C  Define problem on the host (processor 0)
         IF ( mumps_par%MYID .eq. 0 ) THEN
            mumps_par%N=n
            mumps_par%NZ=nz
c           READ(5,*) mumps_par%N
c           READ(5,*) mumps_par%NZ
            ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
            ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
            ALLOCATE( mumps_par%A ( mumps_par%NZ ) )
            ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
            DO I = 1, mumps_par%NZ
               mumps_par%IRN(I)=irn(i)
               mumps_par%JCN(I)=jcn(i)
               mumps_par%A(I)=a(i)
c              WRITE(*,*) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
            END DO
         END IF
C  Call package for solution
         mumps_par%JOB = 4
         CALL DMUMPS(mumps_par)
         
         IF ( mumps_par%MYID .eq. 0 ) THEN
            DO I = 1, mumps_par%NZ
               irn(i)=mumps_par%IRN(I)
               jcn(i)=mumps_par%JCN(I)
               a(i)=mumps_par%A(I)
            END DO
         END IF
      else
         IF ( mumps_par%MYID .eq. 0 ) THEN
            DO I = 1, mumps_par%N
               mumps_par%RHS(I)=rhs(i)
c              READ(5,*) mumps_par%RHS(I)
            END DO
         END IF
         mumps_par%JOB = 3
         CALL DMUMPS(mumps_par)
      endif
C  Solution has been assembled on the host
      IF ( mumps_par%MYID .eq. 0 ) THEN
        DO I = 1, mumps_par%N
           rhs(i)=mumps_par%RHS(I)
c          READ(5,*) mumps_par%RHS(I)
        END DO
      END IF
      RETURN
      END


! =========================================================================
