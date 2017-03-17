      subroutine ICM_CG(x,K,MaxK,ColK,M,MaxM,ColM,r,z,p,neq,nK,nM,
     +                  epsilon,tol,niter)
c     Solve a symmetric system of linear equations by the PCG iterative
c     method with the ICM preconditioner. r already initialized: r=b.
      integer*4 i,neq,nK,nM,niter,
     +          MaxK(neq+1),ColK(nK),MaxM(neq+1),ColM(nM)
      real*8    x(neq),K(nK),r(neq),z(neq),p(neq),
     +          epsilon,tol,alfa,beta,gamma,gamma0,delta,delta0
      real*4    M(nM)
C      DIMENSION NZAD(3,*),ID(1,*),ZADVRE(*)
c ----------------------------------------------------------------------
c     PCG-initialization.
      niter=1
      call PREC(neq,nM,z,M,MaxM,ColM,r)
      gamma0=0.0D0
      delta0=0.0D0
      do i=1,neq
         x(i)=0.0D0
         p(i)=z(i)
         gamma0=gamma0+r(i)*r(i)
         delta0=delta0+r(i)*z(i)
      enddo

C PRESCRIBED VALUES:
C      CALL PRESCR(X,NZAD,ZADVRE,ID,ITER)

      gamma0=DSQRT(gamma0)
c ----------------------------------------------------------------------
c     PCG-iterations.
 1000 continue
      call COMVEC(z,K,MaxK,ColK,p,neq,nK)
      alfa=0.0D0
      do 100 i=1,neq
         alfa=alfa+p(i)*z(i)
  100 continue
      alfa=delta0/alfa
      gamma=0.0D0
      do 200 i=1,neq
         x(i)=x(i)+alfa*p(i)
         r(i)=r(i)-alfa*z(i)
         gamma=gamma+r(i)*r(i)
  200 continue
C PRESCRIBED VALUES:
C      CALL PRESCR(X,NZAD,ZADVRE,ID,ITER)
      tol=DSQRT(gamma)/gamma0
      write(*,'(1x,a,i8,6x,a,d20.14)') 'iter=',niter,'tol=',tol
      if (tol.le.epsilon) goto 2000
      niter=niter+1
      call PREC(neq,nM,z,M,MaxM,ColM,r)
      delta=0.0D0
      do i=1,neq
         delta=delta+r(i)*z(i)
      enddo
      beta=delta/delta0
      delta0=delta
      do 300 i=1,neq
         p(i)=z(i)+beta*p(i)
  300 continue
      goto 1000
c ----------------------------------------------------------------------
 2000 return
      end
c **********************************************************************

      subroutine ICM(neq,nK,mM,nM,nMrej,psi,
     +               K,MaxK,ColK,M,MaxM,ColM,Sky,Mi,DiagM)
c     Creates the ICM (Incomplete Cholesky by Magnitude) preconditioner.
      integer*4 neq,nK,mM,nM,nMrej,MaxK(neq+1),ColK(nK),MaxM(neq+1),
     +          ColM(mM),Sky(neq),i,ii,j,ii_end,j_end,colj,mColM,isky
      real*8    K(nK),Mi(neq),DiagM(neq),psi,x,mii,mjj,msky,s
      real*4    M(mM)
c ----------------------------------------------------------------------
      nM=0
      nMrej=0
      mColM=0
      do 10 i=1,neq
	   if(dabs(K(MaxK(i))).lt.1.d-15) K(MaxK(i))=1.
         DiagM(i)=K(MaxK(i))
         Sky(i)=i
         Mi(i)=0.0D0
   10 continue
      do 1000 i=1,neq
         nM=nM+1
         MaxM(i)=nM
         ColM(nM)=i
         mii=DiagM(i)
         j_end=MaxK(i+1)-1
         colj=ColK(j_end)
         if (colj.gt.mColM) mColM=colj
         do 100 j=MaxK(i)+1,j_end
            colj=ColK(j)
            Mi(colj)=K(j)
            if (Sky(colj).eq.colj) Sky(colj)=i
  100    continue
         ii_end=i-1
         do 200 ii=Sky(i),ii_end
            isky=Sky(ii)
            if (ColM(isky).eq.i) then
               msky=M(isky)*M(MaxM(ii))
               mii=mii-msky*M(isky)
               Sky(ii)=isky+1
               j_end=MaxM(ii+1)-1
               do 250 j=Sky(ii),j_end
                  colj=ColM(j)
                  Mi(colj)=Mi(colj)-msky*M(j)
  250          continue
            endif
  200    continue
         do 300 j=i+1,mColM
            if (Sky(j).eq.j) goto 300
            x=Mi(j)
            if (x.eq.0.0D0) goto 300
            Mi(j)=0.0D0
            mjj=DiagM(j)
            if (x*x.lt.psi*mii*mjj) then
               x=DABS(x)
               s=DSQRT(mii/mjj)
               mii=mii+x*s
               DiagM(j)=mjj+x/s
               nMrej=nMrej+1
            else
               nM=nM+1
               if (nM.gt.mM) then
                  write(*,*) 'Memory for preconditioner exceeded!'
                 stop 
               endif
               M(nM)=x
               ColM(nM)=j
            endif
  300    continue
         if (mii.le.0.0D0) then
            write(*,*) 'Diagonal term <=0 in ICM-precond. - row',i
C            stop
         endif
         M(MaxM(i))=mii
         do 400 j=MaxM(i)+1,nM
            M(j)=M(j)/mii
  400    continue
         Sky(i)=MaxM(i)+1
         if (Sky(i).gt.nM) Sky(i)=MaxM(i)
 1000 continue
c      write(*,*) 'mm,nm',mm,nm
c      write(3,*) 'mm,nm',mm,nm
      MaxM(neq+1)=nM+1
c ----------------------------------------------------------------------
      return
      end
c **********************************************************************

      subroutine PREC(neq,nM,z,M,MaxM,ColM,r)
      integer*4 i,neq,nM,MaxM(neq),ColM(nM)
      real*8    z(neq),r(neq)
      real*4    M(nM)
c ----------------------------------------------------------------------
      do 10 i=1,neq
         z(i)=r(i)
   10 continue
      call FORS1R(neq,nM,z,M,MaxM,ColM)
      do 100 i=1,neq
         z(i)=z(i)/M(MaxM(i))
  100 continue
      call BAKS1R(neq,nM,z,M,MaxM,ColM)
c ----------------------------------------------------------------------
      return
      end
c **********************************************************************

      subroutine COMVEC(c,A,MaxA,ColA,b,n,nA)
      integer*4 n,nA,MaxA(n),ColA(nA),i,j,colj,n_end,j_end
      real*8    c(n),A(nA),b(n),ci,aj,bi
c ----------------------------------------------------------------------
      do 10 i=1,n
   10 c(i)=A(MaxA(i))*b(i)
      n_end=n-1
      do 100 i=1,n_end
         bi=b(i)
         ci=c(i)
         j_end=MaxA(i+1)-1
         do 110 j=MaxA(i)+1,j_end
            aj=A(j)
            colj=ColA(j)
            c(colj)=c(colj)+aj*bi
            ci=ci+aj*b(colj)
  110    continue
         c(i)=ci
  100 continue
c ----------------------------------------------------------------------
      return
      end
c **********************************************************************

      subroutine FORS1R(n,nA,x,A,MaxA,ColA)
      integer*4 n,nA,MaxA(n),ColA(nA),i,j,j_end,n_end,colj
      real*8    x(n),xi
      real*4    A(nA)
c ----------------------------------------------------------------------
      n_end=n-1
      do 100 i=1,n_end
         xi=x(i)
         j_end=MaxA(i+1)-1
         do 110 j=MaxA(i)+1,j_end
            colj=ColA(j)
            x(colj)=x(colj)-A(j)*xi
  110    continue
  100 continue
c ----------------------------------------------------------------------
      return
      end
c **********************************************************************

      subroutine BAKS1R(n,nA,x,A,MaxA,ColA)
      integer*4 n,nA,MaxA(n),ColA(nA),i,j,j_end
      real*8    x(n),xi
      real*4    A(nA)
c ----------------------------------------------------------------------
      do 100 i=n-1,1,-1
         xi=x(i)
         j_end=MaxA(i+1)-1
         do 110 j=MaxA(i)+1,j_end
            xi=xi-A(j)*x(ColA(j))
  110    continue
         x(i)=xi
  100 continue
c ----------------------------------------------------------------------
      return
      end

