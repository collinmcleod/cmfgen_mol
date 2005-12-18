      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      IMPLICIT NONE
      INTEGER iter,mp,ndim,np
      REAL*8 ftol
      REAL*8 p(mp,np)
      REAL*8 y(mp)
      REAL*8 funk
      INTEGER, PARAMETER :: NMAX=20
      INTEGER, PARAMETER :: ITMAX=20000
      EXTERNAL funk
!
! Subroutine uses amotry,funk
!
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
!
      iter=0
!      do j=1,ndim
!	do i=1,ndim+1
!	   write(15,*)p(i,j)
!        end do
!      end do
!
! Main iteration loop.
!
1     continue
      do n=1,ndim
        sum=0.0d0
        do m=1,ndim+1
          sum=sum+p(m,n)
        end do
        psum(n)=sum
      end do
!
2     continue
      ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
!
      do i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
      end do
!
! Check whether desired tolerance has been obtained. If so,
! update and return.
!
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
!
      if (iter.ge.ITMAX)then
        write(6,*)'ITMAX exceeded in amoeba'
        write(6,*)'returning'
        return
      end if
!
! Continue iteration procdure.
!
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)
      write(6,*)'done one'
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)
        if (ytry.ge.ysave) then
          do i=1,ndim+1
            if(i.ne.ilo)then
              do j=1,ndim
                psum(j)=0.5D0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
              end do
              y(i)=funk(psum)
            end if
          end do
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 
