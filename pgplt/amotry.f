      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      IMPLICIT NONE
      INTEGER ihi,mp,ndim,np
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),funk
      INTEGER, PARAMETER :: NMAX=20
      EXTERNAL funk
!
! Uses funk
!
      INTEGER j
      REAL*8 fac1,fac2,ytry,ptry(NMAX)
!
      fac1=(1.D0-fac)/ndim
      fac2=fac1-fac
      write(6,'(3ES14.4,I4)')fac,fac1,fac2,ndim
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
!
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5
