!
! Subroutine to compute dMU. Designed to compute dMU of higher accuracy
! when MU is close to 1.
!
	SUBROUTINE SET_ACC_dMU(X,dX,P,RVAL,N)
	IMPLICIT NONE
!
! Created 14-Jun-2023
!
	INTEGER N
	REAL(10) X(N),dX(N)
	REAL(10) P(N),RVAL
!
	REAL(10) A(8)
	REAL(10) T1,T2
	INTEGER I
!
! Compute the coeficients in the binomial expansion of 1-sqrt(1-x^2)
!
	A(1)=0.5D0
	A(2)=1.0D0/8.0D0
	A(3)=1.0D0/16.0D0
	A(4)=5.0D0/128.0D0
	A(5)=7.0D0/256.0D0
	A(6)=21.0D0/1024.0D0
!
! Compute dMU. We check the size of X(I+1) since it has the larger P/RVAL, 
! and hence is more influenced by the accuracy Taylor series expansion. 
!
	dX(N)=0.0D0
	DO I=1,N-1
          IF(X(I+1) .LT. 0.999D0)THEN
	    dX(I)=X(I)-X(I+1)
	  ELSE
	    T1=(P(I)/RVAL)**2	
	    T2=(P(I+1)/RVAL)**2
	    T1=T1*(A(1)+T1*(A(2)+T1*(A(3)+T1*(A(4)+T1*(A(5)+T1*A(6))))))
	    T2=T2*(A(1)+T2*(A(2)+T2*(A(3)+T2*(A(4)+T2*(A(5)+T2*A(6))))))
	    dX(I)=T2-T1
	  END IF
	END DO 
!
	RETURN
	END
