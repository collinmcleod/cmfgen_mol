!
! Given the lower and upper limits of integration X1 and X2, and given N, this
! routine returns X and W of LENGTH N containing the abscissas and weights of
! the Gauss-Legendre N-point quadrature formulae.
!
! Adapted from Numerical recipes.
!
	SUBROUTINE GAULEG(X1,X2,X,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 10-Apr-1995
!
	INTEGER N
	REAL(KIND=LDP) X1,X2,X(N),W(N)
!
! Local variables.
!
	REAL(KIND=LDP), PARAMETER :: EPS=3.0E-14_LDP
!
	REAL(KIND=LDP) XM,XL
	REAL(KIND=LDP) P1,P2,P3,dP
	REAL(KIND=LDP) Z,Z1
	REAL(KIND=LDP) ANS,SUM
	INTEGER I,J,M
!
! Check sensible limits.
!
	IF(ABS(X2-X1) .LT. EPS)THEN
	  WRITE(5,*)'Invalid intgration limits in GAULEG'
	  WRITE(5,*)'X2=',X2
	  WRITE(5,*)'X1=',X1
	END IF
!
! XM and XL are used to convert the abscissa and weights to the desired
! integration interval.
!
	XM=0.5_LDP*(X2+X1)
	XL=0.5_LDP*(X2-X1)
!
	M=(N+1)/2
	DO I=1,M			!Loop to find desired roots.
	  Z=COS(3.141592654_LDP*(I-.25_LDP)/(N+.5_LDP))
	  Z1=100.0
!
! starting with the approximation Z given above we reine the Roots using
! Newton's method.
!
	  DO WHILE( ABS(Z-Z1).GT.EPS)
	    P1=1._LDP
	    P2=0._LDP
!
! Loop to evaluate the recurrance relation to get the Legendre polynomial
! evaluated at Z.
!
! P1 is the desired polynomial.
! P2 is the polynomial of one lower order.
!
	    DO J=1,N
	      P3=P2
	      P2=P1
	      P1=((2._LDP*J-1._LDP)*Z*P2-(J-1._LDP)*P3)/J
	    END DO
!
! Copute the derivative, dP, of P1
!
	    dP=N*(Z*P1-P2)/(Z*Z-1._LDP)
	    Z1=Z
	    Z=Z1-P1/dP
	  END DO
	  X(I)=XM-XL*Z
	  X(N+1-I)=XM+XL*Z
	  W(I)=2._LDP*XL/((1._LDP-Z*Z)*dP*dP)
	  W(N+1-I)=W(I)
	END DO
!
! Check that evrything has worked. Qudrature should be "exact" for all
! polynomials of degree 2N-1 or less.
!
! We compare the error to the absolute largest term in the formulae for ANS,
! since it is possible for the integration to be identically zero. We use P1
! as a temporary variable.
!
	DO J=0,2*N-1
	  SUM=0.0_LDP
	  DO I=1,N
	    SUM=SUM+W(I)*(X(I)**J)
	  END DO
	  ANS=(X2**(J+1)-X1**(J+1))/(J+1)
	  P1=MAX( ABS(X2**(J+1)), ABS(X1**(J+1)) )/(J+1)
	  IF( ABS(ANS-SUM)/P1 .GT. 1.0E-12_LDP)THEN
	    WRITE(6,*)'Error in GAULEG --- bad quadrature computation.'
	    WRITE(6,*)'J=',J,'ANS=',ANS,'SUM=',SUM
	  END IF
	END DO
!
	RETURN
	END
