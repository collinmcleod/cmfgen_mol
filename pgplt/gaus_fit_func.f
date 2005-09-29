!
! Program to evaluate a multiple Gaussian fit to a set of
! data points passed by the module GAUS_FIT_DATA.
!
! The function has the form
!
!   Y =P1 + P2*(X-X(1))+ P5*EXP( -((X-P3)/P4)^2) + ... 
!
! The routine returns the fit (YFIT) and the SQUARED error.
! between the fit and data.
!
! Created 21-Jul-2005
!
	REAL*8 FUNCTION GAUS_FIT_FUNC(PARAMS)
	USE GAUS_FIT_DATA
	IMPLICIT NONE
	REAL*8 PARAMS(NG_PAR)
!
	REAL SUM
	INTEGER I,J,K
!
	GAUS_FIT_FUNC=0.0D0
	DO J=1,NG_DATA
	  SUM=PARAMS(1)+PARAMS(2)*(X_GAUS(J)-X_GAUS(1))
	  DO K=3,NG_PAR,3
	   SUM=SUM+PARAMS(K+2)*EXP(-((X_GAUS(J)-PARAMS(K))/PARAMS(K+1))**2)
	  END DO
	  YFIT(J)=SUM
	  GAUS_FIT_FUNC=GAUS_FIT_FUNC+(Y_GAUS(J)-SUM)**2
	END DO
!
	RETURN
	END
