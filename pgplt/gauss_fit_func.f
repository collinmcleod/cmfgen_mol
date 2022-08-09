!
! Program to evaluate a multiple pseudo-Gaussian fit to a set of
! data points passed by the module GAUSS_FIT_DATA.
!
! The function has the form
!
!   Y =P1 + P2*(X-X(1))+ P5*EXP( -((X-P3)/P4)^P6) + ... 
!
! The routine returns the fit (YFIT) and the SQUARED error.
! between the fit and data.
!
! Altered 09-Aug-2022 : To get consistency inthe different routines  changed to use Gauss.
! Altered   -Sep-2007 : Use of alternative exponent to 2 installed.
! Created 21-Jul-2005
!
	REAL*8 FUNCTION GAUSS_FIT_FUNC(PARAMS)
	USE GAUSS_FIT_DATA
	IMPLICIT NONE
	REAL*8 PARAMS(NG_PAR)
!
	REAL*8 SUM
	REAL*8 T1
	INTEGER I,J,K
!
	GAUSS_FIT_FUNC=0.0D0
	DO J=1,NG_DATA
	  SUM=PARAMS(1)+PARAMS(2)*(X_GAUSS(J)-X_GAUSS(1))
	  DO K=3,NG_PAR,4
	   T1=PARAMS(K+2)*EXP(-(ABS((X_GAUSS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	   IF(ABS(PARAMS(K+3)) .LT. 1.0D0)T1=100.0*T1
	   SUM=SUM+T1
	  END DO
	  YFIT(J)=SUM
	  GAUSS_FIT_FUNC=GAUSS_FIT_FUNC+(Y_GAUSS(J)-SUM)**2
	END DO
!
	RETURN
	END
