!
! Function to calclate the exponential integral for
! n=0, 1, 2, 3 etc. Accuracy has been set to 1.0D-10
!
!
	FUNCTION EXPN(N,X)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 05-Apr-2006
!
! Referemces: Abramowitz and Stegun
!             Numerical recipes.
!
	INTEGER N
	REAL(KIND=LDP) X
	REAL(KIND=LDP) EXPN
!
	REAL(KIND=LDP), PARAMETER :: EULER=0.5772156649015328606D0
	REAL(KIND=LDP), PARAMETER :: EPS=1.0D-10
	INTEGER, PARAMETER :: MAX_STEPS=50
!
	REAL(KIND=LDP) Y,FAC,ANS,VAL,T1,OLD_ANS
	REAL(KIND=LDP) A0,A1,A2
	REAL(KIND=LDP) B0,B1,B2
	INTEGER I,J
!
	IF(N .EQ. 0)THEN
	  EXPN=EXP(-X)/X
	ELSE IF(X .EQ. 0.0D0 .AND. N .EQ .0)THEN
	  WRITE(6,*)'EXP1(X) is undefined for X=0'
	  STOP
	ELSE IF(X .EQ. 0.0D0)THEN
	  EXPN=EXP(-X)/(N+X)
	ELSE IF(X .LE. 1.5D0)THEN
	  ANS=0.0D0
	  FAC=1.0D0
	  DO I=1,N-1
	    ANS=ANS+1.0D0/I
	    FAC=FAC*I
	  END DO
	  ANS=(ANS-EULER-LOG(X))*(-X)**(N-1)/FAC
	  Y=1.0D0
	  VAL=0.0D0
	  IF(N .NE. 1)VAL=Y/(1.0D0-N)
	  DO J=1,MAX_STEPS
	    Y=-Y*X/J
	    IF(J .NE. N-1)THEN
	      T1=Y/(J-N+1)
	      VAL=VAL+T1
	      IF(ABS(T1) .LT. EPS)EXIT
	    END IF
	  END DO
	  EXPN=ANS-VAL
	ELSE
!
! In some cases, A0, A1, B0, and B1 mght need to be rescaled.
!
!
	  A0=0.0D0; B0=1.0D0
	  A1=1.0D0; B1=X+N
	  ANS=0.0D0
	  DO I=1,MAX_STEPS
	    A2=A1*(X+N+2*I)-A0*I*(N+I-1)
	    B2=B1*(X+N+2*I)-B0*I*(N+I-1)
	    OLD_ANS=ANS
	    ANS=A2/B2
	    IF(ABS(ANS-OLD_ANS) .LT. EPS)EXIT
	    A0=A1; B0=B1
	    A1=A2; B1=B2
	  END DO
	  EXPN=ANS*EXP(-X)
	END IF
!
	RETURN
	END
