	PROGRAM TST_Z
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER, PARAMETER :: N=10000000
	REAL(KIND=LDP) A(N)
	REAL(KIND=LDP) B(N)
	REAL(KIND=LDP) ANS(N)
C
	INTEGER I,J,ISEED
C
	call random_seed
	call random_number(A)
	call random_number(B)
	CALL ZERO(B,N)
C
	CALL TUNE(1,'TST')
	DO J=1,10
	  DO I=1,N
	    IF(B(I) .NE. 0)ANS(I)=A(I)*B(I)
	  END DO
	END DO
	CALL TUNE(2,'TST')
	CALL TUNE(3,'TST')
C
	STOP
	END
