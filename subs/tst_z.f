	PROGRAM TST_Z
	IMPLICIT NONE
C
	INTEGER, PARAMETER :: N=10000000
	REAL(10) A(N)
	REAL(10) B(N)
	REAL(10) ANS(N)
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
