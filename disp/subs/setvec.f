	SUBROUTINE SETVEC(X,Y,ND)
	IMPLICIT NONE
	INTEGER*4 ND,I
	REAL*8 X(ND),Y(ND)
C
	DO I=1,ND
	  X(I)=Y(I)
	END DO
C
	RETURN
	END
