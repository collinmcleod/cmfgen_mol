C
C Set the SAVE vector for an NG acceleration.
C
	SUBROUTINE SETFORNG(RJ,PREVRJ,ITNUM,ND)
	IMPLICIT NONE
	INTEGER*4 I,ITNUM,ND
	REAL*8 RJ(ND),PREVRJ(4,ND)
C
	DO I=1,ND
	  PREVRJ(ITNUM,I)=RJ(I)
	END DO
C
	RETURN
	END	
