C
C Subroutine to compute the Z values at NI depth points for a given
C imact parameter P .
C
	SUBROUTINE ZALONGP(R,X,P,NI)
	IMPLICIT NONE
C
C Altered 09-Dec-2001 : PP removed and repplaced directly in SQRT by P*P
C Altered 28-May-1996 : IMPLCIT NONE installed.
C
	INTEGER*4 NI
	REAL*8 X(NI),R(NI),P
C
	INTEGER*4 I
C
	DO 10 I=1,NI
	  X(I)=SQRT(R(I)*R(I)-P*P)
10	CONTINUE
C
	RETURN
	END
