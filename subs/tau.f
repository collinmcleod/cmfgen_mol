C
C Subroutine to compute the increment in optical depth along a given ray.
C
C Altered 11-Apr-1989 - Implicit NONE installed.
C
	SUBROUTINE  TAU(DTAU,RKI,Z,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER NI,I
	REAL(KIND=LDP) DTAU(NI),RKI(NI),Z(NI)
C
	DO 10 I=1,NI-1
	  DTAU(I)=0.5_LDP*(RKI(I)+RKI(I+1))*(Z(I)-Z(I+1))
10	CONTINUE
C
	RETURN
	END
