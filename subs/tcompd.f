C
C Subroutine to compute T matrix for continuum solution of
C the co-moving radiative transfer equation. For the continuum
C the effect of the velocity gradient can be ignored.
C Uses a Schuster condition or the diffusion approximation
C for the lower boundary condition.
C
	SUBROUTINE TCOMPD(TA,TB,TC,DTAU,DIFF,LS,NC,ND,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 27-May-1996 : IMPLICIT NONE installed.
C Altered 07-Jan-91 - TA,TB,TC loop split so CRAY vectorizes loop.
C Altered 28-JUL-82
C
	INTEGER LS,NC,ND,NI
	REAL(KIND=LDP) TA(NI),TB(NI),TC(NI),DTAU(NI)
	LOGICAL DIFF
C
	INTEGER I
C
	TA(1)=0.0_LDP
	TC(1)=1.0_LDP/DTAU(1)
	TB(1)=-1.0_LDP-TC(1)
C
	DO I=2,NI-1
	  TC(I)=1.0_LDP/DTAU(I)
	END DO
	DO I=2,NI-1
	  TA(I)=TC(I-1)
	  TB(I)=-0.5_LDP*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	END DO
C
	IF(LS .LE. NC .AND. DIFF)THEN
		TA(NI)=-TC(NI-1)
		TB(NI)=TC(NI-1)
	ELSE IF(LS .GT. NC)THEN
		TA(NI)=-TC(NI-1)
		TB(NI)=-TA(NI)+DTAU(NI-1)/2.0_LDP
	ELSE
		TA(NI)=-TC(NI-1)
		TB(NI)=1.0_LDP-TA(NI)
	END IF
	TC(NI)=0.0_LDP
C
	RETURN
	END
