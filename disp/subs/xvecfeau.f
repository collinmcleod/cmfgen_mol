C
C Subroutine to compute the X vector (see notes) . This routine
C uses either the diffusion approximation or a Schuster condition
C for the inner boundary condition.
C
C Created 28-JUL-82
C Changed 31-Oct-86  (Schuster b.c. condition inserted)
C
	SUBROUTINE XVECFEAU(X,R,Q,SOURCE,DIFF,DBB,INBC,IC,CHI,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER ND,I
	REAL(KIND=LDP) X(ND),R(ND),Q(ND),SOURCE(ND),DBB,CHI,INBC,IC
	LOGICAL DIFF
C
	X(1)=0.0
	DO I=2,ND-1
	  X(I)=R(I)*R(I)*SOURCE(I)/Q(I)
	END DO
C
C Note well - DBB =dB/dR (and Q(ND)=1.0 by definition)
C
	IF(DIFF)THEN
	  X(ND)=R(ND)*R(ND)*DBB/CHI/3.0_LDP
	ELSE
	  X(ND)=R(ND)*R(ND)*IC*(0.25_LDP+0.5_LDP*INBC)
	END IF
C
	RETURN
	END
