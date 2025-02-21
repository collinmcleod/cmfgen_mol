	PROGRAM TST_GEN_IN
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
	INTEGER, PARAMETER :: NMAX=10
	REAL(KIND=LDP) VEC(NMAX)
	REAL*4 SP(NMAX)
	INTEGER INTV(NMAX)
	INTEGER NSP,NDP,NI,I
C
	DO I=1,NMAX
	  SP(I)=I*I+0.1
	  VEC(I)=I*I+0.3
	  INTV(I)=I*I
	END DO
C
	CALL GEN_IN(SP,NSP,NMAX,'Single precision Values')
	CALL GEN_IN(VEC,NDP,NMAX,'Double precision Values')
	CALL GEN_IN(INTV,NI,NMAX,'Integer Values')
C
	DO I=1,NSP
	  WRITE(6,*)I,SP(I)
	END DO
	DO I=1,NDP
	  WRITE(6,*)I,VEC(I)
	END DO
	DO I=1,NI
	  WRITE(6,*)I,INTV(I)
	END DO
C
	STOP
	END
