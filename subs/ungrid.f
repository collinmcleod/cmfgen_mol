C
C Subroutine to place desired quntity on old radius grid.
C
	SUBROUTINE UNGRID(RJ,ND,RJEXT,NX,GRID)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER ND,NX,GRID(ND),I
	REAL(KIND=LDP) RJ(ND),RJEXT(NX)
C
	DO I=1,ND
	  RJ(I)=RJEXT(GRID(I))
	END DO
C
	RETURN
	END
