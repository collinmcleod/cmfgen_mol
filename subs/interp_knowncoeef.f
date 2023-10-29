C
C Created 10-Feb-1988.
C
	SUBROUTINE INTERP_KNOWNCOEEF(CHIEXT,COEF,INDX,NX,CHI,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER NX,ND,INDX(NX),I,J
	REAL(KIND=LDP) CHIEXT(NX),CHI(ND),COEF(0:3,NX)
C
	DO I=1,NX
	  CHIEXT(I)=0.0D0
	  DO J=0,3
	    CHIEXT(I)=CHIEXT(I)+COEF(J,I)*LOG( CHI(J+INDX(I)) )
	  END DO
	  CHIEXT(I)=EXP(CHIEXT(I))
	END DO
C
	RETURN
	END
