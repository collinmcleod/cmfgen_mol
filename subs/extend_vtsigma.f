C
C Routine to interpolate V, T and SIGMA onto a new radius
C
	SUBROUTINE EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NX,
	1                     V,T,SIGMA,ND)
	USE SET_KIND_MODULE
C
	IMPLICIT NONE
	INTEGER NX,ND,INDX(NX)
	REAL(KIND=LDP) VEXT(NX),SIGMAEXT(NX),TEXT(NX),COEF(0:3,NX)
	REAL(KIND=LDP) V(ND),SIGMA(ND),T(ND)
C
	REAL(KIND=LDP) T1
	INTEGER I,J
C
C Do intepolation in the log plane - need to add 1 to Sigma as it can be
C negative.
C
	DO I=1,NX
	  VEXT(I)=0.0D0
	  SIGMAEXT(I)=0.0D0
	  TEXT(I)=0.0D0
	  DO J=0,3
	    VEXT(I)=VEXT(I)+COEF(J,I)*LOG( V(J+INDX(I)) )
	    TEXT(I)=TEXT(I)+COEF(J,I)*LOG( T(J+INDX(I)) )
	    T1= SIGMA(J+INDX(I))+1.0D0
	    IF(T1 .LE. 0.0D0)T1=0.0001D0
	    SIGMAEXT(I)=SIGMAEXT(I)+COEF(J,I)*LOG(T1)
	  END DO
	  VEXT(I)=EXP(VEXT(I))
	  SIGMAEXT(I)=EXP(SIGMAEXT(I))-1.0D0
	  IF(SIGMAEXT(I) .LT. -0.9999D0)SIGMAEXT(I)=-0.9999D0
	  TEXT(I)=EXP(TEXT(I))
	END DO
C
	RETURN	
	END
