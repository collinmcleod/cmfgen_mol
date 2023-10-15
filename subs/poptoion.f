C
C Routines to store the populations in the general POPS array, or
C to take them from POPS and put them in the individual arrays
C for each species.
C
	SUBROUTINE POPTOION(POPS,HYD,DHYD,ED,T,
	1                      EQHYD,NLEV,NT,ND,HYD_PRES)
	IMPLICIT NONE
C
	LOGICAL HYD_PRES
	INTEGER EQHYD,NLEV,NT,ND
	REAL(10) POPS(NT,ND)
	REAL(10) HYD(NLEV,ND),DHYD(ND),ED(ND),T(ND)
C
	INTEGER I,J,K
C
	IF(HYD_PRES)THEN
	  DO K=1,ND
	    DO J=1,NLEV
	      I=EQHYD+J-1
	      HYD(J,K)=POPS(I,K)
	    END DO
	    DHYD(K)=POPS(EQHYD+NLEV,K)
	    ED(K)=POPS(NT-1,K)
	    T(K)=POPS(NT,K)
	  END DO
	END IF
C
	RETURN
	END
C
C 
C
	SUBROUTINE IONTOPOP(POPS,HYD,DHYD,ED,T,
	1                      EQHYD,NLEV,NT,ND,HYD_PRES)
	IMPLICIT NONE
C
	LOGICAL HYD_PRES
	INTEGER EQHYD,NLEV,NT,ND
	REAL(10) POPS(NT,ND)
	REAL(10) HYD(NLEV,ND),DHYD(ND),ED(ND),T(ND)
C
	INTEGER I,J,K
C
	IF(HYD_PRES)THEN
	  DO K=1,ND
	    DO J=1,NLEV
	      I=EQHYD+J-1
	      POPS(I,K)=HYD(J,K)
	    END DO
	    POPS(EQHYD+NLEV,K)=DHYD(K)
	    POPS(NT-1,K)=ED(K)
	    POPS(NT,K)=T(K)
	  END DO
	END IF
C
	RETURN
	END
