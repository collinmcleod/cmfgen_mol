	FUNCTION LAMBERT_WM_FUN(Z,LNZ)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) LAMBERT_WM_FUN
!
! Function to compute the lambet function W_{-1}(Z).
!
	REAL(KIND=LDP) Z
	REAL(KIND=LDP) LNZ
	REAL(KIND=LDP) W
!
	REAL(KIND=LDP) EN
	REAL(KIND=LDP) QN
	REAL(KIND=LDP) ZN
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) WOLD
!
	INTEGER J
!
! Make an initial guess -- these are god enough for the
! iteration procdures to converge in a few iterations.
!
	T1=-EXP(-1.0D0)
	IF(LNZ .LT. -100)THEN
	  W=LNZ
	  W=W-LOG(-W)*(1.0D0-1.0D0/W)
	ELSE IF(Z .LT. T1)THEN
	  WRITE(6,*)'Invalid value of z in LAMBERT_WM_FUN. Z=',Z
	  WRITE(6,*)'Argument must be greater (or =) to -EXP(-1.0D0)'
	  STOP
	ELSE IF(Z .EQ. T1)THEN
	  W=-1.0D0
	ELSE
	  IF(Z .GT. -0.15D0)THEN
	    W=LOG(-Z)
	    W=W-LOG(-W)*(1.0D0-1.0D0/W)
	  ELSE
	    W=-1.0D0-9.1D0*(EXP(-1.0D0)+Z)
	  END IF
!
	  WOLD=1
	  DO J=1,10
	    ZN=LNZ-LOG(-W)-W
	    QN=2*(1+W)*(1+W+2*ZN/3.0D0)
	    EN=ZN/(1+W)*(QN-ZN)/(QN-2*ZN)
	    W=W*(1.0D0+EN)
	    IF(W .GT. 0.0D0)W=-W
	    IF( ABS( (W-WOLD)/WOLD ) .LT. 1.0D-08)EXIT
	    WOLD=W
	  END DO
	END IF
!
	LAMBERT_WM_FUN=W
!
	RETURN
	END
