	REAL(10) FUNCTION INTSIGC(ENR,XION_POT,EMIN,EMAX)
	IMPLICIT NONE
!
	REAL(10) ENR              ! Energy at which we compute the integral of the cross section
	REAL(10) XION_POT         ! (first) ionization potential of species under consideration
	REAL(10) EMIN,EMAX        ! Bounds of integration
	REAL(10) XJ,ONEOVERJ
!
	INTSIGC = 0.0D0
	IF (EMIN.GE.EMAX) THEN
	   WRITE(6,*) 'Emin >= emax in intsigc - we stop ',emin,emax,xion_pot
	   STOP
	ELSE
	   XJ = 0.6D0 * XION_POT
	   ONEOVERJ =  1.0D0/XJ
	   INTSIGC = (ATAN((EMAX-XION_POT)*ONEOVERJ)-ATAN((EMIN-XION_POT)*ONEOVERJ) ) &
	          / ATAN((ENR-XION_POT)*ONEOVERJ/2.0D0)
	ENDIF
	
	RETURN
	END FUNCTION INTSIGC
