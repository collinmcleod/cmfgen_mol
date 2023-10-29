	FUNCTION INTSIGC(ENR,XION_POT,EMIN,EMAX)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) INTSIGC 
!
	REAL(KIND=LDP) ENR              ! Energy at which we compute the integral of the cross section
	REAL(KIND=LDP) XION_POT         ! (first) ionization potential of species under consideration
	REAL(KIND=LDP) EMIN,EMAX        ! Bounds of integration
	REAL(KIND=LDP) XJ,ONEOVERJ
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
