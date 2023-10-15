!
! The rates in this routine must be kept consistent with TOTAL_BETHE_RATE_V2(
!
	SUBROUTINE BETHE_APPROX_V3(Q,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
!
! Altered 10-Feb-2012: Improved computation of gbar. Expression now works for
!                        very low energies, andhigh energies.
! Altered 06-Nov-2011: Cleaned.
!                        GBAR changed from vector to scaler
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	REAL(10) Q(NKT)
	REAL(10) XKT(NKT)
	REAL(10) dXKT(NKT)
!
	LOGICAL FAST_METHOD
!
	REAL(10), PARAMETER :: PI=3.141592653589793238462643D0
	REAL(10), PARAMETER :: A0 = 0.529189379D-8    		!Bohr radius in cm
	REAL(10), PARAMETER :: Hz_TO_EV=4.1356691D0
	REAL(10), PARAMETER :: COL_CONST=13.6D0*8.0D0*PI*PI*A0*A0/1.732D0
	REAL(10), PARAMETER :: COEF0=-0.0745397d0
	REAL(10), PARAMETER :: COEF1=0.232715d0
	REAL(10), PARAMETER :: COEF2=-0.00647558d0
	REAL(10), PARAMETER :: CONNECT_POINT=1.2212243D0
!
	REAL(10) GBAR
	REAL(10) X
	REAL(10) T1,T2
	REAL(10) dE
	REAL(10) dE_eV
!
	INTEGER IKT
!
	Q(:)=0.0D0
	GBAR=0.0d0
!
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0)RETURN
	dE_eV=Hz_to_eV*dE
	T2=3.28978D0/dE
!
	T1=3.28978D0*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)*ATM(ID)%XzV_F(NL,DPTH_INDX)/dE*ATM(ID)%CROSEC_NTFAC
	DO IKT=1,NKT
	  IF(XKT(IKT) .GE. dE_eV)THEN
	    X = SQRT(XKT(IKT)/dE_eV-1.0d0)
	    IF((ATM(ID)%ZXzV .NE. 1).AND.(X .LE. CONNECT_POINT))THEN
	      GBAR = 0.2D0
	    ELSE
	      GBAR = COEF0 + COEF1*X + COEF2*X*X
	    END IF
	    IF(GBAR .GE. 0.0D0)THEN
	      Q(IKT)=T1*GBAR*dXKT(IKT)/XKT(IKT)
	    ENDIF
	  END IF
	END DO
!
	RETURN
	END
