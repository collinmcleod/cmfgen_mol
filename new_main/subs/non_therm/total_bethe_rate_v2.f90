!
! The rates in this routine must be kept consistent with BETHE_APPROX.
!
	SUBROUTINE TOTAL_BETHE_RATE_V2(RATE,NL,NUP,YE,XKT,dXKT,NKT,ID,DPTH_INDX,ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered: 06-Nov-2011 : Inserted ATM(ID)%CROSEC_NTFAC into cross-section.
!                        Some cleaning done.
!			 GBAR changed from vector to scaler.
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	INTEGER ND
	REAL(KIND=LDP) RATE
	REAL(KIND=LDP) XKT(NKT)
	REAL(KIND=LDP) dXKT(NKT)
	REAL(KIND=LDP) YE(NKT,ND)
!
	REAL(KIND=LDP), PARAMETER :: PI=3.141592653589793238462643D0
	REAL(KIND=LDP), PARAMETER :: A0 = 0.529189379D-8    		!Bohr radius in cm
	REAL(KIND=LDP), PARAMETER :: Hz_TO_EV=4.1356691D0
	REAL(KIND=LDP), PARAMETER :: COL_CONST=13.6D0*8.0D0*PI*PI*A0*A0/1.732D0
	REAL(KIND=LDP), PARAMETER :: COEF0=-0.0745397d0
	REAL(KIND=LDP), PARAMETER :: COEF1=0.232715d0
	REAL(KIND=LDP), PARAMETER :: COEF2=-0.00647558d0
	REAL(KIND=LDP), PARAMETER :: CONNECT_POINT=1.2212243D0
!
	REAL(KIND=LDP) GBAR
	REAL(KIND=LDP) X
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) dE
	REAL(KIND=LDP) dE_eV
!
	INTEGER IKT
!
	RATE=0.0D0
	GBAR=0.0d0
!
	IF(ATM(ID)%AXzV_F(NL,NUP) .EQ. 0.0D0)RETURN
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0)RETURN
!
	dE_eV=Hz_to_eV*dE
	T1=3.28978D0*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)*ATM(ID)%CROSEC_NTFAC/dE
	DO IKT=1,NKT
	  IF(XKT(IKT) .GE. dE_eV)THEN
	    X = SQRT(XKT(IKT)/dE_eV-1.0d0)
	    IF((ATM(ID)%ZXzV .NE. 1).AND.(X .LE. CONNECT_POINT))THEN
	      GBAR = 0.2D0
	    ELSE
	      GBAR = COEF0 + COEF1*X + COEF2*X*X
	    END IF
	    IF(GBAR .GT. 0.0D0)THEN
	      RATE=RATE+T1*GBAR*YE(IKT,DPTH_INDX)*dXKT(IKT)/XKT(IKT)
	    end if
	  END IF
	END DO
!
	RETURN
	END
