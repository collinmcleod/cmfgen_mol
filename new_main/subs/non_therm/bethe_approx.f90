	SUBROUTINE BETHE_APPROX(Q,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX,FAST_METHOD)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	REAL(KIND=LDP) Q(NKT)
	REAL(KIND=LDP) XKT(NKT)
	REAL(KIND=LDP) dXKT(NKT)
!
	LOGICAL FAST_METHOD
!
	REAL(KIND=LDP), PARAMETER :: PI=3.141592653589793238462643_LDP
	REAL(KIND=LDP), PARAMETER :: A0 = 0.529189379E-8_LDP    		!Bohr radius in cm
	REAL(KIND=LDP), PARAMETER :: Hz_TO_EV=4.1356691_LDP
	REAL(KIND=LDP), PARAMETER :: COL_CONST=13.6_LDP*8.0_LDP*PI*PI*A0*A0/1.732_LDP
!
	REAL(KIND=LDP) GBAR
	REAL(KIND=LDP) X
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) dE
	REAL(KIND=LDP) dE_eV
!
	INTEGER IKT
!
	Q=0.0_LDP
	GBAR=0.2_LDP
!
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0)RETURN
	dE_eV=Hz_to_eV*dE
	T2=3.28978_LDP/dE
	IF(FAST_METHOD)THEN
	  GBAR=0.2
	  T1=3.28978_LDP*COL_CONST*GBAR*ATM(ID)%AXzV_F(NL,NUP)*ATM(ID)%XzV_F(NL,DPTH_INDX)/dE
	  DO IKT=1,NKT
	    IF(XKT(IKT) .GE. dE_eV)THEN
	      Q(IKT)=T1*dXKT(IKT)/XKT(IKT)
	    END IF
	  END DO
	ELSE
	  GBAR=0.2
	  T1=3.28978_LDP*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)*ATM(ID)%XzV_F(NL,DPTH_INDX)/dE
	  DO IKT=1,NKT
	    IF(XKT(IKT) .GE. dE_eV)THEN
	      Q(IKT)=T1*GBAR*dXKT(IKT)/XKT(IKT)
	    END IF
	  END DO
	END IF
!
	RETURN
	END
