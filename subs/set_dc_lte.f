	SUBROUTINE SET_DC_LTE(XzVLTE,DXzV,EDGE,NXzV,T,TMIN,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NXzV
	INTEGER ND
	REAL(KIND=LDP) XzVLTE(NXzV,ND)
	REAL(KIND=LDP) DXzV(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) EDGE(NXzV)
	REAL(KIND=LDP) TMIN
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	INTEGER I,J
!
	DO I=1,ND
	  DXzV(I)=1.0_LDP
	  IF(T(I) .GT. TMIN)THEN
	    XzVLTE(:,I)=1.0_LDP
	  ELSE
	    DO J=1,NXzV
	      XzVLTE(J,I)=EXP(HDKT*EDGE(J)*(1.0_LDP/TMIN-1.0_LDP/T(I)))
	    END DO
	  END IF
	END DO
!
	RETURN
	END
