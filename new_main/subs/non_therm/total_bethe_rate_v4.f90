	SUBROUTINE TOTAL_BETHE_RATE_V4(RATE,NL,NUP,MOD_YE,XKT,NKT,ID,DPTH_INDX,ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	INTEGER ND
	REAL(KIND=LDP) RATE
	REAL(KIND=LDP) XKT(NKT)
	REAL(KIND=LDP) MOD_YE(NKT,ND)
!
	REAL(KIND=LDP), PARAMETER :: PI=3.141592653589793238462643_LDP
	REAL(KIND=LDP), PARAMETER :: A0 = 0.529189379E-8_LDP    		!Bohr radius in cm
	REAL(KIND=LDP), PARAMETER :: Hz_TO_EV=4.1356691_LDP
	REAL(KIND=LDP), PARAMETER :: COL_CONST=13.6_LDP*8.0_LDP*PI*PI*A0*A0/1.732_LDP
	REAL(KIND=LDP), PARAMETER :: COEF0=-0.0745397_LDP
	REAL(KIND=LDP), PARAMETER :: COEF1=0.232715_LDP
	REAL(KIND=LDP), PARAMETER :: COEF2=-0.00647558_LDP
	REAL(KIND=LDP), PARAMETER :: CONNECT_POINT=1.2212243_LDP
	REAL(KIND=LDP), PARAMETER :: CONNECT_POINT_SQ=CONNECT_POINT*CONNECT_POINT
!
	REAL(KIND=LDP) GBAR
	REAL(KIND=LDP) X,XSQ
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) dE
	REAL(KIND=LDP) dE_eV
!
	INTEGER IKT
!
	RATE=0.0_LDP
	GBAR=0.0_LDP
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0.0_LDP)RETURN
	dE_eV=Hz_to_eV*dE
!
	IF(ATM(ID)%AXzV_F(NUP,NL) .GT. 1.0E+03_LDP .OR. ATM(ID)%NT_OMEGA(NL,NUP) .EQ. 0.0_LDP)THEN
	  T1=3.28978_LDP*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)/dE
	  DO IKT=1,NKT
	    IF(XKT(IKT) .GE. dE_eV)THEN
	      XSQ = XKT(IKT)/dE_eV-1.0_LDP
              IF(( NINT(ATM(ID)%ZXzV) .NE. 1) .AND. (XSQ .LE. CONNECT_POINT_SQ) )THEN
                GBAR = 0.2_LDP
              ELSE IF(XSQ .LE. 0.64_LDP)THEN			!X .LE. 0.8D0
                X=SQRT(XSQ)
	        GBAR = 0.074_LDP*X*(1.0_LDP+X)
              ELSE IF(XSQ .LE. 36.0_LDP)THEN				!X .LE. 6.0D0
                X=SQRT(XSQ)
                GBAR = COEF0 + X*(COEF1 + COEF2*X)
              ELSE
                GBAR = +0.105_LDP+LOG(XSQ)/3.6276_LDP                 !Log(X)/1.8138D0
              END IF
	      IF(GBAR .GT. 0.0_LDP)THEN
	        RATE=RATE+GBAR*MOD_YE(IKT,DPTH_INDX)
	      END IF
	    END IF
	  END DO
	  RATE=RATE*T1
	ELSE IF(ATM(ID)%NT_OMEGA(NL,NUP) .NE. 0.0_LDP)THEN
	  T1=13.6_LDP*PI*A0*A0
	  T1=T1*ATM(ID)%NT_OMEGA(NL,NUP)/ATM(ID)%GXzV_F(NL)
	  DO IKT=1,NKT
	    IF(XKT(IKT) .GE. dE_eV)THEN
	      RATE=RATE+MOD_YE(IKT,DPTH_INDX)
	    END IF
	  END DO
	  RATE=RATE*T1
!	  WRITE(147,'(3A,2ES12.4)')ATM(ID)%XzVLEVNAME_F(NL),'->',ATM(ID)%XzVLEVNAME_F(NUP),ATM(ID)%NT_OMEGA(NL,NUP),RATE
!	  flush(147)
!	  RATE=0.0D0
	ELSE
	END IF
!
	RETURN
	END
