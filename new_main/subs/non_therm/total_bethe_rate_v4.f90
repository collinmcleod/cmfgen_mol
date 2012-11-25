	SUBROUTINE TOTAL_BETHE_RATE_V4(RATE,NL,NUP,MOD_YE,XKT,NKT,ID,DPTH_INDX,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	INTEGER ND
	REAL*8 RATE
	REAL*8 XKT(NKT)
	REAL*8 MOD_YE(NKT,ND)
!
	REAL*8, PARAMETER :: PI=3.141592653589793238462643D0
	REAL*8, PARAMETER :: A0 = 0.529189379D-8    		!Bohr radius in cm
	REAL*8, PARAMETER :: Hz_TO_EV=4.1356691D0
	REAL*8, PARAMETER :: COL_CONST=13.6D0*8.0D0*PI*PI*A0*A0/1.732D0
	REAL*8, PARAMETER :: COEF0=-0.0745397d0
	REAL*8, PARAMETER :: COEF1=0.232715d0
	REAL*8, PARAMETER :: COEF2=-0.00647558d0
	REAL*8, PARAMETER :: CONNECT_POINT=1.2212243D0
	REAL*8, PARAMETER :: CONNECT_POINT_SQ=CONNECT_POINT*CONNECT_POINT
!
	REAL*8 GBAR
	REAL*8 X,XSQ
	REAL*8 T1,T2
	REAL*8 dE
	REAL*8 dE_eV
!
	INTEGER IKT
!
	RATE=0.0D0
	GBAR=0.0D0
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0)RETURN
	dE_eV=Hz_to_eV*dE
!
	IF(ATM(ID)%AXzV_F(NUP,NL) .GT. 1.0D+03 .OR. ATM(ID)%NT_OMEGA(NL,NUP) .EQ. 0.0D0)THEN
	  T1=3.28978D0*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)/dE
	  DO IKT=1,NKT
	    IF(XKT(IKT) .GE. dE_eV)THEN
	      XSQ = XKT(IKT)/dE_eV-1.0D0
              IF((ATM(ID)%ZXzV .NE. 1) .AND. (XSQ .LE. CONNECT_POINT_SQ) )THEN
                GBAR = 0.2D0
              ELSE IF(XSQ .LE. 0.64D0)THEN			!X .LE. 0.8D0
                X=SQRT(XSQ)
	        GBAR = 0.074*X*(1.0D0+X)
              ELSE IF(XSQ.LE. 36)THEN				!X .LE. 6.0D0
                X=SQRT(XSQ)
                GBAR = COEF0 + X*(COEF1 + COEF2*X)
              ELSE
                GBAR = +0.105D0+LOG(XSQ)/3.6276D0                 !Log(X)/1.8138D0
              END IF
	      IF(GBAR .GT. 0.0D0)THEN
	        RATE=RATE+GBAR*MOD_YE(IKT,DPTH_INDX)
	      END IF
	    END IF
	  END DO
	  RATE=RATE*T1
	ELSE IF(ATM(ID)%NT_OMEGA(NL,NUP) .NE. 0.0D0)THEN
	  T1=13.6d0*PI*A0*A0
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
