	SUBROUTINE DISP_BETHE_APPROX_V5(Q,NL,NUP,XKT,dXKT_ON_XKT,NKT,ID,DPTH_INDX)
	USE MOD_DISP
	IMPLICIT NONE
!
! Altered 08-Feb-2014: Added ATM(ID)%CROSEC_NTFAC to rate computed with collision strngth.
! Altered 15-Nov-2012: Changed to V5:
!                        Call changed -- dXKT_ON_XKT passed in call instead of dXKT.
!                        Test on XSQ rather than X
!                        OMP instructions inserted.
! Altered 10-Feb-2012: Improved computation of gbar. Expression now works for
!                        very low energies, andhigh energies.
!
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER NKT
	INTEGER NL,NUP
	REAL*8 Q(NKT)
	REAL*8 XKT(NKT)
	REAL*8 dXKT_ON_XKT(NKT)
!
	LOGICAL FAST_METHOD
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
	INTEGER NIN,NOUT
	INTEGER J
	INTEGER IKT
!
! NIN & NOUT are used to make the loop over IKT parallizable.
!
	NOUT=16
	NIN=1+(NKT-1)/NOUT
!
!$OMP PARALLEL DO
	DO J=1,NOUT
	  DO IKT=(J-1)*NIN+1,MIN(J*NIN,NKT)
	    Q(IKT)=0.0D0
	  END DO
	END DO
!
	dE=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(dE .LE. 0)RETURN
	dE_eV=Hz_to_eV*dE
!	T2=3.28978D0/dE
!       3.28978 ~ 13.6/Hz_to_eV
!
	IF(ATM(ID)%AXzV_F(NUP,NL) .GT. 1.0D+03 .OR. ATM(ID)%NT_OMEGA(NL,NUP) .EQ. 0.0D0)THEN
	  T1=3.28978D0*COL_CONST*ATM(ID)%AXzV_F(NL,NUP)*ATM(ID)%XzV_F(NL,DPTH_INDX)/dE*ATM(ID)%CROSEC_NTFAC
!
!$OMP PARALLEL DO PRIVATE(GBAR,IKT,X,XSQ)
	  DO J=1,NOUT
	    DO IKT=(J-1)*NIN+1,MIN(J*NIN,NKT)
	      IF(XKT(IKT) .GE. dE_eV)THEN
	        XSQ = XKT(IKT)/dE_eV-1.0D0
	        IF((ATM(ID)%ZXzV .NE. 1) .AND. (XSQ .LE. CONNECT_POINT_SQ) )THEN
	          GBAR = 0.2D0
	        ELSE IF(XSQ .LE. 0.64D0)THEN			!X .LE. 0.8D0
	          X=SQRT(XSQ)
	          GBAR = 0.074D0*X*(1.0D0+X)
	        ELSE IF(XSQ .LE. 36)THEN		  	!X .LE. 6
	          X=SQRT(XSQ)
	          GBAR = COEF0 + X*(COEF1 + COEF2*X)
	        ELSE
	          GBAR = +0.105D0+LOG(XSQ)/3.6276D0     	!Was LOG(X)/1.8138D0 
	        END IF
	        IF(GBAR .GT. 0.0d0)THEN
	          Q(IKT)=T1*GBAR*dXKT_ON_XKT(IKT)
	        ENDIF
	      END IF
	    END DO
	  END DO
	ELSE IF(ATM(ID)%NT_OMEGA(NL,NUP) .NE. 0.0D0)THEN
	  T1=PI*A0*A0*13.6d0*ATM(ID)%CROSEC_NTFAC
	  T1=T1*ATM(ID)%NT_OMEGA(NL,NUP)/ATM(ID)%GXzV_F(NL)*ATM(ID)%XzV_F(NL,DPTH_INDX)
!$OMP PARALLEL DO  
	  DO J=1,NOUT
	    DO IKT=(J-1)*NIN+1,MIN(J*NIN,NKT)
	      IF(XKT(IKT) .GE. dE_eV)THEN
	        Q(IKT)=T1*dXKT_ON_XKT(IKT)
	      END IF
	    END DO
	  END DO
	ELSE
!
	END IF
!
	RETURN
	END      
