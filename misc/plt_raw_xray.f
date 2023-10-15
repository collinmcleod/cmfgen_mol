	PROGRAM PLT_RAW_XRAY
	USE MOD_XRAY_FLUXES
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL(10) TSHOCK_1
	REAL(10) TSHOCK_2
	REAL(10) CONV_FACT
	INTEGER, PARAMETER :: LUIN=7
!
	CONV_FACT=4.0D-10*ACOS(-1.0D0)
	TSHOCK_1=100.0D0
	TSHOCK_2=200.0D0
	CALL RD_XRAY_SPEC(TSHOCK_1,TSHOCK_2,LUIN)
	X_EMISS1=LOG10(CONV_FACT*X_EMISS1+1.0D-300)
	X_EMISS2=LOG10(CONV_FACT*X_EMISS2+1.0D-300)
	CALL DP_CURVE(N_BINS,X_NU,X_EMISS1)
	CALL GRAMON_PGPLOT('Freq','Log \gL(T=10\u6\dK)',' ',' ')
!
	CALL DP_CURVE(N_BINS,X_NU,X_EMISS2)
	CALL GRAMON_PGPLOT('Freq','Log \gL(T=2 x 10\u6\dK)',' ',' ')
!
	TSHOCK_1=400.0D0
	TSHOCK_2=1000.0D0
	CALL RD_XRAY_SPEC(TSHOCK_1,TSHOCK_2,LUIN)
	X_EMISS1=LOG10(CONV_FACT*X_EMISS1+1.0D-300)
	X_EMISS2=LOG10(CONV_FACT*X_EMISS2+1.0D-300)
	CALL DP_CURVE(N_BINS,X_NU,X_EMISS2)
	CALL GRAMON_PGPLOT('Freq','Log \gL(T=4 x 10\u6\dK)',' ',' ')
!
	CALL DP_CURVE(N_BINS,X_NU,X_EMISS2)
	CALL GRAMON_PGPLOT('Freq','Log \gL(T=10\u7\dK)',' ',' ')
!
	STOP
	END 
