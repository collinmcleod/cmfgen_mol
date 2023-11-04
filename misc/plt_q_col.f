	PROGRAM PLT_Q_COL
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NKT=5000
	REAL(KIND=LDP) Q(NKT)
	REAL(KIND=LDP) COL_STRENGTH(NKT)
	REAL(KIND=LDP) XKT(NKT)
	REAL(KIND=LDP) dXKT(NKT)
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
	INTEGER NCOL
	REAL(KIND=LDP) ECOL(100),YCOL(100)
!
	REAL(KIND=LDP) GBAR(NKT)
	REAL(KIND=LDP) X
	REAL(KIND=LDP) U1
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) dE
	REAL(KIND=LDP) dE_eV
!
	REAL(KIND=LDP) ZION
	REAL(KIND=LDP) GLOW
	REAL(KIND=LDP) FVAL
	REAL(KIND=LDP) OMEGA
!
	REAL(KIND=LDP) ARN_HeI(5)
	DATA ARN_HeI/24.60,17.8,-11.000,7.00,-23.20/
!
	INTEGER IKT
!
	DO IKT=1,NKT
	  XKT(IKT)=IKT
	  dXKT(IKT)=1
	END DO
	Q(:)=0.0D0
	gbar(:)=0.0d0
!
	ZION=1.0D0; GLOW=1.0D0
	dE=5.13; FVAL=0.274	
	GLOW=1.0D0
	OMEGA=0.05D0
	CALL GEN_IN(ZION,'Ion charge e.g., 1 for H)')
	CALL GEN_IN(dE,'Energy of transition in units of 10^15 Hz')
	CALL GEN_IN(FVAL,'Oscilator strength')
	CALL GEN_IN(GLOW,'Statistical weight of lower level')
	CALL GEN_IN(OMEGA,'Collsion strength')
!
	dE_eV=Hz_to_eV*dE
	T1=3.28978D0*COL_CONST*FVAL/dE
	DO IKT=1,NKT
	  IF(XKT(IKT) .GE. dE_eV)THEN
	    X = sqrt(xkt(ikt)/dE_eV-1.0d0)
	    IF( ZION .NE. 1 .AND. X .LE. CONNECT_POINT)THEN
	      GBAR(IKT) = 0.2D0
	    ELSE IF(X .LE. 6)THEN
	      GBAR(IKT) = COEF0 + COEF1*X + COEF2*X*X
	    ELSE
	      GBAR(IKT) = 1.077*LOG(X)/1.7725
	    END IF
	    if(gbar(ikt) .ge. 0.0d0)then
	      Q(IKT)=T1*GBAR(ikt)*dXKT(IKT)/XKT(IKT)
	    endif
	  END IF
	END DO
!
	T1=PI*A0*A0*13.6D0/GLOW
	DO IKT=1,NKT
	  COL_STRENGTH(IKT)=Q(IKT)*XKT(IKT)/dXKT(IKT)/T1
	END DO
!
	CALL DP_CURVE(NKT,XKT,GBAR)
	CALL GRAMON_PGPLOT('E(ev)','Gbar',' ',' ')
!
	CALL DP_CURVE(NKT,XKT,COL_STRENGTH)
	OPEN(UNIT=10,FILE='2p.dat',STATUS='OLD',ACTION='READ')
	  READ(10,*)NCOL
	  DO IKT=1,NCOL
	    READ(10,*)ECOL(IKT),YCOL(IKT)
	    ECOL(IKT)=ECOL(IKT)*dE_eV
	  END DO
	CALL DP_CURVE(NCOL,ECOL,YCOL)
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')PG_PEN(2)//'Bethe approximation'//DEF_PEN
	WRITE(6,'(A)')PG_PEN(3)//'Los Alamos collsion strength'//DEF_PEN
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' '
	CALL GRAMON_PGPLOT('E(ev)','OMEGA',' ',' ')
!
	WRITE(6,*)'Plotting excitation cross-section'
	Q=Q*1.0D+18
	CALL DP_CURVE(NKT,XKT,Q)
!
	T1=PI*A0*A0*13.6d0
	T1=T1*OMEGA/GLOW
	DO IKT=1,NKT
	  IF(XKT(IKT) .GE. dE_eV)THEN
	    Q(IKT)=T1*dXKT(IKT)/XKT(IKT)
	  ELSE
	    Q(IKT)=0.0D0
	  END IF
	END DO
	Q=Q*1.0D+18
	CALL DP_CURVE(NKT,XKT,Q)
!
	DO IKT=1,NKT
	  IF(XKT(IKT) .GT. Arn_HeI(1))THEN
	    U1=XKT(IKT)/ARN_HeI(1)
	    T1=(1.0D0-1.0D0/U1)
	    T2=LOG(U1)
	    Q(IKT) = 1.0D-14*( ARN_HeI(2)*T1 + ARN_HeI(3)*T1*T1 +
	1          ARN_HeI(4)*T2 + Arn_HeI(5)*T2/U1 ) / U1 / (ARN_HeI(1)**2)
	  ELSE
	    Q(IKT)=0.0D0
	  END IF
	END DO
	Q=Q*1.0D+18
	CALL DP_CURVE(NKT,XKT,Q)
	CALL GRAMON_PGPLOT('E(ev)','Q(MB)',' ',' ')
!
	STOP
	END
