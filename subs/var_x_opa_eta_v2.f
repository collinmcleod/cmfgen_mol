!
! Subroutine to compute the contribution to the opacity AND emissivity
! by K shell ionization.
!
	SUBROUTINE VAR_X_OPA_ETA_V2(VCHI,VETA,
	1              HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1              HN_B,HNST_B,dlnHNST_BdlnT,N_B,
	1              ED,DI,T,EQ_A,EQION,AT_NO,Z_A,
	1              NU,EMHNUKT,NT,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 26-Oct-1995 : dlnHNST_... passed in call instead of edge.
!                         Made version 2.
! Altered 22-Jul-1994 : Extensive modifications and testing.
! Created 20-Jul-1993
!
	EXTERNAL XCROSS_V2
!
	INTEGER ND,NT,N_A,N_B,EQ_A,EQION
	REAL(KIND=LDP) VCHI(NT,ND),VETA(NT,ND)
	REAL(KIND=LDP) HN_A(N_A,ND),HNST_A(N_A,ND),dlnHNST_AdlnT(N_A,ND)
	REAL(KIND=LDP) HN_B(N_B,ND),HNST_B(N_B,ND),dlnHNST_BdlnT(N_B,ND)
	REAL(KIND=LDP) ED(ND),DI(ND),T(ND),EMHNUKT(ND),NU
	REAL(KIND=LDP) AT_NO,Z_A
!
! Functions called.
!
	REAL(KIND=LDP) XCROSS_V2
	INTEGER ERROR_LU
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	INTEGER I,J,LEV
	REAL(KIND=LDP) ALPHA,LTE_POP,NO_ELEC
	REAL(KIND=LDP) TCHI1,TETA2,TETA3
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL ERROR_OUTPUT
	DATA ERROR_OUTPUT/.FALSE./
!
	NO_ELEC=AT_NO-Z_A+1
!
! Add in BOUND-FREE contributions (if any). XCROSS_V2 must return 0
! if frequency is to low to cause ionizations. The cross-section
! is assumed to be independent of the level of the valence electron.
!
	ALPHA=XCROSS_V2(NU,AT_NO,NO_ELEC,IZERO,IZERO,L_FALSE,L_FALSE)
	IF(ALPHA .LE. 0.0_LDP)RETURN
	TETA2=ALPHA*TWOHCSQ*(NU**3)
	IF(NO_ELEC .GT. 3)THEN
	  DO I=1,N_A
	    LEV=EQ_A+I-1
	    DO J=1,ND
	      LTE_POP=HNST_A(I,J)*HNST_B(1,J)/HN_B(1,J)
	      VCHI(LEV,J)=VCHI(LEV,J)+ALPHA
	      TCHI1=ALPHA*LTE_POP*EMHNUKT(J)
	      VCHI(EQION,J)=VCHI(EQION,J)-TCHI1/DI(J)
	      VCHI(NT-1,J)=VCHI(NT-1,J)-2.0_LDP*TCHI1/ED(J)
	      VCHI(NT,J)=VCHI(NT,J)-TCHI1*
	1        (HDKT*NU/T(J)+dlnHNST_AdlnT(I,J)+dlnHNST_BdlnT(1,J))/T(J)
!
	      TETA3=TETA2*LTE_POP*EMHNUKT(J)
	      VETA(EQION,J)=VETA(EQION,J)+TETA3/DI(J)
	      VETA(NT-1,J)=VETA(NT-1,J)+2.0_LDP*TETA3/ED(J)
	      VETA(NT,J)=VETA(NT,J)+TETA3*
	1        (HDKT*NU/T(J)+dlnHNST_AdlnT(I,J)+dlnHNST_BdlnT(1,J))/T(J)
	    END DO
	  END DO
	ELSE		!Variation for K shell ionization of Li ions.
	  IF(.NOT. ERROR_OUTPUT)THEN
	    WRITE(ERROR_LU(),*)'**************************************'
	    WRITE(ERROR_LU(),*)'General K shell ionization for Lithium'//
	1              'ions is not yet treated in VAR_X_OPA_ETA'
	  END IF
	  ERROR_OUTPUT=.TRUE.
	END IF
!
	RETURN
	END
