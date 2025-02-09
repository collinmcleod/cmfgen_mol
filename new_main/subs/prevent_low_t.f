!
! Prevent T from becoming too small. When T becomes small, we can get floating
! overflow in some sections of the code. To prevent this, we have added an extra
! local heating term, designed to ensure that T > T_MIN (approximately). We
! have used the exponential term to get a continuous heating term, but which
! goes to zero quickly when T> T_MIN. This may have some problems if the
! T wants to become very small, but it is still well above T_MIN where the
! additional heating term is not effective.
!
	SUBROUTINE PREVENT_LOW_T(HEAT,T_MIN,COMPUTE_BA,LAMBDA_ITERATION,
	1                           T_MIN_EXTRAP,DIAG_INDX,NUM_BNDS,ND,NT)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
        USE MOD_CMFGEN
        IMPLICIT NONE
!
! Altered 15-Jul-2-19 : Added BA_T_EHB equation.
! Created 24-Sep-2004
!
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER NT
!
	LOGICAL COMPUTE_BA
	LOGICAL LAMBDA_ITERATION
	LOGICAL T_MIN_EXTRAP
!
	REAL(KIND=LDP) HEAT(ND)           !Full heating term
	REAL(KIND=LDP) T_MIN
!
	INTEGER GET_DIAG
	INTEGER L,K
	REAL(KIND=LDP) TA(ND)             !Temporary work vector
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) EHB_CONSTANT
!
        GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_INDX)+DIAG_INDX
!
! This function computes the index L on BA( , ,?,K) corresponding
! to the local depth variable (i.e that at K). It is equivalent
! to IF(NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
!
! The expression out the front is a crude expression for the
! free-free cooling rate. Remember that it has been multiplied by
! 10^10 because CHI is 10^10 larger to keep r.chi conserved.
! We should also divide the rate by 4Pi, but we just want a term
! similar in magnitude to what appears in the RE Eqn.
!
	HEAT(1:ND)=0.0_LDP
	T_MIN_EXTRAP=.FALSE.
	EHB_CONSTANT=1.25663706E-09_LDP		!4P x E-10 (for electron heating equation).
!
	IF(T_MIN .GT. 0.0_LDP)THEN
!
	  DO K=1,ND
	    TA(K) = 1.0E-10_LDP*ED(K)*POPION(K)
	    HEAT(K)= TA(K)*EXP(1000.0_LDP*(T_MIN-T(K)))
	  END DO
!
! If the added heating term is much greater than what is needed, we
! take into account the exponential nature of the added term so that
! we can get to the solution quicker.
!
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	     DO K=1,ND
	       WRITE(133,*)TA(K),HEAT(K),STEQ_T(K)
	       L=GET_DIAG(K)
	       IF(HEAT(K) .GT. -100.0_LDP*STEQ_T(K) .AND. STEQ_T(K) .LT. 0.0_LDP)THEN
	         T1=T_MIN-T(K)-LOG(ABS(STEQ_T(K))/TA(K))/1000.0_LDP
	         BA_T(NT,L,K)=BA_T(NT,L,K) - (HEAT(K)-STEQ_T(K))/T1
	         BA_T_EHB(NT,L,K)=BA_T_EHB(NT,L,K) - (EHB_CONSTANT*HEAT(K)-STEQ_T_EHB(K))/T1
	         T_MIN_EXTRAP=.TRUE.
	       ELSE
	         BA_T(NT,L,K)=BA_T(NT,L,K) - 1000.0_LDP*HEAT(K)
	         BA_T_EHB(NT,L,K)=BA_T_EHB(NT,L,K) - 1000.0_LDP*EHB_CONSTANT*HEAT(K)
	       END IF
	       BA_T(NT-1,L,K)=BA_T(NT-1,L,K) + 2.0_LDP*HEAT(K)/ED(K)
	       BA_T_EHB(NT-1,L,K)=BA_T_EHB(NT-1,L,K) + 2.0_LDP*EHB_CONSTANT*HEAT(K)/ED(K)
	     END DO
	  END IF
!
! Add in extra heating to radiative equilibrium equation.
!
	  DO K=1,ND
	    STEQ_T(K)=STEQ_T(K) + HEAT(K)
	  END DO
	  HEAT=HEAT*EHB_CONSTANT
	  STEQ_T_EHB=STEQ_T_EHB+HEAT
!
        END IF
!
	RETURN
	END
