	SUBROUTINE CHK_COOL_RATE_CONSIST(AVE_ENERGY, STEQ_T_NO_LINES, NU, ML, NCF, ND, LU)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NCF
	INTEGER ML
	INTEGER LU
	REAL(KIND=LDP) NU(NCF)
	REAL(KIND=LDP) AVE_ENERGY(ND)
	REAL(KIND=LDP) STEQ_T_NO_LINES(ND)
!
! Local data
!
	INTEGER I,J,L,ID
	REAL(KIND=LDP) PC,T1
	REAL(KIND=LDP) STEQ_SUM(ND)
	REAL(KIND=LDP) COOL_EB_EST(ND)
	REAL(KIND=LDP) COOL_EB_EST2(ND)
!
	EXTERNAL PLANCKS_CONSTANT
	REAL(KIND=LDP) PLANCKS_CONSTANT
!
	PC=1.0D+15*PLANCKS_CONSTANT()
	STEQ_SUM=0.0D0
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO L=1,ND
	      DO I=1,ATM(ID)%NXzV
	        J=ATM(ID)%EQXzV+I-1
	        STEQ_SUM(L)=STEQ_SUM(L)+AVE_ENERGY(J)*SE(ID)%STEQ(I,L)
	      END DO
	    END DO
	  END IF
	END DO
!
	T1=1.0D-10*16*ATAN(1.0D0)
	DO L=1,ND
	  COOL_EB_EST(L)=T1*STEQ_T_NO_LINES(L)+PC*STEQ_SUM(L)
	  COOL_EB_EST2(L)=T1*STEQ_T(L)+PC*STEQ_SUM(L)
	END DO
!
	WRITE(LU,'(I8,3X,ES16.8,6ES17.7)')ML, NU(ML), STEQ_T_EHB(5), COOL_EB_EST(5), COOL_EB_EST2(5),
	1             PC*STEQ_SUM(5), T1*STEQ_T_NO_LINES(5), T1*STEQ_T(5)
!
	RETURN
	END
