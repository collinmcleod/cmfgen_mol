	SUBROUTINE GET_RAD_DECAY_ENERGY(E_RAD_DECAY,ND)
	USE NUC_ISO_MOD
	USE SHOCK_POWER_MOD 
	IMPLICIT NONE
!
! Altered: 17-Nov-2022: Check that SHOCK_POWER is allocated.
! Altered: 14-Aug-2022: Include SHOCK_POWER -- based on LUC.
!
	INTEGER ND
	REAL(10) E_RAD_DECAY(ND)
!
! No need for IF statement at present, since SHOCK_POWER will be zero
!          if unimportant.
!
!	E_RAD_DECAY(1:ND)=RADIOACTIVE_DECAY_ENERGY(1:ND)
	IF(ALLOCATED(SHOCK_POWER))THEN
	  E_RAD_DECAY(1:ND)=RADIOACTIVE_DECAY_ENERGY(1:ND) + SHOCK_POWER(1:ND)
	ELSE
	  E_RAD_DECAY(1:ND)=RADIOACTIVE_DECAY_ENERGY(1:ND)
	END IF
!
	RETURN
	END
