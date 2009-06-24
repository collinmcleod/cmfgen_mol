!
! File to check the consistency of parameters for CMF_FLUX.
!
! Created: 16-Jun-2009
!
	SUBROUTINE CHECK_CMF_FLUX_PARAM_CONSIS()
	USE CMF_FLUX_CNTRL_VAR_MOD
	IMPLICIT NONE
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	CHARACTER(LEN=10) ER_LAB
!
	ER_LAB='Warning:'
	IF(STOP_IF_BAD_PARAM)ER_LAB='ERROR:'
	LUER=ERROR_LU()
!
! Testing parameters for SN model.
!	
	IF(SN_MODEL)THEN
	  IF(DO_REL_IN_OBSFRAME .NE..TRUE.)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_RELO should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(DO_CMF_REL_OBS .NE..TRUE.)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_CMF_RELO should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(USE_J_REL .AND. INCL_REL_TERMS .NE..TRUE.)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_REL should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(USE_J_REL .AND. INCL_ADVEC_TERMS_IN_TRANS_EQ .NE..TRUE.)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_ADV_TRANS should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(USE_J_REL .AND. .NOT. USE_FORMAL_REL)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)TRIM(ER_LAB),' UE_FRM_REL should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
	  IF(PLANE_PARALLEL  .OR. PLANE_PARALLEL_NO_V)THEN
	    WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
	    WRITE(LUER,*)'PP_NOV and PP_MOD cannot be TRUE for a SN model'
	    STOP
	  END IF
	END IF
!
! Testing parameters for plane-parallel models.
!
	IF(PLANE_PARALLEL  .AND. PLANE_PARALLEL_NO_V)THEN
	  WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
	  WRITE(LUER,*)'PP_NOV and PP_MOD cannot both be TRUE'
	  WRITE(LUER,*)'Use PP_NOV for a plane-parallel model with no velocity field'
	  WRITE(LUER,*)'Use PP_MOD for a plane-parallel model with a velocity field'
	  STOP
	END IF
!
! Testing gneric parameters, valid for all models.
!
	IF(OBS_TAU_MAX .LT. 10.0)THEN
	  WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
	  WRITE(LUER,*)'OBS_TAU_MAX should be > 10.0'
	  IF(STOP_IF_BAD_PARAM)STOP
	END IF
	IF(OBS_ES_DTAU .GT. 0.2)THEN
	  WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
	  WRITE(LUER,*)'OBS_ES_DTAU should be < 0.2 '
	  IF(STOP_IF_BAD_PARAM)STOP
	END IF
!
! These parameters control the insertion of extra grid points into the grid before the
! ray/moment routines are called.
!
	IF(ACCURATE .AND. .NOT. ALL_FREQ)THEN
	  WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
          WRITE(LUER,*)'It is recommended that ALL_FREQ be set to TURE when using the INC_GIRD  option'
	  IF(STOP_IF_BAD_PARAM)STOP
	END IF
	IF(ACCURATE .AND. (NPINS .LT .1 .OR. NPINS .GT. 5))THEN
	  WRITE(LUER,*)'Error in control parameters in CMF_FLUX_PARAM_INIT'
          WRITE(LUER,*)'N_INS should be 1, 2, 3, 4 or 5'
	  IF(STOP_IF_BAD_PARAM)STOP
	END IF
!
	RETURN
	END
