!
! Subroutine to read in a model containing the populations from a
! previously converged model at an earlier time step.
!
! NB: TIME_SEQ_NO refers to the current model. This routine reads the populations
! corresponding to model TIME_SEQ_NO-1
!
        SUBROUTINE GET_POPS_AT_PREV_TIME_STEP(R,V,OLD_POPS,OLD_R,TIME_SEQ_NO,ND,NT,LU)
	IMPLICIT NONE
!
! Created : 12-Dec-2005
!
	INTEGER ND
	INTEGER NT
	INTEGER LU
	INTEGER TIME_SEQ_NO
!
	REAL*8 R(ND)
	REAL*8 V(ND)
!
	REAL*8 OLD_R(ND)
	REAL*8 OLD_SIGMA(ND)
	REAL*8 OLD_POPS(NT,ND)
!
	REAL*8 OLD_V(ND)
	REAL*8 LOG_OLD_V(ND)
	REAL*8 LOG_V(ND)
	REAL*8 NEW_VEC(ND)
	REAL*8 OLD_VEC(ND)
!
	INTEGER IREC_RD
	INTEGER I
	INTEGER, PARAMETER :: IONE=1
	LOGICAL, PARAMETER :: RVSIG_WRITTEN=.TRUE.
!
! Get model from the last time step. TIME_SEQ_NO refers to the CURRENT 
! time model. Therefore we must subtract 1.
!
	IREC_RD=TIME_SEQ_NO-1
	CALL  READ_TIME_MODEL_V1(OLD_R,OLD_V,OLD_SIGMA,OLD_POPS,
	1             IREC_RD,RVSIG_WRITTEN,NT,ND,LU)
!
! As a Hubble law, we can use V to interpolate. Note that 
! V is a comoving variable.
!
	OLD_V(1)=V(1); OLD_V(ND)=V(ND)
	LOG_V=LOG(V)
	LOG_OLD_V=LOG(OLD_V)
	DO I=1,NT
	  OLD_VEC=LOG(OLD_POPS(I,:))
	  CALL MON_INTERP(NEW_VEC,ND,IONE,LOG_V,ND,OLD_VEC,ND,LOG_OLD_V,ND)
	  OLD_POPS(I,:)=EXP(NEW_VEC)
	END DO
!
! Now get the interpolated radius scale in the old model.
! This radius scale will have the same velocity coordinates as the
! current model.
!
	OLD_VEC=LOG(OLD_R)
	CALL MON_INTERP(NEW_VEC,ND,IONE,LOG_V,ND,OLD_VEC,ND,LOG_OLD_V,ND)
	OLD_R(2:ND-1)=EXP(NEW_VEC(2:ND-1))
!
	RETURN
	END
