!
! Subroutine to read in a model containing the populations from a
! previously converged model at an earlier time step.
!
! NB: TIME_SEQ_NO refers to the current model. This routine reads the populations
! corresponding to model TIME_SEQ_NO-1
!
        SUBROUTINE GET_POPS_AT_PREV_TIME_STEP_V2(OLD_POPS,OLD_R,TIME_SEQ_NO,ND,NT,LU)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered : 23-Jun-2006 : Removed R, V from call as in MOD_CMFGEN
! Created : 12-Dec-2005
!
	INTEGER ND
	INTEGER NT
	INTEGER LU
	INTEGER TIME_SEQ_NO
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
	REAL*8 OLD_ED(ND)
!
	REAL*8 T1,T2
	INTEGER IREC_RD
	INTEGER I,J,K,ISPEC,ID
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
! Normalize the populations to ensure continuity equation is satisfied.
!
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    T2=0.0D0
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        T2=T2+OLD_POPS(K,J)
	      END DO
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      T2=T2+OLD_POPS(K,J)
!
! Can now do the normalization. Only valid for a Hubble flow.
!
	      T2=POP_SPECIES(J,ISPEC)/T2*(R(ND)/OLD_R(ND))**3
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        DO I=1,ATM(ID)%NXzV
	          K=ATM(ID)%EQXzV+I-1
	          OLD_POPS(K,J)=OLD_POPS(K,J)*T2
	        END DO
	      END DO
	      ID=SPECIES_END_ID(ISPEC)-1
	      K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      OLD_POPS(K,J)=OLD_POPS(K,J)*T2
	    END IF
	  END DO
	END DO
!
! Now determine the electron density. Since we determine Ne from the scaled
! populatons, no scaling is necessary.
!
	OLD_ED=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        OLD_ED(J)=OLD_ED(J)+(ATM(ID)%ZXzV-1.0D0)*OLD_POPS(K,J)
	      END DO
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      OLD_ED(J)=OLD_ED(J)+ATM(ID)%ZXzV*OLD_POPS(K,J)
	    END IF
	  END DO
	END DO
	DO J=1,ND
	  OLD_POPS(NT-1,J)=OLD_ED(J)
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
