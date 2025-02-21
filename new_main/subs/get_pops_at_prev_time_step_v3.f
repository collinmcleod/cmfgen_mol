!
! Subroutine to read in a model containing the populations from a
! previously converged model at an earlier time step. The new and current
! models must have an identical number of grid points, and super levels.
! Further, the velocity at the inner boundary must be identical. The older
! model may extend to larger velocities.
!
! NB: TIME_SEQ_NO refers to the current model. This routine reads the populations
! corresponding to model TIME_SEQ_NO-1
!
        SUBROUTINE GET_POPS_AT_PREV_TIME_STEP_V3(OLD_POPS,OLD_R,DO_ADVECT,DO_RAD_DECAYS,
	1                      NORMALIZE_POPS,TIME_SEQ_NO,ND,NT,LU)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered : 01-Jul-2006 : Now allow outer velocity of older model to be larger than current
!                           model, although number of grid points still must be identical.
! Altered : 23-Jun-2006 : Removed R, V from call as in MOD_CMFGEN
! Created : 12-Dec-2005
!
	INTEGER ND
	INTEGER NT
	INTEGER LU
	INTEGER TIME_SEQ_NO
	LOGICAL DO_ADVECT
	LOGICAL DO_RAD_DECAYS
	LOGICAL NORMALIZE_POPS
!
	REAL(KIND=LDP) OLD_R(ND)
	REAL(KIND=LDP) OLD_POPS(NT,ND)
!
	REAL(KIND=LDP) OLD_V(ND)
	REAL(KIND=LDP) OLD_SIGMA(ND)
	REAL(KIND=LDP) LOG_OLD_V(ND)
	REAL(KIND=LDP) LOG_V(ND)
	REAL(KIND=LDP) NEW_VEC(ND)
	REAL(KIND=LDP) OLD_VEC(ND)
	REAL(KIND=LDP) OLD_ED(ND)
!
	REAL(KIND=LDP) T1,T2
	INTEGER IREC_RD
	INTEGER I,J,K,ISPEC,ID
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU,EQUAL
	LOGICAL EQUAL
!
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
	T1=1.0E-06_LDP
	IF(EQUAL(OLD_V(ND),V(ND),T1))THEN
	  OLD_V(ND)=V(ND)
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in GET_POPS_AT_PREV_TIME_STEP_V3'
	  WRITE(LUER,*)'Velocities at inner boundary are unequal'
	  WRITE(LUER,*)'V(ND)=',V(ND)
	  WRITE(LUER,*)'OLD_V(ND)=',OLD_V(ND)
	  WRITE(LUER,*)'V(ND)/OLD_V(ND)=',V(ND)/OLD_V(ND)
	  WRITE(LUER,*)'T1=',T1
	  STOP
	END IF
	IF(EQUAL(OLD_V(1),V(1),T1))THEN
	  OLD_V(1)=V(1)
	ELSE IF(OLD_V(1) .LT. V(1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in GET_POPS_AT_PREV_TIME_STEP_V3'
	  WRITE(LUER,*)'Old velocity at outer boundary is too small'
	  WRITE(LUER,*)'V(1)=',V(1)
	  WRITE(LUER,*)'OLD_V(1)=',OLD_V(1)
	  WRITE(LUER,*)'V(1)/OLD_V(1)=',V(1)/OLD_V(1)
	  WRITE(LUER,*)'T1=',T1
	END IF
!
! Perform interpolations in Log-Log plane, with V as the independent variable.
! Since V is proportional to r, this is equivalent to assuming r is the
! independent variable.
!
	LOG_V=LOG(V)
	LOG_OLD_V=LOG(OLD_V)
	DO I=1,NT
	  OLD_VEC=LOG(OLD_POPS(I,:))
	  CALL MON_INTERP(NEW_VEC,ND,IONE,LOG_V,ND,OLD_VEC,ND,LOG_OLD_V,ND)
	  OLD_POPS(I,:)=EXP(NEW_VEC)
	END DO
!
! Get the interpolated radius scale in the old model. This radius scale will
! have the same velocity coordinates as the current model. Since we have a
! Hubble law, linear interpolation is accurate.
!
	CALL MON_INTERP(NEW_VEC,ND,IONE,V,ND,OLD_R,ND,OLD_V,ND)
	OLD_R(1:ND-1)=NEW_VEC(1:ND-1)
!
!
! Advect the populations assumubg a pure Hubble Law. This will need to be changed
! when using more sophisticated expansion laws. A change will also need to be made
! to the nomalization section.
!
	IF(DO_ADVECT)THEN
	  T1=(OLD_R(ND)/R(ND))**3
	  OLD_POPS(1:NT-1,:)=OLD_POPS(1:NT-1,:)*T1
	END IF
!
! Adjust the populations for radioactive decays.
!
	IF(DO_RAD_DECAYS)THEN
	  CALL DO_LEV_POP_DECAYS(OLD_POPS,ND,NT)
	END IF
!
! Normalize the populations to ensure continuity equation is satisfied.
! This removes any slight variations introduced by the non-monotonic
! interpolation. It should not be done if radioactive decayse are important,
! and DO_RAD_DECAYS is false.
!
	IF(NORMALIZE_POPS)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO J=1,ND
	      T2=0.0_LDP
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
	        IF(DO_ADVECT)THEN
	          T2=POP_SPECIES(J,ISPEC)/T2
	        ELSE
	          T2=POP_SPECIES(J,ISPEC)/T2*(R(ND)/OLD_R(ND))**3
	        END IF
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
	END IF
!
! Now determine the electron density. Since we determine Ne from the scaled
! populations, no scaling is necessary.
!
	OLD_ED=0.0_LDP
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        OLD_ED(J)=OLD_ED(J)+(ATM(ID)%ZXzV-1.0_LDP)*OLD_POPS(K,J)
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
	RETURN
	END
