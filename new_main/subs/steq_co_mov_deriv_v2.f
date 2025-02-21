!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! the comoving D/Dt term. Presently designed for SN with a Hubble-like
! velocity law.
!
	SUBROUTINE STEQ_CO_MOV_DERIV_V2(RELAXATION_PARAMETER,LINEAR,
	1             INCL_DT_TERM,LAMBDA_ITERATION,COMPUTE_BA,
	1             TIME_SEQ_NO,NUM_BNDS,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 12-Feb-2008 : Changed to V2: HUBBLE_LAW omitted from call.
!                       Some cleaning.
! Altered 22-Mar-2007 : Call GET_POPS_AT_PREV_TIME_STEP_V4
! Altered 23-Jun-2006 : R, V removed from call to GET_POPS_AT_PREV_TIME_STEP_V2 since
!                         passed by module MOD_CMFGEN.
! Created 12-Dec-2005: Based on STEQ_ADVEC_V4
!
	REAL(KIND=LDP) RELAXATION_PARAMETER
!
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER NT
	INTEGER TIME_SEQ_NO
!
	LOGICAL LAMBDA_ITERATION
	LOGICAL COMPUTE_BA
	LOGICAL LINEAR
	LOGICAL INCL_DT_TERM
!
! Local variables.
!
	REAL(KIND=LDP) OLD_POPS(NT,ND)
	REAL(KIND=LDP) OLD_R(ND)
!
	REAL(KIND=LDP) SUM(NUM_IONS,ND)
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DERIV_CONST
	REAL(KIND=LDP) DELTA_TIME_SECS			!Time step
!
	INTEGER K			!Depth index
	INTEGER M			!Band index
	INTEGER I			!Variable index
	INTEGER J			!Variable index
	INTEGER LU
	INTEGER ION_IVAR
	INTEGER IVAR
	INTEGER ID
	INTEGER ID_STRT,ID_END
	INTEGER ISPEC
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	SUM(:,:)=0.0_LDP
	DO ID=1,NUM_IONS
	  SE(ID)%STEQ_ADV(:)=0.0_LDP
	  SE(ID)%STRT_ADV_ID(:)=0.0_LDP
	  SE(ID)%END_ADV_ID(:)=0.0_LDP
	END DO
	BA_ADV_TERM(:,:)=0.0_LDP
	IF(.NOT. INCL_DT_TERM)RETURN
!
	IF(.NOT. LINEAR)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in STEQ_CO_MOV_DERIV_V2'
	  WRITE(LUER,*)'Only linear option is currently installed.'
	  LINEAR=.TRUE.
	END IF
!
! The three L_TRUE indicate that we correct OLD_POPS for both advection (i.e., the expansions
! of the SN), radioactive decay, and that we normalize the population of each species so
! that the continuity equation is exactly satisfied.
!
	LU=7
	CALL GET_POPS_AT_PREV_TIME_STEP_V4(OLD_POPS,OLD_R,
	1         L_TRUE,L_TRUE,L_TRUE,TIME_SEQ_NO,ND,NT,LU)
!
! The relaxation factor should be < 1, and is used to adjust the importance of the
! advection terms. It should be 1 for the final model. It should only be used to
! help converge a model in which advection terms are very important.
!
	DELTA_TIME_SECS=1.0E+05_LDP*(R(ND)-OLD_R(ND))/V(ND)
	DERIV_CONST=RELAXATION_PARAMETER/DELTA_TIME_SECS
!
! We use backward linear differencing. Since the OLD_POPS were corrected
! for expansion when they were read in, this section is valid for all
! expansion laws.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    IVAR=ATM(ID)%EQXzV
	    ION_IVAR=IVAR+ATM(ID)%NXzV
	    DO K=1,ND
	      DO I=1,ATM(ID)%NXzV                                !Which S.E. equation
	        T1=DERIV_CONST*(ATM(ID)%XzV(I,K)-OLD_POPS(IVAR+I-1,K))
	        SE(ID)%STEQ(I,K)=SE(ID)%STEQ(I,K) - T1
	        SUM(ID,K)=SUM(ID,K) + T1
	        IF(K .EQ. 1)THEN
	          WRITE(131,'(I5,4ES14.6)')I,ATM(ID)%XzV(I,K),OLD_POPS(IVAR+I-1,K),T1,SUM(ID,K)
	        END IF
	      END DO
	      IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
	        SUM(ID+1,K)=SUM(ID+1,K)+DERIV_CONST*(ATM(ID)%DXzV(K)-OLD_POPS(ION_IVAR,K))
	        IF(K .EQ. 1)THEN
	          T1=DERIV_CONST*(ATM(ID)%DXzV(K)-OLD_POPS(ION_IVAR,K))
	          WRITE(131,'(I5,4ES14.6)')I,ATM(ID)%DXzV(K),OLD_POPS(ION_IVAR,K),T1,SUM(ID,K)
	        END IF
	      END IF
	    END DO
	  END DO
	END DO
!
! We now have to compute the advection terms for the ion equations. This is complicated
! because of stability issues. Want to use the smallest terms possible.
! We have this option since Sum[all j] Xj = X (X=species population).
!
! Let Xi refer to the total population of ionization stage i. Then
!
!      Sum[j=1,i]      {vdXj/dr} = RRi  or
!
!      Sum[j=i+1,..]  {-vdXj/dr} = RRi
!
! where RRI refers to the net recombination (phot+col) to ionization stage i.
! We choose the equation for which Sum Xj is smallest.
!
! Since we have corrected old populations for radioactive decays, this works
! even in presence of radioactive decays (which causes changes in species
! populations).
!
	DO ISPEC=1,NUM_SPECIES
	  ID_STRT=SPECIES_BEG_ID(ISPEC)
	  ID_END=SPECIES_END_ID(ISPEC)
	  DO K=1,ND
	    DO ID=ID_STRT,ID_END-1
	      T1=0.0_LDP
	      DO I=ID_STRT,ID
	        DO J=1,ATM(I)%NXzV
	          T1=T1+ATM(I)%XzV(J,K)
	        END DO
	      END DO
	      T2=0.0_LDP
	      DO I=ID+1,ID_END-1
	        DO J=1,ATM(I)%NXzV
	          T2=T2+ATM(I)%XzV(J,K)
	        END DO
	      END DO
	      T2=T2+ATM(ID_END-1)%DXzV(K)
	      IF(T1 .GT. T2)THEN
	        DO I=ID+1,ID_END
		  SE(ID)%STEQ_ADV(K)=SE(ID)%STEQ_ADV(K)+SUM(I,K)
	        END DO
	        SE(ID)%STRT_ADV_ID(K)=ID+1
	        SE(ID)%END_ADV_ID(K)=ID_END-1          !Don't count ion.
	      ELSE
	        DO I=ID_STRT,ID
                  SE(ID)%STEQ_ADV(K)=SE(ID)%STEQ_ADV(K)-SUM(I,K)
	        END DO
	        SE(ID)%STRT_ADV_ID(K)=ID_STRT
	        SE(ID)%END_ADV_ID(K)=ID
	      END IF
	    END DO
	  END DO
	END DO
!
! Only have diagonal terms.
!
	M=(NUM_BNDS/2)+1				!Diagonal index
        DO K=1,ND
          BA_ADV_TERM(M,K)=DERIV_CONST
	END DO
	IF(COMPUTE_BA)THEN
	  DO ISPEC=1,NUM_SPECIES
	    ID_STRT=SPECIES_BEG_ID(ISPEC)
	    ID_END=SPECIES_END_ID(ISPEC)
	    DO ID=ID_STRT,ID_END-1
	      DO K=1,ND
	        DO I=1,ATM(ID)%NXzV
	          SE(ID)%BA(I,I,M,K)=SE(ID)%BA(I,I,M,K)-DERIV_CONST
	        END DO		!loop over S.E. equation
	      END DO		!loop over depth
	    END DO		!loop over ionization stages
	  END DO		!loop over species
	END IF
!
	RETURN
	END
