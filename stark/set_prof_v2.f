!
! Routine to compute:
!
!          (1) Doppler line profiles for any species
!          (2) Approximate STARK profiles for any Hydrogenic species of
!                charge Z. The theory is excellent for high Balmer lines, but
!                  only approximate for Ha.
!
! Uses the GRIEM theory as modified by Auer and Mihalas AP J S 24 1972.
!
! Output:
!        PROF - Profile as a function of depth. Either:
!                    (a) Doppler profile or
!                    (b) STARK profile convolved with a Doppler profile.
!
!                    The Doppler profile can have a turbulent contribution
!                    given by VTURB.
! Input:
!        NU        - Velocity in units of C from line center [= (V-Vo)/Vo ].
!        ED_IN     - Electron density (/cm^3) (Vector, length ND)
!        TEMP_IN   - Temperature (10^4 K)     (Vector, length ND)
!        AMASS_IN  - Atomic mass for Doppler profile in amu's.
!        Z_IN      - Ion charge.
!        NL        - Lower level of transition.
!        NUP       - Upper level of transition.
!        ND        - Nuber of (Ne,T) values profile is to be computed for.
!        VTURB_IN   - Turbulent velocity in km/s (function of depth).
!
	SUBROUTINE SET_PROF_V2(PROF,NU,ML_CUR,ML_ST,ML_END,
	1             ED_IN,TEMP_IN,VTURB_IN,ND,
	1             PROF_TYPE,NU_ZERO,NL,NUP,AMASS_IN,Z_IN,
	1             GAM_NAT,GAM_COL,VTURB_FIX,
	1             END_RES_ZONE,NORM_PROFILE)
	USE SET_KIND_MODULE
!
! Altered 07-Jan-1999: ML_CUR now passed in call.
!                      Profile is now recomputed whenever ML_CUR=ML_ST,
!                      or when profile is unavailable.
!
	USE PROF_MOD
	IMPLICIT NONE
	INTEGER ML_CUR
	INTEGER ML_ST
	INTEGER ML_END
	INTEGER ND
	INTEGER NL
	INTEGER NUP
	REAL(KIND=LDP) PROF(ND)
	REAL(KIND=LDP) NU(ML_END)  		!Can actually be larger
	REAL(KIND=LDP) ED_IN(ND)
	REAL(KIND=LDP) TEMP_IN(ND)
	REAL(KIND=LDP) VTURB_IN(ND)
	REAL(KIND=LDP) AMASS_IN
	REAL(KIND=LDP) Z_IN
	REAL(KIND=LDP) NU_ZERO
	REAL(KIND=LDP) GAM_NAT
	REAL(KIND=LDP) GAM_COL
	REAL(KIND=LDP) VTURB_FIX		!Turbulent velocity (km/s): same at all depths
	CHARACTER*(*) PROF_TYPE
	LOGICAL END_RES_ZONE
	LOGICAL NORM_PROFILE
!
	LOGICAL, PARAMETER :: RET_LOG=.FALSE.
!
! Local variables
!
	INTEGER I,J,ML
	INTEGER LOC_INDX		!Indicates which storage
	INTEGER NF_GR
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TMP_ED,NU_DOP
	REAL(KIND=LDP) TMP_VEC(ND)
	REAL(KIND=LDP) A_VOIGT
	REAL(KIND=LDP) V_VOIGT
!
! External functions
!
	REAL(KIND=LDP) VOIGT
!
! Doppler profile is the same for all species, and is the same at all depths.
!
	IF(PROF_TYPE .EQ. 'DOP_FIX')THEN
          T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
          NU_DOP=NU_ZERO*VTURB_FIX/C_KMS
          DO I=1,ND
            PROF(I)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	  END DO
	  RETURN
!
! Doppler profile varies with species and depth.
!
	ELSE IF(PROF_TYPE .EQ. 'DOPPLER')THEN
          T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
          TMP_VEC(1:ND)=( VTURB_IN(1:ND)/12.85_LDP  )**2
          T2=NU_ZERO*12.85_LDP/C_KMS
          DO I=1,ND
            NU_DOP=T2*SQRT( TEMP_IN(I)/AMASS_IN + TMP_VEC(I) )
            PROF(I)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	  END DO
	  RETURN
!
	ELSE IF(PROF_TYPE .EQ. 'VOIGT')THEN
          TMP_VEC(1:ND)=( VTURB_IN(1:ND)/12.85_LDP  )**2
          T2=NU_ZERO*12.85_LDP/C_KMS
	  DO I=1,ND
            NU_DOP=T2*SQRT( TEMP_IN(I)/AMASS_IN + TMP_VEC(I) )
	    A_VOIGT=0.25E-15_LDP*GAM_NAT/PI/NU_DOP
	    V_VOIGT=(NU(ML_CUR)-NU_ZERO)/NU_DOP
	    PROF(I)=1.0E-15_LDP*VOIGT(A_VOIGT,V_VOIGT)/NU_DOP
	  END DO
	  RETURN
!
!
!
! In the case of a STARK profile it is computationally prohibitive to
!   compute the STARK profile for each frequency. Rather we compute
!   the STARK profile [f(Ne,T,Vturb) ] once, and save the data for
!   subsequent calls.
!
	ELSE IF(PROF_TYPE .EQ. 'HZ_STARK')THEN
!
! Check validity of passed parameters for HI and HeII lines.
!
	  IF(NUP .LE. NL)THEN
	    WRITE(LUER,*)'Error in SET_PROF: NUP < NL'
	    STOP
	  END IF
!
! Check if we have already computed profile at the desired frequency set.
! If so we use the data already computed. If ML_CUR .EQ. ML_ST, we assume
! that the full STARK profile needs to be recompyted.
!
	  DO LOC_INDX=1,NSTORE_MAX
	    IF( .NOT. STORE_AVAIL(LOC_INDX) )THEN
	      IF( NU_ZERO .EQ. NU_ZERO_STORE(LOC_INDX)    .AND.
	1           AMASS_IN .EQ. AMASS_STORE(LOC_INDX)     .AND.
	1           NL .EQ. NL_STORE(LOC_INDX)              .AND.
	1           NUP .EQ. NUP_STORE(LOC_INDX)             )THEN
!
	        IF(ML_CUR .EQ. ML_ST)THEN
!
! As first frequency, need to recompute profile.
!
	          STORE_AVAIL(LOC_INDX)=.TRUE.
	          GOTO 1000
	        ELSE
!
! Check store correct, then set profile data
!
	          LST_FREQ_LOC(LOC_INDX)=LST_FREQ_LOC(LOC_INDX)+1
	          IF(NU_STORE(LST_FREQ_LOC(LOC_INDX), LOC_INDX) .NE. NU(ML_CUR))THEN
	            WRITE(LUER,*)'Error in SET_PROF'
	            WRITE(LUER,*)'Frequencies don''t match'
	            WRITE(LUER,*)NU(ML_CUR),NU_STORE(LST_FREQ_LOC(LOC_INDX),LOC_INDX)
	            STOP
	          END IF
	          PROF(1:ND)=PROF_STORE(1:ND,LST_FREQ_LOC(LOC_INDX),LOC_INDX)
!
! If possible, free up STORE:
!
	          IF( END_RES_ZONE )THEN
	            STORE_AVAIL(LOC_INDX)=.TRUE.
	          END IF
	          RETURN
	        END IF
	      END IF
	    END IF
	  END DO
!
! If we get to here, we need to compute the Stark Profile. First we
! determine the location where the STARK profile can be stored.
!
1000	  CONTINUE
	  LOC_INDX=0
	  DO ML=1,NSTORE_MAX
	    IF( STORE_AVAIL(ML) )THEN
	      LOC_INDX=ML
	      EXIT
	    END IF
	  END DO
	  IF(LOC_INDX .EQ. 0)THEN
	    WRITE(LUER,*)'Error in SET_PROF --- insufficient storage'
	    WRITE(LUER,*)'Maximum number of storage locations is:',NSTORE_MAX
	    STOP
	  END IF
	  STORE_AVAIL(LOC_INDX)=.FALSE.
	  PROF_STORE(:,:,LOC_INDX)=0.0_LDP
!
! Determine how many frequencies we need to compute profile for, and check
! we have sufficent storage.
!
	  IF(ML_END-ML_ST+1 .GT. NFREQ_MAX)THEN
	    WRITE(LUER,*)'Error in SET_PROF --- NFREQ_MAX too small'
	    WRITE(LUER,*)'NFREQ_MAX=',NFREQ_MAX
	    WRITE(LUER,*)'Required NFREQ_MAX is',ML_END-ML_ST+1
	    STOP
	  END IF
!
! Convert from Frequency to Angstrom space, measured from line center.
!
	  DO ML=ML_ST,ML_END
	    DWS_GRIEM(ML-ML_ST+1)=0.01_LDP*C_KMS*(1.0_LDP/NU(ML)-1.0_LDP/NU_ZERO)
	  END DO
!
! Now compute, and store the Doppler broadened Stark profile for each
! depth. We assume that the profile is zero outside the computational range.
!
	  NF_GR=ML_END-ML_ST+1
	  IF(NF_GR .GT. NFREQ_MAX)THEN
	    WRITE(LUER,*)'Error in SET_PROF --- insufficient storage'
	    WRITE(LUER,*)'Max # of frequency storage locations is:',NFREQ_MAX
	    WRITE(LUER,*)'# of storage locatins requires is:',NF_GR
	    STOP
	  END IF
	  PROF_STORE(1:ND,1:NF_GR,LOC_INDX)=0.0_LDP
	  WRITE(77,*)NL,NUP,Z_IN
	  DO I=1,ND
	    TMP_ED=1.0E+16_LDP
	    TMP_ED=MIN(ED_IN(I),TMP_ED)
            CALL GRIEM_V2(PR_GRIEM,DWS_GRIEM,NF_GR,
	1        TMP_ED,TEMP_IN(I),VTURB_IN(I),
	1        NL,NUP,Z_IN,AMASS_IN,RET_LOG)
!
! Normalize profile to an integral of unity assuming a trapazoidal integral.
! Since NU is units of 10^15 Hz, we need to multiply the integral by
! 10^15 before scaling.
!
	    IF(NORM_PROFILE)THEN
	      T1=0.0_LDP
	      DO ML=1,NF_GR-1
	        T1=T1+0.5_LDP*(NU(ML_ST+ML-1)-NU(ML_ST+ML))*
	1                       (PR_GRIEM(ML)+PR_GRIEM(ML+1))
	      END DO
	      T1=ABS(T1)*1.0E+15_LDP
	      PROF_STORE(I,1:NF_GR,LOC_INDX)=PR_GRIEM(1:NF_GR)/T1
	      WRITE(77,*)I,T1
	    ELSE
	      PROF_STORE(I,1:NF_GR,LOC_INDX)=PR_GRIEM(1:NF_GR)
	    END IF
	  END DO
!
! Remeber to return profile.
!
	  PROF(1:ND)=PROF_STORE(1:ND,1,LOC_INDX)
!
! Save information to describe profile so that we don't have to compute it
! all over again.
!
	  NU_ZERO_STORE(LOC_INDX)=NU_ZERO
	  AMASS_STORE(LOC_INDX)=AMASS_IN
	  NL_STORE(LOC_INDX)=NL
	  NUP_STORE(LOC_INDX)=NUP
	  LST_FREQ_LOC(LOC_INDX)=1
	  NF_STORE(LOC_INDX)=NF_GR
	  NU_STORE(1:NF_GR,LOC_INDX)=NU(ML_ST:ML_END)
!
	END IF			!Profile type
!
!	I=ND-7
!	WRITE(99,*)AMASS_IN,NL,NUP
!	WRITE(99,*)ED_IN(I),TEMP_IN(I),VTURB_IN(I)
!	WRITE(99,'(1P,5E15.5)')0.01D0*C_KMS/NU_STORE(1:NF_GR,LOC_INDX)
!	WRITE(99,'(1P,5E15.5)')PROF_STORE(I,1:NF_GR,LOC_INDX)
!
	RETURN
	END
