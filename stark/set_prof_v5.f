!
! Routine to compute:
!
!          (1) Doppler line profiles for any species
!          (2) Stark and Voigt profiles (of varrying accuracy) for individual lines.
!          (3) Approximate STARK profiles for any Hydrogenic species of charge
!                Z. The theory is excellent for high Balmer lines, but only
!                approximate for Ha. Uses the GRIEM theory as modified by Auer
!                and Mihalas AP J S 24 1972.
!
!
! Output:
!        PROF - Profile as a function of depth. Either:
!                    (a) Doppler profile or
!                    (b) STARK profile convolved with a Doppler profile.
!
!                    The Doppler profile can have a turbulent contribution
!                    given by VTURB.
! Input:
!        NU         - Freqency (in units of 10^15 Hz)
!        ED_IN      - Electron density (/cm^3) (Vector, length ND)
!        POP_PROTON - Proton density (/cm^3) (Vector, length ND)
!	 POP_HEPLUS - He+ density (/cm^3) (Vector, length ND)
!        TEMP_IN    - Temperature (10^4 K)     (Vector, length ND)
!        AMASS_IN   - Atomic mass for Doppler profile in amu's.
!        Z_IN       - Ion charge.
!        NL         - Lower level of transition.
!        NUP        - Upper level of transition.
!        ND         - Nuber of (Ne,T) values profile is to be computed for.
!        VTURB_IN   - Turbulent velocity in km/s (function of depth).
!
	SUBROUTINE SET_PROF_V5(PROF,NU,ML_CUR,ML_ST,ML_END,
	1             ED_IN,POP_PROTON,POP_HEPLUS,
	1             TEMP_IN,VTURB_IN,ND,
	1             PROF_TYPE,PROF_LIST_LOCATION,
	1             NU_ZERO,NL,NUP,AMASS_IN,Z_IN,
	1             GAM_RAD,C4_INTER,
	1             TDOP,AMASS_DOP,VTURB,MAX_PROF_ED,
	1             END_RES_ZONE,NORM_PROFILE,LU_STK)
	USE SET_KIND_MODULE
!
! Aleteed 12-Nov-2023: Fixed bug when computing STARK profiles when ND=1
! Altered 05-Jul-2022: Call to GRIEM_V2 now in parallel loop. This paralleization
!                        required some changes to GRIEM_STARK_MOD.
! Altered 18-May-2015: LOC_GAM_COL is set to C4 unless specifically set in STRK_LIST.
! Altered 20-May-2014: VERBOSE option introduced, and sone diagnostic output was modified.
! Altered 14-May-2014: MAX_PROF_ED now applies to all profiles
! Altered 31-Jan-2014: Added MAX_PROF_ED to call (changed to V5).
! Altered 06-Jan-2014: VTURB_FIX replaced by TDOP,AMASS_DOP,VTURB in call (changed to V4).
!                          (taken from cur_cm_25jun13 development version).
! Altered 03-Jan-2001: Check for unrecognized profile type.
! Altered 07-Jan-1999: ML_CUR now passed in call.
!                      Profile is now recomputed whenever ML_CUR=ML_ST,
!                      or when profile is unavailable.
!
	USE PROF_MOD
	USE MOD_STRK_LIST
!
	IMPLICIT NONE
	INTEGER ML_CUR
	INTEGER ML_ST
	INTEGER ML_END
	INTEGER ND
	INTEGER NL
	INTEGER NUP
	INTEGER LU_STK
	REAL(KIND=LDP) PROF(ND)
	REAL(KIND=LDP) NU(ML_END)  		!Can actually be larger
	REAL(KIND=LDP) ED_IN(ND)
	REAL(KIND=LDP) POP_PROTON(ND)
	REAL(KIND=LDP) POP_HEPLUS(ND)
	REAL(KIND=LDP) TEMP_IN(ND)
	REAL(KIND=LDP) VTURB_IN(ND)
	REAL(KIND=LDP) AMASS_IN
	REAL(KIND=LDP) Z_IN
	REAL(KIND=LDP) NU_ZERO
	REAL(KIND=LDP) GAM_RAD
	REAL(KIND=LDP) C4_INTER
	REAL(KIND=LDP) VTURB		!Turbulent velocity (km/s): same at all depths
	REAL(KIND=LDP) TDOP
	REAL(KIND=LDP) AMASS_DOP
	REAL(KIND=LDP) MAX_PROF_ED
	INTEGER PROF_LIST_LOCATION
	CHARACTER*(*) PROF_TYPE
	LOGICAL END_RES_ZONE
	LOGICAL NORM_PROFILE
	LOGICAL VERBOSE
!
	LOGICAL, PARAMETER :: RET_LOG=.FALSE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Local variables
!
	INTEGER I,ML,ID,ITR
	INTEGER LOC_INDX		!Indicates which storage
	INTEGER NF
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TMP_ED,NU_DOP
	REAL(KIND=LDP) TMP_VEC(ND)
	REAL(KIND=LDP) ED_MOD(ND)
	REAL(KIND=LDP) XNU(ML_END-ML_ST+1)
	REAL(KIND=LDP) A_VOIGT
	REAL(KIND=LDP) V_VOIGT
	REAL(KIND=LDP) LOC_GAM_RAD
	REAL(KIND=LDP) LOC_GAM_COL
!
! External functions
!
	REAL(KIND=LDP) VOIGT
!
! Doppler profile is the same for all species, and is the same at all depths.
!
	IF(PROF_TYPE .EQ. 'DOP_FIX')THEN
          T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
	  T2=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
          NU_DOP=NU_ZERO*T2/C_KMS
          DO I=1,ND
            PROF(I)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	  END DO
	  RETURN
!
! Doppler profile is the same at all depth, but not the same for all species.
! A reasonable choice for T_DOP is Teff.
!
	ELSE IF(PROF_TYPE .EQ. 'DOP_SPEC')THEN
          T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
	  T2=12.85_LDP*SQRT( TDOP/AMASS_IN + (VTURB/12.85_LDP)**2 )
          NU_DOP=NU_ZERO*T2/C_KMS
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
	  LOC_GAM_RAD=GAM_RAD
	  LOC_GAM_COL=1.55E+04_LDP*(ABS(C4_INTER)**0.667_LDP)
	  ID=PROF_LIST_LOCATION
	  IF(ID .NE. 0)THEN
	    LOC_GAM_RAD=LST_GAM_RAD(ID)
	    LOC_GAM_COL=LST_GAM_COL(ID)
	  END IF
	  DO I=1,ND
	    TMP_ED=MIN(ED_IN(I),MAX_PROF_ED)
            NU_DOP=T2*SQRT( TEMP_IN(I)/AMASS_IN + TMP_VEC(I) )
	    A_VOIGT=0.25E-15_LDP*(LOC_GAM_RAD+LOC_GAM_COL*TMP_ED)/PI/NU_DOP
	    V_VOIGT=(NU(ML_CUR)-NU_ZERO)/NU_DOP
	    PROF(I)=1.0E-15_LDP*VOIGT(A_VOIGT,V_VOIGT)/NU_DOP
	  END DO
	  RETURN
!
!
!
! If we reach here we have two choices:
!
! (1) Profile has already been computed, hence we can use
!       the tabulated values.
! (2) We need to compute full profile, using a variety
!       of different methods.
!
	ELSE
!
! Check if profile data has already been computed. If ML_CUR .EQ. ML_ST,
! we assume that the full STARK profile needs to be recomputed. Thus
! ML_CUR=ML_ST is like a re-initilaization.
!
	  DO LOC_INDX=1,NSTORE_MAX
	    IF( .NOT. STORE_AVAIL(LOC_INDX) )THEN
	      IF( NU_ZERO .EQ. NU_ZERO_STORE(LOC_INDX)    .AND.
	1           NL .EQ. NL_STORE(LOC_INDX)              .AND.
	1           NUP .EQ. NUP_STORE(LOC_INDX)             )THEN
!
!If first frequency, need to recompute profile.
! We thuse exit to computation section.
!
	        IF(ML_CUR .EQ. ML_ST)THEN
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
	END IF
!
! 
!
1000	CONTINUE
!
! In this section of the routine, we compute general STARK profiles
! over the entire line profile, since it is generally computationally
! prohibitive to compute the STARK profile for each frequency. Rather
! we compute the STARK profile [f(Ne,T,Vturb) ] once, and save the data
! for subsequent calls.
!
! Determine the location where the STARK (i.e. Line) profile can be stored.
!
	  CALL GET_VERBOSE_INFO(VERBOSE)
	  DO I=1,ND
	   ED_MOD(I)=MIN(ED_IN(I),MAX_PROF_ED)
	  END DO
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
	    DO ML=1,NSTORE_MAX
	       WRITE(LUER,'(I4,ES14.5,2I8)')ML,NU_ZERO_STORE(ML),NL_STORE(ML),NUP_STORE(ML)
	    END DO
	    STOP
	  END IF
	  STORE_AVAIL(LOC_INDX)=.FALSE.
	  PROF_STORE(:,:,LOC_INDX)=0.0_LDP
!
! Determine how many frequencies we need to compute profile at, and check
! we have sufficent storage.
!
	  IF(ML_END-ML_ST+1 .GT. NFREQ_MAX)THEN
	    WRITE(LUER,*)'Error in SET_PROF --- NFREQ_MAX too small'
	    WRITE(LUER,*)'NFREQ_MAX=',NFREQ_MAX
	    WRITE(LUER,*)'Required NFREQ_MAX is',ML_END-ML_ST+1
	    STOP
	  END IF
!
! We can now compute the individual line profiles. New methods can
! be added to this section.
!
	NF=ML_END-ML_ST+1
	IF(NF .GT. NFREQ_MAX)THEN
	  WRITE(LUER,*)'Error in SET_PROF_V4 --- insufficient storage'
	  WRITE(LUER,*)'Max # of frequency storage locations is:',NFREQ_MAX
	  WRITE(LUER,*)'# of storage locatins requires is:',NF
	  STOP
	END IF
	ID=PROF_LIST_LOCATION		!Location in arrays
	XNU(1:NF)=NU(ML_ST:ML_END)
900	FORMAT(X,A,T20,A,T26,F20.10,2I6,3I10,F6.1)
	IF(PROF_TYPE .EQ. 'BS_HHE' .OR. PROF_TYPE .EQ. 'LEMKE_HI')THEN
	  IF(PROF_TYPE .EQ. 'BS_HHE')THEN
            IF(VERBOSE)WRITE(LUER,900)'Using BS_HHE: ',LST_SPECIES(ID),LST_WAVE(ID),
	1                LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
          ELSE
            IF(VERBOSE)WRITE(LUER,900)'Using LEMKE:  ',LST_SPECIES(ID),LST_WAVE(ID),
	1                LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	  END IF
	  ITR=LST_PID(ID)		!Location in file.
	  CALL RD_BS_STRK_TAB(LST_SPECIES(ID),NU_ZERO,
	1               LST_NL(ID),LST_NUP(ID),PROF_TYPE,ITR,LU_STK)
	  IF(ND .EQ. 1)THEN
	    CALL STRK_BS_HHE(PR_GRIEM,XNU,NF,
	1               ED_MOD,TEMP_IN,VTURB_IN,ND,
	1               NU_ZERO,AMASS_IN,L_FALSE)
	    PROF_STORE(1,1:NF,LOC_INDX)=PR_GRIEM(1:NF)
	  ELSE
	    CALL STRK_BS_HHE(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1               ED_MOD,TEMP_IN,VTURB_IN,ND,
	1               NU_ZERO,AMASS_IN,L_FALSE)
	  END IF
!
	ELSE IF(PROF_TYPE .EQ. 'HeI_IR_STRK')THEN
          IF(VERBOSE)WRITE(LUER,900)'Using HeI_IR_STRK: ',LST_SPECIES(ID),LST_WAVE(ID),
	1              LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	  IF(ND .EQ. 1)THEN
	    CALL STRK_HEI_IR(PR_GRIEM,XNU,NF,
	1                   ED_MOD,TEMP_IN,VTURB_IN,
	1                   POP_PROTON,POP_HEPLUS,ND,NU_ZERO,
	1                   LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                   AMASS_IN,PROF_TYPE,LU_STK)
	    PROF_STORE(1,1:NF,LOC_INDX)=PR_GRIEM(1:NF)
	  ELSE
	    CALL STRK_HEI_IR(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1                   ED_MOD,TEMP_IN,VTURB_IN,
	1                   POP_PROTON,POP_HEPLUS,ND,NU_ZERO,
	1                   LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                   AMASS_IN,PROF_TYPE,LU_STK)
	   END IF
!
	ELSE IF(PROF_TYPE .EQ. 'DS_STRK')THEN
          IF(VERBOSE)WRITE(LUER,900)'Using DS_STRK: ',LST_SPECIES(ID),LST_WAVE(ID),
	1              LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	  IF(ND .EQ. 1)THEN
	    CALL STRK_DS(PR_GRIEM,XNU,NF,
	1	         ED_MOD,TEMP_IN,VTURB_IN,ND,NU_ZERO,
	1                LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                AMASS_IN,PROF_TYPE,LU_STK)
	    PROF_STORE(1,1:NF,LOC_INDX)=PR_GRIEM(1:NF)
	  ELSE
	    CALL STRK_DS(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1	         ED_MOD,TEMP_IN,VTURB_IN,ND,NU_ZERO,
	1                LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                AMASS_IN,PROF_TYPE,LU_STK)
	  END IF
	ELSE IF(PROF_TYPE .EQ. 'HZ_STARK')THEN
!
! Check validity of passed parameters for HI and HeII lines.
!
	  IF(NUP .LE. NL)THEN
	    WRITE(LUER,*)'Error in SET_PROF: NUP < NL'
	    STOP
	  END IF
          IF(VERBOSE)WRITE(LUER,900)'Using HZ_STARK: ',' ',0.01D0*C_KMS/NU_ZERO,
	1              NL,NUP,ML_CUR,ML_ST,ML_END,Z_IN
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
	  PROF_STORE(1:ND,1:NF,LOC_INDX)=0.0_LDP
!$OMP PARALLEL DO PRIVATE(I,PR_GRIEM,TMP_ED)
	  DO I=1,ND
	    TMP_ED=MIN(ED_IN(I),MAX_PROF_ED)
            CALL GRIEM_V2(PR_GRIEM,DWS_GRIEM,NF,
	1        TMP_ED,TEMP_IN(I),VTURB_IN(I),
	1        NL,NUP,Z_IN,AMASS_IN,RET_LOG)
	    PROF_STORE(I,1:NF,LOC_INDX)=PR_GRIEM(1:NF)
	  END DO
!$OMP END PARALLEL DO
	ELSE
	  WRITE(LUER,*)'Error in SET_PROF_V3'
	  WRITE(LUER,*)'Unrecognized profile type'
	  WRITE(LUER,*)'PROF_TYPE=',PROF_TYPE
	  WRITE(LUER,*)NU_ZERO,NL,NUP,AMASS_IN
	  STOP
	END IF			!Profile type
!
! Normalize profile to an integral of unity assuming a trapazoidal integral.
! Since NU is units of 10^15 Hz, we need to multiply the integral by
! 10^15 before scaling.
!
!	NORM_PROFILE=.FALSE.
	IF(NORM_PROFILE)THEN
	  NF=ML_END-ML_ST+1
!$OMP PARALLEL DO IF(ND > 4) PRIVATE (I,T1,ML)
	  DO I=1,ND
	    PROF_STORE(I,1,LOC_INDX)=0.0_LDP
	    PROF_STORE(I,NF,LOC_INDX)=0.0_LDP
	    T1=0.0_LDP
	    DO ML=1,NF-1
	        T1=T1+0.5_LDP*(NU(ML_ST+ML-1)-NU(ML_ST+ML))*
	1                (PROF_STORE(I,ML,LOC_INDX)+PROF_STORE(I,ML+1,LOC_INDX))
	    END DO
	    T1=ABS(T1)*1.0E+15_LDP
	    PROF_STORE(I,1:NF,LOC_INDX)=PROF_STORE(I,1:NF,LOC_INDX)/T1
	  END DO
	END IF
!
! Remember to return profile.
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
	NF_STORE(LOC_INDX)=NF
	NU_STORE(1:NF,LOC_INDX)=NU(ML_ST:ML_END)
!
!
!	I=ND-7
!	WRITE(99,*)AMASS_IN,NL,NUP
!	WRITE(99,*)ED_MOD(I),TEMP_IN(I),VTURB_IN(I)
!	WRITE(99,'(1P,5E15.5)')0.01D0*C_KMS/NU_STORE(1:NF,LOC_INDX)
!	WRITE(99,'(1P,5E15.5)')PROF_STORE(I,1:NF,LOC_INDX)
!
	RETURN
	END
