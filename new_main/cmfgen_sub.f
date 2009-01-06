!                                                     
! Main Subroutine to self consistently solve the equation of transfer and the
! equations of statistical equilibrium for a spherically extended atmosphere 
! in the presence of outflows.
!                                 
! At present the atmosphere considered to consist of 
!      H, He, plus other species such as C, N, O and Fe.
! Eithe H or He must be present.
!
! Abundances are un-normalized relative fractional abundances (i.e. specified 
! specified with respect to some arbitrary species X. Generally we have used
! He as X. If negative they are interpreted as mass-fractions.
!
! Several options are available for handeling both the lines and continuum.
! Some options have not been recently tested.
!
! Inclusion of new species is generally straightforward. Only the
! calling routine, CMFGEN, needs to modified. The main work is in
! collating the atomic data. Dynamical memory allocation is used.
!
! 
!
	SUBROUTINE CMFGEN_SUB(ND,NC,NP,NT,
	1                     NUM_BNDS,NION,DIAG_INDX,
	1                     NDMAX,NPMAX,NCF_MAX,NLINE_MAX,
	1                     TX_OFFSET,MAX_SIM,NM,NM_KI,NLF)
	USE MOD_CMFGEN
	USE ANG_QW_MOD
	USE CMF_SOB_MOD
	USE CONTROL_VARIABLE_MOD
	USE OPAC_MOD
	USE STEQ_DATA_MOD
	USE MOD_LEV_DIS_BLK
	USE LINE_VEC_MOD
	USE LINE_MOD
	USE RADIATION_MOD
	USE UPDATE_KEYWORD_INTERFACE
	USE VAR_RAD_MOD
	IMPLICIT NONE
!
! Altered 16-Feb-2006 : Changed and modified over 2 month period. Section solving for
!                         populations and performing NG acceleration etc removed to
!                         subroutine (SOLVE_FOR_POPS). Routines added to allow time 
!                         variability of the statistical equilibrium equations. 
!                         Currently these are only for a Hubble law flow. Relativistic 
!                         terms added to COMP_OBS (now COMP_OBS_V2).
! Altered 20-Feb-2005 : Changed to use FLUX_MEAN & ROSS_MEAN which are defined in
!                         MOD_CMFGEN. Previusly used FLUXMEAN & ROSSMEAN defined in
!                         RADIATION_MOD.
! 
	INTEGER ND,NC,NP,NT
	INTEGER NUM_BNDS,NION,DIAG_INDX
	INTEGER NDMAX,NPMAX
	INTEGER NCF_MAX,NLINE_MAX
	INTEGER NM,NM_KI,MAX_SIM,NLF
	INTEGER TX_OFFSET
!
	INTEGER NCF
	LOGICAL, PARAMETER :: IMPURITY_CODE=.FALSE.
!
	CHARACTER(LEN=12), PARAMETER :: PRODATE='23-Nov-2007'		!Must be changed after alterations
!
! 
!
	REAL*8 SOL(NT,ND)		!Temp. stor. area for ST. EQ.
	INTEGER DST,DEND
!
! Constants for opacity etc. These are set in CMFGEN.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
!
! Internally used variables
!
	REAL*8 S1,REPA
	REAL*8 MAXCH,MAXCH_SUM
	REAL*8 T1,T2,T3,T4,SRAT
	REAL*8 FL,AMASS,FL_OLD
	REAL*8 FG_COUNT
!
	LOGICAL LST_DEPTH_ONLY
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
 	INTEGER REC_SIZE
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER N_PER_REC
!
! 
!
! Logical Unit assignments. Those indicated with a # after the ! are open in
!  large sections of the code. Other units generally used temprarily.
!
	INTEGER               LUER        !Output/Error file.
	INTEGER, PARAMETER :: LUIN=7      !General input unit (closed after accesses).
	INTEGER, PARAMETER :: LUMOD=8     !Model Description file.
!
	INTEGER, PARAMETER :: LU_DC=9      	!Departure coefficient Output.
	INTEGER, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER, PARAMETER :: LU_SE=16     	!Statistical equilibrium and Solution Arrays.
	INTEGER, PARAMETER :: LU_NET=17    	!# Line Netrate data.
	INTEGER, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER, PARAMETER :: LU_DR=19     	!# Downward rate (Nu. Z. A).
	INTEGER, PARAMETER :: LU_EW=20     	!# EW data.
	INTEGER, PARAMETER :: LU_REC_CHK=21	!# EW data.
!
! For writing scratch file (SCRTEMP). Also used in reading in MODEL data.
!
	INTEGER, PARAMETER :: LUSCR=26
!
	INTEGER, PARAMETER :: LU_HT=27     !#LINEHEAT (i.e Line heating term in R.E. equation)
!
! Used for RVTJ file and POPCARB, POPNIT etc.
!
	INTEGER, PARAMETER :: LU_POP=30
!
	INTEGER, PARAMETER :: LU_IMP=34       !J and CHI for impurity calculation.
	INTEGER, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER, PARAMETER :: LU_JEW=36       !J for EW computation (JEW)
	INTEGER, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER, PARAMETER :: LU_ES=38        !ES_J_CONV
!
! Following is used output the BA matrix, and its associated
! pointer file.
!
	INTEGER, PARAMETER :: LU_BA=40 
!
! For listing of transitions with TOTAL negative opacity values at some depths.
!
	INTEGER, PARAMETER :: LU_NEG=75
!
! 
!
	INTEGER NNM				!Include cont. var in line var.
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER MNL_F,MNUP_F
	INTEGER PHOT_ID
	INTEGER I,J,K,L,ML,LS,LINE_INDX,NEXT_LOC
	INTEGER IREC,MATELIM
!
	CHARACTER*80 TMP_STRING
	CHARACTER*20 TMP_KEY
	INTEGER ID,LOC_ID,ID_SAV,JJ
        INTEGER  L1,L2,U1,U2
	INTEGER IT,MNT,NIV
	INTEGER ISPEC
!
! Main iteration loop variables.
!
	INTEGER MAIN_COUNTER,NITSF,NUM_ITS_RD,NUM_ITS_TO_DO
	LOGICAL LST_ITERATION
!
! Functions called
!
	INTEGER ICHRLEN,ERROR_LU
	REAL*8 DOP_PRO
	REAL*8 S15ADF
	REAL*8 LAMVACAIR
	REAL*8 ATOMIC_MASS_UNIT
	REAL*8 SPEED_OF_LIGHT
	REAL*8 GRAVITATIONAL_CONSTANT
	REAL*8 RAD_SUN
	REAL*8 TEFF_SUN
	REAL*8 MASS_SUN
	LOGICAL EQUAL
	EXTERNAL ICHRLEN,ERROR_LU,SPEED_OF_LIGHT,GRAVITATIONAL_CONSTANT
	EXTERNAL MASS_SUN,RAD_SUN,TEFF_SUN
!
	INTEGER GET_DIAG
	INTEGER BNDST
	INTEGER BNDEND
	INTEGER BND_TO_FULL
!
! Photoionization cross-section routines.
!
! Collisional routines.
!
	EXTERNAL OMEGA_GEN_V3
!
! Wind variablity arrays.
!
	REAL*8 POPS(NT,ND)		!Population for all species.
	REAL*8 MEAN_ATOMIC_WEIGHT	!Mean atomic weight of atoms  (neutrals
!                          		! and ions) in atomic mass units.
	REAL*8 ABUND_SUM
!
!
! Arrays for improving on the initial T structure --- partition functions. 
! Need one for each atomic species.
!
	REAL*8, ALLOCATABLE :: U_PAR_FN(:,:)
	REAL*8, ALLOCATABLE :: PHI_PAR_FN(:,:)
	REAL*8, ALLOCATABLE :: Z_PAR_FN(:)
!
	REAL*8 TGREY(ND)
	REAL*8 T_SAVE(ND)
!
! Variables for scaling the line cooling rates in oder that the radiative
! equilibrium equation is more consistent with the electron heating/cooling 
! equation. The scaling is done when the line frequency is with a fraction
! of SCL_LINE_HT_FAC of the tmean frequency for the super-level under 
! consideration. 0.5 is presently the prefered value.
!
	REAL*8 AVE_ENERGY(NT)		!Average energy of each super level
	REAL*8 STEQ_T_SCL(ND)
	REAL*8 STEQ_T_NO_SCL(ND)
! 
!
! Dielectronic recombination variables and arrays.
!
	INTEGER NMAXDIE
	PARAMETER (NMAXDIE=500)
!
	REAL*8 EDGEDIE(NMAXDIE)		!Ionization frequency (negative)
	REAL*8 EINADIE(NMAXDIE)		!Einstein A coefficient
	REAL*8 GUPDIE(NMAXDIE)		!Stat. weight of autoionizing level.
!
	INTEGER LEVDIE(NMAXDIE)  	!Indicates MNL of low state
	INTEGER INDXDIE(NMAXDIE)
!
! Used for species identification as INDXDIE is not unique (specied not
! present can have same index as a species thats present.)
!
	CHARACTER*10 SPECDIE(NMAXDIE)
	CHARACTER*35 DIENAME(NMAXDIE)
!
	INTEGER NDIETOT
!
! Arrays and variables used for both Dielectronic recombination, and
! the implicit recombination.
!
	INTEGER EQION,EQSPEC
	REAL*8 GLOW,GION
	REAL*8 NUST(ND)			!LTE autoionizing population.
	REAL*8 DION(ND)			!Ion population
	REAL*8, ALLOCATABLE :: DIERECOM(:,:)   !Chk for all species.
	REAL*8, ALLOCATABLE :: DIECOOL(:,:)    !Dielec. cooling check for all spec.
	REAL*8, ALLOCATABLE :: ADDRECOM(:,:)     
! 
!
! Opacity/emissivity
!
	REAL*8 CHIL(ND)                 !Line opacity (without prof.)
	REAL*8 ETAL(ND)                 !Line emissivity (without prof.)
!
! Quadrature weights.
!
	REAL*8 FQW(NCF_MAX)		!Frequency weights
!
! Transfer equation vectors
	REAL*8 R_OLD(NDMAX)		!Used to store previous R grid in SN models.
!
! Line vectors
	REAL*8 AV(ND)
	REAL*8 VB(NDMAX)		!Used for error calculations
	REAL*8 VC(NDMAX)		!Used for error calculations
	REAL*8 H(ND)
	REAL*8 Q(ND)			!FREQ DEPENDENT.
	REAL*8 QH(ND)			!  "      "
	REAL*8 GAM(ND)			!FREQ INDEPENDENT
	REAL*8 GAMH(ND)			!  "      "
! 
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	LOGICAL EDDINGTON
!
! Variables for EW's and LINE blanketing.
!
	REAL*8 CONT_INT,EW
	INTEGER ACCESS_JEW
	LOGICAL COMPUTE_EW,COMPUTE_JEW,COMPUTE_LAM,MID,FULL_ES
!
! ACESS_F is the current record we are writing in EDDFACTOR.
! EDD_CONT_REC is the record in EDDFACTOR which points to the first
! record containing the continuum values.
!
	INTEGER ACCESS_F
	INTEGER, PARAMETER :: EDD_CONT_REC=3
!
	INTEGER NDEXT,NCEXT,NPEXT
!
	REAL*8 CNM(NDMAX,NDMAX)		!For collisions cross-section in
	REAL*8 DCNM(NDMAX,NDMAX)	!STEQGEN
!
!
	INTEGER, PARAMETER :: N_FLUX_MEAN_BANDS=12
	REAL*8     LAM_FLUX_MEAN_BAND_END(N_FLUX_MEAN_BANDS)
	REAL*8     BAND_FLUX_MEAN(ND,N_FLUX_MEAN_BANDS)
	REAL*8     BAND_FLUX(ND,N_FLUX_MEAN_BANDS)
	DATA LAM_FLUX_MEAN_BAND_END/100.0D0,150.0D0,200.0D0,227.83D0,258.90D0,300.0D0,504.25D0,911.75D0,
	1                         1200.0D0,1500.0D0,2000.0D0,1.0D+08/
!
!
! Continuum frequency variables and arrays.
!
	REAL*8 NU(NCF_MAX)		!Continuum and line frequencies
	REAL*8 NU_EVAL_CONT(NCF_MAX)	!Frequencies to evaluate continuum
	REAL*8 OBS(NCF_MAX)		!Observers spectrum
!
! Vectors and arrays used for the observed flux.
!
	INTEGER N_OBS
	REAL*8 OBS_FREQ(NCF_MAX)		!Since N_OBS < NCF =< NCF_MAX
	REAL*8 OBS_FLUX(NCF_MAX)
	LOGICAL FIRST_OBS_COMP
!
	CHARACTER TIME*20
	CHARACTER FMT*120
	CHARACTER*20 SECTION,FORMAT_DATE*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*132
	CHARACTER TEMP_CHAR*132
	CHARACTER*2 FORMFEED
!
! Global vectors:
!
	REAL*8 AMASS_ALL(NT)
	INTEGER N_LINE_FREQ
!
	INTEGER LINES_THIS_FREQ(NCF_MAX)
!
	REAL*8 NU_DOP
	REAL*8 NU_MAX_OBS
	REAL*8 NU_MIN_OBS
!
	INTEGER FREQ_INDX
	INTEGER X_INDX
	INTEGER FIRST_LINE
	INTEGER LAST_LINE
!
! Variables to limit the computation of the continuum opacities and
! emissivities.
!
	REAL*8 JREC(ND)
	REAL*8 dJRECdT(ND)
	REAL*8 JPHOT(ND)
	REAL*8 JREC_CR(ND)
	REAL*8 JPHOT_CR(ND)
	REAL*8 BPHOT_CR(ND)
!
	REAL*8 CONT_FREQ
	LOGICAL FINAL_CONSTANT_CROSS
!
! Indicates whether APRXzV, FFXzZ etc should be zeroed.
!
	LOGICAL ZERO_REC_COOL_ARRAYS 
!
! 
!
	REAL*8 Z_POP(NT)		!Ionic charge for each species
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL FIRST_FREQ
	LOGICAL RAT_TOO_BIG
	LOGICAL NEW_FREQ
!
! 
!
! X-ray variables.
! We dimension from 0 so that we can access a Null vector for the 1st included
! ioinization stage of each species.
! 
	REAL*8 X_RECOM(ND,0:NION)			!Next X-ray recombination rate
	REAL*8 X_COOL(ND,0:NION)			!Next X-ray cooling
!
	REAL*8 OBS_XRAY_LUM_0P1
	REAL*8 OBS_XRAY_LUM_1KEV
	REAL*8 GFF,XCROSS_V2
	EXTERNAL GFF,XCROSS_V2
!
	REAL*8 SPEC_DEN(ND,NUM_SPECIES)		!Used by ELEC_PREP
	REAL*8 AT_NO_VEC(ND,NUM_SPECIES)
!
	REAL*8 AD_COOL_V(ND)
	REAL*8 AD_COOL_DT(ND)
	REAL*8 ARTIFICIAL_HEAT_TERM(ND)
	REAL*8 dE_RAD_DECAY(ND)
	REAL*8 dE_WORK(ND)
!
	LOGICAL FIRST
	LOGICAL CHK,SUCCESS
        LOGICAL VAR_SOB_JC
	LOGICAL NEG_OPACITY(ND),FIRST_NEG
	LOGICAL AT_LEAST_ONE_NEG_OPAC
	LOGICAL FILE_OPEN
!
! Inidicates approximate frequencies for which TAU at outer boundary is written
! to OUTGEN on the last iteration.
!
! They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
!
	INTEGER, PARAMETER :: N_TAU_EDGE=5
	REAL*8 TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16D0,11.60D0,5.95D0,3.29D0,0.83D0/
!
!***********************************************************************
!
!*******************FUNCTION DEFINITIONS********************************
!
! This function takes a band-index and converts it the equivalent index
! in the full matrix. L=BND_TO_FULL(J,K) is equivalent to the statements:
!     IF(NUM_BNDS .EQ. ND)THEN L=J ELSE L=K+J-DIAG_INDX END IF
! The second indice is the equation depth.
!
	BND_TO_FULL(J,K)=(NUM_BNDS/ND)*(DIAG_INDX-K)+K+J-DIAG_INDX
!
! This function computes the index L on BA( , ,?,K) corresponding
! to the local depth variable (i.e that at K). It is equivalent
! to IF (NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
!
	GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_INDX)+DIAG_INDX
!
! These two functions compute the start and end indices when updating
! VJ. eg. we do not wish to update VJ( ,1, ) if we are using the banded
! matrix since this refers to a variable beyond the outer atmosphere.
!
	BNDST(K)=MAX( (NUM_BNDS/ND)*(K-DIAG_INDX)+1+DIAG_INDX-K, 1 )
	BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K, 
	1                 NUM_BNDS )
! 
!
!****************************************************************************
!
! Initialization section
!
	LUER=ERROR_LU()
	ACCESS_F=5
	COMPUTE_LAM=.FALSE.
	COMPUTE_EW=.TRUE.
	FULL_ES=.TRUE.
	SN_MODEL=.FALSE.
	ZERO_REC_COOL_ARRAYS=.TRUE.
	I=12
	FORMFEED=' '//CHAR(I)
	CNT_FIX_BA=0
	MAXCH_SUM=0.0D0
	LST_ITERATION=.FALSE.
!
!
! When TRUE, FIXED_T indicated that T is to be heled fixed (at least at some
! depths) in the linearization. This variable is set automatically by the 
! code depending on the magnitude of the corrections.
!
	FIXED_T=.FALSE.
!
! Set NDIETOT to zero in case no dielectronic lines are include in
! the dielectronic section.
!
	NDIETOT=0
!
! Set so that it is defined for the TEST whether to do an accurate flux
! calculation.
!
	MAXCH=100.D0
!
! A value of 1000 is used to indicate that the last change was greater them
! VAL_DO_NG, or a NEW_MODEL.
!
! A value of 2000 indicates a continuing model. In this case LAST_NG 
! take precedence. NB: NEXT_NG is reset to 1000 for a new model in the
! new model section.
!
	NEXT_NG=2000			!Initial value Indicate model just bega
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! 
!
	CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','UNKNOWN',' ',' ',IZERO,IOS)
!
! Open a scratch file to record model parameters. This file will eventually
! be renamed MODEL.
!
	CALL GEN_ASCI_OPEN(LUSCR,'MODEL_SCR','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL_SCR in CMFGEN, IOS=',IOS
	  STOP
	END IF
	WRITE(LUSCR,'()')
!
! Read in parameters which can change during a single model run. These
! parameters contol the number of iterations, and whether we wish to perform
! a LAMBDA iteration.
!
	  CALL GEN_ASCI_OPEN(LUIN,'IN_ITS','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN  
	     WRITE(LUER,*)'Error opening IN_ITS in CMFGEN, IOS=',IOS
	     STOP
	  END IF
	  CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
	  CALL RD_STORE_INT(NUM_ITS_TO_DO,'NUM_ITS',L_TRUE,'Number of iterations to perform')
	  CALL RD_STORE_LOG(RD_LAMBDA,'DO_LAM_IT',L_TRUE,'Do LAMBDA iterations ?')
	  DO_LAMBDA_AUTO=.TRUE.
	  CALL RD_STORE_LOG(DO_LAMBDA_AUTO,'DO_LAM_AUTO',L_FALSE,
	1                         'Start non-lambda iterations automatically?')
	  DO_GREY_T_AUTO=.TRUE.
	  CALL RD_STORE_LOG(DO_GREY_T_AUTO,'DO_GT_AUTO',L_FALSE,
	1                         'Do a grey temperature iteration after revising USE_FIXED_J?')
	  CALL CLEAN_RD_STORE()
	  CLOSE(UNIT=LUIN)
	NUM_ITS_RD=NUM_ITS_TO_DO
!
! These two parameters were originally read in from IN_ITS, however they
! are no longer required by program when using the BAND solution. They
! are still passed, however, to SOLVEBA.
!
	REPA=1.2
	MATELIM=1

!
! This section does the following :-
!	1) If new model read/determine the new radius scale  and
! populations.
!	2) If old model, populations are read in from scratch file
! previous radius scale etc are used. Ponit1 and point2
! point to the input data record (Note : Single Record)
!
	CALL RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
!
! RMDOT is the density at R=10dex10 cm and V=1km/s (atomic mass units)
!
	RMDOT=RMDOT*3.02286E+23
!
! LAMBDA_ITERATION controls whether a LAMBDA iteration is performed.
! A Lambda iteration is forced if RD_LAMDA is true. FIXED_T is set true
! as the Radiative equilibrium equation is not linearized if we are
! performing a LAMBDA iteration (could be changed with effort).
! The FIX_IMPURITY option is only used in non-LAMBDA mode.
!
	LAMBDA_ITERATION=RD_LAMBDA
	IF(LAMBDA_ITERATION)THEN
	  FIX_IMPURITY=.FALSE.
          FIXED_T=.TRUE.
	ELSE
	  FIX_IMPURITY=RD_FIX_IMP
	END IF
!
! 
!
	IF(ACCURATE)THEN
!
! We first verify that the interpolation range is valid.
!
	  IF(END_INTERP_INDX .GT. ND)END_INTERP_INDX=ND
	  IF(DEEP .GT. ND)DEEP=MIN(5,ND)
	  NDEXT=(END_INTERP_INDX-ST_INTERP_INDX)*NPINS+ND
	  IF(NDEXT .GT. NDMAX)THEN
	    WRITE(LUER,*)' Error - NDEXT larger than NDMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NDMAX in CMFGEN'
	    STOP
	  END IF
	  NCEXT=NC
!
! NB: The following expression guarentees that NPEXT has the same relationship
! to NDEXT and NCEXT as does NP to ND and NC.
!
	  NPEXT=NDEXT+NCEXT+(NP-ND-NC)
	  IF(NPEXT .GT. NPMAX)THEN
	    WRITE(LUER,*)' Error - NPEXT larger than NPMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NPMAX in CMFGEN'
	    STOP
	  END IF
	ELSE
	  NDEXT=ND; NCEXT=NC; NPEXT=NP
	END IF
!
	CALL SET_RADIATION_MOD(ND,NDMAX,NPMAX)
	CALL SET_LINE_MOD(ND,MAX_SIM,NM)
        CALL SET_VAR_RAD_MOD_V2(ND,NDEXT,
	1        NT,NUM_BNDS,NM,MAX_SIM,NM_KI,ACCURATE,L_TRUE)
	CALL SET_CMF_SOB_MOD(ND,NUM_BNDS,NT,NM_KI,NLF,LUER)
!
	T1=10.0/(NLF-1)
	DO ML=1,NLF
	  PF(ML)=5.0D0-T1*(ML-1)
	END DO
! 
!
! Read in bound-free gaunt factors for individual n states of hydrogen,
! and hydrogenic cross-sections for individual l states (n =0 to 30,
! l=0 to n-1)
!
	CALL RD_HYD_BF_DATA(LUIN,LUSCR,LUER)
!
! Read in atomic data for 2-photon transitions.
!
	CALL RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
!
! Read in data for charge exchange reactions. 
!
	CALL RD_CHG_EXCH_V3(LUIN,INCL_CHG_EXCH)
!
! Read in X-ray photoionization cross-sections.
!
	CALL RD_XRAY_FITS(LUIN)
!
	IF(XRAYS .AND. .NOT. FF_XRAYS)THEN
	  CALL RD_XRAY_SPEC(T_SHOCK_1,T_SHOCK_2,LUIN)
	END IF
!
	IF(XRAYS .AND. ADD_XRAYS_SLOWLY .AND. RD_LAMBDA)THEN
	   FILL_X1_SAV=FILL_FAC_XRAYS_1
	   FILL_X2_SAV=FILL_FAC_XRAYS_2
	   FILL_FAC_XRAYS_1=FILL_FAC_X1_BEG
	   FILL_FAC_XRAYS_2=FILL_FAC_X2_BEG
	END IF
!
	CALL RD_NUC_DECAY_DATA(INCL_RADIOACTIVE_DECAY,ND,LUIN)
!
! 
!
! Read in oscillator strengths, the photoionization cross section data,
! dielectronic data, and implicit recombination data for carbon. 
! Individual species are grouped together (rather than grouping all the
! oscillator reads) so that the headers of the INPUT files are grouped
! in the MODEL output file.
!
! We do this in reverse order so that GIONXzV can be correctly set.
!
! Note Well - in GENOSICL  - T1 is returned with the ionization energy.
!                            T2 is returned with screened nuclear charge.
!                             I is returned with the number of transitions.
!
! RDGENDIE returns the dielectronic transitions as a line list. These lines
! are treated as individual lines.
!
! RD_XzV_PHOT_DIE associates the dielectronic transitions with the
! photoionization cross-sections. They are then handeled as part of
! the continuum cross-sections.
!
! NB: The passed GF_CUT is set to zero if the Atomic NO. of the species
!      under consideration is less than AT_NO_GF_CUT.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF( ATM(ID)%XzV_PRES)THEN
	      IF( NINT(AT_NO(SPECIES_LNK(ID))) .LT. NINT(AT_NO_GF_CUT) )THEN
	        T2=0.0D0
	      ELSE
	        T2=GF_CUT
	      END IF
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_OSCDAT'
	      CALL GENOSC_V6( ATM(ID)%AXzV_F, ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1                 ATM(ID)%XzVLEVNAME_F, T1, ATM(ID)%ZXzV,
	1                 ATM(ID)%XzV_OSCDATE, ATM(ID)%NXzV_F,I,
	1                 'SET_ZERO',T2,GF_LEV_CUT,MIN_NUM_TRANS,
	1                 LUIN,LUSCR,TMP_STRING)
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_TO_S'
	      CALL RD_F_TO_S_IDS( ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1           ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F, ATM(ID)%NXzV,
	1           LUIN,TMP_STRING)
	      CALL RDPHOT_GEN_V2( ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1           ATM(ID)%GIONXzV_F,AT_NO(SPECIES_LNK(ID)),
	1           ATM(ID)%ZXzV, ATM(ID)%NXzV_F,
	1           ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT,  NPHOT_MAX,
	1           ATM(ID+1)%XzV_PRES,     ATM(ID+1)%EDGEXzV_F, ATM(ID+1)%GXzV_F,
	1           ATM(ID+1)%F_TO_S_XzV,   ATM(ID+1)%XzVLEVNAME_F, ATM(ID+1)%NXzV_F,
	1           SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE_EDGE,
	1           XRAYS,ID,ION_ID(ID),LUIN,LUSCR)   
              IF(ATM(ID+1)%XzV_PRES) ATM(ID)%GIONXzV_F= ATM(ID+1)%GXzV_F(1)
 	      IF(DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV) )THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RDGENDIE_V4( ATM(ID)%XzVLEVNAME_F, ATM(ID)%INDX_XzV,
	1             ATM(ID)%NXzV_F,
	1             EDGEDIE,EINADIE,GUPDIE,
	1             LEVDIE,INDXDIE,SPECDIE,DIENAME, ATM(ID)%GIONXzV_F,
	1             ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,L_TRUE,TMP_STRING,NMAXDIE,NDIETOT)
	      ELSE IF( ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV)THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RD_PHOT_DIE_V1(ID,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,    ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS, ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,TMP_STRING)
	      END IF
	    END IF
	  END DO
	END DO
!
! 
!
! We open a new MODEL file so that the information is at the head of the
! file.
!
	CALL GEN_ASCI_OPEN(LUMOD,'MODEL','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL in CMFGEN, IOS=',IOS
	  STOP
	END IF
!
! Output description of model. This is done after the reading of most data
! since some of the model information is read in.
!
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(//,'' Model Started on:'',15X,(A))')TIME
	  WRITE(LUMOD,
	1       '('' Main program last changed on:'',3X,(A))')PRODATE
	  WRITE(LUMOD,'()')
	  FMT='(5X,I8,5X,''!Number of depth points'')'
	  WRITE(LUMOD,FMT)ND
	  FMT='(5X,I8,5X,''!Number of core rays'')'
	  WRITE(LUMOD,FMT)NC
	  FMT='(5X,I8,5X,''!Total number of rays'')'
	  WRITE(LUMOD,FMT)NP
	  FMT='(5X,I8,5X,''!Total number of variables'')'
	  WRITE(LUMOD,FMT)NT
	  FMT='(5X,I8,5X,''!Maximum number of frequencies'')'
	  WRITE(LUMOD,FMT)NCF_MAX
	  FMT='(5X,I8,5X,''!Number of bands'')'
	  WRITE(LUMOD,FMT)NUM_BNDS
!	  
	  WRITE(LUMOD,'()')
	  CALL RITE_ATMHD_V2(LUMOD)
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      ISPEC=SPECIES_LNK(ID)
	      CALL RITE_ATMDES_V2( ATM(ID)%XzV_PRES, ATM(ID)%NXzV, 
	1          ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1          ATM(ID)%NXzV_F, ATM(ID)%GIONXzV_F, ATM(ID)%N_XzV_PHOT,
	1          AT_NO(ISPEC),AT_MASS(ISPEC),LUMOD,ION_ID(ID))
	    END IF
	  END DO
!
! Append VADAT information and atomic data headers to model file.
! Thus only have 1 model file output.
!
	  REWIND(LUSCR)
	  IOS=0
	  TEMP_CHAR=FORMFEED
	  DO WHILE(IOS .EQ. 0)
	    I=ICHRLEN(TEMP_CHAR)
	    IF(I .GT. 0)THEN
	      WRITE(LUMOD,'(A)')TEMP_CHAR(1:I)
	    ELSE
	      WRITE(LUMOD,'()')
	    END IF
	    READ(LUSCR,'(A)',IOSTAT=IOS)TEMP_CHAR
	 END DO
!
! Finished all data read, so can close LUMOD (output descriptor).
!
	CLOSE(UNIT=LUSCR,STATUS='DELETE')
	CLOSE(UNIT=LUMOD)
! 
!
! Set the vector Z_POP to contain the ionic charge for each species.
!
	DO I=1,NT
	  Z_POP(I)=0.0D0
	END DO
!
	DO ID=1,NUM_IONS-1
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1              ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
!
! Store atomic masses in vector of LENGTH NT for later use by line 
! calculations. G_ALL and LEVEL_ID  are no longer used due to the use
! of super levels.
!
	AMASS_ALL(1:NT)=0.0D0
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)AMASS_ALL( ATM(ID)%EQXzV: ATM(ID)%EQXzV+ATM(ID)%NXzV-1)=
	1         AT_MASS(SPECIES_LNK(ID))
	END DO
!
! 
!
! Compute profile frequencies such that for the adopted doppler
! velocity the profile ranges from 5 to -5 doppler widths.
! This section needs to be rewritten if we want the profile to
! vary with depth.
!
! ERF is used in computing the Sobolev incident intensity at the
! outer boundary. ERF = int from "x" to "inf" of -e(-x^2)/sqrt(pi).
! Note that ERF is not the error function. ERF is related to the
! complementary error function by ERF =-0.5D0 . erfc(X).
! S15ADF is a NAG routine which returns erfc(x).
!
! The incident Sobolev intensity is S[ 1.0-exp(tau(sob)*ERF) ]
! NB -from the definition, -1<erf<0 .
!
	T1=4.286299D-05*SQRT( TDOP/AMASS_DOP + (VTURB/12.85)**2 )
	J=0
	DO I=1,NLF
	  ERF(I)=-0.5D0*S15ADF(PF(I),J)
	  PF(I)=PF(I)*T1
	END DO
        VDOP_VEC(1:ND)=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
!
! Compute the frequency grid for CMFGEN. Routine also allocates the vectors 
! needed for the line data, sets the line data, and puts the line data into 
! numerical order.
!
	CALL SET_FREQUENCY_GRID(NU,FQW,LINES_THIS_FREQ,NU_EVAL_CONT,
	1               NCF,NCF_MAX,N_LINE_FREQ,
	1               OBS_FREQ,OBS,N_OBS,LUIN,IMPURITY_CODE)
!
! Define the average energy of each super level. At present this is
! depth independent, which should be adequate for most models.
! This average energy is used to scale the line cooling rates in
! the radiative equilibrium equation so that is more consistent
! with the electron cooling rate. The need for this scaling
! arises when levels within a super level have a 'relatively large'
! energy separation, and the dominat rates are scattering.
!
	AVE_ENERGY(:)=0.0D0
	DO ID=1,NUM_IONS-1
	   CALL AVE_LEVEL_ENERGY(AVE_ENERGY, ATM(ID)%EDGEXzV_F,
	1         ATM(ID)%GXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%EQXzV, 
	1         ATM(ID)%NXzV,   ATM(ID)%NXzV_F, NT, ATM(ID)%XzV_PRES)
	END DO
!
! 
!
! Check to see if old model. If so, read in R,V, SIGMA and POPS arrays.
! If not, set NEWMOD to .TRUE. We also check the format of the file,
! in case we are revising the R grid.
!
	NLBEGIN=0		! Initialize for lines.
	IREC=0                  ! Get last iteration
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LAST_NG,
	1                 WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
        IF( (REVISE_R_GRID .OR. DO_HYDRO) .AND. NEWMOD)THEN
	  WRITE_RVSIG=.TRUE.
	ELSE IF(NEWMOD)THEN
	  WRITE_RVSIG=.FALSE.
	ELSE IF(REVISE_R_GRID .OR. DO_HYDRO)THEN
	   IF(.NOT. WRITE_RVSIG)THEN
	     WRITE(LUER,*)'Error in CMFGEN_SUB with SCRTEMP'
	     WRITE(LUER,*)'Inconsistent format request: RVSIG must be writeen for each iteration'
	     WRITE(LUER,*)'Restart a fresh model or use REWRITE_SCR to correct file format'
	     STOP
	   END IF
	END IF
	IF(NEWMOD)THEN
	  WRITE(LUER,*)'Starting a new model.'
	  WRITE(LUER,*)'*_IN files will be used to start model'
	END IF
!
! Now does accurate flux calculation for a single iteration provided not a new model.
!
	IF(.NOT. NEWMOD)MAXCH=0.0D0
!
!		' OLD MODEL '
!
	IF(.NOT. NEWMOD)THEN
!
! Convert back from POPS array to individual matrices.
!
	  DO ID=1,NUM_IONS-1
	    CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV,ED,T,
	1         ATM(ID)%EQXzV, ATM(ID)%NXzV,
	1         NT,ND, ATM(ID)%XzV_PRES)
	  END DO
!
! We have now stored the revised populations back in their individual
! storage locations. For some species we have 2 atomic models. For these
! species we need to take the super-level populations and compute:
!
! 1. The LTE population off all level ls in the FULL atom.
! 2. The population off all levels in the FULL atom.
! 3. The LTE population off all super-levels.
!
! This is done by the following include statement, which  is a sequence of
! calls to the routine SUP_TO_FULL.
!
! Compute the ion population at each depth.                          
! These are required when evaluation the occupation probabilities, and hence
! the LTE populations.
!
	  DO J=1,ND
	    POPION(J)=0.0D0
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	  END DO
!
	  CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	  INCLUDE 'SUP_TO_FULL_V4.INC'
!
	ELSE
!
! Compute R, V, and SIGMA separately from DC read so that can evaluate
! POPHE etc, and angle quadrature weights. These are required in NEWMODEL
! section if iterating on the initial temperature structure.
!
	  IF(VELTYPE .EQ. 1)THEN
	    CALL STARNEW(R,V,SIGMA,RMAX,RP,RN,VRP,VINF,EPPS1,GAMMA1
	1    ,RP2,RN2,VRP2,VINF2,EPPS2,GAMMA2,ND,TA,TB,TC)
	  ELSE IF(VELTYPE .EQ. 2)THEN
	    CALL STARFIN(R,V,SIGMA,RMAX,RP,RN,VRP,VINF,EPPS1,GAMMA1
	1   ,RP2,RN2,VRP2,VINF2,EPPS2,GAMMA2,ND,TA,TB,TC)
	  ELSE IF(VELTYPE .EQ. 3 .OR. VELTYPE .EQ. 6)THEN
	    CALL STARPCYG_V3(R,V,SIGMA,RMAX,RP,
	1             SCL_HT,VCORE,VPHOT,VINF1,V_BETA1,V_EPPS1,
	1             VINF2,V_BETA2,V_EPPS2,
	1             NBND_INS,CONS_FOR_R_GRID,EXP_FOR_R_GRID, 
	1             ND,TA,TB,TC,RDINR,LUIN)
          ELSE IF(VELTYPE .EQ. 4)THEN
            CALL STARRAVE(R,V,SIGMA,ND,LUIN,RMAX,RP)
         ELSE IF(VELTYPE .EQ. 7)THEN
	    CALL RD_RV_FILE_V2(R,V,SIGMA,RMAX,RP,VINF,LUIN,ND,VEL_OPTION,NUM_V_OPTS)
         ELSE IF(VELTYPE .EQ. 10)THEN
	    CALL RV_SN_MODEL_V2(R,V,SIGMA,RMAX,RP,VCORE,V_BETA1,RDINR,LUIN,ND)
	 ELSE IF(VELTYPE .EQ. 11)THEN
	   CALL SET_RV_HYDRO_MODEL(R,V,SIGMA,SN_AGE_DAYS,RMAX,RP,RDINR,ND,LUIN)
           VINF=V(1)
	 ELSE
	   WRITE(LUER,*)'Invalid Velocity Law'
	    STOP
	  END IF
	END IF
!
! 
!
! Compute CLUMP_FAC(1:ND) which allow for the possibility that the wind is
! clumped. At the sime time, we compute the vectors which give the density,
! the atom density, and the species density at each depth.
!
	CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
!
! 
!
	IF(ACCURATE)THEN
	  I=ND-DEEP
	  IF(INTERP_TYPE .NE. 'LOG')THEN
	    WRITE(LUER,*)'Error in CMFGEN_SUB'
	    WRITE(LUER,*)'The INTERP_TYPE currently implemented is LOG'
	    STOP
	  END IF
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  TA(1:ND)=1.0D0	!TEXT not required, T currently zero
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,TA,SIGMA,ND)
!
          VDOP_VEC_EXT(1:NDEXT)=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
!
	END IF
!
! Need to calculate impact parameters, and angular quadrature weights here
! as these may be required when setting up the initial temperature
! distribution of the atmosphere (i.e. required by JGREY).
!
	CALL SET_ANG_QW(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,TRAPFORJ,ACCURATE)
!
! Allocate memory for opacities.
!
        CALL INIT_OPAC_MOD(ND,NT,L_TRUE)
! 
!
!		'NEW MODEL'
!
! Read in old estimates for the departure coefficents for the new model.
! T and ED are also estimated.
!
	IF(NEWMOD)THEN
	  NITSF=0
	  IREC=0
	  LAST_NG=-1000  			!Must be -1000
	  NEXT_NG=1000				!Must be initialized to 1000
	  CALL SET_NEW_MODEL_ESTIMATES(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	END IF
!
! VEXT and SIGMAEXT have already been computed. We need TEXT for
! convolving J with the electron scattering redistribution function.
!
	IF(ACCURATE)THEN
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,T,SIGMA,ND)
	ELSE
	  TEXT(1:ND)=T(1:ND)
	END IF
!
! 
!
! Set TSTAR which is required by TOTOPA_JILA. Its precise value is
! irrelevant (provided >0) as we always adopt the diffusion
! approximation. 
!
	TSTAR=T(ND)
!
! 
! Section allows the user to read in a modified solution vector.
! Useful for assisting convergence when the model is experiencing
! difficulty converging.
!
	CALL TUNE(IONE,'GIT')
	IF(NUM_ITS_TO_DO .EQ. 0)THEN
	  IF(RDINSOL)THEN
!
	    CALL LOCSOLUT(POPS,SOL,TA,30,NT,ND)
	    DO ID=1,NUM_IONS-1
	      CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1            ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND, 
	1            ATM(ID)%XzV_PRES)
	    END DO
!
! This include block also compute the LTE populations of the FULL model atom,
! and the SUPER level model atom.
!
	    CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	    INCLUDE 'SUP_TO_FULL_V4.INC'
	  END IF
!
! Write pointer file and output data necessary to begin a new
! iteration.
!
	  IF(RDINSOL .OR. NEWMOD)THEN
	    MAIN_COUNTER=NITSF+1
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                  LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  END IF
	  LST_ITERATION=.TRUE.
	  GOTO 9999			!End (write out POPS.)
	END IF
! 
!
! Associate charge exchange reactions with levels in the model atoms.
!
	CALL SET_CHG_LEV_ID_V4(ND,LUMOD)
	CALL VERIFY_CHG_EXCH_V3()
!
! Determine the number of important variables for each species, and
! set the links.
!
	CALL DETERMINE_NSE(NION,XRAYS)
        CALL CREATE_IV_LINKS_V2(NT,NION)
!
! Allocate memory for STEQ and BA arrays.
!
        CALL SET_BA_STORAGE(NT,NUM_BNDS,ND,NION)
!
! Read in BA and STEQ arrays. We only attempt this if we have an existing
! model.
!
	CHK=.FALSE.
	IF(.NOT. NEWMOD)THEN
          CALL READ_BA_DATA_V2(LU_BA,NION,NUM_BNDS,ND,COMPUTE_BA,CHK,'BAMAT')
	END IF
	IF(.NOT. CHK .OR. LAMBDA_ITERATION)THEN
	  NLBEGIN=0
          COMPUTE_BA=.TRUE.
	  WRBAMAT=.FALSE.
	  IF(N_ITS_TO_FIX_BA .GT. 0)WRBAMAT=.TRUE.
	ELSE IF(NLBEGIN .EQ. -999)THEN		!Indicate completed iteration
	  NLBEGIN=0				!hence BA matrix available.
	  COMPUTE_BA=COMPUTE_BARDIN
	  WRBAMAT=.FALSE.
	ELSE
	  WRBAMAT=WRBAMAT_RDIN
 	END IF
	IF(FLUX_CAL_ONLY)THEN
	   COMPUTE_BA=.FALSE.
	   WRBAMAT=.FALSE.
	   LAMBDA_ITERATION=.FALSE.
	   MAXCH=0.0D0
!
! For coherent electon scattering, only need 1 iteration.
!
           IF(RD_COHERENT_ES)NUM_ITS_TO_DO=1
	END IF
!
! Temporary check
!
!	CALL WRITE_SEQ_TIME_FILE_V1(SN_AGE_DAYS,ND,LUIN)
!	CALL TST_RD_EQ_FILE(POPS,ND,NT,LUIN)
!
! 
!
!**************************************************************************
!**************************************************************************
!
!                    MAIN ITERATION LOOP
!
!**************************************************************************
!**************************************************************************
!
! MAIN_COUNTER is an integer variable which keeps track of the TOTAL number
! of iterations performed. NB - A NG acceleration is counted as a single
! acceleration.
!
! NUM_ITS_TO_DO indicates the number of iterations left to do. For the
! last iteration this will be zero in the "2000" LOOP.
!
! LST_ITERATION is a logical variable which indicate that the current
! iteration is the last one, and hence DEBUGING and INTERPRETATION data
! should be written out. Its equivalent to NUM_ITS_TO_DO=0 in the loop.
!
	MAIN_COUNTER=NITSF			!Initialize main loop counter
	LAST_LAMBDA=NITSF
	LAST_AV=NITSF
	NEXT_AV=0
20000	CONTINUE
	NUM_ITS_TO_DO=NUM_ITS_TO_DO-1
	IF(NUM_ITS_TO_DO .EQ. 0)LST_ITERATION=.TRUE.
	MAIN_COUNTER=MAIN_COUNTER+1
!
	CALL MESS(LUER,'Current great iteration count is')
	WRITE(LUER,*)MAIN_COUNTER
!
! Used as a initializing switch for COMP_OBS.
!
	FIRST_OBS_COMP=.TRUE.
!
	IF(LST_ITERATION)THEN
	  CALL GEN_ASCI_OPEN(LU_NET,'NETRATE','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_DR,'TOTRATE','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_EW,'EWDATA','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_HT,'LINEHEAT','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LU_NEG,'NEG_OPAC','UNKNOWN',' ',' ',IZERO,IOS)
	END IF
!
	IF(IMPURITY_CODE)THEN
	  I=WORD_SIZE*(4*ND+1)/UNIT_SIZE
	  OPEN(UNIT=LU_EDD,FILE='IMPURITYJ',FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening IMPURITYJ in CMFGEN'
	    STOP
	  END IF
	ELSE
	  IF(ACCURATE .OR. EDD_CONT .OR. EDD_LINECONT)THEN
!
! NB: If not ACCURATE, NDEXT was set to ND. The +1 arises since we write
! NU on the same line as RJ. J is used to get the REC_LENGTH, while string
! will contain the date.
!
	    I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	    CALL READ_DIRECT_INFO_V3(K,J,STRING,'EDDFACTOR',LU_EDD,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error --- unable to open EDDFACTOR_INFO'
	      WRITE(LUER,*)'Will compute new f'
	      COMPUTE_EDDFAC=.TRUE.
	      IOS=0
	    ELSE IF(.NOT. COMPUTE_EDDFAC .AND. K .NE. ND)THEN
	      WRITE(LUER,*)'Error with EDDFACTOR_INFO'
	      WRITE(LUER,*)'Incompatible number of depth points'
	      WRITE(LUER,*)'ND is',ND 
	      WRITE(LUER,*)'ND is EDDFACTOR_INFO is',K
	      STOP
	    END IF
	    IF(.NOT. COMPUTE_EDDFAC)THEN
	      OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	      IF(IOS .EQ. 0)THEN
	        READ(LU_EDD,REC=5,IOSTAT=IOS)T1
	        IF(T1 .EQ. 0 .OR. IOS .NE. 0)THEN
	          WRITE(LUER,*)'Error --- All Eddfactors not'//
	1                      ' computed - will compute new F'
	          COMPUTE_EDDFAC=.TRUE.
	        END IF
	      ELSE
	        IF(.NOT. NEWMOD)THEN
	          WRITE(LUER,*)'Error opening EDDFACTOR - will compute new F'
	        END IF
	        COMPUTE_EDDFAC=.TRUE.
	      END IF
	    END IF
	    IF(COMPUTE_EDDFAC)THEN
	      IF(USE_FIXED_J)THEN
	        WRITE(LUER,*)'Error in CMFGEN_SUB'
	        WRITE(LUER,*)'Program will compute EDDFACTOR but this is'//
	1                      ' incompatable with US_FIXED_J=T'
	        STOP
	      END IF
	      CALL WRITE_DIRECT_INFO_V3(NDEXT,I,'20-Aug-2000','EDDFACTOR',LU_EDD)
	      OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='REPLACE',RECL=I)
	      WRITE(LU_EDD,REC=1)0
	      WRITE(LU_EDD,REC=2)0
	      WRITE(LU_EDD,REC=3)0
	      WRITE(LU_EDD,REC=4)0
!
! We set record 5 to zero, to signify that the eddington factors are
! currently being computed. A non zero value signifies that all values
! have successfully been computed. (Consistent with old Eddfactor
! format a EDD_FAC can never be zero : Reason write a real number).
!
	      T1=0.0
	      WRITE(LU_EDD,REC=5)T1
	    END IF
	  END IF
	END IF
!
! Now open file containing the electron scatterin J (i.e. the convolution of
! J with the e.s. redistribution function.)
!
! If we don't have EDDFACTOR file it is assumed that we don't have
! J_CONV also.
!
	IF(COMPUTE_EDDFAC)THEN
	  COHERENT_ES=.TRUE.
	ELSE IF(.NOT. COHERENT_ES)THEN
	   OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	     IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error opening ES_J_CONV'//
	1                    ' - will compute new J'
	        COHERENT_ES=.TRUE.
	     END IF
	END IF
!
	COMPUTE_JEW=.FALSE.
	I=WORD_SIZE*(ND+1)/UNIT_SIZE
	IF(COMPUTE_EW)THEN
	  OPEN(UNIT=LU_JEW,FILE='JEW',FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    IF(.NOT. NEWMOD)THEN
	      WRITE(LUER,*)'Error opening JEW - will compute new JEW'
	    END IF
	    COMPUTE_JEW=.TRUE.
	  END IF
	END IF
	IF(COMPUTE_JEW)THEN
	    OPEN(UNIT=LU_JEW,FILE='JEW',FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='REPLACE',RECL=I)
	END IF
! 
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	 DO J=1,ND
	    POPION(J)=0.0D0
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	 END DO
!
! 
!
! This routine not only evaluates the LTE populations of both model atoms, but
! it also evaluates the dln(LTE Super level Pop)/dT.
!
	CALL EVAL_LTE_V4(DO_LEV_DISSOLUTION,ND)
!
! 
!
! Set 2-photon data with current atomic models and populations.
!
	DO ID=1,NUM_IONS-1
	  ID_SAV=ID
	  CALL SET_TWO_PHOT_V2(ION_ID(ID), ID_SAV,
	1       ATM(ID)%XzVLTE,     ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F,   ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,  ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ND,
	1       ATM(ID)%ZXzV,       ATM(ID)%EQXzV,  ATM(ID)%XzV_PRES)
	END DO
!
! 
!
! Section to compute DT/DR and D(DT/DR)/D? for use in the
! diffusion approximation. DIFFW is used to store the second
! derivatives.
!
! This was based on the subroutine DTSUB. No longer a subroutine as
! too many variables need to be included.
!
! Zero variation of DTDR vector
!
	DO I=1,NT
	  DIFFW(I)=0.0D0
	END DO
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the portion of the code that computes dTdR.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
!
	CALL TUNE(1,'DTDR')
	DTDR=0.0D0
	SECTION='DTDR'
	IF(IMPURITY_CODE .OR. FLUX_CAL_ONLY .OR. (RD_LAMBDA .AND. NEWMOD .AND. .NOT. SN_MODEL))THEN
	  DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
	  DIFFW(1:NT)=0.0D0
	ELSE
!
! We only need to compute the opacity at the innermost depth, but to save
! programing we will compute it at all depths. As this is only done once
! per iteration, not much time will be wasted.
!
! Setting LST_DEPTH_ONLY to true limits the computation of CHI, ETA, and
! dCHI and dETA to the inner boundary only (in some cases).
!
	  LST_DEPTH_ONLY=.TRUE.
!
! RJ is used in VARCONT to compute the varaition of ETA. In this section
! we only want the variation of CHI, so we initialize its value to zero.
! This prevents a floating point exception.
!
	  RJ(1:ND)=0.0D0
!
	  CONT_FREQ=0.0D0
	  DO ML=1,NCF
	    FREQ_INDX=ML
! 
	    FL=NU(ML)
	    IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	      COMPUTE_NEW_CROSS=.TRUE.
	      CONT_FREQ=NU_EVAL_CONT(ML)
	    ELSE
	      COMPUTE_NEW_CROSS=.FALSE.
	    END IF
!
	    CALL TUNE(1,'DTDR_OPAC')
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!	    INCLUDE 'OPACITIES_V4.INC'
	    CALL TUNE(2,'DTDR_OPAC')
!
! 
!
! Compute variation of opacity/emissivity. Store in VCHI and VETA.
!
	    CALL TUNE(1,'DTDR_VOPAC')
!	    INCLUDE 'VAROPAC_V4.INC'
	    CALL COMP_VAR_OPAC(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	    CALL TUNE(2,'DTDR_VOPAC')
! 
!
! Compute contribution to CHI and VCHI by lines.
!
! Section to include lines automatically with the continuum.
! Only computes line opacity at final depth point. This is used in the
! computation of dTdR.
!
! NB: Care must taken to ensure that this section remains consistent
!      with that in continuum calculation section.
!
	    CALL SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1          LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
!
! Add in line opacity.
!
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	       CHI(ND)=CHI(ND)+CHIL_MAT(ND,SIM_INDX)*LINE_PROF_SIM(SIM_INDX)
	    END IF
	  END DO
!
! Now do the line variation. This presently ignores the effect of a 
! temperature variation.
!
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      NL=SIM_NL(SIM_INDX)
	      NUP=SIM_NUP(SIM_INDX)
	      VCHI(NL,ND)=VCHI(NL,ND)+LINE_PROF_SIM(SIM_INDX)*
	1        LINE_OPAC_CON(SIM_INDX)*L_STAR_RATIO(ND,SIM_INDX)
	      VCHI(NUP,ND)=VCHI(NUP,ND)-LINE_PROF_SIM(SIM_INDX)*
	1        LINE_OPAC_CON(SIM_INDX)*U_STAR_RATIO(ND,SIM_INDX)*
	1        GLDGU(SIM_INDX)
	    END IF
	  END DO
!
! 
!
! Set TA = to the variation vector at the inner boundary.
!
	    CALL TUNE(1,'DTDR_VEC')
	    DO I=1,NT
	      TA(I)=VCHI(I,ND)
	    END DO
	    T1=HDKT*NU(ML)/T(ND)
!
! Increment Parameters
!
	    T3=FQW(ML)*TWOHCSQ*( NU(ML)**3 )*T1*EMHNUKT(ND)/
	1         CHI(ND)/T(ND)/(1.0D0-EMHNUKT(ND))**2
	    DTDR=DTDR+T3
	    DO I=1,NT-1
	      DIFFW(I)=DIFFW(I)+T3*TA(I)/CHI(ND)
	    END DO
	    DIFFW(NT)=DIFFW(NT)+T3*(TA(NT)/CHI(ND)-(T1*(1.0D0+EMHNUKT(ND))
	1           /(1.0D0-EMHNUKT(ND))-2.0D0)/T(ND))
	    CALL TUNE(2,'DTDR_VEC')
!
	  END DO
!
! The luminosity of the Sun is 3.826E+33 ergs/sec. For convenience
! DTDR will have the units  (E+04K)/(E+10cm) .
!
	  T1=LUM*7.2685D+11/R(ND)/R(ND)
	  DTDR=T1/DTDR
	  T1=( DTDR**2 )/T1
	  DO I=1,NT
	    DIFFW(I)=DIFFW(I)*T1
	  END DO
	END IF
	CALL TUNE(2,'DTDR')
	CALL TUNE(3,'  ')
!
	LST_DEPTH_ONLY=.FALSE.
	WRITE(LUER,*)'The value of DTDR is :',DTDR
! 
!
! Zero STEQ and BA arrays.
!
	DO ID=1,NION
	  SE(ID)%STEQ   =0.0D0
	  SE(ID)%BA     =0.0D0
	  SE(ID)%BA_PAR =0.0D0
	END DO
	STEQ_ED=0.0D0
	STEQ_T=0.0D0
	STEQ_T_SCL=0.0D0
	STEQ_T_NO_SCL=0.0D0
        BA_ED   = 0.0D0
        BA_T    = 0.0D0
        BA_T_PAR=0.0D0
!
	IF(.NOT .ALLOCATED(DIERECOM))THEN
	  ALLOCATE (DIERECOM(ND,NION))
	  ALLOCATE (ADDRECOM(ND,NION))
	  ALLOCATE (DIECOOL(ND,NION))
	END IF
	DIERECOM(:,:)=0.0D0
	ADDRECOM(:,:)=0.0D0
	DIECOOL(:,:)=0.0D0
	X_RECOM(:,:)=0.0D0
	X_COOL(:,:)=0.0D0
	DIELUM(:)=0.0D0
! 
!
! Compute the value of the S.E. equations and compute the variation 
! matrix for terms that are independent of Jv .
!
! DST and DEND can be adjusted to that we can avoid reading in the entire 
! diagonal of the BA array for each call to STEQ_MULTI.
!
! Assume all BA mtarix is in memory/
	DO K=1,1
	  DST=1       
	  DEND=ND
!
! In this case, each depth is treated fully (i.e. for all species and 
! ionization stages) before we go onto the next depth.
!
!	DO K=1,ND
!	  DST=K
!	  DEND=K
!
	  CALL TUNE(1,'STEQ')
          DO ID=1,NUM_IONS-1
            LOC_ID=ID
	    IF(ATM(ID)%XzV_PRES)THEN
	      TMP_STRING=TRIM(ION_ID(ID))//'_COL_DATA'
              CALL STEQ_MULTI_V7(CNM,DCNM,ED,T,
	1         ATM(ID)%XzV,       ATM(ID)%XzVLTE,   ATM(ID)%dlnXzVLTE_dlnT,
	1         ATM(ID)%NXzV,      ATM(ID)%DXzV,     ATM(ID)%XzV_F, 
	1         ATM(ID)%XzVLTE_F,  ATM(ID)%W_XzV_F,  ATM(ID)%AXzV_F,
	1         ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,   ATM(ID)%XzVLEVNAME_F,
	1         ATM(ID)%NXzV_F,    ATM(ID)%F_TO_S_XzV,
	1         POP_SPECIES(1,SPECIES_LNK(ID)), ATM(ID+1)%XzV_PRES, ATM(ID)%ZXzV,
	1         LOC_ID,TMP_STRING,OMEGA_GEN_V3,
	1         ATM(ID)%EQXzV,NUM_BNDS,ND,NION,COMPUTE_BA,DST,DEND)
!
! Handle states which can partially autoionize.
!
	      TMP_STRING=TRIM(ION_ID(ID))//'_AUTO_DATA'
              CALL STEQ_AUTO_V2(ED,T,
	1         ATM(ID)%XzV,        ATM(ID)%NXzV,         ATM(ID)%DXzV,
	1         ATM(ID)%XzV_F,      ATM(ID)%XzVLTE_F,     ATM(ID)%EDGEXzV_F, 
	1         ATM(ID)%GXzV_F,     ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F,
	1         ATM(ID)%F_TO_S_XzV, LOC_ID,               
	1         DIERECOM(1,ATM(ID)%INDX_XzV),             DIECOOL(1,ATM(ID)%INDX_XzV),
	1         TMP_STRING,         NUM_BNDS,ND,COMPUTE_BA,DST,DEND)
	    END IF
	  END DO
	  CALL TUNE(2,'STEQ')
!        
! Update charge equation. No longer done in STEQHEII
!
          CALL STEQNE_V4(ED,NT,DIAG_INDX,ND,COMPUTE_BA,DST,DEND)
!
	END DO		!K=1,ND  - DST,DEND
!
	IF(LST_ITERATION)
	1     CALL WR_ASCI_STEQ(NION,ND,'STEQ ARRAY- Collisional Terms',19)
! 
!
! Compute the collisional cooling terms for digestion.
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      TMP_STRING=TRIM(ION_ID(ID))//'_COL_DATA'
	      CALL COLCOOL_SL_V4(
	1        ATM(ID)%CPRXzV, ATM(ID)%CRRXzV,  ATM(ID)%COOLXzV,CNM,DCNM,
	1        ATM(ID)%XzV,    ATM(ID)%XzVLTE,  ATM(ID)%dlnXzVLTE_dlnT, 
	1        ATM(ID)%NXzV,   ATM(ID)%XzV_F,   ATM(ID)%XzVLTE_F, 
	1        ATM(ID)%AXzV_F, ATM(ID)%W_XzV_F, ATM(ID)%EDGEXzV_F,
	1        ATM(ID)%GXzV_F, ATM(ID)%XzVLEVNAME_F,
	1        ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%ZXzV,
	1        ID,TMP_STRING,OMEGA_GEN_V3,ED,T,ND)
	    END IF
	  END DO
!
! 
!
! Reread in BA array if we are not to compute it. There should be no problem
! with this read since it has previously been read in. This double reading
! is necessary to save a rewrite of the STEQ*** routines.
!
	IF(.NOT. COMPUTE_BA .AND. .NOT. FLUX_CAL_ONLY)THEN
	  CALL READ_BA_DATA_V2(LU_BA,NION,NUM_BNDS,ND,IOS,CHK,'BAMAT')
	  IF(.NOT. CHK)THEN
	    WRITE(LUER,*)'Major Error - cant read BA File'
	    WRITE(LUER,*)'Previously read successfully - '//
	1             'before continuum loop'
	    STOP
	  END IF
	END IF
! 
!
	EDDINGTON=EDD_LINECONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    ACCESS_F=6  		!First output record changed from
	    WRITE(LU_EDD,REC=1)ACCESS_F     !5 to 6 on 16-Jan-1992.
	  ELSE
	    READ(LU_EDD,REC=1)ACCESS_F
	  END IF
	END IF
!
	DO ML=1,NDIETOT
	  SECTION='DIELECTRONIC'
	  DO ID=1,NUM_IONS-1
	    IF(SPECDIE(ML) .EQ. ION_ID(ID))THEN
	      MNL_F=LEVDIE(ML)			!Level in full atom
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)		!Level in small_atom atom
	      NL=MNL+ATM(ID)%EQXzV-1			!Level in pops
	      GLOW=ATM(ID)%GXzV_F(LEVDIE(ML))
	      GION=ATM(ID)%GIONXzV_F
	      EQBAL=ATM(ID)%EQXzVBAL
	      EQION=ATM(ID+1)%EQXzV
	      EQSPEC=EQ_SPECIES(SPECIES_LNK(ID))
	      FL=ATM(ID)%EDGEXzV_F(LEVDIE(ML))-EDGEDIE(ML)
	      DO K=1,ND
		LOW_OCC_PROB(K)=ATM(ID)%W_XzV_F(MNL_F,K)		!Occupation prob.
		L_STAR_RATIO(K,1)=ATM(ID)%XzVLTE_F(MNL_F,K)/ATM(ID)%XzVLTE(MNL,K)
		dL_RAT_dT(K,1)=L_STAR_RATIO(K,1)*
	1         (-1.5D0-HDKT*ATM(ID)%EDGEXzV_F(MNL_F)/T(K)
	1               -ATM(ID)%dlnXzVLTE_dlnT(MNL,K))/T(K)
	      END DO
	      ID_SAV=ID
	      EXIT
	    END IF
	  END DO
!
	  COMPUTE_NEW_CROSS=.TRUE.
	  CONT_FREQ=FL
!
! Determine which method will be used to compute continuum intensity.
! Present form is temporary measure for consistency with SAO.
!
	  IF(ALL_FREQ)THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
	  GLDGU(1)=GLOW/GUPDIE(ML)
	  EINA(1)=EINADIE(ML)
	  OSCIL(1)=EINA(1)*EMLIN/( GLDGU(1)*OPLIN*TWOHCSQ*(FL**2) )
!
	  DO I=1,ND
	    DION(I)=POPS(EQION,I)
	  END DO
!
! Compute the LTE population for upper autoionizing state. We multiply
! NUST by Occupation probability of the lower state to correct for the
! fact that some transitions effectively keep the atom ionized.
!
	  CALL LTEPOP(NUST,ED,DION,GUPDIE(ML),EDGEDIE(ML),T,
	1               GION,IONE,ND)
	  DO I=1,ND
	    NUST(I)=NUST(I)*LOW_OCC_PROB(I)
	  END DO
!
! Compute line opacity and emissivity.
!
	  T1=OSCIL(1)*OPLIN
	  T2=FL*EINA(1)*EMLIN
	  DO I=1,ND
	    CHIL(I)=T1*( POPS(NL,I)*L_STAR_RATIO(I,1)-GLDGU(1)*NUST(I) )
	    ETAL(I)=T2*NUST(I)
	  END DO
! 
!
	  IF(IMPURITY_CODE)THEN
!
! Obtain previously compute continuum opacities, and mean intensities.
!
	    INCLUDE 'GET_J_CHI.INC'
	  ELSE
!
! Compute continuum opacity and emissivity at the line frequency.
!
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!	    INCLUDE 'OPACITIES_V4.INC'
!
! Solve for the continuous radiation field.
!
!	    INCLUDE 'COMP_JCONT_V4.INC'
            CALL COMP_J_BLANK(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1                              MAXCH,LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
	  END IF
!
! SOURCE is used by SOBJBAR and in VARCONT. Note that SOURCE is corrupted 
! (i.e. set to line source function) in CMFJBAR. Also note that SOURCEEXT 
! has previously been computed (not for impurity species).
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+RJ(I)*THETA(I)
	  END DO
!  
!
! We define VC(I)=N* d(N* Z)/dN* and VB(I)=d(N* Z)/dNL.
!
	  T3=1.0D0/( TWOHCSQ*(FL**3) )
	  T2=T3/GLDGU(1)
	  DO I=1,ND
	    ZNET(I)=1.0D0-RJ(I)*CHIL(I)/ETAL(I)
	    VC(I)=(1.0D0+RJ(I)*T3)*NUST(I)
	    VB(I)=-RJ(I)*T2
	  END DO
!
! Evaluate contribution to statistical equilibrium equation, and
! and increment variation matrices. 
! To convert cooling to ergs/cm**3/s, we use EDGEDIE(ML) and not FL for the 
! electron cooling component since this is the energy spent in exciting the
! autoionizing state. Can show from statistical equilibrium equations.
! Note that EDGEDIE(ML) is negative.
!
	  T2=-1.256637D-09*EINA(1)*EMLIN*EDGEDIE(ML)
	  T3=EINA(1)*FL*EMLIN
!
	  DO K=1,ND
	    SE(ID)%STEQ(NL,K)=SE(ID)%STEQ(NL,K)+EINA(1)*NUST(K)*ZNET(K)
	    STEQ_T(K)=STEQ_T(K)-T3*NUST(K)*ZNET(K)
	    DIELUM(K)=DIELUM(K)+ETAL(K)*ZNET(K)
	    DIERECOM(K,INDXDIE(ML))=DIERECOM(K,INDXDIE(ML))+
	1          EINA(1)*NUST(K)*ZNET(K)
	    DIECOOL(K,INDXDIE(ML))=DIECOOL(K,INDXDIE(ML))+
	1          T2*NUST(K)*ZNET(K)
	  END DO
!
! Update ionization balance equations if required.
!
	  MNUP=ATM(ID)%NXzV+1
	  DO K=1,ND
	    SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K)-EINA(1)*NUST(K)*ZNET(K)
	  END DO
! 
!
! Update BA matrix - this section must be done even if we are
! performing a LAMBDA iteration. We define a LAMBDA iteration by assuming
! that the variation of J is zero.
!
	  IF(COMPUTE_BA)THEN
!
	    MNUP=ATM(ID)%NXzV+1
	    NUP=ATM(ID)%EQXzV+ATM(ID)%NXzV
	    NIV=SE(ID)%N_IV
	    DO K=1,ND
	      L=GET_DIAG(K)
	      SE(ID)%BA(MNL,MNL,L,K) =SE(ID)%BA(MNL,MNL,L,K)  +EINA(1)*VB(K)
	      SE(ID)%BA(MNL,MNUP,L,K)=SE(ID)%BA(MNL,MNUP,L,K) +EINA(1)*VC(K)/DION(K)
	      SE(ID)%BA(MNL,NIV-1,L,K)=SE(ID)%BA(MNL,NIV-1,L,K) +EINA(1)*VC(K)/ED(K)
	      SE(ID)%BA(MNL,NIV,L,K)  =SE(ID)%BA(MNL,NIV,L,K)   -
	1          EINA(1)*VC(K)*(1.5D0+HDKT*EDGEDIE(ML)/T(K))/T(K) +
	1          EINA(1)*VB(K)*POPS(NL,K)*dL_RAT_dT(K,1)/L_STAR_RATIO(K,1)
!
	      BA_T(NL,L,K) =BA_T(NL,L,K)-T3*VB(K)
	      BA_T(NUP,L,K) =BA_T(NUP,L,K)-T3*VC(K)/DION(K)
	      BA_T(NT-1,L,K)=BA_T(NT-1,L,K)-T3*VC(K)/ED(K)
	      BA_T(NT,L,K)  =BA_T(NT,L,K)+
	1                   T3*VC(K)*(1.5D0+HDKT*EDGEDIE(ML)/T(K))/T(K)
	    END DO
!
! Update ionization equation.
!
	    DO K=1,ND
	      L=GET_DIAG(K)
	      SE(ID)%BA(MNUP,MNL,L,K)  =SE(ID)%BA(MNUP,MNL,L,K)   - EINA(1)*VB(K)
	      SE(ID)%BA(MNUP,MNUP,L,K) =SE(ID)%BA(MNUP,MNUP,L,K)  - EINA(1)*VC(K)/DION(K)
	      SE(ID)%BA(MNUP,NIV-1,L,K)=SE(ID)%BA(MNUP,NIV-1,L,K) - EINA(1)*VC(K)/ED(K)
	      SE(ID)%BA(MNUP,NIV,L,K    )=SE(ID)%BA(MNUP,NIV,L,K) +
	1         EINA(1)*VC(K)*(1.5D0+HDKT*EDGEDIE(ML)/T(K))/T(K) -
	1         EINA(1)*VB(K)*POPS(NL,K)*dL_RAT_dT(K,1)/L_STAR_RATIO(K,1)
	    END DO
	  END IF
! 
!
! Allow for the variation of the continuous radiation field.
!
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1                      .NOT. IMPURITY_CODE)THEN
!
	    DO I=1,ND
	      BETAC(I)=CHIL(I)/ETAL(I)
	    END DO
!
!	    INCLUDE 'VARCONT.INC'
            CALL DO_VAR_CONT(POPS,SECTION,EDDINGTON,
	1                  FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                  ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                  NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
!
! Increment the large simultaneous perturbation matrix due to a variation
! in the continuum. This is only incremented if ZNET is not within 1% of
! 1.0, which indicates that the continuum term is important.
!
	    CALL TUNE(IONE,'DIECONTBA')
	    DO L=1,ND	  	  		  	!S.E. equation depth
	      T1=EINA(1)*NUST(L)*BETAC(L)
	      T2=ETAL(L)*BETAC(L)
	      DO K=BNDST(L),BNDEND(L)	  		!Variable depth.
	        LS=BND_TO_FULL(K,L)
   	        DO J=1,SE(ID)%N_IV	 	   		!Variable
	          JJ=SE(ID)%LNK_TO_F(J)
	          SE(ID)%BA(MNL,J,K,L)=SE(ID)%BA(MNL,J,K,L) - T1*VJ(JJ,K,L)
	          BA_T(JJ,K,L)=BA_T(JJ,K,L) + T2*VJ(JJ,K,L)
	        END DO
	      END DO
	    END DO
!
	    DO L=1,ND	  	  		  	!S.E. equation depth
	      T1=EINA(1)*NUST(L)*BETAC(L)
	      DO K=BNDST(L),BNDEND(L)	  		!Variable depth.
	        LS=BND_TO_FULL(K,L)
   	        DO J=1,NT	 	   		!Variable
	          JJ=SE(ID)%LNK_TO_F(J)
	          SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L) + T1*VJ(JJ,K,L)
	        END DO
	      END DO
	    END DO
	    CALL TUNE(ITWO,'DIECONTBA')
	  END IF			!BA Matrix computed (compute_ba).
!
	  IF(LST_ITERATION)THEN
!
! Estimate the line EW using a Modified Sobolev approximation.
!
! We use TA as a temporary vector which indicates the origin
! of the line emission. Not required in this code as used only
! for display purposes. Variable after THK_CONT is true as we
! want to assume the line opacity is zero --- since dielectronic
! transition.
!
	    CALL SOBEW(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1              V,SIGMA,R,P,AQW,HQW,TA,EW,CONT_INT,
	1              FL,DIF,DBB,IC,THK_CONT,L_TRUE,NC,NP,ND,METHOD)
!
	    T1=LAMVACAIR(FL)			!Wavelength(Angstroms)
	    CALL EW_FORMAT(EW_STRING,DIENAME(ML),T1,CONT_INT,EW,L_TRUE)
	    L=ICHRLEN(EW_STRING)
	    WRITE(LU_NET,40002)EW_STRING(1:L)
	    WRITE(LU_DR,40002)EW_STRING(1:L)
	    WRITE(LU_EW,40005)EW_STRING(1:L)
	    WRITE(LU_HT,40002)EW_STRING(1:L)
	    WRITE(LU_NET,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_DR,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_HT,40009)GUPDIE(ML),FL,EINA(1)
	    WRITE(LU_NET,40003)( ZNET(I),I=1,ND )
	    WRITE(LU_DR,40003)(  ( ZNET(I)*NUST(I)*EINA(1) ),I=1,ND  )
	    WRITE(LU_HT,40003)(  ( ZNET(I)*ETAL(I) ),I=1,ND  )
40009	    FORMAT(1X,F5.0,2X,1P,2E12.4)
	    CLOSE(UNIT=LU_NET)
	    CALL GEN_ASCI_OPEN(LU_NET,'NETRATE','OLD','APPEND',' ',IZERO,IOS)
	  END IF
	END DO		!End of Dielectronic section [do ML]
                                !subsequent iterations.
! 
!***************************************************************************
!***************************************************************************
!
!                         CONTINUUM LOOP
!
!***************************************************************************
!***************************************************************************
!
	EDDINGTON=EDD_CONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    WRITE(LU_EDD,REC=EDD_CONT_REC)ACCESS_F,NCF,NDEXT
	  ELSE
	    READ(LU_EDD,REC=EDD_CONT_REC)ACCESS_F
	  END IF
	END IF
!
! Decide whether to use an file with old J values to provide an initial 
! estimate of J with incoherent electron scattering. The options
!
!             .NOT. RD_COHERENT_ES .AND. COHERENT_ES
!
! indicate that no ES_J_CONV file exists already.
! Note that TEXT and NDEXT contain T and ND when ACCURATE is FALSE.
!
	IF(USE_OLDJ_FOR_ES .AND. .NOT. RD_COHERENT_ES .AND. COHERENT_ES)THEN
	  COHERENT_ES=RD_COHERENT_ES
	  I=SIZE(VJ)
	  CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LUIN,'OLD_J_FILE',
	1           EDD_CONT_REC,L_FALSE,L_TRUE,LU_ES,'ES_J_CONV')
!
! Now open the file so it can be read in the CONTINUUM loop (read in 
! COMP_JCONT).
!
	  OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening ES_J_CONV'//
	1                    ' - will compute new J'
	      COHERENT_ES=.TRUE.
	    END IF
	END IF
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
!
! Ensure none of the storage location for the variation of J with CHIL etc
! are being pointed at.
!
	DO SIM_INDX=1,MAX_SIM
	  LOW_POINTER(SIM_INDX)=0
	  UP_POINTER(SIM_INDX)=0
	END DO
!
	DO I=1,NM
	  VAR_IN_USE_CNT(I)=0
	  VAR_LEV_ID(I)=0
	  IMP_TRANS_VEC(I)=.FALSE.
	END DO
!
	NUM_OF_WEAK_LINES=0.0D0
	CONT_FREQ=0.0D0
!                                                                    
! Enter loop for each continuum frequency.
!
	CALL TUNE(IONE,'MLCF')
	CALL TUNE(IONE,'10000')
	DO 10000 ML=1,NCF
	  FREQ_INDX=ML
	  FL=NU(ML)
	  IF(ML .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	  END IF
	  SECTION='CONTINUUM'
!
	  IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	    COMPUTE_NEW_CROSS=.TRUE.
	    CONT_FREQ=NU_EVAL_CONT(ML)
	  ELSE
	    COMPUTE_NEW_CROSS=.FALSE.
	  END IF
	  FINAL_CONSTANT_CROSS=.TRUE.
	  IF(ML .EQ. NCF)THEN
	    FINAL_CONSTANT_CROSS=.TRUE.
	  ELSE
	    IF(NU_EVAL_CONT(ML+1) .EQ. CONT_FREQ)
	1                      FINAL_CONSTANT_CROSS=.FALSE.
	  END IF
!
! Compute quadrature weights for statistical equilibrium equations.
! TA is used as a work vector (dim ND)
!
	  CALL TUNE(IONE,'QUAD')
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES .AND. COMPUTE_NEW_CROSS)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	       PHOT_ID=J
	       CALL QUAD_MULTI_V7(ATM(ID)%WSXzV(1,1,J), ATM(ID)%dWSXzVdT(1,1,J),
	1             ATM(ID)%WCRXzV(1,1,J),
	1             ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1             ATM(ID)%XzVLTE_F, ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F,
	1             ATM(ID)%F_TO_S_XzV,CONT_FREQ,T,ND,
	1             ION_ID(ID), ATM(ID)%ZXzV, PHOT_ID, ID)
	      END DO
	    END IF  
	  END DO
	  CALL TUNE(ITWO,'QUAD')
!
! 
!
	IF(XRAYS .AND. COMPUTE_NEW_CROSS)THEN
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      T1=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV		!Number of electrons
	      CALL QUAD_X_GEN_V4(AT_NO(SPECIES_LNK(ID)),T1,
	1           ATM(ID)%WSE_X_XzV,   ATM(ID)%WCR_X_XzV, CONT_FREQ,
	1           ATM(ID)%XzVLTE,      ATM(ID)%NXzV,
	1           ATM(ID)%XzVLTE_F,    ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%F_TO_S_XzV,  ATM(ID)%NXzV_F,
	1           ATM(ID+1)%EDGEXzV_F, ATM(ID+1)%NXzV_F, ND)
	    END IF
	  END DO
	END IF
!
! 
!
! Include lines 
!
        CALL SET_LINE_OPAC(POPS,NU,ML,LAST_LINE,N_LINE_FREQ,
	1         LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
!
	CALL INIT_LINE_OPAC_VAR(LAST_LINE,LUER,ND,TX_OFFSET,MAX_SIM,NM)
!
! Determine which method will be used to compute continuum intensity.
!
	  IF(ACCURATE .AND. ALL_FREQ)THEN       
	    THIS_FREQ_EXT=.TRUE.
	  ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE        
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
	  IF(IMPURITY_CODE)THEN
!
! Obtain previously compute continuum opacities, and mean intensities.
!
	    INCLUDE 'GET_J_CHI.INC'
	  ELSE
!
! Compute opacity and emissivity.
!
	    CALL TUNE(IONE,'C_OPAC')
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONlY)
!	    INCLUDE 'OPACITIES_V4.INC'
	    CALL TUNE(ITWO,'C_OPAC')
!
! Since resonance zones included, we must add the line opacity and 
! emissivity to the raw continuum values. We first save the pure continuum 
! opacity and emissivity. These are used in carrying the variation of J from 
! one frequency to the next.
!
	    DO SIM_INDX=1,MAX_SIM
	      IF(RESONANCE_ZONE(SIM_INDX))THEN
	        DO I=1,ND
	          CHI(I)=CHI(I)+CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(SIM_INDX)
	          ETA(I)=ETA(I)+ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(SIM_INDX)
	        END DO
	      END IF
	    END DO
!
!	    DO I=1,ND
!	      IF(CHI(I)*R(I) .LT. 1.0D-04)THEN
!	         CHI(I)=CHI(I)-CHI_SCAT(I)
!	         CHI_SCAT(I)=1.0D-04/R(I)
!	         CHI(I)=CHI(I)+CHI_SCAT(I)
!	      END IF
!	    END DO
!
! CHECK for negative line opacities. NEG_OPAC_FAC is the factor we
! multiply the line opacities by so that the total opacity is positive.
! We do not distinguish between lines. Two different options are possible.
! The first, 'ESEC_CHK' was in use for years, and is probably the preferred
! option. The second, 'SRCE_CHK', was introduced to overcome problems in
! O Star models. In particular, in some models a negative optical depth
! could occur on some iteartions at depths where ABS(TAUL) was still very
! large (primraily in far IT transitions [e.g. H(9-8)].
!
	    AT_LEAST_ONE_NEG_OPAC=.FALSE.
	    NEG_OPACITY(1:ND)=.FALSE.
	    NEG_OPAC_FAC(1:ND)=1.0D0
	    IF(NEG_OPAC_OPTION .EQ. 'SRCE_CHK')THEN
	      DO I=1,ND
	        IF(CHI(I) .LT. CHI_CONT(I) .AND.
	1            CHI(I) .LT. 0.1D0*ETA(I)*CHI_NOSCAT(I)/ETA_CONT(I) )THEN
	          CHI(I)=0.1D0*ETA(I)*CHI_NOSCAT(I)/ETA_CONT(I)
	          NEG_OPACITY(I)=.TRUE.
	          NEG_OPAC_FAC(I)=0.0D0
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        ELSE IF(CHI(I) .LT. 0.1D0*CHI_SCAT(I))THEN
	          CHI(I)=0.1D0*CHI_SCAT(I)
	          NEG_OPACITY(I)=.TRUE.
	          NEG_OPAC_FAC(I)=0.0D0
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        END IF
	      END DO
	    ELSE IF(NEG_OPAC_OPTION .EQ. 'ESEC_CHK')THEN
	      DO I=1,ND
	        IF(CHI(I) .LT. 0.1D0*CHI_SCAT(I))THEN
	          T1=CHI(I)
	          CHI(I)=0.1D0*CHI_SCAT(I)
	          NEG_OPACITY(I)=.TRUE.
!	          NEG_OPAC_FAC(I)=(CHI(I)-CHI_CONT(I))/(T1-CHI_CONT(I))
	          NEG_OPAC_FAC(I)=0.0D0
	          AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        END IF
	      END DO
	    END IF
!
	    IF(LST_ITERATION .AND. AT_LEAST_ONE_NEG_OPAC)THEN
	      WRITE(LU_NEG,'(A,1P,E14.6)')
	1         ' Neg opacity for transition for frequency ',FL
	      DO SIM_INDX=1,MAX_SIM
	        IF(RESONANCE_ZONE(SIM_INDX))THEN
	          WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)
	        END IF
	      END DO
	      J=0
	      K=0
	      DO I=1,ND
	       IF(NEG_OPACITY(I) .AND. K .EQ. 0)K=I
	       IF(NEG_OPACITY(I))J=I
	      END DO
	      WRITE(LU_NEG,'(A,2X,I3,5X,A,2XI3)')
	1        ' 1st depth',K,'Last depth',J
	    END IF
!
	    DO I=1,ND
	      ZETA(I)=ETA(I)/CHI(I)
	      THETA(I)=CHI_SCAT(I)/CHI(I)
	    END DO
!
	    IF(LST_ITERATION .AND. ML .NE. NCF)THEN
	      DO I=1,N_TAU_EDGE
	        IF(NU(ML) .GE. TAU_EDGE(I) .AND. 
	1                       NU(ML+1) .LT. TAU_EDGE(I))THEN
	          T1=LOG(CHI_CONT(5)/CHI_CONT(1))/LOG(R(1)/R(5))
	          WRITE(LUER,'(A,1P,E11.4,A,E10.3)')' Tau(Nu=',NU(ML),
	1            ') at outer boundary is:',CHI_CONT(1)*R(1)/MAX(T1-1.0D0,1.0D0)
	        END IF
	      END DO
	    END IF
!
! Compute continuum intensity.
!
	    CALL TUNE(IONE,'COMP_J')
!	    INCLUDE 'COMP_JCONT_V4.INC'	
            CALL COMP_J_BLANK(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1                              MAXCH,LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
	    CALL TUNE(ITWO,'COMP_J')
	  END IF
! 
!
! Increment the RADIATIVE EQUILIBRIUM equation due to radiation field at this
! frequency. The correction for non-coherent electron scattering allows
! for the fact that the electron-scattering emissivity is ESEC*RJ_ES, not
! ESEC*RJ.
!
! ETA_CONT includes all emissivity sources, including X-ray emission produced
! by mechanical or magnetic energy deposition. This should not be included
! in the radiatively equilibrium equation, hence we subtract out the
! emissivity due to mechanical processes. NB: ETA_NOSCAT does not include
! mechanical term.
!
	  DO K=1,ND
	    STEQ_T(K)=STEQ_T(K)+ FQW(ML)*(CHI_NOSCAT(K)*RJ(K) - ETA_NOSCAT(K))
	  END DO
	  IF(.NOT. COHERENT_ES)THEN
	    STEQ_T(:)=STEQ_T(:)+FQW(ML)*ESEC(:)*(RJ(:)-RJ_ES(:))
	  END IF
!
	  CALL COMP_VAR_JREC(JREC,dJRECdT,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1       RJ,EMHNUKT,T,NU(ML),FQW(ML),TWOHCSQ,HDKT,ND,COMPUTE_NEW_CROSS)
!
! Increment the S.E. equations due to radiation field at this
! frequency.
!
! At the same time, we compute the quadrature weights associated with 
! the intensity for each depth point and each equation. QFV must be zeroed 
! before calling EVALSE_QWVJ. QFV is incremented - not set. Allows for 
! bound-free processes to both the ground and excited states (necessary 
! for CIII).
!
	IF(FINAL_CONSTANT_CROSS)THEN
	  DO ID=1,NION
	    IF(SE(ID)%XzV_PRES)THEN
	      SE(ID)%QFV_R(:,:)=0.0D0		!NT,ND
	      SE(ID)%QFV_P(:,:)=0.0D0
	    END IF
	  END DO
	END IF
!
	IF(FINAL_CONSTANT_CROSS)THEN
	  CALL TUNE(IONE,'EVALSE')
	  DO ID=1,NUM_IONS-1
	    ID_SAV=ID
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        CALL EVALSE_QWVJ_V6(ID_SAV,
	1         ATM(ID)%WSXzV(1,1,J), ATM(ID)%XzV, ATM(ID)%XzVLTE,
	1         ATM(ID)%NXzV, ATM(ID)%XzV_ION_LEV_ID(J),
	1         ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE, ATM(ID+1)%NXzV, 
	1         JREC,JPHOT,NT,ND)
	      END DO
	    END IF
	  END DO
	  CALL TUNE(ITWO,'EVALSE')
	END IF
!
! 
! Note that ATM(ID+2)%EQXzV is the ion equation. Since Auger ionization,
! 2 electrons are ejected.
!
	IF(XRAYS .AND. FINAL_CONSTANT_CROSS)THEN
	  DO ID=1,NUM_IONS-1
	    ID_SAV=ID
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      CALL EVALSE_X_QWVJ_V4(ID_SAV,ATM(ID)%WSE_X_XzV,
	1          ATM(ID)%XzV,   ATM(ID)%XzVLTE,   ATM(ID)%NXzV,
	1          ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE, ATM(ID+1)%NXzV, ATM(ID+2)%EQXzV,
	1          JREC,JPHOT,ND,NION)
	    END IF
	  END DO
	END IF
!
! 
!
! Compute the recombination, photoionization and cooling rates.
!
	IF(LST_ITERATION .AND. FINAL_CONSTANT_CROSS)THEN
	  CALL TUNE(IONE,'PRRRCOOL')
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        CALL PRRR_SL_V5(
	1          ATM(ID)%APRXzV,        ATM(ID)%ARRXzV, 
	1          ATM(ID)%BFCRXzV,      ATM(ID)%FFXzV,
	1          ATM(ID)%WSXzV(1,1,J), ATM(ID)%WCRXzV(1,1,J),
	1          ATM(ID)%XzV,          ATM(ID)%XzVLTE, 
	1          ATM(ID)%NXzV,         ATM(ID)%ZXzV,
	1          ATM(ID+1)%XzV,        ATM(ID+1)%XzVLTE, 
	1          ATM(ID+1)%NXzV, J,    ATM(ID)%XzV_ION_LEV_ID(J),
	1          ED,T,JREC,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1          FL,CONT_FREQ,ZERO_REC_COOL_ARRAYS,ND)
	      END DO
	    END IF
!
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES .AND. XRAYS)THEN
	      CALL X_RRR_COOL_V5(X_RECOM(1,ATM(ID)%INDX_XzV),
	1            X_COOL(1,ATM(ID)%INDX_XzV), ATM(ID)%WSE_X_XzV,
	1            ATM(ID)%WCR_X_XzV,
	1            ATM(ID)%XzV,        ATM(ID)%XzVLTE,     ATM(ID)%NXzV,
	1            ATM(ID+1)%XzV_F,    ATM(ID+1)%XzVLTE_F, ATM(ID+1)%NXzV_F,
	1            JREC,JPHOT,JREC_CR,JPHOT_CR,
	1            ZERO_REC_COOL_ARRAYS,ND,L_TRUE)
	    END IF
	  END DO
	  ZERO_REC_COOL_ARRAYS=.FALSE.
!
	CALL TUNE(ITWO,'PRRRCOOL')
	END IF 			!Only evaluate if last iteration.
!
! 
!
! Update line net rates, and the S.E. Eq. IFF we have finished a line 
! transition.
!                      
	DO SIM_INDX=1,MAX_SIM
	  IF(RESONANCE_ZONE(SIM_INDX))THEN
	    DO I=1,ND
	      ZNET_SIM(I,SIM_INDX)=ZNET_SIM(I,SIM_INDX) +
	1          LINE_QW_SIM(SIM_INDX)*
	1          (1.0D0-RJ(I)*CHIL_MAT(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX))
	      JBAR_SIM(I,SIM_INDX)=JBAR_SIM(I,SIM_INDX) +
	1          LINE_QW_SIM(SIM_INDX)*RJ(I)
	    END DO
	  END IF
	END DO
!
! Update the S.E. Eq. IFF we have finished a line transition (i.e. are at
! the final point of the resonance zone.) 
!
! NB: The line term in the RE equations is not needed since it is included 
! directly with continuum integration.
!
	DO SIM_INDX=1,MAX_SIM
	  IF( END_RES_ZONE(SIM_INDX) )THEN
            T1=FL_SIM(SIM_INDX)*EMLIN 
	    NL=SIM_NL(SIM_INDX)
	    NUP=SIM_NUP(SIM_INDX)
	    I=SIM_LINE_POINTER(SIM_INDX)
	    ID=VEC_ID(I)
	    MNL_F=VEC_MNL_F(I);     MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP_F=VEC_MNUP_F(I);   MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    IF(SCL_LINE_COOL_RATES)THEN
	      T3=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	      IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
	    ELSE
	      T3=1.0D0
	    END IF
	    DO K=1,ND					!Equation depth
	      IF(POP_ATOM(K) .GE. SCL_LINE_DENSITY_LIMIT .AND. K .NE. ND+1)T3=1.0D0
	      T2=ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	      SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K) - T2/T1
	      SE(ID)%STEQ(MNL,K) =SE(ID)%STEQ(MNL,K) + T2/T1
	      STEQ_T(K)=STEQ_T(K) - T2*T3
	    END DO
	  END IF    
	END DO
!
!                                                                    
!
! Allow for the variation of the continuous radiation field.
!
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1                     .NOT. IMPURITY_CODE)THEN
!
! Solve for the perturbations to J in terms of the perturbations
! to CHI and ETA. 
!            
	    CALL TUNE(IONE,'C_VARCONT')
!	      INCLUDE 'VARCONT.INC'
              CALL DO_VAR_CONT(POPS,SECTION,EDDINGTON,
	1                    FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                    ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                    NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
	    CALL TUNE(ITWO,'C_VARCONT')
!
! NB: VJ, VCHI, and VETA must not be modified until we have updated the
!     BA array.
!
	  END IF
!
! 
!
! Modify the BA matrix for terms in the statistical equilibrium
! equations which are multiplied by RJ. NB. This is not for the
! variation of RJ - rather the multiplying factors. This section
! must be done for a LAMBDA iteration.
!
! We use DST and DEND to avoid reading in the entire diagonal of the BA
! array for each call to BA. Calling the routines ND times may take
! too long --- therefore try alternative method.
!
! WE have replaced BA and BAION by BA_PAR and BAION_PAR in calls to
! VSEBYJ. Since these are diagonal components, we have had to replace 
! NUM_BANDS by IONE in the calls.
!
	IF(COMPUTE_BA)THEN
	  CALL TUNE(IONE,'COMPUTE_BA')
!
!	  DO K=1,ND
!	    DST=K
!	    DEND=K
	  DO K=1,1
	    DST=1       
	    DEND=ND
!
! NB: In this equation the matrices are NOT passed for WS, dWS and NU.
!
	    IF(FINAL_CONSTANT_CROSS)THEN
	      DO ID=1,NUM_IONS-1
	        ID_SAV=ID
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO J=1,ATM(ID)%N_XzV_PHOT
	            CALL VSEBYJ_MULTI_V6(ID_SAV,
	1             ATM(ID)%WSXzV(1,1,J), ATM(ID)%dWSXzVdT(1,1,J),
	1             ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, 
	1             ATM(ID)%NXzV,
	1             ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE,
	1             ATM(ID+1)%dlnXzVLTE_dlnT, ATM(ID+1)%NXzV,
	1             ATM(ID)%XzV_ION_LEV_ID(J),ED,T,
	1             JREC,dJRECdt,JPHOT,NUM_BNDS,ND,DST,DEND)
	          END DO
	        END IF
	      END DO
	    END IF
!
! 
! Note that ATM(ID+2)%EQXzV is the ion equation. Since Auger ionization,
! 2 electrons are ejected.
!
	    IF(XRAYS .AND. FINAL_CONSTANT_CROSS)THEN
	      DO ID=1,NUM_IONS-1
	        IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	          CALL VSEBYJ_X_V6(ID,ATM(ID)%WSE_X_XzV,
	1             ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1             ATM(ID+1)%XzV_F, ATM(ID+1)%XzVLTE_F, ATM(ID+1)%EDGEXzV_F,
	1             ATM(ID+1)%NXzV_F,ATM(ID+1)%DXzV,    ATM(ID+2)%EQXzV,
	1             ED,T,JREC,dJRECdT,JPHOT,
	1             ND,NION,DST,DEND)
	        END IF
	      END DO
	    END IF		!End X-rays
!
! 
!
! Increment the large simultaneous perturbation matrix due to the
! variation in J. We generally update the _PAR matrices. This is done to
! ensure better cancelation of terms similar in size as occurs when
! the optical depth is large.
!
! Every N_PAR frequencies we update the the full matrices, and zero the
! ?_PAR arrays. Note that the non-diagonal part of BA and BAION
! (i.e BA(I,J,K,L) with L .NE. DIAG_INDX) are presently updated for
! every frequency.
!
! We now pass the continuum emissivity and opacity without any scattering 
! contribution. We use TA as a zeroed vec for ESEC, which was subtracted
! from CHI_CONT inside BA_UPDATE_V7.
!
	    IF( .NOT. LAMBDA_ITERATION)THEN
	      TA(1:ND)=0.0D0
              CALL BA_UPDATE_V7(VJ,VCHI_ALL,VETA_ALL,
	1             ETA_NOSCAT,CHI_NOSCAT,TA,T,POPS,RJ,NU(ML),FQW(ML),
	1             COMPUTE_NEW_CROSS,FINAL_CONSTANT_CROSS,DO_SRCE_VAR_ONLY,
	1             BA_CHK_FAC,NION,NT,NUM_BNDS,ND,DST,DEND)
	    END IF
!
	  END DO		!K=1,ND : DST,DEND
	  CALL TUNE(ITWO,'COMPUTE_BA')
!
	  CALL TUNE(IONE,'ADD_PAR')
	  IF( MOD(FREQ_INDX,N_PAR) .EQ. 0 .OR. FREQ_INDX .EQ. NCF )THEN
            CALL ADD_PAR_TO_FULL_V2(NION,DIAG_INDX)
  	  END IF
	  CALL TUNE(ITWO,'ADD_PAR')
!
	END IF			!End compute_ba
!
! 
!
! We now allow for the variation of the LINE terms in the S.E. equations.
! In this first section we perform a LAMBDA iteration for the LINE terms.
! We do this if we are doing a FULL Lambda iteration, or if the line is
! being treated as a WEAK_LINE.
!
	CALL TUNE(1,'BA_LAM')
	DO SIM_INDX=1,MAX_SIM
	  IF( COMPUTE_BA .AND. (LAMBDA_ITERATION .OR. WEAK_LINE(SIM_INDX)) )THEN
	    IF (END_RES_ZONE(SIM_INDX) )THEN
!
! If we are doing a lambda iteration, we are assuming that JBAR is fixed.
! Thus, dZ/dCHIL and dZ/dETAL is given by the following.
!
! NB: We do not set VB(I) to if NEG_OPACITY(I)=.TRUE. as we have not altered
!     CHIL : We have only changed CHI which effects JBAR only.
!
! Altered 15-Feb-2005: No longer need VC. We now xplicitly cancel NU/ETAL depdence with
! Nu.
!
	      DO I=1,ND
	        VB(I)=-JBAR_SIM(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX)
!	        VC(I)=JBAR_SIM(I,SIM_INDX)*CHIL_MAT(I,SIM_INDX) /
!	1                       ETAL_MAT(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX)
	      END DO
!
	      OPAC_FAC=OSCIL(SIM_INDX)*OPLIN
	      STIM_FAC=-GLDGU(SIM_INDX)*OPAC_FAC
	      EMIS_FAC=EINA(SIM_INDX)*FL_SIM(SIM_INDX)*EMLIN
	      T4=FL_SIM(SIM_INDX)*EMLIN
	      NUP=SIM_NUP(SIM_INDX)
 	      NL=SIM_NL(SIM_INDX)
!
	      I=SIM_LINE_POINTER(SIM_INDX)
	      ID=VEC_ID(I)
	      MNL_F=VEC_MNL_F(I)
	      MNUP_F=VEC_MNUP_F(I)
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
	      IF(SCL_LINE_COOL_RATES)THEN
	        T3=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	        IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
	        T4=T4*T3
	      END IF
	      DO K=1,ND
	        IF(POP_ATOM(K) .GE. SCL_LINE_DENSITY_LIMIT .AND. K .NE. ND+1)T4=FL_SIM(SIM_INDX)*EMLIN
	        dRATE_dUP=EINA(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)*
	1            (1.0D0+U_STAR_RATIO(K,SIM_INDX)*POPS(NUP,K)*STIM_FAC*VB(K))
	        dRATE_dLOW=EINA(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)*
	1              L_STAR_RATIO(K,SIM_INDX)*OPAC_FAC*POPS(NUP,K)*VB(K)
	        L=GET_DIAG(K)
	        SE(ID)%BA(MNUP,MNUP,L,K)=SE(ID)%BA(MNUP,MNUP,L,K)-dRATE_dUP
	        SE(ID)%BA(MNUP,MNL,L,K) =SE(ID)%BA(MNUP,MNL,L,K) -dRATE_dLOW
	        SE(ID)%BA(MNL,MNUP,L,K) =SE(ID)%BA(MNL,MNUP,L,K) +dRATE_dUP
	        SE(ID)%BA(MNL,MNL,L,K)  =SE(ID)%BA(MNL,MNL,L,K)  +dRATE_dLOW
	        BA_T(NL,L,K) =BA_T(NL,L,K) -T4*dRATE_dLOW
	        BA_T(NUP,L,K)=BA_T(NUP,L,K)-T4*dRATE_dUP
	      END DO
	    END IF	!Resonance zone
	  END IF	!Lambda/Weak line
	END DO		!Loop over SIM_INDX
	CALL TUNE(2,'BA_LAM')
!
! 
	CALL TUNE(IONE,'LINE_BA')
	IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	  IF(NUM_BNDS .EQ. ND)THEN
	     WRITE(LUER,*)'Error --- NUM_BNDS .EQ. ND not implemented'
	     STOP
	  END IF
!
! Update the dZ matrix which gives the variation of the Net Rate
! (=1-JBAR/SL) with respect to eta, chi, ED, ... etc.
!
! dZ( X, X depth [1:NUM_BNDS], Z depth, N_SIM)
!
! Scaling MAX_SIM*ND*NUM_BNDS*NM
!
! NB: J goes from 1 to NUM_BNDS (and not BNDST(K),BNDEND(K)) since we have no
! reference to elements in a full matrix (such as FCHI).
!
! NB: OPAC_FAC and STIM_FAC are not multiplied by NEG_OPAC_FAC since we use
!     the uncorrected line opacity in the source function and hence the
!     expression for ZNET.
!
	  CALL TUNE(IONE,'dZ_LINE')
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX) .AND. .NOT. WEAK_LINE(SIM_INDX))THEN
	      DO K=1,ND
	        MUL_FAC=LINE_QW_SIM(SIM_INDX)*CHIL_MAT(K,SIM_INDX)/
	1                   ETAL_MAT(K,SIM_INDX)
	        OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*L_STAR_RATIO(K,SIM_INDX)
	        EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)
	        STIM_FAC=LINE_OPAC_CON(SIM_INDX)*
	1                        U_STAR_RATIO(K,SIM_INDX)*GLDGU(SIM_INDX)
	        DO J=BNDST(K),BNDEND(K)
	          IF(J .NE. DIAG_INDX)THEN
	            DO L=3,NM
	              IF(DO_THIS_TX_MATRIX(L))THEN
	                dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                        MUL_FAC*dJ_LOC(L,J,K)
	              END IF
	            END DO
	          ELSE 
	            DO L=3,NM
	              IF(DO_THIS_TX_MATRIX(L))THEN
	                IF( L .EQ. LOW_POINTER(SIM_INDX) )THEN
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                    MUL_FAC*( dJ_LOC(L,J,K) +
	1                          RJ(K)*OPAC_FAC/CHIL_MAT(K,SIM_INDX) )
	                ELSE IF( L .EQ. UP_POINTER(SIM_INDX) )THEN
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                    MUL_FAC*( dJ_LOC(L,J,K)-
	1                    RJ(K)*STIM_FAC/CHIL_MAT(K,SIM_INDX)-
	1                    RJ(K)*EMIS_FAC/ETAL_MAT(K,SIM_INDX) )
	                ELSE
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                      MUL_FAC*dJ_LOC(L,J,K)
	                END IF
	              END IF
	            END DO
	          END IF
	        END DO
	      END DO
	    END IF 
	  END DO
	  CALL TUNE(ITWO,'dZ_LINE')
!
! 
!
! Check whether any lines have been finalized. If so we can update the 
! variation matrices.
!
! NB: Weak lines are handled in the LAMBDA ITERATION section.
!
! We use ZENT_VAR_LIMIT to reduce the number of operations. When ZNET
! is approximatley unity, there is little point in linearizing ZNET.
! NB: An ideal vale for ZNET_VAR_LIMIT is probably 0.01 or 0.001. If 
! ZNET_VAR_LIMIT is zero, all depths will be included in the linearization, 
! independent of ZNET. A very large value of ZNET (i.e. 10^4), will imply
! an interation on the NET_RATES, with no linearization.
!
	  DO SIM_INDX=1,MAX_SIM
!
	    IF( END_RES_ZONE(SIM_INDX) .AND. .NOT. WEAK_LINE(SIM_INDX))THEN
!
	      I=SIM_LINE_POINTER(SIM_INDX)
	      ID=VEC_ID(I)
!
! Can now update rate equation for variation in Z.
!
	      I=NT*NUM_BNDS*ND
	      CALL DP_ZERO(dZ_POPS,I)
!
! Check which matrices get used to compute dZ_POPS. We include:
!   (i) all transitions designated to be important.
!  (ii) all transitions for levels in the same ionization stage.
!
	      DO I=TX_OFFSET+1,NM
	        USE_THIS_VAR_MAT(I)=.FALSE.
	        IF(VAR_IN_USE_CNT(I) .GT. 0 .AND. IMP_TRANS_VEC(I))THEN
	          USE_THIS_VAR_MAT(I)=.TRUE.
	        ELSE IF(VAR_IN_USE_CNT(I) .GT. 0)THEN
	          NL=VAR_LEV_ID(I)
	          IF(NL .GE. ATM(ID)%EQXzV .AND.
	1                      NL .LE. ATM(ID)%EQXzV+ATM(ID)%NXzV)
	1         USE_THIS_VAR_MAT(I)=.TRUE.
	        END IF
	      END DO
!
! NOPS= ND*NUM_BNDS*( 4NT + 2 + 7NUM_SIM )
!
!
! Set up a temporary vector to handle Rayleigh scattering.
!
	      TA(1:ND)=0.0D0
	      IF(SPECIES_PRES(1))THEN			!If H present!
	        DO L=1,ND
	          TA(L)=CHI_RAY(L)/ATM(1)%XzV_F(1,L)	          
	        END DO
	      END IF
!
	      DO K=1,ND
		T4=ABS(ZNET_SIM(K,SIM_INDX)-1.0D0)
	          IF(T4 .GT. ZNET_VAR_LIMIT)THEN
	          DO J=BNDST(K),BNDEND(K)	  	  !Since refer to VCHI etc.
 	            L=BND_TO_FULL(J,K)
	            DO JJ=1,SE(ID)%N_IV			!Bad notation
	              I=SE(ID)%LNK_TO_F(JJ)
	              dZ_POPS(I,J,K)=dZ_POPS(I,J,K) +
	1                ( VCHI(I,L)*dZ(3,J,K,SIM_INDX) + 
	1                    VETA(I,L)*dZ(4,J,K,SIM_INDX) )
	            END DO
	            dZ_POPS(NT-1,J,K)=dZ_POPS(NT-1,J,K) +
	1                            ESEC(L)*dZ(5,J,K,SIM_INDX)/ED(L)
	            dZ_POPS(1,J,K)=dZ_POPS(1,J,K) +
	1                            TA(L)*dZ(5,J,K,SIM_INDX)
!
! Now must do line terms.
!
	            DO I=TX_OFFSET+1,NM
	              IF(USE_THIS_VAR_MAT(I))THEN
	                NL=VAR_LEV_ID(I)
	                dZ_POPS(NL,J,K)=dZ_POPS(NL,J,K) + dZ(I,J,K,SIM_INDX)
	              END IF
	            END DO
	          END DO	        !Over Variable depth (1:NUM_BNDS)
	        END IF		!ABS|ZNET-1| > 0.01
	      END DO			!Over J depth.
!
! Update variation equations. NB. We do not need to update the Radiative
! Equilibrium equation, since this is automatically updated with the
! continuum.
!
! Note that there is no update of the NT equation since this is 
! automatically included in the direct CHI*J-ETA term.
!
! NOPS = 6*ND*NT*NUM_BNDS
!
	      CALL TUNE(IONE,'dBA_LINE')
	      NL=SIM_NL(SIM_INDX)
	      NUP=SIM_NUP(SIM_INDX)
	      I=SIM_LINE_POINTER(SIM_INDX)
              MNL_F=VEC_MNL_F(I)
	      MNUP_F=VEC_MNUP_F(I)
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
	      IF(SCL_LINE_COOL_RATES)THEN
	        T3=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	        IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
	      ELSE
	        T3=1.0D0
	      END IF
	      DO L=1,ND					!Equation depth
	        IF(POP_ATOM(L) .GE. SCL_LINE_DENSITY_LIMIT .AND. L .NE. ND+1)T3=1.0D0
	        T1=EINA(SIM_INDX)*U_STAR_RATIO(L,SIM_INDX)*POPS(NUP,L)
	        K=GET_DIAG(L)
	        SE(ID)%BA(MNL,MNUP,K,L)=SE(ID)%BA(MNL,MNUP,K,L) +
	1           EINA(SIM_INDX)*U_STAR_RATIO(L,SIM_INDX)*
	1             ZNET_SIM(L,SIM_INDX)
	        SE(ID)%BA(MNUP,MNUP,K,L)=SE(ID)%BA(MNUP,MNUP,K,L) -
	1           EINA(SIM_INDX)*U_STAR_RATIO(L,SIM_INDX)*
	1            ZNET_SIM(L,SIM_INDX)
	        BA_T(NUP,K,L)=BA_T(NUP,K,L) -
	1           T3*ETAL_MAT(L,SIM_INDX)*ZNET_SIM(L,SIM_INDX)/POPS(NUP,L)
		T4=ABS(ZNET_SIM(L,SIM_INDX)-1.0D0)
	        IF(T4 .GT. ZNET_VAR_LIMIT)THEN
	          DO K=1,NUM_BNDS			!Variable depth
!
! Because of our corrections earlier, unimportant variations are also
! omitted from the radiative Equilibrium Equation.
!
	            DO J=1,SE(ID)%N_IV
	              JJ=SE(ID)%LNK_TO_F(J)
	              SE(ID)%BA(MNL,J,K,L) =SE(ID)%BA(MNL,J,K,L) +dZ_POPS(JJ,K,L)*T1
	              SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L)-dZ_POPS(JJ,K,L)*T1
	              BA_T(JJ,K,L)=BA_T(JJ,K,L)-dZ_POPS(JJ,K,L)*
	1                                        T3*ETAL_MAT(L,SIM_INDX)
	            END DO
	          END DO
	        END IF
	      END DO
	      CALL TUNE(ITWO,'dBA_LINE')
!
! Must now zero dZ since next time it is used it will be for a new line.
!
	      I=NM*NUM_BNDS*ND
	      CALL DP_ZERO(dZ(1,1,1,SIM_INDX),I)
!
	    END IF	!End_res_zone .and. .not. weak_line
	  END DO	!SIM_INDX
	END IF		!Compute_ba .and. .not. lambda_iteration
	CALL TUNE(ITWO,'LINE_BA')
!
! 
!
! Free up LINE storage locations. Removal is done in 3 ways, but only 2 here.
! Recall that frequencies are ordered from highest to lowest, and that we
! integrate from blue to red.
!
! 1. Lines interact over at most 2Vinf from last point of resonance zone. Thus
!      when the current frequency is 2Vinf (converted to frequency units)
!      lower than the last frequency in the lines resonance zone it can safely
!      be removed.
!
! 2. Line is removed when the current frequency is lower by EXT_LINE_VAR 
!      (converted to frequency units) than the last frequency in the resonance
!       zone. This is a control parameter, and may be used to speed up the code.
!       NB: EXT_LINE_VAR >= 0.
!
! 3. To make way for another line. This is only done when necessary, and is
!      done elsewhere. Only requirement is that the current frequency
!      is lower that the last frequency of the resonance zone.
!
	CALL TUNE(1,'CHK_L_FIN')
	T1=1.0D0-EXT_LINE_VAR*V(1)/2.998E+05         
	DO SIM_INDX=1,MAX_SIM
	  IF(LINE_STORAGE_USED(SIM_INDX))THEN
!
! Check whether need storage location for net rate etc. We keep the storage
! until the line levels are removed the line variation setcion.
!
	    L=SIM_LINE_POINTER(SIM_INDX)
	    IF(NU(ML) .LT. NU(LINE_END_INDX_IN_NU(L))*T1)THEN
	      SIM_LINE_POINTER(SIM_INDX)=0
              LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	      LINE_LOC(L)=0
	    END IF
!
! Zero variation storage, if in use, and the storage location is not also
! being used by some other line. 
!
	    IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION .AND.
	1        NU(ML) .LT. NU(LINE_END_INDX_IN_NU(L))*T1 .AND.
	1                                 .NOT. WEAK_LINE(SIM_INDX))THEN
	      NL=LOW_POINTER(SIM_INDX)
	      VAR_IN_USE_CNT(NL)=VAR_IN_USE_CNT(NL)-1
	      IF(VAR_IN_USE_CNT(NL) .EQ. 0)THEN
	        TX(:,:,NL)=0.0D0	!ND,ND,NM
	        TVX(:,:,NL)=0.0D0	!ND-1,ND,NM
	        dZ(NL,:,:,:)=0.0D0	!NM,NUM_BNDS,ND,MAX_SIM
	        VAR_LEV_ID(NL)=0
	      END IF
	      LOW_POINTER(SIM_INDX)=0
!
	      NUP=UP_POINTER(SIM_INDX)
	      VAR_IN_USE_CNT(NUP)=VAR_IN_USE_CNT(NUP)-1
	      IF(VAR_IN_USE_CNT(NUP) .EQ. 0)THEN
	        TX(:,:,NUP)=0.0D0	!ND,ND,NM
	        TVX(:,:,NUP)=0.0D0	!ND-1,ND,NM
	        dZ(NUP,:,:,:)=0.0D0	!NM,NUM_BNDS,ND,MAX_SIM
	        VAR_LEV_ID(NUP)=0
	      END IF
	      UP_POINTER(SIM_INDX)=0
	    END IF			!Outside region of influence by line?
	  END IF			!Line is in use.
	END DO				!Loop over line
	CALL TUNE(2,'CHK_L_FIN')
!
! 
!
	CALL TWO_PHOT_RATE(T,RJ,FL,FQW(ML),ND,NT)
!
! 
!
! Compute flux distribution and luminosity (in L(sun)) of star. NB: For
! NORDFLUX we always assume coherent scattering.
!
!
	CALL TUNE(IONE,'FLUX_DIST')
	IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
!
! Since ETAEXT is not required any more, it will be used
! flux.
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCEEXT,ZETAEXT,THETAEXT,RJEXT,NDEXT)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,
	1               SOURCEEXT,CHIEXT,dCHIdR,HQWEXT,ETAEXT,
	1               S1,THK_CONT,DIF,DBB,IC,NCEXT,NDEXT,NPEXT,METHOD)
	  CALL UNGRID(SOB,ND,ETAEXT,NDEXT,POS_IN_NEW_GRID)
	  SOB(2)=ETAEXT(2)				!Special case
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  N_OBS=NCF
	  OBS_FREQ(ML)=FL
	  OBS_FLUX(ML)=6.599341D0*SOB(1)*2.0D0		!2 DUE TO 0.5U
	ELSE IF(CONT_VEL)THEN
	   H_IN=DBB/CHI(ND)/3.0
	   IF(.NOT. DIF)H_IN=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
	   H_OUT=HBC_CMF(1)*RJ(1)
!
! TA is a work vector. TB initially used for extended SOB.
!
	   IF(ACCURATE)THEN
	     CALL REGRID_H(TB,REXT,RSQHNU,H_OUT,H_IN,NDEXT,TA)
	     DO I=1,ND
	       SOB(I)=TB(POS_IN_NEW_GRID(I))
	     END DO
	   ELSE
	     CALL REGRID_H(SOB,R,RSQHNU,H_OUT,H_IN,ND,TA)
	   END IF
	   IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
	     H_OUT=H_OUT-HBC_CMF(2)
	     SOB(1)=H_OUT; SOB(ND)=H_IN
	     SOB(1:ND)=SOB(1:ND)*R(ND)*R(ND)
	   END IF
	   CALL COMP_OBS_V2(IPLUS,FL,
	1           IPLUS_STORE,NU_STORE,NST_CMF,
	1           MU_AT_RMAX,HQW_AT_RMAX,OBS_FREQ,OBS_FLUX,N_OBS,
	1           V_AT_RMAX,RMAX_OBS,'IPLUS','LIN_INT',DO_FULL_REL_OBS,
	1           FIRST_OBS_COMP,NP_OBS)
!
	ELSE                          
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,THETA,HQW,SOB,
	1               S1,THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  N_OBS=NCF
	  OBS_FREQ(ML)=FL
	  OBS_FLUX(ML)=6.599341D0*SOB(1)*2.0D0		!2 DUE TO 0.5U
	END IF
!
! Evaluate (CMF) observd X-ray luminosities.
!
	  IF(ML .EQ. 1)THEN
	    OBS_XRAY_LUM_0P1=0.0D0
	    OBS_XRAY_LUM_1keV=0.0D0
	  END IF
	  T3=4.1274D-12
	  IF(NU(ML) .GT. 241.7988D0)OBS_XRAY_LUM_1keV=OBS_XRAY_LUM_1keV+T3*FQW(ML)*SOB(1)
	  IF(NU(ML) .GT. 24.17988D0)OBS_XRAY_LUM_0P1=OBS_XRAY_LUM_0P1+T3*FQW(ML)*SOB(1)
!
! Compute the luminosity, the FLUX mean opacity, and the ROSSELAND
! mean opacities.
!
	IF(ML .EQ. 1)THEN		!Need to move to main loop imit.
	  DO I=1,ND
	    RLUMST(I)=0.0D0
	    J_INT(I)=0.0D0
	    DJDt_TERM(I)=0.0D0
	    DJDt_FLUX(I)=0.0D0
	    K_INT(I)=0.0D0
	    FLUX_MEAN(I)=0.0D0
	    ROSS_MEAN(I)=0.0D0
	    INT_dBdT(I)=0.0d0
	  END DO
	END IF
	T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
	DO I=1,ND		              !(4*PI)**2*Dex(+20)/L(sun)
	  T2=SOB(I)*FQW(ML)*4.1274D-12
	  RLUMST(I)=RLUMST(I)+T2
	  J_INT(I)=J_INT(I)+RJ(I)*FQW(ML)*4.1274D-12
	  DJDt_FLUX(I)=DJDt_FLUX(I)+DJDt_TERM(I)*FQW(ML)*4.1274D-12
	  K_INT(I)=K_INT(I)+K_MOM(I)*FQW(ML)*4.1274D-12
	  FLUX_MEAN(I)=FLUX_MEAN(I)+T2*CHI(I)
	  T2=T1*EMHNUKT(I)/(  ( (1.0D0-EMHNUKT(I))*T(I) )**2  )
	  INT_dBdT(I)=INT_dBdT(I)+T2
	  ROSS_MEAN(I)=ROSS_MEAN(I)+T2/CHI(I)
	END DO
	T1=SPEED_OF_LIGHT()*1.0D-05
	DO J=1,N_FLUX_MEAN_BANDS
	  IF(0.01D0*T1/FL .LT. LAM_FLUX_MEAN_BAND_END(J))THEN
	     BAND_FLUX_MEAN(1:ND,J)=FLUX_MEAN(1:ND)
	     BAND_FLUX(1:ND,J)=RLUMST(1:ND)
	     EXIT
	  END IF
	END DO
	CALL TUNE(ITWO,'FLUX_DIST')
!
! The current opacities and emissivities are stored for the variation of the
! radiation field at the next frequency.
!                                                  
	DO I=1,ND
	  CHI_PREV(I)=CHI_CONT(I)
	  CHI_NOSCAT_PREV(I)=CHI_NOSCAT(I)
	  CHI_SCAT_PREV(I)=CHI_SCAT(I)
	  ETA_PREV(I)=ETA_CONT(I)
	END DO
!
	IF(LST_ITERATION)THEN
	  DO SIM_INDX=1,MAX_SIM
	    IF(END_RES_ZONE(SIM_INDX))THEN
	      LS=SIM_LINE_POINTER(SIM_INDX)
	      WRITE(LU_NET,'(/,1X,I6,2X,A,2X,F10.6,2X,I6,2X,I6)')
	1         LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1            VEC_NL(LS),VEC_NUP(LS)
	      WRITE(LU_DR,'(/,1X,I6,2X,A,2X,F10.6,2X,I6,2X,I6)')
	1         LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1            VEC_NL(LS),VEC_NUP(LS)
	      T3=(AVE_ENERGY(SIM_NL(SIM_INDX))-
	1            AVE_ENERGY(SIM_NUP(SIM_INDX)))/VEC_FREQ(LS)
	      WRITE(LU_HT,'(/,1X,I6,2X,A,2X,F10.6,2X,I6,2X,I6,ES14.4)')
	1         LS,TRANS_NAME_SIM(SIM_INDX),VEC_FREQ(LS),
	1            VEC_NL(LS),VEC_NUP(LS),T3
	      WRITE(LU_NET,'(1P,5E14.6)')(ZNET_SIM(I,SIM_INDX),I=1,ND)
	      WRITE(LU_DR,40003)((ZNET_SIM(I,SIM_INDX)*
	1                        POPS(SIM_NUP(SIM_INDX),I)*U_STAR_RATIO(I,SIM_INDX)*
	1                        EINA(SIM_INDX)),I=1,ND)
	      IF(SCL_LINE_COOL_RATES)THEN
	        T3=(AVE_ENERGY(SIM_NL(SIM_INDX))-
	1            AVE_ENERGY(SIM_NUP(SIM_INDX)))/VEC_FREQ(LS)
	        IF(ABS(T3-1.0D0) .GT. SCL_LINE_HT_FAC)T3=1.0D0
	      ELSE
	        T3=1.0D0
	      END IF
	      WRITE(LU_HT,'(1X,1P,5E12.4)')(
	1         T3*ZNET_SIM(I,SIM_INDX)*ETAL_MAT(I,SIM_INDX), I=1,ND)
	      DO K=1,ND
	        T2=ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	        STEQ_T_SCL(K)=STEQ_T_SCL(K) - T2*T3
	        STEQ_T_NO_SCL(K)=STEQ_T_NO_SCL(K) - T2
	      END DO
	      WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(STEQ_T_SCL(I), I=1,ND)
	      WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(STEQ_T_NO_SCL(I), I=1,ND)
	    END IF
	  END DO
	END IF
!
	WRITE(199,*)ML,STEQ_T(1)
10000	CONTINUE
	CALL TUNE(ITWO,'10000')
!
	WRITE(LUER,*)' '
	WRITE(LUER,*)'Number of weak lines is',NUM_OF_WEAK_LINES
	WRITE(LUER,*)' '
!
! 
!
	CALL STEQ_BA_TWO_PHOT_RATE_V3(POPS,NT,ND,
	1         DIAG_INDX,COMPUTE_BA,LUMOD,LST_ITERATION)
	WRITE(199,*)ML,STEQ_T(1),'two'
!
! Store radiative equlibrium equation (radiative terms only) so
! we can check influence on radiation field.
!
	DEP_RAD_EQ(1:ND)=STEQ_T(1:ND)
!
! 
!
! Include influence of charge exchange reactions.
!
	DO ID=1,NUM_IONS-1
	  ID_SAV=ID
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL SET_CHG_EXCH_V4(ID_SAV, ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,  ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV, ATM(ID)%GIONXzV_F, 
	1       ATM(ID)%NXzV_F, ATM(ID)%NXzV, ND, 
	1       ATM(ID)%EQXzV, EQ_SPECIES(SPECIES_LNK(ID)), T)
	  END IF
	END DO
!
	CALL STEQ_BA_CHG_EXCH_V3(POPS,T,NT,ND,DIAG_INDX,COMPUTE_BA)
	WRITE(199,*)ML,STEQ_T(1),'chg'
!
	DO ID=1,NUM_IONS-1
	  CALL EVAL_CHG_RATES_V3(ATM(ID)%CHG_PRXzV, ATM(ID)%CHG_RRXzV,
	1           ION_ID(ID),POPS,T,ND,NT)
	END DO
! 
!
! Output errors that have occurred in MOM_J_CMF
!
	CALL WRITE_J_CMF_ERR(MAIN_COUNTER)
!
! Allow for advection terms.
!
	IF(SN_MODEL .AND. DO_CO_MOV_DDT)THEN
          CALL STEQ_CO_MOV_DERIV_V2(ADVEC_RELAX_PARAM,LINEAR_ADV,
	1             DO_CO_MOV_DDT,LAMBDA_ITERATION,COMPUTE_BA,
	1             TIME_SEQ_NO,NUM_BNDS,ND,NT)
	  WRITE(199,*)ML,STEQ_T(1),'SN_DDT'
	ELSE
	  CALL STEQ_ADVEC_V4(ADVEC_RELAX_PARAM,LINEAR_ADV,NUM_BNDS,ND,
	1            INCL_ADVECTION,LAMBDA_ITERATION,COMPUTE_BA)
	END IF
!
! Allow for adiabatic cooling, if requested.
!
	IF(SN_MODEL .AND. DO_CO_MOV_DDT)THEN
	  CALL EVAL_TEMP_DDT_V2(dE_WORK,AD_COOL_V,AD_COOL_DT,
	1                       POPS,AVE_ENERGY,HDKT,
	1                       COMPUTE_BA,INCL_ADIABATIC,
	1                       TIME_SEQ_NO,DIAG_INDX,NUM_BNDS,NT,ND)
	  WRITE(199,*)ML,STEQ_T(1),'SN_ADD'
	ELSE
	  dE_WORK=0.0D0
	  CALL EVAL_ADIABATIC_V3(AD_COOL_V,AD_COOL_DT,
	1                       POPS,AVE_ENERGY,HDKT,
	1                       COMPUTE_BA,INCL_ADIABATIC,
	1                       DIAG_INDX,NUM_BNDS,NT,ND)
	END IF
!
	IF(SN_MODEL .AND. INCL_RADIOACTIVE_DECAY)THEN
	  CALL EVAL_RAD_DECAY_V1(dE_RAD_DECAY,NT,ND)
	  WRITE(199,*)ML,STEQ_T(1),'SN_RAD'
	END IF
!
! Prevent T from becoming too small by adding a extra heating term.
!
	CALL PREVENT_LOW_T(ARTIFICIAL_HEAT_TERM,T_MIN,COMPUTE_BA,LAMBDA_ITERATION,
	1                       T_MIN_BA_EXTRAP,DIAG_INDX,NUM_BNDS,ND,NT)
!
! Write pointer file and store BA, BA_ED and B_T matrices.
!
!	BA_T(:,:,ND)=0.0D0; STEQ_T(ND)=T(ND)-T(ND-1)
!	BA_T(NT,DIAG_INDX,ND)=1.0D0
!	IF(DIAG_INDX .NE. 1)BA_T(NT,DIAG_INDX-1,ND)=-1.0D0
	IF(COMPUTE_BA .AND. WRBAMAT .AND. .NOT. FLUX_CAL_ONLY .AND. .NOT. LAMBDA_ITERATION)THEN
	  CALL STORE_BA_DATA_V2(LU_BA,NION,NUM_BNDS,ND,COMPUTE_BA,'BAMAT')
	END IF
!               
!
! Write out recombination, photoionization and cooling terms for digestion.
!
! Since X_RECOM (and X_COOL) were dimensioned (ND,0:NION), and since it was 
! initialized, we can assume X_RECOM(1,0) is zero, and hence it can be used 
! as a  dummy vector for the special case when ATM(ID-1)%INDX=0.
!
! NB: For the first ionization stage of a species, ATM(ID-1)%INDX will
!     refer to the highest level of the previous ion.  By convention in
!     CMFGEN this refers to a 1 level state with an INDEX of 0 (and
!     XzV_PRES for this species is FALSE).
!
	IF(LST_ITERATION)THEN
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      J=0
	      IF(ID .NE. 1)J=ATM(ID-1)%INDX_XzV
	      TMP_STRING=TRIM(ION_ID(ID))//'PRRR'
	      CALL WRRECOMCHK_V3(ATM(ID)%APRXzV, ATM(ID)%ARRXzV,
	1          ATM(ID)%CPRXzV, ATM(ID)%CRRXzV,
	1          ATM(ID)%CHG_PRXzV, ATM(ID)%CHG_RRXzV, SE(ID)%STEQ_ADV,
	1          DIERECOM(1,ATM(ID)%INDX_XzV),ADDRECOM(1,ATM(ID)%INDX_XzV),
	1          X_RECOM(1,J),X_RECOM(1,ATM(ID)%INDX_XzV),
	1          R,T,ED,ATM(ID)%DXzV,TA,TB, ATM(ID)%NXzV,
	1          ND,LU_REC_CHK,TMP_STRING,ION_ID(ID))
	    END IF
	  END DO
!
! 
!
! We use TA for NETCR (the net cooling rate).
! We use TB for TOTCR (the sum[absolute cooling rates]).
!
	  DO ML=1,(ND+9)/10
	    LS=ML                
	    CALL FSTCOOL(R,T,ED,TA,TB,ML,ND,LU_REC_CHK)
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        CALL WRCOOLGEN(ATM(ID)%BFCRXzV, ATM(ID)%FFXzV, ATM(ID)%COOLXzV,
	1           DIECOOL(1,ATM(ID)%INDX_XzV), X_COOL(1,ATM(ID)%INDX_XzV),
	1           ATM(ID)%XzV_PRES, ATM(ID)%NXzV, ION_ID(ID),
	1           TA,TB,LS,ND,LU_REC_CHK)
	      END IF
	    END DO
!
! Output adiabatic cooling rate.
!
	    CALL WR_AD_COOL(AD_COOL_V,AD_COOL_DT,TA,TB,
	1            INCL_ADIABATIC,LS,ND,LU_REC_CHK)
!
! Output charge exchange cooling rate, radioactive heating team
! (hence option L_FALSE), and artificial heating term.
!
	    CALL WR_CHG_COOL_V3(TA,TB,LS,ND,LU_REC_CHK)
	    CALL WR_COOLING_TERM(dE_RAD_DECAY,TA,TB,L_FALSE,
	1             'Radiative decay heating term',LS,ND,LU_REC_CHK)
	    CALL WR_ART_HEAT(ARTIFICIAL_HEAT_TERM,TA,TB,LS,ND,LU_REC_CHK)
!
	    CALL ENDCOOL(TA,TB,LS,ND,LU_REC_CHK)
	  END DO
	END IF		!Only output if last iteration.
! 
!
	WRITE(LUER,*)'Luminosity of star is :',RLUMST(1),RLUMST(ND)
	IF(RLUMST(1) .LE. 0.0D0)RLUMST(1)=1.0E-20
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	  CALL WRITV(OBS_FREQ,N_OBS,'Continuum Frequencies',LU_FLUX)
	  CALL WRITV(OBS_FLUX,N_OBS,'Observed intensity (Janskys)',LU_FLUX)
	  CALL WRITV(RLUMST,ND,'Luminosity',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
	CALL TUNE(ITWO,'MLCF')
!
! Compute ROSSELAND and FLUX mean opacities. These MEAN opacities DO NOT 
! include the effect of clumping. Compute the respective optical depth scales; 
! TA is used for the FLUX mean optical depth scale, TB for the ROSSELAND mean 
! optical depth scale, and DTAU for the electron scattering optical depth 
! scale. The optical depth scale INCLDUES the effects of clumping. TCHI is used
! as a temporary work vector.
!
! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0E+16 / pi
!
	T1=7.218771D+11
	DO I=1,ND
	  FLUX_MEAN(I)=FLUX_MEAN(I)/RLUMST(I)
	  INT_dBdT(I)=INT_dBdT(I)/ROSS_MEAN(I)		!Program rosseland opac.
	  ROSS_MEAN(I)=T1*( T(I)**3 )/ROSS_MEAN(I)
	END DO
	DO J=1,N_FLUX_MEAN_BANDS
	  BAND_FLUX_MEAN(:,J)=BAND_FLUX_MEAN(:,J)/RLUMST(:)
	END DO
	TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
        CALL NORDTAU(TA,TCHI,R,R,dCHIdR,ND)
	TCHI(1:ND)=FLUX_MEAN(1:ND)*CLUMP_FAC(1:ND)
	IF(MINVAL(TCHI) .GT. 0)THEN
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
	ELSE
	  dCHIdR(1:ND)=0.0D0
	END IF
        CALL NORDTAU(TB,TCHI,R,R,dCHIdR,ND)
	TCHI(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
        CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
!
	TA(ND)=0.0D0
	TB(ND)=0.0D0
	TC(ND)=0.0D0
	DTAU(ND)=0.0D0
!
	CALL GEN_ASCI_OPEN(LU_OPAC,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(LU_OPAC,
	1  '( ''     R        I   Tau(Ross)   /\Tau   Rat(Ross)'//
	1  '  Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1  '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)     Kappa'' )' )
	  IF(R(1) .GE. 1.0D+05)THEN
	    FMT='( 1X,1P,E10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1        '4(2X,E9.3),2(2X,E8.2),3(2X,E8.2) )'
	  ELSE
	    FMT='( 1X,F10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1        '4(2X,E9.3),2(2X,E8.2),3(2X,E8.2) )'
	  END IF
	  DO I=1,ND
	    IF(I .EQ. 1)THEN
	      T1=LOG(ROSS_MEAN(1)*CLUMP_FAC(1)/ROSS_MEAN(4)/CLUMP_FAC(4))/LOG(R(4)/R(1))
	      IF(T1 .LT. 2.0)T1=2.0D0
	      T1=ROSS_MEAN(1)*CLUMP_FAC(1)*R(1)/(T1-1.0D0)		!Rosseland optical depth scale
	      T2=LOG(ABS(FLUX_MEAN(1)*CLUMP_FAC(1)/FLUX_MEAN(4)/CLUMP_FAC(4)))/LOG(R(4)/R(1))
	      IF(T2 .LT. 2.0)T2=2.0D0
	      T2=FLUX_MEAN(1)*CLUMP_FAC(1)*R(1)/(T2-1.0D0)		!Flux optical depth scale
	      T3=LOG(ESEC(1)*CLUMP_FAC(1)/ESEC(4)/CLUMP_FAC(4))/LOG(R(4)/R(1))
	      IF(T3 .LT. 2.0)T3=2.0D0
	      T3=ESEC(1)*CLUMP_FAC(1)*R(1)/(T3-1.0D0)			!Electon scattering optical depth scale
	      TC(1:3)=0.0D0
	    ELSE
	      T1=T1+TA(I-1)
	      T2=T2+TB(I-1)
	      T3=T3+DTAU(I-1)
	      TC(1)=TA(I)/TA(I-1)
	      TC(2)=TB(I)/TB(I-1)
	      TC(3)=DTAU(I)/DTAU(I-1)
	    END IF
	    WRITE(LU_OPAC,FMT)R(I),I,T1,TA(I),TC(1),
	1      ROSS_MEAN(I),INT_dBdT(I),FLUX_MEAN(I),ESEC(I),
	1      T2,T3,TC(2),TC(3),1.0D-10*ROSS_MEAN(I)/DENSITY(I)
	  END DO
	  WRITE(LU_OPAC,'(//,A,A)')
	1     'NB: Mean opacities do not include effect of clumping',
	1     'NB: Optical depth scale includes effect of clumping'
	CLOSE(UNIT=LU_OPAC)
!
	IF(LST_ITERATION)THEN
!
! Compute the grey temperature distribution and the Rosseland optical 
! depth scale (returned in TA).
!
	  CHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL COMP_GREY_V2(POPS,TGREY,TA,CHI,LUER,NC,ND,NP,NT)
!
	  OPEN(UNIT=LUIN,FILE='GREY_SCL_FACOUT',STATUS='UNKNOWN')
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,'(A,8X,A,7X,A,7X,A,6X,A)')'!','Log(Tau)','T/T(grey)','T(10^4 K)','L'
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,*)ND
	    DO I=1,ND
	      IF(TA(I) .GT. 0)THEN
	        WRITE(LUIN,'(2X,3ES16.6,4X,I3)')LOG10(TA(I)),T(I)/TGREY(I),T(I),I
	      ELSE
	        WRITE(LUER,'(A)')' Bad Roseeland optical depth scale for T/TGREY output'
	        WRITE(LUIN,'(A)')' Bad Roseeland optical depth scale for T/TGREY output'
	        EXIT
	      END IF
            END DO
	  CLOSE(LUIN)
	END IF
!
! Output hydrodynamical terms to allow check on radiation driving of the wind.
!
	IF(.NOT. SN_MODEL .AND. .NOT. USE_FIXED_J)THEN
	  I=18
	  CALL HYDRO_TERMS_V3(POP_ATOM,R,V,T,SIGMA,ED,RLUMST,
	1                 STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1		  FLUX_MEAN,ESEC,
	1                 LAM_FLUX_MEAN_BAND_END,BAND_FLUX_MEAN,
	1		  PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
	1                 BAND_FLUX,N_FLUX_MEAN_BANDS,I,ND)
	END IF
!
	IF(LST_ITERATION)
	1     CALL WR_ASCI_STEQ(NION,ND,'STEQ ARRAY- Continuum Terms',19)
!                                   
		CALL TUNE(ITHREE,' ')
! 
	CALL TUNE(IONE,'LINE LOOP')
!
! Zero LLUMST here, since total line luminosity output to OBSFLUX file
! even if the section of code is not executed.
!
	LLUMST(1:ND)=0.0D0
	IF(GLOBAL_LINE_SWITCH(1:5) .EQ. 'BLANK')THEN
!
! These lines have been treated with the continuum.
!
	ELSE IF(FLUX_CAL_ONLY .AND. .NOT. DO_SOBOLEV_LINES)THEN
!
! Don't waste time by computing rates and EW's for lines treated with the
! Sobolev approximation.
!
	ELSE
!
! This section of the program solves the transfer equation for
! the case of lines (SOB or CMF options).
!
	IF(OVERLAP .AND. GLOBAL_LINE_SWITCH(1:3)  .NE. 'SOB')THEN
	  WRITE(LUER,*)'WARNING in CMFGEN_SUB'
	  WRITE(LUER,*)' Overlapping lines not allowed in single line'//
	1              ' CMF calculation.'
	  WRITE(LUER,*)'Use Sobolev approx for overlapping lines.'
	  OVERLAP=.FALSE.
	END IF
!
! Set access counter for Continuum Eddington file, and for EW Jint.
!
	EDDINGTON=EDD_LINECONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    WRITE(LU_EDD,REC=4)ACCESS_F
	  ELSE
	    READ(LU_EDD,REC=4)ACCESS_F
	  END IF
	END IF
	ACCESS_JEW=1
!
! Enter line loop.
!
	IF(NLBEGIN .LE. 0)NLBEGIN=1
	LINE_INDX=NLBEGIN
	DO WHILE (LINE_INDX .LE. N_LINE_FREQ)
!
	  IF(OVERLAP)THEN
	    TMP_MAX_SIM=MAX_SIM
	  ELSE
	    TMP_MAX_SIM=1
	    OVER_FREQ_DIF=0.0D0
	  END IF
!
! Determine next line (or lines), and store line parameters (eg frequency,
! Einstein A coefficient etc) in appropriate locations for the subsequent 
! treatment of line.
!
! The temporary variable J [=MAX(N_LINE_FREQ,LINE_INDX+SIM_INDX) ] is used so that
! VEC_FREQ(N_LINE_FREQ+1) is not accessed (remember that the arguments of
! a "DO WHILE" can be done in any order.
!
	  SIM_INDX=0
	  J=LINE_INDX
	  SOBOLEV=.FALSE.
	  DO WHILE(LINE_INDX+SIM_INDX .LE. N_LINE_FREQ .AND.
	1      (VEC_FREQ(LINE_INDX)-VEC_FREQ(J))/VEC_FREQ(LINE_INDX) 
	1          .LE. OVER_FREQ_DIF 
	1          .AND. SIM_INDX .LT. TMP_MAX_SIM .AND.
	1      (VEC_TRANS_TYPE(J)(1:3) .EQ. 'CMF' .OR.
	1         VEC_TRANS_TYPE(J)(1:3) .EQ. 'SOB')    )
	    SIM_INDX=SIM_INDX+1
	    SIM_NL(SIM_INDX)=VEC_NL(LINE_INDX+SIM_INDX-1)
	    SIM_NUP(SIM_INDX)=VEC_NUP(LINE_INDX+SIM_INDX-1)
	    NL=SIM_NL(SIM_INDX)
	    NUP=SIM_NUP(SIM_INDX)
	    EINA(SIM_INDX)=VEC_EINA(LINE_INDX)
	    OSCIL(SIM_INDX)=VEC_OSCIL(LINE_INDX)
	    FL_SIM(SIM_INDX)=VEC_FREQ(LINE_INDX)
	    SIM_LINE_POINTER(SIM_INDX)=LINE_INDX
	    IF(VEC_TRANS_TYPE(LINE_INDX)(1:3) .EQ. 'SOB')SOBOLEV=.TRUE.
	    TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LINE_INDX))
	    AMASS_SIM(SIM_INDX)=AMASS_DOP
	    J=MIN(LINE_INDX+SIM_INDX,N_LINE_FREQ)
	  END DO
!
! Check to see that we do have a line. If not we go straight to end of
! this section. This may occur if we have treated some lines in the
! blanketing section.
!
	  NUM_SIM_LINES=SIM_INDX
	  IF(NUM_SIM_LINES .EQ. 0)THEN
	    LINE_INDX=LINE_INDX+1
	    GOTO 50000
	  END IF
!
! Determine center of blend. Used as continuum wavelength.
!
	  FL=0.0D0
	  DO SIM_INDX=1,NUM_SIM_LINES
	    FL=FL+FL_SIM(SIM_INDX)
	  END DO
	  FL=FL/NUM_SIM_LINES
!
	  CONT_FREQ=FL
	  COMPUTE_NEW_CROSS=.TRUE.
!
! VAR_SOB_JC indicates that the variation of Jc should be taken into account
!       when computing BA.
! NNM=4 indicates that the variation of Jc, ETAc and CHIc should be taken
!       into account when computing BA.
!
	 VAR_SOB_JC=.TRUE.
	 NNM=4
	 IF(NM_KI .LT. NNM)THEN
	   WRITE(LUER,*)'Error in CMFGEN -- NM_KI must be at least 4'
	   STOP
	 END IF
!
! Check to see if we have already completed this transition due to
! a restart of the program. We do the change here ( rather than
! "do nl=nlbegin,nt" ) so as the correct record is accessed for the
! continuum eddington factor.
!
	IF(LINE_INDX .LT. NLBEGIN)THEN
	  ACCESS_F=ACCESS_F + 1
	  ACCESS_JEW=ACCESS_JEW + 1
	  GO TO 50000
	END IF
!
! As a temporary measure, we set AMASS to AMASS_DOP. This ensures
! good profile coverage for all species.
!
	AMASS=AMASS_DOP
!
! Determine method to compute continuum intensity.
!
	IF(ACCURATE .AND. ALL_FREQ)THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE
	  THIS_FREQ_EXT=.FALSE.
	END IF
! 
!
! Compute LINE profile. This must be done here, so that we can check
! whether the TOTAL opacity (LINE+CONTINUUM) becomes negative.
! The profile is assumed to be depth independent, We normalize the
! frequency quadrature weights so that integral over the profile is one.
! 
	IF(.NOT. SOBOLEV)THEN
	  CALL TRAPUNEQ(PF,LFQW,NLF)
	  T1=0.0
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)*FL*1.0D+15
	    PROF(ML)=DOP_PRO(PF(ML),FL,TDOP,VTURB,AMASS)
	    T1=T1+LFQW(ML)*PROF(ML)
	  END DO
!
! Now normalize frequency weights.
!
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)/T1
	  END DO
	END IF
! 
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS. 
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  I=SIM_LINE_POINTER(SIM_INDX)
	  MNL_F=VEC_MNL_F(I)
	  MNUP_F=VEC_MNUP_F(I)
	  DO K=1,ND
	    L_STAR_RATIO(K,SIM_INDX)=0.0D0
	    U_STAR_RATIO(K,SIM_INDX)=0.0D0
	    dL_RAT_dT(K,SIM_INDX)=0.0D0
	    dU_RAT_dT(K,SIM_INDX)=0.0D0
	  END DO
	  DO ID=1,NUM_IONS-1
	    IF(VEC_SPEC(I) .EQ. ION_ID(ID))THEN
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	      DO K=1,ND
	        T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzVLTE_F(MNL_F,K)) /
	1                (ATM(ID)%XzV(MNL,K)/ATM(ID)%XzVLTE(MNL,K))
	        L_STAR_RATIO(K,SIM_INDX)=T1*ATM(ID)%W_XzV_F(MNUP_F,K)*
	1             ATM(ID)%XzVLTE_F(MNL_F,K)/ATM(ID)%XzVLTE(MNL,K)/
	1             ATM(ID)%W_XzV_F(MNL_F,K)
	        T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzVLTE_F(MNUP_F,K)) /
	1             (ATM(ID)%XzV(MNUP,K)/ATM(ID)%XzVLTE(MNUP,K))
	        U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F(MNUP_F,K)/
	1              ATM(ID)%XzVLTE(MNUP,K)
	        dL_RAT_dT(K,SIM_INDX)=L_STAR_RATIO(K,SIM_INDX)*
	1         (-1.5D0-HDKT*ATM(ID)%EDGEXzV_F(MNL_F)/T(K)
	1               -ATM(ID)%dlnXzVLTE_dlnT(MNL,K))/T(K)
	        dU_RAT_dT(K,SIM_INDX)=U_STAR_RATIO(K,SIM_INDX)*
	1         (-1.5D0-HDKT*ATM(ID)%EDGEXzV_F(MNUP_F)/T(K)
	1               -ATM(ID)%dlnXzVLTE_dlnT(MNUP,K))/T(K)
	      END DO
	      GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	      TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1           '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1                TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	      EXIT
	    END IF
	  END DO
	END DO
!
! Compute line opacity and emissivity. 
!                     
	DO SIM_INDX=1,NUM_SIM_LINES
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    NEG_OPACITY(I)=.FALSE.
	  END DO
	END DO
! 
!
	  SECTION='LINE'
	  IF(IMPURITY_CODE)THEN
!
! Obtain previously computed continuum opacities, and mean intensities.
!
	    INCLUDE 'GET_J_CHI.INC'
	  ELSE
!
! Compute continuum opacity and emissivity at the line frequency.
!
!	    INCLUDE 'OPACITIES_V4.INC'
	    CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)

!
! Compute continuum intensity.
!
!	    INCLUDE 'COMP_JCONT_V4.INC'	
            CALL COMP_J_BLANK(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1                              MAXCH,LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
!
	  END IF
!
! SOURCE is used by SOBJBAR and in VARCONT. Note that SOURCE is corrupted 
! (i.e. set to line source function) in CMFJBAR. Also note that SOURCEEXT 
! has previously been computed.
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+RJ(I)*THETA(I)
	  END DO
! 
!
! Scale emissivity because of different frequencies.
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  DO I=1,ND
	    T1=FL**3/(EXP(HDKT*FL/T(I))-1.0D0)
	    T2=FL_SIM(SIM_INDX)**3/(EXP(HDKT*FL_SIM(SIM_INDX)/T(I))-1.0D0)
	    BB_COR(I,SIM_INDX)=T1/T2
	  END DO
	END DO
!
! At present overlapping lines in the CMF frame is not installed. We still
! use CHIL and ETAL --- not CHIL_MAT etc. 
!
	IF(SOBOLEV .OR. NUM_SIM_LINES .EQ. 1)THEN
	  DO I=1,ND
	    CHIL(I)=0.0D0
	    ETAL(I)=0.0D0
	    DO SIM_INDX=1,NUM_SIM_LINES
	      CHIL(I)=CHIL(I)+CHIL_MAT(I,SIM_INDX)
	      ETAL(I)=ETAL(I)+ETAL_MAT(I,SIM_INDX)*BB_COR(I,SIM_INDX)
	    END DO
	  END DO
	END IF
!
        IF(CHECK_LINE_OPAC .AND. SOBOLEV)THEN
	  T1=3.0D-10/FL
	  FIRST_NEG=.TRUE.
	  DO I=1,ND
	    T2=CHIL(I)*T1*R(I)/V(I)		!Tau(Sobolev) for dlnV/dlnr=1
	    IF(T2 .LT. -0.5)THEN
	      NEG_OPACITY(I)=.TRUE.
	      IF(LST_ITERATION)THEN
	        IF(FIRST_NEG)THEN
	          J=ICHRLEN(TRANS_NAME_SIM(1))
	          WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(1)(1:J)
	          FIRST_NEG=.FALSE.
	        END IF
	        WRITE(LU_NEG,2290)I,T2,ED(I)
2290	        FORMAT(1X,'I= ',I3,'  : TAU(Sob)= ',1P,E10.2,'  : Ne=',E9.2)
	        WRITE(LU_NEG,*)CHIL(I),POPS(NL,I),R(I),V(I),T1
	      END IF
	      CHIL(I)=1.0D0	!Reset after we output its value.
	    END IF
	  END DO
!
! Note that both CHIL and CHIL_MAT are used by SOBJBAR_SIM to compute ZNET.
! The following ensures consistency. Alternative is to make positive
! only those opacities that are negative --- then need  NEG_OPACITY indicator
! for each line?
!
	  DO SIM_INDX=1,NUM_SIM_LINES
	    DO I=1,ND
              IF(NEG_OPACITY(I))CHIL_MAT(I,SIM_INDX)=1.0D0/NUM_SIM_LINES
	    END DO
	  END DO
!
	ELSE IF(CHECK_LINE_OPAC)THEN
	  FIRST_NEG=.TRUE.
	  DO I=1,ND
	    IF( (CHIL(I)*PROF(NLF/2+1)+CHI(I)) 
	1            .LT. 0.2D0*CHI(I) )THEN
	      T2=CHIL(I)*PROF(NLF/2+1)+CHI(I)
	      CHIL(I)=1.0D-04*ESEC(I)/PROF(NLF/2+1)
	      NEG_OPACITY(I)=.TRUE.
	      IF(LST_ITERATION)THEN
	        IF(FIRST_NEG)THEN
	          J=ICHRLEN(TRANS_NAME_SIM(1))
	          WRITE(LU_NEG,'(A,4I6)')TRANS_NAME_SIM(1)(1:J),
	1                                    MNL,NL,MNUP,NUP
	          FIRST_NEG=.FALSE.
	        END IF
	        WRITE(LU_NEG,2295)I,T2/ESEC(I),ESEC(I)/CHI(I),ED(I)
2295	        FORMAT(1X,'I= ',I3,' : CHIL/ESEC=',1P,E9.2,
	1           '  : ESEC/CHI=',E9.2,'  : Ne=',E9.2)
	      END IF
	    END IF
	  END DO
	END IF
!
! 
!
! Determine the net rates and the variation of the net rates with both
! upper and lower levels. NB: ATM(4)%EQXzV refers to EQHE2.
!
	IF(NL .EQ. ATM(4)%EQXzV .AND. SETZERO .AND. .NOT. OVERLAP)THEN
!
! This routine sets the rates in the He2 resonance lines to Zero.
! As the He2 continuum is thick, the n=2 level of He2 will be collisionally
! coupled to the n=1 level. 
!
	  DO I=1,ND
	    ZNET(I)=0.0D0
	    VB(I)=0.0D0
	    VC(I)=0.0D0
	  END DO
	ELSE IF(.NOT. SOBOLEV)THEN
!
	  CALL SUB_CMF_LINE(SECTION,POPS,CHIL,ETAL,NEG_OPACITY,
	1                    FL,CONT_FREQ,AMASS,
	1                    EDDINGTON,IMPURITY_CODE,
	1                    EW,CONT_INT,COMPUTE_EW,
	1                    COMPUTE_JEW,LU_JEW,ACCESS_JEW,
	1                    NL,NUP,NT,ND,NC,NP,
	1                    NDEXT,NCEXT,NPEXT,
	1                    NLF,NNM,
	1                    DIAG_INDX,NUM_BNDS)
!
	ELSE
!
! Use the escape probability approximation for lines originating
! in all levels.
!
	  CALL SUB_SOB_LINE_V2(SECTION,POPS,AVE_ENERGY,CHIL,ETAL,NEG_OPACITY,
	1                    FL,CONT_FREQ,AMASS,
	1                    EDDINGTON,IMPURITY_CODE,VAR_SOB_JC,LST_ITERATION,
	1                    EW,CONT_INT,
	1                    NL,NUP,NT,ND,NC,NP,
	1                    NDEXT,NCEXT,NPEXT,
	1                    NLF,DIAG_INDX,NUM_BNDS)
!
	END IF
!
! Outpute line EW, net rate, total rate, contribution of line to the luminosisty.
!
	IF(LST_ITERATION)THEN
	  T1=LAMVACAIR(FL_SIM(1)) 		!Wavelength(Angstroms)
	  DO SIM_INDX=1,NUM_SIM_LINES
	    NUP=SIM_NUP(SIM_INDX)
	    CALL EW_FORMAT(EW_STRING,TRANS_NAME_SIM(SIM_INDX),T1,
	1                     CONT_INT,EW,SOBOLEV)
	    IF(SIM_INDX .NE. 1)THEN
	      I=INDEX(EW_STRING,'  ')
	      EW_STRING(3:I+1)=EW_STRING(1:I-1)
	      EW_STRING(2:3)='##'
	    END IF
	    L=ICHRLEN(EW_STRING)
	    WRITE(LU_NET,40002)EW_STRING(1:L)
	    WRITE(LU_DR,40002)EW_STRING(1:L)
	    WRITE(LU_EW,40005)EW_STRING(1:L)
	    WRITE(LU_HT,40002)EW_STRING(1:L)
	    WRITE(LU_NET,40003)(ZNET_SIM(I,SIM_INDX),I=1,ND)
	    WRITE(LU_DR,40003)((ZNET_SIM(I,SIM_INDX)*POPS(NUP,I)*
	1                        EINA(SIM_INDX)),I=1,ND)
	    WRITE(LU_HT,40003)
	1        (ZNET_SIM(I,SIM_INDX)*ETAL_MAT(I,SIM_INDX),I=1,ND)
	  END DO
40002	  FORMAT(//,A)
40003	  FORMAT(3X,1P,5E16.5)
40005	  FORMAT(A)
	END IF
!
! Update line luminosity.
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  DO I=1,ND
	    LLUMST(I)=LLUMST(I)+ZNET_SIM(I,SIM_INDX)*ETAL_MAT(I,SIM_INDX)
	  END DO
	END DO
!
50000	CONTINUE
	LINE_INDX=LINE_INDX+NUM_SIM_LINES
	END DO 			!END NL LOOP
	CALL TUNE(ITWO,'LINE LOOP')
!
! Write pointer file and then store BA and STEQ matrices.
!
! These are the final matrices.
!
	IF(COMPUTE_BA
	1          .AND. WRBAMAT
	1          .AND. .NOT. FLUX_CAL_ONLY 
	1          .AND. .NOT. LAMBDA_ITERATION)THEN
	  CALL STORE_BA_DATA_V2(LU_BA,NION,NUM_BNDS,ND,COMPUTE_BA,'BAMAT')
	END IF
!
! Ends check on GLOBAL_LINE.
!
			END IF
!
!
! Compute the the total line luminosity, and the total dielectronic line 
! emission emitted in each shell between i and i+1.
! Not all of this flux will be received by the observer due to continuum
! absorption. The factor of 0.5 arrises from the average of R*R*ZNET in the
! shell. [ 4.1274D-12=(4pi*1.0D+20)/Lsun * 4PI ]
!
	IF(.NOT. XRAYS)THEN
	  XRAY_LUM_TOT(1:ND)=0.0D0
	  XRAY_LUM_0P1(1:ND)=0.0D0
	  XRAY_LUM_1KEV(1:ND)=0.0D0
	END IF
	T1=4.1274D-12
	T2=1.0D+10*T1/4.0D0/ACOS(-1.0D0)
	DO I=1,ND
	  LLUMST(I)=LLUMST(I)*R(I)*R(I)*T1
	  DIELUM(I)=DIELUM(I)*R(I)*R(I)*T1
	  DEP_RAD_EQ(I)=DEP_RAD_EQ(I)*R(I)*R(I)*T1
	  XRAY_LUM_TOT(I)=XRAY_LUM_TOT(I)*R(I)*R(I)*T1
	  XRAY_LUM_0P1(I)=XRAY_LUM_0P1(I)*R(I)*R(I)*T1
	  XRAY_LUM_1KEV(I)=XRAY_LUM_1KEV(I)*R(I)*R(I)*T1
	  dE_WORK(I)=dE_WORK(I)*R(I)*R(I)*T1
	  dE_RAD_DECAY(I)=dE_RAD_DECAY(I)*R(I)*R(I)*T2		!As ergs/cm^3
	END DO
!
	CALL LUM_FROM_ETA(LLUMST,R,ND)
	CALL LUM_FROM_ETA(DIELUM,R,ND)
	CALL LUM_FROM_ETA(DEP_RAD_EQ,R,ND)
	CALL LUM_FROM_ETA(dE_WORK,R,ND)
	CALL LUM_FROM_ETA(dE_RAD_DECAY,R,ND)
	CALL LUM_FROM_ETA(XRAY_LUM_TOT,R,ND)
	CALL LUM_FROM_ETA(XRAY_LUM_0P1,R,ND)
	CALL LUM_FROM_ETA(XRAY_LUM_1KeV,R,ND)
!
	IF(PLANE_PARALLEL_NO_V)THEN
	  MECH_LUM(1:ND)=0.0D0
	ELSE IF(PLANE_PARALLEL)THEN
	  T1=R(ND)*R(ND)*1.0D+05/SPEED_OF_LIGHT()	!As V in km/s, c in cgs units.
	  DO I=1,ND
	    MECH_LUM(I)=T1*V(I)*(1.0D0+SIGMA(I))*K_INT(I)/R(I)
	  END DO
	ELSE IF(USE_J_REL .AND. INCL_REL_TERMS)THEN
	  T1=1.0D+05/SPEED_OF_LIGHT()			!As V in km/s, c in cgs units.
	  DO I=1,ND
	    T2=1.0D0/(1.0D0-(T1*V(I))**2)		!Gamma^2
	    T3=SQRT(T2)					!Gamma
	    WRITE(205,'(3ES16.6)')T1,T2,T3
	    WRITE(205,'(I5,3ES14.4)')I,J_INT(I),K_INT(I),T1*V(I)*RLUMST(I)/R(I)/R(I)/T3
	    MECH_LUM(I)=T1*R(I)*V(I)*T3*( J_INT(I)-K_INT(I) +
	1          T2*(SIGMA(I)+1.0D0)*(K_INT(I)+T1*V(I)*RLUMST(I)/R(I)/R(I)) )
	    RLUMST(I)=T3*(RLUMST(I)+T1*V(I)*J_INT(I)*R(I)*R(I))
	  END DO
	ELSE
	  T1=1.0D+05/SPEED_OF_LIGHT()	!As V in km/s, c in cgs units.
	  DO I=1,ND
	    MECH_LUM(I)=T1*R(I)*V(I)*(J_INT(I)+SIGMA(I)*K_INT(I))
	  END DO
	END IF
	CALL LUM_FROM_ETA(MECH_LUM,R,ND)
	CALL LUM_FROM_ETA(DJDT_FLUX,R,ND)
!
! Increment the continuum luminosity by the total line luminosity.
!
        T1=0.0D0
	T2=0.0D0
	DO I=1,ND-1
	  DO J=I,ND-1
	    RLUMST(I)=RLUMST(I)+LLUMST(J)+DIELUM(J)
	  END DO
          T1=T1+LLUMST(I)
	  T2=T2+DIELUM(I)
	END DO
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','OLD','APPEND',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening OBSFLUX to output Rec emission'
	    WRITE(LUER,*)'IOS=',IOS
	  END IF
	  IF(USE_J_REL)CALL WRITV(RLUMST,ND,'Luminosity [g.r^2.H + beta.g.r^2.J]',LU_FLUX)
	  IF(SUM(DIELUM) .NE. 0.0D0)CALL WRITV(DIELUM,ND,
	1    'Dielectronic and Implicit Recombination Line Emission',LU_FLUX)
	  IF(SUM(LLUMST) .NE. 0.0D0)CALL WRITV(LLUMST,ND,'Line Emission',LU_FLUX)
	  CALL WRITV(MECH_LUM,ND,'Mechanical Luminosity',LU_FLUX)
	  IF(SUM(de_WORK) .NE. 0.0D0)CALL WRITV(dE_WORK,ND,'Internal/adiabatic term',LU_FLUX)
	  IF(SUM(DJDt_FLUX) .NE. 0.0D0)CALL WRITV(DJDt_FLUX,ND,'Flux arrising from Dr^3J/Dt term',LU_FLUX)
	  IF(SUM(dE_RAD_DECAY) .NE. 0.0D0)
	1     CALL WRITV(dE_RAD_DECAY,ND,'Energy deposited locally due to radioactive decay',LU_FLUX)
	  CALL WRITV(DEP_RAD_EQ,ND,'Departure from Rad Equilibrium Correction',LU_FLUX)
	  CALL WRITV(RLUMST,ND,'Total Radiative Luminosity',LU_FLUX)
	  IF(SUM(XRAY_LUM_TOT) .NE. 0.0D0)CALL WRITV(XRAY_LUM_TOT,ND,'Total Shock Luminosity (Lsun)',LU_FLUX)
!
! Include the machanical luminosity imparted to the wind by the radiation 
! field in the total luminosity, and subtract out radioactive energy deposition..
!
	  T3=0.0D0
	  DO I=ND-1,1,-1
	    T3=T3+MECH_LUM(I)+DJDT_FLUX(I)+dE_WORK(I)-dE_RAD_DECAY(J)
	    RLUMST(I)=RLUMST(I) + T3
	  END DO
	  CALL WRITV(RLUMST,ND,'Luminosity Check (not observed luminosity)',LU_FLUX)
!
	  T3=0.0D0
	  DO I=ND-1,1,-1
	    T3=T3+DEP_RAD_EQ(I)
	    RLUMST(I)=RLUMST(I) + T3
	  END DO
	  CALL WRITV(RLUMST,ND,'Consistency check (include dep. from rad. equil.)',LU_FLUX)
!
	  WRITE(LU_FLUX,'(A)')' '
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Line luminosity:',T1
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')
	1   'Total Dielectronic and Implicit Recombination Luminosity:',T2
	  WRITE(LU_FLUX,'(A,T60,1PE12.4)')'Total Mechanical Luminosity:',T3
!
! The seocnd XRAY flux printed is the OBSERVED XRAY luminosity. Its should be very similar
! to the eralier value for optically thin winds.
!
	  WRITE(LU_FLUX,'(A,T60,ES12.4)')'Total Shock Luminosity (Lsun):',SUM(XRAY_LUM_TOT)
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 0.1 keV, Lsun) :',
	1                                     SUM(XRAY_LUM_0P1),OBS_XRAY_LUM_0P1
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 1 keV, Lsun):',
	1                                     SUM(XRAY_LUM_1KEV),OBS_XRAY_LUM_1KEV
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emitted & observed X-ray Luminosity (> 0.1 keV, Lstar) :',
	1                                     SUM(XRAY_LUM_0P1)/LUM,OBS_XRAY_LUM_0P1/LUM
	  WRITE(LU_FLUX,'(A,T60,2ES12.4)')'Emiited & observed X-ray Luminosity (> 1 keV, Lstar):',
	1                                     SUM(XRAY_LUM_1KEV)/LUM,OBS_XRAY_LUM_1KEV/LUM
	CLOSE(UNIT=LU_FLUX)
!
! Quick and dirty way of ensuring people don't take notic of OBSFLUX when USE_FIXED_J is TRUE.
!
	IF(USE_FIXED_J)THEN
	  CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','OLD',' ',' ',IZERO,IOS)
	  WRITE(LU_FLUX,'(///,A)')'THIS FILE IS USELESS WHEN USE_FIXED_J is TRUE'
	  CLOSE(UNIT=LU_FLUX)
	END IF
!
! Because of the use of SL's, it is possible that a error in Luminosity
! can occur at depth, particulary in SN models where line escape is 
! enhanced because of the velocity field. TA provides an estimate of
! this error. This error is not relevant in the outer regions.
!
	IF(LST_ITERATION)THEN
	  DO I=1,ND
	    TA(I)=4.1274D-12*(STEQ_T_NO_SCL(I)-STEQ_T_SCL(I))*R(I)*R(I)
	  END DO
	  CALL LUM_FROM_ETA(TA,R,ND)
	  WRITE(LU_HT,'(//,A)')' Estimated error in L due to use of SLs'
	  WRITE(LU_HT,'(/,(1X,1P,5E12.4))')(TA(I), I=1,ND)
	END IF
!
! Insure Eddington factor file is closed, and indicate all f's successfully
! computed. If COMPUTE_EDDFAC is true we can safely write to record 5
! as file must have new format.
!
	T1=1.0
	IF(COMPUTE_EDDFAC)WRITE(LU_EDD,REC=5)T1
	CLOSE(UNIT=LU_EDD)
	IF(.NOT. COHERENT_ES)CLOSE(UNIT=LU_ES)
	COMPUTE_JEW=.FALSE.
	COMPUTE_EDDFAC=.FALSE.
!
! IF we are doing a FLUX computation only we do not corrupt the SCTRMEP file
! or the BA matrix. We also do not output th populations. This can be done
! quickly by setting FLUX_CAL_ONLY=.FALSE. and putting N_ITS=0.
!
	 IF(FLUX_CAL_ONLY .AND. RD_COHERENT_ES)THEN
	   WRITE(LUER,*)'Stopping CMFGEN as finished FLUX calculation.'
	   WRITE(LUER,*)'For a FLUX calculation we do 1 iteration only'
	   STOP
!
! Compute the convolution of J with the electron redistribution function.
! May need to compute a flux spectrum sveral times in order for e.s.
! redistributon to be correctly allowed for. 
! RD_NU and ALLOW_UNEQUAL_FREQ are both set to FALSE. Note that TEXT and NDEXT
! contain T and ND when ACCURATE is FALSE.
!
	 ELSE IF(FLUX_CAL_ONLY .AND. .NOT. RD_COHERENT_ES)THEN
	   COHERENT_ES=RD_COHERENT_ES
	   I=SIZE(VJ)
	   CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LU_EDD,'EDDFACTOR',
	1             EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
!
! Close units 2 and 16 to force writing of information.
!
	   CLOSE(UNIT=LUER)
	   CLOSE(UNIT=LU_SE)
	   CALL GEN_ASCI_OPEN(LUER,'OUTGEN','OLD','APPEND',' ',IZERO,IOS)
	   CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','OLD','APPEND',' ',IZERO,IOS)
!
! Now do another iteration. We continue to iterate so that the accuracy of
! RJ_ES is improved (i.e. to allow for multiple scattering). The number of
! iterations is set by NUM_ITS_TO_DO in the input file.
!
	   IF(.NOT. LST_ITERATION)GOTO 20000
	   STOP
	END IF
!
	CALL WRITV(STEQ_T,ND,'Radiative Equlibrium Equation',LU_SE)
	
	CALL SOLVE_FOR_POPS(POPS,NT,NION,ND,NC,NP,NUM_BNDS,DIAG_INDX,
	1      MAXCH,MAIN_COUNTER,IREC,LU_SE,LUSCR,LST_ITERATION)
!
! If we have changed the R grid, we need to recomput the angular quadrature weitghts,
! and put the atom density ect on the new radius grid.
!
	IF(REVISE_R_GRID)THEN
	  CALL SET_ANG_QW(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,TRAPFORJ,ACCURATE)
!
! Compute CLUM_FAC(1:ND) which allow for the possibility that the wind is
! clumped. At the sime time, we compute the vectors which give the density,
! the atom density, and the species density at each depth.
! The new call replaces the interpolation done in
!				  CALL ADJUST_DEN_VECS(R_OLD,ND)
!
	  CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	END IF

! 
!
! Store populations back into individual arrays. At present, program stops
! if an error occurs in the NG acceleration. Could be changed by doing
! the conversion below before the NG call. If the NG accelerate worked,
! we would do the conversion again. IF it failed, we would do the
! reverse conversion (as POPS might be corrupted).
!
	DO ID=1,NUM_IONS-1
	  CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV,ED,T,
	1    ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, ND, ATM(ID)%XzV_PRES)
	END DO
!
! We have now store the revised populations back in their individual
! storage locations. For some species we have 2 atomic models. For these
! species we need to take the super-level populations and compute:
!
! 1. The LTE population off all level ls in the FULL atom.
! 2. The population off all levels in the FULL atom.
! 3. The LTE population off all super-levels.
!
! This is done by the include file SUP_TO_FULL which calls the subroutine
! SUP_TO_FULL.FOR
!
! We first, however, need to compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=1,ND
	  POPION(J)=0.0D0
	  DO I=1,NT
	    IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Revise vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
	CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!	INCLUDE 'SUP_TO_FULL_V4.INC'
!
! Perform a revision to th hydrostatic structure of the atmosphere. If we have changed the 
! R grid, we need to recomput the angular quadrature weitghts, and put the atom density etc
! on the new radius grid.
!
	DONE_HYDRO_REVISION=.FALSE.
	IF(DO_HYDRO)THEN
	  K=MAIN_COUNTER
	  CALL DO_CMF_HYDRO_V2(POPS,LUM,TEFF,LOGG,STARS_MASS,RP,RMAX,RMDOT,VINF,V_BETA1,
	1          PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
	1          K,DONE_HYDRO_REVISION,NC,ND,NP,NT)
	  IF(DONE_HYDRO_REVISION)THEN
	    LAMBDA_ITERATION=.TRUE.
            FIXED_T=.TRUE.
	    FIX_IMPURITY=.FALSE.
	    COMPUTE_BA=.TRUE.
	    MAIN_COUNTER=MAIN_COUNTER+1
	    IF(ACCURATE)THEN
	      I=ND-DEEP
	      CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1              ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	      TA(1:ND)=1.0D0	!TEXT not required, T currently zero
	      CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,V,TA,SIGMA,ND)
              VDOP_VEC_EXT(1:NDEXT)=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
	    END IF
	    CALL SET_ANG_QW(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,TRAPFORJ,ACCURATE)
	    CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	    CALL GREY_T_ITERATE(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  END IF
	END IF
!
! Initialize pointer file for storage of BA matrix.
!
	I=-1000
	CALL INIT_BA_DATA_PNT_V2(LU_BA,NION,NUM_BNDS,ND,COMPUTE_BA,'BAMAT')
!
!
! 
!
! Compute the convolution of J with the electron redistribution function.
! Fist neeed to UPDATE TEXT because of the linearization. We need TEXT for
! convolving J with the electron scattering redistribution function.
! VEXT and SIGMAEXT have already been computed. 
!
	TEXT(1:ND)=T(1:ND)
	IF(ACCURATE)THEN
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,T,SIGMA,ND)
	END IF
	COHERENT_ES=RD_COHERENT_ES
	IF(.NOT. COHERENT_ES)THEN
	  I=SIZE(VJ)
	  CALL COMP_J_CONV_V2(VJ,I,NU,TEXT,NDEXT,NCF,LU_EDD,'EDDFACTOR',
	1             EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
	END IF          
!
! 
! Output brief summary of the model. This is to facilate the creation
! of compact model logs.
!
	IF(LST_ITERATION)THEN
!
	  CALL GEN_ASCI_OPEN(LUMOD,'MOD_SUM','UNKNOWN',' ',' ',IZERO,IOS)
!
	  WRITE(LUMOD,'(/,''Model Started on:'',15X,(A))')TIME
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(''Model Finalized on:'',13X,(A))')TIME
	  WRITE(LUMOD,
	1       '(''Main program last changed on:'',3X,(A))')PRODATE
	  WRITE(LUMOD,'()')
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'ND',ND)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NC',NC)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NP',NP)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NT',NT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NUM_BNDS',NUM_BNDS)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NCF',NCF)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NLINES',N_LINE_FREQ)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Output brief summary of atomic models.
!
	  DO ISPEC=1,NUM_SPECIES
	    STRING=' '
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IF(ATM(ID)%XzV_PRES)THEN
	       CALL WR_SL_INFO(STRING,ATM(ID)%NXzV,ATM(ID)%NXzV_F,
	1                        ATM(ID)%ZXzV,ION_ID(ID),LUMOD)
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	  END DO
	  WRITE(LUMOD,'(A)')' '
!
! Output stellar parameters.
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'L*',LUM)
	  T1=RMDOT/3.02286D+23
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Mdot',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*  ',RP)
	  T1=RMAX/RP   ; CALL WR_VAL_INFO(STRING,NEXT_LOC,'RMAX/R*',T1)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
! Compute the Radius and Velocity at Tau=10, and at Tau=2/3.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(1:ND)=0.0D0 ; DO I=2,ND ; TA(I) = TA(I-1)+DTAU(I-1) ; END DO
	  TB(1)=MIN(2.0D0/3.0D0,TA(ND))  ; TB(2)=MIN(10.0D0,TA(ND))
	  TB(3)=MIN(20.0D0,TA(ND))
	  CALL MON_INTERP(TC,ITHREE,IONE,TB,ITHREE,R,ND,TA,ND)
	  CALL MON_INTERP(AV,ITHREE,IONE,TB,ITHREE,V,ND,TA,ND)
!
! Output summary of Teff, R, and V at RSTAR, Tau=10, and TAU=2/3.
! For a plane-parallel atmosphere, R(ND) defines Teff.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TA(ND))
	  T1=1.0D+10*RP/RAD_SUN() ; CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*/Rsun',T1)
	  T1=TEFF_SUN()*(ABS(LUM)/T1**2)**0.25					!ABS for SN
	  IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)THEN
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff  ',T1)
	  ELSE
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'T*  ',T1)
	  END IF
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',V(ND))
	  IF(DO_HYDRO)THEN
	    T1=LOG10(1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/RP/RP)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	  END IF
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  IF(TA(ND) .GT. 20.0D0 .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(3))		!20.0D0
	    T1=1.0D+10*TC(3)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM)/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(3))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(3)/TC(3))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 10.0D0 .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(2))		!10.0D0
	    T1=1.0D+10*TC(2)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM)/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(2))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(2)/TC(2))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 0.67D0 .AND. .NOT. (PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL) )THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(1))		!0.67D0
	    T1=1.0D+10*TC(1)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(ABS(LUM)/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(1))
	    IF(DO_HYDRO)THEN
	      T1=LOG10(1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()/TC(1)/TC(1))
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,'Log g',T1)
	    END IF
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  STRING=' '
	  NEXT_LOC=1
	  T1=4.9376D+07*(RMDOT/3.02286D+23)*V(1)/LUM
	  T2=8.235D+03*(RMDOT/3.02286D+23)*V(1)*V(1)/LUM
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Eta',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Ek/L(%)',T2)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Velocity law information.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf1',VINF1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta1',V_BETA1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'SCL_HT/RP',SCL_HT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VCORE',VCORE)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VPHOT',VPHOT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  IF(VELTYPE .EQ. 6)THEN
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf2',VINF2)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta2',V_BETA2)
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
! Output abundance information.
!
	  DO ISPEC=1,NUM_SPECIES
	    CALL WR_ABUND_INFO_V2(SPECIES(ISPEC),AT_MASS(ISPEC),
	1           AT_ABUND(ISPEC),ABUND_SUM,
	1           MEAN_ATOMIC_WEIGHT,SOL_MASS_FRAC(ISPEC),LUMOD)
	  END DO
!
	  WRITE(LUMOD,'(A)')' '
	  IF(DO_CLUMP_MODEL)THEN
	    WRITE(LUMOD,'(A,A)')
	1            'Running clumped model: ',TRIM(CLUMP_LAW)
	    WRITE(LUMOD,'(A,1PE10.3)')
	1            'Filling factor at boundary is: ',CLUMP_FAC(1)
	    STRING=' '
	    NEXT_LOC=1
	    DO I=1,N_CLUMP_PAR
	      TEMP_CHAR='CL_P_'
	      WRITE(TEMP_CHAR(6:6),'(I1)')I
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,TEMP_CHAR,CLUMP_PAR(I))
	      IF(NEXT_LOC .GT. 80)THEN
	        WRITE(LUMOD,'(A)')TRIM(STRING)
	        STRING=' '
	        NEXT_LOC=1
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
	  WRITE(LUMOD,'(A,1PE10.3)')
	1            'Maximum correcion (%) on last iteration: ',MAXCH
	  CLOSE(LUMOD)
!
	END IF
!
! Check to see if corrections are reasonable.
!
	IF(MAXCH .GT. MAX_CHNG_LIM)THEN
	  WRITE(LUER,*)'Error - bad initial population guesses.'
	  WRITE(LUER,*)'Predicted changes are too large. '
	  WRITE(LUER,*)'New populations written to SCRTEMP file.'
	  WRITE(LUER,*)'Edit POINT1 file to recover older populations.'
	  STOP
	END IF
!
! 
9999	CONTINUE
!
! NB - Need GE because of NG acceleration.
!
	IF(LST_ITERATION)THEN
!
	  CALL TUNE(ITWO,'GIT')
	  CALL TUNE(ITHREE,' ')
!
! This file is a direct access file and contains the models
! output (i.e. T,density,population levels etc). No longer
! assumes that RVTJ exists (28-Sep-1990). Altered (1-Dec-1991) to be
! an asci file for CRAY/VAX compatibility. NB. Obs may be ZERO if program
! was inadvertently started during the last iteration.
!
	  CALL DATE_TIME(TIME)
!
	  CALL GEN_ASCI_OPEN(LU_POP,'RVTJ','UNKNOWN',' ',' ',IZERO,IOS)
	    FORMAT_DATE='15-Jun-2000'
	    WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	    WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	    WRITE(LU_POP,'(1X,A,T30,A)')'Program Date:',PRODATE
!
	    WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NC:',NC
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NP:',NP
	    WRITE(LU_POP,'(1X,A,T30,I6)')'NCF:',N_OBS
!
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'Mdot(Msun/yr):',
	1                                   RMDOT/3.02286D+23
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'L(Lsun):',LUM
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'H/He abundance:',AT_ABUND(1)
	    WRITE(LU_POP,'(1X,A,T30,L1)')'Was T fixed?:',RD_FIX_T
	    WRITE(LU_POP,'(1X,A,T30,A)')'Species naming convention:',
	1                                        NAME_CONVENTION
!
	    WRITE(LU_POP,'(A)')' Radius (10^10 cm)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')R
	    WRITE(LU_POP,'(A)')' Velocity (km/s)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')V
	    WRITE(LU_POP,'(A)')' dlnV/dlnr-1'
	    WRITE(LU_POP,'(1X,1P8E16.7)')SIGMA
	    WRITE(LU_POP,'(A)')' Electron density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')ED
	    WRITE(LU_POP,'(A)')' Temperature (10^4K)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')T
!
! These are written to OBSFLUX, and hence do not need to be output to
! RVTJ. They are not accessed by DISPGEN.
!
!	    WRITE(LU_POP,'(A)')' Continuum Frequencies (10^15 Hz)'
!	    WRITE(LU_POP,'(1X,1P6E20.12)')(OBS_FREQ(I),I=1,N_OBS)
!	    WRITE(LU_POP,'(A)')' Continuum Fluxes (Jy kpc^2)'
!	    WRITE(LU_POP,'(1X,1P8E16.7)')(OBS_FLUX(I),I=1,N_OBS)
!
	    WRITE(LU_POP,'(A)')' Rosseland Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(ROSS_MEAN(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Flux Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(FLUX_MEAN(I),I=1,ND)
!
! Compute the ion population at each depth.                          
! These are required when evaluation the occupation probabilities.
!
	    DO J=1,ND
	      POPION(J)=0.0D0
	      DO I=1,NT
	        IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	      END DO
	    END DO
	    WRITE(LU_POP,'(A)')' Atom Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POP_ATOM(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Ion Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POPION(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Mass Density (gm/cm^3)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(DENSITY(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Clumping Factor'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(CLUMP_FAC(I),I=1,ND)
!
	    WRITE(LU_POP,'(A)')' Hydrogen Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,1)
	    WRITE(LU_POP,'(A)')' Helium Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,2)
!
	    IF(ATM(1)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')TRIM(ION_ID(1))//' populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Hydrogen levels:',ATM(1)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(1))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(1)%XzV_F
	      WRITE(LU_POP,'(A)')' DHYD population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(1)%DXzV_F
	    END IF
!
	    IF(ATM(3)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' HeI populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium I levels:',ATM(3)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(3))//
	1                ' oscillator date:',ATM(3)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(3)%XzV_F
	      WRITE(LU_POP,'(A)')' DHeI population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(3)%DXzV_F
	    END IF
!
	    IF(ATM(4)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' He2 populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium II levels:',ATM(4)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(4))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(4)%XzV_F
	      WRITE(LU_POP,'(A)')' DHe2 population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(4)%DXzV_F
	    END IF
!
	  CLOSE(UNIT=LU_POP)
!
	  IF(SN_HYDRO_MODEL)THEN
	    CALL OUT_SN_POPS_V2('SN_HYDRO_FOR_NEXT_MODEL',SN_AGE_DAYS,ND,LUMOD)
	  END IF
! 
!
	  FORMAT_DATE='27-JAN-1992'
	  DO ISPEC=1,NUM_SPECIES
	    IF(POP_SPECIES(ND,ISPEC) .NE. 0)THEN
	      TMP_STRING='POP'//TRIM(SPECIES(ISPEC))
	      CALL GEN_ASCI_OPEN(LU_POP,TMP_STRING,'UNKNOWN',' ',' ',IZERO,IOS)
	      WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	      WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	      WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	      WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')
	1              TRIM(SPECIES(ISPEC))//'/He abundance:',AT_ABUND(ISPEC)
     	      WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,ISPEC)
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        CALL RITE_ASC( ATM(ID)%XzV_PRES, ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1              ATM(ID)%NXzV_F, ND,
	1              ATM(ID)%XzV_OSCDATE,TRIM(ION_ID(ID)),LU_POP)
	      END DO
	    END IF
	  END DO
! 
!
! Write out departure coefficients to ASCI file. 
! NB - 1 refers to dimension of DHYD (i.e. DHYD(1,nd)
!      1 refers to format for output.
!      1,NHY - For use with HeI.
!
	CALL EVAL_LTE_V4(DO_LEV_DISSOLUTION,ND)
!      
! GAM_SPECIES refers to the number of electrons arising from each species (eg
! carbon).
!
	DO ISPEC=1,NUM_SPECIES
	  GAM_SPECIES(1:ND,ISPEC)=0.0D0
	  FIRST=.TRUE.
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	    IF( ATM(ID)%XzV_PRES)THEN
	      CALL UPDATE_GAM( GAM_SPECIES(1,ISPEC),
	1          ATM(ID)%XzV_F, ATM(ID)%DXzV_F, ATM(ID)%ZXzV,
	1          ATM(ID)%NXzV_F,ND,
	1          ATM(ID+1)%XzV_PRES,FIRST)
	      TMP_STRING=TRIM(ION_ID(ID))//'OUT'
	      CALL WRITEDC_V2( ATM(ID)%XzV_F, ATM(ID)%XzVLTE_F,
	1          ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,IONE,
	1          R,T,ED,V,CLUMP_FAC,LUM,ND,
	1          TRIM(TMP_STRING),'DC',IONE)
	    END IF
	  END DO
	END DO
!
! We only output GAM when it is not zero.
!
	  CALL RITE_GAM_HEAD(R,ED,T,ND,LUIN,'GAMMAS')
	  DO ISPEC=1,NUM_SPECIES
	    CALL RITE_GAM_V2(POP_SPECIES(1,ISPEC),GAM_SPECIES(1,ISPEC),
	1                       AT_NO(ISPEC),SPECIES(ISPEC),ND,LUIN)
	  END DO
	  CLOSE(LUIN)
!
! If we have a time variability model, we output the model so it can be used by the
! next model in the time sequence. IREC is the ouput record, and is thus simply
! the current TIME_SEQ_NO.
!
	  IF(SN_MODEL .AND. TIME_SEQ_NO .NE. 0)THEN
!	    IREC=TIME_SEQ_NO
!	    CALL RITE_TIME_MODEL(R,V,SIGMA,POPS,IREC,L_FALSE,L_TRUE,
!	1               NT,ND,LUSCR,CHK)
	    CALL WRITE_SEQ_TIME_FILE_V1(SN_AGE_DAYS,ND,LUSCR)
	  END IF
!
	  CLOSE(UNIT=LUER)
	  CLOSE(UNIT=LU_SE)
	  STOP
!
! 
	ELSE
!
!*****************************************************************************
!*****************************************************************************
!                END OF MAIN ITERATION LOOP
!*****************************************************************************
!*****************************************************************************
!
! Close units 2 and 16 to force writing of information. We check if OUTGEN
! has been opened --- if not we are writing to the terminal and nothing
! needs to be done.
!
	  INQUIRE(FILE='OUTGEN',OPENED=FILE_OPEN)
	  IF(FILE_OPEN)THEN
	    CLOSE(UNIT=LUER)
	    CALL GEN_ASCI_OPEN(LUER,'OUTGEN','OLD','APPEND',' ',IZERO,IOS)
	  END IF
	  CLOSE(UNIT=LU_SE)
	  CALL GEN_ASCI_OPEN(LU_SE,'STEQ_VALS','OLD','APPEND',' ',IZERO,IOS)
!
! Adjust X-ray filling factors upwards if within a factor of 100 of convergence.
! We adjust MAXCH to ensure that NUM_ITS_TO_DO is not set to 1.
!
	  IF(XRAYS .AND. ADD_XRAYS_SLOWLY .AND. RD_LAMBDA .AND. MAXCH .LT. 100)THEN
	     IF(FILL_FAC_XRAYS_1 .NE. FILL_X1_SAV .OR.  FILL_FAC_XRAYS_2 .NE. FILL_X2_SAV)THEN
	       FILL_FAC_XRAYS_1=MIN(FILL_FAC_XRAYS_1*SLOW_XRAY_SCL_FAC,FILL_X1_SAV)
	       FILL_FAC_XRAYS_2=MIN(FILL_FAC_XRAYS_2*SLOW_XRAY_SCL_FAC,FILL_X2_SAV)
	       WRITE(LUER,*)'Have adjusted X-ray values closer to the desired values'
	       WRITE(LUER,*)'Current filling factor is (1st component)',FILL_FAC_XRAYS_1
	       WRITE(LUER,*)'Current filling factor is (2nd component)',FILL_FAC_XRAYS_2
	       WRITE(LUER,*)'Current on desired (1st component)=',FILL_FAC_XRAYS_1/FILL_X1_SAV
	       WRITE(LUER,*)'Current on desired (2nd component)=',FILL_FAC_XRAYS_2/FILL_X2_SAV
	       CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_1,'[XFI1_BEG]','VADAT',L_TRUE,L_FALSE,LUIN)
	       CALL UPDATE_KEYWORD(FILL_FAC_XRAYS_2,'[XFI2_BEG]','VADAT',L_FALSE,L_TRUE,LUIN)
	       MAXCH=100			!To force run to continue
	    ELSE
	       CALL UPDATE_KEYWORD(L_FALSE,'[XSLOW]','VADAT',L_TRUE,L_TRUE,LUIN)
	       IF(DO_LAMBDA_AUTO)RD_LAMBDA=.FALSE.
	    END IF
	  END IF
	  IF(INCL_ADVECTION .AND. ADVEC_RELAX_PARAM .LT. 1.0D0 .AND. MAXCH .LT. 100)THEN
	    COMPUTE_BA=.TRUE.
	    ADVEC_RELAX_PARAM=MIN(1.0D0,ADVEC_RELAX_PARAM*2.0D0)
	    WRITE(LUER,*)'Have adjuseted advection relaxation parameter to:',ADVEC_RELAX_PARAM
	    MAXCH=100
	  END IF
!
! If we have reached desired convergence, we do one final loop
! so as to write out all relevant model data.  
!
	  IF( ( (RD_LAMBDA .AND. .NOT. DO_LAMBDA_AUTO) .OR. .NOT. LAMBDA_ITERATION) .AND.
	1      MAXCH .LT. EPS .AND. NUM_ITS_TO_DO .NE. 0)THEN
	      NUM_ITS_TO_DO=1
	  ELSE
!
! If we are USING a fixed J, autmatically switched to variable J when convergence achieved.
! We switch when the convergence is 20%.
!
	    IF(USE_FIXED_J .AND. DO_LAMBDA_AUTO .AND. RD_LAMBDA .AND. MAXCH .LT. 50.0D0)THEN
	       USE_FIXED_J=.FALSE.
	       CALL UPDATE_KEYWORD(L_FALSE,'[USE_FIXED_J]','VADAT',L_TRUE,L_TRUE,LUIN)
	       IF(DO_GREY_T_AUTO)THEN
	         CALL GREY_T_ITERATE(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1               LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	         CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1               LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	       END IF
	       I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	       CALL WRITE_DIRECT_INFO_V3(NDEXT,I,'20-Aug-2000','EDDFACTOR',LU_EDD)
	       COMPUTE_EDDFAC=.TRUE.
	       OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',ACCESS='DIRECT',STATUS='REPLACE',RECL=I)
	       WRITE(LU_EDD,REC=1)0; WRITE(LU_EDD,REC=2)0
	       WRITE(LU_EDD,REC=3)0; WRITE(LU_EDD,REC=4)0
	       T1=0.0
	       WRITE(LU_EDD,REC=5)T1
	       COHERENT_ES=.TRUE.
!
! If DONE_HYDRO_REVISION is TRUE, we do a LAMBDA iteration immediately afterwoods.
! LAMBDA_ITERATION has already been set.
!
	    ELSE IF(DONE_HYDRO_REVISION)THEN
!
! If RD_LABDA is TRUE., switch to full iteration when convergence has been achieved.
! We switch when the convergence is 20%.
!
	    ELSE IF(DO_LAMBDA_AUTO .AND. RD_LAMBDA .AND. MAXCH .LT. 50.0D0)THEN
	       CALL UPDATE_KEYWORD(L_FALSE,'[DO_LAM_IT]','IN_ITS',L_TRUE,L_TRUE,LUIN)
	       RD_LAMBDA=.FALSE.
	       LAMBDA_ITERATION=.FALSE.
	       FIXED_T=RD_FIX_T
	    END IF
!
! Check to see if the user has changed IN_ITS to modify the number of iterations being
! undertaken. If the file has not been modified, no action will be taken. The use may
! also change whether LAMBDA iterations are being done. At least one final iteration
! will be undertaken.
!
	    CALL GEN_ASCI_OPEN(LUSCR,'MODEL_SCR','UNKNOWN',' ',' ',IZERO,IOS)
	    IF(IOS .EQ. 0)CALL GEN_ASCI_OPEN(LUIN,'IN_ITS','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN  
	       WRITE(LUER,*)'Error opening IN_ITS or MODEL_SCR in CMFGEN, IOS=',IOS
	       WRITE(LUER,*)'Error occurs at the end of CMFGEN_SUB.'
	       WRITE(LUER,*)'Error will be ignored.'
	       GOTO 20000
	    END IF
	    OLD_RD_LAMBDA=RD_LAMBDA
	    I=NUM_ITS_RD
	    CALL RD_INT(NUM_ITS_RD,'NUM_ITS',LUIN,LUSCR,'Number of iterations to perform')
	    CALL RD_LOG(RD_LAMBDA,'DO_LAM_IT',LUIN,LUSCR,'Do LAMBDA iterations ?')
	    CLOSE(UNIT=LUIN)
	    NUM_ITS_TO_DO=NUM_ITS_TO_DO+(NUM_ITS_RD-I)
	    IF(NUM_ITS_TO_DO .LE. 0)NUM_ITS_TO_DO=1
	    IF(RD_LAMBDA)THEN
	      LAMBDA_ITERATION=.TRUE.
	      FIX_IMPURITY=.FALSE.
	      FIXED_T=.TRUE.
!
! Don't wish to chnage ITERATION cycle values unless we have to.
!
	    ELSE IF(OLD_RD_LAMBDA)THEN
	      FIX_IMPURITY=RD_FIX_IMP
	      FIXED_T=RD_FIX_T
	    END IF
	    CLOSE(UNIT=LUSCR,STATUS='DELETE')
	  END IF
	  GOTO 20000				!Begin another iteration
	END IF
!
	END
