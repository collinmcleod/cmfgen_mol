!	Module for all the vectors and information that will be read out 
!	from the NUC_DECAY_DATA file, which will then be stored into a structure
!
	MODULE GAMMA_NUC_DECAY_V2
!
        INTEGER, PARAMETER :: SPECIES_MAX=50
	INTEGER, PARAMETER :: LINE_MAX=200
	INTEGER :: N_GAMMA, N_SPECIES
	INTEGER :: NGAM  ! Number of gamma-ray lines used. Repeated lines are ignored
	CHARACTER(LEN=4), DIMENSION(LINE_MAX) :: SPECIES_VEC(LINE_MAX)
        REAL*8 :: GAMMA_E_VEC(LINE_MAX), PROB_VEC(LINE_MAX)
        REAL*8 :: E_KIN_VEC(LINE_MAX), ATOM_MASS_VEC(LINE_MAX)
!
	TYPE GAMMA_NUC_ISO
	 CHARACTER(LEN=10) :: SPECIES
	 INTEGER :: NUM_GAMMA 		! # of gamma lines read in from NUC_DECAY_DATA
	 REAL*8 :: ATOMIC_MASS
	 REAL*8, DIMENSION(30) :: E_GAMMA
	 REAL*8, DIMENSION(30) :: PROB
	 REAL*8 :: KIN_ENERGY
	END TYPE
!
	TYPE (GAMMA_NUC_ISO) :: GAM_ISO(SPECIES_MAX)
        END MODULE
