!	Module for all the vectors and information that will be read out
!	from the NUC_DECAY_DATA file, which will then be stored into a structure
!
	MODULE GAMMA_NUC_DECAY_V2
	USE SET_KIND_MODULE
!
        INTEGER, PARAMETER :: SPECIES_MAX=50
	INTEGER, PARAMETER :: LINE_MAX=200
	INTEGER :: N_GAMMA, N_SPECIES
	INTEGER :: NGAM  ! Number of gamma-ray lines used. Repeated lines are ignored
	CHARACTER(LEN=4), DIMENSION(LINE_MAX) :: SPECIES_VEC(LINE_MAX)
        REAL(KIND=LDP) :: GAMMA_E_VEC(LINE_MAX), PROB_VEC(LINE_MAX)
        REAL(KIND=LDP) :: E_KIN_VEC(LINE_MAX), ATOM_MASS_VEC(LINE_MAX)
!
	TYPE GAMMA_NUC_ISO
	 CHARACTER(LEN=10) :: SPECIES
	 INTEGER :: NUM_GAMMA 		! # of gamma lines read in from NUC_DECAY_DATA
	 REAL(KIND=LDP) :: ATOMIC_MASS
	 REAL(KIND=LDP), DIMENSION(30) :: E_GAMMA
	 REAL(KIND=LDP), DIMENSION(30) :: PROB
	 REAL(KIND=LDP) :: KIN_ENERGY
	END TYPE
!
	TYPE (GAMMA_NUC_ISO) :: GAM_ISO(SPECIES_MAX)
        END MODULE
