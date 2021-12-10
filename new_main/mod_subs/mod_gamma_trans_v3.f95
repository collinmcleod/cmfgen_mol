	MODULE MOD_GAMMA_V3
	IMPLICIT NONE
!
	INTEGER :: NF_GRID_PTS
        INTEGER, DIMENSION(:), ALLOCATABLE :: KDW_INDEX
        INTEGER, DIMENSION(:), ALLOCATABLE :: NU_END
!
	REAL*8 :: GAMRAY_LUM
!
	REAL*8, DIMENSION(:), ALLOCATABLE :: NU_GRID_VEC
	REAL*8, DIMENSION(:), ALLOCATABLE :: ETA_NORM_CONST
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAM_ENERGY_DEP
	REAL*8, DIMENSION(:), ALLOCATABLE :: NU_VEC !vector of gamma frequencies used from nuc_decay_data file
        REAL*8, DIMENSION(:), ALLOCATABLE :: WORK1
	REAL*8, DIMENSION(:), ALLOCATABLE :: CROSS_KN
	REAL*8, DIMENSION(:), ALLOCATABLE :: XCHI,E_SCAT_CHI
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAM_OPAC,GAM_OPAC_CLUMP,GAM_OPAC_COPY
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAM_E_ME ! will be h*nu/mec^2 thus unitless
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAM_NU_MAX ! will be h*nu/mec^2 thus unitless
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMRAY_LUM_VEC
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMRAY_EMISS
	REAL*8, DIMENSION(:), ALLOCATABLE :: ED_TOT
	REAL*8, DIMENSION(:), ALLOCATABLE :: DECAY_KIN_E
	REAL*8, DIMENSION(:), ALLOCATABLE :: TOTAL_DECAY_LUM
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMMA_TAU
	REAL*8, DIMENSION(:), ALLOCATABLE :: XRAY_TAU
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMMA_FLUX
	REAL*8, DIMENSION(:), ALLOCATABLE :: NU_VEC_15
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMMA_NU_STORE
	REAL*8, DIMENSION(:), ALLOCATABLE :: GAMMA_LAM
!
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: KLEIN_ARRAY
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAM_ETA !array that gets passed as (mu,depth)
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAM_INT
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: ETA_ISO
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: E_SCAT_ARRAY
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: CHI_ARRAY_ABS
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAMMA_J
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAM_ETA_MUAVG
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAM_ETA_NUAVG
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAMMA_STORE
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: GAM_I_STORE
!
	REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: GAM_ETA_SCAT
	REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: GAM_ETA_SCAT_T  ! Will be used as the transpose array
	REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: GAM_I
	REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: GAM_I_T  	 ! Will be used as the transpose array
!
	REAL*8, PARAMETER :: PLANCK = 4.135668D-21 ! UNITS OF MeV*s SINCE PHOTON ENERGIES IN MeV
!	REAL*8, PARAMETER :: SOL = 299792 ! Units of km/s
	REAL*8, PARAMETER :: home = 8.09330118D-21	! home = h/mc^2 in units of seconds
!
	END MODULE
