	MODULE OPAC_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP), ALLOCATABLE :: ETA(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_CONT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_C_EVAL(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_NOSCAT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_NOSCAT_EVAL(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_CLUMP(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_MECH(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: CHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_CONT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_C_EVAL(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_NOSCAT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_NOSCAT_EVAL(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_CLUMP(:)
	REAL(KIND=LDP), ALLOCATABLE :: TCHI(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: CHI_SCAT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_CLUMP(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_SCAT_CLUMP(:)
!
! Use to relate opacities at the current frequency to that at last frequency.
! They store the opacities at the previous frequency. Note that CHI_SCAT is
! only constant when we have pure electron scattering --- it varies if we allow
! for Rayleigh scattering.
!
	REAL(KIND=LDP), ALLOCATABLE ::  CHI_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE ::  CHI_NOSCAT_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE ::  CHI_SCAT_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE ::  ETA_PREV(:)
!
! To allow the variation of non-coherent electron scattering to be treated
! in a partially coherent approximation.
!
	REAL(KIND=LDP), ALLOCATABLE :: ES_COH_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: THETA(:)
	REAL(KIND=LDP), ALLOCATABLE :: ZETA(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: EMHNUKT(:)
	REAL(KIND=LDP), ALLOCATABLE :: EMHNUKT_CONT(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_0P1(:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_1KEV(:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_LUM_TOT(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: VCHI(:,:)              !Variation of CHI array.
        REAL(KIND=LDP), ALLOCATABLE :: VETA(:,:)              !Variation of ETA array.
        REAL(KIND=LDP), ALLOCATABLE :: VCHI_SAV(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: VETA_SAV(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: VCHI_ALL(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VCHI_ALL_SAV(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: VETA_ALL(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VETA_ALL_SAV(:,:)
!
	END MODULE OPAC_MOD
