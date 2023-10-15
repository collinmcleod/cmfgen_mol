	MODULE OPAC_MOD
	IMPLICIT NONE
!
	REAL(10), ALLOCATABLE :: ETA(:)
	REAL(10), ALLOCATABLE :: ETA_CONT(:)
	REAL(10), ALLOCATABLE :: ETA_C_EVAL(:)
	REAL(10), ALLOCATABLE :: ETA_NOSCAT(:)
	REAL(10), ALLOCATABLE :: ETA_NOSCAT_EVAL(:)
	REAL(10), ALLOCATABLE :: ETA_CLUMP(:)
	REAL(10), ALLOCATABLE :: ETA_MECH(:)
!
	REAL(10), ALLOCATABLE :: CHI(:)
	REAL(10), ALLOCATABLE :: CHI_CONT(:)
	REAL(10), ALLOCATABLE :: CHI_C_EVAL(:)
	REAL(10), ALLOCATABLE :: CHI_NOSCAT(:)
	REAL(10), ALLOCATABLE :: CHI_NOSCAT_EVAL(:)
	REAL(10), ALLOCATABLE :: CHI_CLUMP(:)
	REAL(10), ALLOCATABLE :: TCHI(:)
!
	REAL(10), ALLOCATABLE :: CHI_SCAT(:)
	REAL(10), ALLOCATABLE :: CHI_RAY(:)
	REAL(10), ALLOCATABLE :: ESEC(:)
	REAL(10), ALLOCATABLE :: ESEC_CLUMP(:)
	REAL(10), ALLOCATABLE :: CHI_SCAT_CLUMP(:)
!
! Use to relate opacities at the current frequency to that at last frequency.
! They store the opacities at the previous frequency. Note that CHI_SCAT is
! only constant when we have pure electron scattering --- it varies if we allow
! for Rayleigh scattering.
!
	REAL(10), ALLOCATABLE ::  CHI_PREV(:)
	REAL(10), ALLOCATABLE ::  CHI_NOSCAT_PREV(:)
	REAL(10), ALLOCATABLE ::  CHI_SCAT_PREV(:)
	REAL(10), ALLOCATABLE ::  ETA_PREV(:)
!
! To allow the variation of non-coherent electron scattering to be treated
! in a partially coherent approximation.
!
	REAL(10), ALLOCATABLE :: ES_COH_VEC(:)
	REAL(10), ALLOCATABLE :: THETA(:)
	REAL(10), ALLOCATABLE :: ZETA(:)
	REAL(10), ALLOCATABLE :: SOURCE(:)
!
	REAL(10), ALLOCATABLE :: EMHNUKT(:)
	REAL(10), ALLOCATABLE :: EMHNUKT_CONT(:)
!
	REAL(10), ALLOCATABLE :: XRAY_LUM_0P1(:)
	REAL(10), ALLOCATABLE :: XRAY_LUM_1KEV(:)
	REAL(10), ALLOCATABLE :: XRAY_LUM_TOT(:)
!
	REAL(10), ALLOCATABLE :: VCHI(:,:)              !Variation of CHI array.
        REAL(10), ALLOCATABLE :: VETA(:,:)              !Variation of ETA array.
        REAL(10), ALLOCATABLE :: VCHI_SAV(:,:)
        REAL(10), ALLOCATABLE :: VETA_SAV(:,:)
        REAL(10), ALLOCATABLE :: VCHI_ALL(:,:)
	REAL(10), ALLOCATABLE :: VCHI_ALL_SAV(:,:)
        REAL(10), ALLOCATABLE :: VETA_ALL(:,:)
	REAL(10), ALLOCATABLE :: VETA_ALL_SAV(:,:)
!
	END MODULE OPAC_MOD
