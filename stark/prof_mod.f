!
! Data module for SET_PROF. It is used so that we don't have to compute
! the entire STARK profile (necessary for convolution purposes) on
! each call.
!
	MODULE PROF_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 21-Sep-1999
!
	INTEGER NSTORE_MAX
	INTEGER NFREQ_MAX
	INTEGER LUER
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) PI
!
	REAL(KIND=LDP), ALLOCATABLE :: PROF_STORE(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NU_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AMASS_STORE(:)
	REAL(KIND=LDP), ALLOCATABLE :: NU_ZERO_STORE(:)
!
	INTEGER, ALLOCATABLE :: NL_STORE(:)
	INTEGER, ALLOCATABLE :: NUP_STORE(:)
	INTEGER, ALLOCATABLE :: LST_FREQ_LOC(:)
	INTEGER, ALLOCATABLE :: NF_STORE(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: PR_GRIEM(:)
	REAL(KIND=LDP), ALLOCATABLE :: DWS_GRIEM(:)
!
	LOGICAL, ALLOCATABLE :: STORE_AVAIL(:)
!
	INTEGER NVGT_ST
	TYPE PROF_STORE_VGT
	  REAL(KIND=LDP), ALLOCATABLE :: VGT_PROF(:,:)
	  INTEGER NF
	END TYPE PROF_STORE_VGT
	INTEGER, ALLOCATABLE :: VGT_POINTER(:)
!
	TYPE (PROF_STORE_VGT), ALLOCATABLE :: VGT(:)
	
!
	END MODULE PROF_MOD
