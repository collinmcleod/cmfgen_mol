!
! Data module for SET_PROF. It is used so that we don't have to compute
! the entire STARK profile (necessary for convolution purposes) on
! each call.
!
	MODULE PROF_MOD
	IMPLICIT NONE
!
! Created 21-Sep-1999
!
	INTEGER NSTORE_MAX
	INTEGER NFREQ_MAX
	INTEGER LUER
	REAL*8 C_KMS
	REAL*8 PI
!
	REAL*8, ALLOCATABLE :: PROF_STORE(:,:,:)
	REAL*8, ALLOCATABLE :: NU_STORE(:,:)
	REAL*8, ALLOCATABLE :: AMASS_STORE(:)
	REAL*8, ALLOCATABLE :: NU_ZERO_STORE(:)
!         
	INTEGER, ALLOCATABLE :: NL_STORE(:)
	INTEGER, ALLOCATABLE :: NUP_STORE(:)
	INTEGER, ALLOCATABLE :: LST_FREQ_LOC(:)
	INTEGER, ALLOCATABLE :: NF_STORE(:)
!
	REAL*8, ALLOCATABLE :: PR_GRIEM(:)
	REAL*8, ALLOCATABLE :: DWS_GRIEM(:)
!
	LOGICAL, ALLOCATABLE :: STORE_AVAIL(:)
!
	END MODULE PROF_MOD
