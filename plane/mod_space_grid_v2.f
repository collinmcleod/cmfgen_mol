!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
	MODULE MOD_SPACE_GRID_V2
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module to define the r-p-z ray or characteristic ray variables.
! This save defining them in CMFGEN and having to pass them
! to subroutines.  Variables are made allocatable so array sizes can
! be passed and array sizes created at runtime.
!
! written  4/8/97  DLM  Previously used f77 and varaibles could not be
!                       allocated.
!
! altered 5/22/97  DLM  Added b_p and b_m for relativistic transfer terms
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Variables for both types of grids
!
	TYPE RAY_DATA
!
	  REAL(10), ALLOCATABLE :: MU(:)
	  REAL(10), ALLOCATABLE :: Z(:)
	  REAL(10), ALLOCATABLE :: R_RAY(:)
!
! Variables for characteristics
!
	  REAL(10), ALLOCATABLE :: S_P(:)                    !Path Length in positive direction
	  REAL(10), ALLOCATABLE :: S_M(:)                    !Path Length in negative direction
	  REAL(10), ALLOCATABLE :: MU_P(:)                   !Angle [cos(theta)] in positive direction
	  REAL(10), ALLOCATABLE :: MU_M(:)                   !Angle [cos(theta)] in negative direction
!
! Variables for frequency independent parts of advection and abberation terms
!
	  REAL(10), ALLOCATABLE :: B_P(:)                    !positive direction
	  REAL(10), ALLOCATABLE :: B_M(:)                    !negative direction
!
	  REAL(10), ALLOCATABLE :: I_P(:)
	  REAL(10), ALLOCATABLE :: I_P_PREV(:)
	  REAL(10), ALLOCATABLE :: I_P_SAVE(:)
	  REAL(10), ALLOCATABLE :: I_M(:)
	  REAL(10), ALLOCATABLE :: I_M_PREV(:)
	  REAL(10), ALLOCATABLE :: I_M_SAVE(:)
!
	  INTEGER, ALLOCATABLE :: LNK(:)		!Indicates link to original grid
!
	  REAL(10), ALLOCATABLE :: I_IN_BND_STORE(:) 
	  REAL(10) P_RAY				      	!P for ray
	  REAL(10) FREQ_CONV_FAC				!Frequency conversion factor for hollow core.

	  INTEGER  NZ	 	                    	!number of grid points along a p-ray
!
	END TYPE RAY_DATA
!
	TYPE (RAY_DATA) RAY(500)
!
! Arrays for "optical depth"
!
	REAL(10), ALLOCATABLE :: TAU(:)
	REAL(10), ALLOCATABLE :: DTAU(:)
!
! Defined along a ray and on the CMFGEN grid.
!
	REAL(10), ALLOCATABLE :: I_P_GRID(:)
	REAL(10), ALLOCATABLE :: I_M_GRID(:)
!
! Quadrature weights for characteristics. These are defined on the
! CMFGEN grid.
!
	REAL(10), ALLOCATABLE :: JQW_P(:,:)                  !Positive mu
	REAL(10), ALLOCATABLE :: HQW_P(:,:)
	REAL(10), ALLOCATABLE :: KQW_P(:,:)
	REAL(10), ALLOCATABLE :: NQW_P(:,:)
!
	REAL(10), ALLOCATABLE :: JQW_M(:,:)                  !NEGATIVE MU
	REAL(10), ALLOCATABLE :: HQW_M(:,:)
	REAL(10), ALLOCATABLE :: KQW_M(:,:)
	REAL(10), ALLOCATABLE :: NQW_M(:,:)
!
	REAL(10), ALLOCATABLE :: R_EXT_SAV(:)
	REAL(10), ALLOCATABLE :: FREQ_STORE(:)
!
	INTEGER CUR_LOC
	INTEGER N_STORE
	INTEGER ND_SAV
	INTEGER ND_EXT_SAV
	INTEGER NP_SAV
	LOGICAL, SAVE :: RAY_POINTS_INSERTED=.FALSE.
!
      END MODULE MOD_SPACE_GRID_V2
