! Module to store the created mu grids into a structure for every depth point
!
! Created July 25, 2014
!-----------------------------------------------------------------------------
	MODULE GAM_MU_MOD
	IMPLICIT NONE
!
	TYPE MU_GRID
	  INTEGER :: MU_PTS
	  INTEGER :: MON_MU_PTS
	  REAL*8, DIMENSION(:), ALLOCATABLE :: MU_VECTOR
	  REAL*8, DIMENSION(:), ALLOCATABLE :: MON_MU
	  REAL*8 :: dMU
	END TYPE
!
	TYPE (MU_GRID), ALLOCATABLE :: R_MU(:)
	END MODULE
