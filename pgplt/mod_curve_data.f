	MODULE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Altered 22-Jul-2022: Extra variables added.
! Altered 24-Jan-2013: Vectors in CURVE_DATA in no longer a pointer.
!                      This required ASSOICATED test to be replaced by ALLOCATED
!                         test in several routines.
! Altered 28-Mar-2003: LU_ER changed from 5 to 6.
!
	INTEGER, PARAMETER :: MAX_PLTS=50
	INTEGER, PARAMETER :: LU_ER=6
	INTEGER IOS
!
	INTEGER NPLTS
	INTEGER NPTS(MAX_PLTS)
	LOGICAL   ERR(MAX_PLTS)
!
	TYPE CURVE_DATA
	  REAL*4, ALLOCATABLE :: XVEC(:)
	  REAL*4, ALLOCATABLE :: DATA(:)
	  REAL*4, ALLOCATABLE :: EMIN(:)
	  REAL*4, ALLOCATABLE :: EMAX(:)
	  CHARACTER(LEN=50) CURVE_ID
	END TYPE CURVE_DATA
!
	DATA NPTS/MAX_PLTS*0/
	DATA ERR/MAX_PLTS*.FALSE./
	DATA NPLTS/0/
!
	INTEGER MAX_TIT_LENGTH/60/
	INTEGER N_TIT_SET/0/
	INTEGER N_HEADER/0/			!Number of non curve identifiers
        INTEGER, PARAMETER :: N_TITLE=10
        CHARACTER(LEN=200)  TITLE(N_TITLE)	
	CHARACTER(LEN=200)  STORED_TITLE(N_TITLE)
	CHARACTER(LEN=200)  CURVE_TITLE(N_TITLE)
!
	TYPE (CURVE_DATA) CD(MAX_PLTS)
!
	END MODULE MOD_CURVE_DATA
