!
! Routine to create a LOGICAL vector indicating which variables are
! classified as IMPORTANT VARIABLES (IVs). These variables are considered
! as important for ALL species.
!
! If NIV_XzV is zero, we assume a maximum of 10 levels (for each species)
! are important
!
	SUBROUTINE SET_IMP_VEC(IMP_VAR,NXzV,NIV_XzV,EQxZV,NT,XzV_PRES)
	IMPLICIT NONE
!
! Created 17-Mar-2001
!
	INTEGER*4 NXzV		!Number of levels in species
	INTEGER*4 NIV_XzV	!Number of Important levels in species
	INTEGER*4 EQXzV		!Equation corresponding to ground state of species.
	INTEGER*4 NT		!Total number of varaiables (all ions and species)
	LOGICAL XzV_PRES
	LOGICAL IMP_VAR(NT)
!
	INTEGER*4 K
	LOGICAL, SAVE :: FIRST=.TRUE.
!
! Initialize IMP_VAR. We set it to FALSE for all variables except
! ED and T. 
!
	IF(FIRST)THEN
	  FIRST=.FALSE.
          IMP_VAR(1:NT)=.FALSE.
	  IMP_VAR(NT-1)=.TRUE.
	  IMP_VAR(NT)=.TRUE.
	END IF
!
! Now set it to true for the variables of the current species,
! including the associated ion.
!
	IF(.NOT. XzV_PRES)RETURN
        K=MIN(NIV_XzV,NXzV)
        IF(K .EQ. 0)RETURN
	IMP_VAR(EQXzV:EQXzV+K-1)=.TRUE.
	IMP_VAR(EQXzV+NXzV)=.TRUE.
!
	RETURN
	END
