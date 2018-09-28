!
! Subroutine designed to replace the CMFGEN grid with the new grid.
! We use a subroutine to so we can access the grid in CMFGEN through MOD_CMFGEN.
!
	SUBROUTINE SET_NEW_GRID_V2(REV_R,REV_T,REV_V,REV_SIGMA,REV_ED,REV_CHI_ROSS,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 22-Jun-2016: Changed to V2. Added REV_T to call. T no reset.
!                         Commnets added 16-Aug-2015 (cur_hmi).
! Created 26-Jan-2007
!
	INTEGER ND
	REAL*8 REV_R(ND)
	REAL*8 REV_T(ND)
	REAL*8 REV_V(ND)
	REAL*8 REV_SIGMA(ND)
	REAL*8 REV_ED(ND)
	REAL*8 REV_CHI_ROSS(ND)
!
	R=REV_R
	T=REV_T
	V=REV_V
	SIGMA=REV_SIGMA
	ED=REV_ED
	ROSS_MEAN=REV_CHI_ROSS
!
	RETURN
	END
