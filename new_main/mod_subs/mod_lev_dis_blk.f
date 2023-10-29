!
! Module to store data necssary to compute occupation probabilities
! etc associated with level dissolution.
!
! The vectors should be defined as follows:
!
!  B_LEV_DIS=( 8.589D+14*(POPION**0.333333D0)/ED)**1.5D0
!  A_LEV_DIS=0.0009*(ED**0.1667)/SQRT(T)
!  X_LEV_DIS=(1+A_LEV_DIS)**3.15
!
! NB: B_LEV_DIS is undefined for ED=0, but this should be of no
!       practical concern, if we specifically handle the case when
!       level dissoultion is switched off (corresponds to ED=0)
!
! Based on the include file LEV_DIS_BLK.INC. Has dynamic allocation
! and variable indicating whether level dissolution is switched on
! in now explicitly included in the data block.
!
!
! Altered 23-Oct-2016 :: Inserted PHOT_DIS_PARAMETER. This limits how far from edge
!                          we have non-zero photo-ionization cross-sections.
!
	MODULE MOD_LEV_DIS_BLK
	USE SET_KIND_MODULE
	REAL(KIND=LDP), ALLOCATABLE :: B_LEV_DIS(:)
	REAL(KIND=LDP), ALLOCATABLE :: A_LEV_DIS(:)
	REAL(KIND=LDP), ALLOCATABLE :: X_LEV_DIS(:)
	REAL(KIND=LDP), PARAMETER :: PHOT_DIS_PARAMETER=1.0D-04		!Was 1.0E-06 (29-May-2019)
	LOGICAL MOD_DO_LEV_DIS
!
	END MODULE MOD_LEV_DIS_BLK
!
! Subroutine to:
!          (1) Allocate desired memory if not done already.
!          (2) Compute level dissolution constants.
!
	SUBROUTINE COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) POPION(ND)
	REAL(KIND=LDP) T(ND)
	LOGICAL DO_LEV_DIS
!
	IF(.NOT. ALLOCATED(B_LEV_DIS))THEN
	  ALLOCATE( B_LEV_DIS(ND) )
	  ALLOCATE( A_LEV_DIS(ND) )
	  ALLOCATE( X_LEV_DIS(ND) )
	END IF
!
	MOD_DO_LEV_DIS=DO_LEV_DIS
	IF(DO_LEV_DIS)THEN
	  B_LEV_DIS=( 8.589D+14*(POPION**0.333333D0)/ED)**1.5D0
	  A_LEV_DIS=0.0009D0*(ED**0.1667D0)/SQRT(T)
	  X_LEV_DIS=(1.0D0+A_LEV_DIS)**3.15D0
	ELSE
	  B_LEV_DIS=1.0D+30		!Large number
	  A_LEV_DIS=0.0D0
	  X_LEV_DIS=1.0D0
	END IF
!
	RETURN
	END
