!
! Module containing DATA for charge exchange reactions.
!
	MODULE CHG_EXCH_MOD_V3
	  INTEGER*4 N_CHG_RD
	  INTEGER*4 N_CHG
	  INTEGER*4 LUER
!
! Altered 10-Sep-2000 : Bug fix. THI_CHG was not beeing set is not read in.
! Altered 01-Oct-1998 : Fitting range and FORMAT date is now read in. 
!                       An alternative level name can also be specified
!                       (via {}).
! Altered 20-Aug-1998
! Created 24-Jun-1998
!
! Reaction data
!
          INTEGER*4, ALLOCATABLE :: TYPE_CHG_RD(:)
          REAL*8, ALLOCATABLE :: COEF_CHG_RD(:,:)
          REAL*8, ALLOCATABLE :: TLO_CHG_RD(:)
          REAL*8, ALLOCATABLE :: THI_CHG_RD(:)
          LOGICAL, ALLOCATABLE :: CHG_INCLUDED_RD(:)
	  CHARACTER*8, ALLOCATABLE :: SPEC_ID_CHG_RD(:,:)
	  CHARACTER*30, ALLOCATABLE :: LEV_NAME_CHG_RD(:,:)
	  CHARACTER*30, ALLOCATABLE :: ALT_LEV_NAME_CHG_RD(:,:)
!
          INTEGER*4, ALLOCATABLE :: TYPE_CHG(:)
          REAL*8, ALLOCATABLE :: COEF_CHG(:,:)
          REAL*8, ALLOCATABLE :: G_CHG(:,:)
          REAL*8, ALLOCATABLE :: TLO_CHG(:)
          REAL*8, ALLOCATABLE :: THI_CHG(:)
	  CHARACTER*8, ALLOCATABLE :: SPEC_ID_CHG(:,:)
	  CHARACTER*30, ALLOCATABLE :: LEV_NAME_CHG(:,:)
	  CHARACTER*30, ALLOCATABLE :: ALT_LEV_NAME_CHG(:,:)
!
! Arrays required to update STEQ (statistical equilibrium) and BA (variation
! of STEQ) arrays.
!
          REAL*8, ALLOCATABLE :: Z_CHG(:,:)
          REAL*8, ALLOCATABLE :: AI_AR_CHG(:,:)
          REAL*8, ALLOCATABLE :: dlnAI_AR_CHG_dlnT(:,:)
          REAL*8, ALLOCATABLE :: COOL_CHG(:,:)
!
	  INTEGER*4, ALLOCATABLE ::  ID_ION_CHG(:,:)
	  INTEGER*4, ALLOCATABLE ::  LEV_IN_POPS_CHG(:,:)
	  INTEGER*4, ALLOCATABLE ::  LEV_IN_ION_CHG(:,:)
	  LOGICAL, ALLOCATABLE ::  CHG_REACTION_AVAILABLE(:)
!
	  INTEGER*4, PARAMETER :: N_COEF_MAX=5
!
	  LOGICAL INITIALIZE_ARRAYS
	  LOGICAL DO_CHG_EXCH
!
	END MODULE CHG_EXCH_MOD_V3
