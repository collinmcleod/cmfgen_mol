!
! Subroutine to add BA_PAR to to the full BA matrix. At the same time
! BA_PAR is zeroed.
!
! Routine may be used to update both BAION and BA.
! To update BA, pass NT for NION.
!
! This routine should not be called on very frequency. Rather it should
! be called every 50 or so frequencies. In this way the BA and BAION matrices
! should suffer less cancelation effects due due to the addition of large
! positive and negative terms.
!
! Utilizing fact that consecutive frequency terms should be correlated, and
! hence similar in size. While some minor cancellation, should be much less
! then adding to full BA matrix in which terms have arbitrary size.
!
	SUBROUTINE ADD_PAR_TO_FULL_V2(NION,DIAG_INDX)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered: 09-Apr-2001  Changed to utilize STEQ_DATA_MOD 
! Created: 28-Feb-1995
!
	INTEGER*4 NION
	INTEGER*4 DIAG_INDX
!
	INTEGER*4 ID,I
!
!	WRITE(83,*)' '
!	WRITE(83,*)SE(24)%BA(1,1,DIAG_INDX,1),SE(24)%BA(1,15,DIAG_INDX,1)
!	WRITE(83,*)SE(24)%BA_PAR(1,1,1),SE(24)%BA_PAR(1,15,1)
	DO ID=1,NION
	  SE(ID)%BA(:,:,DIAG_INDX,:)=SE(ID)%BA(:,:,DIAG_INDX,:) + SE(ID)%BA_PAR
          SE(ID)%BA_PAR=0.0D0
	END DO
	BA_T(:,DIAG_INDX,:)=BA_T(:,DIAG_INDX,:)+BA_T_PAR
	BA_T_PAR=0.0D0
!
!        I=SE(4)%LNK_TO_IV(1101)
!	WRITE(99,*)SE(4)%BA(1,I,DIAG_INDX,1),SE(4)%BA(2,I,DIAG_INDX,1)
!
	RETURN
	END
