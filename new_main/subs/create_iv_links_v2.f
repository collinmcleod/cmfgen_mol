!
! Subroutine to create a 2 dimensional array which, for each equation,
! provides a link between a variable, and its location in the "small"
! BA matrix.
!
! This routines should be called successively for each ion of each
! species. It also returns the maximum # of variables (NIV) that have
! been allocated for any equation (including earlier allocations).
!
      SUBROUTINE CREATE_IV_LINKS_V2(NT,NION)
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      IMPLICIT NONE
!
! Created 17-Mar-2001
!
      INTEGER*4 NT
      INTEGER*4 NION            !Total number of species (Phot/recom balance equations)
!
      INTEGER*4 NIV		!Maximum number of IV allocations.
      INTEGER*4 I,J,ID
      LOGICAL LOC_IMP_VAR(NT)
!
	WRITE(93,*)IMP_VAR
!
      DO ID=1,NION
!
	IF(ATM(ID)%XzV_PRES)THEN
          ALLOCATE (SE(ID)%LNK_TO_IV(NT))
          SE(ID)%LNK_TO_IV(:)=0			!1:NT
!
! We first set the IV's for the ground state. For each ion we include all levels,
! plus the first excited state. We set LOC_IMP_VAR to false, so that we know that
! the variable has been included (or does not need to be).
!
          LOC_IMP_VAR(:)=IMP_VAR(:)
          DO I=1,ATM(ID)%NXzV
	    J=ATM(ID)%EQXzV+(I-1)
            SE(ID)%LNK_TO_IV(J)=I
            LOC_IMP_VAR(J)=.FALSE.
          END DO
          NIV=ATM(ID)%NXzV
!
! Include all possible ionization states (includeing those reached by
! X-ray ionization, and charge exchange reactions.
!
          DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	    J=ATM(ID)%EQXzV+ATM(ID)%NXzV+(SE(ID)%EQ_TO_ION_LEV_PNT(I)-1)
            SE(ID)%LNK_TO_IV(J)=I
            LOC_IMP_VAR(J)=.FALSE.
          END DO
	  NIV=SE(ID)%N_SE-1
!
! Now set those variables which are important for ALL species, and which have NOT
! already been included.
!
          DO I=1,NT
            IF(LOC_IMP_VAR(I))THEN
              NIV=NIV+1
              SE(ID)%LNK_TO_IV(I)=NIV
            END IF
          END DO
!
          SE(ID)%N_IV=NIV
!
! Create reverse links.
!
          ALLOCATE (SE(ID)%LNK_TO_F(NIV))
          SE(ID)%LNK_TO_F(:)=0		!1:NIV
	  DO I=1,NT
	    J=SE(ID)%LNK_TO_IV(I)
	    IF(J .NE. 0)SE(ID)%LNK_TO_F(J)=I
	  END DO
!
	  WRITE(93,*)SE(ID)%LNK_TO_F
	  WRITE(93,*)' '
	  WRITE(93,*)SE(ID)%LNK_TO_IV
!
        ELSE
          SE(ID)%XzV_PRES=.FALSE.
        END IF
      END DO
!
      RETURN
      END    
