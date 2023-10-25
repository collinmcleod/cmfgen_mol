!
! Designed to be called by CMFGEN, and outputs:
!            ION_ID, Level name, SL, Level.
!
	SUBROUTINE WR_LEVEL_LINKS
	USE MOD_CMFGEN
	INTEGER LU
!
! Altered: 24-Oct-2023 -- Cleaning; output file name changed.
! Created: 15-Oct-2023
!
	CALL GET_LU(LU,'In WR_LEVEL_LINKS')
	OPEN(UNIT=LU,FILE='LEVEL_SL_STEQ_LINKS',STATUS='UNKNOWN')
	WRITE(LU,'(4X,A,36X,A,4X,A,7X,A,2X,A)')'Ion(ID)','Name','Level','SL','STEQ(lev)'
	DO ID=1,NUM_IONS
           IF(ATM(ID)%XzV_PRES)THEN
             DO J=1,ATM(ID)%NXzV_F
	       NL= ATM(ID)%F_TO_S_XzV(MNL)+ ATM(ID)%EQXzV-1
               WRITE(LU,'(1X,A10,A40,3X,I6,3X,I6,5X,I6)')ION_ID(ID),ATM(ID)%XzVLEVNAME_F(J),
	1                   J,ATM(ID)%F_TO_S_XzV(J),NL
             END DO
           END IF
        END DO
	CLOSE(UNIT=LU)
!
	RETURN
	END
