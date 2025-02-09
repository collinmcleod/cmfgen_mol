	SUBROUTINE GET_TITLES(TITLE_FILE,TITLE,N_TITLE,IOS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 15-Aug-2020: Added IOS to call.
!
	INTEGER N_TITLE
	CHARACTER(LEN=*) TITLE_FILE
	CHARACTER(LEN=*) TITLE(N_TITLE)
!
	INTEGER I,IOS,LU
!
	CALL GET_LU(LU)
	OPEN(UNIT=LU,FILE=TRIM(TITLE_FILE),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  TITLE=' '
	  DO I=1,N_TITLE
	    READ(LU,'(A)',END=100)TITLE(I)
	    TITLE(I)=ADJUSTL(TITLE(I))
	  END DO
100	  CONTINUE
	ELSE
	  WRITE(6,*)'Unable to open file with titles: ',ADJUSTL(TRIM(TITLE_FILE))
	END IF
	CLOSE(LU)
!
	RETURN
	END
