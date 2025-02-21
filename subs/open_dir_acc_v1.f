!
! Simple routine designed to open a direct acces, unformatted, file.
! The file cannot exist --- if it does a _NEW file is created.
! The subroutine also creates the informational file with data on
! the record length.
!
	SUBROUTINE OPEN_DIR_ACC_V1(NV,IRECL,DATE,FILENAME,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 08-MAr-2006
!
	INTEGER NV
	INTEGER IRECL
	INTEGER LU
	CHARACTER(LEN=*) FILENAME
	CHARACTER(LEN=*) DATE
!
	CHARACTER(LEN=132) LOC_FILENAME
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL FILE_EXISTS
!
	LUER=ERROR_LU()
	INQUIRE(FILE=TRIM(FILENAME),EXIST=FILE_EXISTS)
	IF(FILE_EXISTS)THEN
	  LOC_FILENAME=TRIM(FILENAME)//'_NEW'
	  CALL WRITE_DIRECT_INFO_V3(NV,IRECL,DATE,LOC_FILENAME,LU)
	  OPEN(UNIT=LU,FILE=TRIM(LOC_FILENAME),FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=IRECL,IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    WRITE(LUER,*)'Error opening ',TRIM(FILENAME),' as file already exists'
	    WRITE(LUER,*)'Opened ',TRIM(LOC_FILENAME),' instead'
	  ELSE
	    WRITE(LUER,*)'Error opening ',TRIM(FILENAME),'or ',TRIM(LOC_FILENAME)
	    WRITE(LUER,*)'These files may need to be renamed/deleted'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	ELSE
	  CALL WRITE_DIRECT_INFO_V3(NV,IRECL,DATE,FILENAME,LU)
	  OPEN(UNIT=LU,FILE=TRIM(FILENAME),FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=IRECL,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening ',TRIM(FILENAME)
	    WRITE(LUER,*)'These files may need to be renamed/deleted'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	END IF
!
	RETURN
	END
