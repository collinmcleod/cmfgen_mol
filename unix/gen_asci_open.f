C                                               
C Vax routine to open a sequential and formatted file for output.
C Allows
C
C      (i)  A file to be appended to
C      (ii) A file to be opened for input only.
C
C All these asci files are opened as shared files, and with
C CARRIAGECONTROL='LIST'. These two commands are VAX dependent.
C
C The SHARED option enables the file to be copied, typed etc while opened. 
C Useful for examining OUTGEN whilst program is running.
C
	SUBROUTINE GEN_ASCI_OPEN(LU,FILE_NAME,
	1            FILE_STATUS,FILE_POSIT,FILE_ACTION,FILE_RECL,IOS)
	INTEGER*4 LU,FILE_RECL,IOS
	CHARACTER*(*) FILE_NAME,FILE_STATUS,FILE_POSIT,FILE_ACTION
	CHARACTER*20  LOC_NAME,LOC_STATUS,LOC_POSIT,LOC_ACTION
C
	EXTERNAL ERROR_LU
	INTEGER*4 ERROR_LU,LUER,IZERO,LOC_RECL
C
	LUER=ERROR_LU()
	IZERO=0
C
	LOC_POSIT=FILE_POSIT
	IF(LOC_POSIT .EQ. ' ')LOC_POSIT='REWIND'
	CALL SET_CASE_UP(LOC_POSIT,IZERO,IZERO)
	IF(LOC_POSIT .NE. 'REWIND' .AND. 
	1     LOC_POSIT .NE. 'APPEND' .AND. 
	1     LOC_POSIT .NE. 'ASIS')THEN
	  WRITE(LUER,*)'Error in GEN_ASCI_OPEN - invald FILE_POSIT'
	  WRITE(LUER,*)'File name is',FILE_NAME
	  WRITE(LUER,*)'FILE_POSIT is',FILE_POSIT
	  STOP
	END IF
C
	LOC_STATUS=FILE_STATUS
	IF(LOC_STATUS .EQ. ' ')LOC_STATUS='UNKNOWN'
	CALL SET_CASE_UP(LOC_STATUS,IZERO,IZERO)
	IF(LOC_STATUS .NE. 'UNKNOWN' .AND. 
	1     LOC_STATUS .NE. 'REPLACE' .AND. 
	1     LOC_STATUS .NE. 'OLD' .AND.
	1     LOC_STATUS .NE. 'NEW')THEN
	  WRITE(LUER,*)'Error in GEN_ASCI_OPEN - invald FILE_STATUS'
	  WRITE(LUER,*)'File name is',FILE_NAME
	  WRITE(LUER,*)'FILE_STATUS is',FILE_STATUS
	  STOP
	END IF
C
	LOC_ACTION=FILE_ACTION
	IF(LOC_ACTION .EQ. ' ')LOC_ACTION='READWRITE'
	CALL SET_CASE_UP(LOC_ACTION,IZERO,IZERO)
	IF(LOC_ACTION .NE. 'READWRITE' .AND. 
	1     LOC_ACTION .NE. 'READ' .AND. 
	1     LOC_ACTION .NE. 'WRITE')THEN
	  WRITE(LUER,*)'Error in GEN_ASCI_OPEN - invald FILE_ACTION'
	  WRITE(LUER,*)'File name is',FILE_NAME
	  WRITE(LUER,*)'FILE_ACTION is',FILE_ACTION
	  STOP
	END IF
C
C Only the lase line of open statement is VAX dependent.
C Should be deleted for the CRAY.
C
	LOC_RECL=FILE_RECL
	IF(FILE_RECL .EQ. 0)LOC_RECL=132
        OPEN(UNIT=LU,FILE=FILE_NAME,ACTION=LOC_ACTION,STATUS=LOC_STATUS,
	1      POSITION=LOC_POSIT,RECL=LOC_RECL,IOSTAT=IOS)
C	1      SHARED,CARRIAGECONTROL='LIST')
C
	RETURN
	END



