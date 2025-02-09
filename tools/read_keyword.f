!
! Set of subroutines that allows indiviual KEY WORD values in a
! CMFGEN input files to be read. The file is first read into
! STORE, and values are read from the store. Passed parameters
! indicate when the STORE should be updated or deleted.
!
! If the keyword is not present, and MUST_BE_PRES=.FALSE, the passed
! value is returned.
!
! Created 20-Nov-2014:  Basesed on UPDATE_KEYWORD
!
	MODULE READ_KEYWORD_INTERFACE
	INTEGER, PARAMETER :: MAX_RECS=800
	INTEGER, SAVE :: NUM_RECS=0
	CHARACTER(LEN=80), SAVE, ALLOCATABLE :: STORE(:)
	CHARACTER(LEN=20) TMP_STR
	INTERFACE READ_KEYWORD
          MODULE PROCEDURE READ_KEYWORD_DP,
	1                  READ_KEYWORD_LOG,
	1                  READ_KEYWORD_INT,
	1                  READ_KEYWORD_STR
        END INTERFACE
        CONTAINS
!
!
!
	SUBROUTINE READ_KEYWORD_DP(VALUE,KEY_WORD,MUST_BE_PRES,DATA_FILE,RD_FILE,DEL_STORE,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP) VALUE
	CHARACTER(LEN=*) KEY_WORD
	CHARACTER(LEN=*) DATA_FILE
	LOGICAL RD_FILE
	LOGICAL DEL_STORE
	LOGICAL MUST_BE_PRES
	INTEGER LU,LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER I,J,K
!
	CALL RD_KEY_WRD_FILE(DATA_FILE,RD_FILE,LU)
!
	K=0
	DO I=1,NUM_RECS
	  K=INDEX(STORE(I),KEY_WORD)
	  IF(K .NE. 0)THEN
	    TMP_STR=ADJUSTL(STORE(I)(1:K-1))
	    READ(TMP_STR,*)VALUE
	    EXIT
	  END IF
	END DO
!
	IF(K .EQ. 0 .AND. MUST_BE_PRES)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in READ_KEYWORD'
	  WRITE(LUER,*)TRIM(KEY_WORD),' not found in file ',TRIM(DATA_FILE)
	  STOP
	END IF
!
	CALL DELETE_STORAGE(DEL_STORE)
!
	END SUBROUTINE READ_KEYWORD_DP
!
!
	SUBROUTINE READ_KEYWORD_LOG(VALUE,KEY_WORD,MUST_BE_PRES,DATA_FILE,RD_FILE,DEL_STORE,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	LOGICAL VALUE
	CHARACTER(LEN=*) KEY_WORD
	CHARACTER(LEN=*) DATA_FILE
	LOGICAL RD_FILE
	LOGICAL DEL_STORE
	LOGICAL MUST_BE_PRES
	INTEGER LU,LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER I,K
!
	CALL RD_KEY_WRD_FILE(DATA_FILE,RD_FILE,LU)
!
	K=0
	DO I=1,NUM_RECS
	  K=INDEX(STORE(I),KEY_WORD)
	  IF(K .NE. 0)THEN
	    TMP_STR=ADJUSTL(STORE(I)(1:K-1))
	    READ(TMP_STR,*)VALUE
	    EXIT
	  END IF
	END DO
!
	IF(K .EQ. 0 .AND. MUST_BE_PRES)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in READ_KEYWORD'
	  WRITE(LUER,*)TRIM(KEY_WORD),' not found in file ',TRIM(DATA_FILE)
	  STOP
	END IF
!
	CALL DELETE_STORAGE(DEL_STORE)
!
	END SUBROUTINE READ_KEYWORD_LOG
!
!
	SUBROUTINE READ_KEYWORD_INT(VALUE,KEY_WORD,MUST_BE_PRES,DATA_FILE,RD_FILE,DEL_STORE,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER VALUE
	CHARACTER(LEN=*) KEY_WORD
	CHARACTER(LEN=*) DATA_FILE
	LOGICAL RD_FILE
	LOGICAL DEL_STORE
	LOGICAL MUST_BE_PRES
	INTEGER LU,LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER I,J,K
!
	CALL RD_KEY_WRD_FILE(DATA_FILE,RD_FILE,LU)
!
	K=0
	DO I=1,NUM_RECS
	  K=INDEX(STORE(I),KEY_WORD)
	  IF(K .NE. 0)THEN
	    TMP_STR=ADJUSTL(TMP_STR(1:K-1))
	    READ(TMP_STR,*)VALUE
	    EXIT
	  END IF
	END DO
!
	IF(K .EQ. 0 .AND. MUST_BE_PRES)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in READ_KEYWORD'
	  WRITE(LUER,*)TRIM(KEY_WORD),' not found in file ',TRIM(DATA_FILE)
	  STOP
	END IF
!
	CALL DELETE_STORAGE(DEL_STORE)
!
	END SUBROUTINE READ_KEYWORD_INT
!
!
	SUBROUTINE READ_KEYWORD_STR(VALUE,KEY_WORD,MUST_BE_PRES,DATA_FILE,RD_FILE,DEL_STORE,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	CHARACTER(LEN=*) VALUE
	CHARACTER(LEN=*) KEY_WORD
	CHARACTER(LEN=*) DATA_FILE
	LOGICAL RD_FILE
	LOGICAL DEL_STORE
	INTEGER LU,LUER,ERROR_LU
	LOGICAL MUST_BE_PRES
	EXTERNAL ERROR_LU
	INTEGER I,J,K
!
	CALL RD_KEY_WRD_FILE(DATA_FILE,RD_FILE,LU)
!
	K=0
	DO I=1,NUM_RECS
	  K=INDEX(STORE(I),KEY_WORD)
	  IF(K .NE. 0)THEN
	    VALUE=ADJUSTL(STORE(I)(1:K-1))
	    EXIT
	  END IF
	END DO
!
	IF(K .EQ. 0 .AND. MUST_BE_PRES)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in READ_KEYWORD'
	  WRITE(LUER,*)TRIM(KEY_WORD),' not found in file ',TRIM(DATA_FILE)
	  STOP
	END IF
!
	CALL DELETE_STORAGE(DEL_STORE)
!
	END SUBROUTINE READ_KEYWORD_STR
!
!
!
	SUBROUTINE DELETE_STORAGE(DEL_STORE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	LOGICAL DEL_STORE
!
	IF(DEL_STORE)THEN
	  DEALLOCATE (STORE)
	  NUM_RECS=0
	END IF
!
	END SUBROUTINE DELETE_STORAGE
!
!
!
	SUBROUTINE RD_KEY_WRD_FILE(DATA_FILE,RD_FILE,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) DATA_FILE
	LOGICAL RD_FILE
	INTEGER LU,LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER I
!
	LUER=ERROR_LU()
	IF(RD_FILE)THEN
!
	  IF(ALLOCATED(STORE))THEN
	    WRITE(LUER,*)'Error -- storage location associated with READ_KEYWORD is still allocated'
	    WRITE(LUER,*)'A dump of the first 10 keywords in the current STORE follows'
	    DO I=1,MIN(10,NUM_RECS)
	      WRITE(LUER,'(A)')TRIM(STORE(I))
	    END DO
	    STOP
	  END IF
!
	  ALLOCATE (STORE(MAX_RECS))
	  OPEN(UNIT=LU,FILE=TRIM(DATA_FILE),ACTION='READ')
	    DO I=1,MAX_RECS
	       READ(LU,'(A)',END=100)STORE(I)
	       NUM_RECS=I
	    END DO
	    READ(LU,'(A)',END=100)TMP_STR
	    WRITE(LUER,*)'Insufficient storage in READ_KEYWORD'
	    WRITE(LUER,*)'Keyword data file is:',TRIM(DATA_FILE)
	    STOP
100	    CONTINUE
	  CLOSE(LU)
	END IF
!
	IF(NUM_RECS .EQ. 0)THEN
	  WRITE(LUER,*)'No data file with keywords has been opened in READ_KEYWORD'
	  STOP
	END IF
!
	RETURN
	END SUBROUTINE RD_KEY_WRD_FILE
!
	END MODULE READ_KEYWORD_INTERFACE
