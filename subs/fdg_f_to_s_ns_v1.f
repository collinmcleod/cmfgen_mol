!
! Routine to read in the vector describing the matching of actual atomic levels
! to that in the model atom with super levels.
!
	SUBROUTINE FDG_F_TO_S_NS(NF,NS,NV,SL_OPTION,LUIN,FILENAME)
	IMPLICIT NONE
!
	INTEGER NF
	INTEGER NS
	INTEGER NV
	INTEGER LUIN
	CHARACTER(LEN=*) FILENAME
	CHARACTER(LEN=*) SL_OPTION
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER F_TO_S(NF)
	INTEGER INT_SEQ(NF)
!
	INTEGER LUER,I,J,K,IOS,ENTRY_NUM
	INTEGER CNT(NS)
	CHARACTER*132 STRING
!
	IF(NF .EQ. NS .OR. SL_OPTION .EQ. ' ')RETURN
!
	LUER=ERROR_LU()
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS,
	1       ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening ',FILENAME,' in RD_F_TO_S_IDS'
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'Number of energy levels') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
	  READ(LUIN,*)I
	  IF(I .LT. NF)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Insufficient levels in file'
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'Entry number of link') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
	  READ(LUIN,*)ENTRY_NUM
	  IF(ENTRY_NUM .LT. 2 .OR. ENTRY_NUM .GT. 20)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Bad entry number for level link'
	    STOP
	  END IF
!
! NB: All entries must be separated by at LEAST 2 spaces.
!
	  READ(LUIN,'(A)')STRING		!Blankline
	  DO I=1,NF
	    F_TO_S(I)=0.
	    READ(LUIN,'(A)')STRING
	    STRING=ADJUSTL(STRING)
!
! NB: ENTRY_NUM-2 as
!              -1 due to 1st entry being level ID
!              -1 as only need to point to entry.
!
	    J=INDEX(STRING,'  ')
	    DO K=1,ENTRY_NUM-2
	      STRING(1:)=STRING(J:)
	      STRING=ADJUSTL(STRING)
	      J=INDEX(STRING,'  ')
	    END DO
	    READ(STRING(J:),*)F_TO_S(I),INT_SEQ(I)
	  END DO
	CLOSE(LUIN)
!
! Now check that:
! (1) All super-levels have at least 1 level from full atom
! (2) Correct number of super levels.
!
	DO I=1,NF
	  IF(F_TO_S(I) .LE. 0 .OR. F_TO_S(I) .GT. NS)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Incorrect super level ID for level ',I
	    STOP
	  END IF
	END DO
!
! Check that all super levels have at least 1 level.
!
	DO I=1,NS
	  CNT(I)=0
	END DO
	DO J=1,NF
	  I=F_TO_S(J)
	  CNT(I)=CNT(I)+1
	END DO
	DO I=1,NS
	  IF(CNT(I) .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Super level ',I,' has no components'
	    STOP
	  END IF
	END DO
!
	CALL DO_SL_ADJUSTEMENT(F_TO_S,INT_SEQ,NF,NS,NV,SL_OPTION)
!
	RETURN
	END
!
	SUBROUTINE DO_SL_ADJUSTEMENT(F_TO_S,INT_SEQ,NF,NS,NV,SL_OPTION)
!
	INTEGER NF
	INTEGER NS
	INTEGER NV
	INTEGER F_TO_S(NF)
	INTEGER INT_SEQ(NF)
	CHARACTER(LEN=*) SL_OPTION
!
	INTEGER OLD_F_TO_S(NF)
	INTEGER TMP_F_TO_S
!
	INTEGER I, J
	INTEGER ID
	INTEGER COUNT
	INTEGER NS_LOW
!
!***************************************************************************
!***************************************************************************
!
! We now make the correction.
!
	IF(SL_OPTION(1:10) .EQ. ' ')THEN
!
! Nothing to do.
!
	ELSE IF(SL_OPTION(1:10) .EQ. 'SPLIT_LOW_')THEN
	  READ(SL_OPTION(11:),*)NS_LOW
	  COUNT=NF
	  OLD_F_TO_S(1:NF)=F_TO_S(1:NF)
	  DO I=1,NF
	    IF(F_TO_S(I) .LE. NS_LOW)THEN
	      COUNT=COUNT+1
	      F_TO_S(I)=F_TO_S(I)+COUNT
	    END IF
	  END DO
!
	  ID=0
	  F_TO_S(1:NF)=-F_TO_S(1:NF)
	  DO I=1,NF
	    IF(F_TO_S(I) .LT. 0)THEN
	      ID=ID+1
	      TMP_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	      DO J=I+1,NF
	        IF(F_TO_S(J) .EQ. TMP_F_TO_S)F_TO_S(J)=ID
	      END DO
	    ELSE IF(F_TO_S(I) .EQ. 0)THEN
	      ID=ID+1
	      TMP_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	    END IF
	  END DO
!
! Fix sequence number.
!
	  DO I=1,NF
	    IF(INT_SEQ(I) .NE. 0)THEN
	      DO J=1,NF
	        IF(INT_SEQ(I) .EQ. OLD_F_TO_S(J))THEN
	          INT_SEQ(I)=F_TO_S(J)
	          EXIT
	        END IF
	      END DO
	    END IF
	  END DO
!
	  J=0
	  DO I=1,NF
	    J=MAX(J,F_TO_S(I))
	  END DO
	  IF(NV .NE. 0)NV=NV+(NS-J)
	  NS=J
	ELSE
	  WRITE(6,*)'Error in DO_SL_ADJUSTEMENT: SL_OPTION not recognized'
	  WRITE(6,*)'SL_OPTION=',TRIM(SL_OPTION)
	  STOP
	END IF
!
	WRITE(101,*)'Old NS',MAXVAL(OLD_F_TO_S)
	WRITE(101,*)'Revised NS',NS
	DO I=1,NF
	  WRITE(101,'(4I6)')I,F_TO_S(I),OLD_F_TO_S(I),INT_SEQ(I)
	END DO
!
	RETURN
	END
