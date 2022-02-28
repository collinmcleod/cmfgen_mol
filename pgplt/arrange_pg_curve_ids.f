	SUBROUTINE ARRANGE_PG_CURVE_IDS(PEN_COLOR,PEN_OFFSET,NPEN,ADD_COMMA,RESET_TITLES)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
	INTEGER NPEN,PEN_OFFSET
	INTEGER PEN_COLOR(0:NPEN)
	LOGICAL RESET_TITLES
	LOGICAL ADD_COMMA
!
	INTEGER IP,IT,J,K,IPEN
	CHARACTER(LEN=5) CPEN
!
	IF(RESET_TITLES)THEN
	  N_TIT_SET=0
	  TITLE(N_HEADER+1:)=' '
	  CURVE_TITLE=' '
	ELSE
	  DO J=N_HEADER+1,N_TITLE
	    TITLE(J)=CURVE_TITLE(J-N_HEADER)
	  END DO
	  RETURN
	END IF
	N_TIT_SET=N_HEADER
!
! Now arrange the curve titles. We can force a line to end by ending it
! with \\ (as in latex). Maximum title length is 50.
!
	DO IP=1,NPLTS
	  IF(CD(IP)%CURVE_ID .EQ. ' ')THEN
	  ELSE
	    IF(N_TIT_SET .EQ. N_HEADER)THEN
	      TITLE(N_HEADER+1)=' '
	      N_TIT_SET=N_HEADER+1
	    END IF
	    J=INDEX(TITLE(N_TIT_SET),'\\')
	    IF(J .NE. 0)THEN
	      TITLE(N_TIT_SET)(J:)=TITLE(N_TIT_SET)(J+2:)
	      N_TIT_SET=N_TIT_SET+1
	    END IF
	    WRITE(CPEN,'(I2)')PEN_COLOR(IP+PEN_OFFSET)
	    CPEN=' \p'//ADJUSTL(CPEN)
	    IT=N_TIT_SET
	    K=LEN_TRIM(TITLE(IT))
	    IF(LEN_TRIM(TITLE(IT))+LEN_TRIM(CD(IP)%CURVE_ID) .LT. 47)THEN
	      IF(K .EQ. 0)THEN
	        TITLE(IT)=TRIM(CD(IP)%CURVE_ID)//TRIM(CPEN)
	      ELSE IF(ADD_COMMA)THEN
	        TITLE(IT)=TRIM(TITLE(IT))//', '//TRIM(CD(IP)%CURVE_ID)//TRIM(CPEN)
	      ELSE
	        TITLE(IT)=TRIM(TITLE(IT))//' '//TRIM(CD(IP)%CURVE_ID)//TRIM(CPEN)
	      END IF
	    ELSE
	      IT=IT+1; N_TIT_SET=IT
	      TITLE(IT)=TRIM(CD(IP)%CURVE_ID)//TRIM(CPEN)
	    END IF
	    WRITE(6,'(A)')TRIM(TITLE(IT))
	  END IF
	END DO
	CURVE_TITLE(1:N_TIT_SET-N_HEADER)=TITLE(N_HEADER+1:N_TIT_SET)
!
	RETURN
	END
