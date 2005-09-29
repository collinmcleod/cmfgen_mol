	PROGRAM COUNT_PHOT_VALUES
	IMPLICIT NONE
!
	INTEGER COUNT_DATA
	INTEGER COUNT_LEVS
	INTEGER NC
	INTEGER TYPE
	INTEGER I
	REAL*8 FREQ,  OLD_FREQ
	REAL*8 CROSS, OLD_CROSS
!	
	CHARACTER*200 STRING
	CHARACTER*80 FILENAME
!
	WRITE(6,*)'Input PHOT file name'
	READ(5,*)FILENAME
!
	COUNT_LEVS=0
	COUNT_DATA=0
!
	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD',ACTION='READ')
	DO WHILE (1 .EQ. 1)
	  READ(10,FMT='(A)',END=100)STRING
	  IF(INDEX(STRING,'!Type of cross-section') .NE. 0)THEN
	    READ(STRING,*)TYPE
	  END IF
          IF(INDEX(STRING,'!Number of cross-section points') .NE. 0)THEN
	    COUNT_LEVS=COUNT_LEVS+1
	    READ(STRING,*)NC
	    COUNT_DATA=COUNT_DATA+NC
	    IF(TYPE .EQ. 20)THEN
	      READ(10,*)FREQ,CROSS
	      DO I=2,NC
	        OLD_FREQ=FREQ; OLD_CROSS=CROSS
	        READ(10,*)FREQ,CROSS
	        IF(FREQ .LE. OLD_FREQ)THEN
	           WRITE(6,*)'Error in cross-section: frequency grid is non-monotonic'
	           WRITE(6,*)'      TYPE=',TYPE
	           WRITE(6,*)'        NC=',NC
	           WRITE(6,*)' LOW_FREQ=',OLD_FREQ
	           WRITE(6,*)'HIGH_FREQ=',FREQ
	           STOP
	        END IF
	      END DO
	    END IF
	  END IF
	END DO
100	CONTINUE
!
	WRITE(6,*)'Numer of levels is ',COUNT_LEVS
	WRITE(6,*)'Number of data values is',COUNT_DATA
!
	STOP
	END
