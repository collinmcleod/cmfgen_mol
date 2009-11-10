	PROGRAM COUNT_PHOT_VALUES
	IMPLICIT NONE
!
	INTEGER COUNT_DATA
	INTEGER COUNT_LEVS
	INTEGER NC
	INTEGER TYPE
	INTEGER I,K
	REAL*8 FREQ,  OLD_FREQ
	REAL*8 CROSS, OLD_CROSS
!	
	CHARACTER*200 STRING
	CHARACTER*80 FILENAME
	CHARACTER*30 LEVEL_NAME
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
	  IF(INDEX(STRING,'!Configuration name') .NE. 0)THEN
	    K=INDEX(STRING,'  ')
	    LEVEL_NAME=STRING(1:K)
	    READ(10,FMT='(A)',END=100)STRING
	    IF(INDEX(STRING,'!Type of cross-section') .EQ. 0)THEN
	      WRITE(6,*)'Error- tTYPE reconrd does not follow Config record'
	      WRITE(6,*)LEVEL_NAME
	    END IF
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
	        IF(FREQ .LT. OLD_FREQ)THEN
	           WRITE(6,*)'Error in cross-section: frequency grid is non-monotonic'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE
	           WRITE(6,*)'        NC=',NC
	           WRITE(6,*)' LOW_FREQ=',OLD_FREQ
	           WRITE(6,*)'HIGH_FREQ=',FREQ
	           WRITE(6,*)'    CROSS=',CROSS
	           STOP
	        END IF
	        IF(FREQ .EQ. OLD_FREQ)THEN
	           WRITE(6,*)'Warning in cross-section: equal frequencis'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE,'        NC=',NC
	           WRITE(6,*)' LOW_FREQ=',OLD_FREQ,'HIGH_FREQ=',FREQ
	           WRITE(6,*)'    CROSS=',CROSS
	        END IF
	        IF(CROSS .LT. 0)THEN
	           WRITE(6,*)'Warning in cross-section: negatve cross-section'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE,'        NC=',NC
	           WRITE(6,*)'    CROSS=',CROSS
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
