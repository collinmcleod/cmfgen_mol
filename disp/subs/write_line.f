C
C Altered 11-Mar-1991: --- Conversion to microns was in wrong direction.
C
	SUBROUTINE WRITE_LINE(LEV1,LEV2,FREQ,DESC)
	IMPLICIT NONE
	INTEGER LEV1,LEV2
	REAL*8 FREQ,LAMVACAIR,WAVE
	CHARACTER*(*) DESC
C
	WAVE=LAMVACAIR(FREQ)
C
	IF(WAVE .LT. 9999.99)THEN
	  WRITE(6,100)TRIM(DESC),LEV2,LEV1,WAVE
100	  FORMAT(1X,A,'(',I3,'-',I3,') at ', F8.3,' Angstroms')
	ELSE IF(WAVE .LT. 9.99999E+07)then
	  WAVE=WAVE*1.0D-04
	  WRITE(6,200)TRIM(DESC),LEV2,LEV1,WAVE
200	  FORMAT(1X,A,'(',I3,'-',I3,') at ', F8.3,' microns')
	ELSE
	  WAVE=WAVE*1.0D-04
	  WRITE(6,300)TRIM(DESC),LEV2,LEV1,WAVE
300	  FORMAT(1X,A,'(',I3,'-',I3,') at ', 1PE11.4,' microns')
	END IF
C
	RETURN
	END
