	SUBROUTINE RD_MONTE_LINE(FREQ,FLUX,NMAX,NPTS,FILENAME,IOS)
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Routine to plot the polarization across a line computed using the Monte-Carlo
! code.
!
	INTEGER*4 NMAX
	INTEGER NPTS
	INTEGER IOS
!
	REAL(KIND=LDP) FREQ(NMAX)
	REAL(KIND=LDP) FLUX(NMAX)
	CHARACTER(LEN=80) FILENAME
!
	CHARACTER*80 STRING
	REAL(KIND=LDP) V(NMAX),INTEN(NMAX),P(NMAX),THETA(NMAX),Q(NMAX),U(NMAX)
	REAL(KIND=LDP) EI(NMAX),EP(NMAX),ET(NMAX),EQ(NMAX),EU(NMAX)
!
	REAL*4 QI(NMAX)
!
	INTEGER*4 I,J,K,N
	REAL(KIND=LDP) C_KMS,SUM_I,SUM_Q
	REAL(KIND=LDP) LINE_WAVE
	REAL(KIND=LDP) SCALING_FLUX
!
 	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
	CHARACTER(LEN=20) XLAB
!
	LINE_WAVE=0.0_LDP
	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error opening ',TRIM(FILENAME)
	  WRITE(6,*)'IOS=',IOS
	  RETURN
	END IF
!
	STRING=' '
	DO WHILE(1. EQ. 1)
	  READ(10,'(A)')STRING
          IF(INDEX(STRING,'Number of bins:') .NE. 0)THEN
	    EXIT
	  ELSE IF(INDEX(STRING,'Scaling flux') .NE. 0)THEN
	    I=INDEX(STRING,'   ')
	    READ(STRING(I:),*)SCALING_FLUX
	  ELSE IF(INDEX(STRING,'Wavelength') .NE. 0)THEN
	    I=INDEX(STRING,'   ')
	    READ(STRING(I:),*)LINE_WAVE
	  END IF
	END DO
	READ(STRING(20:),*)NPTS
    	IF(NPTS .GT. NMAX)THEN
	  WRITE(6,*)'Error - NMAX in RD_MONTE_LINE is too small'
	  STOP
	END IF
!
	DO I=1,NPTS
	  READ(10,*)V(I),INTEN(I),P(I),THETA(I),Q(I),U(I)
	  READ(10,*)EI(I),EP(I),ET(I),EQ(I),EU(I)
	END DO
!
! Get correct wavelngth of transition.
!
	IF(LINE_WAVE .EQ. 0.0_LDP)THEN
	  CALL GEN_IN(LINE_WAVE,'Lamba (A): ')
	END IF
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
! Freq is returned in units of 10^15 Hz.
!
	DO I=1,NPTS
	  FLUX(I)=INTEN(I)*SCALING_FLUX
	  FREQ(I)=0.01_LDP*C_KMS/(LINE_WAVE*(1.0_LDP+V(I)/C_KMS))
	END DO
	CLOSE(UNIT=10)
!
	RETURN
	END
