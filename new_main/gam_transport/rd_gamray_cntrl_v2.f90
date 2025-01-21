! This subroutine is to have a quick thing that reads in the input parameters
! so that they aren't hardwired into the code
!
! Created Oct. 12th, 2015
!
! Edited 25 Sept, 2017 : Changed USE MOD_GAMMA_V2 to USE
! 			 MOD_GAMRAY_CNTRL_VARIABLES
!
	SUBROUTINE RD_GAMRAY_CNTRL_V2(FILENAME)
	USE SET_KIND_MODULE
	USE MOD_RD_GAMRAY_CNTRL_VARIABLES
	IMPLICIT NONE
!
	INTEGER :: I,J,K,L
	INTEGER :: IOS
	INTEGER :: LUER
	INTEGER :: ERROR_LU
	EXTERNAL :: ERROR_LU
!
	CHARACTER(LEN=30) :: FILENAME
	CHARACTER(LEN=140) :: STRING
	CHARACTER(LEN=40) :: STRING1,STRING2
!
	LUER=ERROR_LU()
!
	OPEN(UNIT=7,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,'(A,A)')'Error reading file ',FILENAME
	  WRITE(LUER,*)'IOSTAT:',IOS
	  STOP
	END IF
	CHEB_ORDER=8
	ANG_MULT=2
	GAM_IT=1
	GRAY_LINE_MIN=50
	GRAY_TAIL_MIN1=2000
	GRAY_TAIL_MIN2=1000
	GRAY_TAIL_MIN3=500
	GRAY_TAIL_MIN4=100
	GRAY_END_MIN=200
	NU_GRID_MAX=20000	
	V_GAUSS=100.0_LDP		! km/s
	V_INF=1.0E+03_LDP		! km/s
	V_TAIL=5.0E+02_LDP		! km/s
	V_LINE=20.0_LDP		! km/s
	V_BETWEEN=2.5E+03_LDP	! km/s
	V_END=2.5E+03_LDP		! km/s
	BLUE_GAUSS=15.0_LDP
	RED_GAUSS=15.0_LDP
	COMP_TAIL_SPLIT_FRAC=0.05_LDP
	NORM_GAM_LINES=.FALSE.
	NORM_LINES_ONLY=.TRUE.
	GRID_MIN_PTS=.FALSE.
	GRID_TST=.FALSE.
	GAM_IB_METH='ZERO_FLUX'
	OBS_REL_FULL=.TRUE.
	EGAM_MAX=3.6_LDP		! MeV
	EGAM_MIN=1.0E-03_LDP	! MeV
	DO_NORM_ETA=.FALSE.
	VERBOSE_GAMMA=.FALSE.
	WRITE_J_DATA=.FALSE.
!
	OPEN(UNIT=8,FILE='GAMMA_MODEL',STATUS='UNKNOWN',ACTION='WRITE')
	STRING=''
	DO WHILE(1 .EQ. 1)
	  STRING1=''
	  STRING2=''
200	  READ(7,'(A)',END=100)STRING
	  K=INDEX(STRING,'[')
	  L=INDEX(STRING,'!')
	  IF(K .EQ. 0 .OR. L .EQ. 1)GOTO 200
	  STRING1=ADJUSTL(STRING(1:K-1))
	  STRING1=TRIM(STRING1)
	  J=INDEX(STRING,']')
	  STRING2=STRING(K+1:J-1)
	  STRING2=ADJUSTL(STRING2)
	  STRING2=TRIM(STRING2)
	  IF(TRIM(STRING2) .EQ. 'CHEB_ORDER')THEN
	    READ(STRING1,*)CHEB_ORDER
	    WRITE(8,'(1X,I2,T25,1X,A)')CHEB_ORDER,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'NU_GRID_MAX')THEN
	    READ(STRING1,*)NU_GRID_MAX
	    WRITE(8,'(1X,I5,T25,1X,A)')NU_GRID_MAX,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_LINE_MIN')THEN
	    READ(STRING1,*)GRAY_LINE_MIN
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_LINE_MIN,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_TAIL_MIN1')THEN
	    READ(STRING1,*)GRAY_TAIL_MIN1
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_TAIL_MIN1,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_TAIL_MIN2')THEN
	    READ(STRING1,*)GRAY_TAIL_MIN2
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_TAIL_MIN2,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_TAIL_MIN3')THEN
	    READ(STRING1,*)GRAY_TAIL_MIN3
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_TAIL_MIN3,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_TAIL_MIN4')THEN
	    READ(STRING1,*)GRAY_TAIL_MIN4
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_TAIL_MIN4,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GRAY_END_MIN')THEN
	    READ(STRING1,*)GRAY_END_MIN
	    WRITE(8,'(1X,I4,T25,1X,A)')GRAY_END_MIN,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_GAUSS')THEN
	    READ(STRING1,*)V_GAUSS
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_GAUSS,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_INF')THEN
	    READ(STRING1,*)V_INF
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_INF,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_LINE')THEN
	    READ(STRING1,*)V_LINE
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_LINE,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_TAIL')THEN
	    READ(STRING1,*)V_TAIL
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_TAIL,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_BETWEEN')THEN
	    READ(STRING1,*)V_BETWEEN
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_BETWEEN,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'V_END')THEN
	    READ(STRING1,*)V_END
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')V_END,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'BLUE_GAUSS')THEN
	    READ(STRING1,*)BLUE_GAUSS
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')BLUE_GAUSS,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'RED_GAUSS')THEN
	    READ(STRING1,*)RED_GAUSS
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')RED_GAUSS,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'NU_MIN')THEN
	    READ(STRING1,*)EGAM_MIN
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')EGAM_MIN,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'NU_MAX')THEN
	    READ(STRING1,*)EGAM_MAX
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')EGAM_MAX,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'ANG_MULT')THEN
	    READ(STRING1,*)ANG_MULT
	    WRITE(8,'(1X,I2,T25,1X,A)')ANG_MULT,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'GAM_IT')THEN
	    READ(STRING1,*)GAM_IT
	    WRITE(8,'(1X,I2,T25,1X,A)')GAM_IT,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'NORM_GAM')THEN
	    READ(STRING1,*)NORM_GAM_LINES
	    WRITE(8,'(1X,L1,T25,1X,A)')NORM_GAM_LINES,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'NORM_LINES_ONLY')THEN
	    READ(STRING1,*)NORM_LINES_ONLY
	    WRITE(8,'(1X,L1,T25,1X,A)')NORM_LINES_ONLY,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'TAIL_SPLIT_FRAC')THEN
	    READ(STRING1,*)COMP_TAIL_SPLIT_FRAC
	    WRITE(8,'(1X,ES16.6,T25,1X,A)')COMP_TAIL_SPLIT_FRAC,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'USE_MIN_PTS_GRAY_LINES')THEN
	    READ(STRING1,*)GRID_MIN_PTS
	    WRITE(8,'(1X,L1,T25,1X,A)')GRID_MIN_PTS,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'TEST_NU_GRID')THEN
	    READ(STRING1,*)GRID_TST
	    WRITE(8,'(1X,L1,T25,1X,A)')GRID_TST,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'IB_METH')THEN
	    READ(STRING1,*)GAM_IB_METH
	    WRITE(8,'(1X,A10,T25,1X,A)')GAM_IB_METH,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'FULL_REL_OBS')THEN
	    READ(STRING1,*)OBS_REL_FULL
	    WRITE(8,'(1X,L1,T25,1X,A)')OBS_REL_FULL,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'DO_NORM_ETA')THEN
	    READ(STRING1,*)DO_NORM_ETA
	    WRITE(8,'(1X,L1,T25,1X,A)')DO_NORM_ETA,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'VERBOSE_GAM')THEN
	    READ(STRING1,*)VERBOSE_GAMMA
	    WRITE(8,'(1X,L1,T25,1X,A)')VERBOSE_GAMMA,TRIM(STRING2)
	  ELSE IF(TRIM(STRING2) .EQ. 'WRITE_J')THEN
	    READ(STRING1,*)WRITE_J_DATA
	    WRITE(8,'(1X,L1,T25,1X,A)')WRITE_J_DATA,TRIM(STRING2)
	  ELSE
	    WRITE(8,*)'Error: Unknown option in ',FILENAME
	    STOP
	  END IF
	END DO
100 CONTINUE
	CLOSE(UNIT=7)
	CLOSE(UNIT=8)
	END SUBROUTINE
