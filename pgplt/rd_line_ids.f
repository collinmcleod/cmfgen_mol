!
! General purpose line plotting routiine.
!
	SUBROUTINE RD_LINE_IDS(ID_FILNAME,AIR_WAVELENGTHS,XPAR)
	USE MOD_COLOR_PEN_DEF
	USE LINE_ID_MOD
	IMPLICIT NONE
!
	REAL*4 XPAR(2)
	LOGICAL AIR_WAVELENGTHS
	CHARACTER(LEN=*) ID_FILNAME
!
	REAL*8 DP_T1
	REAL*8, EXTERNAL :: LAM_AIR
	INTEGER LU_ID
	INTEGER J,K
	INTEGER IOS
	INTEGER IOFF
        CHARACTER(LEN=200) TMP_STR
!
	CALL GET_LU(LU_ID,'In RD_ID_LINES')
	N_LINE_IDS=0
        OPEN(UNIT=LU_ID,FILE=TRIM(ID_FILNAME),STATUS='OLD',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  TMP_STR='!'
	  DO WHILE(TMP_STR(1:1) .EQ. '!')
	    READ(LU_ID,'(A)')TMP_STR
	  END DO
!
! Check which format
!
	  IF(TMP_STR(1:7) .NE. 'Species')BACKSPACE(LU_ID)
	  J=0
	  DO WHILE(J+1 .LE. 5000)
	    READ(LU_ID,'(A)',END=2000,IOSTAT=IOS)TMP_STR
	    READ(TMP_STR,*)LINE_ID(J+1),ID_WAVE(J+1),TAU(J+1),ID_WAVE_OFF(J+1),
	1               ID_Y_BEG(J+1),ID_Y_END(J+1)
!
	    TMP_STR=TMP_STR(20:)
	    K=INDEX(TMP_STR,TRIM(LINE_ID(J+1)))
	    FULL_LINE_ID(J+1)=TMP_STR(K:)
	    IF(ID_WAVE(J+1) .LT. 0.0D0)THEN
	       ID_WAVE(J+1)=ABS(ID_WAVE(J+1))
	       OBSERVED_WAVE(J+1)=.FALSE.
	    ELSE
	       OBSERVED_WAVE(J+1)=.TRUE.
	    END IF
	    IF(AIR_WAVELENGTHS)THEN
	      DP_T1=ID_WAVE(J+1)
	      ID_WAVE(J+1)=LAM_AIR(DP_T1)
	    END IF
	    IF( (ID_WAVE(J+1)-XPAR(1))*(XPAR(2)-ID_WAVE(J+1)) .GT. 0 .AND.
	1                   ABS(TAU(J+1)) .GT. TAU_CUT)THEN
	      J=J+1
	      N_LINE_IDS=J
	      IF(LINE_ID(J)(2:2) .EQ. 'k')LINE_ID(J)(2:2)='i'
	      IF(LINE_ID(J)(2:2) .EQ. '2')LINE_ID(J)(2:)='II'//LINE_ID(J)(3:)
	      IF(LINE_ID(J)(3:3) .EQ. '2')LINE_ID(J)(3:)='II'//LINE_ID(J)(4:)
	      IF(LINE_ID(J)(2:5) .EQ. 'XSIX')LINE_ID(J)(2:)='XVI'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(2:5) .EQ. 'XSEV')LINE_ID(J)(2:)='XVII'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(3:6) .EQ. 'XSIX')LINE_ID(J)(3:)='XVI'//LINE_ID(J)(7:)
	      IF(LINE_ID(J)(3:6) .EQ. 'XSEV')LINE_ID(J)(3:)='XVII'//LINE_ID(J)(7:)
	      IF(LINE_ID(J)(2:4) .EQ. 'SIX')LINE_ID(J)(2:)='VI'//LINE_ID(J)(5:)
	      IF(LINE_ID(J)(2:4) .EQ. 'SEV')LINE_ID(J)(2:)='VII'//LINE_ID(J)(5:)
	      IF(LINE_ID(J)(3:5) .EQ. 'SIX')LINE_ID(J)(3:)='VI'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(3:5) .EQ. 'SEV')LINE_ID(J)(3:)='VII'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(2:2) .GE. 'a' .AND. LINE_ID(J)(2:2) .LE.  'z')THEN
	        LINE_ID(J)(3:)=' '//LINE_ID(J)(3:)
	      ELSE
	        LINE_ID(J)(2:)=' '//LINE_ID(J)(2:)
	      END IF
	    END IF
	  END DO
2000	  CONTINUE
	  CLOSE(UNIT=LU_ID)
	  WRITE(6,*)'Number of lines read in is',N_LINE_IDS
	ELSE
	  WRITE(6,*)'Unable top open file'
	END IF
!
	RETURN
	END
