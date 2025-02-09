!
! Routine to merge 6 color PGPLOTS onto a single
! page. The plots should be done with:
!      EXPAND_CHAR=1.2; EXPAND_TICK=1.2; ASR=0.89; Plot Size=10.0 cm
! The plots are located in 2 columns, 2 per row.
!
! Landscape format, with CPS as printer.
!
	PROGRAM N_MULT_MERGE
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered   07-Apr-2017 -  Allows for multiple files to be created.
!                          The same plot format is used for all files.
! Finalized 29-May-1997
!
	CHARACTER*132 STRING
	CHARACTER*80 FILE1
	CHARACTER*80 OUTF
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: LU_IN=11
	INTEGER, PARAMETER :: LU_OUT=30
	INTEGER, PARAMETER :: LU_TERM=6
!
	INTEGER IOS,I,J,K,IREC,M
	INTEGER NR,NC
	INTEGER NPLTS,NCOLS,NROWS
!
	CHARACTER*10 TMP_STR
	CHARACTER*4 HOR_OFFSET
	CHARACTER*4 VER_OFFSET
	LOGICAL TEST
!
	WRITE(LU_TERM,*)' '
	WRITE(LU_TERM,*)'Program to merge N plots in C columns'
	WRITE(LU_TERM,*)'1x3 plots (cr): ',
	1    'EXPAND_CHAR=2.5; EXPAND_TICK=2.5; ASR=0.45; Plot Size=17.0 cm'
	WRITE(LU_TERM,*)'2x3 plots (cr): ',
	1    'EXPAND_CHAR=1.2; EXPAND_TICK=1.2; ASR=0.89; Plot Size=9.5 cm'
	WRITE(LU_TERM,*)'2x5 plots (cr): ',
	1    'EXPAND_CHAR=2.5; EXPAND_TICK=2.5; ASR=0.45; Plot Size=10.0 cm'
	WRITE(LU_TERM,*)'3x5 plots (cr): ',
	1    'EXPAND_CHAR=1.5; EXPAND_TICK=1.5; ASR=0.7; Plot Size=6.0 cm'
	WRITE(LU_TERM,*)' '
!
	NCOLS=2;  CALL GEN_IN(NCOLS,'Number of columns on page')
	NROWS=3;  CALL GEN_IN(NROWS,'Number of plot rows')
	NPLTS=NCOLS*NROWS
!
	I=9000/NCOLS
	WRITE(HOR_OFFSET,'(I4)')I
	I=7900-1300*NROWS
	WRITE(VER_OFFSET,'(I4)')I
	IF(NCOLS .EQ. 1 .AND. NROWS .EQ. 1)THEN
	  HOR_OFFSET='0'
	  VER_OFFSET='4100'
	ELSE IF(NCOLS .EQ. 3 .AND. NROWS .EQ. 5)THEN
	  HOR_OFFSET='3200'
	  VER_OFFSET='2300'
	ELSE IF(NCOLS .EQ. 2 .AND. NROWS .EQ. 5)THEN
	  HOR_OFFSET='4800'
	  VER_OFFSET='2400'
	ELSE IF(NCOLS .EQ. 1 .AND. NROWS .EQ. 4)THEN
	  HOR_OFFSET='0'
	  VER_OFFSET='3000'
	END IF
	CALL GEN_IN(HOR_OFFSET,'Horizontal offset')
	CALL GEN_IN(VER_OFFSET,'Vertical offset')
!
	TEST=.FALSE.
	FILE1='pgplot_1.ps'
!
	DO WHILE(1 .EQ. 1)
!
	  IREC=0
	  IOS=0
	  OUTF=' '
	  CALL GEN_IN(OUTF,'Output file')
	  IF(OUTF .EQ. ' ')STOP
	  CALL GEN_ASCI_OPEN(LU_OUT,OUTF,'NEW',' ',' ',IREC,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open ',TRIM(OUTF)
	    WRITE(6,*)'Fille cannot exist already'
	    STOP
	  END IF
!
	  WRITE(6,*)' Append (test) to use same pgplot file'
10	  CALL GEN_IN(FILE1,'Top PGPLOT  file')
	  I=INDEX(FILE1,'(test)')
	  IF(I .NE. 0)THEN
	    TEST=.TRUE.
	    FILE1=FILE1(1:I-1)
	  END IF
	  CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_TERM,*)'Error opening 1st input FILE:',FILE1
	    GOTO 10
	  END IF
	  STRING=' '
	  DO WHILE( INDEX(STRING,'0.072 0.072 scale') .EQ. 0)
	    READ(LU_IN,'(A)')STRING
	    IF(INDEX(STRING,'%%Orientation: Landscape').NE. 0)THEN
	      STRING= '%%Orientation: Portrait'
	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
	  READ(LU_IN,'(A)')STRING   !Skip Translate/rotation string
!
	  WRITE(LU_OUT,'(A)')'  0.8 0.8 scale'
	  IF(NROWS .EQ. 2)THEN
	    WRITE(LU_OUT,'(A)')'  500 7100 translate'
	  ELSE IF(NROWS .EQ. 3)THEN
	    WRITE(LU_OUT,'(A)')'  500 8500 translate'
	  ELSE IF(NROWS .EQ. 4 .AND. NCOLS .EQ. 1)THEN
	    WRITE(LU_OUT,'(A)')' -200 9000 translate'
	  ELSE IF(NROWS .EQ. 4)THEN
	    WRITE(LU_OUT,'(A)')'  500 9700 translate'
	  ELSE IF(NROWS .EQ. 5)THEN
	    WRITE(LU_OUT,'(A)')'  -700 9500 translate'
	  ELSE
	    WRITE(LU_OUT,'(A)')'  -100 9700 translate'
	  END IF
	  DO WHILE(1 .EQ. 1)
	    READ(LU_IN,'(A)')STRING
	    IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	      GOTO 100
	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
100	  CONTINUE
	CLOSE(UNIT=LU_IN)
!
	DO NR=1,NROWS
	  DO NC=1,NCOLS
	    IF(NR .EQ. 1 .AND. NC .EQ. 1)THEN
	      GOTO 500			!Done before as special case
	    ELSE IF(NC .EQ. 1)THEN
	      WRITE(LU_OUT,'(5A)')'  -',TRIM(HOR_OFFSET),' -',TRIM(VER_OFFSET),' translate'
	      DO I=1,NCOLS-2
	        WRITE(LU_OUT,'(4A)')'  -',TRIM(HOR_OFFSET),'    0',' translate'
	      END DO
	    ELSE
	      WRITE(LU_OUT,'(4A)')'  ',TRIM(HOR_OFFSET),'    0',' translate'
	    END IF
!
! Update file name in a systematic way to save typing.
!
	    IF(.NOT. TEST)THEN
	      K=INDEX(FILE1,'_')
	      IF(K .NE. 0)K=INDEX(FILE1(K+1:),'_')+K
	      IF(K .NE. 0)THEN
	        J=K+1
	        DO WHILE(FILE1(J:J) .GE. '0' .AND. FILE1(J:J) .LE. '9')
	          J=J+1
	        END DO
	        IF(J .NE. K+1)THEN
	          READ(FILE1(K+1:J-1),*)I
	          I=I+1
	          WRITE(TMP_STR,'(I6)')I
	          TMP_STR=ADJUSTL(TMP_STR)
	          FILE1=FILE1(1:K)//TRIM(TMP_STR)//FILE1(J:)
	        END IF
	      END IF
	      CALL SET_CASE_LOW(FILE1,IZERO,IZERO)
	      IF(FILE1 .EQ. 'pgplot.ps')FILE1='pgplot_2.ps'
	    END IF
!
150	    CALL GEN_IN(FILE1,'Next PGPLOT file')
	    IF(FILE1 .EQ. ' ')GOTO 1000
	    CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LU_TERM,*)'Error opening input FILE:',FILE1
	      GOTO 150
	    END IF
	     STRING=' '
	    DO WHILE(INDEX(STRING,'%%Page: 1 1') .EQ. 0)
	      READ(LU_IN,'(A)')STRING
	    END DO
	    DO WHILE(INDEX(STRING,' EP ') .EQ. 0)
	      READ(LU_IN,'(A)')STRING
	    END DO
	    I=INDEX(STRING,' EP ')
	    WRITE(LU_OUT,'(A)')TRIM(STRING(I+4:))
	    DO WHILE(1 .EQ. 1)
	      READ(LU_IN,'(A)')STRING
	      IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	        GOTO 200
	      END IF
	      WRITE(LU_OUT,'(A)')TRIM(STRING)
	    END DO
200	    CONTINUE
	    CLOSE(UNIT=LU_IN)
500	    CONTINUE
	    END DO
	  END DO
!
1000	  CONTINUE
	  WRITE(LU_OUT,'(A)')STRING
	  WRITE(LU_OUT,'(A)')'%%EOF'
!
! Update file name in a systematic way to save typing.
!
	  IF(.NOT. TEST)THEN
	    K=INDEX(FILE1,'_')
	    IF(K .NE. 0)K=INDEX(FILE1(K+1:),'_')+K
	    IF(K .NE. 0)THEN
	      J=K+1
	      DO WHILE(FILE1(J:J) .GE. '0' .AND. FILE1(J:J) .LE. '9')
	        J=J+1
	      END DO
	      IF(J .NE. K+1)THEN
	        READ(FILE1(K+1:J-1),*)I
	        I=I+1
	        WRITE(TMP_STR,'(I6)')I
	        TMP_STR=ADJUSTL(TMP_STR)
	        FILE1=FILE1(1:K)//TRIM(TMP_STR)//FILE1(J:)
	      END IF
	    END IF
	    CALL SET_CASE_LOW(FILE1,IZERO,IZERO)
	    IF(FILE1 .EQ. 'pgplot.ps')FILE1='pgplot_2.ps'
	  END IF
!
	END DO
!
	STOP
	END
