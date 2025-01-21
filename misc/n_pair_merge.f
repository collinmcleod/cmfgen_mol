!
! Routine to merge 6 color PGPLOTS onto a single
! page. The plots should be done with:
!      EXPAND_CHAR=1.2; EXPAND_TICK=1.2; ASR=0.89; Plot Size=10.0 cm
! The plots are located in 2 columns, 2 per row.
!
! Landscape format, with CPS as printer.
!
	PROGRAM N_PAIR_MERGE
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
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
	INTEGER IOS,I,J,K,IREC,M,NPLTS
	CHARACTER*10 TMP_STR
!
	WRITE(LU_TERM,*)' '
	WRITE(LU_TERM,*)'Program to merge N plots in 2 columns'
	WRITE(LU_TERM,*)'6 plots (LM): ',
	1    'EXPAND_CHAR=1.2; EXPAND_TICK=1.2; ASR=0.89; Plot Size=9.5 cm'
	WRITE(LU_TERM,*)' '
!
	CALL GEN_IN(OUTF,'Output file')
	CALL GEN_IN(NPLTS,'Total number of plts on page')
	CALL GEN_IN(FILE1,'Top PGPLOT  file')
!
	IREC=0
	IOS=0
	CALL GEN_ASCI_OPEN(LU_OUT,OUTF,'NEW',' ',' ',IREC,IOS)
!
10	CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
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
	  WRITE(LU_OUT,'(A)')'  500 8500 translate'
	  DO WHILE(1 .EQ. 1)
	    READ(LU_IN,'(A)')STRING
	    IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	      GOTO 100
	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
100	  CONTINUE
	CLOSE(UNIT=LU_IN)
C
	DO M=1,NPLTS-1
	  IF( MOD(M,2) .EQ. 0)THEN
	    IF(NPLTS .EQ. 4)WRITE(LU_OUT,'(A)')'  -4500 -5300 translate'
	    IF(NPLTS .EQ. 6)WRITE(LU_OUT,'(A)')'  -4500 -4000 translate'
	    IF(NPLTS .EQ. 8)WRITE(LU_OUT,'(A)')'  -4500 -2700 translate'
	  ELSE
	    WRITE(LU_OUT,'(A)')'  4500 0 translate'
	  END IF
!
! Update file name in a systematic way to save typing.
!
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
!	
150	  CALL GEN_IN(FILE1,'Next PGPLOT file')
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
200	  CONTINUE
	  CLOSE(UNIT=LU_IN)
	END DO
C
1000	CONTINUE
	WRITE(LU_OUT,'(A)')STRING
	WRITE(LU_OUT,'(A)')'%%EOF'
C
	STOP
	END
