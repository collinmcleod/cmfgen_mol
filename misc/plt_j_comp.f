	PROGRAM PLT_J_COMP
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NCF_MAX=150000
!
	REAL*8 INDX(NCF_MAX)
	REAL*8 NU(NCF_MAX)
!
	REAL*8 JM_OUT(NCF_MAX)
	REAL*8 JR_OUT(NCF_MAX)
	REAL*8 DIFF_OUT(NCF_MAX)
	REAL*8 HBC_CMF(NCF_MAX)
!
	REAL*8 JM_IN(NCF_MAX)
	REAL*8 JR_IN(NCF_MAX)
	REAL*8 DIFF_IN(NCF_MAX)
!
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
!
	REAL*8 SUM1
	REAL*8 SUM2
	REAL*8 SUM3
!
	INTEGER I,IBEG
	INTEGER NCF
	INTEGER LUM_STAR 
!
1000	CONTINUE
 	FILENAME='J_COMP'
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	IF(FILENAME .EQ. ' ')GOTO 1000
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'HBC_CMF') .EQ. 0)
	    READ(11,'(A)')STRING
	  END DO
!
	  I=0
	  DO WHILE(1. EQ. 1)
	    I=I+1
	    IF(I .GT. NCF_MAX)THEN
	      WRITE(6,*)'************************************************'
	      WRITE(6,*)' '
	      WRITE(6,*)'Error --- insufficent meory to read in all data'
	      WRITE(6,*)' '
	      WRITE(6,*)'************************************************'
	      EXIT
	    END IF
	    READ(11,*,END=100)INDX(I),NU(I), JM_OUT(I),JR_OUT(I),DIFF_OUT(I),HBC_CMF(I),
	1                                      JM_IN(I),JR_IN(I),DIFF_IN(I)
	  END DO
!
100	  CONTINUE
	  NCF=I-1
	CLOSE(UNIT=11)
!
	CALL DP_CURVE(NCF,NU,DIFF_OUT)
	CALL GRAMON_PGPLOT('\gn(10\u15\ \dHz)','% Difference',' ',' ')
!
	DO I=1,NCF
	  NU(I)=DLOG10(NU(I))
	END DO
	CALL DP_CURVE(NCF,NU,DIFF_OUT)
	CALL GRAMON_PGPLOT('\gn(10\u15\ \dHz)','% Difference',' ',' ')
!
	SUM1=0.0D0
	SUM2=0.0D0
	SUM3=0.0D0
	SUM1=SUM(DIFF_OUT(1:NCF))/NCF
	SUM3=SUM(ABS(DIFF_OUT(1:NCF)))/NCF
	SUM2=SUM(DIFF_OUT(1:NCF)*DIFF_OUT(1:NCF))
	SUM2=SQRT(SUM2/NCF)
!
	WRITE(6,*)'Mean difference is ', SUM1
	WRITE(6,*)'Mean absolute difference is ', SUM3
	WRITE(6,*)'Root mean square is',SUM2
	  SUM1=SUM1
!
	STOP
	END
