	PROGRAM PLT_J_COMP
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NCF_MAX=400000
!
	REAL*8 INDX(NCF_MAX)
	REAL*8 NU(NCF_MAX)
	REAL*8 LAM(NCF_MAX)
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
        CHARACTER*30 UC
        EXTERNAL UC
!
	REAL*8 SUM1
	REAL*8 SUM2
	REAL*8 SUM3
	CHARACTER*20 PLT_OPT
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
	LAM(1:NCF)=2.99702458D+03/NU(1:NCF)
!
	SUM1=0.0D0
	SUM2=0.0D0
	SUM3=0.0D0
	SUM1=SUM(DIFF_OUT(1:NCF))/NCF
	SUM3=SUM(ABS(DIFF_OUT(1:NCF)))/NCF
	SUM2=SUM(DIFF_OUT(1:NCF)*DIFF_OUT(1:NCF))
	SUM2=SQRT(SUM2/NCF)
!
!
	WRITE(6,*)' '
	WRITE(6,'(A,F6.2,A)')'Mean difference is:          ',SUM1,'%'
	WRITE(6,'(A,F6.2,A)')'Mean absolute difference is: ',SUM3,'%'
	WRITE(6,'(A,F6.2,A)')'Root mean square is:         ',SUM2,'%'
!
	PLT_OPT='DW'
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  WRITE(6,*)'Plot options are:'
	  WRITE(6,*)' JF:  Plot J against frequency at outer boundary'
	  WRITE(6,*)' JW:  Plot J against wavelength at outer boundary'
	  WRITE(6,*)' DF:  Plot Error against frequency at outer boundary'
	  WRITE(6,*)' DW:  Plot Error against wavelength at outer boundary'
	  WRITE(6,*)' E(X):  Exit routine'
	  CALL GEN_IN(PLT_OPT,'Plot option')
	  IF(UC(PLT_OPT) .EQ. 'JF')THEN
	    CALL DP_CURVE(NCF,NU,JR_OUT)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','J',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'JW')THEN
	    CALL DP_CURVE(NCF,LAM,JR_OUT)
	    CALL GRAMON_PGPLOT('\g(\A)','J',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'DF')THEN
	    CALL DP_CURVE(NCF,NU,DIFF_OUT)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','% Difference',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'DW')THEN
	    CALL DP_CURVE(NCF,LAM,DIFF_OUT)
	    CALL GRAMON_PGPLOT('\gl(\A)','% Difference',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
	END DO
!
	END
