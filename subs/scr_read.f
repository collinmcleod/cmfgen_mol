C
	SUBROUTINE SCR_READ(R,V,SIGMA,POPS,IREC,NITSF,
	1                        NUM_TIMES,LST_NG,NT,ND,LU,NEWMOD)
	IMPLICIT NONE
C
C Altered 05-Dec-1996 - INQUIRE statement used before closing files
C                       ASCI files opened by GEN_ASCI_OPEN.
C Altered 25-Jun-1996 - ACTION='READ' installed on OPEN statements.
C Altered 12-Jan-1991 - By using call to DIR_ACC_PARS this version is now
C                        compatible with both CRAY and VAX fortran.
C Altered  3-Apr-1989 - LST_NG installed. LU now transmitted in call.
C
	LOGICAL NEWMOD
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the nuber of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the numer of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER*4 LOC_NUMTIMES,ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,L,LUER,ERROR_LU
	INTEGER*4 IOS,IST,IEND
	LOGICAL FILE_OPEN
	INTEGER*4, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
C
	NEWMOD=.FALSE.
	LUER=ERROR_LU()
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C These are computer dependent, hence call to DIR_ACC_PARS. NB.
C REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
C passed to the OPEN statement.
 
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length is too small to output R,V and'//
	1              ' sigma in SCR_RITE'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
C
C Read in pointer to data file. If pointer file does not exist, or bad
C data, initialize parameters for a NEW MODEL.
C
	CALL GEN_ASCI_OPEN(LU,'POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening POINT1 in SCR_READ'
	  ELSE
	    READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	    IF(IOS .NE. 0)
	1      WRITE(LUER,*)'Error reading POINT1 in SCR_READ'
	  END IF
	  IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    CALL GEN_ASCI_OPEN(LU,'POINT2','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)WRITE(LUER,*)'Error opening POINT2 in SCR_READ'
	    IF(IOS .EQ.0)
	1        READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	    IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	      WRITE(LUER,*)'Error on reading Pointers in READ_SCRTEMP'
	      NEWMOD=.TRUE.
	      NITSF=0
	      IREC=0
	      LST_NG=-1000
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    END IF
	  END IF
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	IF(LOC_NUMTIMES .NE. NUM_TIMES)THEN
	  WRITE(LUER,'('' Error in SCR_READ - inconsistent NUM_TIMES'')')
	  STOP
	END IF
C
C		' OLD MODEL '
C
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',
	1       RECL=REC_LEN,IOSTAT=IOS,ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    RETURN
	  END IF
C
C Note that NITSF= # of successful iterations so far.
C
	  READ(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)READ(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading R,V, SIGMA vectors in READ_SCRTEMP'
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    RETURN
	  END IF
	  RECS_FOR_RV=2
C
C Read in the population data.
C
500	CONTINUE		!Try to read an earlier record.
C
C IREC ignores the number of records that it takes to write each time.
C Hence in POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
	ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	DO L=1,NUMRECS
	  IST=(L-1)*N_PER_REC+1
	  IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	  READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error on Scratch Read'
	    IREC=IREC-1
C	    IF(IREC .GE. 0)GOTO 500	!Get another record if one is available.
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    RETURN
	  END IF
	END DO
C
C Successful Read !
C
	CLOSE(UNIT=LU)
	RETURN
	END
C
C 
C
C Routine  to save population data. There is no limit on the
C size of the POPS array.
C
	SUBROUTINE SCR_RITE(R,V,SIGMA,POPS,
	1                      IREC,NITSF,NUM_TIMES,LST_NG,
	1                      NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
C
C Aletered 25-Jun-1996 : GEN_ASCI_OPEN installed to OPEN POINT files.
C
	LOGICAL WRITFAIL
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the nuber of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the numer of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER*4 ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,K,L,LUER,ERROR_LU
	INTEGER*4 IOS,IST,IEND
	INTEGER*4, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C These are computer dependent, hence call to DIR_ACC_PARS. NB.
C REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
C passed to the OPEN statement.
C
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length is too small to output R,V and'//
	1              ' sigma in SCR_RITE'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND			
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
C
C*************************************************************************
C
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP in WRITE_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      WRITFAIL=.TRUE.
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    ELSE
	      IF(IOS .EQ. 0)IREC=0		!Since new file.
	    END IF
	  END IF
C
	IF(IREC .EQ. 0)THEN		!Newfile or newmodel
	  WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	    WRITFAIL=.TRUE.
	    RETURN
	  END IF
	END IF
	RECS_FOR_RV=2
C
C WRITE in the population data.
C
	IREC=IREC+1		!Next record output
	DO K=1,NUM_TIMES	!Write out POPS NUM_TIMES for safety.
C
C IREC ignores the number of records that it takes to write each time.
C Hence in POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)THEN
	      WRITFAIL=.TRUE.
	      WRITE(LUER,*)'Error writing SCRTEMP in SCR_RITE'
	      WRITE(LUER,*)'IOS=',IOS
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    END IF
	  END DO
	END DO
C
C Successful write.
C
	WRITFAIL=.FALSE.
1000	CLOSE(UNIT=LU)
C
C Write pointer to data files.
C
	CALL GEN_ASCI_OPEN(LU,'POINT1','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(1X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	CALL GEN_ASCI_OPEN(LU,'POINT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(1X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	RETURN
	END
