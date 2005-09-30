C
	SUBROUTINE SCR_READ_V2(R,V,SIGMA,POPS,
	1             IREC_RD,NITSF,NUM_TIMES,LST_NG,RVSIG_WRITTEN,
	1             NT,ND,LU,NEWMOD)
	IMPLICIT NONE
C
C Altered 02-May-2004 : Changed to V2.
C                          RVSIG_WRITTEN inserted into call.
C                          Only writes R, V & SIGMA when needed, but all
C                          writes to file must be identical.
C Altered 29-Feb-2004 - Adjusted to read new format file for SCRTEMP.
C                          R, V, & SIGMA can be output for every iteration.
C                          IREC_RD must be zero to read last iteration,
C                          other wise iteration is IREC.
C                          NUM_TIMES in no longer used.
C Altered 05-Dec-1996 - INQUIRE statement used before closing files
C                       ASCI files opened by GEN_ASCI_OPEN.
C Altered 25-Jun-1996 - ACTION='READ' installed on OPEN statements.
C Altered 12-Jan-1991 - By using call to DIR_ACC_PARS this version is now
C                        compatible with both CRAY and VAX fortran.
C Altered  3-Apr-1989 - LST_NG installed. LU now transmitted in call.
C
	LOGICAL NEWMOD
!
!                          On input,  IREC_RD is the iteration to be read.
!                          On output, IREC_RD is the actual iteration read.
!
	INTEGER IREC_RD
	INTEGER NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
	LOGICAL RVSIG_WRITTEN
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the number of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the number of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER REC_SIZE,REC_LEN
	INTEGER NUMRECS
	INTEGER N_PER_REC
	INTEGER ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER LOC_NUMTIMES,ST_REC_M1,RECS_FOR_RV
	INTEGER I,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER IREC
	LOGICAL FILE_OPEN
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
	CHARACTER*80 STRING
	CHARACTER*11 FORMAT_DATE
C
	NEWMOD=.FALSE.
	LUER=ERROR_LU()
	FORMAT_DATE=' '
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C These are computer dependent, hence call to DIR_ACC_PARS. NB.
C REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
C passed to the OPEN statement.
 
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length was too small fro R,V and'//
	1              ' sigma in SCR_READ'
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
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
              READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG,RVSIG_WRITTEN
	    ELSE IF(IOS .EQ. 0)THEN
	      RVSIG_WRITTEN=.FALSE.
              READ(STRING,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	    END IF
	    IF(IOS .NE. 0)WRITE(LUER,*)'Error reading POINT1 in SCR_READ'
	  END IF
!
	  IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    CALL GEN_ASCI_OPEN(LU,'POINT2','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening POINT2 in SCR_READ'
	    ELSE
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	        STRING=ADJUSTL(STRING)
	        FORMAT_DATE=STRING(1:11)
                READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG,RVSIG_WRITTEN
	      ELSE IF(IOS .EQ. 0)THEN
	        RVSIG_WRITTEN=.FALSE.
                READ(STRING,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	      END IF
	    END IF
	  END IF
!
	IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	  WRITE(LUER,*)'Error on reading Pointers in READ_SCRTEMP'
	  NEWMOD=.TRUE.
	  NITSF=0
	  LST_NG=-1000
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
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
	    LST_NG=-1000
	    IREC_RD=0
	    RETURN
	  END IF
	  RECS_FOR_RV=2
C
C Read in the population data.
C
500	CONTINUE		!Try to read an earlier record.
!
! IREC ignores the number of records that it takes to write each time.
! Hence in POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
	IF(IREC_RD .NE. 0)THEN
	  IREC=IREC_RD
	ELSE
	  IREC_RD=IREC
	END IF
	IF(RVSIG_WRITTEN)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	  READ(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	  ST_REC_M1=ST_REC_M1+1
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
!
	IF(IOS .EQ. 0)THEN
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)EXIT
	  END DO
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error on Scratch Read'
	  IREC=IREC-1
	  NEWMOD=.TRUE.
	  NITSF=0
	  LST_NG=-1000
	  IREC_RD=0
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
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
	SUBROUTINE SCR_RITE_V2(R,V,SIGMA,POPS,
	1               IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG,
	1               NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
C
C Altered 02-May-2004 : Changed to V2.
C                          RVSIG_WRITTEN inserted into call.
C                          Only writes R, V & SIGMA when needed, but all
C                          writes to file must be identical.
C Altered 29-Feb-2004 - Adjusted to read new format file for SCRTEMP.
C                          R, V, & SIGMA now output for every iteration.
C                          NUM_TIMES in no longer used.
C Altered 27-Feb-2004 : R,V, and SIGMA rewritten on every iteration.
C Altered 25-Jun-1996 : GEN_ASCI_OPEN installed to OPEN POINT files.
C
	LOGICAL WRITFAIL
	LOGICAL WRITE_RVSIG
	INTEGER IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the number of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the number of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER REC_SIZE,REC_LEN
	INTEGER NUMRECS
	INTEGER N_PER_REC
	INTEGER ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER ST_REC_M1,RECS_FOR_RV
	INTEGER I,K,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
	LOGICAL LOC_WR_RVSIG
	EXTERNAL ERROR_LU
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER*80 STRING
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
!
! We check FORMAT associated with SCRATCH file. This will allow us to
! preserve the same format for an existing file.
!
	IF(IREC .NE. 1)THEN
!
! Writing to existing SCRATCH file.
!
	  FORMAT_DATE=' '
	  CALL GEN_ASCI_OPEN(LU,'POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1     CALL GEN_ASCI_OPEN(LU,'POINT2','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
	      READ(LU,*)I,I,I,I,LOC_WR_RVSIG
	      IF(LOC_WR_RVSIG .NEQV. WRITE_RVSIG)THEN
	        WRITE(LUER,*)'Error in SCR_RITE_V2 -- inconsistent WR_RVSIG option'
	        WRITE(LUER,*)'Restart a fresh model by deleting SCRTEMP etc.'
	        STOP
	      END IF
	    END IF
	  END IF
	  IF(IOS .NE. 0)FORMAT_DATE='28-Feb-2004'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	ELSE
!
! Create a new scratch file.
!
	  FORMAT_DATE='28-Feb-2004'
	END IF
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
	      FORMAT_DATE='28-Feb-2004'
	      IREC=0				!Since new file.
	    END IF
	  END IF
!
! We now write out R,V and SIGMA on very iteration. This is also 
! done for SN model where R can change with iteration. R iread here will
! only be correct for the last iteration. FIxewd with a later read.
!
	WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	  WRITFAIL=.TRUE.
	  RETURN
	END IF
	RECS_FOR_RV=2
C
C WRITE in the population data.
C
	IREC=IREC+1		!Next record output
C
C IREC ignores the number of records that it takes to write each time.
C Hence in POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
C NB: We cannot change the type of file that is being written. If we do
C     SCRTEMP will be corrupted.
C
	IF(WRITE_RVSIG)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	  WRITE(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	  ST_REC_M1=ST_REC_M1+1
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
	IF(IOS .EQ. 0)THEN
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)EXIT
	  END DO
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITFAIL=.TRUE.
	  WRITE(LUER,*)'Error writing SCRTEMP in SCR_RITE'
	  WRITE(LUER,*)'IOS=',IOS
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
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
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	CALL GEN_ASCI_OPEN(LU,'POINT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	RETURN
	END
C
C
C Routine  to save population data. There is no limit on the
C size of the POPS array. This useses the prefix NEW_ attached
C to POINT and SCRTEMP. Otherwise it is identical to SCR_RITE_V2
C
	SUBROUTINE SCR_RITE_NAM_V2(R,V,SIGMA,POPS,
	1               IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG,
	1               NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
C
C Altered 02-May-2004 : Changed to V2.
C                          RVSIG_WRITTEN inserted into call.
C                          Only writes R, V & SIGMA when needed, but all
C                          writes to file must be identical.
C Altered 29-Feb-2004 - Adjusted to read new format file for NEW_SCRTEMP.
C                          R, V, & SIGMA now output for every iteration.
C                          NUM_TIMES in no longer used.
C Altered 27-Feb-2004 : R,V, and SIGMA rewritten on every iteration.
C Altered 25-Jun-1996 : GEN_ASCI_OPEN installed to OPEN NEW_POINT files.
C
	LOGICAL WRITFAIL
	LOGICAL WRITE_RVSIG
	INTEGER IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the number of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the number of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER REC_SIZE,REC_LEN
	INTEGER NUMRECS
	INTEGER N_PER_REC
	INTEGER ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER ST_REC_M1,RECS_FOR_RV
	INTEGER I,K,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
	LOGICAL LOC_WR_RVSIG
	EXTERNAL ERROR_LU
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER*80 STRING
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
!
! We check FORMAT associated with SCRATCH file. This will allow us to
! preserve the same format for an existing file.
!
	IF(IREC .NE. 1)THEN
!
! Writing to existing SCRATCH file.
!
	  FORMAT_DATE=' '
	  CALL GEN_ASCI_OPEN(LU,'NEW_POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1     CALL GEN_ASCI_OPEN(LU,'NEW_POINT2','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
	      READ(LU,*)I,I,I,I,LOC_WR_RVSIG
	      IF(LOC_WR_RVSIG .NE. WRITE_RVSIG)THEN
	        WRITE(LUER,*)'Error in SCR_RITE_NAM_V2 -- inconsistent WR_RVSIG option'
	        WRITE(LUER,*)'Restart a fresh model by deleting NEW_SCRTEMP etc.'
	        STOP
	      END IF
	    END IF
	  END IF
	  IF(IOS .NE. 0)FORMAT_DATE='28-Feb-2004'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	ELSE
!
! Create a new scratch file.
!
	  FORMAT_DATE='28-Feb-2004'
	END IF
C
C*************************************************************************
C
	OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening NEW_SCRTEMP in WRITE_NEW_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening NEW_SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      WRITFAIL=.TRUE.
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    ELSE
	      FORMAT_DATE='28-Feb-2004'
	      IREC=0				!Since new file.
	    END IF
	  END IF
!
! We now write out R,V and SIGMA on very iteration. This is also 
! done for SN model where R can change with iteration. R iread here will
! only be correct for the last iteration. FIxewd with a later read.
!
	WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	  WRITFAIL=.TRUE.
	  RETURN
	END IF
	RECS_FOR_RV=2
C
C WRITE in the population data.
C
	IREC=IREC+1		!Next record output
C
C IREC ignores the number of records that it takes to write each time.
C Hence in NEW_POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
C NB: We cannot change the type of file that is being written. If we do
C     NEW_SCRTEMP will be corrupted.
C
	IF(WRITE_RVSIG)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	  WRITE(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	  ST_REC_M1=ST_REC_M1+1
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
	IF(IOS .EQ. 0)THEN
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)EXIT
	  END DO
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITFAIL=.TRUE.
	  WRITE(LUER,*)'Error writing NEW_SCRTEMP in SCR_RITE'
	  WRITE(LUER,*)'IOS=',IOS
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
C
C Successful write.
C
	WRITFAIL=.FALSE.
1000	CLOSE(UNIT=LU)
C
C Write pointer to data files.
C
	CALL GEN_ASCI_OPEN(LU,'NEW_POINT1','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in NEW_POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	CALL GEN_ASCI_OPEN(LU,'NEW_POINT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in NEW_POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
C
	RETURN
	END
