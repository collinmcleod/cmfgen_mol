!
! Subroutine to rewrite the free-format SCRTEMP file. The user
! input how many of the final iterations are to be output.
! NEW_SCRTEMP, NEW_POINT1, and NEW_POINT1 are created. These need
! to be renamed. This allos space to be saved, and at the same time 
! allows a new model to be run with exactly the same input populations.
!
	PROGRAM REWRITE_SCR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Revised: 13-Feb-2002
!
	INTEGER*4 ND,NT,NIT
!
	REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
!
	INTEGER*4, PARAMETER :: T_OUT=6
!
	INTEGER*4 NPLTS
	INTEGER*4 IREC
	INTEGER*4 IREC_WR
	INTEGER*4 NIT_WR
	INTEGER*4 NITSF
	INTEGER*4 LST_NG
	INTEGER*4 WRITEN_N_TIMES
	INTEGER*4 RITE_N_TIMES
	INTEGER*4 LUSCR
	INTEGER*4 IOS
	INTEGER*4 I,K,J,L
!
	CHARACTER*132 STRING
!
	LUSCR=26
	RITE_N_TIMES=1
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       POINT1.DAT'
	WRITE(T_OUT,*)'                                       SCRTEMP.DAT'
	WRITE(T_OUT,*)'                                       MODEL.DAT'
	WRITE(T_OUT,*)' '
!
!
! Get NT (# of levels) and ND (number of depth points) from MODEL_SPEC.
!
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
!
100	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to read MODEL file'
	    CALL GEN_IN(NT,'Total number of levels')
	    CALL GEN_IN(ND,'Number of depth points')
	  END IF 
	CLOSE(UNIT=12)
!
! Read record information from SCRTEMP information file POINT1.
!
	OPEN(UNIT=12,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)READ(12,*,IOSTAT=IOS)K,NIT,WRITEN_N_TIMES
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Possible error reading POINT1'
	    CALL GEN_IN(NIT,'Number of iterations')
	  END IF
	CLOSE(UNIT=12)
!
	ALLOCATE (POPS(NT,ND))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
!
!
! NIT_WR=2 means that we write out the LAST 2 iterayions only.
!
	NIT_WR=2
	CALL GEN_IN(NIT_WR,'Number of (FINAL) iterations to write')
	IF(NIT_WR .GT. NIT)THEN
	   WRITE(6,*)'Error -- asking to write more iterations than available'
	   STOP
	END IF
!
	IREC_WR=0
	DO I=1,NIT_WR
	   IREC=NIT-(NIT_WR-I)
	   CALL READ_SCR_REC(R,V,SIGMA,POPS,IREC,NITSF,
	1              WRITEN_N_TIMES,LST_NG,
	1              NT,ND,LUSCR)
!
	   LST_NG=-1000
	   CALL RE_RITE_SCR(R,V,SIGMA,POPS,IREC_WR,I,
	1              RITE_N_TIMES,LST_NG,
	1              NT,ND,LUSCR)
!
	END DO
!
	STOP
	END
! 
!
! Routine  to save population data. There is no limit on the
! size of the POPS array. Data is output to NEW_SCRTEMP,
! with iteration information in NEW_POINT1 (& NEW_POINT2).
!
	SUBROUTINE RE_RITE_SCR(R,V,SIGMA,POPS,
	1                      IREC,NITSF,NUM_TIMES,LST_NG,
	1                      NT,ND,LU)
	IMPLICIT NONE
!
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
!
! Local variables.
!
! REC_SIZE     is the (maximum) record length in bytes.
! REC_LEN      is the record length in computer units.
! UNIT_SIZE    is the nuber of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the numer of bytes used to represent the number.
! NUMRECS      is the # of records required to output POPS.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
!	
	INTEGER*4 ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,K,L,LUER
	INTEGER*4 IOS,IST,IEND
	LOGICAL FILE_OPEN
!
! Ouput all errors to terminal
!
	LUER=6
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! These are computer dependent, hence call to DIR_ACC_PARS. NB.
! REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
! passed to the OPEN statement.
!
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
!*************************************************************************
!
	OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening NEW_SCRTEMP in WRITE_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening NEW_SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      STOP
	    ELSE
	      IF(IOS .EQ. 0)IREC=0		!Since new file.
	    END IF
	  END IF
!
! If we opened a NEW file, output R, V and SIGMA
!
	IF(IREC .EQ. 0)THEN		!Newfile or newmodel
	  WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END IF
	RECS_FOR_RV=2
!
! WRITE in the population data.
!
	IREC=IREC+1		!Next iteration output
	ST_REC_M1=(IREC-1)*NUMRECS*NUM_TIMES+RECS_FOR_RV
	DO K=1,NUM_TIMES	!Write out POPS NUM_TIMES for safety.
!
! IREC ignores the number of records that it takes to write each time.
! Hence in POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error writing NEW_SCRTEMP in RE_RITE_SCR'
	      WRITE(LUER,*)'IOS=',IOS
	      CLOSE(UNIT=LU)
	      STOP
	    END IF
	  END DO
	  ST_REC_M1=ST_REC_M1+NUMRECS
	END DO
!
! Successful write.
!
1000	CLOSE(UNIT=LU)
!
! Write pointer to data files. 
!
	OPEN(UNIT=LU,FILE='NEW_POINT1',STATUS='UNKNOWN')
	  WRITE(LU,'(X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	CLOSE(UNIT=LU)
	OPEN(UNIT=LU,FILE='NEW_POINT2',STATUS='UNKNOWN')
	  WRITE(LU,'(X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	CLOSE(UNIT=LU)
!
	RETURN
	END
!
! 
!
! Routine to read in iteration IREC from SCRTEMP file. It assumes
! that each iteration is written only once. Note the distinction
! with SCR_READ, which reads in the LAST iteration, and returns IREC.
!
	SUBROUTINE READ_SCR_REC(R,V,SIGMA,POPS,IREC,NITSF,
	1                        NUM_TIMES,LST_NG,NT,ND,LU)
	IMPLICIT NONE
!
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
!
! Local variables.
!
! REC_SIZE     is the (maximum) record length in bytes.
! REC_LEN      is the record length in computer units.
! UNIT_SIZE    is the nuber of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the numer of bytes used to represent the number.
! NUMRECS      is the # of records required to output POPS.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
!
	INTEGER*4 ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,L,LUER
	INTEGER*4 IOS,IST,IEND
!
	LUER=6
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! These are computer dependent, hence call to DIR_ACC_PARS. NB.
! REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
! passed to the OPEN statement.

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
!		' OLD MODEL '
!
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',
	1       RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    CLOSE(UNIT=LU)
	    STOP
	  END IF
!
! Note that NITSF= # of successful iterations so far.
!
	  READ(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)READ(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading R,V, SIGMA vectors in READ_SCRTEMP'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  RECS_FOR_RV=2
!
! Read in the population data.
!
500	CONTINUE		!Try to read an earlier record.
!
! IREC ignores the number of records that it takes to write each time.
! Hence in POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
	ST_REC_M1=(IREC-1)*NUMRECS*NUM_TIMES+RECS_FOR_RV
	DO L=1,NUMRECS
	  IST=(L-1)*N_PER_REC+1
	  IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	  READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error on Scratch Read'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
!
! Successful Read !
!
	CLOSE(UNIT=LU)
	RETURN
!
	END
