!
! Routine performs an AVERAGE of the last 2 iterations for a Comoving-Frame 
! Model. Progam uses the last 2 iterations which are stored in the last "2 records" 
! (effectively) of SCRTEMP.
!
! Output is to the next record (now last) of SCRTEMP.
!
! Input files required:
!		     MODEL(.DAT)
!                    POINT1(.DAT)
!                    SCRTEMP(.DAT)
! Output files:
!                    POINT1(.DAT)
!                    POINT2(.DAT)
!                    SCRTEMP(.DAT)
!
	PROGRAM DO_AVE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 15-Apr-2003 : IFLAG in RD_4_ITS initialized.
! Altered 18-Feb-2001 : Extensive rewrite.
!                       NBAND option installed.
!
! Altered 04-Jan-1998 : Cleaned. ND, NT now read from MODEL file.
!                       Based on a very old version in [JDH.BASOL] and 
!                        REWRITE_SCR.
!
	REAL*8, ALLOCATABLE :: RDPOPS(:,:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	INTEGER*4 ND,NT
	INTEGER*4 NBAND
!
	INTEGER*4 IOS
	INTEGER*4 IREC
	INTEGER*4 NITSF
	INTEGER*4 LST_NG
	INTEGER*4 IFLAG
	INTEGER*4 I,J
!
	INTEGER*4, PARAMETER :: RITE_N_TIMES=1
	INTEGER*4, PARAMETER :: T_OUT=6
	INTEGER*4, PARAMETER :: LUSCR=10
!
	LOGICAL NEWMOD
	LOGICAL NG_DONE
	LOGICAL DO_REGARDLESS
	LOGICAL SCALE_INDIVIDUALLY
	CHARACTER(LEN=132) STRING
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routines will average the last 2 iterations & '
	WRITE(T_OUT,*)'output the result to SCRTEMP'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       MODEL(.DAT)'
	WRITE(T_OUT,*)'                                       POINT1(.DAT)'
	WRITE(T_OUT,*)'                                       SCRTEMP(.DAT)'
	WRITE(T_OUT,*)' '
!
! Get basic Model data (i.e., NT and ND).
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
! Iw we couldn't successfully red the MODEL file, get NT and ND from terminal.
!
100     CONTINUE
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	END IF 
	CLOSE(UNIT=12)
!
	ALLOCATE (RDPOPS(NT,ND,4))
	ALLOCATE (POPS(NT,ND))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
! 
	NEWMOD=.FALSE.
	CALL SCR_READ(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                NT,ND,LUSCR,NEWMOD)
!
! Read in the last 4 estimates of the poplations, as output to SCRTEMP.
!
	CALL RD_4_ITS(RDPOPS,POPS,NT,ND,IFLAG)
	IF(IFLAG .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read scratch file'
	  STOP
	END IF
!
	DO J=1,ND
	  DO I=1,NT
	    POPS(I,J)=0.5D0*(RDPOPS(I,J,1)+RDPOPS(I,J,2))
	  END DO
	END DO
!
	NITSF=NITSF+1
	LST_NG=NITSF
	CALL SCR_RITE(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,
	1                 LST_NG,NT,ND,LUSCR,NEWMOD)
	WRITE(T_OUT,*)'Results of successful AVERAGE',
	1                 'output to SCRTEMP.'
!
	STOP
	END
!
!
!
	SUBROUTINE RD_4_ITS(RDPOPS,POPS,NT,ND,IFLAG)
	IMPLICIT NONE
!
	INTEGER*4 NT
	INTEGER*4 ND
	INTEGER*4 IFLAG
	REAL*8 RDPOPS(NT*ND,4)
	REAL*8 POPS(NT*ND)
!
! Note that REC_SIZE is the size of the output record in bytes.
!           REC_LEN  is the size of the output record in COMPUTER units.
!           N_PER_REC is the number of numbers per record.
!           NUM_TIME is the number of times each iteration has been written
!                to the scratch file.
!           WORD_SIZE is the size of the number to be output in bytes.
!           UNIT_SIZE is the number of bytes per unit used to specify
!           the size of a direct access file.
!
	INTEGER*4 ARRAYSIZE,IST,IEND
	INTEGER*4 REC_SIZE,REC_LEN,NUM_RECS,N_PER_REC,NUM_TIME
	INTEGER*4 UNIT_SIZE,WORD_SIZE
!
	INTEGER*4 I,L,NPREV,NITSF,IT_CNT
	INTEGER*4 SRT_REC_M1
	INTEGER*4 CHK
	INTEGER*4 LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	IFLAG=0
        LUER=ERROR_LU()
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! As this is computer and installation dependent, we call a subroutine
! to return the parameters.
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,
	1                 WORD_SIZE,N_PER_REC)
	ARRAYSIZE=NT*ND
	NUM_RECS=INT( (ARRAYSIZE-1)/N_PER_REC ) + 1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! Read in pointer to data file. NB. Format has changed. IT_CNT is now
! the iteration count (but may differ from NITSF is problems occurred
! read/riting to SCRTEMP) where as previously it referred to the actual 
! record.
!
	OPEN(UNIT=26,FILE='POINT1',IOSTAT=CHK,STATUS='OLD')
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error opening POINT1 in GENACCEL: IOSTAT=',CHK
	    IFLAG=3
	    RETURN
	  END IF
	  READ(26,*,IOSTAT=CHK)IT_CNT,NITSF,NUM_TIME
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error reading POINT1 in GENACCEL: IOSTAT=',CHK
	    IFLAG=4
	    CLOSE(UNIT=26)
	    RETURN
	  END IF
	CLOSE(UNIT=26)
!
!		' OLD MODEL '
!
	OPEN(UNIT=20,FILE='SCRTEMP',FORM='UNFORMATTED',
	1    ACCESS='DIRECT',STATUS='OLD',RECL=REC_LEN,IOSTAT=CHK)
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input in GENACCEL'
	    WRITE(LUER,*)'IOSTAT=',CHK
	    CLOSE(UNIT=20)
	    IFLAG=5
	    RETURN
	  END IF
!
! SRT_REC_M1 + 1 is the first record for the NPREV previous iteration.
! If output more than once, it refers to the first output.
!
	  DO NPREV=1,4
	    SRT_REC_M1=(IT_CNT-NPREV)*NUM_TIME*NUM_RECS+2
	    DO L=1,NUM_RECS
	      IST=(L-1)*N_PER_REC+1
	      IEND=MIN( IST+N_PER_REC-1,ARRAYSIZE )
	      READ(20,REC=SRT_REC_M1+L,IOSTAT=CHK)
	1          (RDPOPS(I,NPREV),I=IST,IEND)
	      IF(CHK .NE. 0)THEN
	        WRITE(LUER,*)'Error on reading scratch file in GENACCEL'
	        WRITE(LUER,200)SRT_REC_M1,CHK
200	        FORMAT(X,'SRT_REC_M1=',I3,5X,'IOSAT=',I7)
	        IFLAG=6
	        CLOSE(UNIT=20)
	        RETURN
	      END IF
	    END DO
	  END DO
	  CLOSE(UNIT=20)
!
	POPS(:)=RDPOPS(:,1)
!
	WRITE(6,*)'SCRTEMP file was successfully read'
	RETURN
	END
