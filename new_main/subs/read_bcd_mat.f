	SUBROUTINE READ_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,
	1            POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYPE)
	IMPLICIT NONE
!
	INTEGER N,NION,DEPTH_INDX
	REAL*8 B(N,N)
	REAL*8 C(N,N)
	REAL*8 D(N,N)
	REAL*8 ROW_SF(N)
	REAL*8 COL_SF(N)
	REAL*8 POPS(N)
	INTEGER IPIVOT(N)
	LOGICAL REPLACE_EQ(NION)
	LOGICAL ZERO_STEQ(N)
	CHARACTER*(*) TYPE
!
	INTEGER, SAVE :: LU=200
	INTEGER  LUER,ERROR_LU,IOS,I
	EXTERNAL ERROR_LU
	LOGICAL LU_USED
	CHARACTER*12 FILENAME
!
! Get an unused UNIT number.
!
	LU=LU-1;  LU_USED=.TRUE.
	DO WHILE(LU_USED)
	  LU=LU+1
	  INQUIRE(UNIT=LU,OPENED=LU_USED)
	END DO
!
	FILENAME=TRIM(TYPE)//'SCRATCH'
	I=LEN_TRIM(FILENAME)+1
	WRITE(FILENAME(I:I+2),'(I3.3)')DEPTH_INDX
        OPEN(UNIT=LU,FORM='UNFORMATTED',FILE=FILENAME,IOSTAT=IOS,
	1             ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in READ_BCD_MAT'
	    WRITE(LUER,*)'Unable to open ',FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	  IF(TYPE .EQ. 'BC')THEN
	    READ(LU,IOSTAT=IOS)B,C,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ
	  ELSE IF(TYPE .EQ. 'C')THEN
	    READ(LU,IOSTAT=IOS)C,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ
	  ELSE IF(TYPE .EQ. 'D')THEN
	    READ(LU,IOSTAT=IOS)D
	  ELSE
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in READ_BCD_MAT'
	    WRITE(LUER,*)'Unrecognized TYPE: TYPE=',TYPE
	  END IF
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in READ_BCD_MAT'
	    WRITE(LUER,*)'Unable to read ',FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	CLOSE(LU)
!
	RETURN
	END
!
!
	SUBROUTINE WRITE_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,
	1             POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYPE)
	IMPLICIT NONE
!
	INTEGER N,NION,DEPTH_INDX
	REAL*8 B(N,N)
	REAL*8 C(N,N)
	REAL*8 D(N,N)
	REAL*8 ROW_SF(N)
	REAL*8 COL_SF(N)
	REAL*8 POPS(N)
	INTEGER IPIVOT(N)
	LOGICAL REPLACE_EQ(NION)
	LOGICAL ZERO_STEQ(N)
	CHARACTER*(*) TYPE
!
	INTEGER, SAVE :: LU=200
	INTEGER  LUER,ERROR_LU,IOS,I
	EXTERNAL ERROR_LU
	LOGICAL LU_USED
	CHARACTER*12 FILENAME
!
! Get an unused UNIT number.
!
	LU=LU-1;  LU_USED=.TRUE.
	DO WHILE(LU_USED)
	  LU=LU+1
	  INQUIRE(UNIT=LU,OPENED=LU_USED)
	END DO
!
	FILENAME=TRIM(TYPE)//'SCRATCH'
	I=LEN_TRIM(FILENAME)+1
	WRITE(FILENAME(I:I+2),'(I3.3)')DEPTH_INDX
        OPEN(UNIT=LU,FORM='UNFORMATTED',FILE=FILENAME,IOSTAT=IOS,
	1             ACCESS='SEQUENTIAL',STATUS='UNKNOWN',ACTION='WRITE')
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in WRITE_BCD_MAT'
	    WRITE(LUER,*)'Unable to open ',FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	  IF(TYPE .EQ. 'BC')THEN
	    WRITE(LU,IOSTAT=IOS)B,C,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ
	  ELSE IF(TYPE .EQ. 'C')THEN
	    WRITE(LU,IOSTAT=IOS)C,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ
	  ELSE IF(TYPE .EQ. 'D')THEN
	    WRITE(LU,IOSTAT=IOS)D
	  ELSE
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in READ_BCD_MAT'
	    WRITE(LUER,*)'Unrecognized TYPE: TYPE=',TYPE
	  END IF
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in WRITE_BCD_MAT'
	    WRITE(LUER,*)'Unable to write ',FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	CLOSE(LU)
!
	RETURN
	END
