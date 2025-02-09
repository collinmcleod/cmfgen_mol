!
! Subroutine to  read the output model data to a sequential,
! unformatted, file
!      JH_AT_OLD_TIME
! written as 8-byte floating point, and output a new file
!     NEW_JH_AT_OLD_TIME
! written as 10-byte floating point.
!
	PROGRAM REWRITE_JH
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 24-Oct-2023
!
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=10
!
	INTEGER NCF
	INTEGER ND
	INTEGER IREC
	INTEGER ST_IREC
	INTEGER REC_LENGTH

	REAL*8, ALLOCATABLE :: OLD_R(:)
	REAL*8, ALLOCATABLE :: OLD_V(:)
	REAL*8, ALLOCATABLE :: OLD_J(:)
	REAL*8, ALLOCATABLE :: OLD_H(:)
	REAL*8 OLD_NU
	REAL*8 OLD_H_INBC
	REAL*8 OLD_H_OUTBC
!
	REAL(KIND=LDP), ALLOCATABLE :: NEW_R(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_V(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_J(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_H(:)
	REAL(KIND=LDP) NEW_NU
	REAL(KIND=LDP) NEW_H_INBC
	REAL(KIND=LDP) NEW_H_OUTBC
!
	INTEGER I,ML,IOS
	CHARACTER*30 FILE_DATE
	CHARACTER*500 ERROR_MSG
!
        WRITE(6,*)' '
        WRITE(6,*)' This program expects a file JH_AT_OLD_TIME written as an 8-byte floating point number '
        WRITE(6,*)' This program will output NEW_JH_AT_OLD_TIME written as a 16-byte floating point number '
        WRITE(6,*)' '
        WRITE(6,*)' If successful, you will need to rename NEW_JH_AT_OLD_TIME to JH_AT_OLD_TIME'
        WRITE(6,*)' You should also edit JH_AT_OLD_TIME_INFO and change the word size from 8 to 16'
        WRITE(6,*)' '
!
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,'JH_AT_OLD_TIME',LU_IN,IOS)
	REC_LENGTH=REC_LENGTH/2
        WRITE(6,*)I,REC_LENGTH
	IF(IOS .NE. 0)THEN
          WRITE(6,*)'Error opening/reading JH_AT_OLD_TIME_INFO file: check format'
          STOP
	END IF
        OPEN(UNIT=LU_IN,FILE='JH_AT_OLD_TIME',STATUS='OLD',ACTION='READ',
	1        RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS,IOMSG=ERROR_MSG)
        IF(IOS .NE. 0)THEN
          WRITE(6,*)'Error opening JH_AT_OLD_TIME-- IOS=',IOS
          WRITE(6,*)TRIM(ERROR_MSG)
          STOP
	END IF
!
	REC_LENGTH=REC_LENGTH*16/8
        WRITE(6,*)I,REC_LENGTH
        OPEN(UNIT=LU_OUT,FILE='NEW_JH_AT_OLD_TIME',STATUS='UNKNOWN',ACTION='WRITE',
	1        RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS,IOMSG=ERROR_MSG)
        IF(IOS .NE. 0)THEN
          WRITE(6,*)'Error opening NEW_JH_AT_OLD_TIME-- IOS=',IOS
          WRITE(6,*)TRIM(ERROR_MSG)
          STOP
	END IF
	READ(LU_IN,REC=3)ST_IREC,NCF,ND
	WRITE(LU_OUT,REC=3)ST_IREC,NCF,ND
	ALLOCATE (OLD_V(ND),NEW_V(ND))
	ALLOCATE (OLD_R(ND),NEW_R(ND))
!
	IREC=ST_IREC
	READ(LU_IN,REC=IREC)OLD_R,OLD_V
	NEW_R=OLD_R; NEW_V=OLD_V
	WRITE(LU_OUT,REC=IREC)NEW_R,NEW_V
	IREC=IREC+1
!
	ALLOCATE (OLD_J(ND), NEW_J(ND))
	ALLOCATE (OLD_H(ND), NEW_H(ND))
!
! Read and write grey data
!
	IREC=ST_IREC+1
	READ(LU_IN,REC=IREC)(OLD_J(I),I=1,ND),(OLD_H(I),I=1,ND-1),OLD_H_INBC,OLD_H_OUTBC
	NEW_J=OLD_J; NEW_H=OLD_H
	NEW_H_INBC=OLD_H_INBC
	NEW_H_OUTBC=OLD_H_OUTBC
	WRITE(LU_OUT,REC=IREC)(NEW_J(I),I=1,ND),(NEW_H(I),I=1,ND-1),NEW_H_INBC,NEW_H_OUTBC
!
! We skip 2 records, because of the R,V & grey J & H output.
!
	IREC=ST_IREC+2
	WRITE(6,*)'Start record and IREC are ',ST_IREC,IREC
!
	DO ML=1,NCF
	  READ(LU_IN,REC=IREC)(OLD_J(I),I=1,ND),(OLD_H(I),I=1,ND-1),OLD_H_INBC,OLD_H_OUTBC,OLD_NU
	  NEW_J=OLD_J; NEW_H=OLD_H
	  NEW_H_INBC=OLD_H_INBC
	  NEW_H_OUTBC=OLD_H_OUTBC
	  NEW_NU=OLD_NU
	  WRITE(LU_OUT,REC=IREC)(NEW_J(I),I=1,ND),(NEW_H(I),I=1,ND-1),NEW_H_INBC,NEW_H_OUTBC,NEW_NU
	  IREC=IREC+1
	END DO
!
	STOP
	END
