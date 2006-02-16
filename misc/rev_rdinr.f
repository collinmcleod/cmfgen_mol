!
! Simple program to modify a file in departure coeficient format, and modify
! for use with RDINR. IF THEN ELSE structure, so additional options can easily
! be included.
! Current options:
!                 FG_OB:       Finer grid near outer boundary.
!                 DOUB:        Doubel radius grid
!                 SCAL_R:      Doubel radius grid
!
	PROGRAM REV_RDINR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created : 23-Jan-2006
!
	INTEGER, PARAMETER :: MAX_ND=500
	REAL*8 R(MAX_ND)
	REAL*8 DI(MAX_ND)
	REAL*8 ED(MAX_ND)
	REAL*8 T(MAX_ND)
	REAL*8 IRAT(MAX_ND)
	REAL*8 VEL(MAX_ND)
	REAL*8 CLUMP_FAC(MAX_ND)
	REAL*8 RTMP(MAX_ND)
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER I,J
	INTEGER NI,NG
	INTEGER N,ND,NEW_ND
	REAL*8 RMIN,LUM
	REAL*8 GRID_FACTOR,GRID_RATIO
	REAL*8 T1,T2,SCALE_FACTOR,DELR
	CHARACTER*132 STRING
	CHARACTER*132 FILE_IN
	CHARACTER*132 FILE_OUT
!
	CHARACTER*20 OPTION
!
	FILE_IN='HIOUT'; FILE_OUT='RDINR'
	CALL GEN_IN(FILE_IN,'Input file: DC file format')
	CALL GEN_IN(FILE_OUT,'Output file: DC file format')
	OPTION='FG_OB'
	CALL GEN_IN(OPTION,'Action to be taken: FG_OB, DOUB, SCALE_R')
	
	OPEN(UNIT=9,FILE=FILE_IN,STATUS='OLD',ACTION='READ')
	OPEN(UNIT=10,FILE=FILE_OUT,STATUS='UNKNOWN',ACTION='WRITE')
	DO I=1,3
	  READ(9,'(A)')STRING
	  WRITE(10,'(A)')TRIM(STRING)
	END DO
	READ(9,'(A)')STRING
	READ(STRING,*)RMIN,LUM,N,ND
	READ(9,'(A)')STRING                 !Final blank line
!
	DO I=1,ND
	  READ(9,'(A)')STRING
	  READ(STRING,*)R(I),DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I)
	  DO WHILE(STRING .NE. ' ')
	    READ(9,'(A)',END=100)STRING
	  END DO
	END DO
100	CONTINUE
!
	CALL SET_CASE_UP(OPTION,IZERO,IZERO)
	IF(OPTION .EQ. 'DOUB')THEN
	  NEW_ND=2*ND-1
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,NEW_ND
	  WRITE(10,'(A)')' '
	  I=1
	  WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	  WRITE(10,'(F7.1)')1.0D0
!
	  DO I=2,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')SQRT(R(I-1)*R(I)),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	ELSE IF(OPTION .EQ. 'SCALE_R')THEN
	  T1=0.0D0
	  CALL GEN_IN(T1,'New RMIN')
	  IF(T1 .NE. 0.0D0)THEN
	    SCALE_FACTOR=T1/RMIN
	  ELSE
	    SCALE_FACTOR=1.0D0
	    CALL GEN_IN(SCALE_FACTOR,'Factor to scale RMIN by')
	  END IF
	  DO I=1,ND
	    R(I)=SCALE_FACTOR*R(I)
	  END DO
	  RMIN=RMIN*SCALE_FACTOR
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND
	  DO I=1,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
! Make a fine grid near the outer boundary.
!
	ELSE IF(OPTION .EQ. 'FG_OB')THEN
	  WRITE(6,*)'Current grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=1,12
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
	  NI=2; NG=3
	  CALL GEN_IN(NI,'Last depth to replace (e.g. 2):'); NI=NI-1 
	  CALL GEN_IN(NG,'Numer of additional points:')
	  CALL GEN_IN(GRID_RATIO,'Ratio of succesive interval sizes: > 1:')
	  GRID_RATIO=1.0D0/GRID_RATIO
	  WRITE(6,*)GRID_RATIO
!	  GRID_FACTOR=(GRID_RATIO-1)/GRID_RATIO
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND+NG
	  I=1
	  WRITE(10,'(A)')' '
	  WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	  WRITE(10,'(F7.1)')1.0D0
!
! Now do the insertion.
!
	  T1=R(NI+2)
	  WRITE(6,'(I4,ES18.8)')1,R(1)
	  DELR=R(NI+2)-R(NI+3)
	  DO I=NI+NG,1,-1
	    T2=(R(1)-T1)*(1.0D0-GRID_RATIO)/(1.0D0-GRID_RATIO**(I+1))
	    DELR=MIN(T2,1.2D0*DELR)
	    T1=T1+DELR
	    RTMP(I)=T1
	    WRITE(6,*)DELR,I,RTMP(I)
	  END DO
	  DO I=1,NI+NG
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')RTMP(I),
	1                DI(2),ED(2),T(2),IRAT(2),VEL(2),CLUMP_FAC(2),I+1
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
! Output remaining grid.
!
	  DO I=NI+2,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I+NG
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	  R(NI+NG+2:ND)=R(NI+2:ND)
	  R(2:NI+NG+1)=RTMP(1:NI+NG)
	  WRITE(6,*)' '
	  WRITE(6,*)'New grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=1,NI+NG+4
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
	ELSE
	  WRITE(6,*)TRIM(OPTION),' not recognized as valid option.'

	END IF
!
	STOP
	END
