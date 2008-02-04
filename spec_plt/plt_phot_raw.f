!
! Routine to plot OPACITY photoionization cross-sections. Data from 2 distinct
! routines may be plotted. Presently assumed level ordering is the same in
! both files. 
!
	PROGRAM PLT_PHOT_RAW
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 15-Jan-2008 : TYPE_1 etc changed to allocatable arrays.
!                       Errors returned if files cannot be opened.
! Created 09-Jun-1999
!
	INTEGER NPAIRS_1
	INTEGER NLEV_1
	INTEGER, ALLOCATABLE :: TYPE_1(:)
	INTEGER, ALLOCATABLE ::  NUM_VALS_1(:)
	INTEGER, ALLOCATABLE ::  LOC_1(:)
	CHARACTER(LEN=30), ALLOCATABLE :: NAME_1(:)
	REAL*8, ALLOCATABLE :: NU_1(:)
	REAL*8, ALLOCATABLE :: CROSS_1(:)
!
	INTEGER NPAIRS_2
	INTEGER NLEV_2
	INTEGER, ALLOCATABLE :: TYPE_2(:)
	INTEGER, ALLOCATABLE ::  NUM_VALS_2(:)
	INTEGER, ALLOCATABLE ::  LOC_2(:)
	CHARACTER(LEN=30), ALLOCATABLE :: NAME_2(:)
	REAL*8, ALLOCATABLE :: NU_2(:)
	REAL*8, ALLOCATABLE :: CROSS_2(:)
!
	CHARACTER*30 LEVEL_NAME1
	CHARACTER*30 LEVEL_NAME2
!
	INTEGER I,J,K
	INTEGER INDX_1,INDX_2,NV
	INTEGER NXT_LOC
	INTEGER IOS
	LOGICAL OKAY
	CHARACTER*132 STRING,FILENAME
!
	REAL*8 XV(10000),YV(10000)
!
	REAL*8 STAT_WT,GION,EDGE,TOTAL_REC,TEMP
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'Program designed to plot/compare opacity project data from photoionization files '
	WRITE(6,'(A)')'Use DISPGEN to plot no tabulated cross-sections '
	WRITE(6,'(A)')' '
!
10	FILENAME='PHOT1'
	CALL GEN_IN(FILENAME,'First photoionization file1:')
	NPAIRS_1=0
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    GOTO 10
	  END IF
	  DO WHILE(NPAIRS_1 .EQ. 0)
	    READ(10,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading file: ',TRIM(FILENAME)
	      WRITE(6,*)'IOSTAT=',IOS
	      STOP
	    END IF
	    IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NLEV_1
	      ALLOCATE (TYPE_1(NLEV_1),NUM_VALS_1(NLEV_1),LOC_1(NLEV_1), NAME_1(NLEV_1),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error in PLT_PHOT_RAW -- unable to allocate storage (1)'
	        WRITE(6,*)'STAT=',IOS
	        WRITE(6,*)'NLEV_1=',NLEV_1
	        STOP
	      END IF
	    END IF
	    IF( INDEX(STRING,'!Total number of data pairs') .NE. 0)THEN
	      READ(STRING,*)NPAIRS_1
	    END IF
	  END DO
!
	  ALLOCATE (NU_1(NPAIRS_1))
	  ALLOCATE (CROSS_1(NPAIRS_1))
!
	  NXT_LOC=1
	  DO J=1,NLEV_1
	    DO WHILE(INDEX(STRING,'!Configuration name') .EQ. 0)
	      READ(10,'(A)')STRING
	    END DO
	    K=INDEX(STRING,'  ')
	    NAME_1(J)=STRING(1:K-1)
	    READ(10,*)TYPE_1(J)
	    READ(10,*)NUM_VALS_1(J)
	    DO I=NXT_LOC,NXT_LOC+NUM_VALS_1(J)-1
	      IF(TYPE_1(J) .EQ. 20 .OR. TYPE_1(J) .EQ. 21)THEN
	        READ(10,*)NU_1(I),CROSS_1(I)
	      ELSE
	        READ(10,*)NU_1(I)
	      END IF
	    END DO
	    LOC_1(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_1(J)
	    STRING=' '
	  END DO
	CLOSE(UNIT=10)
!
!
!
20	NPAIRS_2=0
	FILENAME='PHOT2'
	CALL GEN_IN(FILENAME,'Second photoionization file1 ("" for null):')
	IF(FILENAME .EQ. ' ')GOTO 1000
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    GOTO 20
	  END IF
	  DO WHILE(NPAIRS_2 .EQ. 0)
	    READ(10,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading file: ',TRIM(FILENAME)
	      WRITE(6,*)'IOSTAT=',IOS
	      STOP
	    END IF
	    IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NLEV_2
	      ALLOCATE (TYPE_1(NLEV_2),NUM_VALS_1(NLEV_2),LOC_1(NLEV_2), NAME_1(NLEV_2),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error in PLT_PHOT_RAW -- unable to allocate storage (2)'
	        WRITE(6,*)'STAT=',IOS
	        WRITE(6,*)'NLEV_2=',NLEV_2
	        STOP
	      END IF
	    END IF
	    IF( INDEX(STRING,'!Total number of data pairs') .NE. 0)THEN
	      READ(STRING,*)NPAIRS_2
	    END IF
	  END DO
!
	  ALLOCATE (NU_2(NPAIRS_2))
	  ALLOCATE (CROSS_2(NPAIRS_2))
!
	  NXT_LOC=1
	  DO J=1,NLEV_2
	    DO WHILE(INDEX(STRING,'!Configuration name') .EQ. 0)
	      READ(10,'(A)')STRING
	    END DO
	    K=INDEX(STRING,'  ')
	    NAME_2(J)=STRING(1:K-1)
	    READ(10,*)TYPE_2(J)
	    READ(10,*)NUM_VALS_2(J)
	    DO I=NXT_LOC,NXT_LOC+NUM_VALS_2(J)-1
	      IF(TYPE_2(J) .EQ. 20 .OR. TYPE_2(J) .EQ. 21)THEN
	        READ(10,*)NU_2(I),CROSS_2(I)
	      ELSE
	        READ(10,*)NU_2(I)
	      END IF
	    END DO
	    LOC_2(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_2(J)
	    STRING=' '
	  END DO
	CLOSE(UNIT=10)
1000	CONTINUE
!
! 
!
	DO WHILE(1 .EQ. 1)
!
	  INDX_1=0
	  DO WHILE(INDX_1 .EQ. 0)
	    CALL GEN_IN(LEVEL_NAME1,'Level name for File 1')
	    IF(LEVEL_NAME1 .EQ. ' ')STOP
            DO J=1,NLEV_1
	      IF(LEVEL_NAME1 .EQ. NAME_1(J))THEN
                INDX_1=J
                EXIT
	      END IF
	    END DO
	    IF(INDX_1 .EQ. 0)WRITE(6,*)'Error - level name not found'
	  END DO
	  WRITE(6,*)'Level_1 name is: ',NAME_1(INDX_1)
!
	  NV=NUM_VALS_1(INDX_1)
	  XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	  YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	  CALL DP_CURVE(NV,XV,YV)
!
	  STAT_WT=1
	  GION=1 
	  EDGE=0
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'A non-zero ionization energy is only required to compute'
	  WRITE(6,'(A)')'the recombination rate.'
	  WRITE(6,'(A)')' '
	  CALL GEN_IN(EDGE,'Ionization energy (10^15 Hz)')
	  IF(EDGE .NE. 0)THEN
	    CALL GEN_IN(TEMP,'Temperature in 10^4K')
	    IF(TEMP .NE. 0)THEN
	      CALL RECOM_OPAC(YV,XV,EDGE,STAT_WT,GION,NV,NV,TOTAL_REC,TEMP)
	      WRITE(6,*)'Total Rec=',TOTAL_REC
	    END IF
	  END IF
!
	  IF(NPAIRS_2 .NE. 0)THEN
	    LEVEL_NAME2=LEVEL_NAME1
	    INDX_2=0
	    DO WHILE(INDX_2 .EQ. 0)
	      CALL GEN_IN(LEVEL_NAME2,'Level name for File 2')
              DO J=1,NLEV_2
	        IF(LEVEL_NAME2 .EQ. NAME_2(J))THEN
                  INDX_2=J
                  EXIT
	        END IF
	      END DO
	      IF(INDX_2 .EQ. 0)WRITE(6,*)'Error - level name not found'
	    END DO
	    NV=NUM_VALS_2(INDX_2)
	    XV(1:NV)=NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	    YV(1:NV)=CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	    CALL DP_CURVE(NV,XV,YV)
	    IF(EDGE .NE. 0)THEN
	      IF(TEMP .NE. 0)THEN
	        CALL RECOM_OPAC(YV,XV,EDGE,STAT_WT,GION,NV,NV,TOTAL_REC,TEMP)
	        WRITE(6,*)'Total Rec=',TOTAL_REC
	      END IF
	    END IF
	  END IF
	  CALL GRAMON_PGPLOT('\gn/\gn\do\u ','\gs(Mb)',LEVEL_NAME1,' ')
!
	  WRITE(6,*)'Now doing a LOG plot'
	  NV=NUM_VALS_1(INDX_1)
	  XV(1:NV)=LOG10(NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1))
	  YV(1:NV)=LOG10(CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)+1.0D-100)
	  CALL DP_CURVE(NV,XV,YV)
	  IF(NPAIRS_2 .NE. 0)THEN
	    NV=NUM_VALS_2(INDX_2)
	    XV(1:NV)=DLOG10(NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1))
	    YV(1:NV)=DLOG10(CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)+1.0D-100)
	    CALL DP_CURVE(NV,XV,YV)
	  END IF
	  CALL GRAMON_PGPLOT('Log \gn/\gn\do\u ','Log \gs(Mb)',LEVEL_NAME1,' ')
!
	  IF(EDGE .NE. 0)THEN
	    WRITE(6,*)'Now doing wavelength plot'
	    NV=NUM_VALS_1(INDX_1)
	    XV(1:NV)=0.299794D+04/(NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1))/EDGE
	    YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	    CALL DP_CURVE(NV,XV,YV)
	    IF(NPAIRS_2 .NE. 0)THEN
	      NV=NUM_VALS_2(INDX_2)
	      XV(1:NV)=0.299794D+04/(NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1))/EDGE
	      YV(1:NV)=CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	      CALL DP_CURVE(NV,XV,YV)
	    END IF
	    CALL GRAMON_PGPLOT('\gl(\gA)','\gs(Mb)',LEVEL_NAME1,' ')
	  END IF
	END DO
!
	END
