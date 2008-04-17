!
! Routine to plot OPACITY photoionization cross-sections. Data from 2 distinct
! routines may be plotted. Presently assumed level ordering is the same in
! both files. 
!
	PROGRAM PLT_PHOT_RAW
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 16-Apr-2008 : Read in energy levels form oscilator file (if available).
!                       Deleted log plot section (since can be done in pgplot)
!                       Can now eneter level index for name. If unrecognized level,
!                          all names are dumped to terminal.
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
	REAL*8, ALLOCATABLE :: ENERGY_1(:)
	REAL*8, ALLOCATABLE :: STAT_WT_1(:)
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
! These are used to read the Oscilator file.
!
	INTEGER NELEV
	REAL*8, ALLOCATABLE :: ENERGY(:)
	REAL*8, ALLOCATABLE :: G(:)
	REAL*8 IONIZATION_ENERGY
	CHARACTER(LEN=30), ALLOCATABLE :: E_NAME(:)
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
	REAL*8 STAT_WEIGHT,GION,EDGE,TOTAL_REC,TEMP
	REAL*8 ANG_TO_HZ,SPEED_OF_LIGHT
	REAL*8, PARAMETER :: RONE=1.0D0
	LOGICAL DO_RECOM
	EXTERNAL SPEED_OF_LIGHT
!
	ANG_TO_HZ=1.0D-07*SPEED_OF_LIGHT()
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'Program designed to plot/compare opacity project data from photoionization files '
	WRITE(6,'(A)')'Use DISPGEN to plot tabulated cross-sections '
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
	      ALLOCATE (TYPE_1(NLEV_1),NUM_VALS_1(NLEV_1),LOC_1(NLEV_1),
	1               NAME_1(NLEV_1),ENERGY_1(NLEV_1),STAT_WT_1(NLEV_1),STAT=IOS)
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
	I=INDEX(FILENAME,'_A')
	IF(I .EQ.  0)I=INDEX(FILENAME,'_B')
	IF(I .GT. 6)THEN
	  FILENAME=FILENAME(5:I)//'F_OSCDAT'	
	ELSE
	  FILENAME=' '
	END IF
!
! Read in energy levels from oscilator file, if it is available.
! May need to change if oscilator file format changes.
!
14	CONTINUE
	CALL GEN_IN(FILENAME,'File with oscillator data: "" to type in indivdual level energies')
	IF(FILENAME .NE. " ")THEN
	  OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    GOTO 14
	  END IF
	  STRING=' '
!	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0)
!	    READ(10,'(A)')STRING
!	  END DO
!	  IF(INDEX(STRING,'17-Oct-2000') .EQ. 0)THEN
!	    WRITE(6,*)'Oscilator format date is not recognized by this routine'
!	    GOTO 14
!	  END IF
	  DO WHILE(INDEX(STRING,'!Number of energy levels') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)NELEV
	  DO WHILE(INDEX(STRING,'!Ionization energy') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)IONIZATION_ENERGY
	  DO WHILE(INDEX(STRING,'!Number of transition') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING		!Get blank line
	  ALLOCATE(G(NELEV),E_NAME(NELEV),ENERGY(NELEV))
!
	  G(:)=0.0D0;ENERGY(:)=0.0D0
	  DO I=1,NELEV
	    READ(10,'(A)')STRING
	    K=INDEX(STRING,'  ')
	    E_NAME(I)=STRING(1:K-1)
	    READ(STRING(K:),*)G(I),ENERGY(I)
	    IF(E_NAME(I)(K-1:K-1) .EQ.  ']')THEN
	      K=INDEX(E_NAME(I),'[')
	      E_NAME(I)(K:)=' '
	    END IF
	    WRITE(6,*)G(I),ENERGY(I)
	  END DO
	  CLOSE(UNIT=10)
	  PAUSE
!
! Now need to match names. We assume photoionzation files do not
! have [].
!
	  ENERGY_1(1:NLEV_1)=0.0D0; STAT_WT_1(1:NLEV_1)=0.0D0
	  DO I=1,NLEV_1
	    DO J=1,NELEV
	      IF(E_NAME(J) .EQ. NAME_1(I))THEN
	        ENERGY_1(I)=ENERGY_1(I)+G(J)*ENERGY(J)
	        STAT_WT_1(I)=STAT_WT_1(I)+G(J)
	      END IF
	    END DO
	    IF(STAT_WT_1(I) .NE. 0.0D0)ENERGY_1(I)=1.0D-15*SPEED_OF_LIGHT()*
	1        (IONIZATION_ENERGY-ENERGY_1(I)/STAT_WT_1(I))
	  END DO
	ELSE
	  ENERGY_1(1:NLEV_1)=0.0D0; STAT_WT_1(1:NLEV_1)=0.0D0
	END IF
	WRITE(6,*)ENERGY_1(1),STAT_WT_1(1)
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
! We can now use a level name, or an index (sas stored in PHOT file).
!
	  INDX_1=0
	  LEVEL_NAME1='1'
	  DO WHILE(INDX_1 .EQ. 0)
	    CALL GEN_IN(LEVEL_NAME1,'Level name [or index] for File 1')
	    IF(LEVEL_NAME1 .EQ. ' ')STOP
	    DO J=1,NLEV_1
	      IF(LEVEL_NAME1 .EQ. NAME_1(J))THEN
                INDX_1=J
                EXIT
	      END IF
	    END DO
	    IF(INDX_1 .EQ. 0)THEN
	      READ(LEVEL_NAME1,*,IOSTAT=IOS)INDX_1
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error - level name not found'
	        DO I=1,NLEV_1,5
	          WRITE(6,'(5A14)')(TRIM(NAME_1(J)),J=I,MAX(I+4,NLEV_1))
	        END DO
	        INDX_1=0
	      END IF
	    END IF
	  END DO
	  WRITE(6,*)'Level_1 name is: ',NAME_1(INDX_1)
!
	  NV=NUM_VALS_1(INDX_1)
	  XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	  YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	  CALL DP_CURVE(NV,XV,YV)
!
	  STAT_WEIGHT=MAX(RONE,STAT_WT_1(INDX_1))
	  EDGE=ENERGY_1(INDX_1)
	  GION=1 
	  IF(EDGE .EQ. 0.0D0)THEN
	    WRITE(6,'(A)')' '
	    WRITE(6,'(A)')'A non-zero ionization energy is only required to compute'
	    WRITE(6,'(A)')'the recombination rate.'
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(EDGE,'Ionization energy (10^15 Hz)')
	    DO_RECOM=.TRUE.
	    IF(EDGE .EQ. 0)DO_RECOM=.FALSE.
	  ELSE
	    WRITE(6,*)'Energy of level is',EDGE
	    DO_RECOM=.FALSE.
	    CALL GEN_IN(DO_RECOM,'Compute recombination rate?')
	  END IF
	  IF(DO_RECOM)THEN
	    CALL GEN_IN(TEMP,'Temperature in 10^4K')
	    IF(TEMP .NE. 0)THEN
	      CALL RECOM_OPAC(YV,XV,EDGE,STAT_WEIGHT,GION,NV,NV,TOTAL_REC,TEMP)
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
	    IF(DO_RECOM)THEN
	      IF(TEMP .NE. 0)THEN
	        CALL RECOM_OPAC(YV,XV,EDGE,STAT_WEIGHT,GION,NV,NV,TOTAL_REC,TEMP)
	        WRITE(6,*)'Total Rec=',TOTAL_REC
	      END IF
	    END IF
	  END IF
	  CALL GRAMON_PGPLOT('\gn/\gn\do\u ','\gs(Mb)',LEVEL_NAME1,' ')
!
	  IF(EDGE .NE. 0.0D0 .AND. XV(1) .NE. 0.0D0)THEN
	    WRITE(6,*)'Now doing wavelength plot'
	    NV=NUM_VALS_1(INDX_1)
	    XV(1:NV)=ANG_TO_HZ/(NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1))/EDGE
	    YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	    CALL DP_CURVE(NV,XV,YV)
	    IF(NPAIRS_2 .NE. 0)THEN
	      NV=NUM_VALS_2(INDX_2)
	      XV(1:NV)=ANG_TO_HZ/(NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1))/EDGE
	      YV(1:NV)=CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	      CALL DP_CURVE(NV,XV,YV)
	    END IF
	    CALL GRAMON_PGPLOT('\gl(\gA)','\gs(Mb)',LEVEL_NAME1,' ')
	  END IF
	END DO
!
	END
