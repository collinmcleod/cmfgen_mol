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
	REAL*8 ZION_1
	REAL*8 GION_1 
	REAL*8 EXC_EN_1
	REAL*8 AMASS
!
	INTEGER NPAIRS_2
	INTEGER NLEV_2
	INTEGER, ALLOCATABLE :: TYPE_2(:)
	INTEGER, ALLOCATABLE ::  NUM_VALS_2(:)
	INTEGER, ALLOCATABLE ::  LOC_2(:)
	CHARACTER(LEN=30), ALLOCATABLE :: NAME_2(:)
	REAL*8, ALLOCATABLE :: NU_2(:)
	REAL*8, ALLOCATABLE :: CROSS_2(:)
	REAL*8 ZION_2
	REAL*8 GION_2 
	REAL*8 EXC_EN_2
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
	REAL*8 T1
	REAL*8 FREQ_SCL_FAC
!
	INTEGER I,J,K
	INTEGER ICOUNT
	INTEGER NT
	INTEGER INDX_1,INDX_2,NV
	INTEGER NXT_LOC
	INTEGER IOS
	LOGICAL OKAY
	LOGICAL CREATE_SUMMARY
	CHARACTER*132 STRING,FILENAME
!
	REAL*8, PARAMETER :: NCF_MAX=10000
	REAL*8 XV(NCF_MAX),YV(NCF_MAX)
!
	INTEGER, PARAMETER :: NREC_MAX=5
	REAL*8 TEMP_VEC(NREC_MAX)
	REAL*8 TOTAL_REC_VEC(NREC_MAX)
	REAL*8 LEVEL_REC_VEC(NREC_MAX)
!
	REAL*8 STAT_WEIGHT,GION,EDGE
	REAL*8 LEVEL_REC,TOTAL_REC,TEMP
	REAL*8 ANG_TO_HZ,SPEED_OF_LIGHT
	REAL*8, PARAMETER :: RONE=1.0D0
	LOGICAL DO_WAVE_PLT
	LOGICAL DO_RECOM
	LOGICAL DO_ALL_RECOM
	LOGICAL DO_SEQ_PLTS
	LOGICAL OSCILLATOR_FILE_AVAIL
	EXTERNAL SPEED_OF_LIGHT
!
	CHARACTER(LEN=30) UC; EXTERNAL UC
	CHARACTER(LEN=30) XLAB
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER, PARAMETER :: LUER=6
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUSCR=8
	INTEGER, PARAMETER :: LUOUT=12
!
! Set constants.
!
        CHIBF=2.815D-06
        CHIFF=3.69D-29
        HDKT=4.7994145D0
        TWOHCSQ=0.0147452575D0
        OPLIN=2.6540081D+08
        EMLIN=5.27296D-03
	ANG_TO_HZ=1.0D-07*SPEED_OF_LIGHT()
	AMASS=40.0D0
!
! Read in bound-free gaunt factors for individual n states of hydrogen,
! and hydrogenic cross-sections for individual l states (n =0 to 30,
! l=0 to n-1)
!
        CALL RD_HYD_BF_DATA(LUIN,LUSCR,LUER)
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'Program designed to plot/compare opacity project data from photoionization files '
	WRITE(6,'(A)')'Use DISPGEN to plot tabulated cross-sections '
	WRITE(6,'(A)')' '
!
10	FILENAME='PHOT1'
	CALL GEN_IN(FILENAME,'First photoionization file1:')
	NPAIRS_1=0
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
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
	    IF( INDEX(STRING,'!Screened nuclear charge') .NE. 0)THEN
	      READ(STRING,*)ZION_1
	    ELSE IF( INDEX(STRING,'Statistical weight of ion') .NE. 0)THEN
	      READ(STRING,*)GION_1
	    ELSE IF( INDEX(STRING,'!Excitation energy of final state') .NE. 0)THEN
	      READ(STRING,*)EXC_EN_1
	    ELSE IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
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
	    WRITE(16,*)TRIM(NAME_1(J)),TYPE_1(J),NUM_VALS_1(J)
	    IF(TYPE_1(J) .EQ. 20 .OR. TYPE_1(J) .EQ. 21)THEN
	      READ(10,*)(NU_1(I),CROSS_1(I), I=NXT_LOC,NXT_LOC+NUM_VALS_1(J)-1)
	    ELSE
	      READ(10,*)(CROSS_1(I), I=NXT_LOC,NXT_LOC+NUM_VALS_1(J)-1)
	    END IF
	    LOC_1(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_1(J)
	    STRING=' '
	  END DO
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
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Not all options available if no oscillator file '
	WRITE(6,'(A)')' Can only treat opacity data when no oscillator file'
	WRITE(6,'(A)')' '
!
	CALL GEN_IN(FILENAME,'File with oscillator data: "" for no file')
	IF(FILENAME .NE. " ")THEN
	  OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
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
	  DO WHILE(INDEX(STRING,'!Number of transitions') .EQ. 0)
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
	  END DO
	  CLOSE(UNIT=10)
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
	  OSCILLATOR_FILE_AVAIL=.TRUE.
	ELSE
	  ENERGY_1(1:NLEV_1)=0.0D0; STAT_WT_1(1:NLEV_1)=0.0D0
	  OSCILLATOR_FILE_AVAIL=.FALSE.
	END IF
!
!
!
20	NPAIRS_2=0
	FILENAME='PHOT2'
	CALL GEN_IN(FILENAME,'Second photoionization file1 ("" for null):')
	IF(FILENAME .EQ. ' ')GOTO 1000
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
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
	    IF( INDEX(STRING,'!Screened nuclear charge') .NE. 0)THEN
	      READ(STRING,*)ZION_2
	    ELSE IF( INDEX(STRING,'Statistical weight of ion') .NE. 0)THEN
	      READ(STRING,*)GION_2
	    ELSE IF( INDEX(STRING,'!Excitation energy of final state') .NE. 0)THEN
	      READ(STRING,*)EXC_EN_2
	    ELSE IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NLEV_2
	      ALLOCATE (TYPE_2(NLEV_2),NUM_VALS_2(NLEV_2),LOC_2(NLEV_2), NAME_2(NLEV_2),STAT=IOS)
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
	    IF(TYPE_2(J) .EQ. 20 .OR. TYPE_2(J) .EQ. 21)THEN
	      READ(10,*)(NU_2(I),CROSS_2(I), I=NXT_LOC,NXT_LOC+NUM_VALS_2(J)-1)
	    ELSE
	      READ(10,*)(CROSS_2(I), I=NXT_LOC,NXT_LOC+NUM_VALS_2(J)-1)
	    END IF
	    LOC_2(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_2(J)
	    STRING=' '
	  END DO
	CLOSE(UNIT=10)
1000	CONTINUE
!
! 
!
	CREATE_SUMMARY=.FALSE.
	CALL GEN_IN(CREATE_SUMMARY,'Create a summary of photo-ionization cross sections (1st file only)?')
	IF(CREATE_SUMMARY)THEN
	  OPEN(UNIT=11,FILE='Phot_summary',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(11,'(A,T20,A,6X,A,3(10X,A))')'Level','Type',' Np','X1','X2','X2'
	    DO K=1,40
	      DO J=1,NLEV_1
	        IF(TYPE_1(J) .EQ. K)THEN
	          WRITE(11,'(A,T20,I4,2X,I6,5ES13.4)')TRIM(NAME_1(J)),TYPE_1(J),NUM_VALS_1(J),
	1            CROSS_1(LOC_1(J):LOC_1(J)+MIN(4,NUM_VALS_1(J)-1))
	        END IF
	      END DO
	    END DO
	  CLOSE(UNIT=11)
	END IF
!
	DO_ALL_RECOM=.FALSE.
	IF(OSCILLATOR_FILE_AVAIL)CALL GEN_IN(DO_ALL_RECOM,'Compute recombination rates for all species?')
	IF(DO_ALL_RECOM)THEN
	  OPEN(UNIT=LUOUT,FILE='RECOM_SUM',STATUS='UNKNOWN',ACTION='WRITE')
	  TEMP_VEC(1)=0.5D0; TEMP_VEC(2)=1.0D0; TEMP_VEC(3)=2.0; TEMP_VEC(4)=5.0D0; TEMP_VEC(5)=10.0D0
	  TOTAL_REC_VEC(:)=0.0D0
	  CALL GEN_IN(TEMP_VEC,NT,NREC_MAX,'Temperature in 10^4K (5 values max)')
	  WRITE(6,'(A,T30,5(5X,F6.2)')'Temperature (10^4 K)=',(TEMP_VEC(I),I=1,NT)
	  WRITE(LUOUT,'(A,T30,5(5X,F6.2)')'Level / Temperature (10^4 K)',(TEMP_VEC(I),I=1,NT)
	  DO INDX_1=1,NLEV_1
	    EDGE=ENERGY_1(INDX_1)
	    STAT_WEIGHT=STAT_WT_1(INDX_1)
	    IF(TYPE_1(INDX_1) .EQ. 20 .OR. TYPE_1(INDX_1) .EQ. 21)THEN
	      NV=NUM_VALS_1(INDX_1)
	      XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	      YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	      FREQ_SCL_FAC=EDGE+EXC_EN_1
	    ELSE
	      NV=1000
	      CALL RAW_SUBPHOT(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NV)
	      FREQ_SCL_FAC=XV(1)
	      XV(1:NV)=XV(1:NV)/FREQ_SCL_FAC
	    END IF
	    T1=EDGE+EXC_EN_1
	    WRITE(6,*)NAME_1(INDX_1)
	    DO I=1,NT
	      CALL RECOM_OPAC_V2(YV,XV,T1,FREQ_SCL_FAC,STAT_WEIGHT,GION_1,NV,NV,LEVEL_REC_VEC(I),TEMP_VEC(I))
	      TOTAL_REC_VEC(I)=TOTAL_REC_VEC(I)+LEVEL_REC_VEC(I)
	    END DO
	    WRITE(LUOUT,'(A,T30,5ES11.3)')TRIM(NAME_1(INDX_1)),(LEVEL_REC_VEC(I),I=1,NT)
	  END DO
	  WRITE(6,'(A,T30,5ES11.3)')'Total Recom. Rate=',(TOTAL_REC_VEC(I),I=1,NT)
	  WRITE(LUOUT,'(A,T30,5ES11.3)')'Total Recom. Rate=',(TOTAL_REC_VEC(I),I=1,NT)
	  CLOSE(LUOUT)
	END IF
!
	DO_WAVE_PLT=.FALSE.
	CALL GEN_IN(DO_WAVE_PLT,'Plot versus wavelength?')
	DO_RECOM=.FALSE.
	CALL GEN_IN(DO_RECOM,'Compute recombination rate?')
	LEVEL_NAME1=' '
!
	DO_SEQ_PLTS=.FALSE.
	CALL GEN_IN(DO_SEQ_PLTS,'Plot photoionization cross-section for a sequence: eg., 6d (T or F)')
	IF(DO_SEQ_PLTS)THEN
	  DO WHILE(1 .EQ. 1)
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(LEVEL_NAME1,'Sequence desciptor')
	    IF(LEVEL_NAME1 .EQ. ' ')STOP
	    WRITE(6,'(A)')' '
	    ICOUNT=0
	    DO J=1,NLEV_1
	      IF(INDEX(NAME_1(J),TRIM(LEVEL_NAME1)) .NE. 0)THEN
                ICOUNT=ICOUNT+1
	        INDX_1=J
	        EDGE=ENERGY_1(INDX_1)
	        WRITE(6,'(4X,A,A,T40,ES10.4)')'Level name/energy is: ',TRIM(NAME_1(INDX_1)),EDGE
	        IF(TYPE_1(INDX_1) .EQ. 20 .OR. TYPE_1(INDX_1) .EQ. 21)THEN
	          NV=NUM_VALS_1(INDX_1)
	          XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	          YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	        ELSE
	          NV=1000
	          CALL RAW_SUBPHOT(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NV)
	          XV(1:NV)=XV(1:NV)/(EDGE+EXC_EN_1)
	        END IF
	        IF(DO_WAVE_PLT)THEN
	          XV(1:NV)=ANG_TO_HZ/XV(1:NV)
	        END IF
	        CALL DP_CURVE(NV,XV,YV)
	      END IF
	      IF(ICOUNT .EQ. 50)THEN
	        WRITE(6,*)'Maximum number of plots exceeded'
	        EXIT
	      END IF
	    END DO
	    WRITE(6,'(A)')' '
	   CALL GRAMON_PGPLOT('\gn/\gn\do\u ','\gs(Mb)',TRIM(LEVEL_NAME1),' ')
	  END DO
	END IF
!
	XLAB='\gn/\gn\do\u'
	IF(DO_WAVE_PLT)XLAB='\gl(\A)'
	DO WHILE(1 .EQ. 1)
!
! We can now use a level name, or an index (as stored in PHOT file).
!
	  INDX_1=0
	  LEVEL_NAME1='1'
	  DO WHILE(INDX_1 .EQ. 0)
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(LEVEL_NAME1,'Level name [or index] for File 1 (P to plot, E to exit)')
	    IF(UC(LEVEL_NAME1) .EQ. 'P')THEN
	      CALL GRAMON_PGPLOT(XLAB,'\gs(Mb)',LEVEL_NAME1,' ')
	      CALL GEN_IN(LEVEL_NAME1,'Level name [or index] for File (E to exit)')
	    END IF
	    IF( UC(LEVEL_NAME1) .EQ. 'E' .OR. UC(LEVEL_NAME1(1:2)) .EQ. 'EX' 
	1         .OR. UC(LEVEL_NAME1) .EQ. ' ')STOP
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
	      ELSE
	         LEVEL_NAME1=NAME_1(INDX_1)
	      END IF
	    END IF
	  END DO
	  EDGE=ENERGY_1(INDX_1)
	  WRITE(6,'(/,8X,A,A)')     '          Level_1 name is: ',NAME_1(INDX_1)
	  WRITE(6,'(8X,A,ES15.8,A,')'       Energy of level is: ',EDGE,' (10^15Hz)'
	  WRITE(6,'(8X,A,I2,/)')    ' Type of cross-section is: ',TYPE_1(INDX_1)
!
	  IF(TYPE_1(INDX_1) .EQ. 20 .OR. TYPE_1(INDX_1) .EQ. 21)THEN
	    NV=NUM_VALS_1(INDX_1)
	    XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	    YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	  ELSE
	    NV=1000
	    CALL RAW_SUBPHOT(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NV)
	    XV(1:NV)=XV(1:NV)/(EDGE+EXC_EN_1)
	  END IF
	  IF(DO_WAVE_PLT)THEN
	    XV(1:NV)=ANG_TO_HZ/XV(1:NV)
	  END IF
	  CALL DP_CURVE(NV,XV,YV)
!
	  STAT_WEIGHT=MAX(RONE,STAT_WT_1(INDX_1))
	  GION=1 
	  IF(EDGE .EQ. 0.0D0)THEN
	    WRITE(6,'(A)')' '
	    WRITE(6,'(A)')'A non-zero ionization energy is only required to compute'
	    WRITE(6,'(A)')'the recombination rate.'
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(EDGE,'Ionization energy (10^15 Hz)')
	    DO_RECOM=.TRUE.
	    IF(EDGE .EQ. 0)DO_RECOM=.FALSE.
	  END IF
	  IF(DO_RECOM .AND. EDGE .NE. 0)THEN
	    TEMP=1.0D0
	    CALL GEN_IN(TEMP,'Temperature in 10^4K')
	    IF(TEMP .NE. 0)THEN
	      WRITE(6,'(A,ES10.4,3X,A,F5.1,3X,A,F4.1,3X,A,I6,3X,A,F6.2)')
	1        ' EDGE=',EDGE,'g=',STAT_WEIGHT,'gion=',GION,'NV=',NV,'T(10^K)=',TEMP
	      CALL RECOM_OPAC(YV,XV,EDGE,STAT_WEIGHT,GION,NV,NV,TOTAL_REC,TEMP)
	      WRITE(6,'(A,ES11.4)')' Total Rec=',TOTAL_REC
	    END IF
	  END IF
!
	  IF(NPAIRS_2 .NE. 0)THEN
	    LEVEL_NAME2=NAME_1(INDX_1)
	    INDX_2=0
	    DO WHILE(INDX_2 .EQ. 0)
	      CALL GEN_IN(LEVEL_NAME2,'Level name for File 2')
              DO J=1,NLEV_2
	        IF(LEVEL_NAME2 .EQ. NAME_2(J))THEN
                  INDX_2=J
                  EXIT
	        END IF
	      END DO
	      IF(INDX_2 .EQ. 0)THEN
	        READ(LEVEL_NAME2,*,IOSTAT=IOS)INDX_2
	        IF(IOS .NE. 0)THEN
	          WRITE(6,*)'Error - level name not found'
	          DO I=1,NLEV_2,5
	            WRITE(6,'(5A14)')(TRIM(NAME_2(J)),J=I,MAX(I+4,NLEV_2))
	          END DO
	          INDX_2=0
	        ELSE
	          LEVEL_NAME2=NAME_2(INDX_2)
	        END IF
	      END IF
	    END DO
!
	    EDGE=ENERGY_1(INDX_1)
	    IF(INDX_2 .GT. 0)THEN
	      IF(TYPE_2(INDX_1) .EQ. 20 .OR. TYPE_2(INDX_1) .EQ. 21)THEN
	        NV=NUM_VALS_2(INDX_2)
	        XV(1:NV)=NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	        YV(1:NV)=CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	      ELSE
	        NV=1000
	        CALL RAW_SUBPHOT(YV,XV,CROSS_2(LOC_2(INDX_2)),TYPE_2(INDX_2),NUM_VALS_2(INDX_2),
	1                    EDGE,EXC_EN_2,ZION_2,AMASS,NV)
	         XV(1:NV)=XV(1:NV)/(EDGE+EXC_EN_2)
	      END IF
	      IF(DO_RECOM .AND. OSCILLATOR_FILE_AVAIL)THEN
	        IF(TEMP .NE. 0)THEN
	          CALL RECOM_OPAC(YV,XV,EDGE,STAT_WEIGHT,GION,NV,NV,TOTAL_REC,TEMP)
	          WRITE(6,*)'Total Rec=',TOTAL_REC
	        END IF
	      END IF
	    END IF
	    IF(DO_WAVE_PLT)THEN
	      XV(1:NV)=ANG_TO_HZ/XV(1:NV)
	    END IF
	    CALL DP_CURVE(NV,XV,YV)
	  END IF
!
	END DO		!Loop over levels
!
	END
