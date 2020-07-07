!
! Designed to plot/examine NETRATE or TOTRATE.
! Program requires:
!       NETRATE or TOTRATE
!       MODEL
!       MODEL_SPEC
!       RVTJ
!       XzVPRRR   (XzV is the species being investigated)  
!
! Program will run quicker (since less to read in) if you use a secondary
! file create by:
!                grep XzV NETRATE >> SMALL_FILE
!
	PROGRAM PLT_RATES
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created 24-May-2020
!
	INTEGER, PARAMETER :: NMAX=600000
!
	REAL*8, ALLOCATABLE :: NT_NU(:)
	REAL*8, ALLOCATABLE :: NT_LAM(:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: LAM(:)
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	CHARACTER(LEN=60), ALLOCATABLE :: TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: NT_TRANS_NAME(:)
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: ED(:)
!
	CHARACTER(LEN=10) PR_SPECIES
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=30) LEVEL
!
	INTEGER*4, ALLOCATABLE ::  LINK(:)
	REAL*8, ALLOCATABLE :: RATES(:,:)
	REAL*8, ALLOCATABLE :: NEW_RATES(:,:)
!
	INTEGER*4, ALLOCATABLE ::  NT_LINK(:)
	REAL*8, ALLOCATABLE :: NT_RATES(:,:)
	REAL*8, ALLOCATABLE :: NEW_NT_RATES(:,:)
!
	REAL*8, ALLOCATABLE :: PHOT_RATE(:,:)
	REAL*8, ALLOCATABLE :: REC_RATE(:,:)
!
	REAL*8 T1,T2
	REAL*8 SUM_TO,SUM_FROM
	REAL*8 CUR_SUM
	INTEGER ND
	INTEGER N_LINES
	INTEGER N_TRANS
	INTEGER N_NT_LINES
	INTEGER N_NT_TRANS
	INTEGER RD_COUNT
	INTEGER I,J,K,L,ML,IBEG,IDEPTH
	INTEGER MAX_TRANS_LENGTH
	INTEGER DPTH_INDX
	INTEGER NL,NUP
	INTEGER IOS
	INTEGER NLEV
	INTEGER LU
	LOGICAL FILE_OPEN
	LOGICAL DO_NT_RATES
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=200) STRING
        CHARACTER(LEN=30) UC
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) XLABEL
	CHARACTER(LEN=20) YLABEL
        EXTERNAL UC
!
	NLEV=0
	SPECIES=' '
	PR_SPECIES=' '
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot/examin NETRATE or TOTRATE'
	WRITE(6,'(A)')' '
!
100	CONTINUE
 	FILENAME='TOTRATE'
	CALL GEN_IN(FILENAME,'File with data to be plotted/examined')
	IF(FILENAME .EQ. ' ')GOTO 100
!
! Get number of depth points
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of data points (must be exact)')
	END IF
!
	N_LINES=NMAX
	CALL GEN_IN(N_LINES,'Maximum number of lines to be read')
!
! These are updated when the data file is read in.
!
	ALLOCATE (RATES(ND,N_LINES))
	ALLOCATE (TRANS_NAME(N_LINES))
	ALLOCATE (NU(N_LINES))
	ALLOCATE (LAM(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
! 
	ALLOCATE (NEW_RATES(ND,N_LINES))
	ALLOCATE (LINK(N_LINES))
!
! Used for ploting.
!
	ALLOCATE (XV(ND))
	ALLOCATE (YV(ND))
!
! These are updated when the data file is read in.
!
	DO_NT_RATES=.FALSE.
	CALL GEN_IN(DO_NT_RATES,'Read in non thermal rates')
	IF(DO_NT_RATES)THEN
	  ALLOCATE (NT_RATES(ND,N_LINES))
	  ALLOCATE (NT_TRANS_NAME(N_LINES))
	  ALLOCATE (NT_NU(N_LINES))
	  ALLOCATE (NT_LAM(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
! 
	  ALLOCATE (NEW_NT_RATES(ND,N_LINES))
	  ALLOCATE (NT_LINK(N_LINES))
	END IF
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	RD_COUNT=0
	CALL TUNE(1,'READ')
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	    READ(11,'(A)',END=5000)STRING
	  END DO
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')
	  STRING=STRING(K+2:)
	  K=INDEX(STRING,'  ')
	  TRANS_NAME(ML)=STRING(1:K)
	  READ(STRING(K+1:),*)NU(ML)
	  RD_COUNT=ML
	  READ(11,*)(RATES(I,ML),I=1,ND)
	END DO
!
5000	CONTINUE
	CALL TUNE(2,'READ')
	CALL TUNE(3,'TERMINAL')

	N_LINES=RD_COUNT
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Number of lines read is',N_LINES 
	WRITE(6,'(A)')DEF_PEN
!
        LAM(1:N_LINES)=2.99792458D+03/NU(1:N_LINES)
	LU=10
	FILENAME='RVTJ'
	ALLOCATE (R(ND),V(ND),ED(ND))
	CALL RD_SING_VEC_RVTJ(R,ND,'Radius',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(V,ND,'Velocity',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(ED,ND,'Electron',FILENAME,LU,IOS)
!
	OPEN(UNIT=11,FILE='NT_RATES',STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	RD_COUNT=0
	CALL TUNE(1,'READ_NT')
	DO ML=1,N_LINES             !NT_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	    READ(11,'(A)',END=6000)STRING
	  END DO
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')
	  NT_TRANS_NAME(ML)=STRING(1:K)
	  READ(STRING(K+1:),*)NT_NU(ML)
	  RD_COUNT=ML
	  READ(11,*)(NT_RATES(I,ML),I=1,ND)
	END DO
6000	CONTINUE
	N_NT_LINES=RD_COUNT
	CALL TUNE(2,'READ_NT')
	CALL TUNE(3,'TERMINAL')

	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Number of number of nonthermal lines read is',N_NT_LINES 
	WRITE(6,'(A)')DEF_PEN
!
        NT_LAM(1:N_LINES)=2.99792458D+03/NT_NU(1:N_NT_LINES)
!
! Set def axis.
!

	XV(1:ND)=V(1:ND)
	XLABEL='V(km/s)'
!
! Enter main plotting/examination loop.
!
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  PLT_OPT='P'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'LOGR')THEN
	    T1=R(ND)
	    DO I=1,ND
	      XV(I)=LOG10(R(I)/T1)
	    END DO
	    XLABEL='Log R/R(ND)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'V')THEN
	    XV(1:ND)=V(1:ND)
	    XLABEL='V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LOGV')THEN
	    XV(1:ND)=LOG10(V(1:ND))
	    XLABEL='Log V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XN')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    XLABEL='Depth index'
!
! Set the species, and reads in the SPECIES//'PRRR' file.
!
	  ELSE IF(UC(PLT_OPT(1:4)) .EQ. 'SPEC')THEN
	    CALL GEN_IN(SPECIES,'Species that will be examined')
!
! Get number of super levels. This is need when we read in the the
! recombination and photoioizations rates for each super level.
!
	    FILENAME='MODEL_SPEC'
	    OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	        GOTO 1000
	      END IF
	      DO WHILE(INDEX(STRING,TRIM(SPECIES)) .EQ. 0)
	        READ(11,'(A)')STRING
	      END DO
	      READ(STRING,*)I,NLEV
	    CLOSE(UNIT=11)
	    IF(ALLOCATED(PHOT_RATE))DEALLOCATE(PHOT_RATE,REC_RATE)
	    ALLOCATE (PHOT_RATE(NLEV,ND),REC_RATE(NLEV,ND))
!
! Read in the recombination and photoioizations rates for each super level.
! Note that ND is the second index.
!
	    FILENAME=TRIM(SPECIES)//'PRRR'
	    OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	      WRITE(6,*)'Was you species name correct: Species=',TRIM(SPECIES)
	      GOTO 1000
	    END IF
!
	    DO K=1,ND,10
	      DO WHILE(INDEX(STRING,'Photoionization') .EQ. 0)
	        READ(11,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(11,*)(PHOT_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	      DO WHILE(INDEX(STRING,'Recombination') .EQ. 0)
	        READ(11,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(11,*)(REC_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	    END DO
	    CLOSE(UNIT=11)
	    PR_SPECIES=SPECIES 
!
! Set the level for which rates will be examined, and store the
! data in a new array.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LEVEL')THEN
	    N_TRANS=0
	    IF(SPECIES .EQ. ' ')THEN
	       WRITE(6,*)'Need to use SPECIES option first to set the species'
	       GOTO 1000
	    END IF
	    CALL GEN_IN(LEVEL,'Level -- full name -- case sensitive')
	    MAX_TRANS_LENGTH=0
	    DO ML=1,N_LINES
	       IF(INDEX(TRANS_NAME(ML),TRIM(SPECIES)) .NE. 0 .AND.
	1         ( INDEX(TRANS_NAME(ML),'('//TRIM(LEVEL)) .NE. 0 .OR.
	1           INDEX(TRANS_NAME(ML),'-'//TRIM(LEVEL)) .NE. 0) )THEN
	          N_TRANS=N_TRANS+1
	          NEW_RATES(:,N_TRANS)=RATES(:,ML)
	          LINK(N_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(TRANS_NAME(ML)))
	       END IF
	    END DO
	    WRITE(6,*)'Number of transitions identified is',N_TRANS
!
	    N_NT_TRANS=0
	    IF(DO_NT_RATES)THEN
	      DO ML=1,N_NT_LINES
	        IF(INDEX(NT_TRANS_NAME(ML),TRIM(SPECIES)) .NE. 0 .AND.
	1          ( INDEX(NT_TRANS_NAME(ML),'('//TRIM(LEVEL)) .NE. 0 .OR.
	1            INDEX(NT_TRANS_NAME(ML),'-'//TRIM(LEVEL)) .NE. 0) )THEN
	          N_NT_TRANS=N_NT_TRANS+1
	          NEW_NT_RATES(:,N_NT_TRANS)=NT_RATES(:,ML)
	          NT_LINK(N_NT_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(NT_TRANS_NAME(ML)))
	        END IF
	     END DO
	     WRITE(6,*)'Number of non_termal transitions identified is',N_NT_TRANS
	   END IF
!
! Plot the normalized rates at each depth
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'NORM')THEN
	    DO I=1,ND
	      T1=0.0D0
	      DO ML=1,N_TRANS
	         T1=MAX(T1,ABS(NEW_RATES(I,ML)))
	      END DO
	      DO ML=1,N_NT_TRANS
	         T1=MAX(T1,ABS(NEW_NT_RATES(I,ML)))
	      END DO
	      NEW_RATES(I,:)=NEW_RATES(I,:)/T1
	      NEW_NT_RATES(I,:)=NEW_NT_RATES(I,:)/T1
	    END DO
	    DO ML=1,N_TRANS
	      CALL DP_CURVE(ND,XV,NEW_RATES(:,ML))
	    END DO
!
! Examine the rates at a specific depth.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EXD')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine'
	      WRITE(6,*)'Use SPEC abd LEVELS options'
	      GOTO 1000
	    END IF
!
	    I=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth to be examined')
	    CALL GEN_IN(K,'SL identification for state '//TRIM(LEVEL))
	    T1=0.0D0; I=DPTH_INDX
	    DO ML=1,N_TRANS
	       T1=MAX(T1,ABS(NEW_RATES(I,ML)))
	    END DO
	    DO ML=1,N_NT_TRANS
	      T1=MAX(T1,ABS(NEW_NT_RATES(I,ML)))
	    END DO
	    WRITE(6,*)' '
	    WRITE(6,'(2X,A,ES16.8)'),'   R(I)=',R(I)
	    WRITE(6,'(2X,A,ES16.8)'),'   V(I)=',V(I)
	    WRITE(6,'(2X,A,ES16.8)'),'  ED(I)=',ED(I)
	    WRITE(6,'(2X,A,ES16.8)'),'Sc.Fac.=',T1
	    WRITE(6,*)' '
	    WRITE(TMP_STR,'(I2.2)')MAX_TRANS_LENGTH+6
	    TMP_STR='(A,T'//TMP_STR(1:2)//',F8.5)'
	    DO ML=1,N_TRANS
	      WRITE(6,TMP_STR)TRIM(TRANS_NAME(LINK(ML))),NEW_RATES(I,ML)/T1
	    END DO
	    IF(N_NT_TRANS .NE. 0)WRITE(6,*)'Non thermal transitions'
	    DO ML=1,N_NT_TRANS
	      WRITE(6,TMP_STR)TRIM(NT_TRANS_NAME(NT_LINK(ML))),NEW_NT_RATES(I,ML)/T1
	    END DO
!
! Compute the total rates. 
!    Note: SUM_FROM is the net rate where the level if interest is the UPPER level.
!    Note: SUM_TO   is the net rate where the level if interest is the LOWER level.
!
	    I=DPTH_INDX
	    SUM_TO=0.0D0
	    SUM_FROM=0.0D0
	    DO ML=1,N_TRANS
	      IF(INDEX(TRANS_NAME(LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        SUM_FROM=SUM_FROM+NEW_RATES(I,ML)
	      ELSE
	        SUM_TO=SUM_TO+NEW_RATES(I,ML)
	      END IF
	    END DO
	    DO ML=1,N_NT_TRANS
	      IF(INDEX(NT_TRANS_NAME(NT_LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        SUM_FROM=SUM_FROM+NEW_NT_RATES(I,ML)
	      ELSE
	        SUM_TO=SUM_TO+NEW_NT_RATES(I,ML)
	      END IF
	    END DO
!
	   WRITE(6,'(A)')
	   WRITE(6,'(A,ES14.4,F11.5)')'  SUM_TO=',SUM_TO,SUM_TO/T1
	   WRITE(6,'(A,ES14.4,F11.5)')'SUM_FROM=',SUM_FROM,SUM_FROM/T1
	   IF(PR_SPECIES .EQ. SPECIES)THEN
	     WRITE(6,'(A,ES14.4,F11.5)')'    Phot=',PHOT_RATE(K,I),PHOT_RATE(K,I)/T1
	     WRITE(6,'(A,ES14.4,F11.5)')'     Rec=',REC_RATE(K,I),REC_RATE(K,I)/T1
	     T2=REC_RATE(K,I)-PHOT_RATE(K,I)
	     WRITE(6,'(A,ES14.4,F11.5)')' Net Rec=',T2,T2/T1
	   END IF
	   T2=SUM_TO-SUM_FROM+REC_RATE(K,I)-PHOT_RATE(K,I)
	   WRITE(6,'(A,ES14.4,F11.5)')'     Net=',T2,T2/T1
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:      Set X-axis to depth index'
	    WRITE(6,*)' LOGR:    Set X-axis to Log R/R(ND)'
	    WRITE(6,*)' LOGV:    Set X-axis to Log V'
	    WRITE(6,*)' V:       Set X-axis to V'
	    WRITE(6,*)' '
	    WRITE(6,*)' SPEC:    Set the species and reads in SPCIES//PRRR file'
	    WRITE(6,*)' LEVEL:   Sets the level to be studied'
	    WRITE(6,*)' EXD:     Examine rates at a single depth'
	    WRITE(6,*)' NORM:    Plot rates as a function of depth'
	    WRITE(6,*)' E(X):    Exit routine'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
1000	  CONTINUE
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
