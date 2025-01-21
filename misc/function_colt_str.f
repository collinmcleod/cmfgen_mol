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
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL*8, ALLOCATABLE :: NT_NU(:)
	REAL*8, ALLOCATABLE :: NT_LAM(:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: LAM(:)
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
	CHARACTER(LEN=60), ALLOCATABLE :: TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: NT_TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: AUTO_LEV_NAME(:)
	CHARACTER(LEN=60)  color_lev_name
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: ED(:)
!
	CHARACTER(LEN=10) PR_SPECIES
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=30) LEVEL
!
	CHARACTER(LEN=40), ALLOCATABLE :: LEVEL_NAME(:)
	REAL*8, ALLOCATABLE :: STAT_WT(:)
	REAL*8, ALLOCATABLE :: ENERGY(:)
	REAL*8, ALLOCATABLE :: FEDGE(:)
	INTEGER, ALLOCATABLE :: F_TO_S(:)
	INTEGER, ALLOCATABLE :: INT_SEQ(:)
	REAL*8 IONIZATION_ENERGY,ZION
	CHARACTER(LEN=20) OSCDATE
	CHARACTER(LEN=20) SL_OPTION 
!
	INTEGER, ALLOCATABLE ::  LINK(:)
	REAL*8, ALLOCATABLE :: RATES(:,:)
	REAL*8, ALLOCATABLE :: NEW_RATES(:,:)
	REAL*8, ALLOCATABLE :: SUM_RATES(:)
!
	INTEGER, ALLOCATABLE ::  NT_LINK(:)
	REAL*8, ALLOCATABLE :: NT_RATES(:,:)
	REAL*8, ALLOCATABLE :: NEW_NT_RATES(:,:)
!
	REAL*8, ALLOCATABLE :: PHOT_RATE(:,:)
	REAL*8, ALLOCATABLE :: REC_RATE(:,:)
!
	REAL*8, ALLOCATABLE :: AUTO_RATE(:,:)
	REAL*8, ALLOCATABLE :: AUTO_REC_RATE(:,:)
!
	REAL*8, ALLOCATABLE :: DPTH_VEC(:)
	REAL*8, ALLOCATABLE :: WRK_VEC(:)
	INTEGER, ALLOCATABLE :: INDX(:)
!
	REAL*8 T1,T2
	REAL*8 SUM_TO,SUM_FROM
	REAL*8 NET_AUTO
	REAL*8 MAX_RATE
	REAL*8 CUR_SUM
	INTEGER ND
	INTEGER N_LINES
	INTEGER N_TRANS
	INTEGER N_NT_LINES
	INTEGER N_NT_TRANS
	INTEGER N_AUTO
	INTEGER RD_COUNT
	INTEGER I,J,K,L,ML,IBEG,IDEPTH
	INTEGER SL_INDX
	INTEGER MAX_TRANS_LENGTH
	INTEGER DPTH_INDX
	INTEGER NF            	!Number of full levels
	INTEGER NLEV		!Number of super levels
	INTEGER LU
	LOGICAL FILE_OPEN
	LOGICAL DO_NT_RATES
	LOGICAL DO_AUTO_RATES
	INTEGER IOS
!
	INTEGER, PARAMETER :: LUIN=11
	INTEGER, PARAMETER :: LUOUT=12	
	INTEGER IVEC(2)
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
	SPECIES=' ';  LEVEL=' '; PR_SPECIES=' '
	SL_OPTION=' '
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot/examine NETRATE or TOTRATE'
	WRITE(6,'(A)')' '
!
100	CONTINUE
 	FILENAME='TOTRATE'
	CALL GEN_IN(FILENAME,'File with data to be plotted/examined')
	IF(FILENAME .EQ. ' ')GOTO 100
!
! Get number of depth points
!
	OPEN(UNIT=LUIN,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=LUIN,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LUIN)
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
	ALLOCATE (DPTH_VEC(ND))
	ALLOCATE (RATES(N_LINES,ND))
	ALLOCATE (TRANS_NAME(N_LINES))
	ALLOCATE (NU(N_LINES))
	ALLOCATE (LAM(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
! 
	ALLOCATE (NEW_RATES(N_LINES,ND))
	ALLOCATE (LINK(N_LINES))
	ALLOCATE (SUM_RATES(N_LINES))
!
! Used for ploting.
!
	ALLOCATE (XV(ND),YV(ND),ZV(ND))
!
! These are updated when the data file is read in.
!
	DO_NT_RATES=.FALSE.
	CALL GEN_IN(DO_NT_RATES,'Read in non thermal rates')
	IF(DO_NT_RATES)THEN
	  ALLOCATE (NT_RATES(N_LINES,ND))
	  ALLOCATE (NT_TRANS_NAME(N_LINES))
	  ALLOCATE (NT_NU(N_LINES))
	  ALLOCATE (NT_LAM(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
! 
	  ALLOCATE (NEW_NT_RATES(N_LINES,ND))
	  ALLOCATE (NT_LINK(N_LINES))
	END IF
!
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	RD_COUNT=0
	CALL TUNE(1,'READ')
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	    READ(LUIN,'(A)',END=5000)STRING
	  END DO
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')
	  STRING=STRING(K+2:)
	  K=INDEX(STRING,'  ')
	  TRANS_NAME(ML)=STRING(1:K)
	  READ(STRING(K+1:),*)NU(ML)
	  RD_COUNT=ML
	  READ(LUIN,*)(RATES(ML,I),I=1,ND)
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
	ALLOCATE (R(ND),V(ND),ED(ND),SIGMA(ND))
	CALL RD_SING_VEC_RVTJ(R,ND,'Radius',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(V,ND,'Velocity',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(SIGMA,ND,'dlnV/dlnr-1',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(ED,ND,'Electron',FILENAME,LU,IOS)
!
	IF(DO_NT_RATES)THEN
	  OPEN(UNIT=LUIN,FILE='NT_RATES',STATUS='OLD',ACTION='READ')
!
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')' Reading data -- this may take a while'
	  WRITE(6,'(A)')DEF_PEN
	  RD_COUNT=0
	  CALL TUNE(1,'READ_NT')
	  DO ML=1,N_LINES             !NT_LINES
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	      READ(LUIN,'(A)',END=6000)STRING
	    END DO
	    STRING=ADJUSTL(STRING)
	    K=INDEX(STRING,'  ')
	    NT_TRANS_NAME(ML)=STRING(1:K)
	    READ(STRING(K+1:),*)NT_NU(ML)
	    RD_COUNT=ML
	    READ(LUIN,*)(NT_RATES(ML,I),I=1,ND)
	  END DO
6000	  CONTINUE
	  N_NT_LINES=RD_COUNT
	  CALL TUNE(2,'READ_NT')
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,*)'Number of number of nonthermal lines read is',N_NT_LINES 
	  WRITE(6,'(A)')DEF_PEN
          NT_LAM(1:N_LINES)=2.99792458D+03/NT_NU(1:N_NT_LINES)
	ELSE
	  N_NT_LINES=0
	END IF
	CALL TUNE(3,'TERMINAL')
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
	  IF(UC(PLT_OPT) .EQ. 'XLOGR')THEN
	    T1=R(ND)
	    DO I=1,ND
	      XV(I)=LOG10(R(I)/T1)
	    END DO
	    XLABEL='Log R/R(ND)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XV')THEN
	    XV(1:ND)=V(1:ND)
	    XLABEL='V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XLOGV')THEN
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
	    CALL GEN_IN(SPECIES,'Species that will be examined: e.g. C2, NIII')
	    SPECIES=UC(SPECIES)
!
! Get number of super levels. This is need when we read in the the
! recombination and photoioizations rates for each super level.
!
	    FILENAME='MODEL_SPEC'
	    OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	        GOTO 1000
	      END IF
	      STRING=' '
	      DO WHILE(INDEX(STRING,TRIM(SPECIES)) .EQ. 0)
	        READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      END DO
	      IF(IOS .EQ. 0)THEN
	        READ(STRING,*)I,NLEV,NF
	        WRITE(6,*)NF,NLEV
	      ELSE
	        WRITE(6,*)'Error - species not found: IOS=',IOS
	        CLOSE(UNIT=LUIN)
	        GOTO 1000
	      END IF
	    CLOSE(UNIT=LUIN)
	    IF(ALLOCATED(PHOT_RATE))DEALLOCATE(PHOT_RATE,REC_RATE)
	    ALLOCATE (PHOT_RATE(NLEV,ND),REC_RATE(NLEV,ND))
!
! Read in the recombination and photoioizations rates for each super level.
! Note that ND is the second index.
!
	    FILENAME=TRIM(SPECIES)//'PRRR'
	    OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	      WRITE(6,*)'Was you species name correct: Species=',TRIM(SPECIES)
	      GOTO 1000
	    END IF
!
	    DO K=1,ND,10
	      DO WHILE(INDEX(STRING,'Photoionization') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(LUIN,*)(PHOT_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	      DO WHILE(INDEX(STRING,'Recombination') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(LUIN,*)(REC_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	    END DO
	    CLOSE(UNIT=LUIN)
	    PR_SPECIES=SPECIES 
!
	    DO_AUTO_RATES=.FALSE.
	    CALL GEN_IN(DO_AUTO_RATES,'Read in autoioization rates?')
	    IF(DO_AUTO_RATES)THEN
	      FILENAME=TRIM(SPECIES)//'_AUTO_RATES'
	      CALL GEN_IN(FILENAME,'File with auto rates')
	      OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      STRING=' '
	      DO WHILE(INDEX(STRING,'Number of ') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      READ(STRING,*)N_AUTO
	      READ(LUIN,'(A)')STRING; STRING=' '
	      DO WHILE(INDEX(STRING,'!') .NE. 0 .OR. STRING .EQ. ' ')
	        READ(LUIN,'(A)')STRING
	      END DO
	      BACKSPACE(UNIT=LUIN)
	      IF(ALLOCATED(AUTO_RATE))DEALLOCATE(AUTO_RATE,AUTO_REC_RATE,AUTO_LEV_NAME)
	      ALLOCATE (AUTO_RATE(N_AUTO,ND))
	      ALLOCATE (AUTO_REC_RATE(N_AUTO,ND))
	      ALLOCATE (AUTO_LEV_NAME(N_AUTO))
	      DO I=1,N_AUTO
	        STRING=' '
	        DO WHILE(INDEX(STRING,'!') .NE. 0 .OR. STRING .EQ. ' ')
	          READ(LUIN,'(A)')STRING
	        END DO
	        STRING=ADJUSTL(STRING)
	        K=INDEX(STRING,'  ');  AUTO_LEV_NAME(I)=STRING(1:K-1)
	        READ(LUIN,*)(AUTO_REC_RATE(I,K),AUTO_RATE(I,K),K=1,ND)
	      END DO
	      CLOSE(UNIT=LUIN)
	    END IF
	    LEVEL=' '
	    N_TRANS=0
!
	    IF(ALLOCATED(LEVEL_NAME))DEALLOCATE(LEVEL_NAME,STAT_WT,ENERGY,FEDGE)
	    ALLOCATE (LEVEL_NAME(NF),STAT_WT(NF),ENERGY(NF),FEDGE(NF))
	    FILENAME=TRIM(SPECIES)//'_F_OSCDAT'
	    WRITE(6,*)'Getting energy levels',TRIM(SPECIES),'/',NF
	    CALL RD_ENERGY(LEVEL_NAME,STAT_WT,ENERGY,FEDGE,NF,NF,
	1             IONIZATION_ENERGY,ZION,OSCDATE,FILENAME,LUIN,LUOUT,IOS)
	    WRITE(6,*)'Got energy levels'
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to read in energy level names'
	      GOTO 1000
	    END IF
	    CLOSE(LUIN)
! 
	    WRITE(6,*)'Setting SL assignent'
	    IF(ALLOCATED(F_TO_S))DEALLOCATE(F_TO_S,INT_SEQ)
	    ALLOCATE (F_TO_S(NF),INT_SEQ(NF))
	    FILENAME=TRIM(SPECIES)//'_F_TO_S'
	    WRITE(6,*)TRIM(FILENAME)
	    CALL RD_F_TO_S_IDS_V2(F_TO_S,INT_SEQ,LEVEL_NAME,NF,NLEV,LUIN,.FALSE.,TRIM(FILENAME)) !,SL_OPTION)
	    WRITE(6,*)'Read SLs'
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
	          NEW_RATES(N_TRANS,:)=RATES(ML,:)
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
	          NEW_NT_RATES(N_NT_TRANS,:)=NT_RATES(ML,:)
	          NT_LINK(N_NT_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(NT_TRANS_NAME(ML)))
	        END IF
	     END DO
	     WRITE(6,*)'Number of non_termal transitions identified is',N_NT_TRANS
	   END IF
!
! Get the SL link.
!
	  DO I=1,NF
	    IF(INDEX(LEVEL_NAME(I),TRIM(LEVEL)) .NE. 0)THEN
	      SL_INDX=F_TO_S(I)
	      EXIT 
	    END IF
	  END DO
! 
	  WRITE(6,*)' '
	  WRITE(6,*)' The following levels have the same super level assignment'
	  DO I=1,NF
	    IF(F_TO_S(I) .EQ. SL_INDX)THEN
	      WRITE(6,'(4X,A,T50,A)')TRIM(LEVEL_NAME(I)),TRIM(LEVEL)
	    END IF
	  END DO  
!
! Plot the normalized rates at each depth
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'NORM')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine'
	      WRITE(6,*)'Use SPEC and LEVELS options'
	      GOTO 1000
	    END IF
!
	    DO I=1,ND
	      T1=0.0D0
	      DO ML=1,N_TRANS
	         T1=MAX(T1,ABS(NEW_RATES(ML,I)))
	      END DO
	      DO ML=1,N_NT_TRANS
	         T1=MAX(T1,ABS(NEW_NT_RATES(ML,I)))
	      END DO
	      NEW_RATES(:,I)=NEW_RATES(:,I)/T1
	      IF(DO_NT_RATES)NEW_NT_RATES(:,I)=NEW_NT_RATES(:,I)/T1
	    END DO
	    DO ML=1,N_TRANS
	      DPTH_VEC=NEW_RATES(ML,:)
	      CALL DP_CURVE(ND,XV,DPTH_VEC)
	    END DO
!
! Examine the rates at a specific depth.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'RATE')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine'
	      WRITE(6,*)'Use SPEC and LEVELS options'
	      GOTO 1000
	    END IF
!
	    CALL DERIVCHI(YV,XV,R,ND,'LINEAR')
	    WRITE(6,*)'Done deriv'
	    DO I=1,ND
	      NEW_RATES(1:N_TRANS,I)=NEW_RATES(1:N_TRANS,I)*R(I)*R(I)/XV(I)
	      IF(DO_NT_RATES)NEW_NT_RATES(:,I)=NEW_NT_RATES(:,I)*R(I)*R(I)/XV(I)
	    END DO
	    SUM_RATES=0.0D0
	    DO I=1,ND
	      DO ML=1,N_TRANS
	        SUM_RATES(ML)=SUM_RATES(ML)+ABS(NEW_RATES(ML,I))
	       END DO
	    END DO
	    WRITE(6,*)'Set sum rates'
!
	    DO ML=1,MIN(6,N_TRANS)
	      J=MAXLOC(SUM_RATES(1:N_TRANS),1)
	      DPTH_VEC=NEW_RATES(J,:)
	      CALL DP_CURVE(ND,XV,DPTH_VEC)
	      L=LINK(J)
              WRITE(6,'(1X,A)')TRIM(TRANS_NAME(LINK(L)))
	      SUM_RATES(J)=0.0D0
	    END DO
	    WRITE(6,*)'Done line curves'
!
	    IF(DO_AUTO_RATES)THEN
	      T1=0.0D0
	      YV=0.0D0; ZV=0.0D0
	      DO J=1,N_AUTO
	         IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	           YV=YV+AUTO_REC_RATE(J,:)
	           ZV=ZV+AUTO_RATE(J,:)
	         END IF
	      END DO
	      CALL DP_CURVE(ND,XV,YV)
1	      CALL DP_CURVE(ND,XV,ZV)
	      WRITE(6,'(1X,A,T30,ES12.3)')'Maximum auto/anti autoionization rate',T1
	    END IF
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
	    DPTH_INDX=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth to be examined')
!
! We first deduce the maximum rate to each level.
!
	    MAX_RATE=0.0D0
	    IF(DO_AUTO_RATES)THEN
	      T1=0.0D0
	      DO J=1,N_AUTO
	         IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	           T1=MAX(T1,AUTO_REC_RATE(J,DPTH_INDX))
	           T1=MAX(T1,AUTO_RATE(J,DPTH_INDX))
	         END IF
	      END DO
	      WRITE(6,'(1X,A,T30,ES12.3)')'Maximum auto/anti autoionization rate',T1
	    END IF
	    MAX_RATE=T1
!
	    T1=REC_RATE(SL_INDX,DPTH_INDX)
	    T1=MAX(T1,PHOT_RATE(SL_INDX,DPTH_INDX))
	    WRITE(6,*)'Maximum phot/rec rate:=',T1
	    MAX_RATE=MAX(MAX_RATE,T1)
!
	    I=DPTH_INDX; T1=0D0
	    DO ML=1,N_TRANS
	       T1=MAX(T1,ABS(NEW_RATES(ML,I)))
	    END DO
	    DO ML=1,N_NT_TRANS
	      T1=MAX(T1,ABS(NEW_NT_RATES(ML,I)))
	    END DO
	    WRITE(6,'(1X,A,T30,ES12.3)')'Maximum line rate=',T1
	    MAX_RATE=MAX(MAX_RATE,T1)

!
	    WRITE(6,*)' '
	    WRITE(6,'(2X,A,I5)')'Depth index=',DPTH_INDX
	    WRITE(6,'(2X,A,ES16.8)'),'     R(DPTH_INDX) =',R(DPTH_INDX)
	    WRITE(6,'(2X,A,ES16.8)'),'     V(DPTH_INDX) =',V(DPTH_INDX)
	    WRITE(6,'(2X,A,ES16.8)'),'    ED(DPTH_INDX) =',ED(DPTH_INDX)
	    WRITE(6,'(2X,A,ES16.8)'),'  Scaling Factor  =',MAX_RATE
	    WRITE(6,*)' '
	    WRITE(6,*)' Upper level is listed first in transition name.'
	    WRITE(6,*)' Negative rate imples there is a net flow to the upper state'
	    WRITE(6,*)' '
!
	    ALLOCATE (WRK_VEC(N_TRANS),INDX(N_TRANS))
	    WRK_VEC(:)=ABS(NEW_RATES(1:N_TRANS,DPTH_INDX))
	    CALL INDEXX(N_TRANS,WRK_VEC,INDX,L_TRUE)
!
	    WRITE(TMP_STR,'(I2.2)')MAX_TRANS_LENGTH+10
	    TMP_STR='(1X,A,T'//TMP_STR(1:2)//',F8.5,3X,F12.5,3X,F8.5)'
	    DO ML=1,N_TRANS
	      L=INDX(ML)
	      string=color_lev_name(TRIM(TRANS_NAME(LINK(L))),LEVEL)
	      WRITE(6,TMP_STR)color_lev_name(TRIM(TRANS_NAME(LINK(L))),LEVEL),NEW_RATES(L,DPTH_INDX)/MAX_RATE,LAM(LINK(L))
	    END DO
	    DEALLOCATE(WRK_VEC)
!
	    IF(N_NT_TRANS .NE. 0)WRITE(6,*)'Non thermal transitions'
	    DO ML=1,N_NT_TRANS
	      WRITE(6,TMP_STR)TRIM(NT_TRANS_NAME(NT_LINK(ML))),
	1              NEW_NT_RATES(ML,I)/MAX_RATE,NT_LAM(NT_LINK(ML))
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
	        SUM_FROM=SUM_FROM+NEW_RATES(ML,I)
	      ELSE
	        SUM_TO=SUM_TO+NEW_RATES(ML,I)
	      END IF
	    END DO
	    DO ML=1,N_NT_TRANS
	      IF(INDEX(NT_TRANS_NAME(NT_LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        SUM_FROM=SUM_FROM+NEW_NT_RATES(ML,I)
	      ELSE
	        SUM_TO=SUM_TO+NEW_NT_RATES(ML,I)
	      END IF
	    END DO
!
	    WRITE(6,'(A)')
	    WRITE(6,'(A,ES14.4,F11.5)')' Net rate into from higher levels =',SUM_TO,SUM_TO/MAX_RATE
	    WRITE(6,'(A,ES14.4,F11.5)')' Net rate into lower levels       =',SUM_FROM,SUM_FROM/MAX_RATE
	    IF(PR_SPECIES .EQ. SPECIES)THEN
	      WRITE(6,'(A,ES14.4,F11.5)')' Photoioization rate              =',
	1             PHOT_RATE(SL_INDX,DPTH_INDX),PHOT_RATE(SL_INDX,DPTH_INDX)/MAX_RATE
	      WRITE(6,'(A,ES14.4,F11.5)')' Recombination rate               =',
	1             REC_RATE(SL_INDX,DPTH_INDX),REC_RATE(SL_INDX,DPTH_INDX)/MAX_RATE
	      T2=REC_RATE(SL_INDX,DPTH_INDX)-PHOT_RATE(SL_INDX,DPTH_INDX)
	      WRITE(6,'(A,ES14.4,F11.5)')' Net recombination rate           =',T2,T2/MAX_RATE
	    END IF
!
	    WRITE(6,'(A)')' '
	    T1=0.0D0; I=DPTH_INDX
	    DO ML=1,N_TRANS
	      L=INDX(ML)
	      IF( (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .NE.  0 .AND. NEW_RATES(L,I) .GT. 0) .OR.
	1         (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .EQ.  0 .AND. NEW_RATES(L,I) .LT. 0) )THEN
	        T1=T1+ABS(NEW_RATES(L,I))
                WRITE(6,TMP_STR)TRIM(TRANS_NAME(LINK(L))),NEW_RATES(L,I)/MAX_RATE,LAM(LINK(L)),T1/MAX_RATE
	      END IF
	    END DO
!
	    WRITE(6,'(A)')' '
	    T1=0.0D0
	    DO ML=1,N_TRANS
	      L=INDX(ML)
	      IF( (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .EQ.  0 .AND. NEW_RATES(L,I) .GT. 0) .OR.
	1         (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .NE.  0 .AND. NEW_RATES(L,I) .LT. 0) )THEN
	        T1=T1+ABS(NEW_RATES(L,I))
                WRITE(6,TMP_STR)TRIM(TRANS_NAME(LINK(L))),
	1         NEW_RATES(L,I)/MAX_RATE,LAM(LINK(L)),T1/MAX_RATE
	      END IF
	    END DO
!
	    IF(DO_AUTO_RATES)THEN
	      NET_AUTO=0.0D0
	      DO J=1,N_AUTO
	        IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	          WRITE(6,'(A)')' '
	          WRITE(6,'(A,ES14.4,F11.5,5X,A)')' Autoionizaton rate               =',
	1                AUTO_RATE(J,DPTH_INDX), AUTO_RATE(J,DPTH_INDX)/MAX_RATE, AUTO_LEV_NAME(J) 
	          WRITE(6,'(A,ES14.4,F11.5)')' Rec. autoioization rate          =',
	1                AUTO_REC_RATE(J,DPTH_INDX), AUTO_REC_RATE(J,DPTH_INDX)/MAX_RATE 
	          NET_AUTO=NET_AUTO+(AUTO_REC_RATE(J,DPTH_INDX)-AUTO_RATE(J,DPTH_INDX))
	          WRITE(6,'(A,ES14.4,F11.5)')' Net recom. autoization rate      =',NET_AUTO,NET_AUTO/MAX_RATE
	        END IF
	      END DO
	    END IF
	    WRITE(6,'(A)')' '
	    T2=(SUM_TO-SUM_FROM)+(REC_RATE(SL_INDX,I)-PHOT_RATE(SL_INDX,I))+NET_AUTO
	    WRITE(6,'(A,T35,A,ES14.4,F11.5)')' Net Rate','=',T2,T2/MAX_RATE
!
	    DEALLOCATE(INDX)
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:      Set X-axis to depth index'
	    WRITE(6,*)' XLOGR:   Set X-axis to Log R/R(ND)'
	    WRITE(6,*)' XLOGV:   Set X-axis to Log V'
	    WRITE(6,*)' XV:      Set X-axis to V'
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
cotains
	function color_lev_name(trans_name,level) result(res)
	USE MOD_COLOR_PEN_DEF
	character(len=*) res
	character(len=*), intent(in) :: trans_name
	character(len=*), intent(in) :: level 
	character(len=100) string
	integer ist,iend
	integer lev_len
!
	string=trans_name
	lev_len=len_trim(level)
	ist=index(string,trim(level))
	iend=index(string,'-')
	if(ist .lt.iend)then
	  iend=index(string,'-')
	  string=BLUE_PEN//string(1:iend-1)//DEF_PEN//string(iend:)
	else if(ist .ne. 0)then
	  ist=iend
	  string=string(1:ist)//BLUE_PEN//string(ist+1:)
	  iend=len_trim(string)
	  string=string(1:iend)//DEF_PEN
	end if
	res = trim(string)
	end function
