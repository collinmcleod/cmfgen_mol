!
! Simple program to plot data from the SCRTEMP file. This can be used to check
!   Convergence of a program.
!   Convergence as a function of depth etc.
!
	PROGRAM PLT_SCR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
	INTEGER ND,NT,NIT
!
! Altered 24-Sep-2023: Added SFDG option (needs work) - 7sep23.
!                      Added DREP option -- allows one depth to replace other depths.
!                      Added (modified) DNRG option -- designed to allow testing of subroutine.
!                      Added CNRG option -- designed to allow modification of R grid using cursors.
! Altered 04-Jul-2020: Altered FDG option so harder to put in depth for variable.
! Altered 14-Feb-2019: Updated to 3 digit exponent in 2 places
! Altered 13-Sep-2018: Added options SM and NINT
! Altered 06-Dec-2017: Made compatible with osiris version
!                        For multiple fudges, only one new record is now written.
!                        REP and RAT options inserted                  
! Altered 26-Feb-2015: Fixed error with FDG ouput if called multiple times.
! Altered 30-Mar-2014: Altered FDG option so as to print out adjacent values.
! Altered 11-Mar-2014: Installed INT option.
! Altered 25-Feb-2014: Modified WRST option.
! Altered 12-Dec-2013: Installed LY option so that we can plot very small populations on
!                        a logarithmic scale.
!
	REAL(10), ALLOCATABLE :: POPS(:,:,:)		!NT,ND,NIT
	REAL(10), ALLOCATABLE :: R_MAT(:,:)		!ND,NIT
	REAL(10), ALLOCATABLE :: V_MAT(:,:)		!ND,NIT
	REAL(10), ALLOCATABLE :: RAT(:,:)			!NT,ND
	REAL(10), ALLOCATABLE :: CORRECTIONS(:,:)			!NT,ND
!
	REAL(10), ALLOCATABLE :: R(:)			!ND
	REAL(10), ALLOCATABLE :: V(:)			!ND
	REAL(10), ALLOCATABLE :: SIGMA(:)			!ND
C
	REAL(10), ALLOCATABLE :: X(:)			!NIT
	REAL(10), ALLOCATABLE :: Y(:)			!NIT
	REAL(10), ALLOCATABLE :: Z(:)			!NIT
	REAL(10), ALLOCATABLE :: TA(:)			!
	REAL(10), ALLOCATABLE :: TB(:)			!
!
	INTEGER, ALLOCATABLE :: I_BIG(:)		!NT
	REAL(10), ALLOCATABLE :: Z_BIG(:)			!NT
C
	INTEGER, PARAMETER :: NUM_IONS_MAX=500
	INTEGER ION_INDEX(NUM_IONS_MAX)
	INTEGER NION(NUM_IONS_MAX)
	INTEGER NUM_IONS
	INTEGER CNT,CNT_NEG,CNT_POS
	INTEGER NAN_CNT
	LOGICAL NANS_PRESENT
	LOGICAL NORM_R
	CHARACTER(LEN=10) ION_ID(NUM_IONS_MAX)
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10
	INTEGER, PARAMETER :: LU_OUT=20
C
	INTEGER, SAVE :: FDG_COUNTER=0
	INTEGER, PARAMETER :: NLIM_MAX=10
	INTEGER LIMITS(NLIM_MAX)
	INTEGER NPLTS
	INTEGER IREC
	INTEGER LST_IREC
	INTEGER IVAR
	INTEGER IV1,IV2
	INTEGER ID1,ID2
	INTEGER LOW_ID,UP_ID
	INTEGER I
	INTEGER J
	INTEGER K
	INTEGER ID
	INTEGER IT,IT2
	INTEGER NY
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER RITE_N_TIMES
	INTEGER LUSCR
	INTEGER IOS
!
	INTEGER, PARAMETER :: NLEV_MAX=10
	INTEGER LEVELS(NLEV_MAX)
C
	LOGICAL LOG_Y_AXIS
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
	LOGICAL DO_ABS
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*10 TMP_STR
	CHARACTER*10 PLT_OPT
	CHARACTER*80 XLABEL
	CHARACTER*80 YLABEL
	CHARACTER*200 STRING
C
	INTEGER IMAX,IMIN
	REAL(10) RMAX,RMIN
	REAL(10) T1,T2,T3
	REAL(10) MAX_CHANGE
	REAL(10) AVE_NEG_CHANGE, MIN_NEG_CHANGE
	REAL(10) AVE_POS_CHANGE, MAX_POS_CHANGE
	LOGICAL COR_READ
	LOGICAL READ_IN_STEQ
	LOGICAL READ_AGAIN
	LOGICAL TMP_LOG
C
	LUSCR=26
	RITE_N_TIMES=1
	NEWMOD=.TRUE.
	COR_READ=.FALSE.
	READ_AGAIN=.FALSE.
	READ_IN_STEQ=.TRUE.
!C
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       POINT1.DAT'
	WRITE(T_OUT,*)'                                       SCRTEMP.DAT'
	WRITE(T_OUT,*)'                                       MODEL.DAT'
	WRITE(T_OUT,*)' '
C
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
C
100	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	ELSE
	  CLOSE(UNIT=12)
	  CALL RD_ION_LOCATIONS(ION_ID,NION,ION_INDEX,NUM_IONS,NUM_IONS_MAX,LU_IN)
	END IF 
C
	OPEN(UNIT=12,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)READ(12,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .EQ. 0)THEN
            IF(INDEX(STRING,'!Format date') .EQ. 0)THEN
              READ(STRING,*,IOSTAT=IOS)K,NIT
	    ELSE
              READ(12,*,IOSTAT=IOS)K,NIT
	    END IF
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Possible error reading POINT1'
	    CALL GEN_IN(NIT,'Number of iterations')
	  END IF
	CLOSE(UNIT=12)
C
	ALLOCATE (POPS(NT,ND,NIT))
	ALLOCATE (R_MAT(ND,NIT))
	ALLOCATE (V_MAT(ND,NIT))
	ALLOCATE (CORRECTIONS(NT,ND))
	ALLOCATE (RAT(NT,ND))
!
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
C
	J=MAX(NIT,ND); J=MAX(NT,J)
	ALLOCATE (X(J))
	ALLOCATE (Y(J))
	ALLOCATE (Z(J))
	ALLOCATE (TA(J))
	ALLOCATE (TB(J))
!
	WRITE(6,*)'   Number of depth points is:',ND
	WRITE(6,*)'Number of variables/depth is:',NT
	WRITE(6,*)'     Number of iterations is:',NIT
!
	ALLOCATE (I_BIG(NT))
	ALLOCATE (Z_BIG(NT))
C
	DO IREC=1,NIT
	  CALL SCR_READ_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  R_MAT(:,IREC)=R(:) 
	  V_MAT(:,IREC)=V(:) 
	END DO
!
	NANS_PRESENT=.FALSE.
	DO K=1,NIT
	  DO ID=1,ND
	    CNT=0
	    NAN_CNT=0
	    DO IVAR=1,NT
	      IF(POPS(IVAR,ID,K) .LE. 0.0D0 .AND. CNT .LE.4)THEN
	         WRITE(6,*)'Invalid population at depth ',IVAR,ID,NT
	         CNT=CNT+1
	      END IF
	      IF(POPS(IVAR,ID,K) .NE. POPS(IVAR,ID,K))THEN
	         IF(NAN_CNT .LE. 4)WRITE(6,*)'NaN at depth ',IVAR,ID,K
	         NAN_CNT=NAN_CNT+1
	    END IF
	    END DO
	    IF(CNT .NE. 0)WRITE(6,'(A,I3,A,I5,A,I3)')
	1          'Number of invalid poulations at depth ',ID,' is',CNT,' for it',K
	    IF(NAN_CNT .NE. 0)WRITE(6,'(A,I3,A,I5,A,I3)')
	1          'Number of NaNs at depth ',ID,' is',NAN_CNT,' for it',K
	    IF(NAN_CNT .NE. 0)NANS_PRESENT=.TRUE.
	   END DO
	END DO
!
	IF(NANS_PRESENT)THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Your data contains NaNs'
	  WRITE(6,'(A)')'Depending on the iteration, you may need to reset POINT1 to an eralier iteration'
	  WRITE(6,'(A)',ADVANCE='NO')'Input any character to continue'
	  READ(5,'(A)')STRING
	  WRITE(6,'(A)')DEF_PEN
	END IF 
C
200	CONTINUE
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'ND  :: ',ND
	WRITE(T_OUT,*)'NT  :: ',NT
	WRITE(T_OUT,*)'NIT :: ',NIT
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)' For the next 4 options, we plot versus iteration number '
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'F   :: Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)'
	WRITE(T_OUT,*)'R   :: [Y(K+2)-Y(K+1)]/[Y(K+1)-Y(K)]'
	WRITE(T_OUT,*)'D   :: Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)'
	WRITE(T_OUT,*)'Y   :: Z(K)=Y(K)'
	WRITE(T_OUT,*)' '
        WRITE(T_OUT,*)'PD  :: Plot 100.0D0*(Y(K)-Y(K-1))/Y(K) as a function of depth index.'
	WRITE(T_OUT,*)'PF  :: Plot 100.0D0*(Y(K+1)-Y(K))/Y(K+1) for all variables at a given depth.'
	WRITE(T_OUT,*)' '
        WRITE(T_OUT,*)'PN  :: Plot a variable as a function of depth index.'
        WRITE(T_OUT,*)'PV  :: Plot a variable as a function of velocity.'
        WRITE(T_OUT,*)'PYD :: Plot all variables at a given depth -- change to log space before using'
        WRITE(T_OUT,*)'VR  :: Plot velocity as a function of radius.'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'MED_R  :: Median corection as a function of depth'
	WRITE(T_OUT,*)'MR     :: Z(K)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)'IR     :: Z(ID)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)'WRST   :: Writes fractional corections to file (FRAC_COR -- same format as STEQ_VALS'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)' The next set of options allow corections to the populations to be made.'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'FDG      :: Fudge individual values at a single depth and output to SCRTEMP'
	WRITE(T_OUT,*)'FDGV     :: Fudge values over a ranges of depths (% change) and output to SCRTEMP'
	WRITE(T_OUT,*)'SM       :: Fudge values over a ranges of depths using values at higher/adjacent depth'
	WRITE(T_OUT,*)'INT      :: Interpolate values using adjacent depths whose',
	1                             ' corrections are above a certain % limit'
	WRITE(T_OUT,*)'NINT     :: Interpolate values whose corrections are above a certain % limit'
	WRITE(T_OUT,*)'RAT      :: Compare populations at adjacent depths'
	WRITE(T_OUT,*)'REP      :: Replace populations on one iteration with those of another'
	WRITE(T_OUT,*)'FDG_OSC  ::'
	WRITE(T_OUT,*)'UNDO     :: Undo corrections over a range of depths'
	WRITE(T_OUT,*)' ' 
	WRITE(T_OUT,*)'LY  :: Switch to/from Log(Y) for options where appropriate (not full implemented)'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'E   :: EXIT'
	WRITE(T_OUT,*)' '
	PLT_OPT='R'
	CALL GEN_IN(PLT_OPT,'Plot option: R(atio), F(rac) or D(elta), Y, IR, MR or E(xit)')
	CALL SET_CASE_UP(PLT_OPT,0,0)
!
	IF(PLT_OPT(1:2) .EQ. 'E ' .OR. PLT_OPT(1:2) .EQ. 'EX')STOP
	IF(PLT_OPT(1:2) .EQ. 'LY')THEN
	   WRITE(6,*)RED_PEN
	  IF(LOG_Y_AXIS)WRITE(6,*)'Switching to linear Y axis'
	  IF(.NOT. LOG_Y_AXIS)WRITE(6,*)'Switching to logarithmic Y axis'
	  WRITE(6,'(A)')DEF_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  LOG_Y_AXIS=.NOT. LOG_Y_AXIS
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'CHK_R')THEN
	  DO IT=1,NIT
	    WRITE(6,'(2ES16.6)')R_MAT(1,IT),R_MAT(ND,IT)
	  END DO
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'WRST')THEN
	  I=1; K=NIT
	  CALL GEN_IN(K,'Iteration to check')
	  CALL GEN_IN(I,'Comparison iteration')
	  RAT(:,:)=(POPS(:,:,I)-POPS(:,:,K))/POPS(:,:,I)
	  OPEN(UNIT=11,FILE='FRAC_COR',STATUS='UNKNOWN')
	    WRITE(TMP_STR,'(I4)')K; TMP_STR=ADJUSTL(TMP_STR)
	    STRING='Fractional changes: 1- POP(IT='//TRIM(TMP_STR)
	    WRITE(TMP_STR,'(I4)')I; TMP_STR=ADJUSTL(TMP_STR)
	    STRING=TRIM(STRING)//')/POP(IT='//TRIM(TMP_STR)//')'
	    CALL WR2D_V2(RAT,NT,ND,STRING,'#',L_TRUE,11)
	  CLOSE(UNIT=11)
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'CHNG')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT
	    CALL GEN_IN(IT,'Iteraton # of lower iteration')
	  END DO
	  ID=ND
	  CALL GEN_IN(ID,'Depth index')
	  DO IVAR=1,NT
	    X(IVAR)=IVAR
	    T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	    T1=POPS(IVAR,ID,NIT)-POPS(IVAR,ID,IT-1)
	    IF(ABS(T2). GT. 1.0D-50*ABS(T1))THEN
	     Y(IVAR)=T1/T2
	    ELSE
	     Y(IVAR)=20000
	    END IF
	  END DO
	  CALL DP_CURVE(NT,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID','Ratio',' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:5) .EQ. 'MED_R')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(ABS(T2). GT. 1.0D-50*ABS(T1))THEN
	        Z_BIG(IVAR)=T1/T2
	      ELSE
	        Z_BIG(IVAR)=20000
	      END IF
	    END DO
!
! Now find median value. First order.
!
	    CALL INDEXX(NT,Z_BIG,I_BIG,.TRUE.)
	    J=(NT+1)/2
	    Y(ID)=100.0D0*(Z_BIG(I_BIG(J))-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(r\dmed\u-1)'
!	    IF(ID .EQ. 30)THEN
	      WRITE(74,*)'ID=',ID
	      DO J=1,NT
	        K=I_BIG(J)
!	        WRITE(74,*)J,I_BIG(J),Z_BIG(I_BIG(J))
	        WRITE(74,'(2I4,4ES16.8)')J,K,Z_BIG(K),POPS(K,ID,IT),POPS(K,ID,IT-1),POPS(K,ID,IT-2)
	      END DO
!	    END IF
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'MR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO_ABS=.TRUE.; CALL GEN_IN(DO_ABS,'Use absolute value')
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        T1=T1/T2
	        IF(DO_ABS)T1=ABS(T1)
	        Y(ID)=Y(ID)+ T1
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR
	      END IF
	    END DO
	    Y(ID)=100.0D0*(Y(ID)/K-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'HR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT
	    CALL GEN_IN(IT,'Iteration #')
	  END DO
	  ID=0
	  DO WHILE(ID .LT. 1 .OR. ID .GT. ND)
	    ID=ND/2; CALL GEN_IN(ID,'Depth index')
	  END DO
	  Y(ID)=0.0D0
!
	  K=MIN(201,NT)
	  T3=200.0D0/(K-1)
	  DO J=1,K
	    X(J)=-50.0D0+(J-1)*T3
	  END DO
	  Y(1:K)=0.0D0
!
	  DO IVAR=1,NT
	    T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	    T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	    IF(T2 .NE. 0)THEN
	      T1=(100*(T1/T2-1.0D0)+51)/T3
	      IF(T1 .LT. 1)T1=1
	      IF(T1 .GT. K)T1=K
	      Y(NINT(T1))=Y(NINT(T1))+1
	    ELSE
	      Y(1)=Y(1)+1
	    END IF
	    Ylabel='N(r)'
	  END DO
	  CALL DP_CURVE(K,X,Y)
	  CALL GRAMON_PGPLOT('100(r-1)',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'IR')THEN
	  ID=0
	  DO WHILE(ID .LT. 1 .OR. ID .GT. ND)
	    ID=ND; CALL GEN_IN(ID,'Depth index')
	  END DO
	  DO IT=3,NIT
	    Y(IT)=0.0D0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        Y(IT)=Y(IT)+ T1/T2
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR,'for iteration',IT
	      END IF
	    END DO
	    Y(IT)=100.0D0*(Y(IT)/NT-1.0D0)
	    X(IT)=FLOAT(IT)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  IT=NIT-2
	  CALL DP_CURVE(IT,X(3),Y(3))
	  CALL GRAMON_PGPLOT('Iteration',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PF')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depth (zero to exit)')
	    IF(ID .EQ. 0)EXIT
	    DO IVAR=1,NT
	      Y(IVAR)=100.0D0*(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	      X(IVAR)=IVAR
	    END DO
	    CALL DP_CURVE(NT,X,Y)
	  END DO
	  Ylabel=''
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'PYD')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depths',LOW_LIM=IONE,UP_LIM=ND)
	    IF(ID .EQ. 0)EXIT
	    DO IVAR=1,NT
	      Y(IVAR)=POPS(IVAR,ID,IT)
	      Z(IVAR)=LOG10(POPS(IVAR,ID,IT))
	      X(IVAR)=IVAR
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(NT,X,Z)
	      Ylabel='Log (Pop)'
	    ELSE
	      CALL DP_CURVE(NT,X,Y)
	      Ylabel='Pop'
	    END IF
	  END DO
	  CALL GRAMON_PGPLOT('Variable',Ylabel,' ',' ')
	  GOTO 200
!
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'PD')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    DO ID=1,ND
	      Y(ID)=100.0D0*(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	      X(ID)=ID
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  YLABEL='[Y(I-1)-Y(I)]/Y(I) [%]'
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'PN')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=ID
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'PV')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(V(1) .GT. 10000.0D0)T1=1.0D-03
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=T1*V_MAT(ID,IT)
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  IF(V(1) .GT. 10000.0D0)THEN
	    CALL GRAMON_PGPLOT('V(Mm/s)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('V(km/s)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
!
        ELSE IF(PLT_OPT(1:2) .EQ. 'PT')THEN
          IT=NIT; ID=ND
          DO WHILE(1 .EQ. 1)
            CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
            IF(IT .EQ. 0)EXIT
            IVAR=NT
            CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
            IF(IVAR .EQ. 0)EXIT
            DO ID=1,ND
              Y(ID)=POPS(IVAR,ID,IT)
              Z(ID)=LOG10(POPS(IVAR,ID,IT))
              X(ID)=POPS(NT,ID,IT)
            END DO
            IF(LOG_Y_AXIS)THEN
              CALL DP_CURVE(ND,X,Z)
              Ylabel='Log'
            ELSE
              CALL DP_CURVE(ND,X,Y)
              Ylabel=''
            END IF
          END DO
          CALL GRAMON_PGPLOT('T(10\u4\ dK))',Ylabel,' ',' ')
          GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'PR')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(R(1) .GT. 1.0D+04)T1=1.0D-04
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=1.0D-04*R_MAT(ID,IT)
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  IF(R(1) .GT. 1.0D+04)THEN
	    CALL GRAMON_PGPLOT('R(10\u14 \dcm)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('R(10\u10 \dcm)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'VR')THEN
	  IT=NIT; ID=ND
	  NORM_R=.TRUE.
	  CALL GEN_IN(NORM_R,'Normlize R by R(ND)?')
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    T1=1.0D0; T2=1.0D0
	    IF(IT .EQ. 0)EXIT
	    IF(NORM_R)THEN
	      T1=1.0D0/R_MAT(ND,IT)
	      XLABEL='R/R(ND)'
	    ELSE
	      IF(R(1) .GT. 1.0D+04)THEN
	        T1=1.0D-04
	        XLABEL='R(10\u14 \dcm)'
	      ELSE
	        XLABEL='R(10\u10 \dcm)'
	      END IF
	    END IF
	    IF(V(1) .GT. 1.0D+04)T2=1.0D-04
	    DO ID=1,ND
	      X(ID)=T1*R_MAT(ID,IT)
	      Y(ID)=T2*V_MAT(ID,IT)
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  Ylabel='V(km/s)'
	  IF(V(1) .GT. 1.0D+04)Ylabel='V(Mm/s)'
	  CALL GRAMON_PGPLOT(XLABEL,Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'UNDO')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT; LIMITS(:)=0; LIMITS(1)=1; LIMITS(2)=ND
	  CALL GEN_IN(IT,'Primary Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IT2=IT-1
	  CALL GEN_IN(IT2,'Iteration to be merged (replaces values)',LOW_LIM=IZERO,UP_LIM=NIT)
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(LIMITS,J,NLIM_MAX,'Depths (L1 to L2, L3 to L4 etc)')
	    DO I=1,J,2
	      IF(LIMITS(I) .EQ. 0)EXIT
	      DO ID=LIMITS(I),LIMITS(I+1)
	        POPS(:,ID,IT)=POPS(:,ID,IT2)
	      END DO
	    END DO
	    IVAR=0
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'The same record is writted every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'FDGV')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT; LIMITS(:)=0; LIMITS(1)=1; LIMITS(2)=ND
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IF(IT .EQ. 0)GOTO 200
	  T1=0.0D0; T2=0.0D0
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    CALL GEN_IN(T1,'% change in variable')
	    IF(T1 .EQ. 0.0D0)CALL GEN_IN(T2,'Change in variable')
	    CALL GEN_IN(LIMITS,J,NLIM_MAX,'Depths (L1 to L2, L3 to L4 etc)')
	    DO I=1,J,2
	      IF(LIMITS(I) .EQ. 0)EXIT
	      DO ID=LIMITS(I),LIMITS(I+1)
	        T3=POPS(IVAR,ID,IT)
	        POPS(IVAR,ID,IT)=T3*(1.0D0+T1/100.0D0)+T2
	      END DO
	    END DO
	    IVAR=0
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'The same record is written every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'SFDG')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  DO K=1,16
	    DO J=1,26
	       POPS(27,K,IT)=POPS(27,K,IT)+POPS(J,K,IT)
	       POPS(J,K,IT)=POPS(J,K,IT)/100.0D0
	    END DO
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'DREP')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  CALL GEN_IN(ID1,'Depth to be used for replacement')
	  CALL GEN_IN(ID2,'ID2: Depths ID2 to ID1 pm 1 to be replaced')
	  I=1
	  IF(ID2 .LT. ID1)I=-1
	  DO K=ID1+I,ID2,I
	    DO J=1,NT
	       POPS(J,K,IT)=POPS(J,ID1,IT)
	    END DO
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'FDG')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IF(IT .EQ. 0)GOTO 200
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depth of variable')
	    IF(ID.EQ. 0)EXIT
	    WRITE(6,'(7(9X,I5))')(I,I=MAX(ID-3,1),MIN(ID+3,ND))
	    WRITE(6,'(7ES14.4E3)')(POPS(IVAR,I,IT),I=MAX(ID-3,1),MIN(ID+3,ND))
	    STRING=' '
	    DO WHILE(INDEX(STRING,'.') .EQ. 0)
	      WRITE(STRING,'(ES14.4E3)')POPS(IVAR,ID,IT); STRING=ADJUSTL(STRING)
	      CALL GEN_IN(STRING,'New value of variable - must contain decimal point')
	    END DO
	    READ(STRING,*,IOSTAT=IOS)T1
	    IF(IOS .EQ. 0)THEN
	      POPS(IVAR,ID,IT)=T1
	    ELSE
	      WRITE(6,'(A)')'Error reading value -- variable not updated.'
	    END IF
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'A new record is writted every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
! This option flips the sign of the correctons predicted by NG. Developed
! mainly for testing purposes. The flip is done in LOG space.
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'FLIP')THEN
	  IT=NIT; ID=ND; IVAR=NT
	  IREC=NIT			!IREC is updated on write
          NITSF=NITSF+1
	  DO J=1,ND
	    DO I=1,NT
 	      T1=LOG(POPS(I,J,IT-1))-LOG(POPS(I,J,IT)/POPS(I,J,IT-1))
	      POPS(I,J,IT)=EXP(T1)
	    END DO
	  END DO
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'Flip should only be called once -- it make nos sense to call it more than once'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  WRITE(6,'(A)')' Input any character to continue: '
	  READ(5,'(A)')STRING
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'SM')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT
	  ID1=1; ID2=ND/2
!
	  WRITE(6,'(A)')BLUE_PEN
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  CALL GEN_IN(ID1,'Initital depth index (inclusive)')
	  CALL GEN_IN(ID2,'Upper depth index (inclusive)')
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')' Summary of changes written to SUM_SM_CHANGES'//DEF_PEN

	  IF(ID1*ID2.EQ. 0)GOTO 200
	  ID1=MIN(ID1,ND-1); ID2=MIN(ID2,ND-1)
!
	  OPEN(UNIT=LU_OUT,FILE='SUM_SM_CHANGES',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU_OUT,'(1X,A,2X,A,7X,A,8X,A,3X,A,2X,A,2X,A)')'Depth(ID)','IVAR','Average','P(rev)',
	1               'OP(IVAR,ID)','P(IVAR,ID+1)','Log(P[ID]/P[ID+1])'
!
! Compute the average shift (log plane) with next higher depth.
! We treat each ioization stage separatly.
!
	  DO ID=MAX(ID1,ID2),MIN(ID1,ID2),-1
	    WRITE(6,*)'ID=',ID
	    DO I=1,NUM_IONS
	      CNT=0; T1=0.0D0
	      DO IVAR=ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	        T2=LOG10(POPS(IVAR,ID,IT)/POPS(IVAR,ID+1,IT))
!	        IF(ABS(T2) .LT. 5)THEN
	          T1=T1+T2
	          CNT=CNT+1
!	        END IF
	      END DO
	      WRITE(6,*)'T1=',CNT,T1,T2
	      T1=T1/CNT
!
! With the factor of 4, these statements will only change the poulations
! if they change by more than a factor of 5.
!
!	      T2=0.175D0
!	      IF(T1 .LE. 0.0D0)T1=MIN(-T2,T1)
!	      IF(T1 .GE. 0.0D0)T1=MAX(T2,T1)
!	      WRITE(6,*)IVAR,T1,T2
!	      T1=SIGN(MIN(0.175D0,T1),T1)
!
	      T2=0.7D0
	      IF(T1 .LE. 0.0D0)T1=MIN(-T2,T1)
	      IF(T1 .GE. 0.0D0)T1=MAX(T2,T1)
	      WRITE(6,*)IVAR,T1,T2
	      T1=SIGN(MIN(0.7D0,ABS(T1)),T1)
!
! If shift, for any pop, is 4 times large than average shift, replace
! the population.
!
	      WRITE(6,'(A4,4I5)')'LIMS',I,ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	      DO IVAR=ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	        T2=LOG10(POPS(IVAR,ID,IT)/POPS(IVAR,ID+1,IT))
	        WRITE(6,*)'HOPE',IVAR,T2
	        IF(ABS(T2) .GT. ABS(T1))THEN
	          T3=10**(LOG10(POPS(IVAR,ID+1,IT))+T1)
	          WRITE(6,'(I10,I6,4ES14.4E3,ES20.4E3)')ID,IVAR,T1,T3,POPS(IVAR,ID,IT),POPS(IVAR,ID+1,IT),T2
	          FLUSH(LU_OUT)
	          POPS(IVAR,ID,IT)=T3
	        END IF
	      END DO
	      WRITE(6,*)NUM_IONS,IVAR,T2
	    END DO
	  END DO
	  CLOSE(UNIT=LU_OUT)
!
! We only update NITSF for the first correction. 
!
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'The same record is written if SM is called a second time'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'INT')THEN
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; T2=100.0D0
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  CALL GEN_IN(ID,'Depth of variable')
	  CALL GEN_IN(T2,'Interpolate values with correction > >%')
	  DO IVAR=1,NT-1
	    T1=100.0D0*ABS(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	    IF(T1 .GT. T2)THEN
	      WRITE(6,*)'Replacing population for variable',IVAR
	      IF(ID .EQ. 2 .OR. ID .EQ. ND)THEN
	        POPS(IVAR,ID,IT)=POPS(IVAR,ID-1,IT)
	      ELSE IF(ID .EQ. 1)THEN
	        POPS(IVAR,ID,IT)=POPS(IVAR,ID+1,IT)
	      ELSE
	        T1=LOG(R(ID)/R(ID-1))/LOG(R(ID+1)/R(ID-1))
	        POPS(IVAR,ID,IT)=EXP( T1*LOG(POPS(IVAR,ID+1,IT)) +
	1                     (1.0D0-T1)*LOG(POPS(IVAR,ID-1,1)) )
	      END IF
	    END IF
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'FINT')THEN
	  INCLUDE 'fint.inc'
!
! This option intepolates new populations between two deths (not
! inclusive). Interpolation is only peformed for those population where
! very the large corrections were made.
! 
! Preferred option is to read in STEQ_VALS.
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'NINT')THEN
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID1=2; ID2=ND/2-1; MAX_CHANGE=1000.0D0
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  CALL GEN_IN(ID1,'Start of interpolating rangir (exclusive)')
	  CALL GEN_IN(ID2,'End of interpolating range (exclusive)')
	  CALL GEN_IN(MAX_CHANGE,'Interpolate values with FRACTIONAL correction > >')
!
	  CALL READ_CORRECTIONS(CORRECTIONS,POPS,ND,NT,NIT,LU_IN)
!
	  OPEN(UNIT=LU_OUT,FILE='NEW_COR_SUM',STATUS='UNKNOWN',ACTION='WRITE')	
	  WRITE(LU_OUT,'(A,4X,A,2X,A,6X,A,2(2X,A),4(8X,A))')' Depth','Ion','Ion','N(levs)',
	1                 'Neg.','Pos.','Ave(-ve)','Min(-ve)','Ave(+ve)','Max(+ve)'
	  WRITE(LU_OUT,'(34X,A,3X,A)')'Cnt','Cnt'
	  T2=(1.0D0-1.0D0/MAX_CHANGE)
	  DO ID=MIN(ID1+1,ID2+1),MAX(ID1-1,ID2-1)
	    DO I=1,NUM_IONS
	      CNT_NEG=0; CNT_POS=0
	      AVE_NEG_CHANGE=0.0D0; MIN_NEG_CHANGE=0.0D0
	      AVE_POS_CHANGE=0.0D0; MAX_POS_CHANGE=0.0D0
	      DO IVAR=ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	        IF(CORRECTIONS(IVAR,ID) .LT. -MAX_CHANGE)THEN
	          AVE_NEG_CHANGE=AVE_NEG_CHANGE+CORRECTIONS(IVAR,ID)
	          MIN_NEG_CHANGE=MIN(MIN_NEG_CHANGE,CORRECTIONS(IVAR,ID))
	          CNT_NEG=CNT_NEG+1
	        ELSE IF(CORRECTIONS(IVAR,ID) .GT. T2)THEN
	          AVE_POS_CHANGE=AVE_POS_CHANGE+CORRECTIONS(IVAR,ID)
	          MAX_POS_CHANGE=MAX(MAX_POS_CHANGE,CORRECTIONS(IVAR,ID))
	          CNT_POS=CNT_POS+1
	        END IF
	      END DO
	      IF(CNT_NEG .GT. 0)AVE_NEG_CHANGE=AVE_NEG_CHANGE/CNT_NEG
	      IF(CNT_POS .GT. 0)AVE_POS_CHANGE=AVE_POS_CHANGE/CNT_POS
	      WRITE(LU_OUT,'(1X,I5,3X,I4,2X,A10,3I6,4ES16.6)')ID,ION_INDEX(I),ION_ID(I),
	1                  NION(I),CNT_NEG,CNT_POS,
	1                  AVE_NEG_CHANGE,MIN_NEG_CHANGE,AVE_POS_CHANGE,MAX_POS_CHANGE	
!
	      T1=(1.0D0-1.0D0/MAX_CHANGE)
	      IF(CNT_NEG .GT. NION(I)/3 .OR. CNT_POS .GT. NION(I)/3)THEN
	        WRITE(6,'(A,A10,A,I5,A,I3)')'Replacing populations for ion ',ION_ID(I),'-',
	1             ION_INDEX(I),' at depth ',ID
	        DO IVAR=ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	          T1=LOG(R(ID)/R(ID1))/LOG(R(ID2)/R(ID1))
	          POPS(IVAR,ID,IT)=EXP( T1*LOG(POPS(IVAR,ID2,IT)) +
	1                     (1.0D0-T1)*LOG(POPS(IVAR,ID1,IT)) )
	        END DO
	      END IF
	    END DO
	  END DO
	  CLOSE(LU_OUT)
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)' '
	  WRITE(6,*)' Broad summary of changes as function of ION written to NEW_COR_SUM'
	  WRITE(6,*)' Values listed only refer to changes exceeding MAX_CHANGE'
	  WRITE(6,*)' '
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'//RED_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  WRITE(6,*)DEF_PEN
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'FIX_OSC')THEN
          FDG_COUNTER=FDG_COUNTER+1
	  CALL FIX_POP_OSCILLATIONS(POPS(1,1,NIT),R,V,SIGMA,LUSCR,ND,NT)
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'DNRG')THEN
	  LEVELS=0; T1=1.2D0; T2=0.0D0; T3=0.0D0
	  CALL GEN_IN(LEVELS,I,NLEV_MAX,'Levels for defining new r grid')
	  CALL GEN_IN(T1,'Max ratio of level populations between cons. grid points')
	  CALL DEF_NEW_RG_V1(Y,R,POPS(1,1,NIT),LEVELS,T1,T2,T3,TMP_LOG,I,NT,ND)
	  IF(TMP_LOG)THEN
	    DO J=1,I
	      TA(1:ND)=LOG10(POPS(LEVELS(J),1:ND,NIT))
	      WRITE(STRING,*)LEVELS(J); STRING='OLD '//ADJUSTL(STRING)
	      CALL DP_CURVE_LAB(ND,R,TA,STRING)
	      CALL MON_INTERP(TB,ND,IONE,Y,ND,TA,ND,R,ND)
	      WRITE(STRING,*)LEVELS(J); STRING='NEW '//ADJUSTL(STRING)
	      CALL DP_CURVE_LAB(ND,Y,TB,STRING)
	    END DO
	    CALL GRAMON_PGPLOT(' ',' ',' ',' ')
!
	    DO I=2,ND-1
	      TA(I)=(Y(I-1)-Y(I))/(Y(I)-Y(I+1))
	    END DO
	    CALL DP_CURVE(ND-2,Y(2:ND-1),TA(2:ND-1))
	    CALL GRAMON_PGPLOT('R','dR/dR',' ',' ')
	  ELSE
	    WRITE(6,*)'Construction of a new R Grid failed'
	    WRITE(6,*)'Iteration failed to converge to required number of grid points'
	    WRITE(6,'(A)',ADVANCE='NO')' Enter any character to continue: '; READ(5,*)TMP_STR
	  END IF
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'CNRG')THEN
!
	   WRITE(6,*)RED_PEN
	   WRITE(6,*)'This option assumes you have created plots using the PR option'
	   WRITE(6,*)'It also assumes that you used the NOI option so that the plots were not initialized'
	   WRITE(6,*)DEF_PEN
! 
	   I=MIN(3*ND,SIZE(TA))
	   CALL CHANGE_XAXIS_GRIDDING(TA,K,I)
!
! For SN models, we sometimes plot in uts of 10^14 cm.
! We also ensure boundary vales are absolutely correct.
!
	   IF(R(1) .GT. 1.0D+04)TA(1:K)=1.0D+04*TA(1:K)
	   TA(1)=R(1); TA(K)=R(ND)
	   OPEN(UNIT=LU_OUT,FILE='NEW_RDINR',STATUS='UNKNOWN',ACTION='WRITE')
	   WRITE(LU_OUT,'(/,2X,A,10X,A,/)')'24-FEB-2004','!Format date'
	   WRITE(LU_OUT,'(2ES14.5,3X,I5,3X,I5,/)')1.0E+05,1.0D0,1,K
	   DO I=1,K
	     WRITE(LU_OUT,'(ES20.12,6ES14.7,3X,I5)')TA(I),(1.0D0, J=1,6),I
	     WRITE(LU_OUT,'(F12.4,/)')1.0D0
	   END DO
	   CLOSE(LU_OUT)
	   GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'RAT')THEN
	  WRITE(6,*)' '
	  IT=NIT; ID=ND; T2=100.0D0
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  CALL GEN_IN(K,'Maximum depth to consider')
!	  CALL GEN_IN(T2,'Interpolate values with correction > >%')
	  DO ID=2,K
	    RMAX=0.0D0; RMIN=100.0D0
	    DO IVAR=1,NT-1
	      T3=POPS(IVAR,ID-1,IT)/POPS(IVAR,ID+1,IT)
	      IF(T3 .GT. RMAX)THEN
	        RMAX=T3; IMAX=IVAR
	      END IF
	      IF(T3 .LT. RMIN)THEN
	        RMIN=T3; IMIN=IVAR
	      END IF
	    END DO
	    WRITE(6,'(1X,A,I4,2(3X,A,I5,ES11.2E3),A)')'At depth',ID,RED_PEN,IMAX,RMAX,
	1                 BLUE_PEN,IMIN,RMIN,DEF_PEN
	  END DO
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'RNEW')THEN
	  WRITE(6,*)' '
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND
	  CALL GEN_IN(IT,'Iteration to be replaced')
	  IF(IT .GT. NIT)THEN 
	    WRITE(6,*)'Invalid iteration -- maximum is ',NIT
	    GOTO 200
	  END IF
	  CALL GEN_IN(ID,'Depth to to be replaced')
	  IF(ID .LE. 0)GOTO 200
	  IF(ID .GT. ND)THEN
	    WRITE(6,*)'Invalid depth -- maximum is ',ND
	    GOTO 200
	  END IF
	  CALL GEN_IN(STRING,'File with new POP estimates')
	  OPEN(UNIT=10,FILE=STRING,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file. IOS=',IOS
	    GOTO 200
	  END IF
	  K=1; CALL GEN_IN(K,'Column')
!
	  DO J=1,NT
	    READ(10,*)(T1,I=1,K-1),POPS(J,ID,IT)
	  END DO
!	  CALL GEN_IN(LIMITS,I,2,'Variable range')
!	  DO IVAR=LIMITS(1),LIMITS(2)
!	    POPS(IVAR,ID,IT)=POPS(IVAR,ID+1,IT)
!	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  GOTO 200
!
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'DREP')THEN
	  WRITE(6,*)' '
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND
	  CALL GEN_IN(IT,'Iteration to be replaced')
	  IF(IT .GT. NIT)THEN 
	    WRITE(6,*)'Invalid iteration -- maximum is ',NIT
	    GOTO 200
	  END IF
	  CALL GEN_IN(ID,'Depth to to be replaced')
	  IF(ID .LE. 0)GOTO 200
	  IF(ID .GT. ND)THEN
	    WRITE(6,*)'Invalid depth -- maximum is ',ND
	    GOTO 200
	  END IF
	  CALL GEN_IN(LIMITS,I,2,'Variable range')
	  DO IVAR=LIMITS(1),LIMITS(2)
	    POPS(IVAR,ID,IT)=POPS(IVAR,ID+1,IT)
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'REP')THEN
	  WRITE(6,*)' '
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT; J=1; K=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration to be replaced')
	  CALL GEN_IN(IT2,'Iteration with estimates')
	  CALL GEN_IN(J,'Minimum depth to consider')
	  CALL GEN_IN(K,'Maximum depth to consider')
	  CALL GEN_IN(IVAR,'Variable to replace')
	  DO ID=J,K
	   POPS(IVAR,ID,IT)=POPS(IVAR,ID,IT2)
	  END DO
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  GOTO 200
!
	ELSE IF(PLT_OPT .EQ. 'R' .OR.
	1       PLT_OPT .EQ. 'F' .OR.
	1       PLT_OPT .EQ. 'D' .OR.
	1       PLT_OPT .EQ. 'Y')THEN
	  ID=ND
	  DO WHILE(0 .EQ. 0)	
	    NPLTS=0
	    IVAR=1
	    DO WHILE(IVAR .NE. 0)
500	      IVAR=0
	      WRITE(STRING,'(I5,A)')NT,'](0 to plot)'
	      DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:) ; END DO
	      STRING='Variable to be plotted ['//STRING
	      CALL GEN_IN(IVAR,STRING)
	      IF(IVAR .LT. 0 .OR. IVAR .GT. NT)GO TO 500
	      IF(IVAR .EQ. 0)GOTO 1000
C
600	      CONTINUE
	      WRITE(STRING,'(I5,A)')ND,'](0 to plot)'
	      DO WHILE(STRING .EQ. ' ') ; STRING(1:)=STRING(2:); END DO
	      STRING='Depth of variable to be plotted ['//STRING
	      CALL GEN_IN(ID,STRING)
	      IF(ID .LE. 0 .OR. ID .GT. ND)GO TO 600
C
	      Y(1:NIT)=POPS(IVAR,ID,1:NIT)
C
  	      IF(PLT_OPT .EQ. 'F')THEN
  	        DO K=1,NIT-1
  	          Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)
  	          X(K)=FLOAT(K)
  	        END DO
  	        NY=NIT-1
  	        T1=MAXVAL(ABS(Z(1:NY)))
  	        IF(T1 .LT. 1.0D-02)THEN
  	          Z(1:NY)=Z(1:NY)*1.0D+03
  	          YLABEL='\gDY/Y(%)\d \ux10\u3\d'
  	          WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
  	        ELSE
	          YLABEL='\gDY/Y(%)'
	        END IF
!	      
	      ELSE IF(PLT_OPT .EQ. 'R')THEN
	        DO K=1,NIT-2
	          T1=Y(K+2)-Y(K+1)
	          T2=Y(K+1)-Y(K)
	          IF(T2 .NE. 0)THEN
	             Z(K)=T1/T2
	          ELSE
	             Z(K)=10.0
	          END IF
	          X(K)=FLOAT(K)+2
	        END DO
	        NY=NIT-2
	        YLABEL='\gDY(K+1)/\gDY(K)'
	      ELSE IF(PLT_OPT .EQ. 'D')THEN
	        DO K=1,NIT
	          Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)
	          X(K)=FLOAT(K)
	        END DO
	        NY=NIT-1
	        T1=MAXVAL(ABS(Z(1:NY)))
	        IF(T1 .LT. 1.0D-02)THEN
	          Z(1:NY)=Z(1:NY)*1.0D+03
	          YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]\d \ux10\u3\d'
	          WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
	        ELSE
	        YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]'
	        END IF
	      ELSE IF(PLT_OPT .EQ. 'Y')THEN
	        DO K=1,NIT
	          Z(K)=Y(K)
	          IF(LOG_Y_AXIS)Z(K)=LOG10(Z(K))
	          X(K)=FLOAT(K)
	        END DO
	        NY=NIT
	        YLABEL='Y(K)'
	        IF(LOG_Y_AXIS)YLABEL='Log Y(K)'
	      END IF
	      CALL DP_CURVE(NY,X,Z)
	      NPLTS=NPLTS+1
	    END DO
1000	    CONTINUE
	    IF(NPLTS .NE. 0)THEN
	      CALL GRAMON_PGPLOT('Iteration number K',Ylabel,' ',' ')
	    END IF
	    GOTO 200
	  END DO
	ELSE
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'Unrecognized command'
	  WRITE(6,'(A)')DEF_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  GOTO 200
	END IF
!
	END
