	PROGRAM PLT_LN_HEAT
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 30-Aug-2018 - All important arrays now allocatable.
!                          Cleaned up output from option MAIN.
!                          Added option IDI, and options to set X-axis.
!                          Options now done in same way as maingen.
!                          Changes done over a week (23 ? on).
! Altered:  6-Dec-2017 - Made copatible with version from osiris.
! Altered: 14-Jun-2017 - Added option to read LINEHEAT created by SOBOLEV approximation.
!
	INTEGER, PARAMETER :: NMAX=600000
	INTEGER, PARAMETER :: NS_MAX=400
!
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SCALE_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: Y(:)
	INTEGER, ALLOCATABLE:: INDX(:)
	CHARACTER(LEN=60), ALLOCATABLE :: NAME(:)
!
! These are used to determine which species contribute most to the
! difference between SE_SCL and SE_NOSCL.
!
	REAL(KIND=LDP) SUMD(NS_MAX)
	CHARACTER*10 SPEC(NS_MAX)
!
	REAL(KIND=LDP), ALLOCATABLE :: LH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SE_SCL(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SE_NOSCL(:,:)
	CHARACTER(LEN=120), ALLOCATABLE ::  TRANS_INFO(:)
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) CUR_SUM
	INTEGER ND
	INTEGER N_LINES
	INTEGER COUNT
	INTEGER I,J,K,L,ML,IBEG,IDEPTH
	INTEGER NL,NUP
	INTEGER IOS
	INTEGER NSPEC
	LOGICAL FILE_OPEN
	LOGICAL DO_SORT
	LOGICAL SOB_MODEL
	LOGICAL, ALLOCATABLE :: IS_CO_LINE(:)
	INTEGER V1,V2,DV,LOC1,LOC2,LOC3
	REAL(8), ALLOCATABLE :: CO_NET_COOL(:)
	LOGICAL DO_RAT
	LOGICAL IS_OPEN
	INTEGER DEPTH_POINT
	INTEGER LU_CO
	REAL(8), ALLOCATABLE :: CO_SCR_DATA(:,:,:)
	INTEGER MAX_ITERS
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=200) STRING
        CHARACTER(LEN=30) UC
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) CO_OPT
	CHARACTER(LEN=30) XLABEL
	CHARACTER(LEN=30) YLABEL
        EXTERNAL UC
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot the radiative equilibrium equation.'
	WRITE(6,'(A)')' Can be used to show how RE changes as we integrate from blue to red.'
	WRITE(6,'(A)')' Designed to see effect of scaling the heating rates.'
	WRITE(6,'(A)')' Can be used to identify SL assignments that might be changed.'
	WRITE(6,'(A)')' '
!
	LU_CO=42
!
        SOB_MODEL=.FALSE.
        CALL GEN_IN(SOB_MODEL,'Was LINEHEAT created using the SOBOLEV approximation?')
100	CONTINUE
 	FILENAME='LINEHEAT'
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	IF(FILENAME .EQ. ' ')GOTO 100
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
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of data points (must be exact)')
	END IF
	N_LINES=NMAX
	CALL GEN_IN(N_LINES,'Maximum number of lines to be read')
!
	ALLOCATE (LH(ND,N_LINES))
	ALLOCATE (SE_SCL(ND,N_LINES))
	ALLOCATE (SE_NOSCL(ND,N_LINES))
	ALLOCATE (TRANS_INFO(N_LINES))
!
	ALLOCATE (NU(N_LINES))
	ALLOCATE (LAM(N_LINES))
	ALLOCATE (SCALE_FAC(N_LINES))
	ALLOCATE (XV(N_LINES))
	ALLOCATE (Y(N_LINES))
	ALLOCATE (INDX(N_LINES))
	ALLOCATE (NAME(N_LINES))
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	COUNT=0
	CALL TUNE(1,'READ')
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(11,'(A)',END=5000)STRING
	  END DO
	  IF(INDEX(STRING,'error in L due to') .NE. 0)GOTO 5000
	    IF(SOB_MODEL)THEN
	    READ(STRING,*,ERR=200)LAM(ML)
	    READ(11,*)LH(1:ND,ML)
	    COUNT=COUNT+1
	  ELSE
	    STRING=ADJUSTL(STRING)
	    TRANS_INFO(ML)=STRING
	    K=INDEX(STRING,' ')
	    STRING=ADJUSTL(STRING(K:))
	    K=INDEX(STRING,' ')
	    NAME(ML)=STRING(1:K-1)
	    K=INDEX(STRING,')')
	    DO WHILE(K .NE. 0)
	      STRING(1:)=STRING(K+1:)
	      K=INDEX(STRING,')')
	    END DO
	    K=INDEX(STRING,'  ')
	    STRING(1:)=STRING(K:)
	    READ(STRING,*,ERR=200)NU(ML),NL,NUP,SCALE_FAC(ML)
	    READ(11,*)LH(1:ND,ML)
	    READ(11,'(A)')STRING
	    IF(STRING .NE. ' ')THEN
	      WRITE(6,'(A)')' '
	      WRITE(6,*)'Error: invalid data format'
	      WRITE(6,*)'Check ND values'
	      WRITE(6,*)'Current line count is',COUNT
	      WRITE(6,'(A)')' '
	      STOP
	    END IF
	    READ(11,*)SE_SCL(1:ND,ML)
	    READ(11,*)SE_NOSCL(1:ND,ML)
	    COUNT=COUNT+1
	  END IF
	END DO
5000	CONTINUE
	CALL TUNE(2,'READ')
	CALL TUNE(3,'TERMINAL')

	N_LINES=COUNT
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Number of lines read is',N_LINES
	WRITE(6,'(A)')DEF_PEN
!
        IF(SOB_MODEL)THEN
          NU(1:N_LINES)=2.99792458E+03_LDP/LAM(1:N_LINES)
        ELSE
          LAM(1:N_LINES)=2.99792458E+03_LDP/NU(1:N_LINES)
	END IF
!
! Set def axis.
!
	XV(1:N_LINES)=NU(1:N_LINES)
	XLABEL='\gn(10\u15 \dHz)'
!
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  PLT_OPT='P'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'XLAM')THEN
	    DO I=1,N_LINES
	      XV(I)=LAM(I)
	    END DO
	    XLABEL='\gl(\AA)'
	  ELSE IF(UC(PLT_OPT) .EQ. 'XNU')THEN
	    DO I=1,N_LINES
	      XV(I)=NU(I)
	    END DO
	    XLABEL='\gn(10\u15 \dHz)'
	  ELSE IF(UC(PLT_OPT) .EQ. 'XN')THEN
	    DO I=1,N_LINES
	      XV(I)=I
	    END DO
	    XLABEL='I'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LH')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=LH(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='LH'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SS')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_SCL(K,1:N_LINES)
	    WRITE(6,*)Y(N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='STEQ(scaled)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SN')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_NOSCL(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='STEQ(not scaled)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'DIFF')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    DO L=1,N_LINES-1
	      Y(L)=(SE_NOSCL(K,L+1)-SE_NOSCL(K,L))-(SE_SCL(K,L+1)-SE_SCL(K,L))
	    END DO
	    Y(N_LINES)=0
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='Diff'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'MAIN')THEN
!
! Get the influence due to scaling of a given line on the radiative equilibrium
! equation.
!
	    IDEPTH=ND
	    CALL GEN_IN(IDEPTH,'Depth for plotting')
	    SUMD=0.0_LDP; SPEC(1:NS_MAX)=' '
!
	    NSPEC=0
	    DO ML=1,N_LINES
	      Y(ML)=LH(IDEPTH,ML)*(1.0_LDP-SCALE_FAC(ML))
	      I=INDEX(NAME(ML),'(')
	      DO J=1,200
	        IF(SPEC(J) .EQ. ' ')THEN
	          SPEC(J)=NAME(ML)(1:I-1)
	          NSPEC=NSPEC+1
	        END IF
	        IF(NAME(ML)(1:I-1) .EQ. SPEC(J))THEN
	          SUMD(J)=SUMD(J)+Y(ML)
	          EXIT
	        END IF
	      END DO
	    END DO
	    CALL INDEXX(NSPEC,SUMD,INDX,.TRUE.)
!
	    WRITE(6,'(A)')RED_PEN
	    WRITE(6,'(A)')' Output is being written to the file: HEAT_CONTRIBUTIONS'
	    WRITE(6,'(A)')DEF_PEN
	    OPEN(UNIT=27,FILE='HEAT_CONTRIBUTIONS',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	      WRITE(27,'(A)')' '
	      WRITE(27,'(A,2X,I4)')' Contributions to Rad. Equil. Equ at depth',IDEPTH
	      WRITE(27,'(A)')' '
	      WRITE(27,'(X,A,ES14.4,5X,A,ES14.4)')'Scaled sum=',SE_SCL(IDEPTH,N_LINES),'Unscsaled sum',SE_NOSCL(IDEPTH,N_LINES)
	      WRITE(27,'(A)')' '
	      DO I=1,NSPEC
	        J=INDX(I)
	        WRITE(27,'(A10,ES14.4)')TRIM(SPEC(J)),SUMD(J)
	      END DO
!
	      DO_SORT=.TRUE.
	      CALL GEN_IN(DO_SORT,'Create a list of the 30 lines having the largest influence')
	      IF(DO_SORT)THEN
	        T1=SUM(SUMD)
	        Y(1:N_LINES)=ABS(Y(1:N_LINES))/T1
	        CALL INDEXX(N_LINES,Y,INDX,.TRUE.)
!
	        WRITE(27,'(A)')' '
	        WRITE(27,'(A,2X,I4)')' Lines with cooling (+ve) / heating (-ve) contributions'
	        WRITE(27,'(A)')' '
	        WRITE(27,'(1X,A,T71,A,2X,A,6X,A,4X,A)')'Transition','Lam(A)','Frac. Contr.','Cur. Sum','Scale Fac.'
	        WRITE(27,'(A)')' '
	        T1=SUM(SUMD)
	        CUR_SUM=0.0_LDP
	        DO I=N_LINES,N_LINES-29,-1
	          J=INDX(I)
	          T2=LH(IDEPTH,J)*(1.0_LDP-SCALE_FAC(J))/T1
	          CUR_SUM=CUR_SUM+T2
	          WRITE(27,'(1X,A60,ES16.6,2ES14.3,ES14.5)')NAME(J),LAM(J), LH(IDEPTH,J)*(1.0D0-SCALE_FAC(J))/T1,CUR_SUM,SCALE_FAC(J)
	        END DO
	      END IF
	    CLOSE(UNIT=27)
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SF')THEN
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SCALE_FAC(1:N_LINES)-1.0_LDP
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='SF-1.0'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDI')THEN
	    CALL GEN_IN(I,'Line index on plot')
	    WRITE(6,'(A)')TRANS_INFO(I)
!	
	  ELSE IF(UC(PLT_OPT) .EQ. 'FV')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    Y(1:ND)=SE_SCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    Y(1:ND)=SE_NOSCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    CALL GRAMON_PGPLOT('Depth Index','STEQ(sc,scl)',' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'WR')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    DO I=1,N_LINES
	      WRITE(100,'(I7,4ES14.4,3X,A)')I,NU(I),SE_NOSCL(K,I)-SE_SCL(K,I),
	1                SE_NOSCL(K,I),SE_SCL(K,I),TRIM(NAME(I))
	    END DO
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'CO')THEN
	     WRITE(6,*) 'Carbon Monoxide Line Options: '
	     WRITE(6,*) 'DV: plot all lines for a given delta v (vibrational quantum number)'
	     WRITE(6,*) 'NET: plot net CO cooling as a function of depth'
	     WRITE(6,*) 'BAND: plot net cooling from a given Band (delta v) as a function of depth'
	     WRITE(6,*) 'RATIO: plot proportion of total CO cooling coming from a particular band'
	     WRITE(6,*) 'LEVEL: plot cooling as a function of depth from lines with a given upper (v) state'
	     WRITE(6,*) 'SCR: Plot CO cooling properties over several iterations'
!	     WRITE(6,*) 'DIFF_ND: plot difference in LINEHEAT at two depths'
	     WRITE(6,*)
!
	     IF (.NOT. ALLOCATED(IS_CO_LINE)) THEN
		ALLOCATE(IS_CO_LINE(N_LINES),STAT=IOS)
		DO I=1,N_LINES
		   IF (INDEX(NAME(I),'COMI') .NE. 0) THEN
		      IS_CO_LINE(I)=.TRUE.
		   ELSE
		      IS_CO_LINE(I)=.FALSE.
		   END IF
		END DO
	     END IF
!
	     CO_OPT='NET'
	     CALL GEN_IN(CO_OPT,'Carbon Monoxide Plotting Option')
!
	     IF (UC(CO_OPT) .EQ. 'DV') THEN
		DV=1
		DO
		   CALL GEN_IN(DV,'Band (Delta V)')
		   IF (DV .LE. 12) EXIT
		   WRITE(6,*) 'DV must be no greater than 12'
		END DO
		DO I=1,N_LINES
		   IF (IS_CO_LINE(I)) THEN
		      STRING=ADJUSTL(NAME(I))
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V1
		      ELSE IF (LOC2 .EQ. 4) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V1
		      END IF
		      STRING=STRING(LOC1+10:)
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V2
		      ELSE IF (LOC2 .EQ. 4) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V2
		      END IF
		      IF (V1-V2 .EQ. DV) THEN
			 DO J=1,ND
!			    CO_BAND_COOL(J)=CO_BAND_COOL(J)+LH(J,I)
			 END DO
		      END IF
		   END IF
		END DO
!
	     ELSE IF (UC(CO_OPT) .EQ. 'NET') THEN
		Y(1:ND)=0.0D0
		DO I=1,N_LINES
		   IF (IS_CO_LINE(I)) THEN
		      DO J=1,ND
			 Y(J)=Y(J)+LH(J,I)
		      END DO
		   END IF
		END DO
		DO I=1,ND
		   XV(I)=I
		END DO
		CALL DP_CURVE(ND,XV,Y)
		CALL GRAMON_PGPLOT('Depth Index','CO Cooling',' ',' ')
!
	     ELSE IF (UC(CO_OPT) .EQ. 'BAND') THEN
		DV=1
		Y(1:ND)=0.0D0
		DO
		   CALL GEN_IN(DV,'Band (Delta V)')
		   IF (DV .LE. 12) EXIT
		   WRITE(6,*) 'DV must be no greater than 12'
		END DO
		DO I=1,N_LINES
		   IF (IS_CO_LINE(I)) THEN
		      STRING=ADJUSTL(NAME(I))
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 2) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V1
		      ELSE IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V1
		      END IF
		      STRING=STRING(LOC1+10:)
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 2) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V2
		      ELSE IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V2
		      END IF
		      IF (V1-V2 .EQ. DV) THEN
			 DO J=1,ND
			    Y(J)=Y(J)+LH(J,I)
			 END DO
		      END IF
		   END IF
		END DO
		DO I=1,ND
		   XV(I)=I
		END DO
		CALL DP_CURVE(ND,XV,Y)
		WRITE(YLABEL,*)'CO Cooling'
		WRITE(XLABEL,*)'Depth Index'
!
	     ELSE IF(UC(CO_OPT) .EQ. 'RATIO') THEN
		IF (.NOT. ALLOCATED(CO_NET_COOL)) THEN
		   ALLOCATE(CO_NET_COOL(ND),STAT=IOS)
		   CO_NET_COOL(:)=0.0D0
		END IF
		DV=1
		Y(1:ND)=0.0D0
		CO_NET_COOL(:)=0.0D0
		DO
		   CALL GEN_IN(DV,'Band (Delta V)')
		   IF (DV .LE. 12) EXIT
		   WRITE(6,*) 'DV must be no greater than 12'
		END DO
		DO I=1,N_LINES
		   IF (IS_CO_LINE(I)) THEN
		      STRING=ADJUSTL(NAME(I))
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 2) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V1
		      ELSE IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V1
		      END IF
		      STRING=STRING(LOC1+10:)
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 2) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V2
		      ELSE IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V2
		      END IF
		      DO J=1,ND
			 IF (V1-V2 .EQ. DV) THEN
			    Y(J)=Y(J)+LH(J,I)
			 END IF
			 CO_NET_COOL(J)=CO_NET_COOL(J)+LH(J,I)
		      END DO
		   END IF
		END DO
		DO I=1,ND
		   XV(I)=I
		   Y(I)=Y(I)/CO_NET_COOL(I)
		END DO
		CALL DP_CURVE(ND,XV,Y)
		WRITE(XLABEL,*)'Depth Index'
		WRITE(YLABEL,*)'CO Cooling Ratio'
	     ELSE IF(UC(CO_OPT) .EQ. 'LEVEL') THEN
!		
		DV=1
		Y(1:ND)=0.0D0
		IF (.NOT. ALLOCATED(CO_NET_COOL)) THEN
		   ALLOCATE(CO_NET_COOL(ND),STAT=IOS)
		   CO_NET_COOL(:)=0.0D0
		END IF
		CO_NET_COOL(:)=0.0D0
		DO
		   CALL GEN_IN(DV,'Vibrational Level (V)')
		   IF (DV .LE. 12) EXIT
		   WRITE(6,*) 'V must be no greater than 12'
		END DO
!
		DO_RAT=.TRUE.
		CALL GEN_IN(DO_RAT,'Plot as ratio of total CO cooling?')
		DO I=1,N_LINES
		   IF (IS_CO_LINE(I)) THEN
		      STRING=ADJUSTL(NAME(I))
		      LOC1=INDEX(STRING,'|')
		      LOC2=INDEX(STRING(LOC1+1:),'|')
		      IF (LOC2 .EQ. 2) THEN
			 READ(STRING(LOC1+1:LOC1+1),'(I1)') V1
		      ELSE IF (LOC2 .EQ. 3) THEN
			 READ(STRING(LOC1+1:LOC1+2),'(I2)') V1
		      END IF
		      DO J=1,ND
			 IF (V1 .EQ. DV) THEN
			    Y(J)=Y(J)+LH(J,I)
			 END IF
			 CO_NET_COOL(J)=CO_NET_COOL(J)+LH(J,I)
		      END DO
		   END IF
		END DO
		DO I=1,ND
		   XV(I)=I
		   IF (DO_RAT) Y(I)=Y(I)/CO_NET_COOL(I)
		END DO
		CALL DP_CURVE(ND,XV,Y)
		WRITE(XLABEL,*)'Depth Index'
		WRITE(YLABEL,*)'CO Cooling'
!
	     ELSE IF (UC(CO_OPT) .EQ. 'SCR') THEN
		INQUIRE(UNIT=LU_CO,EXIST=IS_OPEN)
		IF (IS_OPEN) THEN
		   OPEN(LU_CO,STATUS='OLD',ACTION='READ')
		ELSE
		   WRITE(*,*) 'fort.42 (CO cooling file) not found'
		END IF
		MAX_ITERS=1
		DO 
		   READ(LU_CO,'(A)',IOSTAT=IOS) STRING
		   IF (IOS .NE. 0) EXIT
		   IF (LEN(TRIM(ADJUSTL(STRING))) .GT. 1 .AND. LEN(TRIM(ADJUSTL(STRING))) .LT. 6) THEN
		      READ(STRING,'(I4)') I
		      IF (I .GT. MAX_ITERS) THEN
			 MAX_ITERS=I
		      END IF
		   END IF
		END DO
		
		ALLOCATE(CO_SCR_DATA(MAX_ITERS,ND,3),STAT=IOS)
		IF (IOS .NE. 0) THEN
		   WRITE(6,*) 'Error allocating space for data array CO_SCR_DATA'
		   EXIT
		END IF
		CO_SCR_DATA(:,:,:)=0.0D0

		REWIND(LU_CO)
		DO I=1,MAX_ITERS
		   READ(LU_CO,'(A)') STRING
		   READ(LU_CO,'(I4)') J
		   READ(LU_CO,'(A)') STRING
		   DO J=1,ND
		      READ(LU_CO,'(I4,ES12.4,ES12.4)') K,CO_SCR_DATA(I,J,2),CO_SCR_DATA(I,J,3)
		   END DO
		   READ(LU_CO,'(A)') STRING
		END DO

		DO
		   WRITE(6,*) 'I: plot CO cooling as a function of iteration at a depth of your choice'
		   WRITE(6,*) 'N: plot CO cooling as a function of depth at an iteration of your choice'
		   WRITE(6,*) 'T: plot CO cooling as a function of temperature (for a given depth)'
		   WRITE(6,*) 'E: Exit'
!
		   CO_OPT='I'
		   CALL GEN_IN(CO_OPT,'CO Scratch option')
		   IF (UC(CO_OPT) .EQ. 'I') THEN
		      I=ND
		      CALL GEN_IN(I,'Depth')
		      XV(:)=0.0D0
		      Y(:) =0.0D0
		      DO J=1,MAX_ITERS
			 XV(J)=J
			 Y(J) =CO_SCR_DATA(J,I,3)
		      END DO
		      CALL DP_CURVE(MAX_ITERS,XV,Y)
		      WRITE(XLABEL,*)'Iteration Number'
		      WRITE(YLABEL,'(A,I4,A)')'CO Cooling(depth ',I,')'
		   ELSE IF (UC(CO_OPT) .EQ. 'N') THEN
		      I=MAX_ITERS
		      CALL GEN_IN(I,'Iteration')
		      XV(:)=0.0D0
		      Y(:) =0.0D0
		      DO J=1,ND
			 XV(J)=J
			 Y(J) =CO_SCR_DATA(I,J,3)
		      END DO
		      CALL DP_CURVE(ND,XV,Y)
		      WRITE(XLABEL,*)'Depth Index'
		      WRITE(YLABEL,'(A,I4,A)')'CO Cooling (iteration ',I,')'
		   ELSE IF (UC(CO_OPT) .EQ. 'T') THEN
		      I=ND
		      CALL GEN_IN(I,'Depth')
		      XV(:)=0.0D0
		      Y(:) =0.0D0
		      DO J=1,MAX_ITERS
			 XV(J)=CO_SCR_DATA(J,I,2)
			 Y(J) =CO_SCR_DATA(J,I,3)
		      END DO
		      DO J=1,MAX_ITERS
			 L=MINLOC(XV(1:MAX_ITERS),1)
			 T1=XV(J)
			 XV(J)=XV(J+L-1)
			 XV(J+L-1)=T1
			 T2=Y(J)
			 Y(J)=Y(J+L-1)
			 Y(J+L-1)=T2
		      END DO
		      CALL DP_CURVE(MAX_ITERS,XV,Y)
		      WRITE(XLABEL,*)'Temperature (10^4 K)'
		      WRITE(YLABEL,'(A,I4,A)')'CO Cooling (depth ',I,')'
		   ELSE IF (UC(CO_OPT) .EQ. 'P') THEN
		      CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
		   ELSE IF (UC(CO_OPT) .EQ. 'E' ) THEN
		      EXIT
		   END IF
		END DO
!
	     END IF
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:    Set X-axis to line index'
	    WRITE(6,*)' XNU:   Set X-axis to frequency'
	    WRITE(6,*)' XLAM:  Set X-axis to wavelength'
	    WRITE(6,*)' '
	    WRITE(6,*)' IDI:   Get transition information for line I'
	    WRITE(6,*)' MAIN:  List most important lines contributing to heatinge error at given depth'
	    WRITE(6,*)'           Also lists contributions by ioization stage'
	    WRITE(6,*)' '
	    WRITE(6,*)' LH:    Plot LH  at given depth'
	    WRITE(6,*)' SS:    Plot STEQ (scaling) at given depth'
	    WRITE(6,*)' SN:    Plot STEQ (no scaling) at a given depth'
	    WRITE(6,*)' FV:    Plot final values (SCL, NO SCL) as a function of depth'
	    WRITE(6,*)' WR:    Write SS & SN data at a single depth to file'
	    WRITE(6,*)' DIFF   Plot contribution by individual lines to the difference between SS and SE.'
	    WRITE(6,*)' CO:    List options for plotting information about Carbon Monoxide cooling'
	    WRITE(6,*)' E(X):  Exit routine'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
