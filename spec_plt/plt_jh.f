!
! General routine for plotting and comparing J (or H) obtained from the 
! EDDFACTOR (or similar) data file.
!
! Various options are available to redden and normalize the model spectra. 
! Several different units can be used for the X and Y axes.
!
	PROGRAM PLT_JH
!
! Altered 16-Jun-2000 : DIRECT_INFO call inserted.
!
! Interface routines for IO routines.
!
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
!
	IMPLICIT NONE
!
	TYPE MODEL_INTENSITY
	  INTEGER*4 NCF
	  INTEGER*4 ND
	  REAL*8, POINTER :: RJ(:,:)
	  REAL*8, POINTER :: NU(:)
	  CHARACTER*10 DATA_TYPE
	  CHARACTER*40 FILE_DATE
	  CHARACTER*80 FILENAME
	END TYPE MODEL_INTENSITY
	TYPE (MODEL_INTENSITY) ZM(5)
	INTEGER*4 ND,NCF
	INTEGER*4 NUM_FILES
	INTEGER*4 ID
!
	INTEGER*4 NCF_B
	INTEGER*4 ND_B
	REAL*8, ALLOCATABLE :: RJ_B(:,:)
	REAL*8, ALLOCATABLE :: NU_B(:)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	INTEGER*4 ND_ATM,NC_ATM,NP_ATM
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)

!
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
!
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Absisca
	CHARACTER*80 YAXIS		!Label for Ordinate
!
	REAL*8 ANG_TO_HZ
	REAL*8 KEV_TO_HZ
	REAL*8 C_CMS
!
	LOGICAL LOG_X,LOG_Y
	CHARACTER*10 Y_PLT_OPT,X_UNIT
!
	CHARACTER*6 METHOD,TYPE_ATM
	CHARACTER*10 NAME_CONVENTION
!
        INTEGER*4 ACCESS_F
        INTEGER*4, PARAMETER :: EDD_CONT_REC=3
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
        INTEGER*4 REC_SIZE
        INTEGER*4 UNIT_SIZE
        INTEGER*4 WORD_SIZE
        INTEGER*4 N_PER_REC

! Miscellaneous variables.
!
	INTEGER*4 IOS			!Used for Input/Output errors.
	INTEGER*4 I,J,K,L,ML,ISAV
	INTEGER*4 ST_REC
	INTEGER*4 REC_LENGTH
	REAL*8 SCALE_FAC
	REAL*8 T1,T2
	REAL*8 LAMC
	REAL*8 T_ELEC
	LOGICAL AIR_LAM
	LOGICAL USE_V
!
	INTEGER*4, PARAMETER :: IZERO=0
	INTEGER*4, PARAMETER :: IONE=1
	INTEGER*4, PARAMETER :: T_IN=5		!For file I/O
	INTEGER*4, PARAMETER :: T_OUT=6
	INTEGER*4, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER*4, PARAMETER :: LU_OUT=11
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER*4 GET_INDX_DP
!
! USR_OPTION variables
!
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main otion
	CHARACTER X*10			!Used for the idividual option
	CHARACTER STRING*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
	CHARACTER*80 RVTJ_FILE_NAME
!
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC
	LOGICAL EQUAL
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP
!
! 
! Set constants.
!
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
!
	C_CMS=SPEED_OF_LIGHT()
!
        CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! Set defaults.
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
	METHOD='LOGMON'
	TYPE_ATM=' '			!i.e. not Exponential at outer boundary
!
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='NAT'
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838E+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
!
!  Read in default model.
!
	ID=1
	NUM_FILES=1
	ZM(ID)%FILENAME='EDDFACTOR'
5	CALL GEN_IN(ZM(ID)%FILENAME,'First data file')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,ZM(ID)%FILE_DATE,ZM(ID)%FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening/reading INFO file: check format'
	  WRITE(T_OUT,*)'Also check eroror file or fort.2'
	  GOTO 5
	END IF
	OPEN(UNIT=LU_IN,FILE=ZM(ID)%FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error opening ',TRIM(ZM(ID)%FILENAME)
	     WRITE(T_OUT,*)'IOS=',IOS
	     GOTO 5
	  END IF
	  READ(LU_IN,REC=3)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	  ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	  ALLOCATE (ZM(ID)%RJ(ND,NCF))
	  ALLOCATE (ZM(ID)%NU(NCF))
	  DO ML=1,ZM(ID)%NCF
	    READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),ZM(ID)%NU(ML)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading all frequencies'
	      ZM(ID)%NCF=ML-1
	      EXIT
	    END IF
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file as MODEL A (default)'
	WRITE(T_OUT,*)'Number of depth points is',ZM(ID)%ND
	WRITE(T_OUT,*)'Number of frequencies is ',ZM(ID)%NCF
!
! Set default data types
!
	STRING=ZM(ID)%FILENAME
	CALL SET_CASE_UP(STRING,IZERO,IZERO)
	IF(INDEX(STRING,'EDDF') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='J'
	ELSE IF(INDEX(STRING,'FLUX') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='H'
	ELSE IF(INDEX(STRING,'FORCE') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='M(t)'
	ELSE IF(INDEX(STRING,'ETA') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='ETA'
	ELSE IF(INDEX(STRING,'CHI') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='CHI'
	ELSE
	   ZM(ID)%DATA_TYPE='UNKNOWN'
	END IF
!
!
!
! *************************************************************************
!
! Read in basic model [i.e. R, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	RVTJ_FILE_NAME='RVTJ'
	CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
	OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    GOTO 10
	  END IF
	CLOSE(LU_IN)
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND_ATM,NC_ATM,NP_ATM,RVTJ_FILE_NAME,LU_IN)
	ALLOCATE (R(ND_ATM))
	ALLOCATE (V(ND_ATM))
	ALLOCATE (SIGMA(ND_ATM))
	ALLOCATE (T(ND_ATM))
	ALLOCATE (ED(ND_ATM))
	ALLOCATE (ROSS_MEAN(ND_ATM))
	ALLOCATE (FLUX_MEAN(ND_ATM))
	ALLOCATE (POP_ATOM(ND_ATM))
	ALLOCATE (MASS_DENSITY(ND_ATM))
	ALLOCATE (POPION(ND_ATM))
	ALLOCATE (CLUMP_FAC(ND_ATM))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND_ATM,LU_IN)
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND_ATM))
	 ALLOCATE (TB(ND_ATM))
	 ALLOCATE (TC(ND_ATM))
	 ALLOCATE (TAU_ROSS(ND_ATM))
	 ALLOCATE (TAU_ES(ND_ATM))
	 IF(ROSS_MEAN(ND_ATM) .NE. 0)THEN
	   CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
	 ELSE
	  TAU_ROSS(1:ND)=0.0D0
	 END IF
	 TA(1:ND_ATM)=6.65D-15*ED(1:ND_ATM)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
!
! 
!
! This message will only be printed once
!
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'main_option.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
!
! This call resets the .sve algorithm.  Specifically it sets the next
! input answer to be a main option, and all subsequent inputs to be
! sub-options.  
!
 3	CALL SVE_FILE('RESET')
!
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	DESCRIPTION=' '					!Obvious main option
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
!
!   If the main option begins with a '.', a previously
!   written .sve file is read.
!
!   If the main option begins with a '#', a previously
!   written .box file is read.
!
!   If sve= is apended to the end of this main option, a new .sve file
!   is opened with the given name and the main option and all subsequent
!   sub-options are written to this file.
!
!   If box= is input then a .box file is created, which contains the name
!   of several .sve files to process.
!
!   If only a main option is given, the option and subsequent sub-options
!   are saved in a file called 'main option.sve'.  All following main
!   options are saved in separate files.
!
	X=UC(MAIN_OPT_STR)
	I=INDEX(X,'(')
	IF(I .NE. 0)X=X(1:I-1)		!Remove line variables
!
! 
!
	IF(X(1:3) .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ',' ')
!                    
! Set X-Ais plotting options.
!
	ELSE IF(X(1:2) .EQ.'LX' .OR. X(1:4) .EQ. 'LOGX' .OR. 
	1                            X(1:4) .EQ. 'LINX')THEN
	  LOG_X=.NOT. LOG_X
	  IF(LOG_X)WRITE(T_OUT,*)'Now using Logarithmic X axis'
	  IF(.NOT. LOG_X)WRITE(T_OUT,*)'Now using Linear X axis'
	ELSE IF(X(1:2) .EQ.'XU' .OR. X(1:6) .EQ. 'XUNITS')THEN
	  CALL USR_OPTION(X_UNIT,'X_UNIT','Ang',
	1                  'Ang, um, eV, keV, Hz, Mm/s, km/s')
	  CALL SET_CASE_UP(X_UNIT,IZERO,IZERO)
	  IF(X_UNIT .NE. 'ANG' .AND.
	1        X_UNIT .NE. 'UM' .AND.
	1        X_UNIT .NE. 'EV' .AND.
	1        X_UNIT .NE. 'KEV' .AND.
	1        X_UNIT .NE. 'HZ' .AND.
	1        X_UNIT .NE. 'MM/S' .AND.
	1        X_UNIT .NE. 'KM/S')THEN
	     WRITE(T_OUT,*)'Invalid X unit: Try again'
	   END IF
!
! NB: We offer the option to use the central frequency to avoid
! air/vacuum confusions. Model data is in vacuum wavelngths, which
! we use in plotting at all wavelngths.
!
	   IF(X_UNIT .EQ. 'MM/S' .OR. X_UNIT .EQ. 'KM/S')THEN
	     CALL USR_OPTION(LAMC,'LAMC','0.0',
	1             'Central Lambda(Ang) [-ve for frequency (10^15 Hz)]')
	     IF(LAMC .LT. 0)THEN
	       LAMC=1.0D-07*C_CMS/ABS(LAMC)
	     ELSE
	       IF(LAMC .GT. 2000)THEN
                 CALL USR_OPTION(AIR_LAM,'AIR','T',
	1                'Air wavelength [only for Lam > 2000A]?')
	       ELSE
	         AIR_LAM=.FALSE.
	       END IF
	     END IF
	     IF(AIR_LAM)LAMC=LAM_VAC(LAMC)
	   END IF
!
! Set Y axis plotting options.
!
	ELSE IF(X(1:2) .EQ. 'LY' .OR. X(1:4) .EQ. 'LOGY' .OR. 
	1                             X(1:4) .EQ. 'LINY')THEN
	  LOG_Y=.NOT. LOG_Y
	  IF(LOG_Y)WRITE(T_OUT,*)'Now using Logarithmic Y axis'
	  IF(.NOT. LOG_Y)WRITE(T_OUT,*)'Now using Linear Y axis'
	ELSE IF(X(1:2) .EQ.'YU' .OR. X(1:6) .EQ. 'YUNITS')THEN
	  CALL USR_OPTION(Y_PLT_OPT,'X_UNIT',' ',
	1          'NAT(URAL), FNU, NU_FNU, FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'NAT' .AND.
	1        Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
!
! 
!
	ELSE IF(X(1:4) .EQ. 'WRID')THEN
	  DO ID=1,NUM_FILES
	    WRITE(T_OUT,'(A,I2,A,A)')' ID=',ID,'          ',TRIM(ZM(ID)%FILENAME)
	  END DO
	ELSE IF(X(1:6) .EQ. 'RD_MOD')THEN
	  NUM_FILES=NUM_FILES+1
	  ID=NUM_FILES
	  ZM(ID)%FILENAME='EDDFACTOR'
50	  CALL GEN_IN(ZM(ID)%FILENAME,'First data file')
	  CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,ZM(ID)%FILE_DATE,ZM(ID)%FILENAME,LU_IN,IOS)
	  IF(IOS .NE. 0)GOTO 50
	  OPEN(UNIT=LU_IN,FILE=ZM(ID)%FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	     IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error opening ',TRIM(ZM(ID)%FILENAME)
	       WRITE(T_OUT,*)'IOS=',IOS
	       GOTO 50
	    END IF
	    READ(LU_IN,REC=3)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    ALLOCATE (ZM(ID)%RJ(ND,NCF))
	    ALLOCATE (ZM(ID)%NU(NCF))
	    DO ML=1,ZM(ID)%NCF
	      READ(LU_IN,REC=ST_REC+ML-1)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),ZM(ID)%NU(ML)
	    END DO
	  CLOSE(LU_IN)
	  WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file'
	  WRITE(T_OUT,*)'Number of depth points is',ZM(ID)%ND
	  WRITE(T_OUT,*)'Number of frequencies is ',ZM(ID)%NCF
!
! Set default data types
!
	  STRING=ZM(ID)%FILENAME
	  CALL SET_CASE_UP(STRING,IZERO,IZERO)
	  IF(INDEX(STRING,'EDDF') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='J'
	  ELSE IF(INDEX(STRING,'FLUX') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='H'
	  ELSE IF(INDEX(STRING,'FORCE') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='M(t)'
	  ELSE IF(INDEX(STRING,'ETA') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='ETA'
	  ELSE IF(INDEX(STRING,'CHI') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='CHI'
	  ELSE
	    ZM(ID)%DATA_TYPE='UNKNOWN'
	  END IF
!
	ELSE IF(X(1:4) .EQ. 'WSMJ')THEN
!
! Quick and dirty option to write out J on the normal size grid.
! For this to work, points must be inserted equally.
!
	  ID=1
	  K=(ZM(ID)%ND-1)/(ND_ATM-1)	!Number of points inserted/interval
	  IF( MOD(ZM(ID)%ND-1,ND_ATM-1) .EQ. 0)THEN
	    ACCESS_F=5
	    I=WORD_SIZE*(ND_ATM+1)/UNIT_SIZE; J=83
	    ZM(1)%FILE_DATE='20-Aug-2000'
	    CALL WRITE_DIRECT_INFO_V3(ND_ATM,I,ZM(1)%FILE_DATE,'J_DATA',J)
	    OPEN(UNIT=83,FILE='J_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	      WRITE(83,REC=EDD_CONT_REC)ACCESS_F,NCF,ND_ATM
	      DO ML=1,NCF
	        WRITE(83,REC=ACCESS_F-1+ML)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND,K),ZM(ID)%NU(ML)
	      END DO
	    CLOSE(UNIT=83)
	  ELSE
	    WRITE(T_OUT,*)'Error - nonuniform grid extension'
	    WRITE(T_OUT,*)'Unable to write out small J file'
	   STOP
	  END IF
!
	ELSE IF(X(1:2) .EQ. 'JD')THEN
!
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0','Scale factor to prevent overflow')
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  ISAV=I
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    I=ISAV
	    IF(ID .NE. 1 .AND. ND .NE. ZM(1)%ND)CALL USR_OPTION(I,'Depth',' ','Depth index')
	    IF(I .GT. ND)THEN
	      WRITE(T_OUT,*)'Invalid depth; maximum value is',ND
	      GOTO 1
	    END IF
	    IF(ID .EQ. 1 .AND. ND .EQ. ND_ATM)THEN
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'    R(I)/R*=',R(I)/R(ND)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       V(I)=',V(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       T(I)=',T(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'      ED(I)=',ED(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'TAU_ROSS(I)=',TAU_ROSS(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'  TAU_ES(I)=',TAU_ES(I)
	    END IF
!
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(NCF))
	    ALLOCATE (YV(NCF))
!
	    XV(1:NCF)=ZM(ID)%NU(1:NCF)
	    YV(1:NCF)=ZM(ID)%RJ(I,1:NCF)*SCALE_FAC
	    CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(NCF,XV,YV)
	  END DO
!
	ELSE IF(X(1:3) .EQ. 'JNU')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Wavelength in Ang')
	  T1=0.299794E+04/T1
!
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE',' ','Scale factor to prevent overflow')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
            I=GET_INDX_DP(T1,ZM(ID)%NU,NCF)
	    WRITE(6,*)'Index=',I,'NCF=',NCF
	    IF(ZM(ID)%NU(I)-T1 .GT. T1-ZM(ID)%NU(I+1))I=I+1
	    WRITE(6,*)'Index=',I,'NCF=',NCF
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
!
	    WRITE(6,*)'R(1)=',R(1)
	    WRITE(6,*)'R(ND)=',R(ND)
	    T2=R(ND)
	    XV(1:ND)=DLOG10(R(1:ND)/T2)
	    DO J=1,ND
	      WRITE(6,*)J,ZM(ID)%RJ(J,I)
	      IF(ZM(ID)%RJ(J,I) .GT. 0)THEN
	        YV(J)=DLOG10(ZM(ID)%RJ(J,I))
	      ELSE
	        YV(J)=-100.0
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
          END DO
!
	ELSE IF(X(1:3) .EQ. 'MT')THEN
	  CALL USR_OPTION(USE_V,'USE_V','T','Use V for x-axis (otherwise R)')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
            I=GET_INDX_DP(T1,ZM(ID)%NU,NCF)
	    IF(ZM(ID)%NU(I)-T1 .GT. T1-ZM(ID)%NU(I+1))I=I+1
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
!
	    IF(USE_V)THEN
	      XAXIS='V(km/s)'
	      XV(1:ND)=V(1:ND)
	    ELSE
	      XAXIS='R/R(ND)'
	      T2=R(ND)
	      XV(1:ND)=R(1:ND)/T2
	    END IF
	    YV(1:ND)=ZM(ID)%RJ(1:ND,NCF)
	    YAXIS='M(t)'
	    CALL DP_CURVE(ND,XV,YV)
          END DO

	ELSE IF(X(1:2) .EQ. 'CF')THEN
!
	  DO ID=1,NUM_FILES
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    IF(ALLOCATED(TA))DEALLOCATE(TA)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
	    ALLOCATE (TA(ND))
!
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    TA(1:ND)=0.0D0
	    DO ML=1,NCF-1
	      T1=0.5D0*(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))
	      DO I=1,ND
	        TA(I)=TA(I)+T1*(ZM(ID)%RJ(I,ML)+ZM(ID)%RJ(I,ML+1))
	      END DO
	    END DO
!
	   WRITE(T_OUT,*)'Boundary Luminosity is:',TA(1)*4.1274D+03
	    DO I=1,ND
	      XV(I)=I
	      YV(I)=TA(I)/TA(1)
	    END DO
	    CALL DP_CNVRT_J_V2(XV,YV,ND,LOG_X,LOG_Y,' ',Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='L/L(d=1)'
	  END DO
!
! 
! Convolve J with Electron scattering redistribution function.
!
	ELSE IF(X(1:2) .EQ. 'ES')THEN
	  CALL USR_OPTION(K,'Depth',' ','Depth index')
	  IF(K .LE. 0 .OR. K .GE. ND)THEN
	    WRITE(T_OUT,*)'Invalid depth index, ND_MAX=',ND
	    GOTO 1
	  END IF
	  ID=1; CALL GEN_IN(ID,'Which data file')
	  ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(ZM(ID)%NCF))
	  ALLOCATE (YV(ZM(ID)%NCF))
!
	  XV(1:NCF)=ZM(ID)%NU(1:NCF)
	  YV(1:NCF)=ZM(ID)%RJ(K,1:NCF)
	  CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
!
! Convolve RJ with electron-scattering redistribution function.
!
!	  DEFAULT=T(K)
!	  CALL USR_OPTION(T_ELEC,'Depth',Default,'Depth index')
	  T_ELEC=T(K)
!
	  IF(ALLOCATED(TB))DEALLOCATE(TB)	!J input
	  ALLOCATE (TB(NCF))
	  IF(ALLOCATED(TC))DEALLOCATE(TC)	!J e.s. output
	  ALLOCATE (TC(NCF))
!
	  TB(1:NCF)=ZM(ID)%RJ(K,1:NCF)
!	  IF(ONE_PAR)THEN
!	    CALL CNVLV_ES_ONE_PAR_V2(NU,TB,TC,T_ELEC,T_OUT,NCF)
!	  ELSE
	    CALL CNVLV_ES_TWO_PAR_V2(ZM(ID)%NU,TB,TC,T_ELEC,T_OUT,NCF)
!	  END IF
!
	  YV(1:NCF)=TC(1:NCF)
	  CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL DP_CURVE(NCF,XV,YV)
!
! 
! Plot section:
!
	ELSE IF(X(1:2) .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
!
	ELSE IF(X(1:4) .EQ.'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
!
! 
!
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR. 
	1                             X(1:2) .EQ. 'HE' 
	1          .OR. X(1:4) .EQ. 'HELP')THEN
	  IF(X(1:2) .EQ. 'LI')THEN
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_JH_OPT_DESC','OLD',' ','READ',
	1                  IZERO,IOS)
	  ELSE 
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_JH_OPTIONS','OLD',' ','READ',
	1                  IZERO,IOS)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening HELP or DESCRIPTER file for PLT_JH'
	    GOTO 1
	  END IF
	  READ(LU_IN,*)I,K			!For page formating (I=22,K=12)
	  DO WHILE(1.EQ. 1)
	    DO J=1,I
	      READ(LU_IN,'(A)',END=700)STRING
	      L=LEN(TRIM(STRING))
	      IF(L .EQ. 0)THEN
	        WRITE(T_OUT,'(X)')
	      ELSE
	        WRITE(T_OUT,'(X,A)')STRING(1:L)
	      END IF
	    END DO
	    READ(T_IN,'(A)')STRING
	    IF(STRING(1:1) .EQ. 'E' .OR.
	1       STRING(1:1) .EQ. 'e')GOTO 700		!Exit from listing.
	    I=K
	  END DO
700	  CONTINUE
	  CLOSE(UNIT=LU_IN)
!
	ELSE IF(X(1:2) .EQ. 'EX') THEN
	  CALL DP_CURVE(0,XV,YV)
	  STOP
	ELSE IF(X(1:3) .EQ. 'BOX') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
!
1	CONTINUE
	GO TO 3
!
	END
!
!
	FUNCTION FAC(N)
	REAL*8 FAC
	INTEGER*4 N
	INTEGER*4, PARAMETER :: T_OUT=5
!
	IF(N .EQ. 0)THEN
	  FAC=1
	ELSE IF(N .LT. 0)THEN
	  WRITE(T_OUT,*)'Error in FAC --- invalid argument'
	  STOP
	ELSE
	  FAC=1
	  DO I=2,N
	    FAC=FAC*I
	  END DO
	END IF
!
	RETURN
   	END
