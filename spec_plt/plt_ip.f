!
! General routine for plotting and comparing I ) obtained from IP_DATA.
!
! Various options are available to redden and normalize the model spectra. 
! Several different units can be used for the X and Y axes.
!
	PROGRAM PLT_IP
!
! Aleterd 18-Aug-2003 : IOS added to DIRECT_INFO call
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
	INTEGER*4 NCF
	INTEGER*4 ND
	INTEGER*4 NC
	INTEGER*4 NP
	REAL*8, ALLOCATABLE :: IP(:,:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: P(:)
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
	INTEGER*4 NC2,NP2
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POPTOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)

!
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: YV(:)
	REAL*4, ALLOCATABLE :: ZV(:)
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
	CHARACTER*80 FILENAME
	CHARACTER*80 FILE_DATE
!
	CHARACTER*6 METHOD,TYPETM
	CHARACTER*10 NAME_CONVENTION
!
	REAL*8 DISTANCE			!kpc
	REAL*8 SLIT_WIDTH		!arcseconds
	REAL*8 PIXEL_LENGTH		!arcseconds
!
! Miscellaneous variables.
!
	INTEGER*4 IOS			!Used for Input/Output errors.
	INTEGER*4 I,J,K,L,ML
	INTEGER*4 ST_REC
	INTEGER*4 REC_LENGTH
	REAL*8 T1,T2
	REAL*8 LAMC
	REAL*8 PI
	REAL*8 T_ELEC
	LOGICAL AIR_LAM
	LOGICAL COMPUTE_P
	LOGICAL USE_ARCSEC
	LOGICAL MULT_BY_PSQ
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
!
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC,PARSEC
	LOGICAL EQUAL
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP,PARSEC
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
	PI=ACOS(-1.0D0)
!
! Set defaults.
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
	METHOD='LOGMON'
	TYPETM=' '			!i.e. not Exponential at outer boundary
!
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='FNU'
!
	DISTANCE=1.0		!kpc
	SLIT_WIDTH=0.1          !arcseconds
	PIXEL_LENGTH=0.0254     !arcseconds
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838E+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
!
!  Read in model.
!
	FILENAME='IP_DATA'
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read IP_DATA_INFO'
	  STOP
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED')
	  READ(LU_IN,REC=3)ST_REC,NCF,NP
	  ALLOCATE (IP(NP,NCF))
	  ALLOCATE (P(NP))
	  ALLOCATE (NU(NCF))
	  IF( INDEX(FILE_DATE,'20-Aug-2000') .NE. 0)THEN
	    READ(LU_IN,REC=ST_REC)(P(I),I=1,NP)
	    ST_REC=ST_REC+1
	    COMPUTE_P=.FALSE.
	  ELSE IF( INDEX(FILE_DATE,'Unavailable') .NE. 0)THEN
	    COMPUTE_P=.TRUE.
	  ELSE
	    WRITE(T_OUT,*)'Unrecognized date when reading IP_DATA' 
	    WRITE(T_OUT,*)'Date=',FILE_DATE
	    STOP
	  END IF
	  DO ML=1,NCF
	    READ(LU_IN,REC=ST_REC+ML-1)(IP(I,ML),I=1,NP),NU(ML)
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in IP_DATA file as MODEL A (default)'
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
10	FILENAME='RVTJ'
	CALL GEN_IN(FILENAME,'File with R, V, T etc (RVTJ)')
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    GOTO 10
	  END IF
	CLOSE(LU_IN)
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND,NC2,NP2,FILENAME,LU_IN)
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
	ALLOCATE (T(ND))
	ALLOCATE (ED(ND))
	ALLOCATE (ROSS_MEAN(ND))
	ALLOCATE (FLUX_MEAN(ND))
	ALLOCATE (POPTOM(ND))
	ALLOCATE (MASS_DENSITY(ND))
	ALLOCATE (POPION(ND))
	ALLOCATE (CLUMP_FAC(ND))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POPTOM,POPION,MASS_DENSITY,CLUMP_FAC,ND,LU_IN)
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND))
	 ALLOCATE (TB(ND))
	 ALLOCATE (TC(ND))
	 ALLOCATE (TAU_ROSS(ND))
	 ALLOCATE (TAU_ES(ND))
	 IF(ROSS_MEAN(ND) .NE. 0)THEN
	   CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND,METHOD,TYPETM)
	 END IF
	 TA(1:ND)=6.65D-15*ED(1:ND)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND,METHOD,TYPETM)
!
	 IF(COMPUTE_P)THEN
	   NC=NP-ND
	   P(NC+1:NP)=R(ND:1:-1)
	   DO I=1,NC
	     P(I)=R(ND)*(I-1)/NC
	   END DO
	   WRITE(T_OUT,*)'P computed internally'
	 END IF
!
	CALL GEN_IN(DISTANCE,'Distance to star (kpc)')
	CALL GEN_IN(SLIT_WIDTH,'Width of spectrograph slit (arcsec)')
	CALL GEN_IN(PIXEL_LENGTH,'Size of pixed in spatial direction (arcsec)')
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
	1          'FNU, NU_FNU, FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
!
! 
!
! Simple option to compute inetgrated spectrum inside impact index I,
! and outside index I.
!
	ELSE IF(X(1:2) .EQ. 'SP')THEN
	  CALL USR_OPTION(I,'P',' ','Impact parameter index cuttoff')
	  IF(I .GT. NP)THEN
	    WRITE(T_OUT,*)'Invalid depth; maximum value is',NP
	    GOTO 1
	  END IF
!
	  T1=P(I)*1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	  WRITE(T_OUT,'(X,A,1P,E14.6,A)')'       P(I)=',P(I)*1.0E+10,' cm'
	  WRITE(T_OUT,'(X,A,1P,E14.6,A)')'       P(I)=',T1,' arcsec'
	  WRITE(T_OUT,'(X,A,1P,E14.6)')'    P(I)/R*=',P(I)/R(ND)
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=0.0D0
	  DO J=1,I-1
	    YV(1:NCF)=YV(1:NCF)+0.5D0*(IP(J,1:NCF)*P(J)+IP(J+1,1:NCF)*P(J+1))*(P(J+1)-P(J))
	  END DO
          T1=DISTANCE*1.0E+03*PARSEC()
	  T1=2.0D0*PI*1.0D+23*(1.0E+10/T1)**2
	  YV(1:NCF)=YV(1:NCF)*T1
!
! NB: J and I have the same units, apart from per steradian.
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  CALL CURVE(NCF,XV,YV)
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=0.0D0
	  DO J=I,NP-1
	    YV(1:NCF)=YV(1:NCF)+0.5D0*(IP(J,1:NCF)*P(J)+IP(J+1,1:NCF)*P(J+1))*(P(J+1)-P(J))
	  END DO
          T1=DISTANCE*1.0E+03*PARSEC()
	  T1=2.0D0*PI*1.0D+23*(1.0E+10/T1)**2
	  YV(1:NCF)=YV(1:NCF)*T1
!
! NB: J and I have the same units, apart from per steradian/
!
	  WRITE(6,*)XV(1),YV(1)
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL CURVE(NCF,XV,YV)
	  WRITE(6,*)XV(1),YV(1)
!
	  YAXIS(J:J)='I'
	  J=INDEX(YAXIS,')')
	  YAXIS(J:)=' Jy'
!
	ELSE IF(X(1:2) .EQ. 'IP' .OR. X(1:2) .EQ. 'FP')THEN
	  CALL USR_OPTION(I,'P',' ','Impact parameter index')
	  IF(I .GT. NP)THEN
	    WRITE(T_OUT,*)'Invalid depth; maximum value is',NP
	    GOTO 1
	  END IF
!
	  T1=P(I)*1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	  WRITE(T_OUT,'(X,A,1P,E14.6,A)')'       P(I)=',P(I)*1.0E+10,' cm'
	  WRITE(T_OUT,'(X,A,1P,E14.6,A)')'       P(I)=',T1,' arcsec'
	  WRITE(T_OUT,'(X,A,1P,E14.6)')'    P(I)/R*=',P(I)/R(ND)
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=IP(I,1:NCF)
!
! Convert to per/pixel
!
	  IF(X(1:2) .EQ. 'FP')THEN
	    YV=YV*SLIT_WIDTH*PIXEL_LENGTH*2.3504D-11
	    YAXIS='F'
	  END IF
!
! NB: J and I have the same units, apart from per steradian/
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  J=INDEX(YAXIS,'J')
	  IF(X(1:2) .EQ. 'FP')THEN
	    YAXIS(J:J)='F'
	    J=INDEX(YAXIS,')')
	    YAXIS(J:)=' pixel\u-1\d)'
	  ELSE
	    YAXIS(J:J)='I'
	    J=INDEX(YAXIS,')')
	    YAXIS(J:)=' \gW\u-1\d)'
	  END IF
	  
!
	  CALL CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:4) .EQ. 'IF2')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'Lambda',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  ALLOCATE (ZV(NP))
!
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-1
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-1)=LOG10(P(2:NP)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  YV(:)=0.0D0
	  DO K=MIN(I,J),MAX(I,J)
	    YV(1:NP)=YV(1:NP)+IP(1:NP,K)*0.5D0*(NU(K-1)-NU(K+1))
	  END DO
	  T1=ABS(NU(I)-NU(J))
	  YV(1:NP)=YV(1:NP)/T1		!Normalize so per Hz
!
	  ZV(1:NP)=0.0D0
	  DO I=1,NP-1
	    T1=0.5D0*(P(I)*YV(I)+P(I+1)*YV(I+1))*(P(I+1)-P(I))
	    ZV(I+1)=ZV(I)+T1
	  END DO
	  T1=ZV(NP)
	  DO I=2,NP
	    ZV(I-1)=ZV(I)/T1
	  END DO
          T2=DISTANCE*1.0E+03*PARSEC()
	  T2=2.0D0*PI*1.0D+23*(1.0E+10/T2)**2
	  WRITE(6,'(A,ES12.3,A)')'The average flux in band is',T1*T2,'Jy'
!
	  YAXIS='F(p)'
	  CALL CURVE(NP-1,XV,ZV)
!
	ELSE IF(X(1:4) .EQ. 'JNU2')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'Lambda',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  CALL USR_OPTION(MULT_BY_PSQ,'PSQ','F','Multiply by P^2?')
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
!
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-2
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-2)=LOG10(P(2:NP-1)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  YV(:)=0.0D0
	  DO K=MIN(I,J),MAX(I,J)
	    YV(1:NP-2)=YV(1:NP-2)+IP(2:NP-1,K)*0.5D0*(NU(K-1)-NU(K+1))
	  END DO
	  T1=ABS(NU(I)-NU(J))
	  YV(1:NP-2)=LOG10(YV(1:NP-2)/T1)
	  IF(MULT_BY_PSQ)THEN
	    DO I=1,NP-2
	     YV(I)=YV(I)+2.0D0*LOG10(P(I+1))+20.0D0
	    END DO
	    YAXIS='p\u2\dI(ergs \u-1\d Hz\u-1\d steradian\u-1\d)' 
	  ELSE
	    YAXIS='I(ergs cm\u-2\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	  END IF
	  CALL CURVE(NP-2,XV,YV)
!
	ELSE IF(X(1:3) .EQ. 'JNU')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Wavelength in Ang')
	  T1=0.299794E+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-2
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-2)=LOG10(P(2:NP-1)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  YV(1:NP-2)=DLOG10(IP(2:NP-1,I))
	  YAXIS='I(ergs cm\u-2\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	  CALL CURVE(NP-2,XV,YV)

	ELSE IF(X(1:2) .EQ. 'CF')THEN
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(TA))DEALLOCATE(TA)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  ALLOCATE (TA(NP))
!
	  TA(1:ND)=0.0D0
	  DO ML=1,NCF-1
	    T1=0.5D0*(NU(ML)-NU(ML+1))
	    DO I=1,NP
	      TA(I)=TA(I)+T1*(IP(I,ML)+IP(I,ML+1))
	    END DO
	  END DO
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-2
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-2)=LOG10(P(2:NP-1)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  DO I=1,NP
	    YV(I)=TA(I)/TA(1)
	  END DO
!	  CALL CNVRT_J(XV,YV,ND,LOG_X,LOG_Y,' ',Y_PLT_OPT,
!	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  YAXIS='I(ergs cm\u-2\d s\u-1\d steradian\u-1\d)' 
	  CALL CURVE(ND,XV,YV)
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
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_IP_OPT_DESC','OLD',' ','READ',
	1                  IZERO,IOS)
	  ELSE 
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_IP_OPTIONS','OLD',' ','READ',
	1                  IZERO,IOS)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening HELP or DESCRIPTER file for PLT_JH'
	    WRITE(6,*)' '
	    WRITE(6,*)' TIT  LX  LY  XU '
	    WRITE(6,*)' CF:   Plot cummulative spectrum as a function of impact parameter'
	    WRITE(6,*)' IP:   Plot spectrum as a function of impact parameter'
	    WRITE(6,*)' SP:   Plot spectrum inside and outside impact parameter p'
	    WRITE(6,*)' JNU:  Plot I(p) for a given frequency'
	    WRITE(6,*)' JNU2: Plot I(p) for a given frequency band'
	    WRITE(6,*)' IF2:  Plot normalize Flux originating inside p for a given frequency band'
	    WRITE(6,*)' '
	    WRITE(6,*)' '
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
