!
! Simple program to plot a set of CMFGEN SN spectra for a model sequence.
! The spectra directory must be specified in the file DIRECTORIES.
!
! The program also plots the comoving and observer's frame luminosity (at the outer boundary).
!
! Program reads:
!              VADAT
!              RVTJ
!              OBSFLUX
! from the specified directory.
!
	PROGRAM PLT_MANY_SN_SPEC
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NMAX=200
	INTEGER, PARAMETER :: NCF_MAX=500000
!
	REAL(10) OBS_LUM(NMAX)
	REAL(10) CMF_LUM(NMAX)
	REAL(10) SN_AGE(NMAX)
	REAL(10) VINF(NMAX)
	INTEGER ND(NMAX)
!
	REAL(10) OBS_FREQ(NCF_MAX)
	REAL(10) OBS_FLUX(NCF_MAX)
	REAL(10) NU_CONT(NCF_MAX)
	REAL(10) FLUX_CONT(NCF_MAX)
	REAL(10) WORK(NCF_MAX)
	REAL(10) LAMC
	REAL(10) T1,T2
	REAL(10) LUM_CONV_FACTOR
	REAL(10) QW
	REAL(10) MAXIMUM,MINIMUM
	REAL(10) XEMIS_MEAN
	REAL(10) HA_WAVE
	REAL(10) C_Mms
	REAL(10) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	REAL(10) HA_VEL_EST
	REAL(10) HB_VEL_EST
	REAL(10) FeII_5170_VEL_EST
	REAL(10) FeII_5020_VEL_EST
!
! For plotting.
!
	REAL*4 XV(NCF_MAX)
	REAL*4 YV(NCF_MAX)
!
	INTEGER NCF
	INTEGER NCF_CONT
	INTEGER NMOD
	INTEGER I,J,K
	INTEGER IOS
	INTEGER IMIN,IMAX
!
	LOGICAL NORM
!
	CHARACTER(LEN=80) DIR_NAME(NMAX)
	CHARACTER(LEN=80) FILE_NAME
	CHARACTER(LEN=50) OBS_SPEC_FILE_NAME
	CHARACTER(LEN=50) OBS_CONT_FILE_NAME
	CHARACTER(LEN=200) STRING
!
	CHARACTER(LEN=30) X_UNIT,Y_PLT_OPT
	CHARACTER(LEN=30) X_LAB,Y_LAB
!
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL LIN_INT
!
	OBS_LUM=0.0D0; CMF_LUM=0.0D0
	SN_AGE=0.0D0; VINF=0.0D0
	ND=0
	X_UNIT='Ang'
	Y_PLT_OPT='FNU'
	HA_WAVE=6564.55326D0
	C_MMS=1.0D-08*SPEED_OF_LIGHT()
	HA_VEL_EST=-10.0D0
!
	FILE_NAME='DIRECTORIES'
	CALL GEN_IN(FILE_NAME,'File with list of directories')
	OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	DO I=1,NMAX
	  READ(10,'(A)',END=100)DIR_NAME(I)
	  NMOD=I
	  K=LEN_TRIM(DIR_NAME(I))
	  IF(DIR_NAME(I)(K:K) .NE. '/')DIR_NAME(I)(K+1:K+1)='/'
	END DO
100	CLOSE(UNIT=10)
	WRITE(6,*)'Number of directory names read is:',NMOD
!
	OBS_SPEC_FILE_NAME='obs_zf/obs_fin'
	CALL GEN_IN(OBS_SPEC_FILE_NAME,'Obs specrum filename: OBSFLUX')
	OBS_CONT_FILE_NAME='obs_cont/obs_cont'
	CALL GEN_IN(OBS_CONT_FILE_NAME,'Obs specrum filename: OBSFLUX')
	LIN_INT=.FALSE.
	CALL GEN_IN(LIN_INT,'Use linear interpolation')
	NORM=.TRUE.
	IF(OBS_CONT_FILE_NAME .EQ. ' ')NORM=.FALSE.
	CALL GEN_IN(NORM,'Normalize by continuum spectrum?')
!
	DO I=1,NMOD
	  FILE_NAME=TRIM(DIR_NAME(I))//'OUTGEN'
	  OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	  DO WHILE(1 .EQ. 1)
	    READ(10,'(A)',END=200)STRING
	    IF(INDEX(STRING,'Luminosity of star') .NE. 0)THEN
	      K=INDEX(STRING,':')+1
	      READ(STRING(K:),*)CMF_LUM(I)
	    END IF
	  END DO
200	  CLOSE(UNIT=10)
	END DO
	WRITE(6,*)'Read in CMF luminosities'
!
	DO I=1,NMOD
	  FILE_NAME=TRIM(DIR_NAME(I))//'VADAT'
	  OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	  DO WHILE(SN_AGE(I) .EQ. 0)
	    READ(10,'(A)')STRING
	    IF(INDEX(STRING,'[SN_AGE]') .NE. 0)THEN
	      READ(STRING,*)SN_AGE(I)
	      CLOSE(UNIT=10)
	    END IF
	  END DO
	END DO
	WRITE(6,*)'Read in SN ages'
!
	DO I=1,NMOD
	  WRITE(6,*)TRIM(FILE_NAME)
	  FILE_NAME=TRIM(DIR_NAME(I))//'RVTJ'
	  OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	  DO WHILE(ND(I) .EQ. 0)
	    READ(10,'(A)')STRING
	    IF(INDEX(STRING,'ND:') .NE. 0)THEN
	      K=INDEX(STRING,':')+1
	      READ(STRING(K:),*)ND(I)
	    END IF
	  END DO
	  DO WHILE(VINF(I) .EQ. 0)
	    READ(10,'(A)')STRING
	    IF(INDEX(STRING,'Velocity') .NE. 0)THEN
	      READ(10,*)VINF(I)
	      CLOSE(UNIT=10)
	    END IF
	  END DO
	END DO
	WRITE(6,*)'Read in SN terminal velocities'
	WRITE(6,'(A,T20,A,4X,A,3X,A)')'Model','Age(days)','dlogt','Vinf(kms)'
	DO I=1,NMOD
	  WRITE(6,'(A,T20,F9.3,3X,F6.3,4X,F8.1)')TRIM(DIR_NAME(I)),SN_AGE(I),
	1            SN_AGE(I)/SN_AGE(MAX(I-1,1)),VINF(I)
	END DO
	WRITE(6,'(A)',ADVANCE='NO')'Input any character to continue:'
	READ(5,'(A)')STRING
!
	CALL GEN_IN(X_UNIT,'Ang, Hz,  keV, um (not case sensitive)')
	CALL SET_CASE_UP(X_UNIT,IZERO,IZERO)
	CALL GEN_IN(Y_PLT_OPT,'FNU, NUFNU, FLAM')
	CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	LUM_CONV_FACTOR=16*ATAN(1.0D0)*1.0D+15*1.0D-23*((3.0856D+21)**2)/3.826D+33
	DO I=1,NMOD
	  FILE_NAME=TRIM(DIR_NAME(I))//TRIM(OBS_SPEC_FILE_NAME)
	  CALL RD_MOD(OBS_FREQ,OBS_FLUX,NCF_MAX,NCF,FILE_NAME,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open ',TRIM(FILE_NAME)
	    STOP
	  END IF
	  OBS_LUM(I)=0.0D0
	  DO K=2,NCF
	    OBS_LUM(I)=OBS_LUM(I)+(OBS_FREQ(K-1)-OBS_FREQ(K))*(OBS_FLUX(K)+OBS_FLUX(K+1))
	  END DO
	  OBS_LUM(I)=LUM_CONV_FACTOR*OBS_LUM(I)*0.5D0
	  IF(NORM)THEN
	    FILE_NAME=TRIM(DIR_NAME(I))//TRIM(OBS_CONT_FILE_NAME)
	    CALL RD_MOD(NU_CONT,FLUX_CONT,NCF_MAX,NCF_CONT,FILE_NAME,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to open ',TRIM(FILE_NAME)
	      STOP
	    END IF
	    CALL DIVIDE_BY_CONT(WORK,OBS_FREQ,OBS_FLUX,NCF,NU_CONT,FLUX_CONT,NCF_CONT,LIN_INT)
	    XV(1:NCF)=OBS_FREQ(1:NCF); YV(1:NCF)=WORK(1:NCF)
	    CALL CNVRT(XV,YV,NCF,L_FALSE,L_FALSE,X_UNIT,' ',LAMC,X_LAB,Y_LAB,L_TRUE)
!
	    IF(I .EQ. 1)OPEN(UNIT=20,FILE='HA_VEL_DATA',STATUS='UNKNOWN')
	    CALL GET_PROF_PARAMS_V2(XV,YV,NCF,HA_WAVE,HA_VEL_EST,L_TRUE,SN_AGE(I),DIR_NAME(I),20)
	    IF(I .EQ.NMOD)CLOSE(UNIT=20)
!
	    T1=4862.691D0
	    IF(I .EQ. 1)HB_VEL_EST=HA_VEL_EST
	    IF(I .EQ. 1)OPEN(UNIT=21,FILE='HB_VEL_DATA',STATUS='UNKNOWN')
	    CALL GET_PROF_PARAMS_V2(XV,YV,NCF,T1,HB_VEL_EST,L_FALSE,SN_AGE(I),DIR_NAME(I),21)
	    IF(I .EQ.NMOD)CLOSE(UNIT=21)
!
	    T1=5170.468D0 
	    IF(I .EQ. 1)FEII_5170_VEL_EST=HA_VEL_EST
	    IF(I .EQ. 1)OPEN(UNIT=22,FILE='FeII_5170_VEL_DATA',STATUS='UNKNOWN')
	    CALL GET_PROF_PARAMS_V2(XV,YV,NCF,T1,FeII_5170_VEL_EST,L_FALSE,SN_AGE(I),DIR_NAME(I),22)
	    IF(I .EQ.NMOD)CLOSE(UNIT=22)
!
	    T1=5019.836D0 
	    IF(I .EQ. 1)FEII_5020_VEL_EST=HA_VEL_EST
	    IF(I .EQ. 1)OPEN(UNIT=23,FILE='FeII_5020_VEL_DATA',STATUS='UNKNOWN')
	    CALL GET_PROF_PARAMS_V2(XV,YV,NCF,T1,FeII_5020_VEL_EST,L_FALSE,SN_AGE(I),DIR_NAME(I),23)
	    IF(I .EQ.NMOD)CLOSE(UNIT=23)
!
	  ELSE
	    XV(1:NCF)=OBS_FREQ(1:NCF); YV(1:NCF)=OBS_FLUX(1:NCF)
	    CALL CNVRT(XV,YV,NCF,L_FALSE,L_FALSE,X_UNIT,Y_PLT_OPT,LAMC,X_LAB,Y_LAB,L_FALSE)
	  END IF
	  CALL CURVE(NCF,XV,YV)
	END DO
	CLOSE(UNIT=20)
	CALL GRAMON_PGPLOT(X_LAB,Y_LAB,' ',' ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')'The red curve show the luminosity in the comoving-frame'//BLUE_PEN
	WRITE(6,'(A)')'The blue curve show the luminosity in the observers'' frame'
	WRITE(6,'(A)')DEF_PEN
!
	CALL DP_CURVE(NMOD,SN_AGE,CMF_LUM)
	CALL DP_CURVE(NMOD,SN_AGE,OBS_LUM)
	CALL GRAMON_PGPLOT('Age(days)','Luminosity(Lsun)',' ',' ')
!
	T1=0.0D0; T2=0.0D0
	DO I=2,NMOD
	  T1=T1+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*OBS_LUM(I)+SN_AGE(I-1)*OBS_LUM(I-1))
	  T2=T2+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*CMF_LUM(I)+SN_AGE(I-1)*CMF_LUM(I-1))
	END DO
	T1=T1*0.5D0*3.826D+33*(24.0D0*3600.0D0)**2
	T2=T2*0.5D0*3.826D+33*(24.0D0*3600.0D0)**2
	WRITE(6,'(A,ES14.4,A)')'Int. t.L(t) dt (Obs frame) is',T1,' s. ergs'
	WRITE(6,'(A,ES14.4,A)')'Int. t.L(t) dt (CMF frame) is',T2,' s. ergs'
	STOP
	END
