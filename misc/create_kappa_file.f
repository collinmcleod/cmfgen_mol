!
! Ths program is designed to output opacities and JH (cmf) for the Ia code comparison.
!
! The directory should be specified in the file DIRECTORIES (column format).
!
! Program reads:
!              MODEL
!              RVTJ
!              VADAT
!              PLANCK_KAPPA_MEAN (if Plank mean is not on RVTJ)
!              JH_AT_CURRENT_TIME
!              JH_AT_CURRENT_TIME_INFO
! from the specified directory.
!
! If PLANCK_KAPPA_MEAN is not found, code currently outputs 1.0D-99 for the PLANCK mean.
!
! Most of the time in this routine is spent reading the JH_AT_CURRENT_TIME file.
!
	PROGRAM CREATE_KAPPA_FILE
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	USE READ_KEYWORD_INTERFACE
 	IMPLICIT NONE
!
! Created : 16-Aug-2019
!
	INTEGER, PARAMETER :: NION_MAX=20
	INTEGER, PARAMETER :: NTS_MAX=200
!
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: HTMP(:)
	REAL(KIND=LDP), ALLOCATABLE :: JNU(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: HNU(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RMID(:)
	INTEGER NCF
!
	TYPE TIME_SEQ
	  INTEGER ND
	  REAL(KIND=LDP) SN_AGE
	  REAL(KIND=LDP), ALLOCATABLE :: V(:)
	  REAL(KIND=LDP), ALLOCATABLE :: R(:)
	  REAL(KIND=LDP), ALLOCATABLE :: T(:)
	  REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	  REAL(KIND=LDP), ALLOCATABLE :: DENSITY(:)		!Mass density
	  REAL(KIND=LDP), ALLOCATABLE :: PLANCK_MEAN(:)
	  REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)
	  REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)
	  REAL(KIND=LDP), ALLOCATABLE :: JCMF(:)
	  REAL(KIND=LDP), ALLOCATABLE :: HCMF(:)
	END TYPE TIME_SEQ
	TYPE(TIME_SEQ) TS(NTS_MAX)
!
! Not currently needed.
!
	REAL(KIND=LDP), ALLOCATABLE :: VGRID(:)
	REAL(KIND=LDP), ALLOCATABLE :: YOUT(:)
	REAL(KIND=LDP), ALLOCATABLE :: XVEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: YVEC(:)
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) CHI_ROSS
	REAL(KIND=LDP) CHI_FLUX
!
	INTEGER NMOD
	INTEGER ID,L
	INTEGER I,J,K,ML
	INTEGER NX
	INTEGER IOS
	INTEGER ST_REC
	INTEGER REC_LENGTH
	INTEGER TMP_NLEV
	INTEGER BIGGEST_ID
!
	LOGICAL NORM
!
	CHARACTER(LEN=80) DIR_NAME(NTS_MAX)
	CHARACTER(LEN=80) FILE_NAME
	CHARACTER(LEN=200) STRING
	CHARACTER(LEN=30) LC,UC
	CHARACTER(LEN=30) FILE_DATE
	EXTERNAL LC,UC
!
	CHARACTER(LEN=30) X_UNIT,Y_PLT_OPT
	CHARACTER(LEN=30) X_LAB,Y_LAB
	CHARACTER(LEN=2)  CHEM_SYMB
	CHARACTER(LEN=30) SPECIES
	CHARACTER(LEN=30) TIME,FORMAT_DATE
	CHARACTER(LEN=6) TMP_ABR
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU=10
	LOGICAL FILE_OPEN
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	CHARACTER(LEN=4) GEN_ION_ID(NION_MAX)
	DATA GEN_ION_ID /'I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
	TS(:)%SN_AGE=0.0D0
	TS(:)%ND=0
	BIGGEST_ID=0
	TIME=' '
!
	FILE_NAME='DIRECTORIES'
	CALL GEN_IN(FILE_NAME,'File with list of directories')
	OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to open file with directories'
	  STOP
	END IF
	DO L=1,NTS_MAX
	  READ(10,'(A)',END=100)DIR_NAME(L)
	  NMOD=L
	  K=LEN_TRIM(DIR_NAME(L))
	  IF(DIR_NAME(L)(K:K) .NE. '/')DIR_NAME(L)(K+1:K+1)='/'
	END DO
100	CLOSE(UNIT=10)
	WRITE(6,*)'Number of directory names read is:',NMOD
!
	DO L=1,NMOD
	  WRITE(6,*)' '
	  WRITE(6,*)RED_PEN//'Reading directory ',TRIM(DIR_NAME(L))//DEF_PEN
!
! Get number of depth points.
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'MODEL'
	  OPEN(UNIT=LU,STATUS='OLD',ACTION='READ',FILE=FILE_NAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file ',TRIM(FILE_NAME)
	    STOP
	  END IF
	  DO WHILE(TS(L)%ND .EQ. 0)
	    READ(LU,'(A)')STRING
	    IF( INDEX(STRING,'Number of depth points') .NE. 0)THEN
	      READ(STRING,*)TS(L)%ND
	      WRITE(6,*)'Number of depth points is',TS(L)%ND
	    END IF
	  END DO
	  CLOSE(UNIT=LU)
!
! Get age of SN
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'VADAT'
	  CALL READ_KEYWORD(TS(L)%SN_AGE,'SN_AGE',L_TRUE,FILE_NAME,L_TRUE,L_TRUE,LU)
!
! Get atmospheric vectors. R and V are both needed. We don't need T and ED
!
	  ALLOCATE (TS(L)%R(TS(L)%ND))
	  ALLOCATE (TS(L)%V(TS(L)%ND))
	  ALLOCATE (TS(L)%T(TS(L)%ND))
	  ALLOCATE (TS(L)%ED(TS(L)%ND))
	  ALLOCATE (TS(L)%FLUX_MEAN(TS(L)%ND))
	  ALLOCATE (TS(L)%ROSS_MEAN(TS(L)%ND))
	  ALLOCATE (TS(L)%PLANCK_MEAN(TS(L)%ND))
	  ALLOCATE (TS(L)%JCMF(TS(L)%ND))
	  ALLOCATE (TS(L)%HCMF(TS(L)%ND))
	  ALLOCATE (TS(L)%DENSITY(TS(L)%ND))
!
! This reading is a little inefficient, but its nice and clean.
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'RVTJ'
	  CALL RD_SING_VEC_RVTJ(TS(L)%R,TS(L)%ND,'Radius',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%V,TS(L)%ND,'Velocity',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%T,TS(L)%ND,'Temperature',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%ED,TS(L)%ND,'Electron',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%FLUX_MEAN,TS(L)%ND,'Flux Mean Opacity',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%ROSS_MEAN,TS(L)%ND,'Rosseland Mean Opacity',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%PLANCK_MEAN,TS(L)%ND,'Planck Mean Opacity',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%DENSITY,TS(L)%ND,'Mass Density',FILE_NAME,LU,IOS)
	  TS(L)%T=1.0D+04*TS(L)%T
!
! Convert to mass-absorption coefficient.
!
	  TS(L)%ROSS_MEAN=1.0D-10*TS(L)%ROSS_MEAN/TS(L)%DENSITY
	  TS(L)%FLUX_MEAN=1.0D-10*TS(L)%FLUX_MEAN/TS(L)%DENSITY
	  TS(L)%PLANCK_MEAN=1.0D-10*TS(L)%PLANCK_MEAN/TS(L)%DENSITY
!
	  IF(TS(L)%PLANCK_MEAN(1) .EQ. 0.0D0)THEN
	    FILE_NAME=TRIM(DIR_NAME(L))//'obs_planck/'//'PLANCK_KAPPA_MEAN'
	    OPEN(UNIT=LU,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      IF(IOS .EQ. 0)THEN
	        STRING='!'
	        DO WHILE(STRING(1:1) .EQ. '!')
	          READ(LU,'(A)')STRING
	        END DO
	        BACKSPACE(LU)
	        DO I=1,TS(L)%ND
	          READ(LU,*)J,T1,T1,TS(L)%PLANCK_MEAN(I)
	        END DO
	      ELSE
	        WRITE(6,*)'Error opening file: IOSTAT=',IOS
	        WRITE(6,*)'File is ',TRIM(FILE_NAME)
	        TS(L)%PLANCK_MEAN(I)=1.0D-99
!               STOP
	      END IF
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(LU)
	  END IF
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'JH_AT_CURRENT_TIME'
	  CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILE_NAME,LU,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error opening/reading INFO file: check format'
	    WRITE(6,*)'Also check error file or fort.2'
	    WRITE(6,*)'File is ',TRIM(FILE_NAME)
	    STOP
	  END IF
	  OPEN(UNIT=LU,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error opening ',FILE_NAME
	    WRITE(6,*)'IOS=',IOS
	    STOP
	  END IF
	  READ(LU,REC=3)ST_REC,NCF,NX
	  ST_REC=ST_REC+2
	  IF(NX .NE. TS(L)%ND)THEN
	    WRITE(6,*)'Inconsistent ND in MODEL and JH_AT_CURRENT_TIME'
	    WRITE(6,*)'ND (JH_AT_CURRENT_TIME)=',NX
	    WRITE(6,*)'             ND (MODEL)=',TS(L)%ND
	    STOP
	  END IF
	  ALLOCATE (NU(NCF))
	  ALLOCATE (JNU(NX,NCF))
	  ALLOCATE (HNU(NX+1,NCF))
	  ALLOCATE (HTMP(NX+1))
	  ALLOCATE (RMID(NX+1))
	  DO ML=1,NCF
            READ(LU,REC=ST_REC+ML-1,IOSTAT=IOS)(JNU(I,ML),I=1,NX),
	1            (HNU(I,ML),I=2,NX),HNU(NX+1,ML),HNU(1,ML),NU(ML)    !H_INBC,H_OUTBC,NU(ML)
	     HNU(1,ML)=HNU(1,ML)*JNU(1,ML)
	  END DO
	  TS(L)%JCMF=0.0D0; HTMP=0.0D0
	  DO ML=1,NCF-1
	    T1=0.5D+15*(NU(ML)-NU(ML+1))
	    DO I=1,TS(L)%ND
	      TS(L)%JCMF(I)=TS(L)%JCMF(I)+T1*(JNU(I,ML)+JNU(I,ML+1))
	    END DO
	    DO I=1,TS(L)%ND+1
	      HTMP(I)=HTMP(I)+T1*(HNU(I,ML)+HNU(I,ML+1))
	    END DO
	  END DO
	  DO I=2,TS(L)%ND
	    RMID(I)=0.5D0*(TS(L)%R(I)+TS(L)%R(I-1))
	  END DO
	  RMID(1)=TS(L)%R(1); RMID(NX+1)=TS(L)%R(NX)
	  K=NX+1
	  CALL MON_INTERP(TS(L)%HCMF,NX,IONE,TS(L)%R,NX,HTMP,K,RMID,K)
!
	  TS(L)%JCMF=TS(L)%JCMF/TS(L)%R/TS(L)%R
	  TS(L)%HCMF=TS(L)%HCMF/TS(L)%R/TS(L)%R
	  DEALLOCATE(JNU,HNU,NU,RMID,HTMP)
!
	END DO
!
	WRITE(6,*)'Beginning output the models'
!
	OPEN(UNIT=30,FILE='opac_cmfgen',STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(30,'(A,I5)')'#NTIMES:  ',NMOD
	WRITE(30,'(A,100(F9.3))')'#NTIMES[d]: ',(TS(L)%SN_AGE,L=1,NMOD)
!
	DO L=1,NMOD
	  WRITE(30,'(A)')'#'
	  WRITE(30,'(A,F9.3)')'#TIME: ',TS(L)%SN_AGE
	  WRITE(30,'(A,I3)')'#NVEL: ',TS(L)%ND
	  WRITE(30,'(A)')'#vel_mid[km/s] Rossmean [cm^2/g] Fluxmean [cm^2/g] Planckmean [cm^2/g] Jcmf [erg/s/cm^2] Hcmf [erg/s/cm^2]'
	  DO I=TS(L)%ND,1,-1
	    WRITE(30,'(ES16.6,ES16.4,4ES18.4)')TS(L)%V(I),TS(L)%ROSS_MEAN(I),TS(L)%FLUX_MEAN(I),TS(L)%PLANCK_MEAN(I),
	1                                    TS(L)%JCMF(I),TS(L)%HCMF(I)
	  END DO
	END DO
	CLOSE(UNIT=30)
!
	STOP
	END
