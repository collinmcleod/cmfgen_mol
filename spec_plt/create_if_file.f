!
! The directory should be specified in the file DIRECTORIES (column format).
!
! Program reads:
!              MODEL
!              RVTJ
!              VADAT
!              POPDUM  (DUM=IRON, SIL etc)
! from the specified directory.
!
! Code outputs ionization file for one species in format sutable for code comparison.
!
	PROGRAM CREATE_IF_FILE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	USE READ_KEYWORD_INTERFACE
 	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NION_MAX=20
	INTEGER, PARAMETER :: NTS_MAX=200
!
	TYPE ION_STAGE
	  REAL*8, ALLOCATABLE :: XzV(:,:)
	  REAL*8, ALLOCATABLE :: DXzV(:)
	  REAL*8, ALLOCATABLE :: ION_POP(:)
	  REAL*8 ZXzV
	  INTEGER NLEV
	  LOGICAL PRES
	  CHARACTER(LEN=6) ABR
	END TYPE ION_STAGE
!
	TYPE TIME_SEQ
	  INTEGER ND
	  INTEGER MAX_ID
	  INTEGER MIN_ID
	  REAL*8 SN_AGE
	  REAL*8, ALLOCATABLE :: INT_POP_SPEC(:)
	  REAL*8, ALLOCATABLE :: POP_SPEC(:)
	  REAL*8, ALLOCATABLE :: FRAC(:)
	  REAL*8, ALLOCATABLE :: V(:)
	  REAL*8, ALLOCATABLE :: R(:)
	  REAL*8, ALLOCATABLE :: T(:)
	  REAL*8, ALLOCATABLE :: ED(:)
	  TYPE(ION_STAGE) ION(NION_MAX)
	END TYPE TIME_SEQ
	TYPE(TIME_SEQ) TS(NTS_MAX)
!
! Not currently needed.
!
	REAL*8, ALLOCATABLE :: VGRID(:)
	REAL*8, ALLOCATABLE :: YOUT(:)
	REAL*8, ALLOCATABLE :: XVEC(:)
	REAL*8, ALLOCATABLE :: YVEC(:)
!
	REAL*8 ABUND
!
	INTEGER NMOD
	INTEGER ID,L
	INTEGER I,J,K,IW
	INTEGER IOS
	INTEGER TMP_NLEV
	INTEGER BIGGEST_ID
!
	LOGICAL NORM
!
	CHARACTER(LEN=80) DIR_NAME(NTS_MAX)
	CHARACTER(LEN=80) FILE_NAME
	CHARACTER(LEN=200) STRING
	CHARACTER(LEN=30) LC,UC
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
	SPECIES='IRON'
	CALL GEN_IN(SPECIES,'Species -- e.g., IRON, COB, NICK etc [not case sensitive]')
	SPECIES=UC(SPECIES)
	CHEM_SYMB='Fe'
	CALL GEN_IN(CHEM_SYMB,'Chemical symbol -- e.g., Fe, H, Co, Sk etc [not case sensitive]')
	CHEM_SYMB(1:1)=UC(CHEM_SYMB(1:1)); CHEM_SYMB(2:2)=LC(CHEM_SYMB(2:2))
	WRITE(6,*)CHEM_SYMB
!
	DO L=1,NMOD
	  WRITE(6,*)' '
	  WRITE(6,*)RED_PEN//'Reading directory ',TRIM(DIR_NAME(L))//DEF_PEN
	  TS(L)%ION(:)%PRES=L_FALSE
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
	  FILE_NAME=TRIM(DIR_NAME(L))//'VADAT' 
	  CALL READ_KEYWORD(TS(L)%SN_AGE,'SN_AGE',L_TRUE,FILE_NAME,L_TRUE,L_TRUE,LU)
!
	  ALLOCATE (TS(L)%R(TS(L)%ND))
	  ALLOCATE (TS(L)%V(TS(L)%ND))
	  ALLOCATE (TS(L)%T(TS(L)%ND))
	  ALLOCATE (TS(L)%ED(TS(L)%ND))
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'RVTJ'
	  CALL RD_SING_VEC_RVTJ(TS(L)%R,TS(L)%ND,'Radius',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%V,TS(L)%ND,'Velocity',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%T,TS(L)%ND,'Temperature',FILE_NAME,LU,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%ED,TS(L)%ND,'Electron',FILE_NAME,LU,IOS)
	  TS(L)%T=1.0D+04*TS(L)%T
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'POP'//TRIM(SPECIES)
	  ALLOCATE (TS(L)%POP_SPEC(TS(L)%ND))
	  CALL OP_SPEC_FILE_V2(FILE_NAME,LU,ABUND,TS(L)%POP_SPEC,TS(L)%ND,FORMAT_DATE,IOS,TIME,TRIM(SPECIES))
	  IF(IOS .NE. 0)STOP
	  DO IW=1,NION_MAX
	    READ(LU,'(A)',END=200)STRING
	    K=INDEX(STRING,':')
	    READ(STRING(K+1:),*)TMP_NLEV
	    J=INDEX(STRING,' of ')
	    TMP_ABR=STRING(J+4:INDEX(STRING,'lev')-1)
	    DO I=1,NION_MAX
	       IF(TMP_ABR .EQ. TRIM(CHEM_SYMB)//TRIM(GEN_ION_ID(I)))THEN
	         ID=I
	         TS(L)%ION(ID)%ZXzV=I
	         EXIT
	       END IF
	    END DO
	    TS(L)%ION(ID)%ABR=TMP_ABR
	    TS(L)%ION(ID)%NLEV=TMP_NLEV
	    TS(L)%ION(ID)%PRES=L_TRUE
	    WRITE(6,'(2A)')' Chem symb=',TS(L)%ION(ID)%ABR
	    READ(LU,'(A)')STRING
	    ALLOCATE (TS(L)%ION(ID)%XzV(TS(L)%ION(ID)%NLEV,TS(L)%ND))
	    ALLOCATE (TS(L)%ION(ID)%DXzV(TS(L)%ND))
	    ALLOCATE (TS(L)%ION(ID)%ION_POP(TS(L)%ND))
	    READ(LU,*)TS(L)%ION(ID)%XzV,TS(L)%ION(ID)%DXzV
	    TS(L)%MAX_ID=ID+1
	  END DO
200	  CONTINUE
	  CLOSE(LU)
	  DO ID=1,TS(L)%MAX_ID
	    IF(TS(L)%ION(ID)%PRES)THEN
	      TS(L)%ION(ID)%ION_POP=SUM(TS(L)%ION(ID)%XzV,1)
	    ELSE IF(TS(L)%ION(ID-1)%PRES .AND. .NOT. TS(L)%ION(ID)%PRES)THEN
	      TS(L)%ION(ID)%ION_POP=TS(L)%ION(ID)%DXzV
	    END IF
	  END DO
	  BIGGEST_ID=MAX(BIGGEST_ID,TS(L)%MAX_ID)
	END DO
!

	
	WRITE(6,*)'Beginning output the models'
!
	CHEM_SYMB=LC(CHEM_SYMB)
	IF(CHEM_SYMB(2:2) .EQ. 'k')CHEM_SYMB(2:2)='i'
	OPEN(UNIT=30,FILE='if_cmgen_'//TRIM(CHEM_SYMB),STATUS='UNKNOWN')
	WRITE(30,'(A,I5)')'#NTIMES:  ',NMOD
	WRITE(30,'(A,I5)')'#NSTAGES: ',BIGGEST_ID
	WRITE(30,'(A,100(F8.2))')'#NTIME[d]: ',(TS(L)%SN_AGE,L=1,NMOD)
!
	DO L=1,NMOD
	  WRITE(30,'(A)')'#'
	  WRITE(30,'(A,F9.2)')'#TIME: ',TS(L)%SN_AGE
	  WRITE(30,'(A,I3)')'#NVEL: ',TS(L)%ND
	  WRITE(30,'(A,20(11X,A,I1))')'#vel_mid[km/s]      temp[K]        ne[/cm^3]',(CHEM_SYMB,ID,ID=1,BIGGEST_ID)
	  DO I=1,TS(L)%ND
	    WRITE(STRING,'(ES16.6,2ES14.4)')TS(L)%V(I),TS(L)%T(I),TS(L)%ED(I)
	    DO ID=1,BIGGEST_ID
	      K=LEN_TRIM(STRING)+1
	      IF(TS(L)%ION(ID)%PRES)THEN
	        WRITE(STRING(K:),'(ES14.4)')TS(L)%ION(ID)%ION_POP(I)
	      ELSE IF(TS(L)%ION(ID-1)%PRES .AND. .NOT. TS(L)%ION(ID)%PRES)THEN
	        WRITE(STRING(K:),'(ES14.4)')TS(L)%ION(ID-1)%DXzV(I)
	      ELSE
	        WRITE(STRING(K:),'(ES14.4)')1.0D-99
	      END IF
	    END DO
	    WRITE(30,'(A)')TRIM(STRING)
	  END DO
	END DO
!
	RETURN
	END
