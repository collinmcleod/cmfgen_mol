!
! The spectra directory must be specified in the file DIRECTORIES.
!
! Program reads:
!              MODEL
!              RVTJ
!              OBSFLUX
! from the specified directory.
!
	PROGRAM CREATE_IF_FILE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NION_MAX=21
	INTEGER, PARAMETER :: NTS_MAX=200
!
	TYPE ION_STAGE
	  REAL*8, ALLOCATABLE :: XzV(:,:)
	  REAL*8, ALLOCATABLE :: DXzV(:)
	  REAL*8, ALLOCATABLE :: ION_POP(:)
	  REAL*8, ALLOCATABLE :: INT_DXzV(:)
	  REAL*8, ALLOCATABLE :: INT_ION_POP(:)
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
	  TYPE(ION_STAGE) ION(NION_MAX)
	END TYPE TIME_SEQ
	TYPE(TIME_SEQ) TS(NTS_MAX)
!
	REAL*8, ALLOCATABLE :: VGRID(:)
	REAL*8, ALLOCATABLE :: YOUT(:)
	REAL*8, ALLOCATABLE :: XVEC(:)
	REAL*8, ALLOCATABLE :: YVEC(:)
	REAL*8 ABUND
!
	INTEGER NMOD
	INTEGER ID,L
	INTEGER I,J,K
	INTEGER IOS
	INTEGER IMIN,IMAX
	INTEGER ND_OUT,ND_OLD
!
	LOGICAL NORM
!
	CHARACTER(LEN=80) DIR_NAME(NTS_MAX)
	CHARACTER(LEN=80) FILE_NAME
	CHARACTER(LEN=200) STRING
!
	CHARACTER(LEN=30) X_UNIT,Y_PLT_OPT
	CHARACTER(LEN=30) X_LAB,Y_LAB
	CHARACTER(LEN=2)  CHEM_SYMB
	CHARACTER(LEN=30) SPECIES
	CHARACTER(LEN=30) TIME,FORMAT_DATE
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU=10
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL LIN_INT
!
	CHARACTER(LEN=4) GEN_ION_ID(NION_MAX)
	DATA GEN_ION_ID /'MI','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
	TS(:)%SN_AGE=0.0D0
	TS(:)%ND=0
!
	FILE_NAME='DIRECTORIES'
	CALL GEN_IN(FILE_NAME,'File with list of directories')
	OPEN(UNIT=10,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
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
	CHEM_SYMB='Fe'
	LIN_INT=.FALSE.
!	CALL GEN_IN(LIN_INT,'Use linear interpolation')
!
	DO L=1,NMOD
	  TS(L)%ION(:)%PRES=L_FALSE
	  OPEN(UNIT=LU,STATUS='OLD',ACTION='READ',FILE=TRIM(DIR_NAME(L))//'MODEL')
	  DO WHILE(TS(L)%ND .EQ. 0)
	    READ(LU,'(A)')STRING
	    IF( INDEX(STRING,'Number of depth points') .NE. 0)THEN
	      READ(STRING,*)TS(L)%ND
	      WRITE(6,*)'Number of depth points is',TS(L)%ND
	    END IF
	  END DO
!
	  ALLOCATE (TS(L)%R(TS(L)%ND))
	  ALLOCATE (TS(L)%V(TS(L)%ND))
!
	  CALL RD_SING_VEC_RVTJ(TS(L)%R,TS(L)%ND,'Radius',TRIM(DIR_NAME(L))//'RVTJ',I,IOS)
	  CALL RD_SING_VEC_RVTJ(TS(L)%V,TS(K)%ND,'Velocity',TRIM(DIR_NAME(L))//'RVTJ',I,IOS)
!
	  FILE_NAME=TRIM(DIR_NAME(L))//'POP'//TRIM(SPECIES)
	  ALLOCATE (TS(L)%POP_SPEC(TS(L)%ND))
	  CALL OP_SPEC_FILE_V2(FILE_NAME,LU,ABUND,TS(L)%POP_SPEC,TS(L)%ND,FORMAT_DATE,IOS,TIME,TRIM(SPECIES))
	  DO ID=1,NION_MAX
	    READ(LU,'(A)',END=200)STRING
	    WRITE(6,'(A)')TRIM(STRING)
	    K=INDEX(STRING,':')
	    READ(STRING(K+1:),*)TS(L)%ION(ID)%NLEV
	    J=INDEX(STRING,' of ')
	    TS(L)%ION(ID)%PRES=.TRUE.
	    TS(L)%ION(ID)%ABR=STRING(J+4:INDEX(STRING,'lev')-1)
	    WRITE(6,'(A)')TS(L)%ION(ID)%ABR
	    DO I=1,NION_MAX
	       IF(TS(L)%ION(ID)%ABR .EQ. TRIM(CHEM_SYMB)//TRIM(GEN_ION_ID(I)))THEN
	         TS(L)%ION(ID)%ZXzV=I-1
	         EXIT
	       END IF
	    END DO
	    WRITE(6,'(A)')TS(L)%ION(ID)%ABR
	    WRITE(6,*)TS(L)%ION(ID)%NLEV,TS(L)%ND,TS(L)%ION(ID)%ZXzV
	    READ(LU,'(A)')STRING
	    WRITE(6,'(A)')TRIM(STRING)
	    ALLOCATE (TS(L)%ION(ID)%XzV(TS(L)%ION(ID)%NLEV,TS(L)%ND))
	    ALLOCATE (TS(L)%ION(ID)%DXzV(TS(L)%ND))
	    ALLOCATE (TS(L)%ION(ID)%ION_POP(TS(L)%ND))
	    READ(LU,*)TS(L)%ION(ID)%XzV,TS(L)%ION(ID)%DXzV
	  END DO
200	  CONTINUE
	END DO
!
	DO L=1,NMOD
	  TS(L)%MIN_ID=0
	  DO ID=1,NION_MAX
	    IF(TS(L)%ION(ID)%PRES)THEN
	       IF(TS(L)%MIN_ID .NE. 0)TS(L)%MIN_ID=ID
	       TS(L)%ION(ID)%ION_POP=SUM(TS(L)%ION(ID)%XzV,1)
	       TS(L)%MAX_ID=ID
	       WRITE(6,*)L,TS(L)%MAX_ID
	    END IF
	  END DO
	END DO
!
	ND_OUT=TS(1)%ND
	ALLOCATE (VGRID(ND_OUT),YOUT(ND_OUT))
	ALLOCATE (TS(L)%FRAC(ND_OUT))
	VGRID(:)=LOG(TS(1)%V(:))
	DO L=1,NMOD
	  IF(L .NE. 1)DEALLOCATE(XVEC,YVEC)
	  ND_OLD=TS(L)%ND
	  ALLOCATE(XVEC(ND_OLD),YVEC(ND_OLD))
	  XVEC=LOG(TS(L)%V)
!
	  DO ID=TS(L)%MIN_ID,TS(L)%MAX_ID
	    YVEC=LOG(TS(L)%ION(ID)%ION_POP)
	    CALL MON_INTERP(YOUT,ND_OUT,IONE,VGRID,ND_OUT,YVEC,ND_OLD,XVEC,ND_OLD)
	    ALLOCATE (TS(L)%ION(ID)%INT_ION_POP(ND_OUT))
	    TS(L)%ION(ID)%INT_ION_POP=EXP(YOUT)
	  END DO
	  YVEC=LOG(TS(L)%POP_SPEC)
	  CALL MON_INTERP(YOUT,ND_OUT,IONE,VGRID,ND_OUT,YVEC,ND_OLD,XVEC,ND_OLD)
	  ALLOCATE(TS(L)%INT_POP_SPEC(ND_OUT))
	  TS(L)%INT_POP_SPEC=EXP(YOUT)
!
	  ID=TS(L)%MAX_ID
	  YVEC=LOG(TS(L)%ION(ID)%DXzV)
	  CALL MON_INTERP(YOUT,ND_OUT,IONE,VGRID,ND_OUT,YVEC,ND_OLD,XVEC,ND_OLD)
	  DEALLOCATE(TS(L)%ION(ID)%DXzV(ND_OUT))
	  TS(L)%ION(ID)%DXzV=EXP(YOUT)
	END DO
	VGRID=EXP(VGRID)
!
	DO L=1,NMOD
	  DO ID=1,8
	      IF(TS(L)%ION(ID)%PRES)THEN
	        WRITE(6,*)TS(L)%ION(ID)%PRES
	        IF(ID .GE. TS(L)%MIN_ID .AND. ID .LE. TS(L)%MAX_ID)THEN
	          DO I=1,ND_OUT
	            TS(L)%FRAC(I)=TS(L)%ION(ID)%INT_ION_POP(I)/TS(L)%INT_POP_SPEC(I)
	          END DO
	        END IF
	      ELSE IF(ID .EQ. TS(L)%MAX_ID+1)THEN
	        DO I=1,ND_OUT
	          TS(L)%FRAC(I)=TS(L)%ION(ID-1)%INT_DXzV(I)/TS(L)%INT_POP_SPEC(I)
	        END DO
	      ELSE
	        TS(L)%FRAC0.0D0
	      END IF
	    END DO
	    WRITE(30,'(12ES13.4)')(TS(L)%FRAC(I),I=1,ND_OUT)
	  END DO
	END DO
!
	WRITE(6,*)'Starting ionization fraction calculaiton'
	DO I=1,ND_OUT
	  DO ID=1,8
	    DO L=1,NMOD
	      WRITE(6,'(5I6)')I,ID,L,TS(L)%ND
	      IF(TS(L)%ION(ID)%PRES)THEN
	        WRITE(6,*)TS(L)%ION(ID)%PRES
	        IF(ID .GE. TS(L)%MIN_ID .AND. ID .LE. TS(L)%MAX_ID)THEN
	          TS(L)%FRAC(I)=TS(L)%ION(ID)%INT_ION_POP(I)/TS(L)%INT_POP_SPEC(I)
	        END IF
	      ELSE IF(ID .EQ. TS(L)%MAX_ID+1)THEN
	        TS(L)%FRAC(I)=TS(L)%ION(ID-1)%INT_DXzV(I)/TS(L)%INT_POP_SPEC(I)
	      ELSE
	        TS(L)%FRAC(I)=0.0D0
	      END IF
	    END DO
!	    IF(ID .EQ. 1)THEN
	      WRITE(30,'(12ES13.4)')VGRID(I),(TS(L)%FRAC(I),L=1,10)   !NMOD)
!	    ELSE
!	    END IF
	  END DO
	  WRITE(30,'(A)')' '
	END DO
!
	STOP
	END
