!
! Program to average departure coefficient's from two different
! input files. Designed to create intermediate SN time dependent
! models.
!
	PROGRAM AVE_DC
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: T_IN=5		!Terminal input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal output
	INTEGER, PARAMETER :: DCIN=10
	INTEGER, PARAMETER :: DCOUT=11
	INTEGER, PARAMETER :: LUIN=15
	INTEGER, PARAMETER :: LUOUT=16
!
	INTEGER N,ND
	REAL*8, ALLOCATABLE :: R_NEW(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: DION(:)
	REAL*8, ALLOCATABLE :: DC(:,:)
!
	INTEGER N1,ND1
	REAL*8, ALLOCATABLE :: R1(:)
	REAL*8, ALLOCATABLE :: V1(:)
	REAL*8, ALLOCATABLE :: ED1(:)
	REAL*8, ALLOCATABLE :: T1(:)
	REAL*8, ALLOCATABLE :: DION1(:)
	REAL*8, ALLOCATABLE :: DC1(:,:)
!
	INTEGER N2,ND2
	REAL*8, ALLOCATABLE :: R2(:)
	REAL*8, ALLOCATABLE :: V2(:)
	REAL*8, ALLOCATABLE :: ED2(:)
	REAL*8, ALLOCATABLE :: T2(:)
	REAL*8, ALLOCATABLE :: DION2(:)
	REAL*8, ALLOCATABLE :: DC2(:,:)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
!
	REAL*8, ALLOCATABLE :: XV1(:)
	REAL*8, ALLOCATABLE :: XV2(:)
	REAL*8, ALLOCATABLE :: XNEW(:)
!
	INTEGER, PARAMETER :: NMAX=1000
	REAL*8 EINA(NMAX,NMAX)
	REAL*8 EDGE(NMAX)
	REAL*8 STAT_WT(NMAX)
	REAL*8 ARAD(NMAX)
	REAL*8 GAM2(NMAX)
	REAL*8 GAM4(NMAX)
        LOGICAL  OBSERVED_LEVEL(NMAX)       !If true, level energy is known.
        CHARACTER(LEN=30)  LEVNAME(NMAX)        !
!
        REAL*8 IONIZATION_ENERGY        !Ionization energy of ion (cm^-1}
        REAL*8 ZION                     !Charge on core (i.e. ion charge +1)
	CHARACTER(LEN=11) OSCDATE       !Date oscillator file was written.
	REAL*8  GF_CUT
	INTEGER NTRET
	INTEGER LEV_CUT
	INTEGER MIN_NUM_TRANS
	CHARACTER(LEN=10) GF_ACTION
	LOGICAL   ONLY_OBS_LINES
!
	REAL*8 LSTAR
	REAL*8 RSTAR
	CHARACTER(LEN=200) STRING
!
	REAL*8 SCALE_FAC
	REAL*8 X1
	REAL*8 DELTA_T
	REAL*8 JNK
	REAL*8 FX
	REAL*8 T_EXCITE
!
	REAL*8, PARAMETER :: HDKT=4.7994145D0
	INTEGER I,J,K
	INTEGER IOS
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	CHARACTER(LEN=80) DIRECTORY_ONE
	CHARACTER(LEN=80) DIRECTORY_TWO
	CHARACTER(LEN=80) DC_FILE
!
	WRITE(6,*)' '
	WRITE(6,*)'This program is still under development, and like all porgrams'
	WRITE(6,*)'should be used with caution.'
	WRITE(6,*)' '
	WRITE(6,*)' Input information for this prorgam should be placed in DC_FILE_INFO.'
	WRITE(6,*)' Format should be as follows (with no blank lines):'
	WRITE(6,*)' '
	WRITE(6,*)'1st Directory'
	WRITE(6,*)'2nd Directory'
	WRITE(6,*)'DC file 1 (assumed same name for both directories)'
	WRITE(6,*)'DC file 2 (assumed same name for both directories)'
	WRITE(6,*)'etc'
	WRITE(6,*)' '

	DC_FILE=' '
	OPEN(UNIT=20,FILE='DC_FILE_INFO',STATUS='OLD',ACTION='READ')
	READ(20,'(A)')DIRECTORY_ONE
	READ(20,'(A)')DIRECTORY_TWO
!
	DO K=1,100
!
	  READ(20,'(A)',END=1000)DC_FILE
!
	  OPEN(UNIT=DCIN,FILE=TRIM(DIRECTORY_ONE)//TRIM(DC_FILE),STATUS='OLD',ACTION='READ')
!
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(DCIN,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)THEN
	     REWIND(DCIN)
	     STRING=' '
	  END IF
	  READ(DCIN,*)RSTAR,LSTAR,N1,ND1
!
	  IF(.NOT. ALLOCATED(R1))THEN
	    ALLOCATE (R1(ND1), DION1(ND1), ED1(ND1), T1(ND1))
	    ALLOCATE (V1(ND1),DC1(N1,ND1))
	  ELSE
	    DEALLOCATE (DC1)
	    ALLOCATE (DC1(N1,ND1))
	  END IF
	  DO I=1,ND1
	    READ(DCIN,*)R1(I),DION1(I),ED1(I),T1(I),JNK,V1(I)
	    READ(DCIN,*)(DC1(J,I),J=1,N1)
	  END DO
!
	CLOSE(UNIT=DCIN)
!
	OPEN(UNIT=DCIN,FILE=TRIM(DIRECTORY_TWO)//TRIM(DC_FILE),STATUS='OLD',ACTION='READ')
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(DCIN,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)THEN
	     REWIND(DCIN)
	     STRING=' '
	  END IF
	  READ(DCIN,*)RSTAR,LSTAR,N2,ND2
!
	  IF(.NOT. ALLOCATED(R2))THEN
	    ALLOCATE (R2(ND2), DION2(ND2), ED2(ND2), T2(ND2))
	    ALLOCATE (V2(ND2),DC2(N2,ND2))
	  ELSE
	    DEALLOCATE (DC2)
	    ALLOCATE (DC2(N2,ND2))
	  END IF
	  DO I=1,ND2
	    READ(DCIN,*)R2(I),DION2(I),ED2(I),T2(I),JNK,V2(I)
	    READ(DCIN,*)(DC2(J,I),J=1,N2)
	  END DO
!
	CLOSE(UNIT=DCIN)
!
	GF_CUT=1.0E+06			!Only want energy levels
	I=INDEX(DC_FILE,'OUT')-1
	CALL GENOSC_V8(EINA,EDGE,STAT_WT,LEVNAME,
	1         ARAD,GAM2,GAM4,OBSERVED_LEVEL,
	1         IONIZATION_ENERGY,ZION,
	1         OSCDATE,N1,NTRET,
	1         GF_ACTION,GF_CUT,LEV_CUT,MIN_NUM_TRANS,ONLY_OBS_LINES,
	1         LUIN,LUOUT,DC_FILE(1:I)//'_F_OSCDAT')
!
	DO I=1,N1
	  DO J=ND1,1,-1
	    DELTA_T=100
	    T_EXCITE=T1(J)
	    DO WHILE(ABS(DELTA_T) .GT. 1.0E-08)
	      FX=DC1(I,J)*EXP(HDKT*EDGE(I)*(1.0D0/T1(J)-1.0D0/T_EXCITE))*
	1           (T_EXCITE/T1(J))**1.5D0
	      DELTA_T=(FX-1.0D0)*T_EXCITE/FX/(1.5D0+HDKT*EDGE(I)/T_EXCITE)
	      IF(DELTA_T .GT.  0.8D0*T_EXCITE)DELTA_T=0.8D0*T_EXCITE
	      IF(DELTA_T .LT. -0.8D0*T_EXCITE)DELTA_T=-0.8D0*T_EXCITE
	      T_EXCITE=T_EXCITE-DELTA_T
	    END DO
	    DC1(I,J)=T_EXCITE
	  END DO
	END DO
!
	DO I=1,N2
	  DO J=ND2,1,-1
	    DELTA_T=100
	    T_EXCITE=T2(J)
	    DO WHILE(ABS(DELTA_T) .GT. 1.0E-08)
	      FX=DC2(I,J)*EXP(HDKT*EDGE(I)*(1.0D0/T2(J)-1.0D0/T_EXCITE))*
	1           (T_EXCITE/T2(J))**1.5D0
	      DELTA_T=(FX-1.0D0)*T_EXCITE/FX/(1.5D0+HDKT*EDGE(I)/T_EXCITE)
	      IF(DELTA_T .GT.  0.8D0*T_EXCITE)DELTA_T=0.8D0*T_EXCITE
	      IF(DELTA_T .LT. -0.8D0*T_EXCITE)DELTA_T=-0.8D0*T_EXCITE
	      T_EXCITE=T_EXCITE-DELTA_T
	    END DO
	    DC2(I,J)=T_EXCITE
	  END DO
	END DO
!
	SCALE_FAC=R2(ND2)/R1(ND1)
	ND=ND1; N=N1
	X1=0.5D0
!
	ALLOCATE (R_NEW(ND1))
	ALLOCATE (XV1(ND1),XV2(ND2),XNEW(ND))
	I=MAX(ND1,ND2)
	ALLOCATE (TA(I),TB(I),TC(I))
	ALLOCATE (ED(ND),V(ND),T(ND),DION(ND))
	ALLOCATE (DC(N,ND))
!
	DELTA_T=0.5D0*(R2(ND2)-R1(ND1))/V1(ND1)
	V(1:ND)=V1(1:ND)
	R_NEW(1:ND)=R1(1:ND)+V(1:ND)*DELTA_T
	XV1(1:ND1)=LOG(V1)
	XV2(1:ND2)=LOG(V2)
	XNEW(1:ND)=XV1(1:ND1)
!
	ED2=LOG(ED2)
	CALL MON_INTERP(TA,ND,IONE,XNEW,ND,ED2,ND2,XV2,ND2)
	ED(1:ND)=EXP( (1.0D0-X1)*LOG(ED1(1:ND))+X1*TA)
	WRITE(6,*)ED1(1),EXP(ED2(1)),ED(1)
!
	DION2=LOG(DION2)
	CALL MON_INTERP(TC,ND,IONE,XNEW,ND,DION2,ND2,XV2,ND2)
	DION(1:ND)=EXP( (1.0D0-X1)*LOG(DION1(1:ND))+X1*TC)
!
	T2=LOG(T2)
	CALL MON_INTERP(TC,ND,IONE,XNEW,ND,T2,ND2,XV2,ND2)
	T(1:ND)=EXP( (1.0D0-X1)*LOG(T1(1:ND))+X1*TC)
!
	DO J=1,N2
	  DO I=1,ND
	    TA(I)=LOG(DC2(J,I))
	  END DO
	  CALL MON_INTERP(TC,ND,IONE,XNEW,ND,TA,ND2,XV2,ND2)
	  DC(J,1:ND)=EXP( (1.0D0-X1)*LOG(DC1(J,1:ND))+X1*TC(1:ND))
	END DO
!
	DO I=1,ND
	  DO J=1,N2
	    DC(J,I)=EXP(HDKT*EDGE(J)*(1.0D0/DC(J,I)-1.0D0/T(I)))*
	1         (T(I)/DC(J,I))**1.5D0
	  END DO
	END DO
!
	WRITE(6,'(A,A)')'Done ',TRIM(DC_FILE)
!
	I=INDEX(DC_FILE,'OUT')
	IF(I .NE. 0)DC_FILE(I:I+2)='_IN'
	OPEN(UNIT=9,FILE=DC_FILE,STATUS='NEW')
	  WRITE(9,'(/,1X,A,T40,A)')'24-FEB-2004','!Format date'
	  WRITE(9,'(/,1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')R_NEW(ND),LSTAR,N,ND
	  DO I=1,ND
	    WRITE(9,'(/,1X,ES15.7,6ES16.7,2X,I4,A1)')R_NEW(I),DION(I),ED(I),T(I),0.0D0,V(I),1.0D0,I,' '
	    WRITE(9,'(5ES16.7)')(DC(J,I),J=1,N)
	  END DO
	CLOSE(UNIT=9)
!
	END DO
1000	CONTINUE
!
	STOP
	END
