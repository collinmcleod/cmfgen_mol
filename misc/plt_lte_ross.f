	MODULE MOD_LTE_ROSS_TAB_PLT
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: CHI(:,:)
	REAL*8, ALLOCATABLE :: ESEC(:,:)
	REAL*8, ALLOCATABLE :: KAP(:,:)
	REAL*8, ALLOCATABLE :: KES(:,:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:,:)
	REAL*8, ALLOCATABLE :: ED(:,:)
	REAL*8, ALLOCATABLE :: RHO(:)
	REAL*8, ALLOCATABLE :: TEMP(:)
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
!
	INTEGER N_T,I_T
	INTEGER N_D,I_D,ID_MIN
	INTEGER LUER
!
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
	END MODULE MOD_LTE_ROSS_TAB_PLT
!
	PROGRAM PLT_LTE_ROSS
	USE MOD_LTE_ROSS_TAB_PLT
	IMPLICIT NONE
!
! Created 30-May-2019 - Still underdevelopment.
!
	REAL*8 ATOM_DEN
	REAL*8 KAP_VAL
	REAL*8 KAP_ES
	REAL*8 LTE_ED
	REAL*8 TVAL
	REAL*8 DENSITY
!
	REAL*8 T1,T2
	REAL*8 D1,D2
	INTEGER LU
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER*132 STRING
!
	LU=10
	IF(FIRST_TIME)THEN
	  LUER=ERROR_LU()
	  OPEN(UNIT=LU,FILE='ROSSELAND_LTE_TAB',STATUS='OLD',ACTION='READ')
	    STRING=' '
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(LU,'(A)')STRING
	    END DO
	    IF(INDEX(STRING,'!Number of temperatures') .NE. 0)THEN
	      READ(STRING,*)N_T
	    ELSE IF(INDEX(STRING,'!Number of densities') .NE. 0)THEN
	      READ(STRING,*)N_D
	    ELSE
	      WRITE(LUER,*)'Error reading ROSSELAND_LTE_TAB'
	      WRITE(LUER,*)'Unrecognized record'
	      WRITE(LUER,'(A)')TRIM(STRING)
	      STOP
	    END IF
	    STRING=' '
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(LU,'(A)')STRING
	    END DO
	    IF(INDEX(STRING,'Number of temperatures') .NE. 0)THEN
	      READ(STRING,*)N_T
	    ELSE IF(INDEX(STRING,'Number of densities') .NE. 0)THEN
	      READ(STRING,*)N_D
	    ELSE
	      WRITE(LUER,*)'Error reading ROSSELAND_LTE_TAB'
	      WRITE(LUER,*)'Unrecognized record'
	      WRITE(LUER,'(A)')TRIM(STRING)
	      STOP
	    END IF
	    STRING=' '
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(LU,'(A)')STRING
	    END DO
	    BACKSPACE(LU)
!
	    IF(N_D .EQ. 0)THEN
	      WRITE(LUER,*)'Error reading N_D from ROSSELAND_LTE_TAB'
	      STOP
	    ELSE IF(N_T .EQ. 0)THEN
	      WRITE(LUER,*)'Error reading N_T from ROSSELAND_LTE_TAB'
	      STOP
	    END IF
!
	    ALLOCATE (CHI(N_T,N_D))
	    ALLOCATE (ESEC(N_T,N_D))
	    ALLOCATE (KAP(N_T,N_D))
	    ALLOCATE (KES(N_T,N_D))
	    ALLOCATE (ED(N_T,N_D))
	    ALLOCATE (POP_ATOM(N_T,N_D))
	    ALLOCATE (RHO(N_D))
	    ALLOCATE (TEMP(N_T))
	    ALLOCATE (XV(MAX(N_T,N_D)))
	    ALLOCATE (YV(MAX(N_T,N_D)))
!
	    DO I_D=1,N_D
	      DO I_T=1,N_T
	        READ(LU,*)TEMP(I_T),RHO(I_D),POP_ATOM(I_T,I_D),ED(I_T,I_D),
	1                 CHI(I_T,I_D),ESEC(I_T,I_D),
	1                 KAP(I_T,I_D),KES(I_T,I_D)
	      END DO
	    END DO
!
	  CLOSE(UNIT=LU)
!
	  FIRST_TIME=.FALSE.
	END IF
!
	WRITE(6,*)' Atom density'
	DO I_D=1,N_D
	  IF(POP_ATOM(1,I_D) .GT. 1.0D+10)EXIT
	END DO
	ID_MIN=I_D 
!
	WRITE(6,'(/,A)')' Atom density'
	WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	WRITE(6,'(/,A,/)')' Plotting Kappa(Ross)/Kappa(es) versus the electron temperature '
	DO I_D=ID_MIN,N_D,2
	   XV(1:N_T)=TEMP(1:N_T)
	   YV(1:N_T)=KAP(1:N_T,I_D)/KES(1:N_T,I_D)
	   CALL DP_CURVE(N_T,XV,YV)
	END DO
	CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K','\gK/\K\des\u',' ',' ')
!
	WRITE(6,'(/,A)')' Atom density'
	WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	WRITE(6,'(/,A,/)')' Plotting Kappa(Ross) versus the electron temperature '
	DO I_D=ID_MIN,N_D,2
	   XV(1:N_T)=TEMP(1:N_T)
	   YV(1:N_T)=KAP(1:N_T,I_D)
	   CALL DP_CURVE(N_T,XV,YV)
	END DO
	CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K','\gK',' ',' ')
!
	WRITE(6,'(/,A)')' Atom density'
	WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	WRITE(6,'(/,A,/)')' Plotting T^{3.5}.Kappa(Ross)/rho versus the electron temperature '
	DO I_D=ID_MIN,N_D,2
	   XV(1:N_T)=TEMP(1:N_T)
	   YV(1:N_T)=(TEMP(1:N_T)**3.5)*(KAP(1:N_T,I_D)/KES(1:N_T,I_D))/RHO(I_D)
	   CALL DP_CURVE(N_T,XV,YV)
	END DO
	CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K','T\u3.5\d(\gK/\K\des\u)/\gr',' ',' ')
!
	WRITE(6,'(/,A)')' Atom density'
	WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	WRITE(6,'(/,A,/)')' Plotting T^{3.5}.(Kappa(Ross)/K(es)-1)/rho versus the electron temperature '
	DO I_D=ID_MIN,N_D,2
	   XV(1:N_T)=TEMP(1:N_T)
	   YV(1:N_T)=(TEMP(1:N_T)**3.5)*(KAP(1:N_T,I_D)/KES(1:N_T,I_D)-1.0D0)/RHO(I_D)
	   CALL DP_CURVE(N_T,XV,YV)
	END DO
	CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K','T\u3.5\d(\gK/\K\des\u-1)/\gr',' ',' ')
!
	STOP
	END
