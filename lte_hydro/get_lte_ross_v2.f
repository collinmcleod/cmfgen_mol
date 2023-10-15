	MODULE MOD_LTE_ROSS_TAB
	IMPLICIT NONE
!
	REAL(10), ALLOCATABLE :: CHI(:,:)
	REAL(10), ALLOCATABLE :: ESEC(:,:)
	REAL(10), ALLOCATABLE :: KAP(:,:)
	REAL(10), ALLOCATABLE :: KES(:,:)
	REAL(10), ALLOCATABLE :: POP_ATOM(:,:)
	REAL(10), ALLOCATABLE :: ED(:,:)
	REAL(10), ALLOCATABLE :: RHO(:)
	REAL(10), ALLOCATABLE :: TEMP(:)
!
	INTEGER N_T,I_T
	INTEGER N_D,I_D
	INTEGER LUER
!
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
	END MODULE MOD_LTE_ROSS_TAB
!
	SUBROUTINE GET_LTE_ROSS_V2(KAP_VAL,KAP_ES,LTE_ED,ATOM_DEN,TVAL)
	USE MOD_LTE_ROSS_TAB
	IMPLICIT NONE
!
	REAL(10) ATOM_DEN
	REAL(10) KAP_VAL
	REAL(10) KAP_ES
	REAL(10) LTE_ED
	REAL(10) TVAL
	REAL(10) DENSITY
!
	REAL(10) T1,T2
	REAL(10) D1,D2
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
	DENSITY=ATOM_DEN*RHO(1)/POP_ATOM(1,1)
	IF(DENSITY .GT. RHO(N_D))THEN
	   WRITE(LUER,*)'Error in GET_LTE_ROSS_V2'
	   WRITE(LUER,*)'Density outside range'
	   WRITE(LUER,*)'Maximum density=',RHO(N_D)
	   WRITE(LUER,*)'Requested =',DENSITY
	   STOP
	END IF
	I_D=1
	DO I_D=2,N_D
	  IF(DENSITY .LT. RHO(I_D))EXIT
	END DO
!
	IF(TVAL .GT. TEMP(N_T))THEN
	   WRITE(LUER,*)'Error in GET_LTE_ROSS_V2'
	   WRITE(LUER,*)'Temperature outside range'
	   WRITE(LUER,*)'Maximum temperature=',TEMP(N_T)
	   WRITE(LUER,*)'Requested temperature=',TVAL
	   STOP
	END IF
	I_T=1
	DO I_T=2,N_T
	  IF(TVAL .LT. TEMP(I_T))EXIT
	END DO
	IF(I_T .GT. N_T)I_T=N_T-1
!
	T1=LOG(TEMP(I_T)/TVAL)/LOG(TEMP(I_T)/TEMP(I_T-1))
	T2=LOG(RHO(I_D)/DENSITY)/LOG(RHO(I_D)/RHO(I_D-1))
!
	D1=T1*LOG(KAP(I_T-1,I_D-1))+(1.0D0-T1)*LOG(KAP(I_T,I_D-1))
	D2=T1*LOG(KAP(I_T-1,I_D))+(1.0D0-T1)*LOG(KAP(I_T,I_D))
	KAP_VAL=EXP(T2*D1+(1.0D0-T2)*D2)
!
	D1=T1*LOG(KES(I_T-1,I_D-1))+(1.0D0-T1)*LOG(KES(I_T,I_D-1))
	D2=T1*LOG(KES(I_T-1,I_D))+(1.0D0-T1)*LOG(KES(I_T,I_D))
	KAP_ES=EXP(T2*D1+(1.0D0-T2)*D2)
!
	D1=T1*LOG(ED(I_T-1,I_D-1))+(1.0D0-T1)*LOG(ED(I_T,I_D-1))
	D2=T1*LOG(ED(I_T-1,I_D))+(1.0D0-T1)*LOG(ED(I_T,I_D))
	LTE_ED=EXP(T2*D1+(1.0D0-T2)*D2)
!
	RETURN
	END
