	MODULE MOD_LTE_ROSS_TAB_PLT
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
	REAL(10), ALLOCATABLE :: XV(:)
	REAL(10), ALLOCATABLE :: YV(:)
!
	INTEGER N_T,I_T,IT_MIN,D_IT
	INTEGER N_D,I_D,ID_MIN,D_ID
	INTEGER LUER
!
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
	END MODULE MOD_LTE_ROSS_TAB_PLT
!
	PROGRAM PLT_LTE_ROSS
	USE GEN_IN_INTERFACE
	USE MOD_LTE_ROSS_TAB_PLT
	IMPLICIT NONE
!
! Altered 20-Sep-2021 - Changed to allow greater option flexibility.
!                          New options added.
! Created 30-May-2019 - Still underdevelopment.
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
	REAL(10) ALPHA
	INTEGER I
	INTEGER LU
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=10) XAXIS
	CHARACTER(LEN=50) YLAB
	CHARACTER(LEN=30) UC
	CHARACTER(LEN=10) XOPT
	EXTERNAL UC
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
	ID_MIN=I_D; D_ID=2
	XAXIS='TEMP'
	IT_MIN=1; D_IT=2
!
	XOPT='KR'
3	CALL GEN_IN(XOPT,'Input plotting option')
	XOPT=UC(XOPT)
!
	IF(XOPT .EQ. 'XATOM')THEN
	   XAXIS='ATOM'
	   WRITE(6,'(4(F12.3,A,I2,A))')(TEMP(I_T),'(',I_T,')',I_T=1,N_T)
	   CALL GEN_IN(IT_MIN,'Mimimum index to be used for plotting')
	   CALL GEN_IN(D_IT,'Step size in temperature index')
!
	ELSE IF(XOPT .EQ. 'XTEMP')THEN
	   XAXIS='TEMP'
	   WRITE(6,'(4(F12.3,A,I2,A))')(LOG10(POP_ATOM(1,I_D)),'(',I_D,')',I_D=1,N_D)
	   CALL GEN_IN(ID_MIN,'Mimimum index to be used for plotting')
	   CALL GEN_IN(D_ID,'Step size in atom index')
!
	ELSE IF(XOPT .EQ. 'KR/KES' .AND. XAXIS .EQ. 'ATOM')THEN
	  WRITE(6,'(/,A,/)')' Plotting Kappa(Ross)/Kappa(es) versus the atom density'
	  WRITE(6,'(5ES12.4)')(TEMP(I_T),I_T=IT_MIN,N_T,2)
	  DO I_T=IT_MIN,N_T,D_IT
	     XV(1:N_D)=LOG10(POP_ATOM(I_T,1:N_D))
	     YV(1:N_D)=KAP(I_T,1:N_D)/KES(I_T,1:N_D)
	     CALL DP_CURVE(N_D,XV,YV)
	  END DO
	  CALL GRAMON_PGPLOT('Log N\dA\u','\gK/\K\des\u',' ',' ')
!
	ELSE IF(XOPT .EQ. 'KR/KES' .AND. XAXIS .EQ. 'TEMP')THEN
	  WRITE(6,'(/,A)')' Atom density'
	  WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	  WRITE(6,'(/,A,/)')' Plotting Kappa(Ross)/Kappa(es) versus the electron temperature '
	  DO I_D=ID_MIN,N_D,D_ID
	     XV(1:N_T)=TEMP(1:N_T)
	     YV(1:N_T)=KAP(1:N_T,I_D)/KES(1:N_T,I_D)
	     CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  CALL GRAMON_PGPLOT('T(10\u4\d)\,K','\gK/\K\des\u',' ',' ')
!
	ELSE IF(XOPT .EQ. 'KR' .AND.  XAXIS .EQ. 'ATOM')THEN
	  WRITE(6,'(/,A,/)')' Plotting Kappa(Ross) versus the atom density'
	  WRITE(6,'(/,A)')' Temperature/10^4 K'
	  WRITE(6,'(5ES12.4)')(TEMP(I_T),I_T=IT_MIN,N_T,D_IT)
	  DO I_T=IT_MIN,N_T,D_IT
	     XV(1:N_D)=LOG10(POP_ATOM(I_T,1:N_D))
	     YV(1:N_D)=KAP(I_T,1:N_D)
	     CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  CALL GRAMON_PGPLOT('Log N\dA\u','\gK',' ',' ')
!
	ELSE IF(XOPT .EQ. 'KR' .AND.  XAXIS .EQ. 'TEMP')THEN
	  WRITE(6,'(/,A)')' Atom density'
	  WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,2)
	  WRITE(6,'(/,A,/)')' Plotting Kappa(Ross) versus the electron temperature '
	  DO I_D=ID_MIN,N_D,D_ID
	     XV(1:N_T)=TEMP(1:N_T)
	     YV(1:N_T)=KAP(1:N_T,I_D)
	     CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  CALL GRAMON_PGPLOT('T(10\u4\d)\,K','\gK',' ',' ')
!
	ELSE IF(XOPT .EQ. 'RK')THEN
	  WRITE(6,'(/,A)')' Atom density'
	  WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,D_ID)
	  WRITE(6,'(/,A,/)')' Plotting T^{3.5}.Kappa(Ross)/rho**ALPHA versus the electron temperature '
          WRITE(6,'(A)',ADVANCE='NO')'Value of alpha -- 0 to exit: '
          READ(5,*)ALPHA
	  DO I_D=ID_MIN,N_D,D_ID
	   XV(1:N_T)=TEMP(1:N_T)
	   YV(1:N_T)=(TEMP(1:N_T)**3.5)*(KAP(1:N_T,I_D)/KES(1:N_T,I_D))/RHO(I_D)**ALPHA
	   CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  WRITE(YLAB,'(F5.2)')ALPHA; YLAB=ADJUSTL(YLAB)
	  I=LEN_TRIM(YLAB); IF(YLAB(I:I) .EQ. '0')I=I-1
	  YLAB='T\u3.5\d(\gK/\K\des\u)/\gr\u'//YLAB(1:I)//'\d'
	  CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K',YLAB,' ',' ')
!
	ELSE IF(XOPT .EQ. 'RK/KES')THEN
	  WRITE(6,'(/,A)')' Atom density'
	  WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,D_ID)
          WRITE(6,'(A)',ADVANCE='NO')'Value of alpha -- 0 to exit: '
          READ(5,*)ALPHA
	  WRITE(6,'(/,A,/)')' Plotting T^{3.5}.(Kappa(Ross)/K(es)-1)/rho**alpha versus the electron temperature '
	  DO I_D=ID_MIN,N_D,D_ID
	     XV(1:N_T)=TEMP(1:N_T)
	     YV(1:N_T)=(TEMP(1:N_T)**3.5)*(KAP(1:N_T,I_D)/KES(1:N_T,I_D)-1.0D0)/RHO(I_D)**ALPHA
	     CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  WRITE(YLAB,'(F5.2)')ALPHA; YLAB=ADJUSTL(YLAB)
	  I=LEN_TRIM(YLAB); IF(YLAB(I:I) .EQ. '0')I=I-1
	  YLAB='T\u3.5\d(\gK/\K\des\u-1)/\gr\u'//YLAB(1:I)//'\d'
	  CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K',YLAB,' ',' ')
!
	ELSE IF(XOPT .EQ. 'ED/NA')THEN
	  WRITE(6,'(/,A)')' Atom density'
	  WRITE(6,'(5ES12.4)')(POP_ATOM(1,I_D),I_D=ID_MIN,N_D,D_ID)
	  DO I_D=ID_MIN,N_D,D_ID
	     XV(1:N_T)=TEMP(1:N_T)
	     YV(1:N_T)=ED(1:N_T,I_D)/POP_ATOM(1:N_T,I_D)
	     CALL DP_CURVE(N_T,XV,YV)
	  END DO
	  CALL GRAMON_PGPLOT('Log T(10\u4\d)\,K','N\de\u/N\dA\u',' ',' ')
!
	ELSE IF(XOPT .EQ. 'HE' .OR. XOPT .EQ. 'HELP') THEN
	  WRITE(6,*)'XATOM:       Set xaxis to atom density'
	  WRITE(6,*)'XTEMP:       Set xaxis to temperature'
	  WRITE(6,*)'KR/KES:      Plot Kappa(Ross)/Kappa(es)'
	  WRITE(6,*)'KR:          Plot Kappa(Ross)'     
	  WRITE(6,*)'RK:          Plot T^{3.5}.Kappa(Ross)/rho**ALPHA '
	  WRITE(6,*)'RK/KES:      Plot T^{3.5}.(Kappa(Ross)/K(es)-1)/rho**alpha '
	  WRITE(6,*)'NE/NA:       Plot Electron density / Atom density'
!
	ELSE IF(XOPT .EQ. 'EX' .OR. XOPT .EQ. 'EXIT') THEN
	  STOP
	ELSE
	  WRITE(6,*)'Option requested does not exist'
	END IF
!
1       GO TO 3
	END
