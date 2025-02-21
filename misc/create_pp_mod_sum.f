!
! Simple program to read in a FILE (MOD_DIRS) containg a list of MODEL directories.
! Program reads the MOD_SUM file in each of these directories, and outputs a summary
! file (MODEL_SUMMARY).
!
	PROGRAM CREATE_PP_MOD_SUM
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL(KIND=LDP) T1,T2,T3
	INTEGER K,IOS
	INTEGER LDIR
	INTEGER I,J
!
	LOGICAL USE_MF
	LOGICAL NO_HYD_RD
	LOGICAL FILE_PRES
!	
	CHARACTER(LEN=200) STRING
	CHARACTER(LEN=100) FRM
	CHARACTER(LEN=100) FILENAME,MOD_FILENAME
	CHARACTER(LEN=10)  LSTAR,MDOT,RSTAR
	CHARACTER(LEN=9)  VINF,VINF2
	CHARACTER(LEN=18) BOTH_VINF
	CHARACTER(LEN=8)  TEFF
	CHARACTER(LEN=6)  LOGG
	CHARACTER(LEN=8)  NH,NHe,XCARB,XNIT,XOXY,XNEON,XIRON
	CHARACTER(LEN=6)  ND
	CHARACTER(LEN=6)  NT
	CHARACTER(LEN=8)  BETA,BETA2
	CHARACTER(LEN=10) BOTH_BETA
	CHARACTER(LEN=6)  CL_F
	CHARACTER(LEN=7)  CL_VF
	CHARACTER(LEN=11) DATE
!
! Initialize header variables. Done this way as change in string length does
! not require a change here.
!
	LSTAR='L/Lsun'; MDOT='Mdot'; BOTH_VINF='Vinf'
	RSTAR='R/Rsun'; LOGG='Logg'; NT='NT'; ND='ND'
	CL_F='f'; CL_VF='V(f)';  BOTH_BETA='Beta'
	XCARB='X(C/s)'; XNIT='X(N/s)'; XOXY='X(O/s)'; XNEON='X(Ne/s)'; XIRON='X(Fe/s)'
	NH='N(H)'; NHe='N(He)'; TEFF='Teff'; DATE='Date'
	FILENAME=' '
!
	NO_HYD_RD=.FALSE.
	USE_MF=.FALSE.
	CALL GEN_IN(NO_HYD_RD,'Are the models devoid of hydrogen')
	CALL GEN_IN(USE_MF,'Output mass fractions')
!
	IF(USE_MF)THEN
	  NH=' X(H)'; NHe=' X(He)'
	  XCARB='  X(C)'; XNIT='  X(N)'; XOXY='  X(O)'; XNEON='  X(Ne)'
	END IF	
!
	TEFF=ADJUSTR(TEFF)
	LOGG=ADJUSTR(LOGG)
	LSTAR=ADJUSTR(LSTAR)
	RSTAR=ADJUSTR(RSTAR)
	NH=ADJUSTR(NH)
	NHe=ADJUSTR(NHe)
	XCARB=ADJUSTR(XCARB)
	XNIT=ADJUSTR(XNIT)
	XOXY=ADJUSTR(XOXY)
	XNEON=ADJUSTR(XNEON)
	XIRON=ADJUSTR(XIRON)
	NT=ADJUSTR(NT)
	ND=ADJUSTR(ND)
	DATE=ADJUSTR(DATE)
!
	OPEN(UNIT=12,FILE='MODEL_SUMMARY',STATUS='UNKNOWN',ACTION='WRITE')
!
! Outut header.
! Note: This write should be identical to that at the end of the program.
!
	OPEN(UNIT=8,FILE='MOD_DIRS',STATUS='OLD',ACTION='READ')
	LDIR=1
	DO WHILE(1 .EQ. 1)
	  READ(8,'(A)',END=500)FILENAME
	  LDIR=MAX(LDIR,LEN_TRIM(FILENAME))
	END DO
500	CONTINUE
	LDIR=LDIR+1
!
	WRITE(FRM,'(I3)')LDIR
	FRM=ADJUSTL(FRM)
	FILENAME=' '
	IF(NO_HYD_RD)THEN
	  FRM='(A,T'//TRIM(FRM)//',11A,2X,A11)'
	  WRITE(6,*)LDIR,FRM
	  WRITE(12,FRM)TRIM(ADJUSTL(FILENAME)),TEFF,LOGG,LSTAR,RSTAR,
	1                          NHe,XCARB,XOXY,XNEON,XIRON,ND,NT,DATE
	ELSE
	  FRM='(A,T'//TRIM(FRM)//',12A,2X,A11)'
	  WRITE(12,FRM)TRIM(ADJUSTL(FILENAME)),TEFF,LOGG,LSTAR,RSTAR,
	1                          NH,NHe,XCARB,XNIT,XOXY,XIRON,ND,NT,DATE
	END IF
!
	REWIND(UNIT=8)	
	DO WHILE(1 .EQ. 1)
	  READ(8,'(A)',END=9999)FILENAME
!
	  MOD_FILENAME=TRIM(FILENAME)//'/MOD_SUM'
	  INQUIRE(FILE=MOD_FILENAME,EXIST=FILE_PRES)
	  IF(.NOT. FILE_PRES)THEN
	    WRITE(12,'(A)')TRIM(FILENAME)
	    GOTO 5000
	  END IF
!
	  OPEN(UNIT=10,FILE=MOD_FILENAME,STATUS='OLD',ACTION='READ')
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'Finalized') .EQ. 0 )
	    READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,'on:')+3
	  STRING=ADJUSTL(STRING(K:))
	  DATE=STRING(1:11)
!
	  DO WHILE( INDEX(STRING,'ND[') .EQ. 0 )
	    READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,'ND[')
	  STRING=STRING(K+3:)
	  K=INDEX(STRING,']')
	  ND=STRING(1:K-1)
!
	  K=INDEX(STRING,'NT[')
	  STRING=STRING(K+3:)
	  K=INDEX(STRING,']')
	  NT=STRING(1:K-1)
	  WRITE(6,*)'NT=',TRIM(NT)
!
	  DO WHILE( INDEX(STRING,'L*') .EQ. 0 )
	    READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,' ')
	  READ(STRING(4:K),*)T1
	  WRITE(LSTAR,'(ES8.2)')T1
	  WRITE(6,*)'Done LSTAR=',TRIM(LSTAR)
!
	  DO WHILE( INDEX(STRING,'Tau') .EQ. 0 )
	    READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,'Rsun')
	  STRING=STRING(K+5:)
	  K=INDEX(STRING,' ')
	  RSTAR=STRING(1:K)
	  READ(RSTAR,*)T1
	  IF(T1 .LT. 1.0_LDP .AND. T1 .GT. 0.1_LDP)THEN
	    WRITE(RSTAR,'(F5.3)')T1
	  END IF
!
	  K=INDEX(STRING,'Teff')
	  IF(INDEX(STRING,'Teff(K)') .NE. 0)THEN
	    STRING=STRING(K+8:)
	  ELSE
	    STRING=STRING(K+4:)
	  END IF
	  READ(STRING,*)T1
	  WRITE(TEFF,'(F6.3)')T1/1.0D+04
!
	  K=INDEX(STRING,'Log g')
	  IF(K .EQ. 0)THEN
	    LOGG=' '
	  ELSE
	    STRING=STRING(K+6:)
	    LOGG=STRING(1:4)
	  END IF
!
	  DO WHILE( INDEX(STRING,'HYD') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,'HYD')+4
	  READ(STRING(K:),*)T2,T3
	  IF(USE_MF)THEN
	    WRITE(NH,'(F6.3)')T3
	  ELSE
	    WRITE(NH,'(F6.3)')T2
	  ENDIF
!
	  READ(10,'(A)')STRING
	  K=INDEX(STRING,'HE')+4
	  IF(USE_MF)THEN
	    READ(STRING(K:),*)T2,T3
	    WRITE(NHe,'(F6.3)')T3
	  ELSE
	    READ(STRING(K:),*)T2
	    WRITE(NHe,'(F5.2)')T2
	  END IF
!
	  READ(10,'(A)')STRING
	  K=INDEX(STRING,'CARB')+4
	  READ(STRING(K:),*)T3,T2,T1
	  IF(USE_MF)THEN
	    WRITE(XCARB,'(F7.3)')T2
	  ELSE
	    WRITE(XCARB,'(F7.3)')T1
	  END IF
!	
	  READ(10,'(A)')STRING
	  K=INDEX(STRING,'NIT')+4
	  READ(STRING(K:),*)T3,T2,T1
	  IF(USE_MF)THEN
	    WRITE(XNIT,'(F6.3)')T2
	  ELSE
	    WRITE(XNIT,'(F6.3)')T1
	  END IF
!	
	  READ(10,'(A)')STRING
	  K=INDEX(STRING,'OXY')+4
	  READ(STRING(K:),*)T3,T2,T1
	  IF(USE_MF)THEN
	    WRITE(XOXY,'(F6.3)')T2
	    IF(T3 .LT. 0.05)WRITE(XOXY,'(F6.4)')T2
	  ELSE
	    WRITE(XOXY,'(F6.3)')T1
	  END IF
!
	  READ(10,'(A)')STRING        !Fluorine
	  IF(INDEX(STRING,'FLU') .NE. 0)READ(10,'(A)')STRING
	  K=INDEX(STRING,'NEON')+5
	  READ(STRING(K:),*)T3,T2,T1
	  IF(USE_MF)THEN
	    WRITE(XNEON,'(F6.3)')T2
	    IF(T2 .LT. 0.05)WRITE(XNEON,'(F6.4)')T2
	  ELSE
	    WRITE(XNEON,'(F6.3)')T1
	  END IF
!
	  DO WHILE(INDEX(STRING,'IRON') .EQ. 0)
	     READ(10,'(A)')STRING
	  END DO
	  K=INDEX(STRING,'IRON')+4
	  READ(STRING(K:),*)T3,T2,T1
	  WRITE(XIRON,'(F6.3)')T1                              !Abundance ratio
!
	  IOS=0
	  DO WHILE(INDEX(STRING,'CL_P') .EQ. 0 .AND. IOS .EQ. 0)
	    READ(10,'(A)',IOSTAT=IOS)STRING
	  END DO
	  IF(IOS .EQ. 0)THEN
	    K=INDEX(STRING,'CL_P_1')+7
	    READ(STRING(K:),*)T1
	    WRITE(CL_F,'(F4.2)')T1
	    K=INDEX(STRING,'CL_P_2')+7
	    READ(STRING(K:),*)T1
	    WRITE(CL_VF,'(F6.1)')T1
	  ELSE
	    WRITE(CL_F,'(F4.2)')1.0D0
	    WRITE(CL_VF,'(F6.1)')0.0D0
	  END IF
!	
	  CLOSE(UNIT=10)
	  WRITE(6,*)'Read ',TRIM(FILENAME)
!
	  MOD_FILENAME=TRIM(FILENAME)//'/HYDRO_DEFAULTS'
	  INQUIRE(FILE=MOD_FILENAME,EXIST=FILE_PRES)
	  IF(FILE_PRES .AND. BETA2 .NE. ' ')THEN
	  OPEN(UNIT=10,FILE=MOD_FILENAME,STATUS='OLD')
	    DO WHILE(1 .EQ. 1)
	      READ(10,'(A)',END=2000)STRING
	      IF(INDEX(STRING,'BETA') .NE. 0)THEN
	         STRING=ADJUSTL(STRING)
	         I=INDEX(STRING,' ')
	         BETA=ADJUSTL(BETA)
	         BOTH_BETA=TRIM(BETA)//'/'//STRING(1:I)
	         CLOSE(UNIT=10)
	      ELSE IF(INDEX(STRING,'ITS_DONE') .NE. 0)THEN
	        READ(STRING,*)I
	        IF(I .NE. 0)BETA=BOTH_BETA
	      END IF
            END DO
	  END IF
2000      CONTINUE
	
	  TEFF=ADJUSTR(TEFF)
	  LOGG=ADJUSTR(LOGG)
	  LSTAR=ADJUSTR(LSTAR)
	  RSTAR=ADJUSTR(RSTAR)
	  MDOT=ADJUSTR(MDOT)
	  CL_F=ADJUSTR(CL_F)
	  CL_VF=ADJUSTR(CL_VF)
	  VINF=ADJUSTR(VINF)
	  BETA=ADJUSTR(BETA)
	  VINF2=ADJUSTR(VINF2)
	  BETA2=ADJUSTR(BETA2)
	  NH=ADJUSTR(NH)
	  NHe=ADJUSTR(NHe)
	  XCARB=ADJUSTR(XCARB)
	  XNIT=ADJUSTR(XNIT)
	  XOXY=ADJUSTR(XOXY)
	  XNEON=ADJUSTR(XNEON)
	  XIRON=ADJUSTR(XIRON)
	  NT=ADJUSTR(NT)
	  ND=ADJUSTR(ND)
!
	  IF(BETA2 .NE. ' ')THEN
	    BOTH_BETA=TRIM(ADJUSTL(BETA))//'/'//ADJUSTL(BETA2)
	  ELSE
	    BOTH_BETA=ADJUSTL(BETA)
	  END IF
	  BOTH_BETA=ADJUSTR(BOTH_BETA)
	  IF(VINF2 .NE. ' ')THEN
	    BOTH_VINF=TRIM(ADJUSTL(VINF))//'/'//ADJUSTL(VINF2)
	  ELSE
	    BOTH_VINF=ADJUSTL(VINF)
	  END IF
	  BOTH_VINF=ADJUSTR(BOTH_VINF)
!
	  WRITE(FRM,'(I5)')LDIR
	  FRM=ADJUSTL(FRM)
	  IF(NO_HYD_RD)THEN
	    FRM='(A,T'//TRIM(FRM)//',11A,2X,A11)'
	    WRITE(12,FRM)TRIM(ADJUSTL(FILENAME)),TEFF,LOGG,LSTAR,RSTAR,
	1                          NHe,XCARB,XOXY,XNEON,XIRON,ND,NT,DATE
	  ELSE
	    FRM='(A,T'//TRIM(FRM)//',12A,2X,A11)'
	    WRITE(12,FRM)TRIM(ADJUSTL(FILENAME)),TEFF,LOGG,LSTAR,RSTAR,
	1                          NH,NHe,XCARB,XNIT,XOXY,XIRON,ND,NT,DATE
	  END IF
	  FLUSH(UNIT=12)
5000	CONTINUE
!
	END DO
9999	CONTINUE
!
	STOP
	END
