!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc. 
!
	PROGRAM MOD_PRRR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 17-Nov-2021: Transfered from OSIRIS.
! Altered: 22-Oct-2021: Major editing to improve optios (and earlier)
! Altered: 27-May-2021: Small bug fixis (transferred from OSIRIS 21-Jul-20201).
! Altered: 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered: 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered: 08-Feb-2008: Extra terms (such as V term) sheck and output.
!
!
	REAL*8, ALLOCATABLE :: PHOT(:,:)
	REAL*8, ALLOCATABLE :: RECOM(:,:)
	REAL*8, ALLOCATABLE :: RECOM_SUM(:)
	REAL*8, ALLOCATABLE :: PHOT_SUM(:)
	REAL*8, ALLOCATABLE :: COL_IR(:)
	REAL*8, ALLOCATABLE :: CHG_IR(:)
	REAL*8, ALLOCATABLE :: NT_IR(:)
	REAL*8, ALLOCATABLE :: COL_RR(:)
	REAL*8, ALLOCATABLE :: XRAY_RR(:)
	REAL*8, ALLOCATABLE :: CHG_RR(:)
	REAL*8, ALLOCATABLE :: ADVEC_RR(:)
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: TEMP(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: DI(:)
!
	REAL*8, ALLOCATABLE :: XVEC(:)
	REAL*8, ALLOCATABLE :: YVEC(:)
	REAL*8, ALLOCATABLE :: TOT_SUM(:)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,L
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	INTEGER IST,IEND,NLEV
	REAL*8 T1
	REAL*8 MIN_VAL
	LOGICAL FILE_OPEN
	LOGICAL NET_RECOM_PER_LEVEL
	LOGICAL DO_INDIV_RATES
	LOGICAL NORM
	LOGICAL ABS_VALUE
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
        CHARACTER(LEN=2), PARAMETER :: FORMFEED=' '//CHAR(12)
!
	CHARACTER(LEN=132) TMP_STR
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) FILE_NAME
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=5) OPTION
	CHARACTER(LEN=30) XLABEL
	CHARACTER(LEN=30) YLABEL
        CHARACTER(LEN=30) UC
        EXTERNAL UC
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	ALLOCATE(PHOT_SUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(RECOM_SUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(COL_IR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(CHG_IR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(NT_IR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(COL_RR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(XRAY_RR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(CHG_RR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(ADVEC_RR(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE(R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(TEMP(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(DI(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE(XVEC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(YVEC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(TOT_SUM(ND),STAT=IOS)
!
	IOS=0
	OPEN(UNIT=20,FILE='RVTJ',STATUS='OLD',IOSTAT=IOS)
	  DO WHILE(INDEX(STRING,'Velocity (km/s)') .EQ. 0)
            READ(20,'(A)',IOSTAT=IOS)STRING
	  END DO
	  READ(20,*,IOSTAT=IOS)(V(I),I=1,ND)
	CLOSE(UNIT=20)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to get V(km/s) from RVTJ'
	  V(1:ND)=1.0
	END IF
!
	NET_RECOM_PER_LEVEL=.FALSE.
	DO_INDIV_RATES=.FALSE.
	CALL GEN_IN(DO_INDIV_RATES,'Ouput recombination rate and photioization rate for each leve')
	CALL GEN_IN(NET_RECOM_PER_LEVEL,'Ouput net recombination rate to each level?')
	SPECIES='FeI'
!
1000	CONTINUE
!
	CHG_IR(1:ND)=0.0D0
	COL_IR(1:ND)=0.0D0
	NT_IR(1:ND)=0.0D0
	CHG_RR(1:ND)=0.0D0
	COL_RR(1:ND)=0.0D0
	XRAY_RR(1:ND)=0.0D0
	ADVEC_RR(1:ND)=0.0D0
	RECOM_SUM(1:ND)=0.0D0
	PHOT_SUM(1:ND)=0.0D0
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(SPECIES,'File is assumed to be SPECIES//PRRR- EX to exit')
	  IF(UC(SPECIES) .EQ. 'EX')STOP
	  FILE_NAME=TRIM(SPECIES)//'PRRR'
	  OPEN(UNIT=20,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,'(2A,I5,A)')' Error -- unable to open file: '//TRIM(FILE_NAME),' (IOS=',IOS,')'
	    WRITE(6,'(A)')' Check capatilzation of ION name'
	  END IF
	END DO
!
	STRING=' '
	DO WHILE(INDEX(STRING,'Photoionization Rate') .EQ. 0)
	  READ(20,'(A)')STRING
	END DO
	NLEV=-1
	DO WHILE(STRING .NE. ' ')
	   READ(20,'(A)')STRING
	   NLEV=NLEV+1
	END DO
	WRITE(6,'(A,I4)')' Number of levels in the model is:',NLEV
	CLOSE(UNIT=20)
!
	ALLOCATE (PHOT(ND,NLEV))
	ALLOCATE (RECOM(ND,NLEV))
	OPEN(UNIT=20,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
!
	FILE_NAME=TRIM(FILE_NAME)//'_SUM'
	OPEN(UNIT=21,FILE=FILE_NAME,STATUS='UNKNOWN',ACTION='WRITE')
!
	STRING=' '
	DO L=1,ND,10
	   IST=L; IEND=MIN(ND,IST+9)
	   WRITE(6,*)IST,IEND
!
	   DO WHILE(INDEX(STRING,'Photoionization Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Photoionization Rate') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Radius') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(R(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(TEMP(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Electron Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(ED(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Ion Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(DI(I),I=IST,IEND)
	     END IF
	   END DO
!
	   DO J=1,NLEV
	     READ(20,*)(PHOT(I,J),I=IST,IEND)
	     PHOT_SUM(IST:IEND)=PHOT_SUM(IST:IEND)+PHOT(IST:IEND,J)
	   END DO
!
	   WRITE(21,'(A,I4,A,I4,A)')'   Total photoionization rate (d=',IST,' to',IEND,'):'
	   WRITE(21,'(X,10ES12.4)')(PHOT_SUM(I),I=IST,IEND)
!
	   IF(.NOT. DO_INDIV_RATES)THEN
	   ELSE
	     WRITE(21,'(A)')
	     WRITE(21,'(A)')TRIM(STRING)
	     DO J=1,NLEV
	       WRITE(21,'(X,10ES12.4)')(PHOT(I,J),I=IST,IEND)
	     END DO
	   END IF
!
	   DO WHILE(INDEX(STRING,'Recombination Rates') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Recombination Rates') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Colisional Ionization Rate') .NE. 0)THEN
	        READ(20,'(A)')STRING
	        WRITE(21,'(A)')TRIM(STRING)
	        READ(STRING,*)(COL_IR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Ionization Rate') .NE. 0)THEN
	        READ(20,'(A)')STRING
	        WRITE(21,'(A)')TRIM(STRING)
	        READ(STRING,*)(CHG_IR(I),I=IST,IEND)
	     END IF
	   END DO
!
	   DO J=1,NLEV
	     READ(20,*)(RECOM(I,J),I=IST,IEND)
	     RECOM_SUM(IST:IEND)=RECOM_SUM(IST:IEND)+RECOM(IST:IEND,J)
	   END DO
!
	   WRITE(21,'(A,I4,A,I4,A)')'   Total recombination rate (d=',IST,' to',IEND,'):'
	   WRITE(21,'(X,10ES12.4)')(RECOM_SUM(I),I=IST,IEND)
!
	   IF(.NOT. DO_INDIV_RATES)THEN
	     WRITE(21,'(A)')' '
	   ELSE IF(NET_RECOM_PER_LEVEL)THEN
	     WRITE(21,'(A)')' '
	     WRITE(21,'(A,A)')'   Net ',TRIM(STRING(4:))
	     DO J=1,NLEV
	       WRITE(21,'(X,10ES12.4)')(RECOM(I,J)-PHOT(I,J),I=IST,IEND)
	     END DO
	   ELSE
	     WRITE(21,'(A)')' '
	     WRITE(21,'(A)')TRIM(STRING)
	     DO J=1,NLEV
	       WRITE(21,'(X,10ES12.4)')(RECOM(I,J),I=IST,IEND)
	     END DO
	   END IF
!
	   DO WHILE(INDEX(STRING,'Net Recombination Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Net Recombination Rate') .EQ. 0)THEN
	       WRITE(21,'(A)')TRIM(STRING)
	     END IF
	     IF(INDEX(STRING,'Colisional Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(COL_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(CHG_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Effective Advection Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(ADVEC_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Non-Thermal Ionization') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(NT_IR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Net X-ray recombination rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,'(A)')TRIM(STRING)
	       READ(STRING,*)(XRAY_RR(I),I=IST,IEND)
	     END IF
	   END DO
!
	   WRITE(21,'(A)')TRIM(STRING)
	   FLUSH(21)
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)',END=200)STRING
	     IF(STRING(1:1) .EQ. '1' .OR. STRING(1:1) .EQ. FORMFEED)THEN
	       EXIT
	       WRITE(21,'(A)')FORMFEED
	     END IF
	     WRITE(21,'(A)')TRIM(STRING)
	   END DO
	   FLUSH(21)
	END DO
200	CONTINUE
!
	WRITE(6,'(8ES14.4)')PHOT_SUM(1),RECOM_SUM(1),
	1             COL_IR(1),CHG_IR(1),COL_RR(1),CHG_RR(1),
	1             ADVEC_RR(1),NT_IR(1)
	WRITE(6,'(8ES14.4)')PHOT_SUM(ND),RECOM_SUM(ND),
	1             COL_IR(ND),CHG_IR(ND),COL_RR(ND),CHG_RR(ND),
	1             ADVEC_RR(ND),NT_IR(ND)
	WRITE(6,'(8ES14.4)')R(1),V(1),TEMP(1),R(ND),V(ND),TEMP(ND)
!
	DO I=1,ND
	   TOT_SUM(I)=( PHOT_SUM(I)+RECOM_SUM(I) +
	1             COL_IR(I)+CHG_IR(I) +
	1             COL_RR(I)+CHG_RR(I) + 
	1             ABS(XRAY_RR(I)) +
	1             ABS(NT_IR(I)) +
	1             ABS(ADVEC_RR(I)) )/2.0D0
	   PHOT_SUM(I)=-PHOT_SUM(I)
	   COL_IR(I)=-COL_IR(I)
	   CHG_IR(I)=-CHG_IR(I)
	   NT_IR(I)=NT_IR(I)
!
	   RECOM_SUM(I)=RECOM_SUM(I)
	   XRAY_RR(I)=XRAY_RR(I)
	   COL_RR(I)=COL_RR(I)
	   CHG_RR(I)=CHG_RR(I)
	   ADVEC_RR(I)=ADVEC_RR(I)
	END DO
	YLABEL='Normalized rate'
	XVEC(1:ND)=V(1:ND)
	XLABEL='V(km/s)'
!
2000	CONTINUE
	OPTION='XN'
	CALL GEN_IN(OPTION,'Plot option: XN, R, V, T, ED, EX(stop), NS (new species),RATES')
	OPTION=UC(OPTION)
	WRITE(6,*)'Option=',OPTION
	IF(OPTION(1:2) .EQ. 'XN')THEN
	  DO I=1,ND
	    XVEC(I)=I
	  END DO
	  XLABEL='Depth index'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'NS')THEN
	  GOTO 1000
	ELSE IF(OPTION .EQ. 'R')THEN
	  XVEC(1:ND)=R(1:ND)
	  XLABEL='Radius(10\u10\d cm)'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'V')THEN
	  XVEC(1:ND)=V(1:ND)
	  XLABEL='V(km/s)'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'T')THEN
	  XVEC(1:ND)=TEMP(1:ND)
	  XLABEL='T(10\u4\dK)'
	ELSE IF(OPTION .EQ. 'ED')THEN
	  XVEC(1:ND)=ED(1:ND)
	  XLABEL='Ne(cm\u-3\d)'
	ELSE IF(OPTION(1:2) .EQ. 'EX')THEN
	  STOP
!
	ELSE IF(OPTION(1:5) .EQ. 'RATES')THEN
!
	  NORM=.TRUE.;  ABS_VALUE=.FALSE.
	  CALL GEN_IN(NORM,'Normalize rates')
	  CALL GEN_IN(ABS_VALUE,'Plot absolute values')
	  WRITE(6,*)' '
	  WRITE(6,*)'Calling DP_CURVE'
	  WRITE(6,*)'                : +ve means recombination'
	  WRITE(6,*)'                : -ve means ionizing'
	  WRITE(6,*)' '
!
	  YLABEL='Rates'; I=1
	  CALL PLT_RATE(XVEC,PHOT_SUM,TOT_SUM,'Photoionization rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,RECOM_SUM,TOT_SUM,'Recombination rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,NT_IR,TOT_SUM,'Non-thermal ionization rat',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,ADVEC_RR,TOT_SUM,'Advection rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,CHG_IR,TOT_SUM,'Charge exch. ionization rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,CHG_RR,TOT_SUM,'Charge exch. recombination rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,COL_IR,TOT_SUM,'Collisional ionization rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,COL_RR,TOT_SUM,'Collisional recombination rate',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,XRAY_RR,TOT_SUM,'Net X-ray recombination rate',NORM,ABS_VALUE,I,ND)
!
	ELSE IF(OPTION .EQ. 'P')THEN
	  CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	ELSE IF(OPTION(1:5) .EQ. 'WEIRD')THEN
	  WRITE(6,*)ED
	  WRITE(6,*)DI
	  J=1
	  DO WHILE(J .NE. 0)
	    CALL GEN_IN(J,'Level - zero to cease')
	    IF(J .NE. 0)THEN
	      YVEC(1:ND)=1.0D+14*PHOT(1:ND,J)/ED/DI
	      CALL DP_CURVE(ND,XVEC,YVEC)
	      YVEC(1:ND)=1.0D+14*RECOM(1:ND,J)/ED/DI
	      CALL DP_CURVE(ND,XVEC,YVEC)
	    END IF
	  END DO
	  CALL GRAMON_PGPLOT(XLABEL,'1.0D+14\ga',' ',' ')
	  GOTO 2000
	ELSE
	  WRITE(6,*)'Unrecognized X-axis option'
	  GOTO 2000
	END IF
	GOTO 2000
!
	END
!
	SUBROUTINE PLT_RATE(XVEC,PHOT_SUM,TOT_SUM,DESC,NORM,ABS_VALUE,IP,ND)
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IP
!
	REAL*8 XVEC(ND)
	REAL*8 PHOT_SUM(ND)
	REAL*8 TOT_SUM(ND)
	REAL*8 YVEC(ND)
	CHARACTER(LEN=*) DESC
!
	INTEGER I
	LOGICAL NORM
	LOGICAL ABS_VALUE
!
	IF(MAXVAL(ABS(PHOT_SUM)) .EQ. 0.0D0)RETURN
!
	IF(NORM)THEN
	  YVEC=PHOT_SUM/TOT_SUM
	  CALL DP_CURVE(ND,XVEC,YVEC)
	ELSE IF(ABS_VALUE)THEN
	  YVEC=-40
	  DO I=1,ND
	    IF(PHOT_SUM(I) .GT. 0.0D0)YVEC(I)=LOG10(PHOT_SUM(I))
	  END DO
	  CALL DP_CURVE(ND,XVEC,YVEC)
	ELSE
	  YVEC=PHOT_SUM
	  CALL DP_CURVE(ND,XVEC,YVEC)
	END IF
	WRITE(6,'(A,I2,2A)')PG_PEN(IP+1)//'   Curve ',IP,TRIM(DESC),DEF_PEN
	IP=IP+1
!
	RETURN
	END
