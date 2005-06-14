!
! Routine to read in analytical fits to partial photoionization
! cross-sections. Data is primarily need to handle ionizations by 
! X-rays, although all fits are included. 
! 
! As a consequnce, need to be carefull that photoionization routes
! are not doble counted.
!
! Fits are from the Verner atomic data site.
!
! Ref:
!     Verner, D., Yakovlev, D.G., A\&A, 1995???
!     Analytic fits for partial photoinization cross-sections.
!
	MODULE XRAY_DATA_MOD
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: ATNO_MAX_X=28
	INTEGER, PARAMETER :: PQN_MAX_X=3
	INTEGER, PARAMETER :: ANG_MAX_X=2
!
! We have added _X to make the names more unique. The notation is that
! of Verney and Yakovlev.
!
	REAL*8 E_THRESH_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL*8 E_0_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL*8 SIG_0_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL*8 y_a_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL*8 y_w_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL*8 P_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
!
	LOGICAL XRAY_PHOT_RD_IN
!
	END MODULE XRAY_DATA_MOD
!
!
!
! Subroutine read in all available data from the file XRAY_PHOT_FITS.
!
	SUBROUTINE RD_XRAY_FITS(LU)
	USE XRAY_DATA_MOD
	IMPLICIT NONE
!
	INTEGER LU
	INTEGER I,J,K,L,IOS
	INTEGER IZ,INE,PQN,ANG
!
	EXTERNAL ERROR_LU
	INTEGER LU_ER,ERROR_LU
	CHARACTER*80 STRING
!
	LU_ER=ERROR_LU()
	OPEN(UNIT=LU,FILE='XRAY_PHOT_FITS',STATUS='OLD',ACTION='READ',
	1             IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error opening XRAY_PHOT_FITS file'
	    WRITE(LU_ER,*)'IOSTAT=',IOS
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	  END DO
!
	  SIG_0_X(:,:,:,:)=0.0D0
	  DO WHILE(1 .EQ. 1)	
	     READ(STRING,*,IOSTAT=IOS)IZ,INE,PQN,ANG
	     IF(IOS .NE. 0)THEN
	       WRITE(LU_ER,*)'Error reading IZ,INE... from XRAY_PHOT_FITS'
	       WRITE(LU_ER,*)'IOSTAT=',IOS
	       STOP
	     END IF
	     IF(IZ .GT. ATNO_MAX_X)THEN
	       XRAY_PHOT_RD_IN=.TRUE.
	       EXIT
	     END IF
	     IF(PQN .GT. PQN_MAX_X)GOTO 50
	     IF(ANG .GT. ANG_MAX_X)GOTO 50
!
	     READ(STRING,*,IOSTAT=IOS)I,J,K,L,E_THRESH_X(IZ,INE,PQN,ANG),
	1                          E_0_X(IZ,INE,PQN,ANG),
	1                          SIG_0_X(IZ,INE,PQN,ANG),
	1                          Y_A_X(IZ,INE,PQN,ANG),
	1                          P_X(IZ,INE,PQN,ANG),
	1                          Y_W_X(IZ,INE,PQN,ANG)
	     IF(IOS .NE. 0)THEN
	       WRITE(LU_ER,*)'Error reading IZ,INE... from XRAY_PHOT_FITS'
	       WRITE(LU_ER,*)'IOSTAT=',IOS
	       STOP
	     END IF
!
50	    CONTINUE
	    READ(LU,'(A)')STRING
	  END DO
!
	RETURN
	END
