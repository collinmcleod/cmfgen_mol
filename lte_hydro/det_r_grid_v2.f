!
! Subroutine to create a grid, in TAU space, for solution of the transfer equation.
! Initial grid is created so that it satisfies:
!       V(I+1) > 0.67 V(I)				(set by V_SCL_FAC)
!       LoG(TAU(I+1)) < 0.25 + LOG(TAU(I))		(set by dLOG_TAU)
!
! This scaling will be modified slightly when the number of grid points is
! adjusted to the desired value.
!
	SUBROUTINE DET_R_GRID_V2(REV_TAU_GRID,NEW_ND,ND_MAX,TAU_MAX,
	1           dLOG_TAU,V_SCL_FAC,OUT_BND_OPT,OBND_PARS,NUM_OBND_PARAMS,
	1           R,V,TAU,ND)
	IMPLICIT NONE
!
! Altered 12-Feb-2007 : Modifications to allow more choice in the outer boundary 
!                         condition. dLOG_TAU and V_SCL_FAC now passed in call.
! Created 12-Aug-2006
!
	INTEGER NEW_ND				!Requested number of grid points in new grid.
	INTEGER ND_MAX				!Maximum number of grid points in grid.
	INTEGER NUM_OBND_PARAMS
	REAL*8 REV_TAU_GRID(NEW_ND)
	REAL*8 OBND_PARS(NUM_OBND_PARAMS)
	REAL*8 TAU_MAX				!Maximum optical depth for new grid.
	REAL*8 V_SCL_FAC			!Maxium values fo V(I+1)/V(I)
	REAL*8 dLOG_TAU				!Maximu value for /\Tau
	CHARACTER(LEN=*) OUT_BND_OPT
!
! These describe the old grid.
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 TAU(ND)
!
! Local arrays.
!
	REAL*8 REV_R(ND_MAX)
	REAL*8 REV_V(ND_MAX)
	REAL*8 REV_TAU(ND_MAX)
	REAL*8 OLD_R(ND_MAX)
	REAL*8 OLD_TAU(ND_MAX)
	REAL*8 LOG_TAU(ND)
!
! Local variables.
!
	REAL*8 LOG_TAU_MAX
	REAL*8 dTAU
	REAL*8 T1,T2
	INTEGER ND_TMP
	INTEGER I,J,K,JST
	INTEGER LUER,ERROR_LU
	INTEGER, PARAMETER :: IONE=1
	EXTERNAL ERROR_LU
!
	LOG_TAU=LOG(TAU)
	LOG_TAU_MAX=LOG(TAU_MAX)
!
	I=1
	REV_R(1)=R(1)
	REV_V(1)=V(1)
	REV_TAU(1)=LOG_TAU(1)
!
! Handle the grid near the oputer boundary, according to the option specified.
! For each option, we irst check parameter validity.
!
	LUER=ERROR_LU()
	WRITE(6,*)'OUT_BND_OPT=',OUT_BND_OPT
	IF(OUT_BND_OPT .EQ. 'POW')THEN
	  IF(NUM_OBND_PARAMS .NE. 2)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V2 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'Exactly 2 parameters must be specified.'
	    STOP
	  ELSE IF(NINT(OBND_PARS(1)) .LT. 1 .OR. OBND_PARS(2) .GT. 10)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V2 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'OBND_PARS(1) must be >=1 and <10'
	    STOP
	  ELSE IF(OBND_PARS(2) .LE. 1.0)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V2 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'OBND_PARS(2) must be >1'
	    STOP
	  END IF
	    WRITE(6,*)OBND_PARS(1)
	    WRITE(6,*)OBND_PARS(2)
	  J=NINT(OBND_PARS(1))
	  DO K=2,J+1
	    REV_TAU(K)=REV_TAU(1)+dLOG_TAU/(OBND_PARS(2)**(J+2-K))
	  END DO
	  REV_TAU(J+2)=REV_TAU(1)+dLOG_TAU
	  K=J+2
	ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	  DO K=1,NUM_OBND_PARAMS-1
	    IF(OBND_PARS(K) .LT. 1.0 .OR. OBND_PARS(K) .LE. OBND_PARS(K+1))THEN
	      WRITE(LUER,*)'Error in DET_R_GRID_V2'
	      WRITE(LUER,*)'OB parameters must be >1 and monotonically decreasing'
	      STOP
	    END IF
	  END DO
	  DO K=1,NUM_OBND_PARAMS
	    REV_TAU(K+1)=REV_TAU(1)+dLOG_TAU/OBND_PARS(K)
	  END DO
	  REV_TAU(NUM_OBND_PARAMS+2)=REV_TAU(1)+dLOG_TAU
	  K=NUM_OBND_PARAMS+2
	ELSE IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	  REV_TAU(2)=REV_TAU(1)+dLOG_TAU/9.0D0
	  REV_TAU(3)=REV_TAU(1)+dLOG_TAU/3.0D0
	  REV_TAU(4)=REV_TAU(1)+dLOG_TAU
	  K=4
	END IF
!
!Estimate V: accuracy at outer boundary not important.	
!
	IF(OUT_BND_OPT .NE. ' ' .AND. OUT_BND_OPT .NE. 'NONE')THEN
	  J=2
	  DO WHILE(REV_TAU(K) .GT. LOG_TAU(J))
	    J=J+1
	  END DO
	  DO I=2,K
	    T1=(REV_TAU(I)-LOG_TAU(1))/(LOG_TAU(J)-LOG_TAU(1))
	    REV_V(I)=(1.0D0-T1)*V(1)+T1*V(J)
	  END DO
	  I=K
	END IF
!
	J=2
 	JST=J
	DO WHILE(REV_TAU(I) .LT. LOG_TAU_MAX)
!
	  I=I+1
	  IF(I+2 .GT. ND_MAX)THEN
	    WRITE(6,*)'Error in DET_R_GIRD_V1 --- ND_MAX too small'
	    WRITE(6,*)'ND_MAX=',ND_MAX
	    STOP
	  END IF
!
! This step ensure TAU spacing criterin is satisfied.
!
	  REV_TAU(I)=REV_TAU(I-1)+dLOG_TAU
!
! Check if we are at inner boundary. If so, we decrease spacing systematically
! to give more accuracy for our first orer boundary conditions. We ignore the 
! velocity check.
!
	  IF(REV_TAU(I)+dLOG_TAU .GE. LOG_TAU_MAX)THEN
	    T1=LOG_TAU_MAX-REV_TAU(I-1)
	    REV_TAU(I)=REV_TAU(I-1)+0.6*T1
	    I=I+1
	    REV_TAU(I)=REV_TAU(I-2)+0.9*T1
	    I=I+1
	    REV_TAU(I)=LOG_TAU_MAX
	  ELSE
!
! Check if V spacing criterion is satisifed. If not, shrink spacing. Linear
! interpolation if fine since we are only constructing the grid, not actual values
! on the grid.
!
	    DO WHILE (REV_TAU(I) .GT. LOG_TAU(J))
	      IF(J .EQ. ND)EXIT
	      J=J+1
	    END DO
	    T1=(REV_TAU(I)-LOG_TAU(J-1))/(LOG_TAU(J)-LOG_TAU(J-1))
	    REV_V(I)=(1.0D0-T1)*V(J-1)+T1*V(J)
	    IF(REV_V(I) .GT. 0.1D0 .AND. REV_V(I) .LT. V_SCL_FAC*REV_V(I-1))THEN
	      J=JST
	      DO WHILE (V_SCL_FAC*REV_V(I-1) .LT. V(J))
	        IF(J .EQ. ND)EXIT
	        J=J+1
	       END DO
	      T1=(V_SCL_FAC*REV_V(I-1)-V(J-1))/(V(J)-V(J-1))
	      REV_TAU(I)=(1.0D0-T1)*LOG_TAU(J-1)+T1*LOG_TAU(J)
	      REV_V(I)=V_SCL_FAC*REV_V(I-1)
	    END IF
 	    J=JST
	  END IF
	  ND_TMP=I
	END DO
!
! The next few line can be uncommented if debuging.
!
!	WRITE(80,'(/,A,/)')' Estimate R grid as determined by DET_R_GRID_V2'
!	WRITE(80,'(A,2(13X,A))')'Index','R','V'
!	DO I=1,ND_TMP
!	  WRITE(80,'(I5,3ES14.4)')I,REV_V(I),REV_TAU(I)
!	END DO
!
! We now create a grid with the requested number of grid points. To do so
! we use the array index as the independent variable. This approach should
! preserve the basic grid spacing.
!
! The scaling for REV_R makes the X-limits identical for the new and old grids.
!
	OLD_TAU(1:ND_TMP)=REV_TAU(1:ND_TMP)
	DO I=1,ND_TMP; OLD_R(I)=I; END DO
	DO I=1,NEW_ND
	  REV_R(I)=1+((I-1.0D0)*(ND_TMP-1.0D0) )/(NEW_ND-1.0D0)
	END DO
	CALL MON_INTERP(REV_TAU,NEW_ND,IONE,REV_R,NEW_ND,OLD_TAU,ND_TMP,OLD_R,ND_TMP)
!
! We finish by remembering that we have been working with LOG(TAU).
!
	REV_TAU_GRID(1:NEW_ND)=EXP(REV_TAU(1:NEW_ND))
!
	RETURN
	END
