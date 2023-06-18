	MODULE GAUSS_FIT_DATA
!
! Altered 06-Mar-2023 : YCONT_FIT added
! Altered 09-Aug-2022 : To get consistency in the different routines changed to use Gauss.
!
	INTEGER NUM_GAUSS			!Number of Gaussians
	INTEGER NG_PAR				!Total number of parameters in fit
	INTEGER NG_PAR_MAX		 	!Maximum total number of parameters in fit
	INTEGER NG_DATA				!Number of data points
!
	REAL*8, ALLOCATABLE :: X_GAUSS(:)	!Absica data as stored in module
	REAL*8, ALLOCATABLE :: Y_GAUSS(:)	!Data to be fitted
!
	REAL*8, ALLOCATABLE :: SIM(:,:)		!Simplex (NG_PAR+1 parameter set estimates)
	REAL*8, ALLOCATABLE :: PAR(:)		!Single parameter set
	REAL*8, ALLOCATABLE :: SUM_SQ(:)
	REAL*8, ALLOCATABLE :: SCALE(:)
	REAL*8, ALLOCATABLE :: EW_TABLE(:,:)
	REAL*8, ALLOCATABLE :: LAM_TABLE(:,:)
	REAL*8, ALLOCATABLE :: EW(:)
	REAL*8, ALLOCATABLE :: EW_CONT(:)
	REAL*8, ALLOCATABLE :: EW_ERROR(:)
	REAL*8, ALLOCATABLE :: ALT_ERROR(:)
	REAL*8, ALLOCATABLE :: MIN_ERROR(:)
	INTEGER, ALLOCATABLE :: INDX_VEC(:)
!
	REAL*4, ALLOCATABLE :: XFIT(:)		!Same as X_GAUSS but for PGPLOT routines
	REAL*4, ALLOCATABLE :: YFIT(:)		!Function fit (evaluated by GAUSS_FIT_FUNC)
	REAL*4, ALLOCATABLE :: YCONT_FIT(:)	!Function fit (evaluated by GAUSS_FIT_FUNC)
	REAL*4, ALLOCATABLE :: FIT_DIF(:)	!Function fit (evaluated by GAUSS_FIT_FUNC)
	SAVE
!
	END MODULE GAUSS_FIT_DATA
!
! This subroutined defines the approriate arrrays and stores the data for the
! Gaussin fitting routine. Note the XVEC and YVEC are assumed to be single
! precision as routines is for use with GRAMON_PGPLOT.
!
	SUBROUTINE SET_GAUSS_DATA(XVEC,YVEC,XST,XEND,NX,YST,YEND)
	USE GAUSS_FIT_DATA
	IMPLICIT NONE
!
	INTEGER NX
	REAL*4 XVEC(NX)
	REAL*4 YVEC(NX)
	REAL*8 XST,XEND
	REAL*8 YST,YEND
!
	INTEGER I
	INTEGER IXST,IXEND
!
	IXST=0
	DO I=1,NX-1
	  IF( (XST-XVEC(I))*(XVEC(I+1)-XST) .GE. 0)THEN
	    IF(ABS(XST-XVEC(I)) .GT. ABS(XVEC(I+1)-XST) )THEN 
	      IXST=I
	      EXIT
	    ELSE
	      IXST=I+1
	      EXIT
	    END IF
	  END IF
	END DO
	IF(IXST .EQ. 0)THEN
	  WRITE(6,*)'Error in SET_GAUSS_DATA: XST out of range'
	  WRITE(6,'(X,3(A,E15.8,3X))')'XST=',XST,'XVEC(1)=',XVEC(1),'XVEC(NX)=',XVEC(NX)
	  RETURN
	END IF
!
	IXEND=0
	DO I=1,NX-1
	  IF( (XEND-XVEC(I))*(XVEC(I+1)-XEND) .GE. 0)THEN
	    IF(ABS(XEND-XVEC(I)) .GT. ABS(XVEC(I+1)-XEND) )THEN 
	      IXEND=I
	      EXIT
	    ELSE
	      IXEND=I+1
	      EXIT
	    END IF
	  END IF
	END DO
	IF(IXEND .EQ. 0)THEN
	  WRITE(6,*)'Error in SET_GAUSS_DATA: XEND out of range'
	  WRITE(6,'(X,3(A,E15.8,3X))')'XEND=',XEND,'XVEC(1)=',XVEC(1),'XVEC(NX)=',XVEC(NX)
	  RETURN
	END IF
!
	IF(IXST .GT. IXEND)THEN
	  I=IXEND
	  IXEND=IXST
	  IXST=I
	END IF
	YST=YVEC(IXST)
	YEND=YVEC(IXEND)
!
	NG_DATA=IXEND-IXST+1
	IF(ALLOCATED(X_GAUSS))DEALLOCATE(X_GAUSS)
	IF(ALLOCATED(Y_GAUSS))DEALLOCATE(Y_GAUSS)
	IF(ALLOCATED(XFIT))DEALLOCATE(XFIT)
	IF(ALLOCATED(YFIT))DEALLOCATE(YFIT)
	IF(ALLOCATED(YCONT_FIT))DEALLOCATE(YCONT_FIT)
	IF(ALLOCATED(FIT_DIF))DEALLOCATE(FIT_DIF)
	ALLOCATE (X_GAUSS(NG_DATA),Y_GAUSS(NG_DATA))
	ALLOCATE (XFIT(NG_DATA),YFIT(NG_DATA))
	ALLOCATE (YCONT_FIT(NG_DATA),FIT_DIF(NG_DATA))
	DO I=IXST,IXEND
	  X_GAUSS(I-IXST+1)=XVEC(I)
	  Y_GAUSS(I-IXST+1)=YVEC(I)
	END DO
!
! Set X-values for later plotting
!
	XFIT(1:NG_DATA)=X_GAUSS(1:NG_DATA)
!
	RETURN
	END
