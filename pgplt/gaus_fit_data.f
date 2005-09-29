	MODULE GAUS_FIT_DATA
!
	INTEGER NG_PAR				!Number of parameters in fit
	INTEGER NG_DATA				!Number of data points
!
	REAL*8, ALLOCATABLE :: X_GAUS(:)	!Absica data as stored in module
	REAL*8, ALLOCATABLE :: Y_GAUS(:)	!Data to be fitted
!
	REAL*8, ALLOCATABLE :: SIM(:,:)		!Simplex (NG_PAR+1 parameter set estimates)
	REAL*8, ALLOCATABLE :: PAR(:)		!Single parameter set
	REAL*8, ALLOCATABLE :: SUM_SQ(:)
!
	REAL*4, ALLOCATABLE :: XFIT(:)		!Same as X_GAUS but for PGPLOT routines
	REAL*4, ALLOCATABLE :: YFIT(:)		!Function fit (evaluated by GAUS_FIT_FUNC)
!
	END MODULE GAUS_FIT_DATA
!
! This subroutined defines the approriate arrrays ans store the data for the
! Gaussin fitting routine. Note the XVEC and YVEC are assumed to be single
! precision as routines is for use with GRAMON_PGPLOT.
!
	SUBROUTINE SET_GAUS_DATA(XVEC,YVEC,XST,XEND,NX,N_PAR,YST,YEND,NDATA)
	USE GAUS_FIT_DATA
	IMPLICIT NONE
!
	INTEGER NX
	INTEGER N_PAR
	INTEGER NDATA
	REAL*4 XVEC(NX)
	REAL*4 YVEC(NX)
	REAL*8 XST,XEND
	REAL*8 YST,YEND
!
	INTEGER I
	INTEGER IXST,IXEND
!
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
!
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
	NDATA=NG_DATA
	NG_PAR=N_PAR
	IF(ALLOCATED(X_GAUS))DEALLOCATE(X_GAUS)
	IF(ALLOCATED(Y_GAUS))DEALLOCATE(Y_GAUS)
	IF(ALLOCATED(XFIT))DEALLOCATE(XFIT)
	IF(ALLOCATED(YFIT))DEALLOCATE(YFIT)
	ALLOCATE (X_GAUS(NG_DATA),Y_GAUS(NG_DATA))
	ALLOCATE (XFIT(NG_DATA),YFIT(NG_DATA))
	DO I=IXST,IXEND
	  X_GAUS(I-IXST+1)=XVEC(I)
	  Y_GAUS(I-IXST+1)=YVEC(I)
	END DO
!
! Set X-values for later plotting
!
	XFIT(1:NG_DATA)=X_GAUS(1:NG_DATA)
!
	RETURN
	END
