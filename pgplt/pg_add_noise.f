!
! Subroutine to interpolate a data vector in PGPLOT to a unirom grid in
! X space or LOG(X) space.
!
	  SUBROUTINE PG_ADD_NOISE(IN1,OUT,XMIN,XMAX,COUNTS,R_SEED)
	  USE MOD_CURVE_DATA
	  IMPLICIT NONE
!
! Created 16-Nov-2020 
!
	  INTEGER IN1
	  INTEGER OUT
	  REAL*4 XMIN,XMAX
	  REAL*4 COUNTS
	  REAL*4 R_SEED
	  REAL(10) LAM
!
	  REAL(10) POIDEV
	  INTEGER GET_INDX_SP
	  EXTERNAL GET_INDX_SP,POIDEV
!         
! Local variables
!
	  REAL*4, ALLOCATABLE :: XV(:)
	  REAL*4, ALLOCATABLE :: YV(:)
	  REAL*4, ALLOCATABLE :: XV_NEW(:)
	  REAL*4, ALLOCATABLE :: YV_NEW(:)
!
	  INTEGER N1
	  INTEGER N2
!
	  REAL*4 T1
	  INTEGER I,J,L
	  INTEGER IRAN
	  INTEGER LOW_LIM,UP_LIM
!
	  INTEGER, PARAMETER :: IONE=1
	  INTEGER, PARAMETER :: T_OUT=6
!
! Check validity of maps:
!
	  IF(NPTS(IN1) .EQ. 0)THEN
	    WRITE(T_OUT,*)'Invalid Plot ID''s in DO_VEC_OP'
	    WRITE(T_OUT,*)'Plot ID=',IN1,'  NPTS=',NPTS(IN1)
	    RETURN
	  END IF
!
	  IF(COUNTS .LE. 0.0D0)THEN
	    WRITE(T_OUT,*)'Error -- count smust be > 0'
	    RETURN
	  END IF
!
	  N1=NPTS(IN1)
	  IF(XMIN .EQ. XMAX)THEN
	    LOW_LIM=1; UP_LIM=N1
	  ELSE
	    LOW_LIM=GET_INDX_SP(XMIN,CD(IN1)%XVEC,N1)
	    UP_LIM=GET_INDX_SP(XMAX,CD(IN1)%XVEC,N1)
	    J=LOW_LIM
	    LOW_LIM=MIN(LOW_LIM,UP_LIM)
	    IF(J .NE. LOW_LIM)UP_LIM=J
	  END IF
	  WRITE(6,'(A,2I8)')'LOW_LIM,UP_LIM',LOW_LIM,UP_LIM
!
	  ALLOCATE(XV(NPTS(IN1)))
	  ALLOCATE(YV(NPTS(IN1)))
	  N2=UP_LIM-LOW_LIM+1
	  XV(1:N2)=CD(IN1)%XVEC(LOW_LIM:UP_LIM)
	  YV(1:N2)=CD(IN1)%DATA(LOW_LIM:UP_LIM)
	  WRITE(6,*)'Set old vectors'
!
	  IRAN=-R_SEED*1234567
	  DO I=1,N2
	    LAM=YV(I)*COUNTS
	    YV(I)=POIDEV(LAM,IRAN)/COUNTS
	  END DO
	  NPTS(OUT)=N2
	  ERR(OUT)=.FALSE.
	  IF(OUT .GT. NPLTS)NPLTS=OUT
!
	  IF(ALLOCATED(CD(OUT)%XVEC))THEN
	    DEALLOCATE(CD(OUT)%DATA,CD(OUT)%XVEC)
	  END IF
	  ALLOCATE(CD(OUT)%DATA(N2),CD(OUT)%XVEC(N2))
!
	  WRITE(6,*)'Number of data points in regridded data',N2
	  CD(OUT)%XVEC(1:N2)=XV(1:N2)
	  CD(OUT)%DATA(1:N2)=YV(1:N2)
	  DEALLOCATE(XV,YV)
!
	  RETURN
	  END
