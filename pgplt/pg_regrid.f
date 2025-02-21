	SUBROUTINE PG_REGRID(IP,OP,XMIN,XMAX)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
	INTEGER IP,OP
	REAL*4 XMIN,XMAX
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: YV(:)
!
	REAL*4 DX
	REAL*4 T1
	INTEGER N
	INTEGER I
	INTEGER IL,IU
	INTEGER, PARAMETER :: IONE=1
!
	IL=1; IU=NPTS(IP)
	IF(CD(IP)%XVEC(1) .LT. CD(IP)%XVEC(2))THEN
	  T1=MIN(XMIN,XMAX)
	  DO I=1,NPTS(IP)
	    IF(T1 .GE. CD(IP)%XVEC(I))THEN
	      IL=I; EXIT
	    END IF
	  END DO
	  T1=MAX(XMIN,XMAX)
	  DO I=IL+1,NPTS(IP)
	    IF(T1 .LT. CD(IP)%XVEC(I))THEN
	      IU=I-1; EXIT
	    END IF
	  END DO
	ELSE
	  T1=MIN(XMIN,XMAX)
	  DO I=NPTS(IP),1,-1
	    IF(T1 .GE. CD(IP)%XVEC(I))THEN
	      IU=I; EXIT
	    END IF
	  END DO
	  T1=MAX(XMIN,XMAX)
	  DO I=IU-1,1,-1
	    IF(T1 .LT. CD(IP)%XVEC(I))THEN
	      IL=I+1; EXIT
	    END IF
	  END DO
	END IF
	N=IU-IL+1
	WRITE(6,*)IL,IU,N
!
	IF(ALLOCATED(XV))DEALLOCATE(XV,YV)
	ALLOCATE(XV(N),YV(N))
	DX=(CD(IP)%XVEC(IU)-CD(IP)%XVEC(IL))/(N-1)
	DO I=1,N
	  XV(I)=CD(IP)%XVEC(IL)+(I-1)*DX
	END DO
!
	CALL SP_MON_INTERP(YV,N,IONE,XV,N,CD(IP)%DATA(IL),N,CD(IP)%XVEC(IL),N)
!
	IF(ALLOCATED(CD(OP)%XVEC))DEALLOCATE(CD(OP)%XVEC,CD(OP)%DATA)
	ALLOCATE(CD(OP)%XVEC(N),CD(OP)%DATA(N))
	CD(OP)%XVEC=XV
	CD(OP)%DATA=YV
	NPTS(OP)=N
	NPLTS=NPLTS+1
!
	RETURN
	END
