!
! Routine to define a continum using a cursor. The cursor is used to define
! both the X and Y values for the continuum. The routine can be called
! multiple times. Options are avilable to add new nodes, replace nodes
! or delete nodes.
!
! On the first call, the memory associated with IP need not be nallocated.
!
! This routine can also be used to update corrupted pixels.
!
	SUBROUTINE PG_MOD_CONT_NODES(IP)
	USE MOD_COLOR_PEN_DEF
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
	INTEGER IP           !Address of vector as in CD(IP)%.. (MOD_CURVE_DATA)
!
! Local vectors
!
	REAL*4, ALLOCATABLE :: TMP_XVEC(:)
	REAL*4, ALLOCATABLE :: XVEC(:)
	REAL*4, ALLOCATABLE :: YVEC(:)
!
	REAL*4 XVAL,YVAL,SYMB_EXP_FAC
	REAL*4 T1
	INTEGER I,J
	INTEGER ND		!Curren number of points in continuum vector
	INTEGER NX		!Mximum number of points in continuum vector
        INTEGER PGCURS
        INTEGER CURSERR
!
	INTEGER, PARAMETER :: IONE=1
	CHARACTER(LEN=1) CURSVAL
!
! Allocate XVEC and YVEC which will be modifeied as we set the
! continuum locations. The two CD vectors will only be modied
! as we exit the routine.
!
	ND=NPTS(IP)
	NX=NPTS(IP)+100
	ALLOCATE (XVEC(NX),YVEC(NX),TMP_XVEC(NX))
	XVEC=0; YVEC=0
	IF(ND .NE. 0)XVEC=CD(IP)%XVEC(1:ND)
	IF(ND .NE. 0)YVEC=CD(IP)%DATA(1:ND)
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')' '
	WRITE(6,*)'Click on top bar of plot window to activate cursor'
	WRITE(6,'(A)')' '
	WRITE(6,*)'Use ''r'' to replace a data point'
	WRITE(6,*)'Use ''a'' to add a data point'
	WRITE(6,*)'Use ''d'' to delete a data point'
	WRITE(6,*)'Use ''e'' to exit'
	WRITE(6,'(A)')DEF_PEN
!
	XVAL=5.0; CALL PGSCH(XVAL)
	DO WHILE(1 .EQ. 1)				!Multiple cursor entries
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  WRITE(6,*)'Cursor values are:',XVAL,YVAL
	  IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
!
	  IF(CURSVAL .EQ. 'r' .OR. CURSVAL .EQ. 'R')THEN
	    T1=XVAL
	    TMP_XVEC(1:ND)=ABS(XVEC(1:ND)-XVAL)
	    I=MINLOC(TMP_XVEC(1:ND),IONE)
	    YVEC(I)=YVAL
	    J=1; CALL PGSCI(J)
!	    SYMB_EXP_FAC=5.0; CALL PGSCH(SYMB_EXP_FAC)
	    CALL PGPT(IONE,XVAL,YVAL,IONE)
	    WRITE(6,*)'Replaced X,Y for X=',XVEC(I)
!
	  ELSE IF(CURSVAL .EQ. 'd' .OR. CURSVAL .EQ. 'D')THEN
	    T1=XVAL
	    TMP_XVEC(1:ND)=ABS(XVEC(1:ND)-T1)
	    I=MINLOC(TMP_XVEC(1:ND),IONE)
	    XVEC(I:ND-1)=XVEC(I+1:ND)
	    YVEC(I:ND-1)=YVEC(I+1:ND)
	    ND=ND-1
	    WRITE(6,*)'Deleted X=',XVEC(I),' from grid'
!
	  ELSE IF(CURSVAL .EQ. 'a' .OR. CURSVAL .EQ. 'A')THEN
	    CALL PGPT(IONE,XVAL,YVAL,IONE)
	    IF(ND+1 .GT. NX)THEN
	      WRITE(6,*)'Insuffient storage to add new data point'
	      WRITE(6,*)'Use e option, and renter routine'
	      WRITE(6,*)'This will allow more nodes to be added'
	      EXIT		!Exit DO_WHILE and got to clean up.
	    END IF
	    IF(ND .EQ. 0)THEN
	      ND=1
	      XVEC(1)=XVAL
	      YVEC(1)=YVAL
	    ELSE IF(XVAL .LT. XVEC(1))THEN
	      DO I=1,ND
	        XVEC(I+1)=XVEC(I-1)
	        YVEC(I+1)=YVEC(I-1)
	      END DO
	      XVEC(1)=XVAL;  YVEC(1)=YVAL
	      ND=ND+1
	    ELSE IF(XVAL .GT. XVEC(ND))THEN
	      XVEC(ND+1)=XVAL;  YVEC(ND+1)=YVAL
	      ND=ND+1
	    ELSE    
	      DO I=1,ND-1
	        IF( (XVEC(I)-XVAL)*(XVEC(I+1)-XVAL) .LT. 0 )THEN
	          DO J=ND,I+1,-1
                        XVEC(J+1)=XVEC(J)
                        YVEC(J+1)=YVEC(J)
	          END DO
	          XVEC(I+1)=XVAL; YVEC(I+1)=YVAL
	          ND=ND+1
	          J=1; CALL PGSCI(J)
	          EXIT
	        END IF
	      END DO
	    END IF
	    WRITE(6,*)'Added X=',XVAL,' to grid'
	  ELSE
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'Error - use r(eplace), a(dd), d(elete), e'
	    WRITE(6,*)DEF_PEN
	  END IF
	END DO
!
	IF(ND .NE. NPTS(IP))THEN
	  IF(ALLOCATED(CD(IP)%XVEC))DEALLOCATE (CD(IP)%XVEC,CD(IP)%DATA)
	  ALLOCATE (CD(IP)%XVEC(ND),CD(IP)%DATA(ND))
	  CD(IP)%XVEC(1:ND)=XVEC(1:ND)
	  CD(IP)%DATA(1:ND)=YVEC(1:ND)
	  NPTS(IP)=ND
	END IF
	DEALLOCATE (XVEC,YVEC,TMP_XVEC)
	
	RETURN
	END
