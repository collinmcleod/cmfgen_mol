!
! Routine reads in a set of nodes and clumping factors at those nodes.
! The routine then uses monotonic interpolation to create the model
! clumping factors. At presnet the clumping factors are assumed to
! be functions of velocity.
!
        SUBROUTINE SPL_CLUMP(CLUMP_FAC,VEL,ND)
	IMPLICIT NONE
!
! Created: 30-Jan-2022 : Osiris
!
	INTEGER ND
	REAL(10) CLUMP_FAC(ND)
	REAL(10) VEL(ND)
!
	INTEGER NPAR
	REAL(10), ALLOCATABLE :: VNODE(:)
	REAL(10), ALLOCATABLE :: FVAL(:)
!
	REAL(10) T1
	INTEGER IOS
	INTEGER LU,I,J
!
	LU=10
	OPEN(UNIT=LU,FILE='CLUMP_NODES',STATUS='OLD',IOSTAT=IOS,ACTION='READ')
	READ(LU,*)NPAR
	NPAR=NPAR+2
	ALLOCATE(VNODE(NPAR),FVAL(NPAR))
	DO I=2,NPAR-1
	  READ(LU,*)VNODE(I),FVAL(I)  
	END DO
	CLOSE(UNIT=10)
!
	IF(VNODE(2) .LT. VNODE(3))THEN
	  DO I=1,NPAR/2
	    J=NPAR-I+1
	    T1=VNODE(I); VNODE(I)=VNODE(J); VNODE(J)=T1
	    T1=FVAL(I);  FVAL(I)=FVAL(J);   FVAL(J)=T1
	  END DO
	END IF
	VNODE(1)=VEL(1); FVAL(1)=FVAL(2)
	VNODE(NPAR)=VEL(ND); FVAL(NPAR)=FVAL(NPAR-1)
!
	I=1
	CALL MON_INTERP(CLUMP_FAC,ND,I,VEL,ND,FVAL,NPAR,VNODE,NPAR) 
	DEALLOCATE(VNODE,FVAL)
!
	RETURN	
	END
