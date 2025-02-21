!
! Program reads in population estimates from a normal DC file,
! and uses linear interpolation in the log plane to deduce T and ED.
! If R(model) > R(model), R is scaled to match.
!
	SUBROUTINE REGRID_T_ED_V2(R,ED,T,POPATOM,VOL_EXP_FAC,ND,FILNAME)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 20-Sep-2006 : Based on REGRIDWSC.
!
	INTEGER ND
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) POPATOM(ND)
	REAL(KIND=LDP) VOL_EXP_FAC(ND)
!
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDED(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDT(:)
	REAL(KIND=LDP), ALLOCATABLE :: RLOG(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDR(:)
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local Variables.
!
	INTEGER I,J,IOS
	INTEGER NOLD,NDOLD
	INTEGER NX,NXST
	REAL(KIND=LDP) RPOLD,T1
	CHARACTER*80 STRING
	CHARACTER*(*) FILNAME
!
	LUER=ERROR_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=8,STATUS='OLD',FILE=FILNAME,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRIDWSC'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. CLUMP_FAC is presently not uses in
! this routine, as we regrid in R.
!
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(8,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(8)
!
	  READ(8,*)RPOLD,T1,NOLD,NDOLD
!
! Allocate required memory.
!
	  I=MAX(NDOLD,NOLD,ND)
	  ALLOCATE (OLDED(I))
	  ALLOCATE (OLDT(I))
	  ALLOCATE (RLOG(I))
	  ALLOCATE (TA(I))
	  ALLOCATE (OLDR(I))
!
	  DO I=1,NDOLD
	    READ(8,*)OLDR(I),TA(I),OLDED(I),OLDT(I)
	    READ(8,*)(TA(J),J=1,NOLD)
	  END DO
	CLOSE(UNIT=8)
!
	IF(ABS(OLDR(NDOLD)/R(ND)-1.0_LDP) .GT. 0.0001_LDP)THEN
	  WRITE(LUER,*)'Warning - core radius not identical in REGRIDWSC'
	  WRITE(LUER,*)'Rescaling to make Rcore identical'
	  DO I=1,NDOLD
	    OLDR(I)=R(ND)*( OLDR(I)/OLDR(NDOLD) )
	  END DO
	  OLDR(NDOLD)=R(ND)
	ELSE
	  OLDR(NDOLD)=R(ND)
	END IF
	IF( ABS(1.0_LDP-OLDR(1)/R(1)) .LE. 1.0E-10_LDP )OLDR(1)=R(1)
	IF(OLDR(2) .GE. OLDR(1))THEN
	  WRITE(LUER,*)'Reset OLDR(1) in REGRID_T_ED but now OLDR(2) .GE. OLDR(1))'
	  STOP
	END IF
!
	NXST=1
	DO WHILE(R(NXST) .GT. OLDR(1))
	  NXST=NXST+1
	END DO
	NX=ND-NXST+1
!
! Interpolations are performed in the log plane.
!
	DO I=1,ND
	  RLOG(I)=LOG(R(I))
	END DO
	DO I=1,NDOLD
	  OLDR(I)=LOG(OLDR(I))
	  OLDT(I)=LOG(OLDT(I))
	  OLDED(I)=LOG(OLDED(I))
	END DO
!
	CALL LINPOP(RLOG(NXST),T(NXST),NX,OLDR,OLDT,NDOLD)
	CALL LINPOP(RLOG(NXST),ED(NXST),NX,OLDR,OLDED,NDOLD)
!
	DO I=NXST,ND
	  T(I)=EXP(T(I))
	  ED(I)=EXP(ED(I))/VOL_EXP_FAC(I)
	END DO
	DO I=1,NXST-1
	  T(I)=T(NXST)
	  ED(I)=ED(NXST)*POPATOM(I)/POPATOM(NXST)
	END DO
!
! Free up memory.
!
	DEALLOCATE (TA)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDT)
	DEALLOCATE (RLOG)
	DEALLOCATE (OLDR)
!
	RETURN
	END
