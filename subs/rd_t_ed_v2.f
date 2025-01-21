!
! Subroutine to read in T and ED from a DC file. This T and E will be used to
! set profile limits. Accuracy is not crucial.
!
	SUBROUTINE RD_T_ED_V2(R,T,ED,ND,LU_IN,FILE_NAME)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created: 04-Nov-2013
!
	INTEGER ND
	INTEGER LU_IN
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) T1,T2
!
	REAL(KIND=LDP), ALLOCATABLE :: OLD_R(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_T(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_ED(:)
!
	INTEGER LUER
	INTEGER ERROR_LU
	INTEGER I,J
	INTEGER IOS
	INTEGER NOLD,OLD_ND
	INTEGER, PARAMETER ::IONE=1
	EXTERNAL ERROR_LU
	CHARACTER(LEN=*) FILE_NAME
	CHARACTER(LEN=80) STRING
!
	LUER=ERROR_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=LU_IN,STATUS='OLD',FILE=FILE_NAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILE_NAME),' in RD_T_ED'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file.
!
        I=0
        STRING=' '
        DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
          I=I+1
          READ(LU_IN,'(A)')STRING
        END DO
        IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU_IN)
!
	READ(LU_IN,*)T1,T2,NOLD,OLD_ND
	IF(ND .NE. OLD_ND)THEN
	  WRITE(6,*)'Warning ND .NE. OLD_ND in RD_T_ED'
	  ALLOCATE(OLD_R(OLD_ND))
	  ALLOCATE(OLD_T(OLD_ND))
	  ALLOCATE(OLD_ED(OLD_ND))
	  DO I=1,OLD_ND
	    READ(LU_IN,*)OLD_R(I),T2,OLD_ED(I),OLD_T(I)
	    READ(LU_IN,*)(T1,J=1,NOLD)
	  END DO
	  OLD_R(1)=R(1); OLD_R(OLD_ND)=R(ND)
	  CALL MON_INTERP(T,ND,IONE,R,ND,OLD_T,OLD_ND,OLD_R,OLD_ND)
	  CALL MON_INTERP(ED,ND,IONE,R,ND,OLD_ED,OLD_ND,OLD_R,OLD_ND)
	  DEALLOCATE(OLD_R,OLD_T,OLD_ED)
	ELSE
	  DO I=1,OLD_ND
	    READ(LU_IN,*)T1,T2,ED(I),T(I)
	    READ(LU_IN,*)(T1,J=1,NOLD)
	  END DO
	END IF
	  CLOSE(LU_IN)
!
	RETURN
	END
