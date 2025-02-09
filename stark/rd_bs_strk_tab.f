	MODULE STRK_MOD_HHE
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER STKTA_NWS		!Number of frequencies prof. tabulated at.
	INTEGER STKTA_NTS		!Number of temperatues prof. tabulated at.
	INTEGER STKTA_NES		!Number of Ne densities prof. tabulated at.
	INTEGER STKTA_NPS		!Number of Ne densities prof. tabulated at.
!
	REAL(KIND=LDP), ALLOCATABLE :: STKTA_DWS(:)	!Offset from line center (Ang)
	REAL(KIND=LDP), ALLOCATABLE :: STKTA_TS(:)	!Temperature (Log T(K))
	REAL(KIND=LDP), ALLOCATABLE :: STKTA_ES(:)	!Log Ne (cgs units)
	REAL(KIND=LDP), ALLOCATABLE :: STKTA_PS(:)	!Phi(v)
!
	INTEGER NL_STRK	!Lower levels
	INTEGER NUP_STRK	!Upper level
!
	LOGICAL STKTA_QHALF	!If TRUE, on half of profile tabulated
	LOGICAL STKTA_WSCA      ! Lambda (in A) o
!
	END MODULE STRK_MOD_HHE
!
!
!
	SUBROUTINE RD_BS_STRK_TAB(SPECIES,NU_ZERO,NL,NUP,PROF_TYPE,PROF_ID,LUSTK)
	USE SET_KIND_MODULE
	USE STRK_MOD_HHE
	IMPLICIT NONE
!
	REAL(KIND=LDP) NU_ZERO
	CHARACTER*(*) SPECIES
	CHARACTER*(*) PROF_TYPE
	INTEGER NL,NUP
	INTEGER PROF_ID
	INTEGER LUSTK
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	CHARACTER*3 LOCAL_SPECIES
	INTEGER LTBL
	INTEGER I
!
! I profile was found, we read in the tabulated data
!
	LOCAL_SPECIES=' '
	OPEN(UNIT=LUSTK,FILE=PROF_TYPE,FORM='FORMATTED',
     &               STATUS='OLD',ACTION='READ')
          DO LTBL = 1,PROF_ID
             READ(LUSTK,*)LOCAL_SPECIES(1:3),
     &               NL_STRK,NUP_STRK,STKTA_QHALF,STKTA_WSCA
            READ(LUSTK,*) STKTA_NWS,STKTA_NTS,
     &               STKTA_NES,STKTA_NPS!
!
	    IF(ALLOCATED(STKTA_PS))THEN
	      IF(SIZE(STKTA_DWS) .LT. STKTA_NPS)DEALLOCATE(STKTA_DWS)
	      IF(SIZE(STKTA_TS) .LT. STKTA_NPS)DEALLOCATE(STKTA_TS)
	      IF(SIZE(STKTA_ES) .LT. STKTA_NPS)DEALLOCATE(STKTA_ES)
	      IF(SIZE(STKTA_PS) .LT. STKTA_NPS)DEALLOCATE(STKTA_PS)
	    END IF
            IF(.NOT. ALLOCATED(STKTA_DWS))ALLOCATE(STKTA_DWS(STKTA_NWS))
            IF(.NOT. ALLOCATED(STKTA_TS))ALLOCATE(STKTA_TS(STKTA_NTS))
            IF(.NOT. ALLOCATED(STKTA_ES))ALLOCATE(STKTA_ES(STKTA_NES))
            IF(.NOT. ALLOCATED(STKTA_PS))ALLOCATE(STKTA_PS(STKTA_NPS))
!
            READ(LUSTK,*) (STKTA_DWS(I), I=1,STKTA_NWS)
            READ(LUSTK,*) (STKTA_TS(I), I=1,STKTA_NTS)
            READ(LUSTK,*) (STKTA_ES(I), I=1,STKTA_NES)
            READ(LUSTK,*) (STKTA_PS(I), I=1,STKTA_NPS)
          END DO
	CLOSE(UNIT=LUSTK)
!
! Check that the data read in is that actually required.
!
	IF(SPECIES .NE. LOCAL_SPECIES .OR. NL .NE. NL_STRK
	1   .OR. NUP .NE. NUP_STRK)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in RD_BS_STRK_TAB'
	  WRITE(I,*)'Lines don''t match'
	  WRITE(I,*)'SPECIES=',TRIM(SPECIES),'NL=',NL,'NUP=',NUP
	  WRITE(I,*)'LST_SPE=',TRIM(LOCAL_SPECIES),'NL=',NL_STRK,'NUP=',NUP_STRK
	END IF
!
	RETURN
	END
