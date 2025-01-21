	PROGRAM TST_WR
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER, PARAMETER :: ND=60
	INTEGER, PARAMETER :: LU=10
	INTEGER, PARAMETER :: ISIX=6
	INTEGER I
	REAL(KIND=LDP) F(ND)
!
	DO I=1,ND
	  F(I)=I
	END DO
!
	CALL WRITV_V2(F,ND,ISIX,'Testing model',LU)
!
	STOP
	END
