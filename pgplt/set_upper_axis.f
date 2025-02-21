!
	SUBROUTINE SET_UPPER_AXIS(ED,XV,ND,TOP_LAB)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) XV(ND)
	REAL(KIND=LDP) ED(ND)
	CHARACTER(LEN=*) TOP_LAB
!
	REAL(KIND=LDP) TA(ND),TB(ND)
	INTEGER I,J
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER, PARAMETER :: NSC=31
	COMMON/TOPBORD/ SCED(NSC),XED(NSC),NXED,TOPLABEL
	REAL(KIND=LDP) SCED,XED
	INTEGER NXED
	CHARACTER(LEN=30) TOPLABEL
	DATA SCED/2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,
	1         10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15,15.5,
	1         16,16.5,17.0/
!
	TA(1:ND)=XV(1:ND)
!
	I=LOG10(ED(1));J=LOG10(ED(ND))
	IF(I .LT. ED(1))I=I+1
	NXED=(J-I)*2+1
	NXED=MIN(NXED,NSC)
	SCED(1)=I
	DO I=2,NXED
	  SCED(I)=SCED(1)+(I-1)*0.5D0
	END DO
	DO I=1,ND
	  TB(I)=LOG10(ED(I))
	END DO
	CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	TOPLABEL=TOP_LAB
!
	RETURN
	END
