	PROGRAM TST_GET_INDX
	USE SET_KIND_MODULE
	USE GET_INDX_INTERFACE
	IMPLICIT NONE
C
	INTEGER, PARAMETER :: N=200
	REAL(KIND=LDP) XD,XVECD(N)
	REAL(KIND=LDP) VALD
	REAL*4 XS,XVECS(N)
	REAL*4 VALS
	INTEGER I
C
	DO I=1,N
	  XVECD(I)=I	
	  XVECS(I)=I*2
	END DO
C
	VALD=100.0_LDP
	VALS=100.0
	WRITE(2,*)GET_INDX(VALD,XVECD,N)
	WRITE(2,*)GET_INDX(VALS,XVECS,N)
C
	STOP
	END
