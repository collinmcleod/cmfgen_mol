	PROGRAM TST_MON
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: ND=60
	INTEGER, PARAMETER :: NQ=5*(ND-1)+1
	INTEGER, PARAMETER :: LIN_END=1
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) CHI(ND)
!
	REAL(KIND=LDP) QZR(NQ)
	REAL(KIND=LDP) NEW_CHI(NQ)
!
	INTEGER I,J
	REAL(KIND=LDP) DELR
!	
        DO I=1, ND ; CHI(I)= ( (I+1.0_LDP)**2 )/ND; END DO
        DO I=1,ND
          R(I)=1.0_LDP+0.1_LDP*(I**2)
        END DO
C
	J=1
        DO I=1,ND-1
	  DELR=R(I+1)-R(I)
          QZR(J)=R(I)
          QZR(J+1)=R(I)+0.04_LDP*DELR
          QZR(J+2)=R(I)+0.16_LDP*DELR
          QZR(J+3)=R(I)+0.40_LDP*DELR
          QZR(J+4)=R(I)+0.70_LDP*DELR
	  J=J+5
        END DO
	QZR(NQ)=R(ND)
!
	DO I=1,ND
	  WRITE(6,*)R(I),CHI(I)
	END DO
	DO I=1,NQ
	  WRITE(6,*)QZR(I)
	END DO
!
	CALL TUNE(1,'MON')
	 DO J=1,10000
	   CALL MON_INTERP_FAST(NEW_CHI,NQ,LIN_END,QZR,NQ,CHI,ND,R,ND)
	 END DO
	CALL TUNE(2,'MON')
	CALL TUNE(3,' ')
!
	STOP
	END	
