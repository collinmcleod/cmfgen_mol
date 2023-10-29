	SUBROUTINE WR_CRIT(OMEGA,TEMP,HDKT,AXzV,EDGEXzV,GXzV,XzVLEVNAME,NXzV,DESC,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NXzV
	INTEGER LU
	REAL(KIND=LDP) TEMP
	REAL(KIND=LDP) HDKT
!
	REAL(KIND=LDP) AXzV(NXzV,NXzV)
	REAL(KIND=LDP) OMEGA(NXzV,NXzV)
	REAL(KIND=LDP) EDGEXzV(NXzV)
	REAL(KIND=LDP) GXzV(NXzV)
!
	CHARACTER(LEN=*) DESC
	CHARACTER(LEN=*) XzVLEVNAME(NXzV)
!
	CHARACTER*80 FORM
	REAL(KIND=LDP) ASUM(NXzV)
	REAL(KIND=LDP) COL_SUM(NXzV)
	REAL(KIND=LDP) T1
	INTEGER I,J,IOS
	INTEGER, PARAMETER :: IZERO=0
!
        FORM=TRIM(DESC)//'_CRIT'
        CALL GEN_ASCI_OPEN(LU,FORM,'UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(6, '(A20,ES14.3)')XzVLEVNAME(1)
	WRITE(LU,'(A20,ES14.3)')XzVLEVNAME(1)
!
	ASUM(:)=0.0D0; COL_SUM(:)=0.0D0
	DO I=2,NXzV
	  DO J=I+1,NXzV
	    T1=HDKT*(EDGEXzV(I)-EDGEXzV(J))/TEMP
	    COL_SUM(I)=COL_SUM(I)+8.63D-08*OMEGA(I,J)/GXzV(I)*EXP(-T1)
	  END DO
	  DO J=1,I-1
	    COL_SUM(I)=COL_SUM(I)+8.63D-08*OMEGA(J,I)/GXzV(I)
	    ASUM(I)=ASUM(I)+AXzV(I,J)
	  END DO
	  COL_SUM(I)=COL_SUM(I)/SQRT(TEMP)
	  IF(COL_SUM(I) .NE. 0.0D0)THEN
	    WRITE(6, '(A20,ES14.3)')XzVLEVNAME(I),ASUM(I)/COL_SUM(I)
	    WRITE(LU,'(A20,ES14.3)')XzVLEVNAME(I),ASUM(I)/COL_SUM(I)
	  ELSE
	    WRITE(6, '(A20,ES14.3)')XzVLEVNAME(I)
	    WRITE(LU,'(A20,ES14.3)')XzVLEVNAME(I)
	  END IF
	END DO
	CLOSE(LU)
!
	RETURN
	END
