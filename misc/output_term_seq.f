	SUBROUTINE OUTPUT_TERM_SEQ(LEVEL_NAME,STAT_WGT,FEDGE,NLEV)
	IMPLICIT NONE
!
	INTEGER NLEV
	REAL*8 STAT_WGT(NLEV)
	REAL*8 FEDGE(NLEV)
	CHARACTER(LEN=*) LEVEL_NAME(NLEV)
!
	INTEGER I,J,K,L
	INTEGER, PARAMETER :: NANG=9
	INTEGER, PARAMETER :: NSPIN=8
	CHARACTER(LEN=2) ANGE(NANG)
	CHARACTER(LEN=2) ANGO(NANG)
	CHARACTER(LEN=1) SPIN(NSPIN)
	DATA ANGE/'Se','Pe','De','Fe','Ge','Ie','He','Ke','Le'/
	DATA ANGO/'So','Po','Do','Fo','Go','Io','Ho','Ko','Lo'/
	DATA SPIN/'1','2','3','4','5','6','7','8'/
!
	DO I=1,NANG
	  DO L=1,NSPIN
	    WRITE(20,'(A)')' '
	    WRITE(20,'(A)')SPIN(L)//ANGE(I)//' sequence'
	    DO J=1,NLEV
	      K=LEN_TRIM(LEVEL_NAME(J))
	      IF(LEVEL_NAME(J)(K-2:K) .EQ. SPIN(L)//ANGE(I))THEN
	       WRITE(20,'(A20,4X,F5.1,F15.5)')TRIM(LEVEL_NAME(J)),STAT_WGT(J),FEDGE(J)
	     END IF
	    END DO
	  END DO
	END DO
!
	DO I=1,NANG
	  DO L=1,NSPIN
	    WRITE(20,'(A)')' '
	    WRITE(20,'(A)')SPIN(L)//ANGO(I)//' sequence'
	    DO J=1,NLEV
	      K=LEN_TRIM(LEVEL_NAME(J))
	      IF(LEVEL_NAME(J)(K-2:K) .EQ. SPIN(L)//ANGO(I))THEN
	         WRITE(20,'(A20,4X,F5.1,F15.5)')TRIM(LEVEL_NAME(J)),STAT_WGT(J),FEDGE(J)
	      END IF
	    END DO
	  END DO
	END DO
!
	RETURN
	END

