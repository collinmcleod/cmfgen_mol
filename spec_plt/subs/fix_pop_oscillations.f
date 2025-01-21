!
	SUBROUTINE FIX_POP_OSCILLATIONS(POPS,R,V,SIGMA,LUSCR,ND,NT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NT
	INTEGER LUSCR
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
!
! Local variables
!
	REAL(KIND=LDP) POPS_MAT(NT,ND,5)
	REAL(KIND=LDP) R_TMP(ND)
	REAL(KIND=LDP) V_TMP(ND)
	REAL(KIND=LDP) SIGMA_TMP(ND)
	REAL(KIND=LDP) D5,D4,D3,D2,D1
!
	INTEGER I,J
	INTEGER IREC
	INTEGER IOS
	INTEGER LST_NG
	INTEGER NITSF
	INTEGER RITE_N_TIMES
	INTEGER, SAVE :: LST_OSC_ADJ=0
!
	LOGICAL NEWMOD
	LOGICAL FILE_OPEN
	LOGICAL WRITE_RVSIG
	CHARACTER(LEN=80) STRING
!
        OPEN(UNIT=LUSCR,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
          IF(IOS .EQ. 0)READ(LUSCR,'(A)',IOSTAT=IOS)STRING
          IF(IOS .EQ. 0)THEN
            IF(INDEX(STRING,'!Format date') .EQ. 0)THEN
              READ(STRING,*,IOSTAT=IOS)IREC,NITSF,RITE_N_TIMES,LST_NG
            ELSE
              READ(LUSCR,*,IOSTAT=IOS)IREC,NITSF,RITE_N_TIMES,LST_NG
            END IF
          END IF
	  INQUIRE(UNIT=LUSCR,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LUSCR)
        IF(IOS .NE. 0)RETURN
!
! Check if we have enough iterations to check on oscilatory behaviour.
!
	IF(NITSF  .LT. LST_OSC_ADJ+4)RETURN
	IF(LST_NG .LT. LST_OSC_ADJ+4)RETURN
	IF(IREC   .NE.         NITSF)RETURN
!
! We skip check if we have recently performed a hydro calculations.
!
	DO I=1,4
	  IREC=NITSF-I+1
          CALL SCR_READ_V2(R_TMP,V_TMP,SIGMA_TMP,POPS_MAT(1,1,I),IREC,
	1                 NITSF,RITE_N_TIMES,LST_NG,
	1                 WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  DO J=1,ND
	    IF(R_TMP(J) .NE. R(J))RETURN
	  END DO
	END DO
	POPS_MAT(:,:,5)=POPS
!
	DO J=1,ND
	  DO I=1,NT
	    D5=POPS_MAT(I,J,5)-POPS_MAT(I,J,4)
	    D4=POPS_MAT(I,J,4)-POPS_MAT(I,J,3)
	    D3=POPS_MAT(I,J,3)-POPS_MAT(I,J,2)
	    D2=POPS_MAT(I,J,2)-POPS_MAT(I,J,1)
	    IF(D5*D4 .LT. 0.0_LDP .AND. D4*D3 .LT. 0.0_LDP .AND. D3*D2 .LT. 0.0_LDP)THEN
	      IF(ABS(D5) .GT. ABS(D3) .AND. ABS(D4) .GT. ABS(D2))THEN
	        POPS(I,J)=(POPS_MAT(I,J,2)+POPS_MAT(I,J,1))/2.0_LDP
	      ELSE
	        POPS(I,J)=(POPS_MAT(I,J,5)+POPS_MAT(I,J,4))/2.0_LDP
	      END IF
	      LST_OSC_ADJ=NITSF+1
	    END IF
	  END DO
	END DO
!
	RETURN
	END
