	SUBROUTINE GET_LINE_ID_PG(TRANS_NAME,LINE_WAVE,EW,LINE_CENTER,FWHM_KMS)
	USE LINE_ID_MOD
	IMPLICIT NONE
!
	CHARACTER(LEN=*) TRANS_NAME
	REAL*4 EW
	REAL*4 LINE_CENTER
	REAL*4 FWHM_KMS
	REAL*4 LINE_WAVE
!
! Altered 02-Jul-2022: Now use TAU_VAL rather than 3 times, and use ABS value. 
! Altered 23-Jul-2022: Added LINE_WAVE to call.
! Altered 22-Jul-2022: Better check if ID exits.
!
	INTEGER,PARAMETER :: NMATCH=5
	REAL*4 OFF_STORE(NMATCH)
	INTEGER PNT_STORE(NMATCH)
	INTEGER PNT
!
	REAL*4 T1
	REAL*4 SHIFT
	REAL*4 TAU_VAL
!
	INTEGER I,K,LOC
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
!
	TRANS_NAME=' '; LINE_WAVE=0.0; TAU_VAL=0.01
	IF(N_LINE_IDS .EQ. 0)THEN
	  RETURN
	END IF
	IF( (LINE_CENTER-ID_WAVE(1))*(ID_WAVE(N_LINE_IDS)-LINE_CENTER) .LT. 0)THEN
	  RETURN
	END IF
!
	WRITE(6,*)'Hope',LINE_CENTER,FWHM_KMS
	OFF_STORE=1.0D+10
	PNT_STORE=0
	LOC=GET_INDX_SP(LINE_CENTER,ID_WAVE,N_LINE_IDS)
	DO I=MAX(1,LOC-10),MIN(LOC+10,N_LINE_IDS)
	  T1=2.998D+05*ABS(ID_WAVE(I)/LINE_CENTER-1.0D0)/FWHM_KMS
	  WRITE(70,'(4ES16.6,4X,A)')LINE_CENTER,FWHM_KMS,ID_WAVE(I),T1,TRIM(FULL_LINE_ID(I))
	  FLUSH(UNIT=70)
	  DO K=1,NMATCH
	    IF(T1 .LT. OFF_STORE(K))THEN
	      OFF_STORE(K+1:NMATCH)=OFF_STORE(K:NMATCH-1)
              PNT_STORE(K+1:NMATCH)=PNT_STORE(K:NMATCH-1)
	      OFF_STORE(K)=T1
              PNT_STORE(K)=I
	      EXIT
	    END IF
	  END DO
	END DO
!
! Now decide on best line. We use ABS values since values from the
! file created with PLNID can be negative.
!
	IF(PNT_STORE(1) .EQ. 0)THEN
	ELSE
	  DO K=1,NMATCH
	    WRITE(6,*)K,OFF_STORE(K)
	    IF(OFF_STORE(K) .LT. 0.5)THEN
	      I=PNT_STORE(K)
	      IF(TAU(I) .GT. ABS(TAU_VAL))THEN
	        TRANS_NAME=FULL_LINE_ID(I)
	        TAU_VAL=TAU(I)
	        LINE_WAVE=ID_WAVE(I)
	      END IF
	    END IF
	  END DO
	END IF
!
	RETURN
	END
