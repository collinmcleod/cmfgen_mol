!
! Routine to place X-ray bound-free edge frequencies into a vector.
!
	SUBROUTINE SET_X_FREQ_V2(FREQ,NCF,NCF_MAX,
	1                   MAX_CONT_FREQ,ZCORE,ZION,
	1                   NI_PRES,N2_PRES)
	USE XRAY_DATA_MOD
	IMPLICIT NONE
!
! Created  23-Oct-2000: Bsed on SET_X_FRQ
!
	INTEGER*4 NCF,NCF_MAX
	REAL*8 MAX_CONT_FREQ		!Units 10^15 Hz
	REAL*8 FREQ(NCF_MAX)
	REAL*8 ZCORE,ZION
	LOGICAL NI_PRES,N2_PRES 
!
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	REAL*8 CONV_FAC
	INTEGER*4 I,J,LOOP_PQN_MAX
	INTEGER*4 IZ,NE
!
	IF( .NOT. NI_PRES)RETURN
	IF( .NOT. N2_PRES)RETURN
!
	CONV_FAC=0.24191				!ev to 10^15 Hz
	IZ=ZCORE
	NE=ZCORE-ZION+1
!
	IF(IZ .GT. 30)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in XCROSS_V2'
	  WRITE(LUER,*)'No X-ray data available for species'//
	1                     ' with Atomic No.  > 30'
	  STOP
	END IF
!
! We only wish to return photoionization edges for the core
! (i.e. non-valence) electrons. This section of the code should
! be consistent with that in XCROSS_V2,
!
	DO J=0,ANG_MAX_X
	  LOOP_PQN_MAX=3
	  IF(J .GT. 1 .AND. NE .LT. 20)LOOP_PQN_MAX=2
	  IF(NE .LT. 18)LOOP_PQN_MAX=2
	  IF(NE .LE. 10)LOOP_PQN_MAX=1
	  IF(NE .LE. 2)RETURN
	  DO I=1,LOOP_PQN_MAX
	    IF( SIG_0_X(IZ,NE,I,J) .NE. 0)THEN
              IF(NCF+1 .GT. NCF_MAX)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'Error in SET_X_FREQ --- NCF_MAX too small'
	        STOP
	      END IF
	      NCF=NCF+1
	      FREQ(NCF)=E_THRESH_X(IZ,NE,I,J)*CONV_FAC
              IF(FREQ(NCF) .GT. MAX_CONT_FREQ)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'**** Warning -- Warning -- Warning ****'
	        WRITE(LUER,*)'Max. cont. freq may be too small in in SET_X_FREQ'
	        WRITE(LUER,*)'IZ=',IZ,'NE=',NE,'PQN=',I,'ANG=',J
                NCF=NCF-1
	      END IF
	    END IF
	  END DO
	END DO
!
	RETURN
	END
