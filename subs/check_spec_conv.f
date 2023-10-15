!
! Routine to compare changes in successive spectra comuted by CMFGEN. It should not
! be used with the last spectrum as an alternative computation method is used which
! have improved accuracy and less bleeding.
!
! The rountne is still under development.
!
	SUBROUTINE CHECK_SPEC_CONV(OBS_FLUX,OBS_FREQ,LAM_ST,LAM_END,SPECTRUM_CONVERGED,NOBS)
	IMPLICIT NONE
!
! Altered 06-Mar-2023: Fixed error that occured when no change in spectrum at a particular freqency.
! Altered 21-Aug-2020: Changed output format; Dimension LAMS from 0.
! Created 10-Feb-2020
!
	INTEGER NOBS
	REAL(10) OBS_FLUX(NOBS)
	REAL(10) OBS_FREQ(NOBS)
	REAL(10) LAM_ST			!not yet used
	REAL(10) LAM_END
	LOGICAL SPECTRUM_CONVERGED
!
	INTEGER, PARAMETER :: NLAM=8
	INTEGER, PARAMETER :: NHIST=8
	INTEGER, PARAMETER :: NITS=5
	INTEGER, SAVE :: HIST(NHIST,NLAM,NITS)
	REAL(10), SAVE :: LAM_LIMS(0:NLAM)
	REAL(10), SAVE :: HIST_LIMS(NHIST)
	REAL(10), SAVE :: MAX_FRAC_CHNG(NITS)
!
	DATA LAM_LIMS/0,10.0D0,100.0D0,228.0D0,911.0D0,3600.0D0,10000.0D0,1.0D+05,1.0D+10/
	DATA HIST_LIMS/1.0D-04,1.0D-03,1.0D-02,1.0D-01,1.0,10.0D0,100.0D0,1000.0D0/
!
	INTEGER, SAVE :: NIT_CUR 
	REAL(10), ALLOCATABLE, SAVE :: OBS_SAVE(:)
	REAL(10), ALLOCATABLE, SAVE :: WRK_VEC(:)
!
	REAL(10) C_KMS
	REAL(10) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	REAL(10) T1
	REAL(10) LAM
	INTEGER I
	INTEGER IL
	INTEGER IH
	INTEGER IT
!
! On the first call we simply allocate the required storage, and return.
!
	SPECTRUM_CONVERGED=.FALSE.
	IF(.NOT. ALLOCATED(OBS_SAVE))THEN
	  ALLOCATE(OBS_SAVE(NOBS))
	  ALLOCATE(WRK_VEC(NOBS))
	  HIST=0
	  NIT_CUR=0
	  MAX_FRAC_CHNG=0.0D0
	  OBS_SAVE=OBS_FLUX
	  RETURN
	END IF
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
	NIT_CUR=MIN(NIT_CUR+1,NITS)
	DO IT=NIT_CUR-1,1,-1
	  HIST(:,:,IT+1)=HIST(:,:,IT)
	  MAX_FRAC_CHNG(IT+1)=MAX_FRAC_CHNG(IT)	
	END DO
	HIST(:,:,1)=0
	MAX_FRAC_CHNG(1)=0.0D0
!
! The factor of 10^7 arrises as follows:
!   100 to get in percentage terms.
!   10^5 to get the corect index when we take the log.
!	
	IL=1; IT=1
	DO I=1,NOBS
	  T1=OBS_SAVE(I)+OBS_FLUX(I)
	  IF(T1 .NE. 0.0D0)THEN
	    WRK_VEC(I)=ABS( (OBS_SAVE(I)-OBS_FLUX(I))/T1 )
	    MAX_FRAC_CHNG(IT)=MAX(MAX_FRAC_CHNG(IT),WRK_VEC(I))
	    WRK_VEC(I)=1.0D+07*WRK_VEC(I)
	  ELSE
	    WRK_VEC(I)=1.0D-20
	  END IF
	  IF(WRK_VEC(I) .LT. 1.0D-20)WRK_VEC(I)=1.0D-20
	  LAM=0.01D0*C_KMS/OBS_FREQ(I)
	  IF(LAM .GT. LAM_LIMS(IL))IL=IL+1
	  IH=MAX(1,INT(LOG10(WRK_VEC(I))))
	  IH=MIN(IH,NHIST)
	  HIST(IH,IL,IT)=HIST(IH,IL,IT)+1 
	END DO
!
	WRITE(6,'(A)')' '
	WRITE(6,*)'Spectrum convergence information -- % changes'
	WRITE(6,'(17X,8(2X,ES7.1))')(HIST_LIMS(IH),IH=1,NHIST)
	DO IL=1,NLAM
	  IF(IL .LE. 6)THEN
	    WRITE(6,'(2X,F7.1,A,F7.1,8(2X,I7))')LAM_LIMS(IL-1),'-',LAM_LIMS(IL),(HIST(IH,IL,IT),IH=1,NHIST)
	  ELSE
	    WRITE(6,'(2X,ES7.1,A,ES7.1,8(2X,I7))')LAM_LIMS(IL-1),'-',LAM_LIMS(IL),(HIST(IH,IL,IT),IH=1,NHIST)
	  END IF
	END DO
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Largest % change in spectrum for last 5 full iterations - last listed first'
	WRITE(6,'(10ES10.2)')100.0D0*MAX_FRAC_CHNG
	WRITE(6,'(A)')' '
!
	OBS_SAVE=OBS_FLUX
	IF( MAXVAL(MAX_FRAC_CHNG) .LT. 0.002)SPECTRUM_CONVERGED=.TRUE.
	WRITE(6,*)'Exiting CHECK_SPECT'
!
	RETURN
	END
