	MODULE MOD_EW_VARIABLES
	REAL*4 XST,XEND                 !Line limits in world coordinates
	REAL*4 XMEAN
	REAL*4 YCONT
        REAL*4 EW  			!
	REAL*4 EWL,EWH			!Line equivalent width with a % shift in continuum
	REAL*4 SIGMA			!Standard deviation
	REAL*4 SKEWNESS
	REAL*4 KURTOSIS
	LOGICAL :: USE_MILLI_ANG=.TRUE.
	CHARACTER(LEN=80) TRANS_NAME
	SAVE
	END MODULE MOD_EW_VARIABLES
!
	SUBROUTINE WRITE_EW_HEADER(LUOUT)
	USE MOD_EW_VARIABLES
	IMPLICIT NONE
	INTEGER LUOUT
!
	IF(USE_MILLI_ANG)THEN
	  WRITE(LUOUT,'(A,5(3X,A),20X,3A,20X,A)')
	1           '!','    XST',' XEND','   Line Loc','    F(cont)',
	1               '  EW(mA)','FWHM(km/s)',' Sig(km/s)','    Sig(A)','Plot'
	ELSE
	  WRITE(LUOUT,'(A,(6X,A),20X,4A)')'!','    XST',' XEND','Line Loc',' F(cont)',
	1                  '   EW(A)','FWHM(km/s)',' Sig(km/s)',' Sig(A)',' Plot'
	END IF
!
	RETURN
	END

	SUBROUTINE WR_EW_VARIABLES(IP,LUOUT)
	USE MOD_EW_VARIABLES
	IMPLICIT NONE
	INTEGER IP,LUOUT
	REAL*4 T1
!
! The FWHM assume the line is Gaussian in shape.
!
	T1=2.998D+05*SIGMA/XMEAN
	IF(USE_MILLI_ANG)THEN
	  WRITE(LUOUT,'(3F11.3,ES14.4,F11.2,4F10.2,3F10.3,3X,I3,3X,A)')
	1                XST,XEND,XMEAN,YCONT,1000.0*EW,
	1                100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS,IP,TRIM(TRANS_NAME)
	ELSE
	  WRITE(LUOUT,'(3F14.3,2ES14.4,4F10.2,3F10.3,3X,I3,3X,A)')
	1                XST,XEND,XMEAN,YCONT,EW,
	1                100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS,IP,TRIM(TRANS_NAME)
	END IF
	FLUSH(LUOUT)
!
	RETURN
	END
