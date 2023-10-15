!
! Subroutine formulates a string (>= 132 characters) containing the
! transition name, wavelength, continuum flux and line EW, and outputs
! tological unit LU_EW.
!
	SUBROUTINE WRITE_OUT_EW(TRANSITION,LAMBDA_VAC,FLUX,EW,EMISS_EW,SOBOLEV,NL,NUP,LU_EW)
	IMPLICIT NONE
!
! Altered 01-Sep-2018 : Header output by this subroutine.
!                         Lambda should be passed as a vacuum wavelngth
!                         (note the change to earlier versions).
! Altered 01-Mar-2019 : Added EW_EMISS to call, changed to V3
! Altered 17-MAr-2018 : Changed to V2
!                       Line flux and level ID now output.
!                       String length is not checked.
! Altered 01-Aug-1997 : OUTSTR increased to 132 characters.
!                         Transition name now output at end of string, so
!                         as species with long names will be successfully
!                         identified.
! Altered 24-May-1996 : ERROR_LU, LUER inserted
!                       CNT check inserted.
!
! Created 16-Oct-1989
!
	CHARACTER(LEN=*) TRANSITION
	REAL(10) LAMBDA_VAC
	REAL(10) FLUX
	REAL(10) EW
	REAL(10) EMISS_EW
	INTEGER NL,NUP
	INTEGER LU_EW
	LOGICAL SOBOLEV
!
	REAL(10) T1
	INTEGER L,TAB
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL, SAVE :: FIRST=.TRUE.
	REAL(10) LAMVACAIR,LAMBDA_AIR,SPEED_OF_LIGHT
	EXTERNAL LAM_VAC,SPEED_OF_LIGHT
	CHARACTER(LEN=132) OUTSTR
!
	LUER=ERROR_LU()
	TAB=1
	OUTSTR=' '
!
	IF(FIRST)THEN
	  WRITE(LU_EW,'(A)')'!'
	  WRITE(LU_EW,'(A)')'! EW is the classical EW except is positive when the line is in emission'
	  WRITE(LU_EW,'(A)')'! AEW is is simliar to the classical EW except t use the ABS(Fv-Fc) for the'
	  WRITE(LU_EW,'(A)')'!    integral. As a consequence a P Cygni profile will have non-zero EW even'
	  WRITE(LU_EW,'(A)')'!    when the emission and absorption components have equal strength.'
	  WRITE(LU_EW,'(A)')'! Alam give the air wavelength for Lambda > 2000A'
	  WRITE(LU_EW,'(A)')'!'
	  WRITE(LU_EW,'(A)')'!Formate date: 01-Sep-2019'
	  WRITE(LU_EW,'(A)')'!'
	  WRITE(LU_EW,'(3X,A,2(2X,A),5X,A,5X,A,4X,A,3(2X,A),4X,A)')'Lam(Ang)','Alam(Ang)','C.Flux(Jy)',
	1                  'EW(Ang)','AEW(Ang)','L.Flux','Sob','  NL',' NUP','Trans. Name'
	  FIRST=.FALSE.
	END IF
!
	IF(LAMBDA_VAC .LT. 1.0D+05)THEN
	  WRITE(OUTSTR(TAB:),'(3X,F8.2)')LAMBDA_VAC
	ELSE
	  WRITE(OUTSTR(TAB:),'(ES11.4)')LAMBDA_VAC
	END IF
!
	T1=1.0D-07*SPEED_OF_LIGHT()/LAMBDA_VAC		!Frequency in units of 10^15 Hz
	LAMBDA_AIR=LAMVACAIR(T1)
	TAB=TAB+11
	IF(LAMBDA_AIR .LT. 1.0D+05)THEN
	  WRITE(OUTSTR(TAB:),'(3X,F8.2)')LAMBDA_AIR
	ELSE
	  WRITE(OUTSTR(TAB:),'(ES11.4)')LAMBDA_AIR
	END IF

	TAB=TAB+11
	WRITE(OUTSTR(TAB:),'(ES12.3)')FLUX
!
	T1=2.997924D+10/(1.0D-08*LAMBDA_VAC) 	!Frequency(Hz)
	T1=1.0D-23*FLUX*T1/LAMBDA_VAC               !Flux -- ergs/cm^2/sec/Hz
	T1=EW*T1
	TAB=TAB+12
	WRITE(OUTSTR(TAB:),'(3ES12.3,3X,L1,2X)')EW,EMISS_EW,T1,SOBOLEV
!
	TAB=TAB+42
	WRITE(OUTSTR(TAB:),'(2(I4,2X))')NL,NUP
        WRITE(LU_EW,'(A)')TRIM(OUTSTR)//'   '//TRIM(TRANSITION)
!
	RETURN
	END
