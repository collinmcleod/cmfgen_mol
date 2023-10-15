!
! Subroutine formulates a string (>= 132 characters) containing the
! transition name, wavelength, continuum flux and line EW.
!
	SUBROUTINE EW_FORMAT_V2(OUTSTR,TRANSITION,LAMBDA,FLUX,EW,SOBOLEV,NL,NUP)
	IMPLICIT NONE
!
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
	CHARACTER(LEN=*) TRANSITION,OUTSTR
	REAL(10) LAMBDA,FLUX,EW
	INTEGER NL,NUP
	LOGICAL SOBOLEV
!
	REAL(10) T1
	INTEGER L,TAB
	INTEGER STR_LENGTH
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
	STR_LENGTH=LEN(OUTSTR)
	TAB=1
	OUTSTR=' '
!
	IF(LAMBDA .LT. 1.0D+05)THEN
	  WRITE(OUTSTR(TAB:),'(3X,F8.2)')LAMBDA
	ELSE
	  WRITE(OUTSTR(TAB:),'(ES11.4)')LAMBDA
	END IF
!
	TAB=TAB+11
	WRITE(OUTSTR(TAB:),'(ES12.3)')FLUX
!
	T1=2.997924D+10/(1.0D-08*LAMBDA) 	!Frequency(Hz)
	T1=1.0D-23*FLUX*T1/LAMBDA               !Flux -- ergs/cm^2/sec/Hz
	T1=EW*T1
	TAB=TAB+12
	WRITE(OUTSTR(TAB:),'(2ES12.3,3X,L1,2X)')EW,T1,SOBOLEV
!
	TAB=TAB+30
	WRITE(OUTSTR(TAB:),'(2(I4,2X))')NL,NUP
!
! Idelly allow for a transition name up to 85 characters.
!
	TAB=TAB+14			!4 chracter gap
	L=LEN_TRIM(TRANSITION)
	IF(L .GT. STR_LENGTH-TAB+1)L=1+STR_LENGTH-TAB
	OUTSTR(TAB:)=TRIM(TRANSITION(1:L))
! 
	RETURN
	END
