!
! These two functions allow the logical units associated with general
! input/output file (terminal) to be specified. Virtually all programs
! currently use units 5 and 6.
!
! Altered: 20-Jan-2014 - WARN_LU introduced.
!
! The LU for iterminal output and error messages are currently the same.
!
	FUNCTION TERM_IN()
	USE SET_KIND_MODULE
	INTEGER TERM_IN
	TERM_IN=5
	RETURN
	END
!
	FUNCTION TERM_OUT()
	USE SET_KIND_MODULE
	INTEGER TERM_OUT
	TERM_OUT=6
	RETURN
	END
!
	FUNCTION ERROR_LU()
	USE SET_KIND_MODULE
	INTEGER ERROR_LU
	ERROR_LU=6
	RETURN
	END
!
	FUNCTION WARNING_LU()
	USE SET_KIND_MODULE
	INTEGER WARNING_LU
	WARNING_LU=4
	RETURN
	END
