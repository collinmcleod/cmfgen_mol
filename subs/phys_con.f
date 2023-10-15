C
C File contains physical and astrophysical constants of general interest.
C Each constant can be obtained by a FUNCTION call, and all units are CGS
C unless otherwise is specified in the function name.
C
C Altered 07-Jan-2008: Rsun revised very slightly, TEFF_SUN inserted
C Altered 29-Oct-2007: Stefan Boltzman constant added.
C Altered 29-Apr-1991: ALL routines. REAL(10) now specified in separate
C                      statement for CRAY compatibility.
C Altered 8_Oct-1990 : RYDBERG_CARBON (and NITROGEN) wre returning zero.
C                      IMPLCIT NONE placed in all functions.
C
	FUNCTION SPEED_OF_LIGHT()
	IMPLICIT NONE
	REAL(10) SPEED_OF_LIGHT
	SPEED_OF_LIGHT=2.99792458D+10 			!Exact cm/s
	RETURN
	END

	FUNCTION ELECTRON_VOLT()
	IMPLICIT NONE
	REAL(10) ELECTRON_VOLT
	ELECTRON_VOLT=1.60217733D-12                    !ergs
	RETURN
	END

	FUNCTION GRAVITATIONAL_CONSTANT()
	IMPLICIT NONE
	REAL(10) GRAVITATIONAL_CONSTANT
	GRAVITATIONAL_CONSTANT=6.67259D-08		!cm^3/gm/s
	RETURN
	END

	FUNCTION PLANCKS_CONSTANT()
	IMPLICIT NONE
	REAL(10) PLANCKS_CONSTANT
	PLANCKS_CONSTANT=6.626075D-27			!erg sec
	RETURN
	END

	FUNCTION ATOMIC_MASS_UNIT()
	IMPLICIT NONE
	REAL(10) ATOMIC_MASS_UNIT
	ATOMIC_MASS_UNIT=1.660540D-24			!gm
	RETURN
	END

	FUNCTION ELECTRON_MASS()
	IMPLICIT NONE
	REAL(10) ELECTRON_MASS
	ELECTRON_MASS=9.109389D-28			!gm
	RETURN
	END

	FUNCTION BOLTZMANN_CONSTANT()
	IMPLICIT NONE
	REAL(10) BOLTZMANN_CONSTANT
	BOLTZMANN_CONSTANT=1.380658D-16			!erg/K
	RETURN
	END

	FUNCTION STEFAN_BOLTZ()
	IMPLICIT NONE
	REAL(10) STEFAN_BOLTZ
	STEFAN_BOLTZ=5.670400D-05			!ergs/cm^2/K^4
	RETURN
	END

C 
	FUNCTION RYDBERG_INF()
	IMPLICIT NONE
	REAL(10) RYDBERG_INF
	RYDBERG_INF=109737.31534D0			!/cm
	RETURN
	END

	FUNCTION RYDBERG_HYDROGEN()
	IMPLICIT NONE
	REAL(10) RYDBERG_HYDROGEN
	RYDBERG_HYDROGEN=109677.6D0			!/cm
	RETURN
	END

	FUNCTION RYDBERG_HELIUM()
	IMPLICIT NONE
	REAL(10) RYDBERG_HELIUM
	RYDBERG_HELIUM=109722.3D0			!/cm
	RETURN
	END

	FUNCTION RYDBERG_CARBON()
	IMPLICIT NONE
	REAL(10) RYDBERG_CARBON
	RYDBERG_CARBON=109732.3D0			!/cm
	RETURN
	END

	FUNCTION RYDBERG_NITROGEN()
	IMPLICIT NONE
	REAL(10) RYDBERG_NITROGEN
	RYDBERG_NITROGEN=109733.0D0			!/cm
	RETURN
	END

	FUNCTION RYDBERG_OXYGEN()
	IMPLICIT NONE
	REAL(10) RYDBERG_OXYGEN
	RYDBERG_OXYGEN=109733.5D0			!/cm
	RETURN
	END
C 
	FUNCTION MASS_SUN()
	IMPLICIT NONE
	REAL(10) MASS_SUN
	MASS_SUN=1.989D+33  				!gm
	RETURN
	END

	FUNCTION RAD_SUN()
	IMPLICIT NONE
	REAL(10) RAD_SUN
!	RAD_SUN=6.96D+10		!cm
	RAD_SUN=6.9599D+10		!cm (changed 7-Jan-2008)
	RETURN
	END

	FUNCTION LUM_SUN()
	IMPLICIT NONE
	REAL(10) LUM_SUN
	LUM_SUN=3.826D+33				!erg/s
	RETURN
	END

	FUNCTION TEFF_SUN()
	IMPLICIT NONE
	REAL(10) TEFF_SUN
	TEFF_SUN=5770.0D0				!K
	RETURN
	END

	FUNCTION PARSEC()
	IMPLICIT NONE
	REAL(10) PARSEC
	PARSEC=3.0856D+18       			!cm
	RETURN
	END

	FUNCTION ASTRONOMICAL_UNIT()
	IMPLICIT NONE
	REAL(10) ASTRONOMICAL_UNIT
	ASTRONOMICAL_UNIT=1.49597892D+13		!cm
	RETURN
	END

	FUNCTION JANSKY()
	IMPLICIT NONE
	REAL(10) JANSKY
	JANSKY=1.0D-23					!ergs/cm^2/s/Hz
	RETURN
	END

	FUNCTION PI()			!Use FUN_PI instead so PI
	IMPLICIT NONE			!can be variable
	REAL(10) PI			!Use FUN_PI instead so PI
	PI=3.141592653589793238462643D0
	RETURN
	END
C
	FUNCTION FUN_PI()
	IMPLICIT NONE
	REAL(10) FUN_PI
	FUN_PI=3.141592653589793238462643D0
	RETURN
	END

	FUNCTION SECS_IN_YEAR()
	IMPLICIT NONE
	REAL(10) SECS_IN_YEAR
	SECS_IN_YEAR=31557600.0D0     !365.25*24*60*60
	RETURN
	END
