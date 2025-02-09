	MODULE CONSTANTS_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Module contains physical and astrophysical constants of general interest.
! All constants are cgs. Designed to replace function calls in phys_con.f.
!
! Altered 15-Sep-2023 :" Updated values of several constants.
! Created  7-Jan-2008 :: Based on cur_cmf/subs/phys_con.f
!
	REAL(KIND=LDP), PARAMETER :: SPEED_OF_LIGHT         = 2.99792458E+10_LDP      !Exact cm/s
	REAL(KIND=LDP), PARAMETER :: ELECTRON_VOLT          = 1.60217733E-12_LDP      !ergs
	REAL(KIND=LDP), PARAMETER :: GRAVITATIONAL_CONSTANT = 6.67259E-08_LDP         !cm^3/gm/s
	REAL(KIND=LDP), PARAMETER :: PLANCKS_CONSTANT       = 6.62607015E-27_LDP      !erg sec
	REAL(KIND=LDP), PARAMETER :: ATOMIC_MASS_UNIT       = 1.6605390666E-24_LDP    !gm
	REAL(KIND=LDP), PARAMETER :: ELECTRON_MASS          = 9.1093837015E-28_LDP    !gm
	REAL(KIND=LDP), PARAMETER :: ELECTRON_CHARGE        = 4.80320427E-10_LDP      !esu
	REAL(KIND=LDP), PARAMETER :: BOLTZMANN_CONSTANT     = 1.380649E-16_LDP        !erg/K
	REAL(KIND=LDP), PARAMETER :: STEFAN_BOLTZ           = 5.670374419E-05_LDP     !ergs/cm^2/K^4
!
	REAL(KIND=LDP), PARAMETER :: RYDBERG_INF     =109737.31534_LDP		!/cm
	REAL(KIND=LDP), PARAMETER :: RYDBERG_HYDROGEN=109677.6_LDP		!/cm
	REAL(KIND=LDP), PARAMETER :: RYDBERG_HELIUM  =109722.3_LDP		!/cm
	REAL(KIND=LDP), PARAMETER :: RYDBERG_CARBON  =109732.3_LDP		!/cm
	REAL(KIND=LDP), PARAMETER :: RYDBERG_NITROGEN=109733.0_LDP		!/cm
	REAL(KIND=LDP), PARAMETER :: RYDBERG_OXYGEN  =109733.5_LDP		!/cm
!
	REAL(KIND=LDP), PARAMETER :: MASS_SUN=1.989E+33_LDP  			!gm
	REAL(KIND=LDP), PARAMETER :: RAD_SUN =6.9599E+10_LDP			!cm (changed 7-Jan-2008)
	REAL(KIND=LDP), PARAMETER :: LUM_SUN =3.826E+33_LDP				!erg/s
	REAL(KIND=LDP), PARAMETER :: TEFF_SUN=5770.0_LDP				!K
!
	REAL(KIND=LDP), PARAMETER :: ASTRONOMICAL_UNIT=1.49597892E+13_LDP		!cm
	REAL(KIND=LDP), PARAMETER :: PARSEC           =3.0856E+18_LDP       	!cm
	REAL(KIND=LDP), PARAMETER :: JANSKY           =1.0E-23_LDP			!ergs/cm^2/s/Hz
!
	REAL(KIND=LDP), PARAMETER ::           PI=3.141592653589793238462643_LDP
	REAL(KIND=LDP), PARAMETER ::       FUN_PI=3.141592653589793238462643_LDP
	REAL(KIND=LDP), PARAMETER :: SECS_IN_YEAR=31557600.0_LDP     		!365.25*24*60*60
!
	END MODULE CONSTANTS_MOD	
