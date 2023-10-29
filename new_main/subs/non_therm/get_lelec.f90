	SUBROUTINE GET_LELEC(LELEC,XKT,NKT,N_ELEC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NKT
	REAL(KIND=LDP) LELEC(NKT)
	REAL(KIND=LDP) XKT(NKT)
	REAL(KIND=LDP) N_ELEC
!
	real(kind=LDP), parameter :: ELECTRON_VOLT=1.60217733D-12  ! erg
	real(kind=LDP), parameter :: PLANCKS_CONSTANT=6.626075D-27 ! erg sec
	real(kind=LDP), parameter :: SPEED_OF_LIGHT=2.99792458D+10 ! cm / sec
	real(kind=LDP), parameter :: PI = 3.141592653589793238462643D0
	real(kind=LDP), parameter :: ELECTRON_MASS=9.109389D-28    !gm
	real(kind=LDP), parameter :: ELECTRON_CHARGE = 4.803206814D-10 ! esu
	real(kind=LDP), parameter :: a0 = 0.529189379D-8    ! Bohr radius in cm
!
	INTEGER I
	INTEGER IKT
	REAL(KIND=LDP) T1,T2,XV,XE_CGS
	REAL(KIND=LDP) XI_E
	REAL(KIND=LDP) PLASMA_FREQ
	REAL(KIND=LDP) GAMMA_EULER
!
! plasma_freq is in 1/s, i.e. cgs
! xi_e is put in erg !!! [e^2] has dimensions of energy x length
! we use xe_cgs to work with cgs units
!
	  plasma_freq = sqrt( 4.0D0 * PI * n_elec * ELECTRON_CHARGE**2 / ELECTRON_MASS ) ! 1/s
	  xi_e = PLANCKS_CONSTANT/2.0D0/PI * plasma_freq                                 ! erg
!
! write(6,'(A25,ES15.5)') 'Plasma Frequency [Hz]: ',plasma_freq
! write(6,'(A25,ES15.5)') 'xi electron [eV]',xi_e/ELECTRON_VOLT
	
	  gamma_euler = 0.577215665D0
	
	  do ikt=1,nkt
	     xe_cgs = xkt(ikt) * ELECTRON_VOLT
!
! Not sure if this velcoity shoud be that of thermal electrons.
! Here I take that of the non-thermal electrons we are treating at xe_cgs
!
	     xv = sqrt(2.0D0*xe_cgs/ELECTRON_MASS) ! in cm/s!!!
	
	     if (xkt(ikt).lt.14.0D0) then
	        t1 = n_elec * 2.D0 * PI * ELECTRON_CHARGE**4 / xe_cgs * &
	             log( ELECTRON_MASS * xv**3/ gamma_euler / ELECTRON_CHARGE**2 / plasma_freq  )
	        lelec(ikt) = t1
	     else
	        t2 = n_elec * 2.D0 * PI * ELECTRON_CHARGE**4 / xe_cgs * &
	             log(2.0D0 * xe_cgs / xi_e)
	        lelec(ikt) =  t2
	     endif
	  enddo
	  lelec(1:nkt) = lelec(1:nkt) / ELECTRON_VOLT ! now in eV/cm
	
	return
	end subroutine get_lelec
