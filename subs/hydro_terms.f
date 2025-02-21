C
C Subroutine to compare terms in the Hydronamical equations. Routine is
C meant for illustration purposes only.
C
	SUBROUTINE HYDRO_TERMS(POP_ATOM,R,V,T,SIGMA,ED,LUM_STAR,
	1              STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1              FLUX_MEAN_OPAC,ELEC_MEAN_OPAC,LU_OUT,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 17-Mar-2001 : Altered V output format
C Altered 27-Oct-2000 : Percentage error improved.
C                       STARS mass now output.
C Altered 24-May-1996 : L GEN_ASCI_OPEN and ERROR_LU installed.
C Altered 06-Mar-1996 : Description of terms added to output.
C Created 28-Mar-1994
C
	INTEGER ND,LU_OUT
C
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) LUM_STAR(ND)
C
	REAL(KIND=LDP) FLUX_MEAN_OPAC(ND)
	REAL(KIND=LDP) ELEC_MEAN_OPAC(ND)
C
	REAL(KIND=LDP) STARS_MASS
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT
C
C External functions
C
	REAL(KIND=LDP) BOLTZMANN_CONSTANT
	REAL(KIND=LDP) GRAVITATIONAL_CONSTANT
	REAL(KIND=LDP) MASS_SUN
	REAL(KIND=LDP) LUM_SUN
	REAL(KIND=LDP) FUN_PI
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) ATOMIC_MASS_UNIT
C
	EXTERNAL BOLTZMANN_CONSTANT,GRAVITATIONAL_CONSTANT,MASS_SUN
	EXTERNAL LUM_SUN,FUN_PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
C
C Local variables.
C
	REAL(KIND=LDP) VdVdR
	REAL(KIND=LDP) dPdR_ON_ROH
	REAL(KIND=LDP) g_TOT
	REAL(KIND=LDP) g_ELEC
	REAL(KIND=LDP) g_GRAV
	REAL(KIND=LDP) g_RAD
	REAL(KIND=LDP) ERROR
C
	REAL(KIND=LDP) GRAV_CON
	REAL(KIND=LDP) dP_CON
	REAL(KIND=LDP) RAD_CON
C
	CHARACTER*80 FMT
	INTEGER I,IOS,ERROR_LU
	EXTERNAL ERROR_LU
C
	I=0
	CALL GEN_ASCI_OPEN(LU_OUT,'HYDRO','UNKNOWN',' ',' ',I,IOS)
	IF(IOS .NE. 0)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error opening HYDRO in HYDRO_TERMS'
	  WRITE(I,*)'IOS=',IOS
	END IF
C
C Output file header.
C
	WRITE(LU_OUT,'(1X,6X,A,6X, 8X,A,3X, 2X,A, 5(5X,A,1X), 6X,A)')
	1 'R','V','% Error','   VdVdR',
	1                   'dPdR/ROH',
	1                   '   g_TOT',
	1                   '   g_RAD',
	1                    '  g_ELEC','Gamma'
C
C 1.0D-06 = [ 10^(-10) {from 1/r) * 10^4 {from T}
	dP_CON=1.0E-06_LDP*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/
	1               ATOMIC_MASS_UNIT()
C
	GRAV_CON=1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()
C
C       1.0E-10 (CHI) / 1.0E-20 (1/r^2)
	RAD_CON=1.0E-30_LDP*LUM_SUN()/ATOMIC_MASS_UNIT()/4.0_LDP/FUN_PI()/
	1                  SPEED_OF_LIGHT()/MEAN_ATOMIC_WEIGHT
C
	DO I=1,ND
	  VdVdR=V(I)*V(I)*(SIGMA(I)+1.0_LDP)/R(I)
	  IF(I .EQ. 1)THEN
	    dPdR_ON_ROH=dP_CON*( (POP_ATOM(1)+ED(1))*T(1)  -
	1                     (POP_ATOM(2)+ED(2))*T(2) ) /
	1                     (R(1)-R(2))/POP_ATOM(1)
	  ELSE IF(I .EQ. ND)THEN
	   dPdR_ON_ROH=dP_CON*( (POP_ATOM(ND-1)+ED(ND-1))*T(ND-1)  -
	1                     (POP_ATOM(ND)+ED(ND))*T(ND) ) /
	1                     (R(ND-1)-R(ND))/POP_ATOM(ND)
	  ELSE
	    dPdR_ON_ROH=dP_CON*( (POP_ATOM(I-1)+ED(I-1))*T(I-1)  -
	1                    (POP_ATOM(I+1)+ED(I+1))*T(I+1) ) /
	1                    (R(I-1)-R(I+1))/POP_ATOM(I)
	  END IF
C
C Gravitational force per unit mass (in cgs units).
C
	  g_grav=GRAV_CON/R(I)/R(I)
C
C Radiation pressure forces
C
	  g_rad=RAD_CON*LUM_STAR(I)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  g_elec=RAD_CON*LUM_STAR(I)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
C
	  g_TOT= g_RAD-g_GRAV
	  ERROR=200.0_LDP*(VdVdR+dPdR_ON_ROH-g_TOT)/
	1        ( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )
C
	  IF(R(I) .GT. 9.99E+04_LDP)THEN
	    FMT='(1X,ES12.6,ES13.4,F9.2,5(ES14.4),F11.2)'
	  ELSE
	    FMT='(1X,F12.6,ES13.4,F9.2,5(ES14.4),F11.2)'
	  END IF
	  WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV
C
	END DO
C
C Write out a summary of momentum euqation etc to remind the user of the
C meaning of individual terms.
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )'
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,1PE14.4)')
	1          'Surface gravity is: ',GRAV_CON/R(ND)/R(ND)
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',STARS_MASS,' Msun'
C
	CLOSE(LU_OUT)
	RETURN
	END
