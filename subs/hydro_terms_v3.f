!
! Subroutine to compare terms in the Hydronamical equations. Routine is
! meant for illustration purposes only.
!
	SUBROUTINE HYDRO_TERMS_V3(POP_ATOM,R,V,T,SIGMA,ED,LUM_STAR,
	1              STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1              FLUX_MEAN_OPAC,ELEC_MEAN_OPAC,
	1              LAM_FLUXMEAN_BAND_END,BAND_FLUXMEAN,
	1              PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NOV,
	1              BAND_FLUX,N_FLUXMEAN_BANDS,LU_OUT,ND)
	IMPLICIT NONE
!
! Altered 11-Mar-2008 : Changed from V2 to V3
!                       PRESSURE_VTURB,PLANE_PARALLEL & PLANE_PARALLEL_NOV inserted into call.
! Altered 16-Jan-2005 : Accuracy of dPdR/ROH improved. Now use monotonic
!                         interpolation. Depth index added output.
! Altered 17-Mar-2001 : Altered V output format
! Altered 27-Oct-2000 : Percentage error improved.
!                       STARS mass now output.
! Altered 24-May-1996 : L GEN_ASCI_OPEN and ERROR_LU installed.
! Altered 06-Mar-1996 : Description of terms added to output.
! Created 28-Mar-1994
!
	INTEGER ND,LU_OUT
!
	REAL(10) POP_ATOM(ND)
	REAL(10) R(ND)
	REAL(10) V(ND)
	REAL(10) T(ND)
	REAL(10) SIGMA(ND)
	REAL(10) ED(ND)
	REAL(10) LUM_STAR(ND)                        
!
	REAL(10) FLUX_MEAN_OPAC(ND)
	REAL(10) ELEC_MEAN_OPAC(ND)
!
	INTEGER N_FLUXMEAN_BANDS
	REAL(10) LAM_FLUXMEAN_BAND_END(ND)
	REAL(10) BAND_FLUXMEAN(ND,N_FLUXMEAN_BANDS)
	REAL(10) BAND_FLUX(ND,N_FLUXMEAN_BANDS)
!
	REAL(10) STARS_MASS
	REAL(10) MEAN_ATOMIC_WEIGHT
	REAL(10) PRESSURE_VTURB
	LOGICAL PLANE_PARALLEL
	LOGICAL PLANE_PARALLEL_NOV
!
! External functions
!
	REAL(10) BOLTZMANN_CONSTANT
	REAL(10) GRAVITATIONAL_CONSTANT
	REAL(10) MASS_SUN
	REAL(10) LUM_SUN
	REAL(10) FUN_PI
	REAL(10) SPEED_OF_LIGHT
	REAL(10) ATOMIC_MASS_UNIT
!
	EXTERNAL BOLTZMANN_CONSTANT,GRAVITATIONAL_CONSTANT,MASS_SUN
	EXTERNAL LUM_SUN,FUN_PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
!
	REAL(10) MOD_PRESSURE(ND)
	REAL(10) COEF(ND,4)
!
! Local variables.
!
	REAL(10) VdVdR
	REAL(10) dPdR_ON_ROH
	REAL(10) g_TOT
	REAL(10) g_ELEC
	REAL(10) g_GRAV
	REAL(10) g_RAD
	REAL(10) ERROR
	REAL(10) RLUMST(ND)
!
	REAL(10) GRAV_CON
	REAL(10) dP_CON
	REAL(10) RAD_CON
	REAL(10) T1
!
	CHARACTER*80 FMT
	INTEGER I,J,IOS,ERROR_LU
	EXTERNAL ERROR_LU
!
	I=MAX(132,(N_FLUXMEAN_BANDS+1)*12+16)
	CALL GEN_ASCI_OPEN(LU_OUT,'HYDRO','UNKNOWN',' ',' ',I,IOS)
	IF(IOS .NE. 0)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error opening HYDRO in HYDRO_TERMS'
	  WRITE(I,*)'IOS=',IOS
	END IF
!                                              
! Output file header.
!
	WRITE(LU_OUT,'(1X,6X,A,6X, 8X,A,3X, 2X,A, 5(5X,A,1X), 4X,A,2X,A)')
	1 'R','V','% Error','   VdVdR',
	1                   'dPdR/ROH',
	1                   '   g_TOT',
	1                   '   g_RAD',
	1                    '  g_ELEC','Gamma','Depth'
!
! 1.0D-06 = [ 10^(-10) {from 1/r) * 10^4 {from T}
	dP_CON=1.0D-06*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/
	1               ATOMIC_MASS_UNIT()
!
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()
!                               
!       1.0E-10 (CHI) / 1.0E-20 (1/r^2)
	RAD_CON=1.0D-30*LUM_SUN()/ATOMIC_MASS_UNIT()/4.0D0/FUN_PI()/
	1                  SPEED_OF_LIGHT()/MEAN_ATOMIC_WEIGHT
! 
	DO I=1,ND
	  MOD_PRESSURE(I)=(POP_ATOM(I)+ED(I))*T(I)
	END DO
	CALL MON_INT_FUNS_V2(COEF,MOD_PRESSURE,R,ND)
!	                                     
	DO I=1,ND
	  VdVdR=V(I)*V(I)*(SIGMA(I)+1.0D0)/R(I)
	  dPdR_ON_ROH=dP_CON*COEF(I,3)/POP_ATOM(I)
!
! Gravitational force per unit mass (in cgs units) and radiation pressure forces.
!
	  IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)THEN
	    g_grav=GRAV_CON/R(ND)/R(ND)
	    g_rad=RAD_CON*LUM_STAR(ND)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	    g_elec=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	  ELSE
	    g_grav=GRAV_CON/R(I)/R(I)
	    g_rad=RAD_CON*LUM_STAR(I)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	    g_elec=RAD_CON*LUM_STAR(I)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  END IF
!
	  g_TOT= g_RAD-g_GRAV
	  ERROR=200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/
	1        ( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )
!                  
	  IF(R(I) .GT. 9.99D+04)THEN
	    FMT='(1X,ES12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	  ELSE
	    FMT='(1X,F12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	  END IF                    
	  WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I
!
	END DO
!
! Write out a summary of momentum euqation etc to remind the user of the
! meaning of individual terms.
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )'
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,1PE14.4)')
	1          'Surface gravity is: ',GRAV_CON/R(ND)/R(ND)
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',STARS_MASS,' Msun'
!
	RLUMST(:)=BAND_FLUX(:,N_FLUXMEAN_BANDS)
        DO J=N_FLUXMEAN_BANDS,2,-1
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)-BAND_FLUXMEAN(:,J-1)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)-BAND_FLUX(:,J-1)
        END DO
        DO J=1,N_FLUXMEAN_BANDS
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)/FLUX_MEAN_OPAC(:)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)/RLUMST(:)
        END DO
!
	WRITE(LU_OUT,'(A)')CHAR(12)                         !Formfeed
        WRITE(LU_OUT,'(3X,A,5X,A,6X,14F12.2))')'I','V',
	1                      (LAM_FLUXMEAN_BAND_END(J),J=1,N_FLUXMEAN_BANDS)
        DO I=1,ND
          WRITE(LU_OUT,'(X,I3,15ES12.3))')I,V(I),(BAND_FLUXMEAN(I,J),J=1,N_FLUXMEAN_BANDS)
          WRITE(LU_OUT,'(16X,15ES12.3))')(BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
	  T1=FLUX_MEAN_OPAC(I)/ELEC_MEAN_OPAC(I)
	  DO J=1,N_FLUXMEAN_BANDS
	    IF(ABS(BAND_FLUX(I,J)) .LT. 1.0D-80)BAND_FLUX(I,J)=1.0D-80
          END DO 
	  WRITE(LU_OUT,'(16X,15ES12.3))')(T1*BAND_FLUXMEAN(I,J)/BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
        END DO
!
	CLOSE(LU_OUT)
	RETURN
	END
