	MODULE CMF_HYDRO_MODULE
!
	  REAL(10) GAM_EDD
	  REAL(10) GAM_LIM
	  REAL(10) GAM_FULL
	  REAL(10) AMU			!Atomic mass unit
	  REAL(10) MU_ATOM		!Mean atomic mass in amu
	  REAL(10) BC			!Boltzmann constant
	  REAL(10) C_CMS
	  REAL(10) SIGMA_TH		!Thompson cross-section
	  REAL(10) STEFAN_BC		!Stephan-Boltzmann constant
	  REAL(10) LOGG			!Log (surface gravity) (cgs units)
	  REAL(10) MDOT
	  REAL(10) GRAV_CON
	  REAL(10) REFERENCE_RADIUS
	  REAL(10) PREV_REF_RADIUS
	  REAL(10) RADIUS_AT_TAU_23
	  REAL(10) SOUND_SPEED
	  REAL(10) MOD_SOUND_SPEED
	  REAL(10) VTURB
!
! The following quantities are used when integrating the hydrostatic equation.
! KAP is used to refer to mass absorption opacities.
!
	  REAL(10) R_EST
	  REAL(10) T_EST
	  REAL(10) P_EST
	  REAL(10) TAU_EST
	  REAL(10) ED_ON_NA_EST 
	  REAL(10) ATOM_EST
	  REAL(10) TEFF
	  REAL(10) ROSS_ON_ES
	  REAL(10) KAP_ROSS
	  REAL(10) KAP_ES
	  REAL(10) V_EST
	  REAL(10) F_EST
	  REAL(10) T_SCALING_FACTOR
!
	  REAL(10) dPdR
	  REAL(10) dTAUdR
!
	  REAL(10) OLD_TEFF
	  REAL(10) OLD_REF_RADIUS
	  LOGICAL PURE_LTE_EST
	  LOGICAL OLD_MODEL
	  LOGICAL WIND_PRESENT 
	  LOGICAL PLANE_PARALLEL_MOD
	  LOGICAL RESET_REF_RADIUS
	  INTEGER, SAVE :: LAST_HYDRO_ITERATION=0

	END MODULE CMF_HYDRO_MODULE 	
