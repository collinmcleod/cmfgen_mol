	SUBROUTINE RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
	INTEGER LUIN
	INTEGER LUSCR
	INTEGER LUER
	INTEGER NUM_BNDS
!
! Local variables.
!
	CHARACTER(LEN=80)  TMP_STRING
	CHARACTER(LEN=20)  TMP_KEY
	CHARACTER(LEN=132) TEMP_CHAR
!
	INTEGER I
	INTEGER ISPEC
	INTEGER ID
!
	REAL*8 ATOMIC_MASS_UNIT
	EXTERNAL ATOMIC_MASS_UNIT
!
	CALL GEN_ASCI_OPEN(LUIN,'VADAT','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening VADAT in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
C 
C
C Input model parameters and modelling specifications.
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(RP,'RSTAR',L_TRUE,
	1           'Stellar radius (in 10^10 cm)')
	  CALL RD_STORE_DBLE(RMAX,'RMAX',L_TRUE,
	1           'Maximum radius (in R*)')
	  RMAX=RMAX*RP
C
	  CALL RD_STORE_INT(VELTYPE,'VEL_LAW',L_TRUE,
	1           'Velocity Law to be used')
	  IF(VELTYPE .EQ. 1 .OR. VELTYPE .EQ. 2)THEN
	    CALL RD_STORE_DBLE(VRP,'VRP',L_TRUE,'First velocity component')
	    CALL RD_STORE_DBLE(RN,'RN',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF,'VINF',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS1,'EPSS1',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA1,'GAMMA1',L_TRUE,' ')
!
	    CALL RD_STORE_DBLE(RP2,'RP2',L_TRUE,'Second velocity component')
	    CALL RD_STORE_DBLE(VRP2,'VRP2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(RN2,'RN2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS2,'EPPS2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA2,'GAMMA2',L_TRUE,' ')
	    RN=RN*RP
	    RP2=RP*RP2
	    RN2=RN2*RP2
	  ELSE IF(VELTYPE .EQ. 3)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,
	1           'Terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA',L_TRUE,
	1           'Speed of velocity Law')
	    V_EPPS1=1.0D0
	    VINF=VINF1
	    VINF2=VINF1                !i.e. no 2nd component
	    V_BETA2=1.0D0
	    V_EPPS2=1.0D0
	    NBND_INS=1                 !Old default
	    CONS_FOR_R_GRID=1.0D0
	    EXP_FOR_R_GRID=0.0D0
	  ELSE IF(VELTYPE .EQ. 4)THEN
!
! No parameters required
!
	  ELSE IF(VELTYPE .EQ. 5)THEN
	    WRITE(LUER,*)'Velocity law 5 not implemented in this version',
	1                 ' of CMFGEN'
	    STOP
	  ELSE IF(VELTYPE .EQ. 6)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(VINF1,'VINF1',L_TRUE,
	1            'Terminal velocity (km/s) if no 2nd comp.')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,
	1           'Speed of 1st Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS1,'EPPS1',L_TRUE,
	1           'Scale factor for 1s Beta velocity Law')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,
	1           'True terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA2,'BETA2',L_TRUE,
	1           'Speed of 2nd Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS2,'EPPS2',L_TRUE,
	1           'Scale factor for 2nd Beta V law')
!
	    NBND_INS=1          !Old default
	    CONS_FOR_R_GRID=-1.0D0
	    CALL RD_STORE_INT(NBND_INS,'NBND_INS',L_FALSE,
	1           'Number of additional points to insert in radius grid at boundary')
	    CALL RD_STORE_DBLE(CONS_FOR_R_GRID,'C_R_GRID',L_FALSE,
	1           'Constant to allow imprved shoice of R grid')
	    IF(CONS_FOR_R_GRID .GT. 0)THEN
	      CALL RD_STORE_DBLE(EXP_FOR_R_GRID,'E_R_GRID',L_TRUE,
	1           'Constant to allow imprved shoice of R grid')
	    ELSE
	      CONS_FOR_R_GRID=1.0D0
	      EXP_FOR_R_GRID=0.0D0
	    END IF
C
C !Required by routines other than STARPCYG
C
	    VINF=VINF2	
	  ELSE IF(VELTYPE .EQ. 7)THEN
	    CALL RD_STORE_NCHAR(VEL_OPTION,'VEL_OPT',ITEN,L_TRUE,
	1                        'Velocity option: RVSIG_COL or deKOTER')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,
	1           'Terminal velocity (km/s)')
	    VCORE=0.0D0		!Not used but initialized
	    VPHOT=0.0D0
	    SCL_HT=0.0D0
	    V_BETA1=0.0D0
	    V_EPPS1=1.0D0
	    VINF=VINF1
	    VINF2=VINF1		!i.e. no 2nd component
	    V_BETA2=1.0D0
	    V_EPPS2=1.0D0
	  ELSE IF(VELTYPE .EQ. 10)THEN
	    SN_MODEL=.TRUE.
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,'Initial velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,'Power of velocity Law')
	    CALL RD_STORE_DBLE(RHO_ZERO,'RHO_ZERO',L_TRUE,'Initial density (gm/cm^3)')
	    RHO_ZERO=RHO_ZERO/ATOMIC_MASS_UNIT()
	    CALL RD_STORE_DBLE(N_RHO,'N_RHO',L_TRUE,'Density exponent (+ve)')
	    VINF=VCORE*(RMAX/RP)**V_BETA1
	  ELSE
	    WRITE(LUER,*)'Velocity law ',VELTYPE, ' not implemented',
	1                ' in this version of CMFGEN'
	    STOP
	  END IF
C
	  CALL RD_STORE_DBLE(RMDOT,'MDOT',L_TRUE,
	1            'Mass Loss rate (Msun/yr) ')
	  CALL RD_STORE_DBLE(LUM,'LSTAR',L_TRUE,
	1      'Stellar luminosity (Lsun)')
	  CALL RD_STORE_DBLE(STARS_MASS,'MASS',L_TRUE,
	1      'Stellar mass (Msun)')
C
C All clumping parameters are read in, even when CLUMPING is switched off.
C
	  CALL RD_STORE_LOG(DO_CLUMP_MODEL,'DO_CL',L_TRUE,
	1            'Calculate a model with clumping?')
	  CALL RD_STORE_NCHAR(CLUMP_LAW,'CL_LAW',ISIX,L_TRUE,
	1      'Which clumping law is being utilized?')
	  CALL SET_CASE_UP(CLUMP_LAW,IZERO,IZERO)
	  CALL RD_STORE_INT(N_CLUMP_PAR,'N_CL_PAR',L_TRUE,
	1          'Number of clumping parameters')
	  IF(N_CLUMP_PAR .GT. N_CLUMP_PAR_MAX)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)'N_CLUMP_PAR too large: N_CLUMP_PAR=',N_CLUMP_PAR
	    STOP
	  END IF
	  CLUMP_PAR(:)=0.0D0
	  DO I=1,N_CLUMP_PAR			!Should be less than 10
	    TEMP_CHAR='CL_PAR_'
	    WRITE(TEMP_CHAR(8:8),'(I1)')I
	    CALL RD_STORE_DBLE(CLUMP_PAR(I),TEMP_CHAR(1:8),L_TRUE,
	1             'Clumping parameters:')
	  END DO
!
! These calls read in parameters which automatically control the
! revision of the R grid automatically in SN models that have sharp
! ionization fronts.
!
	  REVISE_R_GRID=.FALSE.
	  JGREY_WITH_V_TERMS=.FALSE.
	  N_RG_PAR=0
	  IF(SN_MODEL)THEN
	    CALL RD_STORE_LOG(REVISE_R_GRID,'REV_RGRID',L_TRUE,
	1            'Automatically revise R grid?')
	    CALL RD_STORE_NCHAR(NEW_RGRID_TYPE,'RG_TYPE',ISIX,REVISE_R_GRID,
	1           'Type of new R grid (UNIFORM or FIX_NX)?')
	    CALL SET_CASE_UP(NEW_RGRID_TYPE,IZERO,IZERO)
	    CALL RD_STORE_INT(N_RG_PAR,'N_RG_PAR',L_FALSE,
	1          'Number of parameters used to determinie new R grid')
	    IF(N_RG_PAR .GT. N_RG_PAR_MAX)THEN
	      WRITE(LUER,*)'Error in CMFGEN'
	      WRITE(LUER,*)'N_RG_PAR too large: N_RG_PAR=',N_RG_PAR
	      STOP
	    END IF
	    RG_PAR(:)=0.0D0
	    DO I=1,N_RG_PAR			!Should be less than 10
	      TEMP_CHAR='RG_PAR_'
	      WRITE(TEMP_CHAR(8:8),'(I1)')I
	      CALL RD_STORE_DBLE(RG_PAR(I),TEMP_CHAR(1:8),L_TRUE,'New R grid parameters:')
	    END DO
	    IF(N_RG_PAR .EQ. 0)N_RG_PAR=1
	    JGREY_WITH_V_TERMS=.TRUE.
	    CALL RD_STORE_LOG(JGREY_WITH_V_TERMS,'JG_W_V',L_FALSE,
	1            'Include V terms when computing Grey temperature?')
	  END IF	  
C
C Read in the un-normalized fractional abundances.
C
	  DO ISPEC=1,NUM_SPECIES
	    TMP_KEY=TRIM(SPECIES(ISPEC))//'/X'
	    TMP_STRING=TRIM(SPECIES(ISPEC))//
	1            '/X fractional abundance by number (un-normalized)'
	    CALL RD_STORE_DBLE(AT_ABUND(ISPEC),TMP_KEY,SPECIES_PRES(ISPEC),
	1            TMP_STRING)
	  END DO
	  WRITE(LUSCR,'()')
C
	  CALL RD_STORE_LOG(RD_CONT_FREQ,'RD_CF_FILE',L_TRUE,
	1            'Read in continuum frequencies from file')
	  CALL RD_STORE_DBLE(MIN_CONT_FREQ,'MIN_CF',L_TRUE,
	1            'Minimum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(MAX_CONT_FREQ,'MAX_CF',L_TRUE,
	1            'Maximum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(SMALL_FREQ_RAT,'FRAC_SP',L_TRUE,
	1            'Fractional spacing for small frequencies')
	  CALL RD_STORE_DBLE(BIG_FREQ_AMP,'AMP_FAC',L_TRUE,
	1            'Amplification factor for large frequency ranges')
	  CALL RD_STORE_DBLE(dFREQ_bf_MAX,'MAX_BF',L_TRUE,
	1            'Maximum frequency spacing close to bf edge')
!
! Installed to allow earlier frequency grids to be used. 0 uses the
! the latest default grid.
!
	  FREQ_GRID_OPTION=0
	  CALL RD_STORE_INT(FREQ_GRID_OPTION,'FR_GRID',L_FALSE,
	1            'Which method to compute frequency grid?')
C
	  CALL RD_STORE_LOG(DO_LEV_DISSOLUTION,'DO_DIS',L_TRUE,
	1            'Allow for level dissolution of upper levels?')
	  CALL RD_STORE_DBLE(dV_LEV_DIS,'dV_LEV',L_TRUE,
	1             'Spacing (in km/s) on low side of bf edge for'//
	1             ' level dissolution')
	  CALL RD_STORE_DBLE(AMP_DIS,'AMP_DIS',L_TRUE,
	1            'Amplification factor on low side bf edge')
	  CALL RD_STORE_DBLE(MIN_FREQ_LEV_DIS,'MIN_DIS',L_TRUE,
	1            'Minimum frequency for level dissolution')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_ALL_CROSS,'CROSS',L_TRUE,
	1            'Compute all photoionization cross-sections?')
	  CALL RD_STORE_DBLE(DELV_CONT,'V_CROSS',L_TRUE,
	1            'Max. vel. sep. (km/s) between evaluations of all'//
	1            '  phot. cross-sections?')
	  CALL RD_STORE_LOG(DIE_AS_LINE,'DIE_AS_LINE',L_TRUE,
	1        'Treat the dielectronic transitions as individual lines?')
	  CALL RD_STORE_DBLE(VSM_DIE_KMS,'VSM_DIE',L_TRUE,
	1        'Velocity (km/s) for smoothing dielectronic transitions')
	  CALL RD_STORE_DBLE(SIG_GAU_KMS,'SIG_GAU_KMS',L_TRUE,
	1        'Sigma of Gaussian used to smooth photoionization data')
	  FRAC_SIG_GAU=0.25D0
	  CALL RD_STORE_DBLE(FRAC_SIG_GAU,'FRAC_SIG',L_FALSE,
	1        'Fractional spacing a across smoothing Gauusian (use 0.25)')
	  CUT_ACCURACY=0.02
	  CALL RD_STORE_DBLE(CUT_ACCURACY,'CUT_ACC',L_FALSE,
	1        'Accuracy to retain data when omitting data points to save space (use 0.02)')
	  ABOVE_EDGE=.TRUE.
	  CALL RD_STORE_LOG(ABOVE_EDGE,'ABV_EDGE',L_FALSE,
	1        'Use only data above edge when smoothing (TRUE)')
!
	  CALL RD_STORE_DBLE(EXT_LINE_VAR,'EXT_LINE_VAR',L_TRUE,
	1        'Extent of line variation zone (V/INF) beyond'//
	1        'the resonancze zone')
          IF(EXT_LINE_VAR .LT. 0.0D0 .OR. EXT_LINE_VAR .GT. 2.0D0)THEN
	    WRITE(LUER,*)'Error in CMFGEN --- invalid range for EXT_LINE_VAR'
	    STOP
	  END IF
!
! NB: An ideal vale for ZNET_VAR_LIMIT is probably 0.01 or 0.001. If 
! ZNET_VAR_LIMIT is zero, all depths will be included in the linearization, 
! independent of ZNET. A very large value of ZNET (i.e. 10^4), will imply
! an interation on the NET_RATES, with no linearization.
!
	  CALL RD_STORE_DBLE(ZNET_VAR_LIMIT,'ZNET_VAR_LIM',L_TRUE,
	1            'Include lines in full varaition when '//
	1            ' ABS(ZNET-1) > ZNET_VAR_LIM')
	  CALL RD_STORE_LOG(WEAK_WITH_NET,'WNET',L_TRUE,
	1            'Use Lambda iteration for weak lines?')
	  CALL RD_STORE_DBLE(WEAK_LINE_LIMIT,'WK_LIM',L_TRUE,
	1            'Maximum opacity ratio for weak lines (0.01)?')
C
	  CALL RD_STORE_LOG(DIF,'DIF',L_TRUE,
	1            'Use Diffusion approximation at inner boundary ?')
C
	  CALL RD_STORE_LOG(RD_COHERENT_ES,'COH_ES',L_TRUE,
	1            'Assume coherent electron scattering? ')
	  CALL RD_STORE_LOG(USE_OLDJ_FOR_ES,'OLD_J',L_TRUE,
	1            'Use old file to provide initial estimate of J_ES?')
	  COHERENT_ES=RD_COHERENT_ES
	  CALL RD_STORE_LOG(MIXED_ES_VAR,'MIX_COH',L_TRUE,
	1            'Mix coherent/non-coherent e.s. in linearization?')
	  CALL RD_STORE_DBLE(ES_VAR_FAC,'ES_FAC',L_TRUE,
	1            'Fractional proximity of RJ and RJ_ES for coherent'//
	1            ' variation')
C
	  CALL RD_STORE_NCHAR(METHOD,'METHOD',ISIX,L_TRUE,
	1         'Which method for continuum tau'//
	1         ' loglog, loglin, linear or zero ?')
	  CALL RD_STORE_NCHAR(N_TYPE,'N_TYPE',ISIX,L_TRUE,
	1         'Method for to handle N for MOM_J_CMF -- '//
	1         'N_ON_J, MIXED, or G_ONLY')
	  CALL RD_STORE_NCHAR(FG_SOL_OPTIONS,'FG_OPT',ITEN,L_TRUE,
	1         'Solution options for FG_J_CMF: DIFF/INS and INT/INS')
	  CALL RD_STORE_DBLE(DELV_FRAC_FG,'VFRAC_FG',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in FG_J_CMF_V10')
	  CALL RD_STORE_DBLE(DELV_FRAC_MOM,'VFRAC_MOM',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in MOM_J_CMF_V10')
											
	  CALL RD_STORE_LOG(RDTHK_CONT,'THK_CONT',L_TRUE,
	1           'Use thick boundary condition for continuum ? ')
	  CALL RD_STORE_LOG(TRAPFORJ,'TRAP_J',L_TRUE,
	1           'Use trapazoidal weights to compute J? ')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(TDOP,'TDOP',L_TRUE,
	1      'Temperature to be used in Doppler profile (10^4K)')
	  CALL RD_STORE_DBLE(AMASS_DOP,'AMASS_DOP',L_TRUE,
	1      'Atomic mass to be used in Doppler profile (amu''s)')
	  CALL RD_STORE_DBLE(VTURB,'VTURB',L_TRUE,
	1      'Turbulent velocity to be used in Doppler profile (km/s)')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(MAX_DOP,'MAX_DOP',L_TRUE,
	1      'Maximum half-width of resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(FRAC_DOP,'FRAC_DOP',L_TRUE,
	1      'Spacing in resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(dV_CMF_PROF,'dV_CMF_PROF',L_TRUE,
	1      'Spacing across cmf profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_CMF_WING,'dV_CMF_WING',L_TRUE,
	1      'Spacing across e.s. wings of cmf profile(in km/s)')
	  CALL RD_STORE_DBLE(ES_WING_EXT,'ES_WING_EXT',L_TRUE,
	1      'Extent of BLUE e.s. wings from resonance core (in km/s)')
	  CALL RD_STORE_DBLE(R_CMF_WING_EXT,'R_CMF_WING_EXT',L_TRUE,
	1      'Extent of RED e.s. wings from RESONANCE core (in Vinf)')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(OBS_PRO_EXT_RAT,'OBS_EXT_RAT',L_TRUE,
	1      'Half width of profile in Vinf.')
	  CALL RD_STORE_DBLE(dV_OBS_PROF,'dV_OBS_PROF',L_TRUE,
	1      'Spacing across observed profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_WING,'dV_OBS_WING',L_TRUE,
	1      'Spacing across e.s. wings of observed profile(in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_BIG,'dV_OBS_BIG',L_TRUE,
	1      'Frequency spacing between lines (in km/s)')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(FLUX_CAL_ONLY,'FLUX_CAL_ONLY',L_TRUE,
	1           'Compute the observers frame flux only ?')
	  CALL RD_STORE_LOG(EXTEND_FRM_SOL,'EXT_FRM_SOL',L_TRUE,
	1           'Extrapolate the formal solution to larger radii?')
	  CALL RD_STORE_LOG(INSERT_FREQ_FRM_SOL,'INS_F_FRM_SOL',L_TRUE,
	1           'Insert extra frequencies for formal solution?')
	  CALL RD_STORE_NCHAR(CMF_FORM_OPTIONS,'FRM_OPT',ITEN,L_TRUE,
	1           'Solution options for CMF_FORM_SOL')
	  CALL RD_STORE_LOG(DO_SOBOLEV_LINES,'DO_SOB_LINES',L_TRUE,
	1        'Compute Sobolev rates and EWs for flux calculation?')
	  CALL RD_STORE_LOG(SOB_FREQ_IN_OBS,'SOB_FREQ_IN_OBS',L_TRUE,
	1        ' Allow for SOB & CMF lines in defining observers'//
	1        ' frequencies?')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_SWITCH,'GLOBAL_LINE',ISIX,L_TRUE,
	1            'Global switch to indicate handeling of line')
	  CALL SET_CASE_UP(GLOBAL_LINE_SWITCH,IZERO,IZERO)
	  IF( GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB' .AND.
	1       GLOBAL_LINE_SWITCH(1:3) .NE. 'CMF' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE' .AND.
	1       GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK')THEN
	    WRITE(LUER,*)'Invalid GLOBAL_LINE SWITCH parameter'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(SET_TRANS_TYPE_BY_LAM,'LAM_SET',L_TRUE,
	1         'Set long wavelengths to SOBOLEV approximation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_BEG,'F_LAM_BEG',L_TRUE,
	1         'Inital wavelength (A) for blanketed flux calculation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_END,'F_LAM_END',L_TRUE,
	1         'Final wavelength (A) for blanketed flux calculation')
	  CALL RD_STORE_DBLE(GF_CUT,'GF_CUT',L_TRUE,
	1          'gf value to omit transitions')
	  CALL RD_STORE_DBLE(AT_NO_GF_CUT,'AT_CUT',L_TRUE,
	1          'Only omit transitions if AT_NO >= AT_CUT')
	  CALL RD_STORE_INT(GF_LEV_CUT,'GF_LEV_CUT',L_TRUE,
	1          'Level above whit transitions omitted if gf < GF_CUT')
	  CALL RD_STORE_INT(MIN_NUM_TRANS,'MIN_TRANS',L_TRUE,
	1          'Minimum number of transitions from each level')
C
	  CALL RD_STORE_LOG(THK_LINE,'THK_LINE',L_TRUE,
	1           'Use thick boundary condition for lines?')
	  CALL RD_STORE_LOG(CHECK_LINE_OPAC,'CHK_L_POS',L_TRUE,
	1      'Ensure Line opacity is positive (SOB & CMF modes only)?')
	  CALL RD_STORE_NCHAR(NEG_OPAC_OPTION,'NEG_OPAC_OPT',ITEN,L_TRUE,
	1            'Method for negative opacities in BLANKETING mode')
	  CALL SET_CASE_UP(NEG_OPAC_OPTION,IZERO,IZERO)
	  IF(NEG_OPAC_OPTION .NE. 'SRCE_CHK' .AND. 
	1                           NEG_OPAC_OPTION .NE. 'ESEC_CHK')THEN
	     WRITE(LUER,*)'Error in CMFGEN_SUB'
	     WRITE(LUER,*)'Invalid NEG_OPAC_OPTION'
	     WRITE(LUER,*)'Valid options are SRCE_CHK and ESEC_CHK'
	     STOP
	  END IF
	  CALL RD_STORE_LOG(SETZERO,'He2_RES=0',L_TRUE,
	1           'Set rates in He2 resonance lines to zero ?')
C
	  CALL RD_STORE_LOG(OVERLAP,'ALLOW_OL',L_TRUE,
	1           'Allow for overlap of close lines (SOB only) ?')
	  CALL RD_STORE_DBLE(OVER_FREQ_DIF,'OL_DIF',L_TRUE,
	1           'Max. difference (in km/s) for overlap')
	  OVER_FREQ_DIF=OVER_FREQ_DIF/2.998E+05
C
	  CALL RD_STORE_LOG(INCL_CHG_EXCH,'INC_CHG',L_TRUE,
	1           'Include charge exchange reactions?')
	  CALL RD_STORE_LOG(INCL_TWO_PHOT,'INC_TWO',L_TRUE,
	1           'Include two photon transitions?')
	  CALL RD_STORE_LOG(INCL_RAY_SCAT,'INC_RAY',L_TRUE,
	1           'Include opacity due to Rayleigh scattering?')
	  CALL RD_STORE_LOG(INCL_ADVECTION,'INC_ADV',L_TRUE,
	1           'Include advection terms in rate equations?')
	  CALL RD_STORE_LOG(INCL_ADIABATIC,'INC_AD',L_TRUE,
	1           'Include adiabatic cooling in energy equation')
	  CALL RD_STORE_LOG(SCL_LINE_COOL_RATES,'SCL_LN',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  CALL RD_STORE_DBLE(SCL_LINE_HT_FAC,'SCL_LN_FAC',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  LINEAR_ADV=.TRUE.
	  CALL RD_STORE_LOG(LINEAR_ADV,'LIN_ADV',L_FALSE,
	1           'Compute advection terms using derivatives in linear plane?')
	  ADVEC_RELAX_PARAM=1.0D0
	  CALL RD_STORE_DBLE(ADVEC_RELAX_PARAM,'ADV_RELAX',L_FALSE,
	1           'Parameter to allow advection terms to be included slowly')
!
! Except for the X-ray switch, the X-ray options are only needed if we 
! are including X-rays.
!
	  VSMOOTH_XRAYS=3000.0D0
	  CALL RD_STORE_LOG(XRAYS,'INC_XRAYS',L_TRUE,
	1           'Include X-ray emission')
	  CALL RD_STORE_LOG(FF_XRAYS,'FF_XRAYS',XRAYS,
	1           'Use free-free processes to compute X-ray emission')
	  CALL RD_STORE_LOG(XRAY_SMOOTH_WIND,'X_SM_WIND',XRAYS,
	1           'Ignore clumping when computing X-ray emission')
	  CALL RD_STORE_DBLE(VSMOOTH_XRAYS,'VS_XRAYS',XRAYS,
	1           'X-ray smoothing width for SOB/CMF options')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_1,'FIL_FAC_1',XRAYS,
	1           'Filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(T_SHOCK_1,'T_SHOCK_1',XRAYS,
	1           'Shock T for X-ray emission [1]')
	  CALL RD_STORE_DBLE(V_SHOCK_1,'V_SHOCK_1',XRAYS,
	1           'Cut off velocity for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_2,'FIL_FAC_2',XRAYS,
	1           'Filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(T_SHOCK_2,'T_SHOCK_2',XRAYS,
	1           'Shock T for X-ray emission [2]')
	  CALL RD_STORE_DBLE(V_SHOCK_2,'V_SHOCK_2',XRAYS,
	1           'Cut off velocity for X-ray emission [2]')
	  CALL RD_STORE_LOG(ADD_XRAYS_SLOWLY,'XSLOW',XRAYS,
	1           'Add X-rays by slowly increasing filling factors?')
	  CALL RD_STORE_DBLE(FILL_FAC_X1_BEG,'XFI1_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_X2_BEG,'XFI2_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(SLOW_XRAY_SCL_FAC,'XSCL_FAC',ADD_XRAYS_SLOWLY,
	1           'Rate to increase X-ray filling factor')
!
	  DELV_XRAY=0.5D0*VSMOOTH_XRAYS
	  CALL RD_STORE_DBLE(DELV_XRAY,'V_XRAY',L_FALSE,
	1            'Max. vel. sep. (km/s) between evaluations of '//
	1            '  phot. cross-sections in X-ray region?')
	  NU_XRAY_END=100.0D0
	  CALL RD_STORE_DBLE(NU_XRAY_END,'NU_XRAY',L_FALSE,
	1            'End of X-ray region for continuum definition')
	  
C
	  WRITE(LUSCR,'()')
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG( RDINR,'RD_IN_R_GRID',L_TRUE,
	1        'Read in a predetermined R grid ?')
	  CALL RD_STORE_LOG(GRID,'LIN_INT',L_TRUE,
	1        'Use direct linear interpolation if  new model ?')
	  CALL RD_STORE_LOG(DO_POP_SCALE,'POP_SCALE',L_TRUE,
	1        'Scale populations so that cons. Eq. satisfied ?')
	  CALL RD_STORE_DBLE(T_INIT_TAU,'T_INIT_TAU',L_TRUE,
	1        'Tau above which T is set exactly to T(spherical)')
	  CALL RD_STORE_LOG(ITERATE_INIT_T,'IT_ON_T',L_TRUE,
	1        'Improve initial T estimate by iteration ?')
	  CALL RD_STORE_DBLE(GREY_PAR,'GREY_TAU',L_TRUE,
	1        'SpecifysTau above which T is set to TGREY in iterative process')
C
	  WRITE(LUSCR,'()')
	  DO ID=1,NUM_IONS-1
	    TMP_KEY='TRANS_'//TRIM(ION_ID(ID))
	    TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	    CALL RD_STORE_NCHAR(ATM(ID)%XzV_TRANS_TYPE,TMP_KEY,
	1          ISIX,ATM(ID)%XZV_PRES,TMP_STRING)
	  END DO
C
	  WRITE(LUSCR,'()')
	  DO ID=1,NUM_IONS-1
	    TMP_KEY='DIE_'//TRIM(ION_ID(ID))
	    TMP_STRING='Include (?) LTDR AUOT, WI calc.s for '//TRIM(ION_ID(ID))
	    CALL RD_STORE_2LOG(ATM(ID)%DIE_AUTO_XzV,ATM(ID)%DIE_WI_XzV,
	1         TMP_KEY,ATM(ID)%XZV_PRES,TMP_STRING)
	  END DO
C
	  DO ISPEC=1,NUM_SPECIES
	    WRITE(LUSCR,'()')
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      TMP_KEY='FIX_'//TRIM(ION_ID(ID))
	      TMP_STRING='Fix ? levels of '//TRIM(ION_ID(ID))
	      CALL RD_STORE_INT(ATM(ID)%FIX_NXzV,TMP_KEY,ATM(ID)%XZV_PRES,
	1           TMP_STRING)
	    END DO
	    TMP_KEY='FIX_'//TRIM(SPECIES(ISPEC))
	    TMP_STRING='Fix (?) highest ionization stage in '//TRIM(SPECIES(ISPEC))
	    CALL RD_STORE_INT(FIX_SPECIES(ISPEC),TMP_KEY,SPECIES_PRES(ISPEC),
	1           TMP_STRING)
	  END DO
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(FIXED_NE,'FIX_NE',L_TRUE,
	1                     'Fix the electron density ?')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_IMP,'FIX_IMP',L_TRUE,
	1            'Automatically fix impurity species?')
C
	  WRITE(LUSCR,'()')                             
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_T,'FIX_T',L_TRUE,
	1            'Keep the Temperature fixed ?')
	  CALL RD_STORE_LOG(VARFIXT,'FIX_T_AUTO',L_TRUE,
	1            'Fix the Temperature automatically ?')
	  CALL RD_STORE_DBLE(TAU_SCL_T,'TAU_SCL_T',L_TRUE,
	1      'Electron scattering optical depth from which to fix T')
	  IF(TAU_SCL_T .EQ. 0.0D0)THEN
	    CON_SCL_T=0.0D0
	  ELSE
	    CON_SCL_T=1000.0D0
	  END IF
	  T_MIN=0.0D0
	  CALL RD_STORE_DBLE(T_MIN,'T_MIN',L_FALSE,'Minimum electron temperature')
	  DO_SRCE_VAR_ONLY=.FALSE. 
	  CALL RD_STORE_LOG(DO_SRCE_VAR_ONLY,'SRCE_ONLY',L_FALSE,
	1            'Allow only ths source function to vary?')
C
	  CALL RD_STORE_NCHAR(METH_SOL,'SOL_METH',ISIX,L_TRUE,
	1            'Which Method To solve Matrix Equations'//
	1            ' DIAG, TRI, PEN, GSIT or MIN')
	  IF(NUM_BNDS .EQ. 1 .AND. METH_SOL .NE. 'DIAG')THEN
	    WRITE(LUER,*)'****************************************'
	    WRITE(LUER,*)'******WARNING in CMFGEN*****************'
	    WRITE(LUER,*)'Solution method inconsistent with NUM_BNDS'
	    WRITE(LUER,*)'METH_SOL=',METH_SOL,'NUM_BNDS=',NUM_BNDS
	    METH_SOL='DIAG'
	  END IF
	  CALL RD_STORE_NCHAR(SCALE_OPT,'SCALE_OPT',ISIX,L_TRUE,
	1           'Scale option (LOCAL, NONE or GLOBAL) ? ')
	  CALL RD_STORE_DBLE(EPS,'EPS_TERM',L_TRUE,
	1      'If maximum fractional % change < EPS terminate model ')
	  CALL RD_STORE_DBLE(MAX_LIN_COR,'MAX_LIN',L_TRUE,
	1      'Maximum fractional change for linearization ')
	  CALL RD_STORE_DBLE(MAX_LAM_COR,'MAX_LAM',L_TRUE,
	1      'Maximum fractional change for lambda iteration ')
	  CALL RD_STORE_DBLE(MAX_CHNG_LIM,'MAX_CHNG',L_TRUE,
	1      'If maximum % fractional change > MAX_CHNG terminate model ')
!
	  CALL RD_STORE_LOG(COMPUTE_BARDIN,'COMP_BA',L_TRUE,
	1            'Compute BA matrix ?')
	  CALL RD_STORE_LOG(WRBAMAT_RDIN,'STORE_BA',L_TRUE,
	1      'Store the BA matrix for/during each iteration ? ')
	  CALL RD_STORE_LOG(WR_BA_INV,'STORE_BA_INV',L_TRUE,
	1      'Store the INVERSE of the BA matrix on each iteration ? ')
	  CALL RD_STORE_LOG(WR_PART_OF_INV,'WR_PRT_INV',L_TRUE,
	1      'Store part of the INVERSE to reduce storage TRIDIAG only)?')
	  CALL RD_STORE_INT(N_ITS_TO_FIX_BA,'N_FIX_BA',L_TRUE,
	1      'Number of iterations to hold BA fixed')
	  CALL RD_STORE_DBLE(BA_CHK_FAC,'BA_CHK_FAC',L_TRUE,
	1      'If dJ < BA_CHK_FAC*RJ, ignore correction to BA')
	  CALL RD_STORE_DBLE(VAL_FIX_BA,'FIX_BA',L_TRUE,
	1      'Switch off BA computation if MAXCH< VAL_FIX_BA ')
	  CALL RD_STORE_DBLE(VAL_DO_LAM,'LAM_VAL',L_TRUE,
	1      'Do Lambda iterations if MAXCH > VAL_DO_LAM')
	  CALL RD_STORE_INT(RD_CNT_LAM,'NUM_LAM',L_TRUE,
	1      '# of Lambda iterations if MAXCH > VAL_DO_LAM')
	  CNT_LAM=0
	  CALL RD_STORE_LOG(RDINSOL,'RD_SOL',L_TRUE,
	1            'RD in solution vector to update populations')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(EDD_CONT,'JC_W_EDD',L_TRUE,
	1        'Compute continuum intensity using Eddington factors')
	  CALL RD_STORE_LOG(EDD_LINECONT,'JBAR_W_EDD',L_TRUE,
	1    'Compute line continuum intensity using Eddington factors')
	
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(ACCURATE,'INC_GRID',L_TRUE,
	1          'Increase grid size to improve accuracy? ')
	  CALL RD_STORE_LOG(ALL_FREQ,'ALL_FREQ',L_TRUE,
	1          'Increase accuracy for all frequencies?')
	  CALL RD_STORE_DBLE(ACC_FREQ_END,'ACC_END',L_TRUE,
	1          'Increase accuracy for all frequencies < ACC_END?')
	  CALL RD_STORE_INT(NPINS,'N_INS',L_TRUE,
	1          'Number of points to be inserted in higher'//
	1          ' accuracy grid (1, 2 or 3) ')
	  CALL RD_STORE_INT(ST_INTERP_INDX,'ST_INT',L_TRUE,
	1          'Interpolate from ? ')
	  CALL RD_STORE_INT(END_INTERP_INDX,'END_INT',L_TRUE,
	1          'Interpolate to ? ')
	  CALL RD_STORE_INT(DEEP,'ND_QUAD',L_TRUE,
	1         'Quadratic interpolation from ND-? to ND')
	  CALL RD_STORE_NCHAR(INTERP_TYPE,'INTERP_TYPE',10,L_TRUE,
	1         'Perform interpolations in LOG or LIN plane')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_INT(N_PAR,'N_PAR',L_TRUE,
	1    'Rate of BA incrementation by BA_PAR in cont. loop (# of freq)')
C
C Next two variables apply for both ACCURATE and EDDINGTON.
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_EDDFAC,'COMP_F',L_TRUE,
	1      'Compute new Eddington factors (f)')
	  CALL RD_STORE_DBLE(ACC_EDD_FAC,'ACC_F',L_TRUE,
	1      'Accuracy with which to compute the eddington factor f')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(NG_DO,'DO_NG',L_TRUE,
	1         'Perform NG acceleration when applicable ?')
	  CALL RD_STORE_DBLE(VAL_DO_NG,'BEG_NG',L_TRUE,
	1       'Percentage accuracy at which to begin NG acceleration')
	  CALL RD_STORE_INT(IT_TO_BEG_NG,'IBEG_NG',L_TRUE,
	1       'Iteration at which to begin NG acceleration')
	  CALL RD_STORE_INT(NG_BAND_WIDTH,'BW_NG',L_TRUE,
	1       'Depth band width for NG acceleration')
	  CALL RD_STORE_INT(ITS_PER_NG,'ITS/NG',L_TRUE,
	1         'Number of iterations between NG accelerations (>=4)')
	  IF(ITS_PER_NG .LT. 4)THEN
	     WRITE(LUER,*)'Error in CMFGEN - ITS_PER_NG too small'
	     STOP
	  END IF
	  CALL CLEAN_RD_STORE()
C
	CLOSE(UNIT=7)
	END
