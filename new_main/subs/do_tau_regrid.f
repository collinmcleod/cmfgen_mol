!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE DO_TAU_REGRID(POPS,ESEC,R_OLD,DONE_R_REV,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ND,NT

	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) ESEC(ND)			!Electron scattering opacity
	REAL(KIND=LDP) R_OLD(ND)
	LOGICAL DONE_R_REV
!
! For specifying grid.
!
	CHARACTER(LEN=10) GRID_TYPE
	CHARACTER(LEN=10) OUT_BND_OPT			!Outer boundary option
	CHARACTER(LEN=10) IN_BND_OPT			!Inner boundary option
!
! Local variables.
!
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) LOG_R_OLD(ND)
	REAL(KIND=LDP) NEW_V(ND)
	REAL(KIND=LDP) dTAU_OLD(ND)
	REAL(KIND=LDP) TAU_OLD(ND)
	REAL(KIND=LDP) TAU(ND)
!
	REAL(KIND=LDP) TA(ND)			!Work vectors
	REAL(KIND=LDP) TB(ND)
!
! The fine grid (FG) is chosen to cover the ionization front. The default values are
! -2.0 to 1.0D0 in log(TAU) space.
!
	REAL(KIND=LDP) FG_MIN			!Min Tau for FG
	REAL(KIND=LDP) FG_MAX			!Max Tau for FG
	REAL(KIND=LDP) FG_RANGE
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) LOG_TAU
	REAL(KIND=LDP) STRETCH_POW		!Power law exponent to stretch tau scale about 1
!
	REAL(KIND=LDP) OBND_PARAMS(5)		!Parameters specifying grid placement at outer boundary.
	REAL(KIND=LDP) IBND_PARAMS(5)		!Parameters specifying grid placement at nner boundary.
!
	INTEGER NUM_IBND_PARAMS		!Number of points inserted near inner boundary.
	INTEGER NUM_OBND_PARAMS		!Number of points inserted near outer boundary.
!
! The parameters describe on what iteration we should do the first grid adjustement,
! how often we should do  a grid adjustment, and how many grid adjustemnts should  be
! done. NO_R_REV is adjusted downwards after each iteration, and is output to ADJUST_R_DEFAULTS.
!
	INTEGER STRT_R_REV
	INTEGER FREQ_R_REV
	INTEGER NO_R_REV
	INTEGER LU
!
	INTEGER NX
	INTEGER I,I1,I2,J
	INTEGER IOS
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LUSCR=8
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER(LEN=80) STRING
!
	WRITE(T_OUT,*)'Entering DO_TAU_REGRID'
	DONE_R_REV=.FALSE.
!
	CALL GET_LU(LU)
	OPEN(UNIT=LU,FILE='R_REGRIDDING_LOG',STATUS='UNKNOWN',ACTION='WRITE')
!
! Set default parameters
!
	OUT_BND_OPT='DEFAULT'
	IN_BND_OPT='DEFAULT'
!
! Read in options/parameters describing how the gridding will be performed.
!
	I=10; CALL RD_STORE_NCHAR(GRID_TYPE,'GRID_TYPE',I,L_TRUE,'Regridding method: MODUN, UNIFORM, FIX_NX')
	IF(GRID_TYPE .EQ. 'MODUN')THEN
	  STRETCH_POW=1.5_LDP
	  CALL RD_STORE_NCHAR(STRETCH_POW,'STRETCH',I,L_FALSE,'Exponent to stretch optical depth scale')
	ELSE IF(GRID_TYPE .EQ. 'UNIFORM')THEN
	  STRETCH_POW=1.0_LDP
	ELSE IF(GRID_TYPE .EQ. 'FIX_NX')THEN
	  CALL RD_STORE_INT(NX,'NX',L_TRUE,'Number of grid points in FINE GRID region')
	  FG_MIN=-2.0_LDP; FG_MAX=1.0_LDP
	  CALL RD_STORE_DBLE(FG_MIN,'FG_MIN',L_FALSE,'Minimum tau to start the FINE GRID')
	  CALL RD_STORE_DBLE(FG_MAX,'FG_MAX',L_FALSE,'Maximum tau to start the FINE GRID')
	ELSE
	  WRITE(T_OUT,*)'Unreconized GRID_TYE in ADJUST_R_DEFAULTS read by ADJUST_R_GRID_V3'
	  WRITE(T_OUT,*)'Skipping adjustment of R grid.'
	  WRITE(T_OUT,*)'GRID_TYPE=',TRIM(GRID_TYPE)
	  RETURN
	END IF
!
	I=10; CALL RD_STORE_NCHAR(OUT_BND_OPT,'OB_OPT',I,L_FALSE,'Outer boundary option: SPECIFY or DEFAULT')
	IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	   NUM_OBND_PARAMS=2
	ELSE
	  CALL RD_STORE_INT(NUM_OBND_PARAMS,'NOB_PARS',L_FALSE,'Number of outer boudary parameters')
	  DO I=1,NUM_OBND_PARAMS
	    WRITE(STRING,'(I3)')I
	    STRING='OB_P'//ADJUSTL(STRING)
	    CALL RD_STORE_DBLE(OBND_PARAMS(I),TRIM(STRING),L_TRUE,'Parameters for outer boundary condition')
	  END DO
	END IF
!
	I=10; CALL RD_STORE_NCHAR(IN_BND_OPT,'IB_OPT',I,L_FALSE,'Inner boundary option: SPECIFY or DEFAULT')
	IF(IN_BND_OPT .EQ. 'DEFAULT')THEN
	   NUM_IBND_PARAMS=2
	ELSE
	  CALL RD_STORE_INT(NUM_IBND_PARAMS,'NIB_PARS',L_FALSE,'Number of inner boundary parameters')
	  DO I=1,NUM_IBND_PARAMS
	    WRITE(STRING,'(I3)')I
	    STRING='IB_P'//ADJUSTL(STRING)
	    CALL RD_STORE_DBLE(IBND_PARAMS(I),TRIM(STRING),L_TRUE,'Paremeters for inner boundary condition')
	  END DO
	END IF
!
! Compute optical depth depth increments. We use the FLUX_MEAN optical depth
! scale, except if it has a problem and is less than ESEC.
!
	DO I=1,ND
	  TA(I)=MAX( ABS(FLUX_MEAN(I)) , ESEC(I) )
	  TA(I)=TA(I)*CLUMP_FAC(I)
	END DO
        TB(1:ND)=0.0_LDP                              !Used for dCHIdR
        CALL NORDTAU(dTAU_OLD,TA,R,R,TB,ND)
!
! Compute optical depth scale. Note that we are given the optical depth
! increments, not the optical depth scale. We assume a power law dependence
! for the opacity when evaluating the optical depth at the outer boundary.
!
	T1=TA(1)/TA(5)			!TA takes clumping into account
	IF(T1 .GT. 0.0_LDP)THEN
	  T1=LOG10(T1)/LOG10(R(5)/R(1))
	ELSE
	  T1=8.0_LDP
	END IF
	IF(T1 .LT. 2.0_LDP)T1=2.0_LDP
	TAU_OLD(1)=FLUX_MEAN(1)*R(1)/(T1-1.0_LDP)
	DO I=2,ND
	  TAU_OLD(I)=TAU_OLD(I-1)+dTAU_OLD(I-1)
	END DO
	WRITE(LU,'(A,ES12.4)')'! Tau(min)=',TAU_OLD(1)
	WRITE(LU,'(A,ES12.4)')'! Tau(max)=',TAU_OLD(ND)
	TAU_OLD(1:ND)=LOG10(TAU_OLD(1:ND))
!
! Save existing grid, which will be used for the interplations.
!
	R_OLD(1:ND)=R(1:ND)
	LOG_R_OLD=LOG10(R_OLD)
!
! We modify the tau scale so that extra points are inserted around
! TAU=1.0. An STRETCH_POW of 0.8 is reasonable. STRETCH_POW=1.0 gives a uniform
! scale, exept at the outer boundaries.
!
	IF( TRIM(GRID_TYPE) .EQ. 'MODUN')THEN
	  IF(STRETCH_POW .GT. 1 .OR. STRETCH_POW .LT. 0.2_LDP)THEN
	    WRITE(T_OUT,*)'Error in ADJUST_R_GRID_V3'
	    WRITE(T_OUT,*)'Error --- STRETCH_POW outside expected range'
	    WRITE(T_OUT,*)'STRETCH_POW read=',STRETCH_POW
	    STOP
	  END IF
!
! Recall TAU_OLD is in log space. This if Log(tau) goes from -3 to 3, puting
! STRETCH=2 will cause it to go from -6 to 6.
!
	  TAU_OLD=TAU_OLD*STRETCH_POW
	END IF
!
! Uniform grid in LOG(TAU) [UNIFORM] or uniform in LOG(TAU)**STRETCH_POW [MOD_UN]
!
	IF( TRIM(GRID_TYPE) .EQ. 'UNIFORM' .OR. TRIM(GRID_TYPE) .EQ. 'MODUN')THEN
!
! Set grid away from boundaries.
!
	  LOG_TAU=(TAU_OLD(ND)-TAU_OLD(1))/(ND-NUM_OBND_PARAMS-NUM_IBND_PARAMS-1)
	  TAU(1)=TAU_OLD(1)
	  DO I=NUM_OBND_PARAMS+2,ND-NUM_IBND_PARAMS-1
	    TAU(I)=TAU_OLD(1)+LOG_TAU*(I-NUM_IBND_PARAMS-1)
	  END DO
	  TAU(ND)=TAU_OLD(ND)
!
	ELSE IF(TRIM(GRID_TYPE) .EQ. 'FIX_NX')THEN
!
! This option allows us to specify a particular region in which extra points
! are inserted. The grid could have slight jumps in spacing across the
! specied region.
!
! Adjust FG_MIN and FG_MAX so that we don't mess up near the boundaries.
!
	  T1=MAX(FG_MIN,FG_MAX); FG_MIN=MIN(FG_MIN,FG_MAX); FG_MAX=T1
	  I=1
	  DO WHILE(FG_MIN .GT. TAU_OLD(I) .AND. I .LT. ND)
	    I=I+1
	  END DO
	  IF(I .GT. ND-3)THEN
	    WRITE(T_OUT,*)'Error -- FG_MIN, FG_MAX outside range covered by OLD_TAU'
	    STOP
	    RETURN
	  ELSE IF(I .LE. NUM_OBND_PARAMS)THEN
	    I=NUM_OBND_PARAMS+1; FG_MIN=TAU_OLD(I)
	    WRITE(LU,'(A,ES10.3,A,ES10.3)')'! Resetting FG_MIN to',FG_MIN,
	1                        'since Log(Tau(1))=',TAU_OLD(1)
	  END IF
!
	  DO WHILE(FG_MAX .GT. TAU_OLD(I) .AND. I .LT. ND)
	    I=I+1
	  END DO
	  IF(I .GE. ND-NUM_IBND_PARAMS)THEN
	    I=ND-NUM_IBND_PARAMS-1; FG_MIN=TAU_OLD(I)
	    WRITE(LU,'(A,ES10.3,A,ES10.3)')'! Resetting FG_MAX to',FG_MAX,
	1                        'since Log(Tau(ND))=',TAU_OLD(ND)
	  END IF
	  FG_RANGE=FG_MAX-FG_MIN
	  WRITE(LU,'(A,ES10.2)')'!  Minimum optical depth is',FG_MIN
	  WRITE(LU,'(A,ES10.2)')'!  Maximum optical depth is',FG_MAX
	  WRITE(LU,'(A,ES10.2)')'! Range of optical depth is',FG_RANGE
!
! Compute the new grid. In the present version, NX points are used
! to cover the range log(TAU)=FG_MIN and FG_MAX inclusively. Outside
! this region we use a uniform grid in Log TAU, except we add extra points
! ar either boundary.
!
	  T1=TAU_OLD(ND)-TAU_OLD(1)-FG_RANGE
          LOG_TAU=T1/(ND-NUM_IBND_PARAMS-NUM_OBND_PARAMS-NX)
	  I1=(FG_MIN-TAU_OLD(1))/LOG_TAU
	  LOG_TAU=(FG_MIN-TAU_OLD(1))/I1
          I1=I1+NUM_OBND_PARAMS
	  TAU(1)=TAU_OLD(1)
	  DO I=NUM_OBND_PARAMS+2,I1
	    TAU(I)=TAU_OLD(1)+LOG_TAU*(I-3)
	  END DO
	  WRITE(LU,'(A,I3)')'! Number of points from outer boundary to zone is ',I1
	  WRITE(LU,'(A,I3)')'! Number of points in zone is',NX
!
! Do the insertion in the crtical section.
!
	  LOG_TAU=FG_RANGE/(NX-1)
	  DO I=I1+1,I1+NX
	    TAU(I)=FG_MIN+LOG_TAU*(I-I1-1)
	  END DO
!
! Now do the last section, towards the inner boundary.
!
	  I2=ND-(NX+I1+NUM_IBND_PARAMS)
	  dLOG_TAU=(TAU_OLD(ND)-FG_MAX)/I2
	  DO I=I1+NX+1,ND-NUM_IBND_PARAMS-1
	    TAU(I)=TAU(I-1)+LOG_TAU
	  END DO
	  TAU(ND)=TAU_OLD(ND)
	  WRITE(LU,'(A,I3)')'! Number of points from zone to inner boundary is ',I2
!
	ELSE
	  WRITE(T_OUT,*)'Invalid GRID_TYPE in ADJUST_R_GRID'
	  WRITE(T_OUT,*)'GRID_TYPE=',TRIM(GRID_TYPE)
	  STOP
	END IF
!
! Do grid at inner boundary. We set NUM_IBND_PARAMS extra points.
!
	dLOG_TAU=TAU(ND)-TAU(ND-NUM_IBND_PARAMS-1)
	IF(IN_BND_OPT .EQ. 'DEFAULT')THEN
	  TAU(ND-1)=TAU(ND)-0.1_LDP*dLOG_TAU
	  TAU(ND-2)=TAU(ND)-0.35_LDP*dLOG_TAU
	ELSE IF(IN_BND_OPT .EQ. 'SPECIFY')THEN
	  DO J=1,NUM_IBND_PARAMS
	    TAU(ND-J)=TAU(ND)-dLOG_TAU/IBND_PARAMS(J)
	  END DO
	ELSE
	  WRITE(T_OUT,*)'Invalid inner boundary option in ADJUST_R_GRID_V3'
	  WRITE(T_OUT,*)'IN_BND_OPT = ',TRIM(IN_BND_OPT)
	  RETURN
	END IF
	WRITE(LU,'(A)')'! Done inner boundary in DO_TAU_REGRID'
!
! Now do outer boundary.
!
	dLOG_TAU=TAU(NUM_OBND_PARAMS+2)-TAU(1)
	IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	  TAU(2)=TAU(1)+0.02_LDP*dLOG_TAU
	  TAU(3)=TAU(1)+0.2_LDP*dLOG_TAU
	ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A,A,3(7X,A))')'! ','  J',' TAU(1)','dLOGTAU','OBND(1)'
	  DO J=1,NUM_OBND_PARAMS
	    WRITE(LU,'(A,I3,3ES14.4)'),'! ',J,TAU(1),dLOG_TAU,OBND_PARAMS(J)
	    TAU(1+J)=TAU(1)+dLOG_TAU/OBND_PARAMS(J)
	  END DO
	  WRITE(LU,'(A)')'!'
	ELSE
	  WRITE(T_OUT,*)'Invalid outer boundary option in ADJUST_R_GRID_V3'
	  WRITE(T_OUT,*)'OUT_BND_OPT = ',TRIM(OUT_BND_OPT)
	  RETURN
	END IF
	WRITE(LU,'(A)')'! Done outer boundary in DO_TAU_REGRID'
!
! Compute the new radius grid. Linear interpolation is
! more than adequate, since we're just defining a new grid.
! We make sure boundary valuese are set accurately.
!
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A)')'! All Logs are base 10.'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,6(7X,A))')'!Index','      TAU','     dTAU',' TAU(old)',
	1                               'dTAU(old)','Log(Rold)','   V(old)'
	WRITE(LU,'(A)')'!'
	DO I=1,ND
	  WRITE(LU,'(I6,6ES16.5)')I,TAU(I),TAU(MIN(I+1,ND))-TAU(I),
	1                              TAU_OLD(I),TAU_OLD(MIN(I+1,ND))-TAU_OLD(I),LOG_R_OLD(I),V(I)
	END DO
	CALL LININT(TAU,LOG_R,ND,TAU_OLD,LOG_R_OLD,ND)
	CALL LININT(R,NEW_V,ND,R_OLD,V,ND)
	LOG_R(1)=LOG_R_OLD(1)
	LOG_R(ND)=LOG_R_OLD(ND)
	R=10**(LOG_R)
	R(1)=R_OLD(1)
	R(ND)=R_OLD(ND)
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A)')'! Simple linear interplation used to get V from R.'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,5(7X,A))')'!Index','        R','     Rold','   Log(R)','Log(Rold)','        V'
	WRITE(LU,'(A)')'!'
	DO I=1,ND
	  WRITE(LU,'(I6,5ES16.5)')I,R(I),R_OLD(I),LOG_R(I),LOG_R_OLD(I),NEW_V(I)
	END DO
	CLOSE(LU)
	DONE_R_REV=.TRUE.
!
	RETURN
	END
