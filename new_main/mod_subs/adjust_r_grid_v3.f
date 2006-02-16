!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE ADJUST_R_GRID_V3(POPS,FLUXMEAN,ESEC,
	1             GRID_TYPE,RG_PARS,N_PARS,ND,NT,NC,NP)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! STILL UNDER development      
!
! Altered :  9-Jan-2005 SIMPLE option replaced by UNIFORM. Also
!                          include MODUN option. Boundries handled
!                          differently. Changed to V3 -- parameters
!                          removed from call.
! Altered : 27-Dec-2005 Tau at outer boudary computed differently.
!                          SIMPLE option added.
! Altered : 02-May-2004 Now only change, R, V, SIGMA and POPS.
! Altered : 18-Mar-2004 FLUXMEAN passed insted of dTAU_OLD
!                       ESEC passed in call.
!                       As still developing, subroutine name NOT changed. 
! Created : 25-Feb-2004
!
	INTEGER ND,NT,NC,NP
!
	REAL*8 POPS(NT,ND)
!
! Opacities.
!
	REAL*8 FLUXMEAN(ND)
	REAL*8 ESEC(ND)
!
! For specifying grid.
!
	INTEGER N_PARS
	CHARACTER(LEN=*) GRID_TYPE
	REAL*8 RG_PARS(N_PARS)
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8 R_OLD(ND)
	REAL*8 LOG_R_OLD(ND)
	REAL*8 LOG_R(ND)
	REAL*8 dTAU_OLD(ND)
	REAL*8 TAU_OLD(ND)
	REAL*8 TAU(ND)
	REAL*8 TA(ND)
	REAL*8 TB(ND)
!
! The fine grid (FG) is chosen to cover the ionization front. The default valuse are
! -2.0 to 1.0D0 in log(TAU) space.
!
	REAL*8 FG_MIN                   !Default=-3.0D0
	REAL*8 FG_MAX			!Default=1.5D0
	REAL*8 FG_RANGE			!Default=4.5D0
!
	REAL*8 DLOG_TAU
	REAL*8 DTAU
	REAL*8 T1,T2
	REAL*8 ALPHA
!
	INTEGER N_IB_INS		!Number of points inserted near inner boundary.
	INTEGER N_OB_INS		!Number of points inserted near outer boundary.
	INTEGER I,I1,I2,J
	INTEGER NX
	INTEGER, PARAMETER :: IONE=1
!
	LUER=ERROR_LU()
	WRITE(139,*)'Entering ADJUST_R_GRID'
!
! We set these here, as can easily be adjusted.
!
	N_OB_INS=3			!Must be at least 1.
	N_IB_INS=3
!
! Compute optical depth depth increments. We use the FLUXMEAN optical depth 
! scale, except if it has a problem and is less than ESEC.
!
	DO I=1,ND
	  TA(I)=MAX(FLUXMEAN(I),ESEC(I))
	  TA(I)=TA(I)*CLUMP_FAC(I)
	END DO
        TB(1:ND)=0.0D0                              !Used for dCHIdR
        CALL NORDLOG_TAU(dTAU_OLD,TA,R,R,TB,ND)
!
! Compute optical depth scale. Note that we are given the optical depth 
! increments, not the optical depth scale. We assume a power law dependence
! for the opacity when evaluating the optical depth at the outer boundary.
!
	T1=TA(1)/TA(5)			!TA takes clumping into account
	IF(T1 .GT. 0.0D0)THEN
	  T1=LOG10(T1)/LOG10(R(5)/R(1))
	ELSE
	  T1=8.0D0
	END IF
	IF(T1 .LT. 2.0D0)T1=2.0D0
	TAU_OLD(1:ND)=0.0D0
	TAU_OLD(1)=FLUXMEAN(1)*R(1)/(T1-1.0D0)
	DO I=2,ND
	  TAU_OLD(I)=TAU_OLD(I-1)+dTAU_OLD(I-1)
	END DO
	TAU_OLD(1:ND)=DLOG10(TAU_OLD(1:ND))
!
! Save existing grid, which will be used for the interplations.
!
	R_OLD(1:ND)=R(1:ND)
	LOG_R_OLD=LOG(R_OLD)
!
! We modify the tau scale so that extra points are inserted around
! TAU=1.0. An ALPHA of 0.8 is reasonable. ALPHA=1.0 gives a uniform
! scale, exept at the outer boundaries.
!
	IF( TRIM(GRID_TYPE) .EQ. 'MODUN')THEN
	  IF(N_PARS .NE. 1)THEN
	    WRITE(LUER,*)'Error in ADJUST_R_GRID_V3'
	    WRITE(LUER,*)'Error --- N_PARS must be 1 for MODUN option'
	    WRITE(LUER,*)'N_PARS=',N_PARS
	    STOP
	  END IF
	  ALPHA=RG_PARS(1)
	  IF(ALPHA .GT. 1 .OR. ALPHA .LT. 0.2D0)THEN
	    WRITE(LUER,*)'Error in ADJUST_R_GRID_V3'
	    WRITE(LUER,*)'Error --- ALPHA outside expected range'
	    WRITE(LUER,*)'ALPHA read=',ALPHA
	    STOP
	  END IF
	  TAU_OLD=TAU_OLD*ABS(TAU_OLD)**(ALPHA-1.0D0)
	END IF
!
! Uniform grid in LOG(TAU) [UNIFORM] or uniform in LOG(TAU)**ALPHA [MOD_UN]
!
	IF( TRIM(GRID_TYPE(1:6)) .EQ. 'UNIFORM' .OR. TRIM(GRID_TYPE(1:6)) .EQ. 'MODUN')THEN
!
! Do outer boundary. We insert a very close point at outer boundary.
!
	  DLOG_TAU=(TAU_OLD(ND)-TAU_OLD(1))/(ND-N_OB_INS-N_IB_INS-1)
	  TAU(1)=TAU_OLD(1)
	  T2=0.02D0*ALPHA
	  TAU(2)=TAU(1)+T2*DLOG_TAU
	  TAU(2)=MIN(TAU(2),TAU_OLD(2))
	  T1=(1.0D0-T2)*DLOG_TAU/(2**N_OB_INS-1.0D0)
	  DO J=1,N_OB_INS-1
	    TAU(2+J)=TAU(1+J)+T1
	    T1=T1*2.0D0
	  END DO
!
! Now set grid away from boundaries.
!
	  DO I=N_OB_INS+2,ND-N_IB_INS-1
	    TAU(I)=TAU_OLD(1)+DLOG_TAU*(I-N_OB_INS-1)
	  END DO
!
! Do grid at inner boundary. We set N_IB_INS extra points, such that
! the step size varies by a factor of 2. 
!
	  TAU(ND)=TAU_OLD(ND)
	  T1=10.0D0**(TAU(ND))
	  DTAU=T1*(1.0D0-10.0D0**(-DLOG_TAU))
	  DTAU=DTAU/(2**(N_IB_INS+1)-1.0D0)
	  DO J=1,N_IB_INS
	    TAU(ND-J)=LOG10(T1-DTAU)
	    DTAU=DTAU*(2**(J+1)-1.0D0)
	  END DO
!
	ELSE IF(TRIM(GRID_TYPE) .EQ. 'FIX_NX')THEN
!
! This option allows us to specify a particular region in which extra points 
! are inserted. The grid could have slight jumps in spacing across the
! specied region.
!
	  FG_MIN=-2.0D0			!Default parameters
	  FG_MAX=1.00D0
	  NX=RG_PARS(1)
	  IF(N_PARS .EQ. 3)THEN
	    FG_MIN=RG_PARS(2)
	    FG_MAX=RG_PARS(3)
	  ELSE IF(N_PARS .NE. 1)THEN
	    WRITE(LUER,*)'Error in ADJUST_R_GRID_V3'
	    WRITE(LUER,*)'Error --- N_PARS must be 1 or 3 for FIX_NX option'
	    WRITE(LUER,*)'N_PARS=',N_PARS
	    STOP
	  END IF
	  FG_RANGE=FG_MAX-FG_MIN
	  I=MIN(ND/3,10)
	  IF(FG_MIN .LT. TAU_OLD(1))FG_MIN=TAU_OLD(I)
	  IF(FG_MAX .GT. TAU_OLD(ND))FG_MAX=TAU_OLD(ND-I)
!
! Compute the new grid. In the present version, NX points are used
! to cover the range log(TAU)=FG_MIN and FG_MAX. Outside this region
! we use a uniform grid in Log TAU, except we add extra points
! ar either boundary.
!
	  T1=TAU_OLD(ND)-TAU_OLD(1)-FG_RANGE
          DLOG_TAU=T1/(ND-N_IB_INS-N_OB_INS-NX)
	  I1=(FG_MIN-TAU_OLD(1))/DLOG_TAU
	  DLOG_TAU=(FG_MIN-TAU_OLD(1))/I1
          I1=I1+N_OB_INS
	  TAU(1)=TAU_OLD(1)
	  DO I=N_OB_INS+2,I1
	    TAU(I)=TAU_OLD(1)+DLOG_TAU*(I-3)
	  END DO
          TAU(2)=MIN(TAU_OLD(2),TAU(1)+0.02D0*DLOG_TAU)
!
! Do the insertion in the crtical section.
!
	  DLOG_TAU=FG_RANGE/(NX-1)
	  DO I=I1+1,I1+NX
	    TAU(I)=FG_MIN+DLOG_TAU*(I-I1-1)
	  END DO
!
! Now do the last section, towards the inner boundary.
!
	  I2=ND-(NX+I1+2)
	  T1=(TAU_OLD(ND)-FG_MAX)
	  DLOG_TAU=T1/I2
	  DO I=I1+NX+1,ND-N_IB_INS-1
	    TAU(I)=TAU(I-1)+DLOG_TAU
	  END DO
	  TAU(ND)=TAU_OLD(ND)
	  T2=0.02D0
	  TAU(2)=TAU(1)+T2*DLOG_TAU
	  TAU(2)=MIN(TAU(2),TAU_OLD(2))
	  T1=(1.0D0-T2)*DLOG_TAU/(2**N_OB_INS-1.0D0)
	  DO J=1,N_OB_INS-1
	    TAU(2+J)=TAU(1+J)+T1
	    T1=T1*2.0D0
	  END DO
!
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid GRID_TYPE in ADJUST_R_GRID'
	  WRITE(LUER,*)'GRID_TYPE=',TRIM(GRID_TYPE)
	  STOP
	END IF 
!
! Compute the new radius grid. Linear interpolation is
! more than adequate, since we're just defining a new grid.
! We make sure boundary valuese are set accurately.
!
	DO I=1,ND
	  WRITE(6,'(I4,4E16.8)')I,TAU(I),TAU_OLD(I),LOG_R_OLD(I)
	END DO
	CALL LININT(TAU,LOG_R,ND,TAU_OLD,LOG_R_OLD,ND)
	LOG_R(1)=LOG_R_OLD(1)
	LOG_R(ND)=LOG_R_OLD(ND)
	R=EXP(LOG_R)
	R(1)=R_OLD(1)
	R(ND)=R_OLD(ND)
	DO I=1,ND
	  WRITE(6,'(I4,4ES16.9)')I,R(I),R_OLD(I),LOG_R(I),LOG_R_OLD(I)
	END DO
!
! Before refining the grid, we check to see if it is really necessary.
!
	T1=0.0D0
	DO I=2,ND-1
	  T1=MAX( T1,(R(I)-R_OLD(I))/(R(I-1)-R(I+1)) ) 
	END DO
	T1=T1*2.0D0
	IF(T1 .LT. 1.0D-03)THEN
	  R=R_OLD
	ELSE

! We now need to regrid all the populations. All interpolations (except 
! sigma) are performed in the LOG-LOG plane. For SN this is ideal, since
! the density and velocity are power laws in r. For SIGMA, we do not take
! the log.
!
! We do not need to interpolate T, and ED directly, since these are part
! of POPS.
!
	  TA(1:ND)=LOG(V(1:ND))
	  CALL MON_INTERP(V,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	  V(1:ND)=EXP(V(1:ND))
	  DO I=1,NT
	    TA(1:ND)=LOG(POPS(I,1:ND))
	    CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	    POPS(I,1:ND)=EXP(TB(1:ND))
	  END DO
	END IF
!
	RETURN
	END
