!
! Subroutine to compute CLUMP_FAC, and the abundance vectors for
! each species. R and V must be defined. This subroutine can be
! called again if the R grid is altered. 
!
! Notes:
!     AT_ABUND can be a number abundance or a mass fraction (if -ve).
!              On exit, it is always a number fraction.
!     MEAN_ATOMIC_WEIGHT is the mean atomic/ionic weight in amu. Electrons
!              are not included.
!     ABUND_SUM is the total number of atoms present on the scale where some
!         species (usually H or He) is one.
!
	SUBROUTINE SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	USE CONTROL_VARIABLE_MOD
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 19-Dec-2004
!
	INTEGER ND
	INTEGER LUER
	REAL*8 MEAN_ATOMIC_WEIGHT	!In amu (does not include electrons)
	REAL*8 ABUND_SUM
!
! Local variables.
!
	REAL*8 T1,T2,T3,T4
!
	INTEGER I
	INTEGER ISPEC
	INTEGER K
!
        REAL*8 ATOMIC_MASS_UNIT
	EXTERNAL ATOMIC_MASS_UNIT
!
!
! Allow for the possibility that the material in the wind is clumped. The
! assumption that goes into our treatment are not rigorous --- rather this
! treatment should be used as a gauge to discern the errors arising by treating
! O and W-R winds as homogeneous.
!
	IF(DO_CLUMP_MODEL)THEN
	  IF(GLOBAL_LINE_SWITCH .NE. 'BLANK' .AND.
	1                        .NOT. FLUX_CAL_ONLY)THEN
	    WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Clumping is only treated for a fully blanketed model'
	    STOP
	  ELSE IF(GLOBAL_LINE_SWITCH .NE. 'BLANK')THEN
	    WRITE(LUER,*)'Warning in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Clumping is only treated for a fully blanketed model'
	    WRITE(LUER,*)'EW in EWDATA file are incorrect'
	    WRITE(LUER,*)'Continuum fluxes will be okay'
	  END IF
	  IF(CLUMP_LAW(1:4) .EQ. 'EXPO')THEN
!
! CLUMP_PAR(1) is the clumping factor at infinity.
! CLUMP_PAR(2) is a velocity, and determines how fast the clumping factor
! approach CLUMP_PAR(1).
!
	    IF(CLUMP_PAR(3) .EQ. 0.0D0)CLUMP_PAR(4)=1.0D0
	    DO K=1,ND
	      CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0D0-CLUMP_PAR(1)-CLUMP_PAR(3))*
	1                     EXP(-V(K)/CLUMP_PAR(2))+
	1                     CLUMP_PAR(3)*EXP(-V(K)/CLUMP_PAR(4))
	    END DO
	  ELSE
	    WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Invalid law for computing clumping factor'
	    STOP
	  END IF
	ELSE
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0D0
	  END DO
	END IF
! 
! Compute the atomic number density as function of radius. The abundances 
! can be specified in two ways:
!     (1) As a fractional abundance relative to some arbitrary species X.
!     (2) As a mass fraction: Assumed when the abundance is negative.
!
!      N#= SUM[ Xi/X ]            ; Xi/X > 0
!      mu#= SUM[ Ai Xi/X ] / N#   ; Xi/X > 0
!      F#= SUM[Fi]=SUM[ Xi/X ]    ; Xi/X < 0
!      1/mu = f#/mu# + SUM[Fi/Ai]
!
	T1=0.0D0		!Sum of abundance fractions (-> N#)
	T2=0.0D0		!Sum of A_i * abundance fractions
	T3=0.0D0		!Sum of mass fractions      (-> F#)
	T4=0.0D0                !Sum of f/A
!
	DO ISPEC=1,NUM_SPECIES
	  IF(AT_ABUND(ISPEC) .GT. 0)THEN
	    T1=T1+AT_ABUND(ISPEC)
            T2=T2+AT_MASS(ISPEC)*AT_ABUND(ISPEC)
	  ELSE
	    T3=T3+ABS(AT_ABUND(ISPEC))
	    T4=T4+ABS(AT_ABUND(ISPEC))/AT_MASS(ISPEC)
	  END IF
	END DO
	IF(T1 .NE. 0 .AND. T3 .GT. 1.0D0)THEN
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Sum of mass fractions should be less than unity'
	  WRITE(LUER,*)'Mass fraction=',T3
	  STOP
	END IF
	IF(T1 .EQ. 0.0D0 .AND. ABS(T3-1.0D0) .GT. 0.001D0)THEN
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Sum of mass fractions should be unity'
	  WRITE(LUER,*)'Mass fraction=',T3
	  STOP
	END IF
!
!NB: 1-T3 is the mass fraction of species whose abundances were specified in
!          terms of fractional number densities.
!
	IF(T2 .NE. 0)T2=(1.0D0-T3)*T1/T2		!f#/mu#
	MEAN_ATOMIC_WEIGHT=1.0D0/(T2+T4)
!
! Convert any mass fractions to number fractions. If no number fractions
! are present, the first species with non-zero abundance is arbitrarily set 
! to have an abundance of unity.
!
	ABUND_SUM=0.0D0
	IF(T1 .EQ. 0)T1=1.0D0
	DO ISPEC=1,NUM_SPECIES
	  IF(T2 .EQ. 0)T2=ABS(AT_ABUND(ISPEC))/AT_MASS(ISPEC)
	  IF(AT_ABUND(ISPEC) .LT. 0)AT_ABUND(ISPEC)=
	1         T1*ABS(AT_ABUND(ISPEC))/T2/AT_MASS(ISPEC)
	  ABUND_SUM=ABUND_SUM+AT_ABUND(ISPEC)
	END DO
!
	IF(SN_MODEL)THEN
	  DO K=1,ND
	    POP_ATOM(K)=(RHO_ZERO/MEAN_ATOMIC_WEIGHT)*(R(ND)/R(K))**N_RHO
	    POP_ATOM(K)=POP_ATOM(K)/CLUMP_FAC(K)
	  END DO
	ELSE
	  DO K=1,ND
	    POP_ATOM(K)=RMDOT/MEAN_ATOMIC_WEIGHT/(V(K)*R(K)*R(K)*CLUMP_FAC(K))
	  END DO
	END IF
!
	DO ISPEC=1,NUM_SPECIES
	  DO K=1,ND
	    POP_SPECIES(K,ISPEC)=AT_ABUND(ISPEC)*POP_ATOM(K)/ABUND_SUM
	  END DO
	END DO
	T1=ATOMIC_MASS_UNIT()*MEAN_ATOMIC_WEIGHT
!
	T1=ATOMIC_MASS_UNIT()*MEAN_ATOMIC_WEIGHT
	DO I=1,ND
	  DENSITY(I)=POP_ATOM(I)*T1			!gm/cm^3
	END DO
!
	RETURN
	END
