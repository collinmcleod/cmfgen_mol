C
C Data module for CMF_FORM_SOL. Data placed in this module is automatically
C saved between subroutine calls..
C
	MODULE CMF_FORM_MOD
	USE SET_KIND_MODULE
C
	REAL(KIND=LDP), ALLOCATABLE :: R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: P_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ROH_EXT(:)
C
	REAL(KIND=LDP), ALLOCATABLE :: AV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: R_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: Z(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH(:,:)
C
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
C
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_OLD(:)
C
C Used to compute the boundary iteration factors (HBC and NBC).
C
	INTEGER, ALLOCATABLE :: NI_RAY(:)
C
C Interpolations from the original R grid are initially performed using
C monotonic cubic polynomials. To save time, interpolations in CHI and ETA
C are done onto a fine grid --- once for each frequency. Values of CHI and
C ETA on a RAY are then obtained from the fine grid using linear interpolation.
C
	REAL(KIND=LDP), ALLOCATABLE :: R_FINE(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_FINE(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_FINE(:)
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR_FINE(:)
	INTEGER, ALLOCATABLE :: INDX_FINE(:,:)
C
C INDX(I) gives the location of R_RAY(I,...) in R. It is also used to give the
C location of R_FINE(I) in R.
C
	INTEGER, ALLOCATABLE :: INDX(:)
C
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) FL_PREV
	REAL(KIND=LDP) OLD_CHI_AT_IN_BND
	INTEGER I_START
C
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
C
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NRAY
	INTEGER ND_FINE
C
C Used to check for consistency with initialization call.
C
	INTEGER ND_SAV
	INTEGER NC_SAV
	INTEGER NP_SAV
	LOGICAL EXTEND_SAV
C
	INTEGER NFINE_INS
	PARAMETER (NFINE_INS=10)
C
	END MODULE CMF_FORM_MOD
C
C 
C
C Routine to compute the observed intensity IPLUS as a function of
C impact parameter. Subroutine is designed to compute more accurate
C intensities (than FG_J_CMF_V*) for the computation of the observer's frame
C flux. Calculation should only be performed on the last iteration. No
C Eddington factors are returned.
C
C Based on FG_J_CMF_V6.
C
C First order differencing in frequency is used.
C
C The radiation field at the previous (bluer) frequency is stored locally.
C Coupling to this frequency is controlled by the variable INIT. If INIT is
C true, the atmosphere is treated with the assumption that V=0.
C
C INIT must be TRUE on the very first call.
C
C Because of the dynamic memory allocation, the routine CANNOT be called with
C different values for ND, NC, and NP in the same program.
C
	SUBROUTINE CMF_FORM_SOL(ETA_NEW,CHI_NEW,ESEC_NEW,
	1                  ROH,V,SIGMA,R,P,
	1                  P_OBS, IPLUS_P, NP_OBS, NP_OBS_MAX,
	1                  MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                  FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                  EXTEND,INSERT_ADD_FREQ,
	1                  FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                  INIT,NC,NP,ND)
	USE SET_KIND_MODULE
	USE CMF_FORM_MOD
	IMPLICIT NONE
C
C Altered 10-Oct-2004 : Replaced R*R-P*P by (R-P)*(R+P) in Z sqrt calculation.
C Altered 21-Aug-1997 : New version (V2) on NEW_R_SCALE called. This has a
C                         finer grid spacing a the outer boundary which
C                         overcomes numerical difficulties with P Cygni models.
C                       PNT_FAC was installed which allows a the number of
C                         points allong all rays to be altered.
C                       NFINE_INF for fine grid was increased from 4 to 10.
C                         Probably not necessary, but timing should be less
C                         than LS loop.
C                       Minor bug fixed with computation of DCHIDR in extended
C                         part of atmosphere.
C                       CV now computed directly from AV if AV has been
C                         adjusted because a neagtive intensity was obtained.
C Altered Nov-06-1996 : Correction to computation of DTAU installed to prevent
C                        NEGATIVE dTAU's caused by the derivative correction.
C Altered 18-Jun-1996 : Bug fix. Expression for N_FREQ was not being divided
C                         by FL. Consequently inserting too many frequencies
C                         for FL>1, and too few for FL < 1.
C Altered 20/28-28-May-1996 : NEQV used when comparing logicals.
C                             Gnerica calls used for SQRT etc
C                             ETA_OLD, ESEC_OLD, CHIL_OLD initialized
C                               to NEW values when INITis true.
C Finalized 01-Apr-1996
C
	INTEGER ND
	INTEGER NC
	INTEGER NP
	INTEGER NP_OBS
	INTEGER NP_OBS_MAX
C
C _NEW refer to the opacities/emissivities at the frequency passed in the
C call. They are to be distinguished from the opacities/emissivities on the
C previous call [ _OLD], and those at intermediate frequencies (no appendage).
C
	REAL(KIND=LDP) ETA_NEW(ND)	!Emissivity --- contains e.s. term.
	REAL(KIND=LDP) CHI_NEW(ND)
	REAL(KIND=LDP) ESEC_NEW(ND)
C
C We use the density to determine  the grid spacing along a ray. This should
C always be well behaved.
C
	REAL(KIND=LDP) ROH(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) P(NP)
C
	REAL(KIND=LDP) P_OBS(NP_OBS_MAX)
	REAL(KIND=LDP) IPLUS_P(NP_OBS_MAX)
	REAL(KIND=LDP) MU_AT_RMAX(NP_OBS_MAX)
	REAL(KIND=LDP) HQW_AT_RMAX(NP_OBS_MAX)
	REAL(KIND=LDP) RMAX_OBS,V_AT_RMAX
C
C These parameters are use to estimate how many extra points will be inserted
C for the frequency integration.
C
	REAL(KIND=LDP) FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING
C
	REAL(KIND=LDP) DBB,IC,FL,dLOG_NU
	CHARACTER*6 METHOD
	LOGICAL DIF		!Use diffusion approximation
C
C First frequency -- no frequency coupling.

	LOGICAL INIT
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
C If EXTEND is TRUE, atmosphere is extended by a factor of EXT_FAC in Radius.
C An additional ND_ADD_MAX grid points are added.
C
	LOGICAL EXTEND
	LOGICAL INSERT_ADD_FREQ
C
	REAL(KIND=LDP) EXT_FAC
	INTEGER ND_ADD_MAX
	PARAMETER (EXT_FAC=10.0_LDP)
	PARAMETER (ND_ADD_MAX=10)		!10*LOG10(EXT_FAC))
	INTEGER NP_EXT
C
C Intermediate opacities/emissivities. They are obtained by linear
C extrapolation in frequency from .._NEW and ..._OLD.
C
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) ESEC(ND)
C
C The following arrays do not need to be stored, and hence can be created
C dynamically on each call. They must have a MINIMUM length NRAY where
C ND+ND_ADD_MAX (=NRAY) is the (MAXIMUM) number of points along a ray.
C PNT_FAC can be adjusted to different integer values so that more points
C can be inserted along a ray.
C
	INTEGER, PARAMETER :: PNT_FAC=1
	REAL(KIND=LDP) TA(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) TB(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) TC(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) AV(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) CV(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) DTAU(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) XM(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) SOURCE(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) U(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) VB(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) VC(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) GB(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) H(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) Q(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) QH(PNT_FAC*(ND+ND_ADD_MAX))
C
C Ray quantities. V_RAY and SIGMA ray are only used to compute GAM and GAMH
C in the initialization. CHI_RAY, ETA_RAY and dCHIdR are needed in the
C formal solution loop.
C
	REAL(KIND=LDP) V_RAY(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) SIGMA_RAY(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) CHI_RAY(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) ETA_RAY(PNT_FAC*(ND+ND_ADD_MAX))
	REAL(KIND=LDP) dCHIdR(PNT_FAC*(ND+ND_ADD_MAX))
C
C Used to perform interpolation of CHI and ETA onto the FINE radius grid.
C
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) LOG_V(ND)
	REAL(KIND=LDP) LOG_SIGMA(ND)
	REAL(KIND=LDP) LOG_CHI(ND)
	REAL(KIND=LDP) LOG_ETA(ND)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER LU_ER
	INTEGER NI
	INTEGER NI_ORIG
	INTEGER I,J,K,L,ML,LS
	INTEGER N_FREQ
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) IBOUND			!Incident intensity on outer boundary.
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DELR
C
	REAL(KIND=LDP) NEW_dLOG_NU		!dv
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) CV_BOUND
	REAL(KIND=LDP) MU,dZ
C
C Power law exponents for the extrapolations.
C
	REAL(KIND=LDP) ROH_ALPHA
	REAL(KIND=LDP) ESEC_ALPHA
	REAL(KIND=LDP) CHI_ALPHA
	REAL(KIND=LDP) ETA_ALPHA

	EXTERNAL SPEED_OF_LIGHT
	REAL(KIND=LDP) SPEED_OF_LIGHT
C
	IF(FIRST_TIME)THEN
C
C Save important dimensions for checking on subsequent calls.
C
	  ND_SAV=ND
	  NC_SAV=NC
	  NP_SAV=NP
	  EXTEND_SAV=EXTEND
C
	  C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT( )
	  IF(EXTEND)THEN
	    ND_ADD=10*LOG10(EXT_FAC)
	  ELSE
	    ND_ADD=0
	  END IF
C
C Determine the extended R Grid, and the impact parameters (P_EXT) on this
C grid. ROH is assumed to be the opacity variable, and will be used to
C determine the new radius grid on each ray.
C
	  ND_EXT=ND+ND_ADD
	  NP_EXT=ND_EXT+NC
	  ALLOCATE ( R_EXT(ND_EXT) )
	  ALLOCATE ( ROH_EXT(ND_EXT) )
	  ALLOCATE ( P_EXT(NP_EXT) )
	  IF(EXTEND)THEN
	    RMAX=R(1)*EXT_FAC
	    DELR=EXP( LOG(EXT_FAC)/ND_ADD )
	    R_EXT(1)=RMAX
	    DO I=2,ND_ADD
	      R_EXT(I)=R_EXT(I-1)/DELR
	    END DO
	    R_EXT(ND_ADD+1:ND_EXT)=R(1:ND)
C
	    P_EXT(1:NC)=P(1:NC)
	    P_EXT(NC+1:NP_EXT)=R_EXT(ND_EXT:1:-1)
C
	    ROH_ALPHA=LOG(ROH(4)/ROH(1)) / LOG(R(1)/R(4))
	    IF(ROH_ALPHA .LT. 2._LDP)ROH_ALPHA=2.0
	    DO I=1,ND_ADD
	      ROH_EXT(I)=ROH(1)*(R(1)/R_EXT(I))**ROH_ALPHA
	    END DO
	    ROH_EXT(ND_ADD+1:ND_EXT)=ROH(1:ND)
	  ELSE
C
C Even though no extension is requested, R_EXT, ROH_ET and P_EXT must be
C defined. We can use all rays since we insert points along each-ray.
C For the last ray we simply assume IPLUS=0.
C
	    ND_EXT=ND
	    NP_EXT=ND+NC
	    R_EXT(1:ND)=R(1:ND)
	    ROH_EXT(1:ND)=ROH(1:ND)
	    P_EXT(1:NC)=P(1:NC)
	    P_EXT(NC+1:NP_EXT)=R_EXT(ND_EXT:1:-1)
	  END IF
C
C ND_EXT is the number of points in our "EXTENDED" atmosphere.
C NRAY is the maximum number of points along a ray. In general these could be
C different for each ray --- here we make them the same.
C
	  NRAY=PNT_FAC*ND_EXT
C
C Now allocate all other required storage locations.
C
	  ALLOCATE ( CHI_OLD(ND) )
	  ALLOCATE ( ETA_OLD(ND) )
	  ALLOCATE ( ESEC_OLD(ND) )
C
C AV_PREV and CV_PREV contain the radiation field on the RAY grid at the
C previous frequency.
C
	  ALLOCATE ( AV_PREV(NRAY,NP_EXT) )
	  ALLOCATE ( CV_PREV(NRAY,NP_EXT) )
C
	  ALLOCATE ( R_RAY(NRAY,NP_EXT) )
	  ALLOCATE ( Z(NRAY,NP_EXT) )
C
	  ALLOCATE ( GAM(NRAY,NP_EXT) )
	  ALLOCATE ( GAMH(NRAY,NP_EXT) )
C
C CHI_COEF(I,...) gives the coefficients of the cubic polynomial that fits
C CHI for R between R(I) and R(I+1).
C
	  ALLOCATE ( CHI_COEF(ND,4) )
	  ALLOCATE ( ETA_COEF(ND,4) )
C
	  ALLOCATE ( NI_RAY(NP_EXT) )
C
	  NP_OBS=NP_EXT
	  IF(NP_OBS .GT. NP_OBS_MAX)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error in COMP_OBS_INT'
	    WRITE(LU_ER,*)'NP_OBS_MAX too small'
	    STOP
	  END IF
	  P_OBS(1:NP_OBS)=P_EXT(1:NP_EXT)
C
C Find the functions to interpolate V and SIGMA. Since these are only required
C in the INITIALIZATION section we save space by using CHI_COEF for V
C and ETA_COEF for SIGMA.
C
	  LOG_V(1:ND)=LOG(V(1:ND))
	  LOG_SIGMA(1:ND)=LOG( (SIGMA(1:ND)+1.0_LDP) )
	  LOG_R(1:ND)=LOG(R(1:ND))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_V,LOG_R,ND)
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_SIGMA,LOG_R,ND)
C
C Define parameters for extrapolating the velocity law.
C
	  BETA=(SIGMA(1)+1.0_LDP)*(R(1)/R(ND)-1.0_LDP)
          VINF=V(1)/(1-R(ND)/R(1))**BETA
C
C Define the fine radius grid in which linear interpolations will be used.
C NFINE_INS additional points are inserted into each interval in the EXTENDED
C R grid.
C
	  ND_FINE=(NFINE_INS+1)*(ND_EXT-1) + 1
	  ALLOCATE ( R_FINE(ND_FINE))
	  ALLOCATE ( CHI_FINE(ND_FINE))
	  ALLOCATE ( ETA_FINE(ND_FINE))
	  ALLOCATE ( dCHIdR_FINE(ND_FINE))
	  ALLOCATE ( INDX_FINE(NRAY,NP_EXT) )
	  ALLOCATE ( INDX(ND_FINE) )
C
	  DO I=1,ND_EXT-1
	    K=(I-1)*(NFINE_INS+1)+1
	    R_FINE(K)=R_EXT(I)
	    DELR=(R_EXT(I)-R_EXT(I+1))/(NFINE_INS+1)
	    DO J=1,NFINE_INS
	      R_FINE(K+J)=R_EXT(I)-J*DELR
	    END DO
	  END DO
	  R_FINE(ND_FINE)=R_EXT(ND_EXT)
C
C Perform initializations.
C
C Choose the radius grid for each ray. ROH is assumed to be the opacity
C variable --- ROH is a convenient choice as it is well behaved. Points are
C chosen equally spaced in Log TAU. Quantities that are independent of the
C the frequency (and hence the opacity and emissivity) are evaluated. These
C data are stored in the module, and will be available for subsequent
C calls.
C
C Don't do the very last ray, as single point.
C
	  DO LS=1,MIN(NP_OBS,NC+ND_EXT-1)
C
	    NI_ORIG=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC)NI_ORIG=ND_EXT
	    CALL NEW_R_SCALE_V2(R_RAY(1,LS),NRAY,R_EXT,NI_ORIG,ROH_EXT,
	1                  P_EXT(LS),L_FALSE,L_TRUE)
C
C NI and NI_RAY(LS) will now refer the NEW number of points along the ray.
C
	    NI_RAY(LS)=NRAY
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      Z(I,LS)=SQRT( (R_RAY(I,LS)-P_EXT(LS))*(R_RAY(I,LS)+P_EXT(LS)) )
	    END DO
C
C We now have the new radius grid. We now need to interpolate V and SIGMA onto
C this new grid for each ray. Monotonic cubic interpolation is used in the
C original R grid --- power law extrapolation outside this grid.
C
	    IF(EXTEND)THEN
C
C Only use V and SIGMA to compute GAMMA. Thus don't need to save there values
C
	      IF(R_RAY(NI,LS) .GT. R(1))THEN
	        DO I=1,NI					!Additional rays.
	          V_RAY(I)=VINF*(1.0_LDP-R(ND)/R_RAY(I,LS))**BETA
	          SIGMA_RAY(I)=BETA/(R_RAY(I,LS)/R(ND)-1.0_LDP)-1.0_LDP
	        END DO
	        I_START=NI+1
	      ELSE
	        I=0
	        DO WHILE(R_RAY(I+1,LS) .GT. R(1))
	          I=I+1
	          V_RAY(I)=VINF*(1.0_LDP-R(ND)/R_RAY(I,LS))**BETA
	          SIGMA_RAY(I)=BETA/(R_RAY(I,LS)/R(ND)-1.0_LDP)-1.0_LDP
	        END DO
	        I_START=I+1
	     END IF
	   ELSE
	     I_START=1
	   END IF
C
C Now interpolate V and SIGMA onto the RAY grid. We first locate the interval for
C the interpolation. The interpolation loop will then be vectorizable.
C
	    L=1
	    DO I=I_START,NI
	      DO WHILE(R_RAY(I,LS) .LT. R(L+1))
	        L=L+1
	        IF(L .GE. ND)THEN
	          LU_ER=ERROR_LU()
	          WRITE(LU_ER,*)'Invalid interpolation range in CMF_FORM_SOL'
	          STOP
	        END IF
	      END DO
	      INDX(I)=L
	    END DO
C
	    DO I=I_START,NI
	      L=INDX(I)
	      T1=LOG(R_RAY(I,LS)/R(L))
	      V_RAY(I)=EXP(
	1         ((CHI_COEF(L,1)*T1+CHI_COEF(L,2))*T1+
	1         CHI_COEF(L,3))*T1+CHI_COEF(L,4) )
	      SIGMA_RAY(I)=EXP(
	1         ((ETA_COEF(L,1)*T1+ETA_COEF(L,2))*T1+
	1         ETA_COEF(L,3))*T1+ETA_COEF(L,4) )-1.0_LDP
	    END DO
C
C Compute GAMMA. This section is straight from the subroutine GAMMA, except
C That _EXT has been added to V, SIGMA, and R.
C
C We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	    DO I=1,NI_RAY(LS)-1
	      MU=Z(I,LS)/R_RAY(I,LS)
	      GAM(I,LS)=3.33564E-06_LDP*V_RAY(I)/R_RAY(I,LS)*
	1                   (  1.0_LDP+SIGMA_RAY(I)*(MU**2)  )
	      MU=(Z(I,LS)+Z(I+1,LS))/(R_RAY(I,LS)+R_RAY(I+1,LS))
	      GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564E-06_LDP*
	1                   (  1.0_LDP+0.5_LDP*(MU**2)*
	1                   (SIGMA_RAY(I)+SIGMA_RAY(I+1))  )/
	1                   (R_RAY(I,LS)+R_RAY(I+1,LS))
	    END DO
	    GAMH(NI,LS)=0.0
	    MU=Z(NI,LS)/R_RAY(NI,LS)
	    GAM(NI,LS)=3.33564E-06_LDP*V_RAY(NI)*
	1              ( 1.0_LDP+SIGMA_RAY(NI)*(MU**2) )/R_RAY(NI,LS)
C
C Find location of R_RAY on the fine grid. This will be used a lot, and hence
C is done in the initialization.
C
	    L=1
	    DO I=1,NI
	      DO WHILE(R_RAY(I,LS) .LT. R_FINE(L+1))
	        L=L+1
	        IF(L .GE. ND_FINE)THEN
	          LU_ER=ERROR_LU()
	          WRITE(LU_ER,*)'Invalid FINE interpolation ',
	1                         ' range in CMF_FORM_SOL'
	          STOP
	        END IF
	      END DO
	      INDX_FINE(I,LS)=L
	    END DO
	    IF(LS .EQ. 1)V_AT_RMAX=V_RAY(1)
C
	  END DO		!LS Loop
C
	  RMAX_OBS=R_EXT(1)
	  DO LS=1,NP_EXT
	    MU_AT_RMAX(LS)=SQRT(1.0_LDP-(P_EXT(LS)/R_EXT(1))**2)
	  END DO
	  IF(P_EXT(NP_EXT) .EQ. R_EXT(1))MU_AT_RMAX(NP_EXT)=0.0_LDP
	  CALL HTRPWGT(MU_AT_RMAX,HQW_AT_RMAX,NP_EXT)
C
	  FIRST_TIME=.FALSE.
C
C Find the interval in R in which R_FINE lies. This is used for interpolating
C ETA and CHI.
C
	  IF(EXTEND)THEN
	    I_START=0
	    DO WHILE(R_FINE(I_START+1) .GE. R(1))
	      I_START=I_START+1
	    END DO
	    I_START=I_START+1
	  ELSE
	    I_START=1
	  END IF
C
	  L=1
	  DO I=I_START,ND_FINE
	    DO WHILE(R_FINE(I) .LT. R(L+1))
	      L=L+1
	      IF(L .GE. ND)THEN
	        LU_ER=ERROR_LU()
	        WRITE(LU_ER,*)'Invalid interpolation range ',
	1                        'in CMF_FORM_SOL [2]'
	        STOP
	      END IF
	    END DO
	    INDX(I)=L
	  END DO
C
	ELSE
	  IF( ND_SAV .NE. ND .OR. NC_SAV .NE. NC .OR. NP_SAV .NE. NP .OR.
	1     EXTEND_SAV .NEQV. EXTEND)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error in CMF_FORM_SOL -- inconsistent dimensions'
	    WRITE(LU_ER,*)'ND_SAV=',ND_SAV,'  ND=',ND
	    WRITE(LU_ER,*)'NC_SAV=',NC_SAV,'  NC=',NC
	    WRITE(LU_ER,*)'NP_SAV=',NP_SAV,'  NP=',NP
	    WRITE(LU_ER,*)'EXTEND_SAV=',EXTEND_SAV,'  EXTEND=',EXTEND
	    STOP
	   END IF
C
	END IF			!First.
C
	IF(INIT)THEN
	  AV_PREV(:,:)=0.0_LDP
	  CV_PREV(:,:)=0.0_LDP
C
C Required for the diffusion approximation.
C
	  OLD_CHI_AT_IN_BND=CHI_NEW(ND)
C
C The old vales are actually multiplied by zero if INIT is true, but not
C initializing them may cause problems on some compilers.
C
	  ETA_OLD(1:ND)=ETA_NEW(1:ND)
	  CHI_OLD(1:ND)=CHI_NEW(1:ND)
	  ESEC_OLD(1:ND)=ESEC_NEW(1:ND)
	END IF
C
C 
C
	IPLUS_P(1:NP_OBS_MAX)=0.0_LDP
C
C Determine how many intermediate frequencies will be inserted. These
C extra frequencies are mainly inserted to assist in the propagation of the
C radiation field from grid point to gird point --- it is to improve the
C treatment of the d/dv term in the transfer equation.
C
C At present we only insert points in the profile region --- not the
C wings or resonance core.
C
	N_FREQ=1
	IF(INSERT_ADD_FREQ .AND. .NOT. INIT)THEN
       	  T1=C_KMS*(FL_PREV-FL)/FL
	  IF( T1 .GT. 1.05_LDP*FRAC_DOP*V_DOP .AND. T1 .LT. 1.05_LDP*dV_CMF_PROF)THEN
	     N_FREQ=NINT(C_KMS*(FL_PREV-FL)/(V_DOP*FRAC_DOP)/FL+0.7_LDP)
	  END IF
	END IF
	IF(.NOT. INIT)NEW_dLOG_NU=dLOG_NU/N_FREQ
C
C	CALL TUNE(1,'FRM_ML')
	DO ML=1,N_FREQ
C
C Use a linear interpolation to get opacities/emissivities at intermediate
C frequencies.
C
	  DO I=1,ND
	    CHI(I)=( CHI_OLD(I)*(N_FREQ-ML)+CHI_NEW(I)*ML )/N_FREQ
	    ETA(I)=( ETA_OLD(I)*(N_FREQ-ML)+ETA_NEW(I)*ML )/N_FREQ
	    ESEC(I)=( ESEC_OLD(I)*(N_FREQ-ML)+ESEC_NEW(I)*ML )/N_FREQ
	  END DO
C
C
C Put CHI, ETA, and dCHIdR onto the FINE radius grid.
C
C	  CALL TUNE(1,'FRM_EXP')
C
C We first do the extrapolation of CHI and ETA.
C
C We also extrapolate in ESEC, since ESEC (in the absence of negative
C absorption) provides a lower bound to the opacity. NB: When CHI is much
C larger then ESEC its variation with r dominates, and it is possible to
C extrapolate CHI below ESEC.
C
	  IF(EXTEND)THEN
	    ESEC_ALPHA=LOG(ESEC(4)/ESEC(1)) / LOG(R(1)/R(4))
	    ETA_ALPHA=LOG(ETA(4)/ETA(1))/LOG(R(1)/R(4))
	    IF(ETA_ALPHA .LT. 3.5_LDP)ETA_ALPHA=3.5
	    IF(CHI(1) .LE. ESEC(1) .OR. CHI(4) .LE. ESEC(4))THEN
	      CHI_ALPHA=ESEC_ALPHA
	    ELSE
	      CHI_ALPHA=LOG( (CHI(4)-ESEC(4)) / (CHI(1)-ESEC(1)) )
	1          / LOG(R(1)/R(4))
	      IF(CHI_ALPHA .LT. 2._LDP)CHI_ALPHA=2.0
	    END IF
	    DO I=1,I_START-1
	      T1=(CHI(1)-ESEC(1))*(R(1)/R_FINE(I))**CHI_ALPHA
	      T2=ESEC(1)*(R(1)/R_FINE(I))**ESEC_ALPHA
	      CHI_FINE(I)=T1+T2
	      dCHIdR_FINE(I)=-(T1*CHI_ALPHA+T2*ESEC_ALPHA)/R_FINE(I)
	      ETA_FINE(I)=ETA(1)*(R(1)/R_FINE(I))**ETA_ALPHA
	    END DO
	  END IF
C
C Determine the interpolation functions for CHI and ETA. Interpolations
C are done in the LOG plane.
C
	  LOG_CHI(1:ND)=LOG(CHI(1:ND))
	  LOG_ETA(1:ND)=LOG(ETA(1:ND))
	  LOG_R(1:ND)=LOG(R(1:ND))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI,LOG_R,ND)
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA,LOG_R,ND)
C
C This loop should be vectorizable.
C
	  DO I=I_START,ND_FINE
	    L=INDX(I)
	    T1=LOG(R_FINE(I)/R(L))
	    CHI_FINE(I)=EXP(
	1         ((CHI_COEF(L,1)*T1+CHI_COEF(L,2))*T1+
	1         CHI_COEF(L,3))*T1+CHI_COEF(L,4) )
	    ETA_FINE(I)=EXP(
	1         ((ETA_COEF(L,1)*T1+ETA_COEF(L,2))*T1+
	1         ETA_COEF(L,3))*T1+ETA_COEF(L,4) )
	    dCHIDR_FINE(I)=CHI_FINE(I)/R_FINE(I)*
	1         ( (3.0_LDP*CHI_COEF(L,1)*T1+2.0_LDP*CHI_COEF(L,2))*T1+
	1         CHI_COEF(L,3) )
	  END DO
C
C 
C
C Enter loop to perform integration along each ray.
C
C	  CALL TUNE(1,'FRM_LS')
	  DO LS=1,MIN(NP_OBS,ND_EXT+NC-1)
	    NI=NI_RAY(LS)
C
C Interpolate, using linear interpolation, from fine to RAY GRID.
C
	    DO I=1,NI
	      L=INDX_FINE(I,LS)
	      T1=(R_RAY(I,LS)-R_FINE(L))/(R_FINE(L+1)-R_FINE(L))
	      T2=1.0_LDP-T1
	      CHI_RAY(I)=T2*CHI_FINE(L)+T1*CHI_FINE(L+1)
	      ETA_RAY(I)=T2*ETA_FINE(L)+T1*ETA_FINE(L+1)
	      dCHIdR(I)=T2*dCHIdR_FINE(L)+T1*dCHIdR_FINE(L+1)
	    END DO
C	
	    SOURCE(1:NI)=ETA_RAY(1:NI) / CHI_RAY(1:NI)
C	    CALL TUNE(2,'FRM_EXP')
C
C Zero AV and CV vectors.
C
	    AV(:)=0.0_LDP
	    CV(:)=0.0_LDP
C
C Incident intensity at outer boundary is assumed to be zero.
C
	    IBOUND=0.0_LDP
C
C 
C
C By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
C calculation for the first frequency.
C
	    IF(INIT)THEN
	      Q(:)=0.0
	      QH(:)=0.0
	    ELSE
	      DO I=1,NI-1
	        QH(I)=GAMH(I,LS)*2.0_LDP/((CHI_RAY(I)+CHI_RAY(I+1))*NEW_dLOG_NU)
	        Q(I)=GAM(I,LS)/(CHI_RAY(I)*NEW_dLOG_NU)
	      END DO
	      QH(NI)=0.0_LDP
	      Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*NEW_dLOG_NU)
	    END IF
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P_EXT(LS)*P_EXT(LS))/
	1        R(ND)/CHI_RAY(NI)
	1        *(1.0_LDP+Q(NI)*(1.0_LDP-CHI_RAY(NI)/OLD_CHI_AT_IN_BND))
	    END IF
C
C Compute the optical depth increments. This code is from TAU, and NORDTAU.
C
	    IF(METHOD .EQ. 'ZERO')THEN
	      DO I=1,NI-1
	        DTAU(I)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*(Z(I,LS)-Z(I+1,LS))
	      END DO
	    ELSE
	       DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I)=0.5_LDP*dZ*(CHI_RAY(I)+CHI_RAY(I+1))
	          T1=dZ*dZ*( dCHIdR(I+1)*Z(I+1,LS)/R_RAY(I+1,LS)-
	1              dCHIdR(I)*Z(I,LS)/R_RAY(I,LS) )/12.0_LDP
	          T1=SIGN( MIN(0.8_LDP*DTAU(I),ABS(T1)),T1 )
	          DTAU(I)=DTAU(I)+T1
	       END DO
	    END IF
C
	    CALL TUVGHD_RH(TA,TB,TC,U,VB,VC,GB,H,XM,
	1             Q,QH,DTAU,SOURCE,DIF,DBC,IC,LS,NC,NI)
	    XM(1)=-IBOUND
C
C Update AV matrix.
C
	    AV(1)=XM(1)+U(1)*AV_PREV(1,LS)
	    DO I=2,NI-1
	        AV(I)=XM(I)+( U(I)*AV_PREV(I,LS)-
	1             (VB(I)*CV_PREV(I-1,LS)+VC(I)*CV_PREV(I,LS)) )
	    END DO
	    AV(NI)=XM(NI)+( U(NI)*AV_PREV(NI,LS)-VB(NI)*CV_PREV(NI-1,LS) )
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS_RH(TA,TB,TC,AV,NI,1)
C
C
C Check for erroneous intensities. These can be obtained because of the
C numerical differencing scheme. CV is computed with the revised AV
C (as in FG_J_CMF_V7).
C
	    DO I=1,NI
	      IF(AV(I) .LE. 0)THEN
	        AV(I)=0.1_LDP*ABS(AV(I))
	      END IF
	    END DO
C
C Update C vector (i.e. flux variable).
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV_PREV(I,LS)
	    END DO
C
C The quadrature weights used to compute HBC and NBC are now reliable.
C Note that V=AV(1)-IBOUND.
C
	    CV_BOUND=AV(1)-IBOUND
C
	    AV_PREV(1:NI,LS)=AV(1:NI)
	    CV_PREV(1:NI,LS)=CV(1:NI)
C
C Because of the frequency insertions, IPLUS_P will be continually overwritten
C until we get to the final frequency.
C
	    IPLUS_P(LS)=2.0_LDP*CV_BOUND          !To compute observed flux
C
	  END DO			!End do LS
C	  CALL TUNE(2,'FRM_LS')
	  OLD_CHI_AT_IN_BND=CHI(ND)
	END DO				!ML Loop
C	CALL TUNE(2,'FRM_ML')
C
C Save opacities/emissivities required for next frequency/
C
	CHI_OLD(1:ND)=CHI_NEW(1:ND)
	ESEC_OLD(1:ND)=ESEC_NEW(1:ND)
	ETA_OLD(1:ND)=ETA_NEW(1:ND)
	FL_PREV=FL
C
	RETURN
	END
