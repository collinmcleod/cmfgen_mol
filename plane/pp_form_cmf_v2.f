!
! Data module for PP_FORM_CMF. Data placed in this module is automatically
! saved between subroutine calls.
!
	MODULE PP_FORM_CMF_MOD_V2
	IMPLICIT NONE
!
! The *_STORE locations are used to store the radiation field, as computed
! on the previous call to PP_FORM. The *_PREV locations are used to store 
! the radiation field, as computed for the previous frequency. They are
! are updated from the *_STORE routines when NEW_FREQ is .TRUE.
!
!***************************************************************************
!***************************************************************************
!
! Data variables required by both routines.
!
! Dimensionded NRAY_MAX,NP
!
	REAL(10), ALLOCATABLE :: GAM(:,:)
	REAL(10), ALLOCATABLE :: DTAU(:,:)
	REAL(10), ALLOCATABLE :: AV(:,:)
	REAL(10), ALLOCATABLE :: CV(:,:)
!
! Dimensioned ND_EXT,4
!
	REAL(10), ALLOCATABLE :: V_COEF(:,:)
	REAL(10), ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL(10), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(10), ALLOCATABLE :: ETA_COEF(:,:)
!
! Dimensioned ND
!
	INTEGER, ALLOCATABLE :: J_PNT(:)
	INTEGER, ALLOCATABLE :: H_PNT(:)
!
! Dimensioned ND+ND_ADD=ND_EXT
!
	REAL(10), ALLOCATABLE :: R_EXT(:)
	REAL(10), ALLOCATABLE :: LOG_R_EXT(:)
	REAL(10), ALLOCATABLE :: V_EXT(:)
	REAL(10), ALLOCATABLE :: SIGMA_EXT(:)
	REAL(10), ALLOCATABLE :: ETA_EXT(:)
	REAL(10), ALLOCATABLE :: CHI_EXT(:)
	REAL(10), ALLOCATABLE :: LOG_ETA_EXT(:)
	REAL(10), ALLOCATABLE :: LOG_CHI_EXT(:)
!
! Dimensioned NRAY_MAX
!
	INTEGER, ALLOCATABLE :: RAY_PNT(:)
	REAL(10), ALLOCATABLE :: R_RAY(:)
	REAL(10), ALLOCATABLE :: dCHIdR(:)
	REAL(10), ALLOCATABLE :: Q(:)
	REAL(10), ALLOCATABLE :: QH(:)
	REAL(10), ALLOCATABLE :: V_RAY(:)
	REAL(10), ALLOCATABLE :: TMPV(:)
!
	REAL(10), ALLOCATABLE :: SIGMA_RAY(:)
	REAL(10), ALLOCATABLE :: ETA_RAY(:)
	REAL(10), ALLOCATABLE :: CHI_RAY(:)
	REAL(10), ALLOCATABLE :: dCHIdR_RAY(:)
	REAL(10), ALLOCATABLE :: SOURCE_RAY(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to the DIFFERENCE equation approach for  solving the
! transfer equation.
!	
	REAL(10), ALLOCATABLE :: AV_PREV(:,:)
	REAL(10), ALLOCATABLE :: AV_STORE(:,:)
	REAL(10), ALLOCATABLE :: CV_PREV(:,:)
	REAL(10), ALLOCATABLE :: CV_STORE(:,:)
	REAL(10), ALLOCATABLE :: PAR_AV(:,:)
!
	REAL(10), ALLOCATABLE :: GAMH(:,:)
	REAL(10), ALLOCATABLE :: TA(:,:)
	REAL(10), ALLOCATABLE :: TB(:,:)
	REAL(10), ALLOCATABLE :: TC(:,:)
	REAL(10), ALLOCATABLE :: GB(:,:)
	REAL(10), ALLOCATABLE :: H(:,:)
!
	REAL(10), ALLOCATABLE :: XM(:)
	REAL(10), ALLOCATABLE :: U(:)
	REAL(10), ALLOCATABLE :: VB(:)
	REAL(10), ALLOCATABLE :: VC(:)
	REAL(10), ALLOCATABLE :: DIV(:)
!
	REAL(10), ALLOCATABLE :: OLDCHI(:)
	REAL(10), ALLOCATABLE :: OLDCHI_STORE(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
! Dimensioned NRAY_MAX,NP
!
	REAL(10), ALLOCATABLE :: I_P_PREV(:,:)
	REAL(10), ALLOCATABLE :: I_P_STORE(:,:)
	REAL(10), ALLOCATABLE :: I_M_PREV(:,:)
	REAL(10), ALLOCATABLE :: I_M_STORE(:,:)
!
	REAL(10), ALLOCATABLE :: A0(:,:)
	REAL(10), ALLOCATABLE :: A1(:,:)
	REAL(10), ALLOCATABLE :: A2(:,:)
	REAL(10), ALLOCATABLE :: A3(:,:)
	REAL(10), ALLOCATABLE :: A4(:,:)
!
	REAL(10), ALLOCATABLE :: I_P(:,:)
	REAL(10), ALLOCATABLE :: I_M(:,:)
!
! Dimensioned NRAY_MAX
!
	REAL(10), ALLOCATABLE :: EE(:)
	REAL(10), ALLOCATABLE :: E0(:)
	REAL(10), ALLOCATABLE :: E1(:)
	REAL(10), ALLOCATABLE :: E2(:)
	REAL(10), ALLOCATABLE :: E3(:)
!
	REAL(10), ALLOCATABLE :: SOURCE_PRIME(:)
	REAL(10), ALLOCATABLE :: S(:)
	REAL(10), ALLOCATABLE :: dS(:)
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!     QW denotes quadrature weight.
!
	REAL(10), ALLOCATABLE :: MU(:)
	REAL(10), ALLOCATABLE :: JQW(:)
	REAL(10), ALLOCATABLE :: HQW(:)
	REAL(10), ALLOCATABLE :: KQW(:)
	REAL(10), ALLOCATABLE :: NQW(:)
!
	REAL(10) PREVIOUS_FREQ
	REAL(10) VDOP_FRAC_SAV
!
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NP_SAV
	INTEGER NRAY_MAX
	INTEGER NI_RAY
	INTEGER NI
	INTEGER IDMIN,IDMAX
!
	CHARACTER*10 OLD_SOLUTION_OPTIONS
	CHARACTER*10 SOLUTION_METHOD
	LOGICAL INSERT
	LOGICAL FIRST_TIME
	LOGICAL NEW_R_GRID
	LOGICAL WRITE_IP
!
	DATA VDOP_FRAC_SAV/-1001.0D0/ 	!Absurd value.
	DATA FIRST_TIME/.TRUE./
	DATA PREVIOUS_FREQ/0.0D0/
!
	END MODULE PP_FORM_CMF_MOD_V2
!
! 
!
! Routine to compute the Eddington F, G and N_ON_J Eddington factors 
! for a single frequency. The transfer is done in the comoving-frame with
! first order frequency differencing. A DIFFERENCE or INTEGRAL equation 
! approach can be used. 
!
! The DIFFERENCE equation approach is at least 30% faster. However,
! the INTEGRAL equation approach is more stable. Because monotonic cubic
! interpolation is used for the Source function in the INTEGRAL equation
! approach, the intensities are guaranteed positive (or zero).
! 
! The F, G, and N_ON_J Eddington factors are computed:
!
! NB:
!	F = K / J
!	G=  N / H
!	N_ON_J(I) = N(I)/( J(I)+ J(I+1))
!
! NB: Only one of G and N_ON_J is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' G is defined at all depths, and N_ON_J=0 at all depths.
!     IF N_TYPE='N_ON_J' N_ON_J is defined at all depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' one of G or N_ON_J is 
!       non-zero, and is the value to be used in MOM_J_CMF
!
! Routine also returns I+, so that observers flux can be computed. Note because
! the outer boundary can be thick, I+ is actually I+ - I- and hence is a flux
! like variable (This I+ = 2v at outer boundary).
!
! The radiation field at the previous (bluer) frequency is stored locally.
! Coupling to this frequency is controlled by the variable INIT. If INIT is
! true, the atmosphere is treated with the assumption that V=0.
!
! INIT must be TRUE on the very first call.
!
!
	SUBROUTINE PP_FORM_CMF_V2(ETA,CHI,ESEC,V,SIGMA,R,
	1               JNU,HNU,KNU,NNU,N_ON_J,
	1               IN_HBC,HBC,HBC_INCID,NBC,NBC_INCID,
	1               IPLUS_P,FREQ,LOG_NU,DIF,DBB,IC,
	1               VDOP_VEC,VDOP_FRAC,
	1               METHOD,SOLUTION_OPTIONS,THK,INCL_INCID_RAD,
	1               INIT,NEW_FREQ,N_TYPE,NP,ND)
	USE PP_FORM_CMF_MOD_V2
	IMPLICIT NONE
!
! Altered 14-May-2009: Altered value of IDMAX (for ETA and CHI interpolaton).
!                        IDMAX & IDMIN now stored in FG_J_CMF_MOD_V11.
! Altered 19-Jan-2009: Changed IP_DATA to IP_PP_DATA which must exist for I(p)
!                        to be output.
! Altered 08-Feb-2008 : Bug fix: NI_SMALL not set, and is no longer required as all rays have
!                                    same number of points.
! Altered 17-Feb-2007 : Bug fix. NI was not in module block.
!                          Two irrelevant lines left over from FG_CMF_V10 deleted.
! Created 01-Mar-2006 : Based on FG_J_CMF_V10
!
	INTEGER ND				!Number of depth points
	INTEGER NP				!number of anlges
	REAL(10) R(ND)
	REAL(10) V(ND)				!In km/s
	REAL(10) SIGMA(ND)			!dlnV/dlnR
	REAL(10) ETA(ND)
	REAL(10) CHI(ND)
	REAL(10) ESEC(ND)
!
	REAL(10) JNU(ND),HNU(ND),KNU(ND),NNU(ND),N_ON_J(ND)
	REAL(10) IN_HBC,HBC,HBC_INCID,NBC,NBC_INCID
	REAL(10) IPLUS_P(NP)
!
! VDOP_VEC(I) is the minimum DOPPLER width for all species at depth I. It will include
! both a turbulent, and and thermal contribution for the ionization species with the 
! highest mass. VDOP_FRAC is used to set the minimum velocity step size along a ray.
!
	REAL(10) VDOP_VEC(ND)
	REAL(10) VDOP_FRAC
!
	REAL(10) DBB,IC,FREQ,LOG_NU
	CHARACTER*(*) SOLUTION_OPTIONS
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
	LOGICAL DIF		!Use diffusion approximation
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INIT 
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INIT is set to TRUE again).
!
	LOGICAL THK
!
! First frequency -- no frequency coupling.

	LOGICAL INIT
	LOGICAL INCL_INCID_RAD
!
! Upon leaving this routine the radiation field along each ray is stored. This
! will provide the blue wing information necessary for the next frequency.
! This routine may, however, be used in an iterative loop. In this case the
! "blue wing" information should remain unaltered between calls.
! NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
! wing" information should be updated.
!
	LOGICAL NEW_FREQ	
!
	INTEGER ND_ADD_MAX
	PARAMETER (ND_ADD_MAX=24)
!
	INTEGER, SAVE :: FREQ_CNT
	INTEGER LU_IP
	INTEGER ACCESS_F
        INTEGER IOS,REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC
!
! The following arrays do not need to be stored, and hence can be created 
! each time.
!
	REAL(10) CV_BOUND(NP)		!Outer boundary V
	REAL(10) I_M_IN_BND(NP)		!Inner boundary
	REAL(10) IBOUND(NP)		!Incident intensity on outer boundary.
!
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(10) FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE,NEG_NNU_VALUE,BAD_NNU_VALUE
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL(10), PARAMETER :: ZERO=0.0D0
	REAL(10), PARAMETER :: ONE=1.0D0
	INTEGER, PARAMETER :: NINS=4
!
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!                
	INTEGER I,J,K,LS
!
	REAL(10) DBC
	REAL(10) I_CORE
	REAL(10) T1,T2
	REAL(10) ALPHA
	REAL(10) ESEC_POW
	REAL(10) BETA
	REAL(10) VINF
	REAL(10) RMAX,DEL_R_FAC
	REAL(10) dR
	REAL(10) DEL_R
!
! Change the following statement to TRUE if running on a VECTOR machine.
!
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
!
!
!
	IF(INIT .OR. NP .NE. NP_SAV)THEN
	  IF(ALLOCATED(MU))DEALLOCATE(MU,JQW,HQW,KQW,NQW)
	  ALLOCATE (MU(NP),JQW(NP),HQW(NP),KQW(NP),NQW(NP))
	  CALL GAULEG(ZERO,ONE,MU,JQW,NP)
	  HQW=JQW*MU
	  KQW=HQW*MU
	  NQW=KQW*MU
	  NP_SAV=NP
	END IF
! Check to see whether we have a new R grid, or the solution options have
! changd. This can only happen when INIT is TRUE.
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
	    IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	      NEW_R_GRID=.TRUE.
	      J=ERROR_LU()
	      WRITE(J,*)'Warning: Updating RGRID in PP_FORM_CMF_V2'
	      EXIT
	    END IF
	  END DO
	  IF(VDOP_FRAC .NE. VDOP_FRAC_SAV
	1       .OR.  OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)NEW_R_GRID=.TRUE.
	END IF
	
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF( ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT )
	  DEALLOCATE ( LOG_R_EXT )
	  DEALLOCATE ( V_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  DEALLOCATE ( LOG_ETA_EXT )
	  DEALLOCATE ( LOG_CHI_EXT )
!
	  DEALLOCATE ( RAY_PNT )
	  DEALLOCATE ( R_RAY )
	  DEALLOCATE ( J_PNT )
	  DEALLOCATE ( H_PNT )
	  DEALLOCATE ( V_COEF )
	  DEALLOCATE ( SIGMA_COEF )
	  DEALLOCATE ( ETA_COEF )
	  DEALLOCATE ( CHI_COEF )
	  DEALLOCATE ( V_RAY )
	  DEALLOCATE ( SIGMA_RAY )
	  DEALLOCATE ( ETA_RAY )
	  DEALLOCATE ( CHI_RAY )
	  DEALLOCATE ( dCHIdR_RAY )
	  DEALLOCATE ( SOURCE_RAY )
	  DEALLOCATE ( dCHIdR )
	  DEALLOCATE ( Q )
	  DEALLOCATE ( AV )
	  DEALLOCATE ( CV )
	  DEALLOCATE ( DTAU )
	  DEALLOCATE ( GAM )
!
! These arrays are only used when using the DIFFERENCE apporach.
!
	  IF( ALLOCATED(AV_PREV) )THEN
	    DEALLOCATE ( AV_PREV )
	    DEALLOCATE ( AV_STORE )
	    DEALLOCATE ( PAR_AV )
	    DEALLOCATE ( CV_PREV )
	    DEALLOCATE ( CV_STORE )
	    DEALLOCATE ( GAMH )
!
	    DEALLOCATE ( TA )
	    DEALLOCATE ( TB )
	    DEALLOCATE ( TC )
	    DEALLOCATE ( GB )
	    DEALLOCATE ( H )
	    DEALLOCATE ( OLDCHI )
	    DEALLOCATE ( OLDCHI_STORE )
!
	    DEALLOCATE ( XM )
	    DEALLOCATE ( U )
	    DEALLOCATE ( VB )
	    DEALLOCATE ( VC )
	    DEALLOCATE ( QH )
	    DEALLOCATE ( DIV )
	  END IF
!
! These arrays are only used when using the INTEGERAL apporach.
!
	  IF( ALLOCATED(I_P) )THEN
	    DEALLOCATE ( I_P )
	    DEALLOCATE ( I_M )
!
	    DEALLOCATE ( I_P_PREV )
	    DEALLOCATE ( I_P_STORE )
	    DEALLOCATE ( I_M_PREV )
	    DEALLOCATE ( I_M_STORE )
!
	    DEALLOCATE ( A0 )
	    DEALLOCATE ( A1 )
	    DEALLOCATE ( A2 )
	    DEALLOCATE ( A3 )
	    DEALLOCATE ( A4 )
!
	    DEALLOCATE (EE)
	    DEALLOCATE (E0)
	    DEALLOCATE (E1)
	    DEALLOCATE (E2)
	    DEALLOCATE (E3)
!
	    DEALLOCATE (SOURCE_PRIME)
	    DEALLOCATE (S)
	    DEALLOCATE (dS)
	  END IF
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! 
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
	  OLD_SOLUTION_OPTIONS=SOLUTION_OPTIONS
	  IF(SOLUTION_OPTIONS .EQ. 'DIFF/INS')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:4) .EQ. 'DIFF')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.FALSE.
	  ELSE IF(SOLUTION_OPTIONS .EQ. 'INT/INS')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:3) .EQ. 'INT')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.FALSE.
	  ELSE
	    J=ERROR_LU()
	    WRITE(J,*)'Error in PP_FORM_CMF_V1: Invalid solution option'
	    WRITE(J,*)SOLUTION_OPTIONS
	    STOP
	  END IF
!
          ND_ADD=0
          IF(THK)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT) )
	  ALLOCATE ( LOG_R_EXT(ND_EXT) )
	  ALLOCATE ( V_EXT(ND_EXT) )
	  ALLOCATE ( SIGMA_EXT(ND_EXT) )
	  ALLOCATE ( ETA_EXT(ND_EXT) )
	  ALLOCATE ( CHI_EXT(ND_EXT) )
	  ALLOCATE ( LOG_ETA_EXT(ND_EXT) )
	  ALLOCATE ( LOG_CHI_EXT(ND_EXT) )
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    IF(V(ND) .LT. 10.0D0 .AND. R(1)/R(ND) .GE. 9.99D0)THEN
	      RMAX=10.0D0*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0D0 .OR. V(1) .GT. 2.0D+04)THEN
	      RMAX=1.5D0*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0D0,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))                               
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-3) )
	    R_EXT(1)=RMAX
	    R_EXT(4)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.1D0*(R_EXT(1)-R_EXT(4))
	    R_EXT(3)=R_EXT(1)-0.4D0*(R_EXT(1)-R_EXT(4))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
	  END IF
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THK)THEN
	    IF(V(1) .EQ. V(5))THEN
	      V_EXT(1:ND_ADD)=V(1)
	      SIGMA_EXT(1:ND_ADD)=SIGMA(1)
	    ELSE
	      BETA=(SIGMA(1)+1.0D0)*(R(1)/R(ND)-1.0D0)
              VINF=V(1)/(1-R(ND)/R(1))**BETA
	      DO I=1,ND_ADD
	        V_EXT(I)=VINF*(1.0D0-R_EXT(ND_EXT)/R_EXT(I))**BETA
	        SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0D0)-1.0D0
	      END DO
	    END IF
	  END IF
!
! Define zone used to extrapolate opacities and emissivities.
!
	IF(ND_ADD .NE. 0)THEN
	  IDMIN=1
	  T1=R_EXT(1)/R(1)
	  IF(T1 .LT. R(1)/R(ND/3))THEN
	    T1=MIN(3.0D0,T1)
	    DO I=1,ND
	      IF(R(1)/R(I) .GT. T1)THEN
	        IDMAX=MAX(4,I)
	        EXIT
	      END IF
	    END DO
	    IDMAX=MIN(IDMAX,ND/6)
	  ELSE
	    IDMAX=ND/6
	  END IF
	  WRITE(6,'(A,I4,A)')' Using depth 1 and',IDMAX,' to extrapolate opacities'
	END IF
!
! Work out the maximim number of points per ray. This will allow us to allocate
! the required memory.
!
	  NRAY_MAX=0
	  T2=VDOP_FRAC*MINVAL(VDOP_VEC)
	  K=1
	  DO I=1,ND_EXT-1
	    T1=(V_EXT(I)-V_EXT(I+1))/T2
	    IF(T1 .GT. 1.0D0)K=K+INT(T1)
	    K=K+1
	  END DO
	  NRAY_MAX=MAX(NRAY_MAX,K)
	  ALLOCATE ( RAY_PNT(NRAY_MAX))
	  ALLOCATE ( R_RAY(NRAY_MAX) )
!
! 
!
! Define the grid along the ray. In our plane-parallel grid, this is the
! same for all rays.
!
	  K=1
	  R_RAY(K)=R_EXT(1)
	  RAY_PNT(1)=1
	  T2=VDOP_FRAC*MINVAL(VDOP_VEC)
	  DO I=1,ND_EXT-1
	    T1=(V_EXT(I)-V_EXT(I+1))/T2
	    IF(T1 .GT. 1.0D0)THEN
	      dR=(R_EXT(I+1)-R_EXT(I))/(INT(T1)+1)
	      DO J=1,INT(T1)
	        K=K+1
	        R_RAY(K)=R_RAY(K-1)+dR
	        RAY_PNT(K)=I
	      END DO
	    END IF
	    K=K+1
	    IF(K .GT. NRAY_MAX)THEN
	      J=ERROR_LU()
	      WRITE(J,*)'Error in PP_FORM_CMF_V2'
	      WRITE(J,*)K,NRAY_MAX,ND
	      STOP
	    END IF
	    R_RAY(K)=R_EXT(I+1)
	    RAY_PNT(K)=I
	  END DO
	  NI_RAY=K
	  NI=NI_RAY
!
! 
!
! J_PNT and H_PNT are used to indicate the postion of J(I) amd H(I) along the ray
! so that J and H can be computed on our regular radius grid.
!
	  ALLOCATE ( J_PNT(ND) ); J_PNT(:)=0
	  ALLOCATE ( H_PNT(ND) ); H_PNT(:)=0
!
	  K=1
	  DO I=1,ND
	    DO WHILE(J_PNT(I) .EQ. 0)
	      IF(R(I) .LE. R_RAY(K) .AND. R(I) .GE. R_RAY(K+1))THEN
	        IF( (R_RAY(K)-R(I)) .LT. (R(I)-R_RAY(K+1)) )THEN
	          J_PNT(I)=K
	        ELSE
	          J_PNT(I)=K+1
	        END IF
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! CV is defined on the nodes. Need to interpolate to the midpoint of R.
!
	    K=1
	    DO I=1,ND-1
	      DO WHILE(H_PNT(I) .EQ. 0)
	        T1=0.5D0*(R(I)+R(I+1))
	        IF(T1 .LE. R_RAY(K) .AND. T1 .GE. R_RAY(K+1))THEN
	          H_PNT(I)=K
	        ELSE
	          K=K+1
	        END IF
	      END DO
	    END DO
	  ELSE
!
! CV is defined at the mid-points, but these don't, in general, correspond to
! the midpoints of the R grid.
!
	    K=1
	    DO I=1,ND-1
	      DO WHILE(H_PNT(I) .EQ. 0)
	        T1=0.5D0*(R(I)+R(I+1))
	        IF(T1 .LE. 0.5D0*(R_RAY(K)+R_RAY(K+1)) .AND. T1 .GE. 0.5D0*(R_RAY(K+1)+R_RAY(K+2)) )THEN
	          H_PNT(I)=K
	        ELSE
	          K=K+1
	        END IF
	      END DO
	    END DO
	  END IF
!
! 
!
	  ALLOCATE ( V_COEF(ND_EXT,4) )
	  ALLOCATE ( SIGMA_COEF(ND_EXT,4) )
	  ALLOCATE ( ETA_COEF(ND_EXT,4) )
	  ALLOCATE ( CHI_COEF(ND_EXT,4) )
!
!***************************************************************************
!***************************************************************************
!
	  ALLOCATE ( V_RAY(NRAY_MAX) )
	  ALLOCATE ( TMPV(NRAY_MAX) )
	  ALLOCATE ( SIGMA_RAY(NRAY_MAX) )
	  ALLOCATE ( ETA_RAY(NRAY_MAX) )
	  ALLOCATE ( CHI_RAY(NRAY_MAX) )
	  ALLOCATE ( dCHIdR_RAY(NRAY_MAX) )
	  ALLOCATE ( SOURCE_RAY(NRAY_MAX) )
	  ALLOCATE ( dCHIdR(NRAY_MAX) )
	  ALLOCATE ( Q(NRAY_MAX))
!
	  ALLOCATE ( AV(NRAY_MAX,NP) )
	  ALLOCATE ( CV(NRAY_MAX,NP) )
	  ALLOCATE ( DTAU(NRAY_MAX,NP) )
	  ALLOCATE ( GAM(NRAY_MAX,NP) )
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by both DIFFERENCE solution approach.
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    ALLOCATE ( AV_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( AV_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( PAR_AV(NRAY_MAX,NP) )
	    ALLOCATE ( CV_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( CV_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( GAMH(NRAY_MAX,NP) )
!
	    ALLOCATE ( TA(NRAY_MAX,NP) )
	    ALLOCATE ( TB(NRAY_MAX,NP) )
	    ALLOCATE ( TC(NRAY_MAX,NP) )
	    ALLOCATE ( GB(NRAY_MAX,NP) )
	    ALLOCATE ( H(NRAY_MAX,NP) )
	    ALLOCATE ( OLDCHI(NP) )
	    ALLOCATE ( OLDCHI_STORE(NP) )
!
	    ALLOCATE ( XM(NRAY_MAX) )
	    ALLOCATE ( U(NRAY_MAX) )
	    ALLOCATE ( VB(NRAY_MAX) )
	    ALLOCATE ( VC(NRAY_MAX) )
	    ALLOCATE ( QH(NRAY_MAX) )
	    ALLOCATE ( DIV(NRAY_MAX) )
	  END IF
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by INTEGRAL solution.
!
	  IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
	    ALLOCATE (I_P(NRAY_MAX,NP))
	    ALLOCATE (I_M(NRAY_MAX,NP))
!
	    ALLOCATE ( I_P_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( I_P_STORE(NRAY_MAX,NP) )
	    ALLOCATE ( I_M_PREV(NRAY_MAX,NP) )
	    ALLOCATE ( I_M_STORE(NRAY_MAX,NP) )
!
	    ALLOCATE ( A0(NRAY_MAX,NP) )
	    ALLOCATE ( A1(NRAY_MAX,NP) )
	    ALLOCATE ( A2(NRAY_MAX,NP) )
	    ALLOCATE ( A3(NRAY_MAX,NP) )
	    ALLOCATE ( A4(NRAY_MAX,NP) )
!
	    ALLOCATE (EE(NRAY_MAX))
	    ALLOCATE (E0(NRAY_MAX))
	    ALLOCATE (E1(NRAY_MAX))
	    ALLOCATE (E2(NRAY_MAX))
	    ALLOCATE (E3(NRAY_MAX))
!
	    ALLOCATE (SOURCE_PRIME(NRAY_MAX))
	    ALLOCATE (S(NRAY_MAX))
	    ALLOCATE (dS(NRAY_MAX))
!
	  END IF
!
	  FIRST_TIME=.FALSE.
	ELSE
	  IF(OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in PP_FORM_CMF_V1'
	    WRITE(J,*)'Can''t switch SOLUTION_OPTIONS while runing code'
	    WRITE(J,*)'New setting:',SOLUTION_OPTIONS
	    WRITE(J,*)'Old setting:',OLD_SOLUTION_OPTIONS
	    STOP
	  END IF
	END IF
!
!
	NEG_AV_VALUE=.FALSE.
	NEG_NNU_VALUE=.FALSE.
	BAD_NNU_VALUE=.FALSE.
!         
! Perform initializations. 
!
	IF(INIT)THEN
!
	  FG_ERR_ON_FREQ(:)=0.0D0
	  FG_ERR_TYPE(:)=0
	  FG _ERR_CNT=0
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    AV_PREV(:,:)=0.0D0
	    CV_PREV(:,:)=0.0D0
          ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=0.0D0
	    I_M_PREV(:,:)=0.0D0
	  END IF
!
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(V_COEF,V_EXT,LOG_R_EXT,ND_EXT)
	  CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_EXT,LOG_R_EXT,ND_EXT)
!
	  DO I=1,NI_RAY
	    K=RAY_PNT(I)
	    T1=LOG(R_RAY(I)/R_EXT(K))
	    V_RAY(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
	    SIGMA_RAY(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
	  END DO
!
! Compute GAMMA. This section is straight from the subroutine GAMMA, except
! That _EXT has been added to V, SIGMA, and R.
!
! We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
!  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  DO LS=1,NP
	    DO I=1,NI_RAY
	      T1=3.33564D-06*V_RAY(I)/R_RAY(I)
	      GAM(I,LS)=T1*MU(LS)*MU(LS)*(SIGMA_RAY(I)+1.0D0)
	    END DO
	  END DO		!LS Loop
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    DO LS=1,NP
	      DO I=1,NI_RAY-1
	        GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564D-06*
	1                   0.5D0*MU(LS)*MU(LS)*
	1                   (SIGMA_RAY(I)+SIGMA_RAY(I+1)+2.0D0)/
	1                   (R_RAY(I)+R_RAY(I+1))
	      END DO
	      GAMH(NI,LS)=0.0D0
	    END DO		!LS Loop
	  END IF
	ELSE IF(NEW_FREQ)THEN
          IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    OLDCHI(1:NP)=OLDCHI_STORE(1:NP)
	    AV_PREV(:,:)=AV_STORE(:,:)
	    CV_PREV(:,:)=CV_STORE(:,:)
!	  ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL' .AND. I_P_PREV(1,1) .EQ. 0.0D0)THEN
	  ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=I_P_STORE(:,:)
	    I_M_PREV(:,:)=I_M_STORE(:,:)
	  END IF
	  IF(FREQ .GE. PREVIOUS_FREQ)THEN
	    WRITE(ERROR_LU(),*)'Error in PP_FORM_CMF_V1'
	    WRITE(ERROR_LU(),*)'Frequencies must be monotonically decreasng'
	    STOP
	  END IF
	ELSE IF(FREQ .NE. PREVIOUS_FREQ)THEN
	  WRITE(ERROR_LU(),*)'Error in PP_FORM_CMF_V1'
	  WRITE(ERROR_LU(),*)'Frequencies must not change if NEW_FREQ=.FALSE.'
	  STOP
	END IF
	PREVIOUS_FREQ=FREQ
!
! 
!
	CALL TUNE(1,'FG_BIG_SEC')
!
! Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
! during the iteration procedure. ETA does however change (since it depends
! of J) and thus ETA_EXT must be re-computed on each entry.
!
! The first checks whether we may have negative line opacities due
! to stimulated emission. In such a case we simply assume an 1/r^2
! extrapolation.
!
! We also interpolate in ESEC, since ESEC (in the absence of negative 
! absorption) provides a lower bound to the opacity. NB: When CHI is much
! larger then ESEC its variation with r dominates, and it is possible to
! extrapolate CHI below ESEC.
!
	CALL TUNE(1,'FG_CHI_BEG')
	IF(ND_ADD .NE. 0)THEN
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0D0)ESEC_POW=2.0D0
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.0D0)ALPHA=2.0D0
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0D0)ESEC_POW=2.0D0
   	    DO I=1,ND_ADD
	      T1=(CHI(IDMIN)-ESEC(IDMIN))*(R(IDMIN)/R_EXT(I))**ALPHA
	      T2=ESEC(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	      dCHIdR(I)=(-ALPHA*T1-ESEC_POW*T2)/R_EXT(I)
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
!
! We limit alpha to 3.5 to avoid excess envelope emission. If alpha were
! 3 we would get a logarithmic flux divergence as we increase the volume.
!
	  ALPHA=LOG(ETA(IDMAX)/ETA(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	  IF(ALPHA .LT. 3.5D0)ALPHA=3.5D0
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(IDMIN)*(R(IDMIN)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0D-280)ETA_EXT(I)=1.0D-280
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
        CALL DERIVCHI(dCHIdr,CHI_EXT,R_EXT,ND_EXT,METHOD)
!
	IF(NI_RAY .NE. ND_EXT)THEN
	  LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	  LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	END IF
	CALL TUNE(2,'FG_CHI_BEG')
!
	IF(NI_RAY .EQ. ND_EXT)THEN
	  CHI_RAY(1:ND_EXT)=CHI_EXT(1:ND_EXT)
	  ETA_RAY(1:ND_EXT)=ETA_EXT(1:ND_EXT)
	  dCHIdR_RAY(1:ND_EXT)=dCHIdR(1:ND_EXT)
	ELSE
	  DO I=1,NI
	    K=RAY_PNT(I)
	    IF(R_RAY(I) .EQ. R_EXT(K))THEN
	      CHI_RAY(I)=CHI_EXT(K)
	      ETA_RAY(I)=ETA_EXT(K)
	      dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I)
	    ELSE
	      T1=LOG(R_RAY(I)/R_EXT(K))
	      T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	      CHI_RAY(I)=EXP(T2)
	      T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	      ETA_RAY(I)=EXP(T2)
	      T2=(3.0D0*CHI_COEF(K,1)*T1+2.0D0*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	      dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I)
	    END IF
	  END DO
	END IF 
	SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!
	IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is not necessary when
! we are using monotonic interpolation. 
!
	  NI=NI_RAY
	  IF(METHOD .EQ. 'ZERO')THEN
	    DO I=1,NI-1
	      dR=R_RAY(I)-R_RAY(I+1)
	      TMPV(I)=0.5D0*(CHI_RAY(I)+CHI_RAY(I+1))*dR
	    END DO  
	  ELSE IF(INSERT .OR. METHOD(4:6) .EQ. 'MON')THEN
	    DO I=1,NI-1
	      dR=R_RAY(I)-R_RAY(I+1)
	      TMPV(I)=0.5D0*dR*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dR*(dCHIdR_RAY(I+1) - dCHIdR_RAY(I))/6.0D0 )
	    END DO
	  ELSE
	    DO I=1,NI-1
	      dR=R_RAY(I)-R_RAY(I+1)
	      TMPV(I)=0.5D0*dR*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dR*(dCHIdR_RAY(I+1)-dCHIdR_RAY(I))/6.0D0 )
	      IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	        TMPV(I)=MAX(CHI_RAY(I)*dR,TMPV(I))
	        TMPV(I)=MIN(CHI_RAY(I+1)*dR,TMPV(I))
	      ELSE
	        TMPV(I)=MIN(CHI_RAY(I)*dR,TMPV(I))
	        TMPV(I)=MAX(CHI_RAY(I+1)*dR,TMPV(I))
	      END IF
	    END DO
	  END IF
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE' .OR. INIT)THEN
	    DO LS=1,NP
	      DO I=1,NI-1
	        DTAU(I,LS)=TMPV(I)/MU(LS)
	      END DO
	    END DO
	  ELSE
	    DO LS=1,NP
	      DO I=1,NI-1
	        DTAU(I,LS)=TMPV(I)/MU(LS)+ 3.33564D-06*MU(LS)*(V_RAY(I)-V_RAY(I+1))/LOG_NU
	      END DO
	    END DO
	  END IF
!
	END IF
!
! Even in the THK case we assume that the intensity incident on the
! outer boundary is zero. In THK case we effectively get IBOUND not
! equal zero at the model atmosphere boundary as a consequence of the
! extension.
!
	IF(INCL_INCID_RAD)THEN
	  CALL GET_IBOUND(IBOUND,MU,NP,FREQ,INIT)
	ELSE
	  IBOUND(:)=0.0D0
	END IF
!
! 
!***************************************************************************
!***************************************************************************
!
	IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
!
! Zero AV and CV matrices.
!
	  AV(:,:)=0.0D0
	  CV(:,:)=0.0D0
!
!	  CALL TUNE(1,'LS_LOOP')
!
	  DO LS=1,NP
	    NI=NI_RAY
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum 
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      DO I=1,NI
	        Q(I)=0.0D0
	        QH(I)=0.0D0
	      END DO                
	      OLDCHI(LS)=CHI_RAY(NI)
	    ELSE IF(NEW_FREQ)THEN
	      DO I=1,NI-1
	        QH(I)=GAMH(I,LS)*2.0D0/((CHI_RAY(I)+CHI_RAY(I+1))*LOG_NU)
	        Q(I)=GAM(I,LS)/(CHI_RAY(I)*LOG_NU)
	      END DO
	      QH(NI)=0.0D0
	      Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*LOG_NU)
	    END IF
!
	    IF(DIF)THEN
	      T1=0.0D0
	      IF(.NOT. INIT)T1=GAM(NI,LS)/(CHI_RAY(NI)*LOG_NU)	     !Q(NI)
	      DBC=DBB*MU(LS)/CHI_RAY(NI)
	1            *(1.0D0+T1*(1.0D0-CHI_RAY(NI)/OLDCHI(LS)))
	    END IF
	    OLDCHI_STORE(LS)=CHI_RAY(NI)
!
	    IF(NEW_FREQ)THEN
!
	      CALL TUVGHD_RH(TA(1,LS),TB(1,LS),TC(1,LS),
	1            U,VB,VC,GB(1,LS),H(1,LS),XM,
	1            Q,QH,DTAU(1,LS),SOURCE_RAY,DIF,DBC,IC,LS,NP,NI)
	      XM(1)=-IBOUND(LS)
!
	      DO I=1,ND
	        WRITE(35,'(5ES16.8)')TA(I,LS),TB(I,LS),TC(I,LS),XM(I),DTAU(I,LS)
	      END DO
!
! Update PAR_AV matrix. U, VC, VC, AV_PREV, and CV_PREV do not change on
! successive calls with the same frequency. Thus PAR_AV can be used on
! successive calls with the same frequency.
!
	      PAR_AV(1,LS)=U(1)*AV_PREV(1,LS)
	      DO I=2,NI-1
	         PAR_AV(I,LS)=U(I)*AV_PREV(I,LS)-
	1             (VB(I)*CV_PREV(I-1,LS)+VC(I)*CV_PREV(I,LS)) 
	      END DO
	      PAR_AV(NI,LS)=U(NI)*AV_PREV(NI,LS)-VB(NI)*CV_PREV(NI-1,LS)
!
	      TC(NI,LS)=0.0D0				!As used.
	      DIV(1)=1.0D0/(TC(1,LS)+TB(1,LS))
	      TC(1,LS)=TC(1,LS)*DIV(1)
	      TB(1,LS)=TB(1,LS)*DIV(1)
	      DO I=2,NI
	        DIV(I)=1.0D0/(TA(I,LS)*TB(I-1,LS)+TB(I,LS)+TC(I,LS))
	        TB(I,LS)=(TA(I,LS)*TB(I-1,LS)+TB(I,LS))*DIV(I)
	        TC(I,LS)=TC(I,LS)*DIV(I)
	      END DO
	      TB(1:NI,LS)=DIV(1:NI)
	    ELSE
!
	      XM(1)=0.0D0
	      DO I=2,NI-1
	        XM(I)=-0.5D0*SOURCE_RAY(I)*(DTAU(I-1,LS)+DTAU(I,LS))
	      END DO
!
	      IF(DIF)THEN
	        XM(NI)=DBC
	      ELSE
	        XM(NI)=IC
	     END IF
!
	    END IF
!
! Update AV matrix.
!
	    DO I=1,NI
	      AV(I,LS)=XM(I)+PAR_AV(I,LS)
	    END DO
!
	  END DO				!LS
!
! 
! 
! Solve for the radiation field along ray for this frequency.
!
	  IF(VECTOR_MACHINE)THEN
!
! This section is for a VECTOR machine. Loop over inner index is outer
! loop to remove a dependency. This is inefficient on scaler machines
! as array is not accessed sequentially.
!
!
! Forward substitution.
!
	    AV(1,1:NP)=AV(1,1:NP)*TB(1,1:NP)
	    DO I=2,NI_RAY
	      DO LS=1,NP
	        AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY,LS)=-AV(NI_RAY,LS)
	    END DO
	    DO I=NI_RAY,1,-1
	      DO LS=1,NP
	          AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  ELSE
C
C This section of code is for a scaler machine where the dependence
C on a previous computation does not matter. More efficient than previous
C code as array is accessed in correct manner.
C
C Forward substitution.
C
	    DO LS=1,NP
	      AV(1,LS)=AV(1,LS)*TB(1,LS)
	      DO I=2,NI_RAY
	        AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY,LS)=-AV(NI_RAY,LS)
	      DO I=NI_RAY-1,1,-1
	        AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  END IF
!	  CALL TUNE(2,'LS_LOOP')
!
!
!
! Verify validity of AV values. If these go negative, we set them +ve 
! but a factor of 10 smaller. We also ensure that the AV are not 
! extremely close to zero by comparing with neighboring AV values.
! The CV (fluxes) are computed using the revised values.
!
	  DO LS=1,NP
	   NI=NI_RAY
!
	    K=2
	    DO I=1,NI
	      T1=1.0D-08*ABS(AV(K,LS))
	      IF(AV(I,LS) .LE. T1)THEN
	        AV(I,LS)=MAX(0.01D0*ABS(AV(I,LS)),T1)
	        NEG_AV_VALUE=.TRUE.
	      ELSE
	        K=K+1
	      END IF
	      IF(I .EQ. 1)K=1
	    END DO
!
! Update C vector (i.e. flux variable).
!                          
	    DO I=1,NI_RAY-1
	      CV(I,LS)=GB(I,LS)*(AV(I,LS)-AV(I+1,LS))+H(I,LS)*CV_PREV(I,LS)
	    END DO
!
! Note that V=AV(1)-IBOUND.
!
	    IF(ND_ADD .EQ. 0)THEN
	      CV_BOUND(LS)=AV(1,LS)-IBOUND(LS)
	      IPLUS_P(LS)=AV(K,LS)+CV_BOUND(LS)
	    ELSE IF(LS .EQ. NP)THEN
	      CV_BOUND(LS)=0.0D0
	      IPLUS_P(LS)=0.0D0
	    ELSE
	      K=J_PNT(1)
	      CV_BOUND(LS)=0.5D0*(CV(K+1,LS)+CV(K-1,LS))
	      IPLUS_P(LS)=AV(K,LS)+CV_BOUND(LS)
	    END IF
	    I_M_IN_BND(LS)=AV(NI,LS)-(AV(NI,LS)-AV(NI-1,LS))/DTAU(NI-1,LS)
	  END DO
!
	  AV_STORE(:,:)=AV(:,:)
	  CV_STORE(:,:)=CV(:,:)
!
	END IF			!Solution method
!
!************************************************************************
!************************************************************************
!
	CALL TUNE(1,'FG_SOL_INT')
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! Enter loop to perform integration along each ray.
!
	  NI=NI_RAY
	  DO LS=1,NP
!
	    IF(DIF)THEN
              I_CORE=(ETA_RAY(NI)+ MU(LS)*DBB)/CHI_RAY(NI)
	    ELSE
	      I_CORE=IC
	    END IF
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum 
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      Q(1:NI)=0.0D0
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	    ELSE
	      Q(1:NI)=GAM(1:NI,LS)/LOG_NU
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/(Q(1:NI)+CHI_RAY(1:NI))
	      Q(1:NI)=Q(1:NI)/(Q(1:NI)+CHI_RAY(1:NI))
	    END IF
!	    CALL TUNE(2,'FG_CHI_Q')
!
!
!	    CALL TUNE(1,'FG_NEW_F')
	    IF(NEW_FREQ)THEN
!                       
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	      DO I=1,NI-1
	        T1=DTAU(I,LS)
	        EE(I)=0.0D0
	        IF(T1 .LT. 700.0D0)EE(I)=EXP(-T1)
	        IF(T1 .GT. 0.5D0)THEN
	          E0(I)=1.0D0-EE(I)
	          E1(I)=1.0D0-E0(I)/T1
	          E2(I)=1.0D0-2.0D0*E1(I)/T1
	          E3(I)=1.0D0-3.0D0*E2(I)/T1
	        ELSE IF(T1 .GT. 0.1D0)THEN
	          E3(I)=0.25D0*T1*( 1.0D0-0.20*T1*
	1               (1.0D0-T1/6.0D0*(1.0D0-T1/7.0D0*
	1               (1.0D0-T1/8.0D0*(1.0D0-T1/9.0D0*
	1               (1.0D0-T1/10.0D0*(1.0D0-T1/11.0D0*
	1               (1.0D0-T1/12.0D0*(1.0D0-T1/13.0D0)))))))) )
	          E2(I)=T1*( 1.0D0-E3(I) )/3.0D0
	          E1(I)=T1*( 1.0D0-E2(I) )/2.0D0
	          E0(I)=T1*( 1.0D0-E1(I) )
	        ELSE
	          E3(I)=0.25D0*T1*( 1.0D0-0.20*T1*
	1               (1.0D0-T1/6.0D0*(1.0D0-T1/7.0D0*
	1               (1.0D0-T1/8.0D0*(1.0D0-T1/9.0D0) ))))
	          E2(I)=T1*( 1.0D0-E3(I) )/3.0D0
	          E1(I)=T1*( 1.0D0-E2(I) )/2.0d0
	          E0(I)=T1*( 1.0D0-E1(I) )
	        END IF
	      END DO
!
              DO I=1,NI-1
	        A0(I,LS)=EE(I)
	        A1(I,LS)=E0(I)-3.0D0*E2(I)+2.0D0*E3(I)
	        A2(I,LS)=3.0D0*E2(I)-2.0D0*E3(I)
	        A3(I,LS)=DTAU(I,LS)*(E1(I)-2.0D0*E2(I)+E3(I))
	        A4(I,LS)=DTAU(I,LS)*(E3(I)-E2(I))
	      END DO
	    END IF
!	    CALL TUNE(2,'FG_NEW_F')
!
! ******************* INWARD DIRECTED RAYS *********************************
!
! Compute the Source function for inward directed rays, and find the
! monotonic interpolating polynomial.
!
!	    CALL TUNE(1,'FG_INT_INT')
!
	    SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*I_M_PREV(1:NI,LS)
	    DO I=1,NI-1
	      S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	    END DO
!
! Now compute the derivatives node I. 
!
	    dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/(DTAU(1,LS)+DTAU(2,LS))
	    DO I=2,NI-1
	      dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/
	1                       (DTAU(I-1,LS)+DTAU(I,LS))
	    END DO
	    dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/
	1                       (DTAU(NI-2,LS)+DTAU(NI-1,LS))
!             
! Adjust first derivatives so that function is monotonic  in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5D0*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1               MIN(ABS(S(I-1)),ABS(S(I)),0.5D0*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1               MIN(ABS(S(NI-1)),0.5D0*ABS(dS(NI)))
!
            I_M(1,LS)=IBOUND(LS)
	    DO I=1,NI-1
	      I_M(I+1,LS)=I_M(I,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I)*A1(I,LS)
	1        +    SOURCE_PRIME(I+1)*A2(I,LS)
	1        +    dS(I)*A3(I,LS)
	1        +    dS(I+1)*A4(I,LS) )
	    END DO
!
! ******************* OUTWARD DIRECTED RAYS *********************************
!
! Compute the Source function for outward directed rays, and find the
! monotonic interpolating polynomial.
!
	    SOURCE_PRIME(1:NI)=SOURCE_RAY(1:NI)+Q(1:NI)*I_P_PREV(1:NI,LS)
	    DO I=1,NI-1
	      S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	    END DO
!
! Now compute the derivatives at node I. 
!
	    dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/(DTAU(1,LS)+DTAU(2,LS))
	    DO I=2,NI-1
	      dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/
	1                  (DTAU(I-1,LS)+DTAU(I,LS))
	    END DO
	    dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/
	1                  (DTAU(NI-2,LS)+DTAU(NI-1,LS))
!             
! Adjust the first derivatives so that function is monotonic in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5D0*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5D0*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1            MIN(ABS(S(NI-1)),0.5D0*ABS(dS(NI)))
!
	    I_P(NI,LS)=I_CORE
	    DO I=NI-1,1,-1
	      I_P(I,LS)=I_P(I+1,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I+1)*A1(I,LS)
	1        +    SOURCE_PRIME(I)*A2(I,LS)
	1        -    dS(I+1)*A3(I,LS)
	1        -    dS(I)*A4(I,LS) )
	    END DO
C
C Note that V=AV(1)-IBOUND.
C
1000	    CONTINUE
	    K=J_PNT(1)
	    CV_BOUND(LS)=0.5D0*(I_P(K,LS)-I_M(K,LS))
	    IPLUS_P(LS)=I_P(K,LS)
	    I_M_IN_BND(LS)=I_M(NI,LS)
!
!	    CALL TUNE(2,'FG_INT_INT')
	  END DO			!LS
!
! Compute the mean intensity like variable U at each grid point. Used to
! compute J.
!
	  DO LS=1,NP
	    DO I=1,NI_RAY
	      AV(I,LS)=0.5D0*( I_P(I,LS)+I_M(I,LS) )
	    END DO
	  END DO
!
! We define CV on the gregular grid. We perform the interpolation onto
! the original R grid as we compute H.
!
	  DO LS=1,NP
	    DO I=1,NI_RAY
	      CV(I,LS)=0.5D0*( I_P(I,LS)-I_M(I,LS) )
	    END DO
	  END DO
!
	  I_P_STORE(:,:)=I_P(:,:)
	  I_M_STORE(:,:)=I_M(:,:)
!
	END IF				!Solution method
	CALL TUNE(2,'FG_SOL_INT')
!
!
!***************************************************************************
!***************************************************************************
!
!
! Zero boundary conditions.
!
	HBC=0.0D0			!uI+/J or H/J at model outer boundary.
	HBC_INCID=0.0D0			!iI-/J at model outer boundary.
	NBC=0.0D0			!N/J at model outer boundary.
	NBC_INCID=0.0D0
	IN_HBC=0.0D0
!
! Initialize intensity matrices.
!
	JNU(:)=0.0D0			!1:ND
	HNU(:)=0.0D0
	KNU(:)=0.0D0
	NNU(:)=0.0D0
	N_ON_J(:)=0.0D0
!
	CALL TUNE(1,'J_LOOP')
	DO LS=1,NP
	  DO I=1,ND
	    K=J_PNT(I)
	    JNU(I)=JNU(I)+JQW(LS)*AV(K,LS)
	    KNU(I)=KNU(I)+KQW(LS)*AV(K,LS)
	  END DO
	END DO
!
! Perform the interpolation of CV onto the mid points of the R grid, and
! compute H & N on this grid. With the INTEGRAL method, we will always
! need to interpolate, since, V is always defined on the grid points.
!
! For the DIFFERNECE method, we will only need to interpolate when we have
! included additional points along each ray.
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!	  Q(1)=0.0D0
!	  DO I=2,ND
!	     Q(I)=Q(I-1)+TMPV(I-1)
!	  END DO
!	  DO I=1,ND-1
!	    TMPV(I)=0.5D0*(Q(I)+Q(I+1))
!	  END DO
!	  DO LS=1,NP
!	    S(1:NI_RAY)=CV(1:NI_RAY,LS); 
!	    CALL MON_INTERP(dS,ND-1,1,TMPV,ND-1,S,1,Q,NI_RAY)
!	    DO I=1,ND-1
!	      HNU(I)=HNU(I)+HQW(LS)*dS(I)
!	      NNU(I)=NNU(I)+NQW(LS)*dS(I)
!	    END DO
!	  END DO
	  DO I=1,ND-1
	     TMPV(I)=0.5D0*(R(I)+R(I+1))
	  END DO
	  DO LS=1,NP
	    S(1:NI_RAY)=CV(1:NI_RAY,LS); 
	    CALL MON_INTERP(dS,ND-1,1,TMPV,ND-1,S,1,R_RAY,NI_RAY)
	    DO I=1,ND-1
	      HNU(I)=HNU(I)+HQW(LS)*dS(I)
	      NNU(I)=NNU(I)+NQW(LS)*dS(I)
	    END DO
	  END DO
!	  DO LS=1,NP
!	    DO I=1,ND-1
!	      K=H_PNT(I)
!	      T2=0.5D0*(R(I)+R(I+1))
!	      T1=(T2-R_RAY(K))/(R_RAY(K+1)-R_RAY(K))
!	      HNU(I)=HNU(I)+HQW(LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
!	      NNU(I)=NNU(I)+NQW(LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
!	    END DO
!	  END DO			!End do LS
	ELSE
	  DO LS=1,NP
	    DO I=1,ND-1
	      K=H_PNT(I)
	      T2=0.5D0*(R(I)+R(I+1))
	      T1=( 2.0D0*T2-(R_RAY(K)+R_RAY(K+1)) )/(R_RAY(K+2)-R_RAY(K))
	      HNU(I)=HNU(I)+HQW(LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	      NNU(I)=NNU(I)+NQW(LS)*((1.0D0-T1)*CV(K,LS)+T1*CV(K+1,LS))
	    END DO
	  END DO			!End do LS
	END IF
!
	DO LS=1,NP
!	  IPLUS_P(LS)=2.0D0*CV_BOUND(LS)
	  IN_HBC=IN_HBC + HQW(LS)*I_M_IN_BND(LS)
	END DO
!
! NB: J=Int 0.5(I+ + I-) & H = Int mu 0.5(I+ - I-). With the factor of 0.5 included,
! H = HBC . J - HBC_INCID. In the second case we simply have H = HBC.J.
!
	IF(INCL_INCID_RAD)THEN
	  K=J_PNT(1)
	  DO LS=1,NP
	    HBC=HBC+0.5D0*(AV(K,LS)+CV_BOUND(LS))*HQW(LS)
	    HBC_INCID=HBC_INCID+0.5D0*(AV(K,LS)-CV_BOUND(LS))*HQW(LS)
	    NBC=NBC+0.5D0*(AV(K,LS)+CV_BOUND(LS))*NQW(LS)
	    NBC_INCID=NBC_INCID+0.5D0*(AV(K,LS)-CV_BOUND(LS))*NQW(LS)
	  END DO
	ELSE
	  DO LS=1,NP
	    HBC=HBC+CV_BOUND(LS)*HQW(LS)
	    NBC=NBC+CV_BOUND(LS)*NQW(LS)
	  END DO
	END IF
!
	CALL TUNE(2,'J_LOOP')
!
!
!
! Compute the boundary Eddington factors.
!
	HBC=HBC/JNU(1)
	NBC=NBC/JNU(1)
	IN_HBC=IN_HBC/(2.0D0*JNU(ND)-IC)
!
! Compute the Eddington F and G factors, which are defined by
! F=K/J and G=N/H. The F and G factors are returned in KNU and NNU
! respectively.
!
	DO I=1,ND
	  KNU(I)=KNU(I)/JNU(I)
	END DO
!
! Now compute the "G" Eddington factors. Because H may be zero we have to
! be careful. Depending on N_TYPE we can specify N in terms of J, N in
! terms of H, or N in terms of J and H.
!
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  DO I=1,ND-1
	    N_ON_J(I)=NNU(I)/(JNU(I)+JNU(I+1))
	    NNU(I)=0.0D0
	  END DO
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
	  N_ON_J(1:ND-1)=0.0D0
	  NNU(1)=NNU(1)/HNU(1)
	  DO I=2,ND-1
	    IF(HNU(I) .NE. 0.0D0)THEN
	      T1=NNU(I)/HNU(I)
	    ELSE
	      T1=100.0D0
	    END IF
	    IF(T1 .GT. 1.1D0 .OR. T1 .LT. 0.05D0)THEN
	      N_ON_J(I)=NNU(I)/(JNU(I)+JNU(I+1))
	      NNU(I)=0.0D0
	    ELSE
	      NNU(I)=T1
	    END IF
	  END DO
	ELSE IF(N_TYPE .EQ. 'G_ONLY')THEN
!
! Compute G Eddington factor storing in N. We also check the validity of 
! the Eddington factor in case strange values are occurring because H 
! is near zero (switching sign?).
!
	  DO I=1,ND-1
	    IF(HNU(I) .NE. 0)THEN
	      NNU(I)=NNU(I)/HNU(I)
	      IF(NNU(I) .GT. 1.0D0)THEN
                BAD_NNU_VALUE=.TRUE.
	        NNU(I)=1.0D0
	      ELSE IF( NNU(I) .LT. 0.01D0) THEN
	        NEG_NNU_VALUE=.TRUE.
	        NNU(I)=0.01D0
	      END IF
	    ELSE
	      NNU(I)=0.0D0		!Later replaced by average.
	      J=ERROR_LU()
	      WRITE(J,'(1X,A,1PE16.8)')'HNU zero for frequency:',FREQ
	      WRITE(J,'(1X,A,I4)')'Error occurred at depth:',I
	    END IF
	  END DO
!
	ELSE
	  J=ERROR_LU()
	  WRITE(J,*)'Unrecognized N_TYPE in PP_FORM_CMF_V1'
	  WRITE(J,*)'N_TYPE=',N_TYPE
	  STOP
	END IF
!
! Store frequencies at which errors occurred, and give an indication of the
! error.
!
	J=1
	DO WHILE (J .LE. N_ERR_MAX .AND. 
	1            (NEG_AV_VALUE .OR. BAD_NNU_VALUE .OR. NEG_NNU_VALUE) )
	  IF(FG_ERR_ON_FREQ(J) .EQ. FREQ)THEN
	    FG_ERR_TYPE(J)=0
	    IF(NEG_AV_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    IF(NEG_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+5
	    IF(BAD_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+10
	    NEG_AV_VALUE=.FALSE.
	    BAD_NNU_VALUE=.FALSE.
	    NEG_NNU_VALUE=.FALSE.
	  ELSE IF(J .EQ. FG_ERR_CNT+1)THEN
	    FG_ERR_CNT=J
	    FG_ERR_TYPE(J)=0
	    FG_ERR_ON_FREQ(J)=FREQ
	    IF(NEG_AV_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    IF(NEG_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+5
	    IF(BAD_NNU_VALUE)FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+10
	    NEG_AV_VALUE=.FALSE.
	    BAD_NNU_VALUE=.FALSE.
	    NEG_NNU_VALUE=.FALSE.
	  END IF
	  J=J+1
	END DO
!
	DO I=1,ND
	  IF(JNU(I) .LT. 0.0D0)THEN
	    WRITE(111,*)'FREQ=',FREQ
	    WRITE(111,'(3X,A,5X,4(5X,A,5X))')'I','JNU','ETA','CHI','ESEC'
	    DO J=1,ND
	      WRITE(111,'(X,I5,4ES16.6)')J,JNU(J),ETA(J),CHI(J),ESEC(J)
	    END DO
	    J=ERROR_LU()
	    WRITE(J,*)'Error on PP_FORM_CMF_V1--- negative mean intensities.'
	    WRITE(J,*)'Check out file FOR111 for aditional information.'
	    WRITE(J,*)'Halting code execution.'
	    STOP
	  END IF
	END DO
!
        ACCESS_F=5
	LU_IP=56
	IF(INIT)THEN
	  INQUIRE(FILE='IP_PP_DATA',EXIST=WRITE_IP)
          IF(WRITE_IP)THEN
            CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
            I=WORD_SIZE*(NP+1)/UNIT_SIZE
            CALL WRITE_DIRECT_INFO_V3(NP,I,'20-Aug-2000','IP_PP_DATA',LU_IP)
            OPEN(UNIT=LU_IP,FILE='IP_PP_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I,IOSTAT=IOS)
	    FREQ_CNT=0
            WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
            WRITE(LU_IP,REC=ACCESS_F)(MU(I),I=1,NP)
	  END IF
	END IF
        IF(WRITE_IP)THEN
	  T1=0.0D0
          IF(FREQ_CNT .NE. 0)READ(LU_IP,REC=ACCESS_F+FREQ_CNT)(IBOUND(LS),LS=1,NP),T1
	  IF(T1 .NE. FREQ)FREQ_CNT=FREQ_CNT+1
          WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
          WRITE(LU_IP,REC=ACCESS_F+FREQ_CNT)(CV_BOUND(LS),LS=1,NP),FREQ
        END IF
	CALL TUNE(2,'FG_BIG_SEC')
!
	RETURN
	END
