!
! Data module for FG_J_CMF. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE FG_J_CMF_MOD_V11
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! The *_STORE routines are used to store the radiation field, as computed
! on the previous call to FG_J_CMF.
! The *_PREV routines are used to store the radiation field, as computed
! for the previous frequency.
!
! The *_PREV routines are updated from the *_STORE routines when NEW_FREQ
! is .TRUE.
!
!***************************************************************************
!***************************************************************************
!
! Data variables required by both routines.
!
! Dimensionded NRAY_MAX,NP
!
	REAL(KIND=LDP), ALLOCATABLE :: R_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: Z(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV(:,:)
!
! Dimensiond ND_EXT,4
!
	REAL(KIND=LDP), ALLOCATABLE :: V_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
!
	INTEGER, ALLOCATABLE :: J_PNT(:,:)
	INTEGER, ALLOCATABLE :: H_PNT(:,:)
!
! Dimensioned ND+ND_ADD=ND_EXT
!
	REAL(KIND=LDP), ALLOCATABLE :: R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: Z_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_CHI_EXT(:)
!
! Dimensioned NRAY_MAX,NP
!
	INTEGER, ALLOCATABLE :: RAY_PNT(:,:)
!
	INTEGER, ALLOCATABLE :: NI_RAY(:)
	INTEGER, ALLOCATABLE :: MAX_LS(:)
!
! Dimensioned NRAY_MAX
!
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR(:)
	REAL(KIND=LDP), ALLOCATABLE :: Q(:)
	REAL(KIND=LDP), ALLOCATABLE :: QH(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_RAY(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_RAY(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to the DIFFERENCE equation approach for  solving the
! transfer equation.
!	
	REAL(KIND=LDP), ALLOCATABLE :: AV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AV_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: PAR_AV(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: GAMH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TC(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GB(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: H(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: XM(:)
	REAL(KIND=LDP), ALLOCATABLE :: U(:)
	REAL(KIND=LDP), ALLOCATABLE :: VB(:)
	REAL(KIND=LDP), ALLOCATABLE :: VC(:)
	REAL(KIND=LDP), ALLOCATABLE :: DIV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI_STORE(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
! Dimensioned NRAY_MAX,NP
!
	REAL(KIND=LDP), ALLOCATABLE :: I_P_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_P_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: dGAMdR(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: A0(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A1(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A2(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A3(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A4(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: I_P(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M(:,:)
!
! Dimensioned NRAY_MAX
!
	REAL(KIND=LDP), ALLOCATABLE :: EE(:)
	REAL(KIND=LDP), ALLOCATABLE :: E0(:)
	REAL(KIND=LDP), ALLOCATABLE :: E1(:)
	REAL(KIND=LDP), ALLOCATABLE :: E2(:)
	REAL(KIND=LDP), ALLOCATABLE :: E3(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_PRIME(:)
	REAL(KIND=LDP), ALLOCATABLE :: S(:)
	REAL(KIND=LDP), ALLOCATABLE :: dS(:)
!
	REAL(KIND=LDP) PREVIOUS_FREQ
	REAL(KIND=LDP) VDOP_FRAC_SAV
!
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NP_MAX
	INTEGER NRAY_MAX
	INTEGER IDMIN,IDMAX
!
	CHARACTER*10 OLD_SOLUTION_OPTIONS
	CHARACTER*10 SOLUTION_METHOD
	LOGICAL INSERT
	LOGICAL FIRST_TIME
	LOGICAL NEW_R_GRID
	LOGICAL WRITE_IP
!
	DATA VDOP_FRAC_SAV/-1001.0_LDP/ 	!Absurd value.
	DATA FIRST_TIME/.TRUE./
	DATA PREVIOUS_FREQ/0.0_LDP/
!
	END MODULE FG_J_CMF_MOD_V11
!
! 
!
! Routine to compute the Eddington F, G and RSQN_ON_RSQJ Eddington factors
! for a single frequency. The transfer is done in the comoving-frame with
! first order frequency differencing. A DIFFERENCE or INTEGRAL equation
! approach can be used.
!
!The DIFFERENCE equation approach is at least 30% faster. However,
! the INTEGRAL equation approach is more stable. Because monotonic cubic
! interpolation is used for the Source function in the INTEGRAL equation
! approach, the intensities are guaranteed positive (or zero).
!
! The F, G, and RSQN_ON_RSQJ Eddington factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' G is defined at all depths, and
!       RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' RSQN_ON_RSQJ is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' one of G or RSQN_ON_RSQJ is
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
	SUBROUTINE FG_J_CMF_V11(ETA,CHI,ESEC,V,SIGMA,R,P,
	1                  JNU,FEDD,JQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  RETURNED_IN_HBC,RETURNED_OUT_HBC,IPLUS_P,
	1                  FREQ,dLOG_NU,DIF,DBB,IC,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  METHOD,SOLUTION_OPTIONS,
	1                  THK,INIT,NEW_FREQ,NC,NP,ND)
	USE SET_KIND_MODULE
	USE FG_J_CMF_MOD_V11
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Altered 07-Jul-2011: Changed extrapolation region I to max-value of 10 (rathee than ND/6).
!                         This was done for shell model. Mayu need further changes.
! Altered 21-Apr-2011: Bug fix for HNU_AT_OB when using DIFF option (not generally sed anyway).
! Altered 14-Jan-2010: Added JPLUS_IB, JMIN_IB etc to improve computation of
!                        boundary conditions.
! Altered 14-May-2009: Altered value of IDMAX (for ETA and CHI interpolaton).
!                        IDMAX & IDMIN now stored in FG_J_CMF_MOD_V11.
! Altered 19-Jan-2009: Changed IP_DATA to IP_FG_DATA which must exist for I(p)
!                        to be output.
! Created 17-Dec-2008: Based on FG_J_CMF_V10
!
	INTEGER NC,NP,ND
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP)
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!     QW denotes quadrature weight.
!
	REAL(KIND=LDP) JQW(ND,NP)
	REAL(KIND=LDP) HQW(ND,NP)
	REAL(KIND=LDP) KQW(ND,NP)
	REAL(KIND=LDP) NQW(ND,NP)
!
	REAL(KIND=LDP) HMIDQW(ND,NP)
	REAL(KIND=LDP) NMIDQW(ND,NP)
!
	REAL(KIND=LDP) JNU(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) IPLUS_P(NP)
!
! VDOP_VEC(I) is the minimum DOPPLER width for all species at depth I. It will include
! both a turbulent, and and thermal contribution for the ionization species with the
! highest mass. VDOP_FRAC is used to set the minimum velocity step size along a ray.
!
	REAL(KIND=LDP) VDOP_VEC(ND)
	REAL(KIND=LDP) VDOP_FRAC
	REAL(KIND=LDP) RETURNED_IN_HBC
	REAL(KIND=LDP) RETURNED_OUT_HBC
!
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) IC
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) dLOG_NU
!
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
	INTEGER ACCESS_F
        INTEGER LU_IP
	INTEGER IOS,REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC
!
! The following arrays do not need to be stored, and hence can be created
! each time.
!
	REAL(KIND=LDP) CV_BOUND(NP)		!Outer boundary V
	REAL(KIND=LDP) I_M_IN_BND(NP)		!Inner boundary
	REAL(KIND=LDP) IBOUND(NP)		!Incident intensity on outer boundary.
!
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(KIND=LDP) FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL(KIND=LDP), PARAMETER :: ONE=1
	INTEGER, PARAMETER :: NINS=4
!
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!
	INTEGER NI_SMALL
	INTEGER I,J,K,LS
	INTEGER NI
!
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) I_CORE
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DELTA_Z
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) ESEC_POW
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) RMAX,DEL_R_FAC
	REAL(KIND=LDP) MU,dZ,PSQ
	REAL(KIND=LDP) DEL_R
!
! Change the following statement to TRUE if running on a VECTOR machine.
!
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
	LOGICAL, PARAMETER :: DEFINE_AT_MID_POINTS=.FALSE.          !TRUE.
!
!
!
! Allocate data for moments which will be used to construct the Eddington
! factors.
!
	LUER=ERROR_LU()
	IF(.NOT. ALLOCATED(JNU_STORE))THEN
	  ND_STORE=ND
	  ALLOCATE (JNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (R_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RMID_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EXT_RMID_STORE(ND+1),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating JNU_STORE block in FG_J_CMF_V11: Status=',IOS
	    STOP	
          END IF
	END IF
!
	IF(INIT)THEN
	  R_STORE(1:ND)=R(1:ND)
	  GAM_REL_STORE(1:ND)=1.0_LDP/SQRT(1.0_LDP-(V(1:ND)/2.99792458E+10_LDP)**2)
	  DO I=1,ND-1
	    RMID_STORE(I)=0.5_LDP*(R(I)+R(I+1))
	    EXT_RMID_STORE(I+1)=0.5_LDP*(R(I)+R(I+1))
	  END DO
	  EXT_RMID_STORE(1)=R(1); EXT_RMID_STORE(ND+1)=R(ND)
	END IF
!
! Check to see whether we have a new R grid, or the solution options have
! changd. This can only happen when INIT is TRUE.
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
	    IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	      NEW_R_GRID=.TRUE.
	      J=ERROR_LU()
	      WRITE(J,*)'Warning: Updating RGRID in FG_J_CMF_V11'
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
	  DEALLOCATE ( Z_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  DEALLOCATE ( LOG_ETA_EXT )
	  DEALLOCATE ( LOG_CHI_EXT )
!
	  DEALLOCATE ( Z )
	  DEALLOCATE ( RAY_PNT )
	  DEALLOCATE ( R_RAY )
	  DEALLOCATE ( NI_RAY )
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
	  DEALLOCATE ( MAX_LS )
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
	    DEALLOCATE ( dGAMdR )
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
	    WRITE(J,*)'Error in FG_J_CMF_V11: Invalid solution option'
	    WRITE(J,*)SOLUTION_OPTIONS
	    STOP
	  END IF
!
          ND_ADD=0
          IF(THK)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( Z_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating R_EXT block in FG_J_CMF_V11: Status=',IOS
	    STOP	
          END IF
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    IF(V(ND) .LT. 10.0_LDP .AND. R(1)/R(ND) .GE. 9.99_LDP)THEN
	      RMAX=10.0_LDP*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0_LDP .OR. V(1) .GT. 2.0E+04_LDP)THEN
	      RMAX=1.5_LDP*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0_LDP,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-3) )
	    R_EXT(1)=RMAX
	    R_EXT(4)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.1_LDP*(R_EXT(1)-R_EXT(4))
	    R_EXT(3)=R_EXT(1)-0.4_LDP*(R_EXT(1)-R_EXT(4))
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
	    BETA=(SIGMA(1)+1.0_LDP)*(R(1)/R(ND)-1.0_LDP)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0_LDP-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0_LDP)-1.0_LDP
	    END DO
	    J=ERROR_LU()
	    WRITE(J,'(A)')' '
	    WRITE(J,*)'Using thick boundary condition in FG_J_CMF_V11'
	    WRITE(J,'(2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	    WRITE(J,'(2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	  END IF
!
! Define zone used to extrapolate opacities and emissivities.
!
	  IF(ND_ADD .NE. 0)THEN
	    IDMIN=1
	    T1=R_EXT(1)/R(1)
	    IF(T1 .LT. R(1)/R(ND/3))THEN
	      T1=MIN(3.0_LDP,T1)
	      DO I=1,ND
	        IF(R(1)/R(I) .GT. T1)THEN
	          IDMAX=MAX(4,I)
	          EXIT
	        END IF
	      END DO
	      IDMAX=MIN(IDMAX,ND/6)
	    ELSE
!	      IDMAX=ND/6
	      IDMAX=MIN(10,ND/6)
	    END IF
	    WRITE(6,'(A,I4,A)')' Using depth 1 and',IDMAX,' to extrapolate opacities'
	  END IF
!
! Work out the maximim number of points per ray. This will allow us to allocate
! the require memory.
!
	  NRAY_MAX=0.0
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0_LDP
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	   END DO
!
	    K=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. 1)K=K+INT(T1)
	      K=K+1
	    END DO
            NRAY_MAX=MAX(NRAY_MAX,K)
	  END DO
	  J=ERROR_LU()
	  WRITE(J,'(A,I6)')' Maximum number of points along ray in FG_J_CMF_V11 is',NRAY_MAX
!
	  ALLOCATE ( Z(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( RAY_PNT(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( R_RAY(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( NI_RAY(NP),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating Z block in FG_J_CMF_V11: Status=',IOS
	    STOP	
          END IF
!
! 
!
! Define the grid along each ray. NB: RAY_PNT is used to indicate the
! interpolation interval on the original (EXT) grid.
!
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0_LDP
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	    END DO
!
	    K=1
	    Z(1,LS)=Z_EXT(1)
	    RAY_PNT(1,LS)=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    DO I=1,NI_SMALL-1
	      T1=(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. 1)THEN
	        DELTA_Z=(Z_EXT(I+1)-Z_EXT(I))/(INT(T1)+1)
	        DO J=1,INT(T1)
	          K=K+1
	          Z(K,LS)=Z(K-1,LS)+DELTA_Z
	          RAY_PNT(K,LS)=I
	        END DO
	      END IF
	      K=K+1
	      Z(K,LS)=Z_EXT(I+1)
	      RAY_PNT(K,LS)=I
	    END DO
	    NI_RAY(LS)=K
!
	    PSQ=P(LS)*P(LS)
	    R_RAY(1,LS)=R_EXT(1)
	    DO I=2,NI_RAY(LS)
	      R_RAY(I,LS)=SQRT(Z(I,LS)*Z(I,LS)+PSQ)
	    END DO
	    R_RAY(NI_RAY(LS),LS)=R_EXT(NI_SMALL)
!
	  END DO
!
! 
!
! J_PNT and H_PNT are used to indicate the postion of J(I) amd H(I) along the ray
! so that J and H can be computed on our regular radius grid.
!
	  ALLOCATE ( J_PNT(ND,NP),STAT=IOS); J_PNT(:,:)=0
	  IF(IOS .EQ. 0)ALLOCATE ( H_PNT(ND,NP),STAT=IOS); H_PNT(:,:)=0
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating J_PNT block in FG_J_CMF_V11: Status=',IOS
	    STOP	
          END IF
!
	  DO LS=1,NP
	    NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	    K=1
	    DO I=1,NI_SMALL
	      DO WHILE(J_PNT(I,LS) .EQ. 0)
	        IF(R(I) .LE. R_RAY(K,LS) .AND. R(I) .GE. R_RAY(K+1,LS))THEN
	          IF( (R_RAY(K,LS)-R(I)) .LT. (R(I)-R_RAY(K+1,LS)) )THEN
	            J_PNT(I,LS)=K
	          ELSE
	            J_PNT(I,LS)=K+1
	          END IF
	        ELSE
	          K=K+1
	        END IF
	      END DO
	    END DO
	  END DO
!
	  DO I=1,ND
	    DO LS=1,NC+(ND-I+1)
	      K=J_PNT(I,LS)
	      IF(J_PNT(I,LS) .LT. 1 .OR. J_PNT(I,LS) .GT. NI_RAY(LS))THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_J_CMF_V11 -- invalid values'
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',J_PNT(I,LS)
	        STOP
	      ELSE IF( ABS(R_RAY(K,LS)-R(I))/R(I) .GT. 1.0E-12_LDP)THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_JCMF_V11 -- invalid values'
	        WRITE(LUER,*)'Fractional difference is',ABS(R_RAY(K,LS)-R(I))/R(I)
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',J_PNT(I,LS)
	        STOP
	      END IF
	    END DO
	  END DO
!
	  IF( SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! CV is defined on the nodes. Need to interpolate to the midpoint of R.
!
	    DO LS=1,NP
	      NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	      K=1
	      DO I=1,NI_SMALL-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5_LDP*(R(I)+R(I+1))
	          IF(T1 .LE. R_RAY(K,LS) .AND. T1 .GE. R_RAY(K+1,LS))THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	          END IF
	        END DO
	      END DO
	    END DO
	  ELSE
!
! CV is defined at the mid-points, but these don't, in general, correspond to
! the midpoints of the R grid.
!
	    DO LS=1,NP
	      K=1
	      DO I=1,ND-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5_LDP*(R(I)+R(I+1))
	          IF(T1 .LE. 0.5_LDP*(R_RAY(K,LS)+R_RAY(K+1,LS)) .AND. T1 .GE. 0.5_LDP*(R_RAY(K+1,LS)+R_RAY(K+2,LS)) )THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	          END IF
	        END DO
	      END DO
	    END DO
	  END IF
!
! 
!
	  ALLOCATE ( V_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_EXT,4),STAT=IOS )
!
!***************************************************************************
!***************************************************************************
!
	  IF(IOS .EQ. 0)ALLOCATE ( V_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SOURCE_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( Q(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( MAX_LS(NRAY_MAX),STAT=IOS )
!
	  IF(IOS .EQ. 0)ALLOCATE ( AV(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CV(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( DTAU(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( GAM(NRAY_MAX,NP),STAT=IOS )
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by both DIFFERENCE solution approach.
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    ALLOCATE ( AV_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( AV_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( PAR_AV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( CV_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( CV_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( GAMH(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( TA(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( TB(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( TC(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( GB(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( H(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( OLDCHI(NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( OLDCHI_STORE(NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( XM(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( U(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( VB(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( VC(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( QH(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( DIV(NRAY_MAX),STAT=IOS )
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating DIFFERENCE block in FG_J_CMF_V11: Status=',IOS
	      STOP	
            END IF
	  END IF
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by INTEGRAL solution.
!
	  IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
	    ALLOCATE (I_P(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (I_M(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( I_P_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_P_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_M_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_M_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( dGAMdR(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( A0(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A1(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A2(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A3(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A4(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE (EE(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E0(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E1(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E2(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E3(NRAY_MAX),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE (SOURCE_PRIME(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (S(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (dS(NRAY_MAX),STAT=IOS )
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating INTEGRAL block in FG_J_CMF_V11: Status=',IOS
	      STOP	
            END IF
!
	  END IF
!
	  FIRST_TIME=.FALSE.
	ELSE
	  IF(OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_V11'
	    WRITE(J,*)'Can''t switch SOLUTION_OPTIONS while runing code'
	    WRITE(J,*)'New setting:',SOLUTION_OPTIONS
	    WRITE(J,*)'Old setting:',OLD_SOLUTION_OPTIONS
	    STOP
	  END IF
	END IF
!
!
!
	NEG_AV_VALUE=.FALSE.
!
! Perform initializations.
!
	IF(INIT)THEN
!
! Insert extra points into radius grid. Not all points will be used along a ray.
! We only insert additional points in the interval between Z(NI)=0 and Z(NI-1),
! and between Z(NI-1) and  Z(NI-2).
!
	  FG_ERR_ON_FREQ(:)=0.0_LDP
	  FG_ERR_TYPE(:)=0
	  FG _ERR_CNT=0
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    AV_PREV(:,:)=0.0_LDP
	    CV_PREV(:,:)=0.0_LDP
          ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=0.0_LDP
	    I_M_PREV(:,:)=0.0_LDP
	  END IF
!
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(V_COEF,V_EXT,LOG_R_EXT,ND_EXT)
	  CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_EXT,LOG_R_EXT,ND_EXT)
!
	  DO LS=1,NP
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      T1=LOG(R_RAY(I,LS)/R_EXT(K))
	      V_RAY(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
	      SIGMA_RAY(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
	    END DO
!
!
! Compute GAMMA. This section is straight from the subroutine GAMMA, except
! That _EXT has been added to V, SIGMA, and R.
!
! We assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
!  	    (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      MU=Z(I,LS)/R_RAY(I,LS)
	      T1=3.33564E-06_LDP*V_RAY(I)/R_RAY(I,LS)
	      GAM(I,LS)=T1*( 1.0_LDP+SIGMA_RAY(I)*(MU**2) )
	    END DO
!
	    IF(SOLUTION_METHOD .EQ. 'INTEGRAL' .AND. NI_RAY(LS) .NE. 1)THEN
	      DO I=1,NI_RAY(LS)
	        MU=Z(I,LS)/R_RAY(I,LS)
	        T1=3.33564E-06_LDP*V_RAY(I)/R_RAY(I,LS)
	        dGAMdR(I,LS)=GAM(I,LS)*SIGMA_RAY(I)/R_RAY(I,LS)
	        J=MAX(I-1,1); K=MIN(NI,I+1)
	        dGAMdR(I,LS)=dGAMdR(I,LS)+T1*(MU**2)*
	1           (SIGMA_RAY(J)-SIGMA_RAY(K))/(R_RAY(J,LS)-R_RAY(K,LS))
	      END DO
	    ELSE IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	      DO I=1,NI_RAY(LS)-1
	        MU=(Z(I,LS)+Z(I+1,LS))/(R_RAY(I,LS)+R_RAY(I+1,LS))
	        GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564E-06_LDP*
	1                   (  1.0_LDP+0.5_LDP*(MU**2)*
	1                   (SIGMA_RAY(I)+SIGMA_RAY(I+1))  )/
	1                   (R_RAY(I,LS)+R_RAY(I+1,LS))
	     END DO
	     GAMH(NI,LS)=0.0
	    END IF
!
	  END DO		!LS Loop
!
! Determine maximum ray index having I points. Used in the
! vector section of the DIFFERENCE approach.
!
	  DO I=1,NRAY_MAX
	    DO LS=1,NP
	      IF(NI_RAY(LS) .GE. I)MAX_LS(I)=LS
	    END DO
	  END DO
!
	ELSE IF(NEW_FREQ)THEN
          IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    OLDCHI(1:NP)=OLDCHI_STORE(1:NP)
	    AV_PREV(:,:)=AV_STORE(:,:)
	    CV_PREV(:,:)=CV_STORE(:,:)
	  ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=I_P_STORE(:,:)
	    I_M_PREV(:,:)=I_M_STORE(:,:)
	  END IF
	  IF(FREQ .GE. PREVIOUS_FREQ)THEN
	    WRITE(ERROR_LU(),*)'Error in FG_J_CMF_V11'
	    WRITE(ERROR_LU(),*)'Frequencies must be monotonically decreasng'
	    STOP
	  END IF
	ELSE IF(FREQ .NE. PREVIOUS_FREQ)THEN
	   WRITE(ERROR_LU(),*)'Error in FG_J_CMF_V11'
	   WRITE(ERROR_LU(),*)
	1      'Frequencies must not change if NEW_FREQ=.FALSE.'
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
	    IF(ESEC_POW .LT. 2._LDP)ESEC_POW=2
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2._LDP)ALPHA=2.0
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2._LDP)ESEC_POW=2
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
	  IF(ALPHA .LT. 3.5_LDP)ALPHA=3.5
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(IDMIN)*(R(IDMIN)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0E-280_LDP)ETA_EXT(I)=1.0E-280_LDP
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
!
	LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	CALL TUNE(2,'FG_CHI_BEG')
!
	JPLUS_IB=0.0_LDP; HPLUS_IB=0.0_LDP; KPLUS_IB=0.0_LDP; NPLUS_IB=0.0_LDP
	JMIN_IB=0.0_LDP;  HMIN_IB=0.0_LDP; KMIN_IB=0.0_LDP;   NMIN_IB=0.0_LDP
	JPLUS_OB=0.0_LDP;  HPLUS_OB=0.0_LDP; KPLUS_OB=0.0_LDP; NPLUS_OB=0.0_LDP
	JMIN_OB=0.0_LDP;   HMIN_OB=0.0_LDP;  KMIN_OB=0.0_LDP;  NMIN_OB=0.0_LDP
!
! 
!***************************************************************************
!***************************************************************************
!
	IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
!
! Zero AV and CV matrices.
!
	  AV(:,:)=0.0_LDP
	  CV(:,:)=0.0_LDP
!
!	  CALL TUNE(1,'LS_LOOP')
!
	  DO LS=1,NP
	    NI=NI_RAY(LS)
!
! NB: dCHIdR = LOG(CHI)/LOG(R) * CHI/R
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0_LDP*CHI_COEF(K,1)*T1+2.0_LDP*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	       END IF
	     END DO
!
	    SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!
! Even in the THK case we assume that the intensity incident on the
! outer boundary is zero. In THK case we effectively get IBOUND not
! equal zero at the model atmosphere boundary as a consequence of the
! extension.
!
	    IBOUND(LS)=0.0_LDP
!
! 
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      DO I=1,NI
	        Q(I)=0.0_LDP
	        QH(I)=0.0_LDP
	      END DO
	      OLDCHI(LS)=CHI_RAY(NI)
	    ELSE IF(NEW_FREQ)THEN
	      DO I=1,NI-1
	        QH(I)=GAMH(I,LS)*2.0_LDP/((CHI_RAY(I)+CHI_RAY(I+1))*dLOG_NU)
	        Q(I)=GAM(I,LS)/(CHI_RAY(I)*dLOG_NU)
	      END DO
	      QH(NI)=0.0_LDP
	      Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)
	    END IF
!
	    IF(DIF .AND. LS .LE. NC)THEN
	      T1=0.0_LDP
	      IF(.NOT. INIT)T1=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)	     !Q(NI)
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI_RAY(NI)
	1            *(1.0_LDP+T1*(1.0_LDP-CHI_RAY(NI)/OLDCHI(LS)))
	    END IF
	    OLDCHI_STORE(LS)=CHI_RAY(NI)
!
	    IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is not necessary when
! we re using monotonic interpolation.
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO
	      ELSE IF(INSERT .OR. METHOD(4:6) .EQ. 'MON')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	       END DO
	      END IF
!
	      CALL TUVGHD_RH(TA(1,LS),TB(1,LS),TC(1,LS),
	1            U,VB,VC,GB(1,LS),H(1,LS),XM,
	1            Q,QH,DTAU(1,LS),SOURCE_RAY,DIF,DBC,IC,LS,NC,NI)
	      XM(1)=-IBOUND(LS)
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
	      TC(NI,LS)=0				!As used.
	      DIV(1)=1.0_LDP/(TC(1,LS)+TB(1,LS))
	      TC(1,LS)=TC(1,LS)*DIV(1)
	      TB(1,LS)=TB(1,LS)*DIV(1)
	      DO I=2,NI
	        DIV(I)=1.0_LDP/(TA(I,LS)*TB(I-1,LS)+TB(I,LS)+TC(I,LS))
	        TB(I,LS)=(TA(I,LS)*TB(I-1,LS)+TB(I,LS))*DIV(I)
	        TC(I,LS)=TC(I,LS)*DIV(I)
	      END DO
	      TB(1:NI,LS)=DIV(1:NI)
	    ELSE
!
	      XM(1)=0.0_LDP
	      DO I=2,NI-1
	        XM(I)=-0.5_LDP*SOURCE_RAY(I)*(DTAU(I-1,LS)+DTAU(I,LS))
	      END DO
!
	      IF(LS .GT. NC)THEN
	        XM(NI)=0.5_LDP*DTAU(NI-1,LS)*SOURCE_RAY(NI)
	      ELSE IF(DIF)THEN
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
	    DO I=2,NRAY_MAX
	      DO LS=1,MAX_LS(I)
	        IF(I .LE. NI_RAY(LS))
	1           AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	    END DO
	    DO I=NRAY_MAX,1,-1
	      DO LS=1,MAX_LS(I)
	        IF(I .LT. NI_RAY(LS))
	1          AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
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
	      DO I=2,NI_RAY(LS)
	        AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	      DO I=NI_RAY(LS)-1,1,-1
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
	   NI=NI_RAY(LS)
!
	    K=2
	    DO I=1,NI
	      T1=1.0E-08_LDP*ABS(AV(K,LS))
	      IF(AV(I,LS) .LE. T1)THEN
	        AV(I,LS)=MAX(0.01_LDP*ABS(AV(I,LS)),T1)
	        NEG_AV_VALUE=.TRUE.
	      ELSE
	        K=K+1
	      END IF
	      IF(I .EQ. 1)K=1
	    END DO
!
! Update C vector (i.e. flux variable).
!
	    DO I=1,NI_RAY(LS)-1
	      CV(I,LS)=GB(I,LS)*(AV(I,LS)-AV(I+1,LS))+H(I,LS)*CV_PREV(I,LS)
	    END DO
!
! Note that V=AV(1)-IBOUND.
!
	    IF(ND_ADD .EQ. 0)THEN
	      CV_BOUND(LS)=AV(1,LS)-IBOUND(LS)
	    ELSE IF(LS .EQ. NP)THEN
	      CV_BOUND(LS)=0.0_LDP
	    ELSE
	      K=J_PNT(1,LS)
	      CV_BOUND(LS)=0.5_LDP*(CV(K+1,LS)+CV(K-1,LS))
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
!$OMP PARALLEL DO PRIVATE(SOURCE_PRIME,SOURCE_RAY,CHI_RAY,ETA_RAY,dCHIdR_RAY,Q,EE,E0,E1,E2,E3,S,dS,T1,T2,I_CORE,dZ,NI,I,K)
!
	  DO LS=1,NP
!
	    IF(NI_RAY(LS) .EQ. 1)THEN
	      I_P(1,LS)=0.0_LDP
	      I_M(1,LS)=0.0_LDP
	      GOTO 1000
	    END IF
!
! NB: dCHIdR = LOG(CHI)/LOG(R) * CHI/R
!
!	    CALL TUNE(1,'FG_CHI_Q')
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0_LDP*CHI_COEF(K,1)*T1+2.0_LDP*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	      END IF
	    END DO
!
	    IF(DIF .AND. LS .LE. NC)THEN
	      I_CORE=( ETA_RAY(NI)+
	1         DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND) )/CHI_RAY(NI)
	    ELSE
	      I_CORE=IC
	    END IF
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      Q(1:NI)=0.0_LDP
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	    ELSE
	      Q(1:NI)=GAM(1:NI,LS)/dLOG_NU
	      CHI_RAY(1:NI)=CHI_RAY(1:NI)+Q(1:NI)
	      dCHIdR_RAY(1:NI)=dCHIdR_RAY(1:NI)+dGAMdR(1:NI,LS)/dLOG_NU
	      Q(1:NI)=Q(1:NI)/CHI_RAY(1:NI)
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
	    END IF
!	    CALL TUNE(2,'FG_CHI_Q')
!
!
!
!	    CALL TUNE(1,'FG_NEW_F')
	    IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is still necessary when
! we re using monotonic interpolation, since the monotonic interpolation only
! applies to CHI. CHI_RAY contains an additional term due to the frequency
! derivative. As non relativistic, DTAU is the same for both directions.
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	        END DO
	      END IF
!
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	      DO I=1,NI-1
	        T1=DTAU(I,LS)
	        EE(I)=0.0_LDP
	        IF(T1 .LT. 700.0_LDP)EE(I)=EXP(-T1)
	        IF(T1 .GT. 0.5_LDP)THEN
	          E0(I)=1.0_LDP-EE(I)
	          E1(I)=1.0_LDP-E0(I)/T1
	          E2(I)=1.0_LDP-2.0_LDP*E1(I)/T1
	          E3(I)=1.0_LDP-3.0_LDP*E2(I)/T1
	        ELSE IF(T1 .GT. 0.1_LDP)THEN
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP*
	1               (1.0_LDP-T1/10.0_LDP*(1.0_LDP-T1/11.0_LDP*
	1               (1.0_LDP-T1/12.0_LDP*(1.0_LDP-T1/13.0_LDP)))))))) )
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	        ELSE
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP) ))))
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	        END IF
	      END DO
!
              DO I=1,NI-1
	        A0(I,LS)=EE(I)
	        A1(I,LS)=E0(I)-3.0_LDP*E2(I)+2.0_LDP*E3(I)
	        A2(I,LS)=3.0_LDP*E2(I)-2.0_LDP*E3(I)
	        A3(I,LS)=DTAU(I,LS)*(E1(I)-2.0_LDP*E2(I)+E3(I))
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
	1                      MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1               MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1               MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
            I_M(1,LS)=0.0_LDP
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
	1                      MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1            MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
	    IF(LS .LE. NC)THEN
	      I_P(NI,LS)=I_CORE
	    ELSE
	      I_P(NI,LS)=I_M(NI,LS)
	    END IF
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
	    K=J_PNT(1,LS)
	    CV_BOUND(LS)=0.5_LDP*(I_P(K,LS)-I_M(K,LS))
	    I_M_IN_BND(LS)=I_M(NI,LS)
!
!	    CALL TUNE(2,'FG_INT_INT')
	  END DO			!LS
!$OMP END PARALLEL DO
!
! Compute the mean intensity like variable U at each grid point. Used to
! compute J.
!
!$OMP PARALLEL DO
!
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      AV(I,LS)=0.5_LDP*( I_P(I,LS)+I_M(I,LS) )
	    END DO
	  END DO
!
!$OMP END PARALLEL DO
!
! We define CV on the gregular grid. We perform the interpolation onto
! the original R grid as we compute H.
!
!$OMP PARALLEL DO
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      CV(I,LS)= 0.5_LDP*( I_P(I,LS)-I_M(I,LS) )
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	  I_P_STORE(:,:)=I_P(:,:)
	  I_M_STORE(:,:)=I_M(:,:)
!
!	  IF(NEW_FREQ)THEN
!	     WRITE(117,*)FREQ
!	     DO LS=1,NP
!	       WRITE(117,'(5ES15.6)')P(LS),I_M_PREV(1,LS),I_P_PREV(1,LS), I_M_PREV(2,LS),I_P_PREV(2,LS)
!	     END DO
!	  END IF
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
	HBC=0.0_LDP			!H/J at model outer boundary.
	NBC=0.0_LDP			!N/J at model outer boundary.
	IN_HBC=0.0_LDP
!
! Initialize intensity matrices.
!
	JNU_STORE(:)=0.0_LDP			!1:ND
	HNU_STORE(:)=0.0_LDP
	KNU_STORE(:)=0.0_LDP
	NNU_STORE(:)=0.0_LDP
!
! J and K are always evaluated on the nodes.
!
	CALL TUNE(1,'J_LOOP')
	DO LS=1,NP
	  DO I=1,MIN(ND,ND-(LS-NC-1))
	    K=J_PNT(I,LS)
	    JNU_STORE(I)=JNU_STORE(I)+JQW(I,LS)*AV(K,LS)
	    KNU_STORE(I)=KNU_STORE(I)+KQW(I,LS)*AV(K,LS)
	  END DO
	END DO
!
! This procedure works for the INTEGRAL or DIFFERENCE approach.
!
	DO LS=1,NP
	  K=J_PNT(1,LS)
	  T1=AV(K,LS)+CV_BOUND(LS)
	  T2=0.0_LDP; T2=MAX(T2,AV(K,LS)-CV_BOUND(LS))
	  JPLUS_OB=JPLUS_OB+JQW(1,LS)*T1
	  HPLUS_OB=HPLUS_OB+HQW(1,LS)*T2
	  KPLUS_OB=KPLUS_OB+KQW(1,LS)*T1
	  NPLUS_OB=NPLUS_OB+NQW(1,LS)*T1
	  JMIN_OB=JMIN_OB+JQW(1,LS)*T2
	  HMIN_OB=HMIN_OB+HQW(1,LS)*T2
	  KMIN_OB=KMIN_OB+KQW(1,LS)*T2
	  NMIN_OB=NMIN_OB+NQW(1,LS)*T2
	END DO
	DO LS=1,NC+1
	  K=J_PNT(ND,LS)
	  T1=AV(K,LS)-0.5_LDP*I_M_IN_BND(LS)
	  T2=I_M_IN_BND(LS)
	  JPLUS_IB=JPLUS_IB+JQW(ND,LS)*T1
	  HPLUS_IB=HPLUS_IB+HQW(ND,LS)*T2
	  KPLUS_IB=KPLUS_IB+KQW(ND,LS)*T1
	  NPLUS_IB=NPLUS_IB+NQW(ND,LS)*T1
	  JMIN_IB=JMIN_IB+JQW(ND,LS)*T2
	  HMIN_IB=HMIN_IB+HQW(ND,LS)*T2
	  KMIN_IB=KMIN_IB+KQW(ND,LS)*T2
	  NMIN_IB=NMIN_IB+NQW(ND,LS)*T2
	END DO
!
! Evaluate the H & N moments. These may be evaluated at the nodes,
! or at the node midpoints, depending on the solution technique.
!
! With the integral technique, H & N are evaluated at the nodes.
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
          IF(DEFINE_AT_MID_POINTS)THEN
	    DO LS=1,NP
	      DO I=1,MIN(ND,ND-(LS-NC-1))-1
                K=H_PNT(I,LS)
                T2=0.5_LDP*(R(I)+R(I+1))
                T1=(T2-R_RAY(K,LS))/(R_RAY(K+1,LS)-R_RAY(K,LS))
                HNU_STORE(I)=HNU_STORE(I)+HMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
                NNU_STORE(I)=NNU_STORE(I)+NMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
              END DO
            END DO
	    HNU_AT_OB=0.0_LDP; NNU_AT_OB=0.0_LDP
	    DO LS=1,NP
	      K=H_PNT(1,LS)
	      HNU_AT_OB=HNU_AT_OB+HQW(1,LS)*CV(K,LS)
	      NNU_AT_OB=NNU_AT_OB+NQW(1,LS)*CV(K,LS)
	    END DO
	    HN_DEF_ON_NODES=.FALSE.
	  ELSE
	    HN_DEF_ON_NODES=.TRUE.
	    DO LS=1,NP
	      DO I=1,MIN(ND,ND-(LS-NC-1))
	        K=J_PNT(I,LS)
	        HNU_STORE(I)=HNU_STORE(I)+HQW(I,LS)*CV(K,LS)
	        NNU_STORE(I)=NNU_STORE(I)+NQW(I,LS)*CV(K,LS)
	      END DO
	      IPLUS_P(LS)=2.0_LDP*CV_BOUND(LS)
	    END DO
	    HNU_AT_OB=HNU_STORE(1); NNU_AT_OB=NNU_STORE(1)
	    HNU_AT_IB=HNU_STORE(ND); NNU_AT_IB=NNU_STORE(ND)
	  END IF
	ELSE
!
! With the DIFFERENCE technique, H & N are evaluated at the midpoints
! of the nodes. We only need to interpolate when we have included
! additional points along each ray.
!
	  HN_DEF_ON_NODES=.FALSE.
	  DO LS=1,NP
	    DO I=1,MIN(ND,ND-(LS-NC-1))-1
	      K=H_PNT(I,LS)
	      T2=0.5_LDP*(R(I)+R(I+1))
	      T1=( 2.0_LDP*T2-(R_RAY(K,LS)+R_RAY(K+1,LS)) )/(R_RAY(K+2,LS)-R_RAY(K,LS))
	      HNU_STORE(I)=HNU_STORE(I)+HMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
	      NNU_STORE(I)=NNU_STORE(I)+NMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
	    END DO
	  END DO			!End do LS
!
	  HNU_AT_OB=0.0_LDP; NNU_AT_OB=0.0_LDP
	  DO LS=1,NP
	    HNU_AT_OB=HNU_AT_OB+CV_BOUND(LS)*HQW(I,LS)
	    NNU_AT_OB=NNU_AT_OB+CV_BOUND(LS)*NQW(I,LS)
	    IPLUS_P(LS)=2.0_LDP*CV_BOUND(LS)
	  END DO
	  IF(DIF)THEN
	    HNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/3.0_LDP
	    NNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/5.0_LDP
	  ELSE
	    HNU_AT_IB=0.0_LDP; NNU_AT_IB=0.0_LDP
	    DO LS=1,NC+1
	      HNU_AT_IB=HNU_AT_IB+0.5_LDP*(IC-I_M_IN_BND(LS))*HQW(ND,LS)
	      NNU_AT_IB=NNU_AT_IB+0.5_LDP*(IC-I_M_IN_BND(LS))*NQW(ND,LS)
	    END DO
	  END IF
	END IF
!
! Get boundary conditons for moment calculations.
!
	IN_HBC=0.0_LDP
	DO LS=1,NC+1
	   IN_HBC=IN_HBC + HQW(ND,LS)*I_M_IN_BND(LS)
	END DO
!
	CALL TUNE(2,'J_LOOP')
!
!
!
! Compute the boundary Eddington factors.
!
	HBC=HNU_AT_OB/JNU_STORE(1)
	NBC=NNU_AT_OB/JNU_STORE(1)
	IN_HBC=IN_HBC/(2.0_LDP*JNU_STORE(ND)-IC)
	RETURNED_OUT_HBC=HBC
	RETURNED_IN_HBC=IN_HBC
!
! Compute the Eddington factor, F. This required in CMFGEN to
! compute the K moment (from J computed by MOM_J_CMF).
!
	DO I=1,ND
	  JNU(I)=JNU_STORE(I)
	  FEDD(I)=KNU_STORE(I)/JNU_STORE(I)
	END DO
!
! Store frequencies at which errors occurred, and give an indication of the
! error.
!
	DO J=1,N_ERR_MAX
	  IF(.NOT. NEG_AV_VALUE)EXIT
	  IF(FG_ERR_ON_FREQ(J) .EQ. FREQ)THEN
	    FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    NEG_AV_VALUE=.FALSE.
	  ELSE IF(J .EQ. FG_ERR_CNT+1)THEN
	    FG_ERR_CNT=J
	    FG_ERR_TYPE(J)=1
	    FG_ERR_ON_FREQ(J)=FREQ
	    NEG_AV_VALUE=.FALSE.
	  END IF
	END DO
!
	DO I=1,ND
	  IF(JNU_STORE(I) .LT. 0)THEN
	    OPEN(UNIT=7,FILE='FG_J_CMF_V11_ERRORS',STATUS='UNKNOWN')
	      WRITE(7,*)'FREQ=',FREQ
	      WRITE(7,'(3X,A,5X,4(5X,A,5X))')'I','JNU_STORE','ETA','CHI','ESEC'
	      DO J=1,ND
	        WRITE(7,'(X,I5,4ES16.6)')J,JNU_STORE(J),ETA(J),CHI(J),ESEC(J)
	      END DO
	    CLOSE(UNIT=7)
	    J=ERROR_LU()
	    WRITE(J,*)'Error on FG_J_CMF_V11 --- negative mean intensities.'
	    WRITE(J,*)'Check out file FG_J_CMF_V11_ERRORS for aditional information.'
	    WRITE(J,*)'Halting code execution.'
	    STOP
	  END IF
	END DO
!
! This will only output IP_FG_DATA if a file IP_FG_DATA exists.
! This is only for debugging purposes.
!
        ACCESS_F=5
	LU_IP=56
        IF(INIT)THEN
	  INQUIRE(FILE='IP_FG_DATA',EXIST=WRITE_IP)
	  IF(WRITE_IP)THEN
	    CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
            I=WORD_SIZE*(NP+1)/UNIT_SIZE
            CALL WRITE_DIRECT_INFO_V3(NP,I,'20-Aug-2000','IP_FG_DATA',LU_IP)
            OPEN(UNIT=LU_IP,FILE='IP_FG_DATA',FORM='UNFORMATTED',
	1         ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I,IOSTAT=IOS)
	    FREQ_CNT=0
            WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
            WRITE(LU_IP,REC=ACCESS_F)(P(I),I=1,NP)
	  END IF
	END IF
        IF(WRITE_IP)THEN
	  T1=0.0_LDP
          IF(FREQ_CNT .NE. 0)READ(LU_IP,REC=ACCESS_F+FREQ_CNT)(IBOUND(LS),LS=1,NP),T1
	  IF(T1 .NE. FREQ)FREQ_CNT=FREQ_CNT+1
          WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
          WRITE(LU_IP,REC=ACCESS_F+FREQ_CNT)(CV_BOUND(LS),LS=1,NP),FREQ
        END IF
	CALL TUNE(2,'FG_BIG_SEC')
!
	RETURN
	END
