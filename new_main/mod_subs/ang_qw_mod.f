!
! Data module defining Quadrature Weights, and related variables,  
! for CMFGEN. 
!
	MODULE ANG_QW_MOD
!
! Created 02-May-2004
!
! Angular qudrature weights for J (AQW), H, and K. These are
! frequency indepoendent.
!
	REAL(10), ALLOCATABLE :: P(:)              !NP
	REAL(10), ALLOCATABLE :: AQW(:,:)          !ND,NP
	REAL(10), ALLOCATABLE :: HQW(:,:)          !ND,NP
	REAL(10), ALLOCATABLE :: KQW(:,:)          !ND,NP
	REAL(10), ALLOCATABLE :: NQW(:,:)          !ND,NP
!
! Defined at the mid points of the mesh.
!
	REAL(10), ALLOCATABLE :: HMIDQW(:,:)       !ND,NP
	REAL(10), ALLOCATABLE :: NMIDQW(:,:)       !ND,NP
!
! If required, these arrays shouLd have size NDEXT*NPEXT
!
	REAL(10), ALLOCATABLE :: PEXT(:)            !NPEXT
	REAL(10), ALLOCATABLE :: AQWEXT(:,:)
	REAL(10), ALLOCATABLE :: HQWEXT(:,:)
	REAL(10), ALLOCATABLE :: KQWEXT(:,:)
	REAL(10), ALLOCATABLE :: NQWEXT(:,:)
	REAL(10), ALLOCATABLE :: HMIDQWEXT(:,:)
	REAL(10), ALLOCATABLE :: NMIDQWEXT(:,:)
!
! Parameters, vectors, and arrays for computing the observed flux.
! Was 2000 (changed 20-Apr-2009).
!
        INTEGER, PARAMETER :: NST_CMF=6000
        REAL(10)  NU_STORE(NST_CMF)
!
        INTEGER NP_OBS_MAX
        INTEGER NP_OBS
        REAL(10) V_AT_RMAX                !Used if we extend the atmosphere.
        REAL(10) RMAX_OBS
        REAL(10) HFLUX_AT_OB		!In comoving frame
        REAL(10) HFLUX_AT_IB
!
! We allocate memory for the following vectors as we use them for the regular
! flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
!
        REAL(10), ALLOCATABLE :: IPLUS_STORE(:,:)
        REAL(10), ALLOCATABLE :: P_OBS(:)
        REAL(10), ALLOCATABLE :: IPLUS(:)
        REAL(10), ALLOCATABLE :: MU_AT_RMAX(:)
        REAL(10), ALLOCATABLE :: HQW_AT_RMAX(:)
!
	END MODULE ANG_QW_MOD
