!
! Subroutine to compute the opacity and emissivity  variation due to FREE-FREE
! processes as a function of the level populations for a general ion.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine can handle ionizations to differnt super levels.
!
! Notation:
!
!         We use _F to denote poupulations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote poupulations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE VAR_FREE_FREE_V1(VCHI,VETA,VCHI_NOT_IMP,VETA_NOT_IMP,
	1             HN_S,HNST_S,LOG_HNST_S,dlnHNST_S_dlnT,N_S,
	1	      HNST_F_ON_S,EDGE_F,N_F,F_TO_S_MAPPING,
	1             DI_S,LOG_DIST_S,dlnDIST_S_dlnT,N_DI,
	1             PHOT_ID,ION_LEV,ED,T,EMHNUKT,IMP_VAR,
	1             NU,Z,ID,IONFF,
	1             EQHN,GS_ION_EQ,NT,ND,LST_DEPTH_ONLY)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
! Created 11-Jul-2019 - Basd on VAR_OP_V10
C
	INTEGER ID
	INTEGER N_S
	INTEGER N_F
	INTEGER N_DI		!Number of levels in ion
	INTEGER EQHN
	INTEGER GS_ION_EQ
	INTEGER NT,ND
	LOGICAL IONFF			!Include free-free opacity for level?
	LOGICAL LST_DEPTH_ONLY	!for computing dTdR
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	REAL*8, TARGET :: VCHI(NT,ND)		!VCHI(I,K)=dCHI(K)/dN(I,K)
	REAL*8, TARGET :: VETA(NT,ND)		!VETA(I,K)=dETA(K)/dN(I,K)
	REAL*8, TARGET :: VCHI_NOT_IMP(NT,ND)
	REAL*8, TARGET :: VETA_NOT_IMP(NT,ND)
C
	REAL*8 HN_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 LOG_HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S_MAPPING(N_F)
C
C Ion population information.
C
	REAL*8 DI_S(N_DI,ND)
	REAL*8 LOG_DIST_S(N_DI,ND)
	REAL*8 dlnDIST_S_dlnT(N_DI,ND)
C
	LOGICAL IMP_VAR(NT)
C
	INTEGER PHOT_ID		!Photoionization ID
	INTEGER ION_LEV		!target level in ION for ionizations.
C
	REAL*8 NU			!Frequency (10^15 Hz)
	REAL*8 Z			!Charge on ion
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperatuure (10^4 K)
	REAL*8 EMHNUKT(ND)		!exp(-hv/kT)
C
C Local variables.
C
	REAL*8, POINTER :: PCHI(:,:)
	REAL*8, POINTER :: PETA(:,:)
C
	REAL*8 GFF_VAL(ND)		!Used as work vector
	REAL*8 LOG_DI_RAT(ND)		!Used as work vector
	REAL*8 DT_TERM(ND)		!Used as work vector
	REAL*8 HDKT_ON_T(ND)		!Used as work vector
	REAL*8 POP_SUM(ND)
C
	REAL*8 YDIS(ND)			!Constant for computing level dissolution/
	REAL*8 XDIS(ND)			!Constant for computing level dissolution/
	REAL*8 DIS_CONST(N_F)		!Constant appearing in dissolution formula.
	REAL*8 ALPHA_VEC(N_F)		!Photionization cross-section
	REAL*8 VCHI_TMP(N_F,ND)
	REAL*8 SUM_ION; REAL*8 SUM_T1; REAL*8 SUM_T2
	REAL*8 SUM_ION_NOT_IMP; REAL*8 SUM_T1_NOT_IMP; REAL*8 SUM_T2_NOT_IMP
	REAL*8 NEFF,ZION_CUBED,T1,T2
C
	INTEGER ND_LOC
	INTEGER I                     !Used as level index (same) in atom.
	INTEGER L			!index of level in full atom.
	INTEGER K_ST,K		!Used as depth index.
	INTEGER GENLEV		!Level index in VCHI, VETA
	INTEGER EQION			!Ion variable in VCHI,VETA
	INTEGER NO_NON_ZERO_PHOT
C
	REAL*8 TCHI1,TCHI2,TETA1,TETA2,TETA3
	REAL*8 HNUONK,ALPHA
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	HNUONK=HDKT*NU
	EQION=GS_ION_EQ+(ION_LEV-1)
!
! If EQION is not important, all free-free and bound-free processes
! are added to the IMP variation.
!
	IF( IMP_VAR(EQION) )THEN
	  PCHI=>VCHI
	  PETA=>VETA
	ELSE
	  PCHI=>VCHI_NOT_IMP
	  PETA=>VETA_NOT_IMP
	END IF
C
C ND_LOC indicates the number of depth points we are going to compute the
C opacity at.
C
C K_ST indicates the depth point to start, and is either 1, or ND.
C
	IF(LST_DEPTH_ONLY)THEN
	  ND_LOC=1
	  K_ST=ND
	ELSE
	  ND_LOC=ND
	  K_ST=1
	END IF
C
C Free-free processes
C
	IF( IONFF .AND. PHOT_ID .EQ. 1)THEN
C
C Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
C
	  IF(LST_DEPTH_ONLY)THEN
	    CALL GFF_VEC(GFF_VAL(ND),NU,T(ND),Z,ND_LOC)
	    POP_SUM(ND)=SUM(DI_S(:,ND))
	  ELSE
	    CALL GFF_VEC(GFF_VAL,NU,T,Z,ND_LOC)
	    POP_SUM=SUM(DI_S,1)
	  END IF
C
	  TCHI1=CHIFF*Z*Z/( NU**3 )
	  TETA1=CHIFF*Z*Z*TWOHCSQ
!
!$OMP PARALLEL DO PRIVATE(ALPHA,TCHI2,TETA2,I,K,L)
	  DO K=K_ST,ND
	    ALPHA=GFF_VAL(K)/SQRT(T(K))
C
	    TCHI2=TCHI1*ALPHA
	    PCHI(NT-1,K)=PCHI(NT-1,K)+POP_SUM(K)*TCHI2*(1.0D0-EMHNUKT(K))
	    PCHI(NT,K)=PCHI(NT,K)+ED(K)*POP_SUM(K)*TCHI2/T(K)*( -0.5D0+(0.5D0-HNUONK/T(K))*EMHNUKT(K) )
C
	    TETA2=TETA1*ALPHA*EMHNUKT(K)
	    PETA(NT-1,K)=PETA(NT-1,K)+TETA2*POP_SUM(K)
	    PETA(NT,K)  =PETA(NT,K)  +TETA2*POP_SUM(K)*ED(K)*(HNUONK/T(K)-0.5D0)/T(K)
!
	    TCHI2=TCHI2*ED(K)*(1.0D0-EMHNUKT(K))
	    TETA2=TETA2*ED(K)
	    DO I=1,N_DI
	      L=EQION+I-1
	      PCHI(L,K)=PCHI(L,K)+TCHI2
	      PETA(L,K)=PETA(L,K)+TETA2
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	END IF
!
! Nullify pointer assignments.
!
	NULLIFY (PCHI)
	NULLIFY (PETA)
!
	RETURN
	END
