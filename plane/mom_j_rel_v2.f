	MODULE MOD_J_REL
	IMPLICIT NONE
!
	INTEGER MOM_ERR_CNT
	INTEGER, PARAMETER :: N_ERR_MAX=1000
	REAL(KIND=LDP) MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! These vectors must be saved as they will be used on subsequent iterations.
!
	REAL(KIND=LDP), ALLOCATABLE :: AV_SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: BETA(:)
	REAL(KIND=LDP), ALLOCATABLE :: BETA_FREQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_REL(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_REL_SQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_DELTA(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_DELTAH(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_dKdNUH(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_dNdNUH(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_dKdNU(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_dHdNU(:)
!
        REAL(KIND=LDP), ALLOCATABLE :: FEDD_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: GEDD_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: H_ON_J_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: RSQN_ON_RSQJ_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: KMID_ON_J_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: JNU_PREV(:)
        REAL(KIND=LDP), ALLOCATABLE :: GAM_RSQHNU_PREV(:)
!
        REAL(KIND=LDP), ALLOCATABLE :: FEDD_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: GEDD_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: H_ON_J_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: RSQN_ON_RSQJ_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: KMID_ON_J_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: JNU_SAVE(:)
        REAL(KIND=LDP), ALLOCATABLE :: GAM_RSQHNU_SAVE(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: TA(:),TB(:),TC(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_H(:),CHI_J(:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU_H(:),DTAU_J(:),DTAUONQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: Q(:),XM(:),SOURCE(:)
	REAL(KIND=LDP), ALLOCATABLE :: VB(:),VC(:),COH_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: HU(:),HL(:),HS(:),HD(:)
	REAL(KIND=LDP), ALLOCATABLE :: P_H(:),P_J(:),JOLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: VdJdR_TERM(:),VdHdR_TERM(:)
	REAL(KIND=LDP), ALLOCATABLE :: DELTA(:),DELTAH(:),W(:),WPREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: PSI(:),PSIPREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS(:),EPS_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_RSQJNU_PREV(:)
!
	REAL(KIND=LDP) HBC_PREV,HBC_SAVE
	REAL(KIND=LDP) NBC_PREV,NBC_SAVE
	REAL(KIND=LDP) IN_HBC_PREV,IN_HBC_SAVE
!
	REAL(KIND=LDP) FREQ_SAVE
	INTEGER ND_SAV
!
	END MODULE MOD_J_REL
!
!
! Subroutine to allocate the required vectors.
!
	SUBROUTINE ALLOC_MOD_J_REL(ND)
	USE SET_KIND_MODULE
	USE MOD_J_REL
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
!
	IF(.NOT. ALLOCATED(AV_SIGMA))THEN
!
	  FREQ_SAVE=0.0_LDP
	  ND_SAV=ND
!
	  ALLOCATE (AV_SIGMA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (BETA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (BETA_FREQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL_SQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_DELTA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_DELTAH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dKdNUH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dNdNUH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dKdNU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dHdNU(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(1) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
          IF(IOS .EQ. 0)ALLOCATE (FEDD_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GEDD_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (H_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (N_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RSQN_ON_RSQJ_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KMID_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (JNU_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_PREV(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(2) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
          IF(IOS .EQ. 0)ALLOCATE (FEDD_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GEDD_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (H_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (N_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RSQN_ON_RSQJ_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KMID_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (JNU_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_SAVE(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(3) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
	  ALLOCATE (TA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TB(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAUONQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (Q(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SOURCE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VB(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (COH_VEC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HL(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HS(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HD(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (P_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (P_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (JOLD(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VdJdR_TERM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VdHdR_TERM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DELTA(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (DELTAH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (W(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (WPREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSI(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSIPREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EPS(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EPS_PREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJNU_PREV(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(3) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
	END IF
!
        FEDD_SAVE(:)=0.0_LDP
        GEDD_SAVE(:)=0.0_LDP
        H_ON_J_SAVE(:)=0.0_LDP
        N_ON_J_SAVE(:)=0.0_LDP
        RSQN_ON_RSQJ_SAVE(:)=0.0_LDP
        KMID_ON_J_SAVE(:)=0.0_LDP
        JNU_SAVE(:)=0.0_LDP
        GAM_RSQHNU_SAVE(:)=0.0_LDP
	HBC_SAVE=0.0_LDP
	NBC_SAVE=0.0_LDP
	IN_HBC_SAVE=0.0_LDP
!
	IF(ND_SAV .EQ. ND)THEN
	  RETURN
	ELSE
	  WRITE(LUER,*)'Error in ALLOC_MOD_J_REL'
	  WRITE(LUER,*)'At present allocation size cannot be changed'
	  STOP
	END IF
!
	RETURN
	END
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame for an EXPANDING atmosphere with a monotonic velocity law.
! The computed intensity thus depends on the intensity computed for the
! previous (bluer) frequency.
!
! FULLY RELATIVISTIC SOLUTION:
!          Ref: Mihalas, ApJ, 237, 574
!               Solution of the Comoving-Frame Equation of Transfer in
!               Spherically Symmetric Flows. VI Relativistic flows.
!               See notes for implementation.
!
! The solution involves the use of several moment ratios. These moment
! ratios should be computed using a formal ray by ray solution. Routine
! should be called in a loop to converge the "MOMENT" ratios (or
! EDDINGTON) factors.
!
! Required MOMENT ratios:
!
!                	     F = K / J    (d=1,2,..., N)
!	                H_ON_J = H/J      (d=1,2,...,N)
!	                N_ON_J = N/J      (d=1,2,...,N)
!
!	                     G =  N / H (d=1.5, 2.5, ..., N-0.5)
!	       RSQN_ON_RSQJ(I) = GAM RSQ_N(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!	          KMID_ON_J(I) = GAM RSQ_K(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!	RSQ_K=0.25*K(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!                           RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) RSQN_ON_RSQJ is defined at all
!                           depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or RSQN_ON_RSQJ is non-zero,
!                           and is the value to be used in MOM_J_CMF
!
	SUBROUTINE MOM_J_REL_V2(ETA,CHI,ESEC,V,SIGMA,R,
	1		   H_ON_J,FEDD,N_ON_J,GEDD,RSQN_ON_RSQJ,KMID_ON_J,
	1                  JNU,GAM_RSQHNU,dlnJdlnR,
	1                  HBC,IN_HBC,NBC,FREQ,dLOG_NU,
	1                  DIF,DBB,IC,METHOD,COHERENT,
	1                  INCL_ADVEC_TERMS,INCL_REL_TERMS,INIT,ND)
	USE SET_KIND_MODULE
	USE MOD_J_REL
	IMPLICIT NONE
!
! Created: 31-Dec-2004 : Based on MOM_J_REL_V1
!                        Removed all *PREV quantities from call and placed
!                          in module.
!
	INTEGER NC
	INTEGER NP
	INTEGER ND
!
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) V(ND)			!in km/s
	REAL(KIND=LDP) SIGMA(ND)		!dlnV/dlnR
	REAL(KIND=LDP) R(ND)			!in units of 10^10 cm
!
! Moment ratio variables. All must be supplied.
!
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) GEDD(ND)
	REAL(KIND=LDP) H_ON_J(ND)
	REAL(KIND=LDP) N_ON_J(ND)
	REAL(KIND=LDP) RSQN_ON_RSQJ(ND)
	REAL(KIND=LDP) KMID_ON_J(ND)
!
! These values are computed, and returned.
!
	REAL(KIND=LDP) JNU(ND)
	REAL(KIND=LDP) GAM_RSQHNU(ND)
	REAL(KIND=LDP) dlnJdlnR(ND)
!
! Boundary conditions: Must be supplied.
!
	REAL(KIND=LDP) HBC,NBC,IN_HBC
!
	REAL(KIND=LDP) DBB,IC
	REAL(KIND=LDP) FREQ,dLOG_NU
	CHARACTER*6 METHOD
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields). When INIT is true, we also calcualte quantities
! that will be used for other frequencies.
!
	LOGICAL INIT
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
	LOGICAL COHERENT
	LOGICAL DIF			!Use diffusion approximation?
!
! If INCL_REL_TERMS is true, terms of order v/c, which are usually
!         are included. These terms are usually excluded from stellar
!         wind models, but included in SN models.
!
	LOGICAL INCL_REL_TERMS
!
! If INCL_ADVEC_TERMS is true, the advection terms (dJ/dr and dH/dr)
!   are included. For WIND models these terms are usually neglected.
!   For supernova models, the ADVECTION terms are sometimes neglected.
!
	LOGICAL INCL_ADVEC_TERMS
!
! Local variables.
!
	REAL(KIND=LDP) DAMP_FAC
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) MAX_ER
	INTEGER COUNT
	INTEGER IFAIL
	INTEGER I,J
	INTEGER LUER,ERROR_LU
	LOGICAL ACCURATE
!
! 
!
	LUER=ERROR_LU()
	IF(INIT)THEN
!
! Allocate storage
!
	  CALL ALLOC_MOD_J_REL(ND)
!
! Zero all vectors.
!
	  DELTAH(:)=0.0_LDP
	  DELTA(:)=0.0_LDP
	  W(:)=0.0_LDP
	  WPREV(:)=0.0_LDP
	  PSI(:)=0.0_LDP
	  PSIPREV(:)=0.0_LDP
	  JNU_PREV(:)=0.0_LDP
	  GAM_RSQHNU_PREV(:)=0.0_LDP
	  EPS(:)=0.0_LDP
	  EPS_PREV(:)=0.0_LDP
!
	  H_ON_J_PREV(:)=0.0_LDP
	  FEDD_PREV(:)=0.0_LDP
	  N_ON_J_PREV(:)=0.0_LDP
	  GEDD_PREV(:)=0.0_LDP
	  RSQN_ON_RSQJ_PREV(:)=0.0_LDP
	  KMID_ON_J_PREV(:)=0.0_LDP
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	
	  BETA_FREQ(1:ND)=V(1:ND)/2.99794E+05_LDP
	  IF(INCL_ADVEC_TERMS .OR. INCL_REL_TERMS)THEN
	    BETA(1:ND)=BETA_FREQ(1:ND)
	  ELSE
	    BETA(1:ND)=0.0_LDP
	  END IF
!
	  GAM_REL_SQ(1:ND)=1.0_LDP/(1.0_LDP-BETA(1:ND)*BETA(1:ND))
	  GAM_REL(1:ND)=SQRT(GAM_REL_SQ(1:ND))
	  CON_DELTA(1:ND)=BETA_FREQ(1:ND)/R(1:ND)
	  CON_dKdNU(1:ND)=GAM_REL_SQ(1:ND)*(SIGMA(1:ND)+1.0_LDP)-1.0_LDP
	  CON_dHdNU(1:ND)=BETA(1:ND)*GAM_REL_SQ(1:ND)*(SIGMA(1:ND)+1.0_LDP)
!
	  DO I=1,ND-1
	    T1=0.5_LDP*(BETA(I)+BETA(I+1))
	    AV_SIGMA(I)=0.5_LDP*(SIGMA(I)+SIGMA(I+1))
	    CON_DELTAH(I)=2.0_LDP*(BETA_FREQ(I)+BETA_FREQ(I+1))/(R(I)+R(I+1))
	    CON_dKdNUH(I)=T1*(AV_SIGMA(I)+1.0_LDP)/(1.0_LDP-T1*T1)
	    CON_dNdNUH(I)=(AV_SIGMA(I)+1.0_LDP)/(1.0_LDP-T1*T1) -1.0_LDP
	  END DO
!
! These are set to zero to insure all velocity terms are neglected.
!
	  H_ON_J(1:ND)=0.0_LDP
	  N_ON_J(1:ND)=0.0_LDP
	  KMID_ON_J(1:ND)=0.0_LDP
          RSQN_ON_RSQJ(1:ND)=0.0_LDP
!
	END IF
!
! If new frequency, we need to update PREV vectors which refer to the previous
! frequency. We can't update these on exit, since we iterate for each freqyency
! in this routine.
!
	IF(FREQ .NE. FREQ_SAVE)THEN
          H_ON_J_PREV(:)=H_ON_J_SAVE(:)
          FEDD_PREV(:)=FEDD_SAVE(:)
          N_ON_J_PREV(:)=N_ON_J_SAVE(:)
          GEDD_PREV(:)=GEDD_SAVE(:)
          RSQN_ON_RSQJ_PREV(:)=RSQN_ON_RSQJ_SAVE(:)
          KMID_ON_J_PREV(:)=KMID_ON_J_SAVE(:)
          JNU_PREV(:)=JNU_SAVE(:)
          GAM_RSQHNU_PREV(:)=GAM_RSQHNU_SAVE(:)
	  HBC_PREV=HBC_SAVE
	  NBC_PREV=NBC_SAVE
	  IN_HBC_PREV=IN_HBC_SAVE
	END IF
!
! 
!
! Zero relevant vectors and matrices.
!
	JNU(:)=0.0_LDP
	GAM_RSQHNU(:)=0.0_LDP
!
	IF(INCL_ADVEC_TERMS)THEN
	  VdHdR_TERM(1:ND)=H_ON_J(1:ND)*BETA(1:ND)
	ELSE
	  VdHdR_TERM(1:ND)=0.0_LDP
	END IF
!
!*****************************************************************************
!
! CHI_H refers to the modified CHI term that multiplies H in the 1ST
! moment equation. We evaluate it at the grid points, and then use the
! averaging procedure used in the non-relativistic case.
!
! CHI_J refers the the modified CHI term that multiplies J in the 0th
! moment equation.
!
	IF(INCL_REL_TERMS)THEN
	  DO I=1,ND
	    CHI_H(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(
	1       1.0_LDP+2.0_LDP*GAM_REL_SQ(I)*(SIGMA(I)+1.0_LDP) )
	    P_H(I)=1.0_LDP
	  END DO
	  IF(INCL_ADVEC_TERMS)THEN
	    DO I=1,ND
	      CHI_J(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(1.0_LDP+SIGMA(I))
	    END DO
	  ELSE
	    DO I=1,ND
	      CHI_J(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(
	1                   2.0_LDP+GAM_REL_SQ(I)*(1.0_LDP+SIGMA(I)) )
	    END DO
	  END IF
	ELSE
	  DO I=1,ND
	    CHI_H(I)=CHI(I)
	    P_H(I)=1.0_LDP
	    CHI_J(I)=CHI(I)
	  END DO
	END IF
!
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0_LDP
	  END DO
	  MOM_ERR_CNT=0
	END IF
!
	SOURCE(1:ND)=ETA(1:ND)/CHI_J(1:ND)
	IF(COHERENT)THEN
	  COH_VEC(1:ND)=ESEC(1:ND)/CHI_J(1:ND)/GAM_REL(1:ND)
	ELSE
	  COH_VEC(1:ND)=0.0_LDP
	END IF
!
! NB: We actually solve for r^2 J, not J.
!
! Compute the Q factors from F. The Q is used to make the J and K terms
! appearing in the first moment equation a perfect differential.
! Note that we integrate from the core outwards, and normalize Q to the
! value of unity at the core.
!
	DO I=1,ND
	  TA(ND-I+1)=( 3.0_LDP*FEDD(I)-1.0_LDP+
	1           BETA(I)*N_ON_J(I)-(SIGMA(I)+1.0_LDP)*VdHdR_TERM(I)+
	1      GAM_REL_SQ(I)*BETA(I)*(SIGMA(I)+1.0_LDP)*
	1       (BETA(I)*(1.0_LDP-FEDD(I)-VdHdR_TERM(I))-N_ON_J(I))
	1            )/(FEDD(I)+VdHdR_TERM(I))/R(I)
	  TB(I)=R(ND-I+1)
	END DO
	CALL INTEGRATE(TB,TA,Q,IFAIL,ND)
!
	DO I=1,ND-1		!Q(ND) undefined after exiting INTEGRATE
	  TB(I)=EXP(Q(ND-I))
	END DO
!
! Scale by r^2
!
	DO I=1,ND-1
	  Q(I)=TB(I)/(R(I)/R(ND))**2
	END DO
	Q(ND)=1.0_LDP
!
! Compute optical depth scales.
!
	TA(1:ND)=CHI_H(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_H,TA,R,R,TB,ND)
!
	TA(1:ND)=CHI_J(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_J,TA,R,R,TB,ND)
!
	IF(.NOT. INIT)THEN
!
! We are integrating from blue to red. dLOG_NU is define as vd / dv which is
! the same as d / d ln v.
!
! EPS is used if we define N in terms of J rather than H, This is sometimes
! useful as H can approach zero, and hence N/H is undefined.
!
	  DO I=1,ND-1
	    DELTAH(I)=CON_DELTAH(I)/dLOG_NU/(CHI_H(I)+CHI_H(I+1))
	    W(I)=DELTAH(I)*(1.0_LDP+CON_dNdNUH(I)*GEDD(I))
	    WPREV(I)=DELTAH(I)*(1.0_LDP+CON_dNdNUH(I)*GEDD_PREV(I))
	    EPS(I)=DELTAH(I)*(CON_dNdNUH(I)*RSQN_ON_RSQJ(I)+
	1            CON_dKdNUH(I)*KMID_ON_J(I))/(P_H(I)+W(I))
	    EPS_PREV(I)=DELTAH(I)*(CON_dNdNUH(I)*RSQN_ON_RSQJ_PREV(I)+
	1            CON_dKdNUH(I)*KMID_ON_J_PREV(I))/(P_H(I)+W(I))
	  END DO
!
	  DO I=2,ND
	    DELTA(I)=CON_DELTA(I)/CHI_J(I)/dLOG_NU
	  END DO
	  DELTA(1)=CON_DELTA(1)/CHI_H(1)/dLOG_NU
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  PSI(1)=DELTA(1)*( HBC-NBC+(NBC+BETA(1)*FEDD(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0_LDP) )
	  PSIPREV(1)=DELTA(1)*( HBC_PREV-NBC_PREV+(NBC_PREV+
	1             BETA(1)*FEDD_PREV(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0_LDP) )
	END IF
!
	DO I=2,ND
	  DTAUONQ(I)=0.5_LDP*(DTAU_J(I)+DTAU_J(I-1))/Q(I)
	  PSI(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0_LDP+CON_dKdNU(I)*FEDD(I)+CON_dHdNU(I)*H_ON_J(I) )
	  PSIPREV(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0_LDP+CON_dKdNU(I)*FEDD_PREV(I)+
	1                       CON_dHdNU(I)*H_ON_J_PREV(I) )
	END DO
!
! NB: We are initially computing GAM_REL R^2 J. We need to multiply the
!     original JNU_PREV by GAM_REL R^2, since is was divided by
!     GAM_REL R^2 before it was stored.
!
	DO I=1,ND
	  GAM_RSQJNU_PREV(I)=GAM_REL(I)*R(I)*R(I)*JNU_PREV(I)
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=(FEDD(I+1)+VdHdR_TERM(I+1))*Q(I+1)/
	1                (P_H(I)+W(I))/DTAU_H(I)
	  HL(I)=(FEDD(I)+VdHdR_TERM(I))*Q(I)/
	1                (P_H(I)+W(I))/DTAU_H(I)
	  HS(I)=WPREV(I)/(P_H(I)+W(I))
	END DO
!
! 
!
! As we don't know dlnJdlnR and dHdlnR we need to iterate.
!
	ACCURATE=.FALSE.
	COUNT=0
	IF(INIT)dlnJdlnR(1:ND)=0.0_LDP
	JOLD(1:ND)=0.0_LDP
!
	DO WHILE(.NOT. ACCURATE)
!
	  COUNT=COUNT+1
!
	  IF(INCL_ADVEC_TERMS)THEN
	    VdJdR_TERM(1:ND)=CON_DELTA(1:ND)*dlnJdlnR(1:ND)/CHI_J(1:ND)
	    P_J(1:ND)=1.0_LDP+VdJdR_TERM(1:ND)
	  ELSE
	    VdJdR_TERM(1:ND)=0.0_LDP
	    P_J(1:ND)=1.0_LDP
	  END IF
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vectors either depend  on dlnJdlnR etc, or are corrupted in the
! solution.
!
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)-EPS(I-1)
	    TC(I)=-HU(I)+EPS(I)
	    TB(I)=DTAUONQ(I)*(P_J(I)-COH_VEC(I)) + PSI(I) + HL(I) +
	1             HU(I-1)-EPS(I-1)+EPS(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	  TC(1)=-( FEDD(2)+VdHdR_TERM(2) )*Q(2)/DTAU_H(1)
	  TB(1)= ( FEDD(1)+VdHdR_TERM(1) )*Q(1)/DTAU_H(1) + PSI(1) + HBC*P_H(1)
	  XM(1)=0.0_LDP
	  TA(1)=0.0_LDP
	  VB(1)=0.0_LDP
	  VC(1)=0.0_LDP
!
! Need to include relativistic terms.
!
	  TA(ND)=-Q(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	  IF(DIF)THEN
	    TB(ND)=(FEDD(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	    XM(ND)=GAM_REL(ND)*DBB*R(ND)*R(ND)/3.0_LDP/CHI(ND)
	  ELSE
	    TB(ND)=(FEDD(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)+IN_HBC*GAM_REL(ND)
	    XM(ND)=GAM_REL(ND)*R(ND)*R(ND)*IC*(0.25_LDP+0.5_LDP*IN_HBC)
	  END IF
	  TC(ND)=0.0_LDP
	  VB(ND)=0.0_LDP
	  VC(ND)=0.0_LDP
	  PSIPREV(ND)=0.0_LDP
!
	  XM(1)=XM(1) + PSIPREV(1)*GAM_RSQJNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*GAM_RSQHNU_PREV(I-1) + VC(I)*GAM_RSQHNU_PREV(I)
	1          + PSIPREV(I)*GAM_RSQJNU_PREV(I)
	1          - EPS_PREV(I-1)*(GAM_RSQJNU_PREV(I-1)+GAM_RSQJNU_PREV(I))
	1          + EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1))
	  END DO
	  XM(ND)=XM(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	  DO I=1,ND
	    IF(XM(I) .LT. 0)THEN
	      XM(I)=ABS(XM(I))/10.0_LDP
	      RECORDED_ERROR=.FALSE.
	      J=1
	      DO WHILE (J .LE. MOM_ERR_CNT .AND. .NOT. RECORDED_ERROR)
	        IF(MOM_ERR_ON_FREQ(J) .EQ. FREQ)RECORDED_ERROR=.TRUE.
	        J=J+1
	      END DO
	      IF(.NOT. RECORDED_ERROR .AND.
	1                     MOM_ERR_CNT .LT. N_ERR_MAX)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      END IF	
	    END IF
	  END DO
!
! Store J, correcting for the fact that we actually compute gamma r^2 J
!
	  DO I=1,ND
	    JNU(I)=XM(I)/R(I)/R(I)/GAM_REL(I)
	  END DO
!
	  DO I=1,ND-1
	    GAM_RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*GAM_RSQHNU_PREV(I) +
	1        ( EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1)) -
	1          EPS(I)*(XM(I)+XM(I+1)) )
	  END DO
!
	  IF(.NOT. INCL_ADVEC_TERMS)THEN
	     ACCURATE=.TRUE.
	  ELSE
!
! Compute new derivatives, and decide whether to iterate further.
!
	    MAX_ER=0.0_LDP
	    DO I=1,ND
	      MAX_ER=MAX(MAX_ER, ABS((JOLD(I)-JNU(I))/JNU(I)))
	    END DO
	    IF(MAX_ER .LT. 1.0E-08_LDP)ACCURATE=.TRUE.
	    JOLD(1:ND)=JNU(1:ND)
!
	    DAMP_FAC=0.8_LDP
	    IF(.NOT. ACCURATE .AND. INCL_ADVEC_TERMS)THEN
!	      CALL DERIVCHI(TB,JNU,R,ND,'LINMON')
!	      TB(1:ND)=R(1:ND)*TB(1:ND)/JNU(1:ND)
	      CALL DERIVCHI(TB,XM,R,ND,'LOGMON')
	      TB(1:ND)=R(1:ND)*TB(1:ND)/XM(1:ND)
	      IF(MAX_ER .LT. 0.01_LDP)THEN
	        dlnJdlnR(1:ND)=DAMP_FAC*TB(1:ND)+(1.0_LDP-DAMP_FAC)*dlnJdlnR(1:ND)
	      ELSE
	        dlnJdlnR(1:ND)=0.1_LDP*TB(1:ND)+0.9_LDP*dlnJdlnR(1:ND)
	      END IF
	    END IF
	    IF(COUNT .EQ.  100)THEN
	      WRITE(LUER,*)'Error in MOM_J_REL_V2: excessive iteration count.'
	      WRITE(LUER,*)'FREQ=',FREQ
	      WRITE(LUER,*)'Maximum error =',MAX_ER
	      ACCURATE=.TRUE.
	    END IF
	  END IF
!
	END DO
!
! Save variables for next frequency
!
	FREQ_SAVE=FREQ
        FEDD_SAVE(:)=FEDD(:)
        GEDD_SAVE(:)=GEDD(:)
        H_ON_J_SAVE(:)=H_ON_J(:)
        N_ON_J_SAVE(:)=N_ON_J(:)
        RSQN_ON_RSQJ_SAVE(:)=RSQN_ON_RSQJ(:)
        KMID_ON_J_SAVE(:)=KMID_ON_J(:)
        JNU_SAVE(:)=JNU(:)
        GAM_RSQHNU_SAVE(:)=GAM_RSQHNU(:)
	HBC_SAVE=HBC
	NBC_SAVE=NBC
	IN_HBC_SAVE=IN_HBC
!
	RETURN
	END
