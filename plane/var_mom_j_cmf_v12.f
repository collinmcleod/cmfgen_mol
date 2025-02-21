!
! This subroutine computes the variation of J as a function of the emissivity
! ETA and the opacity CHI.
!
! TX has dimensions (ND,ND,NM) and TVX has dimension (ND-1,ND,NM).
!
! NM reflects the total number of matrices required to compute the variation
!    of J. All NM matrices are modified by a call to THOMAS in UP_TX_TVX. The
!    first 2 should represent dJ/dCHI and the second dJ/dETA --- the others
!    can be in any order.
!
! Typically they will be CHI_C, and ETA_C, CHIL_, AND ETAL of othe lines.
!
! NB TX(i,m, ) =dJ(i)/d(CHI(m),ETA(m),...)
! NB TVX(i,m, ) =dRSQH(i)/d(CHI(m),ETA(m),...)
!
! NM_KI reflects the 3rd dimension of KI. For this routine only the first 2
!   are important and are used to compute the variation of J with respect to
!   the CURRENT opacity and the CURRENT emissivity.
!
! TX_DIF_d_T and TX_DIFF_d_dTdR are used to describe the variations in J
! caused by the DIFFUSION approximation at the inner boundary.
!
! The particular choice of the outer boundary condition adopted is irrelevant
! for this routine. Such information is incorporated by the outer boundary
! Eddington factors HBC, and NBC.
!
! Note TA, TBC, TC, HU etc have been defined so that the program computes
! JNU and r^2 HNU.
!
	SUBROUTINE VAR_MOM_J_CMF_V12(ETA,CHI,ESEC,THETA,V,SIGMA,R,
	1               TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1               TVX_DIF_d_T,TVX_DIF_d_dTdR,
	1               KI,WORKMAT,RHS_dHdCHI,
	1               INIT,dLOG_NU,
	1               INNER_BND_METH,dTdR,DBB,dDBBdT,IC,IB_STAB_FACTOR,
	1               FREQ,H_CHK_OPTION,OUT_BC_TYPE,
	1               DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
	USE SET_KIND_MODULE
	USE MOD_RAY_MOM_STORE
	USE MOD_VAR_MOM_J_CMF
	IMPLICIT NONE
!
! Altered:   11-Jun-2109 : Added point source option.
! Altered:   17-Oct-2016 : Changed to V12: CHECK_H_ON J replaced by H_CHK_OPTION.
! Altered:   17-Feb-2016 : Changed to V11. DIF replaced by INNER_BND_METH and
!                            IB_STAB_FACTOR added.
! Altered:   24-Feb-2015 : Added CHECK_H_ON_J to call.
! Altered:   06-Aug-2007 : Installed integer OUT_BC_TYPE
!                            so that could have choice of boundary conditions
!                            This allows other options to be included.
! Altered:   04-Feb-1997 : THETA now passed. Passed to allow coherent
!                            scattering at some depths, and non-coherent
!                            scattering at other depths. Only introduced
!                            for the variation calculation. Logical variable
!                            COHERENT no longer required.
! Altered:   12-Dec-1996 : NM_KI installed. Changed to version V6.
! Altered:   05-Dec-1996 : PROGDESC set to REAL(KIND=LDP) value, PROG_ID installed
!                             ERROR_LU installed. TUNE installed.
! Altered:   25-Jan-1996 : HBC, NBC (and HBC_PREV, NBC_PREV) are now
!                            scalers. Other quantities were not needed
!                            with THK option processed by EXTENSION.
!                            Several lines deleted with HBC(2) etc.
!                          Changed to V5
!
! Altered:   11-Jan-1996 : Bug in calculation of DT_DIF_d_dTDR when using
!                             RSQN_ON_RSQJ. Initiliazation improved.

! Altered:   11-Mar-1995 : RSQN_ON_RSQJ installed so that N may be written in
!                          terms of H (using G) or J.
!                          Call modified, as were subroutines UP_TX_TVX and
!                          EDD_J_VAR.
!                          _V4 append to name.
!
! Finalized: 04-Nov-1994
! Created:   27-Sep-1995 : Diffusion approximation not tested yet.
!
	INTEGER ND
	INTEGER NM
	INTEGER NM_KI
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND),THETA(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND)
!
! Variation arrays and vectors.
!
	REAL(KIND=LDP) TX(ND,ND,NM),TVX(ND-1,ND,NM)
	REAL(KIND=LDP) KI(ND,ND,NM_KI)
	REAL(KIND=LDP) WORKMAT(ND,ND),RHS_dHdCHI(ND-1,ND)
	REAL(KIND=LDP) TX_DIF_d_T(ND),TX_DIF_d_dTdR(ND)
	REAL(KIND=LDP) TVX_DIF_d_T(ND),TVX_DIF_d_dTdR(ND)
!
	LOGICAL DO_THIS_TX_MATRIX(NM)
!
	REAL(KIND=LDP) dLOG_NU,dTdR,DBB,dDBBdT,IC
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) IB_STAB_FACTOR
	CHARACTER*6 METHOD
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
	LOGICAL INIT
	INTEGER OUT_BC_TYPE
	CHARACTER(LEN=*) INNER_BND_METH
	CHARACTER(LEN=*) H_CHK_OPTION
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) DTAU_BND
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I
	REAL(KIND=LDP) AV_SIGMA
	REAL(KIND=LDP) LOCAL_DBB
	LOGICAL DIF_OR_ZF
!
! 
!
! This call allocates the vectors, and initialzes vectors such as TA etc.
! It only allocates TA, if it is not already allocated.
!
	CALL MOD_VAR_MOM_ALLOC(ND)
!
	IF(INIT)THEN
	  CALL DP_ZERO(TX, ND*ND*NM )
	  CALL DP_ZERO(TVX, (ND-1)*ND*NM )
	  IF(ALLOCATED(TA))THEN
	    IF(ND .NE. VEC_LENGTH)THEN
	      I=ERROR_LU()
	      WRITE(I,*)'Problem in VAR_MOM_J_CMF_V9'
	      WRITE(I,*)'Incompatible vector size'
	      WRITE(I,*)'ND=',ND,'ND(vec_size)=',VEC_LENGTH
	      STOP
	    END IF
	  END IF
	  DO I=1,ND
	    JNUM1(I)=0.0_LDP
	    RSQ_HNUM1(I)=0.0_LDP
 	  END DO
	END IF
!
! Zero relevant vectors and matrices.
!
	DO I=1,ND
	  JNU(I)=0.0_LDP
	  RSQ_HNU(I)=0.0_LDP
	END DO
!
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
!
! Compute the Q factors from F. Then compute optical depth scale.
!
	CALL QFROMF(K_ON_J,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
!
! We need to call d_DERIVCHI_dCHI to set the TRAP derivatives.
!
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
	DO I=2,ND
	  RSQ_DTAUONQ(I)=0.5_LDP*R(I)*R(I)*(DTAU(I)+DTAU(I-1))/Q(I)
	END DO
!
! 
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
! NB - By definition, G is defined at the mesh midpoints.
!
! We evaluate and store the constant terms in the computation of GAMH and
! GAM, since number of operations only proportional to ND. Later on scaling is
! proportional to NM*ND*ND.
!
	IF(.NOT. INIT)THEN
	  DO I=1,ND-1
	    AV_SIGMA=0.5_LDP*(SIGMA(I)+SIGMA(I+1))
	    GAMH(I)=2.0_LDP*3.33564E-06_LDP*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /dLOG_NU/( CHI(I)+CHI(I+1) )
	    W(I)=GAMH(I)*( 1.0_LDP+AV_SIGMA*NMID_ON_HMID(I) )
	    WPREV(I)=GAMH(I)*( 1.0_LDP+AV_SIGMA*NMID_ON_HMID_PREV(I) )
	    EPS_A(I)=GAMH(I)*AV_SIGMA*NMID_ON_J(I)/(1.0_LDP+W(I))
	    EPS_B(I)=EPS_A(I)*R(I+1)*R(I+1)
	    EPS_A(I)=EPS_A(I)*R(I)*R(I)
	    EPS_PREV_A(I)=GAMH(I)*AV_SIGMA*NMID_ON_J_PREV(I)/(1.0_LDP+W(I))
	    EPS_PREV_B(I)=EPS_PREV_A(I)*R(I+1)*R(I+1)
	    EPS_PREV_A(I)=EPS_PREV_A(I)*R(I)*R(I)
	  END DO
	  DO I=1,ND
	    GAM(I)=3.33564E-06_LDP*V(I)/R(I)/CHI(I)/dLOG_NU
	  END DO
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  DO I=2,ND
	    PSI(I)=RSQ_DTAUONQ(I)*GAM(I)*( 1.0_LDP+SIGMA(I)*K_ON_J(I) )
	    PSIPREV(I)=RSQ_DTAUONQ(I)*GAM(I)*( 1.0_LDP+SIGMA(I)*K_ON_J_PREV(I) )
	  END DO
	END IF
!
! 
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=R(I+1)*R(I+1)*K_ON_J(I+1)*Q(I+1)/(1.0_LDP+W(I))/DTAU(I)
	  HL(I)=R(I)*R(I)*K_ON_J(I)*Q(I)/(1.0_LDP+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0_LDP+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS_A(I-1)
	  TC(I)=-HU(I)+EPS_B(I)
	  TB(I)=RSQ_DTAUONQ(I)*(1.0_LDP-THETA(I)) + PSI(I) +HU(I-1) +HL(I)
	1             -EPS_B(I-1)+EPS_A(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=RSQ_DTAUONQ(I)*SOURCE(I)
	END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	IF(OUT_BC_TYPE .LE. 1)THEN
	  PSI(1)=R(1)*R(1)*GAM(1)*( HBC+NBC*SIGMA(1) )
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	  TC(1)=-R(2)*R(2)*K_ON_J(2)*Q(2)/DTAU(1)
	  TB(1)=R(1)*R(1)*( K_ON_J(1)*Q(1)/DTAU(1) + HBC ) + PSI(1)
	  XM(1)=XM(1) + PSIPREV(1)*JNUM1(1)
	  TA(1)=0.0_LDP
	  VB(1)=0.0_LDP
	  VC(1)=0.0_LDP
	ELSE
	  PSI(1)=R(1)*R(1)*GAM(1)*(1.0_LDP+SIGMA(1)*K_ON_J(1))
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*(1.0_LDP+SIGMA(1)*K_ON_J_PREV(1))
	  T1=0.25_LDP*(CHI(2)+CHI(1))*(R(2)-R(1))
	  DTAU_BND=T1
	  TC(1)=(HU(1)-EPS_B(1))/T1
	  TB(1)=-(HL(1)+R(1)*R(1)*HBC+EPS_A(1))/T1-PSI(1)-R(1)*R(1)*(1.0_LDP-THETA(1))
	  XM(1)=-SOURCE(1)*R(1)*R(1)- HS(1)*RSQ_HNUM1(1)/T1- PSIPREV(1)*JNUM1(1)
	  XM(1)=XM(1)-(EPS_PREV_A(1)*JNUM1(1)+EPS_PREV_B(1)*JNUM1(2))/T1
	END IF
!
	TA(ND)=-R(ND-1)*R(ND-1)*K_ON_J(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0_LDP/CHI(ND)
	  LOCAL_DBB=DBB
	  DIF_OR_ZF=.TRUE.
	ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=HNU_AT_IB*R(ND)*R(ND)
	  LOCAL_DBB=HNU_AT_IB                     !Not used
	  DIF_OR_ZF=.FALSE.
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=IB_STAB_FACTOR*TB(ND)*R(ND)*R(ND)*(JPLUS_IB+JMIN_IB)
	  TB(ND)=(1.0_LDP+IB_STAB_FACTOR)*TB(ND)
	  LOCAL_DBB=0.0_LDP
	  DIF_OR_ZF=.FALSE.
	ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25_LDP+0.5_LDP*IN_HBC)
	  DIF_OR_ZF=.FALSE.
	END IF
	TC(ND)=0.0_LDP
	VB(ND)=0.0_LDP
	VC(ND)=0.0_LDP
	PSIPREV(ND)=0.0_LDP
!
! We create PSIPREV_MOD to save multiplications in the UP_TX_TVX routine/
! It is only different from PSIPREV when N_ON_J is non zero.
!
	PSIPREV_MOD(1)=PSIPREV(1)
	PSIPREV_MOD(ND)=PSIPREV(ND)
	DO I=2,ND-1
	  PSIPREV_MOD(I)=(EPS_PREV_A(I)-EPS_PREV_B(I-1)) + PSIPREV(I)
	END DO
!
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*RSQ_HNUM1(I-1) + VC(I)*RSQ_HNUM1(I)
	1          + PSIPREV_MOD(I)*JNUM1(I)
	1          + ( EPS_PREV_B(I)*JNUM1(I+1) - EPS_PREV_A(I-1)*JNUM1(I-1) )
	END DO
	XM(ND)=XM(ND) + PSIPREV_MOD(ND)*JNUM1(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
	DO I=1,ND
	  JNU(I)=XM(I)
	END DO
!
	DO I=1,ND-1
	  RSQ_HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQ_HNUM1(I) +
	1              (EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*XM(I)) +
	1              (EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*XM(I+1))
	END DO
!
! Make sure H satisfies the basic requirement that it is less than J.
!
	IF(H_CHK_OPTION .EQ. 'AV_VAL')THEN
	  DO I=1,ND-1
	    T1=(R(I)*R(I)*XM(I)+R(I+1)*R(I+1)*XM(I+1))/2.0_LDP
	    IF(RSQ_HNU(I) .GT. T1)THEN
	      RSQ_HNU(I)=0.99999_LDP*T1
	    ELSE IF(RSQ_HNU(I) .LT. -T1)THEN
	      RSQ_HNU(I)=-0.99999_LDP*T1
	    END IF
	  END DO
	ELSE IF(H_CHK_OPTION .EQ. 'MAX_VAL')THEN
	  DO I=1,ND-1
	    T1=MAX(R(I)*R(I)*XM(I),R(I+1)*R(I+1)*XM(I+1))
	    IF(RSQ_HNU(I) .GT. T1)THEN
	      RSQ_HNU(I)=0.99999_LDP*T1
	    ELSE IF(RSQ_HNU(I) .LT. -T1)THEN
	      RSQ_HNU(I)=-0.99999_LDP*T1
	    END IF
	  END DO
	END IF
!
! 
!
! The J and H components of the radiation field have been found. We can thus
! begin the variation computation.
!
! Compute d{non-radiation field}/dchi matrix.
!
	CALL TUNE(1,'MOM_EDD')
	CALL EDD_J_VAR_V6(KI,RHS_dHdCHI,WORKMAT,
	1           SOURCE,CHI,ESEC,THETA,DTAU,R,SIGMA,
	1           K_ON_J,Q,HU,HL,HS,RSQ_DTAUONQ,
	1           W,WPREV,PSI,PSIPREV,
	1           EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1           JNU,JNUM1,RSQ_HNUM1,
	1           LOCAL_DBB,DIF_OR_ZF,HBC,OUT_BC_TYPE,
	1           ND,NM_KI)
	CALL TUNE(2,'MOM_EDD')
!
! Evaluate the intensity variations.
!                                  TX=dJ/d(chi,eta,....)
!                          and     TVX=dRSQH/d(chi,eta,...)
!
! WORKMAT is dimension (ND,ND) is is used to temporarily save TX( , ,K) for
! each K.
!
	CALL TUNE(1,'UP_TX')
	CALL UP_TX_TVX_V2(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       WORKMAT,ND,NM,NM_KI,
	1                       DTAU_BND,OUT_BC_TYPE,
	1                       INIT,DO_THIS_TX_MATRIX)
	CALL TUNE(2,'UP_TX')
!
! 
!
! Evaluate diffusion variation. TXD is initially the value from the
! previous frequency.
!
	DO I=1,ND
	  TX_OLD_d_T(I)=TX_DIF_d_T(I)
	  TX_OLD_d_dTdR(I)=TX_DIF_d_dTDR(I)
	END DO
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  TX_DIF_d_T(1)=PSIPREV(1)*TX_DIF_d_T(1)
	  TX_DIF_d_dTdR(1)=PSIPREV(1)*TX_DIF_d_dTdR(1)
	  DO I=2,ND-1
	    TX_DIF_d_T(I)= PSIPREV_MOD(I)*TX_DIF_d_T(I)
	1             + VB(I)*TVX_DIF_d_T(I-1) + VC(I)*TVX_DIF_d_T(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_T(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_T(I-1) )
	    TX_DIF_d_dTdR(I)= PSIPREV_MOD(I)*TX_DIF_d_dTdR(I)
	1             + VB(I)*TVX_DIF_d_dTdR(I-1) + VC(I)*TVX_DIF_d_dTdR(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_dTdR(I-1) )
	  END DO
	  TX_DIF_d_T(ND)=dDBBdT*R(ND)*R(ND)/3.0_LDP/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=DBB/dTdR*R(ND)*R(ND)/3.0_LDP/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_dTdR(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
!
	  DO I=1,ND-1
	    TVX_DIF_d_T(I)=HU(I)*TX_DIF_d_T(I+1) - HL(I)*TX_DIF_d_T(I) +
	1        HS(I)*TVX_DIF_d_T(I) +
	1        (EPS_PREV_A(I)*TX_OLD_d_T(I)-EPS_A(I)*TX_DIF_d_T(I)) +
	1        (EPS_PREV_B(I)*TX_OLD_d_T(I+1)-EPS_B(I)*TX_DIF_d_T(I+1))
	    TVX_DIF_d_dTdR(I)=HU(I)*TX_DIF_d_dTdR(I+1) -
	1     HL(I)*TX_DIF_d_dTdR(I) +
	1     HS(I)*TVX_DIF_d_dTdR(I) +
	1     (EPS_PREV_A(I)*TX_OLD_d_dTdR(I)-EPS_A(I)*TX_DIF_d_dTdR(I)) +
	1     (EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)-EPS_B(I)*TX_DIF_d_dTdR(I+1))
	  END DO
!
	END IF	    	    !DIF END
!
! 
!
! Save JNU and RSQ_HNU for next frequency integration.
!
	DO I=1,ND
	  JNUM1(I)=JNU(I)
	  RSQ_HNUM1(I)=RSQ_HNU(I)
	END DO
!
	RETURN
	END
