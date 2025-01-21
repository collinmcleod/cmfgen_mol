!
! Subroutine to compute the quadrature weights for the statistical
! equilibrium equations. These quadrature weight now have to be multplied
! by FQW/NU before use. This change was made to allow for a fixed continuum
! photioization cross-section.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine also handles level dissolution.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE QUAD_MULTI_V10(WSE_S,dWSE_SdT,WCR,dWCRdT,
	1                       HNST_S,dlnHNST_S_dlnT,N_S,
	1                       HNST_F_ON_S,EDGE_F,N_F,
	1                       F_TO_S_MAPPING,NU_CONT,T,ND,
	1                       COMPUTE_BA,FIXED_T,LAST_ITERATION,
	1                       DESC,ZION,PHOT_ID,ID)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
	EXTERNAL SUB_PHOT_GEN
!
! Altered 01-Oct-2023 - Bug fix. T2 instead of T3 was being compared to PHOT_DIS_PARAMETER.
! Altered 26-Jul-2019 - Added dWCRdT to call.
!                          Changes made for electron energy balance addition.
!                          Incorporated into IBIS -- 17-Aug-2019
! Altered 23-Oct-2016 - Bug fix (collision ionization with dissolution) and introduced
!                          PHOT_DIS_PARAMETER.
! Altered 04-Oct-2016 - Changed to V9
!                       Only compute WCR(I_S,J) and dWSE_SdT on LAST iteration, or
!                          when T is variable.
! Altered 05-Apr-2011 - Changed to V8.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to faciliate
!                         modifications allowing lower temperaturs.
!                       Most of editing done early 2011
! Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
!                         dissolution can be switched off completely.
! Altered 05-Sep-1997 - Option to assume that continuum cross-sections have
!                         not altered since the last call. Should result
!                         in CPU time. The definitions of WSE etc have changed
!                         which means that the following routines also need
!                         changing:
!                                   EVALSE
!                         As call changed, now V5.
!
! Altered 20-Sep-1996 - Extensive changes to allow SUB_PHOT to be called.
!                       Changes designed to improve speed and vectorization.
!                       Level dissolution effects directy incorporated.
!                       (As extensive changes called _V4, 13-Dec-1996)
! Altered 28-May-1996 - Now call PHOT_GEN_BLEND_V2
! Altered 08-Jun-1995 - WCR installed.
! Created 15-May-1995 - Based on QUADGEN_V4
!
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL(KIND=LDP) WSE_S(N_S,ND)
	REAL(KIND=LDP) dWSE_SdT(N_S,ND)
	REAL(KIND=LDP) WCR(N_S,ND)
	REAL(KIND=LDP) dWCRdT(N_S,ND)
	REAL(KIND=LDP) HNST_S(N_S,ND)
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,ND)
!
	REAL(KIND=LDP) HNST_F_ON_S(N_F,ND)
	REAL(KIND=LDP) EDGE_F(N_F)			!In 10^15 Hz
	INTEGER F_TO_S_MAPPING(N_F)
	REAL(KIND=LDP) T(ND)
!
	REAL(KIND=LDP) NU_CONT
	REAL(KIND=LDP) ZION
	CHARACTER*(*) DESC
	INTEGER PHOT_ID
!
	REAL(KIND=LDP) YDIS(ND)		!Constant for computing level dissolution/
	REAL(KIND=LDP) XDIS(ND)		!Constant for computing level dissolution/
	REAL(KIND=LDP) DIS_CONST(N_F)	!Constant appearing in dissolution formula.
	REAL(KIND=LDP) ALPHA_VEC(N_F)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION
!
! Local Variables,
!
	INTEGER I_S,I_F,J
	REAL(KIND=LDP) T1,T2,T3,ZION_CUBED,NEFF,FOUR_PI_D_H
	LOGICAL DO_ALL
!
! NB: WSE_OLD=WSE*FQW/NU
!     dWSEdT_OLD=dWSEdT*FQW/NU
!     WCR_OLD=(NU*WSE+WCR)*FQW/NU
!
! The factor of DEX(-10) in FOUR_PI_D_H is due to the definition of the
! cross-section in SUB_GEN_PHOT
! which is DEX(10) times the photoionization cross section so that
! CHI*R is constant.
!
! Note FOUR_PI_D_H differs by 10^-15 from original constant in QUADGEN because
! FQW has C units of Hz, not 10^15 Hz.
!
	FOUR_PI_D_H=1.0_LDP/5.27296E-03_LDP                    !1.8965D+02		!4*PI/H*DEX(-10)*DEX(-15)
!
	WSE_S(:,:)=0.0_LDP
	WCR(:,:)=0.0_LDP
	dWSE_SdT(:,:)=0.0_LDP
	dWCRdT(:,:)=0.0_LDP
!
! Get photoionization cross-sections for all levels. The first call returns
! the threshold cross-section when NU < EDGE.
!
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_TRUE)
	ELSE
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_FALSE)
	END IF
!
! DIS_CONST is the constant K appearing in the expression for level dissolution.
! A negative value for DIS_CONST implies that the cross-section is zero.
!
	DIS_CONST(1:N_F)=-1.0_LDP
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  ZION_CUBED=ZION*ZION*ZION
	  DO I_F=1,N_F
	    IF(NU_CONT .LT. EDGE_F(I_F) .AND. ALPHA_VEC(I_F) .NE. 0 .AND. NU_CONT .GT. 0.8_LDP*EDGE_F(I_F))THEN
	      NEFF=SQRT(3.289395_LDP*ZION*ZION/(EDGE_F(I_F)-NU_CONT))
	      IF(NEFF .GT. 2*ZION)THEN
	        T1=MIN(1.0_LDP,16.0_LDP*NEFF/(1+NEFF)/(1+NEFF)/3.0_LDP)
	         DIS_CONST(I_F)=( T1*ZION_CUBED/(NEFF**4) )**1.5_LDP
	      END IF
	    END IF
	  END DO
	END IF
!
! Compute dissolution vectors that are independent of level.
!
	IF(MOD_DO_LEV_DIS)THEN
	  DO J=1,ND
	    YDIS(J)=1.091_LDP*(X_LEV_DIS(J)+4.0_LDP*(ZION-1)*A_LEV_DIS(J))*
	1               B_LEV_DIS(J)*B_LEV_DIS(J)
	    XDIS(J)=B_LEV_DIS(J)*X_LEV_DIS(J)
	  END DO
	END IF
!
	DO_ALL=.FALSE.
	IF(COMPUTE_BA)DO_ALL=.TRUE.
	IF(FIXED_T)DO_ALL=.FALSE.
	IF(LAST_ITERATION)DO_ALL=.TRUE.
!
! We have to loop over depth (rather than frequency) because of the
! FULL to SUPER level mapping.
!
! Note that WSE_S and WCR have alternate signs.
!
	IF(DO_ALL)THEN
	  DO I_F=1,N_F
	    I_S=F_TO_S_MAPPING(I_F)
	    IF(NU_CONT .GE. EDGE_F(I_F))THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        WSE_S(I_S,J)=WSE_S(I_S,J) + T1*HNST_F_ON_S(I_F,J)
	        WCR(I_S,J)=WCR(I_S,J) - EDGE_F(I_F)*T1*HNST_F_ON_S(I_F,J)
	        T2=T1*HNST_F_ON_S(I_F,J)*(dlnHNST_S_dlnT(I_S,J)+1.5_LDP+HDKT*EDGE_F(I_F)/T(J))/T(J)
	        dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) - T2
	        dWCRdT(I_S,J)=dWCRdT(I_S,J) + EDGE_F(I_F)*T2
	      END DO
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	    ELSE IF(DIS_CONST(I_F) .GE. 0.0_LDP)THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        T2=7.782_LDP+XDIS(J)*DIS_CONST(I_F)
	        T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	        IF(T3 .GT. PHOT_DIS_PARAMETER)THEN
	          T3=T1*T3
	          WSE_S(I_S,J)=WSE_S(I_S,J) + T3*HNST_F_ON_S(I_F,J)
	          WCR(I_S,J)=WCR(I_S,J) - EDGE_F(I_F)*T3*HNST_F_ON_S(I_F,J)
	          T2=T3*HNST_F_ON_S(I_F,J)*(dlnHNST_S_dlnT(I_S,J)+1.5_LDP+HDKT*EDGE_F(I_F)/T(J))/T(J)
	          dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) - T2
	          dWCRdT(I_S,J)=dWCRdT(I_S,J) + EDGE_F(I_F)*T2
	        END IF
	      END DO
	    END IF
	  END DO
	ELSE
	  DO I_F=1,N_F
	    I_S=F_TO_S_MAPPING(I_F)
	    IF(NU_CONT .GE. EDGE_F(I_F))THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        WSE_S(I_S,J)=WSE_S(I_S,J) + T1*HNST_F_ON_S(I_F,J)
	        WCR(I_S,J)=WCR(I_S,J) - EDGE_F(I_F)*T1*HNST_F_ON_S(I_F,J)
	      END DO
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	    ELSE IF(DIS_CONST(I_F) .GE. 0.0_LDP)THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        T2=7.782_LDP+XDIS(J)*DIS_CONST(I_F)
	        T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	        IF(T3 .GT. PHOT_DIS_PARAMETER)THEN
	          T3=T1*T3
	          WSE_S(I_S,J)=WSE_S(I_S,J) + T3*HNST_F_ON_S(I_F,J)
	          WCR(I_S,J)=WCR(I_S,J) - EDGE_F(I_F)*T3*HNST_F_ON_S(I_F,J)
	        END IF
	      END DO
	    END IF
	  END DO
	END IF
!
	RETURN
	END
