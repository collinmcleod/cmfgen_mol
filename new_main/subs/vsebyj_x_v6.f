!
! Subroutine to increment the variation matrix due to terms which
! depend directly on the intensity J. The Radiative equilibrium equation
! is not altered. This routine is for X-ray ionization only, where
! 2 electrons are ejected.
!
! Routine also increments the ionization equilibrium equations.
!
	SUBROUTINE VSEBYJ_X_V6(ID,WSE_X,
	1             HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1             HN_B,HNST_B,EDGE_B,N_B,DI,ION_EQ_IN_BA,
	1             ED,T,JREC,dJRECdT,JPHOT,
	1             ND,NION,DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 15-May-2001 : Changed to V5; Uses JREC, JPHOT and dJRECdT
! Altered 12-Apr-2001 : Changed to use STEQ_DATA_MOD.
!                       Changed to V5.
! Altered 27-OCt-1995 : Changed to be compatible with super levels.
!                        Call chnaged --- Now _V3.
! Altered 06-Mar-1995 : Dimensions of WSE_X changed to (N,ND) from (N,NCF).
!                       _V2 append to call.
! Testing 22-Jul-1994 : Minor mods.
! Created 19-Jul-1993 : Based on VSEBYJ_COM and EVALSE_QWVJ.
!
	INTEGER*4 ID		!Ion identification
	INTEGER*4 N_A		!Number of levels in ionizations state i
	INTEGER*4 N_B		!Number of levels in ionizations state i+1
	INTEGER*4 ML		!Indicates current freqency in vector NU
	INTEGER*4 ND		!Number of depth points
        INTEGER*4 NION		!Total Number of ions in model.
        INTEGER*4 ION_EQ_IN_BA  !
!
	INTEGER*4 DST,DEND
!
	REAL*8 WSE_X(N_A,ND)		!Quadrature weights (incl. cross. sec.)
!
! _A refers to quantities associated with the atom WITH super levels.
!
	REAL*8 HN_A(N_A,ND)		!Pops. of ith ionzation stage
	REAL*8 HNST_A(N_A,ND)		!LTE   "    "  "      "       "
	REAL*8 dlnHNST_AdlnT(N_A,ND)	
!
! _B refers to quantities associated with the FULL atom of the next 
!  ionization stage.
!
	REAL*8 HN_B(N_B,ND)		!Pops. of (i+1)th ionization stage
	REAL*8 HNST_B(N_B,ND)		!LTE   "   "    "       "        "
	REAL*8 EDGE_B(N_B)
!
	REAL*8 DI(ND)			!Ion density for B levels.
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature in 10^4K.
	REAL*8 JREC(ND)			! (2hv^3/c/c+J).exp()*FQW/v
	REAL*8 dJRECdT(ND)		!
	REAL*8 JPHOT(ND)	        ! J*FQW/v	
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables
!
	INTEGER*4 NIV
	INTEGER*4 ION_EQ		!Ion eqation in SE(ID)%BA matrix.
	INTEGER*4 ION_V                 !Location of ion variable in SE(ID)%BA.
	INTEGER*4 I,J
	REAL*8 T3,T4
	REAL*8 RECIP_B_ION,BSTIM
	REAL*8 WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
!
	NIV=SE(ID)%N_IV
	ION_EQ=SE(ID)%XRAY_EQ
	ION_V=SE(ID)%LNK_TO_IV(ION_EQ_IN_BA)
!
	DO I=DST,DEND			!Which depth point.
	  RECIP_B_ION=HNST_B(1,I)/HN_B(1,I)
	  BSTIM=JREC(I)*RECIP_B_ION
	  DO J=1,N_A			!Which equation (for S.E. only)
	    IF(WSE_X(J,I) .NE. 0)THEN
	      WSE_BY_RJ=WSE_X(J,I)*JPHOT(I)
	      SE(ID)%BA_PAR(J,J,I)=SE(ID)%BA_PAR(J,J,I)-WSE_BY_RJ
!
! NB: In the following, the factor of 2.0 for ED_FAC arrises because 
! the rate is prop. to 
!
!        HNST_A(J,I)*HNST_B(1,I) .
!
! The variation with respect to HN_B(1,I) is intrinsically zero, since 
! HNST_A(J,I)/HN_B(1,I) is independent of HN_B(1,I)
!
	      T3=HNST_A(J,I)*WSE_X(J,I)*BSTIM
	      T4=HNST_A(J,I)*WSE_X(J,I)*RECIP_B_ION
	      DI_FAC=T3/DI(I)
	      ED_FAC=2.0D0*T3/ED(I)
	      T_FAC=T3*( dlnHNST_AdlnT(J,I) - 
	1             HDKT*(EDGE_B(1)+1.5D0)/T(I) )/T(I) + T4*dJRECdT(I)
!
	      SE(ID)%BA_PAR(J,ION_V,I)=SE(ID)%BA_PAR(J,ION_V,I) +DI_FAC
	      SE(ID)%BA_PAR(J,NIV-1,I)=SE(ID)%BA_PAR(J,NIV-1,I) +ED_FAC
	      SE(ID)%BA_PAR(J,NIV,I)  =SE(ID)%BA_PAR(J,NIV,I)   +T_FAC
!
! Include ionizations/recombinations explicitly in the rate equation
! of the target ion (eg He++(gs) for He+ ion/recoms). 
!
	      SE(ID)%BA_PAR(ION_EQ,J,I)    =SE(ID)%BA_PAR(ION_EQ,J,I)    +WSE_BY_RJ
	      SE(ID)%BA_PAR(ION_EQ,ION_V,I)=SE(ID)%BA_PAR(ION_EQ,ION_V,I)-DI_FAC
	      SE(ID)%BA_PAR(ION_EQ,NIV-1,I)=SE(ID)%BA_PAR(ION_EQ,NIV-1,I)-ED_FAC
	      SE(ID)%BA_PAR(ION_EQ,NIV,I)  =SE(ID)%BA_PAR(ION_EQ,NIV,I)  -T_FAC
!
	    END IF		!WSE(J,ML) .NE. 0
	  END DO
	END DO
!
	RETURN
	END
