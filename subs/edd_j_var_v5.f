C
C Routine to compute the opacity & emissivity variation matrices for
C the case with lines. Also computes the matrix dRHS_dCHI which
C multiply's %KI in the V equation. X is the line profile. It is assumed
C that :-
C				VK( , ,1)=dCHI
C				VK( , ,2)=dETA
C
C  				dRHS_dCHI( , ,)=dCHI
C
	SUBROUTINE EDD_J_VAR_V5(VK,RHS_dHdCHI,dTAUdCHI,
	1                  SOURCE,CHI,ESEC,ES_COH_VEC,DTAU,R,SIGMA,
	1                  MIDF,Q,HU,HL,HS,RSQ_DTAUONQ,
	1                  W,WPREV,PSI,PSIPREV,
	1                  EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                  JNU,JNUM1,RSQ_HNUM1,
	1                  DBB,DIF,ND,NM)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 12-Mar-2003 : Evaluation of terms rewritten for better numerical accuracy.
C Altered 05-Feb-1998 : ES_COH_VEC now passed. Logical variable coherent
C                         no longer required. Introduced to allow
C                         coherent scattering at some depths, and
C                         incoherent scattering at other depths.
C                         Equivalent to THETA in many programs.
C                         Changed to version 5.
C Altered 16-Aug-1996 : Loop splitting to improve vectorization.
C Altered 24-May-1996 : NV removed (Max valued for ND) using F90
C                       CALL to DP_ZERO replaced.
C Created 16-Sep-1994 : Based on EDDLINE_VAR
C
	INTEGER ND,NM
	REAL(KIND=LDP) VK(ND,ND,NM),RHS_dHdCHI(ND-1,ND)
	REAL(KIND=LDP) dTAUdCHI(ND,ND)
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),ESEC(ND),ES_COH_VEC(ND)
	REAL(KIND=LDP) DTAU(ND),R(ND),SIGMA(ND)
	REAL(KIND=LDP) MIDF(ND),Q(ND),HU(ND),HL(ND),HS(ND),RSQ_DTAUONQ(ND)
	REAL(KIND=LDP) W(ND),WPREV(ND),PSI(ND),PSIPREV(ND)
	REAL(KIND=LDP) EPS_A(ND),EPS_B(ND)
	REAL(KIND=LDP) EPS_PREV_A(ND),EPS_PREV_B(ND)
	REAL(KIND=LDP) JNU(ND),JNUM1(ND),RSQ_HNUM1(ND)
	REAL(KIND=LDP) DBB
	LOGICAL DIF
C
C Local variables.
C
	REAL(KIND=LDP) dHUdCHI(ND),dHLdCHI(ND),dHSdCHI(ND)
	REAL(KIND=LDP) dHUdTAU(ND),dHLdTAU(ND),EPS_FAC(ND)
	REAL(KIND=LDP) dRHSdI(ND),dRHSdJ(ND)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	INTEGER I,J,K,L
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) dUdCHI
	REAL(KIND=LDP) dTAdCHI_J,dTAdCHI_I
	REAL(KIND=LDP) dTCdCHI_I,dTCdCHI_K
	REAL(KIND=LDP) dTBdCHI_J,dTBdCHI_I,dTBdCHI_K,dTBdCHI
	REAL(KIND=LDP) dXM_EPS_J,dXM_EPS_I,dXM_EPS_K
C
C 
C
	IF(NM .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in EDD_J_VAR_V4 - NM_KI too small'
	  WRITE(I,*)'NM_KI='
	  STOP
	END IF
	VK(:,:,:)=0.0D0
	RHS_dHdCHI(:,:)=0.0D0
C
C Compute the dTAUdCHI matrix.
C
	CALL dSPHEREdCHI(dTAUdCHI,DTAU,R,Q,ND)
C
C The following derivatives are valid for all ML.
C
	DO I=1,ND-1
	  T1=(1.0D0+W(I))*(CHI(I)+CHI(I+1))
	  dHUdCHI(I)=HU(I)*W(I)/T1
	  dHUdTAU(I)=-HU(I)/DTAU(I)
	  dHLdCHI(I)=HL(I)*W(I)/T1
	  dHLdTAU(I)=-HL(I)/DTAU(I)
	  dHSdCHI(I)=-HS(I)/T1
	  EPS_FAC(I)=-1.0D0/T1
	END DO
C
C 
C
C Firstly we compute the variation of the elements with respect to
C DTAU. We then multiply by dTAUdCHI matrix.
C
C We have 3 separate loops over I to allow vectorization.
C To improve rounding error, we note that dTBdCHI needs to be added to
C both dTBdCHI_I and dTBdCHI_J.
C
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  dTAdCHI_J=-dHLdTAU(J)
	  dTCdCHI_I=-dHUdTAU(I)
	  T1=0.5D0*R(I)*R(I)/Q(I)
	  dTBdCHI_I=dHLdTAU(I)
	  dTBdCHI_J=dHUdTAU(J)
	  dTBdCHI=PSI(I)/(DTAU(J)+DTAU(I))+0.5D0*(1.0D0-ES_COH_VEC(I))*R(I)*R(I)/Q(I)
C
C dDELUB is use as correction because UB(I)=-TB(I)-PSI(I)-PSIPREV(I)
C
	  dUdCHI=PSIPREV(I)/(DTAU(J)+DTAU(I))
C
	  dRHSdJ(I)=  T1*SOURCE(I)
	1         - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1         - dTBdCHI*JNU(I)
	1         + dUdCHI*JNUM1(I)
C
	  dRHSdI(I)= T1*SOURCE(I)
	1         - (dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1         + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
C
	END DO
C
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  DO L=1,ND
	    VK(I,L,1)=VK(I,L,1)+dRHSdJ(I)*dTAUdCHI(J,L)
	    VK(I,L,1)=VK(I,L,1)+dRHSdI(I)*dTAUdCHI(I,L)
	  END DO
	END DO
C
C Can now update VK for direct opacity variation.
C
	DO I=2,ND-1
	  J=I-1
	  K=I+1
	  T1=0.5D0*R(I)*R(I)/Q(I)
C
	  dTAdCHI_J=-dHLdCHI(J)
	  dTAdCHI_I=-dHLdCHI(J)
	  dTCdCHI_I=-dHUdCHI(I)
	  dTCdCHI_K=-dHUdCHI(I)
C
	  dTBdCHI_J=dHUdCHI(J)
	  dTBdCHI_I=dHLdCHI(I)+dHUdCHI(J)
	  dTBdCHI= -PSI(I)/CHI(I)+RSQ_DTAUONQ(I)*ES_COH_VEC(I)/CHI(I)
	  dTBdCHI_K=dHLdCHI(I)
C
	  dUdCHI=-PSIPREV(I)/CHI(I)
C
	  dXM_EPS_J=( EPS_A(J)*JNU(J)-EPS_PREV_A(J)*JNUM1(J) +
	1             EPS_B(J)*JNU(I)-EPS_PREV_B(J)*JNUM1(I) )*EPS_FAC(J)
	  dXM_EPS_K=( EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*JNU(I) +
	1             EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*JNU(I+1) )*EPS_FAC(I)
	  dXM_EPS_I=dXM_EPS_J+dXM_EPS_K
C
C NB  :  VB(I)=-HS(J) and VC(I)=HS(I)
C
	  VK(I,J,1)=VK(I,J,1)
	1             - (dTAdCHI_J*JNU(J)+dTBdCHI_J*JNU(I))
	1             - dHSDCHI(J)*RSQ_HNUM1(J)
	1             + dXM_EPS_J
C
	  VK(I,K,1)=VK(I,K,1)
	1             - (dTCdCHI_K*JNU(K)+dTBdCHI_K*JNU(I))
	1             + dHSDCHI(I)*RSQ_HNUM1(I)
	1             + dXM_EPS_K
C
	  T1=T1*(DTAU(J)+DTAU(I))/CHI(I)
	  VK(I,I,1)=VK(I,I,1)
	1             - (dTAdCHI_I*JNU(J)+dTCdCHI_I*JNU(K)+dTBdCHI_I*JNU(I))
	1             + (dUdCHI*JNUM1(I)-dTBdCHI*JNU(I))
	1             + (dHSDCHI(I)*RSQ_HNUM1(I)-dHSDCHI(J)*RSQ_HNUM1(J))
	1             - T1*SOURCE(I)
	1             + dXM_EPS_I
C
	  VK(I,I,2)=T1
	END DO
C
C Now do the boundary conditions.
C
	T1=  ( MIDF(1)*Q(1)*JNU(1)*R(1)*R(1)
	1      - MIDF(2)*Q(2)*JNU(2)*R(2)*R(2) )/DTAU(1)/DTAU(1)
	DO L=1,ND
	  VK(1,L,1)=VK(1,L,1)+T1*dTAUdCHI(1,L)
	END DO
	VK(1,1,1)=VK(1,1,1) +
	1          ( PSI(1)*JNU(1)- PSIPREV(1)*JNUM1(1) )/CHI(1)
C
	IF(DIF)THEN
	  T1= ( R(ND)*R(ND)*MIDF(ND)*JNU(ND) -
	1          R(ND-1)*R(ND-1)*MIDF(ND-1)*Q(ND-1)*JNU(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	  VK(ND,ND,1)=VK(ND,ND,1)-
	1               DBB*R(ND)*R(ND)/3.0D0/CHI(ND)/CHI(ND)
	ELSE
	  T1= ( R(ND)*R(ND)*MIDF(ND)*JNU(ND) -
	1           R(ND-1)*R(ND-1)*MIDF(ND-1)*Q(ND-1)*JNU(ND-1) )
	1           / DTAU(ND-1)/DTAU(ND-1)
	  DO L=1,ND
	    VK(ND,L,1)=VK(ND,L,1)+T1*dTAUdCHI(ND-1,L)
	  END DO
	END IF
C
C 
C
C Now we compute the varaition of the equation which updates the
C flux variation. We don't correct for the line profile, since
C we would the require two matrices. Note that HU(I), HL(I) and
C HS(I) depend directly on CHI(I) and CHI(I+1).
C
	DO I=1,ND-1
	  T1=dHUdTAU(I)*JNU(I+1)-dHLdTAU(I)*JNU(I)
	  DO L=1,ND
	    RHS_dHdCHI(I,L)=RHS_dHdCHI(I,L)+T1*dTAUdCHI(I,L)
	  END DO
	  T1=dHUdCHI(I)*JNU(I+1) - dHLdCHI(I)*JNU(I)
	1                       + dHSdCHI(I)*RSQ_HNUM1(I) +
	1     EPS_FAC(I)*( EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*JNU(I) +
	1                  EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*JNU(I+1) )
	  RHS_dHdCHI(I,I)=RHS_dHdCHI(I,I) + T1
	  RHS_dHdCHI(I,I+1)=RHS_dHdCHI(I,I+1) + T1
	END DO
C
	RETURN
	END
