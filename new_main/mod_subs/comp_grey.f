!
! Subroutine to compute the temperature structure for a grey atmosphere
! with a Rosseland mean opacity of ROSSMEAN. When passed, ROSSMEAN must include
! the effects (if any) of clumping.
!
! At present, two different options are available.
!
! (a) Spherical atmosphere with out velocity terms.
! (b) If JGREY_WITH_V_TERMS =.TRUE. we use a spherical atmosphere,
!       and include the effect of the first order Doppler shift.
!
	SUBROUTINE COMP_GREY(TGREY,TAU_ROSS,ROSSMEAN,LUER,NC,ND,NP)
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Created 21-Dec-2004
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER LUER
!
	REAL*8 TGREY(ND)
	REAL*8 TAU_ROSS(ND)
	REAL*8 ROSSMEAN(ND)
!
	REAL*8 RJ(ND)
	REAL*8 DTAU(ND)
	REAL*8 Z(ND)
	REAL*8 CHI(ND)
	REAL*8 dCHIdR(ND)
!
	REAL*8 TA(ND)
	REAL*8 TB(ND)
	REAL*8 TC(ND)
	REAL*8 Q(ND)
	REAL*8 GAM(ND)
	REAL*8 GAMH(ND)
	REAL*8 SOB(ND)
	REAL*8 XM(ND)
	REAL*8 FEDD(ND)
!
	REAL*8 T1,T2
	REAL*8 HBC_J
!
	INTEGER I,J,K,L
!
	LOGICAL LST_DEPTH_ONLY
	LOGICAL FIRST
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	CHI(1:ND)=ROSSMEAN(1:ND)
	IF(JGREY_WITH_V_TERMS)THEN
!
! This routine will supersede the one above, and included to zero
! order the effect of the velocity field.
!
	   T2=1.0D-05		!Accuracy to converge f
	   CALL JGREY_WITH_FVT(RJ,SOB,CHI,R,V,SIGMA,
	1           P,AQW,HMIDQW,KQW,NMIDQW,
	1           LUM,METHOD,DIF,IC,T2,ND,NC,NP)
	 ELSE
!
! Will use FEDD for F, GAM for NEWRJ, GAMH for NEWRK, and T2 for NEWHBC.
! Will use HBC_J for HBC. No need to modify JGREY, as outer boundary
! will always be optically thin.
!
	    DO I=1,ND
	      FEDD(I)=1.0D0/3.0D0
	    END DO
	    HBC_J=1.0D0
	    T1=1000.0D0
	    DO WHILE(T1 .GT. 1.0D-05)
	      CALL JGREY(TA,TB,TC,XM,DTAU,R,Z,P,RJ,
	1        GAM,GAMH,Q,FEDD,CHI,dCHIdR,
	1        AQW,KQW,LUM,HBC_J,T2,NC,ND,NP,METHOD)
	      T1=0.0D0
	      DO I=1,ND
	        T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	        FEDD(I)=GAMH(I)
	      END DO
	      T1=MAX(ABS(HBC_J-T2),T1)
	      HBC_J=T2
!	      WRITE(LUER,'('' Maximum change in Grey F is '',1P,E11.4)')T1
	    END DO
	  END IF
!
! Compute the temperature distribution, and the Rosseland optical depth scale.
! NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert T from units of 
! K to units of 10^4 K. The ' ' in TORSCL indicates TYPE of atmosphere,
! and here is set to ' ' so as TORSCL assumes a 1/r^2 density dependence
! at boundary.
! 
	DO I=1,ND
	  TGREY(I)=((3.14159265D0/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	END DO
!
! Compute the Rosseland optical depth scale. The ' ' in TORSCL indicates TYPE of atmosphere,
! and here is set to ' ' so that TORSCL assumes a power law dependence. If CHI(3) .LT. CHI(1),
! an 1/r^2 density dependence at the outer boundary is assumed.
!
        CALL TORSCL(TAU_ROSS,CHI,R,TB,TC,ND,METHOD,' ')
!
	RETURN
	END
