!
! Subroutine to compute the opacity due to Rayleigh scattering at
! ND depth points. This routine is only accurate away from the
! Lyman lines. The maximum cross-section is limited to the 
! electron-scattering cross-section.
!
	SUBROUTINE RAYLEIGH_SCAT(RKI,HYD,AHYD,EDGE_HYD,NHYD,FREQ,ND)
	IMPLICIT NONE
!
! Created 17-Mar-2004
!
	INTEGER NHYD
	INTEGER ND
!
	REAL*8 AHYD(NHYD,NHYD)
	REAL*8 HYD(NHYD,ND)
	REAL*8 EDGE_HYD(NHYD)
	REAL*8 RKI(ND)
	REAL*8 FREQ
!
! Local variables:
!
	REAL*8 T1,T2
	REAL*8 GAM_ON_W0
	INTEGER I
!
! GAM_ON_W0= 2e^2 W0 / (3mc^3) where I have used w(1-2) for simplicity.
!	
	GAM_ON_W0=1.027D-07
	T1=0.0D0
	DO I=2,NHYD
	  T2=( FREQ/(EDGE_HYD(1)-EDGE_HYD(I)) )**2
	  T1=T1+AHYD(1,I)*T2/( (1.0D0-T2)**2 + GAM_ON_W0*T2 )
	END DO
	IF(T1 .GT. 1.0D0)T1=1.0D0
!
	T1=6.65D-15*T1
	RKI(1:ND)=RKI(1:ND)+T1*HYD(1,1:ND)
!
	RETURN
	END
