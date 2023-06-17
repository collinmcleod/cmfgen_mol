!
! Subroutine to compute the quadrature weights for a modified cubic rule.
! Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
! The routine makes specific assumptions concerning the behaviour of the
! mean intensity at the boundary, and would need to be moified to obtain
! other quantities other than 2nd moment of the intensity. The program assumes
! that the first point corresponds to mu=1.0 .
!
! The first derivatives are estimated by 
!                               dv(i)= (V(I-1)-V(I))/(MU(I-1)-MU(I))
!
	SUBROUTINE NWEIGHT_V2(X,dX,W,N)
	IMPLICIT NONE
!
! Altered 14-Jun-2023 - Added dX to call. Changed to V2
! Altered 25-May-1996 - CALL to DP_ZERO deleted
!                       ERROR_LU inserted.     
! Created 11-JAN-1988 - Based on KWEIGHT.
!
	INTEGER N
	REAL*8 X(N),dX(N),W(N)
!
	REAL*8 H,R,RE,RF,T1,SUM
	INTEGER I,ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	W(:)=0.0D0
!
! Firstly evaluate weights for integration from I=1 to I=2.
! R, RF and T1 are used in calculating the corection due
! the derivative of V at d=1.
!
	H=dX(1)
	R=dX(2)
!
	RF=( (X(1)**3)-0.3D0*X(1)*H*H ) * H*(R+2*H)/(R+H)/12.0D0
	T1=( (X(1)**3)-0.3D0*X(1)*H*H ) * (H**3)/R/(R+H)/12.0D0
	RE=H*H*( (X(2)**3)-0.3D0*H*H*X(2) )/(X(1)-X(3))/12.0D0
!
	W(1)=W(1)+0.5D0*H*(  (X(1)**3) -0.5D0*H*(X(1)**2)
	1         +(H**3)/60.0D0  )+RE-RF
	W(2)=W(2)+0.5D0*H*(  (X(2)**3) +0.5D0*H*(X(2)**2)
	1         -(H**3)/60.0D0  )+RF+T1
	W(3)=W(3)-RE-T1
!
	DO I=3,N-1
	  H=dX(I-1)
	  RF=dX(I-2)+dX(I-1)
	  RE=dX(I-1)+dX(I)
	  RF=H*H*(X(I-1)**3-0.3D0*H*H*X(I-1))/RF/12.0D0
	  RE=H*H*(X(I)**3-0.3D0*H*H*X(I))/RE/12.0D0
	  W(I-2)=W(I-2)-RF
	  W(I-1)=W(I-1)+0.5D0*H*(  (X(I-1)**3) -0.5D0*H*(X(I-1)**2)
	1         +(H**3)/60.0D0  )+RE
	  W(I)=W(I)+0.5D0*H*(  (X(I)**3) +0.5D0*H*(X(I)**2)
	1         -(H**3)/60.0D0  )+RF
	  W(I+1)=W(I+1)-RE
	END DO
!
! Assume that V(mu)=0 at mu=0, and that it depends linearly on mu.
! When X(N) .ne. 0 we assume that v'(n)=v(n)/mu*
! Subtracting the X to get dX is fine for high N.
!
	IF(X(N) .EQ. 0)THEN
	  W(N-1)=W(N-1)+0.2D0*(X(N-1)**4)
	ELSE
	  H=X(N-1)-X(N)
	  RF=H*H*(X(N-1)**3-0.3D0*H*H*X(N-1))/(X(N-2)-X(N))/12.0D0
	  RE=H*H*(X(N)**3-0.3D0*H*H*X(N))/X(N)/12.0D0
	  W(N-2)=W(N-2)-RF
	  W(N-1)=W(N-1)+0.5D0*H*(  (X(N-1)**3) -0.5D0*H*(X(N-1)**2)
	1         +(H**3)/60.0D0  )
	  W(N)=W(N)+0.5D0*H*(  (X(N)**3) +0.5D0*H*(X(N)**2)
	1         -(H**3)/60.0D0  )+RF+RE
!
! The above was for the integral from n-1 to n. Now need to extend the
! integral to MU=0.
!
	  W(N)=W(N)+0.2D0*(X(N)**4)
	END IF
!
! Ensure that the weights have the correct normalization (shouldnt be 
! necessary. We assume v=mu for the normalization.
!
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=0.2D0/SUM
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights required normalization in NWEIGHT'
	  WRITE(LUER,*)'SUM=',SUM
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
!
	RETURN
	END
