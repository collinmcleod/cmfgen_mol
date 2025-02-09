C
C Subroutine to compute the quadrature weights for a modified cubic rule.
C Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
C The routine makes specific assumptions concerning the behaviour of the
C mean intensity at the boundary, and would need to be moified to obtain
C other quantities other than the mean intensity. The program assumes that
C the first point corresponds to mu=1.0 .
C
C D1, and R2 are work vectors.
C
	SUBROUTINE JWEIGHT(X,W,N)
	USE SET_KIND_MODULE
C
C Altered 24-May-1996 - DP_ZERO deleted
C                       ERROR_LU etc installed.
C Altered 18-Feb-1987 - Quadratic extrapolation to zero if X(N) .NE. 0
C Altered 25-Feb-1987 - Normalization implemented.
C
	IMPLICIT NONE
	INTEGER N,I
	REAL(KIND=LDP) X(N),W(N),H,HN,RE,RF,SUM
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0_LDP
C
	H=X(1)-X(2)
	HN=X(2)-X(3)
	RF=X(1)-X(3)
	W(1)=0.5_LDP*H-H*H/12.0_LDP/RF-HN*H/RF/12.0_LDP
	W(2)=0.5_LDP*H+H*H/RF/6.0_LDP*(1.0_LDP+0.5_LDP*H/HN)+HN*H/RF/12.0_LDP
	W(3)=H*H/RF/12.0_LDP*(-1.0_LDP-H/HN)
C
	DO I=3,N-2
	  H=X(I-1)-X(I)
	  RF=X(I-2)-X(I)
	  RE=X(I-1)-X(I+1)
	  W(I-2)=W(I-2)-H*H/12.0_LDP/RF	
	  W(I-1)=W(I-1)+0.5_LDP*H+H*H/12.0_LDP/RE
	  W(I)=W(I)+0.5_LDP*H+H*H/12.0_LDP/RF	
	  W(I+1)=W(I+1)-H*H/12.0_LDP/RE
	END DO
C
	IF(X(N) .EQ. 0)THEN
C
C Assumes a quadratic variation for U. The is the integral from X(N-2) to zero.
C
	  RF=X(N-3)-X(N-1)
	  H=X(N-2)-X(N-1)
	  HN=X(N-1)-X(N)
	  W(N-3)=W(N-3)-H*H/RF/12.0_LDP
	  W(N-2)=W(N-2)+0.5_LDP*H
	  W(N-1)=W(N-1)+0.5_LDP*(H+2.0_LDP*HN/3.0_LDP)+H*H/HN/6.0_LDP+H*H/12.0_LDP/RF
	  W(N)=W(N)+2.0_LDP*HN/3.0_LDP-H*H/HN/6.0_LDP
	ELSE
C
C For outer boundary. Extrapolate to zero assuming a quadratic
C variation. i.e. U=A + B*mu*mu. Issue a warning indicating
C this. This integration is from X(N-1) to zero. The quadrature
C weight for X(N-1) will be negative. Due to the limit on the
C early DO loop, we also need to compute the integral form X(N-2) to X(N-1).
C
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - Angle extrapolation in JWEIGHT required'
C
C Inegral from X(N-2) to X(N-1)
C
	  H=X(N-2)-X(N-1)
	  RF=X(N-3)-X(N-1)
	  RE=X(N-2)-X(N)
	  W(N-3)=W(N-3)-H*H/12.0_LDP/RF	
	  W(N-2)=W(N-2)+0.5_LDP*H+H*H/12.0_LDP/RE
	  W(N-1)=W(N-1)+0.5_LDP*H+H*H/12.0_LDP/RF	
	  W(N)=W(N)-H*H/12.0_LDP/RE
C
C Integral from X(N-1) to zero.
C
	  RF=(X(N-1)-X(N))*(X(N-1)+X(N))
	  W(N-1)=W(N-1)-X(N-1)*(X(N)*X(N)-X(N-1)*X(N-1)/3.0_LDP)/RF
	  W(N)=W(N)+(X(N-1)**3)*2.0_LDP/RF/3.0_LDP
	END IF
C
C Ensure that the weights have the correct normalization (shouldnt be
C necessary.
C
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	SUM=1.0_LDP/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	   LUER=ERROR_LU()
	   WRITE(LUER,*)' Warning - weights were normalized in JWEIGHT'
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
C
	RETURN
	END
