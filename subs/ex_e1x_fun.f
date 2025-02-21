C
C Subroutine to compute the EXP(X)*E1(X) where E1(X) is the exponential
C integral and X > 0.
C
C The fits are from Abrmowitz and Stegun (p231).
C
	FUNCTION EX_E1X_FUN(X,IFAIL)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Created 14-Sep-1995.   Routine was verified against S13AAF using COMP_EX.
C
	REAL(KIND=LDP) X,EX_E1X_FUN
	INTEGER IFAIL
C
	REAL(KIND=LDP) W(6),A(4),B(4),NUM,DENOM
	INTEGER I
C
	DATA W/-0.57721566_LDP,0.99999193_LDP,-0.24991055_LDP,
	1       0.05519968_LDP,-0.00976004_LDP,0.00107857_LDP/
	DATA A/0.2677737343_LDP,8.6347608925_LDP,18.0590169730_LDP,8.5733287401_LDP/
	DATA B/3.9584969228_LDP,21.0996530827_LDP,25.6329561486_LDP,9.5733223454_LDP/
C
	IF(X .LE. 0.0_LDP)THEN
	 IFAIL=IFAIL+1
	ELSE IF(X .LE. 1)THEN
C
C Error in E1(x)+ln(x) < 2.0E-07
C
	  EX_E1X_FUN=W(6)
	  DO I=5,1,-1
	    EX_E1X_FUN=EX_E1X_FUN*X+W(I)
	  END DO
	  EX_E1X_FUN=EXP(X)*(EX_E1X_FUN-LOG(X))
	ELSE
C
C Error in X.EXP(X). E1(x) < 2.0D-08
C
	  NUM=1.0_LDP
	  DENOM=1.0_LDP
	  DO I=4,1,-1
	    NUM=NUM*X+A(I)
	    DENOM=DENOM*X+B(I)
	  END DO
	  EX_E1X_FUN=NUM/DENOM/X
	END IF
C
	RETURN
	END
