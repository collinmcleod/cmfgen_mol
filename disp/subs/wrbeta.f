	SUBROUTINE WRBETA(CHIL,SIGMA,BETA,R,V,FREQ,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) CHIL(ND),BETA(ND),R(ND),V(ND),SIGMA(ND)
	REAL(KIND=LDP) FREQ
!
! Local variables
!
	INTEGER I,J,K
	REAL(KIND=LDP) Y(1001)
	REAL(KIND=LDP) VV,H,TAU,SIG1,X,S
	REAL(KIND=LDP) T1,T2,T3
!
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0
	REAL(KIND=LDP), PARAMETER :: RONE=1.0
	REAL(KIND=LDP), PARAMETER :: RTHREE=3.0
	REAL(KIND=LDP), PARAMETER :: RFOUR=4.0
!
	H=0.01_LDP
	DO 1000 I=1,ND
	  TAU=3.0E-10_LDP*CHIL(I)*R(I)/V(I)/FREQ
	  SIG1=RONE+SIGMA(I)
!
	  DO 10 K=1,101
	    X=RONE+SIGMA(I)*H*H*(K-1)*(K-1)
	    Y(K)=(RONE-EXP(-TAU/X))*X/TAU
10	  CONTINUE
!
	  S=SQRT(ABS(SIGMA(I)))
	  IF(SIGMA(I) .EQ. 0.0_LDP)THEN
			VV=RONE
	  ELSE IF(SIGMA(I) .GT. 0.0_LDP)THEN
			VV=ATAN(S)/S
	  ELSE
			VV=0.5_LDP*LOG((RONE+S)/(RONE-S))/S
	  END IF
!
	  IF(TAU .LT. 0.1_LDP .AND. TAU .LT. 0.1*SIG1)THEN
		T1=0.5_LDP*TAU
		T2=TAU*TAU/12.0_LDP
		T3=T2*TAU*0.1875_LDP
		BETA(I)=RONE-(T1+T3-T2)*VV+(T2-T3)/SIG1
	  ELSE IF(TAU .GT. 10.0_LDP .AND. TAU .GT. 10.0_LDP*SIG1)THEN
		BETA(I)=(RONE+SIGMA(I)/RTHREE)/TAU
	  ELSE
		BETA(I)=RZERO
		DO 200 J=1,99,2
		  BETA(I)=BETA(I)+Y(J)+RFOUR*Y(J+1)+Y(J+2)
200		CONTINUE
		BETA(I)=BETA(I)*H/RTHREE
	  END IF
!
1000	CONTINUE
!
	RETURN
	END
