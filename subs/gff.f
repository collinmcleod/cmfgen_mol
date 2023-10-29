C
C Function to compute the free-free gaunt factors for a
C hydrogenic ion of charge Z. The gaunt facors are derived
C from a table in ATLAS (1970) which in turn was derived
C from the results of Karzas and Latter. The first 6 values
C for U=0.5 have been correctd by the addition of 0.1 . The
C frequency should be given in units of 1.0E+15 Hz and the
C temperature in units of 1.0E+04 K .
C
	FUNCTION GFF(RNU,T,ZION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 21-Jun-1998 : E0_DP's inserted.
C Altered 07-Dec-1991 : Ifs replaced by MIN and MAX functions.
C                       GFF_VEC subroutine included.
C
	REAL(KIND=LDP) RNU,T,ZION,GFF
C
	REAL(KIND=LDP) U,GB,DU,DGB
	INTEGER IU,IGB
C
	INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,3000)
	REAL(KIND=LDP), SAVE :: A(12,11)
	DATA A
	1/ 5.53E0_DP,4.91E0_DP,4.29E0_DP,3.64E0_DP,3.00E0_DP,2.41E0_DP,1.87E0_DP,1.33E0_DP,0.90E0_DP,0.55E0_DP,0.33E0_DP,0.19E0_DP
	1, 5.49E0_DP,4.87E0_DP,4.25E0_DP,3.61E0_DP,2.98E0_DP,2.41E0_DP,1.89E0_DP,1.39E0_DP,0.95E0_DP,0.59E0_DP,0.36E0_DP,0.21E0_DP
	1, 5.46E0_DP,4.84E0_DP,4.22E0_DP,3.59E0_DP,2.97E0_DP,2.41E0_DP,1.91E0_DP,1.44E0_DP,1.00E0_DP,0.63E0_DP,0.39E0_DP,0.24E0_DP
	1, 5.43E0_DP,4.80E0_DP,4.18E0_DP,3.56E0_DP,2.95E0_DP,2.41E0_DP,1.93E0_DP,1.49E0_DP,1.08E0_DP,0.72E0_DP,0.46E0_DP,0.28E0_DP
	1, 5.40E0_DP,4.77E0_DP,4.15E0_DP,3.54E0_DP,2.94E0_DP,2.41E0_DP,1.95E0_DP,1.55E0_DP,1.17E0_DP,0.86E0_DP,0.59E0_DP,0.38E0_DP
	1, 5.25E0_DP,4.63E0_DP,4.02E0_DP,3.41E0_DP,2.81E0_DP,2.32E0_DP,1.90E0_DP,1.56E0_DP,1.30E0_DP,1.01E0_DP,0.76E0_DP,0.53E0_DP
	1, 5.00E0_DP,4.40E0_DP,3.80E0_DP,3.22E0_DP,2.65E0_DP,2.19E0_DP,1.80E0_DP,1.51E0_DP,1.32E0_DP,1.14E0_DP,0.97E0_DP,0.76E0_DP
	1, 4.69E0_DP,4.13E0_DP,3.57E0_DP,2.97E0_DP,2.44E0_DP,2.02E0_DP,1.68E0_DP,1.42E0_DP,1.30E0_DP,1.18E0_DP,1.09E0_DP,0.96E0_DP
	1, 4.48E0_DP,3.87E0_DP,3.27E0_DP,2.70E0_DP,2.21E0_DP,1.84E0_DP,1.52E0_DP,1.33E0_DP,1.20E0_DP,1.15E0_DP,1.13E0_DP,1.08E0_DP
	1, 4.16E0_DP,3.52E0_DP,2.98E0_DP,2.45E0_DP,2.01E0_DP,1.67E0_DP,1.41E0_DP,1.25E0_DP,1.15E0_DP,1.11E0_DP,1.10E0_DP,1.09E0_DP
	1, 3.85E0_DP,3.27E0_DP,2.70E0_DP,2.20E0_DP,1.81E0_DP,1.50E0_DP,1.30E0_DP,1.17E0_DP,1.11E0_DP,1.08E0_DP,1.08E0_DP,1.09E0_DP/
C
	U=4.7994E0_DP; U=LOG10(U*RNU/T)
	U=MAX(-3.99999999E0_DP,U)
	U=MIN(1.499999E0_DP,U)
	GB=15.789E0_DP; GB=LOG10(GB*ZION*ZION/T)
	GB=MAX(-2.99999999E0_DP,GB)
	GB=MIN(1.999999E0_DP,GB)
	IU=INT(U*2+9.0E0_DP)
	IGB=INT(GB*2+7.0E0_DP)
	DU=2*U+9.0E0_DP-IU
	DGB=2*GB+7.0E0_DP-IGB
C
	GFF=((1.0E0_DP-DU)*A(IU,IGB)+DU*A(IU+1,IGB))*(1.0E0_DP-DGB)
	1 +((1.0E0_DP-DU)*A(IU,IGB+1)+DU*A(IU+1,IGB+1))*DGB
C
	RETURN
	END
C
C 
C
	SUBROUTINE GFF_VEC(GFF_VAL,RNU,T,ZION,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER ND
	REAL(KIND=LDP) ZION  		!Ion charge
	REAL(KIND=LDP) RNU		!Frequency in units of 10^15 Hz
	REAL(KIND=LDP) T(ND)		!T in units 10^4 K
	REAL(KIND=LDP) GFF_VAL(ND)	!Returned free-free Gaunt factors.
C
C Important data
C
	INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,3000)
	REAL(KIND=LDP), SAVE :: A(12,11)
	DATA A
	1/ 5.53E0_DP,4.91E0_DP,4.29E0_DP,3.64E0_DP,3.00E0_DP,2.41E0_DP,1.87E0_DP,1.33E0_DP,0.90E0_DP,0.55E0_DP,0.33E0_DP,0.19E0_DP
	1, 5.49E0_DP,4.87E0_DP,4.25E0_DP,3.61E0_DP,2.98E0_DP,2.41E0_DP,1.89E0_DP,1.39E0_DP,0.95E0_DP,0.59E0_DP,0.36E0_DP,0.21E0_DP
	1, 5.46E0_DP,4.84E0_DP,4.22E0_DP,3.59E0_DP,2.97E0_DP,2.41E0_DP,1.91E0_DP,1.44E0_DP,1.00E0_DP,0.63E0_DP,0.39E0_DP,0.24E0_DP
	1, 5.43E0_DP,4.80E0_DP,4.18E0_DP,3.56E0_DP,2.95E0_DP,2.41E0_DP,1.93E0_DP,1.49E0_DP,1.08E0_DP,0.72E0_DP,0.46E0_DP,0.28E0_DP
	1, 5.40E0_DP,4.77E0_DP,4.15E0_DP,3.54E0_DP,2.94E0_DP,2.41E0_DP,1.95E0_DP,1.55E0_DP,1.17E0_DP,0.86E0_DP,0.59E0_DP,0.38E0_DP
	1, 5.25E0_DP,4.63E0_DP,4.02E0_DP,3.41E0_DP,2.81E0_DP,2.32E0_DP,1.90E0_DP,1.56E0_DP,1.30E0_DP,1.01E0_DP,0.76E0_DP,0.53E0_DP
	1, 5.00E0_DP,4.40E0_DP,3.80E0_DP,3.22E0_DP,2.65E0_DP,2.19E0_DP,1.80E0_DP,1.51E0_DP,1.32E0_DP,1.14E0_DP,0.97E0_DP,0.76E0_DP
	1, 4.69E0_DP,4.13E0_DP,3.57E0_DP,2.97E0_DP,2.44E0_DP,2.02E0_DP,1.68E0_DP,1.42E0_DP,1.30E0_DP,1.18E0_DP,1.09E0_DP,0.96E0_DP
	1, 4.48E0_DP,3.87E0_DP,3.27E0_DP,2.70E0_DP,2.21E0_DP,1.84E0_DP,1.52E0_DP,1.33E0_DP,1.20E0_DP,1.15E0_DP,1.13E0_DP,1.08E0_DP
	1, 4.16E0_DP,3.52E0_DP,2.98E0_DP,2.45E0_DP,2.01E0_DP,1.67E0_DP,1.41E0_DP,1.25E0_DP,1.15E0_DP,1.11E0_DP,1.10E0_DP,1.09E0_DP
	1, 3.85E0_DP,3.27E0_DP,2.70E0_DP,2.20E0_DP,1.81E0_DP,1.50E0_DP,1.30E0_DP,1.17E0_DP,1.11E0_DP,1.08E0_DP,1.08E0_DP,1.09E0_DP/
C
C Local Variables.
C
	INTEGER I,IU,IGB
	REAL(KIND=LDP) U,GB,DU,DGB,A1,A2,LOGT
C
	A1=4.7994E0_DP; A1=LOG10(A1*RNU)
	A2=15.789E0_DP; A2=LOG10(A2*ZION*ZION)
	DO I=1,ND
	  LOGT=LOG10(T(I))
	  U=A1-LOGT
	  U=MAX(-3.99999999E0_DP,U)
	  U=MIN(1.499999E0_DP,U)
	  GB=A2-LOGT
	  GB=MAX(-2.99999999E0_DP,GB)
	  GB=MIN(1.999999E0_DP,GB)
	
C
	  IU=INT(U*2+9.0E0_DP)
	  IGB=INT(GB*2+7.0E0_DP)
	  DU=2*U+9.0E0_DP-IU
	  DGB=2*GB+7.0E0_DP-IGB
	  GFF_VAL(I)=((1.0E0_DP-DU)*A(IU,IGB)+DU*A(IU+1,IGB))*(1.0E0_DP-DGB)
	1               +((1.0E0_DP-DU)*A(IU,IGB+1)+DU*A(IU+1,IGB+1))*DGB
	END DO
C
	RETURN
	END
