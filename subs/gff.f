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
	IMPLICIT NONE
C
C Altered 21-Jun-1998 : D0's inserted.
C Altered 07-Dec-1991 : Ifs replaced by MIN and MAX functions.
C                       GFF_VEC subroutine included.
C
	REAL*8 RNU,T,ZION,GFF
C
	REAL*8 U,GB,DU,DGB
	INTEGER IU,IGB
C
	REAL*8, SAVE :: A(12,11)
	DATA A
	1/5.53,4.91,4.29,3.64,3.00,2.41,1.87,1.33,0.90,0.55,0.33,0.19
	1,5.49,4.87,4.25,3.61,2.98,2.41,1.89,1.39,0.95,0.59,0.36,0.21
	1,5.46,4.84,4.22,3.59,2.97,2.41,1.91,1.44,1.00,0.63,0.39,0.24
	1,5.43,4.80,4.18,3.56,2.95,2.41,1.93,1.49,1.08,0.72,0.46,0.28
	1,5.40,4.77,4.15,3.54,2.94,2.41,1.95,1.55,1.17,0.86,0.59,0.38
	1,5.25,4.63,4.02,3.41,2.81,2.32,1.90,1.56,1.30,1.01,0.76,0.53
	1,5.00,4.40,3.80,3.22,2.65,2.19,1.80,1.51,1.32,1.14,0.97,0.76
	1,4.69,4.13,3.57,2.97,2.44,2.02,1.68,1.42,1.30,1.18,1.09,0.96
	1,4.48,3.87,3.27,2.70,2.21,1.84,1.52,1.33,1.20,1.15,1.13,1.08
	1,4.16,3.52,2.98,2.45,2.01,1.67,1.41,1.25,1.15,1.11,1.10,1.09
	1,3.85,3.27,2.70,2.20,1.81,1.50,1.30,1.17,1.11,1.08,1.08,1.09/
C
	U=LOG10(4.7994D0*RNU/T)
	U=MAX(-3.99999999D0,U)
	U=MIN(1.499999D0,U)
	GB=LOG10(15.789D0*ZION*ZION/T)
	GB=MAX(-2.99999999D0,GB)
	GB=MIN(1.999999D0,GB)
	IU=INT(U*2.0+9.0)
	IGB=INT(GB*2.0+7.0)
	DU=2.0*U+9.0-IU
	DGB=2.0*GB+7.0-IGB
C
	GFF=((1.0-DU)*A(IU,IGB)+DU*A(IU+1,IGB))*(1.0-DGB)
	1 +((1.0-DU)*A(IU,IGB+1)+DU*A(IU+1,IGB+1))*DGB
C
	RETURN
	END
C
C 
C
	SUBROUTINE GFF_VEC(GFF_VAL,RNU,T,ZION,ND)
	IMPLICIT NONE
	INTEGER ND
	REAL*8 ZION  		!Ion charge
	REAL*8 RNU		!Frequency in units of 10^15 Hz
	REAL*8 T(ND)		!T in units 10^4 K
	REAL*8 GFF_VAL(ND)	!Returned free-free Gaunt factors.
C
C Important data
C
	REAL*8, SAVE :: A(12,11)
	DATA A
	1/5.53,4.91,4.29,3.64,3.00,2.41,1.87,1.33,0.90,0.55,0.33,0.19
	1,5.49,4.87,4.25,3.61,2.98,2.41,1.89,1.39,0.95,0.59,0.36,0.21
	1,5.46,4.84,4.22,3.59,2.97,2.41,1.91,1.44,1.00,0.63,0.39,0.24
	1,5.43,4.80,4.18,3.56,2.95,2.41,1.93,1.49,1.08,0.72,0.46,0.28
	1,5.40,4.77,4.15,3.54,2.94,2.41,1.95,1.55,1.17,0.86,0.59,0.38
	1,5.25,4.63,4.02,3.41,2.81,2.32,1.90,1.56,1.30,1.01,0.76,0.53
	1,5.00,4.40,3.80,3.22,2.65,2.19,1.80,1.51,1.32,1.14,0.97,0.76
	1,4.69,4.13,3.57,2.97,2.44,2.02,1.68,1.42,1.30,1.18,1.09,0.96
	1,4.48,3.87,3.27,2.70,2.21,1.84,1.52,1.33,1.20,1.15,1.13,1.08
	1,4.16,3.52,2.98,2.45,2.01,1.67,1.41,1.25,1.15,1.11,1.10,1.09
	1,3.85,3.27,2.70,2.20,1.81,1.50,1.30,1.17,1.11,1.08,1.08,1.09/
C
C
C Local Variables.
C
	INTEGER I,IU,IGB
	REAL*8 U,GB,DU,DGB,A1,A2,LOGT
C
	A1=LOG10(4.7994D0*RNU)
	A2=LOG10(15.789D0*ZION*ZION)
	DO I=1,ND
	  LOGT=LOG10(T(I))
	  U=A1-LOGT
	  U=MAX(-3.99999999D0,U)
	  U=MIN(1.499999D0,U)
	  GB=A2-LOGT
	  GB=MAX(-2.99999999D0,GB)
	  GB=MIN(1.999999D0,GB)
C
	  IU=INT(U*2.0+9.0)
	  IGB=INT(GB*2.0+7.0)
	  DU=2.0*U+9.0-IU
	  DGB=2.0*GB+7.0-IGB
	  GFF_VAL(I)=((1.0-DU)*A(IU,IGB)+DU*A(IU+1,IGB))*(1.0-DGB)
	1               +((1.0-DU)*A(IU,IGB+1)+DU*A(IU+1,IGB+1))*DGB
	END DO
C
	RETURN
	END
