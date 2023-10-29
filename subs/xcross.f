C
C Routine returns the K Shell photionization cross section for C,N and O
C ions, From Daltabuit and Cox, 1972, ApJ, 177, 855
C
C FREQ     : 10^15 Hz
C XCROSS   : Returned cross section in units of 10^-10 cm^2
C ZCORE    : Charge of nucleus.
C NUM_ELEC : Number of electrons in system that is being photionized.
C
	FUNCTION XCROSS(FREQ,ZCORE,NUM_ELEC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 17-Oct-2000 : ZCORE_SAVE introduced to limit number of error
C                          messages.
C Altered 02-Jul-1998 : No error message output for H and He
C                       LUER and ERROR_LU installed.
C
	REAL(KIND=LDP) XCROSS,FREQ,ZCORE,NUM_ELEC
C
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER J
	REAL(KIND=LDP) U,CON_FAC
C
	LOGICAL ZCORE_SAVE(92)
	SAVE ZCORE_SAVE
	DATA ZCORE_SAVE/92*.FALSE./
C
C
C
C Data variables. C, N, AND O
C
	REAL(KIND=LDP) SIG_C(6),SIG_N(7),SIG_O(8)
	REAL(KIND=LDP) EDGE_C(6),EDGE_N(7),EDGE_O(8)
	REAL(KIND=LDP) S_C(6),S_N(7),S_O(8)
	REAL(KIND=LDP) ALPHA_C(6),ALPHA_N(7),ALPHA_O(8)
C
	DATA SIG_C/0.194,0.526,0.850,0.930,0.997,1.06/
	DATA SIG_N/0.142,0.371,0.595,0.643,0.683,0.717,0.747/
	DATA SIG_O/0.109,0.275,0.439,0.470,0.496,0.518,0.537,0.554/
C
	DATA EDGE_C/490,392,347,317,296,280/
	DATA EDGE_N/666,552,496,459,432,412,395/
	DATA EDGE_O/870,739,672,627,595,570,550,533/
C
	DATA S_C/2.95,2.76,2.51,2.49,2.48,2.47/
	DATA S_N/2.95,2.79,2.57,2.55,2.54,2.54,2.53/
	DATA S_O/2.95,2.81,2.62,2.61,2.60,2.59,2.59,2.58/
C
	DATA ALPHA_C/1.287,1.325,1.0,1.0,1.0,1.0/
	DATA ALPHA_N/1.287,1.314,1.0,1.0,1.0,1.0,1.0/
	DATA ALPHA_O/1.287,1.308,1.0,1.0,1.0,1.0,1.0,1.0/
C
C 
C
C Cross sections for Ne, Mg, Si, and S
C
	REAL(KIND=LDP) SIG_NE(10),SIG_MG(12),SIG_SI(14),SIG_S(16)
	REAL(KIND=LDP) EDGE_NE(10),EDGE_MG(12),EDGE_SI(14),EDGE_S(16)
	REAL(KIND=LDP) S_NE(10),S_MG(12),S_SI(14),S_S(16)
	REAL(KIND=LDP) ALPHA_NE(10),ALPHA_MG(12),ALPHA_SI(14),ALPHA_S(16)
C
	DATA SIG_NE/0.075,0.18,0.267,0.282,0.295,0.305,0.314,
	1           0.322,0.329,0.336/
	DATA SIG_MG/0.055,0.13,0.179,0.188,0.194,0.200,0.205,0.210,
	1           0.213,0.217,0.220,0.223/
	DATA SIG_SI/0.044,0.10,0.128,0.134,0.138,0.141,0.144,0.147,
	1           0.149,0.151,0.153,0.155,0.157,0.158/
	DATA SIG_S/0.035,0.08,0.096,0.099,0.102,0.105,0.106,0.108,
	1           0.110,0.111,0.112,0.113,0.114,0.115,0.116,0.117/
C
	DATA EDGE_NE/1360,1195,1100,1050,1000,968,940,916,896,878/
	DATA EDGE_MG/1960,1761,1650,1570,1520,1470,1440,1410,1380,
	1            1360,1340,1320/
	DATA EDGE_SI/2670,2430,2300,2210,2140,2090,2050,2010,1980,
	1            1950,1930,1910,1880,1870/
	DATA EDGE_S/3480,3210,3060,2960,2880,2820,2770,2730,2690,
	1           2660,2630,2600,2580,2560,2540,2520/
C
	DATA S_NE/2.90,2.95,2.72,2.71,2.70,2.69,2.69,2.68,2.68,2.67/
	DATA S_MG/2.90,2.90,2.79,2.78,2.78,2.77,2.77,2.76,2.76,
	1         2.76,2.75,2.75/
	DATA S_SI/2.90,3.00,2.86,2.85,2.85,2.84,2.84,2.83,2.83,2.83,
	1         2.82,2.82,2.82,2.82/
	DATA S_S/2.90,3.00,2.92,2.91,2.91,2.90,2.90,2.89,2.89,2.89,
	1        2.89,2.88,2.88,2.88,2.88,2.88/
C
	DATA ALPHA_NE/1.25,1.28,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
	DATA ALPHA_MG/1.22,1.25,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	1             1.0,1.0/
	DATA ALPHA_SI/1.20,1.20,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	1             1.0,1.0,1.0,1.0/
	DATA ALPHA_S/1.19,1.20,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	1             1.0,1.0,1.0,1.0,1.0,1.0/
C
C 
C
	CON_FAC=1.0/0.24191			!10^15 Hz to ev
	XCROSS=0.0D0
	IF(ZCORE .EQ. 6.0)THEN			!Carbon
	  J=NUM_ELEC
	  U=EDGE_C(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_C(J)*(U**S_C(J))*(ALPHA_C(J)+ (1.0-ALPHA_C(J))*U )
	ELSE IF(ZCORE .EQ. 7.0)THEN 		!Nitrogen
	  J=NUM_ELEC
	  U=EDGE_N(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_N(J)*(U**S_N(J))*(ALPHA_N(J)+ (1.0-ALPHA_N(J))*U )
	ELSE IF(ZCORE .EQ. 8.0)THEN 		!Oxygen
	  J=NUM_ELEC
	  U=EDGE_O(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_O(J)*(U**S_O(J))*(ALPHA_O(J)+ (1.0-ALPHA_O(J))*U )
	ELSE IF(ZCORE .EQ. 10.0)THEN 		!Neon
	  J=NUM_ELEC
	  U=EDGE_NE(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_NE(J)*(U**S_NE(J))*(ALPHA_NE(J)+ (1.0-ALPHA_NE(J))*U )
	ELSE IF(ZCORE .EQ. 12.0)THEN 		!Magnesium
	  J=NUM_ELEC
	  U=EDGE_MG(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_MG(J)*(U**S_MG(J))*(ALPHA_MG(J)+ (1.0-ALPHA_MG(J))*U )
	ELSE IF(ZCORE .EQ. 14.0)THEN 		!Silicon
	  J=NUM_ELEC
	  U=EDGE_SI(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_SI(J)*(U**S_SI(J))*(ALPHA_SI(J)+ (1.0-ALPHA_SI(J))*U )
	ELSE IF(ZCORE .EQ. 16.0)THEN 		!Sulpher
	  J=NUM_ELEC
	  U=EDGE_S(J)/(CON_FAC*FREQ)
	  IF(U .GT. 1.0D0)RETURN
	  XCROSS=SIG_S(J)*(U**S_S(J))*(ALPHA_S(J)+ (1.0-ALPHA_S(J))*U )
	ELSE IF(ZCORE .EQ. 1.0 .OR. ZCORE .EQ. 2.0)THEN
C
C Do nothing for H and He
C
	ELSE
	  IF(.NOT. ZCORE_SAVE(NINT(ZCORE)))THEN
	    ZCORE_SAVE(NINT(ZCORE))=.TRUE.
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error : This ZCORE not yet implemented in XCROSS'
	    WRITE(LUER,*)'ZCORE=',ZCORE
	  END IF
	END IF
C
	XCROSS=XCROSS*1.0D-08
C
	RETURN
	END

