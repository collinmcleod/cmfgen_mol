	SUBROUTINE WRITE_LINE_LAMBDAS(OPTION)
	IMPLICIT NONE
	CHARACTER(LEN=*) OPTION
	LOGICAL DONE
!
	DONE=.FALSE.	
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'HI' .OR. OPTION .EQ. 'HYD')THEN
	  WRITE(6,*)'H I:     3971.204   4102.900   4341.692   4862.691   6564.60'
	  WRITE(6,*)'H I:    18756.1    12821.6    10941.1    40522.6    26257.7    21661.2'
	  DONE=.TRUE.	
	END IF
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'HEI' .OR. OPTION .EQ. 'HE')THEN
	  WRITE(6,*)'He I:    4027.35    4472.76    4923.30    5017.08    5049.15    5877.29'
	  WRITE(6,*)'He I:    6679.99    7067.20    7283.36   10833.1    20586.9'
	  DONE=.TRUE.	
	END IF
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'HE2' .OR. OPTION .EQ. 'HE')THEN
	  WRITE(6,*)'He II:   1640.42    4026.75    4201.02    4542.86    4687.01    5413.02   10125.36'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'C2' .OR. OPTION .EQ. 'CARB')THEN
	  WRITE(6,*)'C II:    1334.532   1335.663   1335.708   6579.87    6584.70'
	  DONE=.TRUE.	
	END IF
! 
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'CIII' .OR. OPTION .EQ. 'CARB')THEN
	  WRITE(6,*)'C III:   1174.933   1175.263   1174.933   1175.711   1175.987   1176.370'
          WRITE(6,*)'C III:   1247.383   1908.73    2297.58'
	  WRITE(6,*)'C III:   4648.719   4651.548   4652.775   M4650.11'
	  WRITE(6,*)'C III:   5697.50    8502.66' 
	  WRITE(6,*)'C III:   6729.34    6732.90    6744.01    6746.04    6746.23    M6742.5'
	  WRITE(6,*)'C III:   9703.76    9708.07    9709.10    9717.75    9720.42    M9713.3'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'CIV' .OR. OPTION .EQ. 'CARB')THEN
	  WRITE(6,*)'C IV:    1548.187   1550.177   4658.88     5471.21 5802.92    5813.58'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NIV' .OR. OPTION .EQ. 'NIT')THEN
          WRITE(6,*)'N IV:    1718.550   3479.72    3484.00    3485.93    4058.91'
	  DONE=.TRUE.	
	END IF
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NV' .OR. OPTION .EQ. 'NIT')THEN
	  WRITE(6,*)'N V:     1238.821   1242.804   4605.02    4621.27    4945.94    7620.56'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'OIV' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O IV:    1338.615   1342.990   1343.514   3063.426   3071.59'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'OV' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O V:     1218.34    1371.296'
	  WRITE(6,*)'O V:     2781.83    2787.81    2790.67    M2784.8'
	  WRITE(6,*)'O V:     3145.57    5115.482'
	  WRITE(6,*)'O V:     5573.36    5581.67    5584.78    5599.44    5605.83    M5599.46'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'OSIX' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O VI:    1031.912   1037.613'
	  WRITE(6,*)'O VI:    2071.01    3434.25    5292.12    7738.74'
	  WRITE(6,*)'O VI:    3812.43    3835.33    M3820.0'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NAI' .OR. OPTION .EQ. 'SOD')THEN
	  WRITE(6,*)'Na I:    5889.951   5895.924'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'PV' .OR. OPTION .EQ. 'PHOS')THEN
	  WRITE(6,*)'P V:     1117.977   1128.008'
	  DONE=.TRUE.	
	END IF
!  
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'SK2' .OR. OPTION .EQ. 'SIL')THEN
	  WRITE(6,*)'Si II:   1190.416   1193.290  1194.500    1197.394'
	  WRITE(6,*)'Si II:   1230.749   1231.658'
	  WRITE(6,*)'Si II:   1260.42    1260.738  1265.002    1304.360   1309.276'
	  WRITE(6,*)'Si II:   1526.72    1533.45    6348.85    6373.12'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'SKIV' .OR. OPTION .EQ. 'SIL')THEN
	  WRITE(6,*)'Si IV:   1128.325   1128.340'
	  WRITE(6,*)'Si IV:   1393.755   1402.770   4090.016   4117.264'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'SV' .OR. OPTION .EQ. 'SUL')THEN
	  WRITE(6,*)'S V:     1122.031   1128.666   1128.779   1133.901   1133.97    1501.763'    
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'CA2' .OR. OPTION .EQ. 'CAL')THEN
	  WRITE(6,*)'Ca II:   3934.777   3969.592   7293.48    7325.91    8500.35    8544.44    8664.52'
	  DONE=.TRUE.	
	END IF
!
	IF(.NOT. DONE)THEN
	  WRITE(6,*)'Species not recognized'
	  WRITE(6,*)'Current ionization stages are:'
	  WRITE(6,*)'HI    HeI    He2    C2   CII    CIV' 
	  WRITE(6,*)'NIV   OV'
	  WRITE(6,*)'OIV   OV     OSIX'
	  WRITE(6,*)'Sk2   SkIV'
	  WRITE(6,*)'SV    NaI    CaII'
	  WRITE(6,*)' '
	  WRITE(6,*)'HYD   HE  NIT  CARB  OXY  SIL   SUL   SOD  CAL'
	  WRITE(6,*)' '
	END IF
!
	RETURN
	END
