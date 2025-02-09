	SUBROUTINE WRITE_LINE_LAMBDAS(OPTION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) OPTION
	LOGICAL DONE
!
! Altered: Fixed NaI wavelengths.

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
	  WRITE(6,*)'C II:    1063.28    1063.31   1141.63     1141.66   1141.74     1334.532   1335.663    1335.708'
	  WRITE(6,*)'C II:    4268.2 4    268.5    6579.87     6584.70   7233.33     7238.41    7239.16'
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
	  WRITE(6,*)'C IV:    1548.187   1550.177   4658.88    5471.21    5802.92    5813.58'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NIII' .OR. OPTION .EQ. 'NIT')THEN
          WRITE(6,*)'N III:    1182.97   1183.03    1184.51    1184.57    1316.44    1318.41'
          WRITE(6,*)'N III:    1746.82   1748.65    1751.22    1751.66    M1749.44'
	  WRITE(6,*)'N III:    4098.51   4104.55    M4100.53'
          WRITE(6,*)'N III:    4635.42   4641.94    4643.15    M4639.9'
	  DONE=.TRUE.	
	END IF
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NIV' .OR. OPTION .EQ. 'NIT')THEN
          WRITE(6,*)'N IV:    1485.50    1718.550   4058.91'
          WRITE(6,*)'N IV:    3479.72    3484.00    3485.93    M3481.83'
          WRITE(6,*)'N IV:    7113.24    7124.24    7129.21    7131.14     M7118.7'
	  DONE=.TRUE.	
	END IF
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'NV' .OR. OPTION .EQ. 'NIT')THEN
	  WRITE(6,*)'N V:     1238.821   1242.804   4605.02    4621.27    4945.94    7620.56'
	  DONE=.TRUE.	
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'OI' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O I:     1302.17    1306.03'
	  WRITE(6,*)'O I:     [6302.5]   [6365.5]'
	  WRITE(6,*)'O I:     7774.08    7776.31   7777.53                           2s2_2p3(4So)3s_5So -2s2_2p3(4So)3p_5Pe'
	  DONE=.TRUE.
	END IF
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'O2' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O II:     [3727.09]    [3729.86]                                     2s2_2p3_4So-2s2_2p3_2Do'
	  WRITE(6,*)'O II:     [7320.94]    [7322.01]    [7331.68]    [7332.75]           2s2_2p3_2Do-2s2_2p3_2Po'
	  DONE=.TRUE.
	END IF	
!
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'OIV' .OR. OPTION .EQ. 'OXY')THEN
	  WRITE(6,*)'O IV:    1338.615   1342.990   1343.514   '
	  WRITE(6,*)'O IV:    1397.23    1399.78    1401.16    1404.81'
	  WRITE(6,*)'O IV:    2916.31    2921.46    2926.18                         2s_2p(3Po)3d_2Do-2s_2p(3Po)3p_2Pe'
	  WRITE(6,*)'O IV:    3063.426   3071.59'
	  WRITE(6,*)'O IV:    3404.50    3412.67    3414.62    M3410.08             2s2_3d_2De-2s2_3p_2Po'
	  WRITE(6,*)'O IV:    3382.17    3386.49    3397.76    3410.64    M3389.72  2s_2p(3Po)3p_4De-2s_2p(3Po)3s_4Po'
	  WRITE(6,*)'O IV:    7032.34    7053.62'
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
	  WRITE(6,*)'Na I:    2853.649   2853.851   3303.320   3303.930'
	  WRITE(6,*)'Na I:    5891.583   5897.558'
	  WRITE(6,*)'Na I:   10749.38   10752.23   22062.42    22089.69'
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
	  WRITE(6,*)'Si II:   1526.72    1533.45   1808.013    1816.928   1817.451'
	  WRITE(6,*)'Si II:   4129.23    4132.05   4132.06                                3s2_3d_2De-3s2_4f_2Do'
	  WRITE(6,*)'Si II:   5042.44    5057.394  5057.726                               3s2_4p_2Po-3s2_4d_2De'
	  WRITE(6,*)'Si II:   6348.85    6373.12                                          3s2_4s_2Se-3s2_4p_2Po'
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
	IF(OPTION .EQ. ' ' .OR. OPTION .EQ. 'FE2' .OR. OPTION .EQ. 'IRON')THEN
	  WRITE(6,*)'[Fe II]:   4640.972   4665.753   4729.392   4799.617   4890.988  4959.611   5008.027'
	  WRITE(6,*)'Fe II:     4925.296   5109.835   5170.468'
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
	  WRITE(6,*)'Fe2   FeIII'
	  WRITE(6,*)' '
	  WRITE(6,*)'HYD   HE  NIT  CARB  OXY  SIL   SUL   SOD  CAL'
	  WRITE(6,*)' '
	END IF
!
	RETURN
	END
