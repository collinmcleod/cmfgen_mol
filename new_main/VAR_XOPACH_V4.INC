!
! Altered 21-Jun-2002:  CHanged to V4 and inserted LST_DEPTH_ONLY in call.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	    CALL VAR_X_OPA_ETA_V4(VCHI,VETA,
	1          ATM(ID)%XzV,ATM(ID)%XzVLTE,ATM(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,
	1          ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE,
	1          ATM(ID+1)%dlnXzVLTE_dlnT,
	1          ATM(ID+1)%NXzV,ED, ATM(ID+1)%DXzV,
	1          T,IMP_VAR,ATM(ID)%EQXzV, ATM(ID+2)%EQXzv,
	1          AT_NO(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1          CONT_FREQ,EMHNUKT_CONT,NT,ND,LST_DEPTH_ONLY)
	  END IF
	END DO
