	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	    CALL VAR_X_OPA_ETA_V2(VCHI,VETA,
	1          ATM(ID)%XzV,ATM(ID)%XzVLTE,ATM(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,
	1          ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE,
	1          ATM(ID+1)%dlnXzVLTE_dlnT,
	1          ATM(ID+1)%NXzV,ED, ATM(ID+1)%DXzV,
	1          T,ATM(ID)%EQXzV, ATM(ID+2)%EQXzv,
	1          AT_NO(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1          CONT_FREQ,EMHNUKT_CONT,NT,ND)
	  END IF
	END DO
