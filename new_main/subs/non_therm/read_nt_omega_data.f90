	SUBROUTINE READ_NT_OMEGA_DATA()
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	CHARACTER*80 FILE_NAME
	INTEGER ID
!
	INTEGER MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE
	PARAMETER (MAX_TRANS=30000)
	PARAMETER (MAX_TVALS=40)
	PARAMETER (MAX_TAB_SIZE=MAX_TVALS*MAX_TRANS)
	INTEGER ID_LOW(MAX_TRANS)
	INTEGER ID_UP(MAX_TRANS)
	INTEGER ID_INDX(MAX_TRANS)
	REAL*8 T_TABLE(MAX_TVALS)
	REAL*8 OMEGA_TABLE(MAX_TAB_SIZE)
	REAL*8 OMEGA_SET,OMEGA_SCALE
	INTEGER NUM_TRANS,NUM_TVALS
	INTEGER I,J,K
	INTEGER NF
!
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    FILE_NAME=TRIM(ION_ID(ID))//'_COL_DATA'
	    CALL GEN_OMEGA_RD_V2(OMEGA_TABLE,T_TABLE,ID_LOW,ID_UP,ID_INDX, &
	                         OMEGA_SCALE,OMEGA_SET,                    &
	                         ATM(ID)%XzVLEVNAME_F,ATM(ID)%GXzV_F,      &
	                         ATM(ID)%NXzV_F,FILE_NAME,                 &
	                         NUM_TRANS,NUM_TVALS,                      &
	                         MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE)
	    NF=ATM(ID)%NXzV_F
	    IF(.NOT. ALLOCATED(ATM(ID)%NT_OMEGA))ALLOCATE(ATM(ID)%NT_OMEGA(NF,NF))
	    ATM(ID)%NT_OMEGA=0.0D0
	    DO K=1,NUM_TRANS
	      I=ID_LOW(K)
	      J=ID_UP(K)
	      ATM(ID)%NT_OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K)+NUM_TVALS-1)
	    END DO
	  END IF
	END DO

	RETURN
	END
