C
C We have now store the revised popuations back in their individual
C storage locations. For some species we have 2 atomic models. For these
C species we need to take the super-level populations and compute:
C
C 1. The LTE population off all leve ls in the FULL atom.
C 2. The population off all levels in the FULL atom.
C 3. The LTE population off all super-levels.
C
C This is done by the routine SUP_TO_FULL.
C
C Compute vector constants for evaluating the level disolution. These
C constants are the same for all species. These are stored in a module,
C and are required by SUP_TO_FULL. We first need to revise POPION.
C
	  DO J=1,ND
	    POPION(J)=0.0D0
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	  END DO
C
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
C
	  DO ID=1,NUM_IONS-1
	    CALL SUP_TO_FULL_V3(
	1       ATM(ID)%XzV_F,       ATM(ID)%XzVLTE_F, ATM(ID)%W_XzV_F,
	1       ATM(ID)%EDGEXzV_F,   ATM(ID)%GXzV_F,   ATM(ID)%F_TO_S_XzV,
	1       ATM(ID)%INT_SEQ_XzV, ATM(ID)%NXzV_F,   ATM(ID)%DXzV_F,
	1       ATM(ID)%GIONXzV_F,   ATM(ID)%XzV,      ATM(ID)%XzVLTE,
	1       ATM(ID)%NXzV,        ATM(ID)%DXzV,     ATM(ID)%ZXzV,
	1       ATM(ID)%XzV_PRES,
	1       ATM(ID+1)%EDGEXzV_F,  ATM(ID+1)%GXzV_F,
	1       ATM(ID+1)%F_TO_S_XzV, ATM(ID+1)%NXzV_F,
	1       ATM(ID+1)%XzV_PRES,   T,ED,ND)
	  END DO
C
