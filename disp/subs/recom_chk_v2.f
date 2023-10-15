 
c
	SUBROUTINE RECOM_CHK_V2(RECOM,EDGE_GS,GC2,GIONC2,NC2,
	1            EXC_EN,PHOT_ID,SUB_PHOT_C2,SPECIES_ID,T)
	IMPLICIT NONE
	INTEGER NC2
	INTEGER PHOT_ID
	INTEGER SPECIES_ID
C
	REAL(10) RECOM(NC2)
	REAL(10) EDGE_GS(NC2)
	REAL(10) GC2(NC2)
	REAL(10) GIONC2
	REAL(10) EXC_EN
	REAL(10) T
C
	EXTERNAL SUB_PHOT_C2
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(10) CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	REAL(10) LOW_CROSS(NC2)
	REAL(10) HIGH_CROSS(NC2)
	REAL(10) LTE_POP(NC2)
C
	REAL(10) T1
	REAL(10) RGU
	REAL(10) DEL_NU
	REAL(10) MAX_dNU
	REAL(10) LOW_NU,HIGH_NU
	REAL(10) FOUR_PI_ON_H
C
	INTEGER LEV
	INTEGER CHK_LEV
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
	FOUR_PI_ON_H=189.65D+15
	RECOM(:)=0.0D0
C
C Compute the normalized LTE population (for ED=1, DI=1, and no level
C dissolution.
C
	RGU=LOG(2.07078D-22)
	DO LEV=1,NC2
	  LTE_POP(LEV)=GC2(LEV)*EXP(HDKT*(EDGE_GS(LEV)+EXC_EN)/T+RGU) /
	1                 (T**1.5)/GIONC2
	END DO
C
	MAX_dNU=0.1D0*T/HDKT
	HIGH_NU=EDGE_GS(NC2)
	CALL SUB_PHOT_C2(SPECIES_ID,HIGH_CROSS,HIGH_NU,EDGE_GS,NC2,PHOT_ID,L_FALSE)
	T1=EXP(-HDKT*HIGH_NU/T) * HIGH_NU * HIGH_NU
	HIGH_CROSS(:)=HIGH_CROSS(:)*T1
	CHK_LEV=NC2-1
	DO WHILE(HIGH_NU .LT. 50.0)
C
	  LOW_NU=HIGH_NU
	  LOW_CROSS(:)=HIGH_CROSS(:)
	  HIGH_NU=MIN(LOW_NU*1.1,LOW_NU+MAX_dNU)
	  IF(CHK_LEV .GT. 0)THEN
	    IF(HIGH_NU .GT. EDGE_GS(CHK_LEV))THEN
	      HIGH_NU=EDGE_GS(CHK_LEV)
	      CHK_LEV=CHK_LEV-1
	    END IF
	  END IF
	  DEL_NU=HIGH_NU-LOW_NU
	  CALL SUB_PHOT_C2(SPECIES_ID,HIGH_CROSS,HIGH_NU,EDGE_GS,NC2,PHOT_ID,L_FALSE)
	  T1=EXP(-HDKT*HIGH_NU/T) * HIGH_NU * HIGH_NU
	  HIGH_CROSS(:)=HIGH_CROSS(:)*T1
C
C To avoid problems with weights near the bound-free edge, we being summing
C the recombination rates at the edge.
C
	  DO LEV=1,NC2
	   IF(LOW_NU .GE. EDGE_GS(LEV))THEN
	      RECOM(LEV)=RECOM(LEV)+DEL_NU*(LOW_CROSS(LEV)+HIGH_CROSS(LEV))
	   END IF
	  END DO
	END DO
	T1=FOUR_PI_ON_H*TWOHCSQ*0.5D0
	RECOM(:)=T1*LTE_POP(:)*RECOM(:)
C
	RETURN
	END
