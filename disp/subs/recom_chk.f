
c
	SUBROUTINE RECOM_CHK(RECOM,EDGE_GS,GC2,GIONC2,NC2,
	1            EXC_EN,PHOT_ID,SUB_PHOT_C2,T)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER NC2,PHOT_ID
C
	REAL(KIND=LDP) RECOM(NC2)
	REAL(KIND=LDP) EDGE_GS(NC2)
	REAL(KIND=LDP) GC2(NC2)
	REAL(KIND=LDP) GIONC2
	REAL(KIND=LDP) EXC_EN
	REAL(KIND=LDP) T
C
	EXTERNAL SUB_PHOT_C2
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	REAL(KIND=LDP) LOW_CROSS(NC2)
	REAL(KIND=LDP) HIGH_CROSS(NC2)
	REAL(KIND=LDP) LTE_POP(NC2)
C
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) RGU
	REAL(KIND=LDP) DEL_NU
	REAL(KIND=LDP) MAX_dNU
	REAL(KIND=LDP) LOW_NU,HIGH_NU
	REAL(KIND=LDP) FOUR_PI_ON_H
C
	INTEGER LEV
	INTEGER CHK_LEV
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
	FOUR_PI_ON_H=189.65E+15_LDP
	RECOM(:)=0.0_LDP
C
C Compute the normalized LTE population (for ED=1, DI=1, and no level
C dissolution.
C
	RGU=LOG(2.07078E-22_LDP)
	DO LEV=1,NC2
	  LTE_POP(LEV)=GC2(LEV)*EXP(HDKT*(EDGE_GS(LEV)+EXC_EN)/T+RGU) /
	1                 (T**1.5_LDP)/GIONC2
	END DO
C
	MAX_dNU=0.1_LDP*T/HDKT
	HIGH_NU=EDGE_GS(NC2)
	CALL SUB_PHOT_C2(HIGH_CROSS,HIGH_NU,EDGE_GS,NC2,PHOT_ID,L_FALSE)
	T1=EXP(-HDKT*HIGH_NU/T) * HIGH_NU * HIGH_NU
	HIGH_CROSS(:)=HIGH_CROSS(:)*T1
	CHK_LEV=NC2-1
	DO WHILE(HIGH_NU .LT. 50.0_LDP)
C
	  LOW_NU=HIGH_NU
	  LOW_CROSS(:)=HIGH_CROSS(:)
	  HIGH_NU=MIN(LOW_NU*1.1_LDP,LOW_NU+MAX_dNU)
	  IF(CHK_LEV .GT. 0)THEN
	    IF(HIGH_NU .GT. EDGE_GS(CHK_LEV))THEN
	      HIGH_NU=EDGE_GS(CHK_LEV)
	      CHK_LEV=CHK_LEV-1
	    END IF
	  END IF
	  DEL_NU=HIGH_NU-LOW_NU
	  CALL SUB_PHOT_C2(HIGH_CROSS,HIGH_NU,EDGE_GS,NC2,PHOT_ID,L_FALSE)
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
	T1=FOUR_PI_ON_H*TWOHCSQ*0.5_LDP
	RECOM(:)=T1*LTE_POP(:)*RECOM(:)
C
	RETURN
	END
